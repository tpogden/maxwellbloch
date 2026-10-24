"""Spatial propagation: ``lax.scan`` over z of (slice solve, source, field update).

Reproduces ``PropagationSolver.solve`` (Adams-Bashforth by default, or Euler):

* ``Omega`` is carried in angular units and is complex and drives the Hamiltonian
  in full (``Omega`` on ``rho_ab``'s element, ``conj(Omega)`` on the transpose), as
  in ``PropagationSolver``; only the deprecated ``real_drive=True`` uses ``Re(Omega)``.
* The first inner step is driven by the analytic input pulse (clerq builds the
  Hamiltonian from the ``rabi_freq_t_func`` there) and is an Euler step; later
  steps are driven by the cubic not-a-knot spline of the propagated field.
* The source is ``2i N(z_next) (2 pi g_k) sum_pairs factor * rho_ab``, averaged over
  velocity classes with clerq's Maxwell-Boltzmann weights; ``N`` is evaluated at
  the *next* z in both Adams-Bashforth terms.
* Velocity class ``v`` shifts the detuning of field ``k`` by
  ``thermal_delta_v / (2 pi) * factor_doppler_shift_k`` (cyclic).
* The very first inner step is the Euler bootstrap (from the coherence at ``z_min``);
  every later inner step is Adams-Bashforth, and the density of inner step ``i`` is
  ``N(z_min + (i + 1) h)``. ``PropagationSolver`` did this wrongly for
  ``z_steps_inner > 1`` (slices one inner step early, density one step ahead) until
  clerq #328; the two now agree for every ``z_steps_inner``.
"""

from __future__ import annotations

import math
import numbers
from dataclasses import dataclass

import diffrax as dfx
import jax
import jax.numpy as jnp
import numpy as np

from clerq.autodiff._config import require_float64
from clerq.autodiff._drive import CubicDrive
from clerq.autodiff._lindblad import solve_slice
from clerq.autodiff._operators import SystemOperators, build_operators
from clerq.autodiff._pulses import (
    jump_times,
    pulse_function,
    pulse_name,
    strip_index,
)

TWO_PI = 2.0 * math.pi


@dataclass(frozen=True)
class PropagationSetup:
    """Static description of a propagation problem, extracted from a solver.

    Attributes:
        ops: operators of the system.
        tlist: time grid, shape ``(n_t,)``.
        h: inner spatial step.
        z_steps: number of stored spatial steps.
        z_steps_inner: finite-difference steps per stored step.
        N_z: number density at ``z_min + (i + 1) h``, shape ``(z_steps * inner,)``.
        g: angular interaction strength ``2 pi g`` of each field.
        shifts: cyclic Doppler detuning ``thermal_delta / (2 pi)`` of each class.
        weights: normalised Maxwell-Boltzmann weight of each class, times the
            detuning interval it represents.
        pulses: per field ``(name, args, rabi_freq)`` of the analytic input pulse.
        real_drive: drive the Hamiltonian with ``Re(Omega)`` only (the deprecated
            behaviour of ``PropagationSolver(real_drive=True)``).
    """

    ops: SystemOperators
    tlist: np.ndarray
    h: float
    z_steps: int
    z_steps_inner: int
    N_z: np.ndarray
    g: np.ndarray
    shifts: np.ndarray
    weights: np.ndarray
    pulses: tuple
    real_drive: bool = False


def setup_from_solver(solver) -> PropagationSetup:
    """Extract the static problem data of a :class:`~clerq.PropagationSolver`.

    Raises:
        NotImplementedError: if a field is counter-propagating or its pulse
            function is not supported.
        ValueError: if the system has no fields, or the number of interaction
            strengths does not match the number of fields.
    """
    solver.check()
    if not solver.system.fields:
        raise ValueError(
            "to_autodiff needs a system with at least one coupled field to "
            "propagate; this system has none."
        )
    ops = build_operators(solver.system)
    h = solver.z_step_inner()
    n_steps = solver.z_steps * solver.z_steps_inner
    N_z = np.array(
        [
            solver.num_density_z_func(
                solver.z_min + (i + 1) * h, solver.num_density_z_args
            )
            for i in range(n_steps)
        ],
        dtype=float,
    )
    weights = np.asarray(solver.thermal_quadrature_weights, dtype=float)
    pulses = tuple(
        (
            pulse_name(f.rabi_freq_t_func),
            strip_index(f.rabi_freq_t_args),
            float(f.rabi_freq),
        )
        for f in solver.system.fields
    )
    # fail early, with a clear message, for pulses that cannot be traced
    for name, args, _ in pulses:
        pulse_function(name, args)
    return PropagationSetup(
        ops=ops,
        tlist=np.asarray(solver.tlist, dtype=float),
        h=float(h),
        z_steps=int(solver.z_steps),
        z_steps_inner=int(solver.z_steps_inner),
        N_z=N_z,
        g=np.asarray(solver.g, dtype=float),
        shifts=np.asarray(solver.thermal_delta_list, dtype=float) / TWO_PI,
        weights=weights / weights.sum(),
        pulses=pulses,
        real_drive=bool(getattr(solver, "real_drive", False)),
    )


def validate_options(step, strategy, checkpoints) -> None:
    """Check the scheme options; called when the problem is built and when it runs.

    Raises:
        ValueError: for an unknown ``step`` or ``strategy``, ``checkpoints``
            missing or not a positive integer for ``"s4"``, or given for any
            other strategy.
    """
    if step not in ("ab", "euler"):
        raise ValueError(f"Unknown step {step!r}; use 'ab' or 'euler'.")
    if strategy not in ("s2", "s3", "s4"):
        raise ValueError(f"Unknown strategy {strategy!r}; use 's2', 's3' or 's4'.")
    if strategy == "s4":
        ok = (
            isinstance(checkpoints, numbers.Integral)
            and not isinstance(checkpoints, bool)
            and checkpoints >= 1
        )
        if not ok:
            raise ValueError(
                "strategy='s4' needs `checkpoints`, a positive integer; "
                f"got {checkpoints!r}."
            )
    elif checkpoints is not None:
        raise ValueError(
            f"`checkpoints` only applies to strategy='s4', not {strategy!r}."
        )


def propagate_fields(
    setup: PropagationSetup,
    *,
    omega0=None,
    pulse_args=None,
    rabi_freqs=None,
    g=None,
    decay_scale=1.0,
    decay_rates=None,
    detunings=None,
    b_field=None,
    step: str = "ab",
    strategy: str = "s3",
    checkpoints: int | None = None,
    rtol: float = 1e-8,
    atol: float = 1e-10,
    n_steps: int | None = None,
    max_steps: int = 100_000,
    solver=None,
    adjoint=None,
):
    """Propagate the fields through the medium.

    Args:
        setup: the problem, from :func:`setup_from_solver`.
        omega0: input field samples at ``z_min``, shape ``(n_fields, n_t)``, in
            *angular* units. If given, the first slice is driven by their spline;
            default: the analytic input pulses of the fields.
        pulse_args: per field, the (suffix-free) arguments of its analytic input
            pulse (``ampl``, ``centre``, ``fwhm``, ...); default: the problem's.
            Values may be traced, so pulse parameters are differentiable. Unused
            when ``omega0`` is given.
        rabi_freqs: cyclic ``rabi_freq`` multiplying each analytic pulse, shape
            ``(n_fields,)``; default: the problem's.
        g: cyclic interaction strengths, shape ``(n_fields,)``; default: the
            problem's.
        decay_scale: multiplies every decay rate.
        decay_rates: cyclic ``rate`` of each entry of ``system.decays``, shape
            ``(n_decays,)``; default: the system's.
        detunings: cyclic fixed detunings of the fields, shape ``(n_fields,)``;
            default: the problem's.
        b_field: cyclic DC B-field in natural units; default: the problem's.
        step: ``"ab"`` (Adams-Bashforth, default) or ``"euler"``.
        strategy: ``"s2"`` stores every slice's residuals; ``"s3"`` (default)
            checkpoints the scan body (recomputes a slice on the backward pass);
            ``"s4"`` additionally limits the solver's inner checkpoints to
            ``checkpoints``.
        checkpoints: inner checkpoint count, required (a positive integer) for
            ``strategy="s4"`` and not allowed with the other strategies.
        rtol, atol: adaptive step-size tolerances.
        max_steps: the most steps one slice solve may take; the solve errors if
            it needs more. Reverse mode also sizes its default checkpoint
            count from it (about its base-2 logarithm), so a smaller value
            lowers memory use for large problems.
        n_steps: fixed-step solver with this many steps per slice (reference
            gradients); must not be commensurate with the time grid.
        solver: a diffrax solver, default ``Tsit5``.
        adjoint: a diffrax adjoint for the slice solves, overriding ``strategy``'s
            (for example ``ForwardMode()`` to use ``jax.jacfwd``).

    Returns:
        ``Omega(z, t)`` of shape ``(n_fields, z_steps + 1, n_t)``, complex,
        angular units, the engine's internal convention. ``AutodiffPropagation``
        divides by ``2 * pi`` to give ``PropagationSolver.rabi_freqs_zt``.

    Raises:
        RuntimeError: if JAX is not in 64-bit mode.
        ValueError: for an unknown ``step`` or ``strategy``.
    """
    require_float64()
    validate_options(step, strategy, checkpoints)

    ops = setup.ops
    n_f = ops.num_fields
    ts = jnp.asarray(setup.tlist)
    h = setup.h
    inner = setup.z_steps_inner
    N_z = jnp.asarray(setup.N_z)

    g_ang = jnp.asarray(setup.g) if g is None else TWO_PI * jnp.asarray(g)
    det0 = jnp.asarray(ops.detunings if detunings is None else detunings)
    det_classes = det0[None, :] + (
        jnp.asarray(setup.shifts)[:, None] * jnp.asarray(ops.doppler_factors)[None, :]
    )
    weights = jnp.asarray(setup.weights)

    # source: sum_pairs factor * rho_ab for each field, class-averaged
    rows = np.concatenate([p[0] for p in ops.pairs]) if n_f else np.zeros(0, int)
    cols = np.concatenate([p[1] for p in ops.pairs]) if n_f else np.zeros(0, int)
    fac = np.concatenate([p[2] for p in ops.pairs]) if n_f else np.zeros(0)
    owner = np.concatenate([np.full(len(p[0]), k) for k, p in enumerate(ops.pairs)])
    onehot = jnp.asarray(np.eye(n_f)[owner]) if n_f else jnp.zeros((0, 0))
    rows, cols, fac = jnp.asarray(rows), jnp.asarray(cols), jnp.asarray(fac)

    def observe(re, im):
        vr = (re[:, rows, cols] * fac) @ onehot
        vi = (im[:, rows, cols] * fac) @ onehot
        return jnp.stack([weights @ vr, weights @ vi])  # (2, n_f)

    if adjoint is None and strategy == "s4":
        adjoint = dfx.RecursiveCheckpointAdjoint(checkpoints=checkpoints)

    def coherence(rabi, jump_ts=None, real_input=False):
        ys = solve_slice(
            ops,
            rabi,
            det_classes,
            ts,
            decay_scale=decay_scale,
            decay_rates=decay_rates,
            observe=observe,
            rtol=rtol,
            atol=atol,
            n_steps=n_steps,
            max_steps=max_steps,
            solver=solver,
            adjoint=adjoint,
            jump_ts=jump_ts,
            b_field=b_field,
            real_drive=setup.real_drive or real_input,
        )  # (n_t, 2, n_f)
        return (ys[:, 0, :] + 1j * ys[:, 1, :]).T  # (n_f, n_t)

    spline = CubicDrive(setup.tlist)

    def spline_rabi(omega):
        cyc = omega / TWO_PI
        re = spline(jnp.real(cyc))
        if setup.real_drive:  # the imaginary part is never used
            return re
        im = spline(jnp.imag(cyc))
        return lambda t: re(t) + 1j * im(t)

    # input field and the drive of the first slice
    if omega0 is None:
        args = (
            list(pulse_args)
            if pulse_args is not None
            else [a for _, a, _ in setup.pulses]
        )
        if len(args) != n_f:
            raise ValueError(f"pulse_args has {len(args)} entries for {n_f} field(s).")
        pulses = [
            pulse_function(name, a)
            for (name, _, _), a in zip(setup.pulses, args, strict=True)
        ]
        if rabi_freqs is None:
            r0 = jnp.asarray([r for _, _, r in setup.pulses])
        else:
            r0 = jnp.asarray(rabi_freqs)

        def rabi0(t):
            return jnp.stack([r0[k] * pulses[k](t) for k in range(n_f)])

        omega_init = (TWO_PI * jax.vmap(rabi0)(ts).T).astype(jnp.complex128)
        jumps = [
            x
            for (name, _, _), a in zip(setup.pulses, args, strict=True)
            for x in jump_times(name, a)
        ]
        jump_ts = jnp.sort(jnp.asarray(jumps, dtype=float)) if jumps else None
    else:
        omega_init = jnp.asarray(omega0, dtype=jnp.complex128)
        rabi0 = spline_rabi(omega_init)
        jump_ts = None

    def advance(omega, coh, coh_prev, N):
        src = 2.0j * N * g_ang[:, None]
        if step == "euler":
            return omega + h * src * coh
        return omega + h * src * (1.5 * coh - 0.5 * coh_prev)

    # first inner step: analytic drive, Euler bootstrap (as clerq)
    # the analytic input pulses are real, so the first slice needs no imaginary term
    coh0 = coherence(rabi0, jump_ts, real_input=omega0 is None)
    omega1 = omega_init + h * 2.0j * N_z[0] * g_ang[:, None] * coh0

    def inner_step(carry, N):
        omega, coh_prev = carry
        coh = coherence(spline_rabi(omega))
        return (advance(omega, coh, coh_prev, N), coh), None

    carry = (omega1, coh0)
    if inner > 1:
        carry, _ = jax.lax.scan(inner_step, carry, N_z[1:inner])
    first = carry[0]

    def outer_step(carry, N_block):
        carry, _ = jax.lax.scan(inner_step, carry, N_block)
        return carry, carry[0]

    if strategy in ("s3", "s4"):
        outer_step = jax.checkpoint(outer_step)
    rest = N_z[inner:].reshape(setup.z_steps - 1, inner)
    _, outs = jax.lax.scan(outer_step, carry, rest)

    omega_zt = jnp.concatenate([omega_init[None], first[None], outs])  # (z+1, f, t)
    return jnp.transpose(omega_zt, (1, 0, 2))
