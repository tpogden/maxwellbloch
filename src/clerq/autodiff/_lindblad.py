"""Single-slice Lindblad solve, differentiable (the building block of propagation).

Solves, for each of ``n_v`` independent copies (velocity classes) of the system,

    d rho/dt = -i [H(t), rho] + D(rho),
    H(t) = H_0 + sum_k [ pi * Omega_k(t) X_k  -  2 pi s_k Delta_{v,k} N_k ]

with the operators of a :class:`~clerq.autodiff._operators.SystemOperators`. The
Rabi frequencies ``Omega_k(t)`` are cyclic and may be analytic or interpolated
from samples; they and the detunings, decay scale and initial state are all
differentiable.

``H`` is real, so ``rho = A + iB`` is integrated as the real pair ``(A, B)``:

    dA = [H, B] + Re D(rho),      dB = -[H, A] + Im D(rho).

Nothing complex enters diffrax, which flags complex dtypes as unreliable.
"""

from __future__ import annotations

from collections.abc import Callable

import diffrax as dfx
import jax.numpy as jnp
import numpy as np

from clerq.autodiff._config import require_float64
from clerq.autodiff._operators import SystemOperators


def check_fixed_step(ts, n_steps: int | None) -> None:
    """Raise if fixed-step boundaries coincide with interior save times.

    With ``RecursiveCheckpointAdjoint``, ``SaveAt(ts=...)`` and a
    ``ConstantStepSize`` whose steps land exactly on the output times, diffrax
    (0.7.x) returns *wrong* reverse-mode gradients while the forward values are
    right. Adaptive stepping is not affected in practice.

    Raises:
        ValueError: if ``n_steps`` puts an integer number of steps between all
            output times.
    """
    try:
        t = np.asarray(ts, dtype=float)
    except Exception:  # traced; cannot check
        return
    if not n_steps or len(t) < 3:
        return
    idx = (t[1:-1] - t[0]) * n_steps / (t[-1] - t[0])
    if np.any(np.abs(idx - np.round(idx)) < 1e-6):
        raise ValueError(
            f"n_steps={n_steps} puts fixed steps exactly on the output times, which "
            "makes reverse-mode gradients unreliable. Choose a step count that is "
            "not commensurate with the output grid."
        )


def solve_slice(
    ops: SystemOperators,
    rabi: Callable,
    detunings,
    ts,
    *,
    decay_scale=1.0,
    rho0=None,
    observe: Callable | None = None,
    rtol: float = 1e-8,
    atol: float = 1e-10,
    n_steps: int | None = None,
    max_steps: int = 100_000,
    solver=None,
    adjoint=None,
):
    """Solve the master equation at the times ``ts``.

    Args:
        ops: operators of the system.
        rabi: ``rabi(t) -> (n_fields,)`` cyclic Rabi frequencies at scalar time ``t``.
        detunings: cyclic detunings, shape ``(n_fields,)`` for one copy or
            ``(n_v, n_fields)`` for a batch (e.g. velocity classes).
        ts: output times (also the integration span).
        decay_scale: multiplies every decay rate.
        rho0: initial density matrices, ``(n, n)`` or ``(n_v, n, n)``; defaults to
            the system's initial state.
        observe: optional ``observe(re, im) -> real array`` applied at each output
            time to the real and imaginary parts of the batch ``(n_v, n, n)``, so
            only a reduction is stored.
        rtol, atol: adaptive step-size tolerances (ignored with ``n_steps``).
        n_steps: use a fixed-step solver with this many steps (for finite-difference
            reference gradients); must not be commensurate with ``ts``.
        max_steps: step budget of the adaptive solver.
        solver: a diffrax solver, default ``Tsit5``.
        adjoint: a diffrax adjoint, default ``RecursiveCheckpointAdjoint``.

    Returns:
        Complex ``rho(t)`` of shape ``(len(ts), n, n)`` for 1-D ``detunings`` (or
        ``(len(ts), n_v, n, n)`` for a batch), or the stacked ``observe`` output.

    Raises:
        RuntimeError: if JAX is not in 64-bit mode.
        ValueError: if ``n_steps`` is commensurate with ``ts``.
    """
    require_float64()
    check_fixed_step(ts, n_steps)

    ts = jnp.asarray(ts, dtype=jnp.float64)
    det = jnp.asarray(detunings, dtype=jnp.float64)
    single = det.ndim == 1
    if single:
        det = det[None, :]
    n_v, n = det.shape[0], ops.num_states

    X = jnp.asarray(ops.X)
    H_static = jnp.asarray(ops.H0)[None] - 2.0 * jnp.pi * jnp.einsum(
        "vk,kij->vij", det * jnp.asarray(ops.detuning_sign), jnp.asarray(ops.N)
    )
    eye = jnp.eye(n)

    sparse = ops.sparse
    if sparse is not None:
        a, b = jnp.asarray(sparse.a), jnp.asarray(sparse.b)
        w, r = jnp.asarray(sparse.w), jnp.asarray(sparse.r)
    elif len(ops.c_ops):
        # c = P + iQ; sum c^+ c = R + iS
        P, Q = jnp.asarray(ops.c_ops.real), jnp.asarray(ops.c_ops.imag)
        R = jnp.einsum("kji,kjl->il", P, P) + jnp.einsum("kji,kjl->il", Q, Q)
        S = jnp.einsum("kji,kjl->il", P, Q) - jnp.einsum("kji,kjl->il", Q, P)

    def rhs(t, y, args):
        rabi_fn, scale, H_s = args
        A, B = y[0], y[1]
        H = H_s + jnp.pi * jnp.einsum("k,kij->ij", rabi_fn(t), X)[None]
        dA = H @ B - B @ H
        dB = A @ H - H @ A
        if sparse is not None:
            zero = jnp.zeros((n_v, n))
            gain_a = zero.at[:, a].add(scale * w * A[:, b, b])
            gain_b = zero.at[:, a].add(scale * w * B[:, b, b])
            rate = scale * r
            damp = -0.5 * (rate[:, None] + rate[None, :])
            dA = dA + gain_a[:, :, None] * eye + damp * A
            dB = dB + gain_b[:, :, None] * eye + damp * B
        elif len(ops.c_ops):
            cr = jnp.einsum("kij,vjl->kvil", P, A) - jnp.einsum("kij,vjl->kvil", Q, B)
            ci = jnp.einsum("kij,vjl->kvil", P, B) + jnp.einsum("kij,vjl->kvil", Q, A)
            out_r = jnp.einsum("kvil,kml->vim", cr, P) + jnp.einsum(
                "kvil,kml->vim", ci, Q
            )
            out_i = jnp.einsum("kvil,kml->vim", ci, P) - jnp.einsum(
                "kvil,kml->vim", cr, Q
            )
            ac_r = R @ A - S @ B + A @ R - B @ S
            ac_i = R @ B + S @ A + A @ S + B @ R
            dA = dA + scale * (out_r - 0.5 * ac_r)
            dB = dB + scale * (out_i - 0.5 * ac_i)
        return jnp.stack([dA, dB])

    if rho0 is None:
        rho0 = ops.rho0
    rho0 = jnp.broadcast_to(jnp.asarray(rho0, dtype=jnp.complex128), (n_v, n, n))
    y0 = jnp.stack([jnp.real(rho0), jnp.imag(rho0)])

    save = (
        dfx.SaveAt(ts=ts)
        if observe is None
        else dfx.SaveAt(ts=ts, fn=lambda t, y, args: observe(y[0], y[1]))
    )
    sol = dfx.diffeqsolve(
        dfx.ODETerm(rhs),
        solver if solver is not None else dfx.Tsit5(),
        ts[0],
        ts[-1],
        (ts[-1] - ts[0]) / (n_steps or 1000),
        y0,
        args=(rabi, jnp.asarray(decay_scale, dtype=jnp.float64), H_static),
        saveat=save,
        stepsize_controller=(
            dfx.ConstantStepSize()
            if n_steps
            else dfx.PIDController(rtol=rtol, atol=atol)
        ),
        adjoint=adjoint if adjoint is not None else dfx.RecursiveCheckpointAdjoint(),
        max_steps=(n_steps + 2) if n_steps else max_steps,
    )
    if observe is not None:
        return sol.ys
    rho = sol.ys[:, 0] + 1j * sol.ys[:, 1]  # (n_t, n_v, n, n)
    return rho[:, 0] if single else rho
