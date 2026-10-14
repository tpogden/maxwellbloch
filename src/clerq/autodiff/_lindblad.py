"""Single-slice Lindblad solve, differentiable (the building block of propagation).

Solves, for each of ``n_v`` independent copies (velocity classes) of the system,

    d rho/dt = -i [H(t), rho] + D(rho),
    H(t) = H_0 + sum_k [ pi * Omega_k(t) X_k  -  2 pi s_k Delta_{v,k} N_k ]

with the operators of a :class:`~clerq.autodiff._operators.SystemOperators`. The
Rabi frequencies ``Omega_k(t)`` are cyclic and may be analytic or interpolated
from samples; they and the detunings, decay scale and initial state are all
differentiable.
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
        observe: optional ``observe(rho) -> array`` applied to the batch
            ``(n_v, n, n)`` at each output time, so only a reduction is stored.
        rtol, atol: adaptive step-size tolerances (ignored with ``n_steps``).
        n_steps: use a fixed-step solver with this many steps (for finite-difference
            reference gradients); must not be commensurate with ``ts``.
        max_steps: step budget of the adaptive solver.
        solver: a diffrax solver, default ``Tsit5``.
        adjoint: a diffrax adjoint, default ``RecursiveCheckpointAdjoint``.

    Returns:
        ``rho(t)`` of shape ``(len(ts), n, n)`` for 1-D ``detunings`` (or
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
    eye = jnp.eye(n, dtype=jnp.complex128)

    sparse = ops.sparse
    if sparse is not None:
        a, b = jnp.asarray(sparse.a), jnp.asarray(sparse.b)
        w, r = jnp.asarray(sparse.w), jnp.asarray(sparse.r)
    c = jnp.asarray(ops.c_ops)
    cdc = jnp.einsum("kji,kjl->il", c.conj(), c) if len(ops.c_ops) else None

    def rhs(t, rho, args):
        rabi_fn, scale, H_s = args
        H = H_s + jnp.pi * jnp.einsum("k,kij->ij", rabi_fn(t), X)[None]
        d = -1j * (H @ rho - rho @ H)
        if sparse is not None:
            gain = jnp.zeros((n_v, n), rho.dtype).at[:, a].add(scale * w * rho[:, b, b])
            d = d + gain[:, :, None] * eye
            rate = scale * r
            return d - 0.5 * (rate[:, None] + rate[None, :]) * rho
        if cdc is not None:
            d = d + scale * jnp.einsum("kij,vjl,kml->vim", c, rho, c.conj())
            d = d - 0.5 * scale * (cdc @ rho + rho @ cdc)
        return d

    if rho0 is None:
        rho0 = jnp.asarray(ops.rho0)
    rho0 = jnp.broadcast_to(jnp.asarray(rho0, dtype=jnp.complex128), (n_v, n, n))

    save = (
        dfx.SaveAt(ts=ts)
        if observe is None
        else dfx.SaveAt(ts=ts, fn=lambda t, y, args: observe(y))
    )
    sol = dfx.diffeqsolve(
        dfx.ODETerm(rhs),
        solver if solver is not None else dfx.Tsit5(),
        ts[0],
        ts[-1],
        (ts[-1] - ts[0]) / (n_steps or 1000),
        rho0,
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
    return sol.ys[:, 0] if (single and observe is None) else sol.ys
