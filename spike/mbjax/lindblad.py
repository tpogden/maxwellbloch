"""Phase 1: Lindblad master equation for one slice, in JAX/diffrax.

Conventions are clerq's (see systems.py): frequencies are cyclic,

    H(t) = sum_k [ pi * Omega_k(t) * X_k  -  2*pi * Delta_k * N_k ]
    d rho/dt = -i [H, rho] + sum_j ( c_j rho c_j^+ - 1/2 {c_j^+ c_j, rho} )

rho is a complex (n, n) matrix carried directly through diffrax. The Rabi frequency is a
callable ``rabi(t) -> (n_fields,)`` in cyclic units, so a drive can be analytic (Phase 1)
or an interpolated array (Phase 2).
"""

from __future__ import annotations

import diffrax as dfx
import jax
import jax.numpy as jnp

from mbjax.systems import operators

jax.config.update("jax_enable_x64", True)


def build_model(spec: dict) -> dict:
    """Operator pieces of a system spec as JAX arrays."""
    ops = operators(spec)
    n = ops["n"]
    c = jnp.stack([jnp.asarray(c) for c in ops["c_ops"]]) if ops["c_ops"] else jnp.zeros((0, n, n), complex)
    return {
        "n": n,
        "X": jnp.stack([jnp.asarray(x) for x in ops["X"]]),
        "N": jnp.stack([jnp.asarray(x) for x in ops["N"]]),
        "c": c,
        "cdc": jnp.einsum("kji,kjl->il", c.conj(), c) if c.shape[0] else jnp.zeros((n, n), complex),
        "rho0": jnp.asarray(ops["rho0"]),
    }


def gaussian(ampl, centre, fwhm):
    """clerq's gaussian t_func in cyclic units: ampl * exp(-4 ln2 ((t-c)/fwhm)^2)."""
    k = 4.0 * jnp.log(2.0) / fwhm**2
    return lambda t: ampl * jnp.exp(-k * (t - centre) ** 2)


def rhs(t, rho, args):
    model, rabi, detunings = args
    H = jnp.pi * jnp.einsum("k,kij->ij", rabi(t), model["X"]) - 2.0 * jnp.pi * jnp.einsum(
        "k,kij->ij", detunings, model["N"]
    )
    drho = -1j * (H @ rho - rho @ H)
    drho += jnp.einsum("kij,jl,kml->im", model["c"], rho, model["c"].conj())
    drho -= 0.5 * (model["cdc"] @ rho + rho @ model["cdc"])
    return drho


def solve_slice(
    model,
    rabi,
    detunings,
    ts,
    *,
    rho0=None,
    rtol=1e-10,
    atol=1e-12,
    adjoint=None,
    max_steps=1_000_000,
):
    """Solve rho(t) at times ``ts`` for one slice. Returns array (len(ts), n, n)."""
    sol = dfx.diffeqsolve(
        dfx.ODETerm(rhs),
        dfx.Tsit5(),
        ts[0],
        ts[-1],
        (ts[-1] - ts[0]) / 1000.0,
        model["rho0"] if rho0 is None else rho0,
        args=(model, rabi, jnp.atleast_1d(jnp.asarray(detunings))),
        saveat=dfx.SaveAt(ts=ts),
        stepsize_controller=dfx.PIDController(rtol=rtol, atol=atol),
        adjoint=adjoint or dfx.RecursiveCheckpointAdjoint(),
        max_steps=max_steps,
    )
    return sol.ys
