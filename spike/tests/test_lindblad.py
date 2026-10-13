"""Phase 1: single-slice Lindblad solve in JAX vs clerq (T1.2) and gradients (T1.3)."""

import copy

import diffrax as dfx
import jax
import jax.numpy as jnp
import numpy as np
import pytest
from clerq import LindbladSolver

from mbjax.fd import central_fd, rel_err
from mbjax.lindblad import build_model, gaussian, solve_slice
from mbjax.systems import DIMS, make_system

jax.config.update("jax_enable_x64", True)

NAMES = list(DIMS)
AMPL, CENTRE, FWHM, DELTA = 0.8, 0.0, 1.0, 0.3
T_MIN, T_MAX, T_STEPS = -2.0, 3.0, 250


def _spec(name):
    spec = make_system(name)
    f = spec["fields"][0]
    f.update(
        rabi_freq=1.0,
        detuning=DELTA,
        rabi_freq_t_func="gaussian",
        rabi_freq_t_args={"ampl": AMPL, "centre": CENTRE, "fwhm": FWHM},
    )
    return spec


def _clerq_states(name):
    s = LindbladSolver(
        system=copy.deepcopy(_spec(name)),
        t_min=T_MIN,
        t_max=T_MAX,
        t_steps=T_STEPS,
        opts={"atol": 1e-12, "rtol": 1e-10, "method": "adams", "nsteps": 1_000_000},
    )
    s.solve(save=False)
    return s.tlist, s.states_t()


@pytest.mark.parametrize("name", NAMES)
def test_matches_clerq(name):
    tlist, ref = _clerq_states(name)
    model = build_model(_spec(name))
    rho = solve_slice(
        model, lambda t: jnp.array([gaussian(AMPL, CENTRE, FWHM)(t)]), [DELTA], jnp.asarray(tlist)
    )
    err = np.abs(np.asarray(rho) - ref).max()
    assert err < 1e-8, err
    # physical sanity
    tr = np.trace(np.asarray(rho), axis1=1, axis2=2)
    assert np.allclose(tr, 1.0, atol=1e-8)
    assert np.allclose(rho, np.conj(np.swapaxes(rho, 1, 2)), atol=1e-8)


# ---- T1.3 gradients -------------------------------------------------------------------

def _loss_factory(name, adjoint=None):
    model = build_model(_spec(name))
    n = model["n"]
    ts = jnp.linspace(T_MIN, T_MAX, 120)
    rng = np.random.default_rng(1)
    W = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    W = jnp.asarray(W + W.conj().T)

    def loss(theta):  # theta = [ampl, detuning]
        rho = solve_slice(
            model,
            lambda t: jnp.array([gaussian(theta[0], CENTRE, FWHM)(t)]),
            theta[1:2],
            ts,
            rtol=1e-11,
            atol=1e-13,
            adjoint=adjoint,
        )
        return jnp.sum(jnp.real(jnp.einsum("ij,tji->t", W, rho)))

    return loss


@pytest.mark.parametrize("name", NAMES)
def test_grad_vs_fd_and_forward(name):
    theta = jnp.array([AMPL, DELTA])
    g = np.asarray(jax.jit(jax.grad(_loss_factory(name)))(theta))
    fd = central_fd(_loss_factory(name), theta, [0, 1], 1e-5)
    assert rel_err(g, fd) < 1e-6, (g, fd)
    gf = np.asarray(jax.jit(jax.jacfwd(_loss_factory(name, dfx.ForwardMode())))(theta))
    assert rel_err(g, gf) < 1e-8, (g, gf)
