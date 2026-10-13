"""Phase 2: array-driven slice. T2.1 forward vs clerq's interpolated path; T2.2 gradients
w.r.t. the field samples."""

import copy

import diffrax as dfx
import jax
import jax.numpy as jnp
import numpy as np
import pytest
import qutip as qu
from clerq import PropagationSolver
from scipy.interpolate import make_interp_spline

from mbjax.drive import CubicDrive, linear_drive
from mbjax.fd import central_fd, rel_err
from mbjax.lindblad import build_model, solve_slice
from mbjax.systems import DIMS, make_system

jax.config.update("jax_enable_x64", True)

NAMES = list(DIMS)
T_MIN, T_MAX, N_T = -2.0, 3.0, 201
TS = np.linspace(T_MIN, T_MAX, N_T)
DELTA = 0.3


def _field(seed=0):
    """Smooth non-Gaussian cyclic-unit field samples (pulse + chirp-ish ripple)."""
    rng = np.random.default_rng(seed)
    base = 0.8 * np.exp(-4 * np.log(2) * TS**2) + 0.25 * np.sin(2.0 * TS) * np.exp(-0.3 * TS**2)
    return base + 0.02 * rng.normal(size=N_T)  # rough enough to stress the spline


def _spec(name):
    spec = make_system(name)
    spec["fields"][0].update(rabi_freq=1.0, detuning=DELTA, rabi_freq_t_func="gaussian",
                             rabi_freq_t_args={"ampl": 1.0, "centre": 0.0, "fwhm": 1.0})
    return spec


# ---- the spline itself ----------------------------------------------------------------

def test_cubic_matches_scipy():
    y = _field()
    f = CubicDrive(TS)(jnp.asarray(y))
    tt = np.random.default_rng(3).uniform(T_MIN, T_MAX, 500)
    ref = make_interp_spline(TS, y, k=3, bc_type=None)(tt)
    assert np.abs(np.array([float(f(t)) for t in tt[:100]]) - ref[:100]).max() < 1e-11
    assert np.abs(np.asarray(jax.vmap(f)(jnp.asarray(tt))) - ref).max() < 1e-11
    assert np.abs(np.asarray(jax.vmap(f)(jnp.asarray(TS))) - y).max() < 1e-12  # interpolates


# ---- T2.1: forward vs clerq's actual interpolated drive --------------------------------

def _clerq_rho(name, field_cyc):
    spec = _spec(name)
    ps = PropagationSolver(system=copy.deepcopy(spec), t_min=T_MIN, t_max=T_MAX,
                           t_steps=N_T - 1, z_steps=1,
                           opts={"atol": 1e-12, "rtol": 1e-10, "method": "adams", "nsteps": 1_000_000})
    Om = np.zeros((1, N_T), dtype=complex)
    Om[0] = 2 * np.pi * field_cyc  # clerq stores angular Rabi frequency
    ps.system.H_Omega_list = ps._build_intp_H_Omega_list(Om)
    ps.system.mesolve(ps.tlist, options=ps.opts, recalc=True)
    return ps.tlist, ps.system.states_t()


@pytest.mark.parametrize("name", NAMES)
def test_array_drive_matches_clerq(name):
    y = _field()
    tl, ref = _clerq_rho(name, y)
    model = build_model(_spec(name))
    drive = CubicDrive(TS)(jnp.asarray(y))
    rho = solve_slice(model, lambda t: jnp.array([drive(t)]), [DELTA], jnp.asarray(tl))
    err = np.abs(np.asarray(rho) - ref).max()
    assert err < 1e-8, err


def test_linear_drive_is_actually_different():
    """Guard: matching clerq is a statement about the *cubic* spline, not any smooth drive."""
    y = _field()
    model = build_model(_spec("two_level"))
    cub = CubicDrive(TS)(jnp.asarray(y))
    lin = linear_drive(TS)(jnp.asarray(y))
    a = solve_slice(model, lambda t: jnp.array([cub(t)]), [DELTA], jnp.asarray(TS))
    b = solve_slice(model, lambda t: jnp.array([lin(t)]), [DELTA], jnp.asarray(TS))
    assert np.abs(np.asarray(a) - np.asarray(b)).max() > 1e-5


# ---- T2.2: gradients w.r.t. the field samples ------------------------------------------

def _loss(name, kind, adjoint=None, n_steps=None):
    model = build_model(_spec(name))
    n = model["n"]
    rng = np.random.default_rng(1)
    W = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    W = jnp.asarray(W + W.conj().T)
    ts_out = jnp.linspace(T_MIN, T_MAX, 80)
    mk = CubicDrive(TS) if kind == "cubic" else linear_drive(TS)

    def loss(y):
        drive = mk(y)
        rho = solve_slice(model, lambda t: jnp.array([drive(t)]), [DELTA], ts_out,
                          rtol=1e-11, atol=1e-13, adjoint=adjoint, n_steps=n_steps)
        return jnp.sum(jnp.real(jnp.einsum("ij,tji->t", W, rho)))

    return loss


IDX = [0, 1, 2, 50, 100, 101, 150, 197, 198, 199, 200]  # edges stress not-a-knot
N_FIXED = 20000


@pytest.mark.parametrize("kind", ["cubic", "linear"])
@pytest.mark.parametrize("name", ["two_level", "j12_12"])
def test_field_grad_vs_fd(name, kind):
    """AD vs central FD on a fixed-step solve (smooth map), per amendment A5; then the
    production path (adaptive PID, rtol 1e-11) must give the same gradient as fixed-step AD."""
    y = jnp.asarray(_field())
    L_fixed = _loss(name, kind, n_steps=N_FIXED)
    g_fixed = np.asarray(jax.jit(jax.grad(L_fixed))(y))
    fd = central_fd(L_fixed, y, IDX, 1e-3)
    assert rel_err(g_fixed[IDX], fd) < 1e-6, rel_err(g_fixed[IDX], fd)
    g_pid = np.asarray(jax.jit(jax.grad(_loss(name, kind)))(y))
    assert rel_err(g_pid, g_fixed) < 1e-4, rel_err(g_pid, g_fixed)


def test_field_grad_reverse_vs_forward_mode():
    y = jnp.asarray(_field())
    g = np.asarray(jax.jit(jax.grad(_loss("j1_2_res", "cubic")))(y))
    gf = np.asarray(jax.jit(jax.jacfwd(_loss("j1_2_res", "cubic", dfx.ForwardMode())))(y))
    assert rel_err(g, gf) < 1e-8
