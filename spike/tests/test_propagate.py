"""Phase 3: z-scan forward match vs clerq (T3.2, T3.4) and chain gradients (T3.3)."""

import copy
import json

import diffrax as dfx
import jax
import jax.numpy as jnp
import numpy as np
import pytest
from clerq import PropagationSolver

from mbjax.fd import central_fd, rel_err
from mbjax.lindblad import gaussian
from mbjax.propagate import build_propagation_model, propagate, velocity_classes
from mbjax.systems import make_system

jax.config.update("jax_enable_x64", True)
TWO_PI = 2 * np.pi


def _spec(name, ampl, fwhm=1.0, gamma=1.0, detuning=0.0):
    spec = make_system(name, gamma=gamma)
    spec["fields"][0].update(rabi_freq=1.0, detuning=detuning, rabi_freq_t_func="gaussian",
                             rabi_freq_t_args={"ampl": ampl, "centre": 0.0, "fwhm": fwhm})
    return spec


def _run_jax(spec, *, t_min, t_max, t_steps, z_steps, g, vc, ampl, fwhm=1.0, density=1.0,
             convention="clerq", **kw):
    model = build_propagation_model(spec)
    ts = jnp.linspace(t_min, t_max, t_steps + 1)
    h = 1.0 / z_steps
    z_next = (np.arange(z_steps) + 1) * h  # z_min = 0; square medium on [0, 1]
    N_z = jnp.asarray(np.where((z_next >= 0.0) & (z_next <= 1.0), density, 0.0))
    deltas, w = velocity_classes(vc, convention)
    pulse = gaussian(ampl, 0.0, fwhm)
    omega0 = TWO_PI * pulse(ts)
    return propagate(model, omega0, ts, N_z, g, h, deltas, w,
                     drive0=lambda t: jnp.array([pulse(t)]), **kw)


def _run_clerq(spec, *, t_min, t_max, t_steps, z_steps, g, vc, step="ab"):
    ps = PropagationSolver(
        system=copy.deepcopy(spec), t_min=t_min, t_max=t_max, t_steps=t_steps, z_min=0.0,
        z_max=1.0, z_steps=z_steps, z_steps_inner=1, interaction_strengths=[g],
        velocity_classes=vc,
        opts={"atol": 1e-12, "rtol": 1e-10, "method": "adams", "nsteps": 1_000_000})
    ps.solve(step=step, recalc=True, progress=False)
    return ps.Omegas_zt[0]


# ---- T3.2: the frozen reference problem -------------------------------------------------

def test_reference_problem_matches_saved_clerq_output():
    ref = np.load("refs/two_level_ref.npz")
    p = json.loads(str(ref["problem"]))
    spec = _spec(p["system"], ampl=p["ampl"], fwhm=p["fwhm"], gamma=p["gamma"])
    om = _run_jax(spec, t_min=p["t_min"], t_max=p["t_max"], t_steps=1000, z_steps=200,
                  g=p["g"], vc=None, ampl=p["ampl"], fwhm=p["fwhm"], rtol=1e-9, atol=1e-11)
    e = rel_err(np.asarray(om.real), ref["omega_zt"])
    assert e < 1e-6, e
    assert np.abs(np.asarray(om.imag)).max() < 1e-9  # resonant two level: real field


# ---- T3.4: multi-level systems and velocity classes -----------------------------------

VC = {"thermal_delta_min": -1.0, "thermal_delta_max": 1.0, "thermal_delta_steps": 4,
      "thermal_width": 0.7}
CASES = [("two_level", 1.0, 0.4), ("j12_12", 1.5, 0.4), ("j1_2_res", 2.5, 0.4),
         ("j72_72", 3.0, 1.5)]  # (system, pulse ampl, interaction strength g)


@pytest.mark.parametrize("name,ampl,g", CASES)
def test_velocity_classes_match_clerq(name, ampl, g):
    spec = _spec(name, ampl, detuning=0.2)
    kw = dict(t_min=-2.0, t_max=6.0, t_steps=160, z_steps=12, g=g, vc=VC)
    ref = _run_clerq(spec, **kw)
    # guard against a vacuous match: the medium must actually change the field in clerq
    assert np.abs(ref[-1] - ref[0]).max() > 1e-2 * np.abs(ref[0]).max(), name
    got = _run_jax(spec, ampl=ampl, **kw)
    e_r = rel_err(np.asarray(got.real), ref.real)
    assert e_r < 1e-6, (name, e_r)
    # clerq drops Im(Omega) from H but keeps it in the array: ours must too
    assert rel_err(np.asarray(got.imag), ref.imag) < 1e-5 or np.abs(ref.imag).max() < 1e-9


def test_physical_convention_differs_from_clerq():
    """Guard: the clerq Doppler convention is not an invisible choice."""
    spec = _spec("two_level", 1.0)
    kw = dict(t_min=-2.0, t_max=6.0, t_steps=120, z_steps=8, g=0.4, vc=VC, ampl=1.0)
    a = _run_jax(spec, convention="clerq", **kw)
    b = _run_jax(spec, convention="physical", **kw)
    assert rel_err(np.asarray(a.real), np.asarray(b.real)) > 1e-3


# ---- T3.3: gradients through the chain ------------------------------------------------

N_Z, N_T = 6, 60
N_FIXED = 1537  # prime: not commensurate with the 0.1 output spacing (see check_fixed_step)
BASIS_C = np.linspace(-1.5, 2.5, 4)


def _chain_loss(name, *, strategy="s3", step="ab", n_steps=None, adjoint=None, n_v=3):
    spec = _spec(name, 1.0)
    model = build_propagation_model(spec)
    ts = jnp.linspace(-2.0, 4.0, N_T + 1)
    h = 1.0 / N_Z
    N_z = jnp.ones(N_Z)
    deltas, w = velocity_classes(
        {"thermal_delta_min": -0.8, "thermal_delta_max": 0.8,
         "thermal_delta_steps": n_v - 1, "thermal_width": 0.6}, "clerq")
    basis = jnp.stack([jnp.exp(-((ts - c) ** 2) / 0.8) for c in BASIS_C])
    target = 0.5 * jnp.exp(-((ts - 1.0) ** 2) / 1.5)
    rng = np.random.default_rng(0)
    probe = jnp.asarray(rng.normal(size=ts.shape))

    def loss(theta):
        # theta = [g, gamma_scale, detuning, c_1..c_4]
        g, gs, det, coef = theta[0], theta[1], theta[2], theta[3:]
        om0 = TWO_PI * (coef @ basis)
        m = dict(model, detuning=det)
        om = propagate(m, om0, ts, N_z, g, h, deltas, w, strategy=strategy, step=step,
                       gamma_scale=gs, n_steps=n_steps, adjoint=adjoint, rtol=1e-10, atol=1e-12)
        out = om[-1].real / TWO_PI
        return jnp.sum((out - target) ** 2) * (6.0 / N_T) + 0.1 * jnp.vdot(probe, om[-1].imag).real

    theta0 = jnp.array([0.5, 1.1, 0.15, 0.4, 0.3, -0.2, 0.5])
    return loss, theta0


@pytest.mark.parametrize("name", ["two_level", "j12_12"])
def test_chain_grad_vs_fd(name):
    L_fixed, th = _chain_loss(name, n_steps=N_FIXED)
    g_fixed = np.asarray(jax.jit(jax.grad(L_fixed))(th))
    fd = central_fd(L_fixed, th, list(range(len(th))), 1e-4)
    assert rel_err(g_fixed, fd) < 1e-5, (rel_err(g_fixed, fd), g_fixed, fd)
    L_pid, _ = _chain_loss(name)
    g_pid = np.asarray(jax.jit(jax.grad(L_pid))(th))
    assert rel_err(g_pid, g_fixed) < 1e-4, rel_err(g_pid, g_fixed)


def test_chain_forward_mode_and_strategies_agree():
    th = _chain_loss("j12_12")[1]
    gs = {s: np.asarray(jax.jit(jax.grad(_chain_loss("j12_12", strategy=s)[0]))(th))
          for s in ("s2", "s3", "s5")}
    assert rel_err(gs["s3"], gs["s2"]) < 1e-9 and rel_err(gs["s5"], gs["s2"]) < 1e-9
    Lf, _ = _chain_loss("j12_12", strategy="s2", adjoint=dfx.ForwardMode())
    gf = np.asarray(jax.jit(jax.jacfwd(Lf))(th))
    assert rel_err(gs["s2"], gf) < 1e-8


def test_commensurate_fixed_step_is_rejected():
    th = _chain_loss("two_level", n_steps=1500)  # 25.0 steps per output interval
    with pytest.raises(ValueError, match="non-commensurate"):
        th[0](th[1])
