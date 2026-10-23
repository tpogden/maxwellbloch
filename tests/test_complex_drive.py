"""The complex-envelope (full RWA) drive: detuned fields give physical results.

An independent NumPy/SciPy solver (``refsolve.py``) is the oracle. The deprecated
``real_drive=True`` drives with ``Re(Omega)`` only, which is wrong for a single
detuned field.
"""

import copy
import json
import os
import sys
import warnings

import numpy as np
import pytest

jax = pytest.importorskip("jax")
pytest.importorskip("diffrax")
jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
import refsolve  # noqa: E402

from clerq import PropagationSolver  # noqa: E402
from clerq.autodiff._propagate import (  # noqa: E402
    propagate_fields,
    setup_from_solver,
)

pytestmark = pytest.mark.jax

JSON_DIR = os.path.join(os.path.dirname(__file__), "json")
TIGHT = {"atol": 1e-12, "rtol": 1e-10, "method": "adams", "nsteps": 1_000_000}
TWO_PI = 2 * np.pi
FWHM = 4.0


def two_level(detuning, real_drive=False, z_steps_inner=8):
    """Weak narrowband probe through an optical depth of 6 (rate 2, g = 3)."""
    return PropagationSolver(
        system={
            "num_states": 2,
            "decays": [{"channels": [[0, 1]], "rate": 2.0}],
            "fields": [
                {
                    "coupled_levels": [[0, 1]],
                    "detuning": detuning,
                    "rabi_freq": 1e-3,
                    "rabi_freq_t_func": "gaussian",
                    "rabi_freq_t_args": {"ampl": 1.0, "centre": 0.0, "fwhm": FWHM},
                }
            ],
        },
        t_min=-6.0,
        t_max=6.0,
        t_steps=120,
        z_min=0.0,
        z_max=1.0,
        z_steps=20,
        z_steps_inner=z_steps_inner,
        interaction_strengths=[3.0],
        opts=TIGHT,
        real_drive=real_drive,
    )


def transmission(fields_zt, tlist, field=0):
    """Intensity transmission of a field through the whole medium."""
    power = np.trapezoid(np.abs(fields_zt[field]) ** 2, tlist, axis=-1)
    return power[-1] / power[0]


def probe_input(solver):
    """The analytic Gaussian input of ``two_level`` as cyclic samples, shape (1, n_t)."""
    return (1e-3 * np.exp(-4 * np.log(2) * (solver.tlist / FWHM) ** 2))[None]


def solve(solver):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        solver.solve(recalc=True, progress=False)
    return solver


@pytest.mark.parametrize("detuning", [0.0, 0.5, 1.0, 2.0])
def test_single_detuned_class_is_lorentzian(detuning):
    """A weak narrowband probe sees exp(-OD / (1 + d^2)), OD = 2 g, d in half-widths."""
    s = solve(two_level(detuning))
    t = transmission(s.rabi_freqs_zt, s.tlist)
    assert t == pytest.approx(np.exp(-6.0 / (1 + detuning**2)), rel=0.05)
    assert t < 1.0


@pytest.mark.parametrize("detuning", [0.5, 2.0])
def test_propagation_solver_matches_reference(detuning):
    s = solve(two_level(detuning))
    ref = refsolve.reference(
        s, probe_input(s), False, zs=20, zi=8, rtol=1e-10, atol=1e-13
    )
    got = s.rabi_freqs_zt * TWO_PI  # (n_f, z + 1, n_t), angular
    np.testing.assert_allclose(
        np.transpose(ref, (1, 0, 2)), got, atol=1e-6 * np.abs(got).max(), rtol=0
    )


@pytest.mark.parametrize("detuning", [0.5, 2.0])
def test_autodiff_matches_propagation_solver_and_reference(detuning):
    s = solve(two_level(detuning))
    omega = np.asarray(
        propagate_fields(setup_from_solver(s), rtol=1e-10, atol=1e-12)
    )  # angular
    np.testing.assert_allclose(
        omega, s.rabi_freqs_zt * TWO_PI, atol=1e-6 * np.abs(omega).max(), rtol=0
    )
    t_ad = transmission(omega, s.tlist)
    ref = refsolve.reference(
        s, probe_input(s), False, zs=20, zi=8, rtol=1e-10, atol=1e-13
    )
    assert t_ad == pytest.approx(
        transmission(np.transpose(ref, (1, 0, 2)), s.tlist), rel=1e-4
    )
    # a detuned field acquires a phase, so the propagated field is complex
    assert np.abs(omega[0, -1].imag).max() > 1e-3 * np.abs(omega[0, -1]).max()


@pytest.mark.parametrize("detuning", [0.5, 2.0])
def test_real_drive_reproduces_the_old_wrong_values(detuning):
    """The deprecated Re(Omega) drive is unphysical (above 1 at d = 2)."""
    with pytest.warns(DeprecationWarning, match="real_drive"):
        s = two_level(detuning, real_drive=True)
    solve(s)
    t = transmission(s.rabi_freqs_zt, s.tlist)
    lorentz = np.exp(-6.0 / (1 + detuning**2))
    assert t > 1.5 * lorentz
    omega = np.asarray(propagate_fields(setup_from_solver(s), rtol=1e-10, atol=1e-12))
    np.testing.assert_allclose(
        omega, s.rabi_freqs_zt * TWO_PI, atol=1e-6 * np.abs(omega).max(), rtol=0
    )
    ref = refsolve.reference(
        s, probe_input(s), True, zs=20, zi=8, rtol=1e-10, atol=1e-13
    )
    np.testing.assert_allclose(
        np.transpose(ref, (1, 0, 2)), omega, atol=1e-6 * np.abs(omega).max(), rtol=0
    )


# ---- Lambda-EIT-like with a two-photon offset -------------------------------------


def eit_problem(real_drive=False):
    with open(os.path.join(JSON_DIR, "mbs_lambda_eit_with_coupling.json")) as f:
        d = json.load(f)
    d["system"]["fields"][0]["detuning"] = 0.4  # probe
    d["system"]["fields"][1]["detuning"] = -0.25  # coupling
    d.update(t_steps=80, z_steps=6, z_steps_inner=4, interaction_strengths=[5.0, 5.0])
    d.pop("savefile", None)
    if real_drive:
        with pytest.warns(DeprecationWarning):
            return PropagationSolver(**d, opts=TIGHT, real_drive=True)
    return PropagationSolver(**d, opts=TIGHT)


def eit_reference(s):
    om0 = np.zeros((2, len(s.tlist)))
    om0[:] = np.asarray(s.rabi_freqs_zt[:, 0, :].real)
    return refsolve.reference(s, om0, False, zs=6, zi=4, rtol=1e-10, atol=1e-13)


def test_detuned_lambda_eit_matches_reference():
    s = solve(eit_problem())
    got = s.rabi_freqs_zt * TWO_PI
    ref = np.transpose(eit_reference(s), (1, 0, 2))
    scale = np.abs(got[0]).max()
    np.testing.assert_allclose(ref[0], got[0], atol=1e-4 * scale, rtol=0)
    omega = np.asarray(propagate_fields(setup_from_solver(s), rtol=1e-10, atol=1e-12))
    np.testing.assert_allclose(omega, got, atol=1e-6 * np.abs(got).max(), rtol=0)
    # the old drive gives a different (wrong) answer for this asymmetric problem
    old = solve(eit_problem(real_drive=True))
    assert np.abs(old.rabi_freqs_zt[0] * TWO_PI - got[0]).max() > 1e-2 * scale


# ---- regression: symmetric resonant problems are unchanged ----------------------------


@pytest.mark.parametrize(
    "name",
    [
        "mb_solve_01.json",
        "velocity-classes.json",
        "mbs_lambda_eit_with_coupling.json",
        "mbs_two_sech_2pi.json",
    ],
)
def test_symmetric_resonant_problems_are_unchanged(name):
    """Zero detunings, real fields and a symmetric class grid: complex == real drive."""
    with open(os.path.join(JSON_DIR, name)) as f:
        d = json.load(f)
    d.pop("savefile", None)
    d["t_steps"] = min(d.get("t_steps", 100), 100)
    d["z_steps"] = min(d.get("z_steps", 10), 10)
    out = []
    for rd in (False, True):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            s = PropagationSolver(**copy.deepcopy(d), real_drive=rd)
            s.solve(recalc=True, progress=False)
        out.append(np.array(s.rabi_freqs_zt))
    scale = np.abs(out[1]).max()
    np.testing.assert_allclose(out[0], out[1], atol=1e-10 * scale, rtol=0)
    assert np.abs(out[0].imag).max() < 1e-10 * scale


def test_real_drive_warning_only_when_requested():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        two_level(0.0)


# ---- gradients --------------------------------------------------------------------

N_FIXED = 1507  # not a multiple of the output intervals


def test_gradients_through_the_complex_drive_match_finite_differences():
    """Gradients with respect to a *complex* input and a detuning."""
    small = PropagationSolver(
        system={
            "num_states": 2,
            "decays": [{"channels": [[0, 1]], "rate": 2.0}],
            "fields": [{"coupled_levels": [[0, 1]], "rabi_freq": 1e-3}],
        },
        t_min=-3.0,
        t_max=3.0,
        t_steps=30,
        z_steps=4,
        z_steps_inner=2,
        interaction_strengths=[2.0],
    )
    setup = setup_from_solver(small)
    n_t = len(setup.tlist)
    rng = np.random.default_rng(3)
    env = TWO_PI * 0.3 * np.exp(-0.5 * (setup.tlist / 1.2) ** 2)
    re0 = env * (1 + 0.1 * rng.normal(size=n_t))
    im0 = env * 0.5 * (1 + 0.1 * rng.normal(size=n_t))
    wt = jnp.asarray(np.exp(-0.5 * ((setup.tlist - 0.5) / 1.0) ** 2))

    def loss(x, **kw):
        o = (x[:n_t] + 1j * x[n_t : 2 * n_t])[None]
        out = propagate_fields(setup, omega0=o, detunings=x[2 * n_t :], **kw)
        z = out[0, -1]
        return jnp.sum(jnp.abs(z) ** 2 * wt) + jnp.sum(jnp.imag(z) * wt)

    x0 = np.concatenate([re0, im0, [0.4]])
    f = jax.jit(lambda x: loss(x, n_steps=N_FIXED))
    g = np.asarray(jax.jit(jax.grad(lambda x: loss(x, n_steps=N_FIXED)))(x0))
    idx = [0, 5, 10, 15, 29, n_t, n_t + 5, n_t + 10, n_t + 20, 2 * n_t]
    fd = []
    for i in idx:
        e = np.zeros_like(x0)
        e[i] = 1e-4
        fd.append((float(f(x0 + e)) - float(f(x0 - e))) / 2e-4)
    fd, ga = np.array(fd), g[idx]
    assert np.linalg.norm(ga - fd) / np.linalg.norm(fd) < 1e-5
    assert np.abs(g[n_t : 2 * n_t]).max() > 1e-4  # the imaginary part matters
