"""Per-decay rates as autodiff parameters (``decay_rates``)."""

import copy
import json

import numpy as np
import pytest

jax = pytest.importorskip("jax")
pytest.importorskip("diffrax")
jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402

from clerq import PropagationSolver  # noqa: E402
from clerq.autodiff._operators import build_operators  # noqa: E402

pytestmark = pytest.mark.jax

TIGHT = {"atol": 1e-12, "rtol": 1e-10, "method": "adams", "nsteps": 1_000_000}
N_FIXED = 1507  # not a multiple of the 40 output intervals
TWO_PI = 2 * np.pi

# A Lambda system with spontaneous decay from the excited level (two channels),
# decoherence of the ground-state coherence through a pure-dephasing channel
# ``[2, 2]``, and population relaxation between the ground states.
PROBLEM = {
    "system": {
        "num_states": 3,
        "decays": [
            {"rate": 1.0, "channels": [[0, 1], [2, 1]], "factors": [1.0, 0.7]},
            {"rate": 0.2, "channels": [[2, 2]]},
            {"rate": 0.05, "channels": [[0, 2]]},
        ],
        "fields": [
            {
                "coupled_levels": [[0, 1]],
                "detuning": 0.0,
                "rabi_freq": 0.05,
                "rabi_freq_t_func": "gaussian",
                "rabi_freq_t_args": {"ampl": 1.0, "centre": 1.5, "fwhm": 1.0},
            },
            {
                "coupled_levels": [[1, 2]],
                "detuning": 0.0,
                "detuning_positive": False,
                "rabi_freq": 1.5,
                "rabi_freq_t_func": "gaussian",
                "rabi_freq_t_args": {"ampl": 1.0, "centre": 1.5, "fwhm": 4.0},
            },
        ],
    },
    "t_min": 0.0,
    "t_max": 3.0,
    "t_steps": 40,
    "z_min": 0.0,
    "z_max": 1.0,
    "z_steps": 4,
    "z_steps_inner": 2,
    "interaction_strengths": [4.0, 4.0],
}

_WT = np.exp(-0.5 * ((np.linspace(0.0, 3.0, 41) - 1.8) / 0.9) ** 2)


def loss(omega_zt):
    """A scalar of the transmitted probe (depends on every decay)."""
    return jnp.sum(jnp.abs(omega_zt[0, -1]) ** 2 * _WT)


def _solver(problem=PROBLEM, opts=TIGHT):
    ps = PropagationSolver.from_json_str(json.dumps(problem))
    ps.build_opts(opts)
    return ps


def _with_rates(rates):
    d = copy.deepcopy(PROBLEM)
    for dec, rate in zip(d["system"]["decays"], rates, strict=True):
        dec["rate"] = float(rate)
    return d


def _solve(problem):
    ps = _solver(problem)
    ps.solve(recalc=True, progress=False)
    return ps.rabi_freqs_zt


def test_decay_channels_reproduce_the_collapse_operators():
    ops = build_operators(_solver().system)
    ch = ops.decays
    assert ch is not None
    np.testing.assert_allclose(ch.rates, [1.0, 0.2, 0.05])
    # w_k = 2 pi rate f^2 is what the sparse dissipator holds
    w = TWO_PI * ch.rates[ch.idx] * ch.f2
    np.testing.assert_allclose(w, ops.sparse.w, rtol=1e-12)
    np.testing.assert_array_equal(ch.a, ops.sparse.a)
    np.testing.assert_array_equal(ch.b, ops.sparse.b)


def test_default_decay_rates_are_the_problems_own():
    p = _solver().to_autodiff().default_params()
    np.testing.assert_allclose(p["decay_rates"], [1.0, 0.2, 0.05])
    assert p["decay_scale"] == 1.0


def test_default_rates_give_the_same_field_as_the_problem():
    ad = _solver().to_autodiff(rtol=1e-10, atol=1e-12)
    ref = _solve(PROBLEM)
    got = np.asarray(ad.propagate())
    np.testing.assert_allclose(got, ref, atol=1e-7 * np.abs(ref).max(), rtol=0)


@pytest.mark.parametrize("which", [0, 1, 2], ids=["spontaneous", "dephasing", "ground"])
def test_changed_rate_matches_solving_the_changed_problem(which):
    """The rate of each kind of channel, including the pure-dephasing ``[2, 2]``."""
    rates = np.array([1.0, 0.2, 0.05])
    rates[which] *= 2.5
    ad = _solver().to_autodiff(rtol=1e-10, atol=1e-12)
    ref = _solve(_with_rates(rates))
    base = _solve(PROBLEM)
    got = np.asarray(ad.propagate({"decay_rates": rates}))
    scale = np.abs(ref).max()
    np.testing.assert_allclose(got, ref, atol=1e-7 * scale, rtol=0)
    # non-vacuous: the rate changes the transmitted probe
    assert np.abs(ref[0, -1] - base[0, -1]).max() > 1e-4 * scale


def test_a_zero_rate_is_allowed_and_matches_the_problem_without_it():
    rates = np.array([1.0, 0.0, 0.05])
    ad = _solver().to_autodiff(rtol=1e-10, atol=1e-12)
    got = np.asarray(ad.propagate({"decay_rates": rates}))
    ref = _solve(_with_rates(rates))
    np.testing.assert_allclose(got, ref, atol=1e-7 * np.abs(ref).max(), rtol=0)


def test_decay_scale_multiplies_the_rates():
    ad = _solver().to_autodiff(rtol=1e-10, atol=1e-12)
    a = np.asarray(ad.propagate({"decay_scale": 2.0}))
    b = np.asarray(ad.propagate({"decay_rates": 2.0 * np.array([1.0, 0.2, 0.05])}))
    np.testing.assert_allclose(a, b, atol=1e-12, rtol=0)


def _fd(f, x, h=1e-5, one_sided=False):
    g = np.empty(len(x))
    for i in range(len(x)):
        e = np.zeros(len(x))
        e[i] = h
        if one_sided:
            g[i] = (float(f(jnp.asarray(x + e))) - float(f(jnp.asarray(x)))) / h
        else:
            g[i] = (float(f(jnp.asarray(x + e))) - float(f(jnp.asarray(x - e)))) / (
                2 * h
            )
    return g


@pytest.mark.parametrize("mode", ["reverse", "forward"])
def test_gradient_with_respect_to_the_rates_matches_finite_differences(mode):
    ad = _solver().to_autodiff(n_steps=N_FIXED)
    x0 = np.array([1.0, 0.2, 0.05])

    def f(x):
        return loss(ad.propagate({"decay_rates": x}))

    _, g = ad.value_and_grad(loss, {"decay_rates": x0}, mode=mode)
    g = np.asarray(g["decay_rates"])
    fd = _fd(f, x0)
    assert np.all(np.abs(g) > 1e-6)  # every rate matters
    assert np.linalg.norm(g - fd) / np.linalg.norm(fd) < 1e-5


def test_gradient_at_a_zero_rate_is_finite_and_correct():
    """The dissipator is linear in the rate, so rate = 0 is a regular point."""
    ad = _solver().to_autodiff(n_steps=N_FIXED)
    x0 = np.array([1.0, 0.0, 0.05])

    def f(x):
        return loss(ad.propagate({"decay_rates": x}))

    _, g = ad.value_and_grad(loss, {"decay_rates": x0})
    g = np.asarray(g["decay_rates"])
    assert np.all(np.isfinite(g))
    fd = _fd(f, x0, h=1e-6, one_sided=True)
    assert abs(g[1] - fd[1]) < 1e-4 * abs(fd[1])


def test_decay_rates_can_be_fitted():
    """Gradient descent on the rates fits a transmitted probe made with other rates.

    The fit is not unique: one transmitted probe constrains a combination of the
    spontaneous, dephasing and ground-state rates, so this checks the mismatch and
    the direction of the dephasing rate, not that each rate is recovered.
    """
    ad = _solver().to_autodiff(rtol=1e-10, atol=1e-12)
    target = np.asarray(ad.propagate({"decay_rates": np.array([1.0, 0.6, 0.05])}))

    def mismatch(omega_zt):
        return jnp.sum(jnp.abs(omega_zt[0, -1] - target[0, -1]) ** 2)

    start = {"decay_rates": np.array([1.0, 0.2, 0.05])}
    v0, _ = ad.value_and_grad(mismatch, start)
    res = ad.minimize(mismatch, start, learning_rate=0.05, steps=80)
    assert res.value < 0.01 * float(v0)
    assert res.params["decay_rates"][1] > 0.3  # moved towards 0.6 from 0.2


def test_wrong_shape_is_rejected():
    ad = _solver().to_autodiff()
    with pytest.raises(ValueError, match="decay_rates must have shape"):
        ad.propagate({"decay_rates": np.array([1.0, 0.2])})


def test_system_without_decays_has_no_rates():
    d = copy.deepcopy(PROBLEM)
    d["system"]["decays"] = []
    ad = _solver(d).to_autodiff()
    assert ad.default_params()["decay_rates"].shape == (0,)
    ad.propagate()  # still runs
    with pytest.raises(ValueError, match="decay_rates must have shape"):
        ad.propagate({"decay_rates": np.array([0.1])})
