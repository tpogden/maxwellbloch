"""The public autodiff API: ``PropagationSolver.to_autodiff()`` and its parameters."""

import copy
import json

import numpy as np
import pytest

jax = pytest.importorskip("jax")
pytest.importorskip("diffrax")
jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402
from jax.flatten_util import ravel_pytree  # noqa: E402

from clerq import PropagationSolver  # noqa: E402
from clerq.autodiff import AutodiffPropagation  # noqa: E402
from clerq.autodiff._api import _merge_pulse_args  # noqa: E402

pytestmark = pytest.mark.jax

TIGHT = {"atol": 1e-12, "rtol": 1e-10, "method": "adams", "nsteps": 1_000_000}

PROBLEM = {
    "system": {
        "num_states": 4,
        "energies": [0.0, 0.3, 1.2, 2.0],
        "decays": [
            {"rate": 0.7, "channels": [[0, 2], [1, 2]], "factors": [0.6, 0.8]},
            {"rate": 1.3, "channels": [[2, 3]]},
        ],
        "fields": [
            {
                "coupled_levels": [[0, 2], [1, 2]],
                "factors": [0.5, -0.8],
                "detuning": 0.4,
                "detuning_positive": False,
                "rabi_freq": 1.7,
                "rabi_freq_t_func": "gaussian",
                "rabi_freq_t_args": {"ampl": 1.2, "centre": 1.0, "fwhm": 0.8},
            },
            {
                "coupled_levels": [[2, 3]],
                "detuning": -0.2,
                "factor_doppler_shift": 0.5,
                "rabi_freq": 0.9,
                "rabi_freq_t_func": "sech",
                "rabi_freq_t_args": {"ampl": 0.9, "centre": 1.5, "width": 0.4},
            },
        ],
        "initial_state": [0.5, 0.5, 0.0, 0.0],
    },
    "t_min": 0.0,
    "t_max": 3.0,
    "t_steps": 40,
    "z_min": 0.0,
    "z_max": 1.0,
    "z_steps": 4,
    "z_steps_inner": 1,
    "interaction_strengths": [2.0, 1.0],
    "velocity_classes": {
        "thermal_delta_min": -1.0,
        "thermal_delta_max": 1.0,
        "thermal_delta_steps": 2,
        "thermal_width": 0.7,
    },
}
N_FIXED = 2003  # not a multiple of the 40 output intervals


def _solver(opts=TIGHT):
    ps = PropagationSolver.from_json_str(json.dumps(PROBLEM))
    ps.build_opts(opts)
    return ps


_WT = np.exp(-0.5 * ((np.linspace(0.0, 3.0, 41) - 1.5) / 0.8) ** 2)
_PROBE = np.random.default_rng(3).normal(size=(2, 41))


def loss(omega_zt):
    """A scalar of the output field that also depends on its imaginary part."""
    last = omega_zt[:, -1, :]
    return jnp.sum(jnp.abs(last) ** 2 * _WT) + jnp.sum(jnp.imag(last) * _PROBE)


def _flat_fd(f, x, h=1e-5):
    g = np.empty(len(x))
    for i in range(len(x)):
        e = np.zeros(len(x))
        e[i] = h
        g[i] = (float(f(jnp.asarray(x + e))) - float(f(jnp.asarray(x - e)))) / (2 * h)
    return g


def _rel(a, b):
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300))


# ---- entry point ---------------------------------------------------------------


def test_to_autodiff_returns_the_differentiable_object():
    assert isinstance(_solver().to_autodiff(), AutodiffPropagation)


@pytest.mark.parametrize("step", ["euler", "ab"])
def test_propagate_matches_solve(step):
    ps = _solver()
    ad = ps.to_autodiff(step=step, rtol=1e-11, atol=1e-13)
    ps.solve(step=step, recalc=True, progress=False)
    got = np.asarray(ad.propagate())
    ref = ps.Omegas_zt
    assert got.shape == ref.shape
    assert np.abs(got - ref).max() / np.abs(ref).max() < 1e-8
    assert np.abs(ref[:, -1] - ref[:, 0]).max() > 1e-3  # a real effect


def test_default_params_are_the_problems_own():
    ps = _solver()
    p = ps.to_autodiff().default_params()
    f = ps.system.fields
    assert np.allclose(p["interaction_strengths"], ps.interaction_strengths)
    assert np.allclose(p["detunings"], [x.detuning for x in f])
    assert np.allclose(p["rabi_freq"], [x.rabi_freq for x in f])
    assert p["pulse_args"][0] == {"ampl": 1.2, "centre": 1.0, "fwhm": 0.8}
    assert p["pulse_args"][1] == {"ampl": 0.9, "centre": 1.5, "width": 0.4}
    assert p["decay_scale"] == 1.0


# ---- gradients --------------------------------------------------------------------


def test_gradients_of_every_user_facing_parameter_match_finite_differences():
    """Pulse parameters, rabi frequencies, interaction strengths, detunings and the
    decay scale: AD (reverse) vs central FD on a fixed-step solve."""
    ad = _solver().to_autodiff(n_steps=N_FIXED)
    params = jax.tree_util.tree_map(jnp.asarray, ad.default_params())
    flat, unravel = ravel_pytree(params)
    assert flat.size == 3 + 3 + 2 + 2 + 2 + 1  # 2 pulses (3 args each) + 7 others

    _, grads = ad.value_and_grad(loss, params)
    g_ad = np.asarray(ravel_pytree(grads)[0])
    fd = _flat_fd(lambda x: loss(ad.propagate(unravel(x))), np.asarray(flat))
    assert _rel(g_ad, fd) < 1e-5, (_rel(g_ad, fd), g_ad, fd)
    assert np.abs(g_ad).min() > 0.0 and np.abs(g_ad).max() > 1e-3  # non-vacuous


def test_reverse_and_forward_mode_agree():
    ad = _solver().to_autodiff(rtol=1e-10, atol=1e-12)
    params = ad.default_params()
    v_rev, g_rev = ad.value_and_grad(loss, params, mode="reverse")
    v_fwd, g_fwd = ad.value_and_grad(loss, params, mode="forward")
    assert abs(float(v_rev) - float(v_fwd)) < 1e-10 * max(1.0, abs(float(v_rev)))
    assert (
        _rel(np.asarray(ravel_pytree(g_fwd)[0]), np.asarray(ravel_pytree(g_rev)[0]))
        < 1e-8
    )


def test_gradient_of_a_subset_of_parameters():
    ad = _solver().to_autodiff(rtol=1e-10, atol=1e-12)
    full = ad.value_and_grad(loss, ad.default_params())[1]
    only = ad.value_and_grad(loss, {"interaction_strengths": np.array([2.0, 1.0])})[1]
    assert set(only) == {"interaction_strengths"}
    assert np.allclose(only["interaction_strengths"], full["interaction_strengths"])


def test_input_field_samples_are_a_parameter():
    ad = _solver().to_autodiff(rtol=1e-10, atol=1e-12)
    analytic = np.asarray(ad.propagate())
    # omega0 is a *cyclic* Rabi frequency, like rabi_freq: the output field is
    # angular (the layout and units of Omegas_zt), so divide by 2 pi
    omega0 = analytic[:, 0, :].real / (2.0 * np.pi)
    sampled = np.asarray(ad.propagate({"omega0": omega0}))
    # the first slice is then driven by the spline of the samples, not the analytic pulse
    assert np.abs(sampled - analytic).max() / np.abs(analytic).max() < 1e-4
    _, g = ad.value_and_grad(loss, {"omega0": omega0})
    assert g["omega0"].shape == omega0.shape
    assert np.abs(np.asarray(g["omega0"])).max() > 1e-6


def test_omega0_units_are_cyclic_not_angular():
    """Guard: passing the angular samples (the output's units) is a 2 pi error."""
    ad = _solver().to_autodiff(rtol=1e-10, atol=1e-12)
    analytic = np.asarray(ad.propagate())
    angular = analytic[:, 0, :].real
    wrong = np.asarray(ad.propagate({"omega0": angular}))
    assert np.abs(wrong - analytic).max() / np.abs(analytic).max() > 0.1


# ---- validation and guards ----------------------------------------------------------


def test_unknown_parameter_and_bad_shapes_are_rejected():
    ad = _solver().to_autodiff()
    with pytest.raises(ValueError, match="Unknown parameter"):
        ad.propagate({"interaction_strength": np.array([1.0, 1.0])})  # typo
    with pytest.raises(ValueError, match="shape"):
        ad.propagate({"detunings": np.zeros(3)})
    with pytest.raises(ValueError, match="omega0"):
        ad.propagate({"omega0": np.zeros((2, 5))})
    with pytest.raises(ValueError, match="pulse_args"):
        ad.propagate({"pulse_args": [{"ampl": 1.0}]})
    with pytest.raises(ValueError, match="mode"):
        ad.value_and_grad(loss, mode="sideways")


def test_counter_propagating_problem_is_rejected():
    d = copy.deepcopy(PROBLEM)
    d["system"]["fields"][0]["counter_propagating"] = True
    d["system"]["fields"][0].pop("factor_doppler_shift", None)
    ps = PropagationSolver.from_json_str(json.dumps(d))
    with pytest.raises(NotImplementedError, match="co-propagating"):
        ps.to_autodiff()


def test_requires_float64():
    ad = _solver().to_autodiff()
    original = jax.config.jax_enable_x64
    try:
        jax.config.update("jax_enable_x64", False)
        with pytest.raises(RuntimeError, match="JAX_ENABLE_X64"):
            ad.propagate()
    finally:
        jax.config.update("jax_enable_x64", original)


def test_compiled_function_is_reused_for_the_same_loss():
    ad = _solver().to_autodiff(rtol=1e-8, atol=1e-10)
    ad.value_and_grad(loss)
    assert len(ad._grad_cache) == 1  # non-vacuous: the program was cached
    fn = next(iter(ad._grad_cache.values()))
    ad.value_and_grad(loss, {"decay_scale": 1.1})  # same loss, same compiled program
    assert len(ad._grad_cache) == 1
    assert next(iter(ad._grad_cache.values())) is fn


def test_compiled_gradient_cache_is_bounded_and_least_recently_used():
    ad = _solver().to_autodiff(rtol=1e-6, atol=1e-8)
    ad.GRAD_CACHE_SIZE = 2
    params = {"decay_scale": 1.0}

    def make(k):
        return lambda omega: k * loss(omega)

    a, b, c = make(1.0), make(2.0), make(3.0)
    ad.value_and_grad(a, params)
    ad.value_and_grad(b, params)
    ad.value_and_grad(a, params)  # a is now more recent than b
    ad.value_and_grad(c, params)  # evicts b
    assert [k[0] for k in ad._grad_cache] == [a, c]
    assert len(ad._grad_cache) == 2


def test_clear_cache_releases_the_loss():
    import gc
    import weakref

    ad = _solver().to_autodiff(rtol=1e-6, atol=1e-8)

    def make():
        return lambda omega: loss(omega)

    fn = make()
    ref = weakref.ref(fn)
    ad.value_and_grad(fn, {"decay_scale": 1.0})
    del fn
    gc.collect()
    assert ref() is not None  # the cache keeps it alive
    ad.clear_cache()
    gc.collect()
    assert ref() is None


def test_evicted_loss_still_gives_the_right_answer():
    ad = _solver().to_autodiff(rtol=1e-6, atol=1e-8)
    ad.GRAD_CACHE_SIZE = 1
    params = {"decay_scale": 1.0}
    v1, g1 = ad.value_and_grad(loss, params)
    ad.value_and_grad(lambda omega: 2.0 * loss(omega), params)  # evicts `loss`
    v2, g2 = ad.value_and_grad(loss, params)  # recompiled
    assert float(v1) == float(v2)
    assert float(g1["decay_scale"]) == float(g2["decay_scale"])


def test_partial_pulse_args_are_merged_consistently():
    gaussian_npi = {"n_pi": 2.0, "centre": 0.0, "fwhm": 1.0}
    assert _merge_pulse_args(gaussian_npi, {"ampl": 0.5}) == {
        "ampl": 0.5,
        "centre": 0.0,
        "fwhm": 1.0,
    }
    sech_width = {"ampl": 1.0, "centre": 0.0, "width": 0.3}
    assert _merge_pulse_args(sech_width, {"fwhm": 1.0}) == {
        "ampl": 1.0,
        "centre": 0.0,
        "fwhm": 1.0,
    }
    assert _merge_pulse_args(sech_width, {"centre": 0.5})["width"] == 0.3


class TestMinimize:
    """``AutodiffPropagation.minimize`` wraps an optax optimiser."""

    @pytest.fixture(scope="class")
    def ad(self):
        return PropagationSolver.from_json_str(json.dumps(PROBLEM)).to_autodiff()

    @staticmethod
    def loss(omegas):
        return jnp.abs(omegas[0, -1, -1]) ** 2

    def test_descends_and_matches_final_value(self, ad):
        optax = pytest.importorskip("optax")
        params = {"interaction_strengths": ad.default_params()["interaction_strengths"]}
        res = ad.minimize(self.loss, params, steps=5, learning_rate=1e-2)
        assert res.n_steps == 5 and len(res.losses) == 6
        assert res.losses[-1] < res.losses[0]
        assert res.value == res.losses[-1]
        assert set(res.params) == set(params)
        v, _ = ad.value_and_grad(self.loss, res.params)
        assert float(v) == pytest.approx(res.value)
        custom = ad.minimize(self.loss, params, steps=2, optimizer=optax.sgd(1e-3))
        assert custom.n_steps == 2

    def test_gtol_stops_immediately(self, ad):
        pytest.importorskip("optax")
        res = ad.minimize(self.loss, steps=10, gtol=1e12)
        assert res.converged and res.n_steps == 0 and len(res.losses) == 1

    def test_bad_steps(self, ad):
        pytest.importorskip("optax")
        with pytest.raises(ValueError, match="steps"):
            ad.minimize(self.loss, steps=0)
