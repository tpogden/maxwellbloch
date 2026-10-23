"""A hyperfine Rb87 system (11 Zeeman sublevels, Clebsch-Gordan factors) through autodiff."""

import json

import numpy as np
import pytest

jax = pytest.importorskip("jax")
pytest.importorskip("diffrax")
jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402

from clerq import PropagationSolver, hyperfine  # noqa: E402

pytestmark = pytest.mark.jax

TIGHT = {"atol": 1e-12, "rtol": 1e-10, "method": "adams", "nsteps": 1_000_000}
N_FIXED = 1507  # not a multiple of the output intervals


def _atom():
    atom1e = hyperfine.Atom1e(element="Rb", isotope="87")
    atom1e.add_F_level(hyperfine.LevelF(I=1.5, J=0.5, F=1))  # 5s1/2 F=1
    atom1e.add_F_level(hyperfine.LevelF(I=1.5, J=0.5, F=2))  # 5s1/2 F=2 (decay target)
    atom1e.add_F_level(hyperfine.LevelF(I=1.5, J=0.5, F=1))  # 5p1/2 F=1
    return atom1e


def _problem(detuning=0.0, g=14.0, decay_rate=1.0):
    atom = _atom()
    return {
        "system": {
            "num_states": atom.get_num_mF_levels(),
            "energies": atom.get_energies(),
            "decays": [
                {
                    "rate": decay_rate,
                    "channels": atom.get_coupled_levels(
                        F_level_idxs_a=(0, 1), F_level_idxs_b=(2,)
                    ),
                    "factors": atom.get_decay_factors(
                        F_level_idxs_a=(0, 1), F_level_idxs_b=(2,)
                    ).tolist(),
                }
            ],
            "fields": [
                {
                    "coupled_levels": atom.get_coupled_levels(
                        F_level_idxs_a=(0,), F_level_idxs_b=(2,)
                    ),
                    "factors": atom.get_clebsch_hf_factors(
                        F_level_idxs_a=(0,), F_level_idxs_b=(2,), q=1
                    ).tolist(),
                    "detuning": detuning,
                    "rabi_freq": 1e-3,
                    "rabi_freq_t_func": "gaussian",
                    "rabi_freq_t_args": {"ampl": 1.0, "centre": 0.0, "fwhm": 1.0},
                }
            ],
            "initial_state": [0.5 / 3.0] * 3 + [0.5 / 5.0] * 5 + [0.0] * 3,
        },
        "t_min": -2.0,
        "t_max": 4.0,
        "t_steps": 40,
        "z_min": 0.0,
        "z_max": 1.0,
        "z_steps": 4,
        "z_steps_inner": 2,
        "interaction_strengths": [g],
    }


def _solver(**kw):
    ps = PropagationSolver.from_json_str(json.dumps(_problem(**kw)))
    ps.build_opts(TIGHT)
    return ps


_WT = np.exp(-0.5 * (np.linspace(-2.0, 4.0, 41) / 1.5) ** 2)


def loss(omega_zt):
    return jnp.sum(jnp.abs(omega_zt[0, -1]) ** 2 * _WT)


def test_hyperfine_system_matches_solve():
    ps = _solver(detuning=0.3)
    ps.solve(recalc=True, progress=False)
    ref = ps.rabi_freqs_zt
    ad = ps.to_autodiff(rtol=1e-10, atol=1e-12)
    got = np.asarray(ad.propagate())
    scale = np.abs(ref).max()
    np.testing.assert_allclose(got, ref, atol=1e-6 * scale, rtol=0)
    # the medium changes the field (the test is not vacuous)
    assert np.abs(ref[0, -1] - ref[0, 0]).max() > 1e-2 * scale


def test_hyperfine_system_has_one_decay_with_many_channels():
    ad = _solver().to_autodiff()
    ops = ad.setup.ops
    assert ops.num_states == 11
    assert ops.decays is not None and len(ops.decays.rates) == 1
    assert len(ops.decays.a) > 10  # one channel per (lower, upper) sublevel pair


def test_hyperfine_gradients_match_finite_differences():
    """Gradients with respect to g, the detuning, the decay rate and the pulse width."""
    ad = _solver(detuning=0.3).to_autodiff(n_steps=N_FIXED)
    x0 = np.array([14.0, 0.3, 1.0, 1.0])

    def params(x):
        return {
            "interaction_strengths": x[0:1],
            "detunings": x[1:2],
            "decay_rates": x[2:3],
            "pulse_args": [{"fwhm": x[3]}],
        }

    def f(x):
        return float(loss(ad.propagate(params(x))))

    _, g = ad.value_and_grad(loss, params(jnp.asarray(x0)))
    g = np.array(
        [
            g["interaction_strengths"][0],
            g["detunings"][0],
            g["decay_rates"][0],
            g["pulse_args"][0]["fwhm"],
        ]
    )
    fd = np.empty(4)
    for i in range(4):
        e = np.zeros(4)
        e[i] = 1e-5 * max(1.0, abs(x0[i]))
        fd[i] = (f(x0 + e) - f(x0 - e)) / (2 * e[i])
    assert np.all(np.abs(g) > 1e-8)
    assert np.linalg.norm(g - fd) / np.linalg.norm(fd) < 1e-5
