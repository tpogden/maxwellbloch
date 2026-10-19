"""Regression tests: calling PropagationSolver.solve() repeatedly on one object
gives the same answer as a fresh solver, and edits between solves take effect.

Previously the z-stepping left system.H_Omega_list holding the last slice's
interpolated coefficients, so the next solve drove slice 0 with stale data.
"""

import copy
import warnings

import numpy as np
import pytest

from clerq import PropagationSolver

T_WIDTH = 0.001

BASE = {
    "system": {
        "fields": [
            {
                "coupled_levels": [[0, 1]],
                "rabi_freq_t_func": "sech",
                "rabi_freq_t_args": {"n_pi": 2.329, "centre": 0.0, "width": T_WIDTH},
            }
        ],
        "num_states": 2,
        "decays": [{"channels": [[0, 1]], "rate": 1.0}],
    },
    "t_min": -10 * T_WIDTH,
    "t_max": 40 * T_WIDTH,
    "t_steps": 100,
    "z_min": 0.0,
    "z_max": 1.0,
    "z_steps": 10,
    "interaction_strengths": [0.6 / T_WIDTH],
}


def _two_field():
    d = copy.deepcopy(BASE)
    d["system"]["num_states"] = 3
    d["system"]["decays"] = [{"channels": [[0, 1]], "rate": 1.0}]
    d["system"]["fields"] = [
        d["system"]["fields"][0],
        {
            "coupled_levels": [[1, 2]],
            "rabi_freq_t_func": "gaussian",
            "rabi_freq_t_args": {"ampl": 1.0, "fwhm": 2 * T_WIDTH, "centre": 0.0},
            "rabi_freq": 100.0,
        },
    ]
    d["interaction_strengths"] = [0.6 / T_WIDTH, 0.0]
    return d


def _counter_prop():
    d = _two_field()
    d["system"]["fields"][1]["counter_propagating"] = True
    return d


def _velocity_classes():
    d = copy.deepcopy(BASE)
    d["velocity_classes"] = {
        "thermal_delta_min": -50.0,
        "thermal_delta_max": 50.0,
        "thermal_delta_steps": 2,
        "thermal_width": 50.0,
    }
    return d


def _solve(s, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s.solve(recalc=True, progress=False, **kw)
    return s._Omegas_zt.copy(), s.states_zt.copy()


@pytest.mark.parametrize("step", ["ab", "euler"])
@pytest.mark.parametrize(
    "make",
    [lambda: BASE, _two_field, _counter_prop, _velocity_classes],
    ids=["single", "two-field", "counter-prop", "velocity-classes"],
)
def test_solve_twice_identical(make, step):
    cfg = make()
    fresh_O, fresh_S = _solve(PropagationSolver(**copy.deepcopy(cfg)), step=step)
    s = PropagationSolver(**copy.deepcopy(cfg))
    for _ in range(3):
        om, st = _solve(s, step=step)
        np.testing.assert_array_equal(om, fresh_O)
        np.testing.assert_array_equal(st, fresh_S)


def test_h_omega_list_restored_after_solve():
    s = PropagationSolver(**copy.deepcopy(BASE))
    _solve(s)
    assert s.system.H_Omega_list[0][1] is s.system.fields[0].rabi_freq_t_func


def test_changing_t_args_between_solves_takes_effect():
    s = PropagationSolver(**copy.deepcopy(BASE))
    O1, _ = _solve(s)
    s.system.fields[0].rabi_freq_t_args["n_pi_0"] = 1.0
    O2, _ = _solve(s)
    assert not np.allclose(O1[0, 0], O2[0, 0])

    other = copy.deepcopy(BASE)
    other["system"]["fields"][0]["rabi_freq_t_args"]["n_pi"] = 1.0
    ref, _ = _solve(PropagationSolver(**other))
    np.testing.assert_array_equal(O2, ref)


def test_deprecated_mbsolve_alias_repeatable():
    s = PropagationSolver(**copy.deepcopy(BASE))
    ref, _ = _solve(PropagationSolver(**copy.deepcopy(BASE)))
    with pytest.warns(DeprecationWarning):
        s.mbsolve(recalc=True, progress=False)
    with pytest.warns(DeprecationWarning):
        s.mbsolve(recalc=True, progress=False)
    np.testing.assert_array_equal(s._Omegas_zt, ref)
