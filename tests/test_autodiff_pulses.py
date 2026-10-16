"""Every JAX pulse port equals the clerq pulse function it re-implements.

``clerq.autodiff._pulses`` re-implements ``clerq.t_funcs`` so the pulses are traceable
by JAX. These tests evaluate each port against the original on a dense time grid
(including the exact edge times), so a change to either side fails here instead of
silently diverging in the first slice only.
"""

import numpy as np
import pytest

jax = pytest.importorskip("jax")
pytest.importorskip("diffrax")
jax.config.update("jax_enable_x64", True)

from clerq import t_funcs  # noqa: E402
from clerq.autodiff._pulses import (  # noqa: E402
    NONDIFFERENTIABLE_ARGS,
    SUPPORTED,
    jump_times,
    pulse_function,
    pulse_name,
    strip_index,
)

pytestmark = pytest.mark.jax

GAUSSIAN = {"centre": 0.4, "fwhm": 0.9}
SECH = {"centre": -0.2}
RAMP = {"ampl": 1.3, "fwhm": 0.35, "on": -0.5, "off": 1.6}

# (name, arguments, extra times that must be hit exactly)
CASES = {
    "gaussian-ampl": ("gaussian", {**GAUSSIAN, "ampl": 1.7}, []),
    "gaussian-n_pi": ("gaussian", {**GAUSSIAN, "n_pi": 2.0}, []),
    "sech-width-ampl": ("sech", {**SECH, "width": 0.3, "ampl": 0.8}, []),
    "sech-fwhm-ampl": ("sech", {**SECH, "fwhm": 0.5, "ampl": 0.8}, []),
    "sech-width-n_pi": ("sech", {**SECH, "width": 0.3, "n_pi": 1.5}, []),
    "sech-fwhm-n_pi": ("sech", {**SECH, "fwhm": 0.5, "n_pi": 1.5}, []),
    "square": ("square", {"ampl": 2.1, "on": -0.5, "off": 1.6}, [-0.5, 1.6]),
    "ramp_on": ("ramp_on", {k: RAMP[k] for k in ("ampl", "fwhm", "on")}, [-0.5]),
    "ramp_off": ("ramp_off", {k: RAMP[k] for k in ("ampl", "fwhm", "off")}, [1.6]),
    "ramp_onoff": ("ramp_onoff", dict(RAMP), [-0.5, 1.6]),
    "ramp_offon": ("ramp_offon", dict(RAMP), [-0.5, 1.6]),
}

T = np.concatenate([np.linspace(-3.0, 4.0, 701), [-0.5, 1.6, 0.4, -0.2]])


def test_every_supported_pulse_is_covered():
    assert {name for name, _, _ in CASES.values()} == set(SUPPORTED)


@pytest.mark.parametrize("case", list(CASES))
def test_port_equals_t_func(case):
    name, args, _ = CASES[case]
    reference = getattr(t_funcs, name)(0)(T, {k + "_0": v for k, v in args.items()})
    got = np.asarray(jax.vmap(pulse_function(name, args))(T))
    assert np.abs(got - reference).max() < 1e-13 * max(1.0, np.abs(reference).max())
    assert np.abs(reference).max() > 0.1  # non-vacuous: the pulse is not ~0 on the grid


@pytest.mark.parametrize("case", [c for c in CASES if CASES[c][2]])
def test_edges_are_handled_exactly_as_t_func_does(case):
    """The comparisons at t == on / off (inclusive for square, the ramps' joins)."""
    name, args, edges = CASES[case]
    ref = getattr(t_funcs, name)(0)
    f = pulse_function(name, args)
    for t in edges:
        assert float(f(t)) == pytest.approx(
            float(ref(t, {k + "_0": v for k, v in args.items()})), abs=1e-14
        )


def test_strip_index_and_pulse_name_round_trip():
    f = t_funcs.ramp_onoff(7)
    assert pulse_name(f) == "ramp_onoff"
    assert strip_index({"ampl_7": 1.0, "n_pi_7": 2.0, "on_7": 0.0}) == {
        "ampl": 1.0,
        "n_pi": 2.0,
        "on": 0.0,
    }


def test_nondifferentiable_arguments_are_exactly_the_square_edges():
    """Only a discontinuous pulse's edge times are excluded; ramps stay differentiable."""
    assert NONDIFFERENTIABLE_ARGS == {"square": ("on", "off")}
    assert jump_times("square", {"ampl": 1.0, "on": 0.2, "off": 0.9}) == (0.2, 0.9)
    assert jump_times("ramp_onoff", RAMP) == ()
