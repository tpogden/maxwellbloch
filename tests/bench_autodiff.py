# -*- coding: utf-8 -*-

"""Performance benchmarks for the optional ``clerq.autodiff`` backend.

Each case times the NumPy/QuTiP ``PropagationSolver.solve`` as a baseline, the
compiled ``AutodiffPropagation.propagate`` (forward pass) and the compiled
``value_and_grad`` (reverse and forward mode). Every compiled function is called
once before timing, so the numbers exclude JAX's first-call compilation.

Reverse mode should cost a small constant multiple of one propagation (about 5x),
independent of the number of parameters; forward mode grows with it.

Run with:
    uv run pytest tests/bench_autodiff.py --benchmark-only
"""

import pytest

jax = pytest.importorskip("jax")
pytest.importorskip("diffrax")
jax.config.update("jax_enable_x64", True)

from clerq import PropagationSolver  # noqa: E402

pytestmark = pytest.mark.jax

# ---------------------------------------------------------------------------
# Benchmark configs
# ---------------------------------------------------------------------------

TWO_LEVEL_SECH_2PI = {
    "system": {
        "num_states": 2,
        "fields": [
            {
                "coupled_levels": [[0, 1]],
                "rabi_freq_t_func": "sech",
                "rabi_freq_t_args": {
                    "ampl": 0.8384014365421667,
                    "centre": 0.0,
                    "width": 0.3796628587572578,
                },
            }
        ],
    },
    "t_min": -2.0,
    "t_max": 10.0,
    "t_steps": 120,
    "z_min": 0.0,
    "z_max": 1.0,
    "z_steps": 30,
    "interaction_strengths": [10.0],
}

LAMBDA_EIT = {
    "system": {
        "num_states": 3,
        "decays": [{"channels": [[0, 1], [1, 2]], "rate": 1.0}],
        "fields": [
            {
                "coupled_levels": [[0, 1]],
                "rabi_freq_t_func": "gaussian",
                "rabi_freq_t_args": {"ampl": 0.1, "centre": 0.0, "fwhm": 0.1},
            },
            {
                "coupled_levels": [[1, 2]],
                "rabi_freq_t_func": "ramp_offon",
                "rabi_freq_t_args": {
                    "ampl": 10.0,
                    "fwhm": 0.05,
                    "off": 0.3,
                    "on": 0.7,
                },
            },
        ],
    },
    "t_min": 0.0,
    "t_max": 1.0,
    "t_steps": 100,
    "z_min": -0.2,
    "z_max": 1.2,
    "z_steps": 20,
    "z_steps_inner": 2,
    "interaction_strengths": [1000.0, 1000.0],
}

TWO_LEVEL_VELOCITY_CLASSES = {
    "system": {
        "num_states": 2,
        "decays": [{"channels": [[0, 1]], "rate": 1.0}],
        "fields": [
            {
                "coupled_levels": [[0, 1]],
                "rabi_freq": 1.0e-3,
                "rabi_freq_t_func": "gaussian",
                "rabi_freq_t_args": {"ampl": 1.0, "centre": 0.0, "fwhm": 1.0},
            }
        ],
    },
    "t_min": -2.0,
    "t_max": 10.0,
    "t_steps": 1000,
    "z_min": -0.2,
    "z_max": 1.2,
    "z_steps": 50,
    "z_steps_inner": 1,
    "interaction_strengths": [1.0],
    "velocity_classes": {
        "thermal_delta_min": -5.0,
        "thermal_delta_max": 5.0,
        "thermal_delta_steps": 4,
        "thermal_delta_inner_min": 0.0,
        "thermal_delta_inner_max": 0.0,
        "thermal_delta_inner_steps": 0,
        "thermal_width": 1.0,
    },
}

CASES = {
    "two_level_sech_2pi": TWO_LEVEL_SECH_2PI,
    "lambda_eit": LAMBDA_EIT,
    "two_level_velocity_classes": TWO_LEVEL_VELOCITY_CLASSES,
}


def loss(omegas_zt):
    """Output pulse energy: a scalar function of the whole propagated field."""
    return (abs(omegas_zt[:, -1, :]) ** 2).sum()


@pytest.fixture(scope="module", params=list(CASES))
def case(request):
    """A solver and its autodiff object, with every compiled function warmed up."""
    solver = PropagationSolver(**CASES[request.param])
    ad = solver.to_autodiff()
    jax.block_until_ready(ad.propagate())
    jax.block_until_ready(ad.value_and_grad(loss))
    jax.block_until_ready(ad.value_and_grad(loss, mode="forward"))
    return solver, ad


# ---------------------------------------------------------------------------
# Benchmarks
# ---------------------------------------------------------------------------


def test_bench_solve_baseline(benchmark, case):
    """The NumPy/QuTiP solver on the same problem, for comparison."""
    solver, _ = case
    benchmark(solver.solve, recalc=True)


def test_bench_autodiff_propagate(benchmark, case):
    """Compiled forward pass (the autodiff counterpart of ``solve``)."""
    _, ad = case
    params = ad.default_params()
    benchmark(lambda: jax.block_until_ready(ad.propagate(params)))


def test_bench_autodiff_value_and_grad_reverse(benchmark, case):
    """Compiled reverse-mode value and gradient over all differentiable params."""
    _, ad = case
    benchmark(lambda: jax.block_until_ready(ad.value_and_grad(loss)))


def test_bench_autodiff_value_and_grad_forward(benchmark, case):
    """Compiled forward-mode value and gradient; cost grows with parameter count."""
    _, ad = case
    benchmark(lambda: jax.block_until_ready(ad.value_and_grad(loss, mode="forward")))
