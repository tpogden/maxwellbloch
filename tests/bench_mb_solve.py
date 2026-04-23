# -*- coding: utf-8 -*-

"""Performance benchmarks for the Maxwell-Bloch solver.

Run with:
    uv run pytest maxwellbloch/tests/bench_mb_solve.py --benchmark-only

Or together with regular tests:
    uv run pytest --benchmark-disable        # skip benchmarks
    uv run pytest --benchmark-enable         # include benchmarks
"""

from maxwellbloch import mb_solve

# ---------------------------------------------------------------------------
# Benchmark configs — inline dicts mirror the corresponding doc examples
# ---------------------------------------------------------------------------

TWO_LEVEL_SECH_2PI = {
    "atom": {
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
    "atom": {
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
                "rabi_freq_t_func": "square",
                "rabi_freq_t_args": {"ampl": 10.0, "on": 0.0, "off": 1.0},
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
    "atom": {
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


# ---------------------------------------------------------------------------
# Benchmarks
# ---------------------------------------------------------------------------


def test_bench_two_level_sech_2pi(benchmark):
    """Two-level system, sech 2π pulse, no velocity classes.

    The canonical single-field propagation benchmark.
    """
    mbs = mb_solve.MBSolve(**TWO_LEVEL_SECH_2PI)
    benchmark(mbs.mbsolve, recalc=True)


def test_bench_lambda_eit(benchmark):
    """Lambda (Λ) three-level EIT system, two fields, decays.

    Tests the more complex Hamiltonian construction and two-field propagation.
    """
    mbs = mb_solve.MBSolve(**LAMBDA_EIT)
    benchmark(mbs.mbsolve, recalc=True)


def test_bench_two_level_velocity_classes(benchmark):
    """Two-level system with Doppler broadening via velocity classes.

    The most expensive path: the master equation is solved independently
    for each velocity class and results are averaged.
    """
    mbs = mb_solve.MBSolve(**TWO_LEVEL_VELOCITY_CLASSES)
    benchmark(mbs.mbsolve, recalc=True)
