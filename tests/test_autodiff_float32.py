"""``clerq.autodiff`` in float32 (the JAX default) against float64.

The float32 half runs in a subprocess with ``jax_enable_x64`` off, because the precision is a
process-wide JAX setting that the rest of the test-suite turns on.
"""

import json
import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest

jax = pytest.importorskip("jax")
pytest.importorskip("diffrax")
jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402

from clerq import PropagationSolver  # noqa: E402
from clerq.autodiff._config import (  # noqa: E402
    complex_dtype,
    default_tolerances,
    float_dtype,
)

pytestmark = pytest.mark.jax

# A lambda medium: a weak probe and a coupling field, with spontaneous decay and ground dephasing.
PROBLEM = {
    "system": {
        "num_states": 3,
        "decays": [
            {"rate": 1.0, "channels": [[0, 1], [2, 1]]},
            {"rate": 0.05, "channels": [[2, 2]]},
        ],
        "fields": [
            {
                "coupled_levels": [[0, 1]],
                "rabi_freq": 0.05,
                "rabi_freq_t_func": "gaussian",
                "rabi_freq_t_args": {"ampl": 1.0, "centre": 1.5, "fwhm": 1.0},
            },
            {
                "coupled_levels": [[1, 2]],
                "rabi_freq": 1.5,
                "rabi_freq_t_func": "gaussian",
                "rabi_freq_t_args": {"ampl": 1.0, "centre": 1.5, "fwhm": 4.0},
            },
        ],
    },
    "t_min": 0.0,
    "t_max": 3.0,
    "t_steps": 60,
    "z_min": 0.0,
    "z_max": 1.0,
    "z_steps": 8,
    "z_steps_inner": 2,
    "interaction_strengths": [4.0, 4.0],
}

PARAMS = {
    "rabi_freq": [0.05, 1.5],
    "interaction_strengths": [4.0, 4.0],
    "decay_rates": [1.0, 0.05],
}

CHILD = textwrap.dedent(
    """
    import json, sys, warnings
    warnings.filterwarnings("error", message=".*dtype.*")  # an explicit float64 request warns
    import numpy as np
    import jax
    import jax.numpy as jnp
    from clerq import PropagationSolver
    from clerq.autodiff._config import default_tolerances, float_dtype

    problem, params = json.loads(sys.argv[1]), json.loads(sys.argv[2])
    ad = PropagationSolver.from_json_str(json.dumps(problem)).to_autodiff()
    field = np.asarray(ad.propagate())
    loss = lambda o: jnp.sum(jnp.abs(o[0, -1]) ** 2)
    p = {k: np.asarray(v) for k, v in params.items()}
    v, g = ad.value_and_grad(loss, p)
    vf, gf = ad.value_and_grad(loss, p, mode="forward")
    flat = lambda d: np.concatenate([np.asarray(d[k]).ravel() for k in sorted(d)]).tolist()
    print("RESULT" + json.dumps({
        "float": str(float_dtype()), "complex": str(field.dtype), "tol": default_tolerances(),
        "field_re": field.real.ravel().tolist(), "field_im": field.imag.ravel().tolist(),
        "value": float(v), "grad": flat(g), "grad_forward": flat(gf),
    }))
    """
)


@pytest.fixture(scope="module")
def float32():
    env = {**os.environ, "JAX_ENABLE_X64": "0", "JAX_PLATFORMS": "cpu"}
    r = subprocess.run(
        [sys.executable, "-c", CHILD, json.dumps(PROBLEM), json.dumps(PARAMS)],
        capture_output=True,
        text=True,
        env=env,
        timeout=600,
    )
    assert r.returncode == 0, r.stderr[-3000:]
    line = next(x for x in r.stdout.splitlines() if x.startswith("RESULT"))
    return json.loads(line[len("RESULT") :])


@pytest.fixture(scope="module")
def float64():
    ad = PropagationSolver.from_json_str(json.dumps(PROBLEM)).to_autodiff()
    field = np.asarray(ad.propagate())
    loss = lambda o: jnp.sum(jnp.abs(o[0, -1]) ** 2)  # noqa: E731
    p = {k: np.asarray(v) for k, v in PARAMS.items()}
    v, g = ad.value_and_grad(loss, p)
    flat = lambda d: np.concatenate([np.asarray(d[k]).ravel() for k in sorted(d)])  # noqa: E731
    return {"field": field, "value": float(v), "grad": flat(g)}


def rel(a, b):
    return float(
        np.linalg.norm(np.asarray(a) - np.asarray(b)) / np.linalg.norm(np.asarray(b))
    )


def test_the_precision_follows_the_jax_config(float32):
    assert float32["float"] == "float32"
    assert float32["complex"] == "complex64"
    assert float_dtype() == jnp.float64 and complex_dtype() == jnp.complex128


def test_default_tolerances_depend_on_the_precision(float32):
    assert tuple(float32["tol"]) == (1e-6, 1e-9)
    assert default_tolerances() == (1e-8, 1e-10)
    assert default_tolerances(rtol=1e-4) == (1e-4, 1e-10)  # an explicit value wins
    assert default_tolerances(atol=1e-5) == (1e-8, 1e-5)


def test_float32_field_matches_float64(float32, float64):
    field = np.array(float32["field_re"]) + 1j * np.array(float32["field_im"])
    ref = float64["field"].ravel()
    assert np.abs(field - ref).max() / np.abs(ref).max() < 1e-3
    assert np.abs(ref).max() > 0.0


def test_float32_gradients_match_float64(float32, float64):
    assert abs(float32["value"] - float64["value"]) < 1e-3 * abs(float64["value"])
    assert rel(float32["grad"], float64["grad"]) < 1e-2
    assert np.abs(float64["grad"]).min() > 0.0  # every parameter matters


def test_float32_forward_and_reverse_mode_agree(float32):
    assert rel(float32["grad_forward"], float32["grad"]) < 1e-2


def test_explicit_tolerances_still_reach_the_solver():
    """Loose tolerances take fewer steps than tight ones, so the values differ at the tolerance level."""
    ps = PropagationSolver.from_json_str(json.dumps(PROBLEM))
    tight = np.asarray(ps.to_autodiff(rtol=1e-10, atol=1e-12).propagate())
    loose = np.asarray(ps.to_autodiff(rtol=1e-3, atol=1e-5).propagate())
    default = np.asarray(ps.to_autodiff().propagate())
    err_loose = np.abs(loose - tight).max()
    err_default = np.abs(default - tight).max()
    assert err_default < err_loose
