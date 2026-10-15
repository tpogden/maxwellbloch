"""Regression tests from the code review of the autodiff API.

1. Integer literals in a problem definition used to make ``value_and_grad`` with the
   default parameters fail (jax.grad rejects integer inputs).
2. The edge times of a ``square`` pulse are exposed as parameters, but the field's
   derivative with respect to a discontinuity's location is a delta function:
   autodiff silently returned exactly 0 where finite differences give thousands.
"""

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
from clerq.autodiff._api import _overlay  # noqa: E402

pytestmark = pytest.mark.jax

N_FIXED = 2003  # not a multiple of the 40 output intervals


def _problem(func="gaussian", args=None, literal=int):
    """A two-level problem; ``literal`` is the type of the numeric literals used."""
    args = args or {"ampl": 1, "centre": 1, "fwhm": 1}
    return {
        "system": {
            "num_states": 2,
            "decays": [{"rate": 1.0, "channels": [[0, 1]]}],
            "fields": [
                {
                    "coupled_levels": [[0, 1]],
                    "rabi_freq": literal(1),
                    "detuning": literal(0),
                    "rabi_freq_t_func": func,
                    "rabi_freq_t_args": {
                        k: literal(v) if float(v).is_integer() else v
                        for k, v in args.items()
                    },
                }
            ],
        },
        "t_min": 0.0,
        "t_max": 3.0,
        "t_steps": 40,
        "z_steps": 3,
        "z_steps_inner": 1,
        "interaction_strengths": [1.0],
    }


def _ad(problem, **options):
    ps = PropagationSolver.from_json_str(json.dumps(problem))
    return ps.to_autodiff(**options)


def loss(omega_zt):
    return jnp.sum(jnp.abs(omega_zt[:, -1, :]) ** 2)


def _rel(a, b):
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300))


def _fd_gradient(ad, params, fixed=None, h=1e-5):
    """Central differences of loss(propagate(...)) over the flattened ``params``."""
    flat, unravel = ravel_pytree(jax.tree_util.tree_map(jnp.asarray, params))

    def f(x):
        return float(loss(ad.propagate(_overlay(fixed or {}, unravel(x)))))

    x = np.asarray(flat)
    out = np.empty(len(x))
    for i in range(len(x)):
        e = np.zeros(len(x))
        e[i] = h
        out[i] = (f(jnp.asarray(x + e)) - f(jnp.asarray(x - e))) / (2 * h)
    return out


# ---- 1. integer literals ----------------------------------------------------------


def test_default_params_are_floats_whatever_the_literal():
    p = _ad(_problem(literal=int)).default_params()
    leaves = jax.tree_util.tree_leaves(p)
    assert all(isinstance(x, float) or np.asarray(x).dtype.kind == "f" for x in leaves)


def test_value_and_grad_with_default_params_and_integer_literals():
    ad_int = _ad(_problem(literal=int), rtol=1e-9, atol=1e-11)
    ad_flt = _ad(_problem(literal=float), rtol=1e-9, atol=1e-11)
    v_i, g_i = ad_int.value_and_grad(loss)  # used to raise TypeError
    v_f, g_f = ad_flt.value_and_grad(loss)
    assert np.isfinite(float(v_i))
    assert abs(float(v_i) - float(v_f)) < 1e-12 * max(1.0, abs(float(v_f)))
    assert (
        _rel(np.asarray(ravel_pytree(g_i)[0]), np.asarray(ravel_pytree(g_f)[0])) < 1e-9
    )


def test_user_supplied_integer_parameters_are_accepted():
    ad = _ad(_problem(literal=float), rtol=1e-9, atol=1e-11)
    v, g = ad.value_and_grad(loss, {"decay_scale": 1, "detunings": np.array([0])})
    assert np.isfinite(float(v)) and np.isfinite(float(g["decay_scale"]))


# ---- 2. square-pulse edges --------------------------------------------------------

SQUARE = {"ampl": 1.5, "on": 0.6, "off": 2.1}


def test_square_edges_are_held_fixed_by_default_and_the_rest_is_exact():
    ad = _ad(_problem("square", SQUARE, literal=float), n_steps=N_FIXED)
    _, g = ad.value_and_grad(loss)
    assert set(g["pulse_args"][0]) == {"ampl"}  # edges are not differentiated
    # every returned gradient agrees with finite differences at the problem's edges
    params = ad._differentiable_default_params()
    fd = _fd_gradient(ad, params)
    got = np.asarray(ravel_pytree(g)[0])
    assert _rel(got, fd) < 1e-5, (got, fd)
    assert np.abs(got).max() > 1e-3  # non-vacuous


def test_asking_for_the_gradient_of_square_edges_is_an_error_not_a_silent_zero():
    ad = _ad(_problem("square", SQUARE, literal=float), n_steps=N_FIXED)
    for edge in ("on", "off"):
        with pytest.raises(ValueError, match="square"):
            ad.value_and_grad(loss, {"pulse_args": [{edge: SQUARE[edge]}]})

    # and the underlying reason holds: the true derivative is large, not zero
    def f(on):
        return float(loss(ad.propagate({"pulse_args": [{"on": on}]})))

    h = 5e-3
    fd_on = (f(0.6 + h) - f(0.6 - h)) / (2 * h)
    assert abs(fd_on) > 100.0


def test_square_edges_can_be_changed_via_fixed():
    ad = _ad(_problem("square", SQUARE, literal=float), n_steps=N_FIXED)
    fixed = {"pulse_args": [{"on": 0.7, "off": 2.0}]}
    v, g = ad.value_and_grad(loss, fixed=fixed)
    expected = loss(ad.propagate({"pulse_args": [{"on": 0.7, "off": 2.0}]}))
    assert abs(float(v) - float(expected)) < 1e-10 * max(1.0, abs(float(expected)))
    assert abs(float(v) - float(loss(ad.propagate()))) > 1e-4  # the edges did change
    fd = _fd_gradient(ad, ad._differentiable_default_params(), fixed=fixed)
    assert _rel(np.asarray(ravel_pytree(g)[0]), fd) < 1e-5


def test_continuous_ramp_edges_are_differentiable():
    """The restriction is specific to discontinuous pulses: ramp on/off are fine."""
    args = {"ampl": 1.2, "fwhm": 0.3, "on": 0.7, "off": 2.0}
    ad = _ad(_problem("ramp_onoff", args, literal=float), n_steps=N_FIXED)
    params = {"pulse_args": [dict(args)]}
    _, g = ad.value_and_grad(loss, params)  # no error: edges are allowed
    got = np.asarray(ravel_pytree(g)[0])
    fd = _fd_gradient(ad, params)
    assert _rel(got, fd) < 1e-5, (got, fd)
    assert abs(float(g["pulse_args"][0]["on"])) > 1e-6
    assert abs(float(g["pulse_args"][0]["off"])) > 1e-6
