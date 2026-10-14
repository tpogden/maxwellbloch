"""The sampled-field interpolant of clerq.autodiff is clerq's spline, and differentiable."""

import numpy as np
import pytest

jax = pytest.importorskip("jax")
pytest.importorskip("diffrax")
jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402
from scipy.interpolate import make_interp_spline  # noqa: E402

from clerq.autodiff._drive import CubicDrive  # noqa: E402

pytestmark = pytest.mark.jax

TS = np.linspace(-2.0, 3.0, 61)


def _samples(seed=0, shape=()):
    rng = np.random.default_rng(seed)
    base = np.exp(-4 * np.log(2) * TS**2) + 0.3 * np.sin(2.0 * TS)
    return base + 0.02 * rng.normal(size=shape + TS.shape)


def test_matches_scipy_not_a_knot_spline():
    y = _samples()
    f = CubicDrive(TS)(jnp.asarray(y))
    t = np.random.default_rng(1).uniform(TS[0], TS[-1], 300)
    ref = make_interp_spline(TS, y, k=3, bc_type=None)(t)
    got = np.asarray(jax.vmap(f)(jnp.asarray(t)))
    assert np.abs(got - ref).max() < 1e-11
    assert np.abs(np.asarray(jax.vmap(f)(jnp.asarray(TS))) - y).max() < 1e-12


def test_leading_axes_are_interpolated_together():
    y = _samples(seed=2, shape=(3,))
    f = CubicDrive(TS)(jnp.asarray(y))
    t = np.linspace(TS[0], TS[-1], 40)
    got = np.asarray(jax.vmap(f)(jnp.asarray(t)))  # (40, 3)
    for k in range(3):
        ref = make_interp_spline(TS, y[k], k=3, bc_type=None)(t)
        assert np.abs(got[:, k] - ref).max() < 1e-11


def test_gradient_is_the_cardinal_spline_basis():
    """d f(t0) / d y_j equals the spline through the unit vector e_j at t0."""
    t0 = 0.4321
    y = jnp.asarray(_samples(seed=3))
    drive = CubicDrive(TS)
    g = np.asarray(jax.grad(lambda yy: drive(yy)(t0))(y))
    basis = make_interp_spline(TS, np.eye(len(TS)), k=3, bc_type=None)(t0)
    assert np.abs(g - basis).max() < 1e-11


def test_differs_from_linear_interpolation():
    """Guard: matching clerq is a statement about the cubic spline specifically."""
    y = _samples(seed=4)
    t = np.linspace(TS[0], TS[-1], 200)
    cubic = np.asarray(jax.vmap(CubicDrive(TS)(jnp.asarray(y)))(jnp.asarray(t)))
    assert np.abs(cubic - np.interp(t, TS, y)).max() > 1e-4


@pytest.mark.parametrize(
    "ts", [np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 1.0, 2.0, 3.0])]
)
def test_invalid_grids_are_rejected(ts):
    with pytest.raises(ValueError):
        CubicDrive(ts)
