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


def _grid(n, seed):
    """Strictly increasing, non-uniform grid of ``n`` times."""
    rng = np.random.default_rng(seed)
    return np.cumsum(rng.uniform(0.2, 1.8, n)) - 3.0


@pytest.mark.parametrize("n", [4, 5, 6, 17, 400])
def test_knot_derivatives_equal_the_dense_scipy_operator(n):
    """Banded solve == the dense SciPy reference, on non-uniform grids of any length."""
    from clerq.autodiff._drive import not_a_knot_derivative_matrix

    ts = _grid(n, seed=n)
    y = np.random.default_rng(n).normal(size=(3, n))
    got = np.asarray(CubicDrive(ts).knot_derivatives(jnp.asarray(y)))
    ref = y @ not_a_knot_derivative_matrix(ts).T
    assert np.abs(got - ref).max() < 1e-10 * max(1.0, np.abs(ref).max())
    assert np.abs(ref).max() > 0.1  # non-vacuous


@pytest.mark.parametrize("n", [4, 5, 40])
def test_non_uniform_grid_matches_scipy_spline_and_gradient(n):
    ts = _grid(n, seed=100 + n)
    y = np.random.default_rng(n).normal(size=n)
    drive = CubicDrive(ts)
    t = np.random.default_rng(7).uniform(ts[0], ts[-1], 200)
    ref = make_interp_spline(ts, y, k=3, bc_type=None)(t)
    got = np.asarray(jax.vmap(drive(jnp.asarray(y)))(jnp.asarray(t)))
    assert np.abs(got - ref).max() < 1e-10
    t0 = float(t[0])
    g = np.asarray(jax.grad(lambda yy: drive(yy)(t0))(jnp.asarray(y)))
    basis = make_interp_spline(ts, np.eye(n), k=3, bc_type=None)(t0)
    assert np.abs(g - basis).max() < 1e-10


def test_long_grid_is_linear_in_memory():
    """A dense operator for 20 000 samples would be 3.2 GB; this builds instantly."""
    ts = np.linspace(0.0, 1.0, 20_000)
    y = np.sin(40 * ts)
    drive = CubicDrive(ts)
    t = np.array([0.1234, 0.5, 0.98765])
    got = np.asarray(jax.vmap(drive(jnp.asarray(y)))(jnp.asarray(t)))
    assert np.abs(got - np.sin(40 * t)).max() < 1e-8
