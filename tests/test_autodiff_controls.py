"""``clerq.autodiff.controls``: smooth, bounded parametrisations of free-form controls."""

import numpy as np
import pytest

jax = pytest.importorskip("jax")
pytest.importorskip("diffrax")
jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402

from clerq import PropagationSolver  # noqa: E402
from clerq.autodiff import controls  # noqa: E402

pytestmark = pytest.mark.jax

T = np.linspace(-2.0, 10.0, 121)


def test_gaussian_basis_shape_and_default_centres():
    B = controls.gaussian_basis(T, 8)
    assert B.shape == (121, 8)
    assert B.max() <= 1.0 and B.min() >= 0.0
    assert np.argmax(B[:, 0]) == 0 and np.argmax(B[:, -1]) == 120  # centres at the ends


def test_gaussian_basis_can_be_confined_to_a_window():
    B = controls.gaussian_basis(T, 6, hi=3.0)
    assert (
        np.abs(B[T > 3.0 + 5 * (5.0 / 5)]).max() < 1e-4
    )  # nothing beyond the last centre


def test_basis_rejects_bad_arguments():
    with pytest.raises(ValueError, match="at least 2"):
        controls.gaussian_basis(T, 1)
    with pytest.raises(ValueError, match="greater than"):
        controls.gaussian_basis(T, 4, lo=3.0, hi=3.0)
    with pytest.raises(ValueError, match="width must be positive"):
        controls.gaussian_basis(T, 4, width=0.0)


def test_bounded_and_unbounded_reject_an_empty_range():
    with pytest.raises(ValueError, match="greater than lower"):
        controls.bounded(0.0, 5.0, 5.0)
    with pytest.raises(ValueError, match="greater than lower"):
        controls.unbounded(1.0, 5.0, 1.0)


def test_fit_reproduces_a_smooth_waveform():
    target = np.exp(-0.5 * ((T - 4.0) / 1.5) ** 2)
    B = controls.gaussian_basis(T, 14)
    fit = B @ controls.fit_coefficients(B, target)
    assert np.abs(fit - target).max() < 5e-3


def test_basis_expansions_are_smooth_unlike_free_samples():
    """The point of the helper: no near-Nyquist content, whatever the coefficients."""

    def roughness(x):
        return np.sum(np.diff(x, 2) ** 2) / np.sum(x**2)

    rng = np.random.default_rng(0)
    B = controls.gaussian_basis(T, 12)
    smooth = roughness(B @ rng.normal(size=12))
    free = roughness(rng.normal(size=len(T)))
    assert smooth < 1e-3 * free


def test_bounded_and_unbounded_are_inverses_inside_the_range():
    y = np.linspace(0.3, 4.7, 9)
    x = controls.unbounded(y, 0.0, 5.0, margin=0.0)
    np.testing.assert_allclose(controls.bounded(x, 0.0, 5.0), y, rtol=1e-12)


def test_unbounded_clips_so_the_logit_is_finite():
    x = controls.unbounded(np.array([0.0, 5.0, 7.0, -1.0]), 0.0, 5.0)
    assert np.all(np.isfinite(x))
    assert controls.bounded(x, 0.0, 5.0).min() > 0.0
    assert controls.bounded(x, 0.0, 5.0).max() < 5.0
    with pytest.raises(ValueError, match="margin"):
        controls.unbounded(1.0, 0.0, 5.0, margin=0.6)


def test_bounded_has_a_gradient_everywhere():
    g = jax.grad(lambda x: controls.bounded(x, 0.0, 5.0))(3.0)
    assert float(g) > 0.0


def test_coefficients_are_differentiable_parameters_of_omega0():
    """Gradient of a propagated-field loss with respect to basis coefficients."""
    ps = PropagationSolver.from_json_str(
        """{"system": {"num_states": 2, "decays": [{"rate": 1.0, "channels": [[0, 1]]}],
            "fields": [{"coupled_levels": [[0, 1]], "rabi_freq": 0.1,
                        "rabi_freq_t_func": "gaussian",
                        "rabi_freq_t_args": {"ampl": 1.0, "centre": 1.5, "fwhm": 1.0}}]},
            "t_min": 0.0, "t_max": 3.0, "t_steps": 30, "z_steps": 3, "z_steps_inner": 2,
            "interaction_strengths": [2.0]}"""
    )
    ad = ps.to_autodiff(n_steps=1507)
    t = ad.tlist
    B = controls.gaussian_basis(t, 6)
    a0 = controls.fit_coefficients(
        B, controls.unbounded(0.1 * np.exp(-0.5 * ((t - 1.5) / 0.4) ** 2), 0.0, 0.3)
    )

    def loss_a(a):
        om = controls.bounded(B @ a, 0.0, 0.3)[None]
        out = ad.propagate({"omega0": om})
        return jnp.sum(jnp.abs(out[0, -1]) ** 2)

    g = np.asarray(jax.grad(loss_a)(jnp.asarray(a0)))
    fd = np.empty(6)
    for i in range(6):
        e = np.zeros(6)
        e[i] = 1e-5
        fd[i] = (float(loss_a(a0 + e)) - float(loss_a(a0 - e))) / 2e-5
    assert np.abs(g).max() > 1e-8
    assert np.linalg.norm(g - fd) / np.linalg.norm(fd) < 1e-5
