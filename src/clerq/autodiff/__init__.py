"""Automatic differentiation of field propagation (optional, JAX-based).

Needs the ``jax`` extra: ``pip install 'clerq[jax]'``. Importing this package
without it raises an ImportError that says so; importing ``clerq`` itself never
touches JAX.

Build the differentiable object from an ordinary solver and differentiate::

    import jax
    jax.config.update("jax_enable_x64", True)  # required: clerq.autodiff is float64

    ad = solver.to_autodiff()             # solver is a clerq.PropagationSolver
    params = ad.default_params()
    omegas = ad.propagate(params)         # same layout as solver.Omegas_zt

    def loss(omegas_zt):                  # any scalar function of the field
        return -abs(omegas_zt[0, -1]).max()

    value, grads = ad.value_and_grad(loss, params)   # exact gradients

Or let an optax optimiser (Adam by default) follow the gradient::

    result = ad.minimize(loss, params, steps=50)
    result.params, result.losses

The result reproduces ``PropagationSolver.solve``; see
:class:`~clerq.autodiff.AutodiffPropagation` for the parameters and caveats.
"""

from clerq.autodiff._config import require_float64, require_jax

require_jax()

from clerq.autodiff._api import AutodiffPropagation, OptimizeResult  # noqa: E402

__all__ = ["AutodiffPropagation", "OptimizeResult", "require_float64"]
