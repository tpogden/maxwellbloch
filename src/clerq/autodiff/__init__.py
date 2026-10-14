"""Automatic differentiation of field propagation (optional, JAX-based).

Needs the ``jax`` extra: ``pip install 'clerq[jax]'``. Importing this package
without it raises an ImportError that says so; importing ``clerq`` itself never
touches JAX.

The solver reproduces ``PropagationSolver`` and returns exact gradients of a
scalar loss of the propagated field with respect to the problem parameters.
"""

from clerq.autodiff._config import require_float64, require_jax

require_jax()

__all__ = ["require_float64"]
