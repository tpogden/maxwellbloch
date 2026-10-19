"""Array extraction and unit-conversion helpers for plot primitives."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from clerq._deprecation import renamed_arg

if TYPE_CHECKING:
    from clerq.propagation import PropagationSolver


@renamed_arg("mbs", "solver")
def omega_abs(solver: PropagationSolver, field_idx: int = 0) -> np.ndarray:
    """Return |Ω(z, t)| for *field_idx*.

    Returns:
        Array of shape (z_steps+1, t_steps+1).
    """
    return np.abs(solver._Omegas_zt[field_idx])


# TODO: make this a method of PropagationSolver
@renamed_arg("mbs", "solver")
def pulse_area_z(solver: PropagationSolver, field_idx: int = 0) -> np.ndarray:
    """Return the pulse area ∫|Ω(z,t)|dt / π vs z.

    Returns:
        Array of shape (z_steps+1,) in units of π.
    """
    Omega = omega_abs(solver, field_idx)
    return np.trapezoid(Omega, solver.tlist, axis=1) / np.pi


@renamed_arg("mbs", "solver")
def field_label(solver: PropagationSolver, field_idx: int = 0) -> str:
    """Return a display label for a field, falling back to its index."""
    label = solver.system.fields[field_idx].label
    return label if label else f"field {field_idx}"
