"""Array extraction and unit-conversion helpers for plot primitives."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from maxwellbloch.mb_solve import MBSolve


def omega_abs(mbs: MBSolve, field_idx: int = 0) -> np.ndarray:
    """Return |Ω(z, t)| for *field_idx*.

    Returns:
        Array of shape (z_steps+1, t_steps+1).
    """
    return np.abs(mbs.Omegas_zt[field_idx])


# TODO: make this a method of MBSolve
def pulse_area_z(mbs: MBSolve, field_idx: int = 0) -> np.ndarray:
    """Return the pulse area ∫|Ω(z,t)|dt / π vs z.

    Returns:
        Array of shape (z_steps+1,) in units of π.
    """
    Omega = omega_abs(mbs, field_idx)
    return np.trapezoid(Omega, mbs.tlist, axis=1) / np.pi


def field_label(mbs: MBSolve, field_idx: int = 0) -> str:
    """Return a display label for a field, falling back to its index."""
    label = mbs.atom.fields[field_idx].label
    return label if label else f"field {field_idx}"
