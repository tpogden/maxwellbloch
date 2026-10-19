# -*- coding: utf-8 -*-

"""Save and load PropagationSolver result data to and from files."""

import numpy as np

from clerq._deprecation import renamed_arg


@renamed_arg("mb_solve", "solver")
def save_csv_rabi_freq(solver, field_idx=0, filename=None):
    """Save the complex field Rabi frequency result for a PropagationSolver field.

    Args:
        solver :
            The PropagationSolver object
        field_idx :
            The field to save
        filename :
            The path to save to. The default will use the savefile value.
    """

    if not filename:
        filename = solver.savefile + "_rabi_freq_" + str(field_idx) + ".csv"

    np.savetxt(filename, solver._Omegas_zt[field_idx], delimiter=",")


@renamed_arg("mb_solve", "solver")
def save_csv_rabi_freq_abs(solver, field_idx=0, filename=None):
    """Save the abs value of the complex field Rabi frequency result for an
        PropagationSolver field.

    Args:
        solver :
            The PropagationSolver object
        field_idx :
            The field to save
        filename :
            The path to save to. The default will use the savefile value.
    """

    if not filename:
        filename = solver.savefile + "_rabi_freq_abs_" + str(field_idx) + ".csv"

    np.savetxt(filename, np.abs(solver._Omegas_zt[field_idx]), delimiter=",")
