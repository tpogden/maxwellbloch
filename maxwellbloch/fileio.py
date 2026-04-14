# -*- coding: utf-8 -*-

"""Save and load MBSolve result data to and from files."""

import numpy as np


def save_csv_rabi_freq(mb_solve, field_idx=0, filename=None):
    """Save the complex field Rabi frequency result for an MBSolve field.

    Args:
        mb_solve :
            The MBSolve object
        field_idx :
            The field to save
        filename :
            The path to save to. The default will use the savefile value.
    """

    if not filename:
        filename = mb_solve.savefile + "_rabi_freq_" + str(field_idx) + ".csv"

    np.savetxt(filename, mb_solve.Omegas_zt[field_idx], delimiter=",")


def save_csv_rabi_freq_abs(mb_solve, field_idx=0, filename=None):
    """Save the abs value of the complex field Rabi frequency result for an
        MBSolve field.

    Args:
        mb_solve :
            The MBSolve object
        field_idx :
            The field to save
        filename :
            The path to save to. The default will use the savefile value.
    """

    if not filename:
        filename = mb_solve.savefile + "_rabi_freq_abs_" + str(field_idx) + ".csv"

    np.savetxt(filename, np.abs(mb_solve.Omegas_zt[field_idx]), delimiter=",")
