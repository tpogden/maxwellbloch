# -*- coding: utf-8 -*-

""" Save and load MBSolve result data to and from files. """

from qutip import fileio
import numpy as np

def save_csv_rabi_freq(mb_solve, field_idx=0, filename=None):
    """ Save the complex field Rabi frequency result for an MBSolve field.

    Args:
        mb_solve :
            The MBSolve object
        field_idx :
            The field to save
        filename :
            The path to save to. The default will use the savefile value.
    """

    if not filename:
        filename = mb_solve.savefile + "_rabi_freq_" + str(field_idx) + \
            ".csv"

    fileio.file_data_store(filename=filename, numformat="exp", sep=",",
        data=mb_solve.Omegas_zt[field_idx], numtype="complex")

def save_csv_rabi_freq_abs(mb_solve, field_idx=0, filename=None):
    """ Save the abs value of the complex field Rabi frequency result for an
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
        filename = mb_solve.savefile + "_rabi_freq_abs_" + str(field_idx) + \
            ".csv"

    fileio.file_data_store(filename=filename, numformat="exp", sep=",",
        data=np.abs(mb_solve.Omegas_zt[field_idx]), numtype="real")
