""" Unit tests for the file io module. """

import os

import unittest

from maxwellbloch import mb_solve, fileio

import numpy as np
from qutip import fileio as qufileio

# Absolute path of tests/json directory, so that tests can be called from
# different directories.
JSON_DIR = os.path.abspath(os.path.join(__file__, '../', 'json'))

class TestSaveCSVRabiFreq(unittest.TestCase):

    def tearDown(self):
        os.remove('rabi_freq.csv')
        os.remove('rabi_freq_abs.csv')

    def test_save_load(self):

        json_path = os.path.join(JSON_DIR, 'mb_solve_01.json')
        mb_solve_00 = mb_solve.MBSolve().from_json(json_path)

        # Don't solve, just check initial state

        fileio.save_csv_rabi_freq(
            mb_solve_00, field_idx=0, filename='rabi_freq.csv')
        fileio.save_csv_rabi_freq_abs(
            mb_solve_00, field_idx=0, filename='rabi_freq_abs.csv')

        test_rf = qufileio.file_data_read(filename='rabi_freq.csv')
        test_rf_abs = qufileio.file_data_read(filename='rabi_freq_abs.csv')

        self.assertTrue(np.allclose(test_rf, mb_solve_00.Omegas_zt[0]))
        self.assertTrue(np.allclose(test_rf_abs, mb_solve_00.Omegas_zt[0]))
