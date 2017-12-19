"""
Unit tests for the spectral analysis module.

Thomas Ogden <t@ogden.eu>

"""

import os

import unittest

from maxwellbloch import mb_solve, spectral

# Absolute path of tests/json directory, so that tests can be called from
# different directories.
JSON_DIR = os.path.abspath(os.path.join(__file__, '../', 'json'))

class TestSpectral(unittest.TestCase):
    """ Unit tests of the spectral methods.

        Note: The real test of the spectral methods is comparison with a
        two-level linear system, as we know the analytic lineshapes. A good
        test might be to compare these lineshapes, however to get good
        agreement a lot of timesteps are needed which makes the test too slow.

        See Appendix C in the notebooks-maxwellbloch repo.
    """

    def test_spectral_twolevel(self):
        """ Check the spectral methods for exceptions. """

        json_path = os.path.join(JSON_DIR, "mb_solve_01.json")
        mb_solve_00 = mb_solve.MBSolve().from_json(json_path)

        mb_solve_00.mbsolve()

        freq_list = spectral.freq_list(mb_solve_00)

        rabi_freq_fft = spectral.rabi_freq(mb_solve_00, 0)

        abs = spectral.absorption(mb_solve_00, 0, -1)
        dis = spectral.dispersion(mb_solve_00, 0, -1)
