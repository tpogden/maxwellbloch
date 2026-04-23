"""
Unit tests for the spectral analysis module.

Thomas Ogden <t@ogden.eu>

"""

import os
import unittest

import numpy as np

from maxwellbloch import mb_solve, spectral, utility

# Absolute path of tests/json directory, so that tests can be called from
# different directories.
JSON_DIR = os.path.abspath(os.path.join(__file__, "../", "json"))


class TestSpectral(unittest.TestCase):
    """Unit tests of the spectral methods."""

    def test_two_spectral(self):

        json_path = os.path.join(JSON_DIR, "mbs-two-spectral.json")
        mbs = mb_solve.MBSolve().from_json(json_path)
        mbs.mbsolve()

        # We're not asserting anything about these calls, apart from that they
        # run without exception.
        freq_list = spectral.freq_list(mbs)
        spectral.rabi_freq(mbs, 0)
        spectral.dispersion(mbs, 0, -1)

        abs = spectral.absorption(mbs, 0, -1)
        hm, r1, r2 = utility.half_max_roots(freq_list, abs)

        # The absorption profile should have a peak at 1.0 and have a FWHM of
        # 1.0 centred at 0.0.
        self.assertAlmostEqual(hm, 0.5, places=1)
        self.assertAlmostEqual(r1, -0.5, places=3)
        self.assertAlmostEqual(r2, 0.5, places=3)

        abs_linear_known = spectral.absorption_two_linear_known(
            freq_list, mbs.interaction_strengths[0], mbs.atom.decays[0]["rate"]
        )

        # Assert that the max of the abs residuals between the absorption
        # profile and the known absorption profile for linear two-level systems
        # is below a tolerance
        self.assertTrue(np.max(np.abs(abs - abs_linear_known)) < 0.05)


class TestVoigtProfile(unittest.TestCase):
    """Unit tests for voigt_two_linear_known (GH#198).

    When thermal_width → 0 the Voigt profile is dominated by the Lorentzian
    and should approach the known linear-absorption profile.
    """

    def test_voigt_approaches_lorentzian_for_small_thermal_width(self):
        """Voigt with very small thermal_width ≈ linear absorption."""
        decay_rate = 1.0
        freq_list = np.linspace(-3.0, 3.0, 301)

        lorentzian = spectral.absorption_two_linear_known(
            freq_list,
            interaction_strength=1.0,
            decay_rate=decay_rate,
        )
        voigt = (
            spectral.voigt_two_linear_known(
                freq_list,
                decay_rate=decay_rate,
                thermal_width=1e-4,
            ).imag
            / 2.0
        )

        # Normalise both to their peaks before comparing shape
        lorentzian_norm = lorentzian / np.max(lorentzian)
        voigt_norm = voigt / np.max(voigt)

        self.assertTrue(np.max(np.abs(voigt_norm - lorentzian_norm)) < 0.01)
