"""
Unit tests for the macro module.

Thomas Ogden <t@ogden.eu>
"""

import unittest

import numpy as np

from clerq import macro

# Rb87 D2 line
GAMMA = 6.065
LAM = 780.24


class TestIntensity(unittest.TestCase):
    def test_saturation_intensity_rb87_d2(self):
        """Steck gives I_sat = 1.669 mW/cm^2 for the D2 cycling transition."""
        self.assertAlmostEqual(
            macro.saturation_intensity(GAMMA, LAM), 1.669, delta=0.01
        )

    def test_rabi_freq_at_saturation(self):
        """At I = I_sat, Omega = Gamma / sqrt(2)."""
        i_sat = macro.saturation_intensity(GAMMA, LAM)
        self.assertAlmostEqual(
            macro.intensity_to_rabi_freq(i_sat, GAMMA, LAM), GAMMA / np.sqrt(2)
        )

    def test_round_trip(self):
        omega = macro.intensity_to_rabi_freq(3.0, GAMMA, LAM, factor=0.5)
        self.assertAlmostEqual(
            macro.rabi_freq_to_intensity(omega, GAMMA, LAM, factor=0.5), 3.0
        )

    def test_gaussian_beam_peak(self):
        self.assertAlmostEqual(
            macro.gaussian_beam_peak_intensity(1.0, 0.1), 2 / (np.pi * 0.01)
        )


class TestMedium(unittest.TestCase):
    def test_od_matches_interaction_strength(self):
        """Weak-field amplitude loss is 2 N g / Gamma, half the intensity OD."""
        n, length = 1e10, 1.0
        g = macro.interaction_strength(GAMMA, LAM)
        od = macro.resonant_optical_depth(n, length, LAM)
        self.assertAlmostEqual(od, 2 * (2 * n * g * length / GAMMA), places=6)

    def test_thermal_width_rb87_room_temperature(self):
        """The Doppler FWHM of Rb87 D2 at 300 K is ~ 510 MHz."""
        w = macro.thermal_width(300.0, 86.909, LAM)
        self.assertAlmostEqual(2 * np.sqrt(np.log(2)) * w, 511.0, delta=3.0)


if __name__ == "__main__":
    unittest.main()
