""" 
Unit tests for the macro module. 
"""

import unittest
from numpy import testing

from maxwellbloch import macro

from numpy import pi

class TestConvertDipoleMomentAtomicToSI(unittest.TestCase):

    def test_unit(self):

        # https://en.wikipedia.org/wiki/Hartree_atomic_units
        C_KNOWN = 8.47835363e-30

        d_si = macro.convert_dipole_moment_atomic_to_si(dipole_moment_ea0=1.0)

        testing.assert_approx_equal(d_si, C_KNOWN, significant=7)

class TestPropagationCoefficientMHzCm2(unittest.TestCase):

    def test_from_thesis(self):

        d = 2.53e-29 # [C m]
        k = 2*pi*1.2578948e6 # [m^-1] 

        # g_si = macro.propagation_coefficient_MHz_cm2(dipole_moment=d, 
        #     wavenumber=k)

        g = macro.propagation_coefficient_MHz_cm2(dipole_moment=d, 
            wavenumber=k) # [2π MHz cm^2]
        self.assertAlmostEqual(g, 4.3115e-9)

class TestInteractionStrengthMHzPerCm(unittest.TestCase):

    def test_from_thesis(self):

        # NOTE: I think my thesis example calculation is out by 10^3!

        # https://steck.us/alkalidata/rubidium85numbers.pdf
        d = 2.53e-29 # [C m]
        k = 2*pi*1.2578948e6 # [m^-1] 

        g_MHz_cm2 = macro.propagation_coefficient_MHz_cm2(dipole_moment=d, 
            wavenumber=k)

        # For T = 200C
        N = 9.26e14 # [cm^-3]

        Ng_MHz_per_cm = macro.interaction_strength_MHz_per_cm(dipole_moment=d, 
            wavenumber=k, number_density=N)

        Gamma = 5.75 # [2pi MHz]
        L = 0.1 # [cm]

        Ng = Ng_MHz_per_cm/Gamma*L # [2π Γ /L]

        self.assertAlmostEqual(Ng, 70)

        pass