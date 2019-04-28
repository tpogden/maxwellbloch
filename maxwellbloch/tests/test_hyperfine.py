""" 
Unit tests for the hyperfine module. 

Thomas Ogden <t@ogden.eu>
"""

import sys
import unittest

import numpy as np

from maxwellbloch import hyperfine


class TestLevelJInit(unittest.TestCase):
    """ Unit tests of the LevelJ.__init__ method. """

    def test_Rb_87_5p12(self):
        """ Test building the 5p_{1/2} level of Rubidium 87. """

        I = 1.5 # Rb87 nuclear spin
        J = 0.5 # 5p_{1/2}
        Rb_87_5p12 = hyperfine.LevelJ(I=I, J=J, energy=0.0, 
            F_energies=[10.0, 20.0])

        # F numbers are in the range |J - I| <= F <= J + I, so in this case
        # F=1 and F=2
        self.assertEqual(len(Rb_87_5p12.F_levels), 2)

        # Each F level should have 2*F+1 sublevels.
        self.assertEqual(len(Rb_87_5p12.F_levels[0].mF_levels), 2*1+1)
        self.assertEqual(len(Rb_87_5p12.F_levels[1].mF_levels), 2*2+1)

class TestLevelFInit(unittest.TestCase):
    """ Unit tests of the LevelF.__init__ method. """

    def test_f_2(self):
        """ Tests that without specified energies, the mF sublevel energies
            are set to zero. """
        
        F = 2
        lf = hyperfine.LevelF(F=F, energy=0.0)
        self.assertEqual(len(lf.mF_levels), 2*F+1)
        for i in lf.mF_levels:
            self.assertEqual(i.energy, 0.0)

    def test_f_2_energies(self):
        """ Tests that the mF sublevel energies are set. """

        F = 2
        mF_energies = [10.0, 20.0, 30.0, 40.0, 50.0]
        lf = hyperfine.LevelF(F=F, energy=0.0, mF_energies=mF_energies)

        self.assertEqual(len(lf.mF_levels), 2*F+1)
        for i, e in enumerate(mF_energies):
            self.assertEqual(e, lf.mF_levels[i].energy)

def main():

    unittest.main(verbosity=3) # Run all

if __name__ == "__main__":
    status = main()
    sys.exit(status)