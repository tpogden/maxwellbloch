""" 
Unit tests for the hyperfine module. 

Thomas Ogden <t@ogden.eu>
"""

import unittest

import numpy as np

from maxwellbloch import hyperfine

class TestAtom1eAddJLevel(unittest.TestCase):
    """ Unit tests of the Atom1e.add_J_level method. """

    def test_Rb87_5s_5p(self):

        Rb87_5s12 = hyperfine.LevelJ(I=1.5, J=0.5)

        Rb87_5p12 = hyperfine.LevelJ(I=1.5, J=0.5)
        Rb87_5p32 = hyperfine.LevelJ(I=1.5, J=1.5)

        Rb87_5s_5p = hyperfine.Atom1e(element='Rb', isotope='87')

        Rb87_5s_5p.add_J_level(Rb87_5s12)
        Rb87_5s_5p.add_J_level(Rb87_5p12)
        Rb87_5s_5p.add_J_level(Rb87_5p32)

        self.assertEqual(Rb87_5s_5p.get_num_mF_levels(), 32)

        map = [0]*8
        map.extend([1]*8)
        map.extend([2] * 16)

        self.assertEqual(Rb87_5s_5p.get_J_level_idx_map(), map)

        self.assertEqual(len(Rb87_5s_5p.get_coupled_levels(0, 1)), 8*8)
        self.assertEqual(len(Rb87_5s_5p.get_coupled_levels(0, 2)), 8*16)
        
class TestLevelNLInit(unittest.TestCase):
    """ Unit tests of the LevelNL.__init__ method. """

    def test_Rb_87_5s(self):
        """ Test the structure built for the Rubidium 87 5s level. """

        Rb_87_5s = hyperfine.LevelNL(n=5, I=1.5, L=1, S=0.5)
        self.assertEqual(len(Rb_87_5s.J_levels), 2)

        # Fine structure
        self.assertEqual(len(Rb_87_5s.J_levels[0].F_levels), 2)
        self.assertEqual(len(Rb_87_5s.J_levels[1].F_levels), 4)

        # Hyperfine structure
        self.assertEqual(len(Rb_87_5s.J_levels[0].F_levels[0].mF_levels), 3)
        self.assertEqual(len(Rb_87_5s.J_levels[0].F_levels[1].mF_levels), 5)

        self.assertEqual(len(Rb_87_5s.J_levels[1].F_levels[0].mF_levels), 1)
        self.assertEqual(len(Rb_87_5s.J_levels[1].F_levels[1].mF_levels), 3)
        self.assertEqual(len(Rb_87_5s.J_levels[1].F_levels[2].mF_levels), 5)
        self.assertEqual(len(Rb_87_5s.J_levels[1].F_levels[3].mF_levels), 7)

class TestLevelJInit(unittest.TestCase):
    """ Unit tests of the LevelJ.__init__ method. """

    def test_not_half_ints(self):
        """ Test that a ValueError is raised if I and J are not both half 
            integer. """

        with self.assertRaises(ValueError):
            hyperfine.LevelJ(I=1.6, J=0.5, energy=0.0)
        with self.assertRaises(ValueError):
            hyperfine.LevelJ(I=1.5, J=0.6, energy=0.0)
         
    def test_no_F_energies(self):
        """ Test with no setting of F_energies or mF_energies (so should all be 
            zero) """

        I = 1.5 # Rb87 nuclear spin
        J = 0.5 # 5p_{1/2}
        Rb_87_5p12 = hyperfine.LevelJ(I=I, J=J, energy=0.0, 
            F_energies=None, mF_energies=None)
        self.assertEqual(len(Rb_87_5p12.F_levels), 2)
        self.assertEqual(len(Rb_87_5p12.F_levels[0].mF_levels), 2*1+1)
        self.assertEqual(len(Rb_87_5p12.F_levels[1].mF_levels), 2*2+1)

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
