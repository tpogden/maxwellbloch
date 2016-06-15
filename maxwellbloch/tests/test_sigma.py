"""Unit tests for the sigma module.

Thomas Ogden <t@ogden.eu>
"""

import sys
import unittest
import numpy as np
import qutip as qu

from maxwellbloch import sigma

class TestSigma(unittest.TestCase):
    """ Tests for sigma.sigma. """

    def test_sigma_2_0_0(self):
        """ Test |0><0| for a two-level system. """ 

        sigma_2_0_0 = qu.Qobj([[1.,0.],[0.,0.]])
        self.assertEqual(sigma.sigma(2,0,0), sigma_2_0_0)

    def test_sigma_10_9_9(self):
        """ Test |9><9| for a ten-level system. """ 

        sigma_10_9_9 = np.zeros((10,10))
        sigma_10_9_9[9,9] = 1.
        sigma_10_9_9 = qu.Qobj(sigma_10_9_9)
        self.assertEqual(sigma.sigma(10,9,9), sigma_10_9_9)

class TestSigmaN(unittest.TestCase):
    """ Tests for sigma.sigma_N. """ 

    def test_sigma_N_2_0_0_0_1(self):
        """ Test that sigma_N with 1 subsystem returns same as sigma. """

        self.assertEqual(sigma.sigma_N(2,0,0,0,1),sigma.sigma(2,0,0))

    def test_sigma_N_2_0_1_1_2(self):
        """ Test for |0><1| on the 2nd of 2 interacting two-level systems. """

        sigma_2_0_1 =  qu.Qobj([[0.,1.],[0.,0.]])
        sigma_N_2_0_1_1_2 = qu.tensor(qu.identity(2), sigma_2_0_1)

        self.assertEqual(sigma.sigma_N(2,0,1,1,2),sigma_N_2_0_1_1_2)

def main():
    unittest.main(verbosity=3)

if __name__ == "__main__":
    status = main()
    sys.exit(status)