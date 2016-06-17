""" 
Unit tests for the field module.

Thomas Ogden <t@ogden.eu>

"""

import unittest

from maxwellbloch import ob_atom

class TestInit(unittest.TestCase):

    def test_init_default(self):
        """  Test Default Initialise """ 

        ob_atom_00 = ob_atom.OBAtom()

        self.assertEqual(ob_atom_00.num_states, 0.)
        self.assertEqual(ob_atom_00.energies, [])
        self.assertEqual(ob_atom_00.decays, [])
        self.assertEqual(ob_atom_00.fields, [])

def main():
    unittest.main(verbosity=3)

if __name__ == "__main__":
    status = main()
    sys.exit(status)