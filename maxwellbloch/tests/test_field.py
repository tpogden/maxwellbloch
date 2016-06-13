""" 
Unit tests for the field module.

Thomas Ogden <t@ogden.eu>

"""

import unittest

from maxwellbloch import field

class TestInit(unittest.TestCase):

    def test_init_default(self):
        """  Test Default Initialise """ 

        field_00 = field.Field()

        self.assertEqual(field_00.label, '')
        self.assertEqual(field_00.coupled_levels, [])
        self.assertEqual(field_00.detuning, 0.0)
        self.assertEqual(field_00.detuning_positive, True)
        self.assertEqual(field_00.rabi_freq, 0.0)
        self.assertEqual(field_00.rabi_freq_t_func, None)
        self.assertEqual(field_00.rabi_freq_t_args, {})

def main():
    unittest.main(verbosity=3)

if __name__ == "__main__":
    status = main()
    sys.exit(status)