""" 
Unit tests for the field module.

Thomas Ogden <t@ogden.eu>

"""

import sys
import os

import unittest

from maxwellbloch import ob_atom

# Absolute path of tests/json directory, so that tests can be called from
# different directories.
JSON_DIR = os.path.abspath(os.path.join(__file__, '../', 'json'))

JSON_STR_02 = ('{'
               '  "decays": ['
               '    {'
               '      "channels": ['
               '        ['
               '          0,'
               '          1'
               '        ],'
               '        ['
               '          2,'
               '          1'
               '        ]'
               '      ],'
               '      "rate": 1.0'
               '    },'
               '    {'
               '      "channels": ['
               '        ['
               '          2,'
               '          1'
               '        ]'
               '      ],'
               '      "rate": 2.0'
               '    }'
               '  ],'
               '  "energies": ['
               '    1.0,'
               '    2.0,'
               '    3.0'
               '  ],'
               '  "fields": ['
               '    {'
               '      "coupled_levels": ['
               '        ['
               '         0,'
               '          1'
               '        ]'
               '      ],'
               '      "detuning": 1.0,'
               '      "detuning_positive": true,'
               '      "label": "probe",'
               '      "rabi_freq": 0.001,'
               '      "rabi_freq_t_args": {},'
               '      "rabi_freq_t_func": null'
               '    },'
               '    {'
               '      "coupled_levels": ['
               '        ['
               '          1,'
               '          2'
               '        ]'
               '      ],'
               '      "detuning": 2.0,'
               '      "detuning_positive": false,'
               '      "label": "coupling",'
               '      "rabi_freq": 10.0,'
               '      "rabi_freq_t_args": {'
               '        "ampl": 1.0,'
               '        "off": 0.7,'
               '        "on": 0.3'
               '      },'
               '      "rabi_freq_t_func": "square"'
               '    }'
               '  ],'
               '  "num_states": 3'
               '}')


class TestInit(unittest.TestCase):

    def test_init_default(self):
        """  Test Default Initialise """ 

        ob_atom_00 = ob_atom.OBAtom()

        self.assertEqual(ob_atom_00.num_states, 1)
        self.assertEqual(ob_atom_00.energies, [])
        self.assertEqual(ob_atom_00.decays, [])
        self.assertEqual(ob_atom_00.fields, [])


class TestJSON(unittest.TestCase):

    ob_atom_02 = ob_atom.OBAtom.from_json_str(JSON_STR_02)

    def test_to_from_json_str(self):

        ob_atom_00 = ob_atom.OBAtom()
        ob_atom_01 = ob_atom.OBAtom().from_json_str(ob_atom_00.to_json_str())

        self.assertEqual(ob_atom_00.to_json_str(),
                         ob_atom_01.to_json_str())        

    def test_to_from_json_str_02(self):

        ob_atom_03 = ob_atom.OBAtom().from_json_str(
            self.ob_atom_02.to_json_str())

        self.maxDiff = None
        self.assertEqual(self.ob_atom_02.to_json_str(),
                         ob_atom_03.to_json_str())           

    def test_to_from_json(self):

        filepath = os.path.join(JSON_DIR, "test_ob_atom_02.json")

        self.ob_atom_02.to_json(filepath)

        ob_atom_03 = ob_atom.OBAtom().from_json(filepath)
        os.remove(filepath)

        self.maxDiff = None
        self.assertEqual(self.ob_atom_02.to_json_str(),
                         ob_atom_03.to_json_str())

    def test_from_json_file_02(self):

        json_path = os.path.join(JSON_DIR, "ob_atom_02.json")

        oba = ob_atom.OBAtom().from_json(json_path)

class TestGetFieldSumCoherence(unittest.TestCase):

    @unittest.skip("states_t() is not initialised so this doesn't work.") 
    def test_default(self):

        ob_atom_00 = ob_atom.OBAtom()

        with self.assertRaises(IndexError) as context:
            ob_atom_00.get_field_sum_coherence(0)

        self.assertTrue("list index out of range" in str(context.exception))
    
    @unittest.skip("states_t() is not initialised so this doesn't work.")    
    def test_initial_condition(self):

        ob_atom_02 = ob_atom.OBAtom.from_json_str(JSON_STR_02)        

        self.assertEqual(ob_atom_02.get_field_sum_coherence(0)[0], 0j)

        self.assertEqual(ob_atom_02.get_field_sum_coherence(1)[0], 0j)

def main():
    unittest.main(verbosity=3)

if __name__ == "__main__":
    status = main()
    sys.exit(status)
