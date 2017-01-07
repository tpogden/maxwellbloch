""" 
Unit tests for the field module.

Thomas Ogden <t@ogden.eu>

"""

import unittest

from maxwellbloch import ob_atom

json_str_02 = ('{'
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
               '        "ampl_1": 1.0,'
               '        "off_1": 0.7,'
               '        "on_1": 0.3'
               '      },'
               '      "rabi_freq_t_func": "square_1"'
               '    }'
               '  ],'
               '  "num_states": 3'
               '}')

class TestInit(unittest.TestCase):

    ob_atom_02 = ob_atom.OBAtom.from_json_str(json_str_02)

    def test_init_default(self):
        """  Test Default Initialise """ 

        ob_atom_00 = ob_atom.OBAtom()

        self.assertEqual(ob_atom_00.num_states, 1)
        self.assertEqual(ob_atom_00.energies, [])
        self.assertEqual(ob_atom_00.decays, [])
        self.assertEqual(ob_atom_00.fields, [])

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

        import os

        filepath = "test_ob_atom_02.json"

        self.ob_atom_02.to_json(filepath)

        ob_atom_03 = ob_atom.OBAtom().from_json(filepath)
        os.remove(filepath)

        self.maxDiff = None
        self.assertEqual(self.ob_atom_02.to_json_str(),
                         ob_atom_03.to_json_str())

class TestGetFieldSumCoherence(unittest.TestCase):

    @unittest.skip("states_t() is not initialised so this doesn't work.") 
    def test_default(self):

        ob_atom_00 = ob_atom.OBAtom()

        with self.assertRaises(IndexError) as context:
            ob_atom_00.get_field_sum_coherence(0)

        self.assertTrue("list index out of range" in str(context.exception))
    
    @unittest.skip("states_t() is not initialised so this doesn't work.")    
    def test_initial_condition(self):

        ob_atom_02 = ob_atom.OBAtom.from_json_str(json_str_02)        

        self.assertEqual(ob_atom_02.get_field_sum_coherence(0)[0], 0j)

        self.assertEqual(ob_atom_02.get_field_sum_coherence(1)[0], 0j)

def main():
    unittest.main(verbosity=3)

if __name__ == "__main__":
    status = main()
    sys.exit(status)