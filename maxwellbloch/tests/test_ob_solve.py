# -*- coding: utf-8 -*-

""" 
Unit tests for the OBSolve class. 

Thomas Ogden <t@ogden.eu>

"""

import unittest

from maxwellbloch import ob_solve

class TestJSON(unittest.TestCase):

    json_str_02 = ('{'
                   '  "ob_atom": {'
                   '    "decays": ['
                   '      { "channels": [[0,1], [1,2]], '
                   '        "rate": 1.0'
                   '      }'
                   '    ],'
                   '    "energies": [],'
                   '    "fields": ['
                   '      {'
                   '        "coupled_levels": ['
                   '          [0, 1]'
                   '        ],'
                   '        "detuning": 0.0,'
                   '        "detuning_positive": true,'
                   '        "label": "probe",'
                   '        "rabi_freq": 5.0,'
                   '        "rabi_freq_t_args": {},'
                   '        "rabi_freq_t_func": null'
                   '      },'
                   '      {'
                   '        "coupled_levels": ['
                   '          [1, 2]'
                   '        ],'
                   '        "detuning": 0.0,'
                   '        "detuning_positive": false,'
                   '        "label": "coupling",'
                   '        "rabi_freq": 10.0,'
                   '        "rabi_freq_t_args": {},'
                   '        "rabi_freq_t_func": null'
                   '      }'
                   '    ],'
                   '    "num_states": 3'
                   '  },'
                   '  "t_min": 0.0,'
                   '  "t_max": 1.0,'
                   '  "t_steps": 100,'
                   '  "method": "mesolve",'
                   '  "opts": {}'
                   '}')

    def test_to_from_json_str_00(self):

        ob_solve_00 = ob_solve.OBSolve()
        ob_solve_01 = ob_solve.OBSolve.from_json_str(ob_solve_00.to_json_str())

        self.assertEqual(ob_solve_00.to_json_str.__repr__(),
                         ob_solve_01.to_json_str.__repr__())

    def test_from_json_str(self):

        ob_solve_02 = ob_solve.OBSolve().from_json_str(self.json_str_02)

        self.assertEqual(ob_solve_02.t_min, 0.0)
        self.assertEqual(ob_solve_02.t_max, 1.0)
        self.assertEqual(ob_solve_02.t_steps, 100)
        self.assertEqual(ob_solve_02.method, "mesolve")

    def test_to_from_json_str_02(self):

        ob_solve_02 = ob_solve.OBSolve().from_json_str(self.json_str_02)
        ob_solve_03 = ob_solve.OBSolve.from_json_str(ob_solve_02.to_json_str())

        self.maxDiff = None
        self.assertEqual(ob_solve_02.to_json_str.__repr__(),
                         ob_solve_03.to_json_str.__repr__())

    def test_to_from_json(self):

        import os

        filepath = "test_ob_solve_02.json"

        ob_solve_02 = ob_solve.OBSolve().from_json_str(self.json_str_02)

        ob_solve_02.to_json(filepath)

        ob_solve_03 = ob_solve.OBSolve().from_json(filepath)
        os.remove(filepath)

        self.maxDiff = None
        self.assertEqual(ob_solve_02.to_json_str(),
                         ob_solve_03.to_json_str())

def main():

    unittest.main(verbosity=3)

if __name__ == "__main__":

    status = main()
    sys.exit(status)