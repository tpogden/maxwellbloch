# -*- coding: utf-8 -*-

"""
Unit tests for the OBSolve class.

Thomas Ogden <t@ogden.eu>

"""

import unittest

from maxwellbloch import ob_solve, t_funcs

JSON_STR_02 = (
    '{'
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
    '}'
    )

class TestSetFieldRabiTFunc(unittest.TestCase):
    """ Test setting custom Rabi frequency time functions. """

    def test_set_field_rabi_t_func_1(self):
        """ Test that a custom double pulse Rabi freq time functions can be
            set.
        """

        ob_solve_02 = ob_solve.OBSolve().from_json_str(JSON_STR_02)

        two_pulse_t_func = lambda t, args: (t_funcs.gaussian_1(t, args) +
            t_funcs.gaussian_2(t, args))

        two_pulse_t_args = {"ampl_1": 1.0, "centre_1": 0.0, "fwhm_1": 0.1,
            "ampl_2": 2.0, "centre_2": 0.5, "fwhm_2": 0.1, }

        ob_solve_02.set_field_rabi_freq_t_func(0, two_pulse_t_func)
        ob_solve_02.set_field_rabi_freq_t_args(0, two_pulse_t_args)

        field_0 = ob_solve_02.ob_atom.fields[0]

        self.assertAlmostEqual(field_0.rabi_freq_t_func(0.0,
            field_0.rabi_freq_t_args), 1.0)
        self.assertAlmostEqual(field_0.rabi_freq_t_func(0.5,
            field_0.rabi_freq_t_args), 2.0)
        self.assertAlmostEqual(field_0.rabi_freq_t_func(1.0,
            field_0.rabi_freq_t_args), 0.0)

class TestJSON(unittest.TestCase):

    def test_to_from_json_str_00(self):

        ob_solve_00 = ob_solve.OBSolve()
        ob_solve_01 = ob_solve.OBSolve.from_json_str(ob_solve_00.to_json_str())

        self.assertEqual(ob_solve_00.to_json_str.__repr__(),
                         ob_solve_01.to_json_str.__repr__())

    def test_from_json_str(self):

        ob_solve_02 = ob_solve.OBSolve().from_json_str(JSON_STR_02)

        self.assertEqual(ob_solve_02.t_min, 0.0)
        self.assertEqual(ob_solve_02.t_max, 1.0)
        self.assertEqual(ob_solve_02.t_steps, 100)
        self.assertEqual(ob_solve_02.method, "mesolve")

    def test_to_from_json_str_02(self):

        ob_solve_02 = ob_solve.OBSolve().from_json_str(JSON_STR_02)
        ob_solve_03 = ob_solve.OBSolve.from_json_str(ob_solve_02.to_json_str())

        self.maxDiff = None
        self.assertEqual(ob_solve_02.to_json_str.__repr__(),
                         ob_solve_03.to_json_str.__repr__())

    def test_to_from_json(self):

        import os

        filepath = "test_ob_solve_02.json"

        ob_solve_02 = ob_solve.OBSolve().from_json_str(JSON_STR_02)

        ob_solve_02.to_json(filepath)

        ob_solve_03 = ob_solve.OBSolve().from_json(filepath)
        os.remove(filepath)

        self.maxDiff = None
        self.assertEqual(ob_solve_02.to_json_str(),
                         ob_solve_03.to_json_str())

class TestSaveLoad(unittest.TestCase):
    """ Tests for the OBSolve save and load methods."""

    def test_save_load_01(self):
        """ Solve a basic OBSolve problem. Save the results to file. Set the
            results in the OBSolve object to null. Load the results from file
            and check that they match the original values.
        """

        ob_solve_02 = ob_solve.OBSolve().from_json("maxwellbloch/tests/json/ob_solve_02.json")

        states_t = ob_solve_02.solve()

        states_t_loaded = ob_solve_02.solve(recalc=False)

        self.assertTrue((states_t == states_t_loaded).all())

def main():

    unittest.main(verbosity=3)

if __name__ == "__main__":

    status = main()
    sys.exit(status)
