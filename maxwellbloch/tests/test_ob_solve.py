# -*- coding: utf-8 -*-

"""
Unit tests for the OBSolve class.

Thomas Ogden <t@ogden.eu>

"""

import os
import unittest

import numpy as np

from maxwellbloch import ob_solve, t_funcs

# Absolute path of tests/json directory, so that tests can be called from
# different directories.
JSON_DIR = os.path.abspath(os.path.join(__file__, '../', 'json'))

JSON_STR_02 = (
    '{'
    '  "atom": {'
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

        two_pulse_t_func = lambda t, args: (t_funcs.gaussian(0)(t, args) +
            t_funcs.gaussian(1)(t, args))

        two_pulse_t_args = {"ampl_0": 1.0, "centre_0": 0.0, "fwhm_0": 0.1,
            "ampl_1": 2.0, "centre_1": 0.5, "fwhm_1": 0.1, }

        ob_solve_02.set_field_rabi_freq_t_func(0, two_pulse_t_func)
        ob_solve_02.set_field_rabi_freq_t_args(0, two_pulse_t_args)

        field_0 = ob_solve_02.atom.fields[0]

        self.assertAlmostEqual(field_0.rabi_freq_t_func(0.0,
            field_0.rabi_freq_t_args), 1.0)
        self.assertAlmostEqual(field_0.rabi_freq_t_func(0.5,
            field_0.rabi_freq_t_args), 2.0)
        self.assertAlmostEqual(field_0.rabi_freq_t_func(1.0,
            field_0.rabi_freq_t_args), 0.0)

class TestSolve(unittest.TestCase):

    def test_two_level_rabi_oscillations(self):
        """ Solve the optical Bloch equations for the two-level atom. 
        
            Notes:
                See https://en.wikipedia.org/wiki/Rabi_cycle
        """

        RABI_FREQ = 5.0
        atom_dict = {"fields": [{"coupled_levels": [[0, 1]], 
            "rabi_freq": RABI_FREQ}], "num_states": 2}
        obs = ob_solve.OBSolve(atom=atom_dict, t_min=0.0, t_max=1.0, 
            t_steps=100)
        obs.solve()

        # Get the populations
        pop_0 = np.absolute(obs.states_t()[:, 0, 0])
        pop_1 = np.absolute(obs.states_t()[:, 1, 1])

        # The solution is known, we should have Rabi cycling at the frequency.
        known_0 = np.cos(2.0*np.pi*RABI_FREQ*obs.tlist/2.0)**2
        known_1 = np.sin(2.0*np.pi*RABI_FREQ*obs.tlist/2.0)**2

        self.assertTrue(np.allclose(pop_0, known_0, rtol=1.e-5, atol=1.e-5))
        self.assertTrue(np.allclose(pop_1, known_1, rtol=1.e-5, atol=1.e-5))

        # If you want to take a look
        # import matplotlib.pyplot as plt
        # plt.plot(obs.tlist, pop_0)
        # plt.plot(obs.tlist, known_0, ls='dashed')
        # plt.plot(obs.tlist, pop_1)
        # plt.plot(obs.tlist, known_1, ls='dashed')
        # plt.show()
    
    def test_two_level_with_opts(self):
        """ Same as test_two_level_rabi_oscillations() but with opts set such
            that the tolerances are lower. The results will be less 
            accurate.
        """

        RABI_FREQ = 5.0
        atom_dict = {"fields": [{"coupled_levels": [[0, 1]],
                                 "rabi_freq": RABI_FREQ}], "num_states": 2,
                                 "initial_state": [1., 0.]}
        obs = ob_solve.OBSolve(atom=atom_dict, t_min=0.0, t_max=1.0,
                               t_steps=100, opts={'atol': 1e-6, 'rtol': 1e-4})
        obs.solve()

        # Get the populations
        pop_0 = np.absolute(obs.states_t()[:, 0, 0])
        pop_1 = np.absolute(obs.states_t()[:, 1, 1])

        # The solution is known, we should have Rabi cycling at the frequency.
        known_0 = np.cos(2.0 * np.pi * RABI_FREQ * obs.tlist / 2.0)**2
        known_1 = np.sin(2.0 * np.pi * RABI_FREQ * obs.tlist / 2.0)**2

        # Compared with test_two_level_rabi_oscillations() we can only assert 
        # a lower tolerance to the known solution.
        self.assertTrue(np.allclose(pop_0, known_0, rtol=1.e-3, atol=1.e-3))
        self.assertTrue(np.allclose(pop_1, known_1, rtol=1.e-3, atol=1.e-3))

        # If you want to take a look
        # import matplotlib.pyplot as plt
        # plt.plot(obs.tlist, pop_0)
        # plt.plot(obs.tlist, known_0, ls='dashed')
        # plt.plot(obs.tlist, pop_1)
        # plt.plot(obs.tlist, known_1, ls='dashed')
        # plt.show()

    def test_two_level_with_inital_state(self):
        """ Same as test_two_level_rabi_oscillations() but with the initial
            state set so that the population starts in the upper level.
        """

        RABI_FREQ = 5.0
        atom_dict = {"fields": [{"coupled_levels": [[0, 1]], 
            "rabi_freq": RABI_FREQ}], "num_states": 2, 
            "initial_state": [0., 1.]}
        obs = ob_solve.OBSolve(atom=atom_dict, t_min=0.0, t_max=1.0, 
            t_steps=100)
        obs.solve()

        # Get the populations
        pop_0 = np.absolute(obs.states_t()[:, 0, 0])
        pop_1 = np.absolute(obs.states_t()[:, 1, 1])

        # The solution is as test_two_level_rabi_oscillations() but swapped
        known_0 = np.sin(2.0*np.pi*RABI_FREQ*obs.tlist/2.0)**2
        known_1 = np.cos(2.0*np.pi*RABI_FREQ*obs.tlist/2.0)**2

        self.assertTrue(np.allclose(pop_0, known_0, rtol=1.e-5, atol=1.e-5))
        self.assertTrue(np.allclose(pop_1, known_1, rtol=1.e-5, atol=1.e-5))

class TestJSON(unittest.TestCase):

    def test_to_from_json_str_00(self):

        ob_solve_00 = ob_solve.OBSolve()
        ob_solve_01 = ob_solve.OBSolve.from_json_str(ob_solve_00.to_json_str())

        self.assertEqual(ob_solve_00.to_json_str(), ob_solve_01.to_json_str())

    def test_from_json_str(self):

        ob_solve_02 = ob_solve.OBSolve().from_json_str(JSON_STR_02)

        self.assertEqual(ob_solve_02.t_min, 0.0)
        self.assertEqual(ob_solve_02.t_max, 1.0)
        self.assertEqual(ob_solve_02.t_steps, 100)
        self.assertEqual(ob_solve_02.method, "mesolve")

    def test_to_from_json_str_03(self):

        json_path = os.path.join(JSON_DIR, "ob_solve_03.json")

        obs = ob_solve.OBSolve().from_json(json_path)
        obs_test = ob_solve.OBSolve.from_json_str(obs.to_json_str())

        self.assertEqual(obs.to_json_str(), obs_test.to_json_str())

    def test_to_from_json(self):

        import os

        filepath = "test_ob_solve_02.json"

        ob_solve_02 = ob_solve.OBSolve().from_json_str(JSON_STR_02)

        ob_solve_02.to_json(filepath)

        ob_solve_03 = ob_solve.OBSolve().from_json(filepath)
        os.remove(filepath)

        self.assertEqual(ob_solve_02.to_json_str(),
                         ob_solve_03.to_json_str())

class TestSaveLoad(unittest.TestCase):
    """ Tests for the OBSolve save and load methods."""

    def test_save_load_01(self):
        """ Solve a basic OBSolve problem. Save the results to file. Set the
            results in the OBSolve object to null. Load the results from file
            and check that they match the original values.
        """

        json_path = os.path.join(JSON_DIR, "ob_solve_02.json")
        ob_solve_02 = ob_solve.OBSolve().from_json(json_path)

        states_t = ob_solve_02.solve()

        states_t_loaded = ob_solve_02.solve(recalc=False)

        self.assertTrue((states_t == states_t_loaded).all())
