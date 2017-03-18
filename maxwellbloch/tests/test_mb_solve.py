""" 
Unit tests for the module.

Thomas Ogden <t@ogden.eu>

"""

import sys

import unittest

import numpy as np

from maxwellbloch import mb_solve

json_str_01 = (
'{'
'  "ob_atom": {'
'    "decays": ['
'      { "channels": [[0,1]], '
'        "rate": 0.0'
'      }'
'    ],'
'    "energies": [],'
'    "fields": ['
'      {'
'        "coupled_levels": [[0, 1]],'
'        "detuning": 0.0,'
'        "detuning_positive": true,'
'        "label": "probe",'
'        "rabi_freq": 0.01,'
'        "rabi_freq_t_args": {' 
'           "ampl_1": 1.0,'
'           "centre_1": 0.0, ' 
'           "fwhm_1": 0.1' 
'        },'
'        "rabi_freq_t_func": "gaussian_1"'
'      }'
'    ],'
'    "num_states": 2'
'  },'
''
'  "t_min": 0.0,'
'  "t_max": 1.0,'
'  "t_steps": 100,'
''
'  "z_min": -0.2,'
'  "z_max": 1.2,'
'  "z_steps": 4,'
'  "z_steps_inner": 1,'
''
'  "num_density_z_func": "square_1",'
'  "num_density_z_args": { '
'    "on_1":0.0, '
'    "off_1":1.0,'
'    "ampl_1": 1.0e0'
'  },'
'  "interaction_strengths": [1.0],'
''
'  "velocity_classes": {'
'    "thermal_delta_min": -1.0,'
'    "thermal_delta_max":  1.0,'
'    "thermal_delta_steps": 2,'
'    "thermal_delta_inner_min": 0.0,'
'    "thermal_delta_inner_max": 0.0,'
'    "thermal_delta_inner_steps": 0,'
'    "thermal_width": 1.0'
'  },'
''
'  "method": "mesolve",'
'  "opts": {},'
''
'  "savefile": "json_str_01"'
''
'}')

json_str_lamda = (
'{'
'  "ob_atom": {'
'    "decays": ['
'      {'
'        "channels": ['
'          ['
'            0,'
'            1'
'          ],'
'          ['
'            1,'
'            2'
'          ]'
'        ],'
'        "rate": 1.0'
'      }'
'    ],'
'    "energies": [],'
'    "fields": ['
'      {'
'        "coupled_levels": ['
'          ['
'            0,'
'            1'
'          ]'
'        ],'
'        "detuning": 0.0,'
'        "detuning_positive": true,'
'        "label": "probe",'
'        "rabi_freq": 0.1,'
'        "rabi_freq_t_args": {'
'          "ampl_1": 1.0,'
'          "centre_1": 0.0,'
'          "fwhm_1": 0.1'
'        },'
'        "rabi_freq_t_func": "gaussian_1"'
'      },'
'      {'
'        "coupled_levels": ['
'          ['
'            1,'
'            2'
'          ]'
'        ],'
'        "detuning": 0.0,'
'        "detuning_positive": false,'
'        "label": "coupling",'
'        "rabi_freq": 10.0,'
'        "rabi_freq_t_args": {'
'          "ampl_2": 1.0,'
'          "on_2": 0.0,'
'          "off_2": 1.0'
'        },'
'        "rabi_freq_t_func": "square_2"'
'      }'
'    ],'
'    "num_states": 3'
'  },'
'  "t_min": 0.0,'
'  "t_max": 1.0,'
'  "t_steps": 100,'
'  "z_min": -0.2,'
'  "z_max": 1.2,'
'  "z_steps": 20,'
'  "z_steps_inner": 2,'
'  "num_density_z_func": "square_1",'
'  "num_density_z_args": {'
'    "on_1": 0.0,'
'    "off_1": 1.0,'
'    "ampl_1": "1.0e3"'
'  },'
'  "interaction_strengths": ['
'    1.0'
'  ],'
'  "velocity_classes": {'
'    "thermal_delta_min": -5.0,'
'    "thermal_delta_max": 5.0,'
'    "thermal_delta_steps": 4,'
'    "thermal_delta_inner_min": 0.0,'
'    "thermal_delta_inner_max": 0.0,'
'    "thermal_delta_inner_steps": 0,'
'    "thermal_width": 1.0'
'  },'
'  "method": "mesolve",'
'  "opts": {},'
'  "savefile": null'
'}')

json_no_atoms = (
'{'
''
'  "ob_atom": {'
'    "decays": ['
'      { "channels": [[0,1]], '
'        "rate": 0.0'
'      }'
'    ],'
'    "energies": [],'
'    "fields": ['
'      {'
'        "coupled_levels": [[0, 1]],'
'        "detuning": 0.0,'
'        "detuning_positive": true,'
'        "label": "probe",'
'        "rabi_freq": 5.0,'
'        "rabi_freq_t_args": { '
'          "ampl_1": 1.0,'
'          "centre_1": 0.0,  '
'          "fwhm_1": 0.1 '
'        },'
'        "rabi_freq_t_func": "gaussian_1"'
'      }'
'    ],'
'    "num_states": 2'
'  },'
''
'  "t_min": -0.5,'
'  "t_max": 1.0,'
'  "t_steps": 100,'
''
'  "z_min": -0.2,'
'  "z_max": 1.2,'
'  "z_steps": 4,'
'  "z_steps_inner": 1,'
''
'  "num_density_z_func": "square_1",'
'  "num_density_z_args": { '
'    "on_1": 0.0, '
'    "off_1":1.0,'
'    "ampl_1": 0.0'
'  },'
'  "interaction_strengths": [1.0],'
''
'  "velocity_classes": {'
'    "thermal_delta_min": 0.0,'
'    "thermal_delta_max": 0.0,'
'    "thermal_delta_steps": 0,'
'    "thermal_delta_inner_min": 0.0,'
'    "thermal_delta_inner_max": 0.0,'
'    "thermal_delta_inner_steps": 0,'
'    "thermal_width": 1.0'
'  },'
''
'  "method": "mesolve",'
'  "opts": {},'
''
'  "savefile": null'
'}'
)

class TestInit(unittest.TestCase):

    def test_init_default(self):
        """  Test Default Initialise """

        mb_solve_00 = mb_solve.MBSolve()

        self.assertEqual(mb_solve_00.ob_atom.num_states, 1)

        # TODO: And the rest!

    def test_init_00(self):

        mb_solve_00 = mb_solve.MBSolve().from_json_str(json_str_01)
        # print(mb_solve_00)

@unittest.skip("TODO")
class TestSolveOverThermalDetunings(unittest.TestCase):

    def test_00(self):

        mb_solve_00 = mb_solve.MBSolve().from_json_str(json_str_01)

        result_Delta = mb_solve_00.solve_over_thermal_detunings()

        self.assertEqual(len(result_Delta),
                         len(mb_solve_00.thermal_delta_list))

class TestMBSolve(unittest.TestCase):

    # @unittest.skip("Testing others")
    def test_mb_solve(self):
        """ Basic test of mb_solve method. """

        mb_solve_00 = mb_solve.MBSolve().from_json_str(json_str_01)

        mb_solve_00.mbsolve()

    def test_no_atoms(self):
        """ Setting the number density ampl to 0.0, i.e. no atoms. The end
            pulse should be the same as the start. """

        json_no_atoms = (
            '{'
            ''
            '  "ob_atom": {'
            '    "decays": ['
            '      { "channels": [[0,1]], '
            '        "rate": 0.0'
            '      }'
            '    ],'
            '    "energies": [],'
            '    "fields": ['
            '      {'
            '        "coupled_levels": [[0, 1]],'
            '        "detuning": 0.0,'
            '        "detuning_positive": true,'
            '        "label": "probe",'
            '        "rabi_freq": 5.0,'
            '        "rabi_freq_t_args": { '
            '          "ampl_1": 1.0,'
            '          "centre_1": 0.0,  '
            '          "fwhm_1": 0.1 '
            '        },'
            '        "rabi_freq_t_func": "gaussian_1"'
            '      }'
            '    ],'
            '    "num_states": 2'
            '  },'
            ''
            '  "t_min": -0.5,'
            '  "t_max": 1.0,'
            '  "t_steps": 100,'
            ''
            '  "z_min": -0.2,'
            '  "z_max": 1.2,'
            '  "z_steps": 4,'
            '  "z_steps_inner": 1,'
            ''
            '  "num_density_z_func": "square_1",'
            '  "num_density_z_args": { '
            '    "on_1": 0.0, '
            '    "off_1":1.0,'
            '    "ampl_1": 0.0'
            '  },'
            '  "interaction_strengths": [1.0],'
            ''
            '  "velocity_classes": {'
            '    "thermal_delta_min": 0.0,'
            '    "thermal_delta_max": 0.0,'
            '    "thermal_delta_steps": 0,'
            '    "thermal_delta_inner_min": 0.0,'
            '    "thermal_delta_inner_max": 0.0,'
            '    "thermal_delta_inner_steps": 0,'
            '    "thermal_width": 1.0'
            '  },'
            ''
            '  "method": "mesolve",'
            '  "opts": {},'
            ''
            '  "savefile": "json_no_atoms"'
            '}'
        )

        mbs = \
            mb_solve.MBSolve().from_json_str(json_no_atoms)

        mbs.mbsolve(step='euler')

        self.assertEqual(mbs.Omegas_zt.shape, (1, 5, 101))

        # Check that the field at the end of the medium matches the field
        # at the start of the medium.
        self.assertTrue(np.allclose(mbs.Omegas_zt[:, 0, :],
                                    mbs.Omegas_zt[:, -1, :], rtol=1.0e-6))

        # self.assertEqual(mbs.Omegas_zt[:,0,:], mbs.Omegas_zt[:,-1,:])

class TestSaveLoad(unittest.TestCase):
    """ Tests for the MBSolve save and load methods. """

    def test_save_load_01(self):
        """ Solve a basic MBSolve problem. Save the results to file. Set the 
            results in the MBSolve object to null. Load the results from
            file and check that they equal the original values. 
        """

        mb_solve_01 = mb_solve.MBSolve().from_json_str(json_str_01)

        Omegas_zt, states_zt = mb_solve_01.mbsolve()

        mb_solve_01.save_results()

        mb_solve_01.Omegas_zt = None
        mb_solve_01.states_zt = None

        mb_solve_01.load_results()

        Omegas_zt_loaded = mb_solve_01.Omegas_zt
        states_zt_loaded = mb_solve_01.states_zt

        self.assertTrue((Omegas_zt == Omegas_zt_loaded).all())
        self.assertTrue((states_zt == states_zt_loaded).all())

    def test_save_load_no_recalc(self):
        
        mb_solve_01 = mb_solve.MBSolve().from_json_str(json_str_01)

        Omegas_zt, states_zt = mb_solve_01.mbsolve()

        mb_solve_01.save_results()

        mb_solve_01.Omegas_zt = None
        mb_solve_01.states_zt = None

        Omegas_zt, states_zt = mb_solve_01.mbsolve(recalc=False)

        Omegas_zt_loaded = mb_solve_01.Omegas_zt
        states_zt_loaded = mb_solve_01.states_zt

        self.assertTrue((Omegas_zt == Omegas_zt_loaded).all())
        self.assertTrue((states_zt == states_zt_loaded).all())

    def test_no_atoms_ab(self):

        mbs = \
            mb_solve.MBSolve().from_json_str(json_no_atoms)

        mbs.mbsolve(step='ab') 

class TestBuildZlist(unittest.TestCase):

    def test_00(self):

        mb_solve_00 = mb_solve.MBSolve()        

        zlist = np.array([0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.])

        self.assertTrue(np.allclose(mb_solve_00.zlist, zlist, rtol=1.0e-6))

class TestInsertFirstInnerZStep(unittest.TestCase):

    def test_00(self):

        mb_solve_00 = mb_solve.MBSolve()

        zlist = np.array([0., .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.])
        
        self.assertFalse(len(mb_solve_00.zlist) == len(zlist))

        mb_solve_00.insert_first_inner_z_step()

        self.assertTrue(np.allclose(mb_solve_00.zlist, zlist, rtol=1.0e-6))

class TestGetOmegasIntpTFuncs(unittest.TestCase):

    def test_one_fields(self):

        mb_solve_00 = mb_solve.MBSolve().from_json_str(json_str_01)

        self.assertEqual(mb_solve_00.get_Omegas_intp_t_funcs(),
                         ['intp_1'])

    def test_two_fields(self):

        mb_solve_lamda = mb_solve.MBSolve().from_json_str(json_str_lamda)

        self.assertEqual(mb_solve_lamda.get_Omegas_intp_t_funcs(),
                         ['intp_1', 'intp_2'])

def main():

    suite = unittest.TestSuite()
    # suite.addTest(TestBuildZlist("test_00"))
    # suite.addTest(TestInsertFirstInnerZStep("test_00"))
    suite.addTest(TestMBSolve("test_no_atoms"))
    # suite.addTest(TestMBSolve("test_no_atoms_ab"))
    runner = unittest.TextTestRunner()

    # runner.run(suite) # Run suite

    unittest.main(verbosity=3) # Run all

if __name__ == "__main__":
    status = main()
    sys.exit(status)

