# -*- coding: utf-8 -*-

import os
import sys

import numpy as np

import qutip as qu # TODO: remve if not needed in average_states_Delta

from maxwellbloch import ob_solve, t_funcs

from copy import deepcopy



class MBSolve(ob_solve.OBSolve):

    def __init__(self, ob_atom={}, t_min=0.0, t_max=1.0, t_steps=100, 
                 method='mesolve', opts={}, savefile=None, z_min=0.0, 
                 z_max=1.0, z_steps=5, z_steps_inner=2, 
                 num_density_z_func=None, num_density_z_args={}, 
                 interaction_strengths=[], velocity_classes={}):

        super().__init__(ob_atom, t_min, t_max, t_steps, 
                 method, opts, savefile)

        self.build_zlist(z_min, z_max, z_steps, z_steps_inner)

        self.build_number_density(num_density_z_func, num_density_z_args,
                                  interaction_strengths)

        self.build_velocity_classes(velocity_classes)

    def __repr__(self):
        """ TODO: needs interaction strengths """

        return ("MBSolve(ob_atom={0}, " +
                "t_min={1}, " +
                "t_max={2}, " +
                "t_steps={3}, " +
                "method={4}, " +
                "opts={5}, " +
                "savefile={6}, " +
                "z_min={7}, " +
                "z_max={8}, " +
                "z_steps={9}, " +
                "z_steps_inner={10}, " +
                "num_density_z_func={11}, " +
                "num_density_z_args={12}, " +
                "velocity_classes={13})").format(self.ob_atom,
                                    self.t_min, 
                                    self.t_max, 
                                    self.t_steps, 
                                    self.method, 
                                    self.opts,
                                    self.savefile,
                                    self.z_min, 
                                    self.z_max, 
                                    self.z_steps,
                                    self.z_steps_inner, 
                                    self.num_density_z_func, 
                                    self.num_density_z_args, 
                                    self.velocity_classes)

        ## TODO: move opts and savefile to end

    def build_zlist(self, z_min, z_max, z_steps, z_steps_inner):

        self.z_min = z_min
        self.z_max = z_max
        self.z_steps = z_steps 
        self.z_steps_inner = z_steps_inner

        # TODO: does this even work? Did it ever? Where are inner steps?
        z_inner_stepsize = (z_max - z_min)/(z_steps*z_steps_inner + 1)
        zlist = np.linspace(z_min + z_inner_stepsize, z_max, z_steps + 1)
        zlist = np.insert(zlist, 0, z_min) # One more for first step    

        self.zlist = zlist

        # print(zlist)
        return zlist

    def build_velocity_classes(self, velocity_classes={}):

        self.velocity_classes = {}

        # What if only some elements are defined?
        if velocity_classes:
            self.velocity_classes = velocity_classes
        
        else:
            self.velocity_classes['thermal_delta_min'] = 0.0
            self.velocity_classes['thermal_delta_max'] = 0.0
            self.velocity_classes['thermal_delta_steps'] = 0
            self.velocity_classes['thermal_delta_inner_min'] = 0.0
            self.velocity_classes['thermal_delta_inner_max'] = 0.0
            self.velocity_classes['thermal_delta_inner_steps'] = 0
            self.velocity_classes['thermal_width'] = 1.0

        Delta_min = 2*np.pi*self.velocity_classes['thermal_delta_min']
        Delta_max = 2*np.pi*self.velocity_classes['thermal_delta_max']
        Delta_steps = self.velocity_classes['thermal_delta_steps']
        Delta_range = np.linspace(Delta_min, Delta_max, Delta_steps+1)

        # Narrow set of Delta classes around 0
        Delta2_min = 2*np.pi*self.velocity_classes['thermal_delta_inner_min']
        Delta2_max =  2*np.pi*self.velocity_classes['thermal_delta_inner_max']
        Delta2_steps = self.velocity_classes['thermal_delta_inner_steps']
        Delta2_range = np.linspace(Delta2_min, Delta2_max, Delta2_steps+1)

        # Merge the two ranges
        self.thermal_delta_list = np.unique(np.sort(np.append(Delta_range, 
                                                              Delta2_range)))

        self.thermal_weights = maxwell_boltzmann(self.thermal_delta_list, 
                                        self.velocity_classes['thermal_width'])

        return self.thermal_delta_list, self.thermal_weights

    def build_number_density(self, num_density_z_func, num_density_z_args,
                         interaction_strengths):

        self.interaction_strengths = interaction_strengths
        self.g = np.zeros(len(interaction_strengths))
        for i, g in enumerate(interaction_strengths):
            self.g[i] = 2*np.pi*g

        if num_density_z_func:
            self.num_density_z_func = getattr(t_funcs, num_density_z_func)
        else:
            self.num_density_z_func = t_funcs.square_1

        if num_density_z_args:
            self.num_density_z_args = num_density_z_args
        else:
            self.num_density_z_args = {'on_1': 0.0, 'off_1': 1.0, 'ampl_1':1.0}

    def init_Omegas_zt(self):

        self.Omegas_zt = np.zeros((len(self.ob_atom.fields), len(self.zlist), 
                                  len(self.tlist)), dtype=complex) 

        return self.Omegas_zt

    def init_result_zt(self):

        self.result_zt = [None]*len(self.zlist)     

        return self.result_zt

    def mbsolve(self, rho0=None, recalc=True, show_pbar=False):

        # num_fields = len(self.ob_atom.fields)

        self.init_Omegas_zt()      

        self.init_result_zt()

        result_zt = [None]*len(self.zlist) 

        savefile_exists = os.path.isfile(str(self.savefile) + '.qu')

        if (recalc or not savefile_exists):

            # Set Omega_0 to the initial tfunc of the field
      
            Omega_0_args = {}
            for f_i, f in enumerate(self.ob_atom.fields):
                self.Omegas_zt[f_i][0] = f.rabi_freq_t_func(self.tlist, 
                                                          f.rabi_freq_t_args)
                Omega_0_args.update(f.rabi_freq_t_args)  

            len_z = len(self.zlist)-1 # only for printing progress

    ### First step: Euler

            print('j: 0/{:d}, z = {:.3f}'.format(len_z, self.zlist[0])) # Progress 

            # H_Omega should be as it was defined
            # self.ob_atom.build_H_Omega()

            # TODO: If I did this with a np array instead of an array of 
            # qu results, it might be easier
            result_Delta = self.solve_over_thermal_detunings()

            # Need any one of the result objects to fill with averaged states
            self.result_zt[0] = deepcopy(result_Delta[0])

            for k, t in enumerate(self.tlist):
                self.result_zt[0].states[k] = self.average_states_over_thermal_detunings(result_Delta, k)

            self.z_step_fields_euler(j=1)

        return self.Omegas_zt, self.result_zt                        

    def z_step_fields_euler(self, j):
        """ Should this take and return? """

        h_0 = self.zlist[j] - self.zlist[j-1] # J SHOULD BE 1 for first step

        # Number density
        N = self.num_density_z_func(self.zlist[j], self.num_density_z_args) 

        # Array, one derivative for each field
        # dOmegas_dz =  np.zeros(len(self.ob_atom.fields))

        for k, t in enumerate(self.tlist):

            rho = self.result_zt[j-1].states[k]

            for f_i, f in enumerate(self.ob_atom.fields):

                # print(f.coupled_levels)

                # a, b = f.coupled_levels
                sum_rho_ij = 0
                for cl in f.coupled_levels:
                    sum_rho_ij += rho[cl[0], cl[1]]

                print('sum_rho_ij: ', sum_rho_ij)

                dOmega_f_dz = 1j*N*self.g[f_i]*sum_rho_ij
                self.Omegas_zt[f_i, j, k] = self.Omegas_zt[f_i, j-1, k] + h_0*dOmega_f_dz

    def z_step_fields_adams_bashforth(self):

        pass

    def solve_over_thermal_detunings(self):

        len_Delta = len(self.thermal_delta_list)

        result_Delta = [None]*len_Delta

        for Delta_i, Delta in enumerate(self.thermal_delta_list):
            """ Which field? All of them!"""

            # Progress Indicator
            print('    i: {:d}/{:d}, Delta = {:.3f}'
                  .format(Delta_i, len_Delta-1, Delta))  

            # Shift each detuning by Delta
            self.ob_atom.shift_H_Delta([Delta]*len(self.ob_atom.fields))

            result_Delta[Delta_i] = super().solve()

        return result_Delta

    def average_states_over_thermal_detunings(self, result_Delta, t_step):

        len_Delta = len(self.thermal_delta_list)

        states_Delta = np.zeros((len_Delta, self.ob_atom.num_states, 
                                self.ob_atom.num_states), 
                                dtype=np.complex)

        for Delta_i, Delta in enumerate(self.thermal_delta_list):
            states_Delta[Delta_i] = result_Delta[Delta_i].states[t_step].full()

        # Average of states TODO: Do we need this to be a qu object?
        return qu.Qobj(np.average(states_Delta, axis=0, 
                        weights=self.thermal_weights))

### Helper Functions

# def Omega_intp(t, args):

#     t_range = args['t_range']
#     Omega_range = args['Omega_range']

#     Omega_intp = interp1d(t_range, Omega_range,
#                           bounds_error=False, fill_value=0.0)

#     return Omega_intp(t)

def maxwell_boltzmann(v, fwhm):
    """ Maxwell Boltzmann probability distribution function. """

    # TODO: Allow offset, v_0.

    return 1./(fwhm*np.sqrt(np.pi))*np.exp(-(v/fwhm)**2)

def main():

    # print(MBSolve())

    a = MBSolve()
    print(a)
    # a.mbsolve()

if __name__ == '__main__':
    status = main()
    sys.exit(status)