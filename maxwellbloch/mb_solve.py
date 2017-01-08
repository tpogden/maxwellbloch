# -*- coding: utf-8 -*-

import os
import sys

import numpy as np

import qutip as qu # TODO: remve if not needed in average_states_Delta

from maxwellbloch import ob_solve, t_funcs

from copy import deepcopy

from tqdm import tqdm # progress bar

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

        self.init_Omegas_zt()
        self.init_states_zt()

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
        # z_inner_stepsize = (z_max - z_min)/(z_steps*z_steps_inner + 1)
        # zlist = np.linspace(z_min + z_inner_stepsize, z_max, z_steps + 1)
        # zlist = np.insert(zlist, 0, z_min) # One more for first step    

        z_steps_total = z_steps #*z_steps_inner

        zlist = np.linspace(z_min, z_max, z_steps_total+1)

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
                                  len(self.tlist)), dtype=np.complex) 

        ### Set the initial Omegas to the field time func values
        for f_i, f in enumerate(self.ob_atom.fields):
            self.Omegas_zt[f_i][0] = 2.0*np.pi*f.rabi_freq* \
                f.rabi_freq_t_func(self.tlist, f.rabi_freq_t_args)

        return self.Omegas_zt

    def init_states_zt(self):

        # TODO: Change states_zt to state_zt. Omegas refers to the fact that 
        # there are multiple fields.

        self.states_zt = np.zeros((len(self.zlist), len(self.tlist), 
                            self.ob_atom.num_states, self.ob_atom.num_states), 
                            dtype=np.complex) 

        return self.states_zt

    def mbsolve(self, rho0=None, recalc=True, show_pbar=False):

        self.init_Omegas_zt()
        self.init_states_zt()

        savefile_exists = os.path.isfile(str(self.savefile) + '.qu')

        if (recalc or not savefile_exists):

            len_z = len(self.zlist)-1 # only for printing progress

            rabi_freq_ones = np.ones(len(self.ob_atom.fields))

            ### Set initial states at z=0
            self.states_zt[0, :] = \
                self.solve_and_average_over_thermal_detunings()

        ### All Steps:

            h = (self.zlist[2] - self.zlist[1])/self.z_steps_inner # Stepsize

            # print('z_step = {0}'.format(h))

            for j, z in tqdm(enumerate(self.zlist[:-1]), 
                             total=len(self.zlist)-1):

                # print('j: {:d}/{:d}, z = {:.3f}'.format(j, len_z, z))

                # Set initial fields and state
                Omegas_z_this = self.Omegas_zt[:, j, :]
                states_z_this = self.states_zt[j, :]

                z_this = z

            ### Inner z loop

                for jj in np.arange(self.z_steps_inner):

                    # print('  jj: {:d}/{:d}, zz = {:.3f}'
                          # .format(jj, self.z_steps_inner, z_this))

                    z_next = z_this + h

                    self.solve_and_average_over_thermal_detunings()

                    Omegas_z_next = self.z_step_fields_euler(z_this, z_next, 
                        Omegas_z_this)

                    Omegas_z_next_args = \
                        self.get_Omegas_intp_t_args(Omegas_z_next)
                    self.ob_atom.set_H_Omega(rabi_freq_ones,
                        self.get_Omegas_intp_t_funcs(), Omegas_z_next_args)

                    # Set up for next inner step
                    Omegas_z_this = Omegas_z_next
                    z_this = z_next

                self.states_zt[j+1, :] = self.states_t()
                self.Omegas_zt[:, j+1, :] = Omegas_z_next

        return self.Omegas_zt, self.states_zt                        

    def z_step_fields_euler(self, z_this, z_next, Omegas_z_this):
        """ For the current state of the atom, given field Omegas_z_this.

        """

        h = z_next - z_this

        N = self.num_density_z_func(z_next, self.num_density_z_args)

        Omegas_z_next = np.zeros((len(self.ob_atom.fields), len(self.tlist)),
                               dtype=np.complex)


        for f_i, f in enumerate(self.ob_atom.fields):

            sum_coh = self.ob_atom.get_field_sum_coherence(f_i)

            dOmega_f_dz = 1.0j*N*self.g[f_i]*sum_coh

            Omegas_z_next[f_i, :] = Omegas_z_this[f_i, :] + h*dOmega_f_dz

        return Omegas_z_next

    def z_step_fields_adams_bashforth(self, j):

        pass

    def get_Omegas_intp_t_funcs(self):

        rabi_freq_t_funcs = []
        for f_i, f in enumerate(self.ob_atom.fields, start=1):
            rabi_freq_t_funcs.append('intp_{0}'.format(f_i))
        return rabi_freq_t_funcs

    def get_Omegas_intp_t_args(self, Omegas_z):
        """ Return the values of Omegas at a given point as a list of 
            args for interpolation
            TODO complete doc!  
    
            e.g. [{'tlist_1': [], 'ylist_1': []}, 
                  {'tlist_2': [], 'ylist_2': []}]

        """

        fields_args = [{}]*len(self.ob_atom.fields)

        for f_i, f in enumerate(Omegas_z):
            fields_args[f_i] = {'tlist_{0}'.format(f_i+1): self.tlist, 
                                'ylist_{0}'.format(f_i+1): Omegas_z[f_i]}

        return fields_args

    def solve_and_average_over_thermal_detunings(self):

        states_t_Delta = np.zeros((len(self.thermal_delta_list), 
                                  len(self.tlist), self.ob_atom.num_states, 
                                  self.ob_atom.num_states), dtype=np.complex)

        for Delta_i, Delta in enumerate(self.thermal_delta_list):

            # Progress Indicator
            # print('    i: {:d}/{:d}, Delta = {:.3f}'
                  # .format(Delta_i, len(self.thermal_delta_list)-1, 
                          # Delta/(2*np.pi)))    

            # Shift each detuning by Delta
            self.ob_atom.shift_H_Delta([Delta]*len(self.ob_atom.fields))            

            self.solve(opts=qu.Options(max_step=self.t_step()))

            states_t_Delta[Delta_i] = self.states_t()

        thermal_states_t = np.average(states_t_Delta, axis=0, 
                                      weights=self.thermal_weights)

        return thermal_states_t

    def fields_area(self):
        """ Get the integrated pulse area of each field.

        Returns:
            np.array [num_fields, num_z_steps]: Integrated area of each field 
            over time
        """

        return np.trapz(np.abs(self.Omegas_zt), self.tlist, axis=2)

### Helper Functions

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