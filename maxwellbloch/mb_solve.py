# -*- coding: utf-8 -*-

import os
import sys

import numpy as np

from maxwellbloch import ob_solve, t_funcs

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
        Delta_range = np.unique(np.sort(np.append(Delta_range, Delta2_range)))

        print(Delta_range)

        return Delta_range

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

    def solve(self, rho0=None, recalc=True, show_pbar=False):

        num_fields = len(self.ob_atom.fields)

        Omegas_zt = np.zeros((num_fields, len(self.zlist), len(self.tlist)), 
                             dtype=complex)        

        result_zt = [None]*len(self.zlist) 

        savefile_exists = os.path.isfile(str(self.savefile) + '.qu')

        if (recalc or not savefile_exists):
      
            Omega_0_agrs = {}
            for i, f in enumerate(self.ob_atom.fields):
                Omegas_zt[i][0] = f.rabi_freq_t_func(tlist, rabi_freq_t_args)
                Omega_0_args.update(f.rabi_freq_t_args)  

            len_z = len(self.zlist)-1 # only for printing progress

    ### First step: Euler

            print('j: 0/{:d}, z = {:.3f}'.format(len_z, self.zlist[0])) # Progress 

            h_0 = self.zlist[1] - self.zlist[0]

            # Number density
            N = self.num_density_z_func(self.zlist[1], self.num_density_z_args) 

            pass

        return Omegas_zt, result_zt                        

def main():

    # print(MBSolve())

    a = MBSolve()
    print(a)
    a.solve()

if __name__ == '__main__':
    status = main()
    sys.exit(status)