# -*- coding: utf-8 -*-

import sys

from numpy import linspace, insert

from maxwellbloch import ob_solve

class MBSolve(ob_solve.OBSolve):

    def __init__(self, ob_atom={}, t_min=0.0, t_max=1.0, t_steps=100, 
                 method='mesolve', opts={}, savefile=None, z_min=0.0, 
                 z_max=1.0, z_steps=5, z_steps_inner=2, 
                 num_density_z_func=None, num_density_z_args={}, 
                 velocity_classes={}):

        super().__init__(ob_atom, t_min, t_max, t_steps, 
                 method, opts, savefile)

        self.build_zlist(z_min, z_max, z_steps, z_steps_inner)

        self.num_density_z_func = num_density_z_func

        self.num_density_z_args = num_density_z_args

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
                "velocity_classes={12})").format(self.ob_atom,
                                    self.t_min, 
                                    self.t_max, 
                                    self.t_steps, 
                                    self.method, 
                                    self.opts,
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
        zlist = linspace(z_min + z_inner_stepsize, z_max, z_steps + 1)
        zlist = insert(zlist, 0, z_min) # One more for first step    

        print(zlist)
        return zlist

    def build_velocity_classes(self, velocity_classes):        

        self.velocity_classes = velocity_classes

        # TODO: Build it here



def main():

    print(MBSolve())

if __name__ == '__main__':
    status = main()
    sys.exit(status)