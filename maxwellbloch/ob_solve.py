# -*- coding: utf-8 -*-

import os
import sys

import json
from time import strftime

import qutip as qu

from maxwellbloch import ob_atom


# Main

class OBSolve(object):
    """docstring for OBSolve"""

    def __init__(self, ob_atom={}, t_min=0.0, t_max=1.0, t_steps=100, 
                 method='mesolve', opts={}, savefile=None):
    
        self.build_ob_atom(ob_atom)

        self.build_tlist(t_min, t_max, t_steps)    

        self.method = method
        self.build_opts(opts)

        self.build_savefile(savefile)

    def __repr__(self):
        return ("OBSolve(ob_atom={0}, " +
                "t_min={1}, " +
                "t_max={2}, " +
                "t_steps={3}, " +
                "method={4}, " +
                "opts={5})").format(self.ob_atom,
                                    self.t_min, 
                                    self.t_max, 
                                    self.t_steps, 
                                    self.method, 
                                    self.opts)

    def build_ob_atom(self, ob_atom_dict):
    
        self.ob_atom = ob_atom.OBAtom(**ob_atom_dict)
        return self.ob_atom

    def build_tlist(self, t_min, t_max, t_steps):

        from numpy import linspace

        self.t_min=t_min
        self.t_max=t_max
        self.t_steps=t_steps

        self.tlist = linspace(t_min, t_max, t_steps+1)
        return self.tlist

    def build_opts(self, opts):

        self.opts = qu.Options()
        return self.opts

    def solve(self, rho0=None, e_ops=[], 
              recalc=True, show_pbar=False):

        if self.method == 'mesolve':
            self.ob_atom.mesolve(self.tlist, rho0=rho0, e_ops=e_ops, 
                                opts=qu.Options(), recalc=recalc, 
                                savefile=self.savefile, show_pbar=show_pbar)

        return self.ob_atom.result

    def states_t(self):

        return self.ob_atom.states_t()

    def build_savefile(self, savefile):

        self.savefile = savefile

    @classmethod
    def from_json_str(cls, json_str):
        json_dict = json.loads(json_str)
        return cls(**json_dict)

    @classmethod
    def from_json(cls, file_path):
        with open(file_path) as json_file:
            json_dict = json.load(json_file)

            # This needs to be here to get the savefile name from json
            if 'savefile' not in json_dict:
                s = os.path.splitext(file_path)[0]
                json_dict['savefile'] = s

                print(json_dict['savefile'])

        return cls(**json_dict)

def main():

    print(OBSolve())

if __name__ == '__main__':
    status = main()
    sys.exit(status)