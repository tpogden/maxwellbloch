# -*- coding: utf-8 -*-

import os
import sys
import json
import logging

import qutip as qu

from maxwellbloch import ob_atom


class OBSolve(object):
    """docstring for OBSolve"""

    def __init__(self, atom={}, t_min=0.0, t_max=1.0, t_steps=100,
                 method='mesolve', opts={}, savefile=None):

        self.build_atom(atom)

        self.build_tlist(t_min, t_max, t_steps)

        self.method = method
        self.build_opts(opts)

        self.build_savefile(savefile)

    def __repr__(self):
        return (
            "OBSolve(atom={0}, " +
            "t_min={1}, " +
            "t_max={2}, " +
            "t_steps={3}, " +
            "method={4}, " +
            "opts={5})").format(
                self.atom,
                self.t_min,
                self.t_max,
                self.t_steps,
                self.method,
                self.opts)

    def build_atom(self, atom_dict):

        self.atom = ob_atom.OBAtom(**atom_dict)
        return self.atom

    def build_tlist(self, t_min, t_max, t_steps):

        from numpy import linspace

        self.t_min = t_min
        self.t_max = t_max
        self.t_steps = t_steps

        self.tlist = linspace(t_min, t_max, t_steps + 1)
        return self.tlist

    def build_opts(self, opts={}):
        """ Build the options dict to be passed into the QuTiP solver.

            Any option available to the QuTiP solver is available here, we 
            provide defaults for solving the optical Bloch equations. See [0] 
            for details of all the available options.

            Notes:
                - For a stiff problem, it may help to set 'method' to 
                    'bdf' instead of 'adams'.
                - If the solver times out, try more 'nsteps', though
                    this will take longer. 
                - To speed up the solver, reduce atol and rtol, though this 
                    will reduce accuracy.

            Warning:
                There is no validation checking here. If you pass in an option
                which is not known to QuTiP it will throw an exception.
            
            [0]: http://qutip.org/docs/4.2/guide/dynamics/dynamics-options.html
        """
        self.opts = {
                'atol': 1e-8,
                'rtol': 1e-6, 
                'method': 'adams',
                'order': 12,
                'nsteps': 1000,
                'first_step': 0, 
                'max_step': self.t_step(),
                'min_step': 0,
                # Options below here I do not think will be useful, they are
                # mostly for the MC solver. But they may be set.
                # 'average_expect': True, # What is this
                # 'average_states': False, 
                # 'tidy': True,
                # 'rhs_reuse': False,
                # 'rhs_filename': None,
                # 'rhs_with_state': False,
                # 'store_final_state': False,
                # 'store_states': False,
                # 'steady_state_average': False, 
                # 'normalize_output':True, 
                # 'use_openmp': None
                # 'openmp_threads': None
        }
        if opts:
            self.opts.update(opts) # Update any specified in the parameter
        return self.opts

    def set_field_rabi_freq_t_func(self, field_idx, t_func):
        """ Set the Rabi frequency time function of a field to a new time
            function. This is useful when you want to set a custom function,
            not one available in t_funcs.py

        Args:
            field_idx: The field for which to set the Rabi frequency t_func
            t_func: Rabi frequency as a function of time, f(t, args)
        """

        self.atom.fields[field_idx].rabi_freq_t_func = t_func
        self.atom.build_operators() # Rebuild so H_Omega is updated

    def set_field_rabi_freq_t_args(self, field_idx, t_args):
        """ Set the Rabi frequency time function arguments. To be used with
            set_field_rabi_freq_t_func

        Args:
            field_idx: The field for which to set the Rabi frequency t_args
            t_args: A dict representing the args to go with the t_func.
        """

        self.atom.fields[field_idx].rabi_freq_t_args = t_args

    # TODO: Rename to obsolve for clarity when calling from derived class
    def solve(self, e_ops=[], opts=None, recalc=True,
              show_pbar=False, save=True):

        # When we're calling from MBSolve, we don't want to save each step.
        # So we pass in save=False.
        if save:
            savefile = self.savefile
        else:
            savefile = None

        # Choosing to overwrite opts at the solve stage.
        if opts:
            self.build_opts(opts)
        options = qu.Options(**self.opts)

        if self.method == 'mesolve':
            self.atom.mesolve(self.tlist, e_ops=e_ops, options=options, 
                recalc=recalc, savefile=savefile, show_pbar=show_pbar)

        return self.atom.states_t()  # self.atom.result

    def states_t(self):

        return self.atom.states_t()

    def t_step(self):

        return (self.t_max - self.t_min) / self.t_steps

    def build_savefile(self, savefile):

        self.savefile = savefile

    def savefile_exists(self):
        """ Returns true if savefile (with appended extension .qu) exists. """

        return os.path.isfile(str(self.savefile) + '.qu')

### JSON Serialise and Deserialise

    def get_json_dict(self):

        json_dict = {"atom": self.atom.get_json_dict(),
                     "t_min": self.t_min,
                     "t_max": self.t_max,
                     "t_steps": self.t_steps,
                     "method": self.method,
                     "opts": self.opts
                    }
        return json_dict

    def to_json_str(self):

        return json.dumps(self.get_json_dict(), sort_keys=True)

    def to_json(self, file_path):

        with open(file_path, 'w') as fp:
            json.dump(self.get_json_dict(), fp=fp, indent=2, separators=None,
                      sort_keys=True)

    @classmethod
    def from_json_str(cls, json_str):
        """ Initialise OBSolve from JSON string. """
        json_dict = json.loads(json_str)
        return cls(**json_dict)

    @classmethod
    def from_json(cls, file_path):
        """ Initialise OBSolve from JSON file. """
        with open(file_path) as json_file:
            json_dict = json.load(json_file)

            # This needs to be here to get the savefile name from json
            if 'savefile' not in json_dict:

                logging.debug('The savefile JSON element is missing.')

                savefile = os.path.splitext(file_path)[0]
                json_dict['savefile'] = savefile

            elif not json_dict['savefile']:

                logging.debug('The savefile JSON element is null.')

        return cls(**json_dict)
