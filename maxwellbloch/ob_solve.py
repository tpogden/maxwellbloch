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
        return ("OBSolve(atom={0}, " +
                "t_min={1}, " +
                "t_max={2}, " +
                "t_steps={3}, " +
                "method={4}, " +
                "opts={5})").format(self.atom,
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

    def build_opts(self, opts):
        """ This currently just sets the options to default.
            Issue #96.
        """

        self.opts = qu.Options()
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
    def solve(self, rho0=None, e_ops=[], opts=qu.Options(), recalc=True,
              show_pbar=False, save=True):

        # When we're calling from MBSolve, we don't want to save each step.
        # So we pass in save=False.
        if save:
            savefile = self.savefile
        else:
            savefile = None

        if self.method == 'mesolve':
            self.atom.mesolve(self.tlist, rho0=rho0, e_ops=e_ops,
                                 opts=opts, recalc=recalc,
                                 savefile=savefile, show_pbar=show_pbar)

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
                     "opts": '{}'  #  TODO fix when opts fixed.
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


def parse_args():

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--file',
                        help='path of a JSON file containing an ob_solve \
                              problem definition', required=False)

    args = vars(parser.parse_args())

    if args['file']:

        print('Loading problem definition from file {0}'.format(args['file']))

        ob_solve_obj = OBSolve().from_json(args['file'])
        ob_solve_obj.solve()

def main():

    parse_args()

    return 0

if __name__ == '__main__':

    STATUS = main()
    sys.exit(STATUS)
