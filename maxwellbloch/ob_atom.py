# -*- coding: utf-8 -*-

import json

import numpy as np
from numpy import pi
import qutip as qu

from maxwellbloch import ob_base, field, t_funcs


class OBAtom(ob_base.OBBase):

    def __init__(self, label=None, num_states=1, energies=[], decays=[], 
        fields=[]):
        """ Initialise OBAtom. 

        Args:
            label: string label describing the atom-field object.
            num_states: 
            energies: absolute or relative energy levels of the states.
            decays: list of dicts representing decays.
            e.g.
            [ { "rate": 1.0, "channels": [[0,1]] }
              { "rate": 2.0, "channels": [[2,1], [3,1]] } ]
            fields: list of Field objects that couple atom states.
        """

        super().__init__()

        self.label = label
        self.num_states = num_states
        self.energies = energies
        self.decays = decays

        self.build_fields(fields)
        self.build_operators()
        
        self.init_rho()

    def __repr__(self):
        return ("Atom(label={0}, " + "num_states={1}, " + "energies={2}, " +
            "decays={3}, " + "fields={4})").format(self.label, self.num_states,
            self.energies, self.decays, self.fields)

    def add_field(self, field_dict):
        """ Add a Field to the list given a dict representing the field. """ 

        self.fields.append(field.Field(**field_dict))

    def build_fields(self, field_dicts):
        """ Build the field list given a list of dicts representing fields. """

        self.fields = []
        for f in field_dicts:
            self.add_field(f)
        return self.fields

    def build_operators(self):
        """ Build the quantum operators representing the bare Hamiltonian,
            collapse operators and interaction Hamiltonian.
        """

        self.build_H_0()
        self.build_c_ops()
        self.build_H_Delta()
        self.build_H_Omega()

    def build_H_0(self):
        """ Makes a Bare Hamiltonian with the energies as diagonals.

        Leave the list empty for all zero energies (i.e. if you don't care
        about absolute energies.)

        Returns:
            H_0 = [energies[0]             0  ...]
                  [             0 energies[1] ...]
                  [           ...         ... ...]

        """

        if self.energies:
            H_0 = np.diag(np.array(self.energies))
        else:
            H_0 = np.zeros([self.num_states, self.num_states])

        self.H_0 = qu.Qobj(H_0)
        return self.H_0

    def build_c_ops(self):
        """ Takes the list of decays and makes a list of collapse operators to
            be passed to the solver.

        Notes:
            We want at least one collapse operator to force the master equation
            solver to produce density matrices, not state vectors. In the case
            that self.decays is empty, we'll add a zero collapse operator.
        """

        self.c_ops = []

        if not self.decays:
            self.c_ops.append(qu.Qobj(np.zeros([self.num_states,
                self.num_states])))
        else:
            for d in self.decays:
                r = d["rate"]
                for c in d["channels"]:
                    self.c_ops.append(np.sqrt(2*pi*r)*self.sigma(c[0], c[1]))

        return self.c_ops

    def build_H_Delta(self):

        self.H_Delta = qu.Qobj(np.zeros([self.num_states, self.num_states]))

        for f in self.fields:
            if f.detuning_positive:
                sgn = 1.
            else:
                sgn = -1.
            for c in f.coupled_levels:
                self.H_Delta -= sgn * f.detuning * self.sigma(c[1], c[1])

        return self.H_Delta

    def set_H_Delta(self, detunings):
        """ Set the detuning part of the interaction Hamiltonian, H_Delta,
            given a list of detunings.

        Args: 
            detunings: list of floats: detunings of each field in the list
        """

        assert(len(detunings) == len(self.fields))

        for i, f in enumerate(self.fields):
            f.detuning = detunings[i]

        return self.build_H_Delta()

    def shift_H_Delta(self, detuning_shifts):
        """ TODO: Could make use of set_H_Delta """

        for i, f in enumerate(self.fields):
            f.detuning = f.detuning + detuning_shifts[i]

        return self.build_H_Delta()

    def build_H_Omega(self):

        self.H_Omega_list = []

        H_Omega = qu.Qobj(np.zeros([self.num_states, self.num_states]))

        for f in self.fields:

            H_Omega = qu.Qobj(np.zeros([self.num_states, self.num_states]))

            for c in f.coupled_levels:
                H_Omega += self.sigma(c[0], c[1]) + self.sigma(c[1], c[0])

            H_Omega *= pi * f.rabi_freq  # 2π*rabi_freq/2

            if self.is_field_td():  # time-dependent interaction
                self.H_Omega_list.append([H_Omega, f.rabi_freq_t_func])

            else:  # time-independent
                self.H_Omega_list.append(H_Omega)

        return self.H_Omega_list

    def set_H_Omega(self, rabi_freqs, rabi_freq_t_funcs, rabi_freq_t_args):
        """
        Args:
            rabi_freqs: [floats]
            rabi_freq_t_funcs: [strings]
            rabi_freqs_t_args: [dicts]

        Returns:

        """

        for f_i, f in enumerate(self.fields):
            f.rabi_freq = rabi_freqs[f_i]

            f.build_rabi_freq_t_func(rabi_freq_t_funcs[f_i], f_i)
            f.build_rabi_freq_t_args(rabi_freq_t_args[f_i], f_i)

        return self.build_H_Omega()

    def init_rho(self):

        self.rho = self.ground_state() * self.ground_state().dag()

        return self.rho

    def get_field_args(self):

        args = {}
        for f in self.fields:
            args.update(f.rabi_freq_t_args)
        return args

    def is_field_td(self):

        # Time-dependent if there are any t_funcs specified
        return any(f.rabi_freq_t_func is not None for f in self.fields)

    def get_fields_sum_coherence(self):
        """ Returns the sum coherences of the atom density matrix for each
            field

        Returns:

        """

        sum_coh = np.zeros((len(self.fields), len(self.states_t())),
                           dtype=np.complex)
        for f_i, f in enumerate(self.fields):
            for cl in f.coupled_levels:
                sum_coh[f_i, :] += self.states_t()[:, cl[0], cl[1]]
        return sum_coh

    def mesolve(self, tlist, rho0=None, e_ops=[], opts=qu.Options(),
                recalc=True, savefile=None, show_pbar=False):

        args = self.get_field_args()

        td = self.is_field_td()

        self.result = super().mesolve(tlist=tlist, rho0=rho0, td=td,
                                      e_ops=e_ops, args=args, opts=opts,
                                      recalc=recalc, savefile=savefile,
                                      show_pbar=show_pbar)

        return self.result

    def get_json_dict(self):

        json_dict = {"label": self.label,
                     "num_states": self.num_states,
                     "energies": self.energies,
                     "decays": self.decays,
                     "fields": [f.get_json_dict() for f in self.fields]}
        return json_dict

    def to_json_str(self):
        """ Return a JSON string representation of the Atom object.

        Returns:
            (string) JSON representation of the Atom object.
        """

        return json.dumps(self.get_json_dict(), sort_keys=True)

    def to_json(self, file_path):

        with open(file_path, 'w') as fp:
            json.dump(self.get_json_dict(), fp=fp, indent=2, separators=None,
                      sort_keys=True)

    @classmethod
    def from_json_str(cls, json_str):
        json_dict = json.loads(json_str)
        return cls(**json_dict)

    @classmethod
    def from_json(cls, file_path):
        with open(file_path) as json_file:
            json_dict = json.load(json_file)
        return cls(**json_dict)


def main():

    print(OBAtom())
    return 0

if __name__ == '__main__':
    status = main()
