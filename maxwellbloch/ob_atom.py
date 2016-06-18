# -*- coding: utf-8 -*-

import json

from maxwellbloch import ob_base, field

class OBAtom(ob_base.OBBase):

    def __init__(self, num_states=0, energies=[], decays=[], fields=[]):

        self.num_states = num_states
        self.energies = energies
        self.decays = decays

        self.build_fields(fields)

        self.build_H_0(energies)
        self.build_c_ops(decays)


    def __repr__(self):
        return ("Atom(num_states={0}, " +
                "energies={1}, " +
                "decays={2}, " +
                "fields={3})").format(self.num_states,
                                    self.energies,
                                    self.decays,
                                    self.fields)

    def add_field(self, field_dict):
        self.fields.append(field.Field(**field_dict))

    def build_fields(self, field_dicts):
        self.fields = []
        for f in field_dicts:
            self.add_field(f)
        return self.fields

    def build_H_0(self, energies=[]): 
        """ Takes a list of energies and makes a Bare Hamiltonian with the
        energies as diagonals.

        Leave the list empty for all zero energies (i.e. if you don't care
        about absolute energies.)

        Args:
            energies: list of pre-interaction energies

        Returns:
            H_0 = [energies[0]             0  ...]
                  [             0 energies[1] ...]
                  [           ...         ... ...]

        """

        if energies:
            H_0 = np.diag(np.array(energies))
        else:
            H_0 = np.zeros([self.num_states, self.num_states])

        self.H_0 = qu.Qobj(H_0)
        return self.H_0

    def build_c_ops(self, decays=[]):
        """ Takes a list of spontaneous decay rates and makes a list of
        collapse operators to be passed to the solver.

        Args: 
            decays: list of dicts representing decays.
            e.g.
            [ { "rate": 1.0, "channels": [[0,1]] }
              { "rate": 2.0, "channels": [[2,1], [3,1]] } ]
        """

        self.c_ops = []

        for d in decays:
            r = d["rate"]
            for c in d["channels"]:
                self.c_ops.append(np.sqrt(r)*self.sigma(c[0],c[1]))
        return self.c_ops

    def get_json_dict(self):

        json_dict = { "num_states": self.num_states,
                      "energies": self.energies,
                      "decays": self.decays,
                      "fields": [f.__dict__ for f in self.fields] }
        return json_dict

    def to_json_str(self):
        """ Return a JSON string representation of the Atom object.

        Returns:
            (string) JSON representation of the Atom object.
        """

        return json.dumps(self.get_json_dict())

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

if __name__ == '__main__':
    status = main()