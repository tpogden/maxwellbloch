# -*- coding: utf-8 -*-

import json

from maxwellbloch import ob_base

class OBAtom(ob_base.OBBase):

    def __init__(self, num_states=0, energies=[], decays=[], fields=[]):

        self.num_states = num_states
        self.energies = energies
        self.decays = decays

        self.build_fields(fields)

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