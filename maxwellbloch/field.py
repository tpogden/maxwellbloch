# -*- coding: utf-8 -*-

import json
from maxwellbloch import t_funcs

class Field(object):
    """ Field object to address the OBAtom object, describing the atomic levels
    coupled, detuning and Rabi frequency function.

    Attributes:
        label (string): a name for the field e.g. "probe"
        
        coupled_levels (list): pairs of levels coupled by the field.
            e.g. [[0,1], [0,2]]

        detuning (float): detuning of the fields from resonance with the 
            coupled_levels transitions.
        detuning_positive (bool): is the detuning positive?

        rabi_freq (float): Rabi frequency of the field on the transition.
        rabi_freq_t_func (func): Time-dependency of rabi_freq as function of 
            time f(t, args)
        rabi_freq_t_args (dict): arguments to be passed to rabi_freq_t_func.
    """        

    def __init__(self, label="", coupled_levels=[], detuning=0.0,
                 detuning_positive=True, rabi_freq=0.0,
                 rabi_freq_t_func=None, rabi_freq_t_args={}):
        
        self.label = label
        
        self.coupled_levels = coupled_levels
        
        self.detuning = detuning
        self.detuning_positive = detuning_positive

        self.rabi_freq = rabi_freq
        # self.rabi_freq_t_func = rabi_freq_t_func
        self.rabi_freq_t_args = rabi_freq_t_args

        self.build_rabi_freq_t_func(rabi_freq_t_func)

    def __repr__(self):

        return ("Field(label={0}, " +
                "coupled_levels={1}, " +
                "detuning={2}, " +
                "detuning_positive={3}, "
                "rabi_freq={4}, " +
                "rabi_freq_t_func={5}, " +
                "rabi_freq_t_args={6})").format(self.label, 
                                               self.coupled_levels, 
                                               self.detuning,
                                               self.detuning_positive,
                                               self.rabi_freq,
                                               self.rabi_freq_t_func,
                                               self.rabi_freq_t_args)

    def build_rabi_freq_t_func(self, rabi_freq_t_func):

        if rabi_freq_t_func:
            self.rabi_freq_t_func = getattr(t_funcs, rabi_freq_t_func)
        else:
            self.rabi_freq_t_func = None

        return self.rabi_freq_t_func

    def get_json_dict(self):

        json_dict = { "label": self.label,
                      "coupled_levels": self.coupled_levels,
                      "detuning": self.detuning,
                      "detuning_positive": self.detuning_positive,
                      "rabi_freq": self.rabi_freq,
                      "rabi_freq_t_args": self.rabi_freq_t_args }

        if self.rabi_freq_t_func:
            json_dict.update({"rabi_freq_t_func": 
                             self.rabi_freq_t_func.__name__})
        else:
            json_dict.update({"rabi_freq_t_func": None})

        return json_dict

    def to_json_str(self):
        """ Return a JSON string representation of the Field object.

        Returns:
            (string) JSON representation of the Field object.
        """

        return json.dumps(self.get_json_dict(), indent=2, separators=None,
                          sort_keys=True)

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

    print(Field(rabi_freq_t_func="square_2"))

if __name__ == '__main__':
    status = main()