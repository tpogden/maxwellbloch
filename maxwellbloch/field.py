# -*- coding: utf-8 -*-

import json
from maxwellbloch import t_funcs

class Field(object):
    """ Field object to address the OBAtom object, describing the atomic levels
    coupled, detuning and Rabi frequency function.

    Attributes:
        label (string): a name for the field e.g. "probe"
        TODO: index
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

    def __init__(self, label="", index=0, coupled_levels=[], detuning=0.0,
                 detuning_positive=True, rabi_freq=0.0,
                 rabi_freq_t_func=None, rabi_freq_t_args={}):

        self.label = label
        self.index = index

        self.coupled_levels = coupled_levels

        self.detuning = detuning
        self.detuning_positive = detuning_positive

        self.rabi_freq = rabi_freq
        self.rabi_freq_t_args = rabi_freq_t_args

        self.build_rabi_freq_t_func(rabi_freq_t_func, index)
        self.build_rabi_freq_t_args(rabi_freq_t_args, index)

    def __repr__(self):

        return ("Field(label={0}, " +
                "index={1}, " +
                "coupled_levels={2}, " +
                "detuning={3}, " +
                "detuning_positive={4}, "
                "rabi_freq={5}, " +
                "rabi_freq_t_func={6}, " +
                "rabi_freq_t_args={7})").format(self.label,
                                                self.index,
                                                self.coupled_levels,
                                                self.detuning,
                                                self.detuning_positive,
                                                self.rabi_freq,
                                                self.rabi_freq_t_func,
                                                self.rabi_freq_t_args)

    def build_rabi_freq_t_func(self, rabi_freq_t_func, index=0):

        if rabi_freq_t_func:

            t_func = getattr(t_funcs, rabi_freq_t_func)
            self.rabi_freq_t_func = t_func(index)

        else:
            t_func = t_funcs.square
            self.rabi_freq_t_func = t_func(index)

        return self.rabi_freq_t_func

    def build_rabi_freq_t_args(self, rabi_freq_t_args, index=0):

        self.rabi_freq_t_args = {}

        if rabi_freq_t_args:
            for key, value in rabi_freq_t_args.items():
                self.rabi_freq_t_args[key + '_' + str(index)] = value
        
        else:
            self.rabi_freq_t_args = {'on_' + str(index): 0.0, 
                                     'off_' + str(index): 1.0,
                                     'ampl_' + str(index): 1.0}

        return self.rabi_freq_t_args


    def get_json_dict(self):
        """ Return a dict representation of the Field object to be dumped to
            JSON.

            Note:
                For the rabi_freq_t_func attribute generated with 
                build_rabi_freq_t_func, a suffix for the index will have been 
                added. We remove that. e.g. e.g. ramp_onoff_0 -> ramp_onoff
        """

        json_dict = {"label": self.label,
                     "index": self.index,
                     "coupled_levels": self.coupled_levels,
                     "detuning": self.detuning,
                     "detuning_positive": self.detuning_positive,
                     "rabi_freq": self.rabi_freq,
                     "rabi_freq_t_args": self.rabi_freq_t_args}

        if self.rabi_freq_t_func:
            t_func_name = self.rabi_freq_t_func.__name__
            t_func_name = '_'.join(t_func_name.split('_')[:-1])  # remove index
            json_dict.update({"rabi_freq_t_func": t_func_name})
        else:
            json_dict.update({"rabi_freq_t_func": None})

        if self.rabi_freq_t_args:
            rabi_freq_t_args = {}
            for key, value in self.rabi_freq_t_args.items():
                k = '_'.join(key.split('_')[:-1])  # remove index
                rabi_freq_t_args[k] = value
            json_dict.update({"rabi_freq_t_args": rabi_freq_t_args})
        else:
            json_dict.update({"rabi_freq_t_args": {}})

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
