# -*- coding: utf-8 -*-

import json
from typing import Any

from maxwellbloch import t_funcs
from maxwellbloch.t_funcs import TFunc


class Field(object):
    """Field object to address the OBAtom object, describing the atomic levels
    coupled, detuning and Rabi frequency function.

    Attributes:
        label (string): a name for the field e.g. "probe".
        index (int): index of the field, to reference within an OBAtom.
        coupled_levels (list): pairs of levels coupled by the field.
            e.g. [[0,1], [0,2]]
        factors (list): List of strength factors for each pair in coupled
            levels.
        detuning (float): detuning of the fields from resonance with the
            coupled_levels transitions.
        detuning_positive (bool): is the detuning positive?
        counter_propagating (bool): if True, this field travels in the
            opposite direction to the forward-propagating solver direction.
            Sets factor_doppler_shift = -1.0. Orthogonal to detuning_positive.
        factor_doppler_shift (float): multiplier applied to the thermal
            detuning Delta for each velocity class. Default 1.0 (co-propagating).
            Set to -1.0 for a counter-propagating field, or to an arbitrary
            float for schemes where wavelength differences matter (e.g. a
            three-photon ladder where k-vector magnitudes differ).
        rabi_freq (float): Rabi frequency of the field on the transition.
        rabi_freq_t_func (func): Time-dependency of rabi_freq as function of
            time f(t, args)
        rabi_freq_t_args (dict): arguments to be passed to rabi_freq_t_func.
    """

    def __init__(
        self,
        label: str = "",
        index: int = 0,
        coupled_levels: list[list[int]] | None = None,
        factors: list[float] | None = None,
        detuning: float = 0.0,
        detuning_positive: bool = True,
        counter_propagating: bool = False,
        factor_doppler_shift: float | None = None,
        rabi_freq: float = 1.0,
        rabi_freq_t_func: str | None = None,
        rabi_freq_t_args: dict[str, float] | None = None,
    ) -> None:

        if coupled_levels is None:
            coupled_levels = []
        if rabi_freq_t_args is None:
            rabi_freq_t_args = {}

        self.label = label
        self.index = index

        self.coupled_levels = coupled_levels  # TODO should I convert to array?

        self._build_factors(factors)

        self.detuning = detuning
        self.detuning_positive = detuning_positive

        self._build_doppler_shift(counter_propagating, factor_doppler_shift)

        self.rabi_freq = rabi_freq
        self.rabi_freq_t_args = rabi_freq_t_args

        self.build_rabi_freq_t_func(rabi_freq_t_func=rabi_freq_t_func, index=index)
        self.build_rabi_freq_t_args(rabi_freq_t_args=rabi_freq_t_args, index=index)

    def __repr__(self):

        return (
            "Field(label={0}, "
            + "index={1}, "
            + "coupled_levels={2}, "
            + "factors={3}, "
            + "detuning={4}, "
            + "detuning_positive={5}, "
            + "counter_propagating={6}, "
            + "factor_doppler_shift={7}, "
            + "rabi_freq={8}, "
            + "rabi_freq_t_func={9}, "
            + "rabi_freq_t_args={10})"
        ).format(
            self.label,
            self.index,
            self.coupled_levels,
            self.factors,
            self.detuning,
            self.detuning_positive,
            self.counter_propagating,
            self.factor_doppler_shift,
            self.rabi_freq,
            self.rabi_freq_t_func,
            self.rabi_freq_t_args,
        )

    def _build_doppler_shift(
        self, counter_propagating: bool, factor_doppler_shift: float | None
    ) -> None:
        """Set counter_propagating and factor_doppler_shift with consistency check.

        counter_propagating=True is sugar for factor_doppler_shift=-1.0. Setting
        counter_propagating=True alongside an explicit positive factor_doppler_shift
        is a conflict and raises ValueError.
        """
        if (
            counter_propagating
            and factor_doppler_shift is not None
            and factor_doppler_shift > 0.0
        ):
            raise ValueError(
                "counter_propagating=True implies factor_doppler_shift < 0, "
                "but factor_doppler_shift={} was also supplied. "
                "Either set counter_propagating=True (to use -1.0) or set "
                "factor_doppler_shift directly, not both.".format(factor_doppler_shift)
            )
        self.counter_propagating = counter_propagating
        if counter_propagating:
            self.factor_doppler_shift = -1.0
        elif factor_doppler_shift is None:
            self.factor_doppler_shift = 1.0
        else:
            self.factor_doppler_shift = factor_doppler_shift

    def lower_levels(self) -> list[int]:
        """Return the unique lower (ground) level indices coupled by this field.

        Returns:
            Sorted list of unique indices ``c[0]`` from ``coupled_levels``.
        """
        return sorted(set(c[0] for c in self.coupled_levels))

    def upper_levels(self) -> list[int]:
        """Return the unique upper (excited) level indices coupled by this field.

        Returns:
            Sorted list of unique indices ``c[1]`` from ``coupled_levels``.
        """
        return sorted(set(c[1] for c in self.coupled_levels))

    def _build_factors(self, factors: list[float]) -> list[float]:
        """Builds the factors list.

        Args:
            factors (list). List of strength factors for each pair in
            coupled levels.
        Returns:
            self.factors (list)
        Notes:
            - There are no checks on what these factors are, or if they
                make any sense.
        """
        if not factors:
            factors = [1.0] * len(self.coupled_levels)
        if len(factors) != len(self.coupled_levels):
            raise ValueError(
                "The length of factors must be the same as the "
                "length of coupled_levels."
            )
        else:
            self.factors = factors
        return self.factors

    def build_rabi_freq_t_func(
        self, rabi_freq_t_func: str | None, index: int = 0
    ) -> TFunc:

        if rabi_freq_t_func:
            t_func = getattr(t_funcs, rabi_freq_t_func)
            self.rabi_freq_t_func = t_func(index)

        else:
            t_func = t_funcs.square
            self.rabi_freq_t_func = t_func(index)

        return self.rabi_freq_t_func

    def build_rabi_freq_t_args(
        self, rabi_freq_t_args: dict[str, float], index: int = 0
    ) -> dict[str, float]:

        self.rabi_freq_t_args = {}

        if rabi_freq_t_args:
            for key, value in rabi_freq_t_args.items():
                self.rabi_freq_t_args[key + "_" + str(index)] = value

        else:
            self.rabi_freq_t_args = {
                "on_" + str(index): 0.0,
                "off_" + str(index): 1.0,
                "ampl_" + str(index): 1.0,
            }

        return self.rabi_freq_t_args

    def get_json_dict(self) -> dict[str, Any]:
        """Return a dict representation of the Field object to be dumped to
        JSON.

        Note:
            For the rabi_freq_t_func attribute generated with
            build_rabi_freq_t_func, a suffix for the index will have been
            added. We remove that. e.g. e.g. ramp_onoff_0 -> ramp_onoff
        """

        json_dict = {
            "label": self.label,
            "index": self.index,
            "coupled_levels": self.coupled_levels,
            "factors": self.factors,
            "detuning": self.detuning,
            "detuning_positive": self.detuning_positive,
            "counter_propagating": self.counter_propagating,
            "factor_doppler_shift": self.factor_doppler_shift,
            "rabi_freq": self.rabi_freq,
            "rabi_freq_t_args": self.rabi_freq_t_args,
        }

        if self.rabi_freq_t_func:
            t_func_name = self.rabi_freq_t_func.__name__
            t_func_name = "_".join(t_func_name.split("_")[:-1])  # remove index
            json_dict.update({"rabi_freq_t_func": t_func_name})
        else:
            json_dict.update({"rabi_freq_t_func": None})

        if self.rabi_freq_t_args:
            rabi_freq_t_args = {}
            for key, value in self.rabi_freq_t_args.items():
                k = "_".join(key.split("_")[:-1])  # remove index
                rabi_freq_t_args[k] = value
            json_dict.update({"rabi_freq_t_args": rabi_freq_t_args})
        else:
            json_dict.update({"rabi_freq_t_args": {}})

        return json_dict

    def to_json_str(self) -> str:
        """Return a JSON string representation of the Field object.

        Returns:
            (string) JSON representation of the Field object.
        """

        return json.dumps(
            self.get_json_dict(), indent=2, separators=None, sort_keys=True
        )

    def to_json(self, file_path: str) -> None:

        with open(file_path, "w") as fp:
            json.dump(
                self.get_json_dict(), fp=fp, indent=2, separators=None, sort_keys=True
            )

    @classmethod
    def from_json_str(cls, json_str: str) -> "Field":
        json_dict = json.loads(json_str)
        return cls(**json_dict)

    @classmethod
    def from_json(cls, file_path: str) -> "Field":
        with open(file_path) as json_file:
            json_dict = json.load(json_file)
        return cls(**json_dict)
