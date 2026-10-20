# -*- coding: utf-8 -*-

import copy
import json
from typing import Any

import numpy as np
import qutip as qu
from numpy import pi

from clerq import field, hyperfine, system_base


class System(system_base.SystemBase):
    """An open quantum system: its states, decays and the fields coupling them.

    The system can be an atom or any other few-level quantum system, such as
    a quantum dot or a colour centre.

    Args:
        label: string label describing the system.
        num_states: number of states.
        energies: absolute or relative energy levels of the states.
        decays: list of dicts representing decays.
            e.g.
            [ { "rate": 1.0, "channels": [[0,1]], "factors": [1.0] }
              { "rate": 2.0, "channels": [[2,1], [3,1]],
                "factors": [0.707, 0.707] } ]
        initial_state: initial state, as a list or array
        fields: list of dicts representing the fields that couple the
            system's states.
        b_field: a DC magnetic field along the quantisation axis, as a dict
            with the keys ``magnitude`` (default 0), ``units`` (``"natural"``,
            the default, or ``"gauss"``), and the per-level lists ``g_F`` and
            ``m_F`` (length ``num_states``). It adds ``g_F * m_F * B`` to
            the level energies (``H_0 += 2 pi g_F m_F B``), with the sign
            convention of Steck's alkali data (see ``clerq.hyperfine``).
            ``"natural"`` means ``mu_B * B / h`` in the same frequency unit as
            ``energies``; ``"gauss"`` assumes energies in MHz (time in
            microseconds). For hyperfine systems,
            ``Atom1e.get_zeeman_levels()`` gives ``g_F`` and ``m_F``.
    """

    def __init__(
        self,
        label: str | None = None,
        num_states: int = 1,
        energies: list[float] | None = None,
        decays: list[dict] | None = None,
        initial_state: list[float] | None = None,
        fields: list[dict] | None = None,
        b_field: dict | None = None,
    ) -> None:

        super().__init__()

        self.label = label
        self.num_states = num_states
        self.energies = energies if energies is not None else []
        self.decays = decays if decays is not None else []
        self._b_field = self._check_b_field(b_field)

        self.build_initial_state(initial_state if initial_state is not None else [])

        self.build_fields(fields if fields is not None else [])
        self.build_operators()

    def __repr__(self):
        return (
            "System(label={0}, "
            + "num_states={1}, "
            + "energies={2}, "
            + "decays={3}, "
            + "fields={4})"
        ).format(self.label, self.num_states, self.energies, self.decays, self.fields)

    @property
    def b_field(self) -> dict | None:
        """The validated ``b_field`` dict (a copy: assign a new dict to change it)."""

        return copy.deepcopy(self._b_field)

    @b_field.setter
    def b_field(self, value: dict | None) -> None:
        self._b_field = self._check_b_field(value)
        self._build_H_0()

    def _check_b_field(self, b_field: dict | None) -> dict | None:
        """Validate the ``b_field`` dict and fill in its defaults."""

        if b_field is None:
            return None
        unknown = set(b_field) - {"magnitude", "units", "g_F", "m_F"}
        if unknown:
            raise ValueError(f"Unknown b_field key(s) {sorted(unknown)}.")
        units = b_field.get("units", "natural")
        if units not in ("natural", "gauss"):
            raise ValueError("b_field units must be 'natural' or 'gauss'.")
        for key in ("g_F", "m_F"):
            if key not in b_field:
                raise ValueError(f"b_field requires a per-level '{key}' list.")
            if len(b_field[key]) != self.num_states:
                raise ValueError(f"b_field '{key}' must have num_states elements.")
        return {
            "magnitude": float(b_field.get("magnitude", 0.0)),
            "units": units,
            "g_F": [float(g) for g in b_field["g_F"]],
            "m_F": [float(m) for m in b_field["m_F"]],
        }

    @property
    def b_field_natural(self) -> float:
        """The B-field magnitude in natural units (``mu_B B / h``); 0 if none."""

        if self._b_field is None:
            return 0.0
        if self._b_field["units"] == "gauss":
            return hyperfine.gauss_to_natural(self._b_field["magnitude"])
        return self._b_field["magnitude"]

    @property
    def zeeman_coeffs(self) -> np.ndarray:
        """``g_F * m_F`` of each level: the Zeeman shift per unit natural B."""

        if self._b_field is None:
            return np.zeros(self.num_states)
        return np.array(self._b_field["g_F"]) * np.array(self._b_field["m_F"])

    def build_initial_state(self, initial_state: list[float] | None = None) -> qu.Qobj:
        """Build the initial density matrix for the system.

        The default is for all of the population to be in ``|0><0|``.

        Args:
            initial_state: a list or array of populations, length num_states.
                If not given, all the population starts in level 0.

        Returns:
            The initial density matrix, diagonal with the given populations.
            It is also stored as ``self.initial_state`` and ``self.rho``.

        Raises:
            ValueError: If ``initial_state`` has the wrong length or its
                elements do not sum to 1.

        Notes:
            - The elements must sum to 1 (as this will be the trace of the
                density matrix).

        """

        if initial_state:
            if len(initial_state) != self.num_states:
                raise ValueError("initial_state must have num_states elements.")
            self.initial_state = qu.qzero(self.num_states)
            for i, g in enumerate(initial_state):
                self.initial_state += g * self.sigma(a=i, b=i)
        else:
            self.initial_state = self.sigma(a=0, b=0)
        if not np.isclose(self.initial_state.tr(), 1.0):
            raise ValueError("initial_state must have diagonal elements sum to 1.")

        # Initialise rho
        self.rho = self.initial_state

        return self.initial_state

    # def init_rho(self):

    #     self.rho = self.rho_0
    #     return self.rho

    def add_field(self, field_dict: dict) -> None:
        """Add a Field to the list given a dict representing the field."""

        self.fields.append(field.Field(**field_dict))

    def build_fields(self, field_dicts: list[dict]) -> list[field.Field]:
        """Build the field list given a list of dicts representing fields."""

        self.fields = []
        for f_idx, f in enumerate(field_dicts):
            f["index"] = f_idx
            self.add_field(f)
        return self.fields

    def build_operators(self) -> None:
        """Build the quantum operators representing the bare Hamiltonian,
        collapse operators and interaction Hamiltonian.
        """

        self._build_H_0()
        self._build_c_ops()
        self._build_H_Delta()
        self._build_H_Omega()
        self._build_sum_coherence_arrays()

    def _build_H_0(self) -> qu.Qobj:
        """Makes a Bare Hamiltonian with the energies as diagonals.

        Leave the list empty for all zero energies (i.e. if you don't care
        about absolute energies.)

        Returns:
            H_0 = [energies[0]             0  ...]
                  [             0 energies[1] ...]
                  [           ...         ... ...]

        """

        if self.energies:
            H_0 = np.diag(2 * pi * np.array(self.energies))
        else:
            H_0 = np.zeros([self.num_states, self.num_states])

        if self._b_field is not None:
            H_0 = H_0 + np.diag(2 * pi * self.b_field_natural * self.zeeman_coeffs)

        self.H_0 = qu.Qobj(H_0)
        return self.H_0

    def _build_c_ops(self) -> list[qu.Qobj]:
        """Takes the list of decays and makes a list of collapse operators to
            be passed to the solver.

        Notes:
            We want at least one collapse operator to force the master equation
            solver to produce density matrices, not state vectors. In the case
            that self.decays is empty, we'll add a zero collapse operator.
        """
        self.c_ops = []
        if not self.decays:
            self.c_ops.append(qu.Qobj(np.zeros([self.num_states, self.num_states])))
        else:
            for d in self.decays:
                # NOTE: This means if there are no decays factors in the JSON
                # input, factors will be added to any JSON output.
                if "factors" not in d:
                    d["factors"] = [1.0] * len(d["channels"])
                if len(d["factors"]) != len(d["channels"]):
                    raise ValueError(
                        "Length of factors list ({}) is not the "
                        "same as length of channels list ({})".format(
                            len(d["factors"]), len(d["channels"])
                        )
                    )
                for c_i, c in enumerate(d["channels"]):
                    self.c_ops.append(
                        d["factors"][c_i]
                        * np.sqrt(2 * pi * d["rate"])
                        * self.sigma(a=c[0], b=c[1])
                    )
        return self.c_ops

    def _build_H_Delta(self) -> qu.Qobj:
        """Builds the detuning part of the interaction Hamiltonian."""

        self.H_Delta = qu.Qobj(np.zeros([self.num_states, self.num_states]))
        for f in self.fields:
            if f.detuning_positive:
                sgn = 1.0
            else:
                sgn = -1.0
            for i in f.upper_levels():
                self.H_Delta -= sgn * 2 * pi * f.detuning * self.sigma(a=i, b=i)
        return self.H_Delta

    def set_H_Delta(self, detunings: list[float]) -> qu.Qobj:
        """Set the detuning part of the interaction Hamiltonian, H_Delta,
            given a list of detunings.

        Args:
            detunings: list of floats: detunings of each field in the list
        """
        assert len(detunings) == len(self.fields)
        for i, f in enumerate(self.fields):
            f.detuning = detunings[i]
        return self._build_H_Delta()

    def _build_H_Omega(self) -> list:
        """Builds the Rabi frequency (off-diagonals) part of the interaction
        Hamiltonian.
        """

        self.H_Omega_list = []
        for f in self.fields:
            H_Omega = sum(
                (self.sigma(a=c[0], b=c[1]) + self.sigma(a=c[1], b=c[0])) * fac
                for c, fac in zip(f.coupled_levels, f.factors, strict=True)
            )
            H_Omega *= pi * f.rabi_freq  # 2π*rabi_freq/2
            if self.is_field_td():  # time-dependent interaction
                self.H_Omega_list.append([H_Omega, f.rabi_freq_t_func])
            else:  # time-independent
                self.H_Omega_list.append(H_Omega)
        return self.H_Omega_list

    def set_H_Omega(
        self,
        rabi_freqs: list[float],
        rabi_freq_t_funcs: list,
        rabi_freq_t_args: list[dict],
    ) -> list:
        """Set the Rabi frequency part of the interaction Hamiltonian.

        Updates each field's Rabi frequency and time function, then rebuilds
        the Hamiltonian.

        Args:
            rabi_freqs: Rabi frequency for each field.
            rabi_freq_t_funcs: Name of the time function for each field, as
                in ``clerq.t_funcs`` (e.g. ``"gaussian"``).
            rabi_freq_t_args: Argument dict for each field's time function.

        Returns:
            The list of Rabi frequency Hamiltonian terms, one for each field.
        """

        for f_i, f in enumerate(self.fields):
            f.rabi_freq = rabi_freqs[f_i]

            f.build_rabi_freq_t_func(rabi_freq_t_func=rabi_freq_t_funcs[f_i], index=f_i)
            f.build_rabi_freq_t_args(rabi_freq_t_args=rabi_freq_t_args[f_i], index=f_i)

        return self._build_H_Omega()

    def get_detunings(self) -> list[float]:
        """Returns a list of detunings, one for each field in fields."""
        return [f.detuning for f in self.fields]

    def get_field_args(self) -> dict[str, float]:

        args = {}
        for f in self.fields:
            args.update(f.rabi_freq_t_args)
        return args

    def is_field_td(self) -> bool:

        # Time-dependent if there are any t_funcs specified
        return any(f.rabi_freq_t_func is not None for f in self.fields)

    def _build_sum_coherence_arrays(self) -> None:
        """Precompute index and weight arrays used by get_fields_sum_coherence.

        Called once during build_operators so that get_fields_sum_coherence
        can use vectorised NumPy indexing rather than Python loops at every
        z-step.  For each field we store:
            _sum_coh_rows[f_i]    – row indices into the density matrix
            _sum_coh_cols[f_i]    – column indices into the density matrix
            _sum_coh_weights[f_i] – corresponding factors (signed dipole
                                    matrix elements, the same weights as in
                                    H_Omega)
        """
        self._sum_coh_rows: list[np.ndarray] = []
        self._sum_coh_cols: list[np.ndarray] = []
        self._sum_coh_weights: list[np.ndarray] = []
        for f in self.fields:
            self._sum_coh_rows.append(np.array([cl[0] for cl in f.coupled_levels]))
            self._sum_coh_cols.append(np.array([cl[1] for cl in f.coupled_levels]))
            self._sum_coh_weights.append(np.asarray(f.factors, dtype=float))

    def get_fields_sum_coherence(
        self, states_t: np.ndarray | None = None
    ) -> np.ndarray:
        """Returns the sum coherences of the system density matrix for each
            field, including the relative strength factors.

        Notes:
            The coherence of each coupled pair is weighted by the pair's
            (signed) factor, i.e. by its dipole matrix element, the same
            factor that multiplies the pair in the interaction Hamiltonian.
            The field's polarisation is then
            ``sum_k factor_k * rho_{a_k b_k}``, and in the weak-probe limit
            the response goes as the line strength ``sum_k factor_k**2``.

        Args:
            states_t: (optional) a states_t object. If not provided, use
                self.states_t()

        Returns:
            np.ndarray of shape (num_fields, t_steps+1)
        """
        if len(self._sum_coh_rows) != len(self.fields):
            self._build_sum_coherence_arrays()
        if states_t is None:
            states_t = self.states_t()
        sum_coh = np.zeros((len(self.fields), len(states_t)), dtype=complex)
        for f_i in range(len(self.fields)):
            # states_t[:, rows, cols] → (t_steps+1, n_pairs)
            # @ weights              → (t_steps+1,)
            sum_coh[f_i] = (
                states_t[:, self._sum_coh_rows[f_i], self._sum_coh_cols[f_i]]
                @ self._sum_coh_weights[f_i]
            )
        return sum_coh

    def mesolve(
        self,
        tlist: np.ndarray,
        e_ops: list | None = None,
        options: dict | None = None,
        recalc: bool = True,
        savefile: str | None = None,
        show_pbar: bool = False,
    ) -> Any:

        if e_ops is None:
            e_ops = []
        if options is None:
            options = {}

        args = self.get_field_args()

        td = self.is_field_td()

        self.result = super().mesolve(
            tlist=tlist,
            rho0=self.initial_state,
            td=td,
            e_ops=e_ops,
            args=args,
            options=options,
            recalc=recalc,
            savefile=savefile,
            show_pbar=show_pbar,
        )

        return self.result

    def get_json_dict(self) -> dict[str, Any]:

        json_dict = {
            "label": self.label,
            "num_states": self.num_states,
            "energies": self.energies,
            "decays": self.decays,
            "fields": [f.get_json_dict() for f in self.fields],
        }
        if self.b_field is not None:
            json_dict["b_field"] = self.b_field
        return json_dict

    def to_json_str(self) -> str:
        """Return a JSON string representation of the System object.

        Returns:
            (string) JSON representation of the System object.
        """

        return json.dumps(self.get_json_dict(), sort_keys=True)

    def to_json(self, file_path: str) -> None:

        with open(file_path, "w") as fp:
            json.dump(
                self.get_json_dict(), fp=fp, indent=2, separators=None, sort_keys=True
            )

    @classmethod
    def from_json_str(cls, json_str: str) -> "System":
        json_dict = json.loads(json_str)
        return cls(**json_dict)

    @classmethod
    def from_json(cls, file_path: str) -> "System":
        with open(file_path) as json_file:
            json_dict = json.load(json_file)
        return cls(**json_dict)
