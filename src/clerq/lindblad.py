# -*- coding: utf-8 -*-

import json
import logging
import os
from typing import Any

import numpy as np

from clerq._deprecation import rename_atom_key, resolve_system_arg, warn_deprecated
from clerq.system import System


class LindbladSolver(object):
    """Time-domain master equation solver for a single spatial point.

    Wraps :func:`qutip.mesolve` to evolve the density matrix of a
    :class:`~clerq.system.System` over a user-defined time grid, and can
    save/load results to avoid recomputation.

    Args:
        system: Dict (or empty dict for defaults) describing the quantum
            system; passed to :class:`~clerq.system.System`.
        t_min: Start time of the simulation.
        t_max: End time of the simulation.
        t_steps: Number of time steps.
        method: Solver method — ``'mesolve'`` (Lindblad master equation) is
            the only supported value.
        opts: Dict of QuTiP solver options passed to ``qutip.Options``.
        savefile: Path prefix for caching solved results (``'.qu'`` extension
            appended automatically). Pass ``None`` to disable caching.
        atom: Deprecated name for ``system``.
    """

    def __init__(
        self,
        system: dict | None = None,
        t_min: float = 0.0,
        t_max: float = 1.0,
        t_steps: int = 100,
        method: str = "mesolve",
        opts: dict | None = None,
        savefile: str | None = None,
        *,
        atom: dict | None = None,
    ) -> None:

        system = resolve_system_arg(system, atom)
        if system is None:
            system = {}

        self._build_system(system)

        self._build_tlist(t_min=t_min, t_max=t_max, t_steps=t_steps)

        self.method = method
        self.build_opts(opts)

        self._build_savefile(savefile)

    def __repr__(self):
        return (
            "LindbladSolver(system={0}, "
            + "t_min={1}, "
            + "t_max={2}, "
            + "t_steps={3}, "
            + "method={4}, "
            + "opts={5})"
        ).format(
            self.system, self.t_min, self.t_max, self.t_steps, self.method, self.opts
        )

    def _build_system(self, system_dict: dict) -> System:

        self.system = System(**system_dict)
        return self.system

    @property
    def atom(self) -> System:
        """Deprecated name for ``system``."""
        warn_deprecated("The 'atom' attribute", "'system'")
        return self.system

    @atom.setter
    def atom(self, value: System) -> None:
        warn_deprecated("The 'atom' attribute", "'system'")
        self.system = value

    def _build_tlist(self, t_min: float, t_max: float, t_steps: int) -> np.ndarray:

        from numpy import linspace

        self.t_min = t_min
        self.t_max = t_max
        self.t_steps = t_steps

        self.tlist = linspace(t_min, t_max, t_steps + 1)
        return self.tlist

    def build_opts(self, opts: dict | None = None) -> dict:
        """Build the options dict to be passed into the QuTiP solver.

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
            "atol": 1e-8,
            "rtol": 1e-6,
            "method": "adams",
            "nsteps": 1000,
            "max_step": self.t_step(),
        }
        if opts:
            self.opts.update(opts)  # Update any specified in the parameter
        return self.opts

    def set_field_rabi_freq_t_func(self, field_idx: int, t_func: Any) -> None:
        """Set the Rabi frequency time function of a field to a new time
            function. This is useful when you want to set a custom function,
            not one available in t_funcs.py

        Args:
            field_idx: The field for which to set the Rabi frequency t_func
            t_func: Rabi frequency as a function of time, f(t, args)
        """

        self.system.fields[field_idx].rabi_freq_t_func = t_func
        self.system.build_operators()  # Rebuild so H_Omega is updated

    def set_field_rabi_freq_t_args(self, field_idx: int, t_args: dict) -> None:
        """Set the Rabi frequency time function arguments. To be used with
            set_field_rabi_freq_t_func

        Args:
            field_idx: The field for which to set the Rabi frequency t_args
            t_args: A dict representing the args to go with the t_func.
        """

        self.system.fields[field_idx].rabi_freq_t_args = t_args

    def solve(
        self,
        e_ops: list | None = None,
        opts: dict | None = None,
        recalc: bool = True,
        show_pbar: bool = False,
        save: bool = True,
    ) -> np.ndarray:
        """Solves the Lindblad master equation for the system.

        Args:
            e_ops: Expectation value operators passed to the QuTiP solver.
            opts: QuTiP solver options, overriding those set at
                initialisation.
            recalc: Recalculate the solution even if a savefile exists?
            show_pbar: Show the QuTiP text progress bar.
            save: Save the result to ``savefile``, if one is set.

        Returns:
            The solved density matrix at each time, shape
            (t_steps+1, num_states, num_states).
        """
        return self._solve_system(
            e_ops=e_ops, opts=opts, recalc=recalc, show_pbar=show_pbar, save=save
        )

    def obsolve(self, *args: Any, **kwargs: Any) -> np.ndarray:
        """Deprecated name for :meth:`LindbladSolver.solve`."""
        warn_deprecated("obsolve()", "solve()")
        return self._solve_system(*args, **kwargs)

    def _solve_system(
        self,
        e_ops: list | None = None,
        opts: dict | None = None,
        recalc: bool = True,
        show_pbar: bool = False,
        save: bool = True,
    ) -> np.ndarray:
        """Solves the master equation at a single spatial point.

        Kept separate from :meth:`solve` so that a subclass can override
        ``solve`` and still use this at each step.
        """

        # When we're calling from PropagationSolver, we don't want to save each step.
        # So we pass in save=False.
        if save:
            savefile = self.savefile
        else:
            savefile = None

        # Choosing to overwrite opts at the solve stage.
        if opts:
            self.build_opts(opts)
        options = self.opts

        if self.method == "mesolve":
            self.system.mesolve(
                self.tlist,
                e_ops=e_ops,
                options=options,
                recalc=recalc,
                savefile=savefile,
                show_pbar=show_pbar,
            )

        return self.system.states_t()  # self.system.result

    def states_t(self) -> np.ndarray:

        return self.system.states_t()

    def t_step(self) -> float:

        return (self.t_max - self.t_min) / self.t_steps

    def _build_savefile(self, savefile: str | None) -> None:

        self.savefile = savefile

    def savefile_exists(self) -> bool:
        """Returns true if savefile (with appended extension .qu) exists."""

        return os.path.isfile(str(self.savefile) + ".qu")

    ### JSON Serialise and Deserialise

    def get_json_dict(self) -> dict[str, Any]:

        json_dict = {
            "system": self.system.get_json_dict(),
            "t_min": self.t_min,
            "t_max": self.t_max,
            "t_steps": self.t_steps,
            "method": self.method,
            "opts": self.opts,
        }
        return json_dict

    def to_json_str(self) -> str:

        return json.dumps(self.get_json_dict(), sort_keys=True)

    def to_json(self, file_path: str) -> None:

        with open(file_path, "w") as fp:
            json.dump(
                self.get_json_dict(), fp=fp, indent=2, separators=None, sort_keys=True
            )

    @classmethod
    def from_json_str(cls, json_str: str) -> "LindbladSolver":
        """Initialise LindbladSolver from JSON string."""
        json_dict = rename_atom_key(json.loads(json_str))
        return cls(**json_dict)

    @classmethod
    def from_json(cls, file_path: str) -> "LindbladSolver":
        """Initialise LindbladSolver from JSON file."""
        with open(file_path) as json_file:
            json_dict = json.load(json_file)

            # This needs to be here to get the savefile name from json
            if "savefile" not in json_dict:
                logging.debug("The savefile JSON element is missing.")

                savefile = os.path.splitext(file_path)[0]
                json_dict["savefile"] = savefile

            elif not json_dict["savefile"]:
                logging.debug("The savefile JSON element is null.")

        return cls(**rename_atom_key(json_dict))
