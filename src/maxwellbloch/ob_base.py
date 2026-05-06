# -*- coding: utf-8 -*-

import os
from typing import Any

import numpy as np
import qutip as qu

from maxwellbloch import sigma

# QuTiP 5 changed the default coefficient function style to 'pythonic'
# (f(t, **kwargs)). This codebase uses the QuTiP 4 'dict' style (f(t, args)).
qu.settings.core["function_coefficient_style"] = "dict"


class OBBase(object):
    """Base class providing the QuTiP interface for density-matrix solvers.

    Constructs and stores the Hamiltonian components and collapse operators
    that are shared by :class:`OBSolve` and :class:`MBSolve`. Not intended
    to be instantiated directly.

    Attributes:
        num_states: Number of atomic states in the system.
        rho: Initial density matrix (QuTiP Qobj).
        H_0: Bare (field-free) Hamiltonian.
        c_ops: List of collapse operators representing decay channels.
        H_Delta: Detuning part of the interaction Hamiltonian.
        H_Omega_list: List of Rabi-frequency interaction Hamiltonians, one
            per field.
        result: QuTiP Result object populated after solving.
    """

    def __init__(self):
        """Initialise the OB object. This will be overridden in every child
        class, just here to make it clear what attributes an OB object has."""

        self.num_states = 0
        self.rho = qu.Qobj()

        self.H_0 = qu.Qobj()
        self.c_ops = []

        # These will be set by functions in the child class.
        self.H_Delta = qu.Qobj()
        self.H_Omega_list = [qu.Qobj()]

        self.result = None

    def sigma(self, a: int, b: int) -> qu.Qobj:
        """The transition or 'flip' operator between states |a> and |b>, given
        by |a><b|. Uses the sigma_n helper function below this class.

        Args:
            a: level transitioned from
            b: level transitioned to

        Returns:
            |a><b|
        """

        return sigma.sigma(n=self.num_states, a=a, b=b)

    def set_H_0(self, energies: list[float] | None = None) -> qu.Qobj:
        """Takes a list of energies and makes a Bare Hamiltonian with the
        energies as diagonals. This function can be overridden in a child class
        if you want to make the bare Hamiltonian a different way.

        Leave the list empty for all zero energies (i.e. if you don't care
        about absolute energies.

        Args:
            energies: list of pre-interaction energies

        Returns
            H_0 = [energies[0]             0  ...]
                  [             0 energies[1] ...]
                  [           ...         ... ...]

        """

        # See PEP 8: 'Use the fact that empty sequences are false.'
        if not energies:
            H_0 = np.zeros([self.num_states, self.num_states])
        else:
            H_0 = np.diag(np.array(energies))

        self.H_0 = qu.Qobj(H_0)
        return self.H_0

    def H_I_list(self) -> list:
        """The interaction terms* are specified as a list. This could be a
        single item list in the case of one laser, or one item for each laser
        in the case of multiple lasers, or a laser and a dipole-dipole
        interaction in the case of multiple atoms.

        * excluding the detuning term, H_Delta.

        Q: Why do it with lists this way?
        A: For time-dependent interactions, it allows us to specify different
            time functions for each interaction in the list, e.g. pulse a
            probe laser on while a coupling laser is continous.

        This is the default interaction, when there's just a H_Omega to
        take care of. In other cases, where there of other terms to be added to
        the Hamiltonian, write a function in the OB object to override this.
        """

        return self.H_Omega_list

    def H_I_sum(self) -> qu.Qobj:
        """In the case of a time-independent interaction, we don't need the
        time functions, so just sum the items in H_I_list()

        This only works in the case of time-independent interaction,
        where no time_funcs have been specified."""

        if any(isinstance(h, list) for h in self.H_I_list()):
            raise ValueError(
                "H_I_sum() called on a time-dependent interaction. "
                "Use H_I_list() directly and pass it to the solver."
            )

        H_I = qu.Qobj(np.zeros([self.num_states, self.num_states]))

        for H_i in self.H_I_list():
            H_I += H_i

        return H_I

    def states_t(self) -> np.ndarray:
        """Returns:
            If c_ops: 3D array, first dim is time, the other two are the
        density matrix at each time slice.
            If no c_ops: 3D array, first dim is time, second is state vector,
                third is 1.
        """

        # No collapse operators, so states is a list of state vectors
        if len(self.c_ops) == 0:
            states_t = np.zeros(
                (len(self.result.times), self.num_states, 1), dtype=complex
            )

        # Collapse operators, so states is a list of density matrices
        else:
            states_t = np.zeros(
                (len(self.result.times), self.num_states, self.num_states),
                dtype=complex,
            )

        for t, state_t in enumerate(self.result.states):
            states_t[t] = state_t.full()

        return states_t

    def mesolve(
        self,
        tlist: np.ndarray,
        rho0: qu.Qobj | None = None,
        td: bool = False,
        e_ops: list | None = None,
        args: dict[str, Any] | None = None,
        options: dict[str, Any] | None = None,
        recalc: bool = True,
        savefile: str | None = None,
        show_pbar: bool = False,
    ) -> Any:

        if e_ops is None:
            e_ops = []
        if args is None:
            args = {}
        if options is None:
            options = {}

        savefile_exists = os.path.isfile(str(savefile) + ".qu")

        # Solve if 1) we ask for it to be recalculated or 2) it *must* be
        # calculated because no savefile exists.
        if recalc or not savefile_exists:
            # Is the Hamiltonian time-dependent?
            if td:  # If so H is a list of [H_i, t_func_i] pairs.
                H = [self.H_0, self.H_Delta]
                H.extend(self.H_I_list())
            else:  # If not it's a single QObj
                H = self.H_0 + self.H_Delta + self.H_I_sum()

            if show_pbar:
                options = dict(options)  # don't mutate caller's dict
                options["progress_bar"] = "text"

            self.result = qu.mesolve(
                H, rho0, tlist, self.c_ops, e_ops=e_ops, args=args, options=options
            )

            self.rho = self.result.states[-1]  #  Set rho to the final state.

            # Only save the file if we have a place to save it.
            if savefile:
                os.makedirs(os.path.dirname(savefile) or ".", exist_ok=True)
                print("Saving OBBase to {0}.qu".format(savefile))
                qu.qsave(self.result, savefile)

        # Otherwise load the steady state rho_v_delta from file
        else:
            print("Loading from {0}.qu".format(savefile))

            self.result = qu.qload(savefile)
            self.rho = self.result.states[-1]

        return self.result

    # def essolve(self, tlist, rho0=None, recalc=True, savefile=None):
    #     """
    #     Evolution of the density matrix by
    #     expressing the ODE as an exponential series.

    #     Args:
    #         tlist: The list of times for which to find the density matrix.
    #         rho0: Define an initial density matrix. Default is ground state.
    #         recalc: Rerun the calculation if a saved file exists?
    #         savefile: (string) Save the data to savefile.qu

    #     Returns:
    #         result: qutip result object containing the solved data.

    #     Notes:
    #         QuTiP essolve method doesn't return the states properly so I use
    #         the underlying ode2es method.

    #         Unlike the mesolve method, the tlist here doesn't have any need
    #         for high resolution to solve. So better when the density matrix
    #         is only needed at a few points.

    #     """

    #     savefile_exists = os.path.isfile(str(savefile) + '.qu')

    #     # Solve if 1) we ask for it to be recalculated or 2) it *must* be
    #     # calculated because no savefile exists.
    #     if recalc or not savefile_exists:

    #         H = self.H_0 + self.H_Delta + self.H_I_sum()
    #         L = qu.liouvillian(H, self.c_ops)

    #         es = qu.ode2es(L, rho0)
    #         states = es.value(tlist)

    #         self.result = qu.solver.Result()
    #         self.result.states = states
    #         self.result.solver = "essolve"
    #         self.result.times = tlist

    #         self.rho = self.result.states[-1]  #  Set rho to the final state.

    #         # Only save the file if we have a place to save it.
    #         if savefile:
    #             qu.qsave(self.result, savefile)

    #     # Otherwise load the steady state rho_v_delta from file
    #     else:
    #         self.result = qu.qload(savefile)
    #         self.rho = self.result.states[-1]

    #     return self.result

    # def steadystate(self, **kwargs):
    #     """ Calculates the steady state of the system in the case of a
    #     time-independent interaction, i.e. if we let time go to infinity.

    #     Returns:
    #         rho: The density matrix in the steady state.
    #     """

    #     H = self.H_0 + self.H_Delta + self.H_I_sum()
    #     self.rho = qu.steadystate(H, self.c_ops, **kwargs)

    #     return self.rho
