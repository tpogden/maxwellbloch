# -*- coding: utf-8 -*-

import os
import sys

from maxwellbloch import sigma

import numpy as np
import qutip as qu


class OBBase(object):
    """ TODO: Desc here. Parent class.

    Attributes:
        num_states: the number of states of the system, in this case: 2.
        rho: the density matrix

        H_0: the Bare Hamiltonian (see the OBObj method make_H_0())
        c_ops: a list of collapse operators (here just 1 item, the decay rate)

        H_Delta: the detuning part of the interaction
        H_Omega: the Rabi frequency part of the interaction

        time_range: an array of time points over which to solve the system
        ob_data: states and expectation values as solved at each step in
            time_range (see QuTip's odedata object).
        """

    def __init__(self):
        """ Initialise the OB object. This will be overridden in every child
        class, just here to make it clear what attributes an OB object has. """

        self.num_states = 0
        self.rho = qu.Qobj()

        self.H_0 = qu.Qobj()
        self.c_ops = []

        # These will be set by functions in the child class.
        self.H_Delta = qu.Qobj()
        self.H_Omega_list = [qu.Qobj()]

        self.result = qu.solver.Result()

    def sigma(self, a, b):
        """ The transition or 'flip' operator between states |a> and |b>, given
        by |a><b|. Uses the sigma_n helper function below this class.

        Args:
            a: level transitioned from
            b: level transitioned to

        Returns:
            |a><b|
        """

        return sigma.sigma(self.num_states, a, b)

    def ground_state(self):
        """ Returns a state vector for the ground state. [1, 0, 0, …] """

        return qu.basis(self.num_states, 0)

    def set_H_0(self, energies=[]):
        """ Takes a list of energies and makes a Bare Hamiltonian with the
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

    def H_I_list(self):
        """ The interaction terms* are specified as a list. This could be a
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

    def H_I_sum(self):
        """ In the case of a time-independent interaction, we don't need the
        time functions, so just sum the items in H_I_list()

        This only works in the case of time-independent interaction,
        where no time_funcs have been specified. """

        ## TODO: raise an error here if a time-dependent interaction has been
        ## specified, i.e. if H_I is a list of [H_I_term, time_func] pairs.

        H_I = qu.Qobj(np.zeros([self.num_states, self.num_states]))

        for H_i in self.H_I_list():
            H_I += H_i

        return H_I

    def states_t(self):
        """ Returns:
            If c_ops: 3D array, first dim is time, the other two are the
        density matrix at each time slice.
            If no c_ops: 3D array, first dim is time, second is state vector,
                third is 1.
        """

        # No collapse operators, so states is a list of state vectors
        if len(self.c_ops) == 0:
            states_t = np.zeros((len(self.result.times), self.num_states, 1),
                                dtype=np.complex)

        # Collapse operators, so states is a list of density matrices
        else:
            states_t = np.zeros((len(self.result.times),
                                 self.num_states, self.num_states),
                                dtype=np.complex)

        for t, state_t in enumerate(self.result.states):
            states_t[t] = state_t.full()

        return states_t

    def mesolve(self, tlist, rho0=None, td=False, e_ops=[], args={},
                opts=qu.Options(), recalc=True, savefile=None,
                show_pbar=False):

        if not rho0:
            rho0 = self.ground_state()

        savefile_exists = os.path.isfile(str(savefile) + '.qu')

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
                pbar = qu.ui.progressbar.TextProgressBar()
            else:
                pbar = qu.ui.progressbar.BaseProgressBar()

            self.result = qu.mesolve(H, rho0, tlist,
                                     self.c_ops, e_ops,
                                     args=args, options=opts,
                                     progress_bar=pbar)

            self.rho = self.result.states[-1]  #  Set rho to the final state.

            # Only save the file if we have a place to save it.
            if savefile:

                print('Saving OBBase to {0}.qu'.format(savefile))

                qu.qsave(self.result, savefile)

        # Otherwise load the steady state rho_v_delta from file
        else:

            print('Loading from {0}.qu'.format(savefile))

            self.result = qu.qload(savefile)
            self.rho = self.result.states[-1]

        return self.result

    def essolve(self, tlist, rho0=None, recalc=True, savefile=None):
        """
        Evolution of the density matrix by
        expressing the ODE as an exponential series.

        Args:
            tlist: The list of times for which to find the density matrix.
            rho0: Define an initial density matrix. Default is ground state.
            recalc: Rerun the calculation if a saved file exists?
            savefile: (string) Save the data to savefile.qu

        Returns:
            result: qutip result object containing the solved data.

        Notes:
            QuTiP essolve method doesn't return the states properly so I use
            the underlying ode2es method.

            Unlike the mesolve method, the tlist here doesn't have any need
            for high resolution to solve. So better when the density matrix
            is only needed at a few points.

        """

        if not rho0:
            rho0 = self.ground_state() * self.ground_state().dag()

        savefile_exists = os.path.isfile(str(savefile) + '.qu')

        # Solve if 1) we ask for it to be recalculated or 2) it *must* be
        # calculated because no savefile exists.
        if recalc or not savefile_exists:

            H = self.H_0 + self.H_Delta + self.H_I_sum()
            L = qu.liouvillian(H, self.c_ops)

            es = qu.ode2es(L, rho0)
            states = es.value(tlist)

            self.result = qu.solver.Result()
            self.result.states = states
            self.result.solver = "essolve"
            self.result.times = tlist

            self.rho = self.result.states[-1]  #  Set rho to the final state.

            # Only save the file if we have a place to save it.
            if savefile:
                qu.qsave(self.result, savefile)

        # Otherwise load the steady state rho_v_delta from file
        else:
            self.result = qu.qload(savefile)
            self.rho = self.result.states[-1]

        return self.result

    def steadystate(self, **kwargs):
        """ Calculates the steady state of the system in the case of a
        time-independent interaction, i.e. if we let time go to infinity.

        Returns:
            rho: The density matrix in the steady state.
        """

        H = self.H_0 + self.H_Delta + self.H_I_sum()
        self.rho = qu.steadystate(H, self.c_ops, **kwargs)

        return self.rho

# Main


def main():

    print(OBBase())


if __name__ == '__main__':
    status = main()
    sys.exit(status)
