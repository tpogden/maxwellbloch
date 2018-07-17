# -*- coding: utf-8 -*-

import sys

import numpy as np

import qutip as qu

from maxwellbloch import ob_solve, t_funcs

class MBSolve(ob_solve.OBSolve):

    def __init__(self, atom={}, t_min=0.0, t_max=1.0, t_steps=100,
                 method='mesolve', opts={}, savefile=None, z_min=0.0,
                 z_max=1.0, z_steps=10, z_steps_inner=2,
                 num_density_z_func=None, num_density_z_args={},
                 interaction_strengths=[], velocity_classes={}):

        super().__init__(atom, t_min, t_max, t_steps,
                         method, opts, savefile)

        self.build_zlist(z_min, z_max, z_steps, z_steps_inner)

        self.build_number_density(num_density_z_func, num_density_z_args,
                                  interaction_strengths)

        self.build_velocity_classes(velocity_classes)

        self.init_Omegas_zt()
        self.init_states_zt()

    def __repr__(self):
        """ TODO: needs interaction strengths """

        return ("MBSolve(atom={0}, " +
                "t_min={1}, " +
                "t_max={2}, " +
                "t_steps={3}, " +
                "method={4}, " +
                "z_min={5}, " +
                "z_max={6}, " +
                "z_steps={7}, " +
                "z_steps_inner={8}, " +
                "num_density_z_func={9}, " +
                "num_density_z_args={10}, " +
                "velocity_classes={11}, " +
                "opts={12}, " +
                "savefile={13})").format(self.atom,
                                         self.t_min,
                                         self.t_max,
                                         self.t_steps,
                                         self.method,
                                         self.z_min,
                                         self.z_max,
                                         self.z_steps,
                                         self.z_steps_inner,
                                         self.num_density_z_func,
                                         self.num_density_z_args,
                                         self.velocity_classes,
                                         self.opts,
                                         self.savefile)

    def build_zlist(self, z_min, z_max, z_steps, z_steps_inner):

        self.z_min = z_min
        self.z_max = z_max
        self.z_steps = z_steps
        self.z_steps_inner = z_steps_inner

        self.zlist = np.linspace(z_min, z_max, z_steps+1)

        return self.zlist

    def z_step(self):

        return (self.z_max - self.z_min)/self.z_steps

    def z_step_inner(self):

        return self.z_step()/self.z_steps_inner

    def build_velocity_classes(self, velocity_classes={}):
        """ Build the velocity class structure from dict. """

        self.velocity_classes = {}

        if velocity_classes:

            self.velocity_classes = velocity_classes

            if 'thermal_width' not in velocity_classes:
                self.velocity_classes['thermal_width'] = 1.0

            if velocity_classes['thermal_width'] <= 0.0:
                raise ValueError('Thermal width must be > 0.0.')

            if 'thermal_delta_steps' not in velocity_classes:
                self.velocity_classes['thermal_delta_min'] = 0.0
                self.velocity_classes['thermal_delta_max'] = 0.0
                self.velocity_classes['thermal_delta_steps'] = 0

            # If no inner steps, set to default
            if 'thermal_delta_inner_steps' not in velocity_classes:
                self.velocity_classes['thermal_delta_inner_min'] = 0.0
                self.velocity_classes['thermal_delta_inner_max'] = 0.0
                self.velocity_classes['thermal_delta_inner_steps'] = 0

        else:
            self.velocity_classes['thermal_delta_min'] = 0.0
            self.velocity_classes['thermal_delta_max'] = 0.0
            self.velocity_classes['thermal_delta_steps'] = 0
            self.velocity_classes['thermal_delta_inner_min'] = 0.0
            self.velocity_classes['thermal_delta_inner_max'] = 0.0
            self.velocity_classes['thermal_delta_inner_steps'] = 0
            self.velocity_classes['thermal_width'] = 1.0

        Delta_min = 2*np.pi*self.velocity_classes['thermal_delta_min']
        Delta_max = 2*np.pi*self.velocity_classes['thermal_delta_max']
        Delta_steps = self.velocity_classes['thermal_delta_steps']
        Delta_range = np.linspace(Delta_min, Delta_max, Delta_steps+1)

        # Narrow set of Delta classes around 0
        Delta2_min = 2*np.pi*self.velocity_classes['thermal_delta_inner_min']
        Delta2_max = 2*np.pi*self.velocity_classes['thermal_delta_inner_max']
        Delta2_steps = self.velocity_classes['thermal_delta_inner_steps']
        Delta2_range = np.linspace(Delta2_min, Delta2_max, Delta2_steps+1)

        # Merge the two ranges
        self.thermal_delta_list = np.unique(np.sort(np.append(Delta_range,
                                                              Delta2_range)))

        self.thermal_weights = maxwell_boltzmann(self.thermal_delta_list,
                                        self.velocity_classes['thermal_width'])

        return self.thermal_delta_list, self.thermal_weights

    def build_number_density(self, num_density_z_func, num_density_z_args,
                             interaction_strengths):

        self.interaction_strengths = interaction_strengths
        self.g = np.zeros(len(interaction_strengths))
        for i, g in enumerate(interaction_strengths):
            self.g[i] = 2*np.pi*g

        if num_density_z_func:
            self.num_density_z_func = getattr(t_funcs, num_density_z_func)(0)
        else:
            self.num_density_z_func = t_funcs.square(0)

        if num_density_z_args:
            self.num_density_z_args = {}
            for key, value in num_density_z_args.items():
                self.num_density_z_args[key + '_0'] = value
        else:
            self.num_density_z_args = {'on_0': 0.0, 'off_0': 1.0, 'ampl_0':1.0}

    def init_Omegas_zt(self):

        self.Omegas_zt = np.zeros((len(self.atom.fields), len(self.zlist),
                                   len(self.tlist)), dtype=np.complex)

        ### Set the initial Omegas to the field time func values
        for f_i, f in enumerate(self.atom.fields):
            self.Omegas_zt[f_i][0] = 2.0*np.pi*f.rabi_freq* \
                f.rabi_freq_t_func(self.tlist, f.rabi_freq_t_args)

        return self.Omegas_zt

    def init_states_zt(self):

        # TODO: Change states_zt to state_zt. Omegas refers to the fact that
        # there are multiple fields.

        self.states_zt = np.zeros((len(self.zlist), len(self.tlist),
                                   self.atom.num_states,
                                   self.atom.num_states),
                                  dtype=np.complex)

        return self.states_zt

    def mbsolve(self, step='ab', rho0=None, recalc=True, pbar_chunk_size=10):

        self.init_Omegas_zt()
        self.init_states_zt()

        if recalc or not self.savefile_exists():

            if step == 'euler':
                self.mbsolve_euler(rho0, recalc, pbar_chunk_size)

            elif step == 'ab':
                self.mbsolve_ab(rho0, recalc, pbar_chunk_size)

            self.save_results()

        else:
            self.load_results()

        return self.Omegas_zt, self.states_zt

    def mbsolve_euler(self, rho0=None, recalc=True, pbar_chunk_size=0):

        len_z = len(self.zlist)-1 # only for printing progress

        # For the interpolation
        rabi_freq_ones = np.ones(len(self.atom.fields))

        ### Set initial states at z=0
        self.states_zt[0, :] = \
            self.solve_and_average_over_thermal_detunings()
            # NOTE: This means the first z gets ob solved twice.
            # Can I avoid?

    ### All Steps:

        pbar = qu.ui.TextProgressBar(iterations=self.z_steps,
                                     chunk_size=pbar_chunk_size)

        for j, z in enumerate(self.zlist[:-1]):

            pbar.update(j)

            # Set initial fields and state
            Omegas_z_this = self.Omegas_zt[:, j, :]

            z_this = z

        ### Inner z loop

            for jj in range(self.z_steps_inner):

                z_next = z_this + self.z_step_inner()

                self.solve_and_average_over_thermal_detunings()

                Omegas_z_next = self.z_step_fields_euler(z_this, z_next,
                                                         Omegas_z_this)

                Omegas_z_next_args = \
                    self.get_Omegas_intp_t_args(Omegas_z_next)

                self.atom.set_H_Omega(rabi_freq_ones,
                                         self.get_Omegas_intp_t_funcs(),
                                         Omegas_z_next_args)

                # Set up for next inner step
                Omegas_z_this = Omegas_z_next
                z_this = z_next

            # Once we've been through the inner loops we can set the field
            # and states at the next outer space step and continue.
            self.states_zt[j+1, :] = self.states_t()
            self.Omegas_zt[:, j+1, :] = Omegas_z_next

        pbar.finished()

        return self.Omegas_zt, self.states_zt

    def mbsolve_ab(self, rho0=None, recalc=True, pbar_chunk_size=0):

        # For use in the loop
        rabi_freq_ones = np.ones(len(self.atom.fields))

        ### Set initial states at z=0
        self.states_zt[0, :] = \
            self.solve_and_average_over_thermal_detunings()

        # We won't need these until the first AB step
        sum_coh_prev = self.atom.get_fields_sum_coherence()
        z_prev = self.z_min

    ### First z step, Euler

        j = 0
        z = self.z_min

        # Set initial fields and state
        Omegas_z_this = self.Omegas_zt[:, j, :]

        z_this = z
        z_next = z_this + self.z_step_inner()

        Omegas_z_next = self.z_step_fields_euler(z_this, z_next,
                                                 Omegas_z_this)

        Omegas_z_next_args = \
            self.get_Omegas_intp_t_args(Omegas_z_next)

        self.atom.set_H_Omega(rabi_freq_ones,
                                 self.get_Omegas_intp_t_funcs(),
                                 Omegas_z_next_args)

        self.states_zt[j + 1, :] = self.states_t()
        self.Omegas_zt[:, j + 1, :] = Omegas_z_next

    ### Remaining steps, Adams-Bashforth

        pbar = qu.ui.TextProgressBar(iterations=self.z_steps,
                                     chunk_size=pbar_chunk_size)

        for j, z in enumerate(self.zlist[1:-1], start=1):

            pbar.update(j)

            # Set initial fields and state
            Omegas_z_this = self.Omegas_zt[:, j, :]

            z_this = z

        ### Inner z loop

            for jj in range(self.z_steps_inner):

                z_next = z_this + self.z_step_inner()

                self.solve_and_average_over_thermal_detunings()

                sum_coh_this = self.atom.get_fields_sum_coherence()

                Omegas_z_next = self.z_step_fields_ab(z_prev, z_this,
                                                      z_next, sum_coh_prev,
                                                      sum_coh_this,
                                                      Omegas_z_this)

                Omegas_z_next_args = \
                    self.get_Omegas_intp_t_args(Omegas_z_next)

                self.atom.set_H_Omega(rabi_freq_ones,
                                         self.get_Omegas_intp_t_funcs(),
                                         Omegas_z_next_args)

                # Set up for next inner step
                Omegas_z_this = Omegas_z_next
                z_this = z_next
                z_prev = z_this
                sum_coh_prev = sum_coh_this

            # Once we've been through the inner loops we can set the field
            # and states at the next outer space step and continue.
            self.states_zt[j+1, :] = self.states_t()
            self.Omegas_zt[:, j+1, :] = Omegas_z_next

        pbar.finished()

        return self.Omegas_zt, self.states_zt

    def z_step_fields_euler(self, z_this, z_next, Omegas_z_this):
        """ For the current state of the atom, given fields Omegas_z_this,
            make an Euler step to determine the field at the next space step.

            Args:
                z_this: The current space point
                z_next: the next space point
                Omegas_z_this: np.complex[num_t_steps]
                    The field rabi frequencies at z_step.
        """

        h = z_next - z_this
        N = self.num_density_z_func(z_next, self.num_density_z_args)
        Omegas_z_next = np.zeros((len(self.atom.fields), len(self.tlist)),
                                 dtype=np.complex)
        sum_coh = self.atom.get_fields_sum_coherence()

        for f_i, f in enumerate(self.atom.fields):
            dOmega_f_dz = 1.0j*N*self.g[f_i]*sum_coh[f_i]
            Omegas_z_next[f_i, :] = Omegas_z_this[f_i, :] + h*dOmega_f_dz

        return Omegas_z_next

    def z_step_fields_ab(self, z_prev, z_this, z_next, sum_coh_prev,
        sum_coh_this, Omegas_z_this):
        """ For the current state of the atom, given fields Omegas_z_this,
            make an Euler step to determine the field at the next space step.

            Args:
                z_prev: The previous space point
                z_this: The current space point
                z_next: the next space point
                Omegas_z_this: np.complex[num_t_steps]
                    The field rabi frequencies at z_step.
        """

        # Assumes same step size
        h = z_next - z_this

        N = self.num_density_z_func(z_next, self.num_density_z_args)

        Omegas_z_next = np.zeros((len(self.atom.fields), len(self.tlist)),
                                 dtype=np.complex)

        for f_i, f in enumerate(self.atom.fields):

            sum_coh_this_f = sum_coh_this[f_i]
            sum_coh_prev_f = sum_coh_prev[f_i]

            dOmega_f_dz_this = 1.0j * N * self.g[f_i] * sum_coh_this_f
            dOmega_f_dz_prev = 1.0j * N * self.g[f_i] * sum_coh_prev_f

            Omegas_z_next[f_i, :] = (Omegas_z_this[f_i, :] +
                                     1.5 * h * dOmega_f_dz_this -
                                     0.5 * h * dOmega_f_dz_prev)

        return Omegas_z_next

    def get_Omegas_intp_t_funcs(self):
        """ Gets a list of strings representing the interpolation t_funcs for
            use the MB solver, which needs a function representing the field
            at a z step to perform the next master equation solver.

            Returns: A list of strings ['intp', 'intp', …]
        """

        return ['intp' for f in self.atom.fields]

    def get_Omegas_intp_t_args(self, Omegas_z):
        """ Return the values of Omegas at a given point as a list of
            args for interpolation

            e.g. [{'tlist': [], 'ylist': []},
                  {'tlist': [], 'ylist': []}]

            Note:
                The factor of 1/2pi is needed as we pass Rabi freq functions
                in without the factor of 2pi.

        """

        fields_args = [{}] * len(self.atom.fields)
        for f_i, f in enumerate(Omegas_z):
            fields_args[f_i] = {'tlist': self.tlist,
                                'ylist': Omegas_z[f_i] / (2.0*np.pi)}
        return fields_args

    def solve_and_average_over_thermal_detunings(self):

        states_t_Delta = np.zeros((len(self.thermal_delta_list),
                                   len(self.tlist), self.atom.num_states,
                                   self.atom.num_states), dtype=np.complex)

        for Delta_i, Delta in enumerate(self.thermal_delta_list):

            # Shift each detuning by Delta
            self.atom.shift_H_Delta([Delta] * len(self.atom.fields))

            # We don't want the obsolve to save.
            self.solve(opts=qu.Options(max_step=self.t_step()), save=False)

            states_t_Delta[Delta_i] = self.states_t()

        thermal_states_t = np.average(states_t_Delta, axis=0,
                                      weights=self.thermal_weights)

        return thermal_states_t

    def save_results(self):

        # Only save the file if we have a place to save it.
        if self.savefile:
            print('Saving MBSolve to', self.savefile+'.qu')
            qu.qsave((self.Omegas_zt, self.states_zt), self.savefile)

    def load_results(self):

        self.Omegas_zt, self.states_zt = qu.qload(self.savefile)

    def fields_area(self):
        """ Get the integrated pulse area of each field.

        Returns:
            np.array [num_fields, num_z_steps]: Integrated area of each field
            over time
        """

        return np.trapz(np.real(self.Omegas_zt), self.tlist, axis=2)

### Helper Functions

def maxwell_boltzmann(v, fwhm):
    """ Maxwell Boltzmann probability distribution function. """

    # TODO: Allow offset, v_0.

    return 1./(fwhm*np.sqrt(np.pi))*np.exp(-(v/fwhm)**2)

def parse_args():

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file',
                        help='path of a JSON file containing an mb_solve \
                              problem definition', required=True)
    parser.add_argument('-s', '--step',
                        help="Finite difference step to use. 'ab' or 'euler'",
                        required=False, choices=['ab', 'euler'], default='ab')
    parser.add_argument('-r', '--recalc',
                        help='',
                        default=False, action="store_true", required=False)
    parser.add_argument('-p', '--pbarchunksize',
                        help="How often should the progress bar be updated? \
                              Default: 10 for update every 10%.",
                        required=False, type=int, default=10)

    args = vars(parser.parse_args())

    print('Loading problem definition from file {0}'.format(args['file']))
    mb_solve_obj = MBSolve().from_json(args['file'])

    mb_solve_obj.mbsolve(step=args['step'], rho0=None, recalc=args['recalc'],
                         pbar_chunk_size=args['pbarchunksize'])

def main():

    parse_args()

    return 0

if __name__ == '__main__':
    STATUS = main()
    sys.exit(STATUS)
