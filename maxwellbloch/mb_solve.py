# -*- coding: utf-8 -*-

# import os
import sys

import numpy as np
from scipy import interpolate

import qutip as qu

from maxwellbloch import ob_solve, t_funcs

# from cppy import deepcopy

class MBSolve(ob_solve.OBSolve):

    def __init__(self, ob_atom={}, t_min=0.0, t_max=1.0, t_steps=100,
                 method='mesolve', opts={}, savefile=None, z_min=0.0,
                 z_max=1.0, z_steps=10, z_steps_inner=2,
                 num_density_z_func=None, num_density_z_args={},
                 interaction_strengths=[], velocity_classes={}):

        super().__init__(ob_atom, t_min, t_max, t_steps,
                         method, opts, savefile)

        self.build_zlist(z_min, z_max, z_steps, z_steps_inner)

        self.build_number_density(num_density_z_func, num_density_z_args,
                                  interaction_strengths)

        self.build_velocity_classes(velocity_classes)

        self.init_Omegas_zt()
        self.init_states_zt()

    def __repr__(self):
        """ TODO: needs interaction strengths """

        return ("MBSolve(ob_atom={0}, " +
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
                "savefile={13})").format(self.ob_atom,
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

        self.velocity_classes = {}

        # What if only some elements are defined?
        if velocity_classes:
            self.velocity_classes = velocity_classes

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
            self.num_density_z_func = getattr(t_funcs, num_density_z_func)
        else:
            self.num_density_z_func = t_funcs.square_1

        if num_density_z_args:
            self.num_density_z_args = num_density_z_args
        else:
            self.num_density_z_args = {'on_1': 0.0, 'off_1': 1.0, 'ampl_1':1.0}

    def init_Omegas_zt(self):

        self.Omegas_zt = np.zeros((len(self.ob_atom.fields), len(self.zlist),
                                   len(self.tlist)), dtype=np.complex)

        ### Set the initial Omegas to the field time func values
        for f_i, f in enumerate(self.ob_atom.fields):
            self.Omegas_zt[f_i][0] = 2.0*np.pi*f.rabi_freq* \
                f.rabi_freq_t_func(self.tlist, f.rabi_freq_t_args)

        return self.Omegas_zt

    def init_states_zt(self):

        # TODO: Change states_zt to state_zt. Omegas refers to the fact that
        # there are multiple fields.

        self.states_zt = np.zeros((len(self.zlist), len(self.tlist),
                                   self.ob_atom.num_states,
                                   self.ob_atom.num_states),
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
        rabi_freq_ones = np.ones(len(self.ob_atom.fields))

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

                self.ob_atom.set_H_Omega(rabi_freq_ones,
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
        rabi_freq_ones = np.ones(len(self.ob_atom.fields))

        ### Set initial states at z=0
        self.states_zt[0, :] = \
            self.solve_and_average_over_thermal_detunings()

        # We won't need these until the first AB step
        sum_coh_prev = self.ob_atom.get_fields_sum_coherence()
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

        self.ob_atom.set_H_Omega(rabi_freq_ones,
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

                sum_coh_this = self.ob_atom.get_fields_sum_coherence()

                Omegas_z_next = self.z_step_fields_ab(z_prev, z_this,
                                                      z_next, sum_coh_prev,
                                                      sum_coh_this,
                                                      Omegas_z_this)

                Omegas_z_next_args = \
                    self.get_Omegas_intp_t_args(Omegas_z_next)

                self.ob_atom.set_H_Omega(rabi_freq_ones,
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
        """ For the current state of the atom, given field Omegas_z_this.

        """

        h = z_next - z_this

        N = self.num_density_z_func(z_next, self.num_density_z_args)

        Omegas_z_next = np.zeros((len(self.ob_atom.fields), len(self.tlist)),
                                 dtype=np.complex)


        sum_coh = self.ob_atom.get_fields_sum_coherence()

        for f_i, f in enumerate(self.ob_atom.fields):

            dOmega_f_dz = 1.0j*N*self.g[f_i]*sum_coh[f_i]

            Omegas_z_next[f_i, :] = Omegas_z_this[f_i, :] + h*dOmega_f_dz

        return Omegas_z_next

    def z_step_fields_ab(self, z_prev, z_this, z_next, sum_coh_prev,
        sum_coh_this, Omegas_z_this):

        # Assumes same step size
        h = z_next - z_this

        N = self.num_density_z_func(z_next, self.num_density_z_args)

        Omegas_z_next = np.zeros((len(self.ob_atom.fields), len(self.tlist)),
                                 dtype=np.complex)

        for f_i, f in enumerate(self.ob_atom.fields):

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

            Returns: A list of strings ['intp_1', 'intp_2', …]
        """

        rabi_freq_t_funcs = []
        for f_i, f in enumerate(self.ob_atom.fields, start=1):
            rabi_freq_t_funcs.append('intp_{0}'.format(f_i))
        return rabi_freq_t_funcs

    def get_Omegas_intp_t_args(self, Omegas_z):
        """ Return the values of Omegas at a given point as a list of
            args for interpolation

            e.g. [{'tlist_1': [], 'ylist_1': []},
                  {'tlist_2': [], 'ylist_2': []}]

            Note:
                The factor of 1/2pi is needed as we pass Rabi freq functions
                in without the factor of 2pi.

        """

        fields_args = [{}] * len(self.ob_atom.fields)

        for f_i, f in enumerate(Omegas_z):
            fields_args[f_i] = {'tlist_{0}'.format(f_i + 1): self.tlist,
                                'ylist_{0}'.format(f_i + 1):
                                Omegas_z[f_i] / (2.0 * np.pi)}

        return fields_args

    def solve_and_average_over_thermal_detunings(self):

        states_t_Delta = np.zeros((len(self.thermal_delta_list),
                                   len(self.tlist), self.ob_atom.num_states,
                                   self.ob_atom.num_states), dtype=np.complex)

        for Delta_i, Delta in enumerate(self.thermal_delta_list):

            # Shift each detuning by Delta
            self.ob_atom.shift_H_Delta([Delta] * len(self.ob_atom.fields))

            # We don't want the obsolve to save.
            self.solve(opts=qu.Options(max_step=self.t_step()), save=False)

            states_t_Delta[Delta_i] = self.states_t()

        thermal_states_t = np.average(states_t_Delta, axis=0,
                                      weights=self.thermal_weights)

        return thermal_states_t

    def save_results(self):

        # Only save the file if we have a place to save it.
        if self.savefile:
            print('Saving MBSolve to', self.savefile, '.qu')
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

    def tlist_fixed_frame(self, speed_of_light):
        """ Return the time points shifted to the fixed (lab) frame of
            reference given a speed-of-light.

            Args:
                speed_of_light: The speed of light in the system.

            Returns:
                Array of time values in the fixed frame of reference.
        """

        t_scale = 1.0 + self.z_max/(speed_of_light*self.t_max)
        return self.tlist*t_scale

    def Omegas_fixed_frame(self, field_idx, speed_of_light,
                           interp_kind='linear'):

        """ Return the solved field results shifted to the fixed (lab) frame of
            reference given a speed-of-light by interpolation.

            Args:
                speed_of_light: The speed of light in the system.
                interp_kind: The kind of spline interpolation to use ('linear',
                    'cubic' or 'quintic')

            Returns:
                Array[num_fields, num_space_points, num_time_points] of field
                values in the fixed frame of reference.
        """

        Omegas_intp = interpolate.interp2d(self.tlist, self.zlist,
                                           self.Omegas_zt[field_idx].real,
                                           bounds_error=False, fill_value=0.,
                                           kind=interp_kind)

        Omegas_fixed = np.zeros(self.Omegas_zt[field_idx].shape, dtype=np.float)

        for j, t_j in enumerate(self.tlist_fixed_frame(speed_of_light)):
            for i, z_i in enumerate(self.zlist):
                Omegas_fixed[i, j] = Omegas_intp(t_j - z_i/speed_of_light, z_i)

        return Omegas_fixed

    def freq_list(self):
        """ Fourier transform of the tlist into the frequency domain for
            spectral analysis.

            Returns:
                Array[num_time_points] of frequency values.

        """

        t_step = self.tlist[1] - self.tlist[0]
        freq_list = np.fft.fftfreq(len(self.tlist), t_step)  # FFT Freq
        return np.fft.fftshift(freq_list)

    def Omega_freq(self, field_idx):
        """ Fourier transform of the field result of field index.

        Returns:
            Array[num_z_steps, num_t_steps] Field result in frequency domain.

        """

        Omega_zt = self.Omegas_zt[field_idx]

        Omega_fft = np.zeros(Omega_zt.shape, dtype=np.complex)

        # TODO: I should be able to do this without the loop by specifying axis?
        for i, Omega_z_i in enumerate(Omega_zt):

            Omega_fft[i] = np.fft.fft(Omega_zt[i])
            Omega_fft[i] = np.fft.fftshift(Omega_fft[i])

        return Omega_fft

    def spectral_absorption(self, field_idx, z_idx):
        """ Field absorption in the frequency domain.

        Args:
            field_idx: Field to return spectrum of.
            z_idx: z step at which to return absorption.

        Returns:
            Array[num_freq_points] of absorption values.

        Note:
            In the linear regime this is the imaginary part of the linear
            susceptibility (with a factor k/2).
            See TP Ogden thesis Eqn (2.58)
        """

        Omega_freq_abs = np.abs(self.Omega_freq(field_idx))

        return -np.log(Omega_freq_abs[z_idx]/Omega_freq_abs[0])

    def spectral_dispersion(self, field_idx, z_idx):
        """ Field dispersion in the frequency domain.

        Note:
            In the linear regime this is the real part of the linear
            susceptibility.

            See TP Ogden Thesis Eqn (2.59)
        """

        Omega_freq_angle = np.angle(self.Omega_freq(field_idx))

        return Omega_freq_angle[0] - Omega_freq_angle[z_idx]

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
