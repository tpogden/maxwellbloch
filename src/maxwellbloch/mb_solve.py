# -*- coding: utf-8 -*-


import numpy as np
import qutip as qu

from maxwellbloch import ob_solve, t_funcs


def _print_progress(j: int, total: int, chunk_size: int) -> None:
    """Print a progress update every chunk_size percent of total steps."""
    if chunk_size > 0 and total > 0:
        interval = max(1, total * chunk_size // 100)
        if j % interval == 0:
            pct = 100 * j // total
            print(f"  z-step {j}/{total} ({pct}%)", flush=True)


class MBSolve(ob_solve.OBSolve):
    def __init__(
        self,
        atom: dict | None = None,
        t_min: float = 0.0,
        t_max: float = 1.0,
        t_steps: int = 100,
        method: str = "mesolve",
        opts: dict | None = None,
        savefile: str | None = None,
        z_min: float = 0.0,
        z_max: float = 1.0,
        z_steps: int = 10,
        z_steps_inner: int = 2,
        num_density_z_func: str | None = None,
        num_density_z_args: dict | None = None,
        interaction_strengths: list[float] | None = None,
        velocity_classes: dict | None = None,
    ) -> None:

        if interaction_strengths is None:
            interaction_strengths = []

        super().__init__(
            atom=atom,
            t_min=t_min,
            t_max=t_max,
            t_steps=t_steps,
            method=method,
            opts=opts,
            savefile=savefile,
        )

        self.build_zlist(
            z_min=z_min, z_max=z_max, z_steps=z_steps, z_steps_inner=z_steps_inner
        )

        self._build_number_density(
            interaction_strengths=interaction_strengths,
            num_density_z_func=num_density_z_func,
            num_density_z_args=num_density_z_args,
        )

        self.build_velocity_classes(velocity_classes=velocity_classes)

        self.init_Omegas_zt()
        self.init_states_zt()

    def __repr__(self):
        return (
            f"MBSolve(atom={self.atom}, "
            + f"t_min={self.t_min}, "
            + f"t_max={self.t_max}, "
            + f"t_steps={self.t_steps}, "
            + f"method={self.method}, "
            + f"z_min={self.z_min}, "
            + f"z_max={self.z_max}, "
            + f"z_steps={self.z_steps}, "
            + f"z_steps_inner={self.z_steps_inner}, "
            + f"num_density_z_func={self.num_density_z_func}, "
            + f"num_density_z_args={self.num_density_z_args}, "
            + f"interaction_strengths={self.interaction_strengths}, "
            + f"velocity_classes={self.velocity_classes}, "
            + f"opts={self.opts}, "
            + f"savefile={self.savefile})"
        )

    def build_zlist(
        self, z_min: float, z_max: float, z_steps: int, z_steps_inner: int
    ) -> np.ndarray:
        """Builds the space grid.

        Args:
            z_min: The front of the medium, in your chosen length unit.
            z_max: The back of the medium.
            z_steps: The number of even-spaced steps on which to solve and
                record the solution.
            z_steps_inner: Between each z_step, make this many inner steps of
                the finite-difference solver (for numerical stability).

        Notes:
            - If the problem requires a lot of space steps for stability, but
                you don't need to record the solution at such a high-level of
                resolution, increase z_steps_inner.
        """
        self.z_min = z_min
        self.z_max = z_max
        self.z_steps = z_steps
        self.z_steps_inner = z_steps_inner
        self.zlist = np.linspace(z_min, z_max, z_steps + 1)
        return self.zlist

    def z_step(self) -> float:
        """Returns the distance from one space point to the next."""
        return (self.z_max - self.z_min) / self.z_steps

    def z_step_inner(self) -> float:
        """Returns the distance from one inner space point to the next."""
        return self.z_step() / self.z_steps_inner

    def _normalise_velocity_classes(self, velocity_classes: dict | None) -> dict:
        """Return a copy of velocity_classes with all keys filled with defaults.

        Args:
            velocity_classes: user-supplied dict, or None / empty dict for the
                no-Doppler-broadening case.

        Returns:
            A new dict with every expected key present.

        Raises:
            ValueError: if thermal_width is not positive.
        """
        defaults: dict = {
            "thermal_width": 1.0,
            "thermal_delta_min": 0.0,
            "thermal_delta_max": 0.0,
            "thermal_delta_steps": 0,
            "thermal_delta_inner_min": 0.0,
            "thermal_delta_inner_max": 0.0,
            "thermal_delta_inner_steps": 0,
        }
        vc = {**defaults, **(velocity_classes or {})}
        if vc["thermal_width"] <= 0.0:
            raise ValueError("Thermal width must be > 0.0.")
        return vc

    def _build_thermal_delta_list(self, vc: dict) -> np.ndarray:
        """Merge the outer and inner detuning ranges into one sorted unique array.

        Args:
            vc: normalised velocity-class dict (from _normalise_velocity_classes).

        Returns:
            1-D array of thermal detuning values in angular-frequency units.
        """
        outer = np.linspace(
            2 * np.pi * vc["thermal_delta_min"],
            2 * np.pi * vc["thermal_delta_max"],
            vc["thermal_delta_steps"] + 1,
        )
        inner = np.linspace(
            2 * np.pi * vc["thermal_delta_inner_min"],
            2 * np.pi * vc["thermal_delta_inner_max"],
            vc["thermal_delta_inner_steps"] + 1,
        )
        return np.unique(np.concatenate([outer, inner]))

    def build_velocity_classes(
        self, velocity_classes: dict | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        """Build the velocity-class detuning grid and Boltzmann weights.

        Args:
            velocity_classes: dict of velocity-class parameters, or None for
                the no-Doppler-broadening case (single detuning class at 0).

        Returns:
            (thermal_delta_list, thermal_weights)
        """
        self.velocity_classes = self._normalise_velocity_classes(velocity_classes)
        self.thermal_delta_list = self._build_thermal_delta_list(self.velocity_classes)
        self.thermal_weights = maxwell_boltzmann(
            self.thermal_delta_list,
            2 * np.pi * self.velocity_classes["thermal_width"],
        )
        return self.thermal_delta_list, self.thermal_weights

    def _build_number_density(
        self,
        interaction_strengths: list[float],
        num_density_z_func: str | None = None,
        num_density_z_args: dict | None = None,
    ) -> None:
        """Builds the number density function and interaction strengths.

        Args:
            interaction_strengths: A list of interaction strengths, `g`. These
                map 1-to-1 to the fields, and so must have the same number of
                values as there are fields.
            num_density_z_func: An optioanl function provided just like the time
                funcs (from t_funcs.py).
            num_density_z_args: A dict providing the args to be passed to
                num_density_z_func.

        Notes:
            - A factor of 2pi is applied to the interaction_strengths.
            - By default the num_density function will be square, with density
                1.0, starting at z=0.0 and ending at z=1.0.
        """
        self.interaction_strengths = interaction_strengths
        self.g = np.zeros(len(interaction_strengths))
        for i, g in enumerate(interaction_strengths):
            self.g[i] = 2 * np.pi * g
        # Set the num_density function
        if num_density_z_func:
            self.num_density_z_func = getattr(t_funcs, num_density_z_func)(0)
        else:
            self.num_density_z_func = t_funcs.square(0)
        # Set the num_density args
        if num_density_z_args:
            self.num_density_z_args = {}
            for key, value in num_density_z_args.items():
                self.num_density_z_args[key + "_0"] = value
        else:
            self.num_density_z_args = {"on_0": 0.0, "off_0": 1.0, "ampl_0": 1.0}

    def check(self) -> bool:
        """Validates the MBSolve object."""
        if len(self.interaction_strengths) != len(self.atom.fields):
            raise ValueError(
                "The number of interaction_strengths must match the number of fields."
            )
        return True

    def init_Omegas_zt(self) -> np.ndarray:
        """Inits the Rabi frequency array.

        Omegas_zt shape: (num_fields, z_steps+1, t_steps+1) — field-first.
        states_zt shape: (z_steps+1, t_steps+1, num_states, num_states) — z-first.
        These are inconsistent; correcting either is a breaking API change,
        deferred to a future major version.
        """
        self.Omegas_zt = np.zeros(
            (len(self.atom.fields), len(self.zlist), len(self.tlist)), dtype=complex
        )
        self._Omegas_z_buf = np.zeros(
            (len(self.atom.fields), len(self.tlist)), dtype=complex
        )
        for f_i, f in enumerate(self.atom.fields):
            self.Omegas_zt[f_i][0] = (
                2.0
                * np.pi
                * f.rabi_freq
                * f.rabi_freq_t_func(self.tlist, f.rabi_freq_t_args)
            )
        return self.Omegas_zt

    def init_states_zt(self) -> np.ndarray:
        """Inits the system density matrices."""
        self.states_zt = np.zeros(
            (
                len(self.zlist),
                len(self.tlist),
                self.atom.num_states,
                self.atom.num_states,
            ),
            dtype=complex,
        )
        return self.states_zt

    # TODO(#96) Should we be able to pass in opts here?
    def mbsolve(
        self,
        step: str = "ab",
        rho0: qu.Qobj | None = None,
        recalc: bool = True,
        pbar_chunk_size: int = 10,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Solves the Maxwell-Bloch equations for the system.

        Args:
            step: 'euler (for Euler method) or 'ab' (for Adams-Bashforth)
            rho0 (Qobj): the initial density matrix state
            recalc (bool): Recalculate the solution even if a savefile exists?
            pbar_chunk_size: Every how many % should we log progress to stdout.

        Returns:
            self.Omegas_zt: The solved field complex Rabi frequency at each
                point in space z and time t.
            self.states_zt: The solved density matrix at each point in space z
                and time t.
        """
        self.check()
        self.init_Omegas_zt()
        self.init_states_zt()
        # Should we recalculate or load a savefile?
        if recalc or not self.savefile_exists():
            if step == "euler":
                self.mbsolve_euler(
                    rho0=rho0, recalc=recalc, pbar_chunk_size=pbar_chunk_size
                )
            elif step == "ab":
                self.mbsolve_ab(
                    rho0=rho0, recalc=recalc, pbar_chunk_size=pbar_chunk_size
                )
            self.save_results()
        else:
            self.load_results()
        return self.Omegas_zt, self.states_zt

    def mbsolve_euler(
        self,
        rho0: qu.Qobj | None = None,
        recalc: bool = True,
        pbar_chunk_size: int = 0,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Solves the Maxwell-Bloch equations using a Euler step.

        Args:
            rho0 (Qobj): the initial density matrix state
            recalc (bool): Recalculate the solution even if a savefile exists?
            pbar_chunk_size: Every how many % should we log progress to stdout.

        Returns:
            self.Omegas_zt: The solved field complex Rabi frequency at each
                point in space z and time t.
            self.states_zt: The solved density matrix at each point in space z
                and time t.
        """
        self.states_zt[0, :] = self._solve_and_average_over_thermal_detunings()
        for j, z in enumerate(self.zlist[:-1]):
            _print_progress(j=j, total=self.z_steps, chunk_size=pbar_chunk_size)
            Omegas_z_this = self.Omegas_zt[:, j, :]
            z_cur = z
            for _ in range(self.z_steps_inner):
                z_next = z_cur + self.z_step_inner()
                h = z_next - z_cur
                N = self.num_density_z_func(z_next, self.num_density_z_args)
                thermal_states_t = self._solve_and_average_over_thermal_detunings()
                sum_coh_this = self.atom.get_fields_sum_coherence(
                    states_t=thermal_states_t
                )
                np.copyto(
                    self._Omegas_z_buf,
                    self._z_step_fields_euler(
                        h=h,
                        N=N,
                        Omegas_z_this=Omegas_z_this,
                        sum_coh_this=sum_coh_this,
                    ),
                )
                self.atom.H_Omega_list = self._build_intp_H_Omega_list(
                    self._Omegas_z_buf
                )
                Omegas_z_this = self._Omegas_z_buf
                z_cur = z_next
            self.states_zt[j + 1, :] = self.states_t()
            self.Omegas_zt[:, j + 1, :] = self._Omegas_z_buf
        return self.Omegas_zt, self.states_zt

    def mbsolve_ab(
        self,
        rho0: qu.Qobj | None = None,
        recalc: bool = True,
        pbar_chunk_size: int = 0,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Solves the Maxwell-Bloch equations using an Adams-Bashforth step.

        Args:
            rho0 (Qobj): the initial density matrix state
            recalc (bool): Recalculate the solution even if a savefile exists?
            pbar_chunk_size: Every how many % should we log progress to stdout.

        Returns:
            self.Omegas_zt: The solved field complex Rabi frequency at each
                point in space z and time t.
            self.states_zt: The solved density matrix at each point in space z
                and time t.
        """
        thermal_states_t = self._solve_and_average_over_thermal_detunings()
        self.states_zt[0, :] = thermal_states_t
        sum_coh_prev = self.atom.get_fields_sum_coherence(states_t=thermal_states_t)
        # First z step: Euler bootstrap
        j = 0
        z_cur = self.z_min
        z_next = z_cur + self.z_step_inner()
        h = z_next - z_cur
        N = self.num_density_z_func(z_next, self.num_density_z_args)
        np.copyto(
            self._Omegas_z_buf,
            self._z_step_fields_euler(
                h=h,
                N=N,
                Omegas_z_this=self.Omegas_zt[:, j, :],
                sum_coh_this=sum_coh_prev,
            ),
        )
        self.atom.H_Omega_list = self._build_intp_H_Omega_list(self._Omegas_z_buf)
        self.states_zt[j + 1, :] = self.states_t()
        self.Omegas_zt[:, j + 1, :] = self._Omegas_z_buf
        # Remaining steps: Adams-Bashforth
        for j, z in enumerate(self.zlist[1:-1], start=1):
            _print_progress(j=j, total=self.z_steps, chunk_size=pbar_chunk_size)
            Omegas_z_this = self.Omegas_zt[:, j, :]
            z_cur = z
            for _ in range(self.z_steps_inner):
                z_next = z_cur + self.z_step_inner()
                h = z_next - z_cur
                N = self.num_density_z_func(z_next, self.num_density_z_args)
                thermal_states_t = self._solve_and_average_over_thermal_detunings()
                sum_coh_this = self.atom.get_fields_sum_coherence(
                    states_t=thermal_states_t
                )
                np.copyto(
                    self._Omegas_z_buf,
                    self._z_step_fields_ab(
                        h=h,
                        N=N,
                        sum_coh_prev=sum_coh_prev,
                        sum_coh_this=sum_coh_this,
                        Omegas_z_this=Omegas_z_this,
                    ),
                )
                self.atom.H_Omega_list = self._build_intp_H_Omega_list(
                    self._Omegas_z_buf
                )
                Omegas_z_this = self._Omegas_z_buf
                z_cur = z_next
                sum_coh_prev = sum_coh_this
            self.states_zt[j + 1, :] = self.states_t()
            self.Omegas_zt[:, j + 1, :] = self._Omegas_z_buf
        return self.Omegas_zt, self.states_zt

    def _z_step_fields_euler(
        self,
        h: float,
        N: float,
        Omegas_z_this: np.ndarray,
        sum_coh_this: np.ndarray,
    ) -> np.ndarray:
        """Euler step: advance all field Rabi frequencies by one spatial step.

        Args:
            h: spatial step size (z_next - z_this)
            N: number density at z_next
            Omegas_z_this: shape (num_fields, t_steps+1) — fields at z_this
            sum_coh_this: shape (num_fields, t_steps+1) — coherences at z_this

        Returns:
            Omegas_z_next: shape (num_fields, t_steps+1)
        """
        # self.g shape (num_fields,) broadcasts with sum_coh_this (num_fields, t_steps+1)
        dOmega_dz = 1.0j * N * self.g[:, None] * sum_coh_this
        return Omegas_z_this + h * dOmega_dz

    def _z_step_fields_ab(
        self,
        h: float,
        N: float,
        sum_coh_prev: np.ndarray,
        sum_coh_this: np.ndarray,
        Omegas_z_this: np.ndarray,
    ) -> np.ndarray:
        """Adams-Bashforth step: advance all field Rabi frequencies by one spatial step.

        Args:
            h: spatial step size (z_next - z_this); assumes uniform step size
            N: number density at z_next
            sum_coh_prev: shape (num_fields, t_steps+1) — coherences at z_prev
            sum_coh_this: shape (num_fields, t_steps+1) — coherences at z_this
            Omegas_z_this: shape (num_fields, t_steps+1) — fields at z_this

        Returns:
            Omegas_z_next: shape (num_fields, t_steps+1)
        """
        dOmega_dz_this = 1.0j * N * self.g[:, None] * sum_coh_this
        dOmega_dz_prev = 1.0j * N * self.g[:, None] * sum_coh_prev
        return Omegas_z_this + 1.5 * h * dOmega_dz_this - 0.5 * h * dOmega_dz_prev

    def _build_intp_H_Omega_list(self, Omegas_z: np.ndarray) -> list:
        """Build H_Omega_list for the next z-step using QuTiP InterCoefficient.

        At each spatial step the field Rabi frequency profile is known as an
        array over time.  Rather than wiring this through the intp() closure
        (which recreates a scipy interp1d on every ODE function call), we
        create a QuTiP InterCoefficient once per z-step — a compiled Cython
        interpolator that evaluates cheaply at each ODE step.

        Args:
            Omegas_z: shape (num_fields, t_steps+1) — complex Rabi frequency
                profiles at the next z position, in angular-frequency units.

        Returns:
            List of [H_Omega (Qobj), InterCoefficient] pairs, one per field.
        """
        H_Omega_list = []
        for f_i, f in enumerate(self.atom.fields):
            # Build the sigma structure with rabi_freq=1 (same as _build_H_Omega)
            H_Omega = qu.Qobj(np.zeros([self.atom.num_states, self.atom.num_states]))
            for c_i, c in enumerate(f.coupled_levels):
                H_Omega += (
                    self.atom.sigma(a=c[0], b=c[1]) + self.atom.sigma(a=c[1], b=c[0])
                ) * f.factors[c_i]
            H_Omega = H_Omega * np.pi  # pi * 1.0 (rabi_freq normalised to 1)
            # Divide by 2π: H_Omega * coeff(t) = sigma * pi * Omega_z/(2π)
            #                                   = sigma * Omega_z/2  ✓
            coeff = qu.coefficient(Omegas_z[f_i].real / (2 * np.pi), tlist=self.tlist)
            H_Omega_list.append([H_Omega, coeff])
        return H_Omega_list

    def _solve_and_average_over_thermal_detunings(self) -> np.ndarray:
        """Solves the Lindblad equation for the OBAtom over a range of
            detuning shifts for velocity classes.

        Returns:
            A states_t object that is the Maxwell-Boltzmann weighted average
            over the velocity classes.
        """
        states_t_Delta = np.zeros(
            (
                len(self.thermal_delta_list),
                len(self.tlist),
                self.atom.num_states,
                self.atom.num_states,
            ),
            dtype=complex,
        )
        # The set detunings, without any thermal shifting
        fixed_detunings = self.atom.get_detunings()
        for Delta_i, Delta in enumerate(self.thermal_delta_list):
            # print('Delta_i: {0}, Delta: {1:.2f}'.format(Delta_i, Delta))
            # Shift each detuning by Delta.
            self.atom.set_H_Delta([fd + Delta for fd in fixed_detunings])
            # We don't want the obsolve to save.
            # TODO(#96) If we decide to pass down opts from mbsolve, put here.
            self.solve(opts=None, save=False)
            states_t_Delta[Delta_i] = self.states_t()
        # Restore fixed detunings
        self.atom.set_H_Delta(fixed_detunings)
        thermal_states_t = np.average(
            states_t_Delta, axis=0, weights=self.thermal_weights
        )
        return thermal_states_t

    def save_results(self) -> None:
        """Saves the solution to a QuTiP pickle file.

        Notes:
            - The path to which the results will be saved is taken from
                self.savefile.
        """
        # Only save the file if we have a place to save it.
        if self.savefile:
            print("Saving MBSolve to", self.savefile + ".qu")
            qu.qsave((self.Omegas_zt, self.states_zt), self.savefile)

    def load_results(self) -> None:
        """Loads the solution from a QuTiP pickle file.

        Notes:
            - The path from which the results will be loaded is taken from
                self.savefile.
        """
        self.Omegas_zt, self.states_zt = qu.qload(self.savefile)

    def fields_area(self) -> np.ndarray:
        """Gets the integrated pulse area of each field.

        Returns:
            np.array [num_fields, num_z_steps]: Integrated area of each field
            over time
        """
        return np.trapezoid(np.real(self.Omegas_zt), self.tlist, axis=2)

    def populations(self, levels: list[int]) -> np.ndarray:
        """Gets the sum of populations in a list of levels.

        Args:
            levels: a list of level indexes]

        Returns:
            np.array, shape (z_steps+1, t_steps+1), dtype=np.real
        """
        return np.abs(np.sum(self.states_zt[:, :, levels, levels], axis=2))

    def populations_field(self, field_idx: int, upper: bool = True) -> np.ndarray:
        """Gets the sum of populations for the upper (excited) level coupled by
            a field.

        Args:
            field_idx: index in the list of fields

        Returns:
            np.array, shape (z_steps+1, t_steps+1), dtype=np.real

        Note:
            - Casting upper to int so upper is 1, lower is 0.
        """
        f = self.atom.fields[field_idx]
        levels = f.upper_levels() if upper else f.lower_levels()
        return self.populations(levels)

    def coherences(self, coupled_levels: list[list[int]]) -> np.ndarray:
        """Gets the sum of coherences (off-diagonals) in a list of coupled
            level pairs.

        Args:
            coupled_levels: a list of pairs of level indexes

        Returns:
            np.array, shape (z_steps+1, t_steps+1), dtype=complex

        Note:
            Unlike ``OBAtom.get_fields_sum_coherence``, which operates on a
            single-z time series with per-field weighting factors, this method
            sums over the full (z, t) array with unit weights.
        """
        rows = [cl[0] for cl in coupled_levels]
        cols = [cl[1] for cl in coupled_levels]
        return self.states_zt[:, :, rows, cols].sum(axis=2)

    def coherences_field(self, field_idx: int) -> np.ndarray:
        """Get the sum of coherences (off-diagonals) for the levels coupled by
        a field.

        Args:
            field_idx: index in the list of fields
        Returns:
            np.array, shape (z_steps+1, t_steps+1), dtype=complex
        """
        return self.coherences(self.atom.fields[field_idx].coupled_levels)


### Helper Functions


def maxwell_boltzmann(v: np.ndarray, fwhm: float) -> np.ndarray:
    """Maxwell Boltzmann probability distribution function."""

    # TODO: Allow offset, v_0.
    # TODO: move this to utility.py

    return 1.0 / (fwhm * np.sqrt(np.pi)) * np.exp(-((v / fwhm) ** 2))
