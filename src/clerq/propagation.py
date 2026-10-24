# -*- coding: utf-8 -*-


import datetime
import hashlib
import importlib.metadata
import numbers
import os
import warnings
from typing import Any

import numpy as np
import qutip as qu
import tqdm
from qutip.solver.integrator.integrator import IntegratorException
from scipy.interpolate import CubicSpline

from clerq import t_funcs
from clerq._deprecation import resolve_system_arg, warn_deprecated
from clerq.exceptions import (
    CounterPropagatingDepletionError,
    PropagationInstabilityError,
)
from clerq.lindblad import LindbladSolver
from clerq.utility import maxwell_boltzmann

# Weighting of each coupled pair's coherence in the field source term, stored
# in savefiles. "factor" is the signed factor (as in the Hamiltonian); results
# solved before it was introduced used factor**2.
_SOURCE_WEIGHTS = "factor"

# thermal_delta_spacing="auto" is this fraction of 1 / (t_max - t_min). A grid
# of spacing h echoes at t = 1 / h, so this puts the echo at 1 / 0.75 times the
# window length. Measured on a two-level problem, the error is at its floor
# (below 1e-7) by 0.75 / T and 29% at 1.5 / T.
THERMAL_ECHO_MARGIN = 0.75

# thermal_delta_spacing="auto" lets the echo into the window if its amplitude,
# exp(-decay * t_echo), is below this. Measured on a two-level problem, the
# result error is within about 20 times this.
THERMAL_ECHO_TOLERANCE = 1e-6

# Above this many velocity classes thermal_delta_spacing warns: each class is a
# master-equation solve at every z-step, and a tiny spacing is usually a units
# slip (angular for cyclic).
THERMAL_CLASSES_WARN = 1000


class PropagationSolver(LindbladSolver):
    """Solver for the propagation of fields through a medium of systems.

    Steps the fields through space with a finite-difference method, solving
    the Lindblad master equation for the :class:`~clerq.system.System` at
    each step (once per velocity class if the medium is Doppler broadened).

    Args:
        system: Dict (or empty dict for defaults) describing the quantum
            system; passed to :class:`~clerq.system.System`.
        t_min: Start time of the simulation.
        t_max: End time of the simulation.
        t_steps: Number of time steps.
        method: Solver method — ``'mesolve'`` (Lindblad master equation) is
            the only supported value.
        opts: Dict of QuTiP solver options.
        savefile: Path prefix for caching solved results (``'.qu'`` extension
            appended automatically). Pass ``None`` to disable caching.
        z_min: The front of the medium.
        z_max: The back of the medium.
        z_steps: Number of space steps at which the solution is recorded.
        z_steps_inner: Number of finite-difference steps made within each
            space step. Pass ``"auto"`` to start from 2 and let ``solve()``
            raise it until the estimated step error is below ``step_error_warn``
            (see :meth:`solve`).
        num_density_z_func: Name of a function in ``clerq.t_funcs`` giving
            the number density profile over z. Square by default.
        num_density_z_args: Dict of arguments for ``num_density_z_func``.
        interaction_strengths: Interaction strength of each field with the
            medium, one for each field.
        velocity_classes: Dict defining the velocity classes for Doppler
            broadening, or ``None`` for no broadening. ``thermal_width`` and
            the ``thermal_delta_*`` values are cyclic frequencies, in the
            same units as ``decays`` rates and ``rabi_freq``. Instead of
            ``thermal_delta_steps`` the classes can be given a maximum
            spacing ``thermal_delta_spacing`` (new in 0.15.0), or ``"auto"``
            for the largest spacing whose spurious Doppler echoes fall
            outside the time window or are damped by the decay below
            ``thermal_delta_echo_tolerance`` (see
            ``thermal_delta_auto_spacing()``). The density
            of the classes is Gaussian, ``exp(-(delta / thermal_width)**2)``
            (``thermal_weights``); each class is weighted by this times the
            detuning interval it represents (``thermal_widths``).
        atom: Deprecated name for ``system``.
        real_drive: Deprecated. Drive the Hamiltonian with ``Re(Omega)`` only,
            the behaviour before the complex-envelope drive. That is exact for
            resonant problems with real fields and a symmetric,
            equally weighted velocity-class grid, but wrong for a single
            detuned field (for example a two-level medium, optical depth 6,
            half-width detuning 2: transmission 1.1 instead of 0.30). The
            default drives the coupled pair with ``Omega`` and ``conj(Omega)``.
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
        z_min: float = 0.0,
        z_max: float = 1.0,
        z_steps: int = 10,
        z_steps_inner: int | str = 2,
        num_density_z_func: str | None = None,
        num_density_z_args: dict | None = None,
        interaction_strengths: list[float] | None = None,
        velocity_classes: dict | None = None,
        *,
        atom: dict | None = None,
        real_drive: bool = False,
    ) -> None:

        system = resolve_system_arg(system, atom)
        self.real_drive = bool(real_drive)
        if self.real_drive:
            warn_deprecated(
                "real_drive=True (driving with Re(Omega) only)",
                "the default complex drive",
                stacklevel=3,
            )

        if interaction_strengths is None:
            interaction_strengths = []

        super().__init__(
            system=system,
            t_min=t_min,
            t_max=t_max,
            t_steps=t_steps,
            method=method,
            opts=opts,
            savefile=savefile,
        )

        self._auto_z_steps_inner = z_steps_inner == "auto"
        if self._auto_z_steps_inner:
            z_steps_inner = 2

        self.build_zlist(
            z_min=z_min, z_max=z_max, z_steps=z_steps, z_steps_inner=z_steps_inner
        )

        self._build_number_density(
            interaction_strengths=interaction_strengths,
            num_density_z_func=num_density_z_func,
            num_density_z_args=num_density_z_args,
        )

        self.build_velocity_classes(velocity_classes=velocity_classes)

        self._stability = None
        self.step_error_estimate = 0.0
        self.t_step_error_estimate = 0.0
        self.window_edge_estimate = 0.0

        self.init_Omegas_zt()
        self.init_states_zt()

    def __repr__(self):
        return (
            f"PropagationSolver(system={self.system}, "
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
        # Not in the defaults, so that the problem hash (and so existing
        # savefiles) is unchanged when it is not used.
        spacing = vc.get("thermal_delta_spacing")
        if spacing is not None:
            if vc["thermal_delta_steps"] or vc["thermal_delta_inner_steps"]:
                raise ValueError(
                    "thermal_delta_spacing replaces thermal_delta_steps and "
                    "thermal_delta_inner_steps; set only one of them."
                )
            if vc["thermal_delta_max"] <= vc["thermal_delta_min"]:
                raise ValueError(
                    "thermal_delta_spacing needs thermal_delta_max > thermal_delta_min."
                )
            if spacing != "auto" and not (
                isinstance(spacing, numbers.Real)
                and not isinstance(spacing, bool)
                and spacing > 0.0
            ):
                raise ValueError("thermal_delta_spacing must be > 0.0 or 'auto'.")
            if spacing != "auto":
                vc["thermal_delta_spacing"] = float(spacing)  # JSON-safe for the hash
        tolerance = vc.get("thermal_delta_echo_tolerance")
        if tolerance is not None:
            if spacing != "auto":
                raise ValueError(
                    "thermal_delta_echo_tolerance only applies with "
                    "thermal_delta_spacing='auto'."
                )
            if not (
                isinstance(tolerance, numbers.Real)
                and not isinstance(tolerance, bool)
                and 0.0 < tolerance < 1.0
            ):
                raise ValueError(
                    "thermal_delta_echo_tolerance must be between 0 and 1."
                )
            vc["thermal_delta_echo_tolerance"] = float(tolerance)
        return vc

    def _level_doppler_slopes(self) -> np.ndarray:
        """Cyclic Doppler slope of each level's energy, in the rotating frame.

        Each field adds ``factor_doppler_shift * Delta`` to its detuning, and
        the detuning is added to the energy of the field's upper levels (see
        ``System._build_H_Delta``), so level ``i`` moves at ``-sum(sign * s)``
        over the fields that have ``i`` as an upper level.
        """
        slopes = np.zeros(self.system.num_states)
        for f in self.system.fields:
            sign = 1.0 if f.detuning_positive else -1.0
            for level in f.upper_levels():
                slopes[level] -= sign * f.factor_doppler_shift
        return slopes

    def _level_decay_rates(self) -> np.ndarray:
        """Angular rate at which each level's coherences are damped.

        Half the sum of ``diag(c^dagger c)`` over the collapse operators, so
        the coherence between levels ``i`` and ``j`` decays at the sum of the
        two entries. Pure dephasing counts the same way.
        """
        total = np.zeros(self.system.num_states)
        for c_op in self.system.c_ops:
            total += np.real(np.diag((c_op.dag() * c_op).full()))
        return 0.5 * total

    def thermal_delta_auto_spacing(self, tolerance: float | None = None) -> float:
        """Largest cyclic class spacing free of Doppler echoes.

        Sampling the detuning at spacing ``h`` makes the summed response repeat
        at ``t = 1 / h``. A coherence between levels ``i`` and ``j`` oscillates
        at ``b * Delta``, where ``b`` is the difference of the two levels'
        Doppler slopes (1 for a single co-propagating field), so its echo is at
        ``t = 1 / (b h)``. Per pair, the spacing is the larger of:

        - the one that puts the echo outside the time window:
          ``THERMAL_ECHO_MARGIN / (T * b)``, with ``T = t_max - t_min``;
        - the one that damps the echo, at decay rate ``g`` of the coherence,
          below ``tolerance``: ``g / (b * ln(1 / tolerance))``.

        The result is the smallest over all pairs of levels whose slopes
        differ. If none do, the classes do not matter and the window rule is
        used with ``b = 1``.

        Args:
            tolerance: largest echo amplitude allowed inside the window, as a
                fraction of the coherence. Default ``THERMAL_ECHO_TOLERANCE``.

        Returns:
            The spacing, a cyclic frequency.

        .. versionadded:: 0.15.0
        """
        tolerance = THERMAL_ECHO_TOLERANCE if tolerance is None else tolerance
        window = self.t_max - self.t_min
        slopes = self._level_doppler_slopes()
        decay = self._level_decay_rates()
        log_inv_tol = np.log(1.0 / tolerance)
        spacings = []
        for i in range(len(slopes)):
            for j in range(i + 1, len(slopes)):
                b = abs(slopes[i] - slopes[j])
                if b == 0.0:
                    continue
                by_window = THERMAL_ECHO_MARGIN / (window * b)
                by_decay = (decay[i] + decay[j]) / (b * log_inv_tol)
                spacings.append(max(by_window, by_decay))
        return min(spacings) if spacings else THERMAL_ECHO_MARGIN / window

    def _build_thermal_delta_list(self, vc: dict) -> np.ndarray:
        """Build the sorted array of thermal detunings.

        With ``thermal_delta_spacing`` this is one uniform grid from
        ``thermal_delta_min`` to ``thermal_delta_max`` (the step keys must be
        unset). Otherwise the outer and inner ranges are merged, with
        duplicates counted once.

        Args:
            vc: normalised velocity-class dict (from _normalise_velocity_classes).

        Returns:
            1-D array of thermal detuning values in angular-frequency units.
        """
        spacing = vc.get("thermal_delta_spacing")
        if spacing is not None:
            if spacing == "auto":
                spacing = self.thermal_delta_auto_spacing(
                    vc.get("thermal_delta_echo_tolerance")
                )
            span = vc["thermal_delta_max"] - vc["thermal_delta_min"]
            steps = int(np.ceil(span / spacing * (1.0 - 1e-12)))
            if steps + 1 > THERMAL_CLASSES_WARN:
                warnings.warn(
                    f"thermal_delta_spacing={spacing:.3g} gives {steps + 1} "
                    "velocity classes, each solved at every z-step. Check that "
                    "it is a cyclic frequency, like thermal_width.",
                    UserWarning,
                    stacklevel=3,
                )
            return (
                2
                * np.pi
                * np.linspace(
                    vc["thermal_delta_min"], vc["thermal_delta_max"], steps + 1
                )
            )
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
        self.thermal_widths = self._thermal_spacings(self.thermal_delta_list)
        return self.thermal_delta_list, self.thermal_weights

    @property
    def thermal_quadrature_weights(self) -> np.ndarray:
        """Weight of each velocity class in the Doppler average.

        The Maxwell-Boltzmann density ``thermal_weights`` times the detuning
        interval ``thermal_widths`` each class represents. Not normalised.

        .. versionadded:: 0.15.0
        """
        return self.thermal_weights * self.thermal_widths

    @staticmethod
    def _thermal_spacings(deltas: np.ndarray) -> np.ndarray:
        """Width of the detuning interval each velocity class represents.

        The classes are quadrature nodes for the Doppler integral, so each
        Boltzmann weight is multiplied by its interval width. Interior classes
        take half the distance between their neighbours (trapezoid rule). The
        end classes take their whole adjacent spacing, not half of it (a
        rectangle rule at the ends), so that the widths are constant, and the
        weights unchanged, on a uniform grid.

        Args:
            deltas: sorted 1-D array of thermal detunings.

        Returns:
            Array of widths, the same length as ``deltas``; ones for a single
            class.
        """
        deltas = np.asarray(deltas, dtype=float)
        if deltas.size < 2:
            return np.ones_like(deltas)
        gaps = np.diff(deltas)
        widths = np.empty_like(deltas)
        widths[1:-1] = 0.5 * (gaps[1:] + gaps[:-1])
        widths[0], widths[-1] = gaps[0], gaps[-1]
        return widths

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
        self._num_density_z_func_name = num_density_z_func or "square"
        if num_density_z_func:
            self.num_density_z_func = getattr(t_funcs, num_density_z_func)(0)
        else:
            self.num_density_z_func = t_funcs.square(0)
        # Set the num_density args
        self._num_density_z_args_raw = num_density_z_args or {
            "on": 0.0,
            "off": 1.0,
            "ampl": 1.0,
        }
        if num_density_z_args:
            self.num_density_z_args = {}
            for key, value in num_density_z_args.items():
                self.num_density_z_args[key + "_0"] = value
        else:
            self.num_density_z_args = {"on_0": 0.0, "off_0": 1.0, "ampl_0": 1.0}

    def check(self) -> bool:
        """Validates the PropagationSolver object."""
        if len(self.interaction_strengths) != len(self.system.fields):
            raise ValueError(
                "The number of interaction_strengths must match the number of fields."
            )
        return True

    @property
    def rabi_freqs_zt(self) -> np.ndarray:
        """The solved complex Rabi frequency of each field, as a cyclic frequency.

        Same units as ``rabi_freq`` in the problem definition. Shape
        ``(num_fields, z_steps+1, t_steps+1)``. This is a new array on every access
        (the solver stores the angular field), so assign it to a local variable
        before indexing it in a loop; modifying it does not change the solver.
        """
        return self._Omegas_zt / (2 * np.pi)

    # ``Omegas_zt`` is deprecated. It is handled here, not as a property, so that
    # introspection (``dir``, ``inspect.getmembers``) never evaluates it and warns.
    def __getattr__(self, name: str) -> Any:
        if name == "Omegas_zt":
            warn_deprecated("PropagationSolver.Omegas_zt", "rabi_freqs_zt (cyclic)")
            return self._Omegas_zt  # the angular array
        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute '{name}'"
        )

    def __setattr__(self, name: str, value: Any) -> None:
        if name == "Omegas_zt":
            warn_deprecated("PropagationSolver.Omegas_zt", "rabi_freqs_zt (cyclic)")
            name = "_Omegas_zt"
        super().__setattr__(name, value)

    def init_Omegas_zt(self) -> np.ndarray:
        """Inits the Rabi frequency array.

        _Omegas_zt is the internal angular-frequency field array;
        ``rabi_freqs_zt`` is its public cyclic view.

        _Omegas_zt shape: (num_fields, z_steps+1, t_steps+1) — field-first.
        states_zt shape: (z_steps+1, t_steps+1, num_states, num_states) — z-first.
        These are inconsistent; correcting either is a breaking API change,
        deferred to a future major version.
        """
        self._Omegas_zt = np.zeros(
            (len(self.system.fields), len(self.zlist), len(self.tlist)), dtype=complex
        )
        self._Omegas_z_buf = np.zeros(
            (len(self.system.fields), len(self.tlist)), dtype=complex
        )
        for f_i, f in enumerate(self.system.fields):
            self._Omegas_zt[f_i][0] = (
                2.0
                * np.pi
                * f.rabi_freq
                * f.rabi_freq_t_func(self.tlist, f.rabi_freq_t_args)
            )
        return self._Omegas_zt

    def init_states_zt(self) -> np.ndarray:
        """Inits the system density matrices."""
        self.states_zt = np.zeros(
            (
                len(self.zlist),
                len(self.tlist),
                self.system.num_states,
                self.system.num_states,
            ),
            dtype=complex,
        )
        return self.states_zt

    def _pbar_postfix(self, progress_show_area: bool) -> dict:
        postfix = {"Ω_max": f"{np.abs(self._Omegas_z_buf).max():.3g}"}
        if progress_show_area:
            for i, field in enumerate(self.system.fields):
                area = np.trapezoid(np.abs(self._Omegas_z_buf[i]), self.tlist)
                key = f"A({field.label})" if field.label else f"A{i}"
                postfix[key] = f"{area:.3g}"
        return postfix

    def solve(
        self,
        step: str = "ab",
        rho0: qu.Qobj | None = None,
        recalc: bool = True,
        progress: bool = True,
        progress_show_area: bool = False,
        check_counter_prop_depletion: bool = True,
        depletion_warn: float = 0.01,
        depletion_error: float = 0.10,
        check_stability: bool = True,
        blowup_factor: float = 1e3,
        step_error_warn: float = 0.1,
        max_z_steps_inner: int = 512,
        t_step_error_warn: float = 0.1,
        window_edge_warn: float = 0.2,
        window_edge_frac: float = 0.05,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Solves the propagation of the fields through the medium.

        Args:
            step: 'euler (for Euler method) or 'ab' (for Adams-Bashforth)
            rho0 (Qobj): the initial density matrix state
            recalc (bool): Recalculate the solution even if a savefile exists?
            progress: Show a tqdm progress bar with elapsed time, ETA and
                current max field amplitude.
            progress_show_area: Extend the progress bar with the pulse area ∫|Ω|dt for
                each field at the current z-step.  Requires progress=True.
            check_counter_prop_depletion: Run a post-solve depletion check on
                any counter-propagating fields. Set False to suppress.
            depletion_warn: Depletion fraction that triggers a UserWarning.
            depletion_error: Depletion fraction that raises
                CounterPropagatingDepletionError.
            check_stability: Stop with :class:`PropagationInstabilityError` as
                soon as the solution blows up (non-finite field, field larger
                than ``blowup_factor`` times its initial peak, or a failed
                master-equation integration), and warn if the z-step looks too
                coarse. Set False to suppress.
            blowup_factor: Growth of the peak field, relative to its peak at
                ``z_min``, treated as a blow-up. Raise it for media with
                genuine gain.
            step_error_warn: Estimated accumulated relative error (from the
                differences of dOmega/dz between z-steps, relative to the peak
                field) above which a UserWarning is issued
                suggesting a larger ``z_steps_inner``. The largest value seen
                is kept in ``self.step_error_estimate``. With
                ``z_steps_inner="auto"`` this is also the tolerance: the solve
                is repeated with more inner steps until the estimate is below
                it, or ``max_z_steps_inner`` is reached.
            max_z_steps_inner: Upper limit for ``z_steps_inner="auto"``.
            t_step_error_warn: With ``check_stability``, warn if the time grid
                is too coarse to represent the fields and coherences: the
                relative error (to the peak) of interpolating them from every
                other time sample, ``self.t_step_error_estimate``, exceeds this.
                The error of the full grid is smaller, by about a factor of 16
                once the grid resolves the signal (a cubic spline), and by
                less on a coarser grid.
            window_edge_warn: With ``check_stability``, warn if a field is still
                changing near the edge of the time window (the first or last 5%
                of it) by more than this fraction of its peak, so that a pulse
                or its ringing is probably cut off. Kept in
                ``self.window_edge_estimate``.
            window_edge_frac: Fraction of the time window at each end that
                ``window_edge_warn`` looks at. Reduce it if the window
                deliberately starts or ends in the middle of a pulse.

        Returns:
            Omegas_zt: The solved field complex Rabi frequency at each point in
                space z and time t, as an *angular* frequency (``rabi_freqs_zt``
                times ``2 * pi``).
            self.states_zt: The solved density matrix at each point in space z
                and time t.
        """
        self.check()
        self._stability = (
            {
                "blowup_factor": blowup_factor,
                "tol": step_error_warn,
                "ref": None,
                "err": 0.0,
                "order": 1,
            }
            if check_stability
            else None
        )
        self.step_error_estimate = 0.0
        self.t_step_error_estimate = 0.0
        self.window_edge_estimate = 0.0
        self.init_Omegas_zt()
        self.init_states_zt()
        # The z-stepping below replaces system.H_Omega_list with interpolated
        # coefficients. Start from the analytic pulses (picking up any edits to
        # the fields' t_func or t_args) and restore them afterwards, so that
        # solve() is repeatable.
        self.system._build_H_Omega()

        def _do_solve() -> None:
            if step == "euler":
                self._solve_euler(
                    rho0=rho0,
                    recalc=recalc,
                    progress=progress,
                    progress_show_area=progress_show_area,
                )
            elif step == "ab":
                self._solve_ab(
                    rho0=rho0,
                    recalc=recalc,
                    progress=progress,
                    progress_show_area=progress_show_area,
                )
            self.save_results()
            if self._stability is not None:
                self.step_error_estimate = self._stability["err"]
                self._warn_step_error(step_error_warn)

        try:
            # Should we recalculate or load a savefile?
            if recalc or not self.savefile_exists():
                if self._auto_z_steps_inner and self._stability is not None:
                    self._solve_auto(_do_solve, step_error_warn, max_z_steps_inner)
                else:
                    _do_solve()
            else:
                try:
                    self.load_results()
                except ValueError:
                    self.system._build_H_Omega()
                    _do_solve()
        finally:
            self.system._build_H_Omega()

        if check_stability:
            self._check_time_resolution(
                t_step_error_warn, window_edge_warn, window_edge_frac
            )

        if check_counter_prop_depletion:
            self._check_counter_prop_depletion(
                depletion_warn=depletion_warn, depletion_error=depletion_error
            )

        return self._Omegas_zt, self.states_zt

    def to_autodiff(self, **options: Any) -> Any:
        """Return a differentiable version of this problem.

        .. versionadded:: 0.14.0

        The result reproduces :meth:`solve` as a pure function of a parameter
        pytree (pulse parameters, interaction strengths, detunings, decay scale,
        or sampled input fields), so the propagated field and any loss of it can
        be differentiated exactly with JAX autodiff. Requires the ``jax`` extra
        (``pip install 'clerq[jax]'``) and ``jax_enable_x64``.

        Args:
            **options: forwarded to
                :class:`clerq.autodiff.AutodiffPropagation` (``step``,
                ``strategy``, ``rtol``, ``atol``, ``max_steps``, ...).

        Returns:
            The differentiable problem, an ``AutodiffPropagation``.

        Raises:
            ImportError: if the ``jax`` extra is not installed.
            NotImplementedError: for counter-propagating fields or pulse
                functions the autodiff backend does not support.

        Examples:
            >>> ad = solver.to_autodiff()  # doctest: +SKIP
            >>> value, grads = ad.value_and_grad(loss)  # doctest: +SKIP
        """
        from clerq.autodiff import AutodiffPropagation

        return AutodiffPropagation.from_solver(self, **options)

    def mbsolve(self, *args, **kwargs) -> tuple[np.ndarray, np.ndarray]:
        """Deprecated name for :meth:`PropagationSolver.solve`."""
        warn_deprecated("mbsolve()", "solve()")
        return self.solve(*args, **kwargs)

    def mbsolve_euler(self, *args, **kwargs) -> tuple[np.ndarray, np.ndarray]:
        """Deprecated; use :meth:`PropagationSolver.solve` with ``step="euler"``."""
        warn_deprecated("mbsolve_euler()", 'solve(step="euler")')
        return self._solve_euler(*args, **kwargs)

    def mbsolve_ab(self, *args, **kwargs) -> tuple[np.ndarray, np.ndarray]:
        """Deprecated; use :meth:`PropagationSolver.solve` with ``step="ab"``."""
        warn_deprecated("mbsolve_ab()", 'solve(step="ab")')
        return self._solve_ab(*args, **kwargs)

    def _solve_euler(
        self,
        rho0: qu.Qobj | None = None,
        recalc: bool = True,
        progress: bool = False,
        progress_show_area: bool = False,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Solves the propagation using a Euler step.

        Args:
            rho0 (Qobj): the initial density matrix state
            recalc (bool): Recalculate the solution even if a savefile exists?
            progress: Show a tqdm progress bar.
            progress_show_area: Add pulse area ∫|Ω|dt per field to the progress bar.

        Returns:
            Omegas_zt: The solved field complex Rabi frequency at each point in
                space z and time t, as an *angular* frequency (``rabi_freqs_zt``
                times ``2 * pi``).
            self.states_zt: The solved density matrix at each point in space z
                and time t.
        """
        h = self.z_step_inner()
        N_z = np.array(
            [
                self.num_density_z_func(
                    self.z_min + (i + 1) * h, self.num_density_z_args
                )
                for i in range(self.z_steps * self.z_steps_inner)
            ]
        )
        self.states_zt[0, :] = self._solve_and_average_over_thermal_detunings(
            z=self.z_min
        )
        pbar = tqdm.tqdm(
            enumerate(self.zlist[:-1]),
            total=self.z_steps,
            disable=not progress,
            unit="z",
        )
        f_hist: list[np.ndarray] = []
        for j, _z in pbar:
            Omegas_z_this = self._Omegas_zt[:, j, :]
            for k in range(self.z_steps_inner):
                N = N_z[j * self.z_steps_inner + k]
                thermal_states_t = self._solve_and_average_over_thermal_detunings(
                    z=_z + k * h
                )
                sum_coh_this = self.system.get_fields_sum_coherence(
                    states_t=thermal_states_t
                )
                Omegas_z_next = self._z_step_fields_euler(
                    h=h,
                    N=N,
                    Omegas_z_this=Omegas_z_this,
                    sum_coh_this=sum_coh_this,
                )
                f_hist.insert(0, 2.0j * N * self.g[:, None] * sum_coh_this)
                del f_hist[2:]
                self._check_stability_step(
                    Omegas_z_next, z=_z + (k + 1) * h, f_hist=f_hist, h=h, order=1
                )
                np.copyto(self._Omegas_z_buf, Omegas_z_next)
                self.system.H_Omega_list = self._build_intp_H_Omega_list(
                    self._Omegas_z_buf
                )
                Omegas_z_this = self._Omegas_z_buf
            self.states_zt[j + 1, :] = self.states_t()
            self._Omegas_zt[:, j + 1, :] = self._Omegas_z_buf
            pbar.set_postfix(self._pbar_postfix(progress_show_area))
        return self._Omegas_zt, self.states_zt

    def _solve_ab(
        self,
        rho0: qu.Qobj | None = None,
        recalc: bool = True,
        progress: bool = False,
        progress_show_area: bool = False,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Solves the propagation using an Adams-Bashforth step.

        Args:
            rho0 (Qobj): the initial density matrix state
            recalc (bool): Recalculate the solution even if a savefile exists?
            progress: Show a tqdm progress bar.
            progress_show_area: Add pulse area ∫|Ω|dt per field to the progress bar.

        Returns:
            Omegas_zt: The solved field complex Rabi frequency at each point in
                space z and time t, as an *angular* frequency (``rabi_freqs_zt``
                times ``2 * pi``).
            self.states_zt: The solved density matrix at each point in space z
                and time t.
        """
        h = self.z_step_inner()
        N_z = np.array(
            [
                self.num_density_z_func(
                    self.z_min + (i + 1) * h, self.num_density_z_args
                )
                for i in range(self.z_steps * self.z_steps_inner)
            ]
        )
        thermal_states_t = self._solve_and_average_over_thermal_detunings(z=self.z_min)
        self.states_zt[0, :] = thermal_states_t
        sum_coh_prev = self.system.get_fields_sum_coherence(states_t=thermal_states_t)
        pbar = tqdm.tqdm(total=self.z_steps, disable=not progress, unit="z")
        f_hist: list[np.ndarray] = []
        for j in range(self.z_steps):
            Omegas_z_this = self._Omegas_zt[:, j, :]
            for k in range(self.z_steps_inner):
                N = N_z[j * self.z_steps_inner + k]
                if j == 0 and k == 0:
                    # Very first inner step: Euler bootstrap (no previous
                    # coherence yet), using the coherence at z_min.
                    sum_coh_this = sum_coh_prev
                    Omegas_z_next = self._z_step_fields_euler(
                        h=h,
                        N=N,
                        Omegas_z_this=Omegas_z_this,
                        sum_coh_this=sum_coh_this,
                    )
                else:
                    thermal_states_t = self._solve_and_average_over_thermal_detunings(
                        z=self.zlist[j] + k * h
                    )
                    sum_coh_this = self.system.get_fields_sum_coherence(
                        states_t=thermal_states_t
                    )
                    Omegas_z_next = self._z_step_fields_ab(
                        h=h,
                        N=N,
                        sum_coh_prev=sum_coh_prev,
                        sum_coh_this=sum_coh_this,
                        Omegas_z_this=Omegas_z_this,
                    )
                if not (j == 0 and k == 0):
                    f_hist.insert(0, 2.0j * N * self.g[:, None] * sum_coh_this)
                    del f_hist[3:]
                self._check_stability_step(
                    Omegas_z_next,
                    z=self.zlist[j] + (k + 1) * h,
                    f_hist=f_hist,
                    h=h,
                    order=2,
                )
                np.copyto(self._Omegas_z_buf, Omegas_z_next)
                self.system.H_Omega_list = self._build_intp_H_Omega_list(
                    self._Omegas_z_buf
                )
                Omegas_z_this = self._Omegas_z_buf
                sum_coh_prev = sum_coh_this
            self.states_zt[j + 1, :] = self.states_t()
            self._Omegas_zt[:, j + 1, :] = self._Omegas_z_buf
            pbar.update(1)
            pbar.set_postfix(self._pbar_postfix(progress_show_area))
        pbar.close()
        return self._Omegas_zt, self.states_zt

    def _solve_auto(self, do_solve, tol: float, max_inner: int) -> None:
        """Repeat do_solve with a larger z_steps_inner until it is stable."""
        while True:
            self._stability["err"] = 0.0
            self._stability["ref"] = None
            try:
                do_solve()
            except PropagationInstabilityError as e:
                nxt = e.suggested_z_steps_inner
            else:
                if self._stability["err"] <= tol:
                    return
                nxt = self._suggest_inner(self._stability["err"], tol)
            if self.z_steps_inner >= max_inner:
                raise PropagationInstabilityError(
                    f"z_steps_inner='auto' reached max_z_steps_inner="
                    f"{max_inner} without a stable solution; reduce the "
                    f"z extent or density, or raise t_steps.",
                    z=None,
                    suggested_z_steps_inner=nxt,
                )
            self.z_steps_inner = min(max(nxt, self.z_steps_inner + 1), max_inner)
            self.init_Omegas_zt()
            self.init_states_zt()
            self.system._build_H_Omega()

    def _check_stability_step(
        self,
        Omegas_next: np.ndarray,
        z: float,
        f_hist: list[np.ndarray],
        h: float,
        order: int,
    ) -> None:
        """Stop the solve if the field has blown up; accumulate the step error.

        ``f_hist`` holds dOmega/dz at the current and previous steps, newest
        first. The local error is estimated from the differences of f,
        ``h|f0 - f1|/2`` for Euler (order 1) and ``5h|f0 - 2 f1 + f2|/12`` for
        Adams-Bashforth (order 2), relative to each field's peak amplitude, and
        summed over steps to give the accumulated relative error.
        """
        st = self._stability
        if st is None:
            return
        peaks = np.max(np.abs(Omegas_next), axis=1)
        if not np.all(np.isfinite(peaks)):
            self._raise_instability(z, "the field became non-finite")
        if st["ref"] is None:
            st["ref"] = float(np.max(np.abs(self._Omegas_zt[:, 0, :])))
        ref = st["ref"]
        if ref > 0.0 and peaks.max() > st["blowup_factor"] * ref:
            self._raise_instability(
                z,
                f"the peak field grew to {peaks.max():.3g}, over "
                f"{st['blowup_factor']:g} times its initial peak {ref:.3g}",
            )
        st["order"] = order
        if order == 1 and len(f_hist) >= 2:
            diff = 0.5 * h * np.abs(f_hist[0] - f_hist[1])
        elif order == 2 and len(f_hist) >= 3:
            diff = 5.0 / 12.0 * h * np.abs(f_hist[0] - 2 * f_hist[1] + f_hist[2])
        else:
            return
        scale = np.maximum(peaks, 1e-12 * max(peaks.max(), 1e-300))
        st["err"] += float((np.max(diff, axis=1) / scale).max())

    def _raise_instability(self, z: float, reason: str) -> None:
        err = self._stability["err"]
        suggest = (
            self._suggest_inner(err, self._stability["tol"])
            if err > 0
            else self.z_steps_inner * 2
        )
        raise PropagationInstabilityError(
            f"Propagation became unstable at z={z:.6g}: {reason}. The z-step "
            f"is probably too coarse for this medium; try z_steps_inner="
            f"{suggest} (currently {self.z_steps_inner}), z_steps_inner='auto', "
            f"or more t_steps. Pass check_stability=False to solve() to "
            f"suppress this check.",
            z=z,
            suggested_z_steps_inner=suggest,
        )

    def _suggest_inner(self, err: float, tol: float) -> int:
        """z_steps_inner expected to bring the accumulated error down to tol.

        The accumulated error scales as h**order.
        """
        order = (self._stability or {}).get("order", 1)
        factor = max((max(err, tol) / tol) ** (1.0 / order), 2.0)
        return int(np.ceil(self.z_steps_inner * factor))

    # The time grid is used through a cubic spline (QuTiP's array coefficient),
    # whose interpolation error falls as dt**4 once the grid resolves the
    # signal; coarser than that it falls more slowly, so the order is measured.
    _T_ORDER = 4

    @staticmethod
    def _decimation_error(y: np.ndarray, peak: np.ndarray, stride: int = 2) -> float:
        """Relative error of rebuilding time samples from every ``stride``-th one.

        Interpolates with the same not-a-knot cubic spline that QuTiP uses for
        array coefficients and compares at the samples left out. ``y`` has time
        on the last axis; ``peak`` broadcasts against ``y[..., 0]``. With
        ``stride=2`` the error of the full grid is smaller by ``2**order``,
        with the order up to 4 for a grid that resolves the signal and lower
        for a coarser one. Only samples between the first and last knot are
        compared, so nothing is extrapolated.
        """
        n = y.shape[-1]
        if n < 3 * stride + 1:
            return 0.0
        t = np.arange(n, dtype=float)
        spline = CubicSpline(
            t[::stride], y[..., ::stride], axis=-1, bc_type="not-a-knot"
        )
        # Compare only between the first and last knot: the samples after the
        # last knot (when n - 1 is not a multiple of the stride) would be
        # extrapolated.
        last = (n - 1) // stride * stride
        rest = np.ones(n, dtype=bool)
        rest[::stride] = False
        rest[last + 1 :] = False
        err = np.max(np.abs(spline(t[rest]) - y[..., rest]), axis=-1)
        return float(np.max(err / peak))

    def _check_time_resolution(
        self, tol: float, edge_tol: float, edge_frac: float
    ) -> None:
        """Warn if the time grid is too coarse, or the window cuts off a field.

        Sets ``t_step_error_estimate`` (decimation error of the fields and
        summed coherences, relative to each field's peak) and
        ``window_edge_estimate`` (largest change of a field, relative to its
        peak, over the first or last ``edge_frac`` of the window).
        """
        Omegas = self._Omegas_zt  # (num_fields, z, t)
        peak_O = np.max(np.abs(Omegas), axis=(1, 2))
        if not np.all(np.isfinite(peak_O)) or not np.any(peak_O > 0.0):
            return
        peak_O = np.where(peak_O > 0.0, peak_O, np.inf)
        err = self._decimation_error(Omegas, peak_O[:, None])
        err4 = self._decimation_error(Omegas, peak_O[:, None], stride=4)
        n_z, n_t = self.states_zt.shape[:2]
        sum_coh = self.system.get_fields_sum_coherence(
            states_t=self.states_zt.reshape(n_z * n_t, *self.states_zt.shape[2:])
        ).reshape(-1, n_z, n_t)  # (num_fields, z, t)
        peak_c = np.max(np.abs(sum_coh), axis=(1, 2))
        if np.all(np.isfinite(peak_c)):
            # Coherences near the integrator's absolute tolerance are noise,
            # not structure the grid should resolve: normalise to at least this.
            floor = 1e3 * (self.opts or {}).get("atol", 1e-8)
            peak = np.maximum(peak_c, floor)[:, None]
            err_c = self._decimation_error(sum_coh, peak)
            if err_c > err:
                # Take the order from the same quantity as the error.
                err, err4 = err_c, self._decimation_error(sum_coh, peak, stride=4)
        self.t_step_error_estimate = err

        m = max(1, int((Omegas.shape[-1] - 1) * edge_frac))
        change = np.maximum(
            np.abs(Omegas[..., m] - Omegas[..., 0]),
            np.abs(Omegas[..., -1] - Omegas[..., -1 - m]),
        )
        self.window_edge_estimate = float(np.max(np.max(change, axis=1) / peak_O))

        if err > tol:
            # Order of convergence seen between the two strides, 1 to 4.
            order = self._T_ORDER
            if err > 0.0 and err4 > err:
                order = float(np.clip(np.log2(err4 / err), 1.0, self._T_ORDER))
            factor = max((err / tol) ** (1.0 / order), 1.5)
            warnings.warn(
                f"The time grid may be too coarse: the fields or coherences "
                f"are interpolated from every other time sample with a relative "
                f"error of {err:.1%} (tolerance {tol:.1%}). Try "
                f"t_steps={int(np.ceil(self.t_steps * factor))} (currently "
                f"{self.t_steps}), and repeat if the warning remains; "
                f"suggest_t_steps() gives an estimate from the problem's "
                f"frequency scales.",
                UserWarning,
                stacklevel=3,
            )
        if self.window_edge_estimate > edge_tol:
            warnings.warn(
                f"A field is still changing near the edge of the time window "
                f"(by {self.window_edge_estimate:.1%} of its peak), so a pulse "
                f"or its ringing may be cut off. Try widening t_min and t_max.",
                UserWarning,
                stacklevel=3,
            )

    def suggest_t_steps(
        self, points_per_period: float = 10.0, bandwidth_fraction: float = 0.99
    ) -> int:
        """Estimate the ``t_steps`` needed to resolve the problem's frequencies.

        Finds the largest frequency scale (cyclic, in the units of the problem)
        and returns enough steps over ``[t_min, t_max]`` for
        ``points_per_period`` samples per period of it. The scales are the
        bandwidth of the input fields (``bandwidth_fraction`` of their power), their peak Rabi
        frequencies combined with the detunings and the Doppler width.

        This is an estimate from the input fields: the medium can add structure
        (sidebands, ringing), so check the solved result too; ``solve()`` warns
        if the grid turns out to be too coarse.

        Args:
            points_per_period: Time samples per period of the fastest scale.
            bandwidth_fraction: Fraction of an input field's spectral power
                that its bandwidth must contain. Lower it for sharp-edged
                pulses, whose slowly decaying tails would otherwise dominate.

        Returns:
            The suggested number of time steps; compare it with ``t_steps``.
        """
        # The input fields as they are now, not as last built: the fields may
        # have been edited since. Cyclic, as in the problem definition.
        Omegas = np.array(
            [
                f.rabi_freq * f.rabi_freq_t_func(self.tlist, f.rabi_freq_t_args)
                for f in self.system.fields
            ],
            dtype=complex,
        ).reshape(len(self.system.fields), len(self.tlist))
        dt = self.tlist[1] - self.tlist[0]
        rabi = float(np.max(np.abs(Omegas))) if Omegas.size else 0.0
        detuning = max(
            (abs(d) for d in self.system.get_detunings()), default=0.0
        ) + float(np.max(np.abs(self.thermal_delta_list)) / (2 * np.pi))
        f_max = float(np.hypot(rabi, detuning))
        # The envelope is complex, so its spectrum is two-sided: order by |f|.
        freqs = np.abs(np.fft.fftfreq(Omegas.shape[-1], d=dt))
        order = np.argsort(freqs)
        for field in Omegas:
            power = np.abs(np.fft.fft(field - field[0])) ** 2
            total = power.sum()
            if total > 0.0:
                cum = np.cumsum(power[order]) / total
                f_max = max(
                    f_max, float(freqs[order][np.searchsorted(cum, bandwidth_fraction)])
                )
        span = self.t_max - self.t_min
        return max(int(np.ceil(points_per_period * f_max * span)), 1)

    def _warn_step_error(self, tol: float) -> None:
        err = self.step_error_estimate
        if err > tol and not self._auto_z_steps_inner:
            warnings.warn(
                f"The z-step may be too coarse: estimated accumulated relative "
                f"error is {err:.1%} (tolerance {tol:.1%}). Try z_steps_inner="
                f"{self._suggest_inner(err, tol)} (currently "
                f"{self.z_steps_inner}) or z_steps_inner='auto'.",
                UserWarning,
                stacklevel=4,
            )

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
        # Factor of 2: the Maxwell equation is ∂Ω/∂z = 2iNg·ρ when the Bloch
        # Hamiltonian uses the standard coupling Ω/2 (vs Ω in the thesis convention).
        dOmega_dz = 2.0j * N * self.g[:, None] * sum_coh_this
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
        dOmega_dz_this = 2.0j * N * self.g[:, None] * sum_coh_this
        dOmega_dz_prev = 2.0j * N * self.g[:, None] * sum_coh_prev
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
        for f_i, f in enumerate(self.system.fields):
            # Build the sigma structure with rabi_freq=1 (same as _build_H_Omega)
            H_Omega = qu.Qobj(
                np.zeros([self.system.num_states, self.system.num_states])
            )
            for c_i, c in enumerate(f.coupled_levels):
                H_Omega += (
                    self.system.sigma(a=c[0], b=c[1])
                    + self.system.sigma(a=c[1], b=c[0])
                ) * f.factors[c_i]
            H_Omega = H_Omega * np.pi  # pi * 1.0 (rabi_freq normalised to 1)
            # Divide by 2π: H_Omega * coeff(t) = sigma * pi * Omega_z/(2π)
            #                                   = sigma * Omega_z/2  ✓
            # Omega sigma_ab + conj(Omega) sigma_ba
            #     = Re(Omega) (sigma_ab + sigma_ba) + Im(Omega) i (sigma_ab - sigma_ba)
            coeff = qu.coefficient(Omegas_z[f_i].real / (2 * np.pi), tlist=self.tlist)
            H_Omega_list.append([H_Omega, coeff])
            max_im = np.abs(Omegas_z[f_i].imag).max()
            if not self.real_drive and max_im > 1e-14 * np.abs(Omegas_z).max():
                H_Omega_im = qu.Qobj(
                    np.zeros([self.system.num_states, self.system.num_states])
                )
                for c_i, c in enumerate(f.coupled_levels):
                    H_Omega_im += (
                        self.system.sigma(a=c[0], b=c[1])
                        - self.system.sigma(a=c[1], b=c[0])
                    ) * (1j * f.factors[c_i])
                H_Omega_im = H_Omega_im * np.pi
                coeff_im = qu.coefficient(
                    Omegas_z[f_i].imag / (2 * np.pi), tlist=self.tlist
                )
                H_Omega_list.append([H_Omega_im, coeff_im])
        return H_Omega_list

    def _solve_and_average_over_thermal_detunings(
        self, z: float | None = None
    ) -> np.ndarray:
        """Solves the Lindblad equation for the System over a range of
            detuning shifts for velocity classes.

        Returns:
            A states_t object that is the Maxwell-Boltzmann weighted average
            over the velocity classes.
        """
        states_t_Delta = np.zeros(
            (
                len(self.thermal_delta_list),
                len(self.tlist),
                self.system.num_states,
                self.system.num_states,
            ),
            dtype=complex,
        )
        # The set detunings, without any thermal shifting
        fixed_detunings = self.system.get_detunings()
        for Delta_i, Delta in enumerate(self.thermal_delta_list):
            # Shift each field's detuning by Delta scaled by that field's
            # Doppler sign. Counter-propagating fields get factor_doppler_shift=-1.
            # thermal_delta_list is angular, but field detunings are cyclic
            # (set_H_Delta multiplies by 2*pi), so convert back to cyclic.
            self.system.set_H_Delta(
                [
                    fd + Delta / (2 * np.pi) * field.factor_doppler_shift
                    for fd, field in zip(
                        fixed_detunings, self.system.fields, strict=False
                    )
                ]
            )
            # We don't want the single-point solve to save.
            try:
                self._solve_system(opts=None, save=False)
            except IntegratorException as e:
                if self._stability is None or z is None:
                    raise
                self._raise_instability(
                    z, f"the master-equation integrator failed ({e})"
                )
            states_t_Delta[Delta_i] = self.states_t()
        # Restore fixed detunings
        self.system.set_H_Delta(fixed_detunings)
        thermal_states_t = np.average(
            states_t_Delta, axis=0, weights=self.thermal_quadrature_weights
        )
        return thermal_states_t

    def _check_counter_prop_depletion(
        self,
        depletion_warn: float = 0.01,
        depletion_error: float = 0.10,
    ) -> None:
        """Check counter-propagating fields for significant depletion.

        For each field with factor_doppler_shift < 0, computes the peak-based
        depletion across the cell:

            depletion = 1 - max_t|Ω(z_max, t)| / max_t|Ω(z_min, t)|

        Emits a UserWarning at depletion_warn; raises
        CounterPropagatingDepletionError at depletion_error.

        The forward-propagation solver is only equivalent to true
        counter-propagation when depletion is negligible. Results at high
        depletion levels may be unreliable.
        """
        for field_idx, field in enumerate(self.system.fields):
            if field.factor_doppler_shift >= 0:
                continue
            Omega_field = self._Omegas_zt[field_idx]  # shape (z_steps+1, t_steps+1)
            peak_zmin = np.max(np.abs(Omega_field[0]))
            peak_zmax = np.max(np.abs(Omega_field[-1]))
            if peak_zmin == 0.0:
                continue
            depletion = 1.0 - peak_zmax / peak_zmin
            name = field.label or f"field[{field_idx}]"
            msg = (
                f"Counter-propagating field '{name}' depleted by {depletion:.1%} "
                f"across the cell. The forward-propagation solver is only equivalent "
                f"to true counter-propagation when depletion is negligible. Results "
                f"at this depletion level may be unreliable."
            )
            if depletion >= depletion_error:
                raise CounterPropagatingDepletionError(
                    msg + " Pass check_counter_prop_depletion=False to solve() to "
                    "suppress this check."
                )
            elif depletion >= depletion_warn:
                warnings.warn(msg, UserWarning, stacklevel=3)

    def get_json_dict(self) -> dict:
        """Return the full problem definition as a JSON-serialisable dict."""
        d = super().get_json_dict()
        d.update(
            {
                "z_min": self.z_min,
                "z_max": self.z_max,
                "z_steps": self.z_steps,
                # "auto" is part of the problem definition, the value it
                # resolves to is a result (kept in the savefile metadata).
                "z_steps_inner": (
                    "auto" if self._auto_z_steps_inner else self.z_steps_inner
                ),
                "num_density_z_func": self._num_density_z_func_name,
                "num_density_z_args": self._num_density_z_args_raw,
                "interaction_strengths": self.interaction_strengths,
                "velocity_classes": self.velocity_classes,
            }
        )
        if self.real_drive:  # only when set, so default JSON and hashes are unchanged
            d["real_drive"] = True
        return d

    def _savefile_metadata(self) -> dict:
        """Build the metadata dict stored alongside results in the .qu file."""
        return {
            "hash": hashlib.sha256(self.to_json_str().encode()).hexdigest(),
            "source_weights": _SOURCE_WEIGHTS,
            "thermal_quadrature": True,
            "drive": self._drive_name(),
            "z_steps_inner": self.z_steps_inner,
            "clerq": importlib.metadata.version("clerq"),
            "qutip": importlib.metadata.version("qutip"),
            "timestamp": datetime.datetime.now(datetime.UTC).isoformat(),
        }

    def save_results(self) -> None:
        """Saves the solution to a QuTiP pickle file.

        Notes:
            - The path to which the results will be saved is taken from
                self.savefile.
        """
        if self.savefile:
            os.makedirs(os.path.dirname(self.savefile) or ".", exist_ok=True)
            print("Saving PropagationSolver to", self.savefile + ".qu")
            qu.qsave(
                (self._Omegas_zt, self.states_zt, self._savefile_metadata()),
                self.savefile,
            )

    def _has_nonunit_factors(self) -> bool:
        """True if any field couples a pair with a factor other than 1."""
        return any(f != 1.0 for fld in self.system.fields for f in fld.factors)

    def _check_source_weights(self, saved: str | None) -> None:
        """Rejects results solved with a different source-term weighting.

        Before the source term was weighted by the (signed) factor, as in the
        Hamiltonian, it was weighted by factor**2. The two agree only when all
        the field factors are 1, so older results are stale otherwise.
        """
        if saved != _SOURCE_WEIGHTS and self._has_nonunit_factors():
            raise ValueError(
                f"Savefile {self.savefile}.qu was solved with the old "
                "factor**2 source-term weights, which differ from the "
                "current ones for fields with factors other than 1. "
                "Delete it or set recalc=True."
            )

    def _check_thermal_quadrature(self, saved: bool | None) -> None:
        """Rejects results solved before the classes were weighted by width.

        Those averaged with the Boltzmann density alone, which differs from
        the current weights only on a non-uniform grid.
        """
        widths = np.asarray(self.thermal_widths, dtype=float)
        if not saved and not np.allclose(widths, widths[0], rtol=1e-9):
            raise ValueError(
                f"Savefile {self.savefile}.qu was solved with velocity-class "
                "weights that ignored the class spacing, which differ from "
                "the current ones for a non-uniform grid (e.g. "
                "thermal_delta_inner_*). Delete it or set recalc=True."
            )

    def _drive_name(self) -> str:
        return "real" if self.real_drive else "complex"

    def _drive_can_differ(self) -> bool:
        """True if the real-part and complex drives can give different results.

        They agree for real fields with zero static detunings, zero level
        energies and a mirror-symmetric class grid with mirrored equal weights.
        """
        if any(f.detuning != 0 for f in self.system.fields):
            return True
        if np.abs(self.system.H_0.full()).max() > 0:
            return True
        deltas = np.asarray(self.thermal_delta_list, dtype=float)
        weights = np.asarray(self.thermal_quadrature_weights, dtype=float)
        order = np.argsort(deltas)
        deltas, weights = deltas[order], weights[order]
        return not (
            np.allclose(deltas, -deltas[::-1], atol=1e-12)
            and np.allclose(weights, weights[::-1], rtol=1e-9)
        )

    def _check_drive(self, saved: str | None) -> None:
        """Rejects results solved with the other drive convention.

        Savefiles without a marker were solved with the old Re(Omega) drive.
        """
        saved = saved or "real"
        if saved != self._drive_name() and self._drive_can_differ():
            raise ValueError(
                f"Savefile {self.savefile}.qu was solved with the {saved}-envelope "
                f"drive, but this solver uses the {self._drive_name()} drive, which "
                "gives different results for detuned fields. Delete it or set "
                "recalc=True."
            )

    def load_results(self) -> None:
        """Loads the solution from a QuTiP pickle file.

        Raises:
            ValueError: if the savefile was built from a different problem
                definition (hash mismatch), or was solved with the old
                factor**2 source-term weights and a field has factors other
                than 1, or was solved with velocity-class weights that ignored
                the class spacing on a non-uniform grid, or was solved with
                the other drive convention (savefiles without a drive marker
                used ``Re(Omega)``) for a problem where the two differ (static
                detunings, level energies or an asymmetric class grid).

        Notes:
            - The path from which the results will be loaded is taken from
                self.savefile.
            - Old savefiles without metadata are loaded without integrity
                checking, with a warning, unless a field has factors other
                than 1 (see Raises).
        """
        data = qu.qload(self.savefile)
        if len(data) == 3:
            Omegas_zt, states_zt, meta = data
            current_hash = hashlib.sha256(self.to_json_str().encode()).hexdigest()
            if meta["hash"] != current_hash:
                raise ValueError(
                    f"Savefile {self.savefile}.qu was built from a different "
                    "problem definition. Delete it or set recalc=True."
                )
            self._check_source_weights(meta.get("source_weights"))
            self._check_thermal_quadrature(meta.get("thermal_quadrature"))
            self._check_drive(meta.get("drive"))
            if self._auto_z_steps_inner and meta.get("z_steps_inner"):
                self.z_steps_inner = meta["z_steps_inner"]
            for key in ("clerq", "qutip"):
                saved = meta.get(key)
                current = importlib.metadata.version(key)
                if saved and saved != current:
                    warnings.warn(
                        f"Savefile was built with {key}=={saved}, "
                        f"current version is {current}.",
                        stacklevel=2,
                    )
        else:
            self._check_source_weights(None)
            self._check_thermal_quadrature(None)
            self._check_drive(None)
            warnings.warn(
                f"Savefile {self.savefile}.qu has no integrity metadata "
                "(old format). Loading without hash check.",
                stacklevel=2,
            )
            Omegas_zt, states_zt = data
        self._Omegas_zt, self.states_zt = Omegas_zt, states_zt

    def fields_area(self) -> np.ndarray:
        """Gets the integrated pulse area of each field.

        Returns:
            np.array [num_fields, num_z_steps]: Integrated area of each field
            over time
        """
        return np.trapezoid(np.real(self._Omegas_zt), self.tlist, axis=2)

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
        f = self.system.fields[field_idx]
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
            Unlike ``System.get_fields_sum_coherence``, which operates on a
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
        return self.coherences(self.system.fields[field_idx].coupled_levels)
