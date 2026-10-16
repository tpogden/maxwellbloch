"""The user-facing object returned by ``PropagationSolver.to_autodiff()``."""

from __future__ import annotations

from collections.abc import Callable

import diffrax as dfx
import jax
import numpy as np

from clerq.autodiff._config import require_float64
from clerq.autodiff._propagate import (
    TWO_PI,
    PropagationSetup,
    propagate_fields,
    setup_from_solver,
)
from clerq.autodiff._pulses import NONDIFFERENTIABLE_ARGS

_PARAM_KEYS = (
    "pulse_args",
    "rabi_freq",
    "interaction_strengths",
    "detunings",
    "decay_scale",
    "omega0",
)


def _merge_pulse_args(default: dict, user: dict) -> dict:
    """Overlay ``user`` on ``default``, keeping each alternative pair consistent.

    A pulse takes either ``ampl`` or ``n_pi``, and a sech pulse either ``width`` or
    ``fwhm``; supplying one must not leave the other from the defaults behind.
    """
    merged = {**default, **user}
    for chosen, other in (("ampl", "n_pi"), ("n_pi", "ampl"), ("width", "fwhm")):
        if chosen in user and other not in user:
            merged.pop(other, None)
    if "fwhm" in user and "width" not in user:
        merged.pop("width", None)
    return merged


def _floatify(tree):
    """Convert integer and boolean leaves to float64.

    Problem definitions often use integer literals (``"centre": 1``); ``jax.grad``
    rejects integer inputs, and the parameters are physical real quantities anyway.
    """

    def convert(x):
        dtype = getattr(x, "dtype", None)
        if dtype is None:
            dtype = np.asarray(x).dtype
        if np.issubdtype(dtype, np.integer) or np.issubdtype(dtype, np.bool_):
            return jax.numpy.asarray(x, dtype=jax.numpy.float64)
        return x

    return jax.tree_util.tree_map(convert, tree)


def _overlay(fixed: dict, params: dict) -> dict:
    """``params`` over ``fixed``; per-field pulse arguments are merged."""
    out = {**fixed, **params}
    if "pulse_args" in fixed and "pulse_args" in params:
        out["pulse_args"] = [
            _merge_pulse_args(f, p)
            for f, p in zip(fixed["pulse_args"], params["pulse_args"], strict=True)
        ]
    return out


class AutodiffPropagation:
    """A differentiable version of a :class:`~clerq.PropagationSolver` problem.

    Build it with :meth:`clerq.PropagationSolver.to_autodiff`. It reproduces
    ``PropagationSolver.solve`` as a pure function of a *parameter pytree*, so the
    propagated field and any loss of it can be differentiated with respect to the
    parameters, with reverse- or forward-mode autodiff, and composed with
    ``jax.jit``, ``jax.vmap`` and optimisers.

    Parameters are a ``dict``; every key is optional and defaults to the value
    in the problem (see :meth:`default_params`):

    ``pulse_args``
        list, one dict per field, of the arguments of its analytic input pulse
        (``ampl``, ``centre``, ``fwhm``/``width``, ``n_pi``, ``on``, ``off``).
        The ``on``/``off`` edge times of a ``square`` pulse can be set but not
        differentiated (the field is discontinuous there, so its derivative is a
        delta function that autodiff would return as 0); ``value_and_grad``
        raises if asked, and holds them fixed by default. Calling
        ``jax.grad`` yourself through :meth:`propagate` with those traced
        gives exactly that silent 0, so keep them out of the differentiated
        parameters.
    ``rabi_freq``
        array ``(n_fields,)``, cyclic ``rabi_freq`` multiplying each pulse.
    ``interaction_strengths``
        array ``(n_fields,)``, as the ``interaction_strengths`` argument.
    ``detunings``
        array ``(n_fields,)``, cyclic fixed detuning of each field.
    ``decay_scale``
        scalar multiplying every decay rate.
    ``omega0``
        array ``(n_fields, n_t)``, the input field samples on ``tlist`` as a
        *cyclic* Rabi frequency, the same unit as ``rabi_freq``. If given, the
        first slice is driven by their cubic spline instead of the analytic
        pulse (and ``pulse_args`` / ``rabi_freq`` are unused).

    All parameters are cyclic (like the problem definition). The propagated field
    is returned in the layout **and units** of ``PropagationSolver.Omegas_zt``:
    shape ``(n_fields, z_steps + 1, n_t)``, complex, *angular* frequency; divide
    by ``2 * pi`` for the cyclic Rabi frequency.

    Note:
        Counter-propagating fields are not supported. The result matches
        ``PropagationSolver`` for every ``z_steps_inner``; versions of clerq before
        #328 stored Adams-Bashforth slices one inner step early when
        ``z_steps_inner > 1`` and will not match there.

    Args:
        setup: the problem, from ``setup_from_solver``.
        step: ``"ab"`` (Adams-Bashforth, default) or ``"euler"``.
        strategy: ``"s3"`` (default) recomputes a slice on the backward pass;
            ``"s2"`` stores every slice; ``"s4"`` additionally limits the
            solver's inner checkpoints to ``checkpoints``.
        checkpoints: inner checkpoint count, required (a positive integer) for
            ``strategy="s4"`` and not allowed with the other strategies.
        rtol, atol: adaptive step-size tolerances.
        max_steps: the most steps one slice solve may take (the solve errors if
            it needs more); reverse mode sizes its default checkpoint count from
            it, so a smaller value lowers memory use for large problems.
        n_steps: fixed-step solver with this many steps per slice (useful for
            finite-difference reference gradients); must not be commensurate with
            the time grid.
        solver: a diffrax solver, default ``Tsit5``.
    """

    def __init__(
        self,
        setup: PropagationSetup,
        *,
        step: str = "ab",
        strategy: str = "s3",
        checkpoints: int | None = None,
        rtol: float = 1e-8,
        atol: float = 1e-10,
        max_steps: int = 100_000,
        n_steps: int | None = None,
        solver=None,
    ) -> None:
        self.setup = setup
        self._options = {
            "step": step,
            "strategy": strategy,
            "checkpoints": checkpoints,
            "rtol": rtol,
            "atol": atol,
            "max_steps": max_steps,
            "n_steps": n_steps,
            "solver": solver,
        }
        self._compiled: dict = {}

    @classmethod
    def from_solver(cls, solver, **options) -> AutodiffPropagation:
        """Build from a :class:`~clerq.PropagationSolver` (see :class:`cls`)."""
        return cls(setup_from_solver(solver), **options)

    # ---- parameters ---------------------------------------------------------

    @property
    def num_fields(self) -> int:
        return self.setup.ops.num_fields

    @property
    def tlist(self) -> np.ndarray:
        return self.setup.tlist

    def default_params(self) -> dict:
        """The parameters of the problem this object was built from.

        Every value is a float, whatever literal the problem definition used.
        """
        s = self.setup
        return {
            "pulse_args": [{k: float(v) for k, v in a.items()} for _, a, _ in s.pulses],
            "rabi_freq": np.array([r for _, _, r in s.pulses], dtype=float),
            "interaction_strengths": np.asarray(s.g, dtype=float) / TWO_PI,
            "detunings": np.asarray(s.ops.detunings, dtype=float).copy(),
            "decay_scale": 1.0,
        }

    def _differentiable_default_params(self) -> dict:
        """``default_params()`` without the arguments that cannot be differentiated.

        These stay at the problem's values (``_propagate`` fills them in).
        """
        params = self.default_params()
        params["pulse_args"] = [
            {
                k: v
                for k, v in a.items()
                if k not in NONDIFFERENTIABLE_ARGS.get(name, ())
            }
            for (name, _, _), a in zip(
                self.setup.pulses, params["pulse_args"], strict=True
            )
        ]
        return params

    def _reject_nondifferentiable(self, params: dict) -> None:
        """Raise if ``params`` asks for the gradient of a discontinuity's location."""
        for k, ((name, _, _), a) in enumerate(
            zip(self.setup.pulses, params.get("pulse_args", []), strict=False)
        ):
            bad = [x for x in NONDIFFERENTIABLE_ARGS.get(name, ()) if x in a]
            if bad:
                raise ValueError(
                    f"pulse_args[{k}] contains {bad}, which place the edges of a "
                    f"{name!r} pulse. The field is discontinuous there, so its "
                    "derivative with respect to them is a delta function that "
                    "autodiff cannot represent (it would return 0). Hold them fixed "
                    "by passing them in fixed=..., or use a continuous pulse such "
                    "as 'ramp_onoff'."
                )

    def _validate(self, params: dict) -> None:
        unknown = set(params) - set(_PARAM_KEYS)
        if unknown:
            raise ValueError(
                f"Unknown parameter(s) {sorted(unknown)}; valid: {list(_PARAM_KEYS)}."
            )
        n_f, n_t = self.num_fields, len(self.setup.tlist)
        for key in ("rabi_freq", "interaction_strengths", "detunings"):
            if key in params and np.shape(params[key]) != (n_f,):
                raise ValueError(
                    f"{key} must have shape ({n_f},), got {np.shape(params[key])}."
                )
        if "omega0" in params and np.shape(params["omega0"]) != (n_f, n_t):
            raise ValueError(
                f"omega0 must have shape ({n_f}, {n_t}), got "
                f"{np.shape(params['omega0'])}."
            )
        if "pulse_args" in params and len(params["pulse_args"]) != n_f:
            raise ValueError(
                f"pulse_args must have one dict per field ({n_f}), got "
                f"{len(params['pulse_args'])}."
            )

    # ---- propagation ----------------------------------------------------------

    def _propagate(self, params: dict, **overrides):
        defaults = self.default_params()
        p = {**defaults, **params}
        if "pulse_args" in params:  # per-field partial overrides of the pulse arguments
            p["pulse_args"] = [
                _merge_pulse_args(d, u)
                for d, u in zip(
                    defaults["pulse_args"], params["pulse_args"], strict=True
                )
            ]
        # omega0 is cyclic like every other parameter; the engine works in angular
        # units (the layout of Omegas_zt).
        omega0 = p.get("omega0")
        return propagate_fields(
            self.setup,
            omega0=None if omega0 is None else TWO_PI * omega0,
            pulse_args=p["pulse_args"],
            rabi_freqs=p["rabi_freq"],
            g=p["interaction_strengths"],
            detunings=p["detunings"],
            decay_scale=p["decay_scale"],
            **{**self._options, **overrides},
        )

    def propagate(self, params: dict | None = None):
        """Propagate the fields; a pure, jit-compiled function of ``params``.

        Args:
            params: parameter pytree (see the class docstring); defaults to the
                problem's.

        Returns:
            ``Omega(z, t)``, shape ``(n_fields, z_steps + 1, n_t)``, complex.

        Raises:
            RuntimeError: if JAX is not in 64-bit mode.
            ValueError: for unknown keys or wrong shapes.
        """
        require_float64()
        params = {} if params is None else _floatify(params)
        self._validate(params)
        if "propagate" not in self._compiled:
            self._compiled["propagate"] = jax.jit(self._propagate)
        return self._compiled["propagate"](params)

    def value_and_grad(
        self,
        loss: Callable,
        params: dict | None = None,
        *,
        fixed: dict | None = None,
        mode: str = "reverse",
    ):
        """Value and gradient of ``loss(Omega_zt)`` with respect to ``params``.

        Args:
            loss: scalar-valued function of the propagated field ``Omega_zt``,
                shape ``(n_fields, z_steps + 1, n_t)``, complex. Reuse the same
                function object between calls to reuse the compiled program.
            params: the parameters to differentiate (any subset of the keys; the
                gradient has the same structure). Default: every differentiable
                parameter of the problem, i.e. ``default_params()`` without the
                edge times (``on``, ``off``) of ``square`` pulses, which stay
                fixed because the field's derivative with respect to a
                discontinuity's location is not defined.
            fixed: parameters held constant at the given values, with the same
                structure as ``params`` (for example the ``on``/``off`` of a
                square pulse, or a different target interaction strength). They
                are not differentiated.
            mode: ``"reverse"`` (cost ~5x one propagation regardless of the number
                of parameters) or ``"forward"`` (cost grows with the number of
                parameters; useful as a cross-check or for very few parameters).

        Returns:
            ``(value, grads)`` with ``grads`` of the same structure as ``params``.

        Raises:
            RuntimeError: if JAX is not in 64-bit mode.
            ValueError: for an unknown ``mode``, unknown keys, wrong shapes, or a
                request for the gradient of a ``square`` pulse's edge times.
        """
        require_float64()
        if mode not in ("reverse", "forward"):
            raise ValueError(f"Unknown mode {mode!r}; use 'reverse' or 'forward'.")
        params = (
            self._differentiable_default_params()
            if params is None
            else _floatify(params)
        )
        fixed = {} if fixed is None else _floatify(fixed)
        self._validate(params)
        self._validate(fixed)
        self._reject_nondifferentiable(params)

        key = (loss, mode)
        if key not in self._compiled:
            overrides = {"adjoint": dfx.ForwardMode()} if mode == "forward" else {}

            def objective(p, fx):
                return loss(self._propagate(_overlay(fx, p), **overrides))

            if mode == "reverse":
                fn = jax.jit(jax.value_and_grad(objective))
            else:

                def with_aux(p, fx):
                    v = objective(p, fx)
                    return v, v

                jac = jax.jacfwd(with_aux, has_aux=True)
                fn = jax.jit(lambda p, fx: (lambda g_v: (g_v[1], g_v[0]))(jac(p, fx)))
            self._compiled[key] = fn
        return self._compiled[key](params, fixed)
