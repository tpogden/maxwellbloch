"""The user-facing object returned by ``PropagationSolver.to_autodiff()``."""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import diffrax as dfx
import jax
import numpy as np
from jax.flatten_util import ravel_pytree

from clerq.autodiff._config import require_float64
from clerq.autodiff._propagate import (
    TWO_PI,
    PropagationSetup,
    propagate_fields,
    setup_from_solver,
    validate_options,
)
from clerq.autodiff._pulses import NONDIFFERENTIABLE_ARGS, exclusive_groups

_PARAM_KEYS = (
    "pulse_args",
    "rabi_freq",
    "interaction_strengths",
    "detunings",
    "decay_scale",
    "decay_rates",
    "b_field",
    "omega0",
)


def _merge_pulse_args(default: dict, user: dict, name: str | None = None) -> dict:
    """Overlay ``user`` on ``default``, keeping each alternative pair consistent.

    Of each group of alternatives that pulse ``name`` has in ``EXCLUSIVE_ARGS``
    (``ampl``/``n_pi``, and for a sech ``width``/``fwhm``) a pulse takes one;
    supplying one must not leave the other from the defaults behind. With no
    ``name``, every known group applies.
    """
    merged = {**default, **user}
    for group in exclusive_groups(name):
        if any(k in user for k in group):
            for k in group:
                if k not in user:
                    merged.pop(k, None)
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


def _norm(tree) -> float:
    """The Euclidean norm of all the leaves of ``tree``."""
    return float(jax.numpy.linalg.norm(ravel_pytree(tree)[0]))


def _overlay(fixed: dict, params: dict, names=None) -> dict:
    """``params`` over ``fixed``; per-field pulse arguments are merged.

    ``names`` optionally gives each field's pulse name (see ``_merge_pulse_args``).
    """
    out = {**fixed, **params}
    if "pulse_args" in fixed and "pulse_args" in params:
        names = names or [None] * len(params["pulse_args"])
        out["pulse_args"] = [
            _merge_pulse_args(f, p, n)
            for f, p, n in zip(
                fixed["pulse_args"], params["pulse_args"], names, strict=True
            )
        ]
    return out


@dataclass
class OptimizeResult:
    """The outcome of :meth:`AutodiffPropagation.minimize`.

    Attributes:
        params: the optimised parameters, same structure as the ones passed in.
        value: the loss at ``params``.
        losses: the loss at the starting point and after each update, a 1-D array
            of length ``n_steps + 1``; its last entry is ``value``.
        n_steps: the number of updates taken (fewer than requested if ``gtol``
            was reached).
        converged: whether the gradient norm fell below ``gtol``.
    """

    params: dict
    value: float
    losses: np.ndarray
    n_steps: int
    converged: bool


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
        raises if asked, and holds them fixed by default. To differentiate
        the edge times use a ramp pulse (``ramp_on``, ``ramp_off``,
        ``ramp_onoff``, ``ramp_offon``), whose edges are smooth. Calling
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
    ``decay_rates``
        array ``(n_decays,)``, the cyclic ``rate`` of each entry of the system's
        ``decays`` (spontaneous decay, ground-state decoherence, or pure
        dephasing through a ``[a, a]`` channel). The effective rate of a channel is
        ``decay_scale * decay_rates[i]``. Differentiable at any value, including
        zero, so a rate can be fitted or optimised from nothing; keep it
        non-negative. Empty for a system without decays.
    ``b_field``
        scalar, the DC B-field in natural units (``mu_B B / h`` in the
        frequency unit of the energies), so the Zeeman shift of the levels can
        be differentiated: ``H_0 = H_base + 2 pi B Z``. Requires the system to
        have a ``b_field`` block giving ``g_F`` and ``m_F`` (it is then in the
        default parameters; without one it raises ``ValueError``).
    ``omega0``
        array ``(n_fields, n_t)``, the input field samples on ``tlist`` as a
        *cyclic* Rabi frequency, the same unit as ``rabi_freq``. If given, the
        first slice is driven by their cubic spline instead of the analytic
        pulse (and ``pulse_args`` / ``rabi_freq`` are unused).

    All parameters are cyclic (like the problem definition), and so is the
    propagated field. It is returned in the layout and units of
    ``PropagationSolver.rabi_freqs_zt``: shape ``(n_fields, z_steps + 1, n_t)``,
    complex, *cyclic* Rabi frequency, the same unit as ``rabi_freq`` and
    ``omega0``.

    Attributes:
        GRAD_CACHE_SIZE: how many compiled ``value_and_grad`` programs (one per
            loss function and mode) are kept; the least recently used is dropped.

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
        rtol: relative tolerance of the adaptive step-size control.
        atol: absolute tolerance of the adaptive step-size control.
        max_steps: the most steps one slice solve may take (the solve errors if
            it needs more); reverse mode sizes its default checkpoint count from
            it, so a smaller value lowers memory use for large problems.
        n_steps: fixed-step solver with this many steps per slice (useful for
            finite-difference reference gradients); must not be commensurate with
            the time grid.
        solver: a diffrax solver, default ``Tsit5``.
    """

    GRAD_CACHE_SIZE = 8

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
        solver: Any = None,
    ) -> None:
        validate_options(step, strategy, checkpoints)
        self.setup = setup
        self._pulse_names = [name for name, _, _ in setup.pulses]
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
        # compiled value_and_grad functions, least recently used evicted first
        self._grad_cache: OrderedDict = OrderedDict()
        self._default_optimizers: dict = {}

    def clear_cache(self) -> None:
        """Drop the compiled programs, releasing the loss functions (and the arrays
        they capture) that ``value_and_grad`` keeps alive."""
        self._compiled.clear()
        self._grad_cache.clear()
        self._default_optimizers.clear()

    @classmethod
    def from_solver(cls, solver: Any, **options: Any) -> AutodiffPropagation:
        """Build from a :class:`~clerq.PropagationSolver`.

        Args:
            solver: the problem to differentiate.
            **options: the constructor options of this class (``step``, ``strategy``,
                ``checkpoints``, ``rtol``, ``atol``, ``max_steps``, ``n_steps``,
                ``solver``).
        """
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
        ``b_field`` is included only if the system has a ``b_field`` block.
        """
        s = self.setup
        params = {
            "pulse_args": [{k: float(v) for k, v in a.items()} for _, a, _ in s.pulses],
            "rabi_freq": np.array([r for _, _, r in s.pulses], dtype=float),
            "interaction_strengths": np.asarray(s.g, dtype=float) / TWO_PI,
            "detunings": np.asarray(s.ops.detunings, dtype=float).copy(),
            "decay_scale": 1.0,
            "decay_rates": (
                np.asarray(s.ops.decays.rates, dtype=float).copy()
                if s.ops.decays is not None
                else np.zeros(0)
            ),
        }
        if s.ops.has_b_field:
            params["b_field"] = float(s.ops.b_field)
        return params

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
                    "by passing them in fixed=..., or use a ramp pulse "
                    "('ramp_on', 'ramp_off', 'ramp_onoff' or 'ramp_offon'), a "
                    "square-like pulse with smooth edges whose on/off times are "
                    "differentiable."
                )

    def _validate(self, params: dict) -> None:
        unknown = set(params) - set(_PARAM_KEYS)
        if unknown:
            raise ValueError(
                f"Unknown parameter(s) {sorted(unknown)}; valid: {list(_PARAM_KEYS)}."
            )
        if "b_field" in params and not self.setup.ops.has_b_field:
            raise ValueError(
                "b_field was given but the system has no b_field block, so it "
                "would have no effect; add a b_field block (with g_F and m_F)."
            )
        n_f, n_t = self.num_fields, len(self.setup.tlist)
        if "decay_rates" in params:
            n_d = (
                0 if self.setup.ops.decays is None else len(self.setup.ops.decays.rates)
            )
            if np.shape(params["decay_rates"]) != (n_d,):
                raise ValueError(
                    f"decay_rates must have shape ({n_d},) (one per entry of the "
                    f"system's decays), got {np.shape(params['decay_rates'])}."
                )
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
                _merge_pulse_args(d, u, n)
                for d, u, n in zip(
                    defaults["pulse_args"],
                    params["pulse_args"],
                    self._pulse_names,
                    strict=True,
                )
            ]
        # Everything is cyclic at this interface; the engine works in angular units,
        # so convert omega0 on the way in and the propagated field on the way out.
        omega0 = p.get("omega0")
        omegas = propagate_fields(
            self.setup,
            omega0=None if omega0 is None else TWO_PI * omega0,
            pulse_args=p["pulse_args"],
            rabi_freqs=p["rabi_freq"],
            g=p["interaction_strengths"],
            detunings=p["detunings"],
            decay_scale=p["decay_scale"],
            decay_rates=p["decay_rates"] if np.size(p["decay_rates"]) else None,
            b_field=p.get("b_field"),
            **{**self._options, **overrides},
        )
        return omegas / TWO_PI

    def propagate(self, params: dict | None = None) -> Any:
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
    ) -> tuple[Any, dict]:
        """Value and gradient of ``loss(Omega_zt)`` with respect to ``params``.

        Args:
            loss: scalar-valued function of the propagated field ``Omega_zt``,
                shape ``(n_fields, z_steps + 1, n_t)``, complex. Reuse the same
                function object between calls to reuse the compiled program: it
                is compiled once per function object, and only the last
                ``GRAD_CACHE_SIZE`` (8) are kept, so a new lambda on every call
                recompiles every time.
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
        if key in self._grad_cache:
            self._grad_cache.move_to_end(key)
        else:
            overrides = {"adjoint": dfx.ForwardMode()} if mode == "forward" else {}
            objective = self._objective(loss, overrides)

            if mode == "reverse":
                fn = jax.jit(jax.value_and_grad(objective))
            else:

                def with_aux(p, fx):
                    v = objective(p, fx)
                    return v, v

                jac = jax.jacfwd(with_aux, has_aux=True)
                fn = jax.jit(lambda p, fx: (lambda g_v: (g_v[1], g_v[0]))(jac(p, fx)))
            self._grad_cache[key] = fn
            while len(self._grad_cache) > self.GRAD_CACHE_SIZE:
                self._grad_cache.popitem(last=False)
        return self._grad_cache[key](params, fixed)

    def _objective(self, loss: Callable, overrides: dict) -> Callable:
        """``objective(params, fixed)``: ``loss`` of the propagated field."""

        def objective(p, fx):
            over = _overlay(fx, p, self._pulse_names)
            return loss(self._propagate(over, **overrides))

        return objective

    def _update_fn(self, opt: Any, key: Any, loss: Callable) -> Callable:
        """A compiled optimiser update, cached like the ``value_and_grad`` programs.

        ``opt`` is wrapped so that optimisers which take the loss value, the
        gradient and the loss itself (for example ``optax.lbfgs``, which line
        searches on it) get them; the others ignore them.
        """
        import optax

        cache_key = ("update", key, loss)
        if cache_key in self._grad_cache:
            self._grad_cache.move_to_end(cache_key)
        else:
            opt = optax.with_extra_args_support(opt)
            objective = self._objective(loss, {})

            def update(grads, state, p, value, fx):
                return opt.update(
                    grads,
                    state,
                    p,
                    value=value,
                    grad=grads,
                    value_fn=lambda q: objective(q, fx),
                )

            self._grad_cache[cache_key] = jax.jit(update)
            while len(self._grad_cache) > self.GRAD_CACHE_SIZE:
                self._grad_cache.popitem(last=False)
        return self._grad_cache[cache_key]

    def minimize(
        self,
        loss: Callable,
        params: dict | None = None,
        *,
        fixed: dict | None = None,
        optimizer: Any = None,
        learning_rate: float | None = None,
        steps: int = 100,
        gtol: float | None = None,
    ) -> OptimizeResult:
        """Minimise ``loss(Omega_zt)`` over ``params`` with an optax optimiser.

        Each step is one :meth:`value_and_grad` call and one optimiser update, so
        the loss is compiled once, on the first step (plus one final evaluation). Needs ``optax``, which the
        ``jax`` extra installs.

        Args:
            loss: scalar-valued function of the propagated field, as for
                :meth:`value_and_grad`; reuse one function object.
            params: starting point and the parameters to optimise (any subset of
                the keys). Default: every differentiable parameter of the problem.
            fixed: parameters held constant, as for :meth:`value_and_grad`.
            optimizer: an ``optax.GradientTransformation`` (for example
                ``optax.adam(1e-2)`` or ``optax.lbfgs()``). Pass the same object
                to reuse its compiled update between calls. Default
                ``optax.adam(learning_rate)``. Parameters differ in scale (a
                pulse area is of order 1, an interaction strength of order
                ``1e3``), so the right ``learning_rate`` depends on what is being
                optimised; rescale the parameters inside ``loss`` if needed.
            learning_rate: step size of the default Adam optimiser (0.05 if not
                given); not allowed together with ``optimizer``.
            steps: the most updates to take.
            gtol: stop when the gradient's global norm is below this (default:
                always take ``steps`` updates).

        Returns:
            The An :class:`~clerq.autodiff.OptimizeResult`, at the last point reached. That is the
            best one for a descending optimiser but need not be for one that
            overshoots; ``losses`` shows the path.

        Raises:
            ImportError: if ``optax`` is not installed.
            RuntimeError: if JAX is not in 64-bit mode.
            ValueError: as for :meth:`value_and_grad`, or if ``steps`` is not a
                positive integer, or both ``optimizer`` and ``learning_rate`` are
                given.
        """
        try:
            import optax
        except ImportError as e:
            raise ImportError(
                "minimize() needs optax: pip install 'clerq[jax]' (or pip install "
                "optax)"
            ) from e
        if not isinstance(steps, int) or steps < 1:
            raise ValueError(f"steps must be a positive integer, got {steps!r}.")
        if optimizer is not None and learning_rate is not None:
            raise ValueError(
                "Give either optimizer or learning_rate, not both: learning_rate "
                "only sets the step size of the default Adam optimiser."
            )
        require_float64()
        if optimizer is None:
            lr = 0.05 if learning_rate is None else learning_rate
            opt = self._default_optimizers.setdefault(lr, optax.adam(lr))
        else:
            opt = optimizer
        x = (
            self._differentiable_default_params()
            if params is None
            else _floatify(params)
        )
        fx = {} if fixed is None else _floatify(fixed)
        update = self._update_fn(opt, opt, loss)
        state = opt.init(x)
        values, converged = [], False
        for _ in range(steps):
            value, grads = self.value_and_grad(loss, x, fixed=fx)
            values.append(value)
            if gtol is not None and _norm(grads) < gtol:
                converged = True
                break
            updates, state = update(grads, state, x, value, fx)
            x = optax.apply_updates(x, updates)
        if not converged:
            # The loss at the point reached; the gradient there is only needed to
            # test gtol.
            if gtol is None:
                value = loss(self.propagate(_overlay(fx, x, self._pulse_names)))
            else:
                value, grads = self.value_and_grad(loss, x, fixed=fx)
                converged = _norm(grads) < gtol
            values.append(value)
        losses = np.array(jax.device_get(values), dtype=float)
        return OptimizeResult(
            params=x,
            value=float(losses[-1]),
            losses=losses,
            n_steps=len(losses) - 1,
            converged=converged,
        )
