class CounterPropagatingDepletionError(RuntimeError):
    """Raised when a counter-propagating field depletes beyond the error threshold.

    The forward-propagation solver is only equivalent to true counter-propagation
    when depletion is negligible. Pass ``check_counter_prop_depletion=False`` to
    ``solve()`` to suppress this check if you know what you are doing.
    """


class PropagationInstabilityError(RuntimeError):
    """Raised when a propagation solve has blown up numerically.

    The field became non-finite or grew far beyond its initial amplitude, or the
    master-equation integrator failed, which is the usual signature of a
    z-step (or t-step) too coarse for the finite-difference scheme. Increase
    ``z_steps_inner`` (or ``z_steps``) or ``t_steps``, or pass
    ``z_steps_inner="auto"``. Pass ``check_stability=False`` to ``solve()``
    to suppress the check.

    Attributes:
        z: Position at which the solve was stopped.
        suggested_z_steps_inner: A larger ``z_steps_inner`` to try, or ``None``.
    """

    def __init__(self, message, z=None, suggested_z_steps_inner=None):
        super().__init__(message)
        self.z = z
        self.suggested_z_steps_inner = suggested_z_steps_inner
