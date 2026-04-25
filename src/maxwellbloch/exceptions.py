class CounterPropagatingDepletionError(RuntimeError):
    """Raised when a counter-propagating field depletes beyond the error threshold.

    The forward-propagation solver is only equivalent to true counter-propagation
    when depletion is negligible. Pass ``check_counter_prop_depletion=False`` to
    ``mbsolve()`` to suppress this check if you know what you are doing.
    """
