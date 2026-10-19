# Troubleshooting

## Numerical stability

**`PropagationInstabilityError`, or a result that looks wrong or grows along the medium.**
The $z$-step is probably too coarse for the density and pulse. Raise `z_steps_inner`
(the error message suggests a value), or pass `z_steps_inner="auto"`. If a dense medium
still fails, lower the density or shorten the medium, or use more `z_steps`.

**"The z-step may be too coarse" warning.** The estimated error of the $z$-stepping,
`solver.step_error_estimate`, is above `step_error_warn`. Raise `z_steps_inner` as above.

**"The time grid may be too coarse" warning.** The fields or coherences are not
resolved on the time grid. Raise `t_steps`; `solver.suggest_t_steps()` estimates a value.

**"A field is still changing near the edge of the time window" warning.** A pulse or its
ringing is cut off by `t_min` or `t_max`. Widen the window; remember the medium delays
the pulse.

**A genuine gain medium raises the instability error.** Raise `blowup_factor` in `solve()`.

All of these checks can be turned off with `solve(check_stability=False)`. See
[Numerical Stability and Step Sizes](../usage/numerical-stability.ipynb) for examples and how
to read the estimates.

If you have any other problems with clerq, please raise [an issue on Github](
https://github.com/tpogden/clerq/issues).
