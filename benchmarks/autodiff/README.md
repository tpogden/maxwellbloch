# clerq.autodiff benchmarks

A reproducible harness for the measurements on the docs page
[Autodiff performance and validation](../../docs/usage/autodiff-performance.qmd). It is not run by CI.

| file | experiment |
|-|-|
| `e1_cost.py` | E1: the time of one gradient against the number of parameters N (forward, reverse, forward-mode, finite differences through JAX; finite differences through clerq are estimated from one `solve()`) |
| `e2_memory.py` | E2: peak memory of one gradient against the number of z steps and velocity classes, for strategies `s2`, `s3`, `s4` (one process per cell) |
| `e3_accuracy.py` | E3: gradient accuracy against fixed-step finite differences, forward against reverse mode, field error against `solve()` |
| `e4_compile.py` | E4: compile time and first-call against steady-state time (fresh process per system) |
| `e5_time_to_solution.py` | E5: a Doppler-broadened memory design problem, free-form autodiff design against a 5-parameter baseline |
| `run_tier2.sh` | everything above, one after another (about 5 hours), plus E6: E1 on the commit before PR 361 and after it |
| `analyse.py` | tables (`results/tables.md`) and the figures in `docs/usage/figures` |
| `results/` | the raw JSONL records (one per cell, with machine and version metadata) and the logs |

Run one experiment, for example:

```bash
uv sync --extra jax --extra plot
uv run python benchmarks/autodiff/e1_cost.py --systems lambda_eit --ns 2 10 50 --tag e1_try
uv run python benchmarks/autodiff/analyse.py --no-figures
```

For clean numbers keep the machine on AC power, awake (`caffeinate -i -s -w <pid>`) and otherwise idle; the
load average is stored with every record. `--quick` shrinks the grids for a smoke test and the output tag
`smoke_*` is ignored by `analyse.py`. To compare against an older source tree (E6), run the same script with
`PYTHONPATH=<checkout>/src` and a different `--tag`; the harness does not depend on newer API.

The benchmarks parametrise a free-form probe control as `omega0[0] = basis @ coefficients` and call the
private `AutodiffPropagation._propagate(..., adjoint=ForwardMode())` for forward mode, so they may need a
small update if those internals change.
