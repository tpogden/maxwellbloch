#!/usr/bin/env bash
# Tier 2 of the clerq.autodiff benchmark plan: two-level, lambda-EIT and Rb87 with 11 Doppler classes.
# Runs the experiments one after another so they do not compete for the CPU. About 4 to 5 hours.
#
#   OLD_TREE=<checkout of the commit before PR 361>  bash benchmarks/autodiff/run_tier2.sh
#
# Keep the machine on AC power and awake (for example: caffeinate -i -s -w $$ &), and close other
# heavy applications; the load average is recorded with every result.
set -u
cd "$(dirname "$0")"
mkdir -p results/logs
OLD_TREE="${OLD_TREE:-}"
step() { echo "[$(date +%H:%M:%S)] $*"; }

step "E4 compile time"
uv run python e4_compile.py --tag e4 > results/logs/e4.log 2>&1

step "E3 accuracy"
uv run python e3_accuracy.py --tag e3 > results/logs/e3.log 2>&1

step "E1 cost of a gradient"
uv run python e1_cost.py --systems two_level lambda_eit --ns 2 5 10 50 100 --tag e1 > results/logs/e1_small.log 2>&1
uv run python e1_cost.py --systems rb87 --ns 2 5 10 50 --tag e1 > results/logs/e1_rb87.log 2>&1

if [ -n "$OLD_TREE" ]; then
  step "E6 before and after PR 361 (decay_rates), back to back"
  PYTHONPATH="$OLD_TREE/src" uv run python e1_cost.py --systems two_level lambda_eit rb87 --ns 10 50 --repeats 5 --tag e6_before > results/logs/e6_before.log 2>&1
  uv run python e1_cost.py --systems two_level lambda_eit rb87 --ns 10 50 --repeats 5 --tag e6_after > results/logs/e6_after.log 2>&1
fi

step "E2 memory"
uv run python e2_memory.py --tag e2 > results/logs/e2.log 2>&1

step "E5 time to solution"
uv run python e5_time_to_solution.py --steps 40 --nm-evals 300 --tag e5 > results/logs/e5.log 2>&1

step "done"
touch results/DONE
