#!/usr/bin/env bash
# Re-execute every notebook and store the outputs in the .ipynb files.
#
# Quarto does not execute .ipynb files on a plain `quarto render`; it shows the
# outputs stored in the notebook. So the executed notebooks are the source of
# truth for the site, and `quarto render --execute` writes the fresh outputs
# back into them. Review the notebook diff and commit it.
#
# Solver results (*.qu) are build artifacts; they are deleted first so every
# problem is solved from scratch. The startup script in .ipython silences
# progress bars and makes plotly div ids and warning paths deterministic, so
# two runs give the same outputs.
#
# Quarto renders one file at a time, so the notebooks are executed here in
# parallel (JOBS at once, default: the number of CPUs), then the site is
# rendered once from the stored outputs. Each notebook runs in its own copy of
# docs/, because concurrent renders in one project race on .quarto/ and _site/.
#
# Usage: ./execute.sh [notebook.ipynb ...]   (default: every notebook in the site)
set -euo pipefail
cd "$(dirname "$0")"

JOBS=${JOBS:-$(getconf _NPROCESSORS_ONLN)}

# The slowest notebooks start first, so they do not run on alone at the end.
SLOW="
examples/mbs-two-stimulated-echo.ipynb
examples/mbs-two-photon-echo.ipynb
usage/autodiff-recipes.ipynb
examples/mbs-lambda-eit-slow-light.ipynb
examples/autodiff-pulse-design.ipynb
usage/numerical-stability.ipynb
"

if [ $# -gt 0 ]; then
    notebooks=$(printf '%s\n' "$@")
else
    # The notebooks in the site, as listed by the render globs in _quarto.yml.
    notebooks=$(quarto inspect | python3 -c '
import json, os, sys
for f in json.load(sys.stdin)["files"]["input"]:
    if f.endswith(".ipynb"):
        print(os.path.relpath(f))
')
fi
ordered=$(
    for nb in $SLOW; do printf '%s\n' "$notebooks" | grep -Fx "$nb" || true; done
    printf '%s\n' "$notebooks" | grep -Fvx -f <(printf '%s\n' $SLOW) || true
)

find examples usage -name '*.qu' -delete

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# One BLAS thread per notebook: with several notebooks at once, threaded BLAS
# oversubscribes the CPUs and the hyperfine notebooks run ~50x slower.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export WORK SRC="$PWD"

run_one() {
    nb=$1
    copy="$WORK/$(printf '%s' "$nb" | tr / _)"
    mkdir "$copy"
    tar -C "$SRC" --exclude=_site --exclude=.quarto -cf - . | tar -C "$copy" -xf -
    start=$(date +%s)
    if (cd "$copy" && IPYTHONDIR="$copy/.ipython" quarto render "$nb" --execute >"$copy.log" 2>&1); then
        cp "$copy/$nb" "$SRC/$nb"
        echo "$(($(date +%s) - start))s  $nb"
    else
        cat "$copy.log" >&2
        echo "FAILED  $nb" >&2
        return 1
    fi
}
export -f run_one

printf '%s\n' "$ordered" | xargs -P "$JOBS" -I{} bash -c 'run_one "$1"' _ {}

quarto render
