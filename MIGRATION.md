# clerq migration: progress

Tracks the `maxwellbloch` → `clerq` rename and the Sphinx → Quarto docs move.
Source plan: `clerq-migration-plan.md`. Tick boxes as work lands. Delete this
file at cutover (Phase 4), or move it to `docs/news/`.

**Branch:** `claude/clerq-spike-migration-plan-efba93` (long-lived; PRs target `next`)
**Last updated:** 2026-10-09 (cutover complete)

Legend: `[x]` done · `[ ]` to do · 👤 you must do it (external or credentialed) · 🤖 Claude can do it

---

## Phase 0: Pre-work

- [x] 0.1 Reserve `clerq` on PyPI 👤: `clerq 0.0.0` live, verified 2026-10-09
- [x] 0.2 Decisions locked (see plan): continue `0.x.y`, rename in place, ~~`freeze: auto`~~ (replaced by outputs stored in the notebooks), shim, new Zenodo DOI
- [x] 0.3 `clerq.org` DNS points at GitHub Pages 👤: apex A records `185.199.108–111.153` and `www` CNAME to `tpogden.github.io`, verified 2026-10-09 (MX left as Hover mail forwarding)
- [ ] 0.4 Revoke the account-wide PyPI token used for the stub; create a `clerq`-scoped token for Phase 4 👤
- [ ] 0.5 Decide whether to batch a small API cleanup into the rename release (open question in plan)

## Phase 1: Quarto build (local)

- [x] `docs/_quarto.yml`, `index.qmd`, `installation.qmd`, `examples/index.qmd`
- [x] YAML frontmatter cell on 12 example + 8 usage notebooks (`ef80153`)
- [x] `quarto render` clean: 25 pages, no broken internal links, Plotly renders
- [x] `_site/` gitignored

**Phase 1 complete** apart from the README badge, which moves to Phase 3.
- [x] Draft `mbs-ladder-rydberg-eit-counter.ipynb` stays excluded from the build (same as under Sphinx); finishing it is separate work
- [x] API reference: quartodoc (28 generated pages, committed; CI checks they are up to date)
- [x] Sphinx-only content ported (figure caption, gallery); no `:doc:`/`:ref:` left in notebooks; navbar uses one sidebar per section
- [x] Sphinx removed: `conf.py`, `*.rst`, `make.bat`, `.readthedocs.yml`, `custom.css`; `docs` extra, Makefile targets and CI `docs` job now use Quarto
- [ ] readthedocs badge in README (handled in Phase 3 README rebrand)

Notes: kept the existing `mbs-*` notebook names. The plan's `a0`–`a12` naming does
not exist in this repo.

## Phase 2: Notebook reproducibility audit ✅

Done 2026-10-09 (`95a0e96`). Factor-of-2 fix (PR #289) was merged into `next` first and is in this branch.

- [x] 2.1 Inventory: no random numbers or seeds anywhere; nondeterminism came from `%time` (3 notebooks),
      tqdm progress rates, plotly random div ids, and random nbformat cell ids
- [x] 2.2 Fixed: `%time` removed; `docs/.ipython/.../00-docs-determinism.py` silences tqdm and makes `uuid4`
      deterministic (used only by `docs/execute.sh`); frontmatter cells given a fixed `id`
- [x] 2.3 Two from-scratch executions (no `.qu`) give byte-identical results (7 min each)
- [x] 2.4 Committed `.qu` caches deleted and `*.qu` ignored everywhere

**Correction (2026-10-09).** Phase 2 originally stored the executed outputs in `docs/_freeze/` and claimed the site
was built from it. That was wrong. Quarto never executes a `.ipynb` on a plain `quarto render`
(`execute.enabled` is false for notebooks), so it shows the outputs *stored in the notebook* and ignores `_freeze/`.
The notebooks committed in `95a0e96` did contain outputs, but from early runs made without the determinism script,
so the deployed pages showed tqdm progress-bar spam (1201 cells in 7 notebooks; 124 lines on the `sech-2pi` page) and
the "clean clone renders from `_freeze`" gate had only tested those stored outputs. Fixed in the PR that removed
`_freeze/`: the executed outputs now live in the notebooks, two from-scratch runs give byte-identical notebooks
(0 of 21 differ), and the new outputs match the audited freeze (every stdout line, same figure counts).

**Examples vs the factor-of-2 fix (2026-10-09).** The examples were tuned before PR #289 doubled the effective optical
depth, so after it they showed different physics (e.g. the 0.5π + 1.5π simulton figure overshot to 3.3π and was absorbed
in half the distance). Audit: for each notebook, the problem was solved with the pre-fix solver (`maxwellbloch 0.12.0`,
original strengths) and with `clerq` at half the strengths; the field arrays are bit-identical for all 19 notebooks that
call `MBSolve`. So every `interaction_strengths` was halved (SIT's `interaction_strength` too), which restores each
example exactly. Exceptions where the number is part of the text: `linear-absorption` was relabelled (Ng = 0.05 / 0.5 / 5,
amplitude `exp(-2Ng z/Γ)`, `alpha(0) = 4Ng/Γ`), the adiabatons OD table now says `x 4`, and `velocity-classes` scales
`spectral.voigt_two_linear_known` (a unit-strength profile) by the interaction strength. Plot limits that clipped the data
(`vee-simultons`) were widened. Any future solver change that alters results should be audited the same way.

How to refresh after changing code or a notebook: `make docs_execute` (runs `docs/execute.sh`, ~8 min; deletes `.qu`
caches, then `quarto render --execute` with the determinism startup script), review the notebook diff, commit. Do not
run `quarto render --execute` without the startup script (`IPYTHONDIR=docs/.ipython`): outputs then include progress
bars and random plotly ids.

Things learned:

- 15 of 20 notebooks had no stored outputs (Sphinx executed them at build time), so the Phase 1 site had empty
  pages until now.
- Notebooks used `recalc=False` against committed `.qu` caches, so they never re-solved. Those caches are gone.
- The `mbs-ladder-rydberg-eit-counter` draft is excluded from the build and from re-execution.

## Phase 3: Code rename ✅

Done 2026-10-09 (`61d85d4` rename, `aa6a2b0` regenerated docs). Version stays `0.12.0`.

- [x] `git mv src/maxwellbloch src/clerq`; imports in `src/`, `tests/`, `bin/` and notebook code cells
- [x] `pyproject.toml` name, URL, `dev` extra self-reference; `uv.lock` re-locked
- [x] Plotly template `"maxwellbloch"` → `"clerq"`; `.qu` metadata key now `clerq` (old savefiles still load)
- [x] `.gitignore`, `bin/` scripts, `CLAUDE.md`, Makefile
- [x] README: rename notice, install, `clerq.org`, citation (original MaxwellBloch citation kept; clerq DOI pending); coveralls badge dropped (not used in CI)
- [x] `CHANGELOG.md`: `[Unreleased]` entry for the rename and docs move
- [x] CI: no `maxwellbloch` strings in the workflows (docs job already switched to Quarto in Phase 1)
- [x] Gate: 208 tests pass; `ruff check` and `ruff format --check` clean
- [x] Gate: wheel builds as `clerq-0.12.0`, installs in a clean venv, `import clerq` works, `import maxwellbloch` fails, `mbsolve`/`obsolve`/movie scripts installed
- [x] Gate: re-executed outputs are identical to the pre-rename ones once the name is normalised, so the rename changes no numerics
- [x] Gate: clean-clone `quarto render` 24 s, 0 broken links, all pages have outputs

Changed on purpose vs the plan: the site description in `_quarto.yml` says "Maxwell–Bloch solver", not
"Differentiable Maxwell–Lindblad solver", because nothing is differentiable until the spike passes. The plan's
pyproject description was left as the original.

Deferred to Phase 4/5 (needs the repo rename, PyPI or you):

- PyPI trusted publishing (`publish.yml`, environment `pypi`) is registered for project `MaxwellBloch`;
  register `clerq` as a trusted publisher (or use a project token) before tagging 👤
- GitHub URLs in README, `pyproject.toml` and `_quarto.yml` point at `tpogden/clerq`, which only exists after
  the repo rename (Phase 4.2)
- ~~Zenodo DOI and clerq citation~~ done (concept DOI `10.5281/zenodo.23269881`)
- Existing bug found in passing, unrelated to the rename: `mbsolve --help` crashes under Python 3.14 (unescaped `%`
  in an argparse help string)

## Phase 4: Cutover (single session)

### Prepared (🤖, done 2026-10-09)

- [x] `maxwellbloch` shim built and tested, **not uploaded**: `.../scratchpad/maxwellbloch-shim` (`dist/maxwellbloch-0.12.1*`,
      `twine check` passes). Unlike the plan's `from clerq import *`, it aliases submodules, so `from maxwellbloch import mb_solve`,
      `from maxwellbloch.field import Field` and `maxwellbloch.plot` return the same objects as `clerq`; it warns with
      `DeprecationWarning`; and it keeps the `[plot]`, `[test]`, `[docs]`, `[dev]` extras. It requires `clerq>=0.13.0`.
- [x] `.github/workflows/docs.yml`: publishes the Quarto site to `gh-pages` on push to `master` (and manually). Verified that a
      clean clone renders with no Jupyter installed, from the outputs stored in the notebooks
- [x] `docs/CNAME` (`clerq.org`), included in the site via `resources`
- [x] Existing `publish.yml` (tag `v*` → tests → build → PyPI via OIDC → GitHub release) needs no code change

### Version decision (agreed 2026-10-09: clerq 0.13.0, shim 0.12.1)

Last `maxwellbloch` on PyPI is **0.12.0**. The first clerq release contains the factor-of-2 physics fix (a breaking change to
results), so I recommend **`clerq 0.13.0`** via `uv run bump-my-version bump minor`, and the shim as `maxwellbloch 0.12.1`
(above 0.12.0, so pip picks it, and pinned `==0.12.0` users are untouched). The stub `clerq 0.0.0` is already on PyPI.

### Order of operations (each step depends on the one before)

Pre-flight 👤

- [ ] This branch's PR (to `next`) is green on CI
- [x] DNS for `clerq.org` verified (2026-10-09)
- [x] Tagged the last pre-rename release as `final-maxwellbloch` (on the `v0.12.0` commit `baa79e2`), 2026-10-09. Not named `v-final-…`: `publish.yml` runs on any `v*` tag and would have tried to republish to PyPI

Steps

1. [x] 👤 **Rename the GitHub repo** `maxwellbloch` → `clerq` (Settings → General). Then update local remotes:
       `git remote set-url origin git@github.com:tpogden/clerq.git`
2. [x] 👤 **Register PyPI trusted publishers** (PyPI → project `clerq` → Publishing): owner `tpogden`, repo `clerq`,
       workflow `publish.yml`, environment `pypi`. The old one is registered for `MaxwellBloch` and will not work for `clerq`.
       Check the GitHub environment `pypi` still exists after the rename
3. [x] 🤖/👤 PR #301 merged into `next` (merge commit `af73bc2`)
4. [x] 🤖 `next` bumped to 0.13.0 (`3347b6e`, tag `v0.13.0`). `uv.lock` records the project version and `bump-my-version` does not update it, so: `bump-my-version bump minor --no-commit --no-tag`, `uv lock`, then commit and tag by hand
5. [x] 👤 Release merge (done 2026-10-09, `master` at `efaf96d`; the `pypi` environment has a required reviewer, so approve the run in the Actions UI): `git checkout master && git merge --no-ff next && git push && git push --tags`; the tag runs
       `publish.yml`, which uploads `clerq 0.13.0` to PyPI and makes the GitHub release. The push to `master` also runs `docs.yml`
6. [x] 👤 Pages (HTTPS enforced). **First run needs a `gh-pages` branch to exist**: `quarto publish` cannot create it, so the first `docs.yml` run failed until an empty orphan `gh-pages` branch was pushed and the workflow re-run. GitHub then configured Pages itself from the branch and `CNAME`: Settings → Pages → source branch `gh-pages`, custom domain `clerq.org`, enforce HTTPS (certificate can take ~1 h)
7. [x] 👤 Shim `maxwellbloch 0.12.1` uploaded. Upload the shim (after `clerq 0.13.0` is live, since it depends on it):
       `cd <scratchpad>/maxwellbloch-shim && uvx twine upload dist/*` (use a `maxwellbloch`-scoped token)
8. [x] 🤖 Verification checklist below

Verification

- [x] `https://clerq.org` loads; one example notebook shows plots
- [x] `https://github.com/tpogden/maxwellbloch` redirects to `clerq`
- [x] Fresh venv: `pip install clerq` and `import clerq` work; `clerq.__version__ == "0.13.0"`
- [x] Fresh venv: `pip install maxwellbloch` installs the shim and clerq; `from maxwellbloch import mb_solve` works and warns
- [x] `pip install --no-cache-dir` if PyPI metadata lags (it did, briefly: the `/pypi/<name>/json` summary endpoint showed 0.12.0 for a few minutes after the shim upload while the simple index already had 0.12.1)

Verified 2026-10-09: `clerq 0.13.0` and `maxwellbloch 0.12.1` on PyPI; fresh-venv installs of both; `maxwellbloch==0.12.0` pin
still installs the old release alone; shim solve output is identical to clerq; `https://clerq.org` (and `www`, and http) serve
the site over HTTPS; GitHub URLs redirect; release `v0.13.0` published; CI green on `master`.

### Rollback notes

- PyPI releases cannot be re-uploaded or deleted usefully; if `clerq 0.13.0` is broken, yank it and release `0.13.1`
- Do the shim upload last so a bad clerq release does not leave `maxwellbloch` users pointed at it
- The repo rename is reversible (rename back), but do not create a new `maxwellbloch` repo; that would break the redirect

## Phase 5: Post-cutover

- [x] Rename announcement: `docs/news/rename.qmd`, banner on the site home page and README (includes the factor-of-2 results change)
- [x] Zenodo: `v0.13.1` archived 2026-10-09. Concept DOI (all versions, used in the README and `CITATION.cff`) `10.5281/zenodo.23269881`; version DOI for v0.13.1 `10.5281/zenodo.23269882`. The old MaxwellBloch citation stays for earlier versions 👤
- [x] Archived `notebooks-maxwellbloch` 2026-10-09 with a README pointer to `clerq.org`, repo description and homepage updated, issue #32 closed 👤
- [x] Releases, both 2026-10-09 and verified from fresh venvs:
  - `clerq 0.13.0`: the rename, the factor-of-2 fix (breaking: results change), Quarto docs, `maxwellbloch 0.12.1` shim
  - `clerq 0.13.1`: `mbsolve` and `obsolve` scripts fixed (both crashed in 0.13.0), examples recalibrated for the
    factor-of-2 fix with clean outputs, new description, rename announcement, `CITATION.cff`, stimulated-echo example
- [x] Release checklist for next time (each item cost a detour once):
  1. Everything merged to `next`, CI green. No `Unreleased` leftovers in `CHANGELOG.md`.
  2. `bump-my-version bump <part> --no-commit --no-tag` on a clean tree (or add `--allow-dirty` after editing the
     changelog), then `uv lock` (the bump does not touch the project version in `uv.lock`), then commit and `git tag -a vX.Y.Z`.
  3. Re-run `ruff`, `pytest -n auto`, `uv build` + `twine check`, and a docs render on the final tree. If `next` moves
     after tagging, redo the bump on the new tip and re-tag.
  4. Tag names starting with `v` start `publish.yml`: never use `v-...` for markers.
  5. Optional: enable Zenodo for the repo first (it only archives releases published afterwards).
  6. `git push origin next`, then on `master`: `git merge --no-ff next && git push && git push origin vX.Y.Z`.
  7. Approve the `pypi` environment gate (required reviewer) on the Publish run.
  8. Verify: fresh venv install, the CLI scripts, the shim, clerq.org, the GitHub release. PyPI's `/pypi/<name>/json`
     summary endpoint lags by minutes; trust `/simple/` and `/pypi/<name>/<version>/json`.
- [ ] External references sweep 👤 (done 2026-10-09 except the items marked open):
  - [x] `tpogden/clerq` repo description and homepage (`https://clerq.org`)
  - [x] Site `index.qmd` Projects entry in `tpogden.github.com` (master, `fb70693`) and `tpogden-quarto` (main, `ce6796e`).
        `ogden.eu` is served from the `gh-pages` branch of `tpogden.github.com`, published by hand, so the live site
        changes when you next run `quarto publish gh-pages`
  - [x] `tpogden.github.com/environment.yml`: pip `clerq`, python 3.11, dropped the `qutip=4` / python 3.7 pins (no
        `maxwellbloch` >= 0.8 could install with them). The file is unreferenced; the site is built with uv
  - [x] Read the Docs, `maxwellbloch.readthedocs.io`: exact redirect with *Force redirect* on, verified 2026-10-09:
        `/`, `/en/latest/`, `/en/stable/`, deep pages (`/en/latest/examples/...`), the API index and old tagged versions
        (`/en/v0.8.1/...`) all return 301 to `https://clerq.org/` (deep links land on the new home page, since page names
        differ). `.readthedocs.yml` is gone, so new RTD builds fail, which no longer matters
  - [ ] CV (`tpogden/cv`, `tpo-cv.tex`): project entry still says MaxwellBloch
  - [ ] `paper-maxwellbloch`: title and "Program Title" say MaxwellBloch and the abstract frames the package as solving
        the Maxwell-Bloch equations; a rewrite decision, not a find and replace
  - [ ] Blog posts that `import maxwellbloch` and use parameters tuned before the factor-of-2 fix (figures are baked, but a
        re-execution gives different optical depths); add a note or leave
  - [ ] Google Scholar and ORCID entries (not reachable from the repo)
  - Left alone on purpose: `phd-thesis`, `research-notes` (historical), the CoOMBE paper's link to
    `github.com/tpogden/maxwellbloch` (redirects), `notebooks-maxwellbloch` (archived)

---

## Follow-ups after the cutover

- [x] PR #278 (A5 stimulated echo): recalibrated, merged 2026-10-09 and released in 0.13.1
- [ ] Dependabot: 19 alerts (1 critical), all transitive docs/dev packages in `uv.lock` (see the lock-upgrade PR). Dependabot PRs target `master`; consider `target-branch: next` in `.github/dependabot.yml`
- [x] `mbsolve` and `obsolve` scripts crashed in 0.13.0 (not only `--help`): fixed in #312, with end-to-end tests, shipped in 0.13.1
- [x] Quarto API page titles came out lowercased in the browser tab (`mb_solve.mbsolve`): fixed with `docs/_filters/pagetitle.lua`, which sets `pagetitle` from the first heading for pages without a title
- [x] Release workflow attached `default.gitignore` (uv writes `dist/.gitignore`): `publish.yml` now lists `dist/*.whl`, `dist/*.tar.gz` and `dist/*.publish.attestation`. Takes effect, and can only be confirmed, at the next release
- [x] quartodoc docstring warnings (16): fixed in `hyperfine.py`, `ob_atom.py`, `spectral.py` (docstrings and annotations only); `quartodoc build` is now silent
- [x] CHANGELOG entries for 0.10.0, 0.11.0 and 0.12.0 backfilled from the git history and API diffs; the 0.13.1 entry now says `obsolve` was broken since 0.10.0 and `mbsolve` since 0.11.0
- [x] CI `bench` gate was flaky (50% tolerance vs measured run-to-run variance of up to +127% on identical code); widened to 200% in #306

## Open issues and risks

- `next` and `master` show as diverged by the release merge commits. That is expected: each release merges `next` into `master`.
- The `.qu` cache files are committed under `docs/` (`!docs/examples/*.qu`). Resolved: `.qu` caches are no longer committed.
- The spike's counter-propagation exclusion vs the `usage/counter-propagating` notebook: the notebook documents existing behaviour, so it migrates as is.
