# MaxwellBloch — Examples & Research Expansion Master Plan

## Context

Extend the MaxwellBloch example library beyond the current two-level / Λ / V / Ξ /
Rb-D1 set. Three parallel workstreams:

- **Phase A — Textbook examples.** Ship as Sphinx docs, following the existing
  `docs/examples/mbs-<scheme>-<descriptor>.ipynb` pattern (the `.rst` files
  in `docs/examples/` are category TOCs, not individual pages). Each item
  is one PR.
- **Phase B — Research JSON seeds.** Parameter-ready configs plus short research
  notes, living under a new `research/` folder. Not shipped as docs yet — these
  are for Thomas to develop into papers / blog explainers.
- **Phase C — Code extensions.** Larger functionality lumps that unlock some
  Phase B items. Usually a new field / decay / coupling term plus tests.

## Workstream summary

### Phase A — Textbook (13 items, A0 is a prep PR)

| ID  | Title                                                      | Scheme | Effort |
|-----|------------------------------------------------------------|--------|--------|
| A0  | Counter-propagating field support + depletion check        | infra  | S      |
| A1  | Autler–Townes spectroscopy on two-level                    | 2-lvl  | S      |
| A2  | Coherent Population Trapping (CPT) in Λ                    | Λ      | S      |
| A3  | Adiabatons in Λ (Grobe–Hioe–Eberly)                        | Λ      | M      |
| A4  | Two-pulse photon echo (Hahn, optical)                      | 2-lvl  | M      |
| A5  | Stimulated three-pulse photon echo                         | 2-lvl  | M      |
| A6  | Area-theorem odd-π series (1π, 3π, 5π)                     | 2-lvl  | S      |
| A7  | Double-Λ four-wave mixing                                  | 4-lvl  | L      |
| A8  | N-type scheme Kerr (Schmidt–Imamoğlu)                      | 4-lvl  | M      |
| A9  | Rb 87 D2 hyperfine (5s₁/₂ → 5p₃/₂)                         | hfs    | M      |
| A10 | Rydberg-EIT basic (Ξ ladder, counter-propagating)          | Ξ      | M      |
| A11 | Saturation spectroscopy (Lamb dips, crossovers)            | 2-lvl / hfs | M |
| A12 | Magneto-optical rotation (Faraday + NMOR)                  | hfs + B| L      |

### Phase B — Research seeds (7 items)

Tracks: **research** = work on privately, publication-track;
**textbook** = could graduate to Phase A (docs example) once mature;
**blog** = suits an interactive visual explainer on ogden.eu.

| ID  | Title                                                   | Depends on    | Track            |
|-----|---------------------------------------------------------|---------------|------------------|
| B1  | Rydberg EIT-AT electrometer (MW field readout)          | A10           | textbook + blog  |
| B2  | Superheterodyne Rydberg electrometry (LO + signal MW)   | A10, B1       | research         |
| B3  | Thermal Rydberg polariton storage                       | A10, C2       | research         |
| B4  | Λ-scheme simultons (PRL follow-on)                      | existing V    | research (PUB)   |
| B5  | Hot-atom diamond FWM for 420 nm generation              | A7            | research         |
| B6  | Doppler-free three-photon ladder spectroscopy           | A10           | textbook + blog  |
| B7  | Sodium D2 laser guide star pulse dynamics               | A9 pattern    | textbook         |

**(PUB)** flags items where publication is the likely outcome — keep
results out of public docs until preprint is up.

### Phase C — Code extensions (5 items)

| ID  | Title                                                     | Unlocks          |
|-----|-----------------------------------------------------------|------------------|
| C1  | Density-dependent mean-field terms (cooperative Rydberg)  | B1/B2 advanced   |
| C2  | z-dependent detuning field (longitudinal gradient)        | B3 (GEM/CRIB)    |
| C3  | Quantum noise source for spontaneous FWM                  | B5 advanced      |
| C4  | DC B-field (Zeeman shift) support                         | A12, magnetometry|
| C5  | Plotting module (Plotly, publication-grade + videos)      | **all A items**  |

**C5 is infrastructure with broad impact.** It replaces scattered
matplotlib code across example notebooks with a consistent Plotly-based
API. Every new Phase A example benefits. Candidate for landing right
after A0.

## Dependency graph

```
A0 ──► A10 ─┬─► B1 ─► B2
            ├─► B3 ◄─ C2
            └─► B6

A0 ──► A11 (saturation spectroscopy — needs counter-prop pump-probe)

C4 ──► A12 (magneto-optical rotation — needs B-field)

A7 ──► B5 ──► (C3 optional)
A9 ──► B7

C1 ──► B1/B2 (advanced, cooperative regime)

C5 ──► (soft dep on all A items — better plots for each)
```

**A0 is a small prep PR that unblocks the Rydberg line *and* saturation
spectroscopy.** Land it first.

**C5 is a broad-impact infrastructure change that makes every downstream
example better.** Land it early — ideally second, right after A0 — so all
new Phase A notebooks use the consistent Plotly-based API from day one.
Existing matplotlib code keeps working; migration of old notebooks is
opportunistic.

**C4 unblocks A12 and the whole magnetometry line.** Consider whether to
sequence it early if A12 / atomic magnetometry is a strategic priority.

## Conventions

### Git flow (from existing CLAUDE.md)

- Feature branch off `next`, PR to `next`
- One task = one branch = one PR
- `master` only receives `next` at release time
- Don't bump version per task; bump at release

### Pre-PR checklist

Every task, before opening a PR:

```bash
uv run ruff format .
uv run ruff check .
uv run pytest -n auto
uv run sphinx-build docs docs/_build -b html
```

All four must be clean.

### Per-task content checklist

Every new example notebook / feature PR should also:

- [ ] Add a `"savefile": "..."` to any JSON config so the `.qu` cache key
      is stable (existing examples follow this pattern — see
      `mbs-two-sech-2pi.ipynb` etc.)
- [ ] Update `CHANGELOG.md` under a new or existing unreleased section
- [ ] If the change is user-facing, update `CLAUDE.md` so future Claude
      Code sessions know about the new feature / pattern
- [ ] Add regression tests only for non-example changes (Phase C code
      changes); examples get smoke tests at most
- [ ] Make sure `.qu` cache files from local runs aren't committed
      (existing `.gitignore` should handle this; verify)
- [ ] For Phase A example notebooks: include a pulse-area sanity check
      cell (`np.trapz(mbs.Omegas_zt[0,0,:].real, mbs.tlist)/np.pi`) where
      meaningful — matches existing examples

### When the plan is wrong

If Claude Code discovers during execution that the plan says something
incorrect about the codebase (wrong file paths, missing attributes,
different API conventions), it should **fix the sub-plan file in the
same PR** rather than silently working around it. This keeps the plan
accurate for subsequent tasks. A line in the PR description noting the
plan correction is sufficient.

### File naming

Docs pattern (confirmed from the published site):

- **Individual example pages are Jupyter notebooks.** Path:
  `docs/examples/mbs-<scheme>-<descriptor>.ipynb`. These contain narrative
  (Markdown cells), JSON config (Python string), solve call, and plot cells.
- **Usage guides are also notebooks.** Path: `docs/usage/<topic>.ipynb`.
- **TOC / category pages are RST.** Paths: `docs/examples/<category>.rst`
  (e.g. `two-level.rst`, `three-level.rst`), `docs/index.rst`, and the
  usage index.

When adding a new example, modify the relevant category RST to add the
entry, and create the notebook alongside the existing ones. Match the
`mbs-` prefix convention and scheme labels already in use: `two`, `lambda`,
`vee`, `ladder`, `Rb87_5s12_5p12`, etc.

Test fixtures (optional, for smoke tests): `tests/json/mbs-<scheme>-<descriptor>.json`.

Phase B seeds: `research/<Bx>-<scheme>/` containing `problem.json` and `notes.md`.

Phase C: normal source additions under `src/maxwellbloch/`, tests under `tests/`.

### Related: archiving notebooks-maxwellbloch

The separate `notebooks-maxwellbloch` repo is stale — the notebooks that
used to live there are now the docs examples. Archive it (Settings →
Danger Zone → Archive) at your convenience to prevent confusion; put a
pointer in its README to the main repo first. Not blocking for this plan.

### Physics defaults

- Natural units (γ = 1) unless the example is specifically SI-calibrated
- `interaction_strengths` parameterises optical depth; include a comment with the
  equivalent number density and cell length in the narrative
- Velocity classes: `thermal_std = 25.0` (γ-units) is a reasonable Rb D2, 320 K default
- `inner_steps = 2` for Doppler examples unless convergence dictates otherwise
- Always include a Rabi-area or transmission sanity check cell in the notebook

### Sub-plan template

Every `plans/phase-x/Xn-<slug>.md` follows this structure:

```markdown
# Xn: <Title>

## Goal               — one paragraph, what's demonstrated
## Physics            — scheme, key effect, the equation that matters
## Level structure    — ASCII or mermaid diagram
## Parameters         — numerical values with justification
## Files              — create / modify (exact paths)
## Target figures     — ordered list, headline figure first
## Acceptance criteria — reproducible physics checks + CI gates
## Branch & PR        — branch name, PR title, target branch
## References         — papers, order by relevance
## Gotchas            — known traps, scope boundaries
```

Deviations are fine — the template is a floor, not a ceiling.

## Notes on existing code

These are easy things to get wrong on first reading of the codebase.
Claude Code should know them before editing.

### `detuning_positive` is not about propagation direction

Every existing field JSON has a `detuning_positive: true|false` attribute.
This is a **sign convention for the RWA Hamiltonian** based on which level
in `coupled_levels` is higher in energy (see the Λ example in the
three-level usage doc: for a Λ where the second "excited" state sits below
the first, you set `detuning_positive: false`). It's orthogonal to the
`counter_propagating` flag added in A0:

| Flag                     | Flips sign of         | Physical meaning                    |
|--------------------------|-----------------------|-------------------------------------|
| `detuning_positive`      | $\omega_\text{atom}$  | Level-ordering convention (RWA)     |
| `counter_propagating`    | $kv$ (Doppler)        | Laser k-vector direction in lab     |

The two can be set independently. Don't let Claude Code "simplify" them
together during refactoring — they're answering different questions.

### Polarization selection via `q` is already supported

The hyperfine / angular-momentum module already supports polarization
selection through the `q` parameter on field couplings: `q=0` for π,
`q=±1` for σ±. Used in the existing Rb 87 D1 examples. A12 (magneto-
optical rotation) models linearly polarized light as a coherent
superposition of σ+ and σ-, both of which are already expressible — the
missing ingredient is the DC B-field that Zeeman-shifts the sublevels
asymmetrically, which is C4.

### `factors` attribute on decays/fields

Existing structured examples show `"factors": [...]` on decay channels and
field couplings, used for hyperfine branching ratios and angular momentum
matrix elements. See `docs/usage/structure.ipynb` for the pattern.

## Execution order

Rough sequencing that respects dependencies and minimises context-switching:

1. **Prep 1** (blocks the Rydberg and saturation-spectroscopy lines):
   A0 — counter-propagating field support. Small PR, ship first.
2. **Prep 2a** (de-risks the plotting pipeline): C5 **spike PR** — adds
   `plotly + kaleido` to docs extras, tests one Plotly figure through
   Sphinx build. 30-minute job. Verifies which of Path 1/2/3 (see C5) is
   free.
3. **Prep 2b** (unifies every downstream example): C5 **scaffold PR** —
   module structure, theme, `field_spacetime`, `spectrum`, `pulse_area`,
   one migrated example. Uses the path chosen in the spike. Ship before
   A1 starts so new examples land with the consistent API from day one.
   Full-catalogue PR can follow on its own timeline.
4. **Warm-up pass** (parallel, independent): A1, A6 — trivial extensions of
   existing two-level code, using the new plot API. Good for validating
   the workflow end-to-end.
5. **Three-level completeness** (parallel, independent): A2, A3, A9, A10.
   Use existing Λ/Ξ machinery. A10 is the gateway to the Rydberg line and
   depends on A0.
6. **Coherent transients** (sequential): A4 → A5. Echo machinery is shared.
   A11 (saturation spectroscopy) fits here too — it's also a coherent-
   transients-adjacent pump-probe technique, depends on A0.
7. **Nonlinear optics** (sequential): A7 → A8. Double-Λ structure precedes N-type.
8. **Research priorities** (Thomas-led): B4 first (natural extension of the PRL
   simulton work, publication-potential), then B1 (direct extension of A10).
9. **Code extensions on demand**:
   - C2 when B3 becomes the priority
   - C1 when B1/B2 need to demonstrate cooperative physics at high density
   - C4 when A12 / the magnetometry line is prioritised
10. **Magnetometry line**: A12 after C4 lands.
11. **Long tail**: B2, B5, B6, B7, C3 in any order.

## Using this with Claude Code

Recommended loop for each sub-plan:

1. `/clear` — start fresh so context doesn't leak between unrelated tasks
2. Point Claude Code at the repo root and the specific sub-plan file
3. Ask it to execute the plan, stopping before opening the PR
4. Review the generated notebook / JSON / code. The solve and plotting it can do;
   the physical narrative benefits from a human pass
5. Run the four pre-PR checks yourself (fast); fix anything Claude Code missed
6. Let Claude Code draft the PR description, push, and open the PR
7. Merge to `next` after your own review

A good prompt skeleton per task:

> Read `plans/phase-a/A10-rydberg-eit-basic.md`. Execute the plan to completion,
> stopping before `gh pr create`. Keep parameters natural-units where the plan
> specifies; don't switch to SI. Before finishing, run the four pre-PR checks
> from MASTER_PLAN.md and confirm they pass.

## Out of scope (for now)

Things discussed but deliberately deferred:

- Transverse effects (filamentation, self-focusing, OAM transfer) — the code is 1D
- Spin-squeezing / fluctuation-level quantum effects — mean-field only
- Real-time Floquet / modulated-coupling dynamics — would benefit from a separate
  time-dependent-Hamiltonian refactor, out of scope for the example expansion
- GUI / interactive Observable widgets for the blog — separate track from package
  development; use MaxwellBloch outputs as inputs there

## Status tracking

Keep task status in this table. Update when a PR merges to `next`.

| ID  | Status   | PR / Commit | Notes                                                      |
|-----|----------|-------------|------------------------------------------------------------|
| A0  | ✅ Done   | #A0         | Counter-propagating field + depletion check                |
| A1  | ✅ Done   | #260        | Autler-Townes splitting (Ξ ladder)                         |
| A2  | ✅ Done   | #261        | CPT dark resonance (Λ)                                     |
| A3  | ✅ Done   | #263, #266  | Adiabatons in Λ (Grobe–Hioe–Eberly); pulse-centre fix      |
| A4  | ✅ Done   | b5ef65d     | Two-pulse photon echo; branch not yet PRed to `next`       |
| A5  | Planned   |             | Stimulated three-pulse photon echo (depends on A4)         |
| A6  | ✅ Done   | #259        | Odd-π sech pulse area theorem (1π, 3π, 5π)                 |
| A7  | Planned   |             | Double-Λ four-wave mixing                                  |
| A8  | Planned   |             | N-type scheme Kerr                                         |
| A9  | ✅ Done   | #264        | Rb 87 D2 hyperfine (5s₁/₂ → 5p₃/₂, σ⁺ optical pumping)   |
| A10 | ✅ Done   | #262        | Rydberg EIT (Ξ ladder, counter-propagating Doppler narrow) |
| A11 | Planned   |             | Saturation spectroscopy (depends on A0)                    |
| A12 | Planned   |             | Magneto-optical rotation (depends on C4)                   |
| C5  | ✅ Done   | #258        | Plotly plotting module (field_spacetime, coherence, etc.)  |

### Phase A implementation notes

**A4 (photon echo)** — two important codebase conventions discovered during development:
- Lindblad collapse operator includes `√(2π)` factor: actual population decay rate
  `Γ = 2π × rate`, so `T₂ = 2/(2π × rate)`. Use `rate = 0.1` for `T₂ ≈ 3.18 γ⁻¹`.
- Hamiltonian detuning is `4π² × thermal_delta` (user units), so `thermal_width = 0.15`
  maps to `σ ≈ 4.2 rad/γ`. Hard-pulse condition requires `Ω_peak >> 4π² × thermal_width`.

(Fill in PRs as they open; can live in a GitHub Project instead if preferred.)
