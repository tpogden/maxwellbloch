"""Turn the JSONL results into markdown tables (results/tables.md) and figures.

uv run python benchmarks/autodiff/analyse.py [--figures-dir docs/usage/figures] [--prefix smoke_]
"""

import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
RESULTS = HERE / "results"


def load(prefix: str = ""):
    recs = {}
    for f in sorted(RESULTS.glob(f"{prefix}*.jsonl")):
        if not prefix and f.name.startswith("smoke_"):
            continue
        recs[f.stem[len(prefix) :]] = [
            json.loads(line) for line in f.read_text().splitlines() if line.strip()
        ]
    return recs


def fmt(rec, unit=1.0, digits=3):
    if rec is None or "median" not in rec:
        return "–"
    s = f"{rec['median'] * unit:.{digits}g}"
    if "ci_lo" in rec:
        s += f" ({rec['ci_lo'] * unit:.{digits}g}–{rec['ci_hi'] * unit:.{digits}g})"
    return s


def e1_tables(recs):
    out, by = [], {}
    for r in recs:
        by.setdefault(r["system"], {}).setdefault(r["N"], {})[r["method"]] = r
    for system, cells in by.items():
        solve = cells.get(0, {}).get("clerq_solve")
        any_cell = next(
            (m for n, d in cells.items() if n for m in d.values() if "dim" in m), {}
        )
        out.append(
            f"**{system}** (dim {any_cell.get('dim', '?')}, {any_cell.get('n_classes', '?')} velocity class(es), "
            f"n_t = {any_cell.get('n_t', '?')}, n_z = {any_cell.get('n_z', '?')}); one `solve()` takes {fmt(solve)} s.\n"
        )
        out.append(
            "| N | forward (s) | reverse (s) | reverse / forward | forward-mode (s) | FD through JAX (s) | FD through clerq (s, estimate) | reverse speed-up over FD (JAX) | over FD (clerq) |"
        )
        out.append("|-|-|-|-|-|-|-|-|-|")
        for n in sorted(k for k in cells if k):
            c = cells[n]
            f, rv, fa, fj = (
                c.get("forward"),
                c.get("reverse"),
                c.get("forward_ad"),
                c.get("fd_jax"),
            )
            skip = c.get("forward_ad", {}).get("skipped")
            ratio = rv["median"] / f["median"] if f and rv else None
            est = (n + 1) * solve["median"] if solve else None
            out.append(
                f"| {n} | {fmt(f)} | {fmt(rv)} | {ratio:.1f} | {('skipped' if skip else fmt(fa))} | {fmt(fj)} | "
                f"{est:.3g} | {fj['median'] / rv['median']:.1f}x | {est / rv['median']:.0f}x |"
                if ratio and fj and est
                else f"| {n} | incomplete |"
            )
        out.append("")
    return "\n".join(out)


def e2_table(recs):
    rows = [r for r in recs if "rss_peak_mib" in r]
    out = [
        "| system | sweep | n_z | classes | strategy | peak resident set (MiB) | above baseline (MiB) | XLA temporary (MiB) | steady gradient (s) |",
        "|-|-|-|-|-|-|-|-|-|",
    ]
    for r in rows:
        out.append(
            f"| {r['system']} | {r['sweep']} | {r['n_z']} | {r['n_classes']} | {r['strategy']} | {r['rss_peak_mib']:.0f} | "
            f"{r['rss_peak_mib'] - r['rss_before_mib']:.0f} | {r.get('xla_temp_mib') or 0:.1f} | {r['steady_s']:.2f} |"
        )
    for r in recs:
        if "error" in r or "skipped" in r:
            out.append(
                f"| {r['system']} | {r['sweep']} | {r['n_z']} | – | {r['strategy']} | {r.get('error', r.get('skipped'))} | | | |"
            )
    return "\n".join(out)


def e3_table(recs):
    out = [
        "| system | N | fixed steps | field error vs `solve()` | AD vs finite differences (best step) | forward vs reverse | adaptive vs fixed-step AD |",
        "|-|-|-|-|-|-|-|",
    ]
    for r in recs:
        out.append(
            f"| {r['system']} | {r['N']} | {r['n_steps_fixed']} | {r['field_error_vs_solve']:.1e} | {r['ad_vs_fd_fixed']:.1e} | {r['fwd_vs_rev']:.1e} | {r['adaptive_vs_fixed']:.1e} |"
        )
    return "\n".join(out)


def e4_table(recs):
    out = [
        "| system | N | forward: first call (s) | steady (s) | reverse gradient: first call (s) | steady (s) |",
        "|-|-|-|-|-|-|",
    ]
    for r in recs:
        out.append(
            f"| {r['system']} | {r['N']} | {r['forward_first_s']:.2f} | {r['forward_steady_s']:.3f} | {r['reverse_first_s']:.2f} | {r['reverse_steady_s']:.3f} |"
        )
    return "\n".join(out)


def e6_table(before, after):
    key = lambda r: (r["system"], r["N"], r["method"])  # noqa: E731
    b = {key(r): r for r in before if "median" in r and r["N"]}
    out = [
        "| system | N | method | before PR 361 (s) | after (s) | after / before |",
        "|-|-|-|-|-|-|",
    ]
    for r in after:
        if (
            "median" in r
            and r["N"]
            and key(r) in b
            and r["method"] in ("forward", "reverse", "forward_ad")
        ):
            k = key(r)
            out.append(
                f"| {k[0]} | {k[1]} | {k[2]} | {fmt(b[k])} | {fmt(r)} | {r['median'] / b[k]['median']:.2f} |"
            )
    return "\n".join(out)


def e5_text(recs):
    by = {r["what"]: r for r in recs}
    if not by:
        return ""
    g, nm, ad = (
        by.get("gradient_cost"),
        by.get("baseline_nelder_mead"),
        by.get("free_form_adam"),
    )
    lines = []
    if g:
        lines.append(
            f"- {g['n_par']}-parameter gradient: reverse AD {fmt(g['reverse_ad'])} s, forward differences {fmt(g['fd_forward_differences'])} s, one forward evaluation {fmt(g['forward'])} s."
        )
    if nm:
        lines.append(
            f"- Baseline (5 parameters, Nelder-Mead): {nm['retrieved']:.4f} retrieved after {nm['evals']} evaluations, {nm['wall_s']:.0f} s."
        )
    if ad:
        s = f"- Free-form (44 parameters, Adam): {ad['retrieved_final']:.4f} after {ad['steps']} steps, {ad['wall_s']:.0f} s."
        if ad.get("step_passing_target") is not None:
            s += (
                f" It passes the baseline value ({ad['target']:.4f}) at step {ad['step_passing_target'] + 1}, after {ad['wall_s_to_target']:.0f} s; "
                f"with finite-difference gradients the same steps would take about {ad['fd_gradient_wall_s_to_target_estimate']:.0f} s (estimate from the measured gradient time)."
            )
        lines.append(s)
    return "\n".join(lines)


def e7_table(recs):
    """Median over seeds of the best retrieved fraction at the full allotment, per optimiser."""
    by = {}
    for r in recs:
        by.setdefault(r["optimiser"], []).append(r)
    if not by:
        return ""
    first = next(iter(by.values()))[0]
    out = [
        f"{first['n_par']} parameters, budget {first['budget']:.0f} objective evaluations "
        f"(a gradient call is charged {first['charge_per_gradient']:.1f}), optical depth {2 * first['g']:.0f}"
        f"{', Doppler' if first['doppler'] else ', no Doppler'}.\n",
        "| optimiser | per seed | median | cost used |",
        "|-|-|-|-|",
    ]
    for name in ("adam", "lbfgsb", "de", "cmaes", "lhs"):
        rs = sorted(by.get(name, []), key=lambda r: r["seed"])
        if rs:
            vals = [r["best_retrieved"] for r in rs]
            out.append(
                f"| {name} | {', '.join(f'{v:.3f}' for v in vals)} | {np.median(vals):.3f} | "
                f"{np.median([r['cost_used'] for r in rs]):.0f} |"
            )
    return "\n".join(out)


def e7_figure(recs, outdir: Path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not recs:
        return
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    labels = {
        "adam": "Adam (gradient)",
        "lbfgsb": "L-BFGS-B (gradient)",
        "de": "differential evolution",
        "cmaes": "CMA-ES",
        "lhs": "Latin hypercube",
    }
    for name, label in labels.items():
        rs = [r for r in recs if r["optimiser"] == name]
        if not rs:
            continue
        grid = np.linspace(0, rs[0]["budget"], 200)
        curves = []
        for r in rs:
            c, v = np.array(r["trace"]).T
            curves.append(np.interp(grid, c, v, left=0.0))
        ax.plot(grid, np.median(curves, axis=0), label=label)
        ax.fill_between(
            grid, np.min(curves, axis=0), np.max(curves, axis=0), alpha=0.15
        )
    ax.set(
        xlabel="cost (objective evaluations)",
        ylabel="best retrieved fraction so far",
        ylim=(0, None),
    )
    ax.legend(fontsize=7)
    fig.tight_layout()
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / "autodiff-bench-equal-cost.png", dpi=150)
    plt.close(fig)


def figures(e1, e2, outdir: Path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir.mkdir(parents=True, exist_ok=True)
    systems = sorted({r["system"] for r in e1})
    if systems:
        fig, axes = plt.subplots(
            1,
            len(systems),
            figsize=(4.2 * len(systems), 3.6),
            sharey=False,
            squeeze=False,
        )
        for ax, system in zip(axes[0], systems, strict=True):
            rs = [r for r in e1 if r["system"] == system]
            solve = next((r for r in rs if r["method"] == "clerq_solve"), None)
            for method, label, style in (
                ("reverse", "reverse-mode AD", "o-"),
                ("forward_ad", "forward-mode AD", "s--"),
                ("fd_jax", "finite differences (JAX)", "^:"),
            ):
                pts = sorted(
                    (r["N"], r) for r in rs if r["method"] == method and "median" in r
                )
                if pts:
                    ns = [n for n, _ in pts]
                    ax.plot(ns, [r["median"] for _, r in pts], style, label=label)
                    lo = [r.get("ci_lo", r["min"]) for _, r in pts]
                    hi = [r.get("ci_hi", r["max"]) for _, r in pts]
                    ax.fill_between(ns, lo, hi, alpha=0.2)
            if solve:
                ns = sorted({r["N"] for r in rs if r["N"]})
                ax.plot(
                    ns,
                    [(n + 1) * solve["median"] for n in ns],
                    "x-.",
                    color="grey",
                    label="finite differences (clerq), estimate",
                )
            ax.set(xscale="log", yscale="log", xlabel="parameters N", title=system)
        axes[0][0].set_ylabel("time per gradient (s)")
        axes[0][0].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(outdir / "autodiff-bench-cost.png", dpi=150)
        plt.close(fig)
    rows = [r for r in e2 if "rss_peak_mib" in r and r["sweep"] == "n_z"]
    if rows:
        systems = sorted({r["system"] for r in rows})
        fig, axes = plt.subplots(
            1, len(systems), figsize=(4.4 * len(systems), 3.6), squeeze=False
        )
        for ax, system in zip(axes[0], systems, strict=True):
            for strat in ("s2", "s3", "s4"):
                pts = sorted(
                    (r["n_z"], r["xla_temp_mib"])
                    for r in rows
                    if r["system"] == system
                    and r["strategy"] == strat
                    and r.get("xla_temp_mib")
                )
                if pts:
                    ax.plot(
                        *zip(*pts, strict=True),
                        "o-",
                        label={
                            "s2": "s2: store every slice",
                            "s3": "s3: checkpoint the scan body (default)",
                            "s4": "s4: s3 + 32 inner checkpoints",
                        }[strat],
                    )
            ax.set(xscale="log", yscale="log", xlabel="z steps", title=system)
        axes[0][0].set_ylabel("XLA temporary memory of one gradient (MiB)")
        axes[0][0].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(outdir / "autodiff-bench-memory.png", dpi=150)
        plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--figures-dir", default=str(HERE.parent.parent / "docs" / "usage" / "figures")
    )
    p.add_argument("--no-figures", action="store_true")
    p.add_argument("--prefix", default="")
    a = p.parse_args()
    d = load(a.prefix)
    parts = []
    if "e1" in d:
        parts += ["## E1 cost of a gradient\n", e1_tables(d["e1"])]
    if "e6_before" in d and "e6_after" in d:
        parts += [
            "## E6 before and after PR 361\n",
            e6_table(d["e6_before"], d["e6_after"]),
        ]
    if "e2" in d:
        parts += ["## E2 memory\n", e2_table(d["e2"])]
    if "e3" in d:
        parts += ["## E3 accuracy\n", e3_table(d["e3"])]
    if "e4" in d:
        parts += ["## E4 compile time\n", e4_table(d["e4"])]
    if "e5" in d:
        parts += ["## E5 time to solution\n", e5_text(d["e5"])]
    if "e7" in d:
        parts += [
            "## E7 equal-cost comparison with gradient-free search\n",
            e7_table(d["e7"]),
        ]
    text = "\n\n".join(parts) + "\n"
    (RESULTS / f"{a.prefix}tables.md").write_text(text)
    print(text)
    if not a.no_figures:
        figures(d.get("e1", []), d.get("e2", []), Path(a.figures_dir))
        e7_figure(d.get("e7", []), Path(a.figures_dir))


if __name__ == "__main__":
    main()
