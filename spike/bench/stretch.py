"""Gate 1: the stretch point (j1_2_res dim 9, n_z 500, n_t 2000, n_v 50, max_steps 4096)."""
import json, subprocess, sys
base = ["--system", "j1_2_res", "--n_z", "500", "--n_t", "2000", "--n_v", "50", "--max_steps", "4096"]
mode = sys.argv[1]
cells = [("s1", []), ("s2", []), ("s3", []), ("s4", ["--checkpoints", "32"]), ("s5", [])]
if mode == "run":
    cells = [("s4", ["--checkpoints", "32"]), ("s3", []), ("s1", []), ("s2", [])]
out = f"bench/results/stretch_{mode}.jsonl"
with open(out, "w") as f:
    for st, extra in cells:
        flags = ["--compile_only"] if mode == "compile" else ["--no_warmup", "--repeat", "1"]
        r = subprocess.run([sys.executable, "bench/phase4.py", *base, "--strategy", st, *flags, *extra],
                           capture_output=True, text=True)
        line = [l for l in r.stdout.splitlines() if l.startswith("RESULT ")]
        d = json.loads(line[0][7:]) if line else {"strategy": st, "error": r.stderr[-400:]}
        f.write(json.dumps(d) + "\n"); f.flush()
        print(st, {k: (round(v, 1) if isinstance(v, float) else v) for k, v in d.items()
                   if k in ("temp_MiB", "run_s", "peak_rss_MiB", "grad_norm", "compile_s", "error")}, flush=True)
