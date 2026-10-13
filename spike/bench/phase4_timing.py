"""Executed timing cells at the base config (j1_2_res, n_z 200, n_t 1000, n_v 10, max_steps 4096)."""
import json, subprocess, sys
base = ["--system", "j1_2_res", "--n_z", "200", "--n_t", "1000", "--n_v", "10",
        "--max_steps", "4096", "--repeat", "2"]
cells = [("s1", []), ("s2", []), ("s3", []), ("s4", ["--checkpoints", "32"]), ("s5", []), ("fwd", [])]
with open("bench/results/timing_base.jsonl", "w") as f:
    for st, extra in cells:
        r = subprocess.run([sys.executable, "bench/phase4.py", *base, "--strategy", st, *extra],
                           capture_output=True, text=True)
        line = [l for l in r.stdout.splitlines() if l.startswith("RESULT ")]
        d = json.loads(line[0][7:]) if line else {"strategy": st, "error": r.stderr[-300:]}
        f.write(json.dumps(d) + "\n"); f.flush()
        print(st, d.get("run_s"), d.get("temp_MiB"), d.get("peak_rss_MiB"), flush=True)
