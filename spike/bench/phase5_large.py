"""Gate 4 on a larger system: j1_2_res (dim 9), n_v 50. Sequential cells.

clerq forward costs ~50 min at n_z=200 here (15 s per slice), so clerq-FD validity is checked at
n_z=8 and the clerq forward is extrapolated from the per-slice cost."""
import json, subprocess, sys
S = ["--system", "j1_2_res", "--n_v", "50"]
cells = [(2, "fd_clerq", ["--n_z", "8", "--real_fd_max", "2"])]     # validity + per-slice clerq cost
for n in (1, 10, 50, 100):
    cells += [(n, "fwd", []), (n, "ad_s2", []), (n, "ad_s3", [])]
cells += [(1, "jacfwd", []), (10, "jacfwd", []), (50, "jacfwd", ["--repeat", "1"]),
          (1, "fd_jax", ["--repeat", "1"]), (10, "fd_jax", ["--repeat", "1"])]
with open("bench/results/scaling_large.jsonl", "a") as f:
    for n, m, extra in cells:
        rep = ["--repeat", "2"] if "--repeat" not in extra else []
        try:
            r = subprocess.run([sys.executable, "bench/phase5_cell.py", "--n", str(n), "--method", m,
                                *S, *rep, *extra], capture_output=True, text=True, timeout=7200)
            line = [l for l in r.stdout.splitlines() if l.startswith("RESULT ")]
            d = json.loads(line[0][7:]) if line else {"n": n, "method": m, "error": r.stderr[-300:]}
        except subprocess.TimeoutExpired:
            d = {"n": n, "method": m, "error": "timeout"}
        f.write(json.dumps(d) + "\n"); f.flush()
        print(n, m, round(d.get("s", -1), 2), d.get("grad_rel_err_ad_vs_clerq_fd", ""), d.get("error", ""), flush=True)
