"""T5.1: gradient cost vs n_params on the reference problem (two_level, n_z 200, n_t 1000)."""
import json, subprocess, sys
NS = [1, 10, 50, 100, 500]
METHODS = ["fwd", "ad_s2", "ad_s3", "jacfwd", "fd_jax", "fd_clerq"]
with open("bench/results/scaling.jsonl", "a") as f:
    for n in NS:
        for m in METHODS:
            try:
                r = subprocess.run([sys.executable, "bench/phase5_cell.py", "--n", str(n), "--method", m,
                                    "--repeat", "3"], capture_output=True, text=True, timeout=5400)
                line = [l for l in r.stdout.splitlines() if l.startswith("RESULT ")]
                d = json.loads(line[0][7:]) if line else {"n": n, "method": m, "error": r.stderr[-300:]}
            except subprocess.TimeoutExpired:
                d = {"n": n, "method": m, "error": "timeout"}
            f.write(json.dumps(d) + "\n"); f.flush()
            print(n, m, round(d.get("s", -1), 3), d.get("extrapolated", ""), d.get("error", ""), flush=True)
