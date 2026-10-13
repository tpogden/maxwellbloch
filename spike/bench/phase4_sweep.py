"""T4.3: compile-time memory matrix over the four sweep axes (base: j1_2_res, n_z 200, n_t 1000,
n_v 10, max_steps 16384). Cells run as subprocesses, 4 in parallel, results appended as JSONL."""
import json, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

SYS = {2: "two_level", 4: "j12_12", 9: "j1_2_res", 16: "j72_72"}
BASE = dict(dim=9, n_z=200, n_t=1000, n_v=10)
AXES = {"n_z": [50, 100, 200, 500, 1000, 2000], "n_t": [200, 500, 1000, 2000, 5000],
        "dim": [2, 4, 9, 16], "n_v": [1, 10, 50, 100]}
STRATS = [("s1", []), ("s2", []), ("s3", []), ("s4", ["--checkpoints", "32"]), ("s5", [])]

cells = []
for axis, vals in AXES.items():
    for v in vals:
        cfg = {**BASE, axis: v}
        for st, extra in STRATS:
            cells.append((axis, v, cfg, st, extra))
cells.append(("base", 0, BASE, "direct", []))

def run(cell):
    axis, v, cfg, st, extra = cell
    cmd = [sys.executable, "bench/phase4.py", "--system", SYS[cfg["dim"]], "--n_z", str(cfg["n_z"]),
           "--n_t", str(cfg["n_t"]), "--n_v", str(cfg["n_v"]), "--strategy", st, "--compile_only", *extra]
    r = subprocess.run(cmd, capture_output=True, text=True)
    line = [l for l in r.stdout.splitlines() if l.startswith("RESULT ")]
    d = json.loads(line[0][7:]) if line else {"error": r.stderr[-300:]}
    d.update(axis=axis, value=v, strategy_label=st, dim=cfg["dim"])
    return d

with ThreadPoolExecutor(4) as ex, open("bench/results/memory_matrix.jsonl", "w") as f:
    for d in ex.map(run, cells):
        f.write(json.dumps(d) + "\n"); f.flush()
        print(d["axis"], d["value"], d["strategy_label"], round(d.get("temp_MiB", -1), 1), flush=True)
