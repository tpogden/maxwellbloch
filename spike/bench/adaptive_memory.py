"""Concern #2: does adaptive stepping (PID) change the memory picture? Compile-only."""
import json, subprocess, sys

rows = []
for nz, nt, dim, nv in [(200, 1000, 4, 1), (500, 2000, 9, 50)]:
    for adaptive, ms in [(False, None), (True, None), (True, 4096), (True, 16384), (True, 65536)]:
        for st in ("s2", "s3"):
            cmd = [sys.executable, "bench/probe_memory.py", "--n_z", str(nz), "--n_t", str(nt),
                   "--dim", str(dim), "--n_v", str(nv), "--strategy", st]
            if adaptive:
                cmd.append("--adaptive")
            if ms:
                cmd += ["--max_steps", str(ms)]
            r = subprocess.run(cmd, capture_output=True, text=True)
            line = [l for l in r.stdout.splitlines() if l.startswith("RESULT ")]
            if not line:
                rows.append({"nz": nz, "nt": nt, "dim": dim, "nv": nv, "st": st,
                             "adaptive": adaptive, "max_steps": ms, "error": r.stderr[-200:]})
            else:
                d = json.loads(line[0][7:])
                rows.append({"nz": nz, "nt": nt, "dim": dim, "nv": nv, "st": st,
                             "adaptive": adaptive, "max_steps": ms, "temp_MiB": round(d["temp_MiB"], 1)})
            print(rows[-1], flush=True)
json.dump(rows, open("bench/results/adaptive_memory.json", "w"), indent=1)
