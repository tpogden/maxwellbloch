import json, subprocess, sys
base = ["--system", "j1_2_res", "--n_z", "60", "--n_t", "500", "--n_v", "10", "--repeat", "2"]
cells = [("s2", []), ("s3", [])] \
      + [("s4", ["--checkpoints", str(k)]) for k in (4, 8, 16, 32, 64)] \
      + [("s5", ["--chunk", str(c)]) for c in (4, 8, 16, 60)]
rows = []
for st, extra in cells:
    r = subprocess.run([sys.executable, "bench/phase4.py", *base, "--strategy", st, *extra],
                       capture_output=True, text=True)
    line = [l for l in r.stdout.splitlines() if l.startswith("RESULT ")]
    d = json.loads(line[0][7:]) if line else {"strategy": st, "error": r.stderr[-200:]}
    d["extra"] = extra
    rows.append(d)
    print(f"{st:3} {' '.join(extra):16} grad {d.get('run_s', float('nan')):6.2f} s  temp {d.get('temp_MiB', float('nan')):7.1f} MiB  rss {d.get('peak_rss_MiB', float('nan')):5.0f}  |g| {d.get('grad_norm', float('nan')):.6f}", flush=True)
json.dump(rows, open("bench/results/tune_s4.json", "w"), indent=1)
