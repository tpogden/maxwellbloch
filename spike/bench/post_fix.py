"""After merging clerq master (Doppler + dipole-weight fixes): full tests, dim-9 validity check
against clerq-FD, and a re-timing of the large system under the corrected convention."""
import json, subprocess, sys, time
steps = [
    ("pytest", [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider"]),
    ("validity_dim9", [sys.executable, "bench/phase5_cell.py", "--n", "5", "--method", "fd_clerq",
                       "--system", "j1_2_res", "--n_v", "50", "--n_z", "4", "--real_fd_max", "5"]),
    ("fwd_N10", [sys.executable, "bench/phase5_cell.py", "--n", "10", "--method", "fwd",
                 "--system", "j1_2_res", "--n_v", "50", "--repeat", "2"]),
    ("ad_s2_N10", [sys.executable, "bench/phase5_cell.py", "--n", "10", "--method", "ad_s2",
                   "--system", "j1_2_res", "--n_v", "50", "--repeat", "2"]),
    ("ad_s3_N10", [sys.executable, "bench/phase5_cell.py", "--n", "10", "--method", "ad_s3",
                   "--system", "j1_2_res", "--n_v", "50", "--repeat", "2"]),
]
for name, cmd in steps:
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    tail = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr[-300:]
    line = [l for l in r.stdout.splitlines() if l.startswith("RESULT ")]
    rec = {"step": name, "returncode": r.returncode, "wall_s": time.time() - t0,
           "result": json.loads(line[0][7:]) if line else None, "tail": tail}
    open("bench/results/post_fix.jsonl", "a").write(json.dumps(rec) + "\n")
    print(name, r.returncode, round(time.time() - t0), tail[:200], flush=True)
