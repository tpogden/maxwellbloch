## E1 cost of a gradient


**two_level** (dim 2, 1 velocity class(es), n_t = 1001, n_z = 200); one `solve()` takes 4.02 s.

| N | forward (s) | reverse (s) | reverse / forward | forward-mode (s) | FD through JAX (s) | FD through clerq (s, estimate) | reverse speed-up over FD (JAX) | over FD (clerq) |
|-|-|-|-|-|-|-|-|-|
| 2 | 0.592 (0.585–0.6) | 3.32 (3.32–3.41) | 5.6 | 3.68 (3.64–3.77) | 1.79 (1.79–1.83) | 12.1 | 0.5x | 4x |
| 5 | 0.574 (0.57–0.575) | 3.22 (3.21–3.22) | 5.6 | 4 (3.98–4) | 3.44 (3.44–3.46) | 24.1 | 1.1x | 7x |
| 10 | 0.728 (0.722–0.739) | 4.27 (4.27–4.3) | 5.9 | 5.82 (5.79–5.82) | 8.11 (8.09–8.12) | 44.2 | 1.9x | 10x |
| 50 | 0.676 (0.669–0.678) | 3.87 (3.86–3.9) | 5.7 | 8.35 (8.34–8.56) | 34.6 (34.5–34.6) | 205 | 8.9x | 53x |
| 100 | 0.618 (0.613–0.623) | 3.49 (3.49–3.51) | 5.7 | skipped | 61.7 (61.6–62.1) | 406 | 17.7x | 116x |

**lambda_eit** (dim 3, 1 velocity class(es), n_t = 201, n_z = 60); one `solve()` takes 0.399 s.

| N | forward (s) | reverse (s) | reverse / forward | forward-mode (s) | FD through JAX (s) | FD through clerq (s, estimate) | reverse speed-up over FD (JAX) | over FD (clerq) |
|-|-|-|-|-|-|-|-|-|
| 2 | 0.12 (0.118–0.126) | 0.697 (0.687–0.706) | 5.8 | 0.679 (0.656–0.716) | 0.361 (0.357–0.369) | 1.2 | 0.5x | 2x |
| 5 | 0.125 (0.125–0.127) | 0.728 (0.725–0.73) | 5.8 | 1.21 (1.21–1.21) | 0.75 (0.747–0.756) | 2.4 | 1.0x | 3x |
| 10 | 0.133 (0.132–0.133) | 0.753 (0.744–0.755) | 5.7 | 2.1 (2.1–2.11) | 1.45 (1.45–1.45) | 4.39 | 1.9x | 6x |
| 50 | 0.132 (0.132–0.133) | 0.762 (0.754–0.774) | 5.8 | 4.31 (4.31–4.32) | 6.77 (6.76–6.79) | 20.4 | 8.9x | 27x |
| 100 | 0.145 (0.144–0.15) | 0.832 (0.826–0.858) | 5.7 | 9.05 (9.05–9.32) | 14.7 (14.6–14.8) | 40.3 | 17.6x | 48x |

**rb87** (dim 11, 11 velocity class(es), n_t = 201, n_z = 50); one `solve()` takes 22.1 s.

| N | forward (s) | reverse (s) | reverse / forward | forward-mode (s) | FD through JAX (s) | FD through clerq (s, estimate) | reverse speed-up over FD (JAX) | over FD (clerq) |
|-|-|-|-|-|-|-|-|-|
| 2 | 1.16 (1.16–1.16) | 7.07 (7.06–7.08) | 6.1 | 13.2 (13.1–13.2) | 3.49 (3.49–3.5) | 66.3 | 0.5x | 9x |
| 5 | 1.29 (1.29–1.3) | 7.83 (7.82–7.86) | 6.1 | 20.9 (20.9–21) | 7.78 (7.72–7.79) | 133 | 1.0x | 17x |
| 10 | 1.24 (1.23–1.24) | 7.45 (7.44–7.5) | 6.0 | 34.3 (34.3–34.4) | 13.4 (13.4–13.5) | 243 | 1.8x | 33x |
| 50 | 1.18 (1.17–1.2) | 7.11 (7.1–7.27) | 6.0 | 108 (108–109) | 60 (59.9–60.7) | 1.13e+03 | 8.4x | 158x |


## E6 before and after PR 361


| system | N | method | before PR 361 (s) | after (s) | after / before |
|-|-|-|-|-|-|
| two_level | 10 | forward | 0.723 (0.721–0.742) | 0.74 (0.721–0.744) | 1.02 |
| two_level | 10 | reverse | 4.23 (4.2–4.28) | 4.26 (4.2–4.26) | 1.01 |
| two_level | 10 | forward_ad | 5.62 (5.57–5.69) | 5.92 (5.85–5.93) | 1.05 |
| two_level | 50 | forward | 0.668 (0.658–0.674) | 0.684 (0.67–0.69) | 1.02 |
| two_level | 50 | reverse | 3.79 (3.77–3.96) | 3.84 (3.8–4) | 1.01 |
| two_level | 50 | forward_ad | 8.28 (8.26–8.48) | 8.25 (8.24–8.34) | 1.00 |
| lambda_eit | 10 | forward | 0.131 (0.13–0.135) | 0.131 (0.131–0.136) | 1.00 |
| lambda_eit | 10 | reverse | 0.75 (0.733–0.756) | 0.752 (0.745–0.775) | 1.00 |
| lambda_eit | 10 | forward_ad | 1.96 (1.96–1.97) | 2.12 (2.12–2.13) | 1.08 |
| lambda_eit | 50 | forward | 0.13 (0.129–0.132) | 0.132 (0.131–0.132) | 1.02 |
| lambda_eit | 50 | reverse | 0.759 (0.749–0.765) | 0.763 (0.751–0.767) | 1.01 |
| lambda_eit | 50 | forward_ad | 4.3 (4.29–4.35) | 4.31 (4.3–4.31) | 1.00 |
| rb87 | 10 | forward | 1.22 (1.21–1.23) | 1.24 (1.22–1.27) | 1.01 |
| rb87 | 10 | reverse | 7.4 (7.36–7.49) | 7.44 (7.44–7.65) | 1.01 |
| rb87 | 10 | forward_ad | 33.6 (33.3–35.3) | 34.3 (33.9–34.5) | 1.02 |
| rb87 | 50 | forward | 1.17 (1.17–1.17) | 1.18 (1.18–1.21) | 1.01 |
| rb87 | 50 | reverse | 7.03 (7.02–7.03) | 7.11 (7.1–7.27) | 1.01 |
| rb87 | 50 | forward_ad | 109 (108–109) | 108 (108–110) | 1.00 |

## E2 memory


| system | sweep | n_z | classes | strategy | peak resident set (MiB) | above baseline (MiB) | XLA temporary (MiB) | steady gradient (s) |
|-|-|-|-|-|-|-|-|-|
| lambda_eit | n_z | 50 | 1 | s2 | 875 | 333 | 20.2 | 0.52 |
| lambda_eit | n_z | 50 | 1 | s3 | 921 | 381 | 2.4 | 0.64 |
| lambda_eit | n_z | 50 | 1 | s4 | 930 | 387 | 1.6 | 0.73 |
| lambda_eit | n_z | 100 | 1 | s2 | 891 | 349 | 39.6 | 1.05 |
| lambda_eit | n_z | 100 | 1 | s3 | 929 | 386 | 3.6 | 1.28 |
| lambda_eit | n_z | 100 | 1 | s4 | 950 | 412 | 2.8 | 1.46 |
| lambda_eit | n_z | 200 | 1 | s2 | 924 | 381 | 78.5 | 2.05 |
| lambda_eit | n_z | 200 | 1 | s3 | 927 | 386 | 6.1 | 2.51 |
| lambda_eit | n_z | 200 | 1 | s4 | 931 | 386 | 5.3 | 2.83 |
| lambda_eit | n_z | 400 | 1 | s2 | 1009 | 461 | 156.2 | 4.11 |
| lambda_eit | n_z | 400 | 1 | s3 | 939 | 396 | 11.0 | 4.99 |
| lambda_eit | n_z | 400 | 1 | s4 | 953 | 406 | 10.2 | 5.65 |
| lambda_eit | n_z | 800 | 1 | s2 | 1171 | 621 | 311.7 | 8.64 |
| lambda_eit | n_z | 800 | 1 | s3 | 953 | 399 | 20.8 | 9.89 |
| lambda_eit | n_z | 800 | 1 | s4 | 940 | 390 | 20.0 | 11.40 |
| rb87 | n_z | 50 | 11 | s2 | 2710 | 2187 | 1881.3 | 6.24 |
| rb87 | n_z | 50 | 11 | s3 | 1012 | 488 | 101.2 | 7.37 |
| rb87 | n_z | 50 | 11 | s4 | 910 | 384 | 8.5 | 7.83 |
| rb87 | n_z | 100 | 11 | s2 | 4526 | 4002 | 3698.4 | 12.49 |
| rb87 | n_z | 100 | 11 | s3 | 1022 | 502 | 101.8 | 14.99 |
| rb87 | n_z | 100 | 11 | s4 | 922 | 402 | 9.1 | 15.75 |
| rb87 | n_z | 200 | 11 | s2 | 8236 | 7701 | 7332.6 | 25.38 |
| rb87 | n_z | 200 | 11 | s3 | 1020 | 497 | 103.0 | 29.61 |
| rb87 | n_z | 200 | 11 | s4 | 922 | 393 | 10.3 | 31.34 |
| rb87 | n_z | 400 | 11 | s2 | 15461 | 14934 | 14601.0 | 51.39 |
| rb87 | n_z | 400 | 11 | s3 | 1013 | 485 | 105.5 | 59.27 |
| rb87 | n_z | 400 | 11 | s4 | 921 | 394 | 12.8 | 62.83 |
| rb87 | n_z | 800 | 11 | s3 | 1013 | 482 | 110.4 | 119.06 |
| rb87 | n_z | 800 | 11 | s4 | 930 | 404 | 17.7 | 125.50 |
| rb87 | n_v | 100 | 1 | s2 | 1165 | 648 | 346.7 | 1.15 |
| rb87 | n_v | 100 | 1 | s3 | 902 | 384 | 10.7 | 1.48 |
| rb87 | n_v | 100 | 1 | s4 | 900 | 381 | 2.1 | 1.30 |
| rb87 | n_v | 100 | 11 | s2 | 4528 | 4004 | 3698.4 | 12.52 |
| rb87 | n_v | 100 | 11 | s3 | 1007 | 480 | 101.8 | 14.83 |
| rb87 | n_v | 100 | 11 | s4 | 921 | 395 | 9.1 | 15.63 |
| rb87 | n_v | 100 | 25 | s2 | 9227 | 8698 | 8391.1 | 38.01 |
| rb87 | n_v | 100 | 25 | s3 | 1124 | 599 | 229.3 | 44.32 |
| rb87 | n_v | 100 | 25 | s4 | 934 | 403 | 18.9 | 45.93 |
| rb87 | n_z | 800 | – | s2 | s2 stores every slice; too large to run here | | | |

## E3 accuracy


| system | N | fixed steps | field error vs `solve()` | AD vs finite differences (best step) | forward vs reverse | adaptive vs fixed-step AD |
|-|-|-|-|-|-|-|
| two_level | 10 | 2003 | 2.5e-04 | 8.0e-09 | 8.5e-14 | 2.5e-05 |
| lambda_eit | 10 | 2003 | 8.9e-08 | 2.3e-12 | 1.1e-16 | 1.9e-07 |
| rb87 | 6 | 1009 | 3.6e-04 | 5.4e-13 | 6.5e-16 | 4.8e-08 |

## E4 compile time


| system | N | forward: first call (s) | steady (s) | reverse gradient: first call (s) | steady (s) |
|-|-|-|-|-|-|
| two_level | 10 | 1.83 | 0.732 | 8.60 | 4.272 |
| lambda_eit | 10 | 1.42 | 0.134 | 5.53 | 0.744 |
| rb87 | 10 | 2.44 | 1.242 | 12.13 | 7.472 |

## E5 time to solution


- 44-parameter gradient: reverse AD 36.3 s, forward differences 279 s, one forward evaluation 6.23 s.
- Baseline (5 parameters, Nelder-Mead): 0.0681 retrieved after 296 evaluations, 2261 s.
- Free-form (44 parameters, Adam): 0.1057 after 40 steps, 1809 s. It passes the baseline value (0.0681) at step 23, after 1037 s; with finite-difference gradients the same steps would take about 6425 s (estimate from the measured gradient time).

## E7 equal-cost comparison with gradient-free search


44 parameters, budget 550 objective evaluations (a gradient call is charged 6.5), optical depth 25, no Doppler.

| optimiser | per seed | median | cost used | wall time per run (s) |
|-|-|-|-|-|
| adam | 0.459, 0.439, 0.447 | 0.447 | 551 | 704 |
| lbfgsb | 0.000, 0.462, 0.000 | 0.000 | 7 | 6 |
| de | 0.315, 0.336, 0.299 | 0.315 | 528 | 868 |
| cmaes | 0.399, 0.385, 0.404 | 0.399 | 555 | 804 |
| lhs | 0.250, 0.259, 0.232 | 0.250 | 550 | 828 |
