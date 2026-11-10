# MANIFEST — every headline number, the probe that produced it

Generated from `claims.json`; the paper is checked against it by `check_paper.py` before every build.

| claim | value | source |
|---|---|---|
| `step_ratio_max` | 15.58 | probe_master/step counts, hashmod seed 2026 |
| `step_ratio_min` | 0.385 | same ensemble, NBLK>3 |
| `monotone_mean` | 2n-1 exactly | Theorem thm:collapse; measured 1.9998n at 4096 |
| `pattern_ces` | 2.57 n | probe_ces_exact/ce_patch.py, 200 meeting tapes at 4096 |
| `runs_ces` | 2.77 n | probe_ces_exact/ce_patch.py (w2_256) |
| `canonical_cost` | 147.7 n | pcr_lib.Runner, n=1024 |
| `canonical_c1a` | 116.5 n, the c1a variant | audit2.py |
| `hashmod_cell` | 0.08-0.096 n^2 | probe_master (cell reads; the 0.62 n ln n figure is rows) |
| `longrun_witness` | n^3/21 | probe_g_witness/gw.py, four sizes, coefficient 20.5-21.9 |
| `longrun_separation` | g/e = 86232 at 4096 | probe_g_witness |
| `cx_order` | n^3/log n for full-rank kernels at s=Theta(log n) | Theorem thm:rank |
| `cx_measured` | n^2.93, an upper bound | probe_cx_down/CX_DOWN.md, s-slot accumulator |
| `ballast_mean` | 6.2 n = 5.44 n + n^-1 0.75 n^2 | probe_semilinear/SEMILINEAR.md |
| `ballast_exponents` | e n^0.997, g n^1.995, G n^1.000 | probe_semilinear |
| `ballast_ratios` | e/n 6.11-6.24, g/n^2 0.751-0.759 | probe_semilinear |
| `info_ceiling_sum` | 1.893 reads/query | Theorem thm:info (3-valued alphabet) |
| `info_ceiling_four` | 1.5 reads/query | Theorem thm:info (four-state cell) |
| `certificate_sum` | 1285703 cells vs union 1020 | probe_ceiling |
| `transfer_slope` | f = 0.104 at n=128, f ~ n^-0.99 | probe_r2k |
| `peak_exact` | 1496/1500 at W2=128; 2.0 n | probe_peak_fix/PEAK_FIX.md |
| `n2_rung` | <= 1.33 n^2 (measured 1.12) | probe_class_rerun |
| `defects` | nine defects, encodings cover eight | PC_PUSH 2.2 |
| `kappa_range` | 0.472-0.925 at c=16, n=512 | probe_kappa_inc2/KAPPA_INC2.md |
| `kappa_star` | kappa* >= 0.93 | same |
| `g_over_e_fl` | 0.997-1.000 | probe_r2k/gfl.py |
| `g_over_e_super` | G/e <= 0.26 for the band rules; 0.997-1.000 for the fully linear ones | probe_r2k/gsup.py |
| `gates_cube` | 88560 tapes at n<=10 | probe gates |
| `gates_windows` | 2192103 windows | probe gates |
| `gates_four` | 130500 walks | probe gates |
| `gates_endgame` | 147555 walks | probe gates |
| `gates_epochs` | 27984 epochs | probe gates |
| `class_census` | 21 natural rules: 4 fully linear, 17 superlinear; 1 constructed semi-linear | catalogue table |

## Spelling and class bans

`check_paper.py` also refuses: `quasi-linear`, `weak superlinear`, `strong superlinear`, `pointer class`, `walk class`, `four level-count`, `four classes`, `n^{1/\log n}`, `\\Omega(n^2)$ & \\textbf{open}`.

## Directory map

- `probe_vandermonde/` — the field/order checks behind the rank argument (C2).
- `probe_cost_dist/` — the per-tape cost densities of Figure 1.
- `probe_semilinear/` — the `ballast` construction and its exponents.
- `probe_r2k/` — the bootstrap, `G`, and the `g/e` measurements.
- `probe_audit/` — the two audit reports.
- `probe_g_witness/`, `probe_k300/`, `probe_e_cheap/`, `probe_ces_exact/`, `probe_en_linear/`, `probe_peak_fix/` — the cost measurements behind the tables.

Build the paper with `bash build.sh`; it runs the consistency check first and refuses to
compile on any hit.
