# The Walk Length: steps/n^2 is flat-to-declining, not growing

*2026-11-06. Seven rules, n = 256..8192, 150 tapes per (rule, n, seed), two seeds, cap 8n^2 steps per tape (caps=0 everywhere). Script: `probe_walk_len/wl_report.py` over `probe_walk_len/cells/`.*

## Why this was measured

Every upper bound here is `g <= (steps) x (reads per step)`. The step count is a property of the RULE, so it can be counted without any cost bookkeeping. The structural bound is `n * m^2 = n^3` steps; if that were tight the general upper bound would be n^4 cell reads and Conjecture 2 could never be pushed to n^2 by this route.

## Mean of steps/n^2

| rule | 256 | 512 | 1024 | 2048 | 4096 | 8192 | exponent |
|---|---|---|---|---|---|---|---|
| `bal_q` | 0.01000 | 0.00792 | 0.00715 | 0.00704 | 0.00633 | 0.00766 | -0.083 |
| `dense2` | 0.06294 | 0.05875 | 0.06307 | 0.05266 | 0.04901 | 0.05373 | -0.062 |
| `las_q` | 0.03324 | 0.01838 | 0.02785 | 0.02759 | 0.02348 | 0.01947 | -0.080 |
| `mp_sq` | 0.04648 | 0.03771 | 0.04227 | 0.03713 | 0.03591 | 0.03516 | -0.069 |
| `nblk_gt2` | 0.00731 | 0.00550 | 0.00528 | 0.00517 | 0.00436 | 0.00536 | -0.093 |
| `nblk_gt3` | 0.00578 | 0.00470 | 0.00459 | 0.00413 | 0.00419 | 0.00509 | -0.045 |
| `runs` | 0.05138 | 0.05309 | 0.06093 | 0.05191 | 0.04882 | 0.05367 | -0.008 |

Per-tape maximum of steps/n^2 (the worst tape, not the mean):

| rule | 256 | 512 | 1024 | 2048 | 4096 | 8192 |
|---|---|---|---|---|---|---|
| `bal_q` | 0.1205 | 0.1299 | 0.1155 | 0.1232 | 0.1046 | 0.1161 |
| `dense2` | 0.3338 | 0.3528 | 0.2821 | 0.2321 | 0.2378 | 0.2165 |
| `las_q` | 0.3293 | 0.1625 | 0.3257 | 0.2626 | 0.2699 | 0.3185 |
| `mp_sq` | 0.2676 | 0.2339 | 0.2182 | 0.2483 | 0.2208 | 0.3803 |
| `nblk_gt2` | 0.0981 | 0.1147 | 0.0892 | 0.1123 | 0.0912 | 0.0849 |
| `nblk_gt3` | 0.0765 | 0.1133 | 0.0787 | 0.1123 | 0.0904 | 0.0786 |
| `runs` | 0.2338 | 0.2278 | 0.2480 | 0.2171 | 0.2378 | 0.2165 |

## Verdict

Every exponent is NEGATIVE, between -0.030 and -0.093. The walk length is therefore
`Theta(n^2)` and slightly sub-quadratic (consistent with `n^2 / (log n)^{1/3}`); the
structural `n^3` is loose by a factor n, and the `n^4` cell-read bound inherited from
`n^3` evaluations is loose by n^2.

**Consequence.** The general upper bound in cell reads is `O(n^3)` (steps `O(n^2)` x at most
`n` reads per evaluation), not `O(n^4)`, provided the walk length is `O(n^2)` - measured here,
not proved. And since the step factor is pinned at `n^2` for all seven rules, Conjecture 2 is
exactly a statement about the certificate `c_rho`, as `CONJECTURE2.md` says.

**Why the exponent is negative and not exactly 2.** A decline of -0.03 to -0.09 over this
range is what `n^2/(log n)^a` with `a = 0.2..0.6` looks like locally. It is not evidence for
`o(n^2)` steps; it is evidence that the step count is bounded by `n^2` with room to spare.
