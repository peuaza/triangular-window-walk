# The repair as a real trim: BTHM-C when a repeat removes a tape cell

`tr_lib.tapes(n)`, 300 random tapes per cell, n = 512 / 1024 / 2048. Reads are `pcr_lib`'s INC accounting in both models. A **cycle** is a repeated quadruple (ell, r, L, R) - the state of the walk, since the tape is fixed until a repair removes a cell.

**The repair.** On a repeat the rule removes the boundary cell of the TAPE - `L <- L+1` or `R <- R-1` - on the side its own action table selects for that state (the entry `TABLE[2*b2+b1]` the step block would have used, inverted for the `c1a` rules, exactly as the `trim=5` variant reads it), re-opens the window on the narrower tape, pays the re-span a rule-initiated cut pays, and continues. The tape only shrinks, so a run admits at most `n` repairs - the bound the theorem states. This is the companion paper's trim, not a window retraction.

BTHM-0 has no memory, so a repeated quadruple is a closed period and the walk dies at the step cap `40n^2+100`; its numbers here are a re-run of the same walk and are checked tape by tape against `probe_cycle` (see Verification). `first repeat` is how far the tape got before the period closed.

Over the 18000 measured (rule, n, tape) runs: **BTHM-0 meets 2681 (14.9%)**, dies at the cap on 11146 (61.9%); **BTHM-C with the cut repair meets 6796 (37.8%) and empties the region on 11204 (62.2%)**. Every one of the 11146 cap-diers terminates (4115 in a meet, 7031 in an empty region); **0 of 18000 cut-repair runs hit the step cap** and 0 ran out of oracle table.

The largest number of repairs on any single tape is **1956** (n=512: 485, n=1024: 1002, n=2048: 1956, i.e. 0.96 n); no tape ever needed n repairs, let alone more, as the theorem requires.

## Verification (before any number below)

- **BTHM-0 reproduces exactly.** Every BTHM-0 field of every tape (`end0, steps0, reads0, cyc_start, cyc_len, width, cuts0`) was compared with the stored `probe_cycle/cyc_<rule>_<n>.json` rows: 18000 tape records, **0 mismatches**. The no-memory walk is untouched by the repair. The derived BTHM-0 cell fields (22 per cell, including `cyc_len_med` and `cyc_width_med`) were compared with `probe_cycle/cyc_cells.json` too: 1320 fields, **0 differences**.
- **The retraction run still reproduces.** `cut_engine.py` is `cyc_engine.py` with one branch added; re-running mode 1 (trim=1) on all 300 tapes of all 20 rules at n=512 reproduces the stored retraction fields exactly (`cut_ctl.py`).
- **The walk is still the record's walk.** `cut_verify.py`: 99/99 INC cells exact against `probe_g_witness/gw_cells.json` at n=512 (reads, delta, end) and 120/120 (delta, end) against `tr_lib`.
- **The repair branch has a second implementation.** `cut_indep.py` re-derives the cut run with a dict-based seen-set and its own repair loop in plain Python: 30/30 tape-runs agree on all six fields (end, steps, reads, ntrim, nl, nr).
- **Structural invariants over all 18000 cut runs.** A tape that never closed a period in BTHM-0 must run identically under the bit (`invariant_viol = 0`); every repair is either a left or a right cell (`nl+nr = ntrim`, violations 0); no run exceeds n repairs (violations 0).

## n = 512

| rule | B0 meet | B0 cap | first repeat | cycle len | width | CUT meet | CUT empty | rescued | trims med/max (all 300) | reads 0 (meet) | reads C (meet) | reads C (all) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `canonical` | 0.170 | 0.620 | 871 | 2 | 98 | 0.393 | 0.607 | 0.360 | 4/437 | 23483 | 40668 | 36440 |
| `runs` | 0.237 | 0.447 | 14870 | 2 | 16 | 0.360 | 0.640 | 0.276 | 0/437 | 21610 | 42652 | 48970 |
| `dense2` | 0.187 | 0.503 | 9081 | 14 | 104 | 0.440 | 0.560 | 0.503 | 1/316 | 12540 | 54548 | 54711 |
| `peak` | 0.070 | 0.780 | 925 | 2 | 132 | 0.340 | 0.660 | 0.346 | 36/428 | 3312 | 48311 | 47815 |
| `trip` | 0.080 | 0.643 | 1587 | 2 | 102 | 0.297 | 0.703 | 0.337 | 23/432 | 4730 | 69612 | 29459 |
| `oddrun` | 0.090 | 0.697 | 1180 | 2 | 121 | 0.303 | 0.697 | 0.306 | 21/463 | 4239 | 39945 | 35811 |
| `pattern` | 0.317 | 0.283 | 17593 | 6 | 17 | 0.420 | 0.580 | 0.365 | 0/324 | 25280 | 43399 | 43399 |
| `LAS` | 0.077 | 0.733 | 1217 | 2 | 107 | 0.447 | 0.553 | 0.504 | 31/475 | 7555 | 54918 | 41568 |
| `MAXP` | 0.117 | 0.470 | 5563 | 4 | 106 | 0.297 | 0.703 | 0.383 | 0/453 | 7572 | 36192 | 27860 |
| `KCHG` | 0.057 | 0.777 | 468 | 2 | 138 | 0.333 | 0.667 | 0.356 | 62/447 | 3143 | 46216 | 42217 |
| `BAL>m/2` | 0.057 | 0.840 | 399 | 2 | 190 | 0.393 | 0.607 | 0.401 | 75/485 | 3163 | 36739 | 36552 |
| `BAL>m/4` | 0.040 | 0.850 | 396 | 4 | 189 | 0.357 | 0.643 | 0.372 | 78/448 | 4038 | 64046 | 40818 |
| `NBLK=1` | 0.033 | 0.773 | 3443 | 4 | 168 | 0.263 | 0.737 | 0.297 | 58/466 | 3103 | 37843 | 27590 |
| `NBLK>2` | 0.023 | 0.857 | 333 | 4 | 205 | 0.343 | 0.657 | 0.373 | 87/442 | 2886 | 47501 | 35440 |
| `NBLK>3` | 0.023 | 0.847 | 309 | 4 | 218 | 0.350 | 0.650 | 0.386 | 94/468 | 2636 | 40428 | 37686 |
| `AMAX` | 0.037 | 0.703 | 751 | 4 | 128 | 0.293 | 0.707 | 0.365 | 15/441 | 2726 | 23964 | 23096 |
| `AMIN` | 0.063 | 0.697 | 1149 | 4 | 138 | 0.287 | 0.713 | 0.321 | 18/448 | 5484 | 24707 | 24324 |
| `longrun` | 0.220 | 0.063 | 11062 | 2 | 11 | 0.267 | 0.733 | 0.737 | 0/43 | 21604 | 22344 | 20221 |
| `pal` | 0.143 | 0.260 | 10509 | 2 | 5 | 0.270 | 0.730 | 0.487 | 0/459 | 32303 | 32132 | 22426 |
| `hashmod` | 1.000 | 0.000 | — | — | — | 1.000 | 0.000 | — | 0/0 | 65514 | 65514 | 65514 |
| `monotone` | 1.000 | 0.000 | — | — | — | 1.000 | 0.000 | — | 0/0 | 1023 | 1023 | 1023 |

## n = 1024

| rule | B0 meet | B0 cap | first repeat | cycle len | width | CUT meet | CUT empty | rescued | trims med/max (all 300) | reads 0 (meet) | reads C (meet) | reads C (all) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `canonical` | 0.110 | 0.647 | 5993 | 2 | 186 | 0.397 | 0.603 | 0.443 | 7/763 | 105595 | 185117 | 121726 |
| `runs` | 0.250 | 0.430 | 28021 | 2 | 17 | 0.367 | 0.633 | 0.271 | 0/888 | 38661 | 50892 | 118712 |
| `dense2` | 0.150 | 0.560 | 38396 | 13 | 175 | 0.387 | 0.613 | 0.423 | 2/395 | 10291 | 102728 | 124105 |
| `peak` | 0.073 | 0.803 | 1740 | 2 | 225 | 0.363 | 0.637 | 0.361 | 98/894 | 6672 | 54075 | 116667 |
| `trip` | 0.063 | 0.753 | 5136 | 2 | 186 | 0.310 | 0.690 | 0.327 | 50/894 | 5332 | 95084 | 95448 |
| `oddrun` | 0.080 | 0.700 | 8746 | 2 | 131 | 0.320 | 0.680 | 0.343 | 24/906 | 5992 | 79187 | 78941 |
| `pattern` | 0.337 | 0.227 | 26315 | 6 | 13 | 0.443 | 0.557 | 0.471 | 0/761 | 62189 | 62714 | 116204 |
| `LAS` | 0.063 | 0.817 | 1844 | 2 | 212 | 0.433 | 0.567 | 0.453 | 79/905 | 5480 | 121459 | 122254 |
| `MAXP` | 0.150 | 0.470 | 24878 | 4 | 158 | 0.327 | 0.673 | 0.376 | 0/915 | 9849 | 59907 | 61603 |
| `KCHG` | 0.040 | 0.873 | 993 | 2 | 200 | 0.383 | 0.617 | 0.393 | 109/958 | 5920 | 154443 | 120509 |
| `BAL>m/2` | 0.053 | 0.897 | 687 | 4 | 357 | 0.400 | 0.600 | 0.387 | 122/915 | 5524 | 96487 | 100859 |
| `BAL>m/4` | 0.047 | 0.887 | 775 | 4 | 300 | 0.377 | 0.623 | 0.372 | 136/946 | 5413 | 76611 | 103902 |
| `NBLK=1` | 0.040 | 0.780 | 6223 | 4 | 297 | 0.267 | 0.733 | 0.291 | 73/949 | 5285 | 43888 | 58981 |
| `NBLK>2` | 0.047 | 0.900 | 573 | 4 | 363 | 0.413 | 0.587 | 0.407 | 146/926 | 5173 | 53408 | 103682 |
| `NBLK>3` | 0.050 | 0.890 | 557 | 4 | 395 | 0.350 | 0.650 | 0.337 | 165/958 | 5200 | 52570 | 99905 |
| `AMAX` | 0.047 | 0.737 | 1311 | 4 | 231 | 0.260 | 0.740 | 0.290 | 33/970 | 6284 | 78388 | 61040 |
| `AMIN` | 0.067 | 0.707 | 2302 | 4 | 201 | 0.237 | 0.763 | 0.241 | 21/904 | 5832 | 46337 | 58576 |
| `longrun` | 0.187 | 0.077 | 26673 | 2 | 8 | 0.237 | 0.763 | 0.652 | 0/80 | 28880 | 29200 | 55027 |
| `pal` | 0.117 | 0.270 | 8931 | 2 | 5 | 0.237 | 0.763 | 0.444 | 0/1002 | 125603 | 39628 | 59461 |
| `hashmod` | 1.000 | 0.000 | — | — | — | 1.000 | 0.000 | — | 0/0 | 278335 | 278335 | 278335 |
| `monotone` | 1.000 | 0.000 | — | — | — | 1.000 | 0.000 | — | 0/0 | 2047 | 2047 | 2047 |

## n = 2048

| rule | B0 meet | B0 cap | first repeat | cycle len | width | CUT meet | CUT empty | rescued | trims med/max (all 300) | reads 0 (meet) | reads C (meet) | reads C (all) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `canonical` | 0.090 | 0.730 | 18546 | 2 | 345 | 0.427 | 0.573 | 0.461 | 23/1655 | 505241 | 718237 | 618975 |
| `runs` | 0.223 | 0.420 | 118109 | 2 | 18 | 0.370 | 0.630 | 0.349 | 0/1922 | 314286 | 314286 | 632205 |
| `dense2` | 0.167 | 0.560 | 156701 | 14 | 328 | 0.433 | 0.567 | 0.476 | 2/1377 | 67184 | 606979 | 649732 |
| `peak` | 0.080 | 0.817 | 6811 | 2 | 426 | 0.380 | 0.620 | 0.367 | 174/1837 | 118978 | 509217 | 590183 |
| `trip` | 0.053 | 0.803 | 33080 | 2 | 288 | 0.377 | 0.623 | 0.403 | 141/1845 | 14404 | 691447 | 753237 |
| `oddrun` | 0.060 | 0.687 | 60450 | 2 | 162 | 0.250 | 0.750 | 0.277 | 30/1838 | 27648 | 396847 | 503676 |
| `pattern` | 0.333 | 0.183 | 122343 | 6 | 13 | 0.407 | 0.593 | 0.400 | 0/180 | 341220 | 324318 | 618856 |
| `LAS` | 0.050 | 0.843 | 15765 | 2 | 480 | 0.437 | 0.563 | 0.459 | 144/1721 | 12988 | 513363 | 712815 |
| `MAXP` | 0.127 | 0.567 | 100513 | 4 | 316 | 0.303 | 0.697 | 0.312 | 7/1605 | 48710 | 465360 | 428084 |
| `KCHG` | 0.033 | 0.893 | 2646 | 2 | 499 | 0.367 | 0.633 | 0.373 | 261/1680 | 26569 | 835826 | 802238 |
| `BAL>m/2` | 0.023 | 0.930 | 1504 | 4 | 720 | 0.377 | 0.623 | 0.380 | 307/1835 | 61439 | 543777 | 466442 |
| `BAL>m/4` | 0.037 | 0.917 | 1490 | 4 | 727 | 0.337 | 0.663 | 0.327 | 295/1840 | 10672 | 523116 | 502748 |
| `NBLK=1` | 0.027 | 0.837 | 24073 | 4 | 777 | 0.277 | 0.723 | 0.299 | 307/1956 | 117912 | 282563 | 383603 |
| `NBLK>2` | 0.030 | 0.923 | 1220 | 4 | 749 | 0.347 | 0.653 | 0.343 | 397/1886 | 61439 | 312377 | 420430 |
| `NBLK>3` | 0.037 | 0.917 | 1191 | 4 | 740 | 0.370 | 0.630 | 0.364 | 393/1904 | 11567 | 326960 | 433148 |
| `AMAX` | 0.037 | 0.817 | 7606 | 4 | 417 | 0.303 | 0.697 | 0.327 | 129/1954 | 17424 | 917383 | 497506 |
| `AMIN` | 0.067 | 0.733 | 15926 | 4 | 409 | 0.320 | 0.680 | 0.345 | 73/1662 | 34678 | 567265 | 526504 |
| `longrun` | 0.270 | 0.063 | 16251 | 2 | 11 | 0.323 | 0.677 | 0.842 | 0/77 | 340580 | 301002 | 359648 |
| `pal` | 0.183 | 0.247 | 71580 | 2 | 5 | 0.290 | 0.710 | 0.432 | 0/770 | 340580 | 301002 | 359643 |
| `hashmod` | 1.000 | 0.000 | — | — | — | 1.000 | 0.000 | — | 0/0 | 1030077 | 1030077 | 1030077 |
| `monotone` | 1.000 | 0.000 | — | — | — | 1.000 | 0.000 | — | 0/0 | 4095 | 4095 | 4095 |

## The effect of the repair semantics, n = 512

Left column: the retraction run in `probe_cycle/cyc_cells.json`, where the repair retracted a window end and the walk could re-widen it (up to 3,627,840 repairs on one tape, more than n). Right column: the same walk with the cut repair, where the cell is gone for good.

| rule | meet (retraction) | meet (cut) | delta | rescued (retr.) | rescued (cut) | trims med/max, rescued (retr.) | trims med/max, rescued (cut) |
|---|---|---|---|---|---|---|
| `canonical` | 0.690 | 0.393 | -0.297 | 0.839 | 0.360 | 107/7395 | 36/437 |
| `runs` | 0.670 | 0.360 | -0.310 | 0.970 | 0.276 | 14/633 | 16/437 |
| `dense2` | 0.573 | 0.440 | -0.133 | 0.768 | 0.503 | 146/114778 | 19/316 |
| `peak` | 0.813 | 0.340 | -0.473 | 0.953 | 0.346 | 130/13095 | 73/428 |
| `trip` | 0.703 | 0.297 | -0.407 | 0.969 | 0.337 | 111/38559 | 112/415 |
| `oddrun` | 0.743 | 0.303 | -0.440 | 0.938 | 0.306 | 139/127141 | 66/390 |
| `pattern` | 0.600 | 0.420 | -0.180 | 1.000 | 0.365 | 27/171 | 4/324 |
| `LAS` | 0.777 | 0.447 | -0.330 | 0.955 | 0.504 | 119/192212 | 91/417 |
| `MAXP` | 0.543 | 0.297 | -0.247 | 0.908 | 0.383 | 144/18000 | 108/453 |
| `KCHG` | 0.827 | 0.333 | -0.493 | 0.991 | 0.356 | 181/7619 | 116/447 |
| `BAL>m/2` | 0.887 | 0.393 | -0.493 | 0.988 | 0.401 | 178/72094 | 113/485 |
| `BAL>m/4` | 0.880 | 0.357 | -0.523 | 0.988 | 0.372 | 139/52613 | 124/448 |
| `NBLK=1` | 0.770 | 0.263 | -0.507 | 0.953 | 0.297 | 336/8970 | 151/456 |
| `NBLK>2` | 0.857 | 0.343 | -0.513 | 0.973 | 0.373 | 72/3223 | 141/428 |
| `NBLK>3` | 0.847 | 0.350 | -0.497 | 0.972 | 0.386 | 104/3046 | 117/468 |
| `AMAX` | 0.740 | 0.293 | -0.447 | 1.000 | 0.365 | 120/601 | 45/441 |
| `AMIN` | 0.760 | 0.287 | -0.473 | 1.000 | 0.321 | 116/724 | 54/446 |
| `longrun` | 0.277 | 0.267 | -0.010 | 0.895 | 0.737 | 8/58 | 5/43 |
| `pal` | 0.400 | 0.270 | -0.130 | 0.987 | 0.487 | 2/13 | 29/457 |
| `hashmod` | 1.000 | 1.000 | +0.000 | — | — | —/— | —/— |

Median change in the meeting fraction at n=512: **-0.407** (range -0.523 to +0.000).

Median change in the median all-tape reads at n=512: **3.29x** (range 1.00x to 10.07x).

## Fitted exponent of the reads over 512 / 1024 / 2048

Least squares of log(median reads) on log n, only where all three cells exist and the median is positive.

| rule | BTHM-0 reads (meeting tapes) | CUT reads (meeting tapes) | CUT reads (all 300) |
|---|---|---|---|
| `canonical` | 2.214 | 2.071 | 2.043 |
| `runs` | 1.931 | 1.441 | 1.845 |
| `dense2` | 1.211 | 1.738 | 1.785 |
| `peak` | 2.583 | 1.699 | 1.813 |
| `trip` | 0.803 | 1.656 | 2.338 |
| `oddrun` | 1.353 | 1.656 | 1.907 |
| `pattern` | 1.877 | 1.451 | 1.917 |
| `LAS` | 0.391 | 1.612 | 2.050 |
| `MAXP` | 1.343 | 1.842 | 1.971 |
| `KCHG` | 1.540 | 2.088 | 2.124 |
| `BAL>m/2` | 2.140 | 1.944 | 1.837 |
| `BAL>m/4` | 0.701 | 1.515 | 1.811 |
| `NBLK=1` | 2.624 | 1.450 | 1.899 |
| `NBLK>2` | 2.206 | 1.359 | 1.784 |
| `NBLK>3` | 1.067 | 1.508 | 1.761 |
| `AMAX` | 1.338 | 2.629 | 2.214 |
| `AMIN` | 1.330 | 2.261 | 2.218 |
| `longrun` | 1.989 | 1.876 | 2.076 |
| `pal` | 1.699 | 1.614 | 2.002 |
| `hashmod` | 1.987 | 1.987 | 1.987 |
| `monotone` | 1.000 | 1.000 | 1.000 |

## What the repair does

- **The bound holds, and the data now show it.** Every repeat removes one boundary cell; the region never widens again, so a run admits at most n repairs. The largest count seen on a single tape is 1956 (n=512: 485, n=1024: 1002, n=2048: 1956) - just under the tape length, never over it.
- **Termination is complete but the meet is no longer the typical outcome.** All 11146 cap-diers terminate and none reaches the cap, but only 4115 of them end in a meet: 7031 end with an empty region, because the repair eats the tape from the boundary the rule's action points at and the region is exhausted before the window closes. Under the retraction repair the walk merely narrowed a window it could re-widen, so it kept meeting; with a real trim the same bit mostly destroys the region.
- **Cost.** A repair pays a re-span, not one cell: `abs(y-L) + abs(x-R_)` plus the head re-read, because the window re-opens on the narrower tape. That is the same charge a rule-initiated cut pays. At n=512 the median all-tape reads are 3.29x the retraction run's (range 1.00 to 10.07 over the 20 rules), while the repair count itself is bounded: at most 485 repairs at n=512 and 1956 at n=2048, against the millions the retraction run needed on a single tape.
- **The side is the rule's own.** Over all 18000 runs the repair removed 1151794 left boundary cells and 1698863 right ones; the left share per rule runs from 0.00 (`pal`) to 0.97 (`runs`), so the repair really is read off the action table and not a fixed-side trim in disguise.
- **The cycle shape is untouched.** `cycle len` and `width` come from the BTHM-0 side, which the repair does not enter; they equal the values in `probe_cycle/cyc_cells.json` tape for tape (they are two of the seven fields in the exact comparison above).

## Files

`cut_cells.json` (this table), `cut_<rule>_<n>.json` (per-tape records: BTHM-0 fields, cut fields, repairs split left/right, hash slots used), `cut_engine.py` (mode 4), `cut_run.py`, `cut_all.py`, `cut_agg.py`, `cut_verify.py`, `cut_ctl.py` (the retraction control), `cut_pilot.py`. `probe_cycle/` is untouched.

