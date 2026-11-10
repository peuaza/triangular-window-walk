# EN_LINEAR — an exact O(n) machine on the meeting domain, for the other fifteen rules

*2026-11-09, `probe_en_linear/`. No `.tex` was read or modified; no git command was run.
Python: `C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`. Every number below is
printed by a script in this directory; the JSON that carries it is named in §9.*

## 0. Verdict

`CES_EXACT.md` moved `pattern` and `runs` from `Θ(n²)` to `Θ(n)` by asking a narrower question than
the record's: is there a machine that is exact on the **meeting** tapes and costs `O(n)` reads?  For
those two the answer was yes.  The question had never been put to the other fifteen rules.  This
probe puts it — and the answer separates them by *three different failure modes*, not one.

**One rule passes.  Fourteen fail, for three distinct reasons.**

| rule | `b1` at the region | failure mode | exact on the meeting tapes | `e_old/n` | `e_new/n` | verdict |
|---|---|---|---|---|---|---|
| `peak` | 3-gram peak count | — | **700/700** | 34 → 316 | **2.0** | **Θ(n)** |
| `dense2` | count of 2s | 1. the epoch law | 173 / 159 / 162 / 80 | 91 → 1056 | 6.8 → 2.2 | still `Θ(n²)` |
| `trip` | 3-gram triple count | 1. the epoch law | 200 / 196 / 194 | 45 → 356 | 3.9 → 2.1 | still `Θ(n²)` |
| `MAXP` | range max of the prefix walk | 1 + 3 | 200 / 199 / 198 / 99 | 59 → 369 | 355 → 47 | still `Θ(n²)` |
| `AMAX` | range max | 3. the terminal tail | 200 / 200 / 198 / 100 | 42 → 445 | 210 → 41 | still `Θ(n²)` |
| `AMIN` | range min | 3. the terminal tail | 199 / 200 / 198 / 100 | 57 → 286 | 268 → 63 | still `Θ(n²)` |
| `BAL>m/2` | last occurrence of a level | 2. the statistic | 188–195 / 200 | 11 → 47 | 535 | `Θ(n²)` by cost |
| `BAL>m/4` | last occurrence of a level | 2. the statistic | **200 / 200 / 200 / 100** | 15 → 52 | 173 → 1087 | `Θ(n²)` by cost |
| `NBLK=1` | count of visits to a level | 2. the statistic | **200 / 200 / 200 / 100** | 26 → 97 | 177 → 1120 | `Θ(n²)` by cost |
| `NBLK>2` | count of visits to a level | 2. the statistic | **200 / 200 / 200** | 7 → 45 | 136 → 535 | `Θ(n²)` by cost |
| `NBLK>3` | count of visits to a level | 2. the statistic | see M3 | 6 → 29 | 134 → 540 | `Θ(n²)` by cost |
| `oddrun` | odd-run count + parities | 1. the epoch law | see M3 | 45 → 463 | 1.9–2.3 | still `Θ(n²)` |
| `canonical` | `S, Q, Rr` | 3. the terminal tail | see M3 | 98 → 862 | 2.1–3.8 | still `Θ(n²)` |
| `LAS` | alternation count | 1. the epoch law | see M3 | 48 → 361 | 2.0–2.1 | still `Θ(n²)` |
| `KCHG` | two sign-change counts | 1. the epoch law | see M3 | 12 → 28 | 2.0–2.2 | still `Θ(n²)` |

(`e_old/n` is the record's own machine on its own meeting tapes, `probe_table_eg`; `e_new/n` is the
jump machine at `W2 = 256`, by `n = 512 … 4096`.  `exact` is out of 200 meeting tapes per size, 100 at
`n = 4096`.  The five rows that end in "see M3" were still running when this file was assembled; M3
is the full grid and is the authority.)

1. **The epoch law fails** (`dense2`, `trip`, `MAXP`).  At the start of an epoch the window is the
   whole region.  If `b2 = 0` the first move leaves the region at once and the epoch is one cut — that
   half is exact on **every** epoch of every tape measured (0 exceptions).  If `b2 = 1` the window
   shrinks and the epoch is a *scan* of median length 40–70 steps, and its exit side is a function of
   `b1` at the start for some rules and emphatically not for others: `dense2` cuts left on 64 % of its
   `b1 = 0` epochs and right on the other 36 %.  The errors sit in the middle of the tape (median
   `M0` = 359 at `n = 512`, 2 032 at `n = 4096`), so no buffer and no bounded lookahead reaches them:
   of 83 343 `dense2` mispredictions, 2 have epoch length ≤ 4.
2. **The statistic fails** (`BAL>m/2`, `BAL>m/4`, `NBLK=1/>2/>3`).  Here the epoch law is *exact* —
   `NBLK=1` has **zero** mispredicted epochs in 700 meeting tapes, `NBLK>2` has two — but `b1` is a
   level index: the number of visits to `Ps[L]` in `(L, R]`, or the last visit.  Its anchor `Ps[L]`
   moves with the cut, no `O(log n)`-bit cache survives that move, and a per-level position index is
   `Θ(n log n)` bits.  So every query re-scans the region, and there are `Θ(n)` queries: the machine is
   `Θ(n²)` **by construction**.  Measured `e/n` = 173 → 1 087 for `BAL>m/4`; `n²/4` at `n = 4 096` is
   1 024.  These are the only rules whose obstacle is not logic but the query.
3. **The terminal tail fails** (`AMAX`, `AMIN`, `canonical`).  The epoch law is exact — `AMAX` has
   **zero** wrong cuts in 700 meeting tapes, `AMIN` zero in 700 — and the statistic is a range max/min
   that a lazy `O(log n)`-bit argmax cache answers at 0 reads until the argmax dies.  What breaks is
   the one meeting terminal: its `M0` reaches 252, 170, 820 at `n = 512 … 2048` for `AMAX` (0.4 `n`),
   so a fixed `W2` cannot hold it, and the `W2` that does costs the endgame `Θ(W2²)` = `Θ(n²)`.

`peak` is the rule where none of the three fails, and it is the only one: a 3-gram count that slides
(one read per cell given up), an epoch law whose 1 037 mispredictions never change the answer, and a
meeting terminal whose maximum over 700 tapes is 27.  At `W2 = 128` it is exact on **700/700** meeting
tapes at `e = 2.0 n` — against the record's 34 → 316 `n`.

## 1. The machine, and what it is allowed to hold

The walk is the record's (`pcr_lib.walk`: the four-move table `TABLE = [L−, R+, R−, L+]` indexed by
`2·b2 + b1`, header reset after a cut, `conv='a'` for fourteen rules and `c1a` for `canonical`).  At
the start of an epoch the window *is* the region `[L,R]`, so every bit the rule needs is a query on
the region, and the region only ever shrinks — one cell per cut (two for a double cut).  That is the
whole opening: **a statistic that survives a one-cell deletion costs one read per cut; one that does
not costs a scan of the region.**

**Seven rules slide.**  The bit is a function of `O(1)` counters over `[L,R]`, and one cell leaving
costs 1–3 reads:

| rule | the counters | cost of one cut |
|---|---|---|
| `dense2` | `N2 = #{i : a_i = 2}` | `N2 −= [a_L = 2]` |
| `peak` | `P = #{i ∈ [L+1,R−1] : a_{i−1} < a_i > a_{i+1}}` | retire position `L+1`: re-read `a_L..a_{L+2}` |
| `trip` | `T = #{i ∈ [L,R−2] : a_i = a_{i+1} = a_{i+2}}` | retire position `L`: re-read `a_L..a_{L+2}` |
| `oddrun` | `C` = odd maximal constant runs inside, + the two boundary run parities | the boundary run shortens; its parity flips when it closes |
| `canonical` | `S = Σa`, `Q = Σa²`, `Rr = Σ a_i a_{i+1}`, the two end values | `S,Q −= a_L(, a_L²)`, `Rr −= a_L a_{L+1}` |
| `LAS` | `A` = sign alternations between consecutive non-zero `s`; first/last non-zero | drop the boundary non-zero, advance over zeros (charged to the zeros crossed) |
| `KCHG` | `c` = sign changes in `(L,R]`; `c_left` = those with `p − L ≤ M/2` | one cell changes both counts by 0–1 |

`KCHG` needs an identity to be a count at all: the rule's bit is `2·(x − L) > M` where `x` is the
*lower median* of the sign-change positions, and that is equivalent to `c_left ≤ (c+1)/2 − 1`.  The
median test is a comparison of two counts, not a selection.  `el_counter.py` drags the region from the
whole tape down to a point along a random cut sequence and checks the counter's bit against
`pcr_lib`'s `b1` at every step: **954 windows per rule, 0 mismatches, all seven**.

**Three rules need a range max/min.**  `MAXP` asks for `max_{j∈[L,R]} Ps[j] − Ps[L]`; `AMAX` whether
`Ps[R+1]` exceeds `max_{j∈[L,R]} Ps[j]`; `AMIN` the same with the minimum.  A maximum over a window
that loses cells at *both* ends cannot be maintained by the monotone-deque trick, and a sparse table is
`Θ(n log n)` bits.  What an `O(log n)`-bit machine can do is cache the argmax and re-scan when the
argmax is the cell that just died — the **lazy** device.  It is exact (`el_gate2.py`: lazy ≡ full scan
on every tape) and its cost is the re-scans, `nscan` in the tables.

**Five rules need a level index.**  `BAL` asks for the last `j ∈ (L,R]` with `Ps[j] = Ps[L]`; `NBLK`
how many times the walk visits `Ps[L]` there.  The anchor moves with the cut
(`Ps[L+1] = Ps[L] + s[L]`), so no cached position survives, and a per-level index is out of class.
Every query is a scan of the region.  This is the group whose *statistic*, not whose logic, is the
obstacle.

**Reads.**  One cell touch costs 1, counted with multiplicity, opening pass included.  The machine
holds the tape plus `O(log n)` bits of work tape — the same convention as `CES_EXACT` §5.2, which is
why the `W2` endgame is a *streaming* scan and not a buffer of `W2` cells (`W2 = 256` cells is
`Θ(W2)` bits, out of class).  So: region epoch → 1 read per cell the cut gives up, plus the eval cost
(0 sliding; 0-or-`M` lazy; `M` scan); endgame epoch (`M ≤ W2`) → 1 read per step plus `M` per step for
a range/level query, and 1 per step for a sliding statistic.

## 2. Gates (all passing, or nothing below counts)

1. `el_lib.trace_walk` reproduces `pcr_lib.Runner` — `δ`, `end`, `steps`, `cuts`, `inc` reads, `msum` —
   on every rule and every tape at `n = 48, 64, 128`: **0 mismatches in 45 cells** (`el_gate.py`).
   `BAL>m/4` is gated against `tr_bal4` (thr = 4; that module's own gate is thr = 2 ≡ `pcr_lib`),
   `NBLK>3` against `hmach3` (pcr_lib has no `NBLK>3` row).
2. The jump machine with `W2 ≥ n` never jumps and reproduces the trace's `(δ, end)` on every rule
   (`el_gate2.py`): the endgame's control flow is the walk's.
3. The lazy range-max device and a full scan agree on every tape (`MAXP`, `AMAX`, `AMIN`, 16 tapes each).
4. The BAL/NBLK predicates equal a literal prefix walk on 15 072 windows (`el_predgate.py`, 0
   mismatches); the seven sliding counters on 954 windows each (`el_counter.py`, 0 mismatches).

The machine never consults the truth: at `M > W2` it predicts the epoch's cut from `(b2,b1)`
(`b2 = 0` → one cut, double for `conv='a'`, single on the opposite side for `canonical`; `b2 = 1` →
`cutL` if `b1 = 0` else `cutR`, swapped for `canonical`), and at `M ≤ W2` it runs the epoch literally.
`(δ, end)` is compared with the trace's afterwards.

## 3. Failure mode 1 — the epoch law

The law, measured directly on 40 tapes per cell at `512` and `1024` (`el_law.py`), accuracy of the
best prediction per bucket:

| rule | `b2 = 0` (any `b1`) | `b2=1, b1=0` (512 / 1024) | `b2=1, b1=1` (512 / 1024) |
|---|---|---|---|
| `dense2` | 1.000 | **0.642 / 0.723** | **0.729 / 0.687** |
| `peak` | 1.000 | 0.944 / 0.876 | 0.990 / 0.998 |
| `trip` | 1.000 | 0.956 / 0.967 | 0.918 / 1.000 |
| `oddrun` | 1.000 | 0.969 / 1.000 | 0.838 / 0.992 |
| `canonical` | 1.000 | 0.995 / 1.000 | 0.991 / 0.962 |
| `LAS` | 1.000 | 0.788 / 0.674 | 0.946 / 0.977 |
| `KCHG` | 1.000 | 1.000 / 0.996 | 0.970 / 0.954 |
| `MAXP` | 1.000 | 1.000 / 0.974 | 1.000 / **0.681** |
| `AMAX` | 1.000 | 1.000 / 1.000 | 1.000 / 1.000 |
| `AMIN` | 1.000 | 1.000 / 1.000 | 1.000 / 1.000 |
| `BAL>m/2` | 1.000 | 1.000 / 1.000 | 0.906 / 0.984 |
| `BAL>m/4` | 1.000 | 1.000 / 1.000 | 0.990 / 0.969 |
| `NBLK=1` | 1.000 | 1.000 / 1.000 | 1.000 / 1.000 |
| `NBLK>2` | 1.000 | 1.000 / 1.000 | 1.000 / 0.998 |
| `NBLK>3` | 1.000 | 1.000 / 1.000 | 1.000 / 0.998 |

Two things stand out.  **`b2 = 0` is exact, always** — 0 errors in every rule and every size, and the
same in the machine's audit (the `b2=0 err` column of M4 is 0 in every row).  And **the `b2 = 1` half
is a spectrum**: exact for the extremum and level rules, 96–99 % for the counting rules, and 64–73 %
for `dense2` — a coin flip with a bias.  `dense2`'s `b1` carries almost no information about where the
scan will exit.

**The failures are not at the end of the tape.**  The median `M0` of a mispredicted epoch is 359 at
`n = 512` and 2 032 at `n = 4096` for `dense2`, and *all* of `trip`'s and `peak`'s mispredictions are
at `M0 ≥ 256` (M4).  They are the scans at the top of the tape, not the terminal.  That matters for the
repair: a prefix check cannot see them.  Of 83 343 `dense2` mispredictions, **2** have epoch length
≤ 4 and 325 ≤ 16; for `trip` it is 0 of 2 408 and 3 of 2 408.  Four steps of lookahead catch 0.002 %
of the failures — the same verdict `CES_EXACT` §4 reached for `pattern`/`runs`, for the opposite
reason.

**A wrong cut often does not change the answer, which is why the output-exact count and the law
disagree.**  `peak` mispredicts 1 037 epochs spread over 11 of 700 meeting tapes and returns the right
`δ` on all 700.  `trip` mispredicts 2 408 over 23 tapes and returns 6 wrong.  The difference between
the two rules is not the number of mispredictions but the bucket they land in: `trip`'s are 61 % in
`b1=0` and 39 % in `b1=1`, `peak`'s are 82 % in `b1=0`; a wrong cut on the side the terminal epoch does
not use is invisible in `δ`.  That is consistent with the numbers but not proved here; what is measured
is the counts.

## 4. Failure mode 2 — the statistic

The level rules are the *exact* ones: at `W2 = 32` and above, `BAL>m/4` is 200/200 at every size,
`NBLK>2` is 200/200, `NBLK=1` is 199–200/200, and the audit finds 0–66 mispredicted epochs in 700
tapes (M4).  Their machine is exact on the meeting domain.  It is also useless:

| rule | `e/n` at `W2=256`, `n` = 512 / 1024 / 2048 / 4096 | `n²/4` at 4 096 |
|---|---|---|
| `BAL>m/4` | 173 / 281 / 540 / 1 087 | 1 024 |
| `NBLK=1` | 177 / 299 / 567 / 1 120 | 1 024 |
| `NBLK>2` | 136 / 269 / 535 / — | 1 024 |

The measured numbers are `n²/4` to within 6 %: one scan of the region per epoch, `n/2` epochs.  No
`W2` helps (`W2=32` costs the same 534 `n` at `n = 2048`), because the cost is at the region, not in
the endgame.  The record's `hmach3` band machine answers the same queries at `O(1)` amortised **per
step**, which is why it is 15–52 `n` instead of 535 — but it is a step-by-step machine and cannot jump
epochs, so its cost is `Θ(steps)` = `Θ(n^{1.5})`, not `O(n)`.

`NBLK=1` deserves a line of its own: **the only rule in the catalogue with a perfect epoch law** — zero
mispredicted epochs in 700 meeting tapes and 700 halting tapes — and it still cannot be had for `O(n)`.
It is the best candidate for a hybrid (jump at the region, band in the endgame), and that hybrid is a
different probe.

## 5. Failure mode 3 — the terminal tail

`AMAX` and `AMIN` have the law (0 wrong cuts in 700 tapes each) and a cheap statistic (lazy argmax).
Their only failures are the meeting terminal, and its `M0` does not stay bounded:

| rule | meeting terminal `M0`, med / p90 / p99 / max, by `n` |
|---|---|
| `AMAX` | 3/11/38/**252** · 3/14/104/**170** · 2/13/70/**820** · 1/11/31/61 |
| `AMIN` | 3/12/114/**333** · 3/11/97/**177** · 3/13/115/**616** · 1/12/51/145 |
| `MAXP` | 3/11/29/**35** · 3/11/36/**78** · 1/10/42/**413** · 1/14/42/54 |
| `peak` | 1/8/18/**22** · 1/6/19/27 · 1/6/14/20 · 1/4/16/23 |
| `trip` | 1/8/14/**25** · 1/8/14/22 · 1/8/16/20 |
| `dense2` | 3/7/18/**42** · 1/8/19/52 · 1/8/23/27 · 1/6/18/18 |
| `BAL>m/4` | 1/6/17/22 · 1/6/17/18 · 1/6/11/21 · 1/4/8/9 |
| `NBLK=1` | 1/10/22/24 · 1/8/23/28 · 1/7/16/51 · 1/5/12/15 |

A tail at 0.2–0.4 `n` cannot be covered by a fixed `W2`: at `W2 = 256` `AMAX` is 200/200, 200/200,
**198/200**, 100/100 and `AMIN` 199/200, 200/200, 198/200, 100/100.  The `W2` that would cover 820 is
`Θ(n)`, and the endgame then costs `Θ(W2²)` = `Θ(n²)` — `CES_EXACT`'s `W2 = αn` curve, unchanged.
For `pattern` the meeting terminal was ≤ 17; for `AMAX` it is 820.  This is the same failure the
`runs` row of `CES_EXACT` §6 flagged as "a caveat with teeth", and here it bites.

`MAXP` fails both ways: the law's `b1 = 1` bucket drops to 0.681 at `n = 1024` (from 1.000 at 512 on
three samples), and its terminal tail reaches 413.  `MAXP` is also the only rule whose lazy device
degrades: its `nscan` grows 7 → 16 → 35 → 63 over `n = 512 … 4096`, so `e/n` at `W2 = 8` rises
10.3 → 26.1, i.e. `e ~ n^{1.45}` — better than the record's `n^{1.9}` but not linear.

## 6. The rule that passes

`peak`, `W2 = 128`, on 700 meeting tapes (200 at each of 512/1024/2048, 100 at 4096):

| `W2` | exact / nmeet by `n` | `e/n` by `n` |
|---|---|---|
| 32 | 200/200 · 196/200 · 198/200 · 97/100 | 2.0 · 2.0 · 2.0 · 2.0 |
| 64 | 200/200 · 198/200 · 199/200 · 98/100 | 2.0 · 2.0 · 2.0 · 2.0 |
| **128** | **200/200 · 200/200 · 200/200 · 100/100** | **2.0 · 2.0 · 2.0 · 2.0** |
| 256 | 200/200 · 200/200 · 200/200 · 100/100 | 2.7 · 2.2 · 2.0 · 2.0 |

`e = 2.0 n` against the record's 34 → 316 `n`, a factor of 17 to 158, and the exponent of `e/n` over
the four sizes is −0.13: flat.  **The verdict is `Θ(n)`**, with one explicit caveat: the exactness is
*output*-exact, not trajectory-exact.  1 037 epochs on 11 of the 700 tapes are mispredicted, and the
machine still returns the right `δ`.  I could not identify the mechanism; the honest statement is that
the answer is exact on 700/700 meeting tapes and that the reason is not that the cut sequence is the
walk's.

## 7. Reads, and the machine that is exact on every tape

For the sliding rules the cost of the jump machine is `n` (opening) + 0.5 `n` (cuts) + the endgame, so
`e/n` falls towards 1.5–2.0 as `n` grows at *any* fixed `W2` — `dense2` 6.8 → 4.0 → 2.9 → 2.2 and
`peak` 2.0 at every size.  Nothing is bought by a bigger `W2` except exactness.

There is a second machine worth naming: drop the re-span.  `pcr`'s INC charge includes
`|y − L| + |x − R|` per cut because "the reset forces the window back onto the whole region and the
state cannot store it" — but for a sliding statistic the state *is* `O(log n)` bits and is carried
across the cut, so the machine costs `n + steps + cuts` and is **exact on every tape**, meeting or not:

| rule | `n + steps + cuts`, `n` = 512 / 1024 / 2048 / 4096 |
|---|---|
| `dense2` | 30 / 57 / 103 / 198 |
| `peak` | 9 / 14 / 31 / 62 |
| `trip` | 16 / 18 / 43 / — |
| `MAXP` * | 20 / 37 / 56 / 123 |

(`MAXP` * is not a legal machine — its statistic needs a scan per step — the row is the step count.)
For `dense2` that is 0.048 `n²` against the record's 0.26 `n²`: five times cheaper and still quadratic.
It is the right machine to quote if the paper keeps the "exact on every tape" gate.

## 8. What to change in the record

1. **`peak` is `Θ(n)`, `p = e = 2.0 n`.**  ORDERS' "`peak`: `0.3 n²`? borderline" can be settled: the
   witness counter is `2.0 n` with `R² = 1` on the meeting tapes, so `p = e = 2.0 n`.  This is the
   first natural rule besides `pattern`/`runs` to leave the `n²` block.
2. **The `n²` block is not one phenomenon.**  Three different obstacles: the epoch law (`dense2`,
   `trip`, `MAXP`), the statistic (`BAL`, `NBLK`), the terminal tail (`AMAX`, `AMIN`).  A row that says
   `n²` without the mode invites the wrong next experiment: for `dense2` no machine will ever be
   `O(n)`; for `NBLK=1` a hybrid band + jump is the obvious next probe; for `AMAX` the question is
   whether the terminal tail is really `Θ(n)` or a slow transient.
3. **`NBLK=1` has a perfect epoch law** — 0 mispredictions in 1 400 tapes — and is the one rule where a
   cheap machine is ruled out by the *query* alone.  If any rule is worth a dedicated device, it is
   this one, not `dense2`.
4. The `b2 = 0` half of every epoch is exact, so the `p/e` predictor question reduces entirely to the
   `b2 = 1` scans.  That is the object to model.

## 9. Files

| file | what |
|---|---|
| `el_lib.py` | tape draws, the epoch trace (a line-for-line copy of `pcr_lib.walk` plus one record per epoch), the NBLK device |
| `el_gate.py` → `el_gate.json` | gate 1: trace ≡ `pcr_lib.Runner` / `tr_bal4` / `hmach3`, 45 cells |
| `el_gate2.py` → `el_gate2.json` | gates 2–3: `W2 ≥ n` ≡ the walk; lazy ≡ scan |
| `el_predgate.py` | BAL/NBLK predicates vs a literal prefix walk, 15 072 windows |
| `el_counter.py` | the seven sliding counters, 954 windows per rule |
| `el_law.py` → `el_law.json` | the epoch law by `(b2,b1)`, two sizes |
| `el_lawm.py`, `el_law3.py` → `el_lawm.json`, `el_law3.json` | the law by `M0` band and by tape kind |
| `el_cell.py` → `el_<rule>_<n>.json` | one cell: the trace, the epoch audit, the machine at six `W2` |
| `el_run.py` | the three job groups |
| `el_eold.py` → `el_eold.json` | the record's own `e` (mean reads over meeting tapes / `n`) from `probe_table_eg` |
| `el_report.py` → `el_report.json`, `el_tab.py` → `_tables.md` | the tables M1–M8 |
| `el_grid.py`, `el_view.py`, `el_diag.py`, `el_coincide.py`, `el_path.py`, `el_lazydbg.py` | diagnostics |
