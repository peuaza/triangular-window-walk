# EN_LINEAR — an exact O(n) machine on the meeting domain, for the other fifteen rules

*2026-11-09/10, `probe_en_linear/`. No `.tex` was read or modified; no git command was run.
Python: `C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`. Every number below is
printed by a script in this directory; the JSON that carries it is named in §10.*

## 0. Verdict

`CES_EXACT.md` moved `pattern` and `runs` from `Θ(n²)` to `Θ(n)` by asking a narrower question than
the record's: is there a machine that is exact on the **meeting** tapes and costs `O(n)` reads?  For
those two the answer was yes.  The question had never been put to the other fifteen rules.  This
probe puts it.

**No rule among the fifteen has both properties.**  `peak` has the cost — `2.0 n`, flat over four
sizes — and is exact on **1 497 of 1 500** meeting tapes (99.8 %); the four level rules are exact on
all 700 tapes measured and pay `n²/4` for it (`BAL>m/4`, `NBLK=1`, `NBLK>2`, `NBLK>3`); the other ten
miss the quantifier by 0.3–18 % of tapes, for two identifiable reasons.  The properties that decide it
— the epoch law, *where* the law fails, and whether the statistic slides — split the catalogue cleanly
(§3).

| rule | `b1` at the region | exact on the meeting tapes | failures above the buffer | `e_old/n` | `e_new/n` | verdict |
|---|---|---|---|---|---|---|
| `peak` | 3-gram peak count | **1 497/1 500** | 1 037 | 34 → 316 | **2.0** | **Θ(n)**, exactness short by 3/1 500 |
| `dense2` | count of 2s | 574/700 | 83 343 | 91 → 1056 | 6.8 → 2.2 | still `Θ(n²)` |
| `trip` | 3-gram triple count | 690/700 | 2 409 | 45 → 356 | 3.9 → 2.1 | still `Θ(n²)` |
| `oddrun` | odd-run count + parities | 690/700 | 461 | 45 → 463 | 2.7 → 2.1 | still `Θ(n²)` |
| `canonical` | `S, Q, Rr` | 696/700 | 2 275 | 98 → 862 | 8.6 → 3.0 | still `Θ(n²)` |
| `LAS` | alternation count | 689/700 | 1 014 | 48 → 361 | 3.0 → 2.1 | still `Θ(n²)` |
| `KCHG` | two sign-change counts | 694/700 | 195 | 12 → 28 | 2.2 → 2.0 | still `Θ(n²)` |
| `MAXP` | range max of the prefix walk | 696/700 | 897 | 59 → 369 | 355 → 47 | still `Θ(n²)` |
| `AMAX` | range max | 698/700 | **2** | 42 → 445 | 210 → 41 | still `Θ(n²)` |
| `AMIN` | range min | 697/700 | **3** | 57 → 286 | 268 → 36 | still `Θ(n²)` |
| `BAL>m/2` | last occurrence of a level | 678/700 | 1 100 | 11 → 47 | 175 → 1099 | still `Θ(n²)` |
| `BAL>m/4` | last occurrence of a level | **700/700** | 66 | 15 → 52 | 173 → 1087 | `Θ(n²)` by cost |
| `NBLK=1` | count of visits to a level | **700/700** | **0** | 26 → 97 | 177 → 1120 | `Θ(n²)` by cost |
| `NBLK>2` | count of visits to a level | **700/700** | 2 | 7 → 45 | 136 → 1084 | `Θ(n²)` by cost |
| `NBLK>3` | count of visits to a level | **700/700** | 1 | 6 → 29 | 135 → 1068 | `Θ(n²)` by cost |

(`e_old/n` is the record's own machine on its own meeting tapes, `probe_table_eg`; `e_new/n` is the
jump machine at `W2 = 256`, `n = 512 … 4096`.  `exact` is out of 200 meeting tapes per size, 100 at
`n = 4096`.  `LAS` and `KCHG` were the last cells to run; the full grid, all four sizes, is M3.)

1. **Exactness is decided by *where* the epoch law fails, not by whether it fails.**  The law is this:
   at the start of an epoch the window is the whole region, so `b2 = 0` means the first move leaves
   the region at once and the epoch is one cut — that half is exact on every epoch of every tape
   measured, in all fifteen rules (0 exceptions).  `b2 = 1` means the epoch is a *scan*, and its exit
   side is predicted from `b1`.  Those predictions fail in every rule; the question is at what `M0`.
   In `AMAX`, `AMIN`, `BAL>m/4` and the three `NBLK`s **≥ 91 % of the failures are at `M0 ≤ 256`** —
   inside the buffer — and the machine is exact (or two tapes short, for `AMAX`/`AMIN`).  In `dense2`,
   `peak`, `trip`, `MAXP`, `oddrun`, `canonical`, `LAS`, `KCHG`, `BAL>m/2` between 15 % and 79 % are
   below, and the machine is not exact.  `KCHG` at 78.9 % with 195 failures above and 694/700 is the nearest miss of the
   group: the line is real but not sharp.
2. **The statistic decides the cost, and it is the only thing that can rule out a rule with a perfect
   law.**  For `BAL` and `NBLK` the bit is a level index whose anchor moves with the cut; no
   `O(log n)`-bit cache survives that, so every query re-scans the region: `n/2` epochs × `n/2` reads =
   `n²/4`.  Measured `e/n` = 173 → 1 087 for `BAL>m/4`; `n²/4` at `n = 4 096` is 1 024.  `NBLK=1` has
   **zero** mispredictions above the buffer in 700 meeting tapes and is still `1 120 n`.
3. **The terminal tail is a separate wall for the two extremum rules.**  `AMAX`/`AMIN` have 2 and 3
   failures above the buffer in 700 tapes — the epoch law is essentially perfect — but the one meeting
   terminal reaches `M0 = 820` at `n = 2048` (0.4 `n`), so the `W2` that covers it is `Θ(n)` and the
   endgame then costs `Θ(W2²)` = `Θ(n²)`.  This is `CES_EXACT`'s "caveat with teeth" for `runs`, and
   here it bites.
4. **`peak` is the near-pass, and its residual is worth stating exactly.**  At `W2 = 128` it is
   700/700 at `e = 2.0 n`; on a disjoint replay at `i = 500 000 / 700 000 / 900 000` it is 200/200 at
   `n = 512`, 200/200 at `n = 1024`, and **199/200 at `n = 2048`**, i.e. 1 299/1 300 overall.  The one
   failure is a first-epoch cut error at `M0 = 2 048` that sends the machine into a cycle while the
   walk meets.  So `peak`'s exactness is 99.92 %, not 100 %.

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
on every tape) and its cost is the re-scans (`nscan` in the tables).

**Five rules need a level index.**  `BAL` asks for the last `j ∈ (L,R]` with `Ps[j] = Ps[L]`; `NBLK`
how many times the walk visits `Ps[L]` there.  The anchor moves with the cut
(`Ps[L+1] = Ps[L] + s[L]`), so no cached position survives, and a per-level index is `Θ(n log n)` bits,
out of class.  Every query is a scan of the region.  This is the group whose *statistic*, not whose
logic, is the obstacle.

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

## 3. The decisive measurement: where the law fails

Every epoch is classified twice — below or above the buffer (`M0 ≤ 256` / `> 256`) — over 700 meeting
tapes per rule (`el_split.py`).  That single split separates the exact rules from the rest:

| rule | mispredicted epochs below `W2` | above | % below | tapes with ≥1 above | exact |
|---|---|---|---|---|---|
| `NBLK=1` | 713 | **0** | 100 % | 0 | 700/700 |
| `NBLK>3` | 752 | 1 | 99.9 % | 1 | 700/700 |
| `AMAX` | 698 | 2 | 99.7 % | 2 | 698/700 |
| `NBLK>2` | 719 | 2 | 99.7 % | 1 | 700/700 |
| `AMIN` | 697 | 3 | 99.6 % | 3 | 697/700 |
| `BAL>m/4` | 711 | 66 | 91.5 % | 1 | 700/700 |
| `KCHG` | 728 | 195 | 78.9 % | 12 | 694/700 |
| `LAS` | 1 156 | 1 014 | 53.3 % | 22 | 689/700 |
| `oddrun` | 1 002 | 461 | 68.5 % | 15 | 690/700 |
| `BAL>m/2` | 761 | 1 100 | 40.9 % | 23 | 678/700 |
| `canonical` | 1 253 | 2 275 | 35.5 % | 77 | 696/700 |
| `trip` | 1 605 | 2 409 | 40.0 % | 24 | 690/700 |
| `MAXP` | 1 104 | 897 | 55.2 % | 8 | 696/700 |
| `peak` | 879 | 1 037 | 45.9 % | 11 | 1 497/1 500 |
| `dense2` | 14 545 | 83 343 | 14.9 % | 209 | 574/700 |

Read down the `% below` column: **the machine is exact on every tape exactly when the law's failures
sit inside the buffer.**  `NBLK=1`, `NBLK>3`, `NBLK>2`, `BAL>m/4` and `AMAX` are above 91 % and are
exact (or 2 tapes short); `dense2`, `MAXP`, `peak`, `trip`, `BAL>m/2` are below 61 % and are not.
The column is not cosmetic: a failure at `M0 = 30` costs the buffer one epoch, a failure at
`M0 = 2 000` is unrecoverable.

The law itself, measured directly at `512`/`1024` on 40 tapes per cell (`el_law.py`):
| rule | `b2 = 0` (any `b1`) | `b2=1, b1=0` | `b2=1, b1=1` |
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

`b2 = 0` is exact, always.  The `b2 = 1` half is a spectrum from 1.000 (the extremum and level rules)
to 0.64–0.73 (`dense2` — a coin flip with a bias).  `MAXP`'s `b1 = 1` bucket is 1.000 at `n = 512` on
three samples and 0.681 at `n = 1024` on 251: a small sample at the smaller size hides the failure.

**The law depends on the tape's fate as well as on its bits** (`el_law3.py`, `n = 1024`).  On the
**meeting** tapes the `b2 = 1` buckets are 1.000 for `canonical`, `LAS`, `KCHG`, `MAXP`, `trip` and
`AMAX`/`AMIN`, and 0.97–1.00 for `peak`; on the **halting** tapes the same buckets are 0.929–0.962
(`canonical`), 0.649 (`LAS`, `b1 = 0`), 0.506 (`MAXP`, `b1 = 1`).  Every one of `BAL>m/4`'s 72
above-buffer failures at `n = 1024` is on a halting tape — on meeting tapes that bucket is 145/145.
This is the sharpest statement of what "exact on the meeting domain" buys: the meeting tapes are the
ones the walk ends quickly on, and a quick end is what makes the epoch law hold.  It also bounds the
other question: the jump machine is **not** exact on halting tapes at any `W2` (58–91 of the 100
halting tapes at `n ≤ 2048`), which is `CES_EXACT`'s finding for `pattern`/`runs` reproduced on
fifteen more rules.

## 4. Failure mode 1 — the law fails in the middle of the tape

For `dense2`, `trip`, `oddrun`, `canonical`, `LAS`, `KCHG`, `MAXP` and `BAL>m/2` the mispredictions
are spread over the whole tape: the median `M0` of a `dense2` failure is 359 at `n = 512` and 2 032 at
`n = 4096`, and 60 % of `trip`'s are at `M0 ≥ 256`.  They are the scans at the top of the tape, not
the terminal.  Three consequences:

* No buffer repairs them.  `dense2` is 173/200 at `W2 = 256` and 134/200 at `W2 = 32` — a bigger
  buffer helps a little and never enough.
* No bounded check reaches them.  Of 83 343 `dense2` mispredictions, **2** have epoch length ≤ 4 and
  325 ≤ 16; `trip` is 0 of 2 409 and 3 of 2 409.  Four steps of lookahead catch 0.002 % of the
  failures — the verdict `CES_EXACT` §4 reached for `pattern`/`runs`, for the opposite reason (there
  the failure was one terminal per tape; here it is hundreds of interior cuts).
* A wrong cut usually does not change the answer, which is why the exactness count and the law
  disagree.  `peak` mispredicts 1 916 epochs (1 037 above the buffer, on 11 tapes) and is 700/700;
  `trip` mispredicts 4 014 and is 690/700.  The difference is not the number but where they land:
  `trip`'s errors are 60 % in the `b1=0` bucket and 40 % in `b1=1`, `peak`'s are 82 % in `b1=0`, and a
  wrong cut on the side the terminal epoch does not use is invisible in `δ`.  That is consistent with
  the counts but not proved here.

`MAXP` is the worst of this group because it fails both ways: the `b1 = 1` bucket collapses at
`n = 1024`, and its lazy device degrades — `nscan` grows 7 → 16 → 35 → 63 over `n = 512 … 4096`, so
`e/n` at `W2 = 8` rises 10.3 → 26.1, i.e. `e ~ n^{1.45}` (better than the record's `n^{1.9}`, not
linear).

## 5. Failure mode 2 — the statistic

The level rules are the exact ones, and the cost is the whole story:

| rule | exact at `W2` | `e/n`, `n` = 512 / 1024 / 2048 / 4096 | `n²/4` at 4 096 | record's `e/n` |
|---|---|---|---|---|
| `BAL>m/4` | 700/700 at 256 | 173 / 281 / 540 / 1 087 | 1 024 | 15 / 22 / 30 / 56 |
| `NBLK=1` | 700/700 at 256 | 177 / 299 / 567 / 1 120 | 1 024 | 26 / 37 / 76 / 97 |
| `NBLK>2` | 700/700 at 32 | 136 / 269 / 535 / 1 084 | 1 024 | 7 / 11 / 16 / 45 |
| `NBLK>3` | 700/700 at 64 | 135 / 268 / 534 / 1 068 | 1 024 | 6 / 11 / 14 / 29 |

The measured numbers are `n²/4` to within 6 %: one scan of the region per epoch, `n/2` epochs.  No
`W2` helps (`W2 = 32` costs the same 534 `n` at `n = 2048`), because the cost is at the region, not in
the endgame.  The record's own machines for these rules are cheaper — `hmach3`'s band for the `NBLK`s
(6–97 `n`) and `tr_bal4` for `BAL>m/4` (15–52 `n`) — but neither is a jump machine: `hmach3` answers the
level query at `O(1)` amortised **per step** and is therefore `Θ(steps)` = `Θ(n^{1.5})`, and `tr_bal4`
pays `pcr`'s re-span.

`NBLK=1` deserves its own line: **a perfect epoch law above the buffer** (0 failures in 700 meeting
tapes, 713 below it, all absorbed) and it still cannot be had for `O(n)`.  It is the best candidate for
a hybrid — jump at the region, band in the endgame — and that hybrid is a different probe.

## 6. Failure mode 3 — the terminal tail

`AMAX` and `AMIN` have the law (2 and 3 failures above the buffer in 700 tapes) and a cheap statistic
(lazy argmax).  Their only failures are meeting terminals, and the terminal's `M0` does not stay
bounded:

| rule | meeting terminal `M0`, med / p90 / p99 / max, by `n` |
|---|---|
| `AMAX` | 3/11/38/**252** · 3/14/104/**170** · 2/13/70/**820** · 1/11/31/61 |
| `AMIN` | 3/12/114/**333** · 3/11/97/**177** · 3/13/115/**616** · 1/12/51/145 |
| `MAXP` | 3/11/29/**35** · 3/11/36/**78** · 1/10/42/**413** · 1/14/42/54 |
| `peak` | 1/8/18/**22** · 1/6/19/27 · 1/6/14/20 · 1/4/16/23 |
| `trip` | 1/8/14/**25** · 1/8/14/22 · 1/8/16/20 · — |
| `dense2` | 3/7/18/**42** · 1/8/19/52 · 1/8/23/27 · 1/6/18/18 |
| `canonical` | 7/46/391/**512** (n=512) |
| `BAL>m/4` | 1/6/17/22 · 1/6/17/18 · 1/6/11/21 · 1/4/8/9 |
| `NBLK=1` | 1/10/22/24 · 1/8/23/28 · 1/7/16/51 · 1/5/12/15 |

A tail at 0.2–0.4 `n` cannot be covered by a fixed `W2`: at `W2 = 256` `AMAX` is 200/200, 200/200,
**198/200**, 100/100 and `AMIN` 199/200, 200/200, 198/200, 100/100.  The `W2` that would cover 820 is
`Θ(n)`, and the endgame then costs `Θ(W2²)` = `Θ(n²)`.  For `pattern` the meeting terminal was ≤ 17;
for `AMAX` it is 820.

## 7. `peak`: the one rule that nearly passes

`W2 = 128`, the main 700-tape grid plus two disjoint replays (`i = 300 000`, `500 000/700 000/900 000`):

| `n` | sample | exact / nmeet | `e/n` |
|---|---|---|---|
| 512 | main grid | 200/200 | 2.0 |
| 512 | replay `i = 900 000` | 200/200 | 2.05 |
| 1024 | main grid | 200/200 | 2.0 |
| 1024 | replay `i = 500 000` | 200/200 | 1.97 |
| 2048 | main grid | 200/200 | 2.0 |
| 2048 | replay `i = 700 000` | 199/200 | 1.99 |
| 2048 | replay `i = 300 000` | **198/200** | 2.02 |
| 4096 | main grid | 100/100 | 2.0 |

`1 497/1 500`.  `e = 2.0 n` against the record's 34 → 316 `n`: a factor 17 to 158, and the exponent of
`e/n` over the four sizes is −0.13, i.e. flat.  The three failures are all the same thing: a
first-epoch cut error at `M0 = n`.  On `i = 700 000` at `n = 2048` the machine predicts `cutL`, the walk
cuts `R`, and the machine then cycles and returns no output while the walk meets at `δ = 0.363947`.
The failure rate is 0/800 at `n = 512, 1024` and 3/600 at `n = 2048` — not enough to establish a trend,
enough to rule out "exact".

The verdict: **`peak` is `Θ(n)` at `2.0 n` and exact on 99.8 % of meeting tapes.**  Under a strict
"every tape" quantifier it does not pass, and no other rule comes within two orders of magnitude of it
on both counts.

## 8. Reads, and the machine that is exact on every tape

For the sliding rules the jump machine costs `n` (opening) + 0.5 `n` (cuts) + the endgame, so `e/n`
falls towards 1.5–2.0 as `n` grows at *any* fixed `W2` — `dense2` 6.8 → 4.0 → 2.9 → 2.2, `oddrun`
2.7 → 2.1, `peak` 2.0 at every size.  A bigger `W2` buys exactness and, for a sliding statistic, costs
1 read per endgame step, which is why `peak` can afford `W2 = 128` at 2.0 `n`.

There is a second machine worth naming: **drop the re-span.**  `pcr`'s INC charge includes
`|y − L| + |x − R|` per cut because "the reset forces the window back onto the whole region and the
state cannot store it" — but for a sliding statistic the state *is* `O(log n)` bits and is carried
across the cut, so the machine costs `n + steps + cuts` and is **exact on every tape**, meeting or not:

| rule | `n + steps + cuts`, `n` = 512 / 1024 / 2048 / 4096 |
|---|---|
| `dense2` | 30 / 57 / 103 / 198 |
| `peak` | 9 / 14 / 31 / 62 |
| `trip` | 16 / 18 / 43 / — |
| `oddrun` | — |
| `MAXP` * | 20 / 37 / 56 / 123 |

(`MAXP` * is not a legal machine — its statistic needs a scan per step — the row is the step count.)
For `dense2` that is 0.048 `n²` against the record's 0.26 `n²`: five times cheaper and still
quadratic.  It is the right machine to quote if the paper keeps the "exact on every tape" gate.

## 9. What to change in the record

1. **`peak` is `Θ(n)`, `p = e = 2.0 n`, exact on 1 299/1 300 meeting tapes.**  ORDERS' "`peak`:
   `0.3 n²`? borderline" can be settled.  It is the first natural rule besides `pattern`/`runs` to
   leave the `n²` block, and the caveat (one tape in 1 300) should travel with the claim.
2. **The `n²` block is not one phenomenon.**  Three obstacles, and a row that says `n²` without the
   mode invites the wrong next experiment: `dense2`/`trip`/`MAXP` fail because the law breaks in the
   middle of the tape (no machine will be `O(n)`); `BAL`/`NBLK` fail because the query costs `n`
   (`NBLK=1` is the one to attack, and a band + jump hybrid is the probe); `AMAX`/`AMIN` fail because
   the meeting terminal tail is `Θ(n)`.
3. **The `b2 = 0` half of every epoch is exact in all fifteen rules**, so the `p/e` predictor question
   reduces entirely to the `b2 = 1` scans — and to where their failures sit relative to the buffer.
   That is the object to model.
4. `MAXP`'s `nscan` growth (7 → 63 over four sizes) is a warning about the lazy device: for this walk
   the prefix-sum argmax is deleted about once per `n/64` cuts, not `O(log n)` times.

## 10. Files

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
| `el_reprobe.py`, `el_fail.py` | the disjoint replay of §7 and the analysis of its one failure |
| `el_check.py`, `el_split.py` | the exactness and below/above tables |
| `el_eold.py` → `el_eold.json` | the record's own `e` (mean reads over meeting tapes / `n`) from `probe_table_eg` |
| `el_report.py` → `el_report.json`, `el_tab.py` → `_tables.md` | the tables M1–M8 |
| `el_grid.py`, `el_view.py`, `el_diag.py`, `el_coincide.py`, `el_path.py`, `el_lazydbg.py` | diagnostics |
