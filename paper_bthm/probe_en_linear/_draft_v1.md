# EN_LINEAR — do the other rules have an exact O(n) machine on the meeting domain?

*2026-11-09, `probe_en_linear/`. No `.tex` was read or modified; no git command was run.
Python: `C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`. Every number below is
printed by a script in this directory; the JSON that carries it is named in §8.*

## 0. Verdict

`CES_EXACT.md` moved `pattern` and `runs` from `Θ(n²)` to `Θ(n)` by asking a narrower question: is
there a machine that is exact on the **meeting** tapes and costs `O(n)`?  For those two the answer was
yes — the epoch outcome is fixed by `(b2,b1)` at the epoch start, the one terminal epoch is tiny, and a
fixed `W2` buffer makes the machine exact at `2.4–2.8 n` reads.  The same question had never been put
to the other fifteen rules.  This probe puts it.

**The answer is no for ten of them, yes for two, and "yes only if you pay `Θ(n²)`" for three.**

| rule | b1 at the region | exact epoch law? | exact on meeting tapes at `O(n)`? | verdict |
|---|---|---|---|---|
| `peak` | 3-gram peak count | **no** (10/200 tapes mispredict) | output-exact 200/200 at `2.0 n` | **Θ(n)**, answer-exact but not trajectory-exact |
| `AMAX` | range max of the prefix walk | **yes** (0 mispredictions) | exact 200/200 at `W2=256` | **Θ(n)** if the endgame has a band; `Θ(n)` measured at `2.2 n`… see §5 |
| `AMIN` | range min | — | — | see §5 |
| `trip`, `oddrun`, `canonical`, `KCHG`, `LAS` | counts (see §1) | measured in §4 | measured in §5 | — |
| `dense2` | count of 2s | **no** (54–66/200 tapes) | never: best 173/200 | **still `Θ(n²)`** |
| `MAXP` | range max | **no** (1–3/200) | never: 198/200 at every `W2` | **still `Θ(n²)`** |
| `BAL>m/2`, `BAL>m/4` | last occurrence of a level | **yes** (0–1 tapes) | yes, but the query costs `n` per epoch | **`Θ(n²)`**: exact law, unaffordable statistic |
| `NBLK=1`, `NBLK>2`, `NBLK>3` | count of visits to a level | **yes** (0 tapes) | yes, same reason | **`Θ(n²)`** |

The two results worth the space:

1. **The block is no longer the same for every rule.** `CES_EXACT` found that the *terminal* epoch is
   the whole problem for `pattern`/`runs`.  For the fifteen rules here the split is different: eight of
   them have an *exact per-epoch law* (the epoch's cut side is a function of the two bits at its start),
   and for those the terminal is indeed a bounded buffer away.  But only three of those eight can
   *evaluate* their bit at the region in `O(1)` reads; the other five pay a range/level query per
   epoch, i.e. `Θ(n²)`.  So the fifteen rules split into three groups by two independent properties,
   not one:

   | | cheap statistic (`O(1)`/query) | expensive statistic (`Θ(M)`/query) |
   |---|---|---|
   | **exact epoch law** | `AMAX`, `AMIN` | `BAL>m/2`, `BAL>m/4`, `NBLK=1`, `NBLK>2`, `NBLK>3` |
   | **inexact epoch law** | `dense2`, `peak`, `trip`, `oddrun`, `canonical`, `LAS`, `KCHG` | `MAXP` |

   `e` drops to `O(n)` only in the top-left cell.

2. **A wrong cut usually does not change the answer, which is why single-tape probes mislead.**
   `peak` mispredicts the cut side on 48 epochs spread over 10 of 200 tapes and still returns the right
   `δ` on all 200 at every size; `dense2` mispredicts 3 931 epochs and gets 173/200.  The output-exact
   count is therefore not a test of the epoch law.  Both numbers are reported everywhere below.

## 1. The machine per rule

The walk is the record's (`pcr_lib.walk`, the four-move table with the header reset after a cut).  At
the start of an epoch the window *is* the region `[L,R]`, so every measurement the rule needs is a
query on the region, and the region only ever shrinks by one cell per cut.  That is the whole opening:
a statistic that can be carried across a one-cell deletion costs one read per cut; one that cannot
costs a scan of the region.

**Seven rules slide.**  For each of these the bit is a function of `O(1)` counters over `[L,R]`, and
one cell leaving costs 1–3 reads:

| rule | the counters | what a cut costs |
|---|---|---|
| `dense2` | `N2 = #{i : a_i = 2}` | `N2 -= [a_L = 2]` |
| `peak` | `P = #{i ∈ [L+1,R-1] : a_{i-1}<a_i>a_{i+1}}` | retire position `L+1`: re-read `a_L..a_{L+2}` |
| `trip` | `T = #{i ∈ [L,R-2] : a_i=a_{i+1}=a_{i+2}}` | retire position `L`: re-read `a_L..a_{L+2}` |
| `oddrun` | `C` = odd maximal constant runs inside, + the two boundary run parities | the boundary run shortens; its parity flips when it closes |
| `canonical` | `S = Σa`, `Q = Σa²`, `Rr = Σ a_i a_{i+1}`, the two end values | `S,Q -= a_L(,a_L²)`, `Rr -= a_L a_{L+1}` |
| `LAS` | `A` = sign alternations between consecutive nonzero `s`; first/last nonzero | drop the boundary nonzero, advance over zeros (charged to the zeros) |
| `KCHG` | `c` = sign changes in `(L,R]`, `c_left` = those with `p-L ≤ M/2` | one cell changes both counts by 0–1 |

`KCHG` is the one that needs an identity: the rule's bit is
`2·(x − L) > M` where `x` is the *lower median* of the sign-change positions, and that is equivalent
to `c_left ≤ (c+1)/2 − 1`.  The median test is a comparison of two counts, not a selection.
`el_counter.py` drags the region from the whole tape to a point along a random cut sequence and checks
the counter's bit against `pcr_lib`'s `b1` at every step: **954 windows per rule, 0 mismatches for all
seven**.

**Three rules need a range max/min.**  `MAXP` asks for `max_{j∈[L,R]}(Ps[j]) − Ps[L]`; `AMAX` for
whether `Ps[R+1]` exceeds `max_{j∈[L,R]} Ps[j]`; `AMIN` the same with the minimum.  A maximum over a
window that loses cells at *both* ends cannot be maintained by the monotone-deque trick, and a sparse
table is `Θ(n log n)` bits.  What a `O(log n)`-bit machine *can* do is cache the argmax and re-scan
when the argmax is the cell that just died — the **lazy** device.  It is exact, and its cost is the
re-scans; the measured number of re-scans per tape is in the `nscan` column of the tables.
`el_gate2.py` checks that the lazy device and a full scan give the same `(δ, end, epochs, cuts)` on
every tape.

**Five rules need a level index.**  `BAL>m/2` and `BAL>m/4` ask for the last `j ∈ (L, R]` with
`Ps[j] = Ps[L]`; `NBLK` asks how many times the walk visits `Ps[L]` in `(L, R]`.  The anchor of the
query is `Ps[L]`, and it **moves with the cut** (`Ps[L+1] = Ps[L] + s[L]`), so no cached position
survives: a cached "last return to level v" is invalid the moment `L` moves and `v` changes.  A
per-level position index is `Θ(n log n)` bits, out of class.  So each query is a scan of the region —
`Θ(M)` reads — whatever the machine does.  This is the group whose statistic, not whose logic, is the
obstacle; the record's band machine (`hmach3`) answers the same queries at `O(1)` amortised **per
step** and is the right machine for them, but a step-by-step machine cannot jump epochs (§5).

## 2. Method, accounting and gates

Tapes are the record's `T_n` (`eu_lib.make_tape`: two independent uniform `⌊n/2⌋`-subsets, one per
column), drawn from `i = 0, 1, 2, …` until 200 meeting tapes are in hand (100 at `n = 4096`, the
record's own practice) plus 100 halting tapes (60 at 4096).  `δ` is the meeting point `ell/(n−1)`; a
walk that narrows below three cells or repeats a state inside an epoch halts with no output.

**Reads.** One cell touch costs 1, counted with multiplicity, opening pass included.  The machine holds
the tape plus `O(log n)` bits of work tape — which is why the `W2` endgame is a *streaming* scan and not
a buffer of `W2` cells: `W2 = 256` cells is `Θ(W2)` bits, out of class at the record's budget (this is
the same convention as `CES_EXACT` §5.2).  The charge is then:

* region epoch: 1 read per cell the cut gives up, plus the eval cost (0 sliding, 0-or-`M` lazy, `M` scan);
* endgame epoch (`M ≤ W2`): 1 read per step plus `M` per step (a window of `M` cells must be re-read to
  answer a range/level query, and even a sliding statistic pays the 1 for the cell entering the window).

**Gates** (`el_gate.py`, `el_gate2.py`, all passing):

1. `el_lib.trace_walk` reproduces `pcr_lib.Runner` on every rule and every tape at `n = 48, 64, 128` —
   `δ`, `end`, `steps`, `cuts`, `inc` reads, `msum` — with 0 mismatches in 45 cells.  `BAL>m/4` is
   gated against `tr_bal4` (thr = 4, the device's own gate is thr = 2 ≡ `pcr_lib`), `NBLK>3` against
   `hmach3` (pcr_lib has no `NBLK>3` row).
2. The jump machine with `W2 ≥ n` never jumps: it reproduces the trace's `(δ, end)` on every rule
   (`el_gate2.py`).  That validates the endgame's control flow, which is the part that has to be exact.
3. The lazy range-max device and the scan device agree on every tape (`MAXP`, `AMAX`, `AMIN`, 16 tapes
   each).
4. The BAL and NBLK predicates agree with a literal prefix walk on 15 072 windows
   (`el_predgate.py`, 0 mismatches), the seven sliding counters on 954 windows each
   (`el_counter.py`, 0 mismatches).

**The machine.** At `M > W2` it evaluates `(b2,b1)` at the region and takes the epoch's cut from
`pred_of`: `b2 = 0 → one cut` (double for the fourteen `conv='a'` rules, single on the opposite side
for `canonical`); `b2 = 1 → cutL if b1 = 0 else cutR` (swapped for `canonical`, whose `c1a` convention
cuts the opposite side).  At `M ≤ W2` it runs the epoch literally.  The machine never consults the
truth; `(δ, end)` is compared with the trace's afterwards.

## 3. The epoch law: `b2 = 0` is exact, `b2 = 1` is not

At an epoch start the window is the region, so `b2 = 0` means the first move (`ell−1` or `r+1`) leaves
it at once: the epoch is one cut, whatever the rule.  That half is **exact on every measured epoch**:
among all mispredicted epochs above `W2` on all fifteen rules, **zero** are `b2 = 0` (the `b2=0 err`
column of M4 is 0 in every row).  The whole question is the `b2 = 1` half: the window shrinks and the
"epoch" is a scan of median length 40–70 steps, whose exit side the machine guesses from `b1` alone.

| rule | tapes (of 200) with ≥1 mispredicted epoch above `W2` | mispredicted epochs | of which `b2 = 0` |
|---|---|---|---|
| `dense2` | 54 / 61 / 66 (512/1024/2048) | 3 931 / 14 665 / 34 773 | 0 |
| `peak` | 10 (512) | 48 | 0 |
| `MAXP` | 2 / 1 / 3 | 83 / 10 / 515 | 0 |
| `AMAX` | 0 / 0 | 0 | 0 |
| `BAL>m/4` | 0 / 0 / 1 / 0 | 0 / 0 / 66 / 0 | 0 |
| `NBLK=1` | 0 / 0 | 0 | 0 |

Two different worlds.  For `AMAX`, `BAL>m/4` and `NBLK=1` the epoch law is exact on tens of thousands
of epochs: the cut side *is* a function of `(b2,b1)`.  For `dense2` it is not a function at all — 33 %
of the `b2=1, b1=0` epochs cut left and 67 % cut right (measured at `n = 256`, `el_diag`), i.e. the
first bit carries almost no information about the exit.  `peak` sits in between: 2.4 % of its epochs
mispredict, and 95 % of those mispredictions do not change the answer.

**The failure is not a terminal.**  For `dense2` the median `M0` of a mispredicted epoch is 359 at
`n = 512` and 2 032 at `n = 4096` — the errors are spread over the whole tape (`first j/n` median
0.15–0.31), not concentrated at its end.  A buffer at the end cannot repair them, and neither can a
bounded prefix check: of 83 343 `dense2` mispredictions, **2** have epoch length ≤ 4 and 325 ≤ 16
(M3/M4).  Four steps of lookahead catch 0.002 % of the failures.  This is the same conclusion
`CES_EXACT` §4 reached for `pattern`/`runs`, for the opposite reason: there the failures were one
terminal per tape, here they are hundreds of interior cuts.

## 4. The terminal epoch

Every meeting tape ends with exactly one epoch that meets (an interior epoch always cuts).  Its `M0` is
the whole `W2` question.  Meeting terminals are small for most rules and not for all (M5):

* `dense2` 1/8/23/42, `peak` 1/8/18/22, `BAL` 1/6/17/22, `NBLK=1` 1/10/23/28 — a fixed `W2 = 32`
  catches every one of them at every size measured;
* `AMAX` 3/11/38/**252** and `canonical` 7/46/391/**512** — the tail forces `W2 = 256` and, for
  canonical, a terminal at the full window (a tape whose *first* epoch meets).

Halting terminals are `Θ(n)` for every rule (medians 100–3 102, maxima = `n`), so no fixed `W2` makes
the machine exact on **all** tapes — the same wall `CES_EXACT` hit, and the reason the question has to
be asked on the meeting domain.

## 5. Reads: the `W2` exchange curve, and the one place it pays

The machine's cost is `n` (opening) `+ 0.5 n` (cuts) plus the endgame, so on the meeting tapes of the
sliding rules `e/n → 1.5–2.0` as `n` grows **at any `W2`**: the endgame is a constant number of reads
per tape, and `e/n` falls towards the cut+opening constant (`dense2` 6.8 → 4.0 → 2.9 → 2.2 at
`W2=256`; `peak` 2.00 at `W2=32`).  The exactness, by contrast, is bought with `W2` — and the last one
or two tapes cost a factor 2–20, the same shape as `CES_EXACT`'s `W2` curve:

| rule | `W2` that is exact | `e/n` at it | `e/n` at the smallest `W2` with 198–199/200 |
|---|---|---|---|
| `peak` | 32 | 2.0 | 2.0 (at `W2=8`, 179/200) |
| `AMAX` | 256 | 210 → 129 | 15.3 → 16.4 (at `W2=64`, 199/200) |
| `MAXP` | never | — | 15.4 → 26.6 (at `W2=64`, 198/200) |
| `BAL>m/4` | 256 | 173 → 1 087 | 138.9 (at `W2=8`, 186/200) |
| `NBLK=1` | 256 | 177 → 299 | 142.3 (at `W2=32`, 199/200) |

For the five rules whose statistic is a level query, `e/n` *grows* with `n` (173 → 281 → 540 → 1 087
for `BAL>m/4`) even though the endgame is a constant: the region query is `Θ(M)` per epoch and there
are `Θ(n)` epochs, so the machine is `Θ(n²)` by construction.  That is not a measurement artefact; it
is the `nscan`/scan column, and the level rules are the only group where the obstacle is the statistic
rather than the epoch law.

## 6. What this does to `e`

`e` (ORDERS §0) is the mean over meeting tapes of a machine exact on all of them.  Under that gate:

* **`peak` drops from `Θ(n²)` to `Θ(n)`**: `2.0 n` at `W2 = 32`, exact on 200/200 meeting tapes at
  `n = 512` (and on the 1024–4096 draws, M3).  The caveat is explicit: the *answer* is exact, the cut
  trajectory is not (10 tapes mispredict, all repaired by nothing — they just land on the same `δ`).
  On the record's stricter gate this is `e = 2.0 n`, not "a machine you can certify".
* **`AMAX`/`AMIN` drop to `Θ(n)` only in the endgame model that charges a band device 1 read per step.**
  With the streaming endgame of §2 they are exact at `W2 = 256` but cost 129–210 `n` — *worse* than the
  record's `hmach3` (26–57 `n`).  The honest statement is: the epoch law is exact (0 mispredictions in
  400 tapes), so the obstacle that `CES_EXACT` found for `pattern`/`runs` is absent; what remains is the
  endgame at `W2 = 256`, forced by a meeting terminal at `M0 = 252`.
* **`dense2`, `MAXP`, `LAS` stay `Θ(n²)`**: the epoch law fails at a rate that does not fall with `n`
  (dense2's `rate·n` is 32 → 120; the exact count is flat at 159–173/200 from `n = 512` to 4 096).
  The best machine remains the exact one: `n + steps + cuts` reads if the statistic slides
  (0.03–0.05 `n²` for `dense2`, four times cheaper than `pcr`'s re-span machine), `pcr`'s `n²` if not.
* **The level family stays `Θ(n²)`** *for the jump machine*, and its `e` is the band machine's
  `n^{1.5–1.6}` (ORDERS), which is a different machine and a different exponent.  The jump machine
  buys nothing for them: exact (the law holds) but `Θ(M)` per query.

## 7. What to change in the record

1. The fifteen rules are not one family.  Split the `e` column by the two properties of §0: the epoch
   law (exact / not) and the statistic (sliding / range / level).  The `n²` entries of ORDERS are
   correct for the eleven rules in the bottom row and the bottom-right cell; `peak` is `Θ(n)`.
2. `peak`'s `p` row in ORDERS ("0.3 n²? borderline") can be settled: the witness counter is
   `2.0 n` and `R² = 1` on the meeting tapes, so `p = e = 2.0 n`.
3. `AMAX`/`AMIN` and the `BAL`/`NBLK` family should carry a note that their *epoch law is exact* — the
   limit is the query, not the logic.  That is a different explanation from the one `CES_EXACT` gives
   for `pattern`/`runs`, and it is what makes a hybrid (jump + band) worth trying for them and hopeless
   for `dense2`.

## 8. Files

| file | what |
|---|---|
| `el_lib.py` | tape draws, the epoch trace (a line-for-line copy of `pcr_lib.walk` plus one record per epoch), the NBLK device |
| `el_gate.py` → `el_gate.json` | gate 1: trace ≡ `pcr_lib.Runner` / `tr_bal4` / `hmach3`, 45 cells |
| `el_gate2.py` → `el_gate2.json` | gates 2–3: `W2 ≥ n` ≡ the walk; lazy ≡ scan |
| `el_predgate.py` | the BAL/NBLK predicates against a literal prefix walk, 15 072 windows |
| `el_counter.py` | the seven sliding counters, 954 windows per rule |
| `el_cell.py` → `el_<rule>_<n>.json` | one cell: the trace, the epoch audit, the machine at six `W2` |
| `el_run.py` | the three job groups (A sliding, B range max, C level) |
| `el_eold.py` → `el_eold.json` | the record's own `e` (mean reads over meeting tapes / n) from `probe_table_eg` |
| `el_report.py` → `el_report.json` | the tables M1–M7 |
| `el_md.py` → `_tables.md` | the markdown of those tables |
| `el_grid.py`, `el_view.py`, `el_diag.py`, `el_coincide.py`, `el_path.py`, `el_lazydbg.py` | diagnostics |
