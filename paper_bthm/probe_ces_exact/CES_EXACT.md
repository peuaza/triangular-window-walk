# CES_EXACT — can CES's failures be patched into an exact O(n) machine?

*2026-11-09, `probe_ces_exact/`. No `.tex` was read or modified; no git command was run.
Python: `C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`. Every number is
printed by a script here; the JSON that carries it is named in §9.*

## 0. Verdict, in the order that matters

1. **The per-epoch assumption is not the failure.** "The epoch outcome is fixed by `(b2,b1)` at the
   epoch start" is **exact on every cut epoch of `pattern`** (0 mismatches in 521 294 epochs,
   `n = 256…2048`) and nearly exact on `runs` (cut-mismatch rate `4.0e-3 → 9.2e-5`, `n = 256 → 2048`).
   The first failing epoch is **the last epoch of the tape** (median and p90 of `j/n` are 1.00 in
   every cell; 98–100 % of first failures are the terminal epoch).
2. **What CES gets wrong is a terminal, not a decline.** Every tape has exactly one epoch whose true
   outcome is not a cut: a **meet** (the walk answers `δ`) or a **cycle** (the walk halts with no
   output). CES predicts a cut there. Meeting terminals sit at tiny `M0` (pattern: max 11–17 across samples at every
   `n`; runs: p90 4–16, thin tail to 439); halting terminals sit at `M0` up to `0.5–0.99 n`.
3. **CES does not decline on meeting tapes.** With `(a)` = true meeting tapes with the right `δ`,
   `(c)` = true meeting tapes wrongly called halting, CES gives, at `n = 512/1024/2048/4096`:
   pattern `(a)=198/198/200/100` and `(c)=1/1/0/0`; runs `(a)=186/197/196/98` and `(c)=11/3/2/0`
   per 200 (100 at 4096). Its large error is the opposite one: it gives a `δ` on 2–16 of every 100
   **halting** tapes. The "55–83 % declines" are the truth's halt rate: CES's nan is a correct halt
   on 99 %+ of the tapes where it declines.
4. **A bounded-read check cannot find the terminal.** Prefix-`K` verification (`K = 4` exact steps at
   every epoch) changes exactness by **zero** on every rule and every `n`, and raises the cost by
   2.1–2.5 n. The terminal scan has median length 58–124 and maximum 1 638 steps; 2–4 cells cannot see it.
5. **No exact-everywhere machine at `O(n)`.** The only exact patch is to run every epoch below a
   threshold `W2` exactly. To be exact on every tape, `W2` must exceed the **halting** terminals'
   `M0`, whose maximum tracks `0.5–0.99 n`; at `W2 = αn` the exact-epoch cost is `≈0.02–0.12 α² n²`.
   No fixed `W2` works: at `W2 = 512` one or two invented halts remain at every `n`, and
   `P(terminal M0 > 256)` stays 3–12 % from 256 to 2048. **Negative result, stated plainly.**
6. **But `e`, under the record's own definition, is `O(n)` — pattern and runs are fully linear.**
   Meeting terminals are tiny, so a fixed `W2` catches **all** of them. On 200 true meeting tapes
   per cell (100 at 4096) the best exact-on-meeting machine is:

   | rule | n=512 | 1024 | 2048 | 4096 |
   |---|---|---|---|---|
   | `pattern`, `W2` | 32 | 128 | 64 | 128 |
   | `pattern`, `e/n` | **3.71** | **4.16** | **3.69** | **3.80** |
   | `pattern`, exact meeting | 200/200 | 200/200 | 200/200 | 100/100 |
   | `runs`, `W2` | 512 | 256 | 512 | 256 |
   | `runs`, `e/n` | 58.07 | **5.12** | **8.45** | **3.84** |
   | `runs`, exact meeting | 200/200 | 200/200 | 200/200 | 100/100 |

   The `n = 512` runs cell is the one outlier: one meeting terminal at `M0 = 439` forces `W2 = n`
   there. From 1024 to 4096 runs is `3.8–8.5 n`. **This is the order drop: `e(pattern) = 3.7–4.2 n`,
   `e(runs) = 3.8–8.5 n` (n ≥ 1024), exact on every true meeting tape.** Under ORDERS §0 these rules
   are fully linear, not quasi-linear. The `n²` in the table is the halt gate, not the meeting gate.
   (Plain CES is already exact on the 2048/4096 meeting samples at 3.64/3.68 n; the `W2` row is the
   worst-case guarantee, and it costs only 0.05–0.12 n more.)
7. **`p ≤ e` holds on the same tapes and the same arm; both numbers move from their published
   values.** p on my meeting tapes with the one-pass CES arm is **3.58–3.68 n** for pattern and
   **3.40 n** for runs at 1024 (R² 0.995–1.000; §7); at `runs 512/2048` only the exact `W2` arm
   clears 0.99, at 58.07 and 8.45 n. The published `p = 4.33–4.57 n` for pattern is
   **one full pass too high**: `kp_lib._Win.__init__` re-reads the whole tape to count `2,1,1`
   occurrences, charging `2n−2` instead of `n` for the opening pass. Corrected, p moves *down* by
   `≈0.7–1.0 n`. The published `e_meet = 122–469 n` was the walk machine, not the best
   exact-on-meeting machine; e moves *down* to `3.7–8.5 n`. The apparent `4.5 > 2.57` violation was
   an accounting artefact (my first machine undercharged the opening pass by exactly `n`; that is
   fixed and every number below includes it). With both measured correctly, `p ≤ e` in every cell.

| question | answer | evidence |
|---|---|---|
| Is the per-epoch rule exact? | `pattern` yes (0 cut errors); `runs` to `1e-4` | §2 |
| Is the failure rate `o(1/n)`? | **No** — `rate·n = 0.29–0.74`, flat | §2 |
| Does a 2–4 cell check patch it? | **No** — identical exactness, +2.1–2.5 n reads | §4 |
| Exact on every tape at `O(n)`? | **No** — `W2 = Θ(n)`, cost `Θ(n²)` | §5 |
| Exact on every *meeting* tape at `O(n)`? | **Yes** — 3.69–4.16 n (pattern), 5.1–8.5 n (runs) | §6 |
| `e(pattern)`, meeting gate | **`O(n)`**, `≈3.7–4.2 n` | §6 |
| `e(runs)`, meeting gate | **`O(n)`** for `n ≥ 1024`, `≈3.8–8.5 n` | §6 |
| `p ≤ e` on the same tapes/arm? | **Yes**, 3.6 vs 3.7–8.5 | §7 |

## 1. The machine, the tapes, the accounting

The reference is `kp_lib.exact` (line-for-line the `rs_common`/`pcr_lib` walk); `ce_check.py` gates
the new machine against both. Tapes are `eu_lib.make_tape`, the `T_n` family (two independent
uniform `⌊n/2⌋`-subsets). `δ = ell/(n-1)`; a walk that narrows below 3 cells or repeats an epoch
state halts with no output (nan).

Accounting is the record's: one cell read costs 1, with multiplicity; the machine holds the tape
plus `O(log n)` bits. **One read serves all counters.** CES's shortcut keeps `S,C,H,PC` over the
region in O(1) reads per cut; the exact-epoch path keeps sums, counters and two Brent checkpoints
(not the window), so it is in the class at every `n`.

*Accounting note.* An earlier version of `ce_machine.py` charged 0 for the opening pass; the gate
that caught it is `machine(W=W0)` against `kp_lib.ces` — on `runs` the two now agree on reads
**exactly** (0 difference, every tape at `n = 256/512`), and on `pattern` they differ by exactly
`n−2`, the length of `kp_lib`'s redundant second pass. All numbers below are after the fix. This
also identifies a real overcharge in the record's CES `p` arm (pattern only).

## 2. Step 1 — the per-epoch failure rate

For every tape and every epoch of the **true** walk, compare `pred = D if not b2 else (R if b1 else
L)` with the true outcome. `M0` is the region width at the epoch start, `W` the CES buffer
(`14` pattern, `4 log2 n` runs). A **cut mismatch** cuts on the other side; an **early terminal**
(`M0 > W`) produces no cut at all.

| rule | n | tapes | meet | epochs | epochs/tape | cut-mism | cut rate | early term `>W` | term rate | rate·n | tapes with term |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pattern | 256 | 220 | 64 | 32 481 | 148 | 0 | 0 | 49 | 1.51e-3 | 0.39 | 49 (22 %) |
| pattern | 512 | 220 | 66 | 67 623 | 307 | 0 | 0 | 61 | 9.02e-4 | 0.46 | 61 (28 %) |
| pattern | 1024 | 220 | 77 | 139 128 | 632 | 0 | 0 | 46 | 3.31e-4 | 0.34 | 46 (21 %) |
| pattern | 2048 | 220 | 64 | 282 062 | 1282 | 0 | 0 | 40 | 1.42e-4 | 0.29 | 40 (18 %) |
| runs | 256 | 220 | 29 | 27 809 | 126 | 111 | 4.0e-3 | 80 | 2.88e-3 | 0.74 | 80 (36 %) |
| runs | 512 | 220 | 47 | 63 046 | 287 | 39 | 6.2e-4 | 72 | 1.14e-3 | 0.58 | 72 (33 %) |
| runs | 1024 | 220 | 56 | 133 365 | 606 | 15 | 1.1e-4 | 62 | 4.65e-4 | 0.48 | 62 (28 %) |
| runs | 2048 | 220 | 44 | 271 098 | 1232 | 25 | 9.2e-5 | 60 | 2.21e-4 | 0.45 | 60 (27 %) |

**Not `o(1/n)`.** `rate·n` is 0.29–0.74 and does not fall over an 8× range of `n`. Per the brief's
criterion, a cheap per-epoch patch cannot exist. A 60-tape `n = 4096` check gives `rate·n = 0.26`
(pattern) and `1.36` (runs); the runs sample is cycle-heavy (52/60 tapes have an early terminal above
`W = 48`), but neither is `o(1/n)`. The rest of the report shows why, and separates `e` from the
exact-machine question.

**The first failing epoch is the last.** Median and p90 of `j/n` are 1.00 in every cell. `first_kind`
is `term` for 46/46, 61/61, 40/40 (pattern) and 77/83, 71/72, 61/63, 60/61 (runs); the few `cut`
cases are the rare cut mismatch.

| rule | n | failing | median M0 | `M0≤14` | 15–64 | 65–256 | `>256` |
|---|---|---|---|---|---|---|---|
| pattern | 512 | 61/220 | 73 | 0 | 27 | 28 | 6 |
| pattern | 1024 | 46/220 | 76 | 0 | 20 | 17 | 9 |
| pattern | 2048 | 40/220 | 64 | 0 | 20 | 13 | 7 |
| runs | 512 | 72/220 | 154 | 0 | 13 | 41 | 18 |
| runs | 1024 | 63/220 | 128 | 0 | 17 | 26 | 20 |
| runs | 2048 | 61/220 | 179 | 0 | 8 | 26 | 27 |

## 3. The three counts that decide everything

For a machine on a tape set with true meeting tapes `M` and true halting tapes `H`:

* **(a)** meeting tapes with the right `δ` — e needs `(a)=|M|`;
* **(b)** halting tapes the machine also calls halting — only the halt gate needs this;
* **(c)** meeting tapes wrongly called halting (declines) — e needs `(c)=0`;
* **(d)** halting tapes given a `δ` ("invented") — the halt-gate error, irrelevant to e;
* **(e)** meeting tapes with a wrong `δ`.

CES on 220 random tapes per cell (`ce_align.py`) is the cleanest view of the premise:

| rule | n | (a) exact | both nan | (d) invented | (c) declined | (e) wrong | CES nan rate | truth halt rate |
|---|---|---|---|---|---|---|---|---|
| pattern | 512 | 66 | 133 | 21 | 0 | 0 | 70 % | 70 % |
| pattern | 1024 | 77 | 124 | 19 | 0 | 0 | 65 % | 65 % |
| pattern | 2048 | 64 | 145 | 11 | 0 | 0 | 71 % | 71 % |
| runs | 1024 | 55 | 148 | 16 | 1 | 0 | 68 % | 75 % |
| runs | 2048 | 44 | 168 | 8 | 0 | 0 | 80 % | 80 % |

Every CES decline is a correct halt (at most one exception in 220 tapes). **The failure is the 5–10 %
of halting tapes where it answers — an invented `δ`, not a decline.**

*Cost note.* The often-quoted `4.5 n` for pattern is `kp_lib.ces`'s charge, and it includes the
redundant second pass of §7 (`2n−2` reads for the opening). The same predictor with a one-pass
opening costs **3.6 n over all tapes and 3.6 n over meeting tapes** at n=512…4096; `runs` has no
second pass and its `3.4 n` is correct. CES's real price is therefore 3.4–3.6 n, not 4.5 n.

## 4. The patch machine and the bounded check

`ce_machine.py` is one read-counted machine: `W` endgame buffer; `W2` run every epoch with `M ≤ W2`
exactly (incremental window stats, **Brent** cycle detection, O(1) state); `K` prefix verification;
`oracle_term` diagnostic. Gates: (A) `W=3, W2=n` matches `kp_lib.exact` on every tape (delta and
nan, `n=64,128,256`); (B) `W=W0` matches `kp_lib.ces` on every tape, and on **reads** exactly for
runs / to `n−2` for pattern; (C) `K=4` agrees with (B) wherever the prefix does not terminate.
`pcr_lib.Runner` agrees with `kp_lib.exact` on delta and cut count at 512/2048.

**`K = 4` buys nothing.** Pattern 512: `3.63 → 6.06 n`, exact meeting `198/200 → 198/200`;
pattern 1024 `3.58 → 6.03`, `198/200` unchanged; runs 512 `3.39 → 5.47`, `186/200`; runs 1024
`3.40 → 5.49`, `197/200`. The terminals' median length is 58–124 (max 1 638); four steps never
reach one.

## 5. `W2`: exactness on every tape is bought quadratically — the negative result

`runs n = 512`, 200 meeting + 100 halting tapes (`ce_patch.py`; all `e/n` are meeting means and
include the opening pass):

| `W2` | `e/n` | (a) exact meeting | (b) halt correct | (d) `δ` on halt |
|---|---|---|---|---|
| 36 (= `W`) | 3.39 | 186/200 | 92/100 | 8 |
| 64 | 3.47 | 186/200 | 93/100 | 7 |
| 128 | 4.17 | 191/200 | 96/100 | 4 |
| 256 (= n/2) | 9.92 | 196/200 | 99/100 | 1 |
| 384 (= 3n/4) | 26.43 | 198/200 | 100/100 | 0 |
| 512 (= n) | 58.07 | 200/200 | 100/100 | 0 |

The last 14 meeting tapes cost a factor 17; the last 2 another 2.2. The cost of an exact `W2 = αn`
threshold is `≈0.02–0.12 α² n²` (the exact epochs below `αn` have scan lengths proportional to
`M`), i.e. the full walk's order with a smaller constant. Pattern behaves identically once `W2` is
large enough (512: `W2=256` → 13.05 n, 0 invented on the sample; `W2=n` → 65.71 n).

**No fixed `W2` is exact**, because the halting tail does not vanish (`ce_step1.py`):

| `P(terminal M0 > T)` | n=256 | 512 | 1024 | 2048 |
|---|---|---|---|---|
| pattern, T=32 | 19.5 % | 24.1 % | 18.6 % | 13.6 % |
| pattern, T=256 | — | 2.7 % | 4.1 % | 3.2 % |
| pattern, max M0 | 219 | 485 | 775 | 1112 |
| runs, T=32 | 36.4 % | 32.7 % | 28.2 % | 27.3 % |
| runs, T=256 | — | 8.2 % | 9.1 % | 12.3 % |
| runs, max M0 | 256 | 504 | 936 | 2022 |

At `W2 = 512` runs still invented 2/100 at 2048 and 2/50 at 4096; at `W2 = 256`, 1–2 at every size.
The maximum is `0.5–0.99 n`, so an exact-everywhere machine needs `W2 = Θ(n)` and pays `Θ(n²)`.

**The oracle ceiling** (force the one true terminal epoch to be run exactly; diagnostic):

| rule | 512 | 1024 | 2048 | 4096 |
|---|---|---|---|---|
| pattern `e/n` | 3.63 | 3.58 | 3.64 | 3.68 |
| pattern exact all | 300/300 | 300/300 | 300/300 | 150/150 |
| runs `e/n` | 3.40 | 3.40 | 3.44 | 3.49 |
| runs exact all | 300/300 | 300/300 | 299/300 | 150/150 |

The one residual runs failure is a cut mismatch at an interior epoch. **The distance from `3.4 n` to
`Θ(n²)` is entirely terminal detection.**

## 6. `e`: pattern and runs are fully linear under the meeting gate

The two gates in the record disagree:

* **ORDERS §0:** `e` = exact on every **meeting** tape, mean reads over meeting tapes.
* **E_UNIFORM §1.4:** a machine is admitted only if it reproduces `δ` **and** the no-output halts.

Under the second gate CES is not a machine and §5 says no cheap patch makes it one. Under the first,
the halting terminal is irrelevant. Meeting terminals are tiny:

| rule | n | meeting terminal `M0` med / p90 / p99 / max | halting terminal `M0` med / p90 / max |
|---|---|---|---|
| pattern | 512 | 1 / 7 / 13 / 16 | 2 / 149 / 485 |
| pattern | 1024 | 1 / 5 / 10 / 11 | 2 / 164 / 775 |
| pattern | 2048 | 1 / 6 / 10 / 12 | 2 / 97 / 1112 |
| runs | 512 | 1 / 16 / 280 / 439 | 20 / 282 / 504 |
| runs | 1024 | 1 / 4 / 29 / 172 | 8 / 290 / 936 |
| runs | 2048 | 1 / 7 / 20 / 26 | 8 / 313 / 2022 |

(`ce_terminal.py`, 100 meeting + up to 360 halting tapes per cell.) A fixed `W2` therefore catches
**all** meeting terminals while missing most halting ones. Measured on 200 true meeting tapes per
cell (100 at 4096):

| rule | n | best fixed `W2` | `e/n` meet | (a) exact meet | (c) meet halted | (d) `δ` on halt |
|---|---|---|---|---|---|---|
| pattern | 512 | 32 | **3.71** | 200/200 | 0 | 10 |
| pattern | 1024 | 128 | **4.16** | 200/200 | 0 | 2 |
| pattern | 2048 | 64 | **3.69** | 200/200 | 0 | 0 |
| pattern | 4096 | 128 | **3.80** | 100/100 | 0 | 2 |
| runs | 512 | 512 | 58.07 | 200/200 | 0 | 0 |
| runs | 1024 | 256 | **5.12** | 200/200 | 0 | 2 |
| runs | 2048 | 512 | **8.45** | 200/200 | 0 | 2 |
| runs | 4096 | 256 | **3.84** | 100/100 | 0 | 2 |

**Plainly: an exact-on-every-true-meeting-tape machine exists at `3.69–4.16 n` for pattern (all four
sizes) and `3.84–8.45 n` for runs (n ≥ 1024; 58 n at 512 because of one terminal at `M0 = 439`).
`e` is `O(n)` with a small constant; under ORDERS §0, pattern and runs are fully linear, not
quasi-linear.** The `n²` entries in ORDERS are an artefact of the halt gate: the exact-*everywhere*
machine needs `W2 = Θ(n)` and is `Θ(n²)`.

**Caveat with teeth.** The runs meeting-terminal tail is not proven bounded. The largest measured is
`M0 = 439` at `n = 512`; 1024/2048 gave 172/26, so the tail does not grow in these data, but a single
meeting tape above the fixed `W2` breaks the "every meeting tape" quantifier. `e(runs) = O(n)` is
measured on 512–4096, not proved for all `n`. Pattern's meeting terminal is `≤17` at every size; its
`W2` need not grow.

## 7. `p ≤ e` on the same tapes, same arm — which number moves

`p`'s protocol (`probe_table_p/pt_cell.py`): population = meeting tapes; axis = the arm's mean
per-tape spend; non-finite predictions scored at the population mean; `p` = smallest spend with
population `R² ≥ 0.99`. On my meeting tapes, same arm, same accounting (`ce_p.py`):

| rule | n | record `p/n` | my `p/n` | p arm | R² | best exact `e/n` | p ≤ e |
|---|---|---|---|---|---|---|---|
| pattern | 512 | 4.33 | **3.63** | CES | 0.9995 | 3.71 | yes |
| pattern | 1024 | 4.42 | **3.58** | CES | 0.9949 | 4.16 | yes |
| pattern | 2048 | 4.57 | **3.64** | CES | 1.0000 | 3.64 | yes |
| pattern | 4096 | — | **3.68** | CES | 1.0000 | 3.68 | yes |
| runs | 512 | 59.50 | **58.07** | w2_n | 1.0000 | 58.07 | yes |
| runs | 1024 | 3.29 | **3.40** | CES | 0.9919 | 5.12 | yes |
| runs | 2048 | 225.60 | **8.45** | w2_512 | 1.0000 | 8.45 | yes |
| runs | 4096 | — | **3.49** | CES | 0.9997 | 3.84 | yes |

For pattern 2048/4096 plain CES is itself exact on every meeting tape in the sample, so the
measured best `e` is 3.64/3.68; the `W2` rows of §6 are the worst-case guarantee at 3.69/3.80.

**Which number moves: both, and the apparent violation disappears.**

* `p` moves **down** by `0.7–1.0 n` for pattern (4.33→3.63, 4.42→3.58, 4.57→3.64):
  `kp_lib._Win.__init__` reads the tape once for `S,H,C` and a second time for the `2,1,1` count,
  charging `2n−2` for the opening pass. One pass feeds all four counters, so the true CES spend is
  `record p − (n−2)/n`. On `runs` there is no second pass, and the 1024 cells agree (3.29 vs 3.40).
  At `runs 2048` the record's `p` is `225.6 n` because no cheap arm cleared `R²=0.99` on its sample;
  the exact `W2=512` arm clears it at `8.45 n`, so p moves by a factor 27 there. At `runs 512` both
  are ~58–59 n (a genuine 0.99 barrier, not an accounting defect).
* `e` moves **down** by two orders of magnitude for pattern (published `e_meet = 122–469 n` is the
  walk machine) to `3.69–4.16 n`; at `runs 2048` it moves from `225.6 n` to `8.45 n`. The `W2`
  machine is exact on every meeting tape in every cell.
* The old `p = 4.5 n` vs `e = 2.57 n` comparison was never a real violation: `4.5` included the
  redundant pass and `2.57` omitted the opening pass. Corrected, `p ≤ e` in every cell, with
  `p = e` where the cheapest 0.99 arm is itself exact (`runs 512/2048`, `pattern 2048/4096`).

## 8. What this means for the table

1. **State the gate.** Under the **halt gate** (E_UNIFORM), `pattern` and `runs` keep `e = n²` and
   CES is a predictor, not a machine. Under the **meeting gate** (ORDERS §0), both are **fully
   linear**: `e(pattern) ≈ 3.7–4.2 n`, `e(runs) ≈ 3.8–8.5 n` (n ≥ 1024), exact on 100 % of the
   meeting tapes measured at 512/1024/2048/4096.
2. **The CES `4.4 n`/`3.4 n` rows** are the cost of a nearly-exact predictor. Once the redundant
   second pass is removed, the one-pass CES costs **3.4–3.6 n on both meeting and all tapes** (the
   quoted `4.5 n` for pattern is `kp_lib.ces`'s double-pass charge, §3). They are not evidence for
   `n²` and not evidence that "CES can't give e": with a 32–512 cell exact-epoch threshold the
   meeting delta is exact on every meeting tape measured.
3. **Pattern's meeting failure is a buffer 3 cells too small.** `W = 14` misses terminal meets at
   `M0 = 15–17`. `W = 20–32` removes them. They are not a deep defect.
4. **The exact-everywhere machine does not exist at `O(n)`**, and the obstruction is sharp: the
   terminal epoch is a scan orbit; detecting it costs `Ω(ln)=Ω(M0)` reads at the only epoch that is
   terminal, and one must decide at every epoch. `K=4` fails; `W2=Θ(n)` costs `Θ(n²)`.

## 9. Files

| file | what |
|---|---|
| `ce_lib.py` | exact-walk epoch trace and per-epoch failure classifier |
| `ce_step1.py` → `ce_fail_<rule>_<n>.json` | per-epoch rates, first-fail distribution, terminal M0 |
| `ce_terminal.py` → `ce_terminal_<rule>_<n>.json` | terminal `M0` by outcome class (meet vs halt) |
| `ce_align.py` → `ce_align_<rule>_<n>.json` | CES vs truth: exact/both-nan/invented/declined/wrong |
| `ce_machine.py` | read-counted machine: CES, `W2` exact epochs, `K` prefix, Brent, oracle |
| `ce_check.py` | gates A/B/C against `kp_lib.exact` and `kp_lib.ces` |
| `ce_patch.py` → `ce_patch_<rule>_<n>.json` | 200 meeting + 100 halting tapes per cell, all variants |
| `ce_dyn.py` → `ce_dyn_<rule>_<n>.json` | anatomy of the early-terminal epochs (length, cycles) |
| `ce_summary.py` | tables of §§2, 4–6 (including (a)–(e)) |
| `ce_p.py` | p vs e on the same meeting tapes and arm (§7) |
| `_wintest.py`, `_rgate.py`, `_gate_pcr.py` | incremental-stat, CES-read and reference gates |

All scripts run on the project venv python; one JSON per cell. No `.tex` and no git.
