# PEAK_FIX — the first-epoch repair of `peak`, and why it does not close the last three tapes

*2026-11-10, `probe_peak_fix/`. No `.tex` was read or modified; no git command was run.
Python: `C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`. Every number below is
printed by a script in this directory; the JSONL that carries it is named in §7.*

## 0. Verdict

**The prescribed fix does not work, and the premise it rests on is false for two of the three tapes.**
On the 1 500-tape grid of EN_LINEAR, the `peak` jump machine at `W2 = 256` fails on three meeting
tapes, but only **one** of the three fails in its first epoch (`i = 701 615`, `n = 2048`,
`M0 = 2048`, and it is the first of a 15-epoch run); the other two fail in the middle of the tape, at
`M0 = 576` (`i = 302 049`) and `M0 = 1 529` (`i = 302 287`), with epoch 0 a `b2 = 0` double cut and a
law-correct scan respectively.  A literal first epoch cannot touch either.

The machine with the first epoch run literally (charged its worst case, `n` reads) is exact on the
**same 1 497/1 500** at `W2 = 256` and the same **1 496/1 500** at `W2 = 128`.  It changes no
`(δ, end)` on any of the 1 500 tapes.  Its measured cost is exactly the +1.0 `n` the accounting
predicts — `e/n = 3.049 / 3.029 / 2.999 / 2.997` by size, `3.014` over the 1 500 — and it buys
nothing.  (Measured on the epoch's own length instead of the `n` cap, the surcharge is `0.048 n`, not
`n`: `e/n = 2.095 / 2.076 / 2.057 / 2.061`, `2.064` over the 1 500.)

Two further consequences for the record:

1. **EN_LINEAR's "1 497/1 500" is the `W2 = 256` grid, not the `2.0 n` machine.**  At `W2 = 128`
   there is a fourth failure, `n = 1024`, `i = 1 247`, first bad epoch at `M0 = 149`; the `2.0 n`
   machine is **1 496/1 500**.  The `3/600 at n = 2 048` in EN_LINEAR §7 counts both 2 048 replays;
   the two from `i = 300 000` were never analysed, and they are not first-epoch cuts.  So "exactly
   three, all first-epoch" reproduces neither at 128 (four, one first-epoch) nor at 256 (three, one
   first-epoch).
2. **The residual is interior.**  To cover `M0 = 1 529` at `n = 2 048` a jump machine needs
   `W2 ≥ 1 529` (0.75 `n`); that tape then costs `42.5 n` reads.  The top tape `i = 701 615` is
   repaired by the one-entry cut memory at `W2 = 1 536` for `24.2 n`, or by the plain jump machine at
   `W2 = n` for `223.9 n`.  A first-epoch repair alone reaches neither, and a `W2` that does is no
   longer a constant-read machine.

**One-line verdict: the fix does not work because two of the three failures are not first-epoch cut
errors (`M0 = 576` and `M0 = 1 529` at `n = 2048`, with epoch 0 a double cut and a law-correct scan),
and the one that is (`M0 = 2048`) is the first of a 15-epoch run; a literal first epoch changes no
outcome on any of the 1 500 tapes and adds 1.0 `n` for nothing.**

## 1. Reproduction of the 1 500-tape grid

The machine is `el_cell.jm_walk` (probe_en_linear), copied into `pf_lib.jm_base_audit` and gated
against it: 180 (tape, `W2`) cells at `n = 512, 1024`, `W2 = 8 … 256`, **0 mismatches** in `end`,
`reads`, `epochs`, `cuts`, `eg`.  The tape family is `el_lib.tape`, the same draws as `el_peak_*`:

| segment | `n` | tapes | base `W2=128` | fix `W2=128` | base `W2=256` | fix `W2=256` |
|---|---|---|---|---|---|---|
| main | 512 | 200 | 200/200 · 2.05 | 200/200 · 3.05 | 200/200 · 2.68 | 200/200 · 3.68 |
| main | 1024 | 200 | **199/200** · 2.03 | **199/200** · 3.03 | 200/200 · 2.15 | 200/200 · 3.15 |
| main | 2048 | 200 | 200/200 · 2.00 | 200/200 · 3.00 | 200/200 · 2.01 | 200/200 · 3.01 |
| main | 4096 | 100 | 100/100 · 2.00 | 100/100 · 3.00 | 100/100 · 2.03 | 100/100 · 3.03 |
| replay `i=900 000` | 512 | 200 | 200/200 · 2.05 | 200/200 · 3.05 | 200/200 · 2.33 | 200/200 · 3.33 |
| replay `i=500 000` | 1024 | 200 | 200/200 · 1.97 | 200/200 · 2.97 | 200/200 · 2.13 | 200/200 · 3.13 |
| replay `i=700 000` | 2048 | 200 | **199/200** · 1.99 | **199/200** · 2.99 | **199/200** · 2.07 | **199/200** · 3.07 |
| replay `i=300 000` | 2048 | 200 | **198/200** · 2.02 | **198/200** · 3.02 | **198/200** · 2.07 | **198/200** · 3.07 |
| **total** | | **1 500** | **1 496/1 500** | **1 496/1 500** | **1 497/1 500** | **1 497/1 500** |

Outside the 1 500, the 4096 replay (`i = 900 000`, 100 tapes) is 100/100 at both `W2`
(`e/n` 1.995 at 128, 2.013 at 256; fix 2.995 / 3.013), matching `el_reprobe_peak_4096_900000.json`.

The base columns reproduce `el_peak_*.json` exactly (512 2.049/2.678, 1024 2.029/2.151, 2048
1.999/2.014, 4096 2.000/2.030).  "fix" is `pf_lib.jm_fixed`: the first epoch run literally, charged
`n` reads, then the same jump machine.  It is exact on exactly the tapes the base machine is exact
on; the 4 (at 128) and 3 (at 256) failures are the same tapes before and after.

## 2. The failure set

| # | `n` | `i` | first bad epoch | `M0` | pred → out | len | `nbad` | base end | epoch 0 |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 1024 | 1 247 | 443 | **149** | cutR → cutL | 26 | 24 | 1 (region) | double cut (`b2=0`, out 3) |
| 2 | 2048 | 701 615 | 0 | **2 048** | cutL → cutR | 440 | 15 | 2 (cycle) | `b2=1,b1=0` scan, out cutR |
| 3 | 2048 | 302 049 | 740 | **576** | cutL → cutR | 14 | 1 | 1 (region) | double cut (`b2=0`, out 3) |
| 4 | 2048 | 302 287 | 509 | **1 529** | cutR → cutL | 520 | 10 | 1 (region) | `b2=1,b1=0` scan, out cutL (**law-correct**) |

Failure 1 is the extra tape of the `W2 = 128` grid; it passes at `W2 = 256` because `M0 = 149 ≤ 256`.
Failures 2–4 are the three `W2 = 256` failures.  Only #2 begins at `M0 = n`.

The law is not close to exact above the buffer: at `W2 = 256` **22 tapes carry 346 mispredicted
epochs with `M0 > 256`**, and 19 of the 22 still meet at the right `δ`.  The three that do not are
the three failures above.  Two of those three (and the fourth tape of the `W2 = 128` grid) fail at
`M0` well inside the tape, where no first-epoch device can see them.

`i = 302 287` is the instructive one: epoch 0 is a 481-step `b2=1,b1=0` scan that cuts **left**, and
the law predicts left — correct.  Five hundred epochs later, at `M0 = 1 529`, the same signature
produces a 520-step scan that cuts **right**, the law says left, and the machine dies.  Nothing at
the top of the tape predicts it.

## 3. Why the literal first epoch fixes none of them

`jm_fixed` opens by running epoch 0 as the walk does (`el_lib.trace_walk`'s own control flow) and
then handing the new region to the jump machine.  If epoch 0's cut is the one the law would have
predicted, the machine's subsequent trajectory is *identical* to the base machine's — the literal
epoch only changes the charge — so no exactness can change.  If epoch 0's cut is law-wrong, the
machine reproduces the truth's region and then evaluates the next epoch; on `i = 701 615` the next
epoch is the same `b2=1,b1=0` scan (len 438) and the law is wrong again:

```
fixed machine, i = 701615, W2 = 256, after the literal epoch 0:
  k=1  M=2047  b2=1 b1=0  pred=1 (cutL)     truth epout=2 (cutR)
```

Fifteen consecutive epochs (M0 2 048 → 2 034, lengths 440 → 408) share that signature; the fix
repairs the first and the machine diverges at the second.  On `i = 302 049` and `i = 302 287` epoch 0
is law-correct (and a double cut, for the first), so the literal epoch leaves the region at the same
point the base machine had it; the interior failures at `M0 = 576` and `M0 = 1 529` are untouched.
On `i = 1 247` the first epoch is a double cut, so the `M0 = 149` failure is untouched.

Only four tapes in the 1 500 have a law-wrong epoch 0 at all (`i = 58` and `900 869` at `n = 512`,
`701 615` and `302 737` at `n = 2048`), and three of the four pass anyway.  The literal first epoch
does not repair the one that fails either: epoch 1 repeats the same `b2=1,b1=0` signature and the
machine diverges there.  The other two `W2 = 256` failures are interior and have no law-wrong epoch
at the top at all.

Across the whole grid the fix is a **total no-op on `(δ, end)`**: over all 1 600 tapes (the 1 500
plus the 4096 replay) at each of `W2 = 128` and `256`, `end` differs in 0 tapes, `δ` differs in 0
tapes, and `exact` differs in 0 tapes.  It is not that it repairs the wrong tapes; it repairs nothing
and costs 1.0 `n`.

## 4. The prescribed fix, measured (200 tapes per size)

"measured" = the machine's charge; "actual" = `base + (len(epoch 0) − the jump cut it replaced)/n`.
The `n` cap is the task's convention ("cost at most n reads for that epoch, once").  The 4096 row
uses the 100 main + 100 replay tapes; the other rows are the main grids.

| `n` | tapes | base `e/n` | fix `e/n` (actual) | fix `e/n` (n-cap) | fix exact |
|---|---|---|---|---|---|
| 512 | 200 | 2.049 | 2.095 | **3.049** | 200/200 |
| 1024 | 200 | 2.029 | 2.076 | **3.029** | 199/200 |
| 2048 | 200 | 1.999 | 2.057 | **2.999** | 200/200 |
| 4096 | 200 | 1.997 | 2.061 | **2.997** | 200/200 |

The `n`-cap column is the predicted `2.0 + 1.0 = 3.0`, measured; over the 1 500 tapes it is `3.014`.
The actual column is the same machine charged the epoch's real length: `2.0 + 0.05` (`2.064` over
the 1 500).  Neither column moves the exactness count; the 1 024 failure in the table is `i = 1 247`
(§2), and the other three failures of the 1 500 are in the 2 048 replays.

## 5. A one-entry cut memory (the strongest first-epoch variant), and its ceiling

`pf_sticky.jm_sticky`: literal epoch 0; if its cut disagrees with the law, remember
`(b2,b1) → cut side` and predict the remembered side while the signature repeats.  It repairs
`i = 701 615` (the 15-epoch run) at `e/n = 2.15`, and the three law-wrong-but-passing tapes
(`i = 58`, `900 869`, `302 737`) are unharmed; it cannot reach failures #3 and #4, whose law-wrong
epoch is interior.

| `n` | tapes | sticky exact | `e/n` (actual) | `e/n` (epoch 0 at `n`) |
|---|---|---|---|---|
| 512 | 400 | 400/400 | 2.550 | 3.502 |
| 1024 | 400 | 400/400 | 2.188 | 3.137 |
| 2048 | 600 | 598/600 | 2.101 | 3.049 |
| 4096 | 200 | 200/200 | 2.085 | 3.021 |
| **total** | **1 600** | **1 598/1 600** | **2.233** | **3.180** |

(The 1 500 of §1, without the 4096 replay, is 1 498/1 500 at 2.241 actual / 3.192 capped.)  1 498, not
1 500.  The two short are `i = 302 049` and `i = 302 287`.

A literal *prefix of `n` steps* (epochs run until the cumulative step count reaches `n`, then the
jump machine) does repair `i = 302 049` — its whole walk is 1 108 steps at `n = 2 048` — but neither
`i = 701 615` (the 15-epoch scan run outlasts `n`) nor `i = 302 287` (its bad epoch sits at step
≫ `n`).  It too tops out at 1 498.

**What exactness would cost.**  `W2 = 1 536` (0.75 `n` at `n = 2048`) turns failures #3 and #4
literal, and the sticky epoch 0 turns #2: on the four tapes, `i = 302 049` is exact at 1.67 `n`,
`i = 302 287` at 42.5 `n`, `i = 1 247` at 4.58 `n`, `i = 701 615` at 24.2 `n`.  A 20-tape sample of
the 2 048 main grid at `W2 = 1 536` averages `e/n = 7.06` (min 1.62, max 30.73).  That is the price
of the interior residual: a `Θ(n)`-wide endgame, not a constant factor.

## 6. State size and charge rule (for comparability)

The machine holds `L, R` (2 `⌈log n⌉` bits), the window headers `ell, r` while an epoch is literal
(2 `⌈log n⌉`), the peak counter `P = #{i ∈ [L+1,R−1] : a_{i−1} < a_i > a_{i+1}}` (`≤ log n` bits), and
two cached boundary cells per side (2 bits each) — `O(log n)` bits in all, no array.  The cache is
what makes the counter slide at **1 read per cut**: retiring the peak at `L+1` needs
`a_L, a_{L+1}, a_{L+2}`, and `a_L, a_{L+1}` are in the state, so one read of `a_{L+2}` decides it
(symmetrically on the right).

Charge, the record's: `n` for the opening pass; 1 read per cell a cut gives up (2 for a double cut);
1 read per step in any literal epoch, including the first; 0 for the predicate evaluation (the
counter slides).  A literal first epoch is at most `n` reads, so the machine's bound is
`n + n + n/2 + endgame` = `3 n` + endgame.  The `e/n = 3.0` in §4 is that bound; the machine
actually reads `2.05 n` because the first epoch is usually a one-step double cut (mean length
`0.048 n`).

## 7. Files

| file | what |
|---|---|
| `pf_lib.py` | `jm_base_audit` (a copy of `el_cell.jm_walk`, with the epoch audit), `jm_fixed` (literal epoch 0), `exact_of`, `first_divergence` |
| `pf_var.py` | `jm_sel` (per-epoch literal decision), `jm_adapt` (literal until the law agrees) |
| `pf_sticky.py` | `jm_sticky` (literal epoch 0 + one-entry cut memory) |
| `pf_repro.py` → `seg_<n>_<off>.jsonl` | one segment: trace, base and fixed at `W2 = 128/256`, per-tape top epochs, bad-epoch list |
| `pf_agg.py` | the 1 500-tape table of §1 and the failure list of §2 |
| `pf_sticky_summary.py` | the sticky table of §5 (identity + measured law-wrong tapes) |
| `pf_lawwrong.py` | which tapes have a law-wrong epoch 0, and the charge of a literal epoch 0 |
| `pf_one.py` | one tape, all `W2`, and the truth's mispredicted epochs |
| `pf_gate.py` | the three gates: base ≡ `el_cell.jm_walk`; the fix is a no-op on `(end, δ, exact)` over 1 600 tapes; the 1 500-tape totals |
| `pf_adapt.py` | the adaptive variant on one segment (kept for the record; it overpays on tapes whose long scans are law-wrong, e.g. 17.2 `n` on `i = 900 869`) |

**Gate.** `pf_gate.py`: 180 cells, 0 mismatches; 1 600 tapes, 0 `(end, δ, exact)` differences at
either `W2`; the 1 500-tape totals 1 496 (base and fix, `W2 = 128`) and 1 497 (base and fix,
`W2 = 256`).  The base columns of §1 reproduce `el_peak_*.json` number for number.
