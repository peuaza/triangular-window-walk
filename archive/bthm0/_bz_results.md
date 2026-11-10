# Breaking the ceiling: BTHM-0 is exactly computable in ~5.3 n reads at O(log n) bits

Status: **the ceiling conjecture (and the dichotomy, and the open problem) are false as
stated.**  The canonical rule's meeting position is computable exactly, cell by cell, in
`5.3 n` reads with `~11 log2 n` bits of state.  R² = 1.000, not 0.766.

All numbers below come from runs of the scripts in this directory, with the commands and
seeds quoted.  Nothing here modifies any `.tex`.

## The mechanism the published cost model misses

`_phase.py` / `_ktab.tex` model the cost as `n + cuts × n/k` — i.e. they assume that after
every cut the window must be rebuilt from scratch, at `Θ(m)` reads, and with `Θ(n)` cuts of
size `Θ(n)` that gives `Θ(n²)`.

That assumption is wrong.  A cut fires only in this test:

```python
if ny < L or nx > R:
    if ny < L:  R = r - 1
    else:       L = ell + 1
    ell, r = L, R
```

`ny < L` means the move `y → y-1` left the interval, which can only happen when `y == L`.
Then the new window is `[L, x-1]` — **the old window with its last cell deleted**.  Symmetric-
ally, `nx > R` means `x == R` and the new window is the old window with its first cell deleted.
So every cut is an O(1) update of the accumulators, never a rebuild.

The state needed is constant-size and integral:

| quantity | meaning |
|---|---|
| `S`, `Q` | `Σ a_i`, `Σ a_i²` over the window |
| `Rr` | `Σ a_i a_{i+1}` over the window's adjacent pairs |
| `Sh`, `hy`, `hx` | sum and endpoints of the head window `[y, y+⌊m/2⌋-1]` |
| `cy`, `cx` | cached `a_y`, `a_x` |
| `y`, `x`, `L`, `R` | window and interval |

`b2` becomes `S·k > Sh·m`; `b1` becomes `(m-1)·Rr - A·B > 0` with `A = S - a_x`,
`B = S - a_y`, plus `(m-1)(Q - a_x²) - A² > 0` and the symmetric variance test — all exact
integer arithmetic if the tape is held scaled by 2 (values 0,1,2).

Measured workspace (`_bz_bits.py`): **87 bits at n=256, 97 bits at n=512** — 10.9 and 10.8
times `log2 n`.  O(log n) bits, and inside the merged draft's `O(log² n)` too.

Measured per-step cost (`_bz_final.py`): **the largest number of cell reads in any single
step is 3**, at every n from 32 to 1024.  The only branch that would need a rebuild
(`nx < L or ny > R`) never fires: 0 occurrences on every tape in every run.

## Read cost against the budget

`n²/ln n` is the conjecture's threshold.  Meeting set only (tapes that produce a δ), since
the law is defined there; `τ = n²/ln n`.

| n | tapes | meeting | reads/tape (mean) | reads/n | reads/τ | worst tape /τ | tapes over τ | max step reads |
|---|---|---|---|---|---|---|---|---|
| 32 | 300 | 85 | 168.2 | 5.26 | 0.569 | 0.640 | 0/85 | 3 |
| 64 | 300 | 69 | 340.2 | 5.32 | 0.345 | 0.415 | 0/69 | 3 |
| 128 | 300 | 74 | 685.2 | 5.35 | 0.203 | 0.232 | 0/74 | 3 |
| 256 | 300 | 49 | 1361.8 | 5.32 | 0.115 | 0.132 | 0/49 | 3 |
| 256 | 400 | 70 | 1383.8 | 5.41 | 0.117 | 0.130 | 0/70 | 3 |
| 512 | 150 | 22 | 2733.0 | 5.34 | 0.065 | 0.075 | 0/22 | 3 |
| 1024 | 100 | 20 | 5501.4 | 5.37 | 0.036 | 0.041 | 0/20 | 3 |

Commands (interpreter `C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`,
run in `paper/code/bthm0/`):

```
python _bz_final.py 32 300 97
python _bz_final.py 64 300 193
python _bz_final.py 128 300 385
python _bz_final.py 256 300 769
python _bz_final.py 256 200 777
python _bz_final.py 512 150 1025
python _bz_final.py 1024 100 999
```

`reads/n` is flat at 5.3 over a 32× range of n while the budget grows like `n²/ln n`, so the
fraction of budget used falls like `5.3 ln n / n`: 0.57 at n=32, 0.036 at n=1024.  The walk
is in the *cheap* class, not the expensive one.

For comparison with the published cost figures: the manuscript's `0.55 n²` at n=256 is 36 045
reads, against a measured 1 362 (26× fewer); its `3.6 n²` at n=128 is 58 982, against a
measured 685 (86× fewer).  Those figures were produced by rescanning the window at every step.

Worst-case bound, if one wants one without measurement: ≤ 3 reads per step, steps capped at
`20n`, plus a `3.5n` build ⇒ `T ≤ 63.5 n` for every tape, meeting or not.  That is `o(n²/ln n)`.
(Measured worst tape among meeting tapes is ≈ 6.1 n.)

## Correctness

Three independent checks:

1. `_bz_inc.py check` mode asserts, **at every step of every tape**, that `S`, `Q`, `Rr`, `Sh`
   and both decision bits equal a from-scratch integer computation on the same window, and
   that both the window and the interval equal what a cut says they should be.
2. `walk_exact_bf` is a separate, deliberately naive oracle (recompute the window statistics
   by direct summation each step, integer arithmetic).  Same output on every tape tested at
   n = 32, 64, 128, 256.
3. Against the published implementation (`walk_ref`, copied verbatim from `_onetable.py`).

On the meeting set the simulated output equals the published output:

| n | tapes | both meet | R²(sim, published X) | bit-identical |
|---|---|---|---|---|
| 32 | 300 | 83 | 1.000000 | 1.0000 |
| 64 | 300 | 64 | 0.997498 | 0.9688 |
| 128 | 300 | 70 | 1.000000 | 1.0000 |
| 256 | 300 | 45 | 1.000000 | 1.0000 |
| 256 | 200 | 46 | 0.999990 | 0.9783 |
| 256 | 400 | 68 | 1.000000 | 1.0000 |
| 512 | 150 | 21 | 1.000000 | 1.0000 |
| 1024 | 100 | 20 | 1.000000 | 1.0000 |

## One wrinkle, and why it does not matter

The published code computes `b1` with `np.corrcoef(...) > 0`.  On 6–14% of steps the sample
autocorrelation is **exactly zero** (the tape is integer-valued, so this is common, not
measure-zero), and then the sign of the test is decided by float roundoff: numpy returns
`+2.9e-17` where the exact value is 0, so `corr > 0` comes out `True`.

* `M_exact` = the rule read literally (exact arithmetic, tie ⇒ `False`): the fast path above,
  ≤ 3 reads per step, no scan.
* `M_code` = what the published scripts actually run.

`_bz_inc2.py` reproduces `M_code` bit-for-bit by paying an `O(m)` float scan **only on tie
steps** (3–4 per meeting tape).  Reproduction is exact on every tape tested:
250/250 at n=64, 200/200 at n=128, 120/120 at n=256.  On meeting tapes that costs
mean 1542 reads at n=256 (6.0 n, 0.13 τ), worst 3293 (0.28 τ) — still far under budget.

So both readings of the rule are inside the budget, exactly:

* `M_exact`: R² = 1.000 at 5.3 n reads.
* `M_code`: R² = 1.000 at 6.0 n reads.

## Why the statistical route really is stuck (so the refutation is algorithmic, not statistical)

Cheap single-pass features (~35 of them: multi-scale block means/sds/argmax, block ac₁,
longest runs and their positions, prefix-imbalance profile, …), gradient boosting, 5-fold CV,
n=256, 700 meeting tapes (`_bz_feat.py 256 700 4242`): **CV R² = 0.7189**; greedy nested
selection tops out at 0.7332.  ac₁ alone: 0.7073 (GB), 0.6848 (cubic) — the last reproduces
the manuscript's own 0.683, so the tape ensemble here matches the paper's.

So the 0.766 figure is real for fitted models on cheap statistics, and 0.85 is not reachable
that way.  It is reachable by simulation, because simulation was never expensive.

## What this does to the manuscript's claims

* **Ceiling conjecture** (`_b0_ceiling.tex`, `bthm_merged.tex`): false for every `c < 1`.
  A predictor with `T = 5.3 n ≪ n²/ln n` and `O(log n)` bits attains R² = 1.
* **Dichotomy conjecture**: false.  The canonical rule is *exactly predictable in O(n) reads*,
  so it belongs to the cheap end; there is no quadratic lower bound, and the "0.55 n² / 3.6 n²"
  measurements are an artefact of rescanning the window at every step.
* **Open problem** ("is the meeting-position map computable in O(n) reads with O(log n) bits?"):
  yes.  The answer is yes by the incremental construction above.
* The "why we expect no middle case" paragraph's premise — "with logarithmic workspace the
  window must be rebuilt after each cut" — is the specific false step.
* The `_feat_cost.py` "no middle case" argument survives only in the trivial sense that the
  canonical rule is not in the middle: it is at the cheap end.

## Files

| file | role |
|---|---|
| `_bz_inc.py` | incremental simulator + brute-force integer oracle + `check` mode |
| `_bz_inc2.py` | bit-exact reproduction of the published (float-tie) machine |
| `_bz_inc3.py` | read cost on the meeting set, one n per run |
| `_bz_final.py` | the table above: reads, budget, agreement, R², max step reads |
| `_bz_bits.py` | measured workspace in bits |
| `_bz_base.py` | ac₁ baseline (reproduces the manuscript's 0.683) |
| `_bz_feat.py` | cheap-feature ceiling (0.72), greedy nested selection |

Everything runs on numpy/scikit-learn; no other dependency.
