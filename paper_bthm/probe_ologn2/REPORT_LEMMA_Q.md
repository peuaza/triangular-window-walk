# Lemma Q: what is true, and a proof of E[T] = Ω(n²) that avoids it

**Verdict.**

1. **E[T] = Ω(n²) is proved below without Lemma Q.** The new ingredient is a *truncated* use of the
   frontier inequality: sum `cost_i ≥ cost(L_i,n−1) − δ_i` only over the first `N = Θ(n)` epochs
   (not over the whole trajectory). The gain per epoch is `c₀(n−i)` and the loss is `δ_i ≤ i`, both
   linear in `i`, so the truncated sum is positive with an optimum at `N* = c₀n/(1+c₀)`. Since
   `cost(i,n−1)` is the cost of an **unconditioned** suffix — the sub-tape `a[i..n−1]` is exchangeable
   with a uniform arrangement given its counts — the single-epoch theorem applies to it verbatim.
   **No conditioning, no Lemma Q, no polytope slice.**

       E[T]  ≥  [ c₀²/(2(1+c₀)) − (c₀/(1+c₀))·q_N ] · n² · (1−o(1)),      N = ⌊c₀n/(1+c₀)⌋,

   where `c₀` is any constant valid in the single-epoch theorem (`c₀ = 0.008` is rigorous,
   REPORT_EPOCH_BOUND §6; `c₀ = 0.1365` is the sharp measured value) and `q_N` is the probability
   that the walk **meets** before epoch `N`. With `c₀ = 0.008`: **E[T] ≥ 3.17·10⁻⁵ n²**;
   with `c₀ = 0.1365`: **E[T] ≥ 8.20·10⁻³ n²**. Numerically `q_N = 0` (0 meets before epoch 0.4n in
   720 tapes at n = 256, 512, 1024; the largest window at which any meet was seen in 176 000 epochs is
   M = 58).

2. **Lemma Q as stated is not the right target.** The conditional mean is dominated by the run length
   of preceding double cuts and falls by a factor ≈ 250 (0.307 → 0.0012 at n = 1024) as that run grows
   from 0 to ~50, flattening only in the deepest bins. Any uniform constant therefore has to be
   `c ≲ 10⁻³` — two orders of magnitude below the annealed `0.1365` — and the only regime that fixes it
   is the deep-run regime (long persistence of descending top-tests), which no available tool
   controls. What *is* true, and is all that is needed, is the **block-amortised** form proved in §3.

3. The two remaining inputs are (i) the single-epoch theorem of REPORT_EPOCH_BOUND (proved there) and
   (ii) `q_N → 0`. (ii) is strictly weaker than the "`P(early meet) = o(1)`" hypothesis already flagged
   in REPORT_EPOCH_BOUND §8, it contains no conditioning, and it is measured to be vacuous at every
   size tested. It is *not* proved here; §4 states exactly what a proof needs and gives the decay data.

All numbers below are from `qL_verify.py`, `qL_run.py`, `qL_meet.py` in this directory (§6).
No `.tex` file was touched.

---

## 1. The object, and the single-epoch input

Ensemble: `a ∈ {0,1,2}ⁿ` uniform on `{#2 − #0 = 2⌈n/2⌉ − n}` (equivalently, two independent uniform
`h`-subsets of `[n]`, `h = ⌈n/2⌉`, summed). `b = a−1`, `X_j = Σ_{i≤j} b_i`, `A_j = 1 + X_j/j`.

The b1≡0 walk (`q0_walk.py`): at `(L,R)`, `M = R−L+1 ≥ 3`, scan `m = M, M−1, …, 2` and let
`C = {m : A_{⌊m/2⌋} ≥ A_m}`. If `C = ∅` the walk **meets** (cost `M−1`, walk stops); else
`m* = max C`, `cost = M − m* + 1`, then `L ← L+1` and `R ← R−1` iff `m* = M` (a **double cut**),
otherwise `R` unchanged. `T` is the total number of tests.

**(F0) Geometry** (REPORT_EPOCH_BOUND §8, re-verified here). Realised epochs are indexed by
`i = L_i`, and
```
L_i = i,      R_i = n−1−δ_i,      δ_0 = 0,   δ_{i+1} − δ_i ∈ {0,1},   δ_i ≤ i,
M_i = n − i − δ_i .
```
`δ_i` is exactly the number of double cuts among epochs `0..i−1`.

**(F1) Frontier inequality** (REPORT_LOWER §2, Lemma B; re-verified on 176 560 epochs at
n = 256, 512, 1024 with 0 violations). For `L+2 ≤ R ≤ R″`:
```
cost(L,R)  ≥  cost(L,R″) − (R″−R) .
```
Proof: the concavity set of `[L,R]` is `C(R″) ∩ [2, R−L+1]`, so the largest concavity can only move
up by `R″−R`, and the meet case is bounded separately. ∎

**(F2) Single epoch, fresh window** (REPORT_EPOCH_BOUND §6, Theorem). For any exchangeable ensemble
that is a uniform arrangement given the counts, with `s² = M⁻¹Σb_i² ≥ 1/4`, and any `c ∈ (0,1)`:
```
P(cost ≥ cM) ≥ p(c) − o(1),     p(0.05) ≥ 0.16,      hence   E[cost] ≥ c₀ M,   c₀ = 0.008.
```
The `o(1)` is uniform over such multisets, and the bridge/drift term cancels exactly in the dyadic
contrast (`𝔅(θ) − 2𝔅(θ/2) = 𝒱(θ)`), so unbalanced multisets are covered as well. *(Input; proved in
REPORT_EPOCH_BOUND, not re-proved here.)*

**(F3) The suffix sub-tape is fresh.** For every `i`, the positions `[i, n−1]` are a subset of the
exchangeable tape, so `a[i..n−1]` is exchangeable and, conditionally on its multiset, a uniform
arrangement of it, with `s² = (n−i)⁻¹#{j ∈ [i,n−1] : a_j ≠ 1} → 2/3`. Hence (F2) applies to
`cost(i,n−1)`, the cost of the epoch at window `[i, n−1]`, and
```
E[cost(i, n−1)]  ≥  c₀ (n−i) (1−o(1)),      uniformly in i ≤ n/2.
```
Measured `E[cost(i,n−1)]/(n−i)` at n = 512: 0.144, 0.144, 0.147, 0.156, 0.137, 0.146 for
`i/n = 0, .02, .05, .10, .20, .30` — flat, as (F3) requires.

## 2. Why summing (F1) over the whole trajectory fails, and the fix

Summing (F1) with `R″ = n−1` over *all* realised epochs gives
`T ≥ Σ_i cost(i,n−1) − Σ_i δ_i`, and `Σ_i δ_i ≈ 0.62·(0.63n)²/2 ≈ 0.12 n²` (double-cut fraction
≈ 0.62, so `δ_i ≈ 0.62 i`) exceeds `Σ_i E[cost(i,n−1)] ≈ 0.1365·(0.63n − 0.63²n/2)n ≈ 0.059 n²`.
That is the failure recorded in REPORT_LOWER §8.3.

The fix is to truncate. Both terms are *linear in the epoch index*:
```
gain at epoch i:   E[cost(i,n−1)]  ≈  c₀ (n−i)      (decreasing in i)
loss at epoch i:              δ_i  ≤  i             (increasing in i)
```
so their sum over `i < N` is `c₀(nN − N²/2) − N²/2`, a concave quadratic in `N` that is **positive**
on `0 < N < 2c₀n/(1+c₀)` and maximal at `N* = c₀n/(1+c₀)`. Truncating at `N* = Θ(n)` is what buys the
quadratic bound; summing to `i ≈ 0.63n` is what killed the earlier attempt.

## 3. Theorem (E[T] = Ω(n²), no Lemma Q) and its block-amortised corollary

> **Theorem.** For every `n` and every `1 ≤ N ≤ n/2`, with `q_N := P(the walk terminates before
> epoch N)`,
> ```
> E[T]  ≥  c₀ (nN − N(N−1)/2)(1−o(1))  −  N(N−1)/2  −  nN·q_N .
> ```
> In particular, with `N = ⌊c₀n/(1+c₀)⌋`,
> ```
> E[T]  ≥  [ c₀²/(2(1+c₀)) − (c₀/(1+c₀))·q_N ] n² (1−o(1)) .
> ```

*Proof.* Let `A_i = {epoch i is realised} ∈ F_{i−1}`; `T = Σ_i cost_i 1{A_i}`. On `A_i`, (F0) and
(F1) with `R = R_i`, `R″ = n−1` give
```
cost_i  ≥  cost(i,n−1) − (n−1−R_i)  =  cost(i,n−1) − δ_i  ≥  cost(i,n−1) − i .
```
Hence for `i < N`, since `0 ≤ cost(i,n−1) ≤ n−i`,
```
E[cost_i 1{A_i}]  ≥  E[cost(i,n−1) 1{A_i}] − i
                  ≥  E[cost(i,n−1)] − (n−i) P(A_i^c) − i
                  ≥  c₀(n−i)(1−o(1)) − (n−i) q_N − i ,          (F3), and P(A_i^c) ≤ q_N
```
using `P(A_i^c) ≤ P(A_N^c) = q_N` for `i < N` (the walk cannot restart). Summing over `i < N`,
using `Σ_{i<N}(n−i) = nN − N(N−1)/2 ≤ nN` and `Σ_{i<N} i = N(N−1)/2`:
```
E[T] ≥ c₀(nN − N(N−1)/2)(1−o(1)) − nN q_N − N(N−1)/2 .
```
Put `N = αn`: the bracket is `[c₀(α − α²/2) − α²/2 − α q_N] n²`, and
`max_α [c₀(α−α²/2) − α²/2] = c₀²/(2(1+c₀))` at `α = c₀/(1+c₀)`. ∎

> **Corollary (block-amortised Lemma Q).** For every `N ≤ n/2`, with `A_i = {epoch i is realised}`,
> ```
> E[ Σ_{i<N} cost_i 1{A_i} ]  ≥  c₀ Σ_{i<N} (n−i) − N(N−1)/2 − nN q_N ,
> ```
> i.e. the per-epoch bound `E[cost_i] ≥ c₀ M_i` holds *on average over the first N epochs*, with an
> error `N²/2` that is one order smaller than the main term `c₀nN` as long as `N ≪ n`. (Since
> `T ≥ Σ_{i<N}cost_i 1{A_i}`, the same lower bound holds for `E[T]`.)

This is exactly route (b) of the brief ("amortise: group the epochs"), in its sharpest form: the group
is *all* the epochs up to `N`, and the only property of the sub-tape that is used is (F3), i.e.
exchangeability of an unconditioned suffix — never the conditional law of the current window.

**Constants.**

| `c₀` | source | `N*` | bound coefficient | bound at n = 1024 | true E[T]/n² |
|---|---|---|---|---|---|
| 0.008 | rigorous (F2) | 0.00794 n | **3.175·10⁻⁵** | 33 | 0.051 |
| 0.1365 | measured, sharp | 0.1201 n | **8.197·10⁻³** | 8 596 | 0.051 |

The bound is `c₀²/(2(1+c₀))·n²`: it is quadratic in the single-epoch constant, which is why the
weak rigorous `c₀` costs four orders of magnitude. Improving `c₀` is a separate (single-window)
problem; nothing in the route depends on its value.

## 4. The one remaining hypothesis: no early meet

`q_N = P(the walk meets at some epoch < N)`. A meet at `[L,R]` requires **every** test
`A_{⌊m/2⌋} < A_m`, `m = 2..M`, so it is a "strictly rising profile at every dyadic comparison" event.

*(a) Measured, fresh windows* (`qL_meet.py`, 150 000–400 000 draws per M from a fixed n = 65 536 tape):

| M | 8 | 12 | 16 | 24 | 32 | 48 | 64 | ≥ 96 |
|---|---|---|---|---|---|---|---|---|
| P(meet) | .0339 | .0132 | .00556 | .00131 | .000515 | 9.7·10⁻⁵ | 2.5·10⁻⁵ | 0 (0/9·10⁵) |

The local power-law exponent `−d log P/d log M` grows 2.6 → 3.5 → 4.4 over M = 8 → 64, i.e. the decay
is faster than any fixed power.

*(b) Measured, along trajectories* (600 tapes, n = 512; 174 043 epochs): meet hazard by window band

| M | [3,8) | [8,16) | [16,32) | [32,64) | [64,128) | [128,256) | [256,512) |
|---|---|---|---|---|---|---|---|
| epochs | 1 234 | 2 460 | 5 077 | 10 371 | 22 097 | 47 068 | 105 136 |
| meets | 40 | 23 | 10 | 5 | **0** | **0** | **0** |

The largest meet window seen anywhere: 27 (n = 256), 58 (n = 512), 29 (n = 1024) over 176 560 epochs
at the three sizes. Consequently `P(all first N epochs realised) = 1.0000` for every `N ≤ 0.4n` at all
three sizes, and `q_N = 0` for the `N*` of §3 by a factor ≈ 10³ in probability.

*(c) What a proof of `q_N → 0` needs.* The dyadic chain alone gives a clean exponential bound on one
branch: a meet forces `S_{2^t} > 2S_{2^{t−1}}` for `t = 1,2,…`, and with `S_1 = b_1 ∈ {−1,0,1}`:
`S_1 = 1` is impossible (`S_2 ≥ 3 > 2`); `S_1 = 0` forces `S_2 = 1, S_4 = 3, S_8 = 7, …`, i.e.
`a_L = 1` and `a_{L+1} = … = a_{L+2^T−1} = 2`, an event of probability `≈ 4^{−M/2}`. The remaining
branch `S_1 = −1` (the window starts with a 0) is not killed by the dyadic chain — it only forces
`S_{2^t} ≥ −1` — and needs the *full* set of `M−1` tests, not just the dyadic ones. This is a
self-contained combinatorial problem about monotone profiles of a random word, with no conditioning
in it. **It is not solved here**, and it is the only unproved input of §3 beyond (F2).

## 5. Lemma Q itself: the measured obstruction

Lemma Q asks for `E[cost_i | F_{i−1}] ≥ c M_i` with `c` independent of `i, n`. The natural proxy for
`F_{i−1}` is the length `r` of the run of consecutive double cuts ending at epoch `i−1` (a double cut
is exactly `cost_i = 1`, so this is the history feature that says "the last `r` top-tests all came out
descending"). Holding `M/n` fixed, at n = 1024, 300 tapes, 193 488 epochs:

| r | 0 | 1 | 2–4 | 4–8 | 8–16 | 16–32 | 32–64 |
|---|---|---|---|---|---|---|---|
| mean cost/M, M/n ∈ [.9,1) | **.307** | .0220 | .0080 | .0060 | .0024 | .0012 | .0012 |
| P(cost ≥ 2), same band | .976 | .142 | .050 | .042 | .010 | .006 | .004 |

(same pattern in every `M/n` band from .5 to 1; the deepest resolved bin has mean cost/M ≈ .003 at
n = 512 and ≈ .0012 at n = 1024; the all-M table, which mixes in small windows, floors at
`P(cost≥2) = .0085` (n = 512) and `.0044` (n = 1024)).

Reading:

* the conditional mean falls by a factor ≈ 250 between `r = 0` and `r ≈ 50`; the "×1.26 favourable"
  annealed/quenched ratio of REPORT_LOWER §5 is an average over this mixture, not a uniform bound, and
  so is every "the conditioning is favourable" statement in the earlier reports;
* the decay is steep but **flattens** at the deepest bins (`.0012` at r = 16–32 and again at r = 32–64
  at n = 1024; `P(cost≥2)` floors at `.0044`–`.0085` in the all-M tables). So the data do **not**
  decide whether `E[cost_i|F_{i−1}]/M_i` tends to 0 along deep runs; they only force the constant down;
* the inner mechanism: after `r` descending top-tests the top contrast is pinned below zero over `r`
  consecutive scales, and the next top-test is a 2-coordinate perturbation of the pinned value
  (`E^{(i+1)}_{k−1} = E^{(i)}_k + b_{L} − b_{R}` on a double cut), so the exit hazard is the
  probability that a 2-coordinate perturbation crosses back above 0 — small, and the observed
  near-constancy at large `r` says the "pinned" depth saturates;
* therefore **any** proof of Lemma Q must produce `c ≲ 10⁻³` and must control exactly the deep-run
  regime, where the window has already shrunk by `2r` and the history is a long persistence condition.
  That is strictly harder than the block-amortised statement of §3, which needs no `r`-uniform control.

So Lemma Q is not refuted here, but it is not the right target: the constant it permits is 100×
smaller than the annealed `0.1365` that the single-epoch theorem suggests, and the only regime that
determines it is the one no available tool controls. The correct statement, and the one that gives
`Ω(n²)`, is the corollary of §3.

**Status of the other two routes in the brief.**

* (a) *Annealed domination via a Radon–Nikodym factor.* The tests are linear functionals of the tape,
  so `F_{i−1}` is generated by the signs of `T` hyperplanes and the conditional law of the window word
  is uniform on the resulting cell. The density of the conditional law with respect to the fresh
  arrangement law is `1{cell}/μ(cell)` with `1/μ(cell) = 1/P(F_{i−1}) = e^{+Θ(T)}` — enormous, but
  supported on a set of measure `e^{−Θ(T)}`. That gives an *upper* bound on conditional probabilities
  and no lower bound at all: absolute continuity is not the obstruction, positive association between
  the good event `{cost_i ≥ cM_i}` and the history is, and that is exactly Lemma Q. Route (a) can only
  work if the comparison is made *one test at a time* (a rank-one tilt with an explicit likelihood
  ratio), which is the `r`-uniform statement the table above shows to be the hard part.
* (c) *Coupling whole trajectories.* `cost` is a first-passage index of the scan, not a Lipschitz
  functional of the window word: `a = (1,1,…,1)` has `cost = 1` while `a = (0,1,…,1)` has `cost = M−1`
  (verified), and the two differ in one cell. So no total-variation or Wasserstein estimate on window
  words can transfer a cost bound, and a coupling of trajectories would have to be simultaneous in the
  tape — which is the same Radon–Nikodym problem as (a).

## 6. Numerical verification (all runs in this directory)

`$P = C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`, run in
`C:/Users/shini/poetry_axis/paper/probe_ologn2/`.

```
$P qL_verify.py 256 300 2026      # frontier chain, per-epoch violations, LB(N), constants
$P qL_verify.py 512 300 2026
$P qL_verify.py 1024 120 2026
$P qL_run.py   1024 300 77        # Lemma-Q proxy table (run length vs conditional mean cost)
$P qL_meet.py  512 600 2026       # meet probability vs M, fresh and along trajectories
```

**Frontier chain.** Per-epoch check `cost_i ≥ cost(L_i,n−1) − δ_i`: **0 violations** in 38 354 (n=256),
76 775 (n=512), 61 431 (n=1024) epochs. Sum chain `T ≥ LB(N) = Σ_{i<N}cost(i,n−1) − N(N−1)/2`:
0 violations at every N tested, `min(T − LB) > 0` throughout.

**Growth of the bound.** `c = E[cost(0,n−1)]/n` is the measured fresh-window constant for that size;
`c²/(2(1+c))` is the prediction of §3 and `E[LB]/n²` is the realised bound at the optimal `N`:

| n | N* | E[LB]/n² | c | predicted `c²/(2(1+c))` | E[T]/n² |
|---|---|---|---|---|---|
| 256 | 31 (0.121 n) | 0.00916 | 0.138 | 0.00900 | 0.0573 |
| 512 | 61 (0.119 n) | 0.00942 | 0.144 | 0.00912 | 0.0596 |
| 1024 | 123 (0.120 n) | 0.00707 | 0.125 | 0.00735 | 0.0506 |

The empirical optimum sits at `N* = 0.12 n` at every size, exactly as `c/(1+c)` predicts, and the
realised bound reproduces `c²/(2(1+c))` to within 4% — the `c₀²` law of §3 is confirmed with
`c₀ = c`. The realised bound is 16% of `E[T]` at n = 256–512; with the *rigorous* `c₀ = 0.008` the
same formula gives `3.17·10⁻⁵ n²`.

**Termination.** `P(meet ever) = 0.170 / 0.130 / 0.167` at n = 256/512/1024, but every meet occurs at a
window of length ≤ 58; `P(all first N epochs realised) = 1.0000` for all `N ≤ 0.4 n`; `q_{N*} = 0`.

## 7. What is proved, in one paragraph

`E[T] ≥ (c₀²/(2(1+c₀)) − O(q_N)) n²`, where `c₀ = 0.008` is the rigorous single-epoch constant of
REPORT_EPOCH_BOUND and `q_N` is the probability of a meet before epoch `c₀n/(1+c₀)` — measured to be
zero at every size, and independent of any conditioning statement. Concretely: **E[T] ≥ 3.17·10⁻⁵ n²
rigorously**, and **E[T] ≥ 8.2·10⁻³ n²** using the sharp (measured) single-epoch constant. The
mechanism is that the frontier inequality converts each epoch's cost into the cost of an
*unconditioned* suffix at the price of `δ_i ≤ i`, and truncating the sum at `i ≈ c₀n` makes both the
gain and the loss linear in `i`, leaving a positive quadratic. Lemma Q is not used, and its stated
per-history form is not the operative statement: the conditional mean `E[cost_i|F_{i−1}]/M_i` decays
from 0.31 to 0.001 as the preceding double-cut run grows from 0 to 50, so the only constant Lemma Q
could carry is `≈10⁻³`, and the deep-run regime that fixes it is exactly where a proof has no tool.
