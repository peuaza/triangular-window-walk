# Lower bound for the b1≡0 concavity-scan walk — what is proved, what is only measured

Question: is the number T of b2 tests of the b1≡0 (largest-concavity) rule Ω(n²) in
expectation over the ensemble (n even: #0 = #2; n odd: #0 = #2+1)?

**Verdict.** The bound is **true as far as I can measure** and I believe it
asymptotically, but **it is not proved here**. It is reduced here to one clean
single-statistic statement (§5, Lemma L″) instead of a statement about the entire
trajectory. What is *proved*: the exact algebra of a test, a frontier inequality, an
exact recursion (the trajectory is a sum of nested sub-tape costs), an exact half-swap
identity, and a **linear** lower bound with an explicit constant.

Three findings that matter for the paper and are independent of the proof:

1. **The cost is governed by one global statistic**, the half-difference
   Δ(a) = mean(a[n/2..n−1]) − mean(a[0..n/2−1]). E[T/n² | Δ ≥ 0] = 0.090–0.104,
   E[T/n² | Δ < 0] = 0.013–0.016, corr(T/n², Δ/σ_Δ) = **+0.74**; 88% of E[T] comes
   from the Δ ≥ 0 half.
2. **The trajectory is not Θ(n²) pathwise.** P(T ≤ n) = 0.17 (n = 1024, 2000 tapes);
   on that set 99% of epochs cost exactly 1 (all double cuts), i.e. T ≈ 0.6 n.
3. **T does not concentrate**: at n = 1024, mean 0.054 n², median 0.025 n², sd
   0.061 n² (sd ≈ mean). "The rule asks Θ(n²) questions" is an expectation statement,
   not a typical-case statement.

All numbers below come from the scripts in this directory (§9). No .tex file touched.

---

## 1. Proved: exact algebra of one test

For a window [L,R] put b = a−1, X_j = Σ_{i=L}^{L+j−1} b_i (so A_j := mean(a[L..L+j−1])
= 1 + X_j/j), and

    D_k := X_{2k} − 2 X_k = sum(a[L+k..L+2k−1]) − sum(a[L..L+k−1]).

**Lemma A.** For m = 2k: b2([L,R]) = 1 ⟺ D_k > 0.
For m = 2k+1: b2 = 1 ⟺ k·(X_{2k+1} − 2X_k) > X_k.

Proof: b2 = 1 is the strict form of Lemma A of REPORT.md (A_{⌊m/2⌋} < A_m). For m = 2k
it is X_{2k}/(2k) > X_k/k ⟺ D_k > 0; for m = 2k+1 multiply through by k(2k+1). ∎

Verified exhaustively: all 3^n tapes, n ≤ 10, all prefix windows = **752,868
(tape,m) pairs, 0 mismatches** (`q10_identity.log`, `q7_quench.py` step (1)).

So one epoch is a first-passage problem for the process (D_k) scanned in decreasing k;
the odd-m tests are the same test with an O(1) + X_k/k threshold.

## 2. Proved: frontier inequality

**Lemma B.** For L + 2 ≤ R ≤ R″ ≤ n−1:  cost(L,R) ≥ cost(L,R″) − (R″ − R).

Proof. The concavity set of [L,R] is C(R) = {m ∈ [2, R−L+1] : A_{⌊m/2⌋} ≥ A_m}, and A_j
depends only on (L,j), so C(R) = C(R″) ∩ [2, R−L+1]. If C(R) = ∅ then cost(L,R) = R−L
and cost(L,R″) ≤ R″−L. If C(R) ≠ ∅ let m = max C(R): either max C(R″) = m (equality),
or max C(R″) > R−L+1, in which case cost(L,R″) − (R″−R) = R−L+2−max C(R″) ≤ 2 ≤ cost(L,R). ∎

Verified on 54,400 (L,R,R″) triples at n = 10,11,12, **0 violations** (`q7_quench.py`).

Lemma B bounds the loss caused by the right header having moved; but that loss is
δ_L = #double cuts before L and Σ_L δ_L ≈ n²/16 exceeds the total cost, so Lemma B
alone cannot convert suffix-cost bounds into trajectory bounds.

## 3. Proved: exact recursion, and an exact half-swap identity

**Lemma C.** Let the first epoch of a length-n tape find largest concavity m* (None if
it meets) with cost c₀, and put R₁ = n−2 if m* = n−1 (double cut), else R₁ = n−1. Then

    T(a) = c₀ + T(a[1..R₁]),

the second term being the same b1≡0 walk on the sub-tape a[1..R₁].

Proof: after the epoch the headers restart at (1,R₁); every later test is b2([ℓ,r])
with ℓ ≥ 1, r ≤ R₁, depending only on a[ℓ..r] ⊆ a[1..R₁]. ∎

Verified on 200 tapes at n = 200, 0 violations (`q8_hist.py`). Iterating:
T = Σ_{L<L*} cost_L over windows [L,R_L], R_L = n−1−δ_L; and L* ∈ [(n−1)/2, n−1]
whenever the run never meets (measured L*/n = 0.62–0.63).

**Lemma D (half-swap identity).** For n even let σ exchange the two halves,
(σa)_i = a_{i+n/2} for i < n/2. σ is a permutation of positions, hence preserves the
ensemble, and Δ(σa) = −Δ(a). Therefore exactly

    E[T] = E[ (T(a) + T(σa)) · 1{Δ(a) > 0} ] + E[ T · 1{Δ(a) = 0} ].

Verified **exactly** by full enumeration at n = 6, 8, 10, 12 (n = 6: E[T] = 4.2057,
first term 3.2128, second 0.9929, sum 4.2057; n = 12: 14.9054 = 12.0578 + 2.8476).
P(Δ = 0) = 0.35, 0.33, 0.29, 0.27 at n = 6, 8, 10, 12 and 0.10 at n = 128, 0.028 at
n = 1024 (≈ c n^{−1/2}), so the Δ = 0 term is o(n²) and asymptotically
E[T] ≈ E[(T + Tσ)1{Δ > 0}].

## 4. The reduction that matters: the drift Δ controls the trajectory

Measured, n even, 600–2000 tapes per n (`q11_drift.py`, `q11.log`):

| n | P(Δ≥0) | E[T/n² \| Δ≥0] | E[T/n² \| Δ<0] |
|---|---|---|---|
| 128 | 0.535 | **0.1038** | 0.0156 |
| 256 | 0.572 | **0.0949** | 0.0145 |
| 512 | 0.555 | **0.0952** | 0.0134 |
| 1024 | 0.517 | **0.0904** | 0.0154 |
| 2048 | 0.503 | **0.0910** | 0.0120 |
| 4096 | 0.545 | **0.0943** | 0.0140 |

So 88% of E[T] is carried by the Δ ≥ 0 half, and the conditional value is **stable**
at 0.090–0.104 over a 32× range of n, while the Δ < 0 half sits at ~0.013.
Also corr(T/n², Δ/σ_Δ) = +0.738 (n = 512, 600 tapes) and
P(T ≥ 0.01 n² | Δ ≥ 0) = 0.93, 0.93, 0.89, 0.87, 0.92, 0.88 for n = 128…4096.

(The sampled "identity" column in `q11.log` — E[(T+Tσ)1{Δ≥0}] vs E[T] — shows ratios
0.94–0.98, exactly the omitted Δ = 0 term: P(Δ=0) = 0.10, 0.028, 0.031, 0.016 at
n = 128, 1024, 2048, 4096. The exact identity including that term was verified by full
enumeration at n ≤ 12.)

The mechanism is visible in the trajectory: on Δ ≥ 0 tapes **51.8% of the epochs'
windows have a positive window drift** (mean of the window's right half > left half),
against 20.8% on Δ < 0 tapes; trajectory mean cost/M is 0.189 vs 0.066 (n = 512, 300
tapes). The drift is approximately inherited by the sub-windows, which is why one
global sign controls the whole trajectory. Consistently, corr(T, T∘σ) = −0.54 ≈
−(0.74)²: the two orientations of the same tape are anti-correlated because σ flips Δ,
and corr(T, T∘rev) = −0.54 as well.

## 5. The one missing lemma, and what it would give

**Lemma L″ (open).** E[T | Δ ≥ 0] ≥ c n² for a constant c > 0.

With E[T] ≥ P(Δ ≥ 0)·E[T | Δ ≥ 0] and P(Δ ≥ 0) → 1/2, Lemma L″ gives
E[T] ≥ (c/2) n² = Ω(n²); measured c ≈ 0.09–0.10, so E[T] ≳ 0.045 n².

Equivalently (Lemma C): the sub-tape seen after an epoch must not be worse than a
fresh tape of the same length — the annealed→quenched step. Evidence that it holds,
n = 512, 250 tapes (`q7_quench.py`): comparing each epoch's realised cost with the
fresh-window mean ĉ(M) for the same window length,

    mean realised T = 16,429   vs   mean Σ_L ĉ(M_L) = 13,043   → ratio 1.26

i.e. the conditioning is *favourable* by ~26%. Bucketed by M, the trajectory's mean
cost/M is 0.065 (M 32–64), 0.113 (64–128), 0.137 (128–256), 0.200 (256–512) against
fresh-window values 0.16, 0.16, 0.15, 0.14.

Why it is hard: the conditioning event is a system of linear inequalities over the
whole tape (each early answer compares two block sums of a[0..n−1]), so the conditional
law of the sub-tape is a polytope slice of the uniform ensemble; exchangeability does
not survive, and I found no argument that controls the slice. The Δ-reduction is the
closest I got: it converts "control the whole trajectory" into "control the
positive-drift half", a single linear event.

> **Update (REPORT_LEMMA_Q.md).** Neither Lemma L″ nor Lemma Q is needed. Summing the frontier
> inequality `cost_i ≥ cost(L_i,n−1) − δ_i` only over the first `N = Θ(n)` epochs (the loss `δ_i ≤ i`
> and the gain `c₀(n−i)` are both linear in `i`, so the truncated sum is a concave quadratic in `N`)
> gives `E[T] ≥ [c₀²/(2(1+c₀)) − O(q_N)]n²`, i.e. `3.17·10⁻⁵ n²` rigorously and `8.2·10⁻³ n²` with
> the sharp single-epoch constant. The reason it works is that `cost(i,n−1)` is the cost of an
> **unconditioned** suffix `a[i..n−1]`, which is exchangeable, so the single-epoch theorem applies to
> it directly — the slice is never touched. The only residual hypothesis is `P(meet before epoch N) → 0`
> (measured 0: no meet at any window ≥ 64 in 174 043 epochs).

## 6. Proved modulo a standard FCLT: the per-window law (source of the constant)

**Proposition E.** For a window of length M whose content is a uniform random
arrangement, cost/M ⇒ 1 − 2θ*, where

    θ* = sup{θ ∈ (0,1/2] : Z(θ) ≤ 0},   Z(θ) = W(2θ)/(2θ) − W(θ)/θ,

W a standard Brownian motion; hence E[cost] = Θ(M) and E[1−2θ*] = **0.1365**
(20,000 simulated BM paths, grid 1024, `q5_brownian_quarters.py`). The limit is the
FCLT for the triangular array D_k; I did not write the tightness proof, so this is
recorded as standard-but-not-fully-written and the constant as measured-and-matched.

Measured: first epoch of a length-n tape (= a fresh window of length n),
`q9_first.py`, 2000 tapes per n (600 for n ≥ 4096), seed 77:

| n | 64 | 128 | 256 | 512 | 1024 | 2048 | 4096 | 8192 |
|---|---|---|---|---|---|---|---|---|
| E[cost₀]/n | 0.158 | 0.151 | 0.140 | 0.138 | 0.129 | 0.138 | 0.128 | 0.143 |

Fresh random windows of length M (`q3_fresh.py`, up to 125,000 windows per M):
E[cost]/M = 0.217 (16), 0.179 (32), 0.160 (64), 0.162 (128), 0.156 (256), 0.128 (512),
0.121 (1024); P(cost/M ≥ 0.1) = 0.36–0.44; P(window meets) = 0.0042 (M=16), 0.0007
(M=32), 0.0000 for M ≥ 64. The per-window constant shows no decay with M over the
range that matters, so the only n-dependence in E[T]/n² is how many windows the
trajectory visits.

**Corollary (proved given Proposition E).** E[T] ≥ E[cost₀] ≥ 0.128 n − o(n). This is
the strongest *proved* lower bound: linear. (A cruder c n is elementary without
Proposition E: the first epoch alone, since P(b2 = 1 at m = n−1) = 1/2 − O(n^{−1/2})
by the half-swap symmetry.)

## 7. Measured: the answer, the decomposition, and the caveats

Exact enumeration (`q2_exact.py`, all tapes):

| n | tapes | E[T] | E[T]/n² | P(T=n−1) | P(T≤2n) | mean epochs | double-cut frac |
|---|---|---|---|---|---|---|---|
| 6 | 141 | 4.2 | 0.1168 | 0.1135 | 1.000 | 2.24 | 0.585 |
| 8 | 1,107 | 7.2 | 0.1121 | 0.0623 | 0.936 | 3.31 | 0.583 |
| 10 | 8,953 | 10.8 | 0.1077 | 0.0321 | 0.830 | 4.47 | 0.582 |
| 12 | 73,789 | 14.9 | 0.1035 | 0.0279 | 0.761 | 5.65 | 0.583 |
| 13 | 201,643 | 18.0 | 0.1068 | 0.0253 | 0.713 | 6.55 | 0.581 |

Sampled (`q4_scale.py` 400 tapes/n seed 2026; `q8_hist.py` 2000 tapes at n = 1024;
`q9_first.py` n = 8192):

| n | 256 | 512 | 1024 | 1024 (2000 tapes) | 2048 | 4096 |
|---|---|---|---|---|---|---|
| E[T]/n² | 0.0556 | 0.0556 | 0.0477 | **0.0540 ± 0.0014** | 0.0527 | 0.0591 |
| ±se | .0031 | .0030 | .0028 | .0014 | .0031 | .0032 |
| median T/n² | 0.026 | 0.029 | 0.019 | 0.025 | 0.022 | 0.033 |
| P(T≥0.1n²) | 0.247 | 0.255 | 0.185 | 0.250 | 0.230 | 0.263 |
| P(T≤n) | 0.263 | 0.180 | 0.210 | 0.173 | 0.158 | 0.130 |
| epochs/n | 0.623 | 0.634 | 0.624 | — | 0.634 | 0.642 |
| double-cut/ep | 0.64 | 0.62 | 0.65 | — | 0.63 | 0.61 |

**Decomposition.** With M_L the actual window lengths (300 tapes per n, seed 555),

    E[T]/n² = [Σ_L M_L/n²] · [mean cost/M along the trajectory]
    n=256 : 0.3438 · 0.1845 = 0.0634
    n=1024: 0.3297 · 0.1614 = 0.0532
    n=4096: 0.3312 · 0.1545 = 0.0512

Σ_L M_L/n² ≈ 1/3 is geometric (≈0.63 n epochs, mean window ≈0.53 n). The asymptotics
sit in the second factor: 0.185 → 0.161 → 0.155, heading towards the fresh-window
value 0.1365 of Proposition E, which would give E[T]/n² → 0.0455. The data
(0.048–0.059 for n = 1024…4096) cannot by themselves exclude a very slow decay such
as 1/log n, but the ingredient that would have to decay — the fresh-window constant —
shows none for M = 64…1024.

**Bimodality** (`q8_hist.py`, n = 1024, 2000 tapes, T/n²): 6.8% in [0, 0.0005)
(T ≈ n/2, all double cuts), 11.1% in [0.0005, 0.001), then a broad spread (6–12% per
octave) up to 0.25, with 25% above 0.1. The T ≤ n tapes are *not* "no-concavity"
tapes (that event needs T = n−1 and has P = 0.028 at n = 12 and < 5e-4 for n ≥ 64):
they are tapes on which almost every top test fails, so every epoch double-cuts; the
cost-1 fraction is 0.97–0.99 there versus 0.52–0.56 on the rest.

## 8. Routes that do not work (all checked numerically)

1. **Four-quarter-block sufficient condition** ("quarter means strictly increasing ⟹
   b2 = 1 for all m > M/2") is **false**: exhaustive over 3^M (`q5`, `q6_diag.py`)
   gives counterexamples — M = 10: 120 of 677 quarter-monotone tapes have their
   largest concavity at m ∈ {8,9,10}; M = 11: 577 of 3,018. Reason: a quarter's mean
   does not control a sub-block that straddles the boundary (one cell can flip it).
   It is also far from tight: measured P(cost/M ≥ 0.25) = 0.22–0.35 ≫ 1/24.
2. **Global coarse-monotone structure** ("all block means increasing at scale cn"):
   self-defeating — such a tape has no concavity at all, so the walk meets at L = 0
   with T = n−1 (P(meet) = 0 for M ≥ 64); and requiring it for every visited window
   is a Θ(n)-condition event, hence exponentially unlikely.
3. **Frontier inequality as a bridge to suffix costs**: loses Σ_L δ_L ≈ n²/16 > cost.
4. **Certificate / query-complexity arguments**: capped at n−1 reads (REPORT.md §5);
   they cannot see n².
5. **Reversal or complement symmetry**: corr(T, T∘rev) = −0.54, and min over tapes of
   (T + T∘rev)/2 = 0.005 n² (a monotone tape has T = n−1 and T∘rev = n/2), so no
   two-orientation argument can be uniform.

## 9. Files, commands, seeds

`$P = C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`, run in
`C:/Users/shini/poetry_axis/paper/probe_ologn2/`.

```
$P q0_walk.py             # fast exact b1=0 walk; agreement with the reference on all
                          # tapes n=6..12, deltas AND step counts (0 mismatches)
$P q2_exact.py 6 7 8 9 10 11 12 13      # exhaustive statistics (q2_exact.json)
$P q3_fresh.py                          # fresh-window cost law vs M (q3.log)
$P q4_scale.py 256 512 1024 2048 4096   # E[T]/n^2 scaling (q4.log, q4_scale.json)
$P q5_brownian_quarters.py              # Brownian constant; quarter-lemma counterexamples
$P q6_diag.py                           # where the quarter condition fails
$P q7_quench.py                         # identity, frontier check, quenched vs annealed
$P q8_hist.py                           # recursion check + T/n^2 histogram
$P q9_first.py                          # E[cost_0(n)] and n=8192
$P q11_drift.py 128 256 512 1024 2048   # drift decomposition, half-swap identity
```

Seeds: 2026 (q1/q4 trajectories), 424242 (q8 histogram), 77 (q9 first epoch), 555
(decomposition), 5 (q7 windows), 7 (reference tape for q3/q7), 11 (recursion), 5150
(drift), 29/23/17 (structure correlations). Replication of E[T]/n² at n = 1024:
0.0477 ± 0.0028 (seed 2026, 400 tapes) vs 0.0540 ± 0.0014 (seed 424242, 2000 tapes) —
consistent. Third replication, seed 777, 500 tapes, n = 1024: E[T]/n² = 0.0564,
P(T ≥ 0.1n²) = 0.250, P(T ≤ n) = 0.188, E[T/n² | Δ≥0] = 0.0974,
E[T/n² | Δ<0] = 0.0151, P(Δ≥0) = 0.502.

## 10. Bottom line for the paper

* The O(n²) upper bound stands. The matching lower bound is **not proved**; it is now
  reduced to Lemma L″ (E[T | Δ ≥ 0] ≥ c n²), a statement about the positive-drift half
  of the ensemble, with c measured at 0.09–0.10 and stable in n.
* Safe statements: E[T] ≥ c n (proved); E[T] ≳ 0.045 n² conjectured, with
  E[T]/n² = (1/3)·E[cost/M] and E[cost/M] → 0.1365 measured and matched by the
  Brownian-functional limit.
* Unsafe statements: "the walk asks Θ(n²) questions" as a pathwise/typical claim (17%
  of tapes are linear; median T is half the mean), and any claim that the cost
  concentrates.
