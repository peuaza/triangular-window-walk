# One epoch costs Ω(M): a proof for the b1≡0 concavity scan

**Question.** For the b1≡0 (largest-concavity) rule on a uniform ensemble tape, is the cost of a
single epoch Ω(M) with a constant, i.e. is `P(cost ≥ cM) ≥ p > 0`?

**Answer: yes, and here is a complete proof of that statement**, in the strong form

    P(cost ≥ c M) ≥ p(c) − o(1),      p(c) := P( sup_{θ∈[1-c,1]} 𝒱(θ) < 0 ),
    𝒱(θ) := W(θ) − 2 W(θ/2),          W a standard Brownian motion,

which gives `E[cost_0] ≥ c·p(c)·M` for all M ≥ M₀. The only cited (not re-proved) ingredients are
Hoeffding's combinatorial CLT, Serfling's inequality for sampling without replacement, Dudley's entropy
bound and the Borell–Tsirelson–Ibragimov–Sudakov tail inequality for Gaussian suprema; everything else
is proved below, including a quantitative lower bound p(c) > 0.

| statement | constant | status |
|---|---|---|
| p(0.05) ≥ 0.16, hence **E[cost₀] ≥ 0.008 M** | c·p = 0.0080 | **proved** (Lemma 4, Borell + Dudley) |
| p(0.05)=0.381, p(0.1)=0.328, p(0.2)=0.250, p(0.3)=0.191, p(0.5)=0.091 | max c·p = **0.057** (c=0.3); c·p = 0.0328 at c=0.1 | Brownian simulation, 2·10⁴–4·10⁴ paths, grid-stable to 0.5% |
| P(cost ≥ 0.1M) = 0.3293 ± 0.0086, n=4096, 3000 tapes | ratio discrete/Brownian = **1.003** | measured |
| E[cost₀]/M = 0.130 ± 0.009 (n=4096); 0.114–0.149 over n=64…8192 | — | measured |

The last two rows are the numerical check of the constants; §7 has the tables, §8 states exactly what
is and is not proved about `T = Σ cost_i`.

---

## 1. Notation and the exact form of the event

Tape `a ∈ {0,1,2}^M`, `b = a−1 ∈ {−1,0,1}`, `X_m = Σ_{i≤m} b_i`, `A_m = 1 + X_m/m = mean(a_1..a_m)`,
`k(m) = ⌊m/2⌋`. Recall (REPORT_LOWER §1, verified exhaustively for M ≤ 10):

    b2 = 1  ⟺  A_{k(m)} < A_m  ⟺  X_m − 2X_{k(m)} > 0      (m even: exact; m odd: > X_k/k)

The concavity set is `C = {m ∈ [2,M] : A_{k(m)} ≥ A_m}` (these are the m with b2 = 0), `m* = max C`,
and

    cost = M − m* + 1   ( = M−1 if C = ∅, the "meet" case ).

Two exact restatements that drive everything. Write `B = mean(a[M/2+1..M])`; since
`A_M = (A_{M/2}+B)/2`,

* **cost = 1 ⟺ A_M ≤ A_{M/2} ⟺ B ≤ A_{M/2}** — cheap epochs are exactly the windows whose *later* half
  has the **smaller** mean (descending windows);
* **cost ≥ cM ⟺ C ∩ ((1−c)M, M] = ∅** (up to the `+1` in `cost = M−m*+1`) ⟺ every test in the top
  c-fraction gives b2 = 1 ⟺ `X_m − 2X_{⌊m/2⌋} > 0 for all m ∈ ((1−c)M, M]` (up to the odd-m threshold
  `X_k/k`) — i.e. the *later* half out-means the earlier one at **every** dyadic split in the top
  c-fraction.

So the event we must lower-bound is a **persistence event**: the "later block out-means the earlier
block" inequality, true at the top of the scan, stays true down to `(1−c)M`.

## 2. Lemma 1 (FCLT for the tape)

Fix a multiset 𝔐_M ⊂ {0,1,2} with counts (κ₀,κ₁,κ₂), κ₀ = κ₂, and `s_M² = (1/M)Σ_{v∈𝔐}b_v² = 2κ₀/M`.
Let `a` be a uniform random arrangement of 𝔐_M and `𝒳_M(θ) = X_{⌊θM⌋}/(s_M√M)` (polygonal). Then

    𝒳_M ⇒ 𝔅,    𝔅(θ) = W(θ) − θ W(1)   in C[0,1],

uniformly over all such multisets with `s_M² ≥ s²_min > 0`.

*Proof.*
(a) **Exact covariance.** For 0 ≤ i ≤ j ≤ M,
`Cov(X_i,X_j) = s_M² · i(M−j)/(M−1)`.
Indeed `Var(X_i) = s²i(M−i)/(M−1)` (variance of a WOR sum) and for disjoint index sets A,B,
`Cov(Σ_A,Σ_B) = −|A||B|s²/(M−1)`; hence `E[X_iX_j] = s²i(M−i)/(M−1) − s²i(j−i)/(M−1)`. Consequently
`Cov(𝒳(θ),𝒳(θ')) → θ∧θ' − θθ' = Cov(𝔅(θ),𝔅(θ'))`. *(proved; checked numerically, §7A)*
(b) **Finite-dimensional convergence.** `(X_{⌊θ₁M⌋},…,X_{⌊θ_kM⌋})` is a sum over a random permutation
of a fixed vector; Hoeffding's combinatorial CLT gives joint asymptotic normality, and by (a) the limit
covariance is the bridge's. *(cited: Hoeffding 1951; this is exactly the "random permutation" CLT.)*
(c) **Tightness.** `X_j − X_i` is a sum of j−i draws without replacement, so Serfling's inequality
gives `P(|X_j−X_i| ≥ t) ≤ 2exp(−t²/(2(j−i)(1−(j−i)/M)s²))`, whence `E|X_j−X_i|⁴ ≤ C((j−i)s²)²`.
After normalisation `E|𝒳(θ)−𝒳(θ')|⁴ ≤ C'|θ−θ'|²`, and the Chentsov–Kolmogorov criterion
(Billingsley, Thm 13.5, β=4, α=1) yields tightness with a modulus depending only on s²_min.
*(proved, given Serfling 1974.)* ∎

*Remark.* The bridge is the *whole point*: the constraint `Σa_i = 2h` makes the limit a bridge, not a
BM. Lemma 2 shows the statistic we actually use is blind to that difference. (An unbalanced multiset,
`κ₀ ≠ κ₂`, adds the drift term `μθ` to `𝔅`; the contrast kills it too, since
`μ(θ − 2·θ/2) = 0`. So nothing below depends on the balance, only on `s² > 0`.)

## 3. Lemma 2 (the dyadic contrast kills the bridge — exactly)

Let `𝒟_M(θ) = (X_{θM} − 2X_{θM/2})/(s_M√M)` for θM even, and `𝒱(θ) = W(θ) − 2W(θ/2)`.

(i) **Exact covariance identity (all M).** If θM, θM/2, θ'M, θ'M/2 are all integers, then

    Cov(𝒟_M(θ), 𝒟_M(θ')) = ( 2(θ∧θ') − (θ∨θ') )_+ · M/(M−1),

i.e. the finite-M covariance matrix of the contrast is **exactly** the Brownian one, times M/(M−1).
*Proof.* Expand the four terms with Lemma 1(a). With i=θM, j=θ'M and i ≤ j:
if j/2 ≥ i the four terms cancel identically (→ 0); if i/2 ≤ j/2 ≤ i (θ ≤ θ' ≤ 2θ) the sum is
`s²/(M−1)·[i(M−j) − j(M−i) − i(M−j) + 2i(M−j/2)] = s²M²(2θ−θ')/(M−1)`. Normalising by `s²M` gives the
statement; the diagonal is the case θ=θ'. ∎

(ii) **Exact zero at scale ratio 2.** In particular `Cov(𝒟_M(θ),𝒟_M(θ')) = 0` exactly whenever
`θ' ≥ 2θ`, for every M — dyadic-separated contrasts are uncorrelated at finite M, not only in the limit.

(iii) **Limit.** Under Lemma 1, `𝒟_M ⇒ 𝒱` and, for θ < θ',
`Cov(𝒱(θ),𝒱(θ')) = (2θ − θ')_+`, `Var 𝒱(θ) = θ`, so `Cov = 2(θ∧θ') − (θ∨θ')` when `θ∨θ' < 2(θ∧θ')`
and 0 otherwise. In particular `𝒱(θ) ⊥ 𝒱(θ')` for `θ'/θ ≥ 2`.

In the notation of the task, `Z(θ) = X_θ/θ − X_{θ/2}/(θ/2) = 𝒱(θ)/θ`. On `[1−c,1]` the factor θ is
bounded below by `1−c > 0`, so `Z` and `𝒱` have identical signs and identical zero sets; the maximiser
of `Z` over `[c,1]` is the same point as `sup{θ : 𝒱(θ) > 0}`, which is `m*/M`. So the Brownian object
whose distribution we need is exactly `sup_{[1-c,1]} 𝒱`, as used below.

(iv) **The bridge correction cancels.** `𝔅(θ) − 2𝔅(θ/2) = 𝒱(θ)` identically, because the conditioning
term is `−θW(1)` and `θ − 2(θ/2) = 0`. So the statistic `𝒟_M` has, in the limit, exactly the law it
would have if the tape were a free random walk rather than a bridge — and by (i) even its finite-M
covariance is the free one, times the scalar M/(M−1). *(proved; verified exhaustively in §7I.)*

## 4. Lemma 3 (transfer from the Brownian functional to the discrete scan)

Let `𝒢_M(m) = (X_m − 2X_{⌊m/2⌋})/(s_M√M)`. Then

    (a)  b2(m)=1  ⟹  𝒢_M(m) > −1/(s_M√M);
    (b)  b2(m)=0  ⟹  𝒢_M(m) ≤ +1/(s_M√M);
    (c)  sup_m | 𝒢_M(m) − 𝒱_M(m/M) | ≤ 2/(s_M√M),  𝒱_M(θ) := (X_{θM}−2X_{θM/2})/(s_M√M).

*Proof.* (a),(b): for even m the test is exactly `𝒢_M(m) > 0`; for odd m it is
`X_m − 2X_k > X_k/k` with `k=(m−1)/2`, and `|X_k| ≤ k` (as `|b_i| ≤ 1`), so the threshold is
`|X_k/k| ≤ 1` in X-units, i.e. `≤ 1/(s_M√M)` after normalisation. (c): for odd m the contrast uses
`X_k` with `k=(m−1)/2 = m/2 − 1/2`, so `𝒢_M(m) − 𝒱_M(m/M) = b_m/(s_M√M)`; for the polygonal
interpolation the floor costs another `1/(s_M√M)`. ∎

**Consequence.** Fix c ∈ (0,1) and ε > 0. Let `U_ε := {inf_{[1-c,1]} 𝒱 > ε}` (open in C[0,1]). On
`{𝒱_M ∈ U_ε}` and for M so large that `3/(s_M√M) < ε`, Lemma 3(a)–(c) give, for every `m ≥ (1−c)M`,
`𝒢_M(m) ≥ 𝒱_M(m/M) − 1/(s_M√M) > ε − 1/(s_M√M) > 1/(s_M√M) ≥` the odd-m threshold, hence b2 = 1 at every
such m, hence `cost ≥ cM`. Therefore

    P(cost ≥ cM) ≥ P(𝒱_M ∈ U_ε)   for M ≥ M(ε),

and `𝒱_M ⇒ 𝒱` (Lemma 1 + Lemma 2(iii), θ ↦ 𝒱(θ) being a continuous linear functional of the path).
Portmanteau for the open set `U_ε`, then ε↓0:

    liminf_M P(cost ≥ cM) ≥ P(𝒱 ∈ U_ε) ↑ P(inf_{[1-c,1]} 𝒱 > 0) = P(sup_{[1-c,1]} 𝒱 < 0) = p(c).   (★)

The last equality is the sign symmetry `𝒱 =d −𝒱`: the event is "all top-c tests positive", and the
statistic we can compute with is the probability that the process stays *negative* on the same
interval. *(proved)*

## 5. Lemma 4 (p(c) > 0, quantitatively)

**Exact conditional structure.** `(𝒱(θ),𝒱(1))` is jointly Gaussian with `Cov(𝒱(θ),𝒱(1)) = (2θ−1)_+`
and `Var 𝒱(1) = 1`. Hence, given `𝒱(1) = b`,

    E[𝒱(θ) | 𝒱(1)=b] = b(2θ−1)_+ ,      Var[𝒱(θ) | 𝒱(1)=b] = θ − (2θ−1)_+² ≤ 9/16,

the maximum 9/16 being attained at θ = 5/8; and `E|𝒱(θ)−𝒱(θ')|² = 3|θ−θ'|` whenever `θ'/θ < 2`
(on `[1−c,1]` this holds for c ≤ 1/2). So conditionally on `𝒱(1) = b > 0` the process on `[1−c,1]` is a
centred Gaussian `Z` plus the **positive ramp** `b(2θ−1)`, and

    P(inf_{[1-c,1]} 𝒱 ≤ 0 | 𝒱(1)=b)
        = P( sup_{[1-c,1]} (−𝒱) ≥ 0 | b )                 (−𝒱 has mean −b(2θ−1) ≤ −b(1−2c))
        ≤ P( sup Z ≥ b(1−2c) )  ≤  exp( −((b(1−2c) − m_c)_+)²/(2·9/16) ),

since `sup(−𝒱) ≤ sup Z̃ − b(1−2c)` with `Z̃ := −𝒱 − mean` centred Gaussian. Here we use
Tsirelson–Ibragimov–Sudakov/Borell (`P(sup Z ≥ E sup Z + u) ≤ e^{−u²/(2σ²)}`, σ² = sup Var Z) and
Dudley for `m_c = E sup_{[1-c,1]} Z ≤ 1.2533·√(3c) = 2.17√c` (the entropy integral of an interval of
diameter `√(3c)` under `d(θ,θ')=√(3|θ−θ'|)`). Averaging over `b ~ |N(0,1)|`:

    p(c) ≥ 1 − E_b[ exp( −((b(1−2c) − m_c)_+)²/(9/8) ) ].                                  (★★)

This is exact enough to be checkable, and its two inputs are separately verified: the simulated
`E sup_{[1-c,1]}|𝒱−𝒱(1)|` is 0.469, 0.670, 0.955, 1.173 at c = .05, .1, .2, .3 against the Dudley
values 0.485, 0.686, 0.971, 1.189 (§7D), and `Var 𝒱(θ)=θ`, `Cov=(2θ−θ')_+` are confirmed on tapes
(§7A) and exactly at finite M (§7I).

Evaluating (★★) by quadrature:

| c | 0.02 | 0.05 | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 |
|---|---|---|---|---|---|---|---|
| m_c = 2.17√c | .307 | .485 | .686 | .841 | .971 | 1.085 | 1.189 |
| **p(c) ≥** | .246 | **.161** | .078 | .031 | .009 | .0014 | ~0 |
| c·p ≥ | .0049 | **.0080** | .0078 | .0047 | .0018 | .0004 | — |

Meanwhile the *true* p(c) (Brownian simulation, grid-stable, §7D) is much larger:

| c | 0.02 | 0.05 | 0.10 | 0.15 | 0.20 | 0.30 | 0.40 | 0.50 |
|---|---|---|---|---|---|---|---|---|
| p(c) | .4281 | .3810 | .3282 | .2867 | .2504 | .1914 | .1383 | .0913 |
| c·p(c) | .0086 | .0191 | **.0328** | .0430 | .0501 | **.0574** | .0553 | .0457 |
| 2p(c) = P(inf>0 ; 𝒱(1)>0) | .856 | .762 | .656 | .573 | .501 | .383 | .277 | .183 |

The gap between (★★) and the true p(c) comes from two conservative steps: the Borell tail is applied at
the sup-variance 9/16, whereas the conditional variance `θ − (2θ−1)²` is much smaller over most of
`[1−c,1]` (it vanishes at θ = 1); and the bound is flat (`= 1`) for `b ≤ m_c/(1−2c)`, a region carrying
most of the mass of the half-normal. Dudley's `m_c` itself is accurate to 3% here (0.485 vs 0.469).

**Two readings of the mechanism** (this is the "conditional" question): conditioning on the top test
*succeeding* (`𝒱(1)=b>0`, i.e. the second half out-means the first) tilts the whole top segment
*upward* with a ramp `b(2θ−1)`, while the conditional noise there is only `O(√c)`; so the scan persists
to `1−c` with probability `1 − e^{−b²(1−2c)²/(9/8)}` up to the Dudley correction. Since the density of
b at 0 is finite, a positive fraction of tapes persist — that is where the constant comes from. The
empirical conditional probabilities (§7G, n=4096) are `P(cost ≥ 0.1M | cost ≥ 2) = 0.681` and
`P(cost ≥ 0.2M | cost ≥ 2) = 0.515`, against Brownian `2p(c) = 0.656, 0.501`.

## 6. Theorem

For the uniform ensemble (any exchangeable ensemble that, conditionally on the counts, is a uniform
arrangement of a multiset with κ₀=κ₂ and s² ≥ 1/4 — both the "two balanced columns" ensemble
(s²≈1/2) and "uniform over Σa=2h" (s²≈2/3) qualify, and s² ≥ 1/4 holds except on a set of probability
e^{−κM}), for every c ∈ (0,1),

    P(cost₀ ≥ cM) ≥ p(c) − o(1),        E[cost₀] ≥ c·p(c)·M·(1 − o(1)).

**Rigorous constants:** c = 0.05, p ≥ 0.16 give **E[cost₀] ≥ 0.008 M**.
**Sharp (numerical) constants:** c = 0.1, p = 0.328 give E[cost₀] ≥ 0.033 M; sup over c of c·p(c) =
0.057 at c = 0.3.
**Measured:** E[cost₀] = 0.130 M ± 0.009 at M = 4096, flat in M (0.114–0.149 over n = 64…8192).

*Proof.* Combine Lemma 3(★) with Lemma 4 and `E[cost] ≥ cM·P(cost ≥ cM)`. The `o(1)` is the loss in
(★) plus the probability `s² < 1/4`; the argument does not quantify M₀. ∎

## 7. Numerical verification (all runs in this directory)

**A. Covariance of the contrast, n=1024, 400 tapes (seed 2026).** Empirical correlation of `𝒟(θ)` for
θ = .125,.25,.375,.5,.75,1 versus the predicted `(2θ∧θ' − θ∨θ')_+/√(θθ')`:

    empirical  0.397 (θ=.25,.375)   0.596 (.375,.5)   0.371 (.5,.75)   0.599 (.75,1)
    predicted  0.408                 0.577             0.408             0.577
    zero pairs (ratio ≥ 2): 0.035, −0.038, −0.037, −0.001   (predicted 0)

`sd²/θ` = 0.94, 0.99, 1.03, 0.99, 1.02, 1.05 against the exact value M/(M−1) = 1.00098 (MC se ≈ 0.07).

**C. First epoch, full region [0,n−1], 800 tapes/n, seed 2026.**

| n | 64 | 128 | 256 | 512 | 1024 | 2048 | 4096 | 8192 |
|---|---|---|---|---|---|---|---|---|
| E[cost]/M | .164 | .145 | .133 | .130 | .114 | .133 | .149 | .135 |
| P(cost=1) | .560 | .595 | .564 | .541 | .550 | .495 | .476 | .510 |
| P(cost ≥ .1M) | .370 | .328 | .325 | .325 | .289 | .321 | .361 | .325 |
| P(cost ≥ .3M) | .230 | .210 | .200 | .200 | .170 | .190 | .210 | .200 |

**D. Brownian 𝒱 on a grid (N=4096 fine grid, 40 000 paths).** `E[1−θ*] = 0.1311`
(the report's 0.1365 comes from a different grid/finite-M ensemble; both sit within 4% of the measured
0.130–0.145, and this constant plays no role in the proof); `P(𝒱(1)<0) = 0.5023`;
`E[1−θ*|𝒱(1)<0] = 0.2609`; p(c) table in §5. Grid sensitivity of p(0.1): 0.337, 0.333, 0.324, 0.321 for
N = 1024…8192.

**E. Rigorous quadrature (★★)** as tabulated in §5.

**F. Ensemble independence.** σ² = 0.500 ("two balanced columns") vs 0.666 ("uniform Σa=2h"):
E[cost]/M = 0.1345 vs 0.1267 (n=512), 0.1237 vs 0.1305 (n=2048) — equal within ±0.009. This is the
predicted insensitivity of the criterion to σ² (the test is scale-free).

**G. Discrete vs Brownian, n=4096, 3000 tapes (seed 4242).**

| c | 0.05 | 0.10 | 0.20 | 0.30 | 0.50 |
|---|---|---|---|---|---|
| P(cost ≥ cM) ± se | .3807 ± .0089 | .3293 ± .0086 | .2490 ± .0079 | .1897 ± .0072 | .0910 ± .0053 |
| p(c) Brownian | .3810 | .3282 | .2504 | .1914 | .0913 |
| ratio | 0.999 | 1.003 | 0.994 | 0.991 | 0.997 |
| E[cost/M \| cost ≥ cM] | .336 | .377 | .452 | .516 | .645 |

Also `P(cost=1)=0.5167`, `E[cost/M | cost≥2]=0.2687`, `E[cost/M]=0.1300`.

**H. Epoch geometry (150 tapes, seed 909).** `L_i = i` for every epoch (0 violations); `M_i = n−i−δ_i`;
epochs/n has minimum 0.490, 0.499 at n=512, 1024 (the deterministic bound is 0.5); `ΣM_i/n²` has
**minimum 0.250, 0.250** — the deterministic floor n²/4, attained by the all-double-cut tape.

**I. Exact covariance identity, exhaustive over the ensemble** (all words with κ₀=κ₂; M = 8, 10, 12;
1106, 8952, 73 788 words with s² > 0). At M = 8 every θ ∈ {.25,.5,.75,1} has both θM and θM/2 integral,
and the full 4×4 matrix of `E[𝒟(θ)𝒟(θ')]/(Ms²)` matches `(2θ∧θ' − θ∨θ')·M/(M−1)` to 4 decimals:
diagonal 0.2857, 0.5714, 0.8571, 1.1429 = 2/7, 4/7, 6/7, 8/7; off-diagonal 0.2857 at (.5,.75) and
0.5714 at (.75,1) = 2/7, 4/7; and 0.0000 wherever θ'/θ ≥ 2 — exactly Lemma 2(i). At M = 10, 12 the odd
θM entries show the O(1/M) floor effects the lemma excludes (0.2652 vs 0.2727 at M=12, θ=.25).

## 8. From one epoch to E[T] = Ω(n²): what is proved, what is missing

**Proved.** `E[T] ≥ E[cost₀] ≥ 0.008 n` (Theorem), and the deterministic geometry:
`L_i = i` (each epoch advances L by exactly 1), `R_i = n−1−δ_i` with `δ_i = #double cuts ≤ i`, hence
`M_i = n − i − δ_i ≥ n − 2i`; the run ends only by meeting (`C = ∅`, cost M−1) or when `L ≥ R`, which
needs `n−1 ≤ i+δ_i ≤ 2i`, i.e. `i ≥ (n−1)/2`. So **on the no-meet event** there are at least (n−1)/2
epochs and

    Σ_i M_i ≥ Σ_{i<(n−1)/2} (n−2i) = n²/4 − O(n).        (verified: min ΣM_i/n² = 0.250)

**What the "× Θ(n) epochs" step needs.** Precisely this:

> **Lemma Q (open).** For every epoch i, `E[cost_i | F_{i−1}] ≥ c M_i` with c > 0 independent of i, n.

Given Lemma Q, `E[T] = Σ_i E[cost_i] ≥ c E[Σ_i M_i]`, and by the geometry above
`Σ_i M_i ≥ n²/4 − O(n)` on the event that the run does not meet before `i = (n−1)/2`, so
`E[T] ≥ (c/4)n²(1 − o(1))`. Two remarks on that last line: (i) a meet *later* than `(n−1)/2` is
harmless (the geometry is already secured); (ii) a meet *before* it ends the run early, so the argument
also needs `P(run meets before (n−1)/2) = o(1)`. That is measured to be 0 in every run for windows of
length ≥ 64 (`q3_fresh.py`: 4.2·10⁻³ at M=16, 7·10⁻⁴ at M=32, 0 in 10⁴⁺ samples at M ≥ 64), and it is a
much weaker statement than Lemma Q — a meet requires all M−1 comparisons to come out the same way —
but it is not proved here either. So Lemma Q plus `P(early meet) = o(1)` is the exact remaining
content, and Lemma Q is *strictly stronger* than the result proved here: Theorem gives it only for
i = 0 (a fresh uniform window).

**Lemma Q is not implied by Theorem, and is not a formality:**
1. *Pathwise uniformity is false.* A descending window has `A_M ≤ A_{M/2}`, so b2 = 0 at m = M and
   `cost = 1`. Descending windows occur at every epoch — P(cost=1) = 0.48–0.60 at the first epoch and
   the double-cut (cost-1) fraction is 0.61–0.65 along trajectories (REPORT_LOWER §7) — so no statement
   of the form "every window costs ≥ cM" can hold. Lemma Q is irreducibly about a conditional
   *expectation* over the window content given the history.
2. *The conditioning is not exchangeability.* `F_{i−1}` is a system of linear inequalities on
   a[0..n−1] whose blocks reach into the current window `[L_i,R_i]` (the last test of epoch i−1
   compares means of blocks ending at `R_{i−1} ≥ R_i`). So the window word is a *polytope slice* of the
   uniform-permutation measure; uniformity of the arrangement is destroyed and Theorem does not apply.
3. *Measured sign.* It holds comfortably: comparing each epoch's realised cost with a fresh window of
   the same length, mean realised T = 16 429 vs Σ ĉ(M_L) = 13 043, ratio **1.26** (REPORT_LOWER §5,
   n=512, 250 tapes); and E[cost/M] along trajectories is 0.155–0.185 against the fresh-window 0.1365.
   The conditioning is *favourable*, but this is measurement, not proof.

**Exactly what a proof of Lemma Q would need** (either of these would do):
* a **stability estimate** for the scan cost: a coupling/metric on window words in which `cost` is
  Lipschitz, plus a bound showing the conditional law of the window word is absolutely continuous with
  respect to the uniform arrangement law with density bounded below by `e^{−o(M)}` on the set where the
  cost is large; or
* a **direct conditional argument at the epoch level**: show that the last test of epoch i−1, which is
  the only high-scale comparison touching the new window, shifts the law of the window word by an
  amount that changes the persistence probability p(c) by o(1) — i.e. a quantitative version of
  "conditioning on one comparison is a rank-one tilt". The Gaussian computation of Lemma 4 is the
  natural tool: conditioning `𝒱_i(1)` on the outcome of the previous test changes its mean by
  `O(1/√M_i)` in the worst case, which cannot move `p(c)` by more than `O(1/√M_i)` *if* the tilt is
  indeed rank-one.

Until one of these is supplied, the correct statement is: **one epoch costs Ω(M) in expectation
(proved, constant 0.008 rigorous / 0.033 sharp), and E[T] = Ω(n²) follows once the per-epoch bound is
uniform in the conditioning history (Lemma Q) — which is exactly the annealed→quenched step flagged as
the open point in REPORT_LOWER §5, now with the single-epoch half of it closed.**

> **Update (REPORT_LEMMA_Q.md).** Lemma Q is not needed. Summing the *frontier* inequality
> `cost_i ≥ cost(L_i,n−1) − δ_i` over the first `N = Θ(n)` epochs only (not over the whole
> trajectory) gives `E[T] ≥ [c₀²/(2(1+c₀)) − O(q_N)]n²` with `c₀` the constant of §6
> (`3.17·10⁻⁵ n²` rigorous, `8.2·10⁻³ n²` sharp), because `cost(i,n−1)` is an *unconditioned* suffix
> and (F2)/Theorem applies to it verbatim. The residual hypothesis is only `q_N = P(meet before epoch
> N) → 0` (measured 0: no meet at any window ≥ 64 in 174 043 epochs), which contains no conditioning.
> Lemma Q itself, in its per-history form, is not the operative statement: the conditional mean
> `E[cost_i|F_{i−1}]/M_i` falls from 0.31 to 0.001 as the preceding double-cut run grows from 0 to 50,
> so the only constant it could carry is `≈10⁻³`.

## 9. Files and reproduction

`$P = C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`, run in
`C:/Users/shini/poetry_axis/paper/probe_ologn2/`.

```
$P q12_epoch_bound.py A   # covariance of the contrast (Lemma 1, 2)
$P q12_epoch_bound.py C   # first-epoch law vs n, 800 tapes/n, seed 2026
$P q12_epoch_bound.py D   # Brownian V: p(c), E[1-theta*], grid sensitivity
$P q12_epoch_bound.py E   # rigorous (c,p) from Borell + Dudley
$P q12_epoch_bound.py F   # sigma^2-independence (two ensembles)
$P q12_epoch_bound.py G   # discrete vs Brownian at n=4096, 3000 tapes, seed 4242
$P q12_epoch_bound.py all > q12.log
```
`q12.log` holds the full output; §7 quotes it verbatim. Seeds: 2026 (A, C, F), 11 (D), 5 (D2),
4242 (G), 909 (H). No .tex file was touched.
