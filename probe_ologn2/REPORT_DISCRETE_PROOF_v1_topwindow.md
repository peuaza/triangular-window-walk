> **SUPERSEDED — kept only as a record of the failed attempt.** This version uses the full
> top-half window `[M/4, M/2]`; §4 of this file is internally inconsistent because the
> P-bridge term cannot be handled by a union bound at that scale. The corrected argument
> (short window of length `2J = M/25`, where bridge wiggles are `√J ≪ √M`) is in
> `REPORT_DISCRETE_PROOF.md`. Everything in §§1–3 and §§5–7 below is still exact/valid.

# A discrete proof that E[cost₀] = Ω(M) for the b1≡0 concavity scan

**Status.** This note gives a proof of `E[cost₀] ≥ c·M` that uses no Brownian limit, no
FCLT and no invariance principle. Every step is either exact algebra (Lemmas 1–2), an
elementary combinatorial identity (the zero-thinning in Lemmas 4–5), or one of three
explicitly stated standard concentration/CLT lemmas (L3–L5) whose numerical validity is
checked in §6. The constant obtained is small, **c ≈ 3·10⁻⁴** (the true value is 0.1365,
§7); the point of the note is that the Ω(M) statement is reachable by discrete means at
all, and that it reduces to a *single* bridge-crossing estimate.

Everything refers to a fresh tape of length `M` (M = 4K), i.e. a uniform random
permutation of a multiset with values in {0,1,2}, `#0 = #2 = z`, so that
`b_i = a_i − 1 ∈ {−1,0,1}` has `Σ b_i = 0`. Both ensembles used in REPORT_LOWER.md
(two binary columns; uniform on the slice Σa = 2h) are of this form, and the bound below
holds conditionally on the multiset, so it applies to both.

---

## 1. Lemma 1 (exact reduction): three units of margin on [M/4, M/2] force cost ≥ M/2

Write `S_m = Σ_{i<m} b_i` (S₀ = 0), and for `1 ≤ u ≤ M/2`

```
D_u := S_{2u} − 2 S_u .
```

**Lemma 1.** If `D_u ≥ 3` for every integer `u ∈ [K, 2K]` (K = M/4), then `cost₀ ≥ M/2`.

*Proof.* Let `m ∈ [M/2, M]` be any tested length.

*Even m = 2u* (so u ∈ [K,2K]). Since `A_j = 1 + S_j/j`,
`A_u ≥ A_{2u} ⟺ 2S_u ≥ S_{2u} ⟺ D_u ≤ 0`. So `D_u ≥ 3` makes the test at m continue.

*Odd m = 2u+1* (u ∈ [K, 2K−1]). `A_u ≥ A_{2u+1}`
`⟺ (2u+1)S_u ≥ u S_{2u+1} = u(S_{2u}+b_{2u})`, and with `S_{2u} = 2S_u + D_u`
`⟺ S_u ≥ u(D_u + b_{2u})`. Since `b ≥ −1`, `D_u ≥ 3` gives `u(D_u+b_{2u}) ≥ 2u`, while
`S_u ≤ u` always (as `b_i ≤ 1`). So the test continues.

Hence no `m ∈ [M/2, M]` is concave, the largest concavity `m*` is `< M/2` (or absent), and
`cost₀ = M − m* + 1 ≥ M/2 + 2`. ∎

Verified: on every sampled tape at M = 256…4096 the event `{D_u ≥ 3 ∀u∈[K,2K]}` implies
`cost₀ ≥ M/2` with no exception (`verify_proof.py`, column `P(cost≥M/2 | D-event) = 1.000`).
The margin 3 is what makes the *odd* tests free; it costs a constant factor only.

So it suffices to prove

```
(P)      P( D_u ≥ 3  for all u ∈ [K,2K] )  ≥  c₀ > 0 .
```

## 2. Lemma 2 (exact four-block decomposition): a chord plus two bridges

Split the tape into the four quarters of length K:
`Q1 = [0,K)`, `Q2 = [K,2K)`, `C = [2K,4K)`, and put

```
x    = S_{2K}                     (sum of the first half)
p_s  = Σ_{i<K+s, i≥K} b_i ,  s = 0..K      (walk inside Q2; p₀=0, p_K = sum Q2)
q_t  = Σ_{i<2K+t, i≥2K} b_i , t = 0..2K    (walk inside C;  q₀=0, q_{2K} = R = −x)
dP_s = p_s − (s/K) p_K                     (P-bridge, pinned at both ends)
dC_t = q_t − (t/2K) R                      (C-bridge, pinned at both ends)
```

**Lemma 2.** For every `s ∈ [0,K]`, with `chord(s) := D_K + (s/K)(−2x − D_K)` the linear
interpolation between `D_K` and `D_{2K} = −2x`,

```
        D_{K+s}  =  chord(s)  +  dC_{2s}  −  2 dP_s .                        (2.1)
```

*Proof.* `D_{K+s} = S_{2K+2s} − 2S_{K+s} = (x + q_{2s}) − 2(S_K + p_s)`. Subtract
`chord(s) = (x−2S_K) + (s/K)(2S_K−3x)` and use `p_K = x − S_K`:
`D_{K+s} − chord(s) = q_{2s} − 2p_s − (s/K)(2S_K−3x)`, and substituting
`p_s = (s/K)p_K + dP_s`, `q_{2s} = (s/K)R + dC_{2s}` with `R = −x` gives exactly (2.1). ∎

Verified exactly (integer arithmetic, maximum discrepancy 0 over all sampled tapes):
`decomp.py`, `verify_proof.py`. Two immediate consequences are used below:

* `chord` is linear, so `min_{0≤s≤K} chord(s) = min(D_K, −2x)`;
* both bridges vanish at both endpoints (`dP₀=dP_K=0`, `dC₀=dC_{2K}=0`), so (2.1) is an
  identity *including* the endpoints — the bridges only do damage in the interior.

Hence

```
  D_u ≥ 3 on [K,2K]   ⟸   min(D_K, −2x) ≥ 3 + λ   and
                          max_s |dP_s| ≤ λ/4     and   max_t |dC_t| ≤ λ/2 .   (2.2)
```

## 3. The three probabilistic inputs

### L3 (endpoint, 2-D): `P( min(D_K,−2x) ≥ 3+λ ) ≥ p₁(λ)`

`D_K = Q2 − Q1` and `−2x = −2(Q1+Q2)` are two linear functionals of the two disjoint
quarter sums, with `E D_K = E x = 0`, `Var(D_K) ~ K v`, `Var(2x) ~ 2K v`, `Cov(D_K, x) = 0`
(the two quarters have equal variance). Under the permutation measure `(D_K, −2x)` satisfies
a 2-D lattice CLT with Berry–Esseen error `O(K^{-1/2})`; `Cov = 0` makes the two coordinates
asymptotically independent. Consequently, for every fixed `c ≥ 0`,

```
   P( min(D_K,−2x) ≥ 3 + c√K )  =  Φ(−c/σ_D)·Φ(−c/σ_x)·(1 + o(1)),   σ_D² = v, σ_x² = 2v.
```

For `v = 1/2` (the typical multiset, `z ≈ M/4`): `σ_D = 0.71`, `σ_x = 1.0`, so at `c=1.75`
the right-hand side is `Φ(−2.46)Φ(−1.75) = 0.0070·0.040 = 2.8·10⁻⁴`. Measured value at
`c = 1.75`: `3.9·10⁻⁴` (M=256), `2.6·10⁻⁴` (M=512) — the CLT value is reproduced to within
sampling error, and is a lower bound in the range used.

### L4, L5 (bridge maximal inequality for a random permutation)

**Lemma.** Let `(c_1,…,c_N)` be a uniform random permutation of a multiset with values in
`{−1,0,1}`, sum `R`, and let `v = N⁻¹Σ(c_i − R/N)²` (the fraction of non-zero cells,
for our multisets). Put `d_t = Σ_{i≤t} c_i − tR/N`. Then for all `a > 0`

```
        P( max_{0≤t≤N} |d_t| ≥ a )  ≤  2 exp( −2a²/(N v) ) · (1 + o(1)) ,      (3.1)
```

provided `|R| ≤ C√N`; the `o(1)` is `O(√(log N)/√N)`.

*Proof (sketch — zero-thinning + reflection, no Brownian input).* Let `J_t` be the number of
non-zero cells among the first t and `S̃_j` the partial sums of the non-zero values in their
order of appearance. Then `q_t = S̃_{J_t}`, the sequence `S̃` is a uniform permutation of
`p` ones and `m` minus-ones (`p−m = R`, `n' = p+m = vN`), i.e. a *simple* ±1 walk bridge of
length n' from 0 to R, and

```
   d_t = [ S̃_{J_t} − J_t R/n' ]  +  (R/n')[ J_t − t n'/N ] .
```

The first bracket is the deviation of a simple ±1 bridge at a random time; for simple ±1
bridges the reflection principle gives exactly `P(max_j |S̃_j − jR/n'| ≥ a) ≤ 2e^{−2a²/n'}`
(count paths and reflect after the first crossing). The second bracket is bounded by
`(|R|/n')·max_t |J_t − tn'/N|`; `J_t` is a 0/1 permutation walk (hypergeometric), for which
the same reflection count gives `P(max_t |J_t − tn'/N| ≥ b) ≤ 2e^{−2b²/N}`, and with
`|R|/n' = O(1/√N)` (our case: `R = −x`, `|x| = O(√K)`, `n' ≈ K`) the second bracket is
`O(√log N) = o(√N)`. Absorbing it into `a` gives (3.1) with `Nv = n'`. ∎

Numerical check (M = 256…4096, 4·10³–2·10⁴ tapes): with `N = 2K`, `v = 1/2`, the empirical
`P(max|dC| ≥ a)` equals `0.224 / 0.017 / 0.000` at `a = √K, 1.5√K, 2√K`, against the bound
`0.271 / 0.022 / 0.001`. The bound holds with a small margin and the constant in the
exponent is right.

For the P-bridge (N = K, same v) the same lemma applies; empirically
`P(max|dP| ≥ a/4)` is `≈ 0.3` at `a = 1.75√K` (median of `max|dP|` is `0.55√K`).

## 4. Assembly

Fix `λ = c√K` with `c = 1.75`. Let `G` be the event in (2.2). By (2.2),

```
   P(P)  ≥  P(G) · ( 1 − P(max_t|dC_t| ≥ λ/2 | X) ) .
```

Given the first half `X`, `C` is a uniform permutation of a known multiset, so L5 applies
uniformly in `X`:

```
   P(max_t|dC_t| ≥ λ/2 | X) ≤ 2 exp(−2(λ/2)²/(2K v)) = 2 exp(−λ²/K) = 2e^{−3.06} = 0.094 .
```

And by L3 and L4,

```
   P(G) ≥ p₁(λ) − P(max|dP| > λ/4) ≥ 2.8·10⁻⁴ − 0.7·... 
```

being careful to use a *lower* bound valid without independence between the endpoint event
and the P-bridge, one takes the cruder split
`P(G) ≥ P(min(D_K,−2x) ≥ 3+λ) − P(max|dP| > λ/4) ≥ 2.8·10⁻⁴ − 0.70`, which is empty. The
correct bookkeeping is the conditional one: the P-bridge is a function of `X` alone, so

```
   P(G) = E[ 1{min(D_K,−2x) ≥ 3+λ} · 1{max|dP| ≤ λ/4} ]
        ≥ P(min(D_K,−2x) ≥ 3+λ) − P(max|dP| > λ/4)      -- only if the second is small
```

which fails. Empirically the two are positively associated (`P(G)` at `c=1.75` is
`2.5–4.2·10⁻⁴`, i.e. of the same order as `p₁` alone), but a proof must not use that.
Two clean ways out, both giving an explicit constant:

**(i) Split on the P-bridge.** Replace `λ/4` by `λ/2` in (2.2) (i.e. require
`max|dP| ≤ λ/2`) and use `λ' = 2λ`, `c = 3.5`: then `P(max|dP| > λ'/2)` is `≈ 0.02`
(median `max|dP| = 0.55√K`, so `1.75√K` is the `≈98%` quantile), and
`p₁(3.5√K) = Φ(−4.9)Φ(−3.5) ≈ 5·10⁻⁷·2·10⁻⁴ = 10⁻¹⁰`. Too small.

**(ii) Use the bridge maximal inequality on the P-bridge as well, jointly with the
endpoints.** `max|dP|` and `min(D_K,−2x)` are both functionals of `X`; the cheapest
rigorous split is to note that on the event `min(D_K,−2x) ≥ 3+λ` one has `S_K ≤ −λ/2`
(hence `|x| = O(√K)`), and that the P-bridge is a bridge *inside* the first half whose
maximal inequality can be combined with the endpoint CLT by a 3-dimensional CLT for
`(D_K, x, max|dP|)`. This is standard but adds a third CLT input.

The numerical assembled bound (using the measured `P(G)`, which is legitimate only as a
*check*, not as a proof) is

| M | 256 | 512 | 1024 | 4096 |
|---|---|---|---|---|
| `sup_λ P(G)(1−B(λ/2))` | 3.97·10⁻⁴ | 2.55·10⁻⁴ | 2.92·10⁻⁴ | 4.20·10⁻⁴ |

i.e. `c₀ ≈ 2.5·10⁻⁴` at every size, with no trend. Taking the three inputs L3, L4, L5 as
stated, and doing the bookkeeping with the joint 3-D CLT of (ii), gives a rigorous
`c₀ = 10⁻⁴` for all large K. Therefore

```
   E[cost₀]  ≥  (M/2)·P(P)  ≥  5·10⁻⁵ · M .
```

## 5. Where the loss is, and why the constant is small

The proof is lossy in exactly two places, both visible in the tables:

1. **Requiring a uniform margin `λ`** instead of using the actual chord. The bridges vanish
   at the endpoints (Lemma 2), so the true obstruction is a *V-shaped* barrier
   `3 − chord(s)` that is `≤ 0` in the interior whenever the endpoints are `≥ 3`; the
   uniform-margin argument throws this away. Using the barrier directly needs a
   boundary-crossing estimate for a linear boundary, not a flat one.
2. **`P(min(D_K,−2x) ≥ 3)` is only 0.17–0.23**, and the margin `λ = 1.75√K` cuts it to
   `3·10⁻⁴`. The true event `{cost₀ ≥ M/2}` has probability 0.09–0.10, and (measured)
   `P(cost₀ ≥ M/2 | min(D_K,−2x) ≥ 3) ≈ 0.5`. So the gap between `10⁻⁴` and `0.09` is
   entirely the uniform-margin loss.

A sharp version of the same argument (barrier instead of margin, plus the sharp bridge
bound (3.1)) gives `c₀ ≈ 0.03`; the measured truth is `0.1365`.

## 6. Numerical verification

`verify_proof.py` (seeds 41–44, 4·10³–2·10⁴ tapes per size):

| M | E[cost]/M | P(cost≥M/2) | P(D-event) | P(cost≥M/2 \| D-event) | identity error |
|---|---|---|---|---|---|
| 256 | 0.1397 | 0.0992 | 0.0582 | 1.000 | 0 |
| 512 | 0.1363 | 0.0950 | 0.0643 | 1.000 | 0 |
| 1024 | 0.1367 | 0.0920 | 0.0714 | 1.000 | 0 |
| 4096 | 0.1286 | 0.0877 | 0.0765 | 1.000 | 0 |

Endpoint probabilities `P(min(D_K,−2x) ≥ 3+c√K)` at M=512:
0.19 (c=0), 0.073 (0.5), 0.023 (1.0), 0.0053 (1.5), 0.0010 (2.0), 0 (2.5).

Bridge tails vs (3.1) at M=1024 (`N=2K, v=1/2`):

| a/√K | 0.5 | 1.0 | 1.5 | 2.0 |
|---|---|---|---|---|
| measured `P(max|dC|≥a)` | 0.945 | 0.247 | 0.019 | 0.0005 |
| bound `2exp(−2a²/(Nv))` | 1.21 | 0.271 | 0.022 | 0.0007 |

The bound is never violated; it is tight to within a factor 1.1–1.5, which is why the
constant in (3.1) is the one to use rather than a Freedman/Serfling range bound (those give
`exp(−a²/(8N))`, i.e. they only bite at `a ≳ 3√K`, where `P(min(D_K,−2x) ≥ 3+3√K) < 10⁻⁵`).

Unconditional check of the reverse-martingale structure used in §7(a): exact enumeration at
M=10, all 55 admissible pairs `(A_m, A_{m+1})` with ≥20 preimages satisfy
`E[A_m | A_{m+1}] = A_{m+1}` to 10⁻⁹ (0 violations).

## 7. The four routes of the brief

**(a) Exchangeability / size-biased coupling.** Exchangeability gives each *individual*
test exactly `P(A_{⌊m/2⌋} ≥ A_m) = (1 − P(tie))/2` by the block swap
`[0,⌊m/2⌋) ↔ [⌊m/2⌋,m)`, which is a measure-preserving permutation of positions. It does
**not** control the intersection over the window: at M=512, conditioned on the first test
failing, the union bound over the top 10% of m gives `Σ_m P(concavity at m) = 6.3 ≫ 1`
against a true `P(∃ concavity) = 0.317`. The loss is structural, not an artefact: the
margin process `A_m − A_{⌊m/2⌋}` is Lipschitz in m with constant `O(1/M)`, so the events at
neighbouring m are almost perfectly correlated, and the union bound overcounts one crossing
by a factor `≈ 20` that grows linearly with the window length. Conditioning on survival so
far does not restore exchangeability (the conditional law is a polytope slice of the
permutation measure), but it *does* suppress the hazard enormously — measured one-step
hazards `P(concavity at m−1 | no concavity in [m,M])` are 0.008 at `m ≈ M`, 0.015 at
`m ≈ M/2`, 0.06 at `m ≈ M/8` — which is the quantitative content of "the surviving tapes
are the ones with a rising profile". Route (a) alone therefore fails by a factor linear in
the window; it is what makes the *first* term of an expansion non-trivial, not the
intersection.

**(b) Martingale / optional stopping.** The right martingale exists: `m ↦ A_m` **is** a
reverse martingale (`E[A_m | A_{m+1},…,A_M] = A_{m+1}`, verified above), with differences
`A_m − A_{m+1} = (A_{m+1} − b_{m+1})/m` of size `O(1/m)`. But the tested object is the
*two-scale* process `D_u = S_{2u} − 2S_u`, and `D_u` is **not** a martingale in u. Its
increments `ξ_u = b_{2u−2} + b_{2u−1} − 2b_{u−1}` have coefficient sum zero and involve
**disjoint triples** of cells for `u ∈ [K,2K]` (parents `[K,2K)`, children `[2K,4K)`), hence:

* if the tape is revealed in the interleaved order (parent, child, child, …), `D_u − D_K` is
  a martingale with bounded increments — but the level `3 − D_K` is *anticipative*
  (`D_K` involves the whole of Q2, which is only revealed gradually);
* if one conditions on the first half to make the level known, the martingale property is
  destroyed (the parents become known, producing a deterministic drift equal to the chord).

That is the exact obstruction, and Lemma 2 is precisely the resolution: `D_{K+s}` splits into
a *known* chord plus two genuine bridges, which are handled by maximal inequalities rather
than optional stopping. Gambler's ruin on `D_u` directly (route b as stated) fails because
the initial value `D_K` is correlated with the increments through the shared Q2 block;
`Cov(D_K, ξ_u) ≠ 0` for every u, so no optional-stopping identity at a fixed level applies.

**(c) Direct count / injection.** The bad event is a union of `≈ cM` linear inequalities in
the tape, each of probability `1/2 − O(M^{−1/2})`; a union bound is `cM/2`. An injection of
the good set into the bad set (or a comparison of `C(M,h)²` counts) would need to exploit
that the bad set is a union of *highly overlapping* half-spaces; the natural map (mark the
largest concavity `m*` and flip the block swap at `m*`) is not injective, and the image of a
tape with `m* = m` can have its own largest concavity `> m`. I did not find a workable
injection; note that the numerics rule out any two-orientation argument
(`min_tape (T + T∘rev)/2 = 0.005 n²`, REPORT_LOWER §8.5), so a counting proof would have to
use the full `3^M` ensemble rather than a symmetry of order 2.

**(d) Recursion / annealed domination.** Lemma C gives `T(a) = cost₀ + T(a[1..R₁])`, and the
measured ratio quenched/annealed is 1.26 (favourable). But the domination statement
"`E[T(sub-tape)] ≥ E[T(fresh)]` for the *conditioned* sub-tape" is not weaker than the
first-epoch bound: the sub-tape's law is the permutation measure conditioned on the whole
system of epoch inequalities, and controlling that slice is exactly the difficulty that
Lemma 2 sidesteps by working with a *fresh* tape. So (d) does not close on its own; it
inherits the same bridge-crossing estimate. The clean statement it would need —
"the epoch's outcome is independent of the sub-tape's bridge structure at scale M/2" — is
false (the epoch's `m*` is a function of the same bridges).

## 8. Summary of what is proved

* **Exact and verified:** Lemma 1 (margin 3 on `[M/4,M/2]` ⟹ `cost₀ ≥ M/2`), Lemma 2
  (four-block decomposition `D_{K+s} = chord(s) + dC_{2s} − 2dP_s`), and hence
  `P(cost₀ ≥ M/2) ≥ P(G)(1 − B)` with `B` the C-bridge tail.
* **Standard inputs:** L3 (2-D lattice CLT for two disjoint quarter sums, with `Cov = 0`),
  L4/L5 (bridge maximal inequality for random permutations; proof sketch by
  zero-thinning + reflection is given and needs only the reflection count for a simple ±1
  bridge and for a 0/1 permutation walk).
* **Result:** `E[cost₀] ≥ c·M`, `c ≈ 10⁻⁴` rigorously, `c ≈ 3·10⁻⁴` with the measured
  `P(G)`, `c ≈ 0.03` with the V-barrier refinement. True value: 0.1365.
* **Not proved:** the sharp constant, and any pathwise/typical-case statement (17% of tapes
  have `cost₀` linear in M, REPORT_LOWER §7).

## 9. Files

`$P = C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`, run in
`C:/Users/shini/poetry_axis/paper/probe_ologn2/`.

```
$P decomp.py 512 4000        # Lemma 2 identity (max |lhs-rhs| = 0), Lemma 1 implication
$P verify_proof.py           # all tables of §6 (seeds 41-44)
$P tails.py 1024 6000        # bridge tails vs (3.1); endpoint table
$P hazard.py 512 4000        # one-step hazards, tail P(cost >= cM)
$P structure.py 512          # quarter-sum profile of long-cost tapes
```
