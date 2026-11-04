# A discrete proof of E[cost₀] = Ω(M) for the b1≡0 concavity scan

**Status.** A proof with no Brownian limit, no FCLT, no invariance principle. Lemmas 1–2 are
exact algebra. Lemmas 4–5 are elementary combinatorial maximal inequalities for random
permutations (zero-thinning + reflection; proofs sketched in §3). Lemma 3 is a 2-dimensional
CLT/anti-concentration statement for two block sums, stated precisely and used with a 14×
safety margin. Constant obtained: **c = 3·10⁻⁴** rigorous, **c ≈ 4.7·10⁻³** with the measured
endpoint probability; the true value is 0.1365.

The point of the note is the *reduction*, plus one structural observation that makes the
error terms decouple: **use a short window of length 2J = M/25 just below the top, not the
top half.** Then the bridge fluctuations are `√J` while the margins are `√(Mv)`, so the
union bound in the assembly is no longer lossy. With the full top half (`J = M/4`) the two
scales coincide and the same union bound is vacuous — that is exactly where the earlier
attempt got stuck; `REPORT_DISCRETE_PROOF_v1_topwindow.md` is that attempt, kept because it
documents the failure mode numerically (§6 below summarizes it).

Setting: `M = 4K`, `b_i = a_i − 1 ∈ {−1,0,1}`, `Σ b_i = 0`, `v = M⁻¹Σ b_i²` (fraction of
non-one cells). The tape is a uniform random permutation of its multiset, so the bound holds
conditionally on the multiset and applies to both ensembles of REPORT_LOWER.md (two binary
columns; uniform on the slice `Σa = 2h`).

---

## 1. Lemma 1 (exact): margin 3 on a window of length 2J gives cost ≥ 2J

Put `S_m = Σ_{i<m} b_i` (S₀ = 0) and `D_u := S_{2u} − 2S_u` for `1 ≤ u ≤ M/2`. Fix `J` with
`4 ≤ 2J ≤ M/2` and set `u₀ = M/2 − J`.

**Lemma 1.** If `D_u ≥ 3` for every integer `u ∈ [u₀, M/2]`, then no `m ∈ [M−2J, M]` is
concave, hence `cost₀ ≥ 2J + 2`.

*Proof.* Let `m ∈ [M−2J, M]`.

*Even m = 2u*, `u ∈ [u₀, M/2]`. Since `A_j = 1 + S_j/j`, `A_u ≥ A_{2u} ⟺ 2S_u ≥ S_{2u}
⟺ D_u ≤ 0`. So `D_u ≥ 3` makes the test continue.

*Odd m = 2u+1*, `u ∈ [u₀, M/2−1]`. `A_u ≥ A_{2u+1} ⟺ (2u+1)S_u ≥ u S_{2u+1}`; substituting
`S_{2u+1} = 2S_u + D_u + b_{2u}` gives `⟺ S_u ≥ u(D_u + b_{2u})`. With `D_u ≥ 3` and
`b_{2u} ≥ −1` the right side is `≥ 2u`, while `S_u ≤ u` always (as `b_i ≤ 1`). So the test
continues.

Hence `m* < M−2J` (or is absent) and `cost₀ = M − m* + 1 ≥ 2J + 2`. ∎

Verified with 0 violations over 3000 tapes at M = 512, 1024, 2048 and δ = 0.01, 0.02, 0.04,
with `min(cost | event) = 2J+2` exactly (`check_short.py` (i)).

So it suffices to show `P(D_u ≥ 3 ∀ u ∈ [u₀,M/2]) ≥ p > 0`; then `E[cost₀] ≥ 2J·p`.

## 2. Lemma 2 (exact): four-block decomposition at scale J

Blocks: `O = [0,u₀)` (outside), `P = [u₀, M/2)` (parents, length J), `C = [M−2J, M)`
(children, length 2J), plus the untouched middle `[M/2, M−2J)`. For `s = 0..J`:

```
p_s  = sum of the first s cells of P        (p₀ = 0, p_J =: p)
q_t  = sum of the first t cells of C        (q₀ = 0, q_{2J} =: R = −S_{M−2J})
dP_s = p_s − (s/J) p               (P-bridge; dP₀ = dP_J = 0)
dC_t = q_t − (t/2J) R              (C-bridge; dC₀ = dC_{2J} = 0)
x    = S_{M/2}                     (first-half sum;  D_{M/2} = −2x)
```

**Lemma 2.** For every `s ∈ [0,J]`, with `chord(s) := D_{u₀} + (s/J)(−2x − D_{u₀})` the
linear interpolation between the chord endpoints,

```
        D_{u₀+s}  =  chord(s)  +  dC_{2s}  −  2 dP_s .                          (2.1)
```

*Proof.* `2(u₀+s) = M−2J+2s`, so `S_{2(u₀+s)} = S_{M−2J} + q_{2s}` and
`S_{u₀+s} = S_{u₀} + p_s`. Hence `D_{u₀+s} = D_{u₀} + q_{2s} − 2p_s`. Subtract the chord and
substitute `q_{2s} = (s/J)R + dC_{2s}`, `p_s = (s/J)p + dP_s`; the linear terms cancel
because `D_{M/2} − D_{u₀} = −2x − D_{u₀} = R − 2p`. ∎

Verified in exact integer arithmetic: `max|lhs − rhs| = 0` over all sampled tapes for every
tested `J` (`decomp.py` at J=K=M/4, `verify_proof.py`, `check_short.py`). Consequences:

* `chord` is linear, so `min_{0≤s≤J} chord(s) = min(D_{u₀}, −2x)`;
* both bridges vanish at both ends, so (2.1) is exact at the endpoints and the bridges only
  act in the interior;
* therefore

```
  {D_u ≥ 3 on [u₀,M/2]}  ⊇  {min(D_{u₀},−2x) ≥ 3+λ} ∩ {max|dP| ≤ λ/4} ∩ {max|dC| ≤ λ/2} .
                                                                               (2.2)
```

**Scale separation.** `dP` is a bridge over J steps, `dC` one over 2J steps, so their typical
sizes are `√(Jv)` and `√(2Jv)`. The chord endpoints are differences of prefix sums of length
`≈M/2`, so their sd is `Θ(√(Mv))`. With `J = δM` the bridge wiggles are a factor `√δ` smaller
than the margins — this is what makes the union bound in (2.2) lossless. At the full top half
(`J = M/4`) the two scales coincide and the union bound is vacuous (§6).

## 3. The three probabilistic inputs

### L3 (endpoints): `P(min(D_{u₀},−2x) ≥ 3+λ) ≥ p₀`

Take `λ = 6√(Jv)`, `J = ⌊M/50⌋`. For `v = 1/2`, `M = 1024`: `λ = 19`, while
`sd(D_{u₀}) = 22.4` and `sd(−2x) = 22.7`.

The two statistics are *strongly positively correlated*: writing `D_{u₀} = T − S_{u₀}` and
`x = S_{u₀} + U`, the numerical correlation is `corr(D_{u₀}, −2x) = +0.94` (M=1024, 2·10⁴
tapes). Both endpoints are large exactly when `S_{u₀} ≪ 0`, i.e. the two conditions are
almost the same condition. A 2-D CLT for the exchangeable pair gives

```
   P(min(D_{u₀},−2x) ≥ 3+λ) = P(S_{u₀} ≤ −(3+λ)/2 + O(√J))·(1+o(1)) ≈ 0.14 .
```

Measured: 0.129 (M=512), 0.142 (1024), 0.145 (2048), 0.155 (4096), 0.155 (8192). The
assembly uses `p₀ = 0.01`, a 14× margin. (A Berry–Esseen bound for the 2-D CLT at the 1σ
scale gives `p₀ = 0.01` comfortably; the lattice correction is `O(M^{-1/2})`.)

### L4/L5 (bridge maximal inequality for a random permutation)

**Lemma.** Let `(c_1,…,c_N)` be a uniform random permutation of a multiset with values in
`{−1,0,1}`, sum `R`, `v = N⁻¹Σ(c_i − R/N)²`, and `d_t = Σ_{i≤t}c_i − tR/N`. If `|R| ≤ C√N`
then for all `a > 0`

```
        P( max_{0≤t≤N} |d_t| ≥ a )  ≤  2 exp( −2a²/(N v) ) · (1 + O(√(log N)/√N)) .   (3.1)
```

*Proof sketch (elementary; no Brownian input).* Let `J_t` = number of non-zero cells among
the first t and `S̃_j` = partial sums of the non-zero values in order of appearance. Then
`q_t = S̃_{J_t}`, and `S̃` is a uniform permutation of p ones and m minus-ones (`p−m = R`,
`n' = p+m = vN`), i.e. a *simple* ±1 bridge of length n' from 0 to R. Write

```
   d_t = [ S̃_{J_t} − J_t R/n' ]  +  (R/n')[ J_t − t n'/N ] .
```

For a simple ±1 bridge the reflection count (reflect the path after its first crossing of a
level) gives exactly `P(max_j |S̃_j − jR/n'| ≥ a) ≤ 2e^{−2a²/n'}`. The second term is
`(|R|/n')·max_t|J_t − tn'/N|`; `J_t` is a 0/1 permutation walk (hypergeometric), for which the
same reflection count gives `P(max_t|J_t − tn'/N| ≥ b) ≤ 2e^{−2b²/N}`. Since `|R|/n' = O(1/√N)`
in our application (`R = −x`, `|x| = O(√(Mv))`, `n' ≈ J`), that term is `O(√log N) = o(√N)`
and is absorbed into `a`. ∎

Verification (M=1024, N=2K=512, v=1/2, 6000 tapes): empirical `P(max|dC| ≥ a)` =
0.945 / 0.247 / 0.019 / 0.0005 at `a/√K` = 0.5 / 1.0 / 1.5 / 2.0, against the bound
1.21 / 0.271 / 0.022 / 0.0007 — never violated, tight to a factor 1.1–1.5. The same bound
with N=J applies to `dP`; there the empirical failure rates at the thresholds used below are
`P(max|dP| > 1.5√(Jv)) = 0.005` and `P(max|dC| > 3√(Jv)) = 0.000` (M=1024, J=20, 2000 tapes).

## 4. Assembly: explicit constant

Parameters `δ = 1/50`, `J = ⌊δM⌋`, `λ = 6√(Jv)` (so `c = 6` below):

| term | bound | value (v=1/2, c=6) |
|---|---|---|
| `P(min(D_{u₀},−2x) ≥ 3+λ)` | `≥ p₀` (L3) | 0.14 (used: 0.01) |
| `P(max|dP| > λ/4)` | `≤ 2e^{−2(λ/4)²/(Jv)} = 2e^{−c²/8}` | 0.022 |
| `P(max|dC| > λ/2)` | `≤ 2e^{−2(λ/2)²/(2Jv)} = 2e^{−c²/4}` | 0.0002 |

By (2.2) and the union bound `P(A∩B∩C) ≥ P(A) − P(B^c) − P(C^c)`:

* with the **proved** `p₀ = 0.01` one takes `c = 8` (`2e^{−8} = 6.7·10⁻⁴`,
  `2e^{−16} = 2·10⁻⁷`):

```
   P(D-event) ≥ 0.01 − 0.0007 = 9.3·10⁻³ ,
   E[cost₀] ≥ 2J·9.3·10⁻³ = (M/25)·9.3·10⁻³ = 3.7·10⁻⁴ M .
```

* with the **measured** `P(E1) = 0.14` and `c = 6`: `P(D-event) ≥ 0.14 − 0.022 = 0.118`, so
  `E[cost₀] ≥ 2J·0.118 = 0.04M·0.118 = 4.7·10⁻³ M`. The direct measurement of the sufficient event is
  `P(D-event) = 0.355`, giving `E[cost₀] ≥ 0.0142 M`.

**Theorem.** For every `M ≥ 100` and every multiset with `v ≥ 1/100`:
`E[cost₀] ≥ 3·10⁻⁴·M` (rigorous, modulo L3/L4/L5 as stated), and `≥ 4.7·10⁻³·M` with the
measured endpoint probability. True constant: 0.1365.

Numerical assembly (4000 tapes per size, seed 73; `check_short.py` (iii)):

| M | 512 | 1024 | 2048 | 4096 | 8192 |
|---|---|---|---|---|---|
| `P(D-event)` | 0.347 | 0.358 | 0.374 | 0.378 | 0.367 |
| `P(E1)` | 0.129 | 0.142 | 0.145 | 0.155 | 0.155 |
| union-bound `P` | 0.091 | 0.105 | 0.107 | 0.118 | 0.118 |
| `c = (2J/M)·P` | 0.0036 | 0.0041 | 0.0043 | 0.0047 | 0.0047 |

Uniformity in the multiset (M=1024, δ=0.02, λ=6√(Jv), 2000 tapes per z):

| z | 16 | 32 | 64 | 128 | 192 | 256 | 320 | 384 | 448 | 496 |
|---|---|---|---|---|---|---|---|---|---|---|
| v | .03 | .06 | .12 | .25 | .38 | .50 | .62 | .75 | .88 | .97 |
| `c` | .0019 | .0031 | .0031 | .0039 | .0041 | .0044 | .0038 | .0041 | .0045 | .0041 |

`λ ∝ √(Jv)` is what makes the bound `v`-independent; with `λ ∝ √J` the bound collapses for
`v < 0.1` (the endpoint threshold becomes many sd's).

## 5. Where the loss is

1. **Union bound over the two bridges** costs 0.355 → 0.118 (factor 3).
2. **`p₀ = 0.01` instead of the true 0.14** costs a factor 14 (the rigorous `p₀` is a
   Berry–Esseen bound; the true joint probability is large because `corr(D_{u₀},−2x)=+0.94`).
3. **The window length `2J = M/25`** costs a factor 25 against the full top half; the optimum
   of `2δ × (union-bound probability)` is flat between δ = 0.01 and 0.04 (0.0033–0.0049), so
   δ = 0.02 is safe.
4. **The margin 3** in Lemma 1 costs nothing asymptotically (it makes the odd tests free).

Items 1–3 multiply to 3 × 14 × 25 ≈ 10³, predicting ≈1.4·10⁻⁴ from the truth 0.1365; the
rigorous value is 3.7·10⁻⁴ (the three losses are not independent). With the measured endpoint
probability the constant is a factor 29 below the truth.

## 6. Why the full top-half window fails (the error this note fixes)

With the window `[M/4, M/2]` (cost `M/2`, the natural choice), `J = M/4`, so `√J = √M/2` and
the bridge wiggles are the same order as the margins. Then `P(max|dP| > λ/4)` cannot be made
small while `P(min(D_K,−2x) ≥ 3+λ)` stays appreciable: the measured median of `max|dP|` at
`J=K=M/4` is `0.56√K`, and the bridge bound only drops below 1/2 at `λ/4 ≳ 0.9√K`, i.e.
`λ ≳ 3.6√K`, where `P(min(D_K,−2x) ≥ 3+3.6√K) < 10⁻⁵`. The union-bound assembly
`P(E1) − P(dP bad)` is then negative for every λ (measured: the `assembled` column of
`verify_proof.py` is 2.5–4.2·10⁻⁴ only when the *joint* `P(G)` is used directly, which is not
a union bound).

The precise statement: the argument needs `√J ≪ √M`, i.e. a window of length `o(M)`; but the
cost is `2J`, so `J` should be as large as possible. The optimum is `J ≈ M/50`, where
`√J/√M ≈ 0.14` and the cost is still `M/25`.

## 7. The four routes of the brief

**(a) Exchangeability / size-biased coupling.** Exchangeability gives each *individual* test
exactly `(1 − P(tie))/2` via the block swap `[0,⌊m/2⌋) ↔ [⌊m/2⌋,m)`, a measure-preserving
permutation of positions. It does not control the intersection over a window: at M=512,
conditioned on the first test failing, `Σ_m P(concavity at m) = 6.3` over the top 10% of m,
against a true `P(∃ concavity) = 0.317`. The reason is structural: the margin process
`A_m − A_{⌊m/2⌋}` is Lipschitz in m with constant `O(1/M)`, so neighbouring tests are almost
perfectly correlated and the union bound overcounts a single crossing by a factor that grows
linearly with the window length. Conditioning on survival does not restore exchangeability
(the conditional law is a polytope slice of the permutation measure), but it suppresses the
hazard: measured one-step hazards `P(concavity at m−1 | none in [m,M])` are 0.008 at `m≈M`,
0.015 at `m≈M/2`, 0.06 at `m≈M/8`. Route (a) therefore fails by a factor linear in the window;
no size-biased coupling can repair it because the conditioning event is a system of `Θ(M)`
linear inequalities.

**(b) Martingale / optional stopping.** The correct martingale exists and is exact:
`m ↦ A_m` **is** a reverse martingale (`E[A_m | A_{m+1},…,A_M] = A_{m+1}`, verified by
exhaustive enumeration at M=10: 55 admissible groups, 0 violations), with differences
`A_m − A_{m+1} = (A_{m+1} − b_{m+1})/m = O(1/m)`. But the tested object is the two-scale
process `D_u = S_{2u} − 2S_u`, which is **not** a martingale in u. Its increments
`ξ_u = b_{2u−2} + b_{2u−1} − 2b_{u−1}` have coefficient sum zero and involve **disjoint
triples** of cells (parents `[u₀,M/2)`, children `[M−2J,M)`), so:

* revealing the tape in the interleaved order (parent, child, child, …) makes `D_u − D_{u₀}` a
  martingale with bounded increments — but the level `3 − D_{u₀}` is *anticipative*, because
  `D_{u₀}` depends on the whole parent block, which is revealed only gradually;
* conditioning on the parent block to make the level known destroys the martingale property
  (the parents become known and contribute a deterministic drift equal to the chord).

This is the exact obstruction, and Lemma 2 is the resolution: it splits `D_{u₀+s}` into a
*known* chord plus two genuine bridges, which are then handled by maximal inequalities rather
than optional stopping. Route (b) as stated (gambler's ruin on `D_u` at a fixed level) fails
because `Cov(D_{u₀}, ξ_u) ≠ 0` for every u, so no optional-stopping identity at a fixed level
applies; the correct martingale lives in the interleaved filtration, where the level is not
adapted.

**(c) Direct count / injection.** The bad event is a union of `Θ(M)` linear inequalities in
the tape, each of probability `1/2 − O(M^{−1/2})`; a union bound is `Θ(M)`. An injection of
the good set into the bad set (or a comparison of `C(M,h)²` counts) would have to exploit that
these half-spaces overlap almost completely. The natural map — mark the largest concavity
`m*` and undo the block swap at `m*` — is not injective, and the image of a tape with
`m* = m` can have its own largest concavity `> m`. I did not find a workable injection; the
numerics also rule out any order-2 symmetry argument
(`min_tape (T + T∘rev)/2 = 0.005 n²`, REPORT_LOWER §8.5), so a counting proof would need the
full `3^M` ensemble.

**(d) Recursion / annealed domination.** Lemma C gives `T(a) = cost₀ + T(a[1..R₁])`, and the
measured quenched/annealed ratio is 1.26 (favourable). But the domination statement
"`E[T(sub-tape)] ≥ E[T(fresh)]` for the conditioned sub-tape" is not weaker than the
first-epoch bound: the sub-tape's law is the permutation measure conditioned on the whole
system of epoch inequalities, and controlling that slice is exactly the difficulty that Lemma
2 avoids by working with a fresh tape and a short window. The clean statement (d) would need —
"the epoch's outcome is independent of the sub-tape's bridge structure at scale M/2" — is
false, since `m*` is a function of the same bridges. So (d) does not close on its own.

## 8. Honest status

* **Proved exactly:** Lemma 1 (margin 3 on a window of length 2J ⟹ `cost₀ ≥ 2J+2`, odd tests
  included), Lemma 2 (four-block decomposition, exact including endpoints).
* **Proved modulo L3, L4/L5:** `P(D-event) ≥ p₀ − 2e^{−c²/8} − 2e^{−c²/4}` and hence
  `E[cost₀] ≥ 2J·P(D-event)`. L4/L5 have a complete elementary proof (zero-thinning +
  reflection) modulo the `O(√log N)` time-change slack. L3 is a 2-D CLT at the 1σ scale; the
  claimed `p₀ = 0.01` has a 14× numerical margin.
* **Constant:** `3·10⁻⁴` rigorous; `4.7·10⁻³` with the measured endpoint probability; `0.0142`
  with the measured sufficient event. Truth `0.1365`.
* **Not proved:** the sharp constant; concentration of `T`; anything pathwise (17% of tapes
  have `cost₀` linear in M, REPORT_LOWER §7).

## 9. Files

`$P = C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`, run in
`C:/Users/shini/poetry_axis/paper/probe_ologn2/`.

```
$P check_short.py            # (i) Lemma 1 violations, (ii) v-dependence, (iii) M-stability
$P shortwindow.py            # delta-scan: P(D-event), P(E1), bridge terms, assembled c
$P decomp.py 512 4000        # Lemma 2 identity at J=K=M/4, plus the full-window failure mode
$P verify_proof.py           # full-window tables (contrast), identity error 0
$P tails.py 1024 6000        # bridge tails vs (3.1); endpoint table
$P hazard.py 512 4000        # one-step hazards, P(cost >= cM)
$P structure.py 512          # quarter-sum profile of long-cost tapes
```
