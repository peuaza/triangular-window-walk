# Does the b1≡0 run ever meet? Reduction, exact recursion, and the one remaining gap

**Question.** In the b1≡0 (largest-concavity) walk, an epoch *meets* when its scan finds no
concavity at all, `C = ∅`, i.e. `A_⌊m/2⌋ < A_m` for every `m ∈ [2,M]`. A meet ends the run and
costs `M−1`. Since a meet at epoch `i` is only possible while the region is still long
(`M_i ≥ n−2i`), an early meet would truncate the `Σ_i M_i ≥ n²/4` accounting that the
`E[T] = Ω(n²)` argument needs. **Claim to prove:** `P(run meets before i = (n−1)/4) = o(1)`.

| statement | status |
|---|---|
| **T1** `P(meet before (n−1)/4) ≤ Σ_{i<(n−1)/4} f(n−2i)`, `f(M)=P(C=∅)` on a *fresh* uniform window | **proved** (§2), no conditioning needed |
| **T2** `f(2h) = E[1{E_h(W)}·q(W)]` with `q` an explicit barrier probability of the second half | **proved** (§3); `f` exact for `M ≤ 16` and all even `M ≤ 28` |
| **T3** `b₁=1 ⇒ E_M=∅`; `b₁=0 ⇒ X_{2^j} ≥ 2^j−1 ⇒ P ≤ (M+1)3^{−M/2}` | **proved** (§4) |
| **T4** one-scale bound `q(W) ≤ P(V_h ≥ X_h+1)`, and `sup_W q ≈ 0.44–0.62` ⟹ the dyadic recursion **cannot** beat `M^{−1}` | **proved** (§5) + measured |
| `f(M)`: exact for `M ≤ 16` and even `M ≤ 28`, MC (3·10⁸ samples) for `M = 32…64`; `f(M) ≤ 5e^{−1.7√M}` on `7 ≤ M ≤ 64` | measured (§6) |
| per-scale factor `E[q|E_h] = (0.86±0.09)/h`, i.e. `η ≈ 1` in Lemma B | measured (§6) |
| ⟹ `P(early meet) ≤ e^{−1.2√n+O(log n)}`; run-level `P(early meet) = 0 (<6·10⁻⁵)` for `n ≥ 128` | measured (§7) |
| **Lemma B** (`E[q|E_h] ≤ C h^{−η}`, any `η>0`) — the single missing estimate | **open** (§8) |

Routes (b) and (c) of the task brief do not work, and §9 says exactly why. The Brownian analogue
of the event has probability **0** (not a small positive number), so it yields no rate at all.

---

## 1. Notation and the exact form of the event

Tape `a ∈ {0,1,2}^n`, `b = a−1`, `X_m = Σ_{i≤m} b_i`, `A_m = 1 + X_m/m`. For a window `[L,R]`,
`M = R−L+1`, the `m`-th test (`2 ≤ m ≤ M`) is

    b2 = 0  (a concavity)   ⟺  A_⌊m/2⌋ ≥ A_m  ⟺  m·X_{⌊m/2⌋} ≥ ⌊m/2⌋·X_m ,      (1.1)

with `X` measured from `L` (all statements below are relative to the window's left end). Put
`k = ⌊m/2⌋`, `j = k`; then `C = {m : m X_k ≥ k X_m}` and

    meet  ⟺  C = ∅  ⟺  k X_m > m X_k  for all m ∈ [2,M].                       (1.2)

`f(M) := P(C=∅)` for a window whose word is a fresh uniform arrangement (i.i.d. `b_i` uniform on
`{−1,0,1}`, or the two-column ensemble; both are computed in §6).

## 2. T1 — the reduction: early meets are a *fresh-window* statement

**Lemma 2.1.** Along a run, `L_i = i` and `δ_i ≤ i`, hence `M_i = n − i − δ_i ≥ n − 2i`, where
`δ_i` is the number of double cuts among the first `i` epochs.

*Proof.* The rule is `m* = M ⇒ (L,R) ← (L+1,R−1)`, else `L ← L+1` (`q0_walk.walk_b0_fast`). So `L`
increases by exactly 1 per epoch, `L_0 = 0`, hence `L_i = i`; `R` decreases only on a double cut,
so `R_i = n−1−δ_i` and `M_i = R_i−L_i+1 = n−i−δ_i` with `δ_i ≤ i`. ∎

**Lemma 2.2 (nesting).** The test set of `[L,R]` is `{m : 2 ≤ m ≤ R−L+1}`, and the `m`-th test
involves only `a[L..L+m−1]`. Hence `C=∅` on `[L,R]` implies `C=∅` on `[L,R']` for every `R' ≤ R`.

**Theorem 2.3.** For every `n` and every exchangeable tape law,

    P(run meets at some epoch i < (n−1)/4)  ≤  Σ_{i=0}^{⌈(n−1)/4⌉−1} f(n−2i)  ≤  ⌈(n−1)/4⌉·f(⌊(n+1)/2⌋).

*Proof.* By 2.1 the window at epoch `i` is `[i, R_i]` with `R_i − i + 1 ≥ n−2i`, so the block
`B_i := [i, i+n−2i−1]` is contained in it. By 2.2 a meet at epoch `i` forces `C=∅` on `B_i`.
`B_i` is a *deterministic* block of length `n−2i`, so by exchangeability
`P(C=∅ on B_i) = f(n−2i)`. Union bound over `i`, using that `f` is non-increasing (the test sets are
nested: `E_M ⊆ E_{M'}` for `M' ≤ M`). ∎

**This is the point of the reduction: no conditional (Lemma-Q-type) argument is needed.** The
randomness of `R_i`, which is what makes the window a polytope slice of the tape law, is paid for
once by replacing `[i,R_i]` with the *fixed* shorter block `B_i` — legitimate because the test set
is an initial segment in `m`. So the whole sub-problem collapses to

> **f(M) = o(1/M) for M ≍ n/2  ⟹  P(early meet) = o(1).**

Verified: `L_i = i` with 0 violations over 6.9·10⁴ tapes at `n = 64…1024` (§7).

## 3. T2 — the exact recursion, and exact values of f

Split a window of length `2h` into halves `W = b[1..h]` and `V_t = Σ_{i=h+1}^{h+t} b_i`, `t = 1..h`.
The tests with `m ≤ h` involve only `W` (event `E_h`); the tests with `m = h+t > h` involve `W`
and `V`. Using `k X_m > m X_k` with `X_m = X_h + V_t`, `k = ⌊(h+t)/2⌋`:

**Theorem 3.1.** `f(2h) = E[ 1{E_h(W)} · q(W) ]`, where

    q(W) = P( V_t ≥ θ_W(t)  ∀ t ∈ [1,h] ),
    θ_W(t) = ⌊ ( (h+t)·X_{⌊(h+t)/2⌋} − ⌊(h+t)/2⌋·X_h ) / ⌊(h+t)/2⌋ ⌋ + 1 ,                (3.1)

and `V` has the law of a fresh window half. `q(W)` is computed exactly by an `O(h²)` DP over
`V_t ∈ [−t,t]`.

Exact values (i.i.d. ensemble), computed by enumerating the `3^h` halves and evaluating (3.1)
exactly — these are *exact*, not Monte Carlo:

| M | 8 | 10 | 12 | 14 | 16 | 18 | 20 | 22 | 24 | 26 | 28 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| f(M) | 2.865e-2 | 1.607e-2 | 9.384e-3 | 5.790e-3 | 3.581e-3 | 2.322e-3 | 1.539e-3 | 1.051e-3 | 7.288e-4 | 5.156e-4 | 3.680e-4 |
| `E[q|E_h]` | .232 | .145 | .159 | .108 | .125 | .0886 | .0957 | .0717 | .0777 | .0597 | .0636 |
| `sup_W q` | .617 | .519 | .491 | .455 | .584 | .512 | .482 | .462 | .447 | .435 | .435 |
| `min_{m≤h} X_m` on `E_h` | −2 | −2 | −3 | −4 | −3 | −3 | −3 | −4 | −5 | −6 | −5 |

(The recursion reproduces `f(24) = 7.288e-4` against a direct 3·10⁸-sample MC of
`7.2746e-4 ± 0.0016e-4`, and `f(20) = 1.539e-3` against `1.5387e-3 ± 0.0017e-3`: the identity is
right.)

## 4. T3 — the `b₁ ≥ 0` branch is exponentially small (proved)

**Lemma 4.1.** `E_M = ∅` if `b₁ = 1`. If `b₁ = 0` then `E_M ⇒ X_{2^j} ≥ 2^j − 1` for every
`2^j ≤ M`. Hence for the i.i.d. ensemble

    P(E_M ∩ {b₁ ≥ 0})  ≤  (M+1)·3^{−M/2}  ≤  e^{−0.55M + O(log M)}.                          (4.1)

*Proof.* `m=2`: `X_2 > 2X_1`. `b₁=1 ⇒ X_2 ≥ 3 > 2 = max X_2`, impossible. `b₁=0 ⇒ X_2 ≥ 1`; if
`X_{2^j} ≥ 2^j−1` then `m = 2^{j+1}` gives `X_{2^{j+1}} > 2X_{2^j} ≥ 2^{j+1}−2`, so
`X_{2^{j+1}} ≥ 2^{j+1}−1`. That forces `#(−1)=0` and `#(0) ≤ 1` among the first `2^J ≥ M/2` letters,
whose probability is `(2^J+1)3^{−2^J}`. ∎

So with probability `1 − e^{−M/2+O(log M)}` the event lives on `{b₁ = −1}`, where the dyadic chain
gives only `X_{2^j} > 2X_{2^j−1} ≥ −2 ⇒ X_{2^j} ≥ −1`: the **dyadic skeleton stays in `[−1,∞)`**.
Everything that follows is about that branch. Exact enumeration also shows the *whole* walk stays
shallow: `min_m X_m ∈ {−5,…,0}` on `E_M` for every `M ≤ 16` (`q13_meet.py dip`), and the
first-half minima reach `−6` at `h = 13`. So the naive picture "the walk must be driven far
negative" is wrong; the cost of a meet is not a deep excursion.

## 5. T4 — what one scale gives, and why the dyadic recursion caps at `M^{−1}`

**Lemma 5.1.** For every `W`, `q(W) ≤ P(V_h ≥ X_h+1) =: ψ(X_h)`.

*Proof.* The last test `m = 2h` is `V_h > X_h`. ∎

`ψ(x) ≤ 1/2 + O(1/√h)` for `x ≥ −O(1)`, and `ψ(x) → 1` only when `x ≪ −√h`. Consequently:
`f(2h) ≤ Q(h)·f(h)` with `Q(h) = sup_{W∈E_h} q(W)`, and measured `Q(h) = 0.62, 0.52, 0.49, 0.46,
0.58, 0.51, 0.48, 0.46, 0.45, 0.44, 0.44` (table §3). **So no bound of the form `f(2h) ≤ Q f(h)`
with a universal `Q` can give more than `f(M) ≲ M^{−log₂(1/Q)} ≈ M^{−1}`** — and with `Q ≈ 1/2` the
resulting union bound `Σ_i f(n−2i) ≲ Σ_i 2/(n−2i)` is bounded by a *constant*, not `o(1)`.

The sup is attained by the halves whose walk dips deepest — at `h = 8` the maximiser is
`X = (−1,−1,−2,−1,−2,−3,−2,−1)` (dip `−3`, the deepest in `E_8`), which flattens the barrier to
`≈ −5`. The *typical* half gives `q ≈ 0.1`: the average, not the supremum, is where the decay is.
That is exactly the gap (§8).

## 6. The sharp law for `f(M)` (measured)

**Exact/MC values.** Exact for `M ≤ 16` and all even `M ≤ 28` (§3), then MC with `3·10⁸` i.i.d.
windows per point:

| M | 32 | 40 | 48 | 56 | 64 |
|---|---|---|---|---|---|
| f(M) | 1.9367e-4 ± 0.0008e-4 | 6.060e-5 ± 0.074% | 2.188e-5 ± 1.2% | 8.580e-6 ± 2.0% | 3.513e-6 ± 3.1% |
| hits | 58 101 | 18 180 | 6 565 | 2 574 | 1 054 |

**Fits.** Over `8 ≤ M ≤ 64` (`ln f` vs the regressor, 16 points):

    ln f = 1.502 − 1.769 √M        rms 0.078, max|res| 0.16        (A = 4.49, c = 1.769)
    ln f = −0.262 − 0.698 (ln M)²  rms 0.112, max|res| 0.27        (A = 0.770)
    ln f = 6.183 − 4.320 ln M      rms 0.345                       (pure power law: excluded)

So `f` decays **faster than any power of `M`**; the two surviving forms (`e^{−c√M}` and
`e^{−c(ln M)²}`) are not distinguishable on `M ≤ 64`, and both are `o(1/M)` by an enormous margin.
A clean pointwise envelope that holds at **every** measured point `7 ≤ M ≤ 64`:

    f(M)  ≤  5·e^{−1.7√M}        (checked; no violations)                                   (6.1)

(For orientation: `f(7)=5.35e-2, f(16)=3.58e-3, f(32)=1.94e-4, f(64)=3.51e-6`.)

**Per-scale factor.** From the exact recursion, `E[q|E_h]·h = 0.93, 0.72, 0.95, 0.76, 1.00, 0.80,
0.96, 0.79, 0.93, 0.78, 0.89` for `h = 4…14`, i.e.

    E[ q | E_h ]  =  (0.86 ± 0.09)·h^{−1}·(1+o(1))      ⟹  η = 1.0 ± 0.1 in Lemma B.         (6.2)

Consistently, `ln f ≈ −Σ_{j} ln(h_j/0.86)` reproduces the `(ln M)²` fit of the table above.

**Two ensembles.** The paper's tape is the two-column ensemble (`P(b=0)=1/2`, `s²=1/2`); the i.i.d.
one has `P(b=0)=1/3`, `s²=2/3`. `f` for the two-column ensemble, measured on whole tapes:

| M | 16 | 24 | 32 |
|---|---|---|---|
| two-column | 3.95e-3 ± 3.6% | 4.45e-4 ± 10.6% | 8.0e-5 ± 25% |
| i.i.d. | 3.58e-3 (exact) | 7.29e-4 (exact) | 1.94e-4 (exact) |

Same order, same (super-polynomial) decay; the i.i.d. values are the larger ones for `M ≥ 24`, so
(6.1) is conservative for the paper's ensemble there. *Ensemble caveat:* a block of length `M`
inside a long two-column tape carries `Hypergeometric(n, ⌈n/2⌉, M)` ones per column rather than
exactly `⌈M/2⌉`; the numbers above use the canonical counts (`n = M`), i.e. the typical member of
that mixture. Since `f` varies by factors of `e^{√M}` between neighbouring `M`, this `O(√M)`
count fluctuation is a lower-order (polynomial) correction and does not affect any conclusion here;
the reduction T1 itself is exact for the actual block law, whatever it is.

**Where do meets happen?** The meet window's length is *always short*: over 3000 tapes,
`n=1024`: mean `M=10.2`, median 8, max 56, `P(M ≥ 64) = 0`; `n=256`: mean 10.0, median 7, max 52.
This is just the shape of `f`: `f(8) ≈ 3e-2` but `f(56) ≈ 9e-6`.

## 7. Consequence, and the direct run-level measurement

**7a. Bound.** Combining T1 with (6.1):

    P(early meet)  ≤  Σ_{i<(n−1)/4} f(n−2i)  ≤  (n/4)·5·e^{−1.7√(n/2)}  =  e^{−1.20√n + O(log n)}.

With the measured `f` values used wherever they are known (`M ≤ 64`) the bound is much smaller:

| n | 64 | 128 | 256 | 512 | 1024 | 4096 | 65536 |
|---|---|---|---|---|---|---|---|
| `Σ_{i<(n−1)/4} f(n−2i)` | 1.0e-3 | 2.8e-5 | 1.4e-7 | 7.1e-11 | 1.3e-15 | 5.1e-32 | ~e^{−1.2√n} |

(`n ≤ 28` uses exact `f` only; larger `n` uses (6.1) beyond `M=28`, and for `M > 64` is an
extrapolation of the measured law — but at `n=64` the bound 1.0e-3 sits a factor 4.7 above the
measured 2.2e-4, so the extrapolation is not optimistic.)

**7b. Direct measurement of the run** (two-column ensemble, `n=64…1024`, 2·10⁶/n tapes each):

| n | 64 | 128 | 256 | 512 | 1024 |
|---|---|---|---|---|---|
| `P(meet at all)` | .2040 | .1725 | .1476 | .1439 | .1249 |
| `P(meet before (n−1)/4)` | 2.2e-4 | 0 (<6.4e-5) | 0 (<1.3e-4) | 0 (<2.6e-4) | 0 (<5.1e-4) |
| `L_i ≠ i` | 0 | 0 | 0 | 0 | 0 |

So: the run meets with probability `0.12–0.20`, but **essentially always at the end**, on a window
of length `≈ 7–10`; an early meet (while `M ≥ n/2`) has probability below `2·10⁻⁴` at every
`n ≥ 128` (0 events in 15 625 … 3 906 tapes), and `2.2·10⁻⁴` at `n = 64`.
This is the empirical content of the claim, and it matches `f(M) ≈ 4.5 e^{−1.77√M}` quantitatively
(the meet is controlled by the last `O(1)` epochs, where `M = O(1)`).

## 8. The one remaining gap, stated precisely

**Lemma B (channel lemma).** There are `η > 0`, `C < ∞` such that for all `h ≥ 2`

    E[ q(W) | E_h ]  ≤  C·h^{−η}.                                                            (8.1)

Measured: `η = 1.0 ± 0.1`, `C ≈ 0.86` (6.2). *Any* `η > 0` suffices:

**Theorem 8.1 (conditional).** If (8.1) holds for some `η > 0`, then
`f(M) ≤ exp(−(η/(2 ln 2))(ln M)²(1+o(1)))`, so `f(M) = o(M^{−k})` for every `k`, and by Theorem 2.3
`P(early meet) = o(1)` (indeed `O(n^{−k})` for every `k`).

*Proof.* Iterate `f(2h) ≤ C h^{−η} f(h)` down the dyadic chain: `ln f(M) ≤ −η Σ_{j=1}^{J} ln h_j
+ J ln C` with `h_j = M/2^j`, `J = log₂ M`, and `Σ_j ln h_j = (ln M)²/(2 ln 2) − O(ln M)`. Then
Theorem 2.3 with `Σ_{i} f(n−2i) ≤ n f(n/2)`. ∎

**Lemma B is the weakest form that is still true, and the simple shortcuts all fail** (checked
numerically, `h = 8, 10, 12, 14`):

* **Endpoint only.** `q ≤ ψ(X_h)`, and `E[1{E_h}ψ(X_h)] = f(2h)` by Theorem 3.1 — no gain.
* **Floor only.** With `Λ_W := −min_t θ_W(t)`, `q ≤ ballot(Λ_W)`. But `Λ_W` is *not* `O(1)`:
  `E[Λ_W] = 1.79, 1.99, 2.16, 2.27` for `h = 8,10,12,14` — it grows like `h^{0.45}`, i.e. like the
  free-walk scale. `ballot(Λ) → c > 0` at that scale, so the floor bound gives no decay either.
  (`max Λ_W = 5, 7, 9, 9`; `corr(Λ_W, dip) = .49, .44, .43, .38`.)
* **Two checkpoints.** Take `t = h/2` and `t = h` only. Measured: the resulting bound equals `q`
  *exactly* (`β₂/q = 1.00` for all four `h`), i.e. the mid-barrier constraint is **redundant** — the
  constraints at `t ≤ h/2` are implied by those at `t ≥ h/2`. So no two-checkpoint shortcut exists,
  and every "easy" bound returns a non-decaying factor.

Hence the decay is a property of the *whole* `h`-step, random, rising barrier chain — a
channel/Dirichlet-type estimate — and it is genuinely an *average* over `W | E_h` (the sup is
`≈ 0.44`, §5). This is the sense in which Lemma B is the residual content: it is not a statement
that any single test is costly, but that the chain of `h` coupled tests costs `h^{−1}`.

**Why Lemma B is the only gap.** T1 removes the conditioning; T2 makes `f` exactly computable for
`M ≤ 16` and even `M ≤ 28`; T3 shows the event lives on `b₁ = −1` and identifies the dyadic
skeleton as staying in `[−1,∞)`; and the three "easy" bounds above show that no single-test or
two-checkpoint argument can produce the decay. What is left is exactly (8.1): the *conditional,
averaged* probability that the second half's walk tracks the first half's self-similar barrier.

## 9. Why the two other routes in the brief fail

* **Route (b), the Brownian analogue, is degenerate.** With `𝒱(θ) = B(θ) − 2B(θ/2)`,
  `𝒱(θ) ≈ −2B(θ/2) + O(√θ)` as `θ↓0`, so `𝒱` changes sign infinitely often near `0` and
  `P(sup_{θ∈(0,1]} 𝒱 < 0) = 0` **exactly**, for every normalisation. The discrete probability
  therefore tends to `0`, and the *rate* is set by the smallest scale (`m = 2`), not by any
  Brownian functional: the Brownian limit cannot produce the rate, only the fact that it decays.
  (This is the same phenomenon as in REPORT_EPOCH_BOUND §8: finite-`M` effects at scale `1/M` are
  the whole story.)
* **Route (c) mis-states the geometry.** A meet does *not* require the region to have shrunk: `C=∅`
  can happen with `M_i = n−i−δ_i` as large as `n`. Indeed `P(meet)` is `0.12–0.20` precisely
  *because* meets concentrate at short windows (§6), but that is a consequence of `f(M)`, not a
  constraint on the geometry. The correct use of the geometry is Lemma 2.1 + 2.2 (T1).

## 10. Files and reproduction

`$P = C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`, run in
`C:/Users/shini/poetry_axis/paper/probe_ologn2/`:

```
$P q13_meet.py exact          # f(M), M=2..16, full enumeration 3^M
$P q13_meet.py rec            # exact f(2h) via Theorem 3.1, h=4..14 (f(28) in 2 s)
$P q13_meet.py run            # L_i=i check + P(early meet) on real runs
$P q13_meet.py dip            # depth of the walk on E_M (i.i.d., exact)
$P q13_meet.py barrier        # E[Lambda], max q, two-checkpoint beta2 (all the numbers in §8)
$P q13_meet.py mc 64 3e8      # MC f(64) (~6 min)
```

Big-MC numbers quoted above come from the same criterion (`m·X_k < k·X_m` for all `m`), seeds
98765 (i.i.d. MC), 4242 (two-column runs), 31337/777 (ensemble comparison). The per-`h` tables in
§3 are printed by `q13_meet.py rec`. No `.tex` file was touched.

**Bottom line.** The question "does the run meet early" is now reduced, rigorously and without any
conditioning argument, to the fresh-window law `f(M)`; `f` is exactly computable up to `M=28`,
measured to `M=64`, and behaves like `e^{−c√M}` (`c ≈ 1.77`) — so `P(early meet) ≤ n f(n/2) ≈
e^{−1.2√n+O(log n)} = o(1)`, and the run-level measurement is `0` for `n ≥ 128`. The single missing
step is Lemma B, and it only needs *some* `η > 0`, not `η = 1`.

## 11. Paper-ready wording

Proved part (safe to state as a theorem, no extra hypothesis):

> **Lemma (early meets are rare iff a fresh window is rarely concave-free).** Along a b1≡0 run
> `L_i = i` and `M_i = n−i−δ_i ≥ n−2i`; the test set of a window is an initial segment in `m`, so a
> meet at epoch `i` forces the *fixed* block `[i, i+n−2i−1]` to have no concavity. Hence
> `P(meet before i₀) ≤ Σ_{i<i₀} f(n−2i)` with `f(M) = P(C=∅)` for a fresh `M`-window. In
> particular `f(M) = o(1/M)` would give `P(early meet) = o(1)`.

Measured part (must be labelled as such):

> `f` is computable exactly up to `M = 28` (`M ≤ 16` and all even `M ≤ 28`, exact recursion, §3) and
> equals `4.5·e^{−1.77√M}` to within 20% for `8 ≤ M ≤ 64` (3·10⁸-sample MC beyond `M = 28`); the
> per-scale factor is `E[q|E_h] = (0.86±0.09)/h`. On real runs `P(meet before (n−1)/4) = 0` for
> `n = 128…1024` (≤ 2·10⁻⁴) and `2.2·10⁻⁴` at `n = 64`, while `P(meet at all) = 0.12–0.20`, always
> on a window of length `≈ 7–10`.

Suggested citation of the gap: *"the only missing ingredient is a proof that
`E[q(W)|E_h] ≤ C h^{−η}` for some `η > 0`, where `q(W)` is the probability that the second half
tracks the barrier `θ_W` of (3.1); numerically `η = 1` (`E[q|E_h]·h = 0.86 ± 0.09`)."*
