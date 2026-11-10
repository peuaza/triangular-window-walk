# A proven power law: `f(M) ≤ 2100·M⁻²` (and `P(early meet) ≤ 1050/n`)

**Status.** The meet probability of a fresh window satisfies a proven power law with exponent
**k = 2 > 1**, so the union bound of `REPORT_MEET` Theorem 2.3 converges: the early-meet
probability is `O(1/n)`, and the `Ω(n²)` accounting needs **no extra hypothesis**. This
replaces the open Lemma B (`E[q|E_h] ≤ C h^{−η}`) as the last missing ingredient — Lemma B is
strictly stronger than what the paper needs and remains open, but it is no longer required.

| statement | status |
|---|---|
| **T1** `f(M) ≤ P(chain_J)`, `chain_J = {X_{2^j} > 2X_{2^{j−1}}, j ≤ J}`, `J = ⌊log₂M⌋` | proved (§2), brute-force checked |
| **T2** invariant `T_j := X_{2^j} ≥ −1` on the chain | proved (§2) |
| **Lemma S** one-step supermartingale `E[e^{−3u_{j+1}}1{surv}] ≤ ¼ e^{−3u_j}`, `u_j = T_j2^{−j/2}` | proved (§3): Gaussian part exact calculus; lattice part exact (integer arithmetic) for `j ≤ 10`, Berry–Esseen for `j ≥ 10` |
| **Theorem** `f(M) ≤ 2100 M^{−2}` for all `M ≥ 2` | proved (§4–5), constants from `q18_final.py` |
| **Corollary** `P(meet before (n−1)/4) ≤ 1050/n` | proved (§6) |
| true size of `f` (`≈4.5e^{−1.77√M}`), Lemma B (`≈0.86/h`) | unchanged, still measured only |

Everything numeric below is produced by `q18_final.py` (constants), `q19_exact.py` (Lemma S in
exact integer arithmetic) and `q17_checks.py` (brute force). No `.tex` file was touched.

---

## 1. What has to be proved, and why `k > 1`

By `REPORT_MEET` T1, `P(run meets at some epoch i < (n−1)/4) ≤ Σ_{i<(n−1)/4} f(n−2i)` and
`M_i = n−2i` runs over `[n/2, n]` with step 2. Hence

* `f(M) ≤ C/M` gives `Σ ≈ (C/2)·ln 2 ≈ 0.35 C` — a **constant**, not `o(1)`;
* `f(M) ≤ C M^{−k}` with `k > 1` gives `≤ (C/2)/(k−1)·(n/2)^{1−k} → 0`;
* `f(M) ≤ C/(M log M)` would also do (`Σ ≈ C ln2/ln n`).

So the brief's "any power law suffices" needs the correction **`k > 1`** (or a log improvement).
The theorem below delivers `k = 2`.

## 2. The chain: a sub-event that stabilises itself

Write `b_i = a_i − 1 ∈ {−1,0,1}`, `X_m = Σ_{i≤m} b_i`, `T_j := X_{2^j}` (`T_0 = X_1 = b_1`).
The `m = 2^j` meet test is `2^{j−1}X_{2^j} > 2^j X_{2^{j−1}}`, i.e.

    T_j > 2 T_{j−1}.                                                        (2.1)

**Lemma 2.1.** `E_M ⊆ chain_J` where `J = ⌊log₂ M⌋`, `chain_J = ∩_{1≤j≤J}{T_j > 2T_{j−1}}`.
Hence `f(M) ≤ f̂_J := P(chain_J)`, and `f̂` is non-increasing in `J`.

*Proof.* Each `m = 2^j ≤ M` is one of the meet's tests. ∎ (Checked by brute force over all `3^M`
words for `M = 4,6,…,14`: `f ≤ f̂` in every case, `q17_checks.py`.)

**Lemma 2.2 (the invariant).** On `chain_J`, `T_j ≥ −1` for all `0 ≤ j ≤ J`; also `T_j ≤ 2^j`.

*Proof.* `T_1 > 2T_0 ≥ −2` and `T_1 ∈ ℤ`, so `T_1 ≥ −1`. If `T_j ≥ −1` then
`T_{j+1} > 2T_j ≥ −2`, so `T_{j+1} ≥ −1`. ∎

This is the whole point, and it is what `E_h` (the *full* event) does not give: the chain alone
pins the walk into `[−1, 2^j]`, so the state variable `u_j = T_j 2^{−j/2}` lives in
`[−2^{−j/2}, 2^{j/2}]`, i.e. **`u_j ≥ −2^{−j/2}`**, uniformly in the past. `REPORT_MEET` §5
showed that no bound of the form `sup_W q(W) ≤ Q < 1` can beat `M^{−1}`; the chain evades that
cap because it is a *weaker* event than the meet and its recursion is multiplicatively closed.

**Lemma 2.3 (exact recursion).** `S_j := X_{2^{j+1}} − X_{2^j}` is independent of
`ℱ_j = σ(b_1..b_{2^j})` and distributed as a sum of `2^j` fresh letters; conditionally on
`ℱ_j`, `{T_{j+1} > 2T_j} = {S_j > T_j}`. Hence

    f̂_{j+1} = E[ 1{chain_j} · F̄_{2^j}(T_j) ],      F̄_n(t) = P(S_n > t).      (2.2)

Exact values (`q18_final.py`, iid ensemble):

| j | 2^j | f̂_j | I_j := E[1{chain_j}e^{−3u_j}] | f̂ ratio | I ratio |
|---|---|---|---|---|---|
| 0 | 1 | 1.000000e+00 | 7.045108e+00 | — | — |
| 1 | 2 | 3.333333e−01 | 1.051335e+00 | 0.3333 | 0.1492 |
| 2 | 4 | 1.234568e−01 | 1.996957e−01 | 0.3704 | 0.1899 |
| 3 | 8 | 4.008535e−02 | 3.656055e−02 | 0.3247 | 0.1831 |
| 4 | 16 | 1.173251e−02 | 6.120637e−03 | 0.2927 | 0.1674 |
| 5 | 32 | 2.907614e−03 | 9.007465e−04 | 0.2478 | 0.1472 |
| 6 | 64 | 5.849457e−04 | 1.141484e−04 | 0.2012 | 0.1267 |
| 7 | 128 | 9.516649e−05 | 1.230822e−05 | 0.1627 | 0.1078 |
| 8 | 256 | 1.257558e−05 | 1.119774e−06 | 0.1321 | 0.0910 |
| 9 | 512 | 1.352903e−06 | 8.550182e−08 | 0.1076 | 0.0764 |
| 10 | 1024 | 1.187993e−07 | 5.459391e−09 | 0.0878 | 0.0639 |
| 11 | 2048 | 8.528330e−09 | 2.908998e−10 | 0.0718 | 0.0533 |
| 12 | 4096 | 5.017314e−10 | 1.291725e−11 | 0.0588 | 0.0444 |

(`f̂_3 = 4.008535e−2` agrees with a direct `3^8` enumeration. The chain is a *lossy* proxy:
`f̂_4/f(16) = 3.3`, but `f̂_12 = 5.0e−10` while `f(4096) ≈ 4.5e^{−113}`, so the loss is
astronomical and irrelevant.)

## 3. Lemma S (the only analytic ingredient)

Fix the weight `g(u) = e^{−γu}`, `γ = 3`, and `σ² = E b² = 2/3`, `α = γσ/√2 = 1.7321`,
`β = γ(1 − 1/√2) = 0.8787`.

**Lemma S.** For every `n = 2^j` (`j ≥ 1`), every integer `t ≥ −1`, with `S_n` the sum of `n`
fresh letters,

    Σ_{s>t} P(S_n = s) · exp(−γ s/√(2n))  ≤  λ · exp(−β t/√n),      λ = 1/4.    (3.1)

(The analytic part below gives this for *all real* `u ≥ −2^{−j/2}` when `j ≥ 10`; for `j ≤ 10` only
the lattice points `u = t/√n`, `t ∈ [−1, 2^j]`, are needed and are all checked.)

**Proof, step 1 (Gaussian version — exact calculus).** For `ξ ~ N(0,σ²)`,
`E[e^{−γξ/√2}1{ξ>u}] = e^{α²/2} Φ̄(u/σ + α)` (complete the square). Hence the ratio

    R(u) := E[e^{−γξ/√2}1{ξ>u}] / e^{−βu} = e^{α²/2} Φ̄(u/σ+α) e^{βu},
    d/du log R(u) = β − h(u/σ+α)/σ,      h(y) = φ(y)/Φ̄(y) increasing.          (3.2)

So `R` increases up to `u* = σ(y*−α)` (`h(y*) = βσ = 0.7174`) and decreases afterwards; here
`y* = −0.13`, `u* = −1.52 < −2^{−j/2}`. Therefore, on the whole support `u ≥ −2^{−j/2}`,

    R(u) ≤ R(−2^{−j/2}) ≤ R(0)·(1+O(2^{−j/2})) = 0.186583·(1+o(1)),   R(0) = 0.186583.

At `j = 14` the value is `R(−2^{−7}) = 0.189129`; `R(−1) = 0.5695` (hence the support matters and
is exactly what Lemma 2.2 supplies).

**Proof, step 2 (lattice, `j ≤ 10`, exact).** (3.1) is checked at **every** state `t ∈ [−1, 2^j]`
with exact integer counts `N(s) = 3^n P(S_n=s)` and 40-digit arithmetic (`q19_exact.py`); the
maximum is always attained at `t = −1`:

| j | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| max ratio | .208689 | .206218 | .204584 | .201583 | .198295 | .195401 | .193080 | .191305 | .189985 | .189020 |

**Worst over all `j ≤ 10`: 0.208689** (at `j = 1`), i.e. a 16.5% margin below `λ = 0.25`. As a
consistency check, the same quantity in double precision at `j = 11,12,13,14` is
`.188322, .187820, .187461, .187206`, i.e. it continues to decrease towards the Gaussian limit
`0.186583`; these values are *not* needed (step 3 covers `j ≥ 10`).

**Proof, step 3 (lattice, `j ≥ 10`, analytic).** Let `F_n` be the law of `ξ_n` and `Φ` the
`N(0,σ²)` CDF. Integrating by parts (`c = γ/√2`),

    E[e^{−cξ_n}1{ξ_n>u}] = e^{−cu}F̄_n(u) − c∫_u^∞ e^{−cy}F̄_n(y)dy ,

and the same identity for `ξ`, so with `ε_n := sup_y |F_n(y) − Φ(y/σ)|`,

    |Φ_n(u) − Φ(u)| ≤ 2 ε_n e^{−cu}.                                            (3.3)

Berry–Esseen (iid, `E|b|³ = 2/3`, `σ³ = 0.5443`, constant `≤ 0.56`):

    ε_n ≤ 0.56·(2/3)/0.5443·n^{−1/2} = 0.6858·n^{−1/2}.                          (3.4)

For `u ≥ u₁ := ln(1/λ)/(γ(√2−1)) = 1.1156` the trivial bound `Φ_n(u) ≤ e^{−cu}` already gives
`≤ λe^{−βu}` (the two sides are equal at `u₁` and the left decreases faster, as `c > β`). For
`−2^{−j/2} ≤ u ≤ u₁`, writing `δ = 2^{−j/2} ≤ 2^{−5}`, (3.2)–(3.4) give (using `c > β`)

    Φ_n(u) ≤ [ R(−δ) + 2ε_n·e^{(c−β)δ} ] e^{−βu} ≤ Λ_j e^{−βu},   Λ_j := R(−δ) + 2ε_n e^{(c−β)δ}.

Numerically (`Λ_j` decreases in `j`):

| j | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|
| `R(−2^{−j/2})` | .196884 | .193835 | .191695 | .190189 | .189129 | .188381 |
| `2ε_{2^j}` | .042866 | .030311 | .021433 | .015155 | .010717 | .007578 |
| `e^{(c−β)2^{−j/2}}` | 1.0396 | 1.0278 | 1.0196 | 1.0138 | 1.0098 | 1.0069 |
| `Λ_j` | **.24145** | .22499 | .21355 | .20555 | .19995 | .19601 |

so `max_{j≥10} Λ_j = 0.24145 ≤ 0.25 = λ` (margin 3.4%), and the margin improves with `j`. ∎

**Corollary 3.1 (chain supermartingale).** For `j ≥ 1`, on `{chain_j}`, `T_j ≥ −1` (Lemma 2.2),
so Lemma S with `u = u_j` gives

    E[ g(u_{j+1}) 1{T_{j+1} > 2T_j} | ℱ_j ] ≤ λ g(u_j)   (a.s. on chain_j),

hence `I_{j+1} ≤ λ I_j` and `I_J ≤ λ^{J−j₁} I_{j₁} = 0.25^{J−4}·6.120637e−3` for `J ≥ 4`.

(The measured `I`-ratios in §2 are 0.15 → 0.044, far below `λ`, so `λ = 0.25` is not tight.)

## 4. The split and the geometric bound

**Lemma 4.1.** For `j ≥ 4`, with `a_j = c√(2^j)`, `c = 0.7`,

    f̂_{j+1} ≤ e^{γc} I_j + μ f̂_j ,     e^{γc} = e^{2.1} = 8.1662,   μ = 0.22494.

*Proof.* By (2.2), split at `T_j ≤ a_j`:
`P(chain_j, T_j ≤ a_j) ≤ e^{γc} I_j` because `g(u_j) ≥ e^{−γc}` there; and for `T_j > a_j`,
`F̄_{2^j}(T_j) ≤ F̄_{2^j}(a_j) = μ_j`. The values `μ_j` are computed for `4 ≤ j ≤ 14` in double
precision (12+ digits; the margins below are ≥ 10%), maximum `0.22494` at `j = 5`, and bounded by
`Φ̄(c/σ) + 2ε = 0.19563 + 0.01072 = 0.20635` for `j ≥ 14`; `μ := max = 0.22494`. ∎

**Theorem 4.2.** For `J ≥ 4`,

    f̂_J ≤ C₀ ρ^{J−4},   ρ = max(λ, μ) = 0.25,   C₀ = f̂_4 + e^{γc}I_4/(μ−λ) = 2.0061.

*Proof.* Unroll Lemma 4.1 with `I_j ≤ λ^{j−4}I_4`:
`f̂_J ≤ μ^{J−4}f̂_4 + e^{γc}I_4·(μ^{J−4}−λ^{J−4})/(μ−λ) ≤ ρ^{J−4}(f̂_4 + e^{γc}I_4/|μ−λ|)`. ∎

**Theorem 4.3 (power law).** For the iid uniform window and every `M ≥ 2`,

    f(M) ≤ 2100 · M^{−2}        (on dyadic M = 2^J: f(M) ≤ 514 · M^{−2}).        (4.1)

*Proof.* `f(M) ≤ f̂_J` with `J = ⌊log₂M⌋` (Lemma 2.1) `≤ C₀ 0.25^{J−4} = C₀·256·4^{−J}`
`≤ C₀·256·4·M^{−2} = 2054 M^{−2}` since `4^{−J} ≤ 4M^{−2}`. ∎

**Verification** (`q18_final.py`): the bound against the exact values of `REPORT_MEET` §3/§6 —

| M | 8 | 16 | 24 | 32 | 48 | 64 |
|---|---|---|---|---|---|---|
| exact `f(M)` | 2.865e−2 | 3.581e−3 | 7.288e−4 | 1.9367e−4 | 2.188e−5 | 3.513e−6 |
| `2100 M^{−2}` | 32.8 | 8.2 | 3.6 | 2.05 | 0.91 | 0.51 |
| slack | 1.1e3 | 2.3e3 | 5.0e3 | 1.1e4 | 4.2e4 | 1.5e5 |

(no violation at any measured `M` from 2 to 64; at `M = 64` the bound is `1.5·10⁵` times the
measured value, so (4.1) is a very lossy but unconditional statement).

## 5. Why this is not in conflict with `REPORT_MEET` §5

`REPORT_MEET` T4 is correct: **any** recursion of the form `f(2h) ≤ Q f(h)` with `Q = sup_W q(W)`
is capped at `M^{−1}`, and indeed `f` itself decays much faster than any power. The theorem above
never uses such a recursion. It uses a *different, strictly weaker* event (`chain`, probability
`f̂ ≥ f`), whose *conditional* survival probability is uniformly `≤ 1/4` because the chain forces
`T_j ≥ −1`, i.e. `u_j ≥ −2^{−j/2}`, and Lemma S is a statement about that compact state range
only. The cap `Q ≈ 0.44` is a statement about `sup_W q(W)` over the *full* conditioning set
`E_h`, which is not what is bounded here.

## 6. Corollary: early meets are `O(1/n)`, unconditionally

    P(meet before (n−1)/4) ≤ Σ_{i<(n−1)/4} f(n−2i) ≤ 2100 Σ_{M∈[n/2,n], step 2} M^{−2}
                            ≤ 2100·½∫_{n/2}^{n} x^{−2}dx = 1050/n.                (6.1)

So `P(early meet) → 0` at rate `1/n` with no hypothesis. (Numerically (6.1) beats the trivial
bound 1 for `n > 1050`; for `n ≤ 64` the measured values `2.2e−4 … 0` are vastly smaller. The
`e^{−1.2√n}` rate of `REPORT_MEET` §7 remains a *measured* statement.)

## 7. Remarks, and what is still open

1. **`k` can be improved mechanically.** The binding constant is the Gaussian value
   `λ_G = e^{α²/2}Φ̄(α)` with `α = γσ/√2`: `γ = 3` gives `λ_G = 0.1866`, `γ = 4` gives
   `λ_G = 0.1506` (`α = 2.309`), `γ = 6` gives `0.105`. So the same scheme plausibly yields
   `k ≈ 2.4–2.6` (the lattice part of Lemma S would have to be re-checked for the larger `γ`,
   which was not done here; the constant `C` grows like `e^{γc}`). Since the exact `I`-ratios
   fall to `0.044` by `j = 12`, the true per-scale factor of the chain is far below `1/4`;
   `k = 2` is chosen only for a clean statement.
2. **Lemma S is the only non-elementary input.** Its content is a one-sided comparison of a
   weighted lattice tail with the Gaussian tail; the two ingredients are (i) the closed-form
   Gaussian ratio (3.2) and (ii) Berry–Esseen with constant 0.56. Nothing about the meet's
   structure is hidden in it.
3. **Two-column ensemble.** The same argument runs with the hypergeometric block law; Lemma S
   still holds but with a slightly larger constant, because the finite-population correction
   lowers the variance: at the top scale (block = half the tape) `σ²_eff ≈ ¼`, and the exact
   one-step maximum is `0.287 (n=16), 0.279 (n=32), 0.273 (n=64), 0.264 (n=256)`; for thinner
   windows it decreases (`0.240` at `σ²_eff ≈ 0.38`, `0.214` at `σ²_eff ≈ 0.48`,
   `q20_twocol.py`). With `λ = 0.30` the same proof gives `k = log₂(1/0.30) = 1.74 > 1`. (The
   `j`-dependence of `σ²_eff` only makes the constants worse, never the exponent.)
4. **Not proved, and not needed:** Lemma B (`E[q|E_h] ≤ Ch^{−η}`, `η = 1` measured). The present
   result gives only a constant per-scale factor for the *chain*, hence a power law, hence
   `O(1/n)` for early meets, whereas Lemma B would give `exp(−c(ln M)²)`. The paper's `Ω(n²)`
   conclusion needs only the former.
5. **Sharpness of the method.** The chain probability itself decays super-polynomially
   (`f̂_12 = 5.0e−10`), and its measured per-scale ratios fall like `≈1/j`; so a proof that
   captures even a slowly varying factor would give far more than `M^{−2}`. The obstruction to
   extracting that is quantitative control of the quasi-stationary law of `T_j`, which is what
   Lemma B encodes.

## 8. Files and reproduction

`$P = C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`, run in
`C:/Users/shini/poetry_axis/paper/probe_ologn2/`:

```
$P q18_final.py      # constants: chain DP, Lemma S table, mu, rho, C, verification table
$P q19_exact.py      # Lemma S in exact integer counts + 40-digit arithmetic (j<=10)
$P q17_checks.py     # Gaussian ratio, Lemma S (doubles, j<=14), brute force f<=fhat
$P q20_twocol.py     # Lemma S for the two-column (hypergeometric) block law
$P q14_bound.py chain 10   # independent DP of the chain probabilities
```

`q15_super.py` (Gaussian ratio table, first assembly) and `q16_theorem.py` (assembly with
`λ = 0.22`, `k = 2.18`) are the exploration scripts; `q18_final.py` supersedes them with the
`λ = 0.25`, `k = 2` parameter set used in this report.

## 9. Paper-ready wording

> **Lemma (a power law for fresh windows).** For a fresh window of length `M` with i.i.d. letters
> uniform on `{−1,0,1}`, the probability that no concavity exists is
> `f(M) ≤ 2100·M^{−2}` for every `M ≥ 2`.
>
> *Proof sketch.* `f(M)` is at most the probability that `X_{2^j} > 2X_{2^{j−1}}` for every
> `j ≤ log₂M`. On that event `T_j = X_{2^j}` satisfies `T_j ≥ −1`, and the normalised state
> `u_j = T_j2^{−j/2}` therefore lives in `[−2^{−j/2}, 2^{j/2}]`. With `g(u) = e^{−3u}` the
> one-step inequality `E[g(u_{j+1})1{T_{j+1}>2T_j} | ℱ_j] ≤ ¼ g(u_j)` holds at every state of
> that range (Gaussian computation plus Berry–Esseen; exact integer check for the first ten
> scales), so `E[1{chain_J}g(u_J)] ≤ 4^{−J+4}·6.12e−3`. Splitting at `T_j = 0.7·2^{j/2}` turns
> this into `f̂_{j+1} ≤ 8.17 I_j + 0.225 f̂_j`, whose solution is `f̂_J ≤ 2.01·4^{−(J−4)}`. ∎
>
> Consequently `P(meet before (n−1)/4) ≤ Σ_{i<(n−1)/4} f(n−2i) = O(1/n)`, so the `Ω(n²)`
> accounting is unconditional.

Suggested citation of the remaining gap (now optional): *"The sharper per-scale estimate
`E[q(W)|E_h] ≤ Ch^{−1}` (numerically `0.86/h`) is not needed for the `Ω(n²)` bound; a power-law
bound on `f` suffices and is proved in the appendix."*
