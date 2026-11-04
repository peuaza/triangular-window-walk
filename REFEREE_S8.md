# Referee report on §8 of `bthm_merged.tex` — the Θ(n²) chain

**Recommendation: major revision.** The skeleton of the chain is sound, but as written the section
proves less than it claims and contains one outright numerical error, one unproved extension, and a
paragraph that says the opposite of the paragraph before it.

**What was read.** §8 of `bthm_merged.tex` (lines 610–756), `probe_ologn2/REPORT_EPOCH_BOUND.md`,
`probe_ologn2/REPORT_LEMMA_Q.md`. **What was checked.** The covariance identity was re-derived
symbolically and checked exactly per multiset; the Gaussian constant was re-evaluated by quadrature
and by Brownian simulation; the deterministic chain was checked **exhaustively on every word in the
support of the ensemble (`#2=#0`) for n = 8, 10, 12 (83,849 tapes)**. Scripts: `probe_ologn2/ref_lemma4.py`,
`ref_chain.py`, `ref_suffix.py`.

---

## 0. Verdict table

| link | statement | verdict |
|---|---|---|
| (1) | ensemble, `b1≡0`, `X_m`, `A_m`, `k=⌊m/2⌋`, `b2=1 ⟺ A_k<A_m` | **PROVED** (with a caveat on which ensemble) |
| (2) | `P(cost₀ ≥ cM) ≥ p(c) − o(1)`, `p(c)=P(sup_{[1−c,1]}𝒱<0)`, hence `E[cost₀] ≥ .008M` | **PROVED except the constant — the last step is BROKEN by a factor 2** |
| (3) | frontier inequality `cost(L,R) ≥ cost(L,R″) − (R″−R)` | **PROVED** (but the paper's gloss of it is false) |
| (4) | `cost_i ≥ cost(L_i,n−1) − δ_i`, `δ_i ≤ i`, suffix "unconditioned" | **GAP**: (a) proved; (b) the prose is false but the inequality survives; (c) unbalanced suffix unproved |
| (5) | summation and `N = c₀n/(1+c₀)` giving the bracket | **PROVED as algebra; GAP as a theorem** (`q_N` unproved; constant wrong) |

The chain is **not** broken in its logical skeleton. It is broken in three places: the constant in
(2), the justification (not the inequality) in (4), and the unconditional claim in (5). The
following paragraph of the paper also contradicts (5) outright.

---

## 1. Link (1) — the ensemble and the reduction: PROVED

The reduction `b1≡0` is Theorem 2 and is checked on every admissible sequence for `n=6,…,13`; nothing
to attack. `b2=1 ⟺ A_{⌊m/2⌋} < A_m` is Lemma 1, checked on 20,000 windows. The identities
`A_m = 1 + X_m/m`, `k=⌊m/2⌋`, and the cost `cost = M − m* + 1` (`M−1` on a meet) are consistent with
the implementation.

**Caveat (definitional, not fatal).** §2 of the paper defines the ensemble as *two independent
uniform `⌈n/2⌉`-subsets summed* (`s² ≈ 1/2`); the brief for this report describes the ensemble as
*uniform on the slice `#2−#0 = 2⌈n/2⌉−n`* (`s² ≈ 2/3`). These are **different** ensembles (the
slice weights count vectors by their multinomial coefficient and concentrates at `(⅓,⅓,⅓)`; the
two-column ensemble concentrates at `(¼,½,¼)`). Both are exchangeable with uniform conditional
arrangements and both have `s² ≥ 1/4`, so the single-epoch theorem covers both — but the paper should
say which one it means, and any measurement of `c₀` is ensemble-specific (`.134` vs `.127` in the
report's §7F).

---

## 2. Link (2) — the single-epoch bound

### 2.1 The FCLT, the covariance identity and the perturbation: PROVED

- **Exact covariance.** `Cov(X_i,X_j) = s² i(M−j)/(M−1)` for `i ≤ j` is the sampling-without-replacement
  identity and holds for **any** multiset (the mean cancels). Expanding the four terms gives
  `Cov(𝒟(θ),𝒟(θ′)) = (2θ∧θ′ − θ∨θ′)·M/(M−1)` for `θ′ < 2θ`, and `0` for `θ′ ≥ 2θ`. I checked this
  exactly per multiset: `(κ₀,κ₁,κ₂) = (2,4,2)`, `(1,3,4)`, `(0,4,4)` at `M=8`, max error `3.3·10⁻¹⁵`.
  The diagonal is `θ·M/(M−1)`, so the limit has `Var 𝒱(θ)=θ`, as claimed.
- **The bridge cancels exactly.** `𝔅(θ) − 2𝔅(θ/2) = 𝒱(θ)` because the conditioning term is `−θW(1)`
  and `θ − 2(θ/2) = 0`. The same computation shows the **drift** cancels too: for an unbalanced
  multiset `X_j = μj + Y_j`, the even-`m` contrast is exactly `Y_m − 2Y_{m/2}` (the `μ` term is
  `μ(m−2(m/2))=0`), and for odd `m` the `μ` cancels between the two sides of the test
  (`X_m − 2X_k > X_k/k` becomes `Y_m − 2Y_k > Y_k/k`). So the limit is the same `𝒱` for every
  multiset with `s² > 0`.
- **Odd-m / floor perturbation.** `|X_k/k| ≤ 1`, so the threshold is `≤ 1/(s√M)`; the odd-`m`
  contrast differs from `𝒱_M(m/M)` by `b_m/(s√M)`; the polygonal interpolation costs another
  `1/(s√M)`. `O(1/(s√M))` is right.
- **Transfer to the discrete scan.** Lemma 3(a)–(c) then gives `b2=1` at every `m ≥ (1−c)M` on
  `{inf_{[1−c,1]}𝒱 > ε}`, and Portmanteau on the open set `U_ε` gives
  `liminf P(cost ≥ cM) ≥ P(inf 𝒱 > 0) = P(sup 𝒱 < 0)` by `𝒱 =_d −𝒱`. This is correct.

### 2.2 The constant: **BROKEN — a factor 2 is missing**

The report's equation (★★) says

```
p(c) ≥ 1 − E_{b∼|N(0,1)|}[ exp( −((b(1−2c) − m_c)_+)² / (9/8) ) ].
```

The table in §5 (`.246, .161, .078, .031, .009, .0014, ~0`) is exactly the right-hand side of this
display: I reproduce it to four decimals (`.2458, .1608, .0775, .0312, .0089, .0014, .0001`). **The
display is wrong by a factor 2.** The correct computation is:

`p(c) = P(inf_{[1−c,1]} 𝒱 > 0)`. Since `1 ∈ [1−c,1]`, this event implies `𝒱(1) > 0`, so
`p(c) = P(inf>0, 𝒱(1)>0) = ½ · E_{b∼|N(0,1)|}[ P(inf>0 | 𝒱(1)=b) ]`. The conditional bound
`P(inf>0|b) ≥ 1 − exp(…)` therefore gives

```
p(c) ≥ ½ ( 1 − E_{b∼|N|}[ exp(…) ] ),
```

not `1 − E`. Equivalently: with `b ∼ N(0,1)` unrestricted the integrand is `1` on `b<0`, so
`1 − E_{N(0,1)}[exp] = ½(1 − E_{|N|}[exp])` — the same half. So the display computes `2p(c)`, not
`p(c)`; read as an inequality for `p(c)` it is false by a factor 2. The true `p(.05)` is `.381`–`.388`
(report §5 and my simulation), so the **statement** `p(.05) ≥ .16` happens to be true; the **proof
printed** establishes only `.080`.

Consequences, all of which propagate into §8:

| quantity | paper | corrected |
|---|---|---|
| `p(.05)` | `.16` | **`.080`** |
| `c₀ = .05·p(.05)` | `.008` | **`.00402`** |
| `N*/n = c₀/(1+c₀)` | `.00794` | **`.00400`** |
| `c₀²/(2(1+c₀))` | `3.2·10⁻⁵` | **`8.05·10⁻⁶`** |
| `q_N` threshold `c₀/2` | `.004` | **`.00201`** |

The error is repairable in two ways, and the paper should say which it uses: (i) halve the constant
and accept `E[cost₀] ≥ .004M`; or (ii) do the sharper Gaussian computation. The report's `9/16` is the
sup of `Var(𝒱(θ)|𝒱(1))` over **all** `θ`, attained at `θ=5/8`, which is **outside** `[1−c,1]` for
every `c<3/8`. On `[.95,1]` the sup is `θ−(2θ−1)² = .14`, and the conditional process has
`E sup Z ≈ .28` (simulation), not the Dudley value `.485` computed with the unconditional metric.
With `σ²=.14` and `m=.28` the same Borell step gives `p(.05) ≈ .21`, enough to restore `.16` with room
to spare — but this uses the simulated `E sup Z`, so a proved Dudley bound for the conditional metric
would be needed; that computation is not in the paper. As printed, the table is the un-halved one, and
the paper's `.008` rests on it.

---

## 3. Link (3) — the frontier inequality: PROVED

For `L+2 ≤ R ≤ R″`, `C(R) = C(R″) ∩ [2, R−L+1]` (the tests for `m ≤ M` read only cells `≤ R`), hence
`max C(R) ≤ max C(R″)` when both are non-empty and

`cost(L,R) = M − max C(R) + 1 ≥ M − m*″ + 1 = cost(L,R″) − (R″−R)`.

If `C(R)=∅` the window meets and `cost(L,R)=M−1`; if also `C(R″)=∅` the two sides are equal, and if
not, `M−1 ≥ M+1−m*″` because `m*″≥2`. So the inequality holds in the meet case too. I checked it at
every realised epoch of every tape for `n=8,10,12`: **0 violations in 83,849 tapes**.

**But the paper's gloss is false.** Line 713 says "extending a window can only add tests". Extending
can **remove** tests: on the word `a = [0,2,2,0]`, `cost(0,2) = 2` (the window `[0,2]` meets) while
`cost(0,3) = 1` (the test at `m=4` succeeds). The cost drops by one when the window grows by one. The
frontier inequality is a *bound* on how much it can drop, not a monotonicity statement; the `δ_i` loss
is doing real work precisely because of this. The sentence should be deleted.

---

## 4. Link (4) — the key step

The pointwise inequality `cost_i ≥ cost(i,n−1) − δ_i` follows from link (3) with `R=R_i`, `R″=n−1`.
The interesting questions are (a), (b), (c).

### (a) `δ_i ≤ i`: **PROVED**

`L` advances by exactly 1 per epoch; `R` decreases by 1 on a double cut and 0 otherwise. Hence
`δ_{i+1}−δ_i ∈ {0,1}`, `δ_0=0`, so `δ_i ≤ i` and `R_i = n−1−δ_i ≥ n−1−i`. No tape violates it; the
exhaustive check confirms `L_i=i`, `R_i=n−1−δ_i`, `δ_i≤i` at every epoch for `n=8,10,12`. The bound is
attained by the tape `[0]^k[2]^k[1]^{n−2k}`, which double-cuts every epoch (`δ_i=i`).

### (b) "The suffix is unconditioned": **the prose is false; the inequality survives**

`L_i = i` is deterministic, not an event. The event that epoch `i` exists,
`A_i = {epochs 0..i−1 realised}`, is **not** a function of the prefix `a[0..i−1]`: epoch 0's test at
`m=n` compares `A_{n/2}` and `A_n`, which reads **every** cell of the tape, including all cells above
`n−1−i`. So `A_i` constrains the suffix, and the suffix is *not* unconditioned given the history. The
suspicion that the suffix is conditioned is correct at the level of the event.

It does **not** break the inequality, because the proof never uses the conditional law. What is used is

```
E[cost_i 1{A_i}] ≥ E[cost(i,n−1) 1{A_i}] − i
                  ≥ E[cost(i,n−1)] − (n−i) P(A_i^c) − i
                  ≥ c₀(n−i)(1−o(1)) − (n−i) q_N − i,
```

and the second line is the trivial bound `cost(i,n−1) ≤ n−i`, valid for any event `A_i` whatsoever.
The price of the conditioning is exactly the `(n−i)q_N` term, which the paper's display does contain
(via `nNq_N`). So the step is valid as an inequality; the sentences at lines 720–721 are not.

**Failing tape for the pathwise reading.** Take `n=512`, `a = [0,0,0,2,2,2,1,…,1]` (506 ones).
This tape is in both ensembles (`#2=#0=3`), it double-cuts every epoch (255 epochs, every `cost_i=1`),
so `A_6` holds, and yet `cost(6,511) = 1` while `c₀(512−6) = 4.05` (with the paper's `c₀=.008`;
`2.02` with the corrected `.004`). Hence

```
cost(i,n−1) ≥ c₀(n−i)     is FALSE pathwise,
```

which is exactly why the paper needs an expectation. What it must not do is claim the *conditional*
expectation: the object bounded is `cost_i·1{A_i}`, and the paper's inequality handles the dependence
through the marginal `E[cost(i,n−1)]` plus the `(n−i)q_N` correction. That is legitimate, but it means
the sentence "the object bounded below is an unconditioned suffix" describes the *comparison*
variable, not the object bounded. `cost(i,n−1)` is unconditioned only in the marginal sense, and the
conditioning is paid for, not absent.

### (c) Unbalanced suffixes: **GAP (repairable, but not proved)**

This is the real hole in link (4), and it is not the obvious one.

The single-epoch theorem is stated and proved for a multiset with `κ₀ = κ₂` (a balanced tape). A
suffix of the tape does **not** have balanced counts: `#2−#0` on `a[i..n−1]` equals `−X_i` (plus the
total `0` or `1`), which is generically non-zero. So "the single-epoch bound applies to it unchanged"
is literally false: the hypothesis of the theorem is not satisfied.

The report's patch (REPORT_EPOCH_BOUND §2, Remark) is not a proof:

> "An unbalanced multiset adds the drift term `μθ` to `𝔅`; the contrast kills it too, since
> `μ(θ − 2·θ/2) = 0`."

The drift in the *normalised* process `X_{θM}/(s√M)` is `(μ/s)√M·θ`, not `μθ`; it diverges, and the
displayed cancellation is dimensional nonsense. The **conclusion** is nevertheless correct, for a
different reason:

1. **The test is translation-invariant.** Adding a constant to every cell shifts `A_k` and `A_m`
   equally, so `cost` is unchanged. The law of `cost` for a multiset with mean `μ` is the law for its
   centred version `{−μ, 1−μ, 2−μ}`.
2. **The centred prefix is a bridge.** For the centred population, `Cov(Y_i,Y_j)=s²i(M−j)/(M−1)` and
   the WOR CLT gives `Y_{θM}/(s√M) ⇒ 𝔅`; the contrast of `𝔅` is `𝒱`.
3. **The exact covariance identity of Lemma 2(i) already holds for any multiset** (it only uses the
   WOR covariance formula). I verified it exactly for `(1,3,4)` and `(0,4,4)` as well as `(2,4,2)`.
4. **Numerically the cost law is the same**: at `M=256`, `E[cost]/M = .136` for `(64,128,64)`,
   `.141` for `(0,128,128)`, `.134` for `(128,128,0)`, `.141` for `(32,192,32)`, `.135` for
   `(0,64,192)`. The last has `s²=.188 < 1/4` and still gives the same constant.

So the gap is one lemma wide: add "any multiset with `s² ≥ 1/4`" to the theorem, and prove the
centring/translation invariance. But as submitted, §8's "applies unchanged" is an unproved extension,
and the paper's "it is the one hypothesis left" (line 729) is therefore false: there are two.

A second, smaller unproved input is `s² ≥ 1/4` for the suffix. It fails on the all-ones suffix
(`s²=0`) but with probability `e^{−cM}`, which is inside the `o(1)`; the paper should say this.

---

## 5. Link (5) — the summation and the optimisation

### (d) The truncation: **PROVED, but only with the indicators**

`T = Σ_i cost_i 1{A_i}` is exact. For `i<N`, `P(A_i^c) ≤ q_N` (the walk cannot restart), and on the
meet event the lost terms are at most `n(N−j) ≤ nN`, so the `nNq_N` subtraction is sufficient. The
paper's display hides the indicators; the derivation needs them, but the mathematics is correct.
Nothing goes wrong on the meet event: the terms after it are zero and the `q_N` term pays for them.

### (e) The algebra: **PROVED**

`Σ_{i<N}[c₀(n−i) − i] = c₀(nN − N(N−1)/2) − N(N−1)/2`. With `N=αn`,
`f(α) = c₀(α − α²/2) − α²/2`, `f′(α) = c₀ − (1+c₀)α`, so `α* = c₀/(1+c₀)` and
`f(α*) = c₀²/(2(1+c₀))`. The displayed bracket is correct. (The exact optimum with the `q_N` term is
`(c₀−q_N)²/(2(1+c₀))` at `α=(c₀−q_N)/(1+c₀)`, so the paper's linearisation is *conservative*: the
bracket is positive for `q_N < c₀`, not only `q_N < c₀/2`. This is a point in the paper's favour, not
against it.)

### The theorem: **GAP, because `q_N → 0` is not proved**

`q_N = P(the walk meets before epoch N)`. The paper admits it is unproved, and it is the only input
that is *labelled* unproved. But the heading "The run costs Θ(n²)" and the sentence "this makes the
run quadratic in expectation" assert the conclusion unconditionally. They should not.

Two further arithmetic consequences of §2.2: with the corrected `c₀=.00402`, the threshold is
`q_N < .00201`, not `.004`, and the rigorous constant is `8.05·10⁻⁶`, not `3.2·10⁻⁵`. The claimed
`7.9·10⁻³` uses the *measured* `c₀=.133`, not a proved one; it is a measurement inserted into a
theorem-style display (and `0.133²/(2·1.133)=7.81·10⁻³`, so even the rounding is off).

---

## 6. (f) Every overclaim in §8, with replacements

| line | sentence (abridged) | why it overclaims | suggested replacement |
|---|---|---|---|
| 615 | "three lemmas and a theorem, **all machine-checked**" | They are checked numerically/exhaustively, not by a proof assistant. | "three lemmas and a theorem, checked by exhaustive enumeration (`n=6,…,13`) where stated and by simulation elsewhere" |
| 683–684 | "What is **not proved** is the matching lower bound" | Contradicts line 731 ("makes the run quadratic"). One of the two must go. | "The matching lower bound is proved only modulo `q_N→0` (§8)." |
| 696–700 | "gives `p(.05) ≥ .16`, hence `E[cost₀] ≥ .008M`" | The printed Gaussian computation gives `.080`, not `.16` (§2.2). | "gives `p(.05) ≥ .08`, hence `E[cost₀] ≥ .004M`" — or add the sharper conditional-variance computation that recovers `.16`. |
| 701 | "the only one in the paper that goes below `n`" | Unclear: `.008n` is below `n`, but the sentence does not say why that is a virtue or what it is contrasted with. | "the only nontrivial lower bound in the paper" (if that is the meaning — say it) |
| 703 | "The limit itself is **sharp**" | "Sharp" means equality; what is shown is a 1% match in simulation. | "The limit is attained to within 1% in simulation at `n=4096`" |
| 710 | "**The run costs Θ(n²)**" | Only `O(n²)` is unconditional; `Ω(n²)` is conditional on `q_N→0` and on the unbalanced-suffix extension. | "The run is `O(n²)`; it is `Θ(n²)` in expectation modulo `q_N→0`" |
| 710 | "The two halves do meet" | Unintelligible in context; no referent for "the two halves". | delete |
| 713 | "extending a window **can only add tests**" | False: `[0,2,2,0]` has `cost(0,2)=2`, `cost(0,3)=1`. | "extending a window by `R″−R` cells can lower its cost by at most `R″−R`" |
| 718–720 | "the single-epoch bound … applies to it **unchanged**" | The theorem is stated for balanced multisets; a suffix has `#2−#0 = −X_i`. | "applies to it once the single-epoch theorem is extended from balanced multisets to arbitrary ones with `s²≥1/4` (translation invariance)" |
| 720–721 | "the object bounded below is an **unconditioned** suffix" | `A_i` depends on the whole tape, including cells above `n−1−i`; the bounded object is `cost_i·1{A_i}`. | "the conditioning enters only through `P(A_i) ≥ 1−q_N`; no conditional law of the window is used" |
| 726–729 | "it is **the one hypothesis** left in the lower bound" | Two hypotheses: `q_N→0` and the unbalanced-suffix extension. | "it is one of two hypotheses; the other is the extension of the single-epoch theorem to unbalanced suffixes" |
| 727 | "positive as soon as `q_N < .004`" | Tied to the wrong `c₀`. | "`q_N < .002`" |
| 730 | "the constant is `3.2×10⁻⁵`" | Should be `8.05×10⁻⁶`. | "`8.0×10⁻⁶`" |
| 731 | "this **makes** the run quadratic in expectation" | It makes it quadratic *modulo* two unproved inputs. | "modulo `q_N→0`, this makes the run quadratic in expectation" |
| 734–735 | "Passing from one epoch to the whole run **needs** the same inequality at every epoch, conditioned on the history" | False after the truncated argument: that is exactly what the frontier truncation avoids. | delete, or "the truncated argument avoids the per-epoch conditional bound" |
| 739 | "**The passage is open.**" | Directly contradicts lines 710–731. | delete, or "the remaining passage is `q_N→0` plus the unbalanced-suffix extension" |
| 740 | "What is proved is that the first epoch contributes `Ω(n)`" | Understates and is inconsistent with the preceding paragraph. | "what is proved unconditionally is `Ω(n)`; the `Ω(n²)` is conditional on `q_N`" |
| 740–741 | "the geometry contributes `Σ_i M_i ≥ n²/4−O(n)` over the `Θ(n)` epochs" | Only on the event that the first `(n−1)/2` epochs are realised, i.e. modulo `q_{(n−1)/2}`. | "on the event that the first `(n−1)/2` epochs are realised, `Σ_i M_i ≥ n²/4−O(n)`" |
| 736 | "such windows occur **at every epoch**" | False: a descending window is exactly a double cut; left-cut epochs are not descending (about 62% are). | "such windows occur at every double-cut epoch, about 62% of epochs" |

---

## 7. The internal contradiction (the most damaging defect)

Lines 710–732 claim the passage is made: truncated frontier + unconditioned suffix + `q_N` small.
Lines 734–741 then say:

> "Passing from one epoch to the whole run needs the same inequality at every epoch, conditioned on
> the history, and that is false pathwise … **The passage is open.** What is proved is that the first
> epoch contributes `Ω(n)` …"

Both paragraphs cannot stand. The second one is the *old* text, written before the truncation
argument, and it is now **wrong**: the truncated argument does not need a per-epoch conditional bound,
and its own explanation ("the window word is a polytope slice, so the convergence does not apply") is
irrelevant to the new route, which never applies the convergence to the window word. A referee will
read §8 as claiming both `Ω(n²)` and "the passage is open", and will reject on that alone. Delete
734–741 (or rewrite it as the list of the two remaining hypotheses).

---

## 8. What is actually proved, stated correctly

For the `b1≡0` rule on the paper's (two-column) ensemble, for every `c∈(0,1)` and every `N ≤ n/2`,

```
E[T] ≥ c₀ (nN − N(N−1)/2)(1−o(1)) − N(N−1)/2 − nN q_N,
c₀ = c·p(c),   p(.05) ≥ .080  (not .16),   so c₀ = .00402,
N* = ⌊c₀n/(1+c₀)⌋ = .00400 n,
E[T] ≥ [8.05·10⁻⁶ − .00400 q_N] n² (1−o(1)),   positive for q_N < .00201.
```

Unproved inputs: (i) `q_N → 0`; (ii) the extension of the single-epoch theorem to suffixes with
`#2−#0 ≠ 0` and `s² ≥ 1/4`. The first is a self-contained combinatorial problem about monotone
profiles; the second is a one-paragraph centring argument. Neither is in the paper.

With the *measured* `c₀=.133` the constant would be `7.8·10⁻³ n²`, but that is measurement, not
theorem, and should not appear in the same display as the proved bound.

---

## 9. Verification log

```
$P = C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe
cd C:/Users/shini/poetry_axis/paper/probe_ologn2

$P ref_lemma4.py    # quadrature of (★★) and its halved version; Brownian p(c)
$P ref_chain.py     # exhaustive n=8,10,12: 83,849 tapes, 0 violations of
                    #   L_i=i, R_i=n-1-delta_i, delta_i<=i,
                    #   cost_i >= cost(i,n-1)-delta_i,
                    #   T >= sum_{i<N} cost(i,n-1) - N(N-1)/2
                    # exact covariance identity at M=8
$P ref_suffix.py    # unbalanced multiset cost law; all-1s suffix; conditional means
```

Key outputs quoted above:

- quadrature: `1−E_{|b|}[·] = .1608` at `c=.05`; halved `.0804`; true `p(.05)=.388` (30,000 paths,
  grid 4096; report's `.381`).
- exhaustive chain: `n=8: 1107 tapes`, `n=10: 8953`, `n=12: 73789`; `geometry_bad=0`,
  `frontier_bad=0`, `sum_bad=0`.
- covariance: per-multiset max error `3.3·10⁻¹⁵`; `(2,4,2)` diagonal `2/7, 4/7, 6/7, 8/7`.
- unbalanced cost law at `M=256`: `.136, .141, .134, .141, .135` for the five multisets listed in
  §4(c).
- all-1s suffix at `n=512`: `cost(6,511)=1`, `c₀(n−6)=4.05` (paper) / `2.02` (corrected).
- monotonicity counterexample: `[0,2,2,0]`, `cost(0,2)=2`, `cost(0,3)=1`.

No `.tex` file was modified.
