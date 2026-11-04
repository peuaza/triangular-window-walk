# Route (e): the dyadic-contrast route

Fifth attempt at the O(log n)-bit / n·polylog(n)-read separation for the walk
rules. `LOWER_BOUND_ATTEMPT.md` closed routes (a) crossing, (b) behaviour
counting, (c) adversary; `LOWER_BOUND_ROUTE4.md` closed route (d), the aggregate
over all cuts with jump charging, and proved the conservation law
`T·J̄·s ≥ Σ_k CC_k − ns`. This file does not use cuts at all: it uses what the
four earlier routes left on the table, namely that the `b2` test is a family of
**dyadic contrasts** on the prefix sums.

Everything below is about the same function `delta` on the same ensemble `A_n`
as the earlier files, with one standing caveat stated once and used throughout:

> **Caveat (b1).** The theorems are stated for `delta_0`, the rule with `b1`
> forced to 0. `delta_0 = delta` on the 85–90 % of tapes on which `b1` never
> fires (measured 34–36 of 40 tapes at n = 64…512, `REPORT.md` §3), because on
> those tapes the two trajectories coincide; the bounds therefore transfer to
> `delta` on that sub-ensemble. For the remaining tapes the epoch structure
> survives with `O(1)` extra probes per epoch (b1 can fire only at `m < 4ρ`,
> `ρ ≤ 6`, measured), which changes no order. Nothing below depends on the
> `b1` correction.

New probe: `lb_route5.py`, output `lb_route5.log`.

---

## 0. Verdict

**Route (e) produces the first superlinear bound in this project that is not a
movement bound — but only for machines whose memory is additive, and that
restriction provably cannot be lifted by any argument in the file.** What
survives unconditionally is a normal form and four ceilings, each a theorem.

| item | statement | status |
|---|---|---|
| **L5.1** normal form | `b2 = 1 ⟺ mean_z(second block) > mean_z(first block)`; for `m = 2k` this is the sign of the Haar coefficient `Σ_{[L+k,L+2k)} z − Σ_{[L,L+k)} z` | **proved** + 200k windows checked |
| **L5.1a** shift invariance | `delta` depends only on `z = a−1`, and `A_n ≅ {z ∈ {−1,0,1}^n : Σz = −ε}` | **proved** |
| **L5.2** epoch form | the walk is a sweep: epoch `L` probes `m = M, M−1, …` at the *single* anchor `L` and stops at the largest negative `m*(L)`; halt at `L` if there is none; `(L,R) ← (L+1,R)` or `(L+1,R−1)` | **proved** + all 111,322 tapes n ≤ 12, and all 201,643 at n = 13 |
| **L5.3** anchor distinctness | every contrast `(L,m)` is asked **at most once** in a run; no contrast *answer* is ever reused | **proved** (corollary of L5.2) |
| **L5.4** one pass per epoch | `m*(L)` is computable by a single left-to-right pass over `[L,R]` with **3 counters** | **proved** + every epoch of every tape n ≤ 13 |
| **T5.5** checkpoint theorem | for a machine whose memory is `r` exact interval sums, `T ≥ Σ_e max(0, ℓ_e − 2r)`; measured `Σ_e ℓ_e = 0.056–0.117 n²`, so `T = Ω(n²)` for `r = o(n)` | **proved** (model), measured constants |
| **C5.6** | in that model `s = O(log n) ⟹ T = Ω(n²)`; sub-quadratic reads require `r = Ω(n)`, i.e. `s = Ω(n log n)` | **proved** (model) |
| **C5.7** | hence a `p`-pass machine needs `p = Ω(n)` passes | **proved** (model) |
| **C5.8–C5.11** four ceilings | counting `≤ n−1`; sign-vector counting over-counts by `2^{Θ(n)}`; movement `≤ O(n)`; information `≤ n` | **proved** |
| **C5.13** | any `s = O(log n)` machine with `T = o(n²)` must hold memory that is *not* additive — the analogue of route (d)'s "must teleport" | **proved** (contrapositive of T5.5) |
| the conjecture | `T = ω(n)` for `s = O(log n)` | **open**, and route (e) does not prove it |

The shape of the failure is new and worth stating at the top, because it is
quantitative rather than qualitative:

> The dyadic route's *invalid* form gives, per epoch, `H(m*) ≈ log₂ log₂ n` bits
> (measured, §3), so summing over the `Θ(n)` epochs would give `O(n log log n)`
> reads — superlinear, so enough to refute an `O(n)`-read simulator, though short
> of `n polylog n`. The *valid* form of the same argument gives
> `H(A_1…A_E)/log₂3 ≤ n`, because all epoch answers are functions of the same `n`
> cells. The whole difference between the two is the correlation between epoch
> answers, measured at `H(joint)/Σ H(marginal) = 0.62` after six rounds (§3.3).
> **The proof obligation of route (e) is exactly "the answers are effectively
> independent for a small-state machine", and that is a branching-program
> statement.**

---

## 1. The normal form

### 1.1 Notation

`A_n = {a ∈ {0,1,2}^n : #0(a) − #2(a) = ε}`, `ε = 0` (n even), `1` (n odd).
`z = a − 1`, so `z ∈ {−1,0,1}^n` and `Σ z = #2 − #0 = −ε`. `Z_j = Σ_{i<j} z_i`,
`Z_0 = 0`. A *window* is `[L, L+m)`, `M = R−L+1` the region size, `k = ⌊m/2⌋`.

### 1.2 Lemma 5.1 (block-mean form; proved)

For every window with `m ≥ 2`, writing `P = Σ_{i∈[L,L+k)} z_i` and
`Q = Σ_{i∈[L+k,L+m)} z_i`,

```
b2([L,L+m)) = 1   ⟺   k·Q > (m−k)·P   ⟺   mean_z(second block) > mean_z(first block).
```

*Proof.* `b2 = 1 ⟺ (Σ_W a)·k > (Σ_{[L,L+k)} a)·m`. Now `Σ_W a = m + P + Q` and
`Σ_{[L,L+k)} a = k + P`. So `b2 = 1 ⟺ (m+P+Q)k > (k+P)m ⟺ mk + kP + kQ > km + Pm
⟺ kQ > P(m−k)`. ∎

Consequences, both used below:

* For `m = 2k` the test is exactly **the sign of the Haar coefficient at scale
  `k`**, `H_k(L) := Σ_{i∈[L+k,L+2k)} z_i − Σ_{i∈[L,L+k)} z_i > 0`. The rule is a
  family of dyadic contrasts, as the task statement says; the "correction for
  odd `m`" is only the `(k, k+1)` weighting `kQ > (k+1)P`.
* **Shift invariance (5.1a).** `b2` is unchanged by `a ↦ a+c`, so `delta` is a
  function of `z` alone, and `A_n` is the set of balanced ternary vectors with
  `Σz = −ε`. Verified: 15,231 windows, 0 mismatches.

Verified: 200,000 random windows, 0 mismatches in all three forms (block-mean,
even Haar, odd `kQ > (k+1)P`).

### 1.3 Lemma 5.2 (the walk is a single-anchor sweep; proved)

Let `b1 ≡ 0`. The move table is `b2=1 → (dy,dx) = (0,−1)`, `b2=0 → (−1,0)`, and a
move that leaves the region triggers a cut. Since `ell = L` at the start of every
epochs and `dy = −1` immediately violates `ell+dy ≥ L`, the walk has this exact
form:

```
epoch L:  M = R−L+1
          m*(L) := max{ m ∈ [2,M] : b2(L,m) = 0 }        (∅ allowed)
          if m* = ∅      : halt, delta = L/(n−1)
          elif m* = M    : (L,R) ← (L+1, R−1)            [double cut]
          else           : (L,R) ← (L+1, R)              [left cut]
          terminate: M ≤ 1 → delta = L/(n−1);  M = 2 → ⊥
```

So **the anchor is constant inside an epoch and increases by exactly 1 between
epochs**; the epoch's probes are the chain `m = M, M−1, …, m*(L)` of prefixes of
one window, and `ℓ(L) := M − m*(L) + 1` probes are made (`M−1` in the halt case).

*Proof.* Read off the move table: `dx = −1` while `b2 = 1`, so `r` descends by 1
per probe with `ell = L` fixed; the first `b2 = 0` gives `dy = −1` and
`ell−1 = L−1 < L`, which is a cut, and `m* = M` exactly when the first probe is
already negative (a double cut, `ell == L and r == R`); otherwise a left cut.
`ell` never moves inside an epoch. ∎

Verified tape-for-tape against `walk_b0` (which is `_arr_rules.walk` with `b1≡0`)
on **all** ensemble tapes for `n = 6,…,12` — 111,322 tapes, 111,322 agreements
— and on `n = 13` (201,643 tapes, all agreeing), with anchors verified strictly
increasing on every tape and every halt condition re-checked directly.

### 1.4 Lemma 5.3 (anchor distinctness; proved)

Every contrast query `(L,m)` occurs **at most once per run**. Consequently no
contrast *answer* is ever reused, and answer-caching is worthless: the only
object that is reused across epochs is a prefix sum `Z_p`.

*Proof.* Anchors are distinct by Lemma 5.2 and each `m` is probed once per
epoch. ∎

This is the first place the dyadic form pays: in the cut formulation the machine
is free to re-read anything, and route (d) had to charge movement; here the
*questions* are provably non-repeating, so all reuse is in the *data*, not in the
answers (compare route 4 §3, where "cache the answers" was a live escape).

### 1.5 Lemma 5.4 (one left-to-right pass per epoch; proved)

`m*(L)` is computable with three counters and one monotone pass:

```
T = 0; P = 0; half = 0; best = 0
for m = 1..M:
    T += z[L+m−1]                       # T = Σ_{[L,L+m)} z
    while half < m//2: P += z[L+half]; half += 1     # P = Σ_{[L,L+half)} z
    if m ≥ 2 and (T−P)·(m//2) ≤ P·(m − m//2):  best = m
return best                              # 0 ⇒ no negative contrast ⇒ halt
```

Invariants `T = Z_{L+m} − Z_L`, `P = Z_{L+⌊m/2⌋} − Z_L` hold because both
prefixes advance monotonically as `m` grows; each `m` costs `O(1)` reads, so the
epoch costs `O(M)` reads and `O(log n)` bits. Correctness is Lemma 5.1.
Verified: 0 mismatches for every epoch of every tape at `n ≤ 12` (and `n = 13`).

So the whole of `delta_0` is "one cheap pass per epoch", `Θ(n)` epochs, `Θ(n²)`
reads — the mechanism of the `0.05–0.07 n²` step counts measured in
`LOWER_BOUND_ATTEMPT.md` §7. **Every single-epoch question is free; all the cost
is in the iteration.** This one sentence is the reason routes 5a–5c fail (§2–§4)
and the reason the surviving theorem (§5) has to be about memory.

---

## 2. Route 5a: counting dyadic-contrast vectors — closed, with an exponential ceiling

The proposal: the `2^s` memory states of a machine induce `2^s` classes of tapes;
count the distinct "dyadic contrast sign vectors" and conclude that a machine
with `s` bits can distinguish only `2^s` of them, so some class is large, and a
large class "needs many reads to identify".

The second half of that sentence is the false step, and it is false in the
direction opposite to the one needed. Two theorems.

### 2.1 Ceiling 5.8 (the counting dichotomy; proved)

Every counting argument in the read-only model concludes one of

```
(α)  T ≥ H(delta)/log₂3   ≤ log₂(n+1)/log₂3        ("the state must determine the output")
(β)  T ≥ c(a)             ≤ n−1                    ("the read set must certify the tape")
```

and both are `≤ n−1`.

*Proof.* (α): the final state determines `delta`, the state has `2^s` values, and
each read has 3 outcomes, so `T log₂3 ≥ H(delta)`. `delta` has `n+1` values.
(β): Lemma 1 of `LOWER_BOUND_ATTEMPT.md` (read set is a certificate) and Lemma 2
(every `(n−1)`-set certifies). There is no third form: any invariant that is a
function of the memory content alone is either a statement about which output is
forced (α) or about which tapes remain (β), and (β) saturates at `n−1` by
Lemma 2. ∎

So the counting family cannot exceed `n−1`, and (α) is *achieved* at
`s = log₂(n+1)`, which is exactly the paper's `O(log n)` budget.

### 2.2 Ceiling 5.9 (sign-vector counting is exponentially wrong; proved + measured)

Let `Φ(z)` be the vector of `b2` signs at anchor 0, `m = 2..n` (this is the
dyadic contrast vector of the first epoch). Measured over the full ensemble:

| n | tapes | #distinct `Φ` | log₂#Φ / n | #distinct `delta` | mean #Φ per `delta` | H(`delta`) | H(`Φ`) |
|---|---|---|---|---|---|---|---|
| 6 | 141 | 20 | 0.72 | 5 | 7.6 | 1.56 | 4.32 |
| 8 | 1,107 | 68 | 0.76 | 7 | 21.4 | 1.93 | 6.09 |
| 10 | 8,953 | 213 | 0.77 | 9 | 67.1 | 2.05 | 7.73 |
| 12 | 73,789 | 707 | 0.79 | 11 | 204.5 | 2.17 | 9.47 |

`#Φ = 2^{(0.72…0.79)n}`. So the naive counting argument gives

```
s ≥ log₂ #Φ = 0.79 n  bits at n = 12   (and ≈0.8n in general),
```

which is **false**: `s = O(log n)` bits suffice (`LOWER_BOUND_ROUTE4.md` Prop 4.1,
`O(n⁴)` reads). The machine is required to distinguish `delta`, of which there
are `n+1`; the map `Φ ↦ delta` is exponentially lossy — 204.5 signatures share
each `delta` at `n = 12`, and the entropy drops from 9.47 to 2.17 bits. So the
counting over-counts by `2^{Θ(n)}`, and no repair helps: any refinement that
counts a *finer* dyadic object than `delta` over-counts, and any coarser one is
`delta` itself, i.e. case (α) of Ceiling 5.8.

(The escape hatch "count the contrast vectors that a machine must *reproduce*"
was already closed by route 4 §3, Prop 3.2: pigeonhole produces a **large** class,
and a large class forces a **small** read set — `r ≤ 0.95 log₂ n` at
`s = log₂ n`. Counting always certifies easy tapes, never hard ones.)

---

## 3. Route 5b: round elimination on the halving tree — closed, three obstructions

The proposal: the contrast at scale `k` "depends on" the contrast at scale `2k`,
so recurse on the halving tree `m → ⌊m/2⌋`, eliminate one round at a time, and
sum the per-round costs.

### 3.1 Obstruction 5.10a (the tree is a comb, not a tree; proved)

For a fixed anchor `L` the tested intervals are `[L, L+⌊m/2⌋)` and `[L, L+m)` for
`m = 2..M` — all prefixes of one window, i.e. a **chain**, not a branching tree.
The halving relation `m → ⌊m/2⌋` links scales but the two intervals tested at
scale `m` are not the intervals tested at scale `⌊m/2⌋`; there is no
containment-respecting tree of sub-instances. Across anchors one gets `n` chains
(a "comb") sharing cells but not anchored sums. Round elimination needs a
*shrinking* instance; a chain of prefixes has no sub-instance smaller than the
window, because every prefix starts at the same anchor.

*Proof.* Immediate from the definitions: all intervals in an epoch have left
endpoint `L`. ∎

### 3.2 Obstruction 5.10b (no geometric shrinkage; proved)

After one epoch the region is `[L+1, R]` or `[L+1, R−1]`, so `M` decreases by
exactly 1 or 2 — additively, not multiplicatively. Any round-elimination scheme
therefore needs `Θ(n)` rounds (measured `E/n = 0.37 → 0.65`, §5.1), whereas the
technique works when the instance halves each round so that `log n` rounds
suffice. `Θ(n)` rounds × `Ω(1)` per round is only `Ω(n)`.

### 3.3 Obstruction 5.11 (round additivity fails; proved + measured)

The tempting per-round cost is "the machine must learn `m*(L)`, a
`log₂ n`-bit number, and a read carries `log₂3` bits, so `Ω(log n)` reads per
epoch, hence `Ω(n log n)` total". This is the exact place where the dyadic route
appears to succeed and does not. The sum is invalid because the epoch answers are
functions of the *same* `n` cells.

*Validity.* The only sound form is the joint one:
`T ≥ H(A_1,…,A_E)/log₂3`, and `H(A_1…A_E) ≤ H(z) = log₂|A_n| = n log₂3 − O(log n)`,
hence `T ≥ n − O(log n)`: capped at the certificate bound, no gain.

*And the invalid form is quantitatively weak anyway.* Measured (full ensemble,
`n ≤ 12`), with `A_e = m*(L_e)`:

| n | k | Σ_{e≤k} H(A_e) | H(A_1..A_k) | ratio | "naive reads" = ΣH/log₂3 | cap n |
|---|---|---|---|---|---|---|
| 8 | 1 | 1.838 | 1.838 | 1.000 | 1.16 | 8 |
| 8 | 2 | 3.552 | 2.959 | 0.833 | 2.24 | 8 |
| 8 | 4 | 8.867 | 6.262 | 0.706 | 5.59 | 8 |
| 12 | 1 | 2.160 | 2.160 | 1.000 | 1.36 | 12 |
| 12 | 2 | 4.352 | 3.604 | 0.828 | 2.75 | 12 |
| 12 | 3 | 6.671 | 4.887 | 0.733 | 4.21 | 12 |
| 12 | 4 | 8.926 | 5.899 | 0.661 | 5.63 | 12 |
| 12 | 6 | 16.984 | 10.452 | 0.615 | 10.72 | 12 |
| 12 | ALL (10) | 16.599 | — | — | 10.47 | 12 |

Two independent weaknesses:

1. **Correlation.** `H(joint)/Σ H(marginal)` falls to `0.62` by six rounds — the
   answers are strongly dependent, as they must be, being functions of one tape.
2. **Low entropy per answer.** `H(A_e) ≈ log₂ log₂ n`, not `log₂ n`: measured
   `H(m*/M) = 1.84, 2.02, 2.16, 2.89, 3.05` bits at `n = 8, 10, 12, 32, 128`
   (20-bin discretisation, so these are lower bounds on `H(m*)` but the growth
   rate is the point). The reason is visible in the same measurement: `m*/M` is
   concentrated near 1 — `mean m*/M = 0.85`, `P(m* = M) = 0.58…0.66` — because a
   double cut is by far the most likely outcome of an epoch
   (`P(m* ≤ M/4) = 0.03…0.07`).

So summing marginals would give `O(n log log n)` reads at best — superlinear, and
still short of `n log n` — and the sum is not a lower bound in any case. **The
dyadic route's information form is invalid, and its invalid form is not even
strong enough to reach the target.**

---

## 4. Route 5c: reduction from a known hard problem — closed, the core is one-pass trivial

The proposal: reduce from "majority of dyadic block means", or element
distinctness on a dyadic hierarchy, or the max-descent problem on the halving
tree.

**Theorem 5.12 (no non-circular reduction).** Every single-anchor or single-scale
dyadic statistic of the tape is computable in **one left-to-right pass with
`O(1)` counters** (hence `O(n)` reads, `O(log n)` bits):

| problem | one-pass algorithm | reads | bits |
|---|---|---|---|
| `m*(L)` — max negative contrast at one anchor (Lemma 5.4) | `T, P, best` | `O(M)` | `O(log n)` |
| is `Z` concave on `[L,R]` | second differences, 2 counters | `O(M)` | `O(1)` |
| majority of the `k`-block-mean comparisons at one scale | running block sums | `O(n)` | `O(log n)` |
| argmax / argmin of `Z`, first/last 2, longest 2-run | 1–2 counters | `O(n)` | `O(log n)` |
| all `log n` scales simultaneously | one counter pair per scale | `O(n log n)` | `O(log² n)` |

All of these are in the *cheap* regime (linear reads, logarithmic bits), so none
of them can be the hard core of a superlinear lower bound. Any problem that
encodes the *iteration* over anchors is `delta` itself up to the `b1` correction,
so the reduction is circular. ∎

In particular the reduction cannot even be *stated* without already assuming what
is to be proved: the hardness of `delta_0` is not located in any dyadic query
(§1.5) but in the sequence of `Θ(n)` anchors.

---

## 5. What survives: the checkpoint theorem

The one place where the dyadic structure yields a superlinear bound is the
accounting of *what the machine must know* at each anchor.

### 5.1 Model R (exact interval-sum memory)

A machine is in model R with parameter `r` if its memory, at every moment, is

* a family of at most `r` pairs `(I_j, σ_j)` with `I_j ⊆ [0,n)` an interval and
  `σ_j = Σ_{i∈I_j} z_i` the exact sum (each `I_j` having been read completely at
  some earlier time),
* plus `O(1)` indices and `O(1)` integers, and nothing else about `z`.

Its next read and its output are functions of that memory. Equivalently: the
memory is a graph on the vertices `{0,…,n}` with at most `r` edges `[p,q)` each
labelled `Z_q − Z_p`, so `Z` is known up to one additive constant per component,
and at most `2r` vertices are endpoints of the stored intervals.

The model also contains a second, separate restriction, stated here because it is
part of the theorem and not a detail: **the machine must *determine* each `b2`
from its memory before acting on it** — its decision must be constant over the
tapes consistent with the memory. A machine that instead acts on a guess can be
right by luck on a given tape, and converting a lucky machine into a determined
one is exactly the indistinguishability problem at which route (c) of
`LOWER_BOUND_ATTEMPT.md` stopped at `n−1`. So model R has *two* restrictions, and
§5.4 names both.

A general `s`-bit memory may hold *arbitrary* functions of `z` (an argmax, a hull
vertex, a run length, a cached bit) and need not decide exactly at each step; it is
not in model R.

### 5.2 Theorem 5.5 (checkpoint lower bound; proved)

For every tape and every machine in model R,

```
T  ≥  Σ_e  max(0, ℓ_e − 2r),          ℓ_e = M_e − m*_e + 1  (or M_e − 1 if the epoch halts).
```

*Proof.* In epoch `e` the machine must decide `b2(L_e, m)` for **every**
`m ∈ [m*_e, M_e]`: it must know that `m*_e` is negative and that each larger `m`
is positive, since otherwise it cannot know that `m*_e` is the maximum. By
Lemma 5.1 each such decision is a comparison of `Z` at the three vertices
`L_e`, `L_e+k`, `L_e+m`, so by the exactness requirement of model R it is a
function of the stored interval sums, hence determined only if `L_e` and `L_e+m`
lie in the same component of the stored-interval graph. If they do not, then
perturbing `z` on a `±1` pair inside the gap between the two components
(`a_i: 0→1` together with `a_j: 2→1` keeps `#0 − #2`, hence `Σz = −ε`, and leaves
every stored interval sum — hence the entire memory, hence the decision —
unchanged) flips `b2` for that `m`, contradicting exactness. Therefore all `ℓ_e`
consecutive vertices `{L_e + m : m*_e ≤ m ≤ M_e}` must be in the component of
`L_e` when the epoch ends.

Bringing a new vertex into a component costs at least one read: reading cell `i`
yields `z_i`, which adds exactly one edge `[i, i+1)`, so one read can enlarge a
component's vertex set by at most one. At the start of the epoch at most `2r`
vertices are endpoints of stored intervals at all. Hence epoch `e` costs at least
`ℓ_e − 2r` reads (trivially if negative). Reads are partitioned by epoch, so the
per-epoch bounds add. ∎

The two restrictions are both used: additivity at "hence the entire memory", and
exactness at "contradicting exactness". Neither can be dropped from the proof,
and §5.4 is about that.

### 5.3 Corollaries (proved in model R; constants measured)

`r ≤ s / log₂(2n+1)` because each stored sum is an integer in `[−n,n]`.

**C5.6.** With `s = O(log n)` we get `r = O(1)`, so `T ≥ Σ_e ℓ_e − O(n)`.
Measured `Σ_e ℓ_e / n²` over the full ensemble (`n ≤ 12`) and 60–200 sampled
tapes (`n ≥ 14`, seed 2026):

| n | tapes | E/n | Σℓ/n² | mean ℓ/n | max ℓ/n | r=0 | r=1 | r=2 | r=4 |
|---|---|---|---|---|---|---|---|---|---|
| 6 | 141 | 0.374 | 0.1168 | 0.313 | 0.833 | 0.1168 | 0.0288 | 0.0032 | 0.0000 |
| 8 | 1,107 | 0.414 | 0.1121 | 0.271 | 0.875 | 0.1121 | 0.0387 | 0.0111 | 0.0000 |
| 10 | 8,953 | 0.447 | 0.1077 | 0.241 | 0.900 | 0.1077 | 0.0443 | 0.0182 | 0.0005 |
| 12 | 73,789 | 0.471 | 0.1035 | 0.220 | 0.917 | 0.1035 | 0.0479 | 0.0236 | 0.0023 |
| 16 | 200 | 0.500 | 0.0911 | 0.182 | 0.938 | 0.0911 | 0.0472 | 0.0279 | 0.0070 |
| 24 | 200 | 0.545 | 0.0912 | 0.167 | 0.958 | 0.0912 | 0.0581 | 0.0408 | 0.0169 |
| 32 | 200 | 0.556 | 0.0797 | 0.143 | 0.969 | 0.0797 | 0.0553 | 0.0432 | 0.0243 |
| 64 | 200 | 0.603 | 0.0756 | 0.125 | 0.969 | 0.0756 | 0.0622 | 0.0550 | 0.0425 |
| 128 | 60 | 0.604 | 0.0565 | 0.093 | 0.875 | 0.0565 | 0.0499 | 0.0465 | 0.0404 |
| 256 | 60 | 0.647 | 0.0696 | 0.108 | 0.816 | 0.0696 | 0.0658 | 0.0635 | 0.0592 |

(last four columns: `Σ_e max(0, ℓ_e − 2r)/n²`, i.e. the theorem's bound normalised.)

So `Σ_e ℓ_e ∈ [0.056, 0.117]·n²` throughout. Writing `E ≤ 0.65 n` (measured), the
theorem reads `T ≥ Σℓ_e − 2rE ≥ (0.056 − 1.3 r/n)·n²`, so the bound stays
`Ω(n²)` for `r ≤ 0.02 n` (giving `≥ 0.03 n²`), and in general `T = Ω(n²)` for
`r = o(n)`. In particular:

```
s = O(log n)  ⟹  T = Ω(n²)   in model R,
```

which is *stronger* than the conjecture (`n² ≫ n polylog`). It also beats route
(d)'s pointer bound `0.396 n²/s` — at `s = log₂ n` by the factor
`0.056·log₂n/0.396 = 0.14 log₂ n`, i.e. by `1.4×` at `n = 1024` and `2.8×` at
`n ≈ 10⁶` — because the epoch count is a sharper accounting than the per-cut
communication (one read can cross all `n` cuts, which is what route (d) had to
pay for; one read can serve only one new vertex here).

**C5.7.** A `p`-pass machine with `r` checkpoint counters has `T ≤ pn`, so
`p ≥ (Σℓ_e − 2rE)/n = Ω(n)` at `r = O(1)`: **the walk needs `Ω(n)` passes.** A
one-pass or `O(1)`-pass streaming algorithm cannot compute `delta_0`, which is a
multi-pass streaming lower bound of the strongest kind (linear in the number of
passes).

### 5.4 Why the theorem cannot be lifted — the model boundary

`C5.6` is a theorem about model R, and model R carries **two** restrictions that
the read model does not: additivity of the memory, and exactness of each decision.
Neither is a technical convenience, and here is exactly what each excludes.

*Exactness.* An `s`-bit machine may act on a value its memory does not fix and be
right by luck on a given tape. Restoring the theorem without exactness requires
an indistinguishability argument: two tapes consistent with the memory must be
shown to have different `delta`, not merely different `b2` at one point. That is
the adversary problem, and `LOWER_BOUND_ATTEMPT.md` §5 showed it is defeated after
`n−1` distinct cells. So dropping exactness immediately returns the bound to
`n−1`.

*Additivity.* An `s`-bit memory can hold

* `Θ(s)` **answers** `b2(L,m)` for chosen pairs — but by Lemma 5.3 no answer is
  ever reused, so this is worthless;
* an **argmax / argmin** of `Z` on a block (`O(log n)` bits, not an interval sum);
* a **hull vertex**, a run length, a monotone-chain summary;
* any other `s`-bit function of the tape.

None of these is expressible as `O(1)` interval sums, and any of them could in
principle replace a re-scan. The known machine
(`LOWER_BOUND_ROUTE4.md` Prop 4.1: `s = O(log n)` bits, `O(n⁴)` reads) uses no
summary at all — it re-reads — so it sits *inside* model R with `r = O(1)` and
pays `Ω(n²)`, consistent with `C5.6`. But `C5.6` therefore says nothing about a
machine that carries a summary, and no technique in the deterministic read-only
toolbox forces a memory to be additive.

**C5.13 (the necessary condition on a counterexample; proved).** Contrapositive
of `C5.6`: any machine with `s = O(log n)` and `T = o(n²)` must, at some epoch,
hold a memory state that is **not** equivalent to `O(1)` interval sums — i.e. it
must carry a *non-additive* summary of the tape. This is the exact analogue of
route (d)'s Corollary D6 ("a counterexample must teleport, `J̄ = Ω(n/log³n)`)",
and it is the sharpest form of the route-(e) diagnosis:

```
route (d):  a counterexample must move far     (break the movement accounting);
route (e):  a counterexample must remember non-additively (break the sum accounting).
```

Both are necessary conditions, neither is sufficient, and both point at the same
missing lemma: a lower bound for arbitrary `s`-bit state, i.e. a branching-program
bound at width `2^s = n^{O(1)}`.

---

## 6. The ceilings, in the style of route 4 §2.5

Four statements, each a theorem about a technique rather than about `delta`.
Together they say that route (e) is closed at the model boundary, not at a
numerical bound.

**Ceiling 5.8 (counting).** Every counting/knowledge argument concludes
`T ≥ H(delta)/log₂3 ≤ log₂(n+1)/log₂3` or `T ≥ c(a) ≤ n−1`. Both `≤ n−1`. ∎

**Ceiling 5.9 (sign vectors).** Counting dyadic contrast vectors gives
`s ≥ 0.79 n` bits at `n = 12`, which is false (`s = Θ(log n)` suffices). The
over-count is `2^{Θ(n)}`, and the only lossless alternative is to count `delta`
itself, i.e. Ceiling 5.8. ∎

**Ceiling 5.10 (movement).** Every dyadic argument of the form "the machine must
visit these `ℓ_e` positions" concludes `TV ≥ Σ_e ℓ_e ∈ [0.056,0.117]n²`; and by
route (d)'s identity `TV = T·J̄` with `J̄ ≤ n−1`, that gives
`T ≥ Σℓ_e/(n−1) = O(n)`. So movement-style dyadic arguments cap at `O(n)` —
quantitatively below the certificate bound. ∎

**Ceiling 5.11 (information).** The valid information form is
`T ≥ H(A_1..A_E)/log₂3 ≤ n log₂3/log₂3 = n`. The invalid per-round form is
`O(n log log n)` (measured `H(A_e) ≈ log₂log₂n`) and wrong. So the information
form of route (e), made valid, is exactly the certificate bound. ∎

**Ceiling 5.12 (model).** `C5.6` holds only for additive memory. The restriction
is not implied by the read model, and no known argument forces it; consequently
the only superlinear statement route (e) produces is conditional, and the
condition is precisely the open problem. ∎

Note the parallel with route (d), which is the strongest evidence that the two
routes have hit the same wall by different roads:

```
route (d):  T·J̄·s ≥ Σ_k CC_k − ns,   Σ_k CC_k ≤ 0.396 n²   ⟹  ceiling 0.396n/s;
route (e):  T      ≥ Σ_e ℓ_e − 2rE,   Σ_e ℓ_e ≤ 0.117 n²    ⟹  ceiling 0.117 n²   (model R).
```

Route (d) divides by the jump (`≤ n`), route (e) divides by nothing but must pay
`2r` per epoch; both are "conservation laws" whose right-hand side is bounded by
the tape's content (`0.396n²` bits of cut communication, `0.117n²` contrast
probes), and both are quadratic in `n` with no route to the read count except
through a model of access.

---

## 7. What a successful proof would have to do

Stated as sharply as the four closed routes allow.

1. **It must prove redundancy, not information.** The tape has `n log₂3 = 1.58n`
   bits and `n−1` cells certify `delta` (Lemma 2), so no accounting of *what the
   machine learns* can exceed `n` reads (Ceilings 5.8, 5.11). The target
   `n·polylog n` is larger than the tape's own information by the factor
   `polylog n`; therefore the proof must exhibit cells re-read
   `ω(polylog n)` times. Equivalently `s = ω(log T)`.

2. **It must therefore control the correlation between epochs.** By §3.3 the
   per-round information sum is `O(n log log n)` and the joint information is
   `≤ n`; the gap between them *is* the target. Concretely, a proof must show
   that for every `s = O(log n)`-bit machine and every `k`, the machine's state
   after `k` epochs does not capture the dependence between `A_1..A_k` and
   `A_{k+1}` — i.e. that the machine cannot "skip" an epoch by predicting it.

3. **It must not be a statement about a restricted access pattern.** Route (d)
   needed small `J̄`; route (e) needs additive memory. Both restrictions are
   *sufficient* for a superlinear bound and neither is necessary. By Lemma 5.3
   the third natural restriction (cached answers) is empty. So the restriction
   that remains available is the real one: arbitrary `s`-bit state, i.e. a
   branching program of width `2^s` and length `T` (Lemma 5 of
   `LOWER_BOUND_ATTEMPT.md`) at width `n^{O(1)}`, length `n polylog n`. For that
   model no superlinear lower bound is known for any explicit function, and the
   specific gap — additive vs arbitrary memory — is not addressed by any
   technique in the deterministic read-only toolbox.

4. **Concretely, the smallest open lemma that would finish route (e)** is:

   > **(⋆)** For every `s = O(log n)`-bit machine `M` that computes `delta_0`
   > and every tape, the memory of `M` at the start of some epoch determines
   > fewer than `ℓ_e − o(ℓ_e)` of the `ℓ_e` contrast signs at that anchor,
   > *unless* `M` made `ω(polylog n)` repeat reads.

   (⋆) is exactly `C5.6` without the additivity hypothesis. It is implied by,
   and probably equivalent to, the branching-program conjecture; it is *not*
   implied by any counting, communication, movement or entropy statement, since
   all four cap at `n` (`LOWER_BOUND_ATTEMPT.md` Lemma 6, route 4 §2.5,
   Ceilings 5.8–5.11).

5. **A cheap falsifier for proposed counterexamples**, in the spirit of route (d)'s
   `J̄` statistic: any candidate `n polylog n` simulator must be checked against
   `Σ_e (ℓ_e − 2r)`. If its memory is additive with `r = o(n)` and it beats
   `0.05 n²`, it is wrong (Theorem 5.5) — whatever else it does. Route (e)
   therefore makes the conjecture falsifiable by a second statistic, independent
   of route (d)'s mean jump.

---

## 8. Proved / measured / conjectured

**Proved unconditionally (in the model of `LOWER_BOUND_ATTEMPT.md` §1).**

1. **L5.1** block-mean / Haar normal form, both parities; **L5.1a** shift
   invariance (`delta` is a function of `z = a−1` on the balanced ternary cube).
2. **L5.2** epoch form: single anchor per epoch, anchors strictly increasing,
   `m*` = largest negative contrast, cuts as tabulated.
3. **L5.3** anchor distinctness: no contrast query, hence no contrast answer, is
   ever repeated in a run.
4. **L5.4** the epoch answer is a one-pass, three-counter function of the window.
5. **T5.5** `T ≥ Σ_e max(0, ℓ_e − 2r)` for machines with `r` exact interval sums
   that determine each decision (model R, both restrictions named in §5.1);
   **C5.6** `s = O(log n) ⟹ T = Ω(n²)` in that model; **C5.7** `p`-pass machines
   need `p = Ω(n)`; **C5.13** `s = O(log n)` and `T = o(n²)` force a memory state
   that is not additive.
6. **C5.8** counting dichotomy (`≤ log₂(n+1)/log₂3` or `≤ n−1`);
   **C5.9** sign-vector counting over-counts by `2^{Θ(n)}`;
   **C5.10** movement-form dyadic arguments cap at `O(n)`;
   **C5.11** the valid information form caps at `n`.
7. **T5.12** no non-circular reduction: every single-anchor / single-scale dyadic
   statistic is one-pass with `O(1)` counters.
8. **Obstructions 5.10a–b**: the query family is a comb of chains (no shrinking
   sub-instance) and the region shrinks additively (1–2 cells per epoch), so round
   elimination needs `Θ(n)` rounds for `Ω(1)` each.

**Measured (this run; exhaustive over `A_n` for `n ≤ 12`, 60–200 sampled tapes,
seed 2026, above that).**

9. Normal form: 200,000 random windows, 0 mismatches (all three forms);
   15,231 shift-invariance checks, 0 mismatches.
10. Epoch form vs `walk_b0`: 111,322/111,322 tapes at `n = 6…12`, and
    201,643/201,643 at `n = 13`; anchors strictly increasing on every tape;
    one-pass epoch solver mismatches 0; halt condition re-verified 0 failures.
11. `Σ_e ℓ_e/n² = 0.1168, 0.1121, 0.1077, 0.1035, 0.0911, 0.0912, 0.0797,
    0.0756, 0.0565, 0.0696` at `n = 6, 8, 10, 12, 16, 24, 32, 64, 128, 256`;
    `E/n = 0.374…0.647`; `mean ℓ_e/n = 0.313…0.108`; `max ℓ_e/n ≤ 0.97`.
12. Model-R columns of the §5.3 table; `T ≥ Σℓ_e − 2rE ≥ (0.056 − 1.3r/n)n²`,
    so `Ω(n²)` for `r ≤ 0.02n`.
13. Lossiness: `#Φ = 20, 68, 213, 707` and `#delta = 5, 7, 9, 11` at
    `n = 6, 8, 10, 12`; `204.5` signatures per `delta` at `n = 12`;
    `H(Φ) = 9.47` vs `H(delta) = 2.17` bits.
14. Round additivity: `H(joint)/ΣH(marginal) = 1.000, 0.833, 0.719, 0.706` at
    `n = 8` and `1.000, 0.828, 0.733, 0.661, 0.615` at `n = 12` for
    `k = 1, 2, 3, 4, 6`; full sum `16.599` bits `= 10.47` "reads" at `n = 12`
    against the cap `n = 12`.
15. `H(m*/M) = 1.838, 2.021, 2.160, 2.892, 3.053` bits at `n = 8, 10, 12, 32,
    128`; `mean m*/M = 0.85`, `P(m* = M) = 0.58…0.66`, `P(m* ≤ M/4) = 0.03…0.07`.
16. Contrast gaps: fraction of probes with `m ≥ M/2` `= 0.935, 0.916, 0.903,
    0.878, 0.872`; fraction with `|gap| ≤ √k` `= 0.596, 0.593, 0.587, 0.589,
    0.564`; `mean |gap|/k = 0.656, 0.615, 0.581, 0.352, 0.184` at
    `n = 8, 10, 12, 32, 128`.

**Conjectured.**

17. The paper's conjecture: `delta` has no `n^{O(1)}`-width, `n·polylog n`-length
    branching program, equivalently no `s = O(log n)` machine with `n polylog n`
    reads. Route (e) does not prove it and, by Ceilings 5.8–5.12, cannot: the
    additive-memory theorem (`C5.6`) is stronger than the conjecture but confined
    to model R, and no step from `s` bits to `r` interval sums is available.
18. **Open lemma (⋆)** of §7.4, the exact residue of route (e). It is the
    "non-additive memory is useless" statement; it is the same object as the
    branching-program conjecture, and it is not reachable from any of the four
    ceilings.

---

## 9. Files

```
probe_ologn/lb_route5.py    normal form, epoch decomposition, one-pass epoch
                            solver, model-R accounting, lossiness, round
                            additivity, gap and m*/M distributions
probe_ologn/lb_route5.log   its full output (the numbers quoted above)
```

Run with `$P = C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`,
in `probe_ologn/`:

```
$P lb_route5.py verify  6 7 8 9 10 11 12     # V1-V4, ~30 s
$P lb_route5.py measure 6 8 10 12 16 24 32 64 128 256
$P lb_route5.py hist    8 10 12 32 128       # m*/M distribution
$P lb_route5.py gaps    8 10 12 32 128       # contrast gap distribution
$P lb_route5.py lossy   6 8 10 12            # delta vs first-anchor signatures
$P lb_route5.py rounds  8 10 12              # round additivity
```

`lb_route5.log` is the concatenation of all six.

---

## 10. Consequences for §8 of `bthm_merged.tex` (if route (e) is cited)

1. The rule has an exact normal form worth stating independently of any lower
   bound: `b2` is the sign of a dyadic (Haar) contrast on the prefix sums of
   `z = a−1`, `delta` is shift-invariant, and the `b1≡0` walk is a
   single-anchor sweep whose epoch answer is a **one-pass, three-counter**
   function of the window. This is the cleanest available description of *why*
   the trajectory is long: `Θ(n)` anchors, each needing a pass.
2. The honest ceiling list becomes five: `n−1` (certificate/adversary/
   transcript), `0.79n/log₂n` (one cut), `log₂n` (state counting),
   `0.396n/log₂n` (all cuts + jump), and now `n` (counting), `0.79n` bits
   (sign vectors — false), `O(n)` (movement), `n` (information). Only the
   last route additionally gives a **quadratic** bound, and only for additive
   memory.
3. The two necessary conditions on a counterexample should be stated together:
   it must **teleport** (mean jump `Ω(n/log³n)`, route d) and it must remember
   **non-additively** (route e). Both are cheap to test and both are provable.
4. Do not describe the dyadic structure as "the route that almost works". It
   gives a superlinear bound only inside a model whose restriction is exactly the
   open problem; the strongest unconditional content is the normal form and the
   four ceilings.
