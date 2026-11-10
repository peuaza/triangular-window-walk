# Route (d): aggregate route (a) over *every* cut, and charge the jumps

Fourth attempt at the O(log n)-bit / n·polylog(n)-read separation for the walk
rules, after `LOWER_BOUND_ATTEMPT.md` closed routes (a), (b), (c). Everything
below is about the *same* function `delta` on the same ensemble `A_n` as that
file; §1 fixes notation. New probes: `lb_cuts.py` (output `lb_cuts.log`),
`lb_cyc.py`, `lb_rank14.py`.

---

## 0. Verdict

**Route (d) exists, is not one of the three, produces a superlinear lower bound
for this function in a movement-charged read model, and then dies at a provable
ceiling that is *below* the certificate bound.** What it adds is not a
proof of the conjecture; it is a theorem plus a sharp diagnosis of what the
conjecture has to prove instead.

| item | statement | status |
|---|---|---|
| Theorem D3 | `T·J̄·s ≥ Σ_{k=1}^{n−1} CC_k(delta) − n·s`, `J̄` = mean jump length of the read sequence | **proved** |
| Cor. D4 | `Σ_k CC_k ≤ log2(3)·n²/4 = 0.396 n²`, so `T·J̄ ≥ 0.396n²/s − n`; with `J̄ ≤ n−1` this is `T ≥ 0.396n/log2 n − 1` | **proved** |
| Cor. D5 | `J̄ = 1` (single-head pointer machine): `T ≥ 0.396n²/s − n`, i.e. `0.396 n²/log2 n` at `s = log2 n` | **proved** |
| Cor. D6 | any machine with `s = log2 n + O(1)` that reads `≤ n·polylog(n)` cells must have mean jump `J̄ ≥ 0.396 n /(log2 n · polylog n)` | **proved** |
| ceiling | route (d) cannot beat `0.396n/log2 n` in the random-access model, and every member of the family is of the form `(bound on T·J̄)` with `J̄ ≤ n` | **proved** (§2.5) |
| parent's idea (a) | configuration counting gives `s ≥ log2(n+1)`; its refinement has the sign backwards — pigeonhole yields a *large* class, and a large class needs a *small* read set (`r ≤ 0.95 log2 n`) | **proved** (§3) |
| parent's idea (c) | the space complexity of `delta` is exactly `Θ(log n)` bits, so no `ω(log n)` space lower bound exists | **proved (upper bd by construction)** (§4) |
| parent's idea (b) | no `O(log²n)`-bit one-pass simulation can exist (it would beat a measured `≈0.7n`-bit one-way bound) | **proved mod. measurement** (§5) |
| CC input | `CC_k` two-way at the middle cut: `≤ CC_k^{1way}` and `≥ log2 rank = 0.72·CC_k^{1way}` at `n = 13`; the best block rank is 0.83–0.98 of its width | **measured** (§2.8) |

Net: the whole *cut/communication* family — single cut (route a), all cuts
aggregated (§2), p-pass, pointer, external-memory — is closed in the
random-access model, and it is closed for a structural reason, not for lack of
effort: **the aggregation gains a factor `n` (there are `n` cuts) and the jump
charge loses exactly a factor `n` (one read can cross all `n` cuts).** The
conservation law is `T·J̄·s ≥ Σ_k CC_k − ns`; the conjecture is strictly stronger
than it by `≈ log² n`.

What survives is a *necessary condition* on any counterexample (§2.7): a machine
that beats `n·polylog(n)` with `O(log n)` bits must, on average, **jump
`Ω(n/log³n)` cells per read**. So the missing lemma is not "the machine must
repeat reads" (Lemma 7 of the earlier file) but "**the machine cannot teleport**"
— a lower bound on jump length, and equivalently a genuine branching-program
lower bound at width `n^{O(1)}`.

---

## 1. Notation

`A_n = {a ∈ {0,1,2}^n : #0(a) − #2(a) = ε}`, `ε = 0` (n even), `1` (n odd);
`delta(a) ∈ {0, 1/(n−1), …, 1} ∪ {⊥}` is the walk's meeting position, `⊥` when
it does not halt. A machine `M` has `s` bits of state, `T_M(a)` reads, and its
state determines the next position read and the halting/output decision.
Mean jump and total variation of the read sequence `σ_1,…,σ_T ∈ [0,n)`:

```
jump_t = |σ_{t+1} − σ_t| ,   J̄ = (1/T)Σ_t jump_t ,   TV = Σ_t jump_t = T·J̄ .
                                                                        (1.1)
```

For a cut `k`, `x_k(a) = #{t : σ_t < k ≤ σ_{t+1} or σ_{t+1} < k ≤ σ_t}` is the
number of crossings, and `CC_k = CC_k(delta)` is the deterministic two-way
communication complexity of `delta` on the ensemble across that cut (a promise
problem: the players' halves must sum to the ensemble's counts).

Recall from the earlier file: Lemma 2 (every `(n−1)`-set certifies), Lemma 4
(`T ≤ 2^{s+O(1)}`), and the three ceilings `n−1` (certificate / adversary /
transcript), `0.79n/log2 n` (one cut), `s ≥ log2 n` (state counting).

---

## 2. Route (d)

### 2.1 Lemma D1 (movement = sum of crossings; proved)

`TV = Σ_{k=1}^{n−1} x_k`.

*Proof.* A step from `σ_t` to `σ_{t+1}` crosses exactly the cuts strictly
between the two positions, i.e. `|σ_{t+1} − σ_t|` of them, and the cuts crossed
are exactly those counted by `x_k` at that step. Sum over `t`. ∎

Note the failure mode of route (a) in one line: **one read can cross all `n`
cuts**, so `x_k` can be incremented `n` at a time while `T` grows by 1. Lemma D1
is what converts per-cut statements into a statement about the *path*.

### 2.2 Lemma D2 (per-cut protocol; proved)

For every `k`, `(x_k(a) + 1)·s ≥ CC_k`.

*Proof.* Alice holds `a[0,k)`, Bob `a[k,n)`. Both simulate `M`; the player who
owns the cell to be read holds the state. Whenever the simulated head moves
across the cut, the state (≤ `s` bits) is sent. The number of messages is
`x_k(a)+1`, and the protocol computes `delta` on the ensemble. ∎

(No move is charged: this is exactly why route (a) is weak, and it is also why
it is *correct* — the protocol is a legitimate communication protocol.)

### 2.3 Theorem D3 (the tradeoff; proved)

For every correct machine and every tape,

```
TV(a)  ≥  Σ_{k=1}^{n−1} CC_k / s  −  (n−1) ,      equivalently
T(a) · J̄(a) · s  ≥  Σ_k CC_k  −  n·s .                                    (2.1)
```

*Proof.* Lemma D2 gives `x_k ≥ CC_k/s − 1` for each `k` (trivially if the right
side is negative). Sum over the `n−1` cuts and apply Lemma D1:
`TV = Σ_k x_k ≥ Σ_k CC_k/s − (n−1)`. Multiply by `s` and use `TV = T·J̄`. ∎

This is the new object. It is a three-way tradeoff — **reads × jump × space ≥
(total cut communication)** — and it contains route (a) as the special case in
which one cut is used and no movement is charged. Two independent factors are
being spent: `Σ_k` over `n` cuts (a gain of `n` over route (a)) and `1/J̄`
(a loss of up to `n`).

### 2.4 The only upper bound available on the right-hand side (proved)

`CC_k ≤ min(k, n−k)·log2 3` (one player sends his half), hence

```
Σ_{k=1}^{n−1} CC_k  ≤  log2(3) · ⌊n²/4⌋  =  0.3962 n²  (n even).           (2.2)
```

So (2.1) becomes, at `s = log2 n`,

```
T·J̄  ≥  0.396 n²/log2 n  −  n ,      and with J̄ ≤ n−1:   T ≥ 0.396n/log2 n − 1.  (2.3)
```

### 2.5 The ceiling of route (d), and of the family (proved)

Two statements, both immediate from (2.3) and Lemma D1, both worth separating
because they are the reason this route cannot finish the job.

**(i) Route (d) is dominated by the certificate bound.** At `s = log2 n`,
(2.3) gives `T ≥ 0.396n/log2 n`, which is smaller than `n−1` for every `n ≥ 4`
(by a factor ≈ 25 at `n = 1024`). It is not "slightly" short: it is a
`log n`-sized factor short of the bound the paper already has.

**(ii) The whole family is closed.** Every argument in the cut family concludes
something of the form `T·J̄ ≥ A`; to obtain a bound on `T` one must use
`J̄ ≤ n−1`, which costs a factor `n`. Since `Σ_k CC_k ≤ 0.396n²` by (2.2), every
such argument yields `T ≥ 0.396n²/(s·n) = 0.396n/s` at the very best — **an upper
bound on what the technique can ever prove**, not on what is true. Concretely:

* a single cut: `0.79n/s` (route (a));
* all cuts, unweighted: `0.396n/s` — *worse*, by 2×;
* all cuts beats one cut only if `J̄ < ½n`, i.e. only for machines that do not
  teleport at the scale of the tape. This is the precise sense in which route
  (d) is a strict improvement of route (a): it buys a factor `≈ 0.49n/J̄`.

The inequality (2.2) is where the ceiling comes from, and it cannot be improved:
`CC_k` is *at most* the length of the smaller half in bits, because that half can
simply be transmitted. Any cut-based proof of a superlinear bound would need
`Σ_k CC_k = ω(n²)`, i.e. more communication than the tape contains.

### 2.6 Corollary D5 (superlinear for pointer machines; proved)

A machine whose reads form a path with unit steps (`J̄ = 1` — a single head on a
two-way read-only tape, i.e. exactly the "pointer strategy" class the paper
contrasts with random access) satisfies

```
T  ≥  0.396 n²/s  −  n        = 0.396 n²/log2 n  at s = log2 n .           (2.4)
```

This is a genuine superlinear bound, and within this project the first for this
function in a movement-charged model. It is *tight up to `log n`*: the rule's own
two-header trajectory evaluates `Θ(n²)` window questions and moves a header
`Θ(n²)` times (`steps/n² = 0.048–0.068` in `LOWER_BOUND_ATTEMPT.md` §7), and its
reads move by `O(1)` cells apiece under the two-pointer trick of `REPORT.md`
§6.2, so it is an instance of (2.1) with `J̄ = O(1)` and `T = Θ(n²)`. The bound
`0.396n²/log2 n` therefore sits a factor `log2 n/6.6` *below* the cost the rule
actually pays — the standard `log` gap, and no conjecture about `delta`'s
structure beyond (2.2) is needed for it.

### 2.7 Corollary D6 (what a counterexample must look like; proved)

Combine (2.1) with the *lower* end of what is known about `Σ_k CC_k`. If
`Σ_k CC_k ≥ c n²` for a constant `c > 0` (measured: `c ≈ 0.383` at `n = 6..13`;
see §2.8), then for `s = log2 n + O(1)` and `T ≤ n·polylog(n)`:

```
J̄  ≥  c n² / (s·T)  −  n/T  =  Ω( n / (log2 n · polylog(n)) ) .            (2.5)
```

**This sharpens Lemma 7 of `LOWER_BOUND_ATTEMPT.md`.** That lemma says a proof
of `T ≥ n polylog n` must exhibit `Ω(n polylog n)` *repeat* reads. Corollary D6
says the repeats must additionally be **far apart**: the machine must revisit
cells at mean distance `Ω(n/log³ n)` (for `polylog = log²`). A machine that
re-reads its own neighbourhood — which is what every known simulator does, and
what the "horizon jump" family does — is *provably* not a counterexample: its
`T·J̄` product is too small. So the object to rule out is not "a machine that
re-reads" but "**a machine that teleports**": one whose read sequence looks like
adaptive random sampling at scale `n/polylog`.

Two consequences worth stating in the paper:

1. It makes the conjecture *falsifiable by a cheap statistic*. Any proposed
   `n polylog n` simulator can be rejected by measuring `J̄`; if `J̄ = o(n/log³n)`
   it cannot be correct, whatever else it does.
2. It localises the difficulty in the *order* of reads, not in their number.
   `T` alone cannot separate anything (Lemma 4: `s ≥ log2 T`); `T·J̄` can.
   The conjecture is the statement that the product cannot be beaten, and that
   is a statement about locality of access.

### 2.8 What is measured

`lb_cuts.py` reads the cached deltas (`cert_n*.npy`, row 0, enum order — the
same walk as `lb_cert2.py`) and computes, for every cut, the exact one-way
communication complexity in both directions (`log2 #distinct rows` and
`log2 #distinct columns`), plus two lower bounds on the two-way value at the
middle cut (log-rank of the largest full count-balanced sub-rectangle, and a
greedy fooling set).

**(i) The aggregate is at its ceiling.** `Σ_k min(dir)`, the upper bound on
`Σ_k CC_k` obtained from one-way protocols:

| n | tapes | Σ_k CC_k (upper bd) | 0.396n² | ratio |
|---|---|---|---|---|
| 6 | 141 | 13.60 | 14.26 | 0.954 |
| 8 | 1,107 | 24.53 | 25.36 | 0.967 |
| 10 | 8,953 | 38.61 | 39.62 | 0.974 |
| 11 | 24,068 | 46.19 | 47.95 | 0.963 |
| 12 | 73,789 | 55.78 | 57.06 | 0.978 |
| 13 | 201,643 | 64.92 | 66.96 | 0.969 |

So `Σ_k CC_k ≈ 0.383 n²`, i.e. **97 % of the largest value the bound (2.2)
permits**. There is nothing to be gained by looking at cuts more cleverly: the
tape's information is already almost fully "communicated" in the one-way sense,
which is what makes (2.2) tight and the ceiling in §2.5 real.

One trap, since it is easy to get wrong when reading the tables: the **larger**
of the two one-way directions can *exceed* `min(k,n−k)·log2 3` (at `n = 6`,
`k = 1`: `4.00` bits one way versus `1.58` the other). The valid upper bound on
`CC_k^{2way}` is the *smaller* direction (either protocol may be used), which is
what the `Σ_k min(dir)` column sums; summing the larger direction gives
`1.5–1.6×` the ceiling and is meaningless.

**(ii) The two-way value is Θ(n), empirically.** At the middle cut
(`k = n//2`), in bits:

| n | k | CC_2way ≥ log-rank | block (rows×cols) | rank/min-dim | CC_1way (better dir.) | min-side·log2 3 |
|---|---|---|---|---|---|---|
| 6 | 3 | 1.58 | 3×3 | 1.00 | 4.09 | 4.75 |
| 8 | 4 | 3.32 | 12×12 | 0.83 | 5.61 | 6.34 |
| 10 | 5 | 4.64 | 30×30 | 0.83 | 7.17 | 7.92 |
| 11 | 5 | 4.81 | 30×90 | 0.93 | 7.11 | 7.92 |
| 12 | 6 | 6.44 | 90×90 | 0.97 | 8.72 | 9.51 |
| 13 | 6 | 6.23 | 90×210 | 0.83 | 8.61 | 9.51 |
| 14 | 7 | 7.68 | 210×210 | 0.98 | — | 11.09 |

(The `n = 14` entry is the only one not computed from the full ensemble: `lb_rank14.py`
enumerates just the block. The three largest `210×210` blocks there give ranks
`92, 205, 133`, i.e. `0.44, 0.98, 0.63` of their width — so *not every* profile is
near-full, and Conjecture R below must be read as "some block", which is all the
cut bound needs, since `CC_k` is a maximum over blocks.)

Fooling sets are *not* a usable second tool here: an unrestricted greedy gives
2.8–5.2 bits, but those sets are unsound (the cross pairs leave the promise),
and after restricting to one `d = #0−#2` class — the condition for
`(x_i, y_j)` to be a tape at all — they collapse to `1.00, 1.58, 2.32, 2.81,
2.81, 3.58` bits. So **log-rank is the only two-way lower bound that scales**,
and it is capped by the block: `log2 rank ≤ log2 min(rows, cols)`.

That cap is much less damaging than it looks, and this is the sharpest
formulation of what route (d) still needs. A count-profile block has
`Θ(3^k/poly(k))` rows and columns, so `log2 min-dim = k log2 3 − O(log k)`;
and the measured `rank/min-dim` is `0.83–0.98` throughout. Hence the log-rank
bound is *already of the right order* — it loses only `O(log k)` — and the whole
route becomes unconditional (up to that additive `O(log k)`, which changes
nothing since `Σ_k O(log k) = O(n log n) = o(n²)`) if the following purely
combinatorial statement holds:

> **Conjecture R.** For every `n` there is a count profile whose balanced block
> of the `delta` matrix (rows = left halves with that profile, columns = right
> halves with the complementary one) has rank `≥ c·` its smaller dimension, for
> some absolute constant `c > 0`. Equivalently: for every `n` there is an
> explicit `c·3^{k}/poly(k)` square submatrix that is nonsingular over `ℝ`.

Conjecture R is a statement about the linear algebra of `delta` with no
machine, no read model and no space in it — an object a combinatorialist can
attack by writing down a triangular submatrix (order the halves by a
"first-difference" or "prefix-average" statistic and look for a staircase). It
is verified above (for the largest block at each `n ≤ 13`) to within a constant,
and it is the *only* missing input to route (d): with it, `Σ_k CC_k = (0.396 − o(1))n²` and Corollaries D4–D6
hold unconditionally.

**(iii) The trajectory the bound is calibrated against.** `steps/n² =
0.048–0.068`, `epochs = 0.51–0.65 n`, mean scan `0.10–0.13 n` (earlier file §7).
The number of *moves* is therefore `Θ(0.06 n²)`, and route (d)'s pointer bound
`0.396n²/log2 n` sits a factor `log2 n/6.6` below it — the standard `log` gap, and
the same gap that separates any communication bound from a time bound.

### 2.9 The three ways out, and why each fails

For route (d) to be beaten one must break one of its three inputs:

1. **Beat (2.2)** — impossible: `CC_k` is at most the smaller half's length,
   because one player can send that half. (This is route (a)'s ceiling restated.)
2. **Prove `J̄` small** — the anti-teleport lemma, §2.7. This is a genuine
   branching-program lower bound in disguise: a machine *is* a BP of width `2^s`
   (Lemma 5 of the earlier file), `J̄ ≤ J` is the restriction that consecutive
   queries be `J`-close, and the conjecture is exactly that lifting this
   restriction does not help. Nothing in the deterministic read-only toolbox
   proves anything of the kind.
3. **Use more cuts** — there are only `n−1`, all used, and (2.2) is tight on
   them (97 % measured). No weighting helps, since the bound is a sum of
   non-negative terms each already at its maximum.

So the route is closed at `0.396n/log2 n`, and the closure is a theorem about
cut arguments, not a gap in this write-up.

---

## 3. The configuration-counting route (parent's idea (a)) is refuted

The suggestion: the run is a path of length `T` in a graph of `2^s` states; the
output is constant on the set of tapes reaching a final state; count distinct
deltas against `2^s`, and compare with the number of tapes sharing a final state.

**Proposition 3.1 (what the counting gives).** Let `q(a)` be the final state.
`delta` is a function of `q`, so `#distinct values = n+1 ≤ 2^s`, i.e.
`s ≥ log2(n+1)` — and that is the whole content.
*Proof.* `q` determines the output. ∎

**Proposition 3.2 (the refinement has the inequality backwards; proved).** Say a
class `A_q` "fits a cylinder" if `A_q ⊆ {a : a|_S = w}` for some set `S` of
`r_q` cells, so `|A_q| ≤ 3^{n−r_q}`. By pigeonhole some `q` has
`|A_q| ≥ |A_n|/2^s`, and therefore

```
r_q  ≤  n − log_3|A_q|  ≤  n − log_3|A_n| + s·log_3 2
     =  ½·log_3 n + log_3 C + 0.631 s      (|A_n| = 3^n/(C√n))
     =  0.95 log2 n   at  s = log2 n.                                        (3.1)
```

So the counting bounds the read set of *some* tape from **above** by `O(log n)`
cells. That is not a step towards `T ≥ n polylog n`; it is the opposite
direction, and it is vacuous. The mechanism is worth naming: **reads and class
size are anti-correlated.** A machine with `s` bits has few classes, so some
class is large, and a large class is one the machine got *cheaply* — many tapes
whose behaviour is forced by a few cells. Pigeonhole over final states therefore
always certifies the existence of *easy* tapes, never of hard ones. (Sanity
check at `n = 12`, `s = log2 n`: (3.1) gives `r ≤ 4.06`, against the true
`min_a c(a) = 3` — consistent, and about as far from `n−1 = 11` as it is
possible to be.)

**Proposition 3.3 (and the hypothesis of 3.2 is false anyway).** `A_q` need not
fit any cylinder: two tapes can reach the same state by different read
sequences, and then no single `S` covers both. What *is* true (Lemma 1 of the
earlier file) is that each **trajectory class** is a cylinder; but trajectory
classes *tile* `A_n`, so summing `|class| ≤ 3^{n−r}` over the tiling reproduces
`|A_n| ≤ 3^n` and yields nothing. Counting is blind to `T` because the number of
*leaves* is bounded by `3^T`, not by `2^s`, while the number of *final states* is
bounded by `2^s` but carries no certificate information.

A last nail, which also disposes of every entropy/direct-product variant: the
range of `delta` is `n+1`, so any argument that charges the output's
`log2(n+1)` bits against the `s` bits of state proves at most (3.1)'s
`s ≥ log2(n+1)`. There is no product structure to exploit: the machine is not
required to reproduce the trajectory, only to end in the right state, and
`log2(n+1)` bits of state suffice for that. This is why `s = O(log n)` is not
merely allowed but *optimal* (§4).

---

## 4. Idea (c) is refuted: the space complexity is exactly `Θ(log n)`

Idea (c) wants a reduction to a function with `Ω(n)`-bit one-pass complexity and
then "a machine with `O(log n)` bits cannot do it even with unbounded reads".
The conclusion is false, and provably so.

**Proposition 4.1 (upper bound; design argument, not implemented).** `delta` is
computable with `s = O(log n)` bits and `O(n⁴)` reads.

*Construction.* Simulate the walk literally. The state is `(L, R, ell, r)` — four
positions, `4 log2 n` bits — plus, between steps, one window statistic. Each step
needs `b2` and `b1` for the current window:
* `b2`: compare `Σz` and `Σ z[0..m/2)`; both are `Θ(log n)`-bit numbers and can
  be recomputed from the tape in `O(m)` reads (or maintained in `O(1)` reads per
  step, since the window changes by one cell, at the cost of re-reading the cell
  dropped from the half-prefix).
* `b1`: `longest run of 2s in z > m/4`. Recompute by one scan of `z`, `O(m)`
  reads.
An epoch has at most `(R−L+1)² ≤ n²` states `(ell,r)`, so if it has not halted
after `n²` steps it has repeated a state and will loop for ever — which is
exactly the `⊥` case, so the simulation caps the epoch at `n²` steps and is
still exact. At most `n` epochs, hence at most `n³` steps, each `O(n)` reads:
`O(n⁴)` reads with `s = 4 log2 n` bits. ∎

Combined with `s ≥ log2(n+1)` (3.1) this pins the space complexity at
`Θ(log n)`. **Therefore no reduction of the kind idea (c) contemplates can
exist**: there is no `ω(log n)` space lower bound to reduce to, because a
`Θ(log n)`-bit machine computes the function. The `Ω(n)`-bit one-pass bound of
§5 is a bound on *one pass*, not on space; that distinction is the entire
question, and idea (c) conflates them.

(For the record, `lb_cyc.py` classifies the non-halting cases, because a literal
simulation must reproduce `⊥`: `meet 33–71 %`, `region<3 25–61 %`,
`cycle 4–7 %` with periods `2, 6, 10, 14, 18` — i.e. `2 mod 4`, growing with
`n`. Cycles are cheap to detect in `O(1)` memory by Floyd/Brent re-running, at a
constant factor in reads, so they do not change Proposition 4.1.)

---

## 5. Idea (b) is refuted: no small-space one-pass simulation exists

Idea (b) wants: "machine with `O(log n)` bits, `O(n polylog n)` reads ⇒ one-pass
streaming algorithm with `O(log²n)` bits ⇒ invoke a streaming lower bound".

**Proposition 5.1 (proved modulo the measured one-way value).** Any one-pass
(left-to-right) algorithm computing `delta` on `A_n` needs
`s ≥ CC^{1way}_{n/2}(delta)` bits, and `CC^{1way}` at the middle cut is
`8.61` bits at `n = 13`, i.e. `0.66n` (measured `0.65n–0.73n` over `n = 6..13`;
`(n/2)log2 3 = 9.51` is the trivial maximum).

*Proof.* Cut at `n/2`. A one-pass `s`-bit algorithm's state after the prefix is
an `s`-bit message; the suffix then determines the output. Hence `s` is at least
the one-way communication complexity, which equals `log2 #distinct rows` of the
cut matrix. ∎

So the simulation contemplated in (b) would need `≈0.7n` bits, not `O(log²n)`.
**Idea (b) is not merely unproved; its hypothesis is false**, and necessarily so:
the random-access model is *strictly* stronger than one-pass at these parameters,
which is exactly the phenomenon the paper is trying to explain. What replaces it
is Theorem D3: a machine with small mean jump is communication-bound; a machine
that escapes the communication bound must have large mean jump; there is no
"simulate it by a stream" step, because streams are one-way and the phenomenon
is two-way.

---

## 6. Where this leaves the conjecture (calibration)

The conjecture is: `delta` has no branching program of width `2^{O(log n)} =
n^{O(1)}` and length `n·polylog(n)`. The literature on *general* (non-oblivious)
branching programs is thin exactly there:

* Beame–Jayram–Saks (JCSS 2001) get the first non-trivial tradeoff for general
  BPs, in the shape "exponential size is needed at length `(1+ε)n`" — a
  *linear-length* statement, i.e. width `2^{Ω(n)}`, the opposite corner;
* Beame's multiparty-communication tradeoffs give `T = Ω(n log²n)` at space
  `S = n^{1−ε}` — superlinear, but only when space is nearly linear, and they
  degrade as `S` shrinks;
* the classical two-lower-bounds results for width `(log n)^{O(1)}` give
  `Ω(n log n)`-type bounds only for **symmetric** functions, in the restricted
  models where they are available; `delta` is neither symmetric, nor read-once,
  nor oblivious;
* "making BPs oblivious requires superlogarithmic overhead" gives
  `Ω(n log(n/S)/log log(n/S))` for pointer jumping, i.e. a `log/loglog` factor,
  in the oblivious model.

At width `n^{O(1)}` and length `n·polylog(n)` **nothing is known for any explicit
function**. So the honest sentence for §8 of `bthm_merged.tex` is:

> The classification is a space-versus-reads tradeoff. Its read half is a
> branching-program lower bound at width `n^{O(1)}` and length `n·polylog n`;
> the three standard families of argument (certificates/adversaries, one-way
> cuts, state counting) each have a provable ceiling, and so does the fourth
> (aggregated cuts with jump charging, Theorem D3): `n−1`, `0.79n/log2 n`,
> `log2 n`, and `0.396n/log2 n` respectively. The sharpest consequence available
> is the three-way tradeoff `T·J̄·s ≥ Σ_k CC_k`, which forces any counterexample
> to teleport with mean jump `Ω(n/log³n)`.

---

## 7. Proved / verified / conjectured

**Proved (unconditionally, in the model of §1).**

1. Lemma D1: `TV = Σ_k x_k`.
2. Lemma D2: `(x_k+1)s ≥ CC_k` for every cut.
3. Theorem D3: `T·J̄·s ≥ Σ_k CC_k − n s`.
4. (2.2): `Σ_k CC_k ≤ log2(3)⌊n²/4⌋ = 0.396n²`; ceiling `T ≥ 0.396n/log2 n − 1` in
   the random-access model, `T ≥ 0.396n²/(s) − n` for `J̄ = 1`.
5. Corollary D6: any `n·polylog(n)`-read `O(log n)`-bit machine has mean jump
   `Ω(n/(log n · polylog n))` (given `Σ_k CC_k = Ω(n²)`).
6. Propositions 3.1–3.3: the configuration-counting route gives exactly
   `s ≥ log2(n+1)`, and its refinement runs the wrong way: pigeonhole produces a
   *large* class, and a large class forces a *small* read set
   (`r ≤ 0.95 log2 n` at `s = log2 n`). Final states are not cylinders in any
   case; trajectory classes are, but they tile `A_n`.
7. Proposition 4.1: `s = O(log n)` bits suffice (`O(n⁴)` reads), so the space
   complexity of `delta` is `Θ(log n)`; hence idea (c)'s space lower bound does
   not exist.
8. Proposition 5.1: one-pass needs `≥ CC^{1way}_{n/2}` bits (`0.65n–0.73n` measured
   at `n = 6..13`), so idea (b)'s `O(log²n)`-bit simulation cannot exist.

**Measured (this run, exact over the full ensemble).**

9. `Σ_k CC_k` (better one-way direction) `= 13.60, 24.53, 38.61, 46.19, 55.78,
   64.92` bits at `n = 6, 8, 10, 11, 12, 13`, i.e. `0.954–0.978` of the ceiling
   `0.396n²`.
10. Middle-cut two-way bounds: log-rank `1.58, 3.32, 4.64, 4.81, 6.44, 6.23` bits
    (block rank `0.83–0.97` of the block's smaller dimension); one-way
    `4.09, 5.61, 7.17, 7.11, 8.72, 8.61`. So `CC_2way ∈ [0.72, 1.00]·CC_1way` at
    `n = 12, 13`. Greedy fooling sets give only `1.00–3.58` bits once restricted
    to a `d`-class (the only sound version), so they are not a usable tool.
11. Non-halting census: `meet 33–71 %`, `region<3 25–61 %`, `cycle 4–7 %`,
    cycle periods `≡ 2 (mod 4)`, `≤ 18` at `n = 12`.

**Conjectured.**

12. **Conjecture R** (§2.8): for every `n` *some* balanced block of the `delta`
    cut matrix has rank `≥ c·(its smaller dimension)`, `c > 0` absolute. This is
    the *only* missing input to route (d) — a statement of pure linear algebra,
    no machine, no space, verified (`0.83–0.98`) at `n ≤ 14`. With it,
    `Σ_k CC_k = (0.396 − o(1))n²` and Corollaries D4–D6 hold unconditionally.
13. `delta` admits no `n^{O(1)}`-width, `n·polylog n`-length BP (the paper's
    conjecture). Route (d) does not prove it and provably cannot (§2.5); a proof
    needs the anti-teleport lemma (§2.9, item 2).

---

## 8. Files

```
probe_ologn/lb_cuts.py    per-cut one-way CC (all k, both directions), sum over
                          cuts, log-rank and fooling-set bounds at the middle
                          cut; reads cached deltas from cert_n*.npy
probe_ologn/lb_cuts.log   its output (n = 6, 8, 10, 11, 12, 13)
probe_ologn/lb_cyc.py     non-halting census: meet / region<3 / cycle, periods
probe_ologn/lb_rank14.py  block-only enumeration: rank of the 210x210 blocks at
                          n = 14 (Conjecture R), no full-ensemble walk needed
```

Run with `$P = C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`:

```
$P lb_cuts.py 6 8 10 11 12 13     # ~25 s
$P lb_cyc.py 6 8 10 11 12         # ~40 s
$P lb_rank14.py 14                # ~40 s, block rank at n = 14
```
