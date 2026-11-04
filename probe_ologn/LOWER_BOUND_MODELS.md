# Which machine models admit a separation?  M1–M4

Sixth file in the `probe_ologn` series. `LOWER_BOUND_ATTEMPT.md` closed routes
(a) certificate, (b) counting, (c) adversary; `LOWER_BOUND_ROUTE4.md` closed the
cut family with the conservation law `T·J̄·s ≥ Σ_k CC_k − ns` and the pointer
bound `0.396 n²/s`; `LOWER_BOUND_ROUTE5.md` closed the dyadic family and proved
the additive-memory (model R) bound `T ≥ Σ_e (ℓ_e − 2r)`, hence `Ω(n²)`.

Two model-specific separations exist:

* **single-head pointer machines** (`J̄ = 1`): `T = Ω(n²/log n)` (route 4, Cor. D5);
* **additive-memory exact-decision machines** (model R, `r = O(1)` interval
  sums): `T = Ω(n²)` (route 5, C5.6).

This file asks whether a *third* model gives a clean separation, or whether
these two are all that this family of arguments can reach. Four candidates are
settled, each with a proof, a measurement, or a proof of non-amenability.

New probe: `probe_models.py`, output `probe_models.log` (every number below is
from that log; the run commands are in §8).

---

## 0. Verdict

**One of the four is a third model with a genuine separation, but it is a
separation of the *same* kind: it is the cut argument, applied to a fixed read
order.** The other three yield nothing; M1's general fixed-order case is closed
at `1.5 n` reads by a new prefix-time bound (a factor `1.5` above the certificate
bound, still far below `Ω(n²)`), and its `p`-pass subclass does reach
`Ω(n²/log n)`.

| model | separation? | the statement | status |
|---|---|---|---|
| **M1** offline, one pass | no *read* separation; a **space** separation | `s ≥ m_H(n) = 5.58, 7.32, 8.93` bits at `n = 8,10,12` (`0.70–0.74 n`), min over **all** balanced bipartitions | proved + measured |
| **M1** offline, `p` passes | **yes** (mod. R′) | `T = Ω(n²/s)` = `Ω(n²/log n)`; constant `ρ/(2) ≈ 0.19–0.25` | proved given Conjecture R′ (§2.5), constant measured, asymptotic (vacuous below `n ≈ 20`) |
| **M1** arbitrary fixed order | **linear, with constant > 1** | `T ≥ (n−1) + (m_H − s)/log₂3` = `1.21–1.39 n` measured, `→1.43 n`; the technique caps at `1.5 n` | proved + measured |
| **M2** covering | no | covering is a **normal form**: `C(n) ∈ [G(n), G(n)+n]`; the `O(n)`-read question *is* the open conjecture | proved |
| **M3** random order | no | one pass has `T = n` by definition, and fails on space for **every** order (`s ≥ 8.93` bits at `n = 12`, all 924 cuts); determination is at `Θ(n)` reads, not `Θ(n log n)` | proved + measured |
| **M4** `log² n` bits | no new model | cut ceiling `0.396n²/s` is *decreasing* in `s`; the additive-memory bound *survives* and strengthens in scope to `s = o(n log n)` | proved |

The load-bearing new facts:

1. **The one-way communication complexity of `delta` is essentially the same for
   every balanced bipartition** (exhaustive at `n = 6, 8, 10, 12`, full ensemble:
   minima `4.09, 5.61, 7.17, 8.72` bits, and in each case the *interval* cut is
   the easiest). So a non-adaptive one-pass machine cannot pick a cheap order: it
   needs `≈0.73 n` bits of state, whatever the order. This kills "offline" as a
   counterexample candidate outright.
2. **A `p`-pass fixed order is a hard cut crossed only `2p` times**, so the cut
   lemma gives `p = Ω(n/s)` and `T = Ω(n²/s)` — the third separation.
3. **Every other fixed order is pushed to `1.2–1.4 n` reads, provably, and the
   technique caps at `1.5 n`.** At the time the prefix has covered `n/2` cells,
   the state plus the `|O|` re-read cells determine the whole row function, so
   `s ≥ CC^{1way} − |O|·log₂3` and `|O| ≤ T−(n−1)`: the "repair" of one
   overlapped cell costs exactly `log₂3` bits, giving
   `T ≥ (n−1) + (m_H − s)/log₂3 = 1.21–1.39 n` (measured) and `1.43 n`
   asymptotically, with no bound above `1.5 n` available from either ceiling.
4. **Covering is a normal form** (`C(n) = G(n) ± n`), so it cannot be a third
   model; the asked question "can a covering machine decide in `O(n)` reads" is
   equivalent to the first open rung of the conjecture, and a *yes* would refute
   the conjecture.
5. **At `s = log² n` the cut ceiling is `0.396 n/log² n`, below the certificate
   bound `n−1` for every `n ≥ 4`**, while model R still gives `Ω(n²)`. So more
   memory kills the movement route and leaves the additive route untouched: the
   conservation law of route 4 *cannot* be strengthened in `s`, and the correct
   monotonicity is the opposite of the one hoped for.

The taxonomy in one line: **every argument in this family is a conservation law
`T · (charge) ≥ (content ≤ c n²)`, and the family knows exactly two charges —
geometric (movement/crossings) and algebraic (memory additivity). M1 supplies a
third *setting* for the geometric charge, not a third charge. M2, M3, M4 supply
no charge at all.**

---

## 1. What an argument in this family is

Every inequality proved in the five earlier files has the shape

```
T · (charge per read)  ≥  (content of the tape)  −  (correction),          (1.1)
```

with the right-hand side bounded by the tape's own information. The instances:

| argument | charge | content bound | ceiling |
|---|---|---|---|
| certificate / adversary / transcript | 1 cell per read | `n` cells | `n−1` |
| one cut (route a) | `s` bits per crossing | `min(k,n−k)·log₂3` | `0.79 n/s` |
| all interval cuts + jump (route d) | `J̄·s` | `Σ_k CC_k ≤ 0.396 n²` | `0.396 n/s` (`0.396 n²/s` if `J̄=1`) |
| additive memory (route e) | 1 new endpoint per read | `Σ_e ℓ_e ≤ 0.117 n²` | `0.117 n²` at `r=0` |
| **arbitrary cut, `p` passes (§2.5)** | `s` per crossing, `n` reads per pass | `ρ(n) ≈ 0.5 n` | `0.25 n²/s` |

A **third model** would have to supply a third charge: a structural reason why
one read is worth strictly less than one cell of tape, independent of (i) how far
the head moves and (ii) what the memory can hold. Sections 2–5 test the four
candidates against that requirement; §6 states why the test is exhaustive for
this family.

Two facts used throughout, both from `LOWER_BOUND_ATTEMPT.md`:

* **Lemma A1** (read set = certificate). If two tapes agree on every cell a
  machine reads, it behaves identically on them. Hence the read set of a correct
  run certifies the tape, and `T(a) ≥ c(a)`, the minimum certificate size of `a`.
* **Lemma A2** (every `(n−1)`-set certifies every tape). The missing cell is
  determined by `#0 − #2 = ε`. Hence `c(a) ≤ n−1` and `K_global = n−1`.

`c(a)` and `K_global` are computed exactly here for `n ≤ 12` by enumerating all
`2ⁿ` subsets per tape (`probe_models.py cert`, `owH`); `K_global = n−1` is
re-confirmed by `p_{n−2} = 0`, `p_{n−1} = 1` in the same run.

---

## 2. M1 — the offline (non-adaptive) model

### 2.1 The model, and two corrections to the received picture

A machine is **non-adaptive** if there is a fixed sequence `σ₁,σ₂,…` of positions,
chosen before the tape is seen, and at step `t` the machine reads `a_{σ_t}`. It
may halt at a tape-dependent time `T(a)`, but it cannot choose *where* to look.
This is exactly an **oblivious branching program** of width `2^s` and length `T`
with a fixed variable order (and repeats allowed).

The received picture in the task statement is: "the machine sees a fixed set `S`;
if `|S| ≥ n−1` it can decide, else not; so the read count is `≥ n−1` or the
machine errs." Two corrections, both measured:

**(i) The per-tape bound is `n−2` or `n−1`, not always `n−1`, and the global
certificate bound is the wrong quantity.** Because the read set is a *prefix* of
the fixed order, different tapes may halt at different prefixes; the correct
per-tape statement is `T(a) ≥ c(a)`, so the worst-case read count is at least
`C_max := max_a c(a)`, not `K_global`. Measured (all subsets, exact):

| `n` | 6 | 8 | 10 | 11 | 12 |
|---|---|---|---|---|---|
| `C_max` | 4 | 7 | 9 | 9 | 11 |
| `C_max − (n−1)` | −1 | 0 | 0 | −1 | 0 |
| mean `c(a)` | 2.64 | 3.64 | 4.79 | 5.66 | 5.88 |
| min `c(a)` | 2 | 2 | 3 | 3 | 3 |

So `C_max ∈ {n−2, n−1}`: the received `n−1` is right up to one. What the table
does *not* say is that one fixed set of `n−1` cells is needed: the per-tape minima
are achieved by different sets on different tapes, and a machine with a fixed
order exploits exactly that by halting at different prefixes on different tapes.
For a lower bound the difference is cosmetic; for the protocol argument it is
not, which is why §2.4 restricts to the tapes `H` on which every run is long.

**(ii) Sufficiency is informational, not operational.** It is true that `n−1`
cells determine the tape. It is false that a machine with `O(log n)` bits can
*decide* from them: it cannot store them, and a fixed order cannot go back. This
is the whole content of §2.6–2.7: the model's read complexity is not `n`; §2.6
puts it at `≥ 1.39 n` (measured at `n = 12`, `→ 1.43 n`) and §2.7 exhibits an
`O(n³)` machine.

### 2.2 Lemma M1-1 (the read set is a prefix; proved)

For a non-adaptive machine with order `σ` and halting time `T(a)`, the set of
cells read on `a` is `{σ₁,…,σ_{T(a)}}`. Hence `T(a) ≥ c(a)` and

```
max_a T(a)  ≥  C_max  =  n − 1 or n − 2   (measured, §2.1).
```

*Proof.* Immediate from the definition and Lemma A1. ∎

### 2.3 Lemma M1-2 (cut lemma for an **arbitrary** bipartition; proved)

For every machine (adaptive or not), every tape `a`, and every `A ⊆ [n]`, let
`x_A(a)` be the number of indices `t` with `σ_t ∈ A ≠ σ_{t+1} ∈ Aᶜ`. Then

```
CC_A(delta)  ≤  (x_A(a) + 1)·s ,                                          (2.1)
```

where `CC_A` is deterministic two-way communication complexity across the
partition `(A, Aᶜ)`.

*Proof.* Alice holds `a|_A`, Bob `a|_{Aᶜ}`; both know `σ`. The player owning the
cell to be read holds the simulated state; at each crossing the `s`-bit state is
sent. Messages: `x_A + 1`. Correctness is Lemma A1's induction. ∎

This is route 4's Lemma D2 with the interval cut replaced by an arbitrary
bipartition. The replacement is what makes M1 tractable: a fixed order has a
*cut* whose two sides are arbitrary sets, and the machine cannot dodge it.

### 2.4 Theorem M1-3 (one pass needs linear space; proved + measured)

Let `M` be non-adaptive with `s` bits, whose read sequence is a **permutation**
`σ` of `[n]` (one pass, at most `n` reads per tape). Put `k = ⌈n/2⌉`,
`A = {σ₁,…,σ_k}`, and `H = {a ∈ A_n : c(a) ≥ k}`. Then

```
s  ≥  CC^{1way}_{A → Aᶜ}(delta | H)  ≥  m_H(n) := min_{|A| = ⌈n/2⌉} CC^{1way}(delta|H).
```

Measured (exact one-way CC = `log₂ #distinct rows`; **every** balanced
bipartition — `C(n,⌈n/2⌉) = 20, 70, 252, 924` cuts at `n = 6, 8, 10, 12`):

| `n` | `H`-fraction | `m_H(n)` bits | `m_H/n` | cuts |
|---|---|---|---|---|
| 6 | 0.41 | 3.58 | 0.60 | 20 (all) |
| 8 | 0.54 | 5.58 | 0.70 | 70 (all) |
| 10 | 0.57 | 7.32 | 0.73 | 252 (all) |
| 12 | 0.59 | 8.93 | 0.74 | 924 (all) |

*Proof.* For `a ∈ H` every certificate has size `≥ k`, so by Lemma M1-1
`T(a) ≥ k`: every tape of `H` is read for at least `k` steps, and the state after
`k` reads is defined on all of `H`. Alice simulates those `k` reads (all
positions lie in `A`), sends her `s`-bit state, and Bob — who knows `a|_{Aᶜ}` and
that all remaining positions lie outside `A` — continues and outputs. This is a
one-way protocol for `delta` on the promise `H`, so `CC^{1way}(delta|H) ≤ s`.
Minimise over `A`. ∎

Three remarks.

1. **Why the promise `H` and not `A_n`.** If the machine may halt before `k`, the
   first message would have to carry the answer, costing `log₂(n+1)+1` extra
   bits. Restricting to `H` removes the correction at a price measured to be
   negligible: the one-way CC on `H` is `5.58` vs `5.61` (n=8), `7.32` vs `7.17`
   (n=10), `8.93` vs `8.72` (n=12) — i.e. *the hard tapes carry all the
   communication*. (On `H` the interval cut is again the cheapest.)
2. **The minimum over bipartitions is attained by the interval cut in every
   exhaustive case (full ensemble).** At `n = 10` (all 252 cuts) and `n = 12` (all
   924 cuts) the minimum equals the middle-interval value, `7.17` and `8.72` bits.
   The median over cuts is higher (`7.65`, `9.33`), so a random bipartition is
   *harder* than the interval one. A non-adaptive machine cannot exploit a clever
   order: the order is a permutation, so its `n/2`-prefix is *some* balanced set,
   and every balanced set costs `≈0.7n` bits.
3. **The dependence on the cut size is the trivial one.** At `n = 10`, a cut of
   size `k ≤ 3` has `CC^{1way} = k·log₂3` exactly (`1.58, 3.17, 4.75`): Alice
   simply transmits her `k` cells. The expensive cuts are the balanced ones
   (`7.17–7.84` at `k=5`). A permutation must pass through `k = n/2`, so the
   machine cannot stay on the cheap cuts.

**Consequence.** One pass over the tape — however cleverly ordered — is
impossible with `O(log n)` bits. The offline model at `s = O(log n)` must
therefore re-read, i.e. `T > n`. This is a *space* statement, and it is the
honest answer to "can the machine read a fixed set of `n−1` cells and decide?",
the operational version of the question in the task: no.

### 2.5 Corollary M1-4 (`p` passes: the third separation; proved mod. R′, constant measured)

Let `σ = π^p` (`p` full passes over a permutation `π`, so at most `pn` reads per
tape). With `A = {π₁,…,π_k}`, `k = ⌈n/2⌉`, we have `x_A ≤ 2p−1`, hence by (2.1)

```
CC_A(delta)  ≤  (2p−1)·s + log₂(n+1) + 1,
p            ≥  (CC_A − log₂(n+1) − 1) / (2s),
T = pn       ≥  n·(CC_A − log₂(n+1) − 1) / (2s)   =   Ω(n²/s).            (2.2)
```

The `log₂(n+1)+1` is the early-halt flag of §2.4; it is `o(n)`.

For `CC_A` we need a **two-way** lower bound; the only scalable tool is the
log-rank of a full count-balanced block (route 4 §2.8). For a *promise* problem
the rank argument must be applied to a full rectangle inside the promise: a
`c`-bit protocol cuts such a rectangle into `≤ 2^c` monochromatic sub-rectangles,
so `rank ≤ 2^c`. A count-profile block is exactly such a rectangle (every pair of
an `A`-half and a `B`-half with complementary profiles is a tape of `A_n`), so
the bound is sound. Measured smallest block rank over balanced bipartitions:

| `n` | cuts | `ρ(n) = min log₂ rank` | `ρ/n` |
|---|---|---|---|
| 8 | all 70 | 3.00 | 0.375 |
| 10 | all 252 | 4.25 | 0.425 |
| 12 | 60 sampled | 6.04 | 0.503 |

**The corollary is conditional, and on exactly the generalisation of route 4's
Conjecture R that the arbitrary cut needs.**

> **Conjecture R′.** There is an absolute constant `c > 0` such that for every
> `n` and *every* balanced bipartition `A` of `[n]`, some count-balanced block
> of the `delta` matrix at `A` has rank `≥ c·` its smaller dimension.

With R′, `CC_A ≥ cn` for every `A`, and (2.2) gives

```
p  ≥  (c n − log₂(n+1) − 1) / (2s),      T = p n  =  Ω(n²/s) .
```

The measurements give `c = 0.375, 0.425, 0.503` at `n = 8, 10, 12`, i.e. a
constant `ρ/2 ≈ 0.19–0.25`: `T ≥ 0.19 n²/s − O(n)`, and `0.19 n²/log₂ n` at
`s = log₂ n`. Route 4 needed R only at the middle interval cut; M1 needs it at an
arbitrary balanced cut, and the exhaustive `n = 8, 10` rows plus the 60-cut
sample at `n = 12` are the evidence that the stronger form holds.

**Status of the constant.** At the sizes that can be enumerated the bound is
vacuous (`p ≥ −0.20, −0.03, +0.19` at `n = 8, 10, 12`); it becomes non-trivial at
`n ≈ 20` and is the usual `log`-factor statement beyond that. The theorem is
asymptotic, like route 4's pointer bound `0.396 n²/log n`, which is below `n−1`
until `n ≈ 500`.

**Why this is a genuinely different model, and not route (d).** A `p`-pass
machine over a wild permutation has mean jump `J̄ ≈ n/3`, so route (d) gives only
`T·J̄·s ≥ 0.396n²`, i.e. `p ≥ 1.2/s` — vacuous. The single-cut version (2.1) is
strong exactly because a pass carries `n` reads while crossing the cut once. So
the two bounds are incomparable: route (d) needs small `J̄` (a *pointer*), (2.2)
needs a *fixed order that sweeps*. Neither implies the other.

**The mechanism, stated plainly.** Covering every cell once per pass forces
`x_A ≥ 1` for every cut and every pass. So the read count is `n` times the
crossing number, and the crossing number is bounded by the number of passes.
The charge is still geometric (crossings), but the geometry is now *chosen by the
designer, not the machine* — which is the whole difference from route (d), where
the machine chose `J̄` and the fix was to sum over all `n−1` interval cuts.

### 2.6 Theorem M1-5 (prefix-time bound: `T ≥ (1.2…1.4) n`, `→ 1.43 n`, cap `1.5 n`)

This is the strongest statement the family gives for a general fixed order, and it
replaces the received `n−1`.

Let `σ` be any fixed order, `k = ⌈n/2⌉`, and suppose the first `k` reads are
distinct (`σ₁,…,σ_k` pairwise different — true for every permutation and every
`p`-pass order). Put `S = {σ₁,…,σ_k}`, `O = S ∩ {σ_{k+1},…,σ_T}` (the cells of
`S` read again later), and `H = {a ∈ A_n : c(a) ≥ k}`. Let `m ≥ n−1` be the
number of distinct cells the order covers (a correct machine covers at least
`n−1`: no `(n−2)`-set certifies every tape, `K_global = n−1`, Lemma 3 and its
re-verification to `n = 12`). Then

```
s  ≥  CC^{1way}_{S → Sᶜ}(delta | H)  −  |O|·log₂3 ,     |O| ≤ T − m ≤ T − (n−1),   (2.3)
```

hence

```
T  ≥  (n−1) + (m_H(n) − s)/log₂3  →  1.43·n   (s = o(n)).                        (2.4)
```

*Proof.* For `a ∈ H`, `T(a) ≥ c(a) ≥ k`, so the state after `k` reads exists on
every tape of `H`. Write a tape as `(x, z)` with `x = a|_S`, `z = a|_{Sᶜ}`. The
first `k` reads are all in `S`, so the state `q(x)` depends on `x` alone. The
later reads are at positions in `Sᶜ` and in `O`; their values are `z` and
`x|_O` respectively. Hence the whole run after `k`, and therefore the answer, is
a function of `(q(x), x|_O, z)`. So the row function `z ↦ delta(x,z)` is
determined by `(q(x), x|_O)`, and

```
#distinct rows  ≤  2^s · 3^{|O|} .
```

The left side is `2^{CC^{1way}_{S→Sᶜ}(delta|H)}` (one-way CC = log₂ #distinct
rows). Finally `|O| ≤ T − m`: every cell of `O` is read at least twice, and
`Σ_i (m_i − 1) = T − m` over the `m` covered cells. ∎

(If the order repeats a cell within its first `k` reads, apply the same argument
at the time of the `k`-th distinct cell; the promise `H` shrinks and the
early-halt flag of §2.4 costs `O(log n)` bits, so the asymptotic bound is
unchanged and the measured constants move by `o(n)`.)

**The numbers.** `m_H(n) = 3.58, 5.58, 7.32, 8.93` bits at `n = 6, 8, 10, 12`
(§2.4, all balanced cuts), `m_H/log₂3 = 2.26, 3.52, 4.62, 5.63`:

| `n` | 6 | 8 | 10 | 12 |
|---|---|---|---|---|
| `T ≥ (n−1) + m_H/log₂3` (s = o(n)) | 1.21 n | 1.32 n | 1.36 n | 1.39 n |
| `T ≥ (n−1) + (m_H − log₂n)/log₂3` | 0.94 n | 1.08 n | 1.15 n | 1.20 n |
| technique cap, `(n−1) + 0.79n/log₂3` | 1.33 n | 1.37 n | 1.40 n | 1.42 n |

So a non-adaptive `O(log n)`-bit machine must read **more than `n` cells**, with a
measured constant `1.21–1.39` at `n = 6…12` and limit `1.43`; this is the correct
replacement for the received `T ≥ n−1` (an improvement by a factor `1.5` at
`n = 12`, `1.43` asymptotically).

**Why the technique stops at `1.5 n`.** The bound (2.3) has two ceilings:
`|O| ≤ |S| = n/2`, and `CC^{1way} ≤ min(|S|,|Sᶜ|)·log₂3 = 0.79 n` (Alice
transmits her half; measured maximum over cuts `9.48` bits at `n = 12`). Hence
`T ≥ (n−1) + (0.79n − s)/log₂3 → 1.5 n`, and no choice of `t`, no weighting and no
repetition of the argument can exceed it: the "repair" of an overlapped cell costs
exactly `log₂3` bits of state, and there are never more than `n/2` of them before
the prefix is balanced. **The general non-adaptive model is thus closed at
`T ≥ 1.5 n` for this family**, with `1.43 n` achieved.

**What one would need for more.** To go past `1.5n` the state would have to be
charged against something other than the row count — i.e. one needs the *cut*
argument with many crossings, which is §2.5's `p`-pass setting, or a proof that
`delta` has no small-width BP at all. Note that the row-function count is
unavoidably capped by `3^{|S|} = 3^{n/2}`, since `S` has only `n/2` cells; this is
the same cap as the trivial protocol and is the exact analogue, for M1, of route
4's `Σ_k CC_k ≤ 0.396 n²` ceiling.

### 2.7 Theorem M1-6 (an `O(n³)` non-adaptive machine; constructed and executed)

`delta` is computable by a non-adaptive `O(log n)`-bit machine with
`T(a) = n·Σ_e ℓ_e(a) ≤ n³` reads, on the fixed order `0,1,…,n−1` repeated.

*Construction.* One full pass over the tape computes `b2(L,m)` for a *single*
pair `(L,m)`: the machine keeps `P = Σ_{[L,L+k)} z` and `T = Σ_{[L,L+m)} z` as it
sweeps `0..n−1` (ignoring cells outside `[L,L+m)`), then evaluates
`b2 = [k·T > m·P]` (Lemma 5.1). The state is `(L, R, m, P, T, clock)` —
`O(log n)` bits. The epoch loop of Lemma 5.2 probes `m = M, M−1, …` at anchor
`L`; at most `M ≤ n` probes per epoch and at most `n` epochs, so at most `n²`
passes of `n` reads.

*Verified* (`probe_models.py m1ub`) against `walk_b0` on **all** ensemble tapes
for `n = 6, 8, 10, 11` — 34,269 tapes, 0 mismatches:

| `n` | tapes | agree | mean reads`/n³` | max reads`/n³` |
|---|---|---|---|---|
| 6 | 141 | 141 | 0.1168 | 0.3056 |
| 8 | 1107 | 1107 | 0.1121 | 0.3438 |
| 10 | 8953 | 8953 | 0.1077 | 0.3500 |
| 11 | 24068 | 24068 | 0.1126 | 0.3636 |

The mean is `n·Σ_e ℓ_e/n³ = Σ_e ℓ_e/n² ≈ 0.11`, matching route 5's measured
`Σ_e ℓ_e ∈ [0.056, 0.117] n²`. So the offline model is not vacuous: it has a
polynomial machine with a small constant. The gap `[(1.39…1.43)n, n³]` is the
exact statement of what is not known.

### 2.8 M1 verdict

* **`Ω(n)` or `Ω(n²)`?** `Ω(n)` — with a constant strictly above 1:
  `T ≥ 1.39 n` measured at `n = 12` (`1.21–1.39 n` over `n = 6…12`, `→ 1.43 n`),
  not the received `n−1`. `Ω(n²)` is **not** available for the model as stated;
  the only superlinear statement in the family is the `p`-pass one (§2.5).
* **The third separation is `p`-pass** (§2.5): `T = Ω(n²/s)`, i.e.
  `Ω(n²/log n)` at `s = log n`, with a measured constant `≈0.2`, conditional on
  Conjecture R′. It is a cut argument and therefore *inside* the family, not a
  new mechanism.
* **One pass is impossible** at `O(log n)` bits (`s ≥ 0.70–0.74 n`, measured over
  all balanced cuts), which is the useful negative result: the natural "offline
  counterexample" cannot exist.
* **The general fixed order is closed at `1.5 n`**: the prefix-time bound's two
  ceilings (`|O| ≤ n/2`, `CC ≤ 0.79 n`) meet there. So the family's complete
  answer for M1 is `(1.39…1.43)n ≤ T ≤ n³`, with `T = Ω(n²/log n)` on the
  `p`-pass subclass and no technique for the gap.

---

## 3. M2 — the covering model

### 3.1 Lemma M2-1 (covering is a normal form; proved)

Let `G(n, s)` be the minimum over all `s`-bit machines of `max_a T(a)`, and
`C(n, s)` the same over covering machines (every cell read at least once before
halting). Then

```
G(n, s)  ≤  C(n, s)  ≤  G(n, s + O(log n)) + n .
```

*Proof.* Left: a covering machine is a machine. Right: prefix any machine `M`
with the fixed cover `0,1,…,n−1` (an `O(log n)`-bit counter); the result reads
`T_M(a) + n` cells and covers all of them. ∎

**Consequences.** (i) `C` and `G` differ by at most `n`, so a covering model is
*not* a separate model: a lower bound for `C` is a lower bound for `G` up to
`−n`, and an upper bound for `G` is one for `C` up to `+n`. (ii) Hence
"can a covering machine decide in `O(n)` reads?" is equivalent, up to constants,
to "can any machine?" — the first open rung of the conjecture (`T = ω(n)`), since
the best unconditional bound is `n−1`. (iii) A *yes* would refute the conjecture
outright (an `O(n)`-read covering machine is an `O(n)`-read machine). So no
counterexample is to be expected, and no proof is available without the
conjecture.

### 3.2 Why the suggested first-coverage order does not help

The suggestion: order the cells by *first coverage* `π(a)` and apply a cut
argument at the moment the `k`-th new cell is read. It fails for exactly the
reason of Theorem M1-5, and the failure is structural, not technical:

* the cut at the `k`-th first coverage has `S` = the covered prefix and `U` = all
  reads after it, and a covering machine may re-read the prefix arbitrarily;
* in the natural covering machine (`cover 0..n−1, then simulate`) the overlap is
  the whole prefix, `|O| = n`, so (2.3) gives `s ≥ CC^{1way} − 1.58n < 0`: nothing;
* and the covering constraint adds exactly one read over the certificate bound:
  `T ≥ n`, against `C_max = n−1` (n = 8, 10, 12) measured. Covering buys one read.

### 3.3 M2 verdict

**Decided, negatively.** Covering is a normal form (`C = G ± n`), so it is not a
third model and neither a proof nor a counterexample for it can exist
independently of the conjecture: proving `ω(n)` for covering machines would prove
it for all machines (add `n`), and an `O(n)`-read covering machine would refute
the conjecture. The only quantitative content of the covering restriction is
`T ≥ n`, one read above the certificate bound `C_max = n−1` (measured at
`n = 8, 10, 12`), and the first-coverage order supplies no charge (§3.2).

---

## 4. M3 — the random order

Two readings of the model; both are settled, and neither gives a separation.

### 4.1 Random *permutation* (one pass): `T = n`, and the model fails on space

If `π` is a uniformly random permutation and the machine reads `π₁,…,πₙ`, then
`T = n` on every tape: there is **no read count to lower-bound**. The only
question is whether the model can be correct, and Theorem M1-3 answers it for
*every* order, random or not: `s ≥ m_H(n)`, measured `5.58, 7.32, 8.93` bits at
`n = 8, 10, 12` over **all** balanced bipartitions (`70, 252, 924` cuts), so no
permutation (random or adversarial) does better than `8.93` bits `= 0.74 n` at
`n = 12`. The `n/2`-prefix of a random permutation is a uniformly random balanced
set, and the *median* cut costs `9.33` bits, so a random order is worse than the
best order, never better.

So M3, read as one pass, is not a counterexample candidate: it needs linear
space. Randomness of the order cannot help because the bound is worst-case over
the ensemble and holds for every order.

### 4.2 Random *sequence with replacement*: `Θ(n log n)` to cover, `Θ(n)` to determine

Here the machine reads i.i.d. uniform cells; covering takes `n H_n` reads.
Measured (`probe_models.py cover`, 2000 trials): `E[T_cover]/n H_n = 1.006, 1.016,
0.999, 1.004, 1.002, 1.002` at `n = 8, 12, 16, 32, 64, 128`, so the coupon
collector is the textbook `n H_n = n ln n + γn + 1/2 − O(1/n)`.

But covering is not what determines the answer. The answer is determined by the
*read set* `S` becoming a certificate, and the measured determination curve says
that happens only at `Θ(n)` distinct cells. With `p_k(a)` = fraction of
`k`-subsets that certify `a`:

| `n` | mean `c(a)` | median | `mean_a p_k` at `k=n/4` | `k=n/2` | `k=3n/4` | `k=n−1` |
|---|---|---|---|---|---|---|
| 6 | 2.64 | 2 | 0.10 | 0.31 | 0.61 | 1 |
| 8 | 3.64 | 4 | 0.01 | 0.14 | 0.61 | 1 |
| 10 | 4.79 | 5 | 0.00 | 0.06 | 0.63 | 1 |
| 11 | 5.66 | 6 | 0.00 | 0.03 | 0.25 | 1 |
| 12 | 5.88 | 6 | 0.00 | 0.03 | 0.36 | 1 |

So a random half of the tape determines the answer with probability `0.03`
(`n = 12`), and a random three-quarters with probability `0.36`. The
determination threshold is `n−1`, i.e. the full cover. The `Θ(n log n)` coupon
collector is an artifact of sampling *with* replacement and is not the binding
constraint: distinct cells, not reads, are what matter, and the answer is not
determined much earlier.

### 4.3 M3 verdict

No separation: one pass has `T = n` by definition (nothing to separate) and is
infeasible at `O(log n)` bits; the with-replacement version needs `Θ(n)` distinct
cells and `Θ(n log n)` reads only because of replacement. Randomness of the order
buys nothing, since the space bound holds for every order.

---

## 5. M4 — machines with `log² n` bits

`T ≤ 2^{s+O(1)} = n^{O(log n)}` at `s = log² n`, so no contradiction with
anything is available from state counting. The question is whether the *known*
bounds improve. They do not; one degrades and one strengthens.

### 5.1 The cut ceiling is monotone decreasing in `s` (proved)

Route 4's law is `T·J̄·s ≥ Σ_k CC_k − ns` with `Σ_k CC_k ≤ log₂(3)·n²/4`. Hence

```
T·J̄  ≥  0.396 n²/s − n ,       ceiling:  T ≥ 0.396 n²/(s·(n−1)) ≈ 0.396 n/s.   (5.1)
```

Every factor of `s` is a *loss*: the per-cut protocol `(x_k+1)s ≥ CC_k` charges
`s` bits per crossing, so a bigger state buys fewer crossings, and the ceiling is
`0.396 n/s` for unrestricted access and `0.396 n²/s` for `J̄ = 1`. At `s = log² n`
(verified numerically, `probe_models.py ceil`):

| `n` | `n−1` | `0.396 n/log²n` (random access) | `0.396 n²/log²n` (`J̄=1`) |
|---|---|---|---|
| 16 | 15 | 0.40 | 6.3 |
| 64 | 63 | 0.70 | 45.1 |
| 256 | 255 | 1.58 | 405.8 |
| 1024 | 1023 | 4.06 | 4154.9 |
| 10⁶ | 999999 | 997.4 | 9.97·10⁸ |

So at `s = log² n` the random-access ceiling is `0.396 n/log² n`, which is below
`n−1` for every `n ≥ 4` (since `0.396 < log² n`): the movement route is *dead*,
not merely weaker. The pointer ceiling stays superlinear
(`0.396 n²/log² n`) but is a factor `log n` worse than at `s = log n`
(`0.396 n²/log n`). **The conservation law cannot be strengthened at larger `s`**:
its right-hand side is capped by the tape, and its charge is `s` per crossing by
construction.

### 5.2 The additive-memory bound survives, and strengthens in scope (proved + measured)

Route 5's model R has `r ≤ s/(3 log₂ n)` interval sums, and

```
T  ≥  Σ_e max(0, ℓ_e − 2r)  ≥  Σ_e ℓ_e − 2rE  ≥  (0.0565 − 1.3 r/n)·n² .    (5.2)
```

At `s = log² n`, `r ≤ log₂n/3`, so `1.3 r/n = O(log n/n)` and

```
T  ≥  (0.0565 − O(log n / n))·n²  =  Ω(n²)      for every s = o(n log n).
```

Verified: `(0.0565 − 1.3 log₂n/(3n))` `= 0.0430, 0.0523, 0.0559, 0.0565` at
`n = 256, 1024, 10⁴, 10⁶`, and the measured checkpoint sums
`Σ_e max(0, ℓ_e − 2r)/n²` at `r = log₂ n` are `0.0321` (`n=64`), `0.0326`
(`n=128`), still `Ω(n²)`; they only die at `r ≈ n/4` (`0.0039`, `0.0029`) and
`r = n/2` (`0`), which needs `s ≈ n log n` bits.

So the second separation is **not** a small-memory phenomenon: it holds for all
`s = o(n log n)`, i.e. all the way to just below `n log n` bits. This is the one
place in the whole family where raising the memory budget *strengthens* a
theorem's scope rather than weakening its constant.

### 5.3 No third mechanism at larger memory

A machine with `log² n` bits can hold `Θ(log n)` interval sums (model R, covered),
or an arbitrary `s`-bit summary (argmin, hull vertex, run length, cached bits) —
which is exactly the general model at `s = log² n`, for which no technique
applies. C5.13 becomes weaker as `s` grows: a counterexample at `s = log² n` must
hold a non-additive summary of information content `ω(log n)` interval sums,
rather than `ω(1)`.

### 5.4 M4 verdict

No third model. The correct monotonicity is: **movement degrades as `1/s`,
additivity is unchanged**. Hence at large `s` the additive-memory separation is
the only one left standing — and it is the second known model, not a third.

---

## 6. Why the family reaches exactly these settings

The family's arguments are conservation laws (1.1). For a *read* lower bound one
needs a charge that is `Ω(1)` per read and is not chosen by the machine. The
family knows two such charges, and M1's `p`-pass case is a third *setting* of the
first:

1. **Geometric, machine-chosen (route d).** The machine pays `J̄` per read. Since
   it chooses `J̄`, one sums over all `n−1` interval cuts: `T·J̄·s ≥ ΣCC`, and
   `ΣCC ≤ 0.396n²` caps the technique at `0.396n/s`. Pointer machines (`J̄=1`) get
   `0.396n²/s`.
2. **Algebraic (route e).** An additive memory cannot represent a contrast
   without owning both endpoints; each read attaches one vertex. `T ≥ Σℓ − 2rE`
   with `Σℓ ≤ 0.117n²`, giving `Ω(n²)` for `r = o(n)`.
3. **Geometric, machine-fixed (M1, `p`-pass).** A fixed order that sweeps gives
   `n` reads per crossing, so `T = pn` and `p ≥ CC_A/(2s)` with `CC_A ≥ ρ(n)`;
   `Ω(n²/s)`. This is charge 1 with the geometry removed from the machine's
   control.

Nothing else in the tested set supplies a charge:

* **covering** supplies no restriction on the read *order*: the machine is
  adaptive, so there is no fixed `σ` to apply the prefix-time bound (M1-5) to,
  and it may re-read anything; it is a normal form anyway (§3.1);
* **random order** is a distribution over fixed orders — the bound is worst-case
  over orders, so the distribution is irrelevant;
* **more memory** weakens charge 1 (`1/s`) and leaves charge 2 intact; it cannot
  create a charge.

The two ceilings that make the exhaustion quantitative are Ceiling 5.8 (any
argument whose content is the tape's information caps at `n`) and Ceiling 5.11
(the valid information form of a per-round argument is again `n`). To exceed `n`
one must prove *redundancy*, and redundancy is a charge. So a genuinely third
model would need a third, non-geometric, non-algebraic reason why a read cannot
be worth a bit of the tape. None of M1–M4 has one; M1 has a third geometry.

---

## 7. Proved / measured / open

**Proved.**

1. **M1-1**: the read set of a non-adaptive machine is a prefix of its order, so
   `T(a) ≥ c(a)` and `max_a T(a) ≥ C_max`.
2. **M1-2**: the cut lemma `CC_A ≤ (x_A+1)s` for an arbitrary bipartition `A`.
3. **M1-3**: one pass over a permutation needs `s ≥ CC^{1way}(delta|H)` for the
   `n/2`-prefix, hence `s ≥ m_H(n)`; no one-pass `O(log n)`-bit machine exists.
4. **M1-4**: `p` passes give `CC_A ≤ (2p−1)s + log₂(n+1) + 1`, hence
   `T = Ω(n²/s)` given `CC_A = Ω(n)` — i.e. given Conjecture R′, which is
   verified at every balanced cut for `n = 8, 10` and at 60 cuts for `n = 12`.
5. **M1-5** (prefix-time bound): with `O` = the cells of the `n/2`-prefix re-read
   later, `s ≥ CC^{1way}_{S}(delta|H) − |O|log₂3` and `|O| ≤ T−(n−1)`, hence
   `T ≥ (n−1) + (m_H − s)/log₂3` (`1.21–1.39 n` measured, `→1.43 n`); the
   technique caps at `1.5n` (`|O| ≤ n/2`, `CC ≤ 0.79n`).
6. **M1-6**: an `O(log n)`-bit non-adaptive machine with `n·Σ_e ℓ_e ≤ n³` reads.
7. **M2-1**: `G(n,s) ≤ C(n,s) ≤ G(n,s+O(log n)) + n`; covering is a normal form.
8. **M4-5.1**: the cut ceiling `0.396n²/s` is strictly decreasing in `s`; at
   `s = log²n` the random-access ceiling `0.396n/log²n` is below `n−1` for all
   `n ≥ 4`.
9. **M4-5.2**: model R gives `Ω(n²)` for every `s = o(n log n)`.

**Measured (exact over the full ensemble unless stated).**

10. `C_max = 4, 7, 9, 9, 11` at `n = 6, 8, 10, 11, 12`; mean `c(a) = 2.64, 3.64,
    4.79, 5.66, 5.88`; min `c(a) = 2, 2, 3, 3, 3`; `K_global = n−1` re-confirmed
    by `p_{n−2} = 0`, `p_{n−1} = 1` at every `n ≤ 12`.
11. One-way CC over balanced bipartitions: minima `4.09, 5.61, 7.17, 7.43, 8.72,
    9.13` bits at `n = 6, 8, 10, 11, 12, 13` (exhaustive at 6, 8, 10, 12 — 20,
    70, 252, 924 cuts; 40 sampled at 11, 13). The interval cut attains the
    minimum in every exhaustive case; the medians are `4.28, 5.95, 7.65, 7.86,
    9.33, 9.47`.
12. On `H = {c(a) ≥ n/2}`: `m_H = 3.58, 5.58, 7.32, 8.93` bits at
    `n = 6, 8, 10, 12` (`0.60–0.74 n`), over **all** `C(n,⌈n/2⌉)` balanced cuts;
    `H`-fraction `0.41, 0.54, 0.57, 0.59`. The interval cut attains the minimum at
    `n = 6, 10, 12`; at `n = 8` the minimum is `5.58` against the interval's
    `5.61`.
13. Log-rank at balanced cuts: min `3.00` (all 70, `n=8`), `4.25` (all 252,
    `n=10`), `6.04` (60 sampled, `n=12`); `0.375–0.503 n`.
14. M1 upper bound: `n·Σ_e ℓ_e` reads, `0.1077–0.1168 n³` mean, `≤0.3636 n³`
    max, 34,269/34,269 tapes agreeing with `walk_b0` at `n = 6, 8, 10, 11`.
15. Determination curve: `mean_a p_k = 0.03` at `k = n/2`, `0.36` at `3n/4`,
    `1` at `n−1` (`n = 12`); coupon collector `E[T]/nH_n ∈ [0.999, 1.016]` at
    `n = 8…128`.
16. Model-R checkpoint sums at `r = log₂n`: `0.0321 n²` (`n=64`), `0.0326 n²`
    (`n=128`); zero at `r = n/2`; `0.0039` at `r = n/4` (`n=64`).
17. Constants: `log₂3/4 = 0.39624`; `(0.0565 − 1.3log₂n/(3n)) = 0.0430, 0.0523,
    0.0559, 0.0565` at `n = 256, 1024, 10⁴, 10⁶`; `p`-pass bound `p ≥ −0.195,
    −0.032, +0.187` at `n = 8, 10, 12`.

**Open / conjectural.**

18. M1's read complexity for arbitrary fixed orders: in `[(1.39…1.43)n, n³]`,
    closed at `1.5 n` by the prefix-time technique (§2.6). In particular
    `T = Ω(n²)` for such machines is *not* proved and not reachable by
    certificates, the cut lemma, movement or counting.
19. The `p`-pass constant `ρ(n)/n`: measured `0.375, 0.425, 0.503` at
    `n = 8, 10, 12`; the natural conjecture is `ρ(n) = (1−o(1))·rank` of the
    balanced block, i.e. `ρ ≈ 0.5–0.7 n` (route 4's Conjecture R in a different
    guise: the block rank is `≥ c·` its smaller dimension).
20. Whether a covering machine can decide in `O(n)` reads: equivalent to the
    first rung of the paper's conjecture (§3.1); open, and a *yes* would refute
    the conjecture.

---

## 8. Files

```
probe_ologn/probe_models.py   all four model probes
        ow     n... [--nsub K]   one-way CC over balanced bipartitions (all, or K sampled)
        owint  n                 one-way CC vs partition size (interval vs random sets)
        owH    n...              one-way CC restricted to H = {c(a) >= n/2}
        rank   n... [--nsub K]   log2 rank of the balanced block at balanced cuts
        cert   n...              exact c(a) for every tape (all 2^n subsets) and p_k
        cover  n...              coupon-collector simulation vs n*H_n
        m4     n...              model-R checkpoint accounting for a grid of r
        m1ub   n...              the O(n^3) non-adaptive machine, executed
        ceil                     every closed-form constant in this file
probe_ologn/probe_models.log  the concatenated output of all of the above
```

Run with `$P = C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`,
in `probe_ologn/`:

```
$P probe_models.py ow 6 8 10                      # ~3 s
$P probe_models.py ow 12                          # ~2 min (924 cuts)
$P probe_models.py owH 8 10 12                    # ~2 min
$P probe_models.py rank 10 12 --nsub 60           # ~1 min
$P probe_models.py cert 6 8 10 11 12              # ~2 min
$P probe_models.py m1ub 6 8 10 11                 # ~3 s
$P probe_models.py m4 12 16 24 32 64 128          # ~6 s
$P probe_models.py cover 8 12 16 32 64 128        # ~10 s
$P probe_models.py ceil
```

`cert` and `owH` use `cert_n*.npy` (row 0 = meeting position `L`, `−1` = no
halt) and `p0_tapes.enum_tapes`; `m1ub` and `m4` verify against `walk_b0` of
`lb_route5.py`, which reproduces `_arr_rules.walk` with `b1 ≡ 0` on all 306,965
ensemble tapes for `n ≤ 13`. No `.tex` file was modified.
