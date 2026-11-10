# Can the O(log n)-bit / O(n polylog n)-read separation be proved?

Attempt at the lower bound for the four walk rules (b1 = longest run > m/4, run
count > m/2, fixed 3-cell pattern, palindrome; b2 = window mean vs first-half
mean; move table of `_arr_rules.walk`). See `REPORT.md` and §8 of
`bthm_merged.tex` for the state of the question and the earlier probe results.

---

## 0. Verdict

**Not proved, and not provable by any of the three proposed routes.** Each route
fails for a reason that can be stated as a theorem, not as a gap in an argument:

| route | best bound it can give | why it stops there |
|---|---|---|
| (a) crossing / one-way communication | `T >= CC_k/s`, and `CC_k <= min(k,n-k)·log2(3)`; at `k=n/2`, `s=log2 n` this is `0.79 n/log2 n` | the number of crossings is chosen by the *machine*, not the adversary (a machine may cross once); the bound then degenerates into a **space** bound for the one-pass subclass |
| (b) counting behaviours / states | `s >= log2 T` (state counting) or `T >= K_global = n-1` (transcript counting) | the machine need not distinguish behaviours, only outputs (there are only `n+1` of those); and the transcript count is *exactly* the decision-tree bound `n-1` |
| (c) adversary | `T >= n-1` | the ensemble's count constraint makes **every** `(n-1)`-set a certificate, so the adversary's ambiguity dies after `n-1` distinct cells and repeats cannot be forced |

The positive structural reason why nothing works is Lemma 4 below: the model
satisfies `T <= 2^s` for every tape, i.e. **state caps reads**, so with
`s = O(log n)` the model is confined to `poly(n)` reads. The target
`T = n·polylog(n)` sits at `s = log2 n + O(log2 log2 n)`, only a lower-order
term above the information-theoretic minimum `log2 T`. So the conjecture
is not a coarse space-versus-reads separation; it is the claim that the tradeoff
curve `T(s)` is a step at `s = Θ(n)` rather than the smooth curve `T ≈ 2^s`.
No technique in the deterministic read-only model proves such a thing for a
function whose range has size `n`.

**What is proved here:** Lemmas 1–6 (§2), the ceilings of routes (a), (b), (c)
(§3–§5). **What is verified computationally:** `K_global = n-1` for
`n = 6..11` (§2.3), the `(n-1)`-sufficiency proof, the four rules' step counts
and epoch structure (§7), and the exact one-way communication complexity across
the middle cut (§3.2). **What remains conjectural:** everything about
`T = n polylog n`.

---

## 1. The model, exactly

Fix `n`. A *tape* is `a ∈ {0,1,2}^n`; the ensemble `A_n` is

```
A_n = { a : #0(a) - #2(a) = eps },   eps = 0 (n even), eps = 1 (n odd),
```

which is the sequence set of two binary columns each with `ceil(n/2)` ones
(`p0_tapes.py`). `delta(a) ∈ {0, 1/(n-1), ..., 1} ∪ {⊥}` is the meeting position
of the walk, `⊥` when the walk does not halt (cycle within an epoch, or region
shorter than 3); `⊥` occurs on 19–67 % of `A_n` depending on `n` (measured, §7).

A *machine* `M` for this `n` has `s = s(n)` bits of memory. Its memory content
determines (i) whether it halts and with which output, and (ii) which cell it
reads next. Reads are at unit cost, the position is arbitrary (random access),
the tape is read-only. `T_M(a)` = number of reads on `a`.

Two readings of "s bits of state" are possible and both are used below:

* **S1 (state is everything, including the index):** the configuration is `s`
  bits, so `T_M(a) <= 2^s`.
* **S2 (a free input head):** the configuration is `(s bits, head position)`, so
  `T_M(a) <= n·2^s`.

The paper's model ("O(log n) bits of state, reads a cell at unit cost") is S1.
Every statement below is given for S1 and the difference for S2 is noted.

---

## 2. Proved lemmas

### Lemma 1 (transcript = certificate).
Let `S(a) ⊆ [n]` be the set of *distinct* cells `M` reads on `a`. Then for every
`a'` with `a'|_S(a) = a|_S(a)`, `M` halts on `a'` after reading the same sequence
of positions and outputs the same value.

*Proof.* Induction on the read index. The memory content after `t` reads is a
function of the values read so far, which agree; hence the next position, the
halting decision and the output agree. ∎

Consequently `S(a)` is a certificate of `a` in the sense that `delta` is
constant on the agreement class, and

```
T_M(a) >= |S(a)| >= c(a) := min{ |S| : S is a certificate for a }.
```

### Lemma 2 (`n-1` cells always suffice; proved for all `n`).
Every `S ⊆ [n]` with `|S| = n-1` is a certificate for every tape, and its
agreement class is a single tape.

*Proof.* Let `c0, c2` be the counts on `a|_S` and `x` the missing cell. The three
candidate values of `#0 - #2` are `c0-c2+1` (`x=0`), `c0-c2` (`x=1`),
`c0-c2-1` (`x=2`); exactly one equals `eps`. So the completion is unique. ∎

Hence `c(a) <= n-1` for every `a`, and
`K_global := min{|S| : S certifies every tape} <= n-1`.

### Lemma 3 (`K_global = n-1`; verified, not proved in general).
No `S` with `|S| <= n-2` certifies `A_n`. Monotonicity (a superset of a
certificate is a certificate) reduces this to testing `|S| = n-2`, which was done
exhaustively for `n = 6,7,8,9,10,11` (`lb_cert4.py`): all `C(n,2)` sets fail, in
both readings (treating `⊥` as a value, and restricting to the halting
sub-family). The same script confirms Lemma 2 numerically: across all `(n-1)`
projections there are **0** classes with more than one `delta`.

So the deterministic *query* complexity of `delta` (unlimited memory, adaptive
queries) is exactly `n-1` in the "some fixed set" sense, and the earlier report's
table is reproduced by an independent implementation. This is the well-known cap:
`T >= n-1` and no argument of this shape can give more.

### Lemma 4 (state caps reads: `T <= 2^{s+O(1)}`, hence `s >= log2 T - O(1)`).
In model S1, if `M` halts on `a` then `T_M(a) <= 2^{s+O(1)}` (the `O(1)`
absorbs the finite control); in S2, `T_M(a) <= n·2^{s+O(1)}`. Equivalently every
halting machine satisfies `s >= log2 T_M(a) - O(1)`.

*Proof.* The memory content plus the finite control is the whole configuration. If
the same configuration occurs at two read indices `t < t'`, the machine's
subsequent behaviour is identical from both, so it reads the same positions
forever and never halts. So the `T` configurations on a halting run are distinct,
and there are `2^{s+O(1)}` of them. ∎

Three consequences, all of which constrain what a proof can look like:

1. With `s = O(log n)` the model can make at most **polynomially many** reads, so
   the conjecture is a statement about *which* polynomial, not a coarse
   separation. The read budget `n·polylog(n)` needs `s >= log2 n + O(log2 log2 n)`;
   the conjecture says that is still not enough, for any constant multiple of
   `log2 n`.
2. Any proof must establish `s = ω(log T)`: the state must exceed the
   information-theoretic minimum by a super-constant factor. Nothing in the
   deterministic read-only toolbox produces such a statement.
3. Corollary: `#outputs(M) <= 2^s`, and `delta` has `n+1` values, so `s >=
   log2(n+1)`. This is the only *space* lower bound available and it is tight at
   `O(log n)`.

### Lemma 5 (machine = branching program).
For fixed `n, s`, the set of deterministic `s`-bit machines with `T` reads is
exactly the set of branching programs over variables `x_1..x_n` with alphabet
`{0,1,2}`, of width `2^s` and length `T`, in which the variable queried at a node
is a function of that node (non-oblivious order) and variables may be re-read.
Size is `T·2^s`.

*Proof.* Node = memory content; at node `u` the program names the position
`i(u)`, and the three outgoing edges are labelled by the three values. Conversely
a BP of that shape is a machine whose state is the node. ∎

So the conjecture is exactly: *`delta` has no branching program of width
`n^{O(1)}` and length `n·polylog(n)`* (size `n² polylog n`). For unrestricted
branching programs no super-linear lower bound is known for any explicit
function; the `Ω(n²/log² n)` bounds in the literature are for *read-once*
programs, and re-reading is precisely what the conjecture is about. A proof would
be a breakthrough in branching-program lower bounds, not a routine extension of
the existing argument.

### Lemma 6 (why the measures cap at `n` and `log n`).
The following are the standard measures for this model, and each is bounded by a
function of `n` and `log|range|` only:

* certificate complexity: `<= n-1` (Lemma 2);
* communication-matrix rank: `rank <= |range| = n+1`, so the log-rank bound only
  forces width `>= n`, i.e. `s >= log2 n` — met exactly by `s = O(log n)`;
* subfunction entropy (Nechiporuk-type): a block `B` has at most
  `|range|^{3^{|B|}}` subfunctions, so `log2 N_B <= 3^{|B|}·log2 n` and any
  normalised sum is `O(n log n)` at best, `O(n)` in the standard normalisation.

So every measure certifies one of the two halves — `T >= n` or `s >= log2 n` —
separately. The conjecture is precisely the statement that **the two cannot be
attained simultaneously**, and no known measure expresses non-simultaneity.

---

## 3. Route (a): cutting the tape

### 3.1 The exact statement, and its ceiling
Fix a cut `k`. For a tape `a` let `x_k(a)` be the number of read indices `t` such
that reads `t` and `t+1` are on opposite sides of `k` (a *crossing*). If the
machine's reads `1..t` are on Alice's side `[0,k)` and read `t+1` is on Bob's
side, Alice can send her `s`-bit content; Bob continues and sends back when the
machine returns. This gives a two-player deterministic protocol with
`x_k(a)+1` messages of `s` bits computing `delta`, so

```
(x_k(a) + 1)·s  >=  CC_k(delta),      T_M(a) >= x_k(a) >= CC_k(delta)/s - 1,     (1)
```

where `CC_k` is deterministic two-way communication complexity for the partition
at `k`.

**Ceiling theorem.** `CC_k(delta) <= min(k, n-k)·log2 3` for every `k` (one player
sends his half). Hence (1) can never give more than

```
T >= min(k,n-k)·log2 3 / s - 1   <=   0.79 n / log2 n      (k = n/2, s = log2 n).
```

That is *smaller than the certificate bound* `n-1` for every `n >= 4`, so route
(a) is dominated by the argument the paper already has. No refinement of the
protocol can beat this, because the ceiling is the length of the input: to reach
`T = n polylog n` with `s = log n` one would need
`CC_k(delta) >= s·n·polylog(n) = n log n polylog(n)` bits, i.e. `n log³ n` for
`polylog = log²`, on an input of `n log2 3` bits. More communication than the
input carries is impossible for any protocol, let alone this one.

### 3.2 The measurement
`CC_k` for `k = n/2` is `ceil(log2 #distinct rows)` of the `delta` matrix, which
can be computed exactly at small `n` (`lb_cert2.py`):

| n | left halves | distinct rows | one-way CC (bits) | max `(n/2)log2 3` | bound (1) with `s=log2 n` |
|---|---|---|---|---|---|
| 6 | 27 | 17 | 4.09 | 4.75 | 1.58 reads |
| 7 | 27 | 18 | 4.17 | 4.75 | 1.49 |
| 8 | 81 | 50 | 5.64 | 6.34 | 1.88 |
| 9 | 81 | 51 | 5.67 | 6.34 | 1.79 |

The one-way value is `≈ 0.87 · (n/2)log2 3`, essentially the maximum; it is an
*upper* bound on the two-way complexity that (1) actually needs (`CC_2way <=
CC_1way`), so the true read bound from (1) is at most these numbers — under two
reads. Four orders of magnitude below the
truth at `n = 1024`, where the known construction costs `0.62 n² ≈ 6·10^5` reads.

### 3.3 Why the suggested refinement makes it worse
Two readings of "condition on the window contents at the crossing":

* If Alice also sends the part of the window on her side, each message grows from
  `s` bits to `Θ(m)` cells (up to `n/2` cells, not `Θ(log n)`: the window may be
  entirely on one side of the cut except for the crossing cell). The per-message
  cost rises, so the number of crossings one can force *falls*: (1) becomes
  `T >= CC_k/(s + Θ(m))`. The refinement moves the bound in the wrong direction.
* If instead one wants to prove the state at a crossing carries no information
  beyond the window position (`2 log n` bits), the conclusion is `s >= 2 log n`,
  which is *consistent* with `s = O(log n)` and therefore proves nothing. To get
  a contradiction one would have to show the state must encode the window
  *contents* (`Θ(m)` bits), and that is simply false: the machine can re-read
  them. Re-reading is exactly the phenomenon the bound is supposed to explain.

### 3.4 The structural failure
Crossings are chosen by the machine. A machine may read all of Alice's half and
then all of Bob's half: `x_k = 1` on every tape, and (1) degenerates to
s `>= CC_k`, which is a statement about **one-pass** machines only. To bound reads
one needs a lower bound on the *number* of crossings, and no communication
argument supplies one, because the machine is free to sweep. This is not a
technicality: it is the reason route (a) cannot see repeats at all, and repeats
are the whole question (Lemma 7).

---

## 4. Route (b): counting behaviours

### 4.1 The stated comparison is invalid
"Count how many distinct `(window, next-action)` sequences a machine with `s`
bits can distinguish and compare with the number of distinct behaviours of the
rule": the machine is not required to reproduce the trajectory, only the final
position. There are `n+1` possible outputs (`⊥` included), so the memory need
only carry `log2(n+1)` bits, which `s = O(log n)` does with room to spare. If one
compared `s` with `log2(#behaviours) = Θ(n²)` one would "prove" `s = Ω(n²)`,
which is false: a machine with `n` bits and `n-1` reads computes `delta`
(Lemma 2). The comparison is between quantities that are not required to match.

### 4.2 The sharpened version is exactly the certificate bound
The only correct counting argument is over *transcripts*: the machine's
execution on `a` is a path of `<= T` queries whose leaf is the output, i.e. a
decision tree of depth `T` with unlimited memory. Conversely every such tree is a
machine. So the transcript count is **exactly** the adaptive decision-tree depth
of `delta`, which by the certificate table is `max_a c(a)`, equal to `n-2` or
`n-1` at the sizes tested. Since a depth-`(n-1)` tree exists (Lemma 2), route (b)
made rigorous proves `T >= n-1` and cannot prove more.

### 4.3 The other counting argument proves the wrong direction
Counting *states* gives Lemma 4: `T <= 2^s`, i.e. `s >= log2 T`. That is a lower
bound on space in terms of reads, not on reads in terms of space, and at
`T = n polylog n` it evaluates to `log2 n + O(log2 log2 n)` bits — satisfied
by `s = O(log n)`. There is no way to run it backwards.

---

## 5. Route (c): the adversary

**Ceiling theorem.** Let an adversary answer the machine's reads so as to keep the
family `F` of tapes consistent with the answers containing two tapes of different
output. Its state is a partial assignment on a set `S` of cells. By Lemma 2, once
`|S| = n-1` the family is a *singleton*, so `F` cannot contain two different
outputs; hence the adversary is defeated after at most `n-1` distinct cells.
Repeats of an already-answered cell do not change `F` at all, so an adversary of
this form cannot force a single repeat. Its conclusion is `T >= n-1`, exactly the
certificate bound, and no strengthening of the ambiguity invariant helps, because
the invariant is a function of the *set* of fixed cells and that set saturates.

To break the cap the adversary would have to track what the machine has
*forgotten* — i.e. to argue "the memory at time `t` no longer contains the value
of cell `i`, so the machine must read `i` again". That is a statement about
memory contents, and proving it is Lemma 5's branching-program problem. Note that
such a statement is false in the form usually attempted: a machine is free to
choose a read schedule that never needs a forgotten value, and the only thing
that stops it is the *computation*, not the schedule.

---

## 6. What a proof would have to show

### Lemma 7 (the missing ingredient is repeats).
Write `T = |S(a)| + R(a)` with `R(a)` the number of repeated reads. `|S(a)| <= n`
always. Hence

```
T >= n·polylog(n)   <=>   R(a) >= n·polylog(n) - n  on some tape,
```

i.e. an average of `polylog(n)` reads per cell. Every certificate, communication,
rank, or counting argument bounds `|S|` (or information about it) and is
therefore blind to `R`. **The conjecture is exactly the statement that an
`O(log n)`-bit machine must re-derive cell values `ω(polylog n)` times.**

Three equivalent formulations, all open:

1. **Branching programs:** no width-`n^{O(1)}`, length-`n polylog n` BP computes
   `delta` (Lemma 5).
2. **Fine-grained, inside L:** random access at unit cost is at least as strong
   as a log-space Turing machine with a two-way read-only input (one random
   access costs a walk of length `<= n`), so the conjecture implies
   `delta ∉ DTISP(n² polylog n, O(log n))`. Fine-grained time lower bounds of
   this kind inside L are not known for any natural function.
3. **Tradeoff curve shape:** the conjecture says `T(s)` is a step at `s = Θ(n)`
   (`n-1` reads once the tape fits in memory, `Θ(n²)` reads for every
   `s = O(log n)`), rather than the smooth curve `T ≈ 2^s` that Lemma 4 allows.

### 6.1 The quantitative reason the natural simulators sit at `n²`
For all four rules the trajectory has `Θ(n²)` steps, organised as `Θ(n)` epochs
each re-scanning `Θ(n)` cells (§7). The trajectory therefore asks `Θ(n)`
"largest-concavity"-type questions, each of certificate complexity `Θ(n)`
(the `b2` answer at window `[L,r]` compares two block averages whose cells are up
to `r-L` apart and whose sign pattern in `m` is non-monotone, so no binary search
exists — `REPORT.md` §1). With `s = O(log n)` there is no room to cache the prefix
sums that would answer them, so each epoch is re-derived from the tape. This is
the mechanism, and it is not a proof: it says what an `n²`-read machine does, not
that every machine must.

### 6.2 Why the per-step cost is *not* the obstruction
Each step changes the window by one cell, so a machine that keeps the window sum,
the half-window sum and any "count of adjacent pairs / occurrences" statistic
answers both bits in `O(1)` reads per step; with the region sums also maintained
across cuts, the epoch-initial query (window = region) is free. So for the two
rules whose `b1` is such a statistic (`runs`, `pattern`) there is an
`O(n²)`-read, `O(log n)`-bit simulator at `~2-3` reads per step instead of `5.7`
— a constant-factor improvement over `0.62 n²`, no asymptotic gain. (Design
argument, not implemented here.) The `n²` is in the *number* of steps. Any
`n polylog n` machine must therefore answer fewer than one read per step once
`n ≳ 1800` (solve `log²n/(0.065n) < 1`), i.e. it must **skip blocks of steps and
predict their outcome**, and
the natural such prediction (the "horizon jump" of `REPORT.md` §4.1) is wrong on
about half the tapes with a mismatch rate that does not fall with `n`.

---

## 7. What was measured here

`lb_steps.py`, `lb_epochs.py`, `lb_cert4.py`; ensemble and walk exactly as
`_arr_rules.walk` (verified: step counts reproduce the paper's table, e.g.
`steps/n² = 0.0482` for `longrun` at `n = 128` against `0.0482` in `REPORT.md`
§3). 30–62 tapes per size, seed 2026.

**Step counts (`steps/n²`)** — stable, so no rule is trivially cheap:

| rule | n=64 | n=128 | n=256 | max/`n²` at 256 |
|---|---|---|---|---|
| longrun | .0682 | .0482 | .0645 | .210 |
| runs | .0421 | .0309 | .0498 | .160 |
| pattern | .0573 | .0427 | .0544 | .178 |
| palindrome | .0581 | .0453 | .0641 | .210 |

**Epoch structure** (`n = 256`): epochs `= 0.51–0.65 n`, mean scan length
`= 0.10–0.13 n`; so `Θ(n)` questions of `Θ(n)` read-complexity each.

**Defined fraction** (`delta ≠ ⊥`): 33, 81, 38, 76, 38, 71 % at
`n = 6,7,8,9,10,11` — the function is partial, oscillating with parity. All
certificate statements were checked both ways (with `⊥` as a value and on the
halting sub-family) and are unchanged.

**Certificates:** `K_global = n-1` for `n = 6..11`; every `(n-1)`-projection class
is a single tape, as Lemma 2 predicts.

**One-way communication across the middle cut:** 4.09 / 4.17 / 5.64 / 5.67 bits at
`n = 6,7,8,9`, against a maximum of `(n/2)log2 3 =` 4.75 / 4.75 / 6.34 / 6.34.

---

## 8. Proved versus conjectured

**Proved (unconditionally, in the model of §1):**

1. Lemma 1 — the read set is a certificate; `T >= c(a)`.
2. Lemma 2 — every `(n-1)`-set is a certificate for every tape; `c(a) <= n-1` and
   `K_global <= n-1` for all `n`.
3. Lemma 4 — `T <= 2^s` (S1) or `n·2^s` (S2); equivalently `s >= log2 T`.
4. Lemma 5 — machine `≡` branching program of width `2^s`, length `T`.
5. Route (a) ceiling — `T >= CC_k/s - 1` and `CC_k <= min(k,n-k)log2 3`, so no
   crossing argument can exceed `0.79 n/log2 n`; also the crossing count is the
   machine's choice.
6. Route (b) — the behaviour-counting comparison is invalid; the transcript
   version is equivalent to "decision tree of depth `T`", hence to `T >= n-1`.
7. Route (c) — an adversary that keeps two outputs alive is defeated after `n-1`
   distinct cells; repeats cannot be forced.
8. Lemma 7 — a proof of `T >= n polylog n` is exactly a proof of
   `R(a) >= n polylog n - n` repeat reads, i.e. `s = ω(log T)`.

**Empirical (not proved, reproduced here or in `REPORT.md`):** `K_global = n-1` for
`n <= 13`; `Θ(n²)` steps and `Θ(n)` epochs for all four rules; the one-way
complexity of §3.2; the `52 %` mismatch rate of the horizon jump.

**Conjectured:** that no `O(log n)`-bit machine computes `delta` in
`n·polylog(n)` reads, for any of the four rules. Nothing found here contradicts
it, and §6.2 explains why every natural simulator is stuck at `Θ(n²)`. But the
statement is equivalent to a branching-program lower bound at width `n^{O(1)}`
(or to a fine-grained lower bound inside L), and no current technique reaches
either. It should be reported as open, not as "we could not find the argument".

---

## 9. Consequences for §8 of `bthm_merged.tex`

1. The phrase "the obstruction is the machine's memory and not the data's
   redundancy" is right but incomplete: by Lemma 4 the memory *caps* the reads at
   `2^s`, so the conjecture is a statement about the sharp constant in a
   polynomial, and by Lemma 6 every available measure certifies `T >= n` or
   `s >= log2 n` separately. The honest sentence is: *the classification is a
   space-versus-reads tradeoff, and no technique in the deterministic read-only
   model proves one for a function of range `n`.*
2. "The two standard routes fail on this function" should become "the three
   standard routes fail, each with a provable ceiling: `n-1` for
   certificate/adversary/transcript counting, `0.79 n/log2 n` for
   crossing/one-way communication, and `s >= log2 n` for state counting".
3. The sentence "no set of `n-2` cells determines the output" can be restored —
   it is reproduced here for `n = 6..11` by an independent implementation, and the
   `(n-1)`-sufficiency proof is the count argument of Lemma 2. The earlier
   instability is very likely a **projection-direction bug**: if the code groups
   tapes by the values *outside* the read set instead of *on* it, it tests
   `2`-cell sets and reports a much smaller certificate size (my first attempt
   had exactly this bug and produced the same kind of void numbers, see
   `lb_cert.py` / `lb_cert3.py`). The check to run before trusting any such
   table: every `(n-1)`-projection class must be a singleton, which fails
   immediately under the wrong direction.
4. Do not claim `T >= n polylog n` is "left open" without saying that proving it
   requires `s = ω(log T)` and that this is a branching-program lower bound.
5. If a theorem is wanted rather than a conjecture, restrict the model — e.g. to
   machines that read only inside the current window and determine every step.
   There `T >= #steps = Θ(n²)` is immediate from the step counts (which are
   themselves measured, `REPORT.md` §3).

---

## 10. Files

```
probe_ologn/lb_steps.py     four rules: steps/n^2, steps/(n log n), max/n^2
probe_ologn/lb_epochs.py    epochs, mean scan length, steps
probe_ologn/lb_cert2.py     one-way CC across the middle cut (rows of the delta matrix)
probe_ologn/lb_cert4.py     K_global via all (n-2)-subsets; (n-1)-class check
probe_ologn/lb_cert.py      first attempt, projection bug (kept as a record)
probe_ologn/lb_cert3.py     partial-function variant, same bug
```

Run with `$P = C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`:

```
$P lb_steps.py 64 128 256
$P lb_epochs.py 128 256
$P lb_cert2.py 6 7 8 9
$P lb_cert4.py 6 7 8 9 10 11
```

`lb_cert.py` and `lb_cert3.py` project onto the *complement* of the read set; the
numbers they print ("0 of 15 certify") are about 2-cell sets and are void. They
are kept only so the mistake is on record.
