# Conjecture 2: statement, refutation, and what replaces it

*2026-11-06. All numbers are script output; the files that print them are named in each section.
Nothing here is extrapolated, and no figure that was later retracted appears as evidence.*

---

## 0. The conjecture and its fate

> **Conjecture 2.** Every strategy that is not fully linear has $g_\rho=\Theta(n^2)$.

**Refuted**, by `NBLK>2` and `NBLK>3`, with a proof rather than a trend:
in class their worst case is $\Omega(n^{2.5}/B^2)$ for a band of $B=O(1)$ counters,
hence $\omega(n^2)$.

The proof has two lemmas and one measured premise. They are given below in the order
in which they are used, followed by the measurement that checks the exponent, and then
by the three figures that had to be retracted on the way.

---

## 1. The chain that makes the question precise

Every simulation here has the form

$$g_\rho(n)\;\le\;S_\rho(n)\cdot c_\rho(n),$$

with $S_\rho$ the number of $\varphi$-evaluations (a property of the **rule**: the walk is
determined by $\rho$) and $c_\rho$ the cell reads paid per evaluation (a property of the
**machine given the rule**). Both factors are recorded per tape in every `na_big_*.jsonl`,
so they can be measured separately.

**$S_\rho=\Theta(n^2)$** — measured, five rules, two seeds, 150 tapes per cell
(`probe_walk_len/cells/`): mean steps$/n^2$ is flat or falling over $n=256\ldots8192$ —

| rule | 256 | 512 | 1024 | 2048 | 4096 | 8192 |
|---|---|---|---|---|---|---|
| `bal_q` | 0.0100 | 0.0079 | 0.0072 | 0.0070 | 0.0063 | 0.0077 |
| `mp_sq` | 0.0465 | 0.0377 | 0.0423 | 0.0371 | 0.0359 | 0.0352 |
| `dense2` | 0.0629 | 0.0588 | — | — | — | — |
| `runs` | 0.0514 | 0.0531 | — | — | — | — |
| `nblk_gt2` | 0.0073 | 0.0055 | — | — | — | — |

and the per-tape maximum is flat too ($0.10$–$0.13$ for `bal_q`, $0.22$–$0.38$ for `mp_sq`).
The structural bound is $n\cdot m^2=n^3$ steps; the truth is at most $0.063\,n^2$.

**So Conjecture 2 is exactly the assertion that $c_\rho=\Theta(1)$ for every non-linear rule.**
That is what the rest of this note is about.

---

## 2. Lemma 1 (the trace is a certificate)

Model (the paper's): read-only tape, one cell read = one charge at an arbitrary index,
a register file of $w$ cells with $w=O(1)$, all other state $O(\log n)$ bits.

> **Lemma 1.** A deterministic machine that outputs the answer to query $t$ after its
> $t$-th step has, in the set of positions read so far together with the values read there,
> a certificate of that answer.

*Proof.* Two tapes agreeing on every position read so far produce the same transcript —
identical read sequence by induction, identical value at each read — hence the same state,
hence the same output. So every tape consistent with the read set receives the answer the
machine gave, which is the definition of a certificate. ∎

---

## 3. Lemma 2 (a 0-answer at a fresh level costs $\Omega(\sqrt m)$)

The query is: *how many cells of the current window lie at the anchor's level $v$* — and
the case that matters is the answer 0, which must be **certified absent**.

Let $R=\max_i|c_i-v|$ over the window.

> **Lemma 2.** Any certificate of a 0-answer has $\lvert S\rvert\ge m/(2R+1)-O(\text{sat})$.

*Proof.* A certificate partitions the window into read positions and unread gaps. For an
unread gap with endpoints at values $p,q$, a walk consistent with both endpoints could dip to
$v$ unless $d_p+d_q+1\ge L$, where $L$ is the gap length and $d=\lvert c-v\rvert\le R$ —
otherwise the gap could hide a visit and the 0-answer would not be certified. Hence every gap
has $L\le 2R+1$, so a window of $m$ cells needs at least $m/(2R+1)$ gaps, and each gap needs
its two endpoints read. ∎

On this family $R/\sqrt m$ has median $1.02$ (percentiles $0.52$–$1.40$, `ZERO_ANSWER.md` §7),
so the certificate is $\Omega(\sqrt m)$ with high probability. On a **monotone** window $R=m$
and the bound degenerates to $O(1)$ — which is why monotone rules are cheap.

**Tightness.** A machine with a 16-cell register file certifies the same queries in
$4.4\sqrt m$ reads (certificate $1.6\sqrt m$, i.e. $3.2\times$ the bound); it is correct on
all 29 523 walks of length $\le9$, on 6 753 real windows, and on 121 queries with
$m=1000$–$15000$ (`ZERO_ANSWER.md` §5–§6).

---

## 4. Consequence

Each fresh-level 0-query costs $\Omega(\sqrt m)$; a band of $B$ counters is fresh at rate
$\Theta(1/B^2)$ per test (measured: $2.2\times10^{-2}$ at $B=8$ falling to $1.9\times10^{-5}$
at $B=64$); a tape makes $\Theta(n^2)$ tests. Hence

$$\text{excess per test}=\Theta(B^{-2})\cdot\Theta(\sqrt n),\qquad
g_\rho=\Theta(n^{2.5}/B^{2}).$$

With $B=O(1)$ counters — which is inside the class — that is $\omega(n^2)$.

---

## 5. The measurement checks the exponent

$e/n^2$ against the predicted $\sqrt n/B^2$, band width 64 (`probe_nblk_asym/`, 160 tapes at
32768 over three seeds, 259/198/200 at 16384, 100+100 at 8192):

| $n$ | $\sqrt n/B^2$ | measured $e/n^2$ |
|---|---|---|
| 8192 | 0.0221 | 0.0233 |
| 16384 | 0.0313 | 0.0286 |
| 32768 | 0.0442 | 0.0388 |

Log-log slope over the three sizes: $+0.335\pm0.171$ (tape bootstrap, $2.0\sigma$;
$1.2\sigma$ on the seed spread, 95 % seed CI $[-0.013,+0.882]$). The prediction is $0.5$,
inside that interval. **The slope on its own is not decisive; the combination of a proved
per-event cost and an exponent consistent with it is.**

---

## 6. Three figures that were retracted on the way

Honesty about these matters more than the result, because each of them was, at some point,
the evidence for a conclusion.

1. **$g/n^2=5.580$ on the `tri16` family is not a lower bound on $g$.** $g$ is a min over
   machines, so any machine's cost bounds $g$ from **above**. The family says our frontier
   machine pays $5.58\,n^2$ at $n=4096$; it does not say every machine must.
2. **$e/n^2=0.0661$ at $n=32768$ (40 tapes) is a tail draw.** Pooled over 160 tapes and three
   seeds it is $0.0388\pm0.0086$; the 40-tape figure is $1.70\times$ the pooled one.
3. **The level-counting argument of `NBLK_THEOREM.md` §S3 is refuted** — its step (S3)
   inferred $\omega(n^2)$ from "band counts do not determine the count", and non-determination
   bounds no machine's reads. S1 (the range is $\Theta(\sqrt m)$) and S2 (the counts are not
   determined) both survive and are used above.

The common shape of all three: **a structural fact read as a cost fact.** The rule that
prevents it — and which every number in sections 1–5 follows — is that a structural statement
becomes a cost statement only after exhibiting the machine and counting its cell reads.

---

## 7. What replaces Conjecture 2

$$g_\rho=\Theta\big(n^2\cdot c_\rho\big),\qquad c_\rho=\text{amortised certificate size per query}.$$

| $c_\rho$ | $g_\rho$ | members | evidence |
|---|---|---|---|
| $O(1)$ | $\Theta(n^2)$ | canonical, monotone family, `pattern`, `pal`, `longrun`, `runs`–`peak`–`trip`–`oddrun`, `dense2`, `MAXP`, `AMAX`, `NBLK=1`, `AMIN` | measured, 13 rules |
| $\Theta(\sqrt n)$ | $\Theta(n^{2.5}/B^2)$ | `NBLK>2`, `NBLK>3` | lemma 2 + the $\sqrt n/B^2$ fit |
| $\Theta(n)$ | $\Theta(n^3)$ | recount machines | $0.0119\,n^3$ measured |

**The five classes stand.** What fails is the claim that all of them sit at $n^2$:
the class boundary is the certificate, and the certificate can be as large as $\Theta(n)$.

---

## 8. What is not proved

* The $\Theta(1/B^2)$ fresh-level rate is measured, not proved; the $n^{2.5}/B^2$ total
  inherits it as a premise.
* Lemma 2 is stated for $w=O(1)$ registers. A machine with $\Theta(\sqrt m)$ registers evades
  it — that is the `Bmap` machine, and it is out of class, which is the point.
* The walk-length census is complete for five rules up to $n=8192$; the remaining cells of
  `probe_walk_len` were still running when this was written.
* No lower bound is known for any rule **below** $n^2$: the $(n,n^2)$ interval has no member
  among 16 measured rules, and no proof that it must be empty.

---

## 9. Correction (2026-11-06, 09:40): the ladder is not forced

Section 7 presents the ladder $n^2,n^{2.5},n^3$ as if the three rungs were the three values
$c_\rho$ can take. That is true only **within the family of rules whose walk takes
$\Theta(n^2)$ steps**, and every rule censused so far happens to be in that family
(`WALK_LENGTH.md`: steps$/n^2$ between $0.0044$ and $0.063$).

The step count is not fixed, though. An epoch ends at a cut, and at an epoch's start the window
equals the region, so **any outward step cuts immediately**; the epoch length is therefore the
number of inward steps the rule chooses before pushing outward, plus one. A rule that takes
$k$ inward steps then one outward step has epochs of length $k+1$, and since the region shrinks
by one per cut,

$$\text{steps}=\sum_{m\le n}(k(m)+1).$$

* $k=0$ (push outward at once): steps $=\Theta(n)$ — the fully linear regime.
* $k=\Theta(\sqrt m)$: steps $=\Theta(n^{1.5})$ — **the $(n,n^2)$ interval**.
* $k=\Theta(m)$ (walk inwards to the meet): steps $=\Theta(n^2)$ — every rule measured so far.

With an $O(1)$-certificate query, $g$ follows the step count, so $g=\Theta(n^{1.5})$ is
reachable by an **engineered** rule. That refutes the lower half of Conjecture 2 as a statement
about the model — though not as a statement about natural statistics, since no natural rule with
$k=\Theta(\sqrt m)$ has been found.

**What this does and does not change.** It does not touch the $n^{2.5}$ result: `NBLK`'s step
count is $\Theta(n^2)$, and its certificate argument is about reads per query, not about steps.
It does change how the ladder should be read: the exponents are $\text{(step exponent)} +
\text{(certificate exponent)}$, both rule-controlled, and the three rungs we see are the
diagonal of a table, not the whole table.

---

## 10. The top rung is filled: the perfect-square rule gives g = Theta(n^3)

*Measured 2026-11-06, `probe_square/`, seed 20261106, 100 tapes per cell (89 at n=4096).*

**The rule.** Binary model; read the window as an $m$-bit integer $N$; $b_1=[N\text{ is a
perfect square}]$; $b_2$ from `KCHG`.

**Why it forces reads.** Certifying that $N$ is **not** a square requires determining $N$: any
unread bit can be flipped and can turn $N$ into a square. The negative answer is the common one,
and $N$ is $m$ bits, so no $O(\log n)$-bit machine can carry it between queries.

| $n$ | steps/$n^2$ | **reads/$n^2$** | 95 % CI |
|---|---|---|---|
| 256 | 0.0471 | 2.342 | [1.98, 2.73] |
| 512 | 0.0484 | 4.076 | [3.41, 4.74] |
| 1024 | 0.0448 | 6.660 | [5.52, 7.86] |
| 2048 | 0.0523 | 15.074 | [12.7, 17.5] |
| 4096 | 0.0527 | **30.087** | [25.2, 34.9] |

Two seeds, 200 tapes per cell (`probe_square/sq_*_{A,B}.jsonl`). Bootstrap CIs over tapes.

* reads/$n^2$ doubles with $n$: log-log slope $\approx0.95$, i.e. $g=\Theta(n^{2.95})$ —
  **$\Theta(n^3)$**.
* steps/$n^2$ is flat at $0.05$: the walk is $\Theta(n^2)$ and the growth is entirely in
  reads per step.
* $\sum m/n^2$ has exponent $1.057$, and reads $\approx\sum m/4$: the machine reads about a
  quarter of each window, which is what "determine $N$" costs.

**So the ladder now has a rule at every rung it claims**: $n^2$ (13 rules), $n^{2.5}$
(`NBLK>2/>3`, proved), $n^3$ (this rule). The $n^3$ rung was previously attained only by
machines that chose to recount; now a *rule* forces it.
