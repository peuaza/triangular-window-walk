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

## 10. The top rung is NOT filled: the perfect-square rule fails, and why

*Corrected 2026-11-06 from a subagent report that checked my numbers. The first version of this
section claimed $g\propto n^{2.93}$ for this rule. Both the exponent and the certificate
argument behind it were wrong.*

**The rule.** Binary model; read the window as an $m$-bit integer $N$; $b_1=[N	ext{ is a
perfect square}]$; $b_2$ from `KCHG`.

**The argument I made.** Certifying that $N$ is *not* a square requires determining $N$, since
any unread bit can be flipped and can turn $N$ into a square. Hence a certificate is all $m$
bits.

**Why it is false.** To certify that $N$ is not a square it suffices to exhibit a modulus $q$
with $N \bmod q$ a quadratic non-residue modulo some prime — an $O(1)$-**bit** witness. The
unread-bit argument only shows that the certificate must *determine the answer*, not that it must
determine $N$.

**The machine that exploits it.** Carry $N \bmod Q$, $Q$ the primorial of the first $k$ primes,
in $O(1)$ words. Then each query costs 1 read per step plus 2 per cut, and a full scan only when
the residue is a quadratic residue. Using $8\log_2 n$ bits — the only budget legal at every $n$:

| $c$ | bits at $n=4096$ | reads/$n^2$ at 256…4096 | exponent |
|---|---|---|---|
| 4 | 47 | 0.102 0.135 0.116 0.146 0.138 | 2.10 ± 0.06 |
| **8** | **94** | **0.061 0.057 0.049 0.059 0.055** | **1.97 ± 0.05 (flat)** |

**So $g=\Theta(n^2)$: the square rule is a natural rule at $n^2$, not at $n^3$.** The
$O(\log n)$-bit budget caps $Q$ at $\mathrm{poly}(n)$, so the false-positive rate is a
*constant*; getting it below $1/n$ would need $\log n\cdot\log\log n$ bits, still out of class,
whence only $g\le n^{3-o(1)}$.

**And a measurement error of mine, corrected.** The first version reported
reads/$n^2=2.34,4.08,6.66,15.07,30.09$ (slope $0.925$). Those are **not** the machine's reads:
they sit within 1–2 % of $1.5	imes r_{	ext{low}}$, the 16-bit residue prefilter, and the scan
column is $\sum m/n^2=6.81,13.85,24.63,58.37,118.41$ (slope $1.032$). Summing the prefilter
counters $r_{	ext{low}}+r_{	ext{mod}}+r_{	ext{mod2}}+r_{	ext{mod3}}$ is not a read count.
Correct numbers: `SQUARE_RULE.md` §3.1.

**What this changes.** The ladder's top rung has **no rule instance**. $n^2$: many (this rule
included); $n^{2.5}$: `NBLK>2/>3`, proved; $n^3$: none known.

**The pattern, now five times over.** A statistic that looks like it needs $\Theta(m)$ reads
yields to a *sketch* — a modular residue here, a prefix count there, a level band in `NBLK`.
Each time the sketch costs $O(1)$ or $\Theta(\sqrt n)$ words and collapses the read count. The
top rung may be empty: we know of no rule that forces $\Theta(m)$ reads per query.

### 10.1 The walk verification is valid

`probe_square/sqverify.py 10`, exit 0: the walk on the whole cube $n=3\ldots10$ (2 040 tapes)
against three independent engines, 0 mismatches; 94 192 rule evaluations, 0; 1 908 support
tapes, 0; 4 000 residue-table values and 4 000 gate-soundness values, 0. So the *walk* was
verified; what was wrong was the cost column I read out of it.

## 11. RETRACTED: the phase family does not fill (n, n^2)

*Retracted 2026-11-06 13:00, after a dedicated check (`N15_RIGOR.md`). This section claimed a
continuously tunable worst-case exponent $1+\alpha$. **That claim is false: the worst-case cost
is $\Theta(n)$ for every $\alpha$.***

**What is true.** The *step* count is $\Theta(n^{1+\alpha})$, and that continuum survives even
when the rule is restricted to the window and its length (steps 1.22 / 1.49 / 1.69 at
$\alpha=.25/.5/.75$). So the step-count remark of §9 stands.

**What is false.** The cost does not follow the step count. The phase rule's per-step bit is
always a single **run** of one bit: $b_1=1$ moves $\ell$ right, $b_1=0$ moves $r$ left, and once
a 0 is met every remaining query in that epoch is the same cell. A machine holding three
integers $(L,R,	ext{run length})$ scans each cell once and reproduces the output exactly:

* verified against the engine on 81 840 walks ($n\le12$, rules A/C, five $\alpha$) and
  2 760 large-$n$ walks, **0 mismatches**;
* cost $\le1.75\,n$ on every tape $n\le14$, $1.00\,n$ at $n=65\,536$, **independent of $\alpha$**;
* at $n=1024,\alpha=.5$: **1 027 reads against the claimed 11 245**; at $n=65\,536,\alpha=1$:
  65 539 against $1.07	imes10^9$ — a factor 16 384.

**So the phase family is fully linear, and $(n,n^2)$ stays empty.**

### 11.1 Three defects in the original measurement

1. The "reads" column was a **4-word FIFO cache cost** — 26 % above Belady's optimum, 11× above
   the 3-word machine — on a walk that touches only **17 distinct cells** at $n=1024$. It
   measured a machine, not $g$.
2. The "pure phase $\Rightarrow$ cut" law is **false**: 10 170 of 122 686 pure phases do not
   cut. The correct condition is `run1 >= k+1`.
3. Rule A **never meets** on any tape (0 outputs in 24 528 walks), so its $g$ was vacuous by the
   project's own test.

And the $\alpha$-ladder had its own defect: $k(M)\equiv1$ for $M<2^{1/\alpha}$, so $\alpha=0.1$
and $\alpha=0$ are the same rule up to $n=10^6$. Part of the "nine-point continuum" was one rule
measured nine times.

### 11.2 The pattern, now six times over

A statistic that looks like it needs $\Theta(m)$ reads yields to a **small sufficient summary**:
a level band (`NBLK`), a modular residue (square rule), a run length (here). Every collapse was
found by asking *what the machine actually has to remember*, not by fitting a curve.

### 11.3 What stands after the retraction

| rung | instance | status |
|---|---|---|
| $\Theta(n)$ | canonical, monotone family, `pattern`, `pal`, **and now the phase family** | measured |
| $\Theta(n^2)$ | 14 natural rules | measured |
| $\Theta(n^{2.5})$ | `NBLK>2/>3` | proved |
| $\Theta(n^3)$ | none | no rule known |

**$(n,n^2)$ is empty across every rule we have. Conjecture 2's lower half
(non-fully-linear $\Rightarrow$ $g=\Omega(n^2)$) is therefore *not* refuted — it has no
counterexample, and its only would-be counterexample collapsed to linear.**
