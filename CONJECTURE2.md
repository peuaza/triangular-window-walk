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

## 4. Consequence — **DOWNGRADED 2026-11-06 14:20: the total does not follow**

> **Audit result.** The per-query bound of §3 is a theorem, but it applies only to queries the
> machine has **not** already certified. On the real walk the band certifies almost all of them:
> 3, 3 and 6 uncertified zero-answers per tape at $n=2048,4096,8192$. Summing the honest
> per-query bounds gives a total whose log-log slope against $n$ is **1.051**, where the claim
> needs 2.5 (`probe_n25_audit/n25_scale.py`).
>
> **So $g=\Omega(n^{2.5})$ is NOT established.** What stands is the measurement of §5 - reads
> per test rising $5.08\to6.21\to10.40$ over $n=8192,16384,32768$ - which is a $2\sigma$ trend
> with no lower bound behind it.
>
> The derivation below is kept only to record what was attempted; do not cite it as a proof.

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

### 11.4 Two independent lines agree, and A is vacuous

The verification line (`N15_RIGOR.md`) and the pushdown line (`N15_PUSHDOWN.md`) were run
separately and agree on every number.

| | verification | pushdown |
|---|---|---|
| in-class worst case | $\Theta(n)$, independent of $\alpha$ | $\Theta(n)$, independent of $\alpha$ |
| $\alpha$-ladder cost exponent | 1.001 at all nine points | 1.001, 1.001, 1.001, 1.001, 1.001 |
| machine | three integers $(L,R,\text{run length})$ | one pass, two counters, $O(\log n)$ bits |

**Rule `A` is vacuous.** It never terminates by *meet*: the gap $r-\ell$ falls only to
$M-1-k(M)\ge1$ and the outward move raises it again, so the meet check can never fire — on any
tape, for any $\alpha$. Its output is the constant $\bot$. The reported $0.34\,n^{1.5}$ was the
read pattern of a computation no output depends on, which is the same defect `N15_RULE.md` §7d
flags for `D0`.

**Rule `C` is a first-long-run scan**: $\min\{L:\text{run}(L)\ge k(n-L)+1\}+2$, one pass, two
counters, $O(\log n)$ bits; worst case $n-1$ on the alternating tape, 33 reads on the all-ones
tape at $n=4096$ against the claimed 89 245.

**Rule `D` is "find the first 1"** and its certificate number is exactly $n-2$ (exhaustive,
$n\le14$) — so no device can beat $n-2$, and the family maximum is pinned there for every
$\alpha$.

**And the $\Theta(\sqrt n)$ register file was an artefact of caching cells**: the $k+1$ reads of
an epoch carry a single run length, so a counter to $k+1$ does the work; $W=1$ suffices.

**Rule-`A`/`D` aside, one honest residual:** `C`'s lower bound is only $\Omega(\sqrt n)$, so `C`
itself might be cheaper than linear — but the family maximum is pinned by `D` at $n-2$, so that
cannot rescue the $1+\alpha$ ladder.

---

## 12. The audit of §3-§4: what is theorem, what is gap (2026-11-06 14:45)

`N25_AUDIT.md`, 458 lines, probes in `probe_n25_audit/`.

| # | step | status |
|---|---|---|
| 1 | trace = certificate | **THEOREM**, independent of state size; $w$ never enters |
| 2 | 0-answer certificate $\ge m/(2R+1)-O(\mathrm{sat})$ | **THEOREM** (0-answers only); the printed gap condition is loose, not wrong (exact: $L\le\max(d_p,1)+\max(d_q,1)-1$); 27 990-case brute force, no counterexample |
| 3 | hence reads at that query $\ge(r+1-x)/(2R+1)-w$ | **GAP, false as printed**: 1-answers have 3-position certificates, and Lemma 1 bounds the *cumulative distinct* read set ($\le n$), not the reads at that query; needs an unproved state/re-read lemma (the project's own `FORMAL_MODEL` Remark 2.4 warns of this) |
| 4 | $R/\sqrt m=1.024$ | **MEASURED**, but it is a property of 0-answer windows (0.48–0.62 over all tests) |
| 5 | band exit rate $\Theta(1/B^2)$ | **FALSE as used**: that is the *band-exit* rate, not the *0-answer* rate |
| 6 | $\Theta(n^2)$ steps | **MEASURED, wrong direction**: mean$/n^2$ flat but median$/n^2$ falls 26× (median tape $\approx n^{1.05}$); no lower bound on steps |

### 12.1 The real failure: two anti-correlated marginals

It is **not** a mean-versus-max error — $g\ge$ mean makes averaging sound. It is that the two
rates that were multiplied are anti-correlated:

$$P(\text{0-answer}\mid\text{band exit})=1.1\times10^{-3}\ \longrightarrow\ 1.1\times10^{-4},
\qquad 0\ \text{at}\ B\ge16\ \text{in}\ 2.8\text{M tests}.$$

Replaying the project's own 121 `za_cert3` windows: 3/121 at $B=8$ and **0/121 at $B\ge16$** sit
at uncounted levels. NBLK's 0-answers occur at levels where the anchor **lingers**, which is
exactly what the band covers.

Honest per-test bound: $2.9\times10^{-4}$, $1.1\times10^{-4}$, $5.8\times10^{-5}$ at
$n=2048,4096,8192$; totals 29 / 73 / 125 reads, log-log slope **1.05**; the maximum over 18
tapes is 51.6 reads against $n^{2.5}=6.1\times10^{9}$.

### 12.2 What would repair it

$P(\text{0-answer}\mid\text{fresh level})=\Omega(1)$ — a property of the **rule**, and NBLK does
not have it. Plus the re-read lemma for step 3.

### 12.3 Consequence for the table

**No rule is now known to have a proved super-quadratic worst case.** NBLK keeps only its
measured trend (reads/test $5.08\to6.21\to10.40$ over $n=8192,16384,32768$, about $2\sigma$).

---

## 13. Three results, and one correction between them (2026-11-06 15:00)

### 13.1 The top rung is filled: `CX` reaches n^3 (`ABOVE_N25.md`)

Rule `CX`: $b_1=[\sum_i \mathrm{sg}[i]\,h(i-\ell)\bmod q<q/2]$ with kernel $h(k)=3^{k^2+1}$.
Because $h(k+1)/h(k)=3^{2k+1}$ is **not constant**, no rolling-hash update exists — which is
exactly the escape the square rule lacked.

* certificate: **0.727 / 0.757 / 0.735** of the window is single-flip sensitive at
  $n=256/512/1024$ — exponent 1, not 1/2;
* cost: $0.00226$–$0.00360\,n^3$; **slope $2.932\pm0.010$** over $n=1024\ldots8192$; two-term fit
  $0.00219\,n^3+0.519\,n^2$ (cubic 97 % at $n=8192$);
* gap stated: no lower bound against *all* $O(\log n)$-bit machines; the certificate plus a
  no-recurrence argument is the evidence.

Rule `PC` is the cautionary twin: naive machine $n^{2.991}$, but the pointer chain is a
functional graph with rho length $1.5\sqrt m$, so Brent detection answers the same queries in
$5.0\sqrt m$ reads — slope **2.493**, i.e. $n^{2.5}$, verified on 174 168 windows.

**So $n^{2.5}$ is not a ceiling.** It is the diameter exponent of the two random structures the
model can express (the range of a $\pm1$ walk, the rho length of a functional graph); `CX`
needs neither.

### 13.2 Semi-linear is not empty, and `longrun` is rescued (`SEMILIN_2.md`)

Second member: $b_1=[2\Lambda_2>m]$ ("one run covers more than half the window"): $e/n=5.55\to
6.85$ over $n=128\ldots8192$ (exponent 1.048), $g=0.299\,n^2$ steps, worst family
$0^{3n/8}2^{3n/8}1^{n/4}$. The whole $c$-family qualifies ($e/n=5.6/6.7/7.7/10.3$,
$g=0.30$–$0.41\,n^2$), and $c\to1^+$ degenerates continuously to `b0`.

**Correction to `LONGRUN_FAST.md`:** the per-epoch $4\rho$ prefix sum is a fixed-length window
sum and *is* incrementally maintainable ($\le3$ reads/epoch). Doing so drops `longrun` from
$12n$ to **$8n$** and removes the log factor. So the "longrun leaves semi-linear" conclusion of
the kappa run is **superseded**: $e=\Theta(n)$ stands, the machine was wrong.

### 13.3 What the semi-linear row now is

| member | $e$ | $g$ | status |
|---|---|---|---|
| `longrun` | $8n$ (corrected machine) | $\approx25n^2$ | semi-linear |
| $b_1=[2\Lambda_2>m]$ and its $c$-family | $5.6$–$10.3\,n$ | $0.30$–$0.41\,n^2$ | semi-linear |

Both mechanisms are the same in kind: one passing integer statistic $A$, with $b_1\equiv0$ for
$m\ge cA$. The two escapes tried (local rules, count-type certificates) fail measurably.

---

## 14. The final table after the second round of attacks (2026-11-06 15:40)

### 14.1 NBLK is dead: exponent 2.0, not 2.5

The audit's certificate bound was turned into a *machine*: band + frontier with `zac.cert_scan` on
the residual scans (725 710 scan events verified). On matched tapes at $n=16384$:

| | band | hybrid |
|---|---|---|
| $B=64$ | 4.840 | **4.547** (floor 4.492; excess cut 6.4×) |
| $B=32$ | 9.277 | **5.569** (cut 4.4×) |

Law: $e\approx a n^2+c\,n^{2.5}/B^2$ with $a=0.021$, $c\approx10^{-6}$–$10^{-5}$; measured
exponent **2.0** over $n=1024\ldots16384$; the $n^{2.5}$ term would only take over at
$n\approx10^5$–$10^7$. The certificate total is **0 events in 25 292 tests** at $B=64$. So
$g=\Omega(n^{2.5})$ is false for NBLK.

### 14.2 PC is the one member at n^2.5, and the best-supported result of the night

Endpoint-only Brent: $2.5$–$2.8\sqrt m$ (was 5.0), 0 mismatches on 195 168 windows; cost
$\mathbf{0.0176\,n^{2.5}}$, exponent $2.474\pm0.003$, tape spread 2 %. Certificate census
$1.0$–$1.1\sqrt m$ (71 % of the orbit), so the machine is within $2.6\times$ of the floor.
**Lemma S proved** (a $b$-bit state fixes at most $0.63(b+1)$ cells on half the tapes); the
$\Omega(\sqrt m)$ per-query bound then needs one hypothesis (H1, random-functional-graph
certificate size).

### 14.3 CX: n^3 measured, but the certificate argument is FALSE

The decisive control: the same action map, tape and comparison with kernel $h(k)=3^k$ (which has
a verified one-read rolling update):

| kernel | certificate | reads/query | slope |
|---|---|---|---|
| $3^{k^2+1}$ (CX) | 0.710 | $m$ | **2.932** |
| $3^k$ | **0.835** (larger) | **1** | **1.999** |

A larger certificate costing $1/m$ as much per query. **The certificate census measures chaos,
not cost** — the same error as §10, reused in §13.1. What holds CX up is **shift-rank**:
$\mathrm{ord}_q(3)=715827882$, no period $\le4096$, BM linear complexity $L/M=1.000$, Hankel
rank $=m$, no recurrence of order $A\le16$. No small summary exists (previous collapses had rank
1, 1, $\le2$; CX has $m$), but there is no lower bound either.

### 14.4 longrun stays weakly super-linear

Corrected $e/n$: 8.54, 8.53, 8.48, 8.66, 8.90, 9.33 over $n=512\ldots16384$. Constant fit
residual sd 0.298; **log fit $7.021+0.1494\log_2 n$**, $b$ at $3.3\sigma$, $F(1,4)=10.9$, AIC
−18.4 vs −12.5. Mechanism: `scanlen` $=0.405\cdot4\rho$ ($r=0.9989$), per-epoch $3+2.45\rho$,
irreducible. **The semi-linear row stays empty** — the $2\Lambda_2>m$ family also drifts
(5.55→6.85).

### 14.5 The table

| rung | rule | evidence | proof |
|---|---|---|---|
| $\Theta(n)$ | fully linear family | measured | — |
| $\Theta(n^2)$ | 14 natural rules | measured | — |
| **$\Theta(n^{2.5})$** | **`PC`** | flat cost/$n^{2.5}$, verified machine, certificate census | Lemma S + H1 |
| $\Theta(n^3)$ | `CX` | slope 2.932, no summary (rank $m$) | **none** |
| | NBLK | exponent 2.0 | — |

---

## 15. The rank: proved for CX, but it is a threshold, not a rate (2026-11-06 16:30)

### 15.1 The right definition

Rank = the **Hankel rank of the kernel**: $H_m[d,k]=h(k+d)$ over $\mathbb F_q$; by Kronecker this
is the minimal recurrence order. The rejected candidate — linear complexity of the window or of
the answer sequence — gives $L/\mathrm{len}=0.500$ for **both** CX and its one-read $3^k$ twin,
so it cannot separate them.

**CX's rank $=m$ is now proved**, not measured: $H_m=3\,\mathrm{diag}(3^{d^2})\,[9^{dk}]\,
\mathrm{diag}(3^{k^2})$, a Vandermonde matrix with distinct nodes since
$\mathrm{ord}_q(9)=357\,913\,941>m$. The earlier "recurrence crossover at $A=m/2$" is a
dimension artefact ($A>m/2$ leaves the windowed system underdetermined) and is withdrawn.

### 15.2 The implication is FALSE

"rank $\Omega(m)$ $\Rightarrow$ $\Omega(m)$ reads per query" fails:

| counterexample | rank | reads/query |
|---|---|---|
| constant comparison | $m$ | **0** |
| right-only move mix | $m$ | **1.00** (840 moves, 0 wrong) |

And the counting proof cannot work: the answer is **one bit**, so the only valid information
bound is $0.63$ reads/query — which is tight for the twin ($0.5078$ against $0.3204\,n^2$).

### 15.3 What is true: a threshold

$$\text{reads/query}=\begin{cases}1,&\mathrm{rank}\le s\ m/s,&\mathrm{rank}>s\end{cases}$$

proved and machine-verified at every step. Rank 3 and rank $m$ cost **identically** at $s=2$
(71.667 for both). **Corollary: kernel-correlation rules are $n^2$ or $n^3$ and never between.**
`PC`, the one intermediate rule, is not a kernel correlation — its exponent comes from the rho
length of a functional graph.

| type | exponent | instances |
|---|---|---|
| kernel correlation, rank $\le s$ | $n^2$ | NBLK, square rule, phase family |
| kernel correlation, rank $>s$ | $n^3$ | **CX** |
| not a kernel correlation | set by a random-structure diameter | **PC** ($n^{2.5}$) |

### 15.4 Two corrections carried forward

`CX_ATTACK.md` §4's "rank $r\Rightarrow m/r$ reads" is inverted (at $r=m$ it would predict 1);
the correct form is the piecewise rule above. And §5's "5.3×" improvement needs $A$ words
$=31A$ bits, which exceeds the $O(\log n)$ budget for $A\ge2$.

### 15.5 Open

An unconditional $\Omega(n^3)$ for all $O(\log n)$-bit machines is still open; it is carried as
hypothesis H (no $b$-bit nonlinear sketch advances $A_0$). The escape is nonlinear — a constant
tape costs one read total — which is exactly why the rank is blind to it. The $H1$ line for
`PC` was interrupted and must be redone.

---

## 16. Semi-linear has a non-degenerate member (2026-11-06 16:55)

`probe_semilin_nat/`, rule `sn_core`:

$$T(w)=[\text{no cell of }w\text{ equals }1],\qquad
E(w)=[\textstyle\sum_k \mathrm{sg}[\ell+k]\,h(k)\bmod q<q/2],\qquad b_1=T\wedge E,$$

with $h(k)=3^{k^2+1}$ (the rank-$m$ kernel of §15) and $b_2$ unchanged.

* $T$ is incremented by one cell per move, so the machine certifies $b_1=0$ with **one** read;
  $T$ fails with probability $2^{-m}$ on a random window, so the expensive branch is rare.
* **$e/n=5.31,5.63,5.36,5.44,5.53,5.36$ over $n=256\ldots8192$ — flat at $5.4$.**
* **$g=0.500\,n^2$ on the worst family** (`0/2` block tapes), log-log slope $1.995$; model
  residuals $n^2: 0.008$, $n^{2.5}: 0.668$, $n^3: 1.179$.
* Without the certificate the same rule costs $e/n=25.7\to53.4$ — rising, so the certificate is
  what makes it semi-linear.
* Verification: 86 050 walks, three machines compared pairwise, **0 mismatches**.

**Why this is not $S_n$ in disguise.** $S_n$ is "`b0` plus a rare expensive branch" and agrees
with `b0` almost everywhere. Here $b_1=T\wedge E=0$ whenever $T$ fails, whereas `b0`'s $b_1$ is
the leading bit — so the rule differs from `b0` *everywhere*, the trigger is **content-driven**
(is there a 1 in the window) rather than a geometric clock, and the expensive branch is the
rank-$m$ checksum rather than a plain scan. It is the first non-degenerate member of the row.

**So the semi-linear row is not empty, and the five-class scheme keeps its second row.** The
price is that the rule is a conjunction — natural in its parts, contrived in its shape — and
that is the honest description.

---

## 17. The continuum: four retractions, and what actually stands (2026-11-06 22:20)

### 17.1 What stands

$$\text{steps}(n)=\sum_{W}\left(\frac{W}{64}+1\right)=\frac{n^2}{128}\bigl(1+o(1)\bigr),\qquad
g_\rho=\text{steps}\times c_\rho .$$

* **The step count is a clock and is rule-independent.** Measured identical across five different
  statistics, **including a constant bit**, to residual 0.0000. No rule in the project changes it.
* **The price tracks the certificate**, not the window and not the support: a prefix-sum rule on
  the same window has support $m$ and price 1.
* **Natural certificates take three values** — $1$, $\Theta(\sqrt m)$, $\Theta(m)$ — hence the
  natural exponents are $\{2,\ 2.5,\ 3\}$.
* Two-sided bracket for any chase-like rule: $0.47\,\sigma\le$ reads/query $\le2.8\,\sigma$.

### 17.2 The four retracted versions

| version | killed by |
|---|---|
| 「a continuous spectrum」 | the chase reads only a prefix $[\ell,\ell+L)$; $p\ge L$ never occurs in 16 907 424 queries; a plain scan of $m^\sigma$ cells reproduces the exponent to 0.003 |
| 「two universality-class lines」 | $m^{2/3}$ exists twice (drift $\theta=2/3$, leaky $\alpha=1/3$) |
| 「$m^{2/3}$ gives a third class」 | $L$ takes **one value across 200 tapes** at every $m$: $\theta$ is a literal, not a computed quantity |
| 「a parameterisation: the exponent moves smoothly in $\theta$」 | at fixed $n$, sweeping $\theta$ from .20 to 1.00 in steps of .05 gives a **non-monotone** local exponent with dips of 0.15–0.20, i.e. half the whole predictive span |

**The convergence check did pass** (local exponents reach $2.497$ and $2.991$ against predictions
$2.500$ and $3.000$ at $n=16384$), but it only shows the exponent converges to *a constant that
was written into the rule*.

### 17.3 The defect this exposes

Every read index in the project is $\ell+p$ ($p<m$) or $O(1)$ boundary cells. The only
sub-window producer is the `domain(m,theta)` of the continuum rules. Measured support slopes:
drift $1.000$; $a=2/3/5$ give $0.503/0.485/0.507$; leaky $0.58$–$0.74$; folded $0.36/0.26$.

**So: no natural rule has a domain strictly between $O(1)$ and $\Theta(m)$.** The intermediate
exponents exist only in rules that declare their own domain — which is a free constant, not a
window statistic. **This is a defect of the model's permissiveness, and the paper should say so
rather than present a spectrum.**

### 17.4 Methodological record

Five wrong calls in one evening, all of the same shape — a machine's cost read as the rule's cost:
the phase family ($n^{1.5}\to n$), the square rule ($n^3\to n^2$), the memoised chase
($2.5\to1.7$ with $\Theta(n)$ state), the two-line model (killed by $m^{2/3}$), and the
parameterisation (killed by non-monotonicity). **The recurring question that catches all five:
what does the machine have to *remember*, and is that within $O(\log n)$ bits?**
