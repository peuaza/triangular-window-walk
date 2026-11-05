# The Upper Bound, and What Conjecture 2 Really Says

*2026-11-06, 01:50. Written while the n=32768 error bars are being collected.
All numbers quoted here come from `KAPPA_MEASURED.md` §46–§49 and the files named there.*

---

## 1. The chain

Every simulation in this project has the form

$$g_\rho(n)\;\le\;\underbrace{S_\rho(n)}_{\text{steps}}\times\underbrace{c_\rho(n)}_{\text{reads per step}}.$$

A *step* is one evaluation of $\varphi$ on the current window, at whatever the machine
pays for it. $S_\rho$ is a property of the **rule** (the walk is determined by $\rho$,
not by the machine). $c_\rho$ is a property of the **machine given the rule**.

So the bound factors, and the two factors can be — and now have been — measured
separately, because every `na_big_*.jsonl` row carries both `tests` and `reads`.

## 2. The step count is Theta(n^2)

| $n$ | 1024 | 2048 | 4096 | 8192 | 16384 | 32768 |
|---|---|---|---|---|---|---|
| `NBLK>2` steps/$n^2$ | 0.00579 | 0.00487 | 0.00483 | 0.00458 | 0.00461 | 0.00682 |

Flat to within a factor 1.5 over a factor 32 in $n$ (the last cell is 34 tapes, one seed).
The structural bound is $n\cdot m^2=n^3$ steps; the measurement says the truth is at most
$0.007\,n^2$, i.e. **the structural bound is loose by $n/140$**, and the conjectured
"walk length is $O(n^2)$" survives the first direct test.

**Consequence.** The route "$n^4\to n^2$" is *not* blocked by the walk length.
It is blocked — if it is blocked — entirely by $c_\rho$.

## 3. The certificate is what $c_\rho$ is

A machine that reads $k$ cells before answering sees $k$ values, so its answer is a
function of those $k$ cells: the cells read form a **certificate** for the answer.
Hence

$$c_\rho(n)\;\ge\;\text{(certificate size of one query)},$$

and conversely a machine that carries a certificate forward between queries pays little.

This is the mechanism behind every classification measured so far:

| certificate | rules | $c_\rho$ | $g_\rho$ |
|---|---|---|---|
| one integer | `longrun` ($4\rho-1$), canonical, monotone family | $O(1)$ | $\Theta(n^2)$ |
| $O(1)$ words | `pattern` ($M_{\text{final}}\le14$), `pal`, `runs` | $O(1)$ | $\Theta(n^2)$ |
| $O(\log n)$ bits | hashmod | $O(\log n)$ | $\Theta(n^2)$ |
| $\Theta(\sqrt m)$ counters | `NBLK>2/>3` via `Bmap` | $\Theta(\sqrt m)$ | $\Theta(n^{2.5})$ |
| $\Theta(m)$ | recount machines | $\Theta(m)$ | $\Theta(n^3)$ |

**The $n^3$ bound is attained**, by machines that recount per test: measured
$g=0.0119\,n^3$ for `NBLK>2`, with the predicted 48.8 and 62.1 matching the measured
48.89 and 62.12. So the general upper bound cannot be pushed below $n^3$ without
restricting the rule class.

## 4. What Conjecture 2 really asserts

> *Every non-fully-linear strategy has $g_\rho=\Theta(n^2)$.*

Given the chain and $S=\Theta(n^2)$, this is **exactly** the assertion

> *every non-fully-linear rule has $c_\rho=\Theta(1)$ amortised.*

That is what the conjecture is about, and it is why the whole measurement programme
reduces to one question per rule: **how big is the certificate?**

The four classes then read:

| class | $c_\rho$ | why |
|---|---|---|
| fully linear | $O(1)$ **and** $g=O(n)$ | every tape has an $O(1)$ certificate, including the worst |
| semilinear | $O(1)$ | $e=O(n)$ means the walk itself is short |
| near-linear | $O(1)$ | the certificate exists and is cheap to carry |
| weakly superlinear | $O(1)$ amortised, but the *predictor* cannot carry it | $\kappa\in(0,1)$ |
| strongly superlinear | grows | $\kappa=0$, and $g$ grows with it |

`NBLK>2/>3` are the one family where $c_\rho$ has been observed to grow
($6.9\to8.05\to10.39$ from $n=4096$ to $16384$ to $32768$, band width 64).

## 5. The two questions that remain

**(Q1) Does $c_\rho$ saturate for `NBLK`, or grow?** This is the whole of Conjecture 2.
$10.39$ at $n=32768$ from $6.9$ at $n=1024$ is a factor $1.5$ per quadrupling, i.e.
$c\propto n^{0.29}$; extrapolated that is unbounded, but three points with one seed
cannot separate $n^{0.29}$ from a saturating curve that has not yet saturated.

**(Q2) Is there a class machine with $c_\rho=O(1)$ for `NBLK`?** If yes, the conjecture
survives and `NBLK` is not a counterexample after all; if no, it falls. The surviving
mechanism is the **zero answer**: 82–90 % of band-epochs contain a query whose answer is
0, and certifying absence is the one thing a band of $B$ counters cannot do without a
scan. The lemma to prove or refute is: *a $B$-word machine that meets a 0-answer at the
anchor pays $\Omega(m)$ reads.*

**Both questions are being attacked tonight.** Until they are answered, the honest
statement about Conjecture 2 is: *no counterexample among 16 measured rules, one
candidate with a $2.3\sigma$ trend and a mechanism, and a proof route that was tried
and refuted (level counting — `NBLK_THEOREM.md` §S3).*

## 6. A caution that cost three results

Three times tonight and before, a **structural** fact was mistaken for a **cost** fact:

1. `runs`: "the walk is $0.15n^2$ steps, so $e$ is that" — false, the certificate is a
   single erosion path ($e=\Theta(n)$).
2. `pal`: "the naive machine is $n^{2.73}$" — false, the machine was naive ($g=\Theta(n^2)$).
3. `NBLK`: "band counts do not determine the count" — true, and it bounds **no machine's
   reads**: a $B/3$-witness machine answers covered queries with 0 reads.

The rule that would have caught all three: **a structural statement becomes a cost
statement only after you exhibit the machine and count its cell reads.** Every table in
this note follows that rule.

---

## 7. Correction (02:00): g is a min over machines, so a family bounds g above, not below

$$g_\rho(n)=\min_{\text{machines }M}\ \max_{\text{tapes }\tau}\ \mathrm{reads}(M,\tau).$$

So **any** machine's cost is an **upper** bound on $g$; a hard tape family shows that the
machines we have pay a lot, not that every machine must. The `tri16` family's
$g/n^2=0.943,1.633,3.019,5.580$ ($n=512..4096$, $B=16$, exponent $2.86$; $2.976$ at $B=32$,
exactly halved, consistent with $n^3/B$) is therefore a statement about the frontier machine,
namely $g\le5.58\,n^2$ at $n=4096$ — **not** $g=\omega(n^2)$.

To refute Conjecture 2 one must show that **every** class machine pays $\omega(n^2)$ on some
family. That is the zero-answer lemma, and it is the only route. Until it is settled the
honest position is:

| proposition | status |
|---|---|
| some class machine attains $O(n^2)$ on `NBLK` | **unknown** |
| best known class machine | $5.58\,n^2$ at $n=4096$ ($B=16$), growing like $n^{2.86}$ |
| out-of-class machine ($\Theta(\sqrt n)$ words) | $0.28\,n^2$ (steps $0.063n^2$ x 4.5 reads/test) |

**This correction matters because the same error — reading a machine's cost as a lower bound —
is the one that produced `runs`' "walk length is $e$" and `pal`'s "$n^{2.73}$". It is now the
fourth occurrence.**
