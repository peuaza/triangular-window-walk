# Filling the top rung: can a window function force n^3?

*2026-11-06. The construction proposed for the k=6 rung, the correction that killed its first
form, and what a valid construction has to satisfy.*

---

## 0. Correction before anything else

The rule "$b_1=1$ iff the window's $m$ values are pairwise distinct" is **empty in the
four-valued model**: there the cell values are $\sigma_i\in\{0,1,2\}$ — three values — so the
answer is 0 for every $m\ge4$, and the yes-answer never occurs on an interesting window.

The construction has to be in the **binary** model, and even there it needs a property that a
counter cannot maintain. What survives is the certificate argument, which is correct as far as
it goes and is the template for the replacement.

## 1. The certificate template (correct)

> **Lemma.** If every unread cell of a window can be changed so as to flip the answer without
> touching any read position, then every certificate of that answer reads all $m$ cells.

*Proof.* Let $i$ be unread and change it. Every read position keeps its value, so by the trace
lemma — two tapes agreeing on the read set produce the same transcript — the machine gives the
same answer on both tapes, yet the answers differ. So the read set does not certify. ∎

## 2. Why that is not enough on its own

"The window has at most $m/2$ ones" satisfies the template: an unread cell can be set to 1, and
its certificate is $\Theta(m)$. Yet a machine keeps a counter and pays $O(1)$ per step. The
certificate bound is about **state**, not about the query — in the residual form
$\text{reads}\ge(r+1-x)/(2R+1)-w$ the maintained count gives $x=r+1$ and the bound is vacuous.

**So a valid construction needs a property whose statistic does not fit in $O(\log n)$ bits.**
That requirement is much stronger than "the certificate is large", and it is the one the first
version of this note missed.

## 3. The candidate

In the binary model, read the window as an $m$-bit integer and take

$$b_1(\ell,m)=\big[\text{that integer is prime}\big].$$

* The integer is $m$ bits. No $O(\log n)$-bit machine can carry it forward: each step replaces
  one bit at one end, and the resulting number is not recoverable from any $O(1)$-word summary
  the machine could keep.
* The answer depends on every bit, so §1 applies.
* No counter-style shortcut exists for primality, unlike for a sum or a count.

**This is the candidate to test**, not the distinctness rule.

## 4. What a measurement must show

1. The rule is well defined and the walk terminates — verified against brute force on the whole
   cube $n=3\ldots10$.
2. `steps/n^2` is $\Theta(n^2)$: the census gives $0.004$–$0.06$ for seven rules, so a wildly
   larger value is itself a finding.
3. `reads/n^2` grows with exponent $\approx3$.
4. The yes-rate and its trend.

If instead some $O(1)$-word machine answers the query in $o(m)$ reads on the tapes where a scan
pays $m$, the construction fails — and that outcome should be reported plainly, because it
would say the top rung is not reachable by a window function at all.

## 5. The ladder, with the gap named

| $c_\rho$ | $g_\rho$ | instances |
|---|---|---|
| $O(1)$ | $\Theta(n^2)$ | 13 rules |
| $\Theta(\sqrt n)$ | $\Theta(n^{2.5})$ | `NBLK>2/>3` |
| $\Theta(n)$ | $\Theta(n^3)$ | **not yet constructed as a rule**; machines attain it ($0.0119n^3$) |
| $o(1)$ | — | impossible: a certificate has at least one cell |

So $n^{1.5}$ is excluded by $c\ge1$ unless a rule has $o(n^2)$ steps, which the census has not
found; and the top rung is the one still to be earned.
