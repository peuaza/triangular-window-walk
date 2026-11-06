# ALLDISTINCT: a rule that forces n^3

*2026-11-06. The construction proposed to fill the k=6 rung of the ladder. The certificate
argument is below; the measurement is delegated.*

---

## 1. The rule

$$b_1(\ell,m)=\big[\text{the }m\text{ values in the window are pairwise distinct}\big],
\qquad b_2=\text{as in \texttt{KCHG}}.$$

Everything else as in the four-valued model: the walk, the cuts, the termination.

## 2. Certificate lower bound

> **Lemma.** A certificate of the answer $b_1=1$ has $\lvert S\rvert\ge m$.

*Proof.* Suppose some cell $i$ in the window is not read. Let the tape's values elsewhere be
as in the certified tape. Change cell $i$ to any value already present in the window at a read
position. The answer becomes 0, and every read position is unchanged — so the two tapes are
consistent with the read set but receive different answers. Hence the read set does not
certify the answer, contradiction. ∎

Note what makes this work and why it fails for most statistics: the negative answer has a
constant-size certificate (two equal cells), so the lemma is expensive only on the *yes*
answer — and $g$ is a maximum over tapes, so a tape that makes the walk ask yes-questions
often is enough.

## 3. Why the state cannot rescue it

A machine in class holds $O(\log n)$ bits plus $w=O(1)$ registers. To answer "are these $m$
values distinct" without rereading, the state must in effect hold the set of values present —
$\Theta(m\log n)$ bits for $m=\Theta(n)$, and $\Theta(m)$ bits even for $m=\sqrt n$ with a
range of $n$. That is out of class.

The residual form of the trace lemma then applies with $x=\ell$ (no maintained aggregate covers
the region) and gives reads $\ge\lvert S\rvert-w-O(1)=\Omega(m)$ per yes-query.

## 4. Consequence

$\Theta(n^2)$ steps, a positive fraction of which are yes-queries (a random window of length
$m$ is all-distinct with probability $\approx e^{-m^2/2n}$, which is $\Theta(1)$ for
$m=O(\sqrt n)$), and $\Omega(m)$ reads each:

$$g_{\text{ALLDISTINCT}}=\Theta(n^2\cdot n)=\Theta(n^3).$$

This is the first rule — as opposed to machine — that attains the general upper bound.

## 5. The ladder, with the gap named

| $c_\rho$ | $g_\rho$ | instances |
|---|---|---|
| $O(1)$ | $\Theta(n^2)$ | 13 rules |
| $\Theta(\sqrt n)$ | $\Theta(n^{2.5})$ | `NBLK>2/>3` |
| $\Theta(n)$ | $\Theta(n^3)$ | `ALLDISTINCT` (this note) |
| $o(1)$ | — | impossible: a certificate has at least one cell |

So the half-integer ladder $n^2,n^{2.5},n^3$ is not an accident of sampling: it is the three
values the certificate can take, and $n^{1.5}$ is excluded by $c\ge1$ unless a rule has
$o(n^2)$ steps, which the census has not found.

## 6. What must be measured

1. That the walk is well defined and terminates on the cube for $n\le10$ (verify against a
   brute-force reference, as with every other rule here).
2. `steps/n^2` for `ALLDISTINCT` — it should be $\Theta(n^2)$ like the others, not larger.
3. `reads/n^2` with the best machine available — the prediction is $\Theta(n)$, i.e. a
   log-log exponent of 3.
4. The yes-rate: the fraction of queries whose answer is 1, and its trend in $n$.

If (3) comes out at 3 and (2) is flat, the ladder is filled and Conjecture 2's replacement
$g=\Theta(n^2c_\rho)$ has an instance at every rung.
