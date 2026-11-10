# SEMI-LINEAR, REALISED

**Verdict.** A legal rule with $e_\rho=\Theta(n)$ and $g_\rho=\Theta(n^2)$ exists, and it is
measured here at four sizes on one machine:

> **semi-linear realised: `e ~ n^1.0` (exponent 0.997), `g ~ n^2.0` (exponent 1.995), `G ~ n^1.0`**
> — over $n=512,1024,2048,4096$, 200 meeting tapes per cell, one machine, $e/n$ flat at
> $6.11$–$6.24$, $g/n^2$ flat at $0.751$–$0.759$, $G/n$ flat at $5.469$–$5.476$.

The rule is `ballast`, `sl2_engine.py:walk2` rule `R_BALLAST`. It is not a natural rule: it is
built for this cell, which is what the cell needs — the natural catalogue has no member and the
construction is the only evidence the class is occupied at all.

Two things make it work, and both are the opposite of the obvious design:

1. **the cheap branch is cheap in *reads* while its output stays content-dependent.** The normal
   branch is a monotone collapse: exactly $n-2$ steps, no cut, one read per step. It is not a
   *fixed* action — the side is the canonical bit of the window, so $\delta$ moves by a third of
   the tape from one tape to the next (sd $0.35\,n$ cells). The answer is not a function of the
   geometry, so no 0-read machine exists (§7).
2. **the defect is a one-shot property of the whole tape, read once at the opening, and the
   branch it selects is then carried by the geometry.** It fires with probability $1/n$ and the
   branch it opens costs $0.75\,n^2$ reads. Mean: $5.44\,n + n^{-1}\cdot0.75\,n^2 = 6.2\,n$.

---

## 1. What had to be shown, and what the model allows

The model (bthm_jalc §2; FORMAL_MODEL.md §2–3) is a memoryless strategy

$$f\bigl(\varphi(\tau[\ell..r]),\; n,\; \ell,\; r,\; L,\; R\bigr)\longrightarrow
\{\ell{+}1,\ell{-}1,r{+}1,r{-}1,\text{Lcut},\text{Rcut},\text{DCut},\text{halt}\},$$

with $\varphi$ a fixed function of the window with a finite value set, the geometry part of the
input, reads charged one per cell sensed (counted with multiplicity; the state's contents free,
`space` charged in bits), and $O(\log n)$ bits of state for the class.

| what the class needs | what that means operationally |
|---|---|
| $e_\rho=O(n)$ | a machine, exact on every meeting tape, whose **mean** reads is linear |
| $g_\rho=\omega(n)$ | that same machine's **max** over the meeting tapes (eq. `three`, bthm_jalc:477-490) is superlinear |
| legal | a memoryless $f$ of the shape above; actions from the eight; no move outside the region |
| real | both measured, on one machine, at four sizes |

The three definitions, quoted for the record: $e_\rho=\min_{M\in\mathcal M_\rho(n)}\mathbb E_\tau C_M$,
$g_\rho=\max_\tau C_{M^\star}$ with $M^\star$ the $e$-attaining machine,
$G_\rho=\exp(\mathbb E_\tau\ln C_{M^\star})$. Everything is over the **meeting** tapes of $T_n$;
a tape that halts has no $\delta$ to reproduce and is not in the domain.

## 2. The rule

$\varphi(w)=(c(w),\,z(w))$, two bits.

* $c(w)$ — the canonical bit of the window, verbatim from `pcr_lib.b1_of`'s `R_CANON` branch
  (the normalised lag-one correlation test with its $m=1,2$ cases and its zero on a constant
  half). Content-dependent, $O(1)$ to evaluate from three window aggregates.
* $z(w)=1$ iff $\;h(w)=\sum_i a_i 3^i \bmod (2^{31}-1)\;$ is divisible by $n$, where
  $a_i=\sigma_i=b^{(1)}_i+b^{(2)}_i\in\{0,1,2\}$.

$f$, read off the geometry — the rule never cuts in the normal branch, so $L=0$ *is* the
normal branch and $L\ge1$ *is* the defect branch:

| condition | action | why |
|---|---|---|
| $L=0,\ \ell=0,\ r=n-1,\ z=1$ | **left cut** | the one-shot defect: taken at step 1, $L$ becomes 1 and stays $\ge1$ |
| $L=0$ otherwise | $\ell{+}1$ if $c=1$, else $r{-}1$ | the monotone collapse: $n-2$ steps, no cut |
| $L\ge1,\ W=R-L+1\le4$ | $\ell{+}1$ if $c=1$, else $r{-}1$ | the closing collapse — this is where the run meets |
| $L\ge1,\ 2m>W$ | $r{-}1$ | shrink the window to half the region |
| $L\ge1,\ 2m\le W$ | left cut if $c=0$, right cut if $c=1$ | the epoch's cut, on the side the content names |

**The defect branch costs $0.75\,n^2$ reads** because an epoch of region width $W$ shrinks the
window by $W/2$ ($W/2$ reads) and then cuts — and a cut is a jump, so the $O(\log n)$-bit state
re-acquires the window aggregates at a cost of $W/2$ reads. The region loses one cell per epoch,
so the branch totals $\sum_{W=n}^{4}\Theta(W)=\Theta(n^2)$. **The normal branch costs
$\Theta(n)$**: no cut ever, hence no re-span, hence one read per step.

Both branches meet, always (measured: 200/200 and 4000/4000 at every size, §5). The normal branch
meets when $\ell=r$ after $n-2$ inward moves; the defect branch shrinks to $W=4$ and collapses
from there, and it can never empty the region because it only cuts while $W>4$.

## 3. Legality, clause by clause

1. **Shape.** $f$ takes $(\varphi(w),n,\ell,r,L,R)$ and returns one of the eight actions. Nothing
   else is consulted. $\varphi$'s range is $\{0,1\}^2$, finite. ✔
2. **No memory.** The one bit the defect needs is *not* a register: the branch is identified by
   $L$, which is part of the state the model already gives $f$. Normal $\Rightarrow L=0$ for the
   whole run (the branch has no cut); defect $\Rightarrow L=1$ from step 1 and $L\ge1$ forever
   (the branch only left-cuts). So the rule is a memoryless $f$, not a machine with a hidden bit.
   (`walk2` keeps the opening bit in a local, but it is consulted only at the step that produced
   it; the run's branch is decided by $L$ from step 2 on.)
3. **The defect is a function of the window.** At step 1 the window is $[0,n-1]$, i.e. the whole
   tape, so a full-tape hash *is* a $\varphi$ value. It is evaluated once. $\varphi$ is never
   re-evaluated for $z$ afterwards, and no later window can satisfy the first row ($\ell=0$ and
   $r=n-1$ never recur, since the normal headers only move inward and the defect branch leaves
   $L\ge1$).
4. **Reads for the defect are charged.** Computing $z$ is one pass over the opening window,
   $n$ reads, and the machine pays exactly that pass anyway ($n+2(n-1)+2+\lfloor n/2\rfloor$,
   cr_lib's window build). No free lunch is taken: the charge is in `inc` and is what the gate
   checks.
5. **Every action is legal.** The normal branch moves $\ell{+}1$ or $r{-}1$ with
   $L\le\ell<r\le R$: inside the region, and when the two headers meet the run outputs. The
   defect branch moves $r{-}1$ with $\ell=L<r$, and cuts only while the region keeps $\ge4$
   cells — so no cut empties the region and no move leaves it.
6. **The cut convention is not in play.** This rule never attempts a move that would leave the
   region, so neither the `a` nor the `c1a` convention (RULE_CATALOGUE §0.2) is exercised. The one
   convention it does use, "a cut removes one cell from the side it names, headers restart",
   is the model's (bthm_jalc:185-190), and it is what `pcr_lib.walk` implements.

## 4. The gate

`sl2_engine.walk2` is a re-statement of `pcr_lib.walk` with the same action table, the same cut
semantics, the same INC accounting, the same repeat stamp — and the canonical rule in both. The
gate demands exact agreement, tape by tape, on **six** fields: output, steps, m-sum, INC reads,
cut count, ending; both cut conventions.

```
n=128   40 tapes x 2 conv   mismatch so far 0/80
n=256   40 tapes x 2 conv   mismatch so far 0/160
n=512   40 tapes x 2 conv   mismatch so far 0/240
n=1024  40 tapes x 2 conv   mismatch so far 0/320
GATE PASS  (320 runs checked, 0 field mismatches)
```

(`sl2_gate.log`; 1920 field comparisons.) A hand-rolled walk that does not reproduce `pcr_lib` is
void; this one does.

## 5. The two traps, both measured on this construction

**Trap 1 — an insensitive $\delta$ has $e=O(1)$.** Control `geo`: the same monotone collapse with
the content bit removed. Its $\delta$ takes **1 distinct value** at $n=512,1024,2048$ (sd 0.0
cells) — a machine that outputs it reads nothing, so $e=0$ and the rule is trivially linear, not
semi-linear. `ballast`'s $\delta$ takes 147–193 distinct values on 200 tapes and has sd
$0.35\,n$ cells. No 0-read machine exists for it. Details in §7.

**Trap 2 — a per-query defect fires on every tape.** Control `perquery`: the same epoch branch,
reached whenever $h(w)\bmod m=0$ for the *current* window ($\Pr=1/m$ per query). Hazard
$\sum_m 1/m\approx\ln n$ per walk, so

| $n$ | tapes where the per-query defect fires | first fire at window |
|---|---|---|
| 512 | 199/200 (99.5 %) | $m=282\approx0.55\,n$ |
| 1024 | 200/200 (100 %) | $m=569\approx0.56\,n$ |
| 2048 | 200/200 (100 %) | $m=1068\approx0.52\,n$ |

The rare set is gone: the defect fraction is a constant, every tape pays the $\Theta(n^2)$
branch, and the rule is simply superlinear. This is why the defect must be one-shot and
tape-level, and why the branch is carried by the geometry afterwards.

## 6. The table

Machine: the `walk2` INC machine, one machine for all three costs. Tapes: `pcr_lib.make_tapes`,
$T_n$ (two independent uniform $h$-subsets, $h=n/2$), seed $2026+n$, 200 per cell, **all** of
which meet. $q$ = $\Pr[z=1]$, measured on 200 000 hash-only draws per cell.
$e=(1-q)C_{norm}+qC_{def}$, with $C_{def}$ measured on three defect witnesses per cell — a 200-tape
uniform sample cannot see a $1/n$ event at $n=4096$ (it contains 0.05 of one), so the tail is
measured where it lives and the mixture is reported. The direct 4000-tape means are in the
second table for comparison.

| $n$ | tapes | meet | $e/n$ | $g/n^2$ | $G/n$ | $C_{norm}/n$ | $C_{def}/n^2$ | $q$ | $1/n$ | $q\cdot n$ |
|---|---|---|---|---|---|---|---|---|---|---|
| 512 | 200 | 200/200 | 6.159 | 0.7593 | 5.471 | 5.438 | 0.7569 | 1.88e-03 | 1.95e-03 | 0.965 |
| 1024 | 200 | 200/200 | 6.202 | 0.7551 | 5.469 | 5.454 | 0.7539 | 9.75e-04 | 9.77e-04 | 0.998 |
| 2048 | 200 | 200/200 | 6.239 | 0.7525 | 5.476 | 5.472 | 0.7517 | 5.00e-04 | 4.88e-04 | 1.024 |
| 4096 | 200 | 200/200 | 6.106 | 0.7513 | 5.473 | 5.477 | 0.7510 | 2.05e-04 | 2.44e-04 | 0.840 |

| exponent (log–log OLS, 4 points) | value |
|---|---|
| `e` | **0.997** |
| `g` | **1.995** |
| `G` | **1.000** |
| $C_{norm}$ | 1.004 |
| $C_{def}$ | 1.996 |
| $q$ | −1.057 |
| `e`, 200-tape raw (no tail term) | 1.006 |

The raw 200-tape cell, before the tail is added:

| $n$ | $e_{raw}/n$ | $g_{raw}/n$ | $G_{raw}/n$ | steps$/n$ | $\delta$ distinct | sd($\delta$) in cells |
|---|---|---|---|---|---|---|
| 512 | 5.439 | 5.986 | 5.427 | 0.9980 | 147 | 182.8 |
| 1024 | 5.425 | 5.995 | 5.414 | 0.9990 | 171 | 353.5 |
| 2048 | 5.433 | 5.997 | 5.422 | 0.9995 | 184 | 720.9 |
| 4096 | 5.508 | 5.998 | 5.496 | 0.9998 | 193 | 1451.9 |

and the independent 4000-tape samples (which do contain the defect tapes — 10, 4, 3, 1 of them):

| $n$ | $e_{4000}/n$ | max$/n^2$ | defects in the sample |
|---|---|---|---|
| 512 | 6.393 | 0.7593 | 10 |
| 1024 | 6.220 | 0.7542 | 4 |
| 2048 | 6.623 | 0.7525 | 3 |
| 4096 | 6.244 | 0.7507 | 1 |

The two estimators agree within the sampling error of the tail ($1/\sqrt{k}$ with $k=1..10$ is
30–100 %); the mixture is the one to quote, because its tail term is measured on 3 witnesses
instead of on one lucky draw.

Defect cost, per witness, INC reads in units of $n^2$: 0.7574/0.7579/0.7554 at 512,
0.7538/0.7528/0.7551 at 1024, 0.7520/0.7521/0.7510 at 2048, 0.7509/0.7513/0.7508 at 4096 — flat,
so the defect branch is $\Theta(n^2)$ with a constant that does not drift. Their $\delta$ is
0.458/0.517/0.192 and 0.425/0.840/0.291: **the defect tapes do not all output the same number**,
so $g$ is the max of a hard output and not of a constant.

## 7. Is $\delta$ a function of the tape, and of how much of it?

Four probes, `sl2_content.py`, all measured.

**(A) Spread.** $\delta$'s sd is $0.35\,n$ cells at every size (table above), quartiles
0.096/0.402/0.809 at $n=512$ and 0.139/0.541/0.864 at $n=4096$. A geometry-determined $\delta$ has
sd 0 (control `geo`, §5).

**(B) One flip pair** — `pb_content.py`'s perturbation, one 2-cell $\to1$ and one 0-cell $\to1$,
which keeps $2\#2+\#1=n$ and so keeps the tape in $T_n$:

| $n$ | probes | unchanged | median move (when it moves) | p90 | max |
|---|---|---|---|---|---|
| 512 | 240 | 92.1 % | 19 cells | 505 | 505 |
| 2048 | 120 | 92.5 % | 661 cells (0.32 n) | 1612 (0.79 n) | 1770 (0.86 n) |

So the trap's strong form — "a flip moves $\delta$ by a large fraction of $n$ on at least half the
meeting tapes" — is **not** met, and the reason is worth stating rather than hiding: $c(w)$ is the
sign of a $\Theta(m^2)$ statistic, so a $\pm1$ perturbation of it does not change the sign at most
windows, and the normal branch's $\delta$ is stable under such a perturbation. What the strong
form is a proxy for is "$e$ is not $O(1)$", i.e. "no cheap machine knows $\delta$", and that is
measured directly in (C) and (D).

**(C) The defect flip.** On a defect tape, one flip pair changes $h$, hence $z$, hence the branch,
and $\delta$ jumps by half the tape. Every 2-cell of the tape was probed with its own partner (the pair flip certifies the 0-cell too,
so half the cells are covered this way; the 1-cells are not, and (D) covers all of them):

| $n$ | 2-cells probed | moved $>n/4$ | median move |
|---|---|---|---|
| 512 | 117 | 116 (99.1 %) | 299 cells (0.58 n) |
| 2048 | 520 | 520 (100 %) | 1417 cells (0.69 n) |

Every cell of the tape is in the certificate set of some tape: drop cell $i$ and there is a tape
whose $\delta$ differs by $\Theta(n)$ from its partner differing only at $i$. That is the
exchange argument of FORMAL_MODEL.md:124-130, and it is what makes $e=\Omega(n)$ — the argument
cannot give more than a linear bound, and a linear bound is exactly the claim.

**(D) Erasure** — the direct form of the same question. Erase a fraction $f$ of the cells (impute
1, the modal value) and re-run; "reproduced" means the same $\delta$ as the full tape:

| $f$ erased | $n=512$ reproduced | $n=2048$ reproduced |
|---|---|---|
| 0.001 (1–2 cells) | 38/40 | 28/30 |
| 0.01 | 35/40 | 24/30 |
| 0.05 | 23/40 | 15/30 |
| 0.25 | 5/40 | 3/30 |
| 0.50 | 2/40 | 0/30 |

A machine that reads $0.95\,n$ cells reproduces $\delta$ on 57 % ($n=512$) and 50 % ($n=2048$) of
the tapes; one that reads $0.75\,n$ on 12 % and 10 %. Exactness on **every** tape therefore needs
essentially all $n$ cells: $e\ge n\,(1-o(1))$. Combined with the machine of §6 ($6.1\,n$) this
gives $e=\Theta(n)$ with both bounds measured.

## 8. Accounting, and what survives it

$e$ is machine-relative in this model, so the accounting is named with the number. The INC
accounting above is `pcr_lib`'s (the $O(\log n)$-bit window machine: 1 read per cell taken in or
given up, one re-span per cut). Under the other accounting `pcr_lib` reports, SCAN ($m$ reads per
query, "a pass over a window of length $m$ costs $m$ reads"), the same rule is

| accounting | $e$ | $g$ | $g/e$ |
|---|---|---|---|
| INC (this report) | $6.1\,n^{1.00}$ | $0.75\,n^{2.00}$ | $\Theta(n)$ |
| SCAN | $0.62\,n^{2}$ | $0.125\,n^{3}$ | $\Theta(n)$ |

because the normal branch pays $0.5\,n^2$ per tape (a pass per step over a shrinking window) and
the defect branch $0.125\,n^3$. The **separation** $e\ll g$ is not an artefact of the accounting;
the *exponent* of $e$ is, and the paper's fully-linear cell has the same property (under SCAN,
`pattern` and `runs` are $\Theta(n^2)$ too). The class table is INC-relative and this row is
consistent with it.

## 9. What this shows, and what it does not

**Shown.**

* A rule of the model's own shape (`memoryless f(φ,n,ℓ,r,L,R)`, legal actions, $O(\log n)$ bits)
  with $e/n$ flat at 6.1–6.2 and $g/n^2$ flat at 0.751–0.759 over a factor of 8 in $n$, measured
  on one machine, 200 meeting tapes per cell, gated against `pcr_lib` at 320 runs and 1920 field
  comparisons.
* The mechanism is exactly the one the cell needs and is not a loophole: a $1/n$-probability
  tape-level defect read once at the opening (measured $q\cdot n=0.965,0.998,1.024,0.840$; the
  decay exponent $-1.057$), a $\Theta(n^2)$ branch behind it ($0.75\,n^2$, exponent 1.996), a
  $\Theta(n)$ branch in front of it ($5.44\,n$, exponent 1.004).
* Both traps fire if the design is wrong and do not fire here: a geometry-only $\delta$ has
  $e=0$ (control `geo`), and a per-query defect fires on 199/200–200/200 tapes (control
  `perquery`).
* $e=\Omega(n)$ by measurement as well as by argument: erasing 5 % of the cells breaks $\delta$
  on 43–50 % of the tapes, and every 2-cell and every 0-cell (through its pair partner) appears
  in the certificate set of some tape through the defect bit; the 1-cells are covered by the
  erasure test only.

**Not shown, and not claimable.**

* $6.1\,n$ is an **upper** bound from this machine (3.5 n of window build + 1 n of moves + 0.94 n
  of head drift). A tighter machine would build all three aggregates in one pass and land near
  $2.5\,n$ — still $\Theta(n)$; no exponent claim depends on the constant.
* $g$ is $\max_\tau C_{M^\star}$ on the $e$-best machine *we can build*, which is how eq. (three)
  defines it. The model's only lower-bound tool is indistinguishability, which
  (FORMAL_MODEL.md:128) "can never prove more than a linear bound" — so no $g=\omega(n)$ claim in
  this model can be a proof, for this rule or for any other. What is measured: any machine that
  reproduces the defect branch pays $\Theta(n^2)$, because the branch's $\delta$ is the fixed
  point of a content-driven epoch walk with no closed form, and the branch cannot be identified
  without a full-tape read.
* The rule is **constructed**, not natural. It occupies the semi-linear cell; it does not make
  the cell natural-looking, and it does not touch the paper's claim that no natural rule in the
  catalogue is in it.
* $n\le4096$. The $1/n$ rate is measured over a factor of 8; the $\Theta(n^2)$ tail is measured on
  3 witnesses per cell, and their spread is under 1 % ($0.7508$–$0.7579\,n^2$ across all sizes).

## 10. Files

| file | what |
|---|---|
| `sl2_engine.py` | `walk2`: `pcr_lib.walk` re-stated for two rules, `ballast` and the gate's `canonical`; `tape_hash`, `Runner2` |
| `sl2_gate.py`, `sl2_gate.log` | the gate: 320 runs, six fields, both cut conventions |
| `sl2_measure.py`, `sl2_measure.json`, `sl2_measure.log` | the cells: 200 + 4000 tapes, 3 defect witnesses, $q$ on 200 000 draws, the $e/g/G$ estimators |
| `sl2_table.py`, `sl2_table.log` | renders §6 from the JSON — no number in the table is typed by hand |
| `sl2_content.py`, `sl2_content.log`, `sl2_content_2048.log` | §7: spread, pair flip, defect flip, erasure |
| `sl2_controls.py` | §5: `geo` (trap 1) and `perquery` (trap 2) |
| `sl2_calib.py`, `sl2_smoke.py`, `sl2_smoke2.py`, `sl2_find.py` | calibration against `pcr_lib` and the first defect-tape search |
| `run_sl2.sh` | reproduces everything in order |

No `.tex` was edited; nothing was committed.

The `sl_*.py` files in this directory predate this report (the 2026-11-05 search); nothing here
reads or writes them, and none of their numbers is quoted above. This report's files are the
`sl2_*` set.
