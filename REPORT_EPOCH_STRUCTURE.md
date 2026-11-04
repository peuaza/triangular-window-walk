# Structure of the epoch sequence of the b1≡0 walk — is there a one-pass algorithm?

Question: is the sequence of cut positions m*(L) (and the epoch lengths cost_L = M_L − m*_L + 1)
monotone, locally computable, or otherwise regular enough that the whole run can be produced in
o(n²) instead of by the M, M−1, … scan?

**Verdict: no exploitable structure was found, and the negative is measured, not asserted.**
Every route that would give a one-pass or amortised algorithm fails a specific measurement:

| route to o(n²) | what it needs | measurement that kills it |
|---|---|---|
| two-pointer sweep on m* | m* non-increasing | m* **rises** on 7.9–9.8 % of steps overall, and on **18–22 %** of the steps that are not forced down by a double cut (§3) |
| sweep with warm start | verdicts inherited: (L,m) fails ⇒ (L+1,m) fails | all four cells of the 2×2 table are populated; 4.2–10.7 % of verdicts flip per step and Θ(M) flip per epoch (§3.2) |
| exact O(1) recurrence for m* | cost_L = f(M_L, m*_{L−1}, a few cells) | explicit counterexamples: identical (M, parity of M, previous cost bucket, 4 prefix cells, 4 suffix cells) with cost **1 vs 93** out of M = 98 (§4.3) |
| O(1) certificate that an epoch is long | a few cells imply cost ≥ M/4 | R² of the best linear predictor of cost/M grows **linearly in w with no saturation**: 0.0055 (w=1) → 0.24 (w=64) at n = 512; best single feature lifts P(long) from 0.266 only to 0.38–0.46 (§5.4) |
| warm-start at m*_{L−1}+1 (best heuristic) | the answer moves by O(1) | right on 92.9 % of epochs and 80 % cheaper, but wrong on 7 % with no O(1) way to detect it (§4.4) |
| maintain the concavity set S_L incrementally | \|S_L Δ S_{L+1}\| = O(1) | \|S_L Δ S_{L+1}\|/M = 0.107 at n = 512 (8.9 % lost, 9.1 % gained); S_L is nested in S_{L+1} on 2.1 % of steps (§6.2) |
| skip a run of equal verdicts | the bit field piecewise constant at scale ≈ cost | adjacent-m flip rate 4.3–12.9 %; runs ≈ M^{0.54} (§5.5), i.e. mean run 16 cells at n = 512 against a mean scan of 48 |

What *is* strongly structured is the **cost process in bulk**, not the individual answer:
β_L := 1{cost_L = 1} forms heavy-tailed sojourns (mean 27–32 epochs at n = 512, median 4–7, max 462)
and H(β_L | β_{L−1}) = 0.197 bits. So the walk's *output* is cheap to describe (Θ(n) bits) while its
*verification* is what costs 0.059 n² tests. That gap is the honest statement of the negative:
the trajectory is predictable, the answers are not checkable without being read.

All numbers below are from the scripts in this directory; command and seed with every table.
No .tex file touched.

---

## 1. Setup and validation

Model exactly as `q0_walk.py`: epoch at [L,R], M = R−L+1, m* = largest m ∈ [2,M] with
A_{⌊m/2⌋} ≥ A_m (A_j = mean of a[L..L+j−1]), cost = M−m*+1, then L += 1 and R −= 1 iff m* = M.
Ensemble: two binary columns with ⌈n/2⌉ ones each (`sample_tapes`).

`q21_epoch_structure.py` builds one record per epoch (L, R, M, m*, cost, top-bit, concavity mask).
Two independent checks:

* `assert sum(cost) == q0_walk.walk_b0_fast(a)[1] == scan_naive(a)` — the last being a literal
  `for m in range(M,1,-1)` scan — on 4 hand tapes and 300 random tapes at n = 40: **0 mismatches**
  (printed as the first line of every `q21_*.log`).
* the pipeline reproduces the published cost (`REPORT_LOWER.md` §7 gives 0.0540 ± 0.0014 at
  n = 1024, 0.0556 at 256/512):

```
$P q21_epoch_structure.py mono 64,128,256,512 400 2026      # and q25_nesting.py C
   n    tests/epoch  epochs/n   tests/n^2
 128        12.70      0.612      0.0607
 256        22.78      0.622      0.0553
 512        48.39      0.640      0.0604
```
epochs/n = 0.61–0.64 and cost/M ≈ 0.156–0.173 agree with `REPORT_LOWER.md` §7 (0.623–0.642,
0.155–0.185). So the measurements below are of the same object as the paper's.

**The two exact reformulations used throughout** (both verified, 0 mismatches):

* m = 2k even: concavity ⟺ `P[L] + P[L+2k] ≥ 2 P[L+k]` — the prefix-sum point (L+k, P[L+k]) lies
  on or **below the chord** from (L, P[L]) to (L+2k, P[L+2k]); 0/16700 mismatches (`q27_chord.log`).
* m = 2k+1 odd: concavity ⟺ mean(a[L..L+k−1]) ≥ mean(a[L+k..L+2k]) — first k cells vs last k+1
  cells; 0/16340 mismatches.
* identity checks inside the walk (`q24_hazard.py D`, 112 986 epochs): cost=1 ⟺ concavity at m=M
  (0 mismatches); even-M top test ⟺ mean(second half) ≤ mean(first half) (0/81563); odd-M top test
  ⟺ mean(first k) ≥ mean(last k+1) (0/31423); M_{L+1} = M_L − 1 − 1{cost_L=1} (0/112477).

The walk is therefore a scan over **chord tests with stride k = ⌊m/2⌋ at a fixed left endpoint**,
and the trajectory is driven by one bit per epoch: β_L = 1{cost_L = 1} (double cut).

## 2. What "the cut position" is

Two sequences matter and they behave differently:

* **m\*_L** — the scan's stopping index; determines cost. Range 2..M_L.
* **R\*_L = L + m\*_L − 1** — the right end of the stopping window. Note ΔR\* = Δm\* + 1, so
  "R\* non-increasing" and "m\* strictly decreasing" are the *same* event; only one of them is
  reported below (columns `P(dR*<=0)` and `P(dm<0)` are identical by algebra, as they must be).

## 3. Q1 — is m*(L) monotone?  No.

`$P q21_epoch_structure.py mono 64,128,256,512 400 2026` (`q21_mono.log`), replicated with seed
424242 (`q21_mono.log`, second run, same conclusion):

```
   n   pairs   P(dm<0)  P(dm=0)  P(dm>0)  P(|dm|<=1)  P(dm<0|lc)  P(dm>0|lc)
  64   14741    0.875    0.033    0.092     0.139       0.695       0.224
 128   30853    0.885    0.030    0.085     0.131       0.716       0.209
 256   63364    0.889    0.029    0.081     0.131       0.733       0.196
 512  129318    0.891    0.030    0.079     0.134       0.747       0.182
seed 424242:  0.871/0.879/0.901/0.893 falling, 0.098/0.089/0.072/0.078 rising
dm* quantiles 1/5/50/95/99 % (n=512): -71 -6 -2 +5 +83 ; mean -1.36, sd 24.1
```

Read this carefully, because the naive reading is wrong:

* The 89 % "decrease" is mostly **forced**. A double cut sets m*_L = M_L and M_{L+1} = M_L − 2, so
  m*_{L+1} ≤ m*_L − 2 automatically. On the 62 % of steps that are double cuts, monotonicity is a
  tautology, not evidence.
* On the steps where m* could rise — the **left cuts** (cost > 1), which are exactly the expensive
  epochs — **18–22 % do rise**, and by a lot: the 99th percentile of Δm* is +83 at n = 512 (0.16 n).
* P(|Δm*| ≤ 1) = 0.13. There is no local "the cut moves by one" relation to exploit.

So the fraction of *consecutive pairs* on which m* rises is 7.9–9.8 %, and 18–22 % conditional on the
only case where it is not forced. A two-pointer sweep that assumed monotonicity would be wrong on
about one in five left cuts — and a left cut is where all the work is (§5.1).

### 3.1 Total movement is cheap — the pointer is not the problem

`$P q22_local_structure.py C 64,128,256,512 300 2026` (`q22_C.log`):

```
   n    TV/n   TV+/n   max|dm*|/n     sym/M   P(up)
  64   1.938   1.509      0.276       0.232   0.096
 128   2.303   1.838      0.257       0.188   0.083
 256   3.071   2.589      0.251       0.143   0.083
 512   4.357   3.863      0.272       0.107   0.081
```

TV = Σ|Δm*| per tape divided by n. TV/n itself grows with n (1.94 → 4.36 over an 8× range, ≈ n^{0.39}),
so a bidirectional pointer would cost ≈ n^{1.4}: still sub-quadratic, but the *prize* is finite —
the ratio cost/pointer is 3.4× at n = 128, 4.6× at 256, 7.0× at 512, i.e. growing like n^{0.65}. A
sweep that worked would be worth a factor of this order, which is why the route deserves the
measurements below rather than a single number.

### 3.2 And monotonicity would not be enough: verdicts are not inherited

A sweep must also know that a verdict it does not re-test is unchanged. That fails at the level of
single cells. `$P q22_local_structure.py E 64,128,256,512 300 2026` (`q22_E.log`), 2×2 tables of
(c(L,m), c(L+1,m)) — c = 1 means concavity (b2 = 0, the scan stops):

```
n=512  same-m  [[13351852, 1112272],        same-R  [[13665587, 817084],
                [ 1087331, 10874566]]                [  773596, 11169754]]
```
Both off-diagonal cells hold ~10⁶ events. **8.3 % of the (L,m) verdicts change when L advances**
(4.2 % 0→1 plus 4.1 % 1→0 unconditionally; 7.7 % conditional on c=0, 9.1 % conditional on c=1).
The same-right-end version (diagonal: R fixed, L and m both shifted) is slightly better — 3.0 % /
3.1 % — and still Θ(M) flips per step. Neither direction has an empty cell, so **no implication
holds in either direction**; there is no cell-level inheritance to build a sweep on.

Per-cell agreement improves with n (0.79 → 0.92 between n = 64 and 512) but so does the scan
length; the ratio is in §5.5.

## 4. Q2 — is the epoch length predictable from local data?

### 4.1 Exact relations that do hold (all checked, 0 mismatches)

1. cost_L = 1 ⟺ concavity at m = M_L ⟺ the top test passes.
2. M_L even: cost_L = 1 ⟺ mean(a[right half]) ≤ mean(a[left half]) — the drift sign.
3. M_L odd: cost_L = 1 ⟺ mean(first ⌊M/2⌋) ≥ mean(last ⌈M/2⌉).
4. M_{L+1} = M_L − 1 − 1{cost_L = 1}; and m*_L = M_L whenever cost_L = 1.
5. Consecutive double cuts give Δm* = −2 exactly — which is why P(Δm* = −2) = 0.68
   (`q23_certificate.py C`). This is a tautology given two double cuts, not a predictive relation.

Beyond these, nothing exact exists. In particular **cost_L is not a function of
(M_L, m*_{L−1}, local cells)** — see the counterexamples in §4.3.

### 4.2 Statistical predictability is strong but not usable

Relative cost classes 0 = {c=1}, 1 = (c/M > 1/2), 2 = (1/4,1/2], 3 = (1/8,1/4], 4 = [0,1/8];
transition matrix (rows = current, `$P q21_epoch_structure.py cost 64,128,256,512 400 2026`):

```
n=512       0      1      2      3      4      n
 0       0.968  0.009  0.006  0.004  0.013  37051        H(b)=1.80
 1       0.045  0.844  0.055  0.020  0.035   7242   H(b|b_prev)=0.59
 2       0.022  0.073  0.836  0.040  0.029   8235        MI=1.21 bits
 3       0.032  0.019  0.107  0.765  0.078   4777
 4       0.110  0.020  0.028  0.075  0.767   6059
```

MI(b_t ; b_{t−1}) = 0.93, 1.07, 1.21 bits at n = 64, 128, 512 — the previous epoch's bucket carries
half to two thirds of the entropy; the "predict the same bucket again" rule is right 0.84 / 0.90 of
the time at n = 64 / 512, against a 0.59 best-marginal base rate. But that predictability is a
*description*, not a rule: from the state alone you cannot say what the next answer is, only what it
probably is, and the algorithm must be exact.

Adding cells to the state confirms this (`$P q24_hazard.py E`, `q24_E.log`):

```
n=512   H(b)=1.875  H(b|prev)=0.488   H(b|prev,4 suffix cells)=0.451   H(b|4 suffix)=1.718
                    H(b|4 prefix cells)=1.856   H(b|prev,4 prefix,4 suffix)=0.201 (17230 levels, sparse)
```
The 4 suffix cells carry 0.16 bits, the 4 prefix cells 0.02 bits. The state "previous bucket +
8 cells" has 17 230 levels for 65 000 epochs, so its 0.20 bits is largely overfitting; the honest
figure is H(bucket | previous bucket) ≈ 0.49–0.72 bits.

The best single "local predictor" of the cost is d := M_L − m*_{L−1}, the gap opened since the last
left cut (`$P q24_hazard.py B`, 200 tapes, seed 2026):

```
   n    corr(c,d)  E|c-d|/M  P(|c-d|<=2)  P(c-d>=3)  sd((c-d)/M)  H(b|d)
 128     0.910      0.105      0.148       0.794       0.172      0.937
 512     0.955      0.045      0.150       0.792       0.114      0.817
```
The correlation is real but it *is* an identity: cost_L = d + 1 − Δm*_L, and Δm* is small and
usually negative (median −2), so cost ≈ d + 3. d is a **lower bound with typical slack 3** and sd of
the slack 0.114 M (P(c−d ≥ 3) = 0.79, tail to 0.16 n). It predicts the *scale* of the cost, not its
value, and saves no test: trusting "cost = d" would be wrong by ≥ 3 cells in 79 % of epochs and the
wrongness is exactly the unknown part of m*.

### 4.3 Explicit counterexamples to any O(1)-state rule

`$P q25_nesting.py B 128,256,512 200 2026` (`q25_B.log`) — search for two epochs with identical
state (M, parity of M, previous cost bucket, a[L..L+3], a[R−3..R]) and different cost:

```
n=128  state=(M=98, even, prev bucket 1, prefix (1,0,0,1), suffix (1,1,1,0))
       cost 1 .. 93   -> cost/M 0.01 .. 0.95
n=256  state=(M=138, even, prev bucket 0, prefix (1,1,0,1), suffix (2,1,1,1))
       cost 1 .. 118  -> 0.01 .. 0.86
n=512  state=(M=452, even, prev bucket 4, prefix (1,1,1,2), suffix (2,1,1,1))
       cost 8 .. 282  -> 0.02 .. 0.62
```

The two n=128 windows (cost 1 and cost 93) share their length, their parity, their previous cost
class and both 4-cell ends, and differ only in the interior. No function of that state determines
the epoch length; the interior Θ(M) cells do.

### 4.4 The strongest O(1) shortcut, and exactly how it fails

The best a one-pass heuristic can do with O(1) state is to warm start epoch L's scan at
min(M_L, m*_{L−1}+1) — the previous stop, one cell up. `$P q28_warmstart.py 64,128,256,512 300
2026` (`q28.log`; replicated with seed 424242):

```
   n    P(ok)     se   save|ok  E[save]  cost warm  cost/epoch
  64   0.919  0.0026    0.224    0.595      3.1        7.8
 128   0.928  0.0017    0.245    0.689      3.9       12.6
 256   0.928  0.0012    0.293    0.751      5.8       23.3
 512   0.929  0.0008    0.332    0.796      9.7       47.8
seed 424242, n=512: 0.933 / 0.312 / 0.798 / 8.8 vs 43.6
```

So the shortcut is **right 93 % of the time and removes 80 % of the tests** (mean cost 47.8 → 9.7,
a 5× saving that grows with n). This is the quantitative reason the sequence *looks* exploitable —
and the reason it is not: on the other 7.1 % of epochs the true m* lies **above** the warm start, so
the scan returns a smaller concave m — or, if there is none below the start, a "meet", which
terminates the walk with the wrong δ. All failures follow a left cut (after a double cut the warm
start is exact); of them, 83 % are left cuts whose m* jumped up past the previous stop and 17 % are
double cuts whose top test the shortcut never performed
(`q28_failure_modes.log`, n = 512, 300 tapes, seed 2026). Either way the answer is
wrong, and nothing in the O(1) state reveals it: the only test that would is the one the shortcut
skipped. A 93 % rule is not an algorithm, and the residual 7 % is a subset of the "m* rises" cases
of §3 (it needs a rise of at least 2), which is what makes |S_L Δ S_{L+1}| = 0.09 M (§6.2).

### 4.5 The one strong "O(1) predictor" is a proxy, not a predictor

`$P q22_local_structure.py A 64,128,256,512 300 2026` (`q22_A.log`, replicated with seed 99 at
n = 512, 1024 in `q22_big.log`) splits P(cost = 1) by the parity of M:

```
   n   epochs  P(M even)  P(top|even)  P(top|odd)  fresh even  fresh odd
 128    23334     0.745       0.734        0.231       0.561      0.501
 256    48062     0.733       0.712        0.212       0.520      0.498
 512    97890     0.721       0.695        0.212       0.534      0.513
1024    97319     0.711       0.700        0.264       0.537      0.511   (seed 99)
```

On *fresh* windows the split is 0.5/0.5, so the trajectory split is selection, and
`$P q22_local_structure.py F 128,256,512 300 2026` (`q22_F.log`) shows what selects it:
P(M even | β = 0) = 0.894 and P(M even | β = 1) = 0.500. Parity is *absorbing* under double cuts
and flips on every left cut, so "M is even" is a 90 %-accurate proxy for "the last cut was double".
It is not usable as a rule (0.70/0.21, not 1/0), and it carries no information independent of β.

## 5. Q3 — what the expensive epochs look like, and why no certificate exists

### 5.1 The cost is concentrated in the deep scans

`$P q21_epoch_structure.py long 128,256,512 300 2026` (`q21_long.log`):

```
   n    epochs  E[cost/M]  P(c>=M/4)  P(c>=M/2)  share of tests in c>=M/4  top decile
 128    23334     0.173      0.268      0.136            0.844             0.523
 256    48062     0.159      0.255      0.117            0.837             0.523
 512    97890     0.156      0.258      0.119            0.856             0.516
```

26 % of epochs (cost ≥ M/4) carry 86 % of all tests; the top decile carries 52 %. And these are
exactly the left cuts: a left cut has E[cost/M] = 0.33–0.42 (`q21_cost.log`), so **the whole cost
of the run is one deep scan per left cut.** Any speed-up must speed up *these*.

β-sojourns (`$P q22_local_structure.py F 128,256,512 300 2026`, `q22_F.log`): at n = 512 the
β = 0 runs number 1707 with mean 32.1 (median 7, P(len=1) = 0.20, max 255) and the β = 1 runs 1578
with mean 27.3 (median 4, P(len=1) = 0.20, max 462). The sojourn law is heavy-tailed: the typical
excursion is short, but the cost is carried by the rare long ones. So long epochs arrive in long
*excursions* — the cost process has a two-state, long-memory structure even though the individual
length does not.

### 5.2 The scan has no decaying hazard

h(u) = P(cost > t+1 | cost > t) at t = u·M (`$P q24_hazard.py A`, `q24_A.log`):

```
n=512  u=0.00:0.445  u=0.02:0.856  u=0.05:0.938  u=0.10:0.967  u=0.20:0.982
       u=0.30:0.986  u=0.40:0.983  u=0.50:0.984  u=0.60:0.968  u=0.70:0.966  u=0.80:0.942
```

This is not a geometric scan that dies out: after the first test (which passes with probability
1 − h(0) = 0.555 in that run), the conditional continuation probability rises to **0.98 and stays
there** out to u = 0.8. Once an epoch is deep it stays deep — the concavity bit is positively
correlated along m. That is what makes the cost distribution bimodal and what makes "will this epoch
be long?" a question about the whole window rather than about its beginning.

### 5.3 Window contents by cost class

`q21_long.log`, n = 512, 300 tapes (class = cost/M as in §4.2):

```
class  frac  P(topfail)  E[z-drift]  P(a[R]=2)  P(aR+aR-1>=3)  E[loc2]  E[M/n]
 0     0.560     0.00       -0.78      0.239        0.294        -0.074   0.455
 1     0.118     1.00       +0.95      0.366        0.372        +0.171   0.604
 2     0.137     1.00       +1.08      0.375        0.360        +0.074   0.641
 3     0.071     1.00       +0.78      0.340        0.288        +0.052   0.636
 4     0.114     1.00       +0.41      0.417        0.431        +0.190   0.642
```

Two things to notice. (i) The drift z-score separates cost=1 from cost>1 (tautologically: for even
M the top test *is* the drift sign) but is **flat across the long classes** (0.41–1.08) — it does
not say *how* long. (ii) P(a[R] = 2) moves only 0.24 → 0.42 and is non-monotone in the class; the
2-cell local drift loc2 moves −0.07 → +0.19. There is no "one big cell at the end" signature.

The first test of a long epoch is *not* a near-tie, which is worth saying because it is the opposite
of what one expects: median |Q| (integer margin k·S_m − m·S_k, |Q| ≤ 2 means a 1-cell decision) is
**1272 at n = 512** for left-cut epochs, and only 0.2 % have |Q| ≤ 2 (`q26.log`). The failing tests
are decisive one at a time; the fragility is in *which* m's are tested, i.e. in the shift of the
window, not in the comparison. Example epoch (n = 64, seed 2026, `q26.log`): M = 60, m* = 30,
cost = 31, and all 24 top tests fail:

```
bit for m = M-23..M : 000000000000000000000000      (0 = b2 = 1 = fail)
margins |Q| at m=M..M-7 : 180 207 174 201 168 85 54 109
```

### 5.4 Where the information about "long" lives: spread linearly over the window

`$P q23_certificate.py A 128,256,512 200 2026` (`q23_A.log`) — OLS of cost/M on centred cells,
R² as a function of the number of cells w:

```
n=512   suffix : w=1 0.0060  w=2 0.0066  w=4 0.0148  w=8 0.0180  w=16 0.0231  w=32 0.0421  w=64 0.1100
        prefix : w=1 0.0055  w=2 0.0109  w=4 0.0220  w=8 0.0424  w=16 0.0776  w=32 0.1383  w=64 0.2371
        both   : w=1 0.0114  w=2 0.0174  w=4 0.0363  w=8 0.0590  w=16 0.0977  w=32 0.1708  w=64 0.3044
n=256   both   : w=1 0.0095  w=8 0.0890  w=16 0.1730  w=32 0.2995
n=128   both   : w=1 0.0198  w=8 0.1406  w=16 0.2618
```

R² ≈ (0.004–0.005)·w, i.e. **linear in w with no saturation anywhere up to w = 64** (w is capped at
n/8 by the script, which keeps enough epochs with M ≥ 2w to fit). Each cell carries the same small
amount of information; there is no scale at which the curve bends over. To explain a quarter of the
variance you need ~50 cells; to explain the rest you need the window. This is the quantitative form
of "no O(1) certificate".

Direct certificate test (`$P q23_certificate.py B 128,256,512 200 2026`, `q23_B.log`), base rate
P(cost ≥ M/4) = 0.266 at n = 512:

```
a[R]              0:0.259  1:0.223  2:0.335
a[R]+a[R-1]       0:0.262  1:0.263  2:0.244  3:0.286  4:0.322
a[L]              0:0.301  1:0.268  2:0.228
loc2=last2-first2 -4:0.187 -2:0.223 -1:0.238  0:0.261 +1:0.286 +2:0.317 +3:0.375 +4:0.380
```

Best lift: 0.375–0.380 against a base of 0.266 (1.4×), on the 3–4 % of epochs with the most extreme
value, and every one of those rules is wrong more than 60 % of the time. At n = 128 the extreme cell
reaches 0.458 (1.7×, n = 553 epochs). Nothing here is a certificate.

### 5.5 The bit field has a short correlation length relative to the scan

`$P q22_local_structure.py B`, `q22_B.log` + `q22_big.log`:

```
   n   P(c(L,m)!=c(L,m+1))  runs  runs/M   runs if independent  P(c)
  64        0.129            5.4   0.153        18.1            0.487
 128        0.101            7.8   0.115        34.5            0.484
 256        0.077           11.4   0.084        68.0            0.474
 512        0.058           16.6   0.061       135.1            0.454
1024        0.043           24.1   0.045       268.9            0.456
```

The concavity set is genuinely smoother than noise (3–11× fewer runs than independent bits) and
runs ≈ M^{0.54}: the parallel-scale structure (m = 2k and 2k+1 share the first block) does make
neighbouring scales agree. But the implied correlation length is only ℓ ≈ 1/(2·0.058) = 8.6 cells at
n = 512 (11.6 at n = 1024) against a mean scan of E[cost] = 48.4 cells (ratio 5.6, growing), and the
run structure cannot be *searched*: the bit is not monotone inside a run, so locating the last
concave m still costs one test per m. The runs statistic explains why the cost is Θ(M) rather than
Θ(M log M); it does not remove it.

## 6. Q4 — residual entropy, and the route-by-route verdict

### 6.1 Best attainable state

Taking the state an algorithm actually has for free at epoch L — M_L (hence parity), cost_{L−1},
and the local cells — the residual uncertainty in the cost class is

    H(bucket) = 1.75–1.88 bits  →  H(bucket | cost_{L−1}) = 0.49–0.82 bits
    →  H(bucket | cost_{L−1}, 4 suffix cells) = 0.45 bits   (n = 512)

and none of it is removable by a *rule*: the counterexamples of §4.3 sit inside every cell of that
state. The remaining entropy is 0.45–0.8 bits per epoch over ~0.64 n epochs, i.e. Θ(n) bits of
genuine uncertainty — the same order as the Θ(n) bits that determine the trajectory, and far less
than the 0.059 n² tests spent to resolve it.

### 6.2 Does the set move slowly?

`$P q25_nesting.py A 128,256,512 200 2026` (`q25_A.log`):

```
   n   |S|/M   lost   gained   agree   P(nested)
 128   0.475  0.162   0.170    0.843    0.065
 256   0.470  0.121   0.125    0.885    0.037
 512   0.454  0.089   0.091    0.918    0.021
```

|S_L| ≈ 0.46 M, and 9 % of it is lost plus 9 % gained per step at n = 512: the set is not nested
(P = 0.021) and the symmetric difference is 0.107 M ≈ 29 cells per step, so an incremental
maintenance of S_L costs Θ(M) per epoch = Θ(n²) total, the same as the scan. This is why the
"maintain the answer set and query max ≤ M" route buys nothing.

### 6.3 Verdict per route

1. **Two-pointer sweep.** Fails twice over: m* rises on 18–22 % of the unforced steps, and verdicts
   are not inherited (8.3 % of them change per step). The pointer itself would be cheap (TV = 4.4 n at
   n = 512, growing like n^{1.4}), which is exactly why this had to be measured rather than assumed.
2. **Exact O(1) recurrence.** Refuted by counterexamples: same state, cost 1 vs 93 (§4.3).
3. **Local certificate of a long epoch.** Refuted: R²(w) linear in w with no knee, best feature lift
   1.4×, and the first test's margin is huge (median 1272), so there is no near-tie signature to
   read off a few cells.
4. **Incremental set maintenance.** Refuted: |S_L Δ S_{L+1}| = 0.09–0.11 M per step (§6.2).
5. **Skip equal runs / binary search inside a run.** Refuted: 4.3–7.7 % flip per m-step (n ≥ 256),
   run count M^{0.54}, no monotone structure within a run, and skipping requires reading what is
   skipped.
6. **Predict the trajectory instead of testing it.** β is nearly deterministic per step
   (P(β=1|β=1) = 0.964, P(β=1|β=0) = 0.026, H(β|prev) = 0.197 bits at n = 512, `q22_D.log`), so the
   *sequence of windows* is cheap to guess — but "guess + verify" is the scan, and the guess is wrong
   on 2.6–3.6 % of steps, which corrupts L\* and δ.

### 6.4 What a sub-quadratic algorithm would have to do instead

Not a local or monotone method: the measurements above close locality (R² grows linearly in w),
monotonicity (18–22 % violations on the free steps) and inheritance (8.3 % of verdicts flip per step). What survives
untested is a *global algebraic* method: the test is the orientation predicate
`2P[L+k] ≥ P[L+2k] + P[L]`, i.e. "is the point at L+k below the chord from L to L+2k"; a fast method
would have to preprocess the prefix-sum polygon so that "largest stride k ≤ M/2 whose midpoint point
is below the chord" is answered without examining the strides one by one. These measurements say
nothing against such a method — they say that the answer cannot be read off local data, cannot be
carried over from the previous epoch, and is not monotone in L, m or R, so it would have to be a
genuine range-searching structure on the polygon, not a sweep.

The one framing that the data does support, and that is independent of any algorithmic claim:

> The walk's trajectory is described by Θ(n) bits (0.64 n nearly deterministic bits, H = 0.2 bits
> each) and its cut positions move a total distance of only 4.4 n at n = 512 (n^{1.4} overall), yet
> the rule spends 0.059 n² tests to produce them. The information in the answer is Θ(n); the cost is
> verification, not information.

## 7. Files, commands, seeds

`$P = C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe`, run in
`C:/Users/shini/poetry_axis/paper/probe_ologn2/`. New files only; no .tex touched.

```
q21_epoch_structure.py   per-epoch records + modes mono|persist|cost|long (also ent, unused here)
  $P q21_epoch_structure.py mono    64,128,256,512 400 2026     -> q21_mono.log/.json
  $P q21_epoch_structure.py persist 64,128,256,512 400 2026     -> q21_persist.log
  $P q21_epoch_structure.py cost    64,128,256,512 400 2026     -> q21_cost.log
  $P q21_epoch_structure.py long    128,256,512    300 2026     -> q21_long.log
  $P q21_epoch_structure.py mono    64,128,256,512 400 424242   (replication, same log file)

q22_local_structure.py   modes A|B|C|D|E|F
  $P q22_local_structure.py A 64,128,256,512 300 2026           -> q22_A.log
  $P q22_local_structure.py B 64,128,256,512 300 2026           -> q22_B.log
  $P q22_local_structure.py C 64,128,256,512 300 2026           -> q22_C.log
  $P q22_local_structure.py D 64,128,256,512 300 2026           -> q22_D.log
  $P q22_local_structure.py E 64,128,256,512 300 2026           -> q22_E.log
  $P q22_local_structure.py F 128,256,512    300 2026           -> q22_F.log
  $P q22_local_structure.py A|B 512,1024 150 99                 -> q22_big.log

q23_certificate.py       modes A (R^2 vs w), B (feature lift), C (D = M-m*_prev)
  $P q23_certificate.py A|B|C 128,256,512 200 2026              -> q23_A/B/C.log
q24_hazard.py            modes A (hazard) B (cost~d) C (runs) D (identities) E (entropy)
  $P q24_hazard.py A|B|C|D|E 128,256,512 200 2026               -> q24_*.log
q25_nesting.py           modes A (nesting) B (counterexamples) C (redundancy)
  $P q25_nesting.py A|B|C 128,256,512 200 2026                  -> q25_*.log
q26_margin.py
  $P q26_margin.py 128,256,512 150 2026                         -> q26.log
q28_warmstart.py         the best O(1)-state shortcut and its failure rate
  $P q28_warmstart.py 64,128,256,512 300 2026                   -> q28.log
  $P q28_warmstart.py 128,512 300 424242                       (replication, overwrites q28.log)

q27_chord.log            (chord / block identities, 0/16700 and 0/16340)
```

Seeds: 2026 (all main tables), 424242 (monotonicity replication, identical conclusions), 99 with
n = 512,1024 (parity/flip-rate controls in `q22_big.log`), 77/555/… unchanged from
`REPORT_LOWER.md` (not used here).

## 8. Caveats

* Everything is measured on the ensemble of `REPORT_LOWER.md` (uniform two-column tapes). The
  negative is about *this* ensemble and *these* structural properties; a different tape family, or a
  different rule (b1 ≠ 0), could behave differently.
* "No O(1) certificate" is established as "R²(w) grows linearly, best feature lift 1.4×, and
  counterexamples exist for the O(1) states tested" — not as a theorem about all O(1) functions of
  the window. The counterexample of §4.3 is the strongest form: it holds the whole O(1) state fixed
  and varies only the interior.
* The cost process is *predictable in bulk* (§5.1: sojourns of 27–32 epochs, cost ≈ 0.3 M inside a
  left-cut excursion). Nothing in this report rules out an algorithm that exploits that at the level
  of a block of epochs; what is ruled out is doing it epoch by epoch from local data.
