# The contrast curve: R² against the read budget, cheap rule vs expensive rule

`_contrast.py`, `_contrast_check.py`, `_contrast_plot.py`, `figs/contrast_curve.png`.
`n = 256`, 2000 tapes, seed 2026, one predictor and one read counter for both rules.
No `.tex` was touched.

---

## 1. What the curve is

A predictor is allotted **B cell reads** and spends them on the exact simulator of
the rule:

* if the tape's run costs **≤ B** reads, the output is determined and the prediction
  is that exact output;
* otherwise the predictor stops at the last step whose cumulative read count is ≤ B,
  sees the walk state there (window statistics, counters, both decision bits — all
  available at zero marginal reads from the accumulators, the same family used for
  the four-rule curve), and predicts δ with the same out-of-fold
  `HistGradientBoostingRegressor(max_iter=80, lr=.05, depth=3, leaf=20, l2=1.0)`.

So R²(B) is the fraction of the output variance a machine with B reads recovers. It
reaches 1 exactly when every tape in the sample is decided.

This is deliberately *not* the survival-conditioned curve of `_phase_rule.py` /
`_phase_arr2.py`, which drops the tapes that have already stopped. That curve can
never reach 1 — at any step s the surviving tapes are precisely the ones whose
answer is not yet known — so it cannot answer where a curve meets 1. It also loses
sample as s grows (495 → 393 → 365 meetings in `_phase_lr_fast256.log`), which is
what makes that curve turn down after s = 128. Here the meeting set is constant
(646 and 495) at every budget.

### Costs measured by the same counter

| rule | meetings / tapes | reads mean | median | min | max |
|---|---:|---:|---:|---:|---:|
| canonical | 646 / 2000 (32.3 %) | 1317 = **5.15 n** | 5.05 n | 4.99 n | 6.10 n |
| longest-run | 495 / 2000 (24.8 %) | 59 267 = **231.5 n** = 0.904 n² | 163.3 n = 0.638 n² | 99.2 n | 601.9 n = 2.35 n² |

The canonical figure reproduces Table 1 (5.32 n there, 5.15 n here; the difference
is the meeting subsample). The paper's 0.62 n² for the longest-run rule is
reproduced as the **median** (0.638 n²); the mean over meeting tapes is 0.90 n²
because the cost is heavy-tailed (90th percentile of steps 10 259 = 0.157 n²).

The canonical simulator's opening read of the tape costs 3.5 n reads
(`Win.__init__` reads every cell twice for Rr), so no point exists below 3.6 n: at
B = 1 n the machine has not finished looking at the tape. The longest-run curve
starts at 2 n for the same reason — its first snapshot is one step past an opening
scan of n + n/2 reads.

### Why the step grid had to become a read grid

The obvious extension of the existing longest-run curve — add s = 192, 256, 384,
512 steps at n = 256 — does not reach the exact cost. At n = 256 the fast
longest-run simulator spends about 8.7 reads per step, so 512 steps is 72 n reads,
**45 per cent of 0.62 n² = 159 n**. On the step grid the curve is still climbing at
s = 512 and covers only the left half of the figure. The exact cost is met near
s = 4 100 steps (measured 157.8 n at s = 4 096), which is why the budgets here are
read budgets up to 640 n rather than a step grid ending at 512.

---

## 2. The two curves

| B/n | canonical R² | decided | longest-run R² | decided |
|---:|---:|---:|---:|---:|
| 2 | -- | -- | 0.376 | 0.00 |
| 2.5 | -- | -- | 0.449 | 0.00 |
| 3 | -- | -- | 0.459 | 0.00 |
| 3.5 | -- | -- | 0.444 | 0.00 |
| 3.6 | 0.544 | 0.00 | -- | -- |
| 3.8 | 0.655 | 0.00 | -- | -- |
| 4 | 0.753 | 0.00 | 0.448 | 0.00 |
| 4.2 | 0.834 | 0.00 | -- | -- |
| 4.4 | 0.880 | 0.00 | -- | -- |
| 4.6 | 0.921 | 0.00 | -- | -- |
| **4.8** | **0.959** | 0.00 | -- | -- |
| 5 | 0.981 | 0.35 | 0.400 | 0.00 |
| 5.15 | 0.977 | 0.66 | -- | -- |
| 5.3 | 0.983 | 0.81 | -- | -- |
| 5.5 | 0.989 | 0.92 | -- | -- |
| 5.7 | 0.986 | 0.98 | -- | -- |
| 6 | 0.995 | 1.00 | 0.419 | 0.00 |
| **6.3** | **1.000** | 1.00 | -- | -- |
| 8 | 1.000 | 1.00 | 0.455 | 0.00 |
| 13 | -- | -- | 0.487 | 0.00 |
| 20 | -- | -- | 0.546 | 0.00 |
| 32 | -- | -- | 0.550 | 0.00 |
| 50 | -- | -- | 0.575 | 0.00 |
| 63 | -- | -- | 0.638 | 0.00 |
| 80 | -- | -- | 0.674 | 0.00 |
| 100 | -- | -- | 0.741 | 0.12 |
| 126 | -- | -- | 0.763 | 0.39 |
| **158.7 = 0.62 n** | -- | -- | **0.751** | 0.48 |
| 200 | -- | -- | 0.809 | 0.55 |
| 252 | -- | -- | 0.846 | 0.62 |
| 320 | -- | -- | 0.877 | 0.71 |
| 400 | -- | -- | 0.943 | 0.81 |
| **500** | -- | -- | **0.976** | 0.94 |
| **640** | -- | -- | **1.000** | 1.00 |

Tapes: 2000 for both. Meetings: 646 (canonical), 495 (longest-run), constant at
every budget. `decided` is the fraction of meetings whose run finished inside B.
Out-of-fold standard errors on the canonical points are .004–.081 (largest at
5.7 n, where 12 undecided tapes remain); on the longest-run points .012–.262
(largest on the plateau). The fitted booster on the *undecided* tapes alone
(`R2learn` in `_data/contrast_*.json`) equals the total while nothing is decided
and falls below it once the decided tapes carry the curve. The booster is not
uniformly the better model on the tail: RidgeCV on the same features scores .54
against the booster's .35 at 400 n, so the level of the longest-run plateau is a
property of the feature set, not of the booster.

### Meeting points

| | first B with R² ≥ 0.95 | first B with R² = 1.000 | exact cost |
|---|---:|---:|---:|
| canonical | **4.8 n** (R² = .959) | **6.3 n** | 5.3 n quoted, 5.15 n measured |
| longest-run | **500 n** (R² = .976) | **640 n** | 0.62 n² = 159 n quoted, 163 n median |
| ratio | **104×** | **102×** | 30× |

Tapes and meetings behind every row are 2000 and 646 / 495. The 0.95 crossing is
bracketed: canonical R² = .921 at 4.6 n and .959 at 4.8 n; longest-run .943 at
400 n and .976 at 500 n.

---

## 3. Does the data say what the classification says?

**Yes on the order, no on the constant.**

* The two curves do separate by two orders of magnitude, and each reaches 1.000 at a
  budget that grows with n: the cheap one at 6.3 n, the expensive one at 640 n =
  2.5 n². Read as an order statement — Θ(n) against Θ(n²) — the figure is exactly
  the picture the classification predicts.
* The factor is **102**, not n = 256. In units of n the measured ratio is **0.40 n**,
  so it grows linearly in n but with a constant near 0.4 rather than 1. "A factor of
  about n apart" is right in order and overstated by 2.5 at this size. If the
  sentence is kept, it should read "a factor of about 0.4 n".
* The quoted exact costs (5.3 n and 0.62 n²) are a factor of **30** apart, and at
  its own exact cost the expensive rule is only at R² = **.751**: half its meetings
  (48.5 %) have not finished by 159 n. "The curve reaches 1 at 0.62 n²" is not what
  the data show — 0.62 n² is the *median* cost, and reaching 1 needs the maximum of
  the sample, which is 2.35 n².
* The cheap curve is the opposite case: its cost is tight (4.99 n to 6.10 n, a
  spread of 22 %), so R² = 1 arrives within a fifth of a n of the mean cost, and
  the curve is already at .983 at the quoted 5.3 n.

**The 0.95 points are the asymmetry that matters, and they are not symmetric.**
The canonical curve reaches .95 at **4.8 n = 0.93 × its measured mean cost**
(0.91 × the quoted 5.3 n) — that is, essentially at its exact cost, not "well
before" it. The longest-run curve reaches .95 at **500 n = 3.1 × its quoted exact
cost** (2.2 × its measured mean). The threefold overshoot is the dispersion of the
run: a budget large enough to decide half the tapes leaves the other half — the
long ones, whose δ carries much of the variance — to a predictor that scores only
R² = .27 on them.

**Is the simulator exact at 5.3 n?** Yes, by construction and by test. Every tape
whose run costs ≤ B is scored with its exact output, so at B = 5.3 n the 81.3 % of
meetings that have finished contribute zero error; the residual 1.7 % of variance
comes entirely from the 18.7 % still running, and `R2learn` on that subset is .93.
The simulator itself was checked independently: `_contrast_check.py canon 256 120`
gives **120/120** agreement on output **and on read count** with
`code/bthm0/_bz_canon2.py walk_inc`, and **0** mismatches against a from-scratch
double-cut oracle. (`_bz_canon2.py`'s own `walk_exact_bf` asserts a *single* cut
where the paper specifies a double cut when both headers are on the boundary; that
in-file oracle now fails on the first tape of a 40-tape run, returning 0.0039 where
`walk_inc` returns None. `walk_inc` is the one that matches the paper and the one
reproduced here.) For the longest-run rule,
`_contrast_check.py longrun 256 40` checks both bits at every step against brute
force — **169 122 bit checks, 0 mismatches** — and the output against a
brute-force walk driven by those bits.

### Caveats

* R² = 1.000 is reached when the last *sampled* tape is decided, so 640 n is the
  maximum of 495 meeting costs, not a population constant; with more tapes it would
  grow. The 0.95 point is the stable statistic.
* R² is computed over tapes that produce an output (the paper's convention), so the
  target distribution differs between the two rules; this affects the level of the
  curves, not the budgets at which they saturate.
* One size only (n = 256). The claim that the ratio grows like 0.40 n is an
  extrapolation from the two cost laws, not a second measurement.

---

## 4. Paper-ready wording

> Figure~\ref{fig:contrast} puts the two sides of the classification on one axis.
> Both curves are built the same way: a machine is given $B$ cell reads, runs the
> exact simulator, returns the output when the run costs at most $B$, and otherwise
> predicts $\delta$ from the walk state with the same out-of-fold booster used
> above. At $n=256$ over 2000 tapes the canonical curve crosses $R^2=.95$ at
> $4.8\,n$ reads and reaches $1.000$ at $6.3\,n$, when the last of 646
> output-producing tapes is decided; the longest-run curve crosses $.95$ at
> $500\,n$ and reaches $1.000$ at $640\,n$, on 495 meetings. The two saturation
> budgets differ by a factor of $102$, which is $0.40\,n$: the ratio grows linearly
> in $n$, as the classification requires, with a constant near $0.4$ rather than
> $1$. The exact costs alone differ by a factor of $30$, and the extra factor is
> dispersion --- the canonical run costs between $4.99\,n$ and $6.10\,n$ reads,
> while the longest-run run costs between $99\,n$ and $602\,n$, so a budget that
> decides half of its tapes still leaves the expensive half to a predictor that
> scores $.27$ on them. At its own exact cost of $0.62\,n^2$ the expensive rule is
> at $R^2=.751$, not $1$; $0.62\,n^2$ is the median cost, and $1$ arrives only at
> the tail.

---

## 5. Reproduce

```bash
P=C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe
$P _contrast.py canon   256 2000 2026   # 3 s,  writes _data/contrast_canon_n256_N2000_s2026.json
$P _contrast.py longrun 256 2000 2026   # 33 s, writes _data/contrast_longrun_n256_N2000_s2026.json
$P _contrast_check.py canon   256 120   # simulator vs _bz_canon2 and vs brute force
$P _contrast_check.py longrun 256 40    # both bits at every step vs brute force
$P _contrast_plot.py                    # figs/contrast_curve.png
```

Logs: `_contrast_canon256.log`, `_contrast_longrun256.log`.
