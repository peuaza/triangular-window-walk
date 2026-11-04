# kappa_n: a finite-n statement, with two-sided intervals and an explicit selection account

*Response to the objection that the ceiling is estimated at one n, bounded on one side, and
uncorrected for the search over budgets, folds, rules and feature sets. No .tex file was touched.
Every number below comes from the scripts in §8; the raw outputs are in `_data/kap/`.*

---

## 0. What the objection is right about

The current §7 text says

> "κ ≥ .56 (one-sided 95 per cent, selection-adjusted)"

from `.579`/`.583` at n=256 (8000 tapes, 1919/1939 meetings) and reports `.48` with a bound of
`.39` at n=512 (900 tapes, 198 meetings). Four things are wrong with reading that as a statement
about κ:

1. **One-sided.** A one-sided bound cannot be contradicted from above, so it cannot distinguish
   "κ = .6" from "κ = .9". A referee is entitled to ask what values the data exclude.
2. **One n.** Nothing in the construction says the number at n=256 is the number at n=128 or 512;
   the paper itself reports .48 at 512 and then reads the difference qualitatively.
3. **Selection over the budget grid is handled, selection over everything else is not.** The
   existing bootstrap (`_boot_sel4.py`) does take the maximum over the budget grid inside each
   resample, which is right. It does not account for the fold realisation, and it cannot account
   for rules or feature sets.
4. **The n=512 row is confounded by training size, and this is the largest single error in the
   section.** 198 meetings train a booster on ~158 observations per fold. At the same n with 1867
   meetings the estimate is not .481 but .550 (§5). The sentence "the drop from .56 is about three
   standard errors of that sample" is also not reproducible: the bootstrap sd of that row is .060,
   so .08 is 1.3 sd; the `se` column of the JSON it appears to use is std(per-fold R²)/√5, a
   fold-to-fold spread and not the sampling sd of the estimate (it is .055 at the selected budget,
   and the .02 values in the same column belong to other budgets).

Everything below is my own run at 8000 tapes per size for each of two seeds, so the rows for the
three sizes are comparable in the only way that matters for a trained predictor — the number of
output-producing tapes it is trained on.

## 1. The estimand

Fix n. Let P_n be the ensemble of §7 (two independent uniform ⌈n/2⌉-subsets, one per column, and
the tape their average), Δ the output, ⊥ on tapes that halt or repeat, and P_n^out = P_n(· | Δ ≠ ⊥)
the law of the output conditional on producing one. Write R²(g) = 1 − E[(Δ − g)²]/Var(Δ), both
moments under P_n^out.

Let B be a finite set of step budgets and N_n the number of tapes drawn. Define

  **C_n(B)** = { the predictor returned by the fixed procedure: gradient boosting
  (`HistGradientBoostingRegressor`, max_iter 80, learning rate .05, max depth 3, min samples per
  leaf 20, L2 1.0, `random_state` 0) on the six prefix features (position, window length, region
  length, cuts, step, running mean) read off the first s steps of the trajectory, for some s ∈ B,
  fitted on 4/5 of the N_n tapes, with its population R² estimated by pooling the five held-out
  folds },

  **κ_n(B) = sup { R²(g) : g ∈ C_n(B) }.**

A predictor trained on all N_n tapes rather than 4/5 of them would do slightly better, so the
out-of-fold numbers below are mildly conservative for the deployed predictor; the class they
describe is the 80 per cent one, and that is stated rather than assumed away.

Three things this definition is careful about, each of which the objection turns on:

* **It is in-class.** κ_n(B) is the best population R² *inside the class the paper actually
  searches*, not the best over all predictors and not the limit of anything. It is a property of
  a procedure (including its training-set size), not of the rule alone.
* **It is finite-n.** The population is the ensemble at size n, the predictor is trained on N_n
  tapes, and the budget is a finite set. Nothing is extrapolated in n.
* **It depends on B.** A larger B gives a larger κ_n(B), because the class is bigger. This is not
  bias; it is a different estimand. So every number below carries its budget set and its largest
  budget in reads per cell. **Comparing κ at two n with different budget sets is not a like-for-like
  comparison**, and §4 gives the fixed-B version separately.

What is estimated: the pooled out-of-fold R² at each s, θ̂_s, and θ̂ = max_s θ̂_s. Because the same
m tapes drive every θ̂_s, the θ̂_s are strongly correlated (correlation with the argmax budget is
.78–.81 at the neighbouring budget; the participation ratio of the 14×14 correlation matrix is 1.86,
i.e. the 14 budgets carry about two independent estimates). θ̂ is therefore barely more variable
than a single θ̂_s, but it is biased up by the selection; §6 quantifies both.

## 2. Procedure

* Tapes: 8000 per size and seed, generated exactly as in `_phase_arr2.py`
  (`default_rng(seed)`, two independent `choice(n, n//2, replace=False)` columns). Meeting counts
  reproduce the paper's: 1919 at n=256 seed 2026, 1939 at seed 90210.
* Budget survival. A tape appears in the budget-s sample only if its trajectory reaches step s.
  Budgets are kept only when **every** output-producing tape reaches them, so all budgets score
  one population. This drops s ≥ 96 at n=128, s ≥ 192 at n=256, s ≥ 384 at n=512. The dropped
  budgets all score below the retained maximum, so the rule does not buy the answer for free:
  n=128 gives .505 at s=96 (1687 of 2068 tapes) and .479 at s=128; n=256 gives .486 at s=192
  (1587 of 1919) and .433 at s=256; n=512 gives .414 at s=384 (1522 of 1867) and .356 at s=512.
* Out-of-fold predictions: 5-fold `KFold(shuffle=True, random_state=seed)`, as in the paper.
* Interval: B = 4000 bootstrap resamples of the m output-producing tapes, with replacement, the
  **same resample index used at every budget**; inside each resample the pooled R² is recomputed
  at every budget from the stored out-of-fold predictions and **the maximum is retaken**; the
  interval is the 2.5/97.5 percentile of that maximum. 25 different 5-fold splits are used, one per
  bootstrap replicate (`1000r + 17`), so the fold realisation is inside the interval as well.
* Reported: percentile interval; BCa endpoints as a check; the one-sided 5th percentile as a link
  to the old number.

## 3. The table

**Table 1. κ_n inside the searched grid (8000 tapes per row, two seeds).**

| n | tapes | meetings | budget grid (steps) | best s | reads/cell | κ̂_n | 95% two-sided | seed |
|---|---|---|---|---|---|---|---|---|
| 128 | 8000 | 2068 | 1–64 (12 pts) | 64 | 28.6 n | **.619** | **[.594, .645]** | 2026 |
| 128 | 8000 | 1976 | 1–64 (12 pts) | 48 | 25.0 n | **.617** | **[.592, .643]** | 90210 |
| 256 | 8000 | 1919 | 1–128 (14 pts) | 96 | 45.4 n | **.582** | **[.559, .611]** | 2026 |
| 256 | 8000 | 1939 | 1–128 (14 pts) | 96 | 44.3 n | **.584** | **[.550, .607]** | 90210 |
| 512 | 8000 | 1867 | 1–256 (16 pts) | 192 | 83.4 n | **.550** | **[.523, .579]** | 2026 |
| 512 | 8000 | 1884 | 1–256 (16 pts) | 192 | 79.9 n | **.550** | **[.522, .582]** | 90210 |

Notes. The point is the pooled out-of-fold R² at the best budget, i.e. the maximum over the grid;
the interval is the percentile bootstrap of that maximum. BCa endpoints differ from the reported
ones by at most .012 (.593/.643, .591/.642, .548/.603, .559/.615, .519/.575, .517/.576). The
one-sided 5th percentile, which is what the paper reports, is .598, .596, .563, .555, .528, .527 —
the n=256 seed-2026 value reproduces the paper's .563 exactly, and .555 is within .006 of its .561.

Read as a bound: the lower endpoint of a two-sided 95 per cent interval is also a one-sided
97.5 per cent lower bound, so the strongest defensible statement of the old form is **κ_256 ≥ .559**
and **κ_512 ≥ .523**, conditional on this rule and this feature set.

**These numbers are uncorrected for the search over rules and over feature sets, and no such
correction is possible with the data in hand** (§6c,d): one rule and one feature set were measured,
so there is nothing to maximise over. Every number in this note carries that label.

**Table 2. A like-for-like class: the same step budgets at every n.**

| n | seed | s=32 (read/cell) | 95% | s=48 (read/cell) | 95% |
|---|---|---|---|---|---|
| 128 | 2026 | .544 (20.6 n) | [.516, .576] | .614 (25.7 n) | [.589, .642] |
| 128 | 90210 | .534 (20.0 n) | [.507, .569] | .617 (25.0 n) | [.588, .642] |
| 256 | 2026 | .443 (22.9 n) | [.411, .477] | .488 (30.8 n) | [.455, .519] |
| 256 | 90210 | .408 (22.2 n) | [.371, .438] | .460 (29.9 n) | [.427, .493] |
| 512 | 2026 | .373 (23.6 n) | [.343, .411] | .391 (33.3 n) | [.358, .423] |
| 512 | 90210 | .378 (22.5 n) | [.353, .419] | .397 (31.8 n) | [.367, .435] |

At s=48 the three sizes cost 25.7n, 30.8n and 33.3n reads (s=32: 20.6n, 22.9n, 23.6n), so the
budget is matched to within 30 per cent, and the budget family is the same set {1,…,48} in steps.
Each cell is the sup over that family, which for all three n is attained at the largest member. All
rows are the 8000-tape samples of Table 1, with the meeting counts given there.

## 4. The trend

**At a fixed budget, κ_n falls.** At s=48: .614 → .488 → .391. The two-sample bootstrap
(independent resamples of the two samples' draws) gives

| difference | estimate | 95% interval |
|---|---|---|
| κ_128 − κ_256 (s=48) | +.126 | [+.088, +.171] |
| κ_256 − κ_512 (s=48) | +.097 | [+.050, +.141] |
| κ_128 − κ_512 (s=48) | +.223 | [+.183, +.267] |
| κ_128 − κ_256 (s=32) | +.100 | [+.059, +.148] |
| κ_256 − κ_512 (s=32) | +.070 | [+.019, +.114] |

Neither difference interval contains zero. The decline is not a null result at a fixed budget.

**At each n's own best budget, κ̂_n falls, but the differences are inside the intervals.** .619 →
.582 → .550; κ_128 − κ_256 = +.038 [−.003, +.071] and κ_256 − κ_512 = +.032 [−.005, +.073], both
touching zero; the whole range +.070 [+.031, +.106] is outside. But these three numbers are not a
like-for-like comparison — the class at n=128 stops at 28.6n reads and at n=512 extends to 83.4n —
so the flatness is an artefact of letting the class grow with n. **The like-for-like answer is that
κ falls with n at fixed budget and is never observed to rise.**

One more thing the table shows: at n=128 the maximum sits on the boundary of the surviving grid
(s=64 is the last budget at which all 2068 tapes exist). κ_128 is therefore censored from above by
the survival rule and .619 is a lower bound for a class that extends further; it is not a peak.

## 5. The n=512 row was measuring the training set, not n

Same rule, same n, different numbers of tapes:

| n | tapes | meetings | seed | κ̂ | 95% two-sided | fold-split sd |
|---|---|---|---|---|---|---|
| 512 | 900 | 198 | 2026 | .481 | [.380, .613] | .017 |
| 512 | 4000 | 936 | 2026 | .513 | [.478, .562] | .007 |
| 512 | 4000 | 933 | 90210 | .518 | [.477, .567] | .007 |
| 512 | 8000 | 1867 | 2026 | **.550** | [.523, .579] | .003 |
| 128 | 800 | 213 | 2026 | .515 | [.431, .629] | .020 |
| 128 | 8000 | 2068 | 2026 | **.619** | [.594, .645] | .002 |

Going from 198 to 1867 meetings is worth +.07 at n=512; going from 213 to 2068 is worth +.10 at
n=128. Both are larger than the .03 "drop" the paper attributes to n. At 198 meetings the
two-sided interval is [.380, .613] (BCa [.316, .563]), .23 wide, so the row cannot separate .48
from .35, and the upward selection over 16 budgets at sd .058 is itself of order .05 — the same
size as the gap being discussed. **No comparison of κ across n is interpretable unless the meeting
count is matched**, which is why Table 1 uses 8000 tapes everywhere (1867–2068 meetings).

## 6. Selection: what is corrected, what is not, and what it costs

**(a) Budget grid — corrected.** The maximum is retaken inside every bootstrap resample, so the
interval is a percentile interval for the max-over-grid functional, not for a single budget.
Cost of taking the max rather than reporting the argmax budget: both endpoints move up, by .003 and
.007, and the width is unchanged (.053 vs .056 at n=256) — the budgets are too correlated for the
grid to be an expensive search. The *extent* of the grid is the expensive part, and it is not
selection bias: extending the budget from 31n to 45n reads at n=256 moves the estimate from .488 to
.582 because a predictor with a bigger budget is a different, better predictor. That is why the
class boundary is printed in Table 1.

**(b) Five folds — corrected two ways, and it is not a fold-max.** The reported statistic is the
*pooled* out-of-fold R², not the maximum over the five held-out folds, so no fold selection enters
the point estimate. If a reader assumed a fold-max, the number would be larger and worse: at the
selected budget the fold-max is .634, .654, .610, .621, .581, .564 for the six rows (i.e. .014–.037
above the pooled value) and its bootstrap interval at n=256 is [.580, .656]. We do not use it.
What *is* in the interval is the fold-split realisation: 25 different 5-fold splits, one per
bootstrap replicate, with a spread of .0018–.0039; adding it in quadrature changes the half-width
from .0265 to .0268 at n=256. The estimator's dependence on which tapes land in which fold is
therefore accounted for and immaterial at 8000 tapes. At 198 meetings it is not immaterial (split
sd .017 against a tape sd of .058).

**(c) Rules — not corrected, and not correctable with these data.** One rule was measured, the
longest-run rule (with the O(n) b₁ witness). There is no second rule's out-of-fold predictions to
take a maximum over, so the number is **conditional on that rule** and is labelled as such. It is
not a bound over the four candidates of §7. A Bonferroni-type sensitivity for a search over K
equally good alternatives, which is *not applied* to the numbers above, would subtract about
σ√(2 ln K) from the lower endpoint (σ = .0135 at n=256): .022 for K=4, .029 for K=10, .034 for
K=25, i.e. κ_256 ≥ .537, .530, .525. The ceiling claim survives that; the exact constant does not.

**(d) Feature sets — not corrected, and not correctable with these data.** One feature set (six
hand-chosen prefix features) was measured. The only second learner available is ridge on the same
features, whose best budget gives .542 at n=256; it never wins in any of the 4000 resamples, so the
maximum over {booster, ridge} × grid is identical to the booster's and the correction over that
pair is exactly zero. That is evidence about one alternative, not about the feature-set space.

**(e) Number of prefixes** is the budget grid, item (a).

## 7. Paper-ready wording

> We estimate κ_n, the largest population R² attainable inside the class we search at size n, on
> 8000 tapes per size, of which 2068, 1919 and 1867 produce an output at n = 128, 256, 512. The
> class is the gradient booster of §7 on the six prefix features, over step budgets up to the
> largest at which every output-producing tape still exists (28.6n, 45.4n and 83.4n cell reads).
> Out-of-fold, κ̂_n = .619 [.594, .645], .582 [.559, .611] and .550 [.523, .579], two-sided 95 per
> cent, from a bootstrap that resamples output-producing tapes, recomputes the out-of-fold R² at
> every budget inside each resample, retakes the maximum over the budget grid, and redraws the
> five-fold split. The intervals exclude R² = 1 by a wide margin: the ceiling is not an artefact of
> a starved predictor, and the strongest one-sided statement of the old form is κ_256 ≥ .559
> (97.5 per cent). They do not locate κ. The estimate moves with the training sample — at n = 512,
> 198 meetings give .481 and 1867 give .550 — so no row here should be read as an approach to a
> limit. On a common budget of 48 steps the three values are .614, .488 and .391 and the falls are
> outside the intervals ([+.088, +.171], [+.050, +.141]); the value therefore falls with n at fixed
> budget, and appears flat only when each n is allowed its own largest budget. The numbers are
> conditional on the longest-run rule and on this feature set: no correction over rules or over
> feature sets is possible with the data we have, and none is claimed.

If a shorter replacement for the display at eq. (κ ≥ .56) is wanted:

> κ_256 ∈ [.559, .611] and κ_512 ∈ [.523, .579] (two-sided 95 per cent, selection-adjusted over the
> budget grid and the fold split; conditional on the longest-run rule and the six prefix features).

## 8. Commands, seeds, files

```
P=C:/Users/shini/poetry_axis/intern_decision/.venv/Scripts/python.exe
cd C:/Users/shini/poetry_axis/paper

# budget-curve data (walk + features + tape index), 8 processes
$P _kappa_gen.py 128 8000 2026  _data/kap/n128_s2026
$P _kappa_gen.py 128 8000 90210 _data/kap/n128_s90210
$P _kappa_gen.py 256 8000 2026  _data/kap/n256_s2026
$P _kappa_gen.py 256 8000 90210 _data/kap/n256_s90210
$P _kappa_gen.py 512 8000 2026  _data/kap/n512b_s2026
$P _kappa_gen.py 512 8000 90210 _data/kap/n512b_s90210
# training-size rows
$P _kappa_gen.py 128 800  2026 _data/kap/n128_N800_s2026
$P _kappa_gen.py 512 900  2026 _data/kap/n512_N900_s2026
$P _kappa_gen.py 512 4000 2026 _data/kap/n512_s2026
$P _kappa_gen.py 512 4000 90210 _data/kap/n512_s90210

# point estimates, 25 fold splits, 4000-resample bootstrap
for t in n128_s2026 n128_s90210 n256_s2026 n256_s90210 n512b_s2026 n512b_s90210; do
  $P _kappa_boot.py _data/kap/$t 4000 25 $t; done
for t in n128_N800_s2026 n512_N900_s2026 n512_s2026 n512_s90210; do
  $P _kappa_boot.py _data/kap/$t 2000 15 $t; done

$P _kappa_final.py     # tables 1-3 and the difference intervals
$P _kappa_extra.py     # dropped budgets, variance arithmetic, search sensitivity
$P _kappa_sd.py _data/kap/n256_s2026   # bootstrap vs jackknife vs half-split sd
```

Seeds: tapes 2026 and 90210; canonical fold split `random_state = tape seed`; the other 24 splits
`1000r + 17`; bootstrap RNG 7; permutation null RNG 12345.

Files written: `_kappa_gen.py`, `_kappa_boot.py`, `_kappa_final.py`, `_kappa_extra.py`,
`_kappa_sd.py` (new), `_data/kap/<tag>/{meta,all,s*}.npz`, `_data/kap/<tag>.json`,
`_data/kap/<tag>_draws.npz` (bootstrap draws, so every interval in this note is recomputable).

Reproduction checks against the paper: 1919 meetings at n=256 seed 2026 ✓; 1939 at seed 90210 ✓;
198 at n=512 with 900 tapes ✓; fold-mean R² at n=256 .5794 (s=96) and .5796 (s=128) against the
paper's .579/.580 ✓; one-sided 5th percentile .563 against .563 ✓; pooled n=512 with 900 tapes
.481 against the paper's .48 ✓. The pooled and mean-of-folds functionals differ by .002 at n=256 and
by .011 at 198 meetings; Table 1 uses the pooled one, which is the functional the bootstrap
resamples.

Note on the stored predictions: `_data/s1big` and `_data/s2` (both n=256) each carry an
`oof_s384.npz`/`oof_s512.npz` pair identical to `_data/seed2`'s (156 and 150 rows, reads 33818.4
and 34274.1) — those budgets are outside the n=256 grid and the files are copies of the n=512 run.
The old bootstrap survived this only because it filtered to the modal sample length. That is one
reason the numbers above were regenerated rather than reused; another is that the fold-split refits
need the feature matrix, which the old npz files do not carry.

## 9. What this does not establish

* No limit. Nothing here says κ_n converges; the training-size table says the opposite — the
  estimate is still moving with m at fixed n, and the largest budgets are censored by survival.
* No bound over rules or feature sets. The number is conditional on the longest-run rule and the
  six prefix features, and there is no correction for either.
* No like-for-like statement about n at the *best* budget, because the best budget grows with n.
  The like-for-like statement is the fixed-budget one, and it is a decline.
* The intervals are percentile bootstrap intervals for a max-over-grid functional; BCa shifts them
  by up to .012 and would widen the 198-meeting row substantially. Neither changes a conclusion
  above.
