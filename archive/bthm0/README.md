# BTHM-0: scripts for the merged paper

Each script reproduces one table or paragraph of `bthm_merged.tex`.

| script | produces | paper location |
|---|---|---|
| `_onetable.py` | the headline run: in-sample and 5-fold CV R² for affine / cubic / boosting on ac₁ and on twenty features, 4000 tapes at n=256 | Table 1 |
| `_by_n.py` | the same measurement at n = 128, 256, 512, 1024 (does the constant drift?) | Table 2 |
| `_rules_by_n.py` | four rules × three sizes, with fold-level standard errors | the multi-rule paragraph |
| `_paired.py` | fold-level paired differences between the one-feature and twenty-feature models, with a 95% interval | the paired-comparison paragraph |
| `_phase.py` | reads versus the number of checkpoints k, at n = 128 and 256 | Table 3 |
| `_feat_cost.py` | per-move read cost of six window features (mean, variance, ac₁, mean-vs-half, pattern, palindrome) | the "why no middle case" argument |
| `_kpass.py` | features computable with 1, 2 and 3 passes, all at O(n) reads | the k-pass paragraph |
| `_b0_battery.py` | the eight-rule battery at n=256 | the battery table in earlier drafts |
| `_b0_sec3.py` | writes the battery section into the paper | — |

## Conventions

- The tape is `n` cells, each the average of a pair of bits; each of the two columns has exactly
  `⌈n/2⌉` ones, drawn uniformly.
- BTHM-0 has no memory beyond the header positions, so a repeated state produces no output; runs are
  cut off after `20n` steps and the cap is reported with every number.
- δ = ℓ/(n−1) at the meeting, and every statistic is computed on the meeting set.

## Running

    python _onetable.py 256 4000      # n, number of tapes
    python _by_n.py                   # all four sizes
    python _phase.py 128 12           # n, tapes per row

The scripts need numpy and scikit-learn.
