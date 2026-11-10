# Code for the triangular-window-walk papers

**Paper 1 preprint:** concept DOI 10.5281/zenodo.23804105

This repository holds the code for three papers that share the same tape and the same walk. The
scripts at the root belong to paper 1; the two subdirectories belong to the papers that grew out of
it, and each has its own README.

| folder | paper |
|---|---|
| `.` (root) | *A triangular window walk on random binary matrices: a conjectured arcsine limit* |
| `family/` | *Memory and randomness control the shape: a family of limit laws from a triangular-window walk* |
| `paper_bthm/` | *A memoryless two-header tape machine: what a window costs, and why two lower-bound techniques fail* (the current paper; the arXiv/Zenodo version is `10.5281/zenodo.24115991`) |

Nothing at the root was moved when the subdirectories were added, so every path in the paper-1
README below still resolves. The root's `archive/` holds the superseded drafts (`bthm_v2.tex`,
`bthm_jalc.tex`) and the logs they produced; it is not part of any paper; it holds the earlier BTHM draft (`bthm0/`) as well as the superseded `bthm_v2.tex` and `bthm_jalc.tex`.

## Paper 1

This folder reproduces every number reported in paper 1.

## Reproducing everything

```bash
python run_all.py --quick   # fast subset first, to check the environment
python run_all.py           # full run, 1-2 hours; writes _results/timings.csv
```

Hardware, software versions and measured wall-clock times are in `RUNTIME.md`.

## Requirements

- Python 3.10+
- NumPy, SciPy (`pip install numpy scipy`)

All scripts are single-file, run from this folder, and print their results to stdout.
Monte Carlo runs use `numpy.random.default_rng` with fixed seeds, so results are reproducible.

## The walk

`walk_bin.py` implements the map of Section 2 on an `n x 2` binary matrix.
It is imported by all other scripts. No re-binarisation is performed.

## Scripts by item

Paper 1 has two tables and six figures. This is the mapping, in the order they appear.

| item in paper 1 | script | what it prints |
|---|---|---|
| Table 1: exhaustive check of the reduction, $n\le 12$ | `wm_exact.py` (n=4..11), `wm12.py` (n=12) | matrices, classes, inconsistent classes |
| Table 2: specification sensitivity, n=256 | `_spec_sens.py` | rate, sd, median and KS for the four variants |
| Fig. 2, 4, 5: histogram, trend and Q-Q against the arcsine | `arcsine_conv.py` | sd, median, KS, n=48..512 |
| Fig. 6: reversal rate against $rccos
ho(m)/\pi$ | `ac_corr3.py`, `switch_m.py` | $
ho(m)$, $s(m)$ |
| Fig. 1, 3: flowchart and one walk | — | schematic, no data |
| rate at n=256..2048, convention A | `confirm_c.py`, `quarter.py` | rate |
| rate, convention B (revisit terminates) | `fallback_big.py` | rate |
| conditional spread of the output | `p_cond.py` | conditional mean and sd by ac bin |
| Link 1 | `ac_path.py` | autocorrelation of the local ac along the path |
| Link 2 | `cluster_src.py` | reversal rate against absolute ac |
| Link 3 | `excursion.py`, `fano.py` | excursion lengths; Fano factor of reversal intervals |
| Link 4 | `by_switch.py`, `maxrun_by_k.py` | longest-run share by number of reversals |
| Link 5 | `burst_pos.py` | burst position against fraction |
| Section 2.6 | `segments.py`, `qseg.py` | decomposition into segments; segment statistics |
| — | `run_all.py` | driver: runs everything, records wall-clock time |

> `enum_small.py` is a legacy implementation of Table 1, superseded by `wm_exact.py`.

## Notes

- `wm12.py` enumerates 853,776 matrices and takes about six minutes.
- `fallback_big.py` is the slowest Monte Carlo script (about 30 minutes).
- Some scripts emit `RuntimeWarning: invalid value encountered in divide` from
  `numpy.corrcoef` on constant windows; this is expected and handled by the
  convention described in Remark 1 of the paper (undefined correlation -> b1 = 0).
- Enumerations for `n <= 12` are exhaustive; everything else is Monte Carlo with
  the sample sizes stated in the paper.
