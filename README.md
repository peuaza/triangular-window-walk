# Code for "A triangular window walk on random binary matrices"

This folder reproduces every number reported in the paper.

## Requirements

- Python 3.10+
- NumPy, SciPy (`pip install numpy scipy`)

All scripts are single-file, run from this folder, and print their results to stdout.
Monte Carlo runs use `numpy.random.default_rng` with fixed seeds, so results are reproducible.

## The walk

`walk_bin.py` implements the map of Section 2 on an `n x 2` binary matrix.
It is imported by all other scripts. No re-binarisation is performed.

## Scripts by table

| Paper | Script | Output |
|---|---|---|
| Table 1 (`n <= 12`, exact) | `wm_exact.py` (n=4..11), `wm12.py` (n=12) | matrices / delta count / rate |
| Table 2 (distribution) | `arcsine_conv.py` | sd, median, KS vs arcsine, n=48..512 |
| Table 3 (rate, convention A) | `confirm_c.py`, `quarter.py` | rate at n=256..2048 |
| Table 3 (convention B) | `fallback_big.py` | rate at n=256..2048 |
| Table 4 (conditional spread) | `p_cond.py` | conditional mean/sd of p by ac bin |
| Table 5 (reversal rate) | `ac_corr3.py`, `switch_m.py` | rho(m), s(m) |
| Link 1 | `ac_path.py` | autocorrelation of local ac along the path |
| Link 2 | `cluster_src.py` | reversal rate vs abs(ac) |
| Link 3 | `fano.py` | Fano factor of reversal intervals |
| Link 4 | `maxrun_by_k.py` | longest-run share by number of reversals |
| Link 5 | `burst_pos.py` | burst position vs fraction |

## Notes

- `wm12.py` enumerates 853,776 matrices and takes about six minutes.
- `fallback_big.py` is the slowest Monte Carlo script (about 30 minutes).
- Some scripts emit `RuntimeWarning: invalid value encountered in divide` from
  `numpy.corrcoef` on constant windows; this is expected and handled by the
  convention described in Remark 1 of the paper (undefined correlation -> b1 = 0).
- Enumerations for `n <= 12` are exhaustive; everything else is Monte Carlo with
  the sample sizes stated in the paper.
