# Runtime and hardware

All numbers in the paper were produced on the machine below, with the scripts in this folder,
using `run_all.py`. No GPU is used; everything is single-threaded NumPy/SciPy.

## Machine

```
CPU      AMD Ryzen AI 7 H 350 w/ Radeon 860M
RAM      31.1 GB
OS       Windows 11
Python   3.12.14
NumPy    2.5.2
SciPy    1.18.1
```

## Measured wall-clock (single run, this machine)

| Script | Paper location | What it does | Time |
|---|---|---|---|
| `wm_exact.py` | Table 1 | exhaustive enumeration, n=4..11 | 65 s total (0,0,0,0,1,3,13,48 s per n) |
| `wm12.py` | Table 1 | exhaustive enumeration, n=12 (853,776 matrices) | ≈ 6 min |
| `fano.py` | Sec. 4, Link 3 | Fano factors of reversal intervals | 73 s |
| `confirm_c.py` | Table 3, convention A | occurrence rate, 2,000 δ at n=2048 | 405 s |
| `walk_bin.py` (single walk) | — | one walk at n=8192 | ≈ 25 s |

The remaining scripts were not timed individually; `run_all.py` writes a complete
`_results/timings.csv` when run end to end, and the full run takes roughly 1–2 hours on this
machine.

## How to reproduce

```bash
python run_all.py --quick     # seconds-to-minutes subset, checks the environment
python run_all.py            # everything; writes _results/<script>.out and _results/timings.csv
python run_all.py --only wm12.py
```

Each script prints its results to stdout and is independent of the others, apart from importing
`walk_bin.py`.

## Determinism

All Monte Carlo scripts use `numpy.random.default_rng(<fixed seed>)`, so a rerun on the same
NumPy version reproduces the reported numbers exactly. Exact enumeration is deterministic by
construction. Where a script reports a rate over several seeds, the seed is printed in its output.
