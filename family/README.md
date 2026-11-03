# Code for paper 2: *Memory and randomness control the shape*

Each script reproduces one table or figure of the paper.

| script | produces |
|---|---|
| `_four_components.py` | the exact enumeration at n ≤ 13, split by first cut (no cut / column / row / both); the four-component table |
| `_counts_by_first.py` | the same counts by first event for n ≤ 12, per number of ones k |
| `_centre_components.py`, `_centre_components2.py` | the centre family split into its two components, with the KS tests against Beta and Normal |
| `_centre_mixture.py`, `_centre_mixture2.py` | the two-component Beta mixture fit (six parameters) — the fit degenerates, see the paper |
| `_tau_scaling.py` | the decay of the displacement with n, with standard errors and the log–log slope |
| `_phase_grid.py`, `_phase_grid_hist.py` | the phase diagram: 6 strategies × 3 jump rules, with the density curves |
| `_palindrome_knob.py` | the palindrome-distance knob |
| `_counter2048.py`, `_counter512big.py` | the counter rule at n = 2048 and a 512-tape re-run |
| `_ac_vs_firstcut.py` | the autocorrelation against the first-cut class |

## Conventions

- The tape is `n` cells, each the average of a pair of bits; each of the two columns has exactly
  `⌈n/2⌉` ones.
- The walk follows the four-way table of paper 1.
- Every script uses `numpy.random.default_rng` with a fixed seed and prints to stdout.

## Running

    python _four_components.py
    python _tau_scaling.py

NumPy is required; `_centre_mixture*.py` also uses SciPy.
