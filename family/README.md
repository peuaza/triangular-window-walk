# Code for paper 2: *Memory and randomness control the shape*

Each script reproduces one table or figure of the paper.

| Table in paper 2 | Script | What it computes |
|---|---|---|
| Table 1: memory levels M0–M3 at n=256 | `_memory_fig.py` | rate, sd and kurtosis of δ for the four levels |
| Table 2: six modifications at n=256 | `_new_actions.py` | jump of length two and three, slide, local rewriting, and the b2 variants |
| Table 2, the b2 rows at N=2000 | `_accumulate_mod.py` | large-sample runs of `self`, `b2last`, `b2half` |
| Table 3: three families at three sizes | `_dist_analysis.py` | reads `_dist/<n>.jsonl`; prints N, sd, median, kurtosis and KS |
| Tables 4 and 5: class and matrix counts | `_counts13.py` | exact enumeration of sequences and matrices with δ=k, n=4..12 |
| Tables 4 and 5, second implementation | `_counts_class.py` | the same counts, written independently |
| Four-component split by first cut | `_four_components.py` | exact enumeration at n≤13, split by first cut |
| Counts by first event | `_counts_by_first.py` | the same split for n≤12, per number of ones k |
| The centre family, fitted | `_centre_law.py`, `_centre_components.py`, `_centre_mixture.py` | candidate limit laws; Beta/Normal fits; the two-component mixture (degenerates, see the paper) |
| Phase diagram | `_phase_grid.py`, `_phase_grid_hist.py` | 6 strategies × 3 jump rules, with the density curves |
| Reversal rate against the first-cut class | `_ac_vs_firstcut.py` | |
| The palindrome knob | `_palindrome_knob.py` | |
| Counter rule at n=2048, and a 512-tape re-run | `_counter2048.py`, `_counter512big.py` | |
| Decay of the displacement with n | `_tau_scaling.py` | |

## Conventions

- The tape is `n` cells, each the average of a pair of bits; each of the two columns has exactly
  `⌈n/2⌉` ones.
- The walk follows the four-way table of paper 1.
- Every script uses `numpy.random.default_rng` with a fixed seed and prints to stdout.

## Running

    python _four_components.py
    python _tau_scaling.py

NumPy is required; `_centre_mixture*.py` also uses SciPy.
