"""The failure-set concentration: is there a budget at which eps*M << e?

For a budget B (a multiple of n), let eps(B) be the share of meeting tapes whose cost exceeds B and
M(B) their mean cost.  Patching the budget predictor with a fallback gives e <= B + eps*M, so

    p >= e - eps*M.

That bound is useful only where eps*M is small compared with e.  The question this script answers:
does any budget have eps(B)*M(B) << e?  If none does, the "failures are cheap" route to p << e is
closed at every linear-to-quadratic budget; if one does, that budget is where to look for a predictor.

Reports, per rule at n = 4096: eps, M/n, eps*M/e, and eps*M/B over a grid of B/n.

Usage: python concentration.py
"""
import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BGRID = (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000)
N = 4096


def main():
    cells = {}
    for f in glob.glob(os.path.join(HERE, f'k3_*_{N}.json')):
        j = json.load(open(f))
        cells[j['rule']] = np.array(j['reads'], float)
    print(f"{'rule':<10}{'B/n':>6}{'eps':>8}{'M/n':>9}{'epsM/e':>9}{'epsM/B':>9}")
    for rule in sorted(cells):
        a = cells[rule]
        e = float(np.mean(a))
        for b in BGRID:
            B = b * N
            sel = a > B
            if not sel.any():
                continue
            eps = float(np.mean(sel))
            M = float(np.mean(a[sel]))
            print(f'{rule:<10}{b:>6}{eps:>8.3f}{M / N:>9.1f}{eps * M / e:>9.3f}'
                  f'{eps * M / B:>9.3f}')
        print()


if __name__ == '__main__':
    main()
