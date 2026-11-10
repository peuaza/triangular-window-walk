"""The reliability curve p(alpha): the budget needed to be right on a 1-alpha fraction of tapes.

A predictor that finishes the tapes inside a budget B and declines the rest is right on exactly the
finished fraction, so R^2 = 1 - alpha with B = q_{1-alpha}(cost).  Hence

    p(alpha) = q_{1-alpha}(cost),

the cost quantile.  p(0) is the maximum (the worst tape), and the mean e sits where q = e.

The curve is threshold-free (alpha is the free variable, not a cut), which is why it is the
principled weaker notion than p: it reports the whole trade-off instead of one point.

Usage: python pcurve.py [n]
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
N = int(sys.argv[1]) if len(sys.argv) > 1 else 4096
ALPHAS = (0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50)


def main():
    cells = {}
    for f in glob.glob(os.path.join(HERE, f'k3_*_{N}.json')):
        j = json.load(open(f))
        cells[j['rule']] = np.array(j['reads'], float)
    head = ''.join(f"{'a=' + str(a):>8}" for a in ALPHAS)
    print(f"{'rule':<10}{head}{'e/n':>9}{'ratio':>7}")
    for rule in sorted(cells):
        a = cells[rule]
        e = float(np.mean(a))
        qs = [float(np.quantile(a, 1 - al)) / N for al in ALPHAS]
        # the alpha at which the curve crosses e
        cross = float('nan')
        for al in np.linspace(0.0, 0.95, 96):
            if np.quantile(a, 1 - al) <= e:
                cross = al
                break
        row = ''.join(f'{q:>8.1f}' for q in qs)
        print(f'{rule:<10}{row}{e / N:>9.1f}{cross:>7.2f}')


if __name__ == '__main__':
    main()
