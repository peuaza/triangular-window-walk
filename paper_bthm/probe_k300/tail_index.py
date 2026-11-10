"""The upper tail index of the per-tape cost, by the Hill estimator.

For a Pareto-type upper tail P(cost > t) ~ (t/tau)^-alpha, the Hill estimator on the top k order
statistics is  alpha = 1 / mean_i log(X_(i) / X_(k+1)).

Why it matters: with the maximum at Theta(n^2) and alpha < 1, the mean scales as
e ~ (n^2)^{1-alpha} = n^{2(1-alpha)}, so an intermediate exponent is stable and e << g is generic;
alpha > 1 makes the mean converge to the scale instead.

Reports alpha for k = 5 %, 10 %, 25 % of the sample at each n, and the stability across the sizes.

Usage: python tail_index.py
"""
import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NS = (512, 1024, 2048, 4096)
FRACS = (0.05, 0.10, 0.25)


def hill(a, frac):
    a = np.sort(np.asarray(a, float))[::-1]
    k = max(5, int(frac * a.size))
    k = min(k, a.size - 1)
    top = a[:k]
    ref = a[k]
    return 1.0 / float(np.mean(np.log(top / ref)))


def main():
    cells = {}
    for f in glob.glob(os.path.join(HERE, 'k3_*.json')):
        j = json.load(open(f))
        cells[(j['rule'], j['n'])] = np.array(j['reads'], float)
    rules = sorted({r for r, _ in cells})
    head = ''.join(f"{'a@' + str(int(100 * q)) + '%':>10}" for q in FRACS)
    print(f"{'rule':<10}{'n':>6}{head}{'npts':>7}")
    for rule in rules:
        for n in NS:
            a = cells.get((rule, n))
            if a is None or a.size < 50:
                continue
            row = ''.join(f'{hill(a, q):>10.2f}' for q in FRACS)
            print(f'{rule:<10}{n:>6}{row}{a.size:>7}')
        print()


if __name__ == '__main__':
    main()
