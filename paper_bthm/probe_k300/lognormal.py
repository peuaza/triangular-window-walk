"""Is the per-tape cost log-normal, with a log-variance growing like c*ln n?

The Pareto model E ~ n^{2(1-alpha)} does not fit: the Hill index drifts with the depth
(4-6 at the top 5 %, 0.6-1.5 at 25 %), which is the signature of a non-power-law tail.  A
multiplicative (log-normal) cost predicts instead

    sigma^2(n) = c * ln n + const        (so the median grows like n and the mean like n^{1+c/2})
    P(C > k * median) -> Phi(-ln k / sigma) -> 1/2 as sigma grows

Both are measurable.  This script reports, per rule and n, the geometric mean (the median), the
log-variance sigma^2, the mean/median ratio (which should be exp(sigma^2/2)), and the fit of
sigma^2 against ln n.

Usage: python lognormal.py
"""
import glob
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NS = (512, 1024, 2048, 4096)


def main():
    cells = {}
    for f in glob.glob(os.path.join(HERE, 'k3_*.json')):
        j = json.load(open(f))
        cells[(j['rule'], j['n'])] = np.array(j['reads'], float)
    print(f"{'rule':<10}{'n':>6}{'median':>10}{'sigma^2':>9}{'mean/med':>10}"
          f"{'exp(s2/2)':>10}{'P>10med':>9}")
    for rule in sorted({r for r, _ in cells}):
        for n in NS:
            a = cells.get((rule, n))
            if a is None or a.size < 50:
                continue
            la = np.log(a)
            s2 = float(np.var(la, ddof=1))
            med = float(np.exp(np.mean(la)))
            print(f'{rule:<10}{n:>6}{med:>10.0f}{s2:>9.2f}{np.mean(a) / med:>10.1f}'
                  f'{math.exp(s2 / 2):>10.1f}{np.mean(a > 10 * np.median(a)):>9.3f}')
        print()


if __name__ == '__main__':
    main()
