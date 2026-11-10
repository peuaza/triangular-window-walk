"""The two asymptotic gates, with their current numerical positions.

Semi-linear gate.  If the mean is tail-dominated with the maximum at Theta(n^2), then
e ~ max(n, n^{2(1-alpha)}); the class is filled as soon as alpha >= 1/2 (the floor at n takes over).

Quasi-linear gate.  If the tail's *fraction* f(n) -> 0, a predictor that finishes the typical tapes
and gives up on the rest reaches R^2 = 1 - f -> 1 at a cost p ~ (1-f)*bulk, i.e. p = O(n) while
e = omega(n).

Both gates are reported per rule, with the measured f (share above 10x the median), its exponent a,
the implied R^2 = 1 - f at n = 4096, and the Hill index at the 25 % depth (the depth nearest the
mass that dominates the mean).

Usage: python gates.py
"""
import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NS = (512, 1024, 2048, 4096)


def exp_of(xs, ys):
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    ok = ys > 0
    if ok.sum() < 2:
        return float('nan')
    return float(np.polyfit(np.log(xs[ok]), np.log(ys[ok]), 1)[0])


def hill(a, frac=0.25):
    a = np.sort(np.asarray(a, float))[::-1]
    k = max(5, min(int(frac * a.size), a.size - 1))
    return 1.0 / float(np.mean(np.log(a[:k] / a[k])))


def main():
    cells = {}
    for f in glob.glob(os.path.join(HERE, 'k3_*.json')):
        j = json.load(open(f))
        cells[(j['rule'], j['n'])] = np.array(j['reads'], float)
    print(f"{'rule':<10}{'f@512':>7}{'f@4096':>8}{'a(f)':>6}{'R^2=1-f':>9}"
          f"{'a25':>6}{'e-implied':>10}{'gate':>8}")
    for rule in sorted({r for r, _ in cells}):
        ns = [n for n in NS if (rule, n) in cells]
        if len(ns) < 3:
            continue
        f = []
        for n in ns:
            a = cells[(rule, n)]
            f.append(float(np.mean(a > 10 * float(np.median(a)))))
        af = exp_of(ns, f)
        a25 = hill(cells[(rule, ns[-1])], 0.25)
        impl = max(1.0, 2.0 * (1.0 - a25)) if a25 == a25 else float('nan')
        gate = 'SL if a>=.5' if a25 >= 0.5 else 'no SL'
        print(f'{rule:<10}{f[0]:>7.3f}{f[-1]:>8.3f}{af:>6.2f}{1 - f[-1]:>9.3f}'
              f'{a25:>6.2f}{impl:>10.2f}{gate:>8}')


if __name__ == '__main__':
    main()
