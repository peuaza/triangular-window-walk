"""The tail's weight f(n) and its decay exponent a, from the per-tape reads.

Model: e = L + f*C with L ~ n (the bulk, taken as the median), C ~ n^2 (the tail's mean), and
f ~ n^-a.  Then e ~ n^{max(1, 2-a)}, so a measured exponent lambda in (1,2) implies a = 2 - lambda,
and a stable lambda means a stable a (a genuine intermediate order, not a crossover).

Two tail definitions are reported, because the choice is not canonical:
  tail10  : cost > 10 * median(n)
  tail100n: cost > 100 * n
For each, f, the tail's mean C, the contribution f*C/e, and the exponents of f and C.

Usage: python tail_f.py
"""
import glob
import json
import math
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


def main():
    cells = {}
    for f in glob.glob(os.path.join(HERE, 'k3_*.json')):
        j = json.load(open(f))
        cells[(j['rule'], j['n'])] = np.array(j['reads'], float)
    rules = sorted({r for r, _ in cells})
    print(f"{'rule':<10}{'a (f>10*med)':>14}{'a (f>100n)':>12}{'f*C/e @4096':>13}"
          f"{'pred lam=max(1,2-a)':>21}")
    for rule in rules:
        ns = [n for n in NS if (rule, n) in cells]
        if len(ns) < 3:
            continue
        f10, f100, contrib = [], [], []
        for n in ns:
            a = cells[(rule, n)]
            med = float(np.median(a))
            t10 = a > 10 * med
            t100 = a > 100 * n
            f10.append(float(np.mean(t10)))
            f100.append(float(np.mean(t100)))
            e = float(np.mean(a))
            c = float(np.mean(a[t100])) if t100.any() else 0.0
            contrib.append(float(np.mean(t100)) * c / e if e else float('nan'))
        a10 = exp_of(ns, f10)
        a100 = exp_of(ns, f100)
        pred = max(1.0, 2.0 - a10) if a10 == a10 else float('nan')
        print(f'{rule:<10}{a10:>14.2f}{a100:>12.2f}{contrib[-1]:>13.3f}{pred:>21.2f}')


if __name__ == '__main__':
    main()
