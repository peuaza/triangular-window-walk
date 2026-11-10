"""The geometric mean of the cost: the definitionally-weak cost G <= e.

By AM-GM, exp(E[ln C]) <= E[C] for any positive cost, with no distributional assumption, so

    G_rho <= e_rho <= g_rho

is three cost-type measures on the same meeting tapes, the first definitionally weaker than the
second.  This script reports all three, their exponents over n = 512..4096, and the ratios G/e and
e/g_max, plus the spread of each across rules (the discriminating power).

Usage: python gmean.py
"""
import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NS = (512, 1024, 2048, 4096)


def ex(xs, ys):
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
    print(f"{'rule':<10}{'G exp':>7}{'e exp':>7}{'g exp':>7}{'G/e @4k':>9}{'e/g @4k':>9}"
          f"{'G/n @4k':>9}{'e/n @4k':>9}")
    Gs, es = [], []
    for rule in sorted({r for r, _ in cells}):
        ns = [n for n in NS if (rule, n) in cells]
        if len(ns) < 3:
            continue
        G = [float(np.exp(np.mean(np.log(cells[(rule, n)])))) for n in ns]
        E = [float(np.mean(cells[(rule, n)])) for n in ns]
        Mx = [float(np.max(cells[(rule, n)])) for n in ns]
        n4 = ns[-1]
        print(f'{rule:<10}{ex(ns, G):>7.2f}{ex(ns, E):>7.2f}{ex(ns, Mx):>7.2f}'
              f'{G[-1] / E[-1]:>9.3f}{E[-1] / Mx[-1]:>9.3f}'
              f'{G[-1] / n4:>9.2f}{E[-1] / n4:>9.2f}')
        Gs.append(G[-1] / n4)
        es.append(E[-1] / n4)
    print(f'\nspread across rules at n=4096:  G/n {min(Gs):.1f}..{max(Gs):.1f}'
          f'  ({max(Gs) / min(Gs):.0f}x);  e/n {min(es):.1f}..{max(es):.1f}'
          f'  ({max(es) / min(es):.0f}x)')


if __name__ == '__main__':
    main()
