"""G for peak, on its own best machine (the jump machine at W2 = 128).

peak was the one fully linear rule whose G was missing, because probe_en_linear stored only the
means.  This re-runs its machine per tape and reports G, e, the max, and the ratios.

Usage: python gpeak.py [nt]
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
for d in (os.path.join(PAPER, 'probe_en_linear'), HERE):
    if d not in sys.path:
        sys.path.insert(0, d)

import pf_lib as F                                              # noqa: E402
import el_lib as E                                              # noqa: E402

NS = (512, 1024, 2048, 4096)
NT = int(sys.argv[1]) if len(sys.argv) > 1 else 100


def main():
    print(f"{'n':>6}{'npts':>6}{'exact':>7}{'G/n':>8}{'e/n':>8}{'g/n':>9}{'G/e':>7}{'g/e':>7}")
    for n in NS:
        costs, nex, tot = [], 0, 0
        for i in range(NT * 40):
            if len(costs) >= NT:
                break
            X = E.tape(n, 5000 + n + 100003 * i)
            M = E.Mach(X, n)
            t = M.trace('peak')
            if t['end'] != 0:
                continue
            tot += 1
            f = F.jm_fixed('peak', M.S, n, 128)
            ok = F.exact_of(f, t, n)
            nex += 1 if ok else 0
            costs.append(float(f['reads']))
        if len(costs) < 40:
            print(f'{n:>6}{len(costs):>6}  too few')
            continue
        a = np.array(costs)
        G = float(np.exp(np.mean(np.log(a))))
        e = float(np.mean(a))
        print(f'{n:>6}{a.size:>6}{nex:>4}/{tot:<2}{G / n:>8.3f}{e / n:>8.3f}{a.max() / n:>9.3f}'
              f'{G / e:>7.3f}{a.max() / e:>7.3f}', flush=True)


if __name__ == '__main__':
    main()
