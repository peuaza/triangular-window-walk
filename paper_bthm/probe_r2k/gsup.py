"""G for the six natural superlinear rules whose per-tape costs were never stored.

Uses the record machine (tr_lib.run) on the meeting tapes, which is the machine the `e` column of
TABLE3.md uses for these rules.  Reports G, e, the max, and G/e -- the quantity that decides
perfection.

Usage: python gsup.py [nt]
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
for d in (os.path.join(PAPER, 'probe_triple_nat'), os.path.join(PAPER, 'probe_e_uniform')):
    if d not in sys.path:
        sys.path.insert(0, d)

import tr_lib as T                                              # noqa: E402
import eu_lib as L                                              # noqa: E402

RULES = ('canonical', 'dense2', 'hashmod', 'pal', 'trip', 'longrun')
NS = (512, 1024, 2048)
NT = int(sys.argv[1]) if len(sys.argv) > 1 else 100


def main():
    print(f"{'rule':<10}{'n':>6}{'npts':>6}{'G/n':>9}{'e/n':>9}{'g/n':>9}{'G/e':>7}{'g/e':>8}")
    for rule in RULES:
        for n in NS:
            out = []
            for i in range(NT * 40):
                if len(out) >= NT:
                    break
                a, s, X = L.make_tape(n, 5000 + n + 100003 * i)
                r = T.run(rule, a, s)
                d = r.get('delta')
                if d is None or d != d:
                    continue
                out.append(float(r['reads']))
            if len(out) < 30:
                print(f'{rule:<10}{n:>6}{len(out):>6}  too few')
                continue
            a = np.array(out)
            G = float(np.exp(np.mean(np.log(a))))
            e = float(np.mean(a))
            print(f'{rule:<10}{n:>6}{a.size:>6}{G / n:>9.2f}{e / n:>9.2f}{a.max() / n:>9.1f}'
                  f'{G / e:>7.3f}{a.max() / e:>8.2f}', flush=True)


if __name__ == '__main__':
    main()
