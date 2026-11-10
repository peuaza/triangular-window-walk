"""G for the fully linear rules, measured on their best machine.

pattern and runs: the CES device (kp_lib.ces), which is exact on the meeting tapes.
monotone: the collapse machine (tr_lib), cost 2n-1 by construction.
peak: needs the jump machine of probe_en_linear, whose per-tape costs were not stored; reported
separately when that is re-run.

G = exp(mean log cost), e = mean cost, g = max cost; AM-GM gives G <= e.

Usage: python gfl.py [nt]
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
for d in (os.path.join(PAPER, 'probe_kappa_pred'), os.path.join(PAPER, 'probe_triple_nat'),
          os.path.join(PAPER, 'probe_e_uniform')):
    if d not in sys.path:
        sys.path.insert(0, d)

import kp_lib as K                                              # noqa: E402
import tr_lib as T                                              # noqa: E402
import eu_lib as L                                              # noqa: E402

NS = (512, 1024, 2048, 4096)
NT = int(sys.argv[1]) if len(sys.argv) > 1 else 150


def ex(xs, ys):
    xs = np.asarray(xs, float)
    ys = np.asarray(ys, float)
    ok = ys > 0
    return float(np.polyfit(np.log(xs[ok]), np.log(ys[ok]), 1)[0]) if ok.sum() > 1 else float('nan')


def read_costs(rule, n, nt):
    out = []
    for i in range(nt * 40):
        if len(out) >= nt:
            break
        a, s, X = L.make_tape(n, 5000 + n + 100003 * i)
        if rule == 'monotone':
            r = T.run('monotone', a, s)
            if r.get('delta') is None or r['delta'] != r['delta']:
                continue
            out.append(float(r['reads']))
        else:
            d, reads, info = K.ces(a, rule=rule, W=16)
            if d is None or d != d:
                continue
            out.append(float(reads))
    return np.array(out)


def main():
    print(f"{'rule':<10}{'n':>6}{'npts':>6}{'G/n':>8}{'e/n':>8}{'g/n':>9}{'G/e':>7}{'g/e':>7}")
    for rule in ('monotone', 'pattern', 'runs'):
        for n in NS:
            a = read_costs(rule, n, NT)
            if a.size < 40:
                print(f'{rule:<10}{n:>6}{a.size:>6}  too few')
                continue
            G = float(np.exp(np.mean(np.log(a))))
            e = float(np.mean(a))
            print(f'{rule:<10}{n:>6}{a.size:>6}{G / n:>8.3f}{e / n:>8.3f}{a.max() / n:>9.3f}'
                  f'{G / e:>7.3f}{a.max() / e:>7.3f}', flush=True)


if __name__ == '__main__':
    main()
