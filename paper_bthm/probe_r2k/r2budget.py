"""1/(1-R^2) as a family: the budget needed for 1-R^2 = n^{-theta}.

p is the budget at which 1-R^2 = 0 (theta = infinity).  This measures the same protocol at finite
theta: for a grid of budgets B (as multiples of n) the predictor follows the walk for B reads, then
predicts delta from its state (kp_lib's protocol, out-of-fold R^2, ten folds), and we report the
curve 1-R^2(B).  From it, q(theta) = min{B : 1-R^2(B) <= n^{-theta}}.

The point: q(theta) <= e for every theta (an exact machine attains 1-R^2 = 0), so this family lives
between the cheapest useful budget and e, and the question is how fast it falls as theta grows.

Usage: python r2budget.py [n] [nt]
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

RULES = ('pattern', 'runs', 'pal', 'longrun')
N = int(sys.argv[1]) if len(sys.argv) > 1 else 1024
NT = int(sys.argv[2]) if len(sys.argv) > 2 else 150
BS = (2, 5, 10, 20, 50, 100, 200)


def main():
    print(f"{'rule':<10}{'e/n':>7}" + ''.join(f"{'B=' + str(b) + 'n':>10}" for b in BS))
    for rule in RULES:
        tapes, ys, es = [], [], []
        for i in range(NT * 30):
            if len(tapes) >= NT:
                break
            a, s, X = L.make_tape(N, 5000 + N + 100003 * i)
            r = T.run(rule, a, s)
            d = r.get('delta')
            if d is None or d != d:
                continue
            tapes.append(np.asarray(a, np.int64))
            ys.append(float(d))
            es.append(float(r['reads']))
        if len(tapes) < 40:
            print(f'{rule:<10}  too few ({len(tapes)})')
            continue
        ys = np.array(ys)
        row = []
        for b in BS:
            gf = [K.global_feats(a) for a in tapes]
            Xf = np.array([K.state_feats_at(a, rule, b * N, g)[0] for a, g in zip(tapes, gf)])
            pr = K.cv_pred(lambda Xtr, ytr: (lambda Xte: K.ridge_fit_pred(Xtr, ytr, Xte)),
                           Xf, ys, folds=10, seed=0)
            row.append(1.0 - K.r2(ys, pr))
        print(f'{rule:<10}{np.mean(es) / N:>7.1f}' + ''.join(f'{v:>10.3f}' for v in row),
              flush=True)


if __name__ == '__main__':
    main()
