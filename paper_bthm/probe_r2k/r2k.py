"""R^2(k): how much of delta is learnable from k random cells?

The statistical version of "can the expensive tapes be predicted cheaply".  For each tape we draw a
uniform random k-subset of cells, hand the predictor their values and positions, and ask how much of
delta's variance the prediction explains (out-of-fold R^2).  Reported on all meeting tapes and on
the hard stratum (cost > 10x the median), for k = 8 ... n/2.

If R^2(k) rises to near 1 at k << e on the hard stratum, then a learned predictor gives p << e.  If it
stays flat, the hard tapes' delta carries no partial-information structure, and p = Theta(e) has an
information-theoretic reason.

Usage: python r2k.py [nt] [n]
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

RULES = ('AMAX', 'AMIN', 'BAL>m/2', 'BAL>m/4', 'KCHG', 'LAS', 'MAXP', 'NBLK=1', 'NBLK>2',
         'NBLK>3', 'dense2', 'oddrun', 'trip', 'pattern', 'runs')
KS = (8, 16, 32, 64, 128, 256)
NT = int(sys.argv[1]) if len(sys.argv) > 1 else 150
N = int(sys.argv[2]) if len(sys.argv) > 2 else 512


def collect(rule, nt, n):
    d, c = [], []
    for i in range(nt * 30):
        if len(d) >= nt:
            break
        a, s, X = L.make_tape(n, 5000 + n + 100003 * i)
        r = T.run(rule, a, s)
        dd = r.get('delta')
        if dd is None or dd != dd:
            continue
        d.append(float(dd))
        c.append(float(r['reads']))
    return np.array(d), np.array(c), [L.make_tape(n, 5000 + n + 100003 * i)[0]
                                      for i in range(0)] or None


def feats(a, idx):
    v = a[idx].astype(float)
    p = idx.astype(float) / max(1, a.size - 1)
    return np.concatenate([v, p])


def main():
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import cross_val_predict, KFold

    print(f"{'rule':<10}{'thr':>6}" + ''.join(f"{'k=' + str(k):>9}" for k in KS) + f"{'e/n':>8}")
    for rule in RULES:
        pairs, deltas, costs = [], [], []
        for i in range(NT * 30):
            if len(pairs) >= NT:
                break
            a, s, X = L.make_tape(N, 5000 + N + 100003 * i)
            r = T.run(rule, a, s)
            dd = r.get('delta')
            if dd is None or dd != dd:
                continue
            pairs.append((np.asarray(a, np.int64), np.asarray(s, np.int64)))
            deltas.append(float(dd))
            costs.append(float(r['reads']))
        if len(pairs) < 30:
            print(f'{rule:<10}  too few tapes ({len(pairs)})')
            continue
        tapes = [p[0] for p in pairs]
        costs = np.array(costs)
        deltas = np.array(deltas)
        med = float(np.median(costs))
        hard = costs > 10 * med
        kap = KFold(5, shuffle=True, random_state=0)
        row = []
        for k in KS:
            Xf = np.zeros((len(tapes), 2 * k))
            for j, a in enumerate(tapes):
                rng = np.random.default_rng(1000 * k + j)
                idx = np.sort(rng.choice(a.size, size=k, replace=False))
                Xf[j] = feats(a, idx)
            rf = RandomForestRegressor(n_estimators=120, min_samples_leaf=2,
                                       random_state=0, n_jobs=1)
            pr = cross_val_predict(rf, Xf, deltas, cv=kap, n_jobs=1)
            ss = np.var(deltas) * len(deltas)
            row.append(1.0 - np.sum((pr - deltas) ** 2) / ss)
        print(f'{rule:<10}{np.mean(hard):>6.2f}' + ''.join(f'{v:>9.2f}' for v in row)
              + f'{np.mean(costs) / N:>8.1f}', flush=True)


if __name__ == '__main__':
    main()
