"""R^2(k), v2: compact summary features, and prefix rather than random sampling.

v1 asked a random forest to learn delta from the raw values and positions of k random cells with 40
tapes; every R^2 came out negative, which is as consistent with overfitting as with 'no information'.
This version removes both confounds:

  features  six summaries only (density in three thirds, longest run, number of runs, mean value)
  sampling  two modes: 'rand' (a uniform k-subset) and 'prefix' (the first k cells)
  tapes     150 per rule

If R^2(k) stays near zero with six features and 150 tapes, the partial-information route really is
closed; if it rises, the hard tapes do carry learnable structure and p << e becomes possible.

Usage: python r2k2.py [nt] [n]
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

RULES = ('AMAX', 'BAL>m/2', 'KCHG', 'LAS', 'MAXP', 'NBLK=1', 'dense2', 'oddrun', 'trip',
         'pattern', 'runs')
KS = (16, 64, 256)
NT = int(sys.argv[1]) if len(sys.argv) > 1 else 150
N = int(sys.argv[2]) if len(sys.argv) > 2 else 512


def summary(v, pos, n):
    """Six numbers from a sampled subsequence."""
    v = np.asarray(v, float)
    pos = np.asarray(pos, float)
    thirds = []
    for lo, hi in ((0, 1 / 3), (1 / 3, 2 / 3), (2 / 3, 1.01)):
        m = (pos >= lo) & (pos < hi)
        thirds.append(v[m].mean() if m.any() else np.nan)
    d = np.diff(v)
    runs = 1 + int(np.sum(d != 0)) if v.size else 0
    # longest run of equal values
    best = cur = 1
    for x in d:
        cur = cur + 1 if x == 0 else 1
        best = max(best, cur)
    return np.array(thirds + [best / max(1, v.size), runs / max(1, v.size), v.mean()])


def main():
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import cross_val_predict, KFold

    for mode in ('rand', 'prefix'):
        print(f'=== mode {mode}')
        print(f"{'rule':<10}{'hard':>6}" + ''.join(f"{'k=' + str(k):>9}" for k in KS) + f"{'e/n':>8}")
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
            if len(pairs) < 40:
                print(f'{rule:<10}  too few ({len(pairs)})')
                continue
            costs = np.array(costs)
            deltas = np.array(deltas)
            hard = float(np.mean(costs > 10 * np.median(costs)))
            kap = KFold(5, shuffle=True, random_state=0)
            row = []
            for k in KS:
                F = np.zeros((len(pairs), 6))
                for j, (a, s) in enumerate(pairs):
                    if mode == 'prefix':
                        idx = np.arange(k)
                    else:
                        rng = np.random.default_rng(1000 * k + j)
                        idx = np.sort(rng.choice(a.size, size=k, replace=False))
                    F[j] = summary(a[idx], idx / max(1, N - 1), N)
                F = np.nan_to_num(F)
                rf = RandomForestRegressor(n_estimators=200, min_samples_leaf=3,
                                           random_state=0, n_jobs=1)
                pr = cross_val_predict(rf, F, deltas, cv=kap, n_jobs=1)
                ss = np.var(deltas) * len(deltas)
                row.append(1.0 - np.sum((pr - deltas) ** 2) / ss)
            print(f'{rule:<10}{hard:>6.2f}' + ''.join(f'{v:>9.2f}' for v in row)
                  + f'{np.mean(costs) / N:>8.1f}', flush=True)


if __name__ == '__main__':
    main()
