# -*- coding: utf-8 -*-
"""Sanity baseline: reproduce the published single-statistic R^2 (ac1, cubic fit)
so that the tape generator and the meeting-set convention in this workspace are
known to match the manuscript's (0.683 for the cubic fit at n=256).

Also reports R^2 for ac1 with gradient boosting, for reference.

Usage: python _bz_base.py <n> <tapes_wanted> [seed]
"""
import sys

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import KFold

from _bz_inc import make_tape, Tap, r2
from _bz_final import run_one


def ac1(a):
    u, v = a[:-1], a[1:]
    if u.std() < 1e-12 or v.std() < 1e-12:
        return 0.0
    return float(np.corrcoef(u, v)[0, 1])


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 256
    want = int(sys.argv[2]) if len(sys.argv) > 2 else 1200
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 2026
    rng = np.random.default_rng(seed)
    X, A = [], []
    tried = 0
    while len(X) < want and tried < 60 * want:
        tried += 1
        a = make_tape(rng, n)
        d = run_one(Tap(np.round(a * 2).astype(np.int32)), n)[0]
        if d is None:
            continue
        X.append(d)
        A.append(ac1(a))
    X = np.array(X)
    A = np.array(A)
    print('n=%d  meeting tapes=%d  (generated %d)  sd(X)=%.4f' % (n, len(X), tried, X.std(ddof=1)))
    kf = KFold(5, shuffle=True, random_state=0)
    folds = list(kf.split(A))
    for deg in (1, 3):
        sc = []
        for tr, te in folds:
            c = np.polyfit(A[tr], X[tr], deg)
            p = np.polyval(c, A[te])
            sc.append(1 - ((X[te] - p) ** 2).sum() / ((X[te] - X[tr].mean()) ** 2).sum())
        print('  ac1, polynomial degree %d : CV R2 = %.4f  (folds %s)'
              % (deg, np.mean(sc), ' '.join('%.3f' % s for s in sc)))
    sc = []
    for tr, te in folds:
        m = HistGradientBoostingRegressor(max_iter=250, learning_rate=.06, max_depth=5, random_state=0)
        m.fit(A[tr].reshape(-1, 1), X[tr])
        p = m.predict(A[te].reshape(-1, 1))
        sc.append(1 - ((X[te] - p) ** 2).sum() / ((X[te] - X[tr].mean()) ** 2).sum())
    print('  ac1, gradient boosting  : CV R2 = %.4f  (folds %s)'
          % (np.mean(sc), ' '.join('%.3f' % s for s in sc)))


if __name__ == '__main__':
    main()
