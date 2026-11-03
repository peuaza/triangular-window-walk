# -*- coding: utf-8 -*-
"""Can cheap single-pass features reach R^2 = 0.85 on the canonical rule?

Feature family: ~24 statistics, each computable in one streaming pass over the
tape with O(log n) counters (no window rescans).  Gradient boosting, 5-fold CV,
on the meeting set.  This is the "statistical predictor" reading of the ceiling
conjecture, tested at the same n as the published 0.766 measurement.

Usage: python _bz_feat.py <n> <meeting_tapes> [seed]
"""
import sys
import time

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import KFold

from _bz_inc import make_tape, Tap, r2
from _bz_final import run_one


def corr(u, v):
    if u.std() < 1e-12 or v.std() < 1e-12:
        return 0.0
    return float(np.corrcoef(u, v)[0, 1])


def feats(a):
    n = len(a)
    f = {}
    f['ac1'] = corr(a[:-1], a[1:])
    f['var'] = float(a.var())
    f['runs'] = float((np.abs(np.diff(a)) > 1e-9).sum()) / n
    idx = np.arange(n, dtype=float)
    f['centroid'] = float((a * idx).sum() / (a.sum() + 1e-12)) / n
    ones = a > 0.9
    zeros = a < 0.1
    f['first1'] = float(np.argmax(ones)) / n if ones.any() else 1.0
    f['last0'] = float(n - 1 - np.argmax(zeros[::-1])) / n if zeros.any() else 1.0
    # longest runs
    best1 = best0 = cur1 = cur0 = 0
    pos1 = pos0 = 0
    for i, v in enumerate(a):
        if v > 0.9:
            cur1 += 1; cur0 = 0
            if cur1 > best1:
                best1, pos1 = cur1, i
        elif v < 0.1:
            cur0 += 1; cur1 = 0
            if cur0 > best0:
                best0, pos0 = cur0, i
        else:
            cur1 = cur0 = 0
    f['run1'] = best1 / n
    f['run0'] = best0 / n
    f['run1pos'] = pos1 / n
    f['run0pos'] = pos0 / n
    # block means at four scales
    for sc in (4, 8, 16, 32):
        bl = a[:(n // sc) * sc].reshape(-1, sc).mean(1)
        f['b%d_mean' % sc] = float(bl.mean())
        f['b%d_sd' % sc] = float(bl.std())
        f['b%d_argmax' % sc] = float(np.argmax(bl)) / len(bl)
        f['b%d_argmin' % sc] = float(np.argmin(bl)) / len(bl)
    # block ac1 at scale 32
    sc = 32
    nblk = n // sc
    acs = [corr(a[i * sc:(i + 1) * sc - 1], a[i * sc + 1:(i + 1) * sc]) for i in range(nblk)]
    acs = np.array(acs)
    f['ac_blk_mean'] = float(acs.mean())
    f['ac_blk_sd'] = float(acs.std())
    f['ac_blk_min'] = float(acs.min())
    f['ac_blk_max'] = float(acs.max())
    f['ac_blk_argmax'] = float(np.argmax(acs)) / nblk
    # prefix imbalance of the ternary walk (one pass)
    d = a - 0.5
    cs = np.cumsum(d)
    f['walk_max'] = float(cs.max()) / n
    f['walk_min'] = float(cs.min()) / n
    f['walk_argmax'] = float(np.argmax(cs)) / n
    f['walk_argmin'] = float(np.argmin(cs)) / n
    return f


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 256
    want = int(sys.argv[2]) if len(sys.argv) > 2 else 800
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 2026
    rng = np.random.default_rng(seed)
    X, F = [], []
    tried = 0
    t0 = time.time()
    while len(X) < want and tried < 60 * want:
        tried += 1
        a = make_tape(rng, n)
        d = run_one(Tap(np.round(a * 2).astype(np.int32)), n)[0]
        if d is None:
            continue
        X.append(d)
        F.append(feats(a))
    X = np.array(X)
    keys = list(F[0].keys())
    M = np.array([[r[k] for k in keys] for r in F])
    print('n=%d meeting tapes=%d (generated %d) features=%d  sd(X)=%.4f  gen %.0fs'
          % (n, len(X), tried, len(keys), X.std(ddof=1), time.time() - t0))
    folds = list(KFold(5, shuffle=True, random_state=0).split(M))

    def cv(cols):
        sc = []
        for tr, te in folds:
            m = HistGradientBoostingRegressor(max_iter=250, learning_rate=.06, max_depth=5, random_state=0)
            m.fit(M[np.ix_(tr, cols)], X[tr])
            p = m.predict(M[np.ix_(te, cols)])
            sc.append(1 - ((X[te] - p) ** 2).sum() / ((X[te] - X[tr].mean()) ** 2).sum())
        return np.array(sc)

    sc = cv(list(range(len(keys))))
    print('all %d cheap features, GBM: CV R2 = %.4f  (folds %s)'
          % (len(keys), sc.mean(), ' '.join('%.3f' % s for s in sc)))
    sc = cv([keys.index('ac1')])
    print('ac1 alone, GBM            : CV R2 = %.4f' % sc.mean())
    # greedy nested selection, 6 rounds (same protocol as _nested_feats.py)
    chosen, remaining = [], list(range(len(keys)))
    for step in range(1, 7):
        best, bj = -9, None
        for j in remaining:
            s = cv(chosen + [j]).mean()
            if s > best:
                best, bj = s, j
        chosen.append(bj); remaining.remove(bj)
        print('  nested #%d += %-14s CV R2 = %.4f' % (step, keys[bj], best))


if __name__ == '__main__':
    main()
