# -*- coding: utf-8 -*-
"""p6: exact certificate sizes on the ensemble.

For every subset S of cells, project all tapes onto S and ask whether delta is
constant inside every projection class.
  K_global(n) = smallest |S| for which SOME S is pure for all tapes
                (a fixed set of cells that determines delta for every tape)
  c(a)        = smallest |S| such that S is pure on the class of a
                (a certificate for that one tape; every correct simulator on
                 input a must read at least c(a) distinct cells)

Usage: python p6_cert.py <n> [n2 ...]
"""
import sys, itertools, time, json
import numpy as np
sys.path.insert(0, '..')
from p0_tapes import enum_tapes
from _arr_longrun_fast import walk_naive


def load(n):
    A, D = [], []
    for a in enum_tapes(n):
        d, s, _ = walk_naive(a, n)
        A.append(a)
        D.append(-1 if d is None else int(round(d * (n - 1))))
    return np.array(A, dtype=np.int8), np.array(D, dtype=np.int16)


def run(n):
    t0 = time.time()
    A, D = load(n)
    N = len(A)
    unresolved = np.ones(N, dtype=bool)
    cert = np.full(N, -1, dtype=np.int16)
    K_global = None
    pw = (3 ** np.arange(n)).astype(np.int64)
    stats = []
    for k in range(1, n + 1):
        n_pure_sets = 0
        for S in itertools.combinations(range(n), k):
            key = A[:, S].astype(np.int64) @ pw[:k]
            order = np.argsort(key, kind='stable')
            ks = key[order]
            ds = D[order]
            b = np.empty(N, dtype=bool)
            b[0] = True
            b[1:] = ks[1:] != ks[:-1]
            gid = np.cumsum(b) - 1
            G = int(gid[-1]) + 1
            mn = np.full(G, 32767, dtype=np.int16)
            mx = np.full(G, -32768, dtype=np.int16)
            np.minimum.at(mn, gid, ds)
            np.maximum.at(mx, gid, ds)
            impure_g = mn != mx
            is_imp = impure_g[gid]
            pure = np.empty(N, dtype=bool)
            pure[order] = ~is_imp
            if pure.all():
                n_pure_sets += 1
                if K_global is None:
                    K_global = k
            new = pure & unresolved
            if new.any():
                cert[new] = k
                unresolved[new] = False
        stats.append((k, n_pure_sets))
        print('  n=%d k=%d  globally pure subsets: %d   unresolved tapes: %d   (%.1fs)'
              % (n, k, n_pure_sets, int(unresolved.sum()), time.time() - t0), flush=True)
        if not unresolved.any() and K_global is not None:
            break
    print('n=%d tapes=%d  K_global=%s   c(a): min %d median %d mean %.2f max %d  (c/n: median %.3f max %.3f)  %.1fs'
          % (n, N, K_global, cert.min(), np.median(cert), cert.mean(), cert.max(),
             np.median(cert) / n, cert.max() / n, time.time() - t0))
    np.save('cert_n%d.npy' % n, np.stack([D, cert]))
    return K_global, cert


if __name__ == '__main__':
    out = {}
    for n in (int(x) for x in sys.argv[1:]):
        K, cert = run(n)
        out[n] = (int(K) if K else None, int(cert.min()), float(np.median(cert)), float(cert.mean()),
                  int(cert.max()), len(cert))
    with open('cert_summary.json', 'w') as fh:
        json.dump({str(k): v for k, v in out.items()}, fh, indent=1)
    print(json.dumps({str(k): v for k, v in out.items()}, indent=1))
