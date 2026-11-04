# -*- coding: utf-8 -*-
"""p1: exhaustive enumeration of ensemble tapes; exact delta; test one-pass candidates.

Ensemble: a[i] = sum of two bits, each of two columns has exactly ceil(n/2) ones.
=> multiset has #0 == #2 == t, #1 == n-2t  (n even: sum = n; n odd: sum = n+1)
"""
import sys, itertools, json
import numpy as np
sys.path.insert(0, '..')
from _arr_longrun_fast import walk_naive, bit_b2_ref, bit_b1_ref, longest_run_vec


def enum_tapes(n):
    for t in range(n // 2 + 1):
        c1 = n - 2 * t
        if c1 < 0:
            continue
        for zeros in itertools.combinations(range(n), t):
            zs = set(zeros)
            rest = [i for i in range(n) if i not in zs]
            for twos in itertools.combinations(rest, t):
                a = np.ones(n, dtype=np.int64)
                a[list(zeros)] = 0
                a[list(twos)] = 2
                yield a


def cand_feats(a):
    n = len(a)
    x = a.astype(np.int64) - 1          # balanced walk (sums to 0 for even n, 1 for odd)
    P = np.concatenate(([0], np.cumsum(x)))
    P2 = np.concatenate(([0], np.cumsum((a == 2).astype(np.int64) - 0.25)))
    out = {}
    out['argmax_P_first'] = int(np.argmax(P))
    out['argmax_P_last'] = int(len(P) - 1 - np.argmax(P[::-1]))
    out['argmin_P_first'] = int(np.argmin(P))
    out['argmin_P_last'] = int(len(P) - 1 - np.argmin(P[::-1]))
    out['argmax_P2_first'] = int(np.argmax(P2))
    out['argmax_P2_last'] = int(len(P2) - 1 - np.argmax(P2[::-1]))
    out['argmin_P2_first'] = int(np.argmin(P2))
    out['argmin_P2_last'] = int(len(P2) - 1 - np.argmin(P2[::-1]))
    # centre of mass of 2s / 0s
    tw = np.flatnonzero(a == 2); zr = np.flatnonzero(a == 0)
    out['cm_2'] = int(round(tw.mean())) if tw.size else 0
    out['cm_0'] = int(round(zr.mean())) if zr.size else 0
    out['first_2'] = int(tw[0]) if tw.size else 0
    out['last_2'] = int(tw[-1]) if tw.size else n - 1
    out['first_0'] = int(zr[0]) if zr.size else 0
    out['last_0'] = int(zr[-1]) if zr.size else n - 1
    return out


def main():
    n = int(sys.argv[1])
    tapes = list(enum_tapes(n))
    print('n=%d  ensemble tapes = %d' % (n, len(tapes)))
    deltas = []
    steps = []
    for a in tapes:
        d, s, _ = walk_naive(a, n)
        deltas.append(d)
        steps.append(s)
    nd = sum(1 for d in deltas if d is not None)
    print('  with delta: %d (%.3f)   steps/tape mean %.1f' % (nd, nd / len(tapes), np.mean(steps)))
    # distribution of delta index
    idx = [int(round(d * (n - 1))) for d in deltas if d is not None]
    vals, cnt = np.unique(idx, return_counts=True) if idx else ([], [])
    print('  delta-index histogram: %s' % dict(zip([int(v) for v in vals], [int(c) for c in cnt])))
    # value at the delta cell
    for v in (0, 1, 2):
        c = sum(1 for a, d in zip(tapes, deltas) if d is not None and a[int(round(d * (n - 1)))] == v)
        print('    a[delta]==%d : %d' % (v, c))
    # candidate features: exact match rate on tapes with delta
    keys = list(cand_feats(tapes[0]).keys())
    hits = {k: 0 for k in keys}
    for a, d in zip(tapes, deltas):
        if d is None:
            continue
        f = cand_feats(a)
        for k in keys:
            if f[k] == int(round(d * (n - 1))):
                hits[k] += 1
    print('  one-pass candidates (exact hits / %d):' % nd)
    for k in keys:
        print('    %-18s %d  (%.3f)' % (k, hits[k], hits[k] / max(1, nd)))
    # save
    with open('p1_n%d.json' % n, 'w') as fh:
        json.dump({'n': n, 'tapes': [a.tolist() for a in tapes],
                   'delta': deltas, 'steps': steps}, fh)


if __name__ == '__main__':
    main()
