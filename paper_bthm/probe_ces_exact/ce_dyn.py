# -*- coding: utf-8 -*-
"""ce_dyn -- anatomy of the early-terminal epochs (the ones CES gets wrong).

For every tape whose true terminal epoch has M0 > W, replay that epoch and
record: type, M0, final m, length, number of direction changes, how much of
the epoch is a one-directional shrink, and the final states.

usage: ce_dyn.py <rule> <n> <nt> [W]
"""
import json
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (HERE, os.path.join(ROOT, 'probe_e_uniform'),
           os.path.join(ROOT, 'probe_kappa_pred'),
           os.path.join(ROOT, 'probe_ologn2')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import ce_lib as CL                                              # noqa: E402
import eu_lib as EL                                              # noqa: E402


def scan_full(a, L, R, rule, maxsteps=100000):
    from rs_common import Capsule
    C = Capsule(a)
    ell, r = L, R
    seen = {}
    path = []
    for t in range(maxsteps):
        m = r - ell + 1
        b2 = C.b2(ell, m); b1 = C.b1(ell, m, rule)
        if b2:
            dy, dx = (1, 0) if b1 else (0, -1)
        else:
            dy, dx = (0, 1) if b1 else (-1, 0)
        ev = 'move'
        if ell + dy < L or r + dx > R:
            ev = 'D' if (ell == L and r == R) else ('L' if ell + dy < L else 'R')
        elif ell + dy >= r + dx:
            ev = 'meet'
        elif r + dx < L or ell + dy > R:
            ev = 'D'
        path.append((ell, r, m, int(b2), int(b1), ev))
        if ev != 'move':
            return path, ev, None
        key = (ell, r)
        if key in seen:
            return path, 'cycle', t - seen[key]
        seen[key] = t
        ell, r = ell + dy, r + dx
    return path, 'cap', None


def main():
    rule = sys.argv[1]; n = int(sys.argv[2])
    nt = int(sys.argv[3]) if len(sys.argv) > 3 else 220
    W = int(sys.argv[4]) if len(sys.argv) > 4 else (14 if rule == 'pattern' else 40)
    stats = []
    types = Counter()
    for i in range(nt):
        a, s, X = EL.make_tape(n, 5000 + n + 100003 * i)
        a = np.asarray(a, np.int64)
        d, why, eps = CL.trace(a, rule)
        M0, b2, b1, oc, ln, L, R, ee = eps[-1]
        if oc in ('D', 'R', 'L') or M0 <= W:
            continue
        path, ev, cyc = scan_full(a, L, R, rule)
        moves = [p for p in path if p[5] == 'move']
        dirs = [1 if p[0] != q[0] else (-1 if p[1] != q[1] else 0)
                for p, q in zip(moves, moves[1:])]
        ndir = sum(1 for x, y in zip(dirs, dirs[1:]) if x != y)
        # how long is the final monotone run?
        run = 1
        for x, y in zip(reversed(dirs[:-1]), reversed(dirs[1:])):
            if x == y:
                run += 1
            else:
                break
        stats.append(dict(i=i, why=why, M0=int(M0), b2=int(b2), b1=int(b1),
                          oc=oc, ln=int(ln), mfinal=int(path[-1][2]),
                          ndir=int(ndir), finalrun=int(run), cyc=cyc,
                          path=[(int(p[0]), int(p[1]), int(p[2]), p[3], p[4], p[5]) for p in path]))
        types[(oc, int(b2), int(b1))] += 1
    out = dict(rule=rule, n=n, W=W, nt=nt, types={str(k): v for k, v in types.items()},
               epochs=stats)
    with open(os.path.join(HERE, 'ce_dyn_%s_%d.json' % (rule, n)), 'w') as f:
        json.dump(out, f, indent=1)
    print('%s n=%d: %d early-terminal epochs, types %s' % (rule, n, len(stats), dict(types)))
    if stats:
        import statistics as st
        for k in ['M0', 'ln', 'mfinal', 'ndir', 'finalrun']:
            v = [s[k] for s in stats]
            print('  %-8s min %4d med %6.1f mean %7.1f max %4d'
                  % (k, min(v), st.median(v), st.mean(v), max(v)))
        print('  cycles: %d, of which 2-cycles %d; meet %d'
              % (sum(1 for s in stats if s['oc'] == 'cycle'),
                 sum(1 for s in stats if s['cyc'] == 2),
                 sum(1 for s in stats if s['oc'] == 'meet')))


if __name__ == '__main__':
    main()
