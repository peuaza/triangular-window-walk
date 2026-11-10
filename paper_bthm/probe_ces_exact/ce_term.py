# -*- coding: utf-8 -*-
"""ce_term -- what do the early-terminal (CES-failing) epochs look like?

Reads ce_fail_<rule>_<n>.json (epoch records) and prints, for a sample of
failing tapes, the region/window at the terminal epoch and the exact scan
trajectory inside that epoch.  Also tabulates P(terminal | b2,b1,M0 bucket).

usage: ce_term.py <rule> <n> [maxshow]
"""
import json
import os
import sys
from collections import Counter, defaultdict

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


def scan_trace(a, L, R, rule, maxsteps=2000):
    """Literal epoch trajectory from (L,R): list of (ell,r,m,b2,b1,event)."""
    from rs_common import Capsule
    C = Capsule(a)
    ell, r = L, R
    seen = set()
    out = []
    for _ in range(maxsteps):
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
        out.append((ell, r, m, int(b2), int(b1), ev))
        if ev != 'move':
            return out, ev
        key = (ell, r)
        if key in seen:
            return out, 'cycle'
        seen.add(key)
        ell, r = ell + dy, r + dx
    return out, 'cap'


def main():
    rule = sys.argv[1]; n = int(sys.argv[2])
    maxshow = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    J = json.load(open(os.path.join(HERE, 'ce_fail_%s_%d.json' % (rule, n))))
    W = J['W']
    eps_by_tape = [r for r in J['recs']]
    # conditional terminal rate by (b2,b1) and M0 bucket
    tot = defaultdict(int); term = defaultdict(int)
    ttypes = Counter(); tln = defaultdict(list)
    for r in eps_by_tape:
        pass
    # need epoch records: recompute from tapes (JSON has only summaries in recs)
    shown = 0
    for i in range(J['nt']):
        a, s, X = EL.make_tape(n, 5000 + n + 100003 * i)
        a = np.asarray(a, np.int64)
        d, why, eps = CL.trace(a, rule)
        for e in eps:
            M0, b2, b1, oc, ln, L, R, ee = e
            key = (int(b2), int(b1), min(M0 // 32, 20))
            tot[key] += 1
            if oc not in ('D', 'R', 'L'):
                term[key] += 1
                ttypes[(oc, int(b2), int(b1))] += 1
                tln[oc].append(ln)
        # show the first few tapes whose terminal sits above W
        last = eps[-1]
        M0, b2, b1, oc, ln, L, R, ee = last
        if oc not in ('D', 'R', 'L') and M0 > W and shown < maxshow:
            shown += 1
            win = ''.join(str(int(x)) for x in a[L:R + 1])
            print('\n=== t=%d rule=%s why=%s delta=%s  terminal M0=%d L=%d R=%d '
                  'b2=%d b1=%d oc=%s len=%d' % (i, rule, why, d, M0, L, R, b2, b1, oc, ln))
            print('window: %s' % win)
            tr, ev = scan_trace(a, L, R, rule)
            print('steps=%d end=%s' % (len(tr), ev))
            for st in tr[:min(24, len(tr))]:
                print('   ell=%3d r=%3d m=%3d b2=%d b1=%d %s' % st)
            if len(tr) > 24:
                print('   ...')
                for st in tr[-6:]:
                    print('   ell=%3d r=%3d m=%3d b2=%d b1=%d %s' % st)
    print('\nterminal counts by (oc,b2,b1): %s' % dict(ttypes))
    for oc in tln:
        v = np.array(tln[oc])
        print('  len[%s]: n=%d mean %.1f med %.1f p90 %.0f max %d'
              % (oc, len(v), v.mean(), np.median(v), np.percentile(v, 90), v.max()))
    print('\nP(terminal | b2,b1,M0//32):')
    for key in sorted(tot):
        if tot[key] >= 20:
            print('  b2=%d b1=%d mb=%2d  n=%6d term=%5d p=%.4f'
                  % (key[0], key[1], key[2], tot[key], term[key], term[key] / tot[key]))


if __name__ == '__main__':
    main()
