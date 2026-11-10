# -*- coding: utf-8 -*-
"""ce_patch -- reads and exactness of CES and its cheap patches.

Machines, all read-counted, all O(log n) bits:
  ces        pure CES (W = 14 pattern, 4 log2 n runs)
  k4         CES + prefix-4 verification at every shortcut epoch
  oracle     CES + the terminal epoch run exactly (diagnostic ceiling)
  w2_<T>     CES + every epoch with M <= T run exactly

Tapes: T_n (eu_lib.make_tape).  Truth: pcr_lib.Runner (numba), gated against
kp_lib.exact.  Per cell: up to 200 meeting tapes and 100 halting tapes.

usage: ce_patch.py <rule> <n> [nt_meet] [nt_halt] [variants]
Writes ce_patch_<rule>_<n>.json
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for _p in (HERE, os.path.join(ROOT, 'probe_e_uniform'),
           os.path.join(ROOT, 'probe_kappa_pred'),
           os.path.join(ROOT, 'probe_ologn2'),
           os.path.join(ROOT, 'probe_class_rerun')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import ce_lib as CL                                              # noqa: E402
import ce_machine as CM                                          # noqa: E402
import eu_lib as EL                                              # noqa: E402

_P = None


def pcr():
    global _P
    if _P is None:
        import pcr_lib as P
        _P = P
    return _P


def truth(X, rule):
    P = pcr()
    r = P.Runner(X).run(rule)
    out, steps, msum, inc, cuts, end, snaps = r
    d = None if out < 0 else float(out)
    return d, int(steps), int(cuts), int(end)


def Wof(rule, n):
    return 14 if rule == 'pattern' else max(16, 4 * int(round(np.log2(n))))


def cls(t, g):
    tn = t is None; gn = g is None or g != g
    if tn and gn:
        return 'both_nan'
    if tn:
        return 'invented'
    if gn:
        return 'declined'
    return 'exact' if abs(t - g) < 1e-12 else 'wrong'


def variants_for(n):
    v = [('ces', dict()), ('k4', dict(K=4)), ('oracle', dict())]
    for T in (32, 64, 128, 256, 512, 1024):
        if T < n:
            v.append(('w2_%d' % T, dict(W2=T)))
    for name, T in (('w2_q', n // 4), ('w2_half', n // 2),
                    ('w2_3q', 3 * n // 4), ('w2_n', n)):
        if T >= 3:
            v.append((name, dict(W2=T)))
    return v


def main():
    rule = sys.argv[1]
    n = int(sys.argv[2])
    nt_meet = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    nt_halt = int(sys.argv[4]) if len(sys.argv) > 4 else 100
    only = sys.argv[5].split(',') if len(sys.argv) > 5 else None
    W = Wof(rule, n)
    # ---- collect tapes: first nt_meet meeting, first nt_halt halting
    tapes = []      # (tag, a, X, truth, steps, cuts, end)
    nmeet = nhalt = 0
    i = 0
    cap = 20 * (nt_meet + nt_halt) + 200
    while (nmeet < nt_meet or nhalt < nt_halt) and i < cap:
        a, s, X = EL.make_tape(n, 5000 + n + 100003 * i)
        d, st, cuts, end = truth(X, rule)
        tag = None
        if d is not None and nmeet < nt_meet:
            tag = 'meet'; nmeet += 1
        elif d is None and nhalt < nt_halt:
            tag = 'halt'; nhalt += 1
        if tag:
            tapes.append((tag, np.asarray(a, np.int64), X, d, st, cuts, end))
        i += 1
        if i % 50 == 0:
            print('  scanned %d: meet %d halt %d' % (i, nmeet, nhalt), flush=True)
    print('%s n=%d W=%d: %d meeting, %d halting out of %d scanned'
          % (rule, n, W, nmeet, nhalt, i), flush=True)
    # ---- machines
    res = {}
    for name, kw in variants_for(n):
        if only is not None and name not in only:
            continue
        rows = []
        for (tag, a, X, d, st, cuts, end) in tapes:
            orc = None
            if name == 'oracle':
                orc = {cuts}
            kw2 = dict(kw)
            if name == 'oracle':
                kw2['oracle_term'] = orc
            g, rd, inf = CM.machine(a, rule, W=W, **kw2)
            rows.append(dict(tag=tag, reads=int(rd), got=None if (g is None or g != g) else float(g),
                             truth=d, cls=cls(d, g), term=inf.get('terminal'),
                             epochs=inf.get('epochs'), exact_epochs=inf.get('exact_epochs')))
        meet = [r for r in rows if r['tag'] == 'meet']
        halt = [r for r in rows if r['tag'] == 'halt']
        rm = np.mean([r['reads'] for r in meet]) if meet else float('nan')
        ra = np.mean([r['reads'] for r in rows])
        em = sum(r['cls'] == 'exact' for r in meet)
        ea = sum(r['cls'] in ('exact', 'both_nan') for r in rows)
        res[name] = dict(reads_meet_mean=float(rm), reads_meet_n=float(rm / n),
                         reads_all_mean=float(ra), e_n=float(ra / n),
                         exact_meet=int(em), n_meet=len(meet),
                         exact_all=int(ea), n_all=len(rows),
                         counts_meet={k: sum(r['cls'] == k for r in meet)
                                      for k in ('exact', 'invented', 'declined', 'wrong', 'both_nan')},
                         counts_halt={k: sum(r['cls'] == k for r in halt)
                                      for k in ('exact', 'invented', 'declined', 'wrong', 'both_nan')},
                         rows=rows)
        r = res[name]
        print('  %-8s e/n %6.2f  exact meet %3d/%3d  exact all %3d/%3d  halt-invented %d'
              % (name, r['e_n'], r['exact_meet'], r['n_meet'], r['exact_all'],
                 r['n_all'], r['counts_halt'].get('invented', 0)), flush=True)
    fn = os.path.join(HERE, 'ce_patch_%s_%d.json' % (rule, n))
    prev = {}
    if os.path.exists(fn):
        with open(fn) as f:
            prev = json.load(f).get('machines', {})
    prev.update(res)
    out = dict(rule=rule, n=n, W=W, nt_meet=nmeet, nt_halt=nhalt,
               scanned=i, machines=prev)
    with open(fn, 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
