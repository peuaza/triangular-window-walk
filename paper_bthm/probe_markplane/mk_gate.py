# -*- coding: utf-8 -*-
"""mk_gate -- exactness gate, run before any cell is measured.

For every (rule, n) and every variant:
  * policy 'none' must reproduce the baseline machine's reads, delta, end, steps
    (and cuts where the machine has them) on every tape;
  * every variant must reproduce the baseline's (delta, end) on every tape;
  * any variant that does not is reported and DISCARDED (never tabulated).

usage:  mk_gate.py [ntapes]      (default 8; the full 300-tape gate is mk_run's
                                  own check, which compares every tape it runs)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PAPER = os.path.dirname(HERE)
for d in ('probe_triple_nat', 'probe_e_uniform'):
    p = os.path.join(PAPER, d)
    if p not in sys.path:
        sys.path.insert(0, p)

import tr_lib as T                                              # noqa: E402
import mk_pal                                                   # noqa: E402
import mk_pat                                                   # noqa: E402


def same(x, y):
    if x != y:
        if isinstance(x, float) and isinstance(y, float):
            return (x != x and y != y) or abs(x - y) <= 1e-12
        return False
    return True


def gate_pal(n, nt=8, policies=('cert', 'cert_naive', 'run')):
    rows = []
    for i, a, s, X in T.tapes(n, nt):
        b = mk_pal.run_base(a)
        m = mk_pal.run_pal(a, 'none')
        ok_none = (m['reads'] == b['reads'] and same(m['delta'], b['delta'])
                   and m['end'] == b['end'] and m['steps'] == b['steps'])
        row = dict(i=i, base=b, none_ok=ok_none, var={})
        for p in policies:
            r = mk_pal.run_pal(a, p)
            row['var'][p] = (r, same(r['delta'], b['delta']) and r['end'] == b['end'])
        rows.append(row)
    return rows


def gate_pat(n, nt=8, policies=('adj', 'occ', 'occ_null')):
    rows = []
    for i, a, s, X in T.tapes(n, nt):
        b = mk_pat.run_base(X)
        m = mk_pat.run_pat(X, 'none')
        ok_none = (m['reads'] == b['reads'] and same(m['delta'], b['delta'])
                   and m['end'] == b['end'] and m['steps'] == b['steps']
                   and m['cuts'] == b['cuts'])
        row = dict(i=i, base=b, none_ok=ok_none, var={})
        for p in policies:
            r = mk_pat.run_pat(X, p)
            row['var'][p] = (r, same(r['delta'], b['delta']) and r['end'] == b['end'])
        rows.append(row)
    return rows


if __name__ == '__main__':
    nt = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    for n in (512, 1024):
        for name, rows in (('pal', gate_pal(n, nt)), ('pattern', gate_pat(n, nt))):
            bad_none = [r['i'] for r in rows if not r['none_ok']]
            print('%-8s n=%-5d tapes=%-3d policy-none mismatches: %d %s'
                  % (name, n, len(rows), len(bad_none), bad_none[:5]))
            for p in rows[0]['var']:
                bad = [r['i'] for r in rows if not r['var'][p][1]]
                rd = np.mean([r['var'][p][0]['reads'] for r in rows])
                rd0 = np.mean([r['base']['reads'] for r in rows])
                ink = np.mean([r['var'][p][0]['ink'] for r in rows])
                print('    %-10s reads=%-10.0f (base %-10.0f, x%.3f) ink=%-9.0f '
                      'exact %d/%d %s'
                      % (p, rd, rd0, rd / max(rd0, 1e-9), ink,
                         len(rows) - len(bad), len(rows), bad[:4]))
