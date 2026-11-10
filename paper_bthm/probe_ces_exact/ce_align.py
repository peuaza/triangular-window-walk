# -*- coding: utf-8 -*-
"""ce_align -- CES vs truth vs CES-with-oracle-outcomes, tape by tape.

For each tape:
  truth  : kp_lib.exact (delta or nan)
  cheap  : kp_lib.ces(W) with the (b2,b1) outcome rule
  oracle : kp_lib.ces(W, oracle_oc=true cut outcomes)
Classifies agreement and prints read counts.  Isolates the outcome-rule error
from the endgame-buffer error.

usage: ce_align.py <rule> <n> <nt> [W]
Writes ce_align_<rule>_<n>.json
"""
import json
import os
import sys

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
import kp_lib as KP                                              # noqa: E402


def Wof(rule, n):
    if rule == 'pattern':
        return 14
    return max(16, 4 * int(round(np.log2(n))))


def cls(truth, got):
    tn = truth is None or truth != truth
    gn = got is None or got != got
    if tn and gn:
        return 'both_nan'
    if tn and not gn:
        return 'invented'
    if not tn and gn:
        return 'declined'
    if abs(truth - got) < 1e-12:
        return 'exact'
    return 'wrong_delta'


def one(rule, a, n, W):
    d, why, eps = CL.trace(a, rule)
    truth = None if (d is None or d != d) else float(d)
    oc = [e[3] for e in eps if e[3] in ('D', 'R', 'L')]
    dc, rc, infc = KP.ces(a, rule, W=W)
    do, ro, info = KP.ces(a, rule, W=W, oracle_oc=oc)
    return dict(why=why, truth=truth, cheap=None if dc != dc else float(dc),
                oracle=None if do != do else float(do), reads_cheap=int(rc),
                reads_oracle=int(ro), eps=len(eps),
                end_cheap=infc.get('end'), end_oracle=info.get('end'))


def main():
    rule = sys.argv[1]
    n = int(sys.argv[2])
    nt = int(sys.argv[3]) if len(sys.argv) > 3 else 220
    W = int(sys.argv[4]) if len(sys.argv) > 4 else Wof(rule, n)
    recs = []
    cnt = {'cheap': {}, 'oracle': {}}
    reads_c = []; reads_o = []
    for i in range(nt):
        a, s, X = EL.make_tape(n, 5000 + n + 100003 * i)
        a = np.asarray(a, np.int64)
        r = one(rule, a, n, W)
        recs.append(r)
        if r['truth'] is not None or True:
            cnt['cheap'][cls(r['truth'], r['cheap'])] = cnt['cheap'].get(cls(r['truth'], r['cheap']), 0) + 1
            cnt['oracle'][cls(r['truth'], r['oracle'])] = cnt['oracle'].get(cls(r['truth'], r['oracle']), 0) + 1
        reads_c.append(r['reads_cheap']); reads_o.append(r['reads_oracle'])
    out = dict(rule=rule, n=n, nt=nt, W=W, counts=cnt,
               reads_cheap_mean=float(np.mean(reads_c)),
               reads_oracle_mean=float(np.mean(reads_o)),
               recs=recs)
    with open(os.path.join(HERE, 'ce_align_%s_%d.json' % (rule, n)), 'w') as f:
        json.dump(out, f, indent=1)
    print('%s n=%d nt=%d W=%d' % (rule, n, nt, W))
    print('  cheap : %s  reads/n %.2f' % (cnt['cheap'], np.mean(reads_c) / n))
    print('  oracle: %s  reads/n %.2f' % (cnt['oracle'], np.mean(reads_o) / n))


if __name__ == '__main__':
    main()
