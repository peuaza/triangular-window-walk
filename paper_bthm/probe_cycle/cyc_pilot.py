# -*- coding: utf-8 -*-
"""cyc_pilot -- how expensive is the measurement, and does anything cycle?"""
import sys
import time

import numpy as np

import cyc_engine as E
import cyc_verify as V

n = int(sys.argv[1]) if len(sys.argv) > 1 else 512
NT = int(sys.argv[2]) if len(sys.argv) > 2 else 20
import tr_lib as T

print('%5s %6s %5s %5s %5s %8s %8s %8s %6s' %
      ('rule', 'meet', 'cyc', 'empty', 'cap', 'med_rep', 'max_rep', 'med_cyc', 'secs'))
for rule in T.RULENAME:
    if rule == 'monotone':
        continue
    t0 = time.time()
    c = dict(meet=0, cyc=0, empty=0, cap=0, other=0)
    rep, cyc, wid = [], [], []
    for i, a, s, X in T.tapes(n, NT):
        C = E.Cyc(X)
        r = C.run(V.rrid(rule), V.conv(rule), V.variant(rule), mode=0)
        e = r['end']
        if e == 0:
            c['meet'] += 1
        elif e == 4:
            c['cyc'] += 1
            rep.append(r['steps'])
            cyc.append(r['cyc_len'])
            wid.append(r['width'])
        elif e == 1:
            c['empty'] += 1
        elif e == 3:
            c['cap'] += 1
        else:
            c['other'] += 1
    dt = time.time() - t0
    print('%5s %6d %5d %5d %5d %8.0f %8.0f %8.0f %6.1f' %
          (rule[:5], c['meet'], c['cyc'], c['empty'], c['cap'],
           np.median(rep) if rep else -1, np.max(rep) if rep else -1,
           np.median(cyc) if cyc else -1, dt), flush=True)
