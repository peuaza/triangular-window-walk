# -*- coding: utf-8 -*-
"""cyc_period -- is the detected repeat a real period?

Mode 3 reruns the walk past the first repeat: the second repeat must be of the
state that was current one period after the first occurrence, i.e.
s2start == cstart + clen and s2step == cstart + 2*clen.  A projector that only
compared (ell, r) could not pass this; a genuine quadruple period must.
"""
import sys

import numpy as np

import cyc_engine as E
import cyc_verify as V
import tr_lib as T

n = int(sys.argv[1]) if len(sys.argv) > 1 else 512
NT = int(sys.argv[2]) if len(sys.argv) > 2 else 40
bad = 0
tot = 0
per = {}
for rule in T.RULENAME:
    if rule == 'monotone':
        continue
    ok = 0
    cls = []
    for i, a, s, X in T.tapes(n, NT):
        C = E.Cyc(X)
        r = C.run(V.rrid(rule), V.conv(rule), V.variant(rule), mode=3)
        if r['end'] != 0 and r['end'] != 1:
            tot += 1
            good = (r['s2start'] == r['cyc_start']
                    and r['s2step'] == r['cyc_start'] + 2 * r['cyc_len'])
            if good:
                ok += 1
                cls.append(r['cyc_len'])
            else:
                bad += 1
                print('PERIOD FAIL', rule, i, r['cyc_start'], r['cyc_len'],
                      r['s2start'], r['s2step'], r['end'], flush=True)
        else:
            tot += 1
            good = (r['s2start'] < 0)
            if good:
                ok += 1
            else:
                bad += 1
    per[rule] = (ok, cls)
print('period check n=%d nt=%d: %d/%d cells consistent' % (n, NT, tot - bad, tot))
for rule, (ok, cls) in per.items():
    if cls:
        print('%-10s periods: %s' % (rule, sorted(cls)[:12]))

import io
import json
import os

_out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cyc_period.json")
io.open(_out, "w", encoding="utf-8", newline="\n").write(json.dumps(
    {"n": n, "nt": NT, "consistent": tot - bad, "total": tot, "bad": bad,
     "per_rule": {r: {"ok": ok, "periods": sorted(cls)} for r, (ok, cls) in per.items()}},
    indent=1, ensure_ascii=False) + "\n")
print("wrote cyc_period.json: %d/%d consistent" % (tot - bad, tot))
