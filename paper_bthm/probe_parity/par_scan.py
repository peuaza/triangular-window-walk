# -*- coding: utf-8 -*-
"""par_scan -- direct engine comparison mode 4 (flag) vs mode 5 (parity).

For every rule and every tape, the two runs must agree on every dynamical
field: steps, reads, cuts, end, ntrim, nl, nr, nused, msum.  (The diagnostic
fields cyc_start / cyc_len / width are populated by mode 4 only; mode 5 does
not stop on a repeat, so they are not compared.)  Any tape with a third visit
to one of its window states would make the repair disciplines differ at that
repeat and show up here as a mismatch.

usage: par_scan.py [n ...]
writes par_divscan.json
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import par_engine as E                                           # noqa: E402
import par_verify as V                                           # noqa: E402
import tr_lib as T                                               # noqa: E402

FIELDS = ('steps', 'reads', 'cuts', 'end', 'ntrim', 'nl', 'nr', 'nused',
          'msum')


def main():
    ns = [int(a) for a in sys.argv[1:] if a.isdigit()] or [256, 512]
    nt = 10
    t0 = time.time()
    out = dict(ns=ns, nt=nt, fields=list(FIELDS), rules={}, mismatches=0,
               first=None)
    for n in ns:
        for rule in [r for r, _, _, _, _ in T.SPEC if r != 'monotone']:
            rid, cv, vr = V.rrid(rule), V.conv(rule), V.variant(rule)
            bad = diff_tape = 0
            for i, a, s, X in T.tapes(n, nt):
                C = E.Cyc(X, hlog=20)
                r4 = C.run(rid, cv, vr, mode=4)
                r5 = C.run(rid, cv, vr, mode=5)
                hits = [k for k in FIELDS if r4[k] != r5[k]]
                if hits:
                    diff_tape += 1
                    bad += len(hits)
                    if out['first'] is None:
                        out['first'] = (rule, n, i, hits)
            out['rules']['%s|%d' % (rule, n)] = dict(tapes=nt,
                                                     field_mismatches=bad,
                                                     tapes_diverging=diff_tape)
            out['mismatches'] += bad
            print('%-10s n=%-5d field mismatches %d (tapes diverging %d/10)'
                  % (rule, n, bad, diff_tape), flush=True)
    out['secs'] = round(time.time() - t0, 1)
    json.dump(out, open(os.path.join(HERE, 'par_divscan.json'), 'w'), indent=1)
    print('total field mismatches:', out['mismatches'], '| secs', out['secs'])
    return 0 if out['mismatches'] == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
