# -*- coding: utf-8 -*-
"""_spot -- does mode 0 in cut_engine.py reproduce probe_cycle's stored rows?

Compares every field cyc_run.py stores for BTHM-0 (end0, steps0, reads0,
cyc_start, cyc_len, width, cuts0) tape by tape against
probe_cycle/cyc_<rule>_<n>.json for a few cells.  Nothing is written.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cut_engine as E                                            # noqa: E402
import cut_verify as V                                            # noqa: E402
import tr_lib as T                                                # noqa: E402

OLD = os.path.join(os.path.dirname(HERE), 'probe_cycle')
FIELDS = ('end0', 'steps0', 'reads0', 'cyc_start', 'cyc_len', 'width', 'cuts0')


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 512
    nt = int(sys.argv[2]) if len(sys.argv) > 2 else 25
    rules = sys.argv[3:] or ['canonical', 'dense2', 'NBLK>2']
    for rule in rules:
        f = os.path.join(OLD, 'cyc_%s_%d.json' % (
            rule.replace('>', 'g').replace('=', 'e').replace('/', 'd'), n))
        old = json.load(open(f))['rows']
        bad = 0
        for i, a, s, X in T.tapes(n, nt):
            C = E.Cyc(X)
            r = C.run(V.rrid(rule), V.conv(rule), V.variant(rule), mode=0)
            got = dict(end0=r['end'], steps0=r['steps'], reads0=r['reads'],
                       cyc_start=r['cyc_start'], cyc_len=r['cyc_len'],
                       width=r['width'], cuts0=r['cuts'])
            for k in FIELDS:
                if int(old[i][k]) != int(got[k]):
                    bad += 1
                    print('MISMATCH %s tape %d %s old=%s new=%s'
                          % (rule, i, k, old[i][k], got[k]), flush=True)
        print('%-9s n=%d  %d tapes x %d fields, mismatches=%d'
              % (rule, n, nt, len(FIELDS), bad), flush=True)


if __name__ == '__main__':
    main()
