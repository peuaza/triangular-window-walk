# -*- coding: utf-8 -*-
"""cut_ctl -- control: does the copied engine still reproduce the RETRACTION run?

probe_cycle's BTHM-C numbers were produced with trim=1 (a left window
retraction).  cut_engine.py is cyc_engine.py with one branch added (mode 4);
mode 1 is untouched.  This script re-runs mode 1, trim=1, on all 300 tapes of
every rule at n = 512 and compares (endC, stepsC, readsC, ntrim, cutsC) tape by
tape with probe_cycle/cyc_<rule>_512.json.  Any mismatch is a bug in the copy,
not a property of the repair.

usage: cut_ctl.py [n] [--jobs K]
writes cut_ctl.json
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OLD = os.path.join(os.path.dirname(HERE), 'probe_cycle')
FIELDS = ('endC', 'stepsC', 'readsC', 'ntrim', 'cutsC')


def safe(name):
    return (name.replace('>', 'g').replace('=', 'e').replace('/', 'd')
            .replace('.', 'p'))


def job(args):
    n, rule = args
    import sys
    sys.path.insert(0, HERE)
    import cut_engine as E
    import cut_verify as V
    import tr_lib as T
    old = json.load(open(os.path.join(OLD, 'cyc_%s_%d.json' % (safe(rule), n))))
    rows = old['rows']
    bad = []
    for i, a, s, X in T.tapes(n, len(rows)):
        C = E.Cyc(X, hlog=24)
        r = C.run(V.rrid(rule), V.conv(rule), V.variant(rule), mode=1, trim=1)
        got = dict(endC=r['end'], stepsC=r['steps'], readsC=r['reads'],
                   ntrim=r['ntrim'], cutsC=r['cuts'])
        for k in FIELDS:
            if int(rows[i][k]) != int(got[k]):
                bad.append((rule, i, k, rows[i][k], got[k]))
    return rule, len(rows) * len(FIELDS), bad


def main():
    n = 512
    jobs = 6
    argv = sys.argv[1:]
    if '--jobs' in argv:
        k = argv.index('--jobs')
        jobs = int(argv[k + 1])
        del argv[k:k + 2]
    if argv and argv[0].isdigit():
        n = int(argv[0])
    sys.path.insert(0, HERE)
    import cut_engine                                             # noqa: F401
    from tr_lib import RULENAME
    rules = [r for r in RULENAME if r != 'monotone']
    from concurrent.futures import ProcessPoolExecutor
    tot = 0
    bad = []
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        for rule, nf, b in ex.map(job, [(n, r) for r in rules]):
            tot += nf
            bad += b
            print('%-10s %5d fields, mismatches=%d' % (rule, nf, len(b)),
                  flush=True)
    print('TOTAL %d/%d fields exact' % (tot - len(bad), tot), flush=True)
    for b in bad[:40]:
        print('MISMATCH', b, flush=True)
    json.dump(dict(n=n, fields=tot, bad=bad),
              open(os.path.join(HERE, 'cut_ctl.json'), 'w'), indent=1)
    return 0 if not bad else 1


if __name__ == '__main__':
    sys.exit(main())
