# -*- coding: utf-8 -*-
"""cyc_record_check -- BTHM-0's meeting count against the published tr_ cells.

For every catalogue rule whose tr_lib machine is the four-move walk, the set of
tapes that meet cannot depend on the machine: a meeting tape never repeats a
window, so no memory is needed to get there.  So BTHM-0's meet count must equal
tr_<rule>_<n>.json's.  Any cell where it does not is a disagreement about the
walk itself, not about cycles.
"""
import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
TAN = os.path.join(PAPER, 'probe_triple_nat')
import cyc_engine                                                  # noqa: F401
from tr_lib import RULENAME

SAFE = lambda s: (s.replace('>', 'g').replace('=', 'e').replace('/', 'd')
                  .replace('.', 'p'))


def norm(e):
    if isinstance(e, str):
        return e
    return {0: 'meet', 1: 'empty', 2: 'repeat', 3: 'cap', None: 'none'}.get(e, str(e))


def main():
    out = {}
    for n in (512, 1024, 2048, 4096):
        for rule in RULENAME:
            rec = os.path.join(TAN, 'tr_%s_%d.json' % (SAFE(rule), n))
            cyc = os.path.join(HERE, 'cyc_%s_%d.json' % (SAFE(rule), n))
            if not (os.path.exists(rec) and (os.path.exists(cyc) or rule == 'monotone')):
                continue
            R = json.load(open(rec))
            C = json.load(open(cyc)) if os.path.exists(cyc) else None
            if not R.get('end') or len(R['end']) != 300:
                continue
            met_rec = sum(1 for e in R['end'] if norm(e) == 'meet')
            if rule == 'monotone':
                met_cyc = 300          # analytic: no cut, ell += 1, always meets
            else:
                met_cyc = sum(1 for x in C['rows'] if x['end0'] == 0)
            out['%s@%d' % (rule, n)] = dict(rec=met_rec, cyc=met_cyc,
                                            agree=met_rec == met_cyc)
    bad = {k: v for k, v in out.items() if not v['agree']}
    json.dump(out, open(os.path.join(HERE, 'cyc_record_check.json'), 'w'), indent=1)
    print('cells compared: %d, disagreements: %d' % (len(out), len(bad)))
    for k, v in bad.items():
        print('  %-16s record %3d  bthm0 %3d' % (k, v['rec'], v['cyc']))


if __name__ == '__main__':
    main()
