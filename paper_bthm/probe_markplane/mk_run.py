# -*- coding: utf-8 -*-
"""mk_run -- the 300-tape sweep: pal and pattern, n = 512/1024/2048, every
variant, with the exactness gate applied to every tape.

usage:  mk_run.py [ntapes] [n1,n2,...] [rule1,rule2] [procs]
Writes mk_cells.json next to this file.
"""
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PAPER = os.path.dirname(HERE)
for d in ('probe_triple_nat', 'probe_e_uniform'):
    p = os.path.join(PAPER, d)
    if p not in sys.path:
        sys.path.insert(0, p)

PAL_POL = ('cert', 'cert_naive', 'run')
PAT_POL = ('adj', 'occ', 'occ_null')


def _tape(n, i):
    import tr_lib as T
    for j, a, s, X in T.tapes(n, i + 1):
        if j == i:
            return a, X
    raise IndexError(i)


def work(job):
    rule, n, policy, i = job
    a, X = _tape(n, i)
    if rule == 'pal':
        import mk_pal
        r = mk_pal.run_base(a) if policy == 'base' else mk_pal.run_pal(a, policy)
    else:
        import mk_pat
        r = mk_pat.run_base(X) if policy == 'base' else mk_pat.run_pat(X, policy)
    r['i'] = i
    return r


def run_cell(rule, n, policy, nt, pool):
    t0 = time.time()
    recs = pool.map(work, [(rule, n, policy, i) for i in range(nt)], chunksize=2)
    recs.sort(key=lambda r: r['i'])
    out = dict(rule=rule, n=n, policy=policy, tapes=nt,
               reads=[float(r['reads']) for r in recs],
               ink=[float(r['ink']) for r in recs],
               delta=[r['delta'] for r in recs],
               end=[int(r['end']) for r in recs],
               steps=[int(r['steps']) for r in recs],
               cuts=[int(r.get('cuts') or 0) for r in recs],
               consults=[int(r['consults']) for r in recs],
               secs=round(time.time() - t0, 1))
    out['reads_mean'] = float(np.mean(out['reads']))
    out['ink_mean'] = float(np.mean(out['ink']))
    return out


def main():
    nt = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    ns = ([int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2
          else [512, 1024, 2048])
    rules = sys.argv[3].split(',') if len(sys.argv) > 3 else ['pal', 'pattern']
    procs = int(sys.argv[4]) if len(sys.argv) > 4 else min(12, os.cpu_count())
    cells = []
    with Pool(procs) as pool:
        for rule in rules:
            pols = ('base',) + (PAL_POL if rule == 'pal' else PAT_POL)
            for n in ns:
                base = None
                for pol in pols:
                    c = run_cell(rule, n, pol, nt, pool)
                    if pol == 'base':
                        base = c
                    c['mismatch_delta'] = sum(
                        1 for x, y in zip(c['delta'], base['delta'])
                        if not _same(x, y))
                    c['mismatch_end'] = sum(1 for x, y in zip(c['end'], base['end'])
                                            if x != y)
                    c['mismatch_steps'] = sum(
                        1 for x, y in zip(c['steps'], base['steps']) if x != y)
                    cells.append(c)
                    print('%-8s n=%-5d %-11s reads=%-11.1f ink=%-10.1f '
                          'x%.4f  bad(d/e/s)=%d/%d/%d  [%.0fs]'
                          % (rule, n, pol, c['reads_mean'], c['ink_mean'],
                             c['reads_mean'] / max(base['reads_mean'], 1e-9),
                             c['mismatch_delta'], c['mismatch_end'],
                             c['mismatch_steps'], c['secs']), flush=True)
    doc = dict(protocol=dict(
        tapes='tr_lib.T.tapes(n, %d): rng default_rng(5000+n+100003 i)' % nt,
        nt=nt, ns=ns,
        plane='one writable/overwritable bit per cell, read with the data for '
              '1 read, written for 1 ink, removed by a trim',
        base_pal='palmach2.walk_mach2 (charge: 1 read per comparison)',
        base_pattern="pcr_lib.Runner.run('pattern') (charge: INC, 1 read per "
                     "cell taken in or given up + re-span per cut)"),
        cells=cells)
    with open(os.path.join(HERE, 'mk_cells.json'), 'w') as f:
        json.dump(doc, f)
    print('wrote mk_cells.json')


def _same(x, y):
    if x is None or y is None:
        return x is None and y is None
    return abs(float(x) - float(y)) <= 1e-12


if __name__ == '__main__':
    main()
