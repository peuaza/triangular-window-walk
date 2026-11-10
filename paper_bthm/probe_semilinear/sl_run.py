# -*- coding: utf-8 -*-
"""sl_run.py -- measure the threshold-run machine on the ensemble.

usage: sl_run.py <n list> <tapes> <seed list> [c]
writes _data/sl_n<n>_s<seed>.jsonl  (one row per tape)
prints the e/g table.
"""
import sys
import json
import math
import os
import numpy as np

sys.path.insert(0, '../probe_ologn2')
from rs_common import sample_tapes  # noqa: E402
from sl_machine import solve, t_of  # noqa: E402


def summarize(n, t, rows):
    rd = np.array([r['reads'] for r in rows], float)
    st = np.array([r['steps'] for r in rows], float)
    rho = np.array([r['rho'] for r in rows])
    exc = rho > t
    d = [r['delta'] for r in rows if r['delta'] is not None]
    cnt = {}
    for x in d:
        k = round(x, 12)
        cnt[k] = cnt.get(k, 0) + 1
    H = -sum((v / len(d)) * math.log2(v / len(d)) for v in cnt.values()) if d else 0.0
    return dict(
        n=n, t=t, tapes=len(rows), rare=int(exc.sum()),
        mean=float(rd.mean()), e_over_n=float(rd.mean() / n),
        max=float(rd.max()), g_over_n=float(rd.max() / n),
        q99=float(np.quantile(rd, 0.99)), q999=float(np.quantile(rd, 0.999)),
        med=float(np.median(rd)),
        exc_mean=float(rd[exc].mean()) if exc.any() else None,
        exc_over_n2=float(rd[exc].mean() / n ** 2) if exc.any() else None,
        exc_steps=float(st[exc].mean()) if exc.any() else None,
        exc_steps_max=float(st[exc].max()) if exc.any() else None,
        exc_steps_over_n=float(st[exc].max() / n) if exc.any() else None,
        steps_mean=float(st.mean()), steps_over_n=float(st.mean() / n),
        steps_max=int(st.max()),
        p_none=float(sum(1 for r in rows if r['delta'] is None) / len(rows)),
        distinct=len(cnt), H=H,
        modal=(max(cnt.values()) / len(d)) if d else 0.0,
        maxcache=int(max(r['maxcache'] for r in rows)),
        modes={k: sum(1 for r in rows if r['mode'] == k) for k in ('shortcut', 'literal')},
        scans_mean=float(np.mean([r['scans'] for r in rows])),
        scan_reads_share=float(np.sum([r['scan_reads'] for r in rows]) / rd.sum()),
    )


def main():
    ns = [int(x) for x in sys.argv[1].split(',')]
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    seeds = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else [2026]
    c = int(sys.argv[4]) if len(sys.argv) > 4 else -1
    store = bool(int(sys.argv[5])) if len(sys.argv) > 5 else False
    tag = ('_st' if store else '') + ('' if c == -1 else '_c%d' % c)
    os.makedirs('_data', exist_ok=True)
    print('%5s %4s %6s %6s %10s %8s %10s %8s %8s %8s %10s %8s' %
          ('n', 't', 'tapes', 'rare', 'mean', 'e/n', 'max', 'g/n', 'q99/n',
           'q999/n', 'exc_mean', 'exc/n2'))
    out = []
    for n in ns:
        t = t_of(n, c)
        for seed in seeds:
            rows = []
            path = '_data/sl_n%d_s%d%s.jsonl' % (n, seed, tag)
            with open(path, 'w') as fh:
                for a in sample_tapes(n, N, seed):
                    r = solve(a, t, store=store)
                    r['delta'] = None if r['delta'] is None else float(r['delta'])
                    r['n'] = n
                    r['seed'] = seed
                    rows.append(r)
                    fh.write(json.dumps({k: r[k] for k in
                                         ('n', 'seed', 'reads', 'steps', 'rho', 'mode',
                                          'delta', 'maxcache', 'scans', 'scan_reads',
                                          'nlong', 'overflow')} |
                                        {'cycle': bool(r.get('cycle', False)),
                                         'cap': bool(r.get('cap', False))}) + '\n')
            s = summarize(n, t, rows)
            s['seed'] = seed
            out.append(s)
            print('%5d %4d %6d %6d %10.1f %8.2f %10.1f %8.1f %8.1f %8.1f %10.1f %8.3f' %
                  (n, t, s['tapes'], s['rare'], s['mean'], s['e_over_n'], s['max'],
                   s['g_over_n'], s['q99'] / n, s['q999'] / n,
                   s['exc_mean'] if s['exc_mean'] else float('nan'),
                   s['exc_over_n2'] if s['exc_over_n2'] else float('nan')), flush=True)
            print('        steps %.1f (%.2f n, max %d)  exc_steps mean %.1f max %.1f (%.1f n)  '
                  'P(none)=%.3f  distinct=%d H=%.2f modal=%.3f  maxcache=%d  scans=%.1f  '
                  'scan_share=%.3f  modes=%s'
                  % (s['steps_mean'], s['steps_over_n'], s['steps_max'],
                     s['exc_steps'] if s['exc_steps'] is not None else float('nan'),
                     s['exc_steps_max'] if s['exc_steps_max'] is not None else float('nan'),
                     s['exc_steps_over_n'] if s['exc_steps_over_n'] is not None else float('nan'),
                     s['p_none'], s['distinct'], s['H'], s['modal'], s['maxcache'],
                     s['scans_mean'], s['scan_reads_share'], s['modes']), flush=True)
    with open('_data/sl_summary%s.json' % tag, 'w') as fh:
        json.dump(out, fh, indent=1, default=float)


if __name__ == '__main__':
    main()
