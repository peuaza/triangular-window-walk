# -*- coding: utf-8 -*-
"""sl_analyze.py -- final tables and fits from the _data files."""
import json
import glob
import math
import numpy as np


def load(pattern):
    rows = []
    for p in sorted(glob.glob(pattern)):
        rows += [json.loads(l) for l in open(p)]
    return rows


def fit(xs, ys):
    xs = np.log(np.array(xs, float))
    ys = np.log(np.array(ys, float))
    return float(np.polyfit(xs, ys, 1)[0])


def main():
    for tag, label in (('_st', 'store machine'), ('_st_c-2', 'store, t=c-2')):
        print('\n==', label)
        print('%5s %4s %6s %6s %8s %8s %8s %10s %10s %10s' %
              ('n', 't', 'tapes', 'rare', 'p*n', 'mean', 'e/n', 'max', 'g/n',
               'exc_mean/n'))
        for n in (128, 256, 512, 1024, 2048, 4096):
            files = sorted(glob.glob('_data/sl_n%d_s*%s.jsonl' % (n, tag)))
            files = [f for f in files if '_c-2' not in f or tag == '_st_c-2']
            rows = []
            for f in files:
                rows += [json.loads(l) for l in open(f)]
            if not rows:
                continue
            t = rows[0]['n'] and None
            rd = np.array([r['reads'] for r in rows], float)
            rho = np.array([r['rho'] for r in rows])
            t = int(np.floor(np.log2(n))) - (2 if tag == '_st_c-2' else 1)
            exc = rho > t
            print('%5d %4d %6d %6d %8.2f %8.1f %8.2f %10.0f %8.1f %10.1f' %
                  (n, t, len(rows), exc.sum(), exc.sum() / len(rows) * n,
                   rd.mean(), rd.mean() / n, rd.max(), rd.max() / n,
                   rd[exc].mean() / n if exc.any() else float('nan')))
        ns, es, gs = [], [], []
        for n in (128, 256, 512, 1024, 2048, 4096):
            rows = []
            for f in sorted(glob.glob('_data/sl_n%d_s*%s.jsonl' % (n, tag))):
                rows += [json.loads(l) for l in open(f)]
            if len(rows) >= 2000:
                ns.append(n)
                es.append(np.mean([r['reads'] for r in rows]) / n)
                gs.append(max(r['reads'] for r in rows) / n)
        if len(ns) >= 3:
            print('  fit e ~ n^%.3f   fit g ~ n^%.3f   (e/n %.2f..%.2f, g/n %.1f..%.1f)'
                  % (fit(ns, [e * n for e, n in zip(es, ns)]),
                     fit(ns, [g * n for g, n in zip(gs, ns)]),
                     min(es), max(es), min(gs), max(gs)))

    print('\n== constructed worst case (store machine)')
    for n, st in ((256, 14542), (512, 52434), (1024, 244927), (2048, 812303)):
        print('  n=%5d steps=%8d (%.1f n, %.3f n^2)' % (n, st, st / n, st / n ** 2))

    print('\n== 2000-tape table, both seeds (store)')
    for n in (128, 256, 512):
        for seed in (2026, 777):
            rows = [json.loads(l) for l in
                    open('_data/sl_n%d_s%d_st.jsonl' % (n, seed))]
            rd = np.array([r['reads'] for r in rows], float)
            rho = np.array([r['rho'] for r in rows])
            t = int(np.floor(np.log2(n))) - 1
            exc = rho > t
            print('  n=%4d seed=%4d rare=%3d mean=%8.1f e/n=%.2f max=%9.0f g/n=%7.1f '
                  'Pnone=%.3f distinct=%3d H=%.2f modal=%.3f'
                  % (n, seed, exc.sum(), rd.mean(), rd.mean() / n, rd.max(),
                     rd.max() / n,
                     sum(1 for r in rows if r['delta'] is None) / len(rows),
                     len({round(r['delta'], 12) for r in rows if r['delta'] is not None}),
                     -sum((v / max(1, sum(1 for r in rows if r['delta'] is not None))) *
                          math.log2(v / max(1, sum(1 for r in rows if r['delta'] is not None)))
                          for v in __import__('collections').Counter(
                              round(r['delta'], 12) for r in rows
                              if r['delta'] is not None).values()),
                     max(__import__('collections').Counter(
                         round(r['delta'], 12) for r in rows
                         if r['delta'] is not None).values()) /
                     max(1, sum(1 for r in rows if r['delta'] is not None))))


if __name__ == '__main__':
    main()
