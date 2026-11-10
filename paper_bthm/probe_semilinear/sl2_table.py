# -*- coding: utf-8 -*-
"""sl2_table -- render the numbers from sl2_measure.json.  Nothing here is typed by hand."""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'sl2_measure.json')))
cells = D['cells']
ex = D['ex']
ns = [c['n'] for c in cells]


def fmt(v, k=3):
    return ('%.' + str(k) + 'f') % v


print('| $n$ | tapes | meet | $e/n$ | $g/n^2$ | $G/n$ | $C_{norm}/n$ | '
      '$C_{def}/n^2$ | $q$ | $1/n$ | $q\\cdot n$ |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for c in cells:
    print('| %d | %d | %d/%d | %s | %s | %s | %s | %s | %s | %s | %s |'
          % (c['n'], c['nt'], c['n_meet'], c['nt'],
             fmt(c['e_over_n'], 3), fmt(c['g_over_n2'], 4), fmt(c['G_over_n'], 3),
             fmt(c['c_normal'] / c['n'], 3), fmt(c['c_defect'] / c['n'] ** 2, 4),
             ('%.2e' % c['q']), ('%.2e' % (1.0 / c['n'])), fmt(c['q'] * c['n'], 3)))
print()
print('| exponent | value |')
print('|---|---|')
for k in ('e_mix', 'g_mix', 'G_mix', 'c_normal', 'c_defect', 'q', 'e_raw'):
    print('| `%s` | %.3f |' % (k, ex[k]))
print()
print('per-cell detail')
for c in cells:
    d = np.array(c['delta'])
    print('n=%-5d 200-tape: e/n %.3f  g/n %.3f  G/n %.3f  steps/n %.4f  '
          'delta: distinct %d  sd %.1f cells  q25 %.3f med %.3f q75 %.3f'
          % (c['n'], c['e_raw'] / c['n'], c['g_raw'] / c['n'], c['G_raw'] / c['n'],
             np.mean(c['steps']) / c['n'], len(set(np.round(d, 9))),
             float(np.std(d)) * c['n'], np.percentile(d, 25), np.median(d),
             np.percentile(d, 75)))
    w = c['witness']
    print('        4000-tape: e/n %.3f  max/n^2 %.4f  defects %d   witnesses: '
          'inc/n^2 %s  msum/n^3 %s'
          % (c['e_big'] / c['n'], c['g_big'] / c['n'] ** 2, c['n_defect_big'],
             '/'.join('%.4f' % (r['inc'] / c['n'] ** 2) for r in w),
             '/'.join('%.4f' % (r['msum'] / c['n'] ** 3) for r in w)))
    ms = np.array(c['msum'], float)
    print('        SCAN: normal msum/n^2 %.4f (constant by geometry)  '
          'defect msum/n^3 %.4f  => e_scan/n^2 %.3f  g_scan/n^3 %.4f'
          % (np.mean(ms) / c['n'] ** 2,
             np.mean([r['msum'] for r in w]) / c['n'] ** 3,
             ((1 - c['q']) * np.mean(ms) + c['q'] * np.mean([r['msum'] for r in w]))
             / c['n'] ** 2,
             max(max(r['msum'] for r in w), ms.max()) / c['n'] ** 3))
print()
print('defect fraction vs 1/n')
for c in cells:
    print('  n=%-5d q = %.6f +- %.6f   (%.3f / n)   hits %d of %d'
          % (c['n'], c['q'], c['q_ci'], c['q'] * c['n'], c['q_hit'], c['q_n']))
