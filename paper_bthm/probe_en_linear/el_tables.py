# -*- coding: utf-8 -*-
"""el_md -- emit the report's tables as markdown, straight from el_report.json.

python el_md.py > _tables.md
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, 'el_report.json')))
EOLD = json.load(open(os.path.join(HERE, 'el_eold.json')))
NS = (512, 1024, 2048, 4096)
W2S = (8, 16, 32, 64, 128, 256)

NAME = {'dense2': 'dense2', 'peak': 'peak', 'trip': 'trip', 'oddrun': 'oddrun',
        'canonical': 'canonical', 'LAS': 'LAS', 'KCHG': 'KCHG', 'MAXP': 'MAXP',
        'AMAX': 'AMAX', 'AMIN': 'AMIN', 'BAL>m/2': 'BAL>m/2',
        'BAL>m/4': 'BAL>m/4', 'NBLK=1': 'NBLK=1', 'NBLK>2': 'NBLK>2',
        'NBLK>3': 'NBLK>3'}
DEV = {'dense2': 'count of 2s', 'peak': '3-gram peak count',
       'trip': '3-gram triple count', 'oddrun': 'odd-run count + 2 boundary parities',
       'canonical': 'S, Q, Rr', 'LAS': 'alternation count + first/last nonzero',
       'KCHG': 'two sign-change counts',
       'MAXP': 'range max of the prefix walk', 'AMAX': 'range max',
       'AMIN': 'range min',
       'BAL>m/2': 'last occurrence of the level', 'BAL>m/4': 'last occurrence of the level',
       'NBLK=1': 'count of visits to the level', 'NBLK>2': 'count of visits to the level',
       'NBLK>3': 'count of visits to the level'}
CLS = {'dense2': 'slide', 'peak': 'slide', 'trip': 'slide', 'oddrun': 'slide',
       'canonical': 'slide', 'LAS': 'slide', 'KCHG': 'slide',
       'MAXP': 'lazy', 'AMAX': 'lazy', 'AMIN': 'lazy',
       'BAL>m/2': 'scan', 'BAL>m/4': 'scan', 'NBLK=1': 'scan',
       'NBLK>2': 'scan', 'NBLK>3': 'scan'}


def f(x, d=1):
    if x is None or (isinstance(x, float) and x != x):
        return '—'
    return ('%.' + str(d) + 'f') % x


print('### M1 — the device per rule\n')
print('| rule | b1 is a query on | class | reads per query at the region |')
print('|---|---|---|---|')
for rule in R:
    cls = CLS[rule]
    cost = {'slide': '0 (one counter)', 'lazy': '0, or the region width when the argmax dies',
            'scan': 'the region width, every query'}[cls]
    print('| `%s` | %s | %s | %s |' % (NAME[rule], DEV[rule], cls, cost))
print()

print('### M2 — e_old/n (the record machine, mean over meeting tapes) vs the '
      'exact-on-every-tape counter machine\n')
print('| rule | machine | ' + ' | '.join('n=%d' % n for n in NS) + ' |')
print('|---|---|' + '---|' * len(NS))
for rule in R:
    v = R[rule]
    ns = v['sizes']
    eo = []
    for n in NS:
        k = '%s|%d' % (rule, n)
        eo.append(f(EOLD.get(k, {}).get('e_meet'), 1) if k in EOLD else '—')
    mach = EOLD.get('%s|%d' % (rule, ns[0]), {}).get('machine', '')
    print('| `%s` | %s | %s |' % (NAME[rule], mach, ' | '.join(eo)))
print()

print('### M3 — the W2 x n grid: exact meeting tapes and e/n\n')
for rule in R:
    v = R[rule]
    ns = v['sizes']
    print('`%s` (%s), e_old/n = %s\n' % (
        NAME[rule], v['mode'],
        ' / '.join(f(EOLD.get('%s|%d' % (rule, n), {}).get('e_meet'), 0)
                   for n in ns)))
    print('| n | meet | ' + ' | '.join('W2=%d' % w for w in W2S) + ' |')
    print('|---|' + '---|' * (len(W2S) + 1))
    for n in ns:
        c = v['cells'][str(n)]
        cells = []
        for w in W2S:
            d = c['w2'][str(w)]
            cells.append('%d/%d · %.1f' % (d['exact_meet'], c['nmeet'], d['e']))
        print('| %d | %d | %s |' % (n, c['nmeet'], ' | '.join(cells)))
    print()

print('### M4 — the failure anatomy at the best fixed W2\n')
print('| rule | W2 | exact at every n | interior cut errors | invented terminals '
      '| b2=0 errors | median M0 of a failure (by n) |')
print('|---|---|---|---|---|---|---|')
for rule in R:
    v = R[rule]
    w = str(v['best_w2'])
    cut = sum(v['cells'][str(n)]['w2'][w]['bad_cut'] for n in v['sizes'])
    term = sum(v['cells'][str(n)]['w2'][w]['bad_term'] for n in v['sizes'])
    b2z = sum(v['cells'][str(n)]['w2'][w]['bad_b2z'] for n in v['sizes'])
    m0 = [v['cells'][str(n)]['w2'][w]['bad_m0'][0] for n in v['sizes']]
    print('| `%s` | %d | %s | %d | %d | %d | %s |' % (
        NAME[rule], v['best_w2'], 'yes' if v['exact_all'] else 'no', cut, term,
        b2z, ' / '.join(f(x, 0) for x in m0)))
print()

print('### M5 — meeting terminal M0 (med / p90 / p99 / max)\n')
print('| rule | ' + ' | '.join('n=%d' % n for n in NS) + ' |')
print('|---|' + '---|' * len(NS))
for rule in R:
    v = R[rule]
    cells = []
    for n in v['sizes']:
        t = v['cells'][str(n)]['term_meet']
        cells.append('%d / %d / %d / %d' % tuple(int(x) for x in t))
    print('| `%s` | %s |' % (NAME[rule], ' | '.join(cells)))
print()

print('### M6 — halting terminal M0 (med / p90 / max): why no fixed W2 is '
      'exact on every tape\n')
print('| rule | ' + ' | '.join('n=%d' % n for n in NS) + ' |')
print('|---|' + '---|' * len(NS))
for rule in R:
    v = R[rule]
    cells = []
    for n in v['sizes']:
        t = v['cells'][str(n)]['term_halt']
        cells.append('%d / %d / %d' % tuple(int(x) for x in t))
    print('| `%s` | %s |' % (NAME[rule], ' | '.join(cells)))
print()

print('### M7 — the exact-on-every-tape counter machine (n + steps + cuts), '
      'sliding rules only\n')
print('| rule | ' + ' | '.join('n=%d' % n for n in NS) + ' |')
print('|---|' + '---|' * len(NS))
for rule in R:
    if CLS[rule] != 'slide':
        continue
    v = R[rule]
    vals = [v['cells'][str(n)]['slide'] for n in v['sizes']]
    print('| `%s` | %s |' % (NAME[rule], ' | '.join(f(x, 0) for x in vals)))
print()
