# -*- coding: utf-8 -*-
"""el_tab -- the tables of EN_LINEAR.md as markdown, from the cell JSONs.

python el_tab.py > _tables.md
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import el_cell as C
import el_lib as E

NS = (512, 1024, 2048, 4096)
W2S = (8, 16, 32, 64, 128, 256)
EOLD = {}
_p = os.path.join(HERE, 'el_eold.json')
if os.path.exists(_p):
    EOLD = json.load(open(_p))

DEV = {'dense2': 'count of 2s', 'peak': '3-gram peak count',
       'trip': '3-gram triple count',
       'oddrun': 'odd-run count + two boundary parities',
       'canonical': 'S, Q, Rr and the two end values',
       'LAS': 'alternation count + first/last nonzero',
       'KCHG': 'two sign-change counts',
       'MAXP': 'range max of the prefix walk', 'AMAX': 'range max',
       'AMIN': 'range min',
       'BAL>m/2': 'last occurrence of the level',
       'BAL>m/4': 'last occurrence of the level',
       'NBLK=1': 'count of visits to the level',
       'NBLK>2': 'count of visits to the level',
       'NBLK>3': 'count of visits to the level'}
CLS = {'dense2': 'slide', 'peak': 'slide', 'trip': 'slide', 'oddrun': 'slide',
       'canonical': 'slide', 'LAS': 'slide', 'KCHG': 'slide',
       'MAXP': 'lazy', 'AMAX': 'lazy', 'AMIN': 'lazy',
       'BAL>m/2': 'scan', 'BAL>m/4': 'scan', 'NBLK=1': 'scan',
       'NBLK>2': 'scan', 'NBLK>3': 'scan'}


def cells(rule):
    out = {}
    for n in NS:
        p = os.path.join(HERE, 'el_%s_%d.json' % (C.safe(rule), n))
        if os.path.exists(p):
            out[n] = json.load(open(p))
    return out


def fmt(x, d=1):
    return ('%.' + str(d) + 'f') % x if x is not None else '—'


def main():
    print('#### M1 — the device per rule\n')
    print('| rule | `b1` at the region is | class | cost of one query |')
    print('|---|---|---|---|')
    for rule in E.RULES:
        c = CLS[rule]
        cost = {'slide': '0 (one counter, 1 read per cell given up)',
                'lazy': '0, or the region width when the cached argmax dies',
                'scan': 'the region width, every query'}[c]
        print('| `%s` | %s | %s | %s |' % (rule, DEV[rule], c, cost))
    print()

    print('#### M2 — `e` on the meeting tapes: the record machine, and the '
          'exact-on-every-tape counter machine\n')
    print('| rule | record machine | ' + ' | '.join('%d' % n for n in NS) +
          ' | `n+steps+cuts` |')
    print('|---|---|' + '---|' * len(NS) + '---|')
    for rule in E.RULES:
        cs = cells(rule)
        if not cs:
            continue
        eo = []
        for n in NS:
            k = '%s|%d' % (rule, n)
            eo.append(fmt(EOLD[k]['e_meet'], 0) if k in EOLD else '—')
        mach = EOLD.get('%s|%d' % (rule, min(cs)), {}).get('machine', '')
        mach = mach.replace('pcr_lib.Runner', '`pcr` INC').replace(
            'hmach3.H3 B=64', '`hmach3`').replace(
            'tr_bal4 (BAL device, thr=4)', '`tr_bal4`')
        sl = []
        for n in NS:
            if n in cs:
                m = [t for t in cs[n]['taps'] if t['end'] == 0]
                sl.append(fmt(np.mean([(n + t['steps'] + t['cuts']) / n
                                       for t in m]), 0))
            else:
                sl.append('—')
        print('| `%s` | %s | %s | %s |' % (rule, mach, ' | '.join(eo),
                                           ' | '.join(sl)))
    print()

    print('#### M3 — the `W2` x `n` grid: `exact/nmeet · e/n`\n')
    for rule in E.RULES:
        cs = cells(rule)
        if not cs:
            continue
        print('`%s` (%s)\n' % (rule, cs[min(cs)]['mode']))
        print('| n | meet | ' + ' | '.join('W2=%d' % w for w in W2S) + ' |')
        print('|---|---|' + '---|' * len(W2S))
        for n in sorted(cs):
            c = cs[n]
            meet = [t for t in c['taps'] if t['end'] == 0]
            row = []
            for w in W2S:
                ex = sum(t['mach'][str(w)]['exact'] for t in meet)
                e = np.mean([t['mach'][str(w)]['reads'] for t in meet]) / n
                row.append('%d/%d · %.1f' % (ex, len(meet), e))
            print('| %d | %d | %s |' % (n, len(meet), ' | '.join(row)))
        print()

    print('#### M4 — failure anatomy at `W2 = 256` (meeting tapes)\n')
    print('| rule | exact by n | mispredicted epochs (cut / terminal) | '
          '`b2=0` errors | at `M0 ≥ 256` | epoch len ≤ 4 | tapes with ≥1 |')
    print('|---|---|---|---|---|---|---|')
    for rule in E.RULES:
        cs = cells(rule)
        if not cs:
            continue
        ex, cut, term, b2z, hi, l4, tb = [], 0, 0, 0, 0, 0, 0
        for n in sorted(cs):
            c = cs[n]
            meet = [t for t in c['taps'] if t['end'] == 0]
            ex.append('%d/%d' % (sum(t['mach']['256']['exact'] for t in meet),
                                 len(meet)))
            bad = [b for t in meet for b in t['bad'] if b[0] > 256]
            cut += sum(1 for b in bad if b[2] in (1, 2, 3))
            term += sum(1 for b in bad if b[2] in (0, 4, 5))
            b2z += sum(1 for b in bad if b[1] == 3)
            hi += len(bad)
            l4 += sum(1 for b in bad if b[5] <= 4)
            tb += sum(1 for t in meet if any(b[0] > 256 for b in t['bad']))
        print('| `%s` | %s | %d / %d | %d | %d | %d | %d |'
              % (rule, ' '.join(ex), cut, term, b2z, hi, l4, tb))
    print()

    print('#### M5 — meeting terminal `M0` (med / p90 / p99 / max)\n')
    print('| rule | ' + ' | '.join('%d' % n for n in NS) + ' |')
    print('|---|' + '---|' * len(NS))
    for rule in E.RULES:
        cs = cells(rule)
        if not cs:
            continue
        row = []
        for n in NS:
            if n in cs:
                m = [t['termM0'] for t in cs[n]['taps'] if t['end'] == 0]
                row.append('%d / %d / %d / %d' % (
                    np.median(m), np.percentile(m, 90), np.percentile(m, 99),
                    max(m)))
            else:
                row.append('—')
        print('| `%s` | %s |' % (rule, ' | '.join(row)))
    print()

    print('#### M6 — halting terminal `M0` (med / p90 / max)\n')
    print('| rule | ' + ' | '.join('%d' % n for n in NS) + ' |')
    print('|---|' + '---|' * len(NS))
    for rule in E.RULES:
        cs = cells(rule)
        if not cs:
            continue
        row = []
        for n in NS:
            if n in cs:
                m = [t['termM0'] for t in cs[n]['taps'] if t['end'] != 0]
                row.append('%d / %d / %d' % (np.median(m),
                                             np.percentile(m, 90), max(m)))
            else:
                row.append('—')
        print('| `%s` | %s |' % (rule, ' | '.join(row)))
    print()

    print('#### M7 — the verdict table\n')
    print('| rule | `W2` | exact on the meeting tapes | `e_new/n` by n | '
          'slope of `e/n` | tapes with a mispredicted epoch | verdict |')
    print('|---|---|---|---|---|---|---|')
    for rule in E.RULES:
        cs = cells(rule)
        if not cs:
            continue
        ns = sorted(cs)
        best = None
        for w in W2S:
            if all(sum(t['mach'][str(w)]['exact'] for t in cs[n]['taps']
                       if t['end'] == 0) ==
                   sum(1 for t in cs[n]['taps'] if t['end'] == 0) for n in ns):
                best = w
                break
        exact_all = best is not None
        if best is None:
            tot = {w: sum(sum(t['mach'][str(w)]['exact'] for t in cs[n]['taps']
                              if t['end'] == 0) for n in ns) for w in W2S}
            best = max(tot, key=tot.get)
        ex, en, tb = [], [], 0
        for n in ns:
            meet = [t for t in cs[n]['taps'] if t['end'] == 0]
            ex.append('%d/%d' % (sum(t['mach'][str(best)]['exact'] for t in meet),
                                 len(meet)))
            en.append(np.mean([t['mach'][str(best)]['reads'] for t in meet]) / n)
            tb += sum(1 for t in meet if any(b[0] > best for b in t['bad']))
        slope = float(np.polyfit(np.log(np.array(ns, float)),
                                 np.log(np.array(en, float)), 1)[0])
        if exact_all and slope < 0.35:
            v = '**Θ(n)**'
        elif exact_all:
            v = 'exact but `e/n` rises: Θ(n²)'
        else:
            v = 'not exact: still Θ(n²)'
        print('| `%s` | %d | %s | %s | %+.2f | %d | %s |'
              % (rule, best, ' '.join(ex), ' '.join('%.1f' % x for x in en),
                 slope, tb, v))
    print()

    print('#### M8 — `e/n` if the endgame cost 1 read per step (a band device), '
          'estimate not measurement\n')
    print('| rule | ' + ' | '.join('%d' % n for n in NS) + ' |')
    print('|---|' + '---|' * len(NS))
    for rule in E.RULES:
        cs = cells(rule)
        if not cs:
            continue
        row = []
        for n in NS:
            if n in cs:
                c = cs[n]
                meet = [t for t in c['taps'] if t['end'] == 0]
                eg = np.mean([t['mach']['256']['eg'] for t in meet]) / n
                cu = np.mean([t['cuts'] for t in meet]) / n
                row.append('%.1f' % (1 + cu + eg))
            else:
                row.append('—')
        print('| `%s` | %s |' % (rule, ' | '.join(row)))
    print()


if __name__ == '__main__':
    main()
