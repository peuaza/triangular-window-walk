# -*- coding: utf-8 -*-
"""el_summary -- the tables of EN_LINEAR.md, from the el_<rule>_<n>.json cells.

Prints
  S1  the machine table: per rule and n, e_old/n (pcr), the best fixed W2, its
      e_new/n, the exactness count on the meeting tapes, the halting bonus.
  S2  the prediction audit: per-epoch mismatch rate above the best W2, split
      into wrong cut / invented terminal, and the tapes affected.
  S3  the terminal M0 quantiles by outcome class.
  S4  the exponent of e_new/n over the four sizes (is it flat?).
  S5  the bounded-check column: the share of bad epochs a prefix-K step check
      would have caught (epoch length <= K).
Writes el_summary.json with every row.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import el_cell as C
import el_lib as E

NS = E.NS


def load(rule, n):
    p = os.path.join(HERE, 'el_%s_%d.json' % (C.safe(rule), n))
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        return json.load(fh)


def q(v, p):
    return float(np.percentile(v, p)) if len(v) else float('nan')


def cell_row(d):
    n = d['n']
    meet = [t for t in d['taps'] if t['end'] == 0]
    other = [t for t in d['taps'] if t['end'] != 0]
    row = dict(rule=d['rule'], n=n, mode=d['mode'], nmeet=len(meet),
               nhalt=len(other), ndraw=d['ndraw'], secs=d['secs'],
               eold=float(np.mean([t['inc'] for t in meet]) / n),
               mean_epoch_len=d['mean_epoch_len'],
               term_meet=[q([t['termM0'] for t in meet], 50),
                          q([t['termM0'] for t in meet], 90),
                          q([t['termM0'] for t in meet], 99),
                          max([t['termM0'] for t in meet]) if meet else -1],
               term_halt=[q([t['termM0'] for t in other], 50),
                          q([t['termM0'] for t in other], 90),
                          max([t['termM0'] for t in other]) if other else -1],
               by_w2={})
    for W in d['w2s']:
        w = str(W)
        em = sum(t['mach'][w]['exact'] for t in meet)
        eh = sum(t['mach'][w]['exact'] for t in other)
        rn = float(np.mean([t['mach'][w]['reads'] for t in meet]) / n)
        bad_above = [b for t in meet for b in t['bad'] if b[0] > W]
        tapes_bad = sum(1 for t in meet if any(b[0] > W for b in t['bad']))
        ncut = sum(1 for b in bad_above if b[2] in (1, 2, 3))
        nterm = sum(1 for b in bad_above if b[2] in (0, 4, 5))
        ep_above = sum(1 for t in meet for b in [0] if True)  # placeholder
        row['by_w2'][w] = dict(exact_meet=int(em), exact_halt=int(eh),
                               e=rn, bad_above=len(bad_above),
                               bad_cut=ncut, bad_term=nterm,
                               tapes_bad=int(tapes_bad),
                               nscan=int(np.mean([t['mach'][w]['nscan'] for t in meet])),
                               eg=float(np.mean([t['mach'][w]['eg'] for t in meet]) / n))
    # epochs above a given W2 (denominator for the per-epoch rate)
    return row


def best_w2(d):
    """smallest W2 that is exact on every meeting tape; else the max-exact one."""
    meet = [t for t in d['taps'] if t['end'] == 0]
    best = None
    for W in d['w2s']:
        w = str(W)
        em = sum(t['mach'][w]['exact'] for t in meet)
        if em == len(meet):
            return W
        if best is None or em > best[1]:
            best = (W, em)
    return best[0]


def main():
    rows = []
    miss = []
    for rule in E.RULES:
        for n in NS:
            d = load(rule, n)
            if d is None:
                miss.append((rule, n))
                continue
            r = cell_row(d)
            r['best_w2'] = best_w2(d)
            bw = r['by_w2'][str(r['best_w2'])]
            r['e_new'] = bw['e']
            r['exact_meet'] = bw['exact_meet']
            r['exact_halt'] = bw['exact_halt']
            # per-epoch audit at the best W2
            meet = [t for t in d['taps'] if t['end'] == 0]
            ep_tot = 0
            for t in meet:
                ep_tot += t['nep'] - 1      # the last epoch is the terminal
            bad = [b for t in meet for b in t['bad'] if b[0] > r['best_w2']]
            r['ep_tot'] = ep_tot
            r['rate'] = (len(bad) / ep_tot) if ep_tot else 0.0
            r['rate_n'] = r['rate'] * n
            r['catch4'] = int(sum(1 for b in bad if b[5] <= 4))
            r['catch16'] = int(sum(1 for b in bad if b[5] <= 16))
            rows.append(r)

    print('\nS1  machine table (meeting tapes)')
    print('%-9s %5s %6s %8s %8s %9s %8s %8s' % (
        'rule', 'n', 'W2', 'e_old/n', 'e_new/n', 'exactMe', 'exactH', 'mode'))
    for r in rows:
        print('%-9s %5d %6d %8.1f %8.2f %5d/%-3d %5d/%-3d %8s' % (
            r['rule'], r['n'], r['best_w2'], r['eold'], r['e_new'],
            r['exact_meet'], r['nmeet'], r['exact_halt'], r['nhalt'], r['mode']))

    print('\nS2  prediction audit at the best W2 (interior epochs, meeting tapes)')
    print('%-9s %5s %8s %8s %6s %7s %7s %9s' % (
        'rule', 'n', 'epochs', 'bad', 'rate', 'rate*n', 'tapes', 'cut/term'))
    for r in rows:
        bw = None
        d = load(r['rule'], r['n'])
        for W in d['w2s']:
            if W == r['best_w2']:
                bw = r['by_w2'][str(W)]
        print('%-9s %5d %8d %8d %8.4f %7.3f %5d/%-3d %4d/%-4d' % (
            r['rule'], r['n'], r['ep_tot'], bw['bad_above'], r['rate'],
            r['rate_n'], bw['tapes_bad'], r['nmeet'], bw['bad_cut'],
            bw['bad_term']))

    print('\nS3  terminal M0 by class: meeting (med/p90/p99/max) / halting (med/p90/max)')
    for r in rows:
        print('%-9s %5d  meet %5.0f %5.0f %5.0f %5d   halt %5.0f %5.0f %5d' % (
            r['rule'], r['n'], r['term_meet'][0], r['term_meet'][1],
            r['term_meet'][2], r['term_meet'][3], r['term_halt'][0],
            r['term_halt'][1], r['term_halt'][2]))

    print('\nS4  exponent of e_new/n and of e_old/n over 512..4096')
    for rule in E.RULES:
        rs = [r for r in rows if r['rule'] == rule]
        if len(rs) < 2:
            continue
        ns = np.array([r['n'] for r in rs], float)
        en = np.array([r['e_new'] for r in rs], float)
        eo = np.array([r['eold'] for r in rs], float)
        ex = float(np.polyfit(np.log(ns), np.log(en), 1)[0])
        xo = float(np.polyfit(np.log(ns), np.log(eo), 1)[0])
        ok = all(r['exact_meet'] == r['nmeet'] for r in rs)
        print('%-9s  e_new/n %s  slope %+.3f   e_old/n %s  slope %+.3f  '
              'exact-at-best-W2 %s' % (
                  rule, ' '.join('%.1f' % v for v in en), ex,
                  ' '.join('%.0f' % v for v in eo), xo,
                  'yes' if ok else 'NO'))

    print('\nS5  bounded check: bad epochs with epoch length <= 4 / <= 16')
    for r in rows:
        if r.get('catch4') or r.get('catch16'):
            print('%-9s %5d  bad %5d  catch<=4 %4d  catch<=16 %4d'
                  % (r['rule'], r['n'],
                     r['by_w2'][str(r['best_w2'])]['bad_above'], r['catch4'],
                     r['catch16']))
    if miss:
        print('\nMISSING cells:', miss)
    with open(os.path.join(HERE, 'el_summary.json'), 'w') as fh:
        json.dump(dict(rows=rows, missing=miss), fh, indent=1)
    print('\nwrote el_summary.json  (%d rows, %d missing)' % (len(rows), len(miss)))


if __name__ == '__main__':
    main()
