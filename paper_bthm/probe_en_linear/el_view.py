# -*- coding: utf-8 -*-
"""el_view -- print one cell's table (debug view)."""
import json
import sys

import numpy as np

for f in sys.argv[1:]:
    d = json.load(open(f))
    n = d['n']
    meet = [t for t in d['taps'] if t['end'] == 0]
    other = [t for t in d['taps'] if t['end'] != 0]
    print('==', d['rule'], n, 'mode', d['mode'], 'meet', len(meet),
          'halt', len(other), 'draws', d['ndraw'],
          'e_old/n(meet) = %.1f' % (np.mean([t['inc'] for t in meet]) / n))
    print('   %6s %11s %10s %11s %10s %8s' % ('W2', 'exactMe', 'e/n Me',
                                              'exactH', 'e/n H', 'egSteps'))
    for W in d['w2s']:
        w = str(W)
        em = sum(t['mach'][w]['exact'] for t in meet)
        eh = sum(t['mach'][w]['exact'] for t in other)
        rm = np.mean([t['mach'][w]['reads'] for t in meet]) / n
        rh = np.mean([t['mach'][w]['reads'] for t in other]) / n
        eg = np.mean([t['mach'][w]['eg'] for t in meet]) / n
        print('   %6d %6d/%-4d %10.2f %6d/%-4d %10.2f %8.2f'
              % (W, em, len(meet), rm, eh, len(other), rh, eg))
    nb = sum(len(t['bad']) for t in meet)
    print('   interior mismatches on meeting tapes: %d over %d epochs (%.2f/tape)'
          % (nb, sum(t['nep'] for t in meet), nb / max(1, len(meet))))
    if meet:
        tm = [t['termM0'] for t in meet]
        print('   meeting terminal M0 med/p90/p99/max = %d/%d/%d/%d'
              % (np.median(tm), np.percentile(tm, 90), np.percentile(tm, 99),
                 max(tm)))
    if other:
        tm = [t['termM0'] for t in other]
        print('   halting  terminal M0 med/p90/max = %d/%d/%d'
              % (np.median(tm), np.percentile(tm, 90), max(tm)))
    print('   mean epoch len %.1f' % d['mean_epoch_len'],
          'term-len hist', dict(zip(d['binnames'], d['ep_len_hist'])))
