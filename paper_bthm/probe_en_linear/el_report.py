# -*- coding: utf-8 -*-
"""el_report -- the tables of EN_LINEAR.md from the cell JSONs.

T1  per rule: device class, e_old/n (the record's own machine, mean over meeting
    tapes), e_slide/n (the exact-on-all-tapes counter machine n+steps+cuts),
    and the verdict inputs.
T2  per rule x n: exact/nmeet and e/n at each fixed W2 (the W2 x n grid).
T3  per rule: the smallest fixed W2 that is exact on every meeting tape at
    every size, with its e/n; or the closest grid point.
T4  the failure anatomy at that W2: interior cut errors vs invented terminals,
    the M0 of the first failure, the epoch length of the failures (what a
    bounded prefix check could catch).
T5  terminal M0 by outcome class.
Writes el_report.json.
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


def load(rule, n):
    p = os.path.join(HERE, 'el_%s_%d.json' % (C.safe(rule), n))
    if not os.path.exists(p):
        return None
    return json.load(open(p))


def q(v, p):
    return float(np.percentile(v, p)) if len(v) else float('nan')


def main():
    old = {}
    po = os.path.join(HERE, 'el_eold.json')
    if os.path.exists(po):
        old = json.load(open(po))
    out = {}
    for rule in E.RULES:
        rs = {}
        for n in NS:
            d = load(rule, n)
            if d is None:
                continue
            meet = [t for t in d['taps'] if t['end'] == 0]
            other = [t for t in d['taps'] if t['end'] != 0]
            row = dict(n=n, nmeet=len(meet), nhalt=len(other), ndraw=d['ndraw'],
                       mode=d['mode'], secs=d['secs'],
                       inc=float(np.mean([t['inc'] for t in meet]) / n),
                       steps=float(np.mean([t['steps'] for t in meet]) / n),
                       cuts=float(np.mean([t['cuts'] for t in meet]) / n),
                       slide=float(np.mean([(n + t['steps'] + t['cuts']) / n
                                            for t in meet])),
                       term_meet=[q([t['termM0'] for t in meet], 50),
                                  q([t['termM0'] for t in meet], 90),
                                  q([t['termM0'] for t in meet], 99),
                                  max(t['termM0'] for t in meet)],
                       term_halt=[q([t['termM0'] for t in other], 50),
                                  q([t['termM0'] for t in other], 90),
                                  max(t['termM0'] for t in other)],
                       w2={})
            for W in W2S:
                w = str(W)
                bad = [b for t in meet for b in t['bad'] if b[0] > W]
                tapes_bad = [t for t in meet if any(b[0] > W for b in t['bad'])]
                fb = []
                for t in tapes_bad:
                    b = [x for x in t['bad'] if x[0] > W][0]
                    fb.append((1.0 - b[0] / n, 1 if b[4] else 0))
                row['w2'][w] = dict(
                    exact_meet=int(sum(t['mach'][w]['exact'] for t in meet)),
                    exact_halt=int(sum(t['mach'][w]['exact'] for t in other)),
                    e=float(np.mean([t['mach'][w]['reads'] for t in meet]) / n),
                    eg=float(np.mean([t['mach'][w]['eg'] for t in meet]) / n),
                    nscan=float(np.mean([t['mach'][w]['nscan'] for t in meet])),
                    bad=len(bad),
                    bad_cut=int(sum(1 for b in bad if b[2] in (1, 2, 3))),
                    bad_term=int(sum(1 for b in bad if b[2] in (0, 4, 5))),
                    bad_len4=int(sum(1 for b in bad if b[5] <= 4)),
                    bad_len16=int(sum(1 for b in bad if b[5] <= 16)),
                    bad_b2z=int(sum(1 for b in bad if b[1] == 3)),
                    bad_b2o=int(sum(1 for b in bad if b[1] != 3)),
                    bad_m0=[q([b[0] for b in bad], 50),
                            q([b[0] for b in bad], 90)],
                    tapes_bad=len(tapes_bad),
                    by_bucket={},
                    by_bucket_hi={},
                    first_pos=[q([x[0] for x in fb], 50),
                               float(np.mean([x[1] for x in fb])) if fb else 0.0])
            for W in W2S:
                w = str(W)
                bad = [b for t in meet for b in t['bad'] if b[0] > W]
                d = row['w2'][w]
                for b in bad:
                    k = 'b2=0' if b[1] == 3 else ('b1=0' if b[1] == 1 else 'b1=1')
                    d['by_bucket'][k] = d['by_bucket'].get(k, 0) + 1
                    if b[0] >= 256:
                        d['by_bucket_hi'][k] = d['by_bucket_hi'].get(k, 0) + 1
            rs[n] = row
        if not rs:
            continue
        # smallest W2 exact at every size present
        best = None
        for W in W2S:
            w = str(W)
            if all(rs[n]['w2'][w]['exact_meet'] == rs[n]['nmeet'] for n in rs):
                best = W
                break
        # smallest W2 that is STRUCTURALLY exact: right answer and no
        # mispredicted epoch, i.e. the machine's cut sequence is the walk's
        struct = None
        for W in W2S:
            w = str(W)
            if all(rs[n]['w2'][w]['exact_meet'] == rs[n]['nmeet']
                   and rs[n]['w2'][w]['bad'] == 0 for n in rs):
                struct = W
                break
        # the W2 with the largest total exactness (for the rules with none)
        if best is None:
            tot = {W: sum(rs[n]['w2'][str(W)]['exact_meet'] for n in rs)
                   for W in W2S}
            best = max(tot, key=tot.get)
        out[rule] = dict(sizes=sorted(rs), cells={str(n): rs[n] for n in rs},
                         best_w2=best, struct_w2=struct,
                         exact_all=all(rs[n]['w2'][str(best)]['exact_meet']
                                       == rs[n]['nmeet'] for n in rs),
                         eold=[old.get('%s|%d' % (rule, n), {}).get('e_meet')
                               for n in sorted(rs)],
                         mode=rs[sorted(rs)[0]]['mode'])

    print('T1  e_old/n (record machine, meeting tapes) and the exact-on-all '
          'counter machine')
    print('%-9s %-8s %s | %s | %s' % ('rule', 'mode', 'e_old/n by n',
                                      'e_slide/n by n', 'steps/n by n'))
    for rule, v in out.items():
        ns = v['sizes']
        print('%-9s %-8s %s | %s | %s' % (
            rule, v['mode'],
            ' '.join('%.0f' % v['cells'][str(n)]['inc'] for n in ns),
            ' '.join('%.0f' % v['cells'][str(n)]['slide'] for n in ns),
            ' '.join('%.0f' % v['cells'][str(n)]['steps'] for n in ns)))

    print('\nT3  smallest W2 exact on every meeting tape at every size')
    print('%-9s %5s %8s  %s' % ('rule', 'W2', 'exact', 'e/n by n'))
    for rule, v in out.items():
        ns = v['sizes']
        w = str(v['best_w2'])
        print('%-9s %5d %8s  %s  %s' % (
            rule, v['best_w2'],
            ' '.join('%d/%d' % (v['cells'][str(n)]['w2'][w]['exact_meet'],
                                v['cells'][str(n)]['nmeet']) for n in ns),
            ' '.join('%.1f' % v['cells'][str(n)]['w2'][w]['e'] for n in ns),
            'EXACT' if v['exact_all'] else 'not exact'))

    print('\nT4  failure anatomy at that W2 (meeting tapes)')
    print('%-9s %5s %6s %6s %6s %7s %7s %8s %8s' % (
        'rule', 'W2', 'bad', 'cut', 'term', 'len<=4', 'len<=16', 'M0 med/p90',
        'first j/n'))
    for rule, v in out.items():
        ns = v['sizes']
        w = str(v['best_w2'])
        agg = {k: sum(v['cells'][str(n)]['w2'][w][k] for n in ns)
               for k in ('bad', 'bad_cut', 'bad_term', 'bad_len4', 'bad_len16',
                         'tapes_bad', 'bad_b2z', 'bad_b2o')}
        m0 = [v['cells'][str(n)]['w2'][w]['bad_m0'][0] for n in ns]
        fp = [v['cells'][str(n)]['w2'][w]['first_pos'][0] for n in ns]
        print('%-9s %5d %6d %6d %6d %7d %7d %7d %8s %8s' % (
            rule, v['best_w2'], agg['bad'], agg['bad_cut'], agg['bad_term'],
            agg['bad_len4'], agg['bad_len16'], agg['bad_b2z'],
            '/'.join('%.0f' % x for x in m0),
            '/'.join('%.2f' % x for x in fp)))

    print('\nT5  terminal M0: meeting med/p90/p99/max, halting med/p90/max')
    for rule, v in out.items():
        ns = v['sizes']
        print('%-9s %s | %s' % (
            rule,
            ' ; '.join('%d: %s' % (n, '/'.join('%.0f' % x for x in
                                              v['cells'][str(n)]['term_meet']))
                       for n in ns),
            ' ; '.join('%d: %s' % (n, '/'.join('%.0f' % x for x in
                                              v['cells'][str(n)]['term_halt']))
                       for n in ns)))

    with open(os.path.join(HERE, 'el_report.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    print('\nwrote el_report.json (%d rules)' % len(out))


if __name__ == '__main__':
    main()
