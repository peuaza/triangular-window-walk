# -*- coding: utf-8 -*-
"""pf_sticky_summary -- the per-size table of the sticky machine.

Identity used: if the trace's epoch 0 is law-correct (pred_of(b2,b1) == out),
a literal epoch 0 ends with exactly the cut the base machine's jump takes, so
the sticky machine follows the base machine's W2 = 256 trajectory from the same
region.  Its charge differs by (len(epoch 0) - the jump cut's 1 or 2 reads).
On a law-wrong epoch 0 the sticky machine is re-run on the tape and measured.

Columns: e/n actual (epoch 0 charged its measured length) and e/n capped
(epoch 0 charged n, the machine's worst-case price for that epoch).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'probe_en_linear'))

import pf_lib as F
import pf_sticky as S
import el_lib as E

SEGS = [
    ('512 main', 512, 0), ('1024 main', 1024, 0), ('2048 main', 2048, 0),
    ('4096 main', 4096, 0), ('512 repl', 512, 900000),
    ('1024 repl', 1024, 500000), ('2048 repl A', 2048, 700000),
    ('2048 repl B', 2048, 300000), ('4096 repl', 4096, 900000),
]
W = 256


def lawwrong(r):
    top = r['top'][0]
    b2, b1, out = int(top[3]), int(top[4]), int(top[2])
    pred = 3 if b2 == 0 else (2 if b1 else 1)
    return pred != out


def main():
    per_n = {}
    for name, n, off in SEGS:
        p = os.path.join(HERE, 'seg_%d_%d.jsonl' % (n, off))
        if not os.path.exists(p):
            print('MISSING', p)
            continue
        for line in open(p):
            r = json.loads(line)
            b = r['base%d' % W]
            top = r['top'][0]
            jump = 2 if int(top[3]) == 0 else 1
            if not lawwrong(r):
                exact = b['exact']
                cost = b['reads'] + top[1] - jump
                capped = b['reads'] + n - jump
            else:
                X = E.tape(n, r['i'])
                M = E.Mach(X, n)
                t = M.trace('peak')
                m = S.jm_sticky('peak', M.S, n, W)
                exact = int(F.exact_of(m, t, n))
                cost = m['reads']
                capped = cost + n - top[1]
                print('  law-wrong %-10s i=%-7d base256 exact=%d sticky exact=%d '
                      'sticky e/n=%.3f mode=%s'
                      % (name, r['i'], b['exact'], exact, cost / n, m['mode']))
            d = per_n.setdefault(n, [0, 0.0, 0.0, 0])
            d[0] += exact
            d[1] += cost / n
            d[2] += capped / n
            d[3] += 1
    print('\n%-5s %5s %9s %10s %10s' % ('n', 'tapes', 'exact', 'e/n actual', 'e/n capped'))
    tot = [0, 0.0, 0.0, 0]
    for n in sorted(per_n):
        ex, ea, ec, k = per_n[n]
        print('%-5d %5d %5d/%d %10.3f %10.3f'
              % (n, k, ex, k, ea / k, ec / k))
        tot[0] += ex
        tot[1] += ea
        tot[2] += ec
        tot[3] += k
    print('%-5s %5d %5d/%d %10.3f %10.3f'
          % ('all', tot[3], tot[0], tot[3], tot[1] / tot[3], tot[2] / tot[3]))


if __name__ == '__main__':
    main()
