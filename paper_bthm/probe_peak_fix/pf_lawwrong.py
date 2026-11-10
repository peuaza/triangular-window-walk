# -*- coding: utf-8 -*-
"""pf_lawwrong -- which tapes can the literal-first-epoch device even touch?

If the trace's epoch 0 satisfies pred_of(b2,b1) == out, then a literal epoch 0
ends with exactly the cut the jump machine would have taken, so every machine
that runs epoch 0 literally (the fixed one, the adaptive one, the sticky one)
has the same region after it as the base machine and follows the same jump
path: only the charge differs.  The machines can differ in (delta, end) only on
tapes where epoch 0 is a law-wrong epoch.

This script lists those tapes across the eight segments of the 1 500-tape grid
and prices the literal epoch 0 against the jump cut it replaces.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
SEGS = [
    ('512 main', 512, 0), ('1024 main', 1024, 0), ('2048 main', 2048, 0),
    ('4096 main', 4096, 0), ('512 repl', 512, 900000),
    ('1024 repl', 1024, 500000), ('2048 repl A', 2048, 700000),
    ('2048 repl B', 2048, 300000), ('4096 repl', 4096, 900000),
]


def main():
    nwrong = 0
    rows = []
    delta = []
    for name, n, off in SEGS:
        p = os.path.join(HERE, 'seg_%d_%d.jsonl' % (n, off))
        if not os.path.exists(p):
            print('MISSING', p)
            continue
        for line in open(p):
            r = json.loads(line)
            top = r['top'][0]
            b2, b1, out = int(top[3]), int(top[4]), int(top[2])
            pred = 3 if b2 == 0 else (2 if b1 else 1)
            jump = 2 if b2 == 0 else 1
            delta.append((top[1] - jump) / n)
            if pred != out:
                nwrong += 1
                rows.append((name, r['i'], n, top, pred))
    print('tapes with a law-wrong epoch 0: %d' % nwrong)
    for x in rows:
        print('  %-11s i=%-7d n=%4d epoch0=(M0=%d,len=%d,out=%d,b2=%d,b1=%d) pred=%d'
              % (x[0], x[1], x[2], *x[3], x[4]))
    d = np.array(delta)
    print('\ncharge of literal epoch 0 vs the jump cut it replaces, steps/n:')
    print('  mean %.5f  median %.5f  p90 %.4f  p99 %.4f  max %.4f'
          % (d.mean(), np.median(d), np.percentile(d, 90),
             np.percentile(d, 99), d.max()))


if __name__ == '__main__':
    main()
