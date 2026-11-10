# -*- coding: utf-8 -*-
"""el_run -- run a group of rules over the four sizes, one JSON per cell.

Usage: python el_run.py GROUP [--force]
Groups: A (the seven sliding rules), B (range max/min + BAL>m/2),
        C (BAL>m/4 + the three NBLK), or a comma-separated rule list.

Tape budget per cell: 200 meeting + 100 halting at n <= 2048, 100 + 60 at 4096
(the record's own practice at 4096 is 100 meeting tapes).
A cell whose JSON already exists is skipped unless --force.
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import el_cell as C
import el_lib as E

GROUPS = {
    'A': ['dense2', 'peak', 'trip', 'oddrun', 'canonical', 'LAS', 'KCHG'],
    'B': ['MAXP', 'AMAX', 'AMIN', 'BAL>m/2'],
    'C': ['BAL>m/4', 'NBLK=1', 'NBLK>2', 'NBLK>3'],
}


def budget(n):
    if n >= 4096:
        return 100, 60, 40000
    return 200, 100, 60000


def main():
    g = sys.argv[1]
    force = '--force' in sys.argv
    rules = GROUPS.get(g, g.split(','))
    t0 = time.time()
    for rule in rules:
        for n in E.NS:
            p = os.path.join(HERE, 'el_%s_%d.json' % (C.safe(rule), n))
            if os.path.exists(p) and not force:
                print('skip', rule, n)
                continue
            nmeet, nhalt, cap = budget(n)
            t1 = time.time()
            d = C.run(rule, n, nmeet, nhalt, cap)
            import numpy as np
            meet = [x for x in d['taps'] if x['end'] == 0]
            eold = np.mean([x['inc'] for x in meet]) / n
            print('%-9s n=%4d %-7s meet=%3d/%3d halt=%3d draws=%5d '
                  'secs=%5.0f e_old/n=%7.1f  [group %s, %.0f s elapsed]'
                  % (rule, n, d['mode'], d['nmeet'], nmeet, d['nhalt'],
                     d['ndraw'], time.time() - t1, eold, g, time.time() - t0),
                  flush=True)


if __name__ == '__main__':
    main()
