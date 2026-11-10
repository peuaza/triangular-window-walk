# -*- coding: utf-8 -*-
"""el_predgate -- the BAL/NBLK predicates against a literal prefix walk.

bal_b1 (tr_bal4) and nblk_b1 (el_lib) must equal the literal definitions on
every window of every tape at n = 32, 64.
"""
import sys

import numpy as np

sys.path.insert(0, '.')
import el_lib as E
import pcr_lib as P
import tr_bal4


def lit_bal(s, ell, r, thr):
    m = r - ell + 1
    acc = 0
    best = 0
    for i in range(ell, r + 1):
        acc += int(s[i])
        if acc == 0:
            best = i - ell + 1
    return 1 if thr * best > m else 0


def lit_nblk(s, ell, r, sat):
    acc = 0
    cnt = 0
    for i in range(ell + 1, r + 2):
        acc += int(s[i - 1])
        if acc == 0:
            cnt += 1
    if sat == 1:
        return 1 if cnt == 1 else 0
    return 1 if cnt >= sat else 0


bad = 0
tot = 0
for n in (32, 64):
    for i in range(6):
        X = E.tape(n, i)
        S = P.St(X)
        s = S.s
        for ell in range(n):
            for r in range(ell + 1, n):
                tot += 1
                for thr in (2, 4):
                    a1 = int(tr_bal4.bal_b1(S.a, S.Pa, S.Ps, S.psmin, S.bstart,
                                            S.bpos, thr, ell, r))
                    b1 = lit_bal(s, ell, r, thr)
                    if a1 != b1:
                        bad += 1
                        if bad < 8:
                            print('BAL thr=%d n=%d tape=%d [%d,%d] mach %d lit %d'
                                  % (thr, n, i, ell, r, a1, b1))
                for sat in (1, 3, 4):
                    a1 = int(E.nblk_b1(S.Pa, S.Ps, S.psmin, S.bstart, S.bpos,
                                       sat, ell, r))
                    b1 = lit_nblk(s, ell, r, sat)
                    if a1 != b1:
                        bad += 1
                        if bad < 8:
                            print('NBLK sat=%d n=%d tape=%d [%d,%d] mach %d lit %d'
                                  % (sat, n, i, ell, r, a1, b1))
print('windows', tot, 'predicate mismatches', bad)
