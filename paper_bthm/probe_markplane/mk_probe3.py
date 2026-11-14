# -*- coding: utf-8 -*-
"""mk_probe3 -- exact decomposition of pcr_lib.walk's `inc` for rule pattern,
plus a census of the cells it is charged for (duplicate reads = the only thing
a per-cell memo could ever save)."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
for d in ('probe_triple_nat', 'probe_class_rerun', 'probe_e_uniform'):
    p = os.path.join(PAPER, d)
    if p not in sys.path:
        sys.path.insert(0, p)
sys.path.insert(0, HERE)

import tr_lib as T                                              # noqa: E402
import pcr_lib as P                                             # noqa: E402

TABLE = P.TABLE


def walk_trace(a, Pa, Ppat, n, cap):
    """pcr_lib.walk's control flow for rule 'pattern', with the charge split."""
    L, R_ = 0, n - 1
    ell, r = L, R_
    y, x = L, R_
    hy, hx = 0, max(1, n // 2) - 1
    steps = 0
    seen = set()
    d = dict(base=n + 2 * (n - 1) + 2 + max(1, n // 2),
             move=0, head=0, cut_span=0, cut_head=0, cuts=0,
             cells=0, dup=0, headreads=0, cutreads=0)
    seen_cells = set()
    while True:
        if ell >= r:
            end = 0
            break
        if R_ - L + 1 < 3:
            end = 1
            break
        steps += 1
        if steps > cap:
            end = 3
            break
        key = ell * (n + 1) + r
        if key in seen:
            end = 2
            break
        seen.add(key)
        m = r - ell + 1
        b1 = 0 if m < 3 else (1 if (Ppat[r - 1] - Ppat[ell]) > 0 else 0)
        k = m >> 1
        if k == 0:
            k = 1
        b2 = 1 if (Pa[r + 1] - Pa[ell]) * k > (Pa[ell + k] - Pa[ell]) * m else 0
        act = TABLE[2 * b2 + b1]
        if act >= P.CUTL:
            d['cuts'] += 1
            if act == P.CUTL:
                L += 1
            elif act == P.CUTR:
                R_ -= 1
            else:
                L += 1
                R_ -= 1
            if L > R_:
                end = 1
                break
            ell, r = L, R_
            c = abs(y - L) + abs(x - R_)
            d['cut_span'] += c
            d['cutreads'] += c
            y, x = L, R_
            k = (x - y + 1) >> 1
            if k == 0:
                k = 1
            ty, tx = y, y + k - 1
            c = abs(hy - ty) + abs(hx - tx)
            d['cut_head'] += c
            d['cutreads'] += c
            hy, hx = ty, tx
            continue
        dy = 1 if act == P.LPLUS else (-1 if act == P.LMINUS else 0)
        dx = 1 if act == P.RPLUS else (-1 if act == P.RMINUS else 0)
        if ell + dy < L or r + dx > R_:
            d['cuts'] += 1
            if ell == L and r == R_:
                L += 1
                R_ -= 1
                y += 1
                x -= 1
                d['cut_span'] += 2
                d['cutreads'] += 2
            elif ell + dy < L:
                L += 1
            else:
                R_ -= 1
            if L > R_:
                end = 1
                break
            ell, r = L, R_
            c = abs(y - L) + abs(x - R_)
            d['cut_span'] += c
            d['cutreads'] += c
            y, x = L, R_
            k = (x - y + 1) >> 1
            if k == 0:
                k = 1
            ty, tx = y, y + k - 1
            c = abs(hy - ty) + abs(hx - tx)
            d['cut_head'] += c
            d['cutreads'] += c
            hy, hx = ty, tx
            continue
        if ell + dy >= r + dx:
            end = 0
            break
        ell += dy
        r += dx
        y, x = ell, r
        d['move'] += 1
        # the cell taken in
        newc = ell if dy else r
        d['cells'] += 1
        if newc in seen_cells:
            d['dup'] += 1
        else:
            seen_cells.add(newc)
        k = (x - y + 1) >> 1
        if k == 0:
            k = 1
        ty, tx = y, y + k - 1
        c = abs(hy - ty) + abs(hx - tx)
        d['head'] += c
        d['headreads'] += c
        hy, hx = ty, tx
    d['steps'] = steps
    d['inc'] = d['base'] + d['move'] + d['head'] + d['cut_span'] + d['cut_head']
    d['end'] = end
    d['ndistinct'] = len(seen_cells)
    return d


if __name__ == '__main__':
    for n in (512, 1024):
        print('n =', n)
        for i, a, s, X in T.tapes(n, 4):
            S = P.St(X)
            d = walk_trace(S.a, S.Pa, S.Ppat, n, 40 * n * n + 100)
            chk = P.Runner(X).run('pattern')
            print('  tape %d inc=%-8d (pcr %-8d) base=%-6d move=%-7d head=%-7d '
                  'cutspan=%-5d cuthead=%-5d steps=%-7d cuts=%-4d dup=%-7d '
                  'distinct=%-5d end=%d/%d'
                  % (i, d['inc'], chk[3], d['base'], d['move'], d['head'],
                     d['cut_span'], d['cut_head'], d['steps'], d['cuts'],
                     d['dup'], d['ndistinct'], d['end'], chk[5]))
