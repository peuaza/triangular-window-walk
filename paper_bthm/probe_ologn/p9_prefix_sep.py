# -*- coding: utf-8 -*-
"""p9: one-pass streaming lower bound test.

If two prefixes p != p' (of length n/2, same count profile) have, for EVERY
admissible suffix s, delta(p.s) == delta(p'.s), then a one-pass streaming
machine may identify them and save bits.  Otherwise every one-pass machine
with O(log n) bits must keep the prefix (=> Omega(n) bits).

Exhaustive over all valid tapes for n = 6, 8, 10.
"""
import sys, itertools
from collections import defaultdict
import numpy as np
sys.path.insert(0, '..')
from p0_tapes import enum_tapes
from _arr_longrun_fast import walk_naive


def run(n):
    h = n // 2
    table = defaultdict(dict)          # prefix tuple -> {suffix tuple: delta}
    pref_counts = defaultdict(int)
    for a in enum_tapes(n):
        d, s, _ = walk_naive(a, n)
        key = -1 if d is None else int(round(d * (n - 1)))
        p = tuple(int(x) for x in a[:h])
        q = tuple(int(x) for x in a[h:])
        table[p][q] = key
        pref_counts[tuple(sorted(p))]  # noqa
    prefs = sorted(table.keys())
    print('n=%d: extendable prefixes=%d' % (n, len(prefs)))
    # group by count profile
    byprof = defaultdict(list)
    for p in prefs:
        byprof[tuple(p.count(v) for v in (0, 1, 2))].append(p)
    bad = 0
    tot = 0
    bad_ex = []
    for prof, ps in byprof.items():
        for i in range(len(ps)):
            for j in range(i + 1, len(ps)):
                p, p2 = ps[i], ps[j]
                common = set(table[p]) & set(table[p2])
                sep = any(table[p][q] != table[p2][q] for q in common)
                tot += 1
                if not sep:
                    bad += 1
                    if len(bad_ex) < 3:
                        bad_ex.append((p, p2))
    print('  same-profile prefix pairs=%d  NOT separated by any common suffix: %d (%.4f)'
          % (tot, bad, bad / max(1, tot)))
    for e in bad_ex:
        print('    example', e)


if __name__ == '__main__':
    for n in (6, 8, 10):
        run(n)
