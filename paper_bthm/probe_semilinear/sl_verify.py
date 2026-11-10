# -*- coding: utf-8 -*-
"""sl_verify.py -- exhaustive / random agreement of the in-class machine with the
literal reference engine, for the threshold-run rule.

n <= 10: every ensemble tape (two independent uniform h-subsets).
n = 11, 12: 4000 random ensemble tapes each.
Also checks the b1=0 shortcut against rs_common.ref_walk('b0').
"""
import sys
import itertools
import numpy as np

sys.path.insert(0, '../probe_ologn2')
from rs_common import ref_walk, sample_tapes  # noqa: E402
from sl_ref import ref_rule  # noqa: E402
from sl_machine import solve, t_of  # noqa: E402


def enum_ens(n):
    h = (n + 1) // 2
    subs = list(itertools.combinations(range(n), h))
    for x in subs:
        for y in subs:
            a = [0] * n
            for i in x:
                a[i] += 1
            for i in y:
                a[i] += 1
            yield a


def main():
    ns = [int(x) for x in sys.argv[1].split(',')]
    c = int(sys.argv[2]) if len(sys.argv) > 2 else -1
    bad = 0
    tot = 0
    bad_short = 0
    for n in ns:
        t = t_of(n, c)
        if n <= 10:
            tapes = list(enum_ens(n))
            kind = 'exhaustive'
        else:
            tapes = [list(a) for a in sample_tapes(n, 4000, 2026)]
            kind = 'random'
        for a in tapes:
            d0, _ = ref_rule(a, n, t)
            for store in (False, True):
                r = solve(a, t, store=store)
                tot += 1
                if (d0 is None) != (r['delta'] is None) or (
                        d0 is not None and abs(d0 - r['delta']) > 1e-12):
                    bad += 1
                    if bad <= 3:
                        print('  MISMATCH n=%d t=%d store=%s tape=%s ref=%s mach=%s mode=%s'
                              % (n, t, store, a, d0, r['delta'], r['mode']))
        # b0 shortcut against the graded engine
        for a in tapes[:len(tapes)]:
            d0, _ = ref_walk(np.asarray(a, dtype=np.int64), n, 'b0')
            r = solve(a, 10 ** 9)     # rho <= t always -> shortcut path
            if (d0 is None) != (r['delta'] is None) or (
                    d0 is not None and abs(d0 - r['delta']) > 1e-12):
                bad_short += 1
                if bad_short <= 3:
                    print('  SHORTCUT MISMATCH n=%d tape=%s ref=%s mach=%s'
                          % (n, a, d0, r['delta']))
        print('  n=%-3d %-10s t=%d tapes=%d rule-mismatch=%d shortcut-mismatch=%d'
              % (n, kind, t, len(tapes), bad, bad_short), flush=True)
    print('TOTAL tapes=%d rule mismatches=%d shortcut mismatches=%d'
          % (tot, bad, bad_short))


if __name__ == '__main__':
    main()
