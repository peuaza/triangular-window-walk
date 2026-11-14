# -*- coding: utf-8 -*-
"""cut_pilot -- cost and shape of the cut repair on a handful of tapes.

usage: cut_pilot.py [n] [ntapes] [rule ...]
prints one line per (rule, tape): mode-0 end/steps, mode-4 end/steps/ntrim,
the reads, the repairs split left/right and the seconds.  Nothing is written.
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cut_engine as E                                            # noqa: E402
import cut_verify as V                                            # noqa: E402
import tr_lib as T                                                # noqa: E402
import numpy as np                                                # noqa: E402


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 512
    nt = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    rules = sys.argv[3:] or ['canonical', 'dense2', 'NBLK>2', 'peak', 'oddrun']
    for rule in rules:
        t0 = time.time()
        ntr = []
        st = []
        ends = []
        for i, a, s, X in T.tapes(n, nt):
            C = E.Cyc(X)
            r4 = C.run(V.rrid(rule), V.conv(rule), V.variant(rule), mode=4)
            ntr.append(r4['ntrim'])
            st.append(r4['steps'])
            ends.append(r4['end'])
            print('%-9s n=%-5d tape %3d  end4=%d steps4=%-10d ntrim=%-6d '
                  'L/R=%d/%d reads=%-9d nused=%d'
                  % (rule, n, i, r4['end'], r4['steps'], r4['ntrim'],
                     r4['nl'], r4['nr'], r4['reads'], r4['nused']), flush=True)
        print('== %-9s n=%d  ntapes=%d  mean ntrim=%.1f max=%-6d mean steps=%.0f '
              'ends=%s  [%.0fs]'
              % (rule, n, nt, np.mean(ntr), max(ntr), np.mean(st),
                 {e: ends.count(e) for e in sorted(set(ends))},
                 time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
