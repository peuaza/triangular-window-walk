"""`canonical`'s worst tape: is it structured or generic?

Backs the sentence in bthm_jalc.tex §classes: the witness panel does not find canonical's worst
tape because that tape is generic (311 runs of length 1-8, counts 132/248/132 at n = 512).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
for d in (os.path.join(PAPER, 'probe_triple_nat'), os.path.join(PAPER, 'probe_e_uniform')):
    if d not in sys.path:
        sys.path.insert(0, d)

import tr_lib as T                                              # noqa: E402
import eu_lib as L                                              # noqa: E402

N = 512
NT = 200


def main():
    best = (-1, None)
    for i in range(NT):
        a, s, X = L.make_tape(N, 5000 + N + 100003 * i)
        r = T.run('canonical', a, s)
        d = r.get('delta')
        if d is None or d != d:
            continue
        if r['reads'] > best[0]:
            best = (r['reads'], i)
    reads, i = best
    a, s, X = L.make_tape(N, 5000 + N + 100003 * i)
    v, c = np.unique(a, return_counts=True)
    seg = np.split(a, np.flatnonzero(np.diff(a) != 0) + 1)
    print(f'canonical n={N}: worst of {NT} tapes = tape {i}, reads {reads} = {reads / N ** 2:.3f} n^2')
    print('sigma counts:', dict(zip(v.tolist(), c.tolist())))
    print('runs:', len(seg), 'longest:', sorted(len(x) for x in seg)[-5:])
    print('first eight run lengths:', [len(x) for x in seg[:8]])


if __name__ == '__main__':
    main()
