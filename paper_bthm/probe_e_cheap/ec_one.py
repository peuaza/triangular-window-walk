"""One cell of the e re-measurement at large n (the decisive 2048 -> 8192 test).

Usage: python ec_one.py <rule> <n> <nt>

Writes ec_<rule>_<n>.json: the record machine's mean reads and the cheap machine's cost
(n + steps + cuts), over the first <nt> meeting tapes.  The cheap machine is only charged where
EN_LINEAR §8 says it is legal (a sliding statistic); for the level rules the re-span may be
necessary, so both numbers are reported and the caller decides.
"""
import json
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


def main(rule, n, nt):
    old, new, drawn = [], [], 0
    cap = max(4000, 200 * nt)
    for i in range(cap):
        if len(new) >= nt:
            break
        drawn = i + 1
        a, s, X = L.make_tape(n, 5000 + n + 100003 * i)
        r = T.run(rule, a, s)
        d = r.get('delta')
        if d is None or d != d:
            continue
        st = r.get('steps')
        if st is None:
            continue
        cu = 0 if r.get('cuts') is None else int(r['cuts'])
        old.append(float(r['reads']))
        new.append(float(n + int(st) + cu))
    out = dict(rule=rule, n=n, nt=len(new), drawn=drawn,
               e_old=float(np.mean(old)) if old else None,
               e_new=float(np.mean(new)) if new else None)
    if out['e_old']:
        out['e_old_n'] = out['e_old'] / n
        out['e_new_n'] = out['e_new'] / n
    with open(os.path.join(HERE, f"ec_{rule.replace('>', 'g').replace('/', 'd')}_{n}.json"), 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out), flush=True)


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
