# -*- coding: utf-8 -*-
"""el_coincide -- tapes whose prediction failed above W2 but whose answer was
still right.  Prints the evidence so the coincidence can be judged."""
import json
import sys

f = sys.argv[1]
W = int(sys.argv[2])
d = json.load(open(f))
n = d['n']
meet = [t for t in d['taps'] if t['end'] == 0]
print('==', d['rule'], n, 'W2 =', W)
cnt = 0
for t in meet:
    above = [b for b in t['bad'] if b[0] > W]
    m = t['mach'][str(W)]
    if above and m['exact']:
        cnt += 1
        if cnt <= 12:
            print('  tape %4d  nep %4d termM0 %4d  bad-above %s  truth d=%.6f '
                  'mach d=%.6f  mach end %d'
                  % (t['i'], t['nep'], t['termM0'],
                     [(b[0], b[1], b[2], b[3]) for b in above[:3]],
                     t['delta'], m['delta'], m['end']))
print('  such tapes: %d / %d' % (cnt, len(meet)))
