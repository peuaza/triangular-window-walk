# -*- coding: utf-8 -*-
"""改装族的大样本：python _accumulate_mod.py <variant> <n> <target> <seconds>
variant: self | b2last | b2half
"""
import io
import json
import os
import sys
import time
import numpy as np

VARIANT = sys.argv[1]
N = int(sys.argv[2])
TARGET = int(sys.argv[3])
BUDGET = float(sys.argv[4])

def walk(a0, variant):
    a = np.array(a0, dtype=float, copy=True)
    n = len(a)
    y0, x1 = 0, n - 1
    while x1 - y0 >= 3:
        y, x = y0, x1
        seen = set()
        while True:
            if (y, x) in seen:
                y0, x1 = y0 + 1, x1 - 1
                break
            seen.add((y, x))
            m = x - y + 1
            s = a[y:x + 1]
            if variant == 'b2half':
                b2 = s.mean() > 0.5
            elif variant == 'b2last':
                b2 = s.mean() > s[m // 2:].mean()
            else:
                b2 = s.mean() > s[:max(1, m // 2)].mean()
            if m == 1:
                b1 = a[y] > 0.5
            elif m == 2:
                b1 = a[x] > a[y]
            else:
                u, v = s[:-1], s[1:]
                b1 = (np.corrcoef(u, v)[0, 1] > 0) if (u.std() > 1e-12 and v.std() > 1e-12) else False
            if variant == 'self':
                a[y] = 1.0 - a[y]
            o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
            ny, nx = {"A": (y - 1, x), "B": (y + 1, x), "C": (y, x + 1), "D": (y, x - 1)}[o]
            if ny < y0 or nx > x1:
                if ny < y0:
                    x1 = x - 1
                else:
                    y0 = y + 1
                break
            if ny >= nx:
                return ny
            if nx < y0 or ny > x1:
                y0, x1 = y0 + 1, x1 - 1
                break
            y, x = ny, nx
    return None

path = '_mod_%s_%d.jsonl' % (VARIANT, N)
done = 0
if os.path.exists(path):
    with io.open(path, encoding='utf-8') as f:
        for ln in f:
            try:
                if json.loads(ln)['x'] is not None:
                    done += 1
            except Exception:
                pass

h = (N + 1) // 2
rng = np.random.default_rng(abs(hash((VARIANT, N))) % (2 ** 31))
t0 = time.time()
walks = 0
f = io.open(path, 'a', encoding='utf-8', buffering=1)
while done < TARGET and time.time() - t0 < BUDGET:
    B = np.zeros((N, 2), dtype=np.int8)
    for c in range(2):
        B[rng.choice(N, h, replace=False), c] = 1
    d = walk(B.mean(1), VARIANT)
    f.write(json.dumps({'x': (None if d is None else float(d) / (N - 1))}) + '\n')
    walks += 1
    if d is not None:
        done += 1
f.close()
print('DONE %s n=%d: %d/%d ｜ %d walks ｜ %.0fs' % (VARIANT, N, done, TARGET, walks, time.time() - t0))
