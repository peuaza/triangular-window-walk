# -*- coding: utf-8 -*-
"""BTHM-0 电池：成环 = 退出（无记忆）

python _b0_battery.py strat <config> <tapes>
python _b0_battery.py nscan <n> <tapes>
python _b0_battery.py hybrid <mode> <tapes>
结果逐行追加到 _b0_battery.jsonl
"""
import io
import json
import sys

import numpy as np

TASK = sys.argv[1]
ARG = sys.argv[2]
NW = int(sys.argv[3]) if len(sys.argv) > 3 else 1500
OUT = '_b0_battery.jsonl'


def walk(a, cfg):
    n = len(a)
    y0, x1 = 0, n - 1
    kind, _, arg = cfg.partition(':')
    sub = ''
    if '+' in arg:
        arg, sub = arg.split('+', 1)
    while x1 - y0 >= 3:
        y, x = y0, x1
        seen = set()
        while True:
            if (y, x) in seen:
                return None                       # 成环 = 退出
            seen.add((y, x))
            m = x - y + 1
            s = a[y:x + 1]
            k = max(1, m // 2)
            b2 = s.mean() > s[:k].mean()
            if kind == 'pos':
                b2 = (y0 + x1) / 2 / n > float(arg)
            if m == 1:
                b1 = a[y] > 0.5
            elif m == 2:
                b1 = a[x] > a[y]
            else:
                u, v = s[:-1], s[1:]
                b1 = (np.corrcoef(u, v)[0, 1] > 0) if (u.std() > 1e-12 and v.std() > 1e-12) else False
            o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
            step = 2 if (kind == 'jump2' or sub == 'jump2') else 1
            ny, nx = {"A": (y - step, x), "B": (y + step, x),
                      "C": (y, x + step), "D": (y, x - step)}[o]
            if kind == 'slide':
                ny, nx = y + 1, x - 1
            if ny < y0 or nx > x1:
                if kind == 'centre' and abs((y0 + x1) / 2 / n - float(arg)) < 0.10:
                    y0, x1 = y0 + 1, x1 - 1
                elif sub == 'alt':
                    if (x1 - y0) % 2 == 0:
                        x1 = x - 1
                    else:
                        y0 = y + 1
                elif ny < y0:
                    x1 = x - 1
                else:
                    y0 = y + 1
                break
            if ny >= nx:
                return ny / (n - 1)
            if nx < y0 or ny > x1:
                y0, x1 = y0 + 1, x1 - 1
                break
            y, x = ny, nx
    return None


def ac1(a):
    u, v = a[:-1], a[1:]
    return 0.0 if u.std() < 1e-12 or v.std() < 1e-12 else float(np.corrcoef(u, v)[0, 1])


if TASK == 'strat':
    n, cfg = 256, ARG
elif TASK == 'nscan':
    n, cfg = int(ARG), 'base'
elif TASK == 'hybrid':
    n, cfg = 256, 'base'
else:
    raise SystemExit('unknown task')

h = (n + 1) // 2
rng = np.random.default_rng(90210)
X, A = [], []
for _ in range(NW):
    B = np.zeros((n, 2), dtype=np.int8)
    for c in range(2):
        B[rng.choice(n, h, replace=False), c] = 1
    a = B.mean(1)
    d = walk(a, cfg)
    if d is None:
        continue
    X.append(d)
    A.append(ac1(a))
X = np.array(X)
A = np.array(A)
rec = {'task': TASK, 'arg': ARG, 'n': n, 'tapes': NW, 'meet': int(len(X)),
       'rate': float(len(X) / NW), 'sd': float(X.std(ddof=1)) if len(X) > 1 else None,
       'support': int(len(np.unique(np.round(X, 6)))) if len(X) else 0, 'r2': None}
if len(X) > 30:
    b = np.polyfit(A, X, 1)
    rec['r2'] = float(1 - (X - np.polyval(b, A)).var() / X.var())
print(json.dumps(rec))
with io.open(OUT, 'a', encoding='utf-8') as f:
    f.write(json.dumps(rec) + '\n')
