# -*- coding: utf-8 -*-
"""
walk_spec.py —— 按《TWW-猜想-v3-离散版》§1.4「五种结局」实现的参考游走

事件与处理（与文档逐条对应）：
  出 y < y0        → x1 ← x-1，重启
  出 x > x1        → y0 ← y+1，重启
  到达 y = x       → 终止，输出 δ = y
  走进已访问的格    → 终止（cycle），无 δ
  起点即出界        → y0+1, x1-1，继续（两侧各弃一条边）
  m = 1            → ε，无 δ
  m = 2 且起步即出界 → ζ，无 δ

返回值：X = δ/(n-1) ∈ [0,1]，或 None。
"""
import numpy as np


def _bits(a, y, x):
    m = x - y + 1
    b2 = a[y:x + 1].mean() > a[y:y + max(1, m // 2)].mean()
    if m == 1:
        b1 = a[y] > 0.5
    elif m == 2:
        b1 = a[x] > a[y]
    else:
        sg = a[y:x + 1]
        b1 = (np.corrcoef(sg[:-1], sg[1:])[0, 1] > 0) if sg.std() > 1e-9 else (a[x] > a[y])
    return ('B' if b1 else 'D') if b2 else ('C' if b1 else 'A')


def _step(y, x, o):
    return {'A': (y - 1, x), 'B': (y + 1, x), 'C': (y, x + 1), 'D': (y, x - 1)}[o]


def walk_seq(a):
    n = len(a)
    y0, x1 = 0, n - 1
    while True:
        m = x1 - y0 + 1
        if m == 1:
            return None                      # ε
        if m == 2:
            o0 = _bits(a, y0, x1)
            ny, nx = _step(y0, x1, o0)
            if ny < y0 or nx > x1:
                return None                  # ζ
        y, x = y0, x1
        seen = set()
        first = True
        while True:
            if (y, x) in seen:
                return None                  # cycle → 终止
            seen.add((y, x))
            o = _bits(a, y, x)
            ny, nx = _step(y, x, o)
            if first and (ny < y0 or nx > x1):
                y0, x1 = y0 + 1, x1 - 1      # 起点即出界
                break
            if ny < y0:
                x1 = x - 1; break
            if nx > x1:
                y0 = y + 1; break
            if ny >= nx:
                return ny / (n - 1)          # δ
            y, x = ny, nx
            first = False


def walk_spec(B):
    """B: n×2 的 0/1 矩阵，每列恰 h 个 1"""
    return walk_seq(B.mean(1))
