# -*- coding: utf-8 -*-
"""特征更新代价谱

对每种窗口特征，实现「O(log n) 状态下的增量更新」，在一个真实游走上
量平均每步读多少格、状态多少位。

预期：双峰 —— O(1) 或 Θ(m)。若出现 Θ(log n)，D 的机制前提被打破。
"""
import sys

import numpy as np

N = int(sys.argv[1]) if len(sys.argv) > 1 else 256
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 60
h = N // 2


class Mean:
    """S1 增量"""
    name = 'mean'
    def __init__(self):
        self.s = 0.0; self.bits = 0
    def build(self, a, lo, hi):
        self.s = float(a[lo:hi + 1].sum()); return hi - lo + 1
    def update(self, a, olo, ohi, lo, hi):
        # 四种情形：左端右移(减)、左端左移(加)、右端左移(减)、右端右移(加)
        if lo > olo:
            self.s -= a[olo]; return 1
        if lo < olo:
            self.s += a[lo]; return 1
        if hi < ohi:
            self.s -= a[ohi]; return 1
        self.s += a[hi]; return 1
    def value(self, lo, hi):
        return self.s / (hi - lo + 1)


class Var(Mean):
    name = 'variance'
    def __init__(self):
        self.s = 0.0; self.q = 0.0
    def build(self, a, lo, hi):
        self.s = float(a[lo:hi + 1].sum()); self.q = float((a[lo:hi + 1] ** 2).sum())
        return 2 * (hi - lo + 1)
    def update(self, a, olo, ohi, lo, hi):
        if lo > olo: x = a[olo]; self.s -= x; self.q -= x * x
        elif lo < olo: x = a[lo]; self.s += x; self.q += x * x
        elif hi < ohi: x = a[ohi]; self.s -= x; self.q -= x * x
        else: x = a[hi]; self.s += x; self.q += x * x
        return 1
    def value(self, lo, hi):
        m = hi - lo + 1
        return self.q / m - (self.s / m) ** 2


class AC(Mean):
    name = 'ac1'
    def __init__(self):
        self.s = self.q = self.r = 0.0
    def build(self, a, lo, hi):
        self.s = float(a[lo:hi + 1].sum()); self.q = float((a[lo:hi + 1] ** 2).sum())
        self.r = float((a[lo:hi] * a[lo + 1:hi + 1]).sum())
        return 3 * (hi - lo + 1)
    def update(self, a, olo, ohi, lo, hi):
        # r 项：收缩时减去涉及的乘积，扩张时加上
        if lo > olo:
            x = a[olo]; self.s -= x; self.q -= x * x; self.r -= x * a[olo + 1]
        elif lo < olo:
            x = a[lo]; self.s += x; self.q += x * x; self.r += x * a[lo + 1]
        elif hi < ohi:
            x = a[ohi]; self.s -= x; self.q -= x * x; self.r -= a[ohi - 1] * x
        else:
            x = a[hi]; self.s += x; self.q += x * x; self.r += a[hi - 1] * x
        return 2
    def value(self, lo, hi):
        m = hi - lo + 1; nu = m - 1
        if nu < 1:
            return 0.0
        mu = (self.s - a_glob[hi]) / nu
        mv = (self.s - a_glob[lo]) / nu
        vu = (self.q - a_glob[hi] ** 2) / nu - mu * mu
        vv = (self.q - a_glob[lo] ** 2) / nu - mv * mv
        if vu <= 1e-24 or vv <= 1e-24:
            return 0.0
        return (self.r / nu - mu * mv) / np.sqrt(vu * vv)


class Half(Mean):
    """窗口均值 vs 前半均值：需要同时维护窗口和与前一半和"""
    name = 'mean-vs-half'
    def __init__(self):
        self.s = 0.0; self.sh = 0.0; self.k = 0
    def build(self, a, lo, hi):
        m = hi - lo + 1; k = max(1, m // 2)
        self.s = float(a[lo:hi + 1].sum()); self.sh = float(a[lo:lo + k].sum()); self.k = k
        return m + k
    def update(self, a, olo, ohi, lo, hi):
        m = hi - lo + 1; k = max(1, m // 2)
        r = 0
        if lo > olo: self.s -= a[olo]; r += 1
        elif lo < olo: self.s += a[lo]; r += 1
        elif hi < ohi: self.s -= a[ohi]; r += 1
        else: self.s += a[hi]; r += 1
        # 前半区间滑动：读变化的 0..2 格
        while self.k > k:
            self.sh -= a[lo + self.k - 1]; self.k -= 1; r += 1
        while self.k < k:
            self.sh += a[lo + self.k]; self.k += 1; r += 1
        return r
    def value(self, lo, hi):
        m = hi - lo + 1; k = max(1, m // 2)
        return (self.s / m) - (self.sh / k)


class Pattern(Mean):
    """窗口开头 L 连 1（L = log4 n）"""
    name = 'pattern-L'
    def __init__(self):
        self.run = 0; self.L = 0
    def build(self, a, lo, hi):
        self.L = max(2, int(round(np.log(max(2, hi - lo + 1)) / np.log(4))))
        self.run = 0
        for i in range(lo, min(hi, lo + self.L - 1) + 1):
            if a[i] > 0.9:
                self.run += 1
            else:
                break
        return min(self.L, hi - lo + 1)
    def update(self, a, olo, ohi, lo, hi):
        # 左端变化 ⇒ 重算 run（O(L)）；右端变化 ⇒ run 不变
        if lo == olo and hi != ohi:
            return 0
        self.run = 0
        r = 0
        for i in range(lo, min(hi, lo + self.L - 1) + 1):
            r += 1
            if a[i] > 0.9:
                self.run += 1
            else:
                break
        return r
    def value(self, lo, hi):
        return 1.0 if self.run >= self.L else 0.0


class Palindrome(Mean):
    """窗口是否回文（不可增量 ⇒ 每步 Θ(m)）"""
    name = 'palindrome'
    def __init__(self):
        self.flag = 0.0
    def build(self, a, lo, hi):
        self.flag = 1.0 if np.allclose(a[lo:hi + 1], a[lo:hi + 1][::-1]) else 0.0
        return hi - lo + 1
    def update(self, a, olo, ohi, lo, hi):
        self.flag = 1.0 if np.allclose(a[lo:hi + 1], a[lo:hi + 1][::-1]) else 0.0
        return hi - lo + 1
    def value(self, lo, hi):
        return self.flag


def walk_with(a, n, feat):
    """跑规范规则；每步调用 feat.update 并记录读次数"""
    global a_glob
    a_glob = a
    L, R = 0, n - 1
    ell, r = 0, n - 1
    reads = feat.build(a, ell, r)
    steps = 0
    while True:
        if ell >= r:
            return reads, steps
        if R - L + 1 < 3:
            return reads, steps
        steps += 1
        if steps > 20 * n:
            return reads, steps
        m = r - ell + 1
        k = max(1, m // 2)
        if m == 1:
            b2 = True
        else:
            b2 = (feat.value(ell, r) > 0) if feat.name == 'mean-vs-half' else (
                a[ell:r + 1].mean() > a[ell:ell + k].mean())
        if m <= 2:
            b1 = a[ell] > 0.5 if m == 1 else a[r] > a[ell]
        else:
            u, v = a[ell:r], a[ell + 1:r + 1]
            b1 = (np.corrcoef(u, v)[0, 1] > 0) if (u.std() > 1e-12 and v.std() > 1e-12) else False
        o = ("B" if b1 else "D") if b2 else ("C" if b1 else "A")
        ny, nx = {"A": (ell - 1, r), "B": (ell + 1, r), "C": (ell, r + 1), "D": (ell, r - 1)}[o]
        if ny < L or nx > R:
            if ny < L:
                R = r - 1
            else:
                L = ell + 1
            ell, r = L, R
            continue
        if ny >= nx:
            return reads, steps
        if nx < L or ny > R:
            L, R = L + 1, R - 1
            ell, r = L, R
            continue
        old_ell, old_r = ell, r
        ell, r = ny, nx
        reads += feat.update(a, old_ell, old_r, ell, r)
        # 说明：update 的 lo/hi 用新窗口；被移除的格在旧窗口端点上，这里用近似记账


rng = np.random.default_rng(7)
print('n=%d ｜ 每种特征 %d 条纸带' % (N, NW))
print('%-14s %14s %14s' % ('特征', '平均读/步', '读/n'))
for cls in (Mean, Var, AC, Half, Pattern, Palindrome):
    f = cls()
    tot_r = tot_s = 0
    for _ in range(NW):
        B = np.zeros((N, 2), dtype=np.int8)
        for c in range(2):
            B[rng.choice(N, h, replace=False), c] = 1
        rd, st = walk_with(B.mean(1), N, f)
        tot_r += rd; tot_s += st
        f = cls()   # 复位
    per = tot_r / max(1, tot_s)
    print('%-14s %14.2f %14.3f' % (f.name, per, per / N))
