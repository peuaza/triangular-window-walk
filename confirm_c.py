# -*- coding: utf-8 -*-
"""confirm_c.py —— 每档 4000 个 δ，把 P(δ) 的 SE 压到 0.007"""
import io, os, sys, time
import numpy as np

os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ".")
from fast_delta2 import walk_fast

rng = np.random.default_rng(7777)
print("%-6s %-8s %-9s %-8s %-9s" % ("n", "尝试", "δ", "δ率", "±SE"))
for n in (256, 512, 1024, 2048):
    target = 4000
    hit = 0; tries = 0
    t0 = time.time()
    while hit < target and tries < 40000 and time.time() - t0 < 3000:
        tries += 1
        if walk_fast(rng.random((n, 2))) is not None:
            hit += 1
    p = hit / tries
    se = np.sqrt(p * (1 - p) / tries)
    print("%-6d %-8d %-9d %-8.4f %-9.4f ｜距 1/3 %+.4f（%.1f SE）｜%.0fs"
          % (n, tries, hit, p, se, p - 1/3, (p - 1/3) / se, time.time() - t0), flush=True)
