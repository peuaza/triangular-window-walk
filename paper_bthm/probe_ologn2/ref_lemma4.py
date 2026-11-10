# -*- coding: utf-8 -*-
"""Referee check 1: the constant p(c) of Lemma 4 (REPORT_EPOCH_BOUND).
Evaluates (##) exactly as stated and with the factor 1/2 that conditioning on b>0 requires.
"""
import numpy as np
from scipy.integrate import quad
from scipy.stats import norm

def m_dudley(c):
    return 1.2533 * np.sqrt(3.0 * c)

def E_absb_exp(c):
    """E_{b~|N(0,1)|}[ exp(-((b(1-2c)-m_c)_+)^2/(9/8)) ]"""
    m = m_dudley(c)
    s2 = 9.0 / 8.0
    f = lambda b: 2.0 * norm.pdf(b) * np.exp(-max(b * (1 - 2 * c) - m, 0.0) ** 2 / s2)
    val, err = quad(f, 0, 20, limit=400)
    return val, err

print(" c      m_c    1-E[.]  0.5*(1-E[.])   true p(c) (sim, below)")
for c in (0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30):
    e, err = E_absb_exp(c)
    print("%.2f  %.4f  %.4f   %.4f" % (c, m_dudley(c), 1 - e, 0.5 * (1 - e)))

# ---- true p(c) by Brownian simulation of V(theta)=W(theta)-2W(theta/2)
rng = np.random.default_rng(12345)
N = 4096          # grid on [0,1]
npath = 30000
dW = rng.normal(0.0, np.sqrt(1.0 / N), size=(npath, N))
W = np.concatenate([np.zeros((npath, 1)), np.cumsum(dW, axis=1)], axis=1)  # W[k]=W(k/N)
th = np.arange(N + 1) / N
# V(theta) = W(theta) - 2 W(theta/2); theta/2 lies on the grid when theta=k/N, k even;
# for odd k interpolate linearly.
idx = np.arange(N + 1)
half = idx / 2.0
lo = np.floor(half).astype(int); hi = np.minimum(lo + 1, N); w = half - lo
W2 = W[:, lo] * (1 - w) + W[:, hi] * w
V = W - 2.0 * W2
print("\n True p(c) by simulation (%d paths, grid %d):" % (npath, N))
for c in (0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50):
    k0 = int(np.ceil((1 - c) * N))
    ok = np.all(V[:, k0:] < 0, axis=1)
    print("  c=%.2f  p=%.4f" % (c, ok.mean()))
