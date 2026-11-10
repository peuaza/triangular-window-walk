# -*- coding: utf-8 -*-
"""check_short.py -- (i) the short-window D-event really implies cost >= 2J;
(ii) the assembled constant as a function of the multiset parameter z (v = 2z/M);
(iii) stability of the assembled constant with M.
"""
import sys
import numpy as np

sys.path.insert(0, '.')
from q0_walk import sample_tapes, prefix, epoch


def cost_of(A):
    nt = len(A)
    out = np.empty(nt, dtype=np.int64)
    for i, a in enumerate(A):
        P = prefix(a)
        _, c = epoch(P, 0, len(a) - 1)
        out[i] = c
    return out


def window_event(A, delta):
    M = A.shape[1]
    B = A - 1
    S = np.zeros((A.shape[0], M + 1), dtype=np.int64)
    np.cumsum(B, axis=1, out=S[:, 1:])
    J = int(round(delta * M))
    u = np.arange(M // 2 - J, M // 2 + 1)
    D = S[:, 2 * u] - 2 * S[:, u]
    return (D >= 3).all(1), J


def assembled(A, delta):
    """union-bound lower bound for P(D_u>=3 on window), with lambda = 4 sqrt J."""
    M = A.shape[1]
    B = A - 1
    S = np.zeros((A.shape[0], M + 1), dtype=np.int64)
    np.cumsum(B, axis=1, out=S[:, 1:])
    J = int(round(delta * M))
    u0 = M // 2 - J
    p = S[:, u0:M // 2 + 1] - S[:, u0:u0 + 1]
    pJ = p[:, -1]
    s = np.arange(J + 1)
    dP = p - (s[None, :] / J) * pJ[:, None]
    q = S[:, M - 2 * J:] - S[:, M - 2 * J:M - 2 * J + 1]
    t = np.arange(2 * J + 1)
    dC = q - (t[None, :] / (2 * J)) * q[:, -1][:, None]
    lam = 4.0 * np.sqrt(J)
    Du0 = S[:, 2 * u0] - 2 * S[:, u0]
    x = S[:, M // 2]
    e1 = (Du0 >= 3 + lam) & (-2 * x >= 3 + lam)
    okP = np.abs(dP).max(1) <= lam / 4
    okC = np.abs(dC).max(1) <= lam / 2
    v = float(np.mean(B != 0))
    bp = min(1.0, 2 * np.exp(-2 * (lam / 4) ** 2 / (J * v)))
    bc = min(1.0, 2 * np.exp(-2 * (lam / 2) ** 2 / (2 * J * v)))
    return max(0.0, e1.mean() - bp - bc), e1.mean(), bp, bc, J


def main():
    # (i) implication
    print('(i) short-window D-event => cost >= 2J')
    for M in [512, 1024, 2048]:
        for delta in [0.01, 0.02, 0.04]:
            A = np.array(list(sample_tapes(M, 3000, 71)))
            ev, J = window_event(A, delta)
            c = cost_of(A)
            bad = np.sum(ev & (c < 2 * J))
            print('   M=%4d delta=%.2f J=%3d  P(ev)=%.3f  violations=%d  min(cost|ev)=%d  2J=%d'
                  % (M, delta, J, ev.mean(), bad, c[ev].min() if ev.any() else -1, 2 * J))
    # (ii) multiset dependence
    print('\n(ii) assembled bound vs z  (M=1024, delta=0.02)')
    rng = np.random.default_rng(5)
    M = 1024
    for z in [8, 32, 64, 128, 192, 256, 320, 384, 448, 500]:
        A = np.empty((1500, M), dtype=np.int64)
        for i in range(1500):
            a = np.ones(M, dtype=np.int64)
            pos = rng.choice(M, 2 * z, replace=False)
            a[pos[:z]] = 0
            a[pos[z:]] = 2
            rng.shuffle(a)
            A[i] = a
        val, pe1, bp, bc, J = assembled(A, 0.02)
        c = cost_of(A)
        print('   z=%3d v=%.2f  E[cost]/M=%.4f  P(ev-ish)=%.4f  bound=%.4e  ->  c=%.4f'
              % (z, 2 * z / M, c.mean() / M, pe1, val, 2 * J / M * val))
    # (iii) M stability at delta=0.02, standard ensemble
    print('\n(iii) M stability (delta=0.02, two-column ensemble)')
    for M in [512, 1024, 2048, 4096, 8192]:
        A = np.array(list(sample_tapes(M, 4000 if M <= 4096 else 1200, 73)))
        val, pe1, bp, bc, J = assembled(A, 0.02)
        c = cost_of(A)
        print('   M=%5d J=%4d  E[cost]/M=%.4f  P(ev)=%.3f  P(E1)=%.4f  bound=%.4f  c=%.4f'
              % (M, J, c.mean() / M, window_event(A, 0.02)[0].mean(), pe1, val, 2 * J / M * val))


if __name__ == '__main__':
    main()
