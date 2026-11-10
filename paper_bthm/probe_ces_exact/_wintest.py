import sys, os, random
import numpy as np
sys.path.insert(0, os.getcwd())
import ce_machine as CM

def brute(a, ell, r, rule):
    M = r-ell+1; k = max(1, M>>1)
    S = sum(int(a[i]) for i in range(ell, r+1))
    C = sum(1 for i in range(ell, r) if a[i] != a[i+1])
    H = sum(int(a[i]) for i in range(ell, ell+k))
    PC = sum(1 for j in range(ell, r-1) if a[j]==2 and a[j+1]==1 and a[j+2]==1)
    return S, C, H, PC, M, k

random.seed(1)
for trial in range(30):
    n = 60
    a = np.random.default_rng(trial).integers(0, 3, n)
    rule = 'pattern'
    L = random.randint(0, 20); R = random.randint(L+5, n-1)
    w = CM.Win(a, L, R, rule)
    for step in range(200):
        S,C,H,PC,M,k = brute(a, w.ell, w.r, rule)
        if (w.S,w.C,w.H,w.PC,w.M,w.k) != (S,C,H,PC,M,k):
            print('MISMATCH at trial',trial,'step',step,'state',w.ell,w.r)
            print(' got ',(w.S,w.C,w.H,w.PC,w.M,w.k))
            print(' want',(S,C,H,PC,M,k))
            raise SystemExit
        ops = []
        if w.M > 3:
            ops += ['tl','tr']
        if w.ell > 0: ops += ['el']
        if w.r < n-1: ops += ['er']
        op = random.choice(ops)
        if op=='tl': w.trim_left()
        elif op=='tr': w.trim_right()
        elif op=='el': w.ext_left()
        else: w.ext_right()
print('window stats OK')
