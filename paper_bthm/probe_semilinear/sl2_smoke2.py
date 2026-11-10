# -*- coding: utf-8 -*-
"""sl2_smoke2 -- delta distribution and one-flip-pair sensitivity of the ballast rule."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); PAPER = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(PAPER, 'probe_class_rerun')); sys.path.insert(0, HERE)
import pcr_lib as P, sl2_engine as E

n = int(sys.argv[1]) if len(sys.argv) > 1 else 512
nt = int(sys.argv[2]) if len(sys.argv) > 2 else 60
npair = int(sys.argv[3]) if len(sys.argv) > 3 else 6
X = P.make_tapes(n, nt, 424242 + n)
D = []
base = []
for t in range(nt):
    S = P.St(X[t]); out, steps, msum, inc, cuts, end = E.Runner2(S).run('ballast')
    if end == 0: base.append(out)
base = np.array(base)
print('n=%d  meeting %d/%d   delta: min %.3f q25 %.3f med %.3f q75 %.3f max %.3f  sd(cells) %.1f'
      % (n, len(base), nt, base.min(), np.percentile(base,25), np.median(base),
         np.percentile(base,75), base.max(), base.std()*n))
print('  frac delta<0.1: %.2f   >0.9: %.2f   in [0.3,0.7]: %.2f'
      % ((base<0.1).mean(), (base>0.9).mean(), ((base>=0.3)&(base<=0.7)).mean()))
# one flip pair (2->1, 0->1), cell 0 fixed: keeps 2*#2 + #1 = n, i.e. the T_n multiset
rng = np.random.default_rng(11)
dd = []
for t in range(nt):
    tau = (X[t,:,0] + X[t,:,1]).astype(np.int64)
    S = P.St(X[t]); out, *_ = E.Runner2(S).run('ballast')
    if out < 0: continue
    ins = [i for i,v in enumerate(tau) if v==2 and i!=0]
    outs= [j for j,v in enumerate(tau) if v==0 and j!=0]
    rng.shuffle(ins); rng.shuffle(outs)
    for i,j in list(zip(ins,outs))[:npair]:
        b = tau.copy(); b[i]=1; b[j]=1
        Y = np.zeros((n,2), np.int8)
        # rebuild a 4-valued tape with the same sigma and a fixed label pattern
        sig = b
        Y[:,0] = (sig>0).astype(np.int8); Y[:,1] = (sig==2).astype(np.int8)
        # keep labels: only where sig==1 may the label differ; use the original
        Y[:,0] = np.where(sig==1, X[t,:,0], Y[:,0])
        Y[:,1] = np.where(sig==1, X[t,:,1], Y[:,1])
        Sb = P.St(Y); o2, *_ = E.Runner2(Sb).run('ballast')
        dd.append(n if o2 < 0 else abs(o2*(n-1) - out*(n-1)))
v = np.array(dd, float)
print('flip pairs: %d probes  zero %.1f%%  med %.1f  med(nz) %.1f  p90(nz) %.1f  max %.1f'
      % (len(v), 100*(v==0).mean(), np.median(v), np.median(v[v>0]) if (v>0).any() else 0,
         np.percentile(v[v>0],90) if (v>0).any() else 0, v.max()))
