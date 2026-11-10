import sys, os
import numpy as np
sys.path.insert(0, os.getcwd()); sys.path.insert(0,'../probe_e_uniform'); sys.path.insert(0,'../probe_kappa_pred'); sys.path.insert(0,'../probe_ologn2')
import eu_lib as EL, ce_lib as CL
n = int(sys.argv[1]); rule = sys.argv[2]; nt = int(sys.argv[3]) if len(sys.argv)>3 else 200
M0s=[]; i=0; nmeet=0
while nmeet < nt and i < 30*nt:
    a,s,X = EL.make_tape(n, 5000+n+100003*i)
    a=np.asarray(a,np.int64)
    d,why,eps = CL.trace(a,rule)
    if why=='meet' and eps:
        M0,b2,b1,oc,ln,L,R,ee = eps[-1]
        M0s.append(int(M0)); nmeet+=1
    i+=1
M0s=np.array(sorted(M0s))
print('%s n=%d meeting tapes=%d scanned=%d'%(rule,n,len(M0s),i))
print('  terminal M0: min %d med %.0f p90 %.0f p99 %.0f max %d'%(M0s.min(),np.median(M0s),np.percentile(M0s,90),np.percentile(M0s,99),M0s.max()))
for T in (14,32,64,128,256,512,1024):
    if T<n: print('  frac terminal M0 > %4d: %.3f'%(T,(M0s>T).mean()))
