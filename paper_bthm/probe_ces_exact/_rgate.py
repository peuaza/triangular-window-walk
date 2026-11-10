import sys, os
import numpy as np
sys.path.insert(0, os.getcwd()); sys.path.insert(0,'../probe_e_uniform'); sys.path.insert(0,'../probe_kappa_pred'); sys.path.insert(0,'../probe_ologn2')
import eu_lib as EL, kp_lib as KP, ce_machine as CM
n=int(sys.argv[1]); nt=int(sys.argv[2]) if len(sys.argv)>2 else 5
for rule,W in (('runs',max(16,4*int(round(np.log2(n))))),('pattern',14)):
    bad=0; diffs=[]
    for i in range(nt):
        a,s,X=EL.make_tape(n,5000+n+100003*i); a=np.asarray(a,np.int64)
        dk,rk,ik=KP.ces(a,rule,W=W)
        dm,rm,im=CM.machine(a,rule,W=W)
        ok = ((dk!=dk) and (dm!=dm)) or (dk==dk and dm==dm and abs(dk-dm)<1e-12)
        if not ok: bad+=1
        diffs.append(rm-rk)
    print('%s n=%d: output mismatches %d/%d ; reads(mine)-reads(kp): %s (n-2=%d)'%(rule,n,bad,nt,sorted(set(diffs)),n-2))
