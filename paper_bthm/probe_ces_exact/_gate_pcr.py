import sys, os
import numpy as np
sys.path.insert(0, os.getcwd())
sys.path.insert(0, '../probe_e_uniform'); sys.path.insert(0, '../probe_class_rerun'); sys.path.insert(0, '../probe_kappa_pred')
import eu_lib as EL, ce_lib as CL, pcr_lib as P
n = int(sys.argv[1]) if len(sys.argv)>1 else 512
nt = int(sys.argv[2]) if len(sys.argv)>2 else 10
bad = 0
for i in range(nt):
    a,s,X = EL.make_tape(n, 5000+n+100003*i)
    a = np.asarray(a, np.int64)
    out, steps, msum, inc, cuts, end, snaps = P.Runner(X).run('pattern')
    d2, why, eps = CL.trace(a,'pattern')
    ncut = sum(1 for e in eps if e[3] in ('D','R','L'))
    okd = (out < 0 and (d2 is None or d2 != d2)) or (out >= 0 and d2 == d2 and abs(out-d2)<1e-12)
    okc = (cuts == ncut)
    if not (okd and okc):
        bad += 1
        print('i',i,'pcr',out,cuts,'kp',d2,'ncut',ncut,'why',why,'okd',okd,'okc',okc)
print('gate n=%d nt=%d bad=%d' % (n,nt,bad))
