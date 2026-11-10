import sys, os, json
import numpy as np
sys.path.insert(0, os.getcwd()); sys.path.insert(0,'../probe_e_uniform'); sys.path.insert(0,'../probe_kappa_pred'); sys.path.insert(0,'../probe_ologn2')
import eu_lib as EL, ce_machine as CM, ce_lib as CL
rule=sys.argv[1]; n=int(sys.argv[2]); T=int(sys.argv[3]) if len(sys.argv)>3 else 200
J=json.load(open('ce_fail_%s_%d.json'%(rule,n)))
idx=[]
for i,r in enumerate(J['recs']):
    ff=r.get('first_fail')
    if ff and ff['M0']>T: idx.append((i,ff['M0'],ff['true'],r['why'],r['nfail']))
print('%s n=%d tapes with first-fail M0>%d: %d'%(rule,n,T,len(idx)))
for (i,M0,oc,why,nf) in idx[:12]:
    a,s,X=EL.make_tape(n,5000+n+100003*i); a=np.asarray(a,np.int64)
    d,why2,eps=CL.trace(a,rule)
    truth=None if d is None or d!=d else float(d)
    out=[]
    for W2 in (64,128,256,512):
        if W2>=n: continue
        g,rd,inf=CM.machine(a,rule,W=(14 if rule=='pattern' else 40),W2=W2)
        gn = g is None or g!=g
        tag=('both_nan' if (truth is None and gn) else 'INVENTED' if truth is None else 'declined' if gn else 'exact' if abs(truth-g)<1e-12 else 'WRONG')
        out.append('W2=%d:%s'%(W2,tag))
    g,rd,inf=CM.machine(a,rule,W=(14 if rule=='pattern' else 40),oracle_term={len([e for e in eps if e[3] in ('D','R','L')])})
    gn = g is None or g!=g
    tag=('both_nan' if (truth is None and gn) else 'INVENTED' if truth is None else 'declined' if gn else 'exact' if abs(truth-g)<1e-12 else 'WRONG')
    print('  i=%3d termM0=%4d term=%s why=%s nfail=%d  %s  oracle:%s'%(i,M0,oc,why,nf,' '.join(out),tag))
