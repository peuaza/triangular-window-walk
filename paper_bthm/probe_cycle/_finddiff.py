import numpy as np, sys
import cyc_engine as E, cyc_verify as V, cyc_indep as I, cyc_dump as D, tr_lib as T, pcr_lib as P
n=512
rule='hashmod'
tapes=list(T.tapes(n,5))
i,a,s,X = tapes[3]
S=P.St(X)
# python trace, recording
L,R=0,n-1; ell,r=L,R; y,x=L,R; hy,hx=0,max(1,n//2)-1
steps=0; py=[]
while True:
    if ell>=r: break
    if R-L+1<3: break
    steps+=1
    if steps>40*n*n+100: break
    py.append((ell,r,L,R,int(I.b1(rule,S,ell,r))))
    v1=py[-1][4]
    act=4 if v1==1 else 3   # CUTL if b1 else RMINUS
    if act>=4:
        L+=1
        if L>R: break
        ell,r=L,R; y,x=L,R
        k=max(1,(x-y+1)>>1); hy,hx=y,y+k-1
        continue
    dx=-1
    if r+dx>R:
        if ell==L and r==R:
            L+=1; R-=1; y+=1; x-=1
        elif ell+dx<L:
            L+=1
        else:
            R-=1
        if L>R: break
        ell,r=L,R; y,x=L,R
        k=max(1,(x-y+1)>>1); hy,hx=y,y+k-1
        continue
    if ell>=r+dx: break
    r+=dx; y,x=ell,r
    k=max(1,(x-y+1)>>1); hy,hx=y,y+k-1
print("python steps",steps,"states",len(py))
traj=np.zeros((20000,4),np.int32); b1=np.zeros(20000,np.int8)
st,en,out,nrec = D.dump(3, P.R_HASHMOD, 0, S.a,S.Pa,S.Paa,S.PR,S.Pd,S.Pn2,S.Ppk,S.Ptr,S.Ps,S.cpos,S.pnz,S.C2,S.bstart,S.bpos,S.psmin,S.SPs,S.SNs,S.LOG2,S.r2s,S.r2e,S.SP2,S.rcs,S.rce,S.rcodd,S.SPC,S.R,S.IP,S.HX,S.Ppat,S.d1,S.d2,n,40*n*n+100,traj,b1)
print("numba steps",st,"end",en,"nrec",nrec)
m=min(len(py),nrec)
for k in range(m):
    nb=(int(traj[k,0]),int(traj[k,1]),int(traj[k,2]),int(traj[k,3]))
    if py[k][:4]!=nb or py[k][4]!=int(b1[k]):
        print("first diff at k=%d"%k)
        for j in range(max(0,k-3),min(m,k+4)):
            print("  k=%d python %s b1=%d | numba %s b1=%d"%(j,py[j][:4],py[j][4],tuple(int(v) for v in traj[j]),int(b1[j])))
        break
else:
    print("no diff in first",m,"states")
