import itertools, mpmath as mp, sys
sys.path.insert(0,'.')
from num import FAM, pmf, moments, div, nbmatched, KAP
mp.mp.dps=60; X=150
res={}
for mu in ['0.05','0.005']:
    mu=mp.mpf(mu); P={f:pmf(f,mu,X) for f in FAM}; D={}
    for f in FAM: m=moments(f,mu,3); D[f]=1+(m[2]-m[1]**2)/m[1]
    rat=[];hk=[];ok=[]
    for F,G in itertools.combinations(FAM,2):
        KL,H2,om=div(P[F],P[G]); rat.append(KL/((D[F]-D[G])**2/8)); hk.append(H2/KL); ok.append(om/KL)
    print(float(mu),'KL/law range',mp.nstr(min(rat),4),mp.nstr(max(rat),4),'H2/KL',mp.nstr(min(hk),4),mp.nstr(max(hk),4),'om/KL',mp.nstr(min(ok),4),mp.nstr(max(ok),4))
