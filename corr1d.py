import numpy as np, sys
sys.path.insert(0,'.')
from gauss_pred import predict, CANDS
KAP={'PGD1':0,'PL':-2,'PSh':0,'PAk':4,'PS':2,'PA':2}
def oned(truth,t,B=400000,rng=np.random.default_rng(2)):
    W=rng.standard_normal(B); one=[c for c in CANDS if c!='NB']
    H=np.array([t*(KAP[c]-KAP[truth])/2 for c in one])
    crit=(W[:,None]-H[None,:])**2
    best=crit.min(1); NBsel=best>2
    # group probabilities by kappa
    res={}
    for k in sorted(set(KAP.values())):
        g=[i for i,c in enumerate(one) if KAP[c]==k]
        res[k]=np.mean((~NBsel)&np.isin(crit.argmin(1),g))
    res['NB']=NBsel.mean(); return res
mu=0.05
for truth in ['PL','PSh','PAk']:
    for t in [0.5,2.0]:
        n=(t/mu**3)**2
        pr=predict(truth,mu,n,B=200000)
        grp={k:sum(pr[i] for i,c in enumerate(CANDS[:-1]) if KAP[c]==k) for k in sorted(set(KAP.values()))}; grp['NB']=pr[-1]
        o=oned(truth,t)
        print(truth,t,'full:',{k:round(v,3) for k,v in grp.items()},'\n        1D :',{k:round(v,3) for k,v in o.items()})
