import numpy as np, sys
sys.path.insert(0,'.')
import gp2; gp2.MATCH=True
from corr1d import oned, KAP
CANDS=gp2.CANDS
for mu in [0.02,0.01,0.005]:
  for truth in ['PL','PSh']:
    t=1.0; n=(t/mu**3)**2
    pr=gp2.predict(truth,mu,n,B=200000)
    grp={k:round(sum(pr[i] for i,c in enumerate(CANDS[:-1]) if KAP[c]==k),3) for k in sorted(set(KAP.values()))}; grp['NB']=round(pr[-1],3)
    print(mu,truth,'full',grp,' 1D',{k:round(v,3) for k,v in oned(truth,t).items()})
