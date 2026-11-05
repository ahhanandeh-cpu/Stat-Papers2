import numpy as np, rdata, sys
from scipy.stats import nbinom
from scipy.optimize import minimize_scalar, minimize, brentq
from scipy.special import gammaln
X=np.arange(0,400)
FAM={'PGD1':lambda t:[1],'PL':lambda t:[1,1],'PSh':lambda t:[t,1],'PAk':lambda t:[1,0,1],
     'PS':lambda t:[1,1,1],'PA':lambda t:[1,1,1,1]}
def pmf1(f,lt):
    t=np.exp(lt); c=FAM[f](t); K=len(c)
    lv=np.array([np.log(c[k]) - (k+1)*lt + gammaln(k+1) if c[k]>0 else -np.inf for k in range(K)])
    w=np.exp(lv-lv.max()); w/=w.sum(); p=t/(1+t)
    return sum(w[k]*nbinom.pmf(X,k+1,p) for k in range(K) if w[k]>0)
def pmfnb(par): a=np.exp(par[0]); p=1/(1+np.exp(-par[1])); return nbinom.pmf(X,a,p)
def mean1(f,lt): return (X*pmf1(f,lt)).sum()
def ltheta(f,mu): return brentq(lambda lt: mean1(f,lt)-mu, -3, 8)
CANDS=['PGD1','PL','PSh','PAk','PS','PA','NB']
def fam_pmf(g,par): return pmfnb(par) if g=='NB' else pmf1(g,par[0])
def project(P,g,init):
    obj=lambda par: -(P*np.log(np.maximum(fam_pmf(g,par),1e-300))).sum()
    r=minimize(obj,init,method='Nelder-Mead',options={'xatol':1e-10,'fatol':1e-14,'maxiter':4000})
    return r.x
def predict(truth,mu,n,B=200000,rng=np.random.default_rng(1)):
    lt0=ltheta(truth,mu); P=pmf1(truth,lt0); m=P>1e-15; Pm=P[m]
    feats=[]; info=[]
    for g in CANDS:
        if g==truth: par=np.array([lt0])
        elif g=='NB':
            v=(X**2*P).sum()-mu**2; b=mu/(v-mu); a=mu*b; par=project(P,g,np.array([np.log(a),np.log(b)]))
        else: par=project(P,g,np.array([ltheta(g,mu)]))
        q=fam_pmf(g,par)[m]; L=np.log(q/Pm)
        h=1e-4; d=len(par); S=[]; H=np.zeros((d,d))
        lq=lambda pp: np.log(np.maximum(fam_pmf(g,pp)[m],1e-300))
        for i in range(d):
            e=np.zeros(d); e[i]=h; S.append((lq(par+e)-lq(par-e))/(2*h))
            for j in range(d):
                f=np.zeros(d); f[j]=h
                H[i,j]=-(Pm*(lq(par+e+f)-lq(par+e-f)-lq(par-e+f)+lq(par-e-f))).sum()/(4*h*h)
        feats.append((L,np.array(S),H,d))
    # joint covariance of [L_j (all j), scores]
    rows=[]; idx=[]
    for j,(L,S,H,d) in enumerate(feats):
        rows.append(L); 
    for j,(L,S,H,d) in enumerate(feats):
        for s in S: rows.append(s)
    M=np.array(rows); mean=(M*Pm).sum(1); C=((M-mean[:,None])*(M-mean[:,None])*1)@np.diag(Pm)@(M-mean[:,None]).T if False else ((M-mean[:,None])*Pm)@(M-mean[:,None]).T
    Z=rng.multivariate_normal(np.zeros(len(rows)),C,size=B,method='eigh')
    K=len(CANDS); aic=np.zeros((B,K)); off=K
    for j,(L,S,H,d) in enumerate(feats):
        zs=Z[:,off:off+d]; off+=d
        quad=np.einsum('bi,ij,bj->b',zs,np.linalg.inv(H),zs)
        ll=n*mean[j]+np.sqrt(n)*Z[:,j]+0.5*quad
        aic[:,j]=-2*ll+2*d
    sel=aic.argmin(1)
    return np.bincount(sel,minlength=K)/B
if __name__=='__main__':
    d=rdata.read_rds(sys.argv[1])
    out=[]
    for _,row in d.iterrows():
        pr=predict(row['truth'],float(row['mu']),float(row['n']),B=100000)
        k=CANDS.index(row['truth'])
        out.append((row['truth'],row['mu'],int(row['n']),row['p_correct'],pr[k]))
        print(*out[-1], flush=True)
    import pickle; pickle.dump(out,open('pred.pkl','wb'))
