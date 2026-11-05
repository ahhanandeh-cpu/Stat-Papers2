import mpmath as mp
mp.mp.dps=60
FAM={'geo':lambda t:[1],'PL':lambda t:[1,1],'PSh':lambda t:[t,1],'PAk':lambda t:[1,0,1],
     'PSu':lambda t:[1,1,1],'PAm':lambda t:[1,1,1,1],'PAr':lambda t:[1,2,1],'PIs':lambda t:[t,0,1]}
KAP={'geo':0,'PL':-2,'PSh':0,'PAk':4,'PSu':2,'PAm':2,'PAr':-4,'PIs':0}
def weights(f,t):
    c=FAM[f](t); v=[c[k]*mp.factorial(k)*t**(-(k+1)) for k in range(len(c))]; Z=sum(v); return [x/Z for x in v]
def mean(f,t): return sum(w*(k+1) for k,w in enumerate(weights(f,t)))/t
def theta(f,mu): return mp.findroot(lambda t: mean(f,t)-mu, 1/mu)
def nb(x,a,p): return mp.exp(mp.loggamma(x+a)-mp.loggamma(a)-mp.loggamma(x+1))*p**a*(1-p)**x
def pmf(f,mu,X):
    t=theta(f,mu); w=weights(f,t); p=t/(1+t)
    return [sum(w[k]*nb(x,k+1,p) for k in range(len(w))) for x in range(X)]
def moments(f,mu,R):
    t=theta(f,mu); w=weights(f,t)
    return [sum(w[k]*mp.factorial(k+r)/mp.factorial(k) for k in range(len(w)))/t**r for r in range(R)]
def nbmatched(f,mu,X):
    m=moments(f,mu,3); var=m[2]-m[1]**2  # mixing var
    b=m[1]/var; a=m[1]*b; p=b/(1+b)
    return [nb(x,a,p) for x in range(X)]
def div(P,Q):
    KL=sum(p*mp.log(p/q) for p,q in zip(P,Q)); H2=1-sum(mp.sqrt(p*q) for p,q in zip(P,Q))
    L=[mp.log(p/q) for p,q in zip(P,Q)]; om=sum(p*l*l for p,l in zip(P,L))-KL**2
    return KL,H2,om
X=200
print("second-moment law, ratio KL/((dD)^2/8), H2/KL, om/KL")
for F,G in [('PL','PSu'),('PL','PAk'),('PSh','PAk'),('PSu','PAm'),('PSh','PIs')]:
    for mu in [mp.mpf('0.05'),mp.mpf('0.005')]:
        P=pmf(F,mu,X);Q=pmf(G,mu,X)
        def D(f): m=moments(f,mu,3); return 1+(m[2]-m[1]**2)/m[1]
        dD=D(F)-D(G); KL,H2,om=div(P,Q)
        print(F,G,float(mu),mp.nstr(KL/(dD**2/8),5),mp.nstr(H2/KL,5),mp.nstr(om/KL,5), mp.nstr(dD,5))
print("third-moment law vs matched NB: KL/(kappa^2 mu^7/18)")
for F in ['PL','PAk','PSu','PAr']:
    for mu in [mp.mpf('0.05'),mp.mpf('0.01'),mp.mpf('0.002')]:
        P=pmf(F,mu,X);Q=nbmatched(F,mu,X);KL,H2,om=div(P,Q)
        print(F,float(mu),mp.nstr(KL/(KAP[F]**2*mu**7/18),6),mp.nstr(H2/KL,5),mp.nstr(om/KL,5))
