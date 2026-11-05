import sympy as sp
s,mu,r=sp.symbols('s mu r')
N=6
FAM={'geo':[1],'PL':[1,1],'PSh':[1,s],'PAk':[1,0,1],'PSu':[1,1,1],'PAm':[1,1,1,1],'PAr':[1,2,1],'PIs':[1,0,s]}  # c_k/c_0 as fn of s
def binom(k): return sp.expand(sp.prod([(r+j) for j in range(1,k+1)])/sp.factorial(k))
out={}
for f,rk in FAM.items():
    v=[sp.factorial(k)*s**k*rk[k] for k in range(len(rk))]; Z=sum(v)
    A=sp.series(sum(v[k]*(k+1) for k in range(len(v)))/Z,s,0,N).removeO()
    B=sp.series(sum(v[k]*binom(k) for k in range(len(v)))/Z,s,0,N).removeO()
    logR=sp.series(sp.log(B)-r*sp.log(A),s,0,N).removeO()
    c=sp.symbols('c1:6'); sg=mu+sum(c[i]*mu**(i+2) for i in range(4))
    eq=sp.expand(sp.series(sp.expand(s*A).subs(s,sg)-mu,mu,0,N).removeO())
    sol=sp.solve([eq.coeff(mu,i) for i in range(2,N)],c[:4],dict=True)[0]
    R=sp.expand(sp.series(sp.exp(logR.subs(s,sg.subs(sol))),mu,0,5).removeO())
    out[f]=R
    print(f,[sp.factor(R.coeff(mu,i)) for i in range(5)])
import itertools
print('pairs: order q of min ord DeltaR_r and whether DeltaR_2 achieves it')
for F,G in itertools.combinations(FAM,2):
    d=sp.expand(out[F]-out[G])
    q=next(i for i in range(5) if sp.simplify(d.coeff(mu,i))!=0)
    co=sp.factor(d.coeff(mu,q)); print(F,G,'q=',q,'coef(r)=',co,' at r=2:',co.subs(r,2))
