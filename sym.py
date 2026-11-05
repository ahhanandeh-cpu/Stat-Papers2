import sympy as sp
s,mu,r=sp.symbols('s mu r')
a1,b1,e1,a2,b2,a3=sp.symbols('a1 b1 e1 a2 b2 a3')
N=5
rk=[1, a1+b1*s+e1*s**2, a2+b2*s, a3]
v=[sp.factorial(k)*s**k*rk[k] for k in range(4)]
Z=sum(v)
def binom(k): # binom(k+r,k) as poly in r
    return sp.expand(sp.prod([(r+j) for j in range(1,k+1)])/sp.factorial(k))
A=sp.series(sum(v[k]*(k+1) for k in range(4))/Z,s,0,N).removeO()
B=sp.series(sum(v[k]*binom(k) for k in range(4))/Z,s,0,N).removeO()
# log R_r = log B - r log A
logR=sp.series(sp.log(B)-r*sp.log(A),s,0,N).removeO()
# invert mu = s A(s)
muexpr=sp.expand(s*A)
# series reversion
c=sp.symbols('c1:6')
sguess=mu+sum(c[i]*mu**(i+2) for i in range(3))
eq=sp.series(muexpr.subs(s,sguess)-mu,mu,0,N).removeO()
sol=sp.solve([sp.expand(eq).coeff(mu,i) for i in range(2,N)],c[:3],dict=True)[0]
sm=sguess.subs(sol)
logRmu=sp.expand(sp.series(logR.subs(s,sm),mu,0,4).removeO())
for i in range(4):
    print(i, sp.factor(logRmu.coeff(mu,i)))
