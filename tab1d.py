import numpy as np
from scipy.stats import norm
KAPG={-2:'PL',0:'geometric, PSh',2:'PS, PA',4:'PAk'}
ks=sorted(KAPG); r2=np.sqrt(2)
def probs(k1,t):
    H=[t*(k-k1)/2 for k in ks]; out={}
    for i,h in enumerate(H):
        a=-np.inf if i==0 else (H[i-1]+h)/2; b=np.inf if i==len(H)-1 else (h+H[i+1])/2
        lo=max(a,h-r2); hi=min(b,h+r2); out[ks[i]]=max(0,norm.cdf(hi)-norm.cdf(lo))
    out['NB']=1-sum(out.values()); return out
rows=[]
for tr in ks:
    for t in [0,0.5,1,2,4]:
        p=probs(tr,t); rows.append((tr,t,p)); print(KAPG[tr],t,{KAPG.get(k,k):round(v,3) for k,v in p.items()})
# latex
L=[r'\begin{tabular}{llrrrrr}',r'\toprule',r'True family ($\kappa$) & $t=\sqrt n\,\mu^3$ & PL & Geo/PSh & PS/PA & PAk & NB\\',r' & & ($-2$) & ($0$) & ($2$) & ($4$) & \\',r'\midrule']
name={-2:'Poisson--Lindley ($-2$)',0:'Poisson--Shanker ($0$)',2:'Poisson--Sujatha ($2$)',4:'Poisson--Akash ($4$)'}
for tr in ks:
    first=True
    for (a,t,p) in [r for r in rows if r[0]==tr]:
        cells=[('\\textbf{%.3f}'%p[k]) if k==tr else '%.3f'%p[k] for k in ks]+['%.3f'%p['NB']]
        L.append((name[tr] if first else '')+' & %g & '%t+' & '.join(cells)+r'\\'); first=False
    if tr!=ks[-1]: L.append(r'\addlinespace')
L+= [r'\bottomrule',r'\end{tabular}']
open('tab_selection_limit.tex','w').write('\n'.join(L)+'\n')
