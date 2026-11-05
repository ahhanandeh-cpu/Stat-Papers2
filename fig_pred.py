import pickle, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
out=pickle.load(open('pred.pkl','rb'))
sim=np.array([o[3] for o in out]); pr=np.array([o[4] for o in out]); ns=np.array([o[2] for o in out])
e=pr-sim; se=np.sqrt(sim*(1-sim)/1000)
print('MAE',np.abs(e).mean(),'max',np.abs(e).max(),'corr',np.corrcoef(sim,pr)[0,1],'within2se',np.mean(np.abs(e)<=2*np.sqrt(se**2+1e-6)))
for nn in sorted(set(ns)): print(nn, np.abs(e[ns==nn]).mean().round(4), np.abs(e[ns==nn]).max().round(3))
i=np.argmax(np.abs(e)); print(out[i])
fam=[o[0] for o in out]; mk={'PL':'o','PSh':'s','PAk':'^','PS':'D','PA':'v'}
fig,ax=plt.subplots(figsize=(4.6,4.4))
for f in mk:
    idx=[k for k in range(len(out)) if fam[k]==f]
    ax.scatter(pr[idx],sim[idx],marker=mk[f],s=22,facecolors='none',edgecolors='k',linewidths=0.8,label={'PL':'Poisson–Lindley','PSh':'Poisson–Shanker','PAk':'Poisson–Akash','PS':'Poisson–Sujatha','PA':'Poisson–Amarendra'}[f])
ax.plot([0,1],[0,1],'k-',lw=0.6); ax.axhline(1/7,ls=':',c='k',lw=0.8)
ax.set_xlabel('Predicted by the Gaussian limit'); ax.set_ylabel('Simulated (1,000 replicates)')
ax.set_xlim(0,1); ax.set_ylim(0,1); ax.legend(frameon=False,fontsize=7,loc='upper left')
fig.tight_layout(); fig.savefig('fig7_selection_limit.pdf'); fig.savefig('fig7.png',dpi=110)
