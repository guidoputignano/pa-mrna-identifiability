import sys, numpy as np
sys.path.insert(0,'src')
import design.run as D, model_A.run as LOC
P=D.P; sigma=LOC.sigma_total(D.LOAD_FRAC,D.ASSAY_CV); Lb=D.load(0.0)

def filt(lv, dd, mode, kap_prior=None, ke_prior=None):
    """mode: 'pub' (no Delta), 'del' (Delta known), 'kap' (Delta scale free)."""
    names=['ln th','ln Km','b']; pri=[0.8**2,0.8**2,(D.LOAD_FRAC*Lb)**2]
    if ke_prior is not None: names.append('ln ke'); pri.append(ke_prior**2)
    if mode=='kap': names.append('ln kap'); pri.append((kap_prior if kap_prior else 1e6)**2)
    Finv=np.linalg.inv(np.diag(pri))
    for m,dd_i in zip(np.atleast_1d(lv),np.atleast_1d(dd)):
        g=m/(P.Km+m); dg=-P.Km*m/(P.Km+m)**2
        eff=P.theta+(0.0 if mode=='pub' else dd_i)
        q=[P.theta*g, eff*dg, 1.0]
        if ke_prior is not None: q.append(P.ke*m)
        if mode=='kap': q.append(dd_i*g)
        q=np.array(q); Finv=Finv+np.outer(q,q)/sigma**2
    return np.linalg.inv(Finv)

T=730.0
t=np.arange(0.0,T,2.0); lv=D.quasi_steady_level(t); dd=np.asarray(D.delta(t),float)
print("DENSE two-daily x 2 years, n=365")
print(f"  published (Delta dropped from the Km column)      sd(ln th) = {np.sqrt(filt(lv,dd,'pub')[0,0]):.4f}")
print(f"  Delta known exactly                               sd(ln th) = {np.sqrt(filt(lv,dd,'del')[0,0]):.4f}")
for kp in (0.05,0.10,0.20,0.30,0.50,1.0,None):
    lab = f"{kp:.0%}" if kp else "flat (unknown)"
    print(f"  Delta scale free, prior {lab:>14}         sd(ln th) = {np.sqrt(filt(lv,dd,'kap',kap_prior=kp)[0,0]):.4f}")
print()
print("  same, with k_e ALSO fitted (0.8 log-prior, the A1 recommended footing):")
for mode,kp,lab in [('pub',None,'published'),('del',None,'Delta known'),
                    ('kap',0.20,'Delta scale free, 20% prior'),('kap',None,'Delta scale free, flat')]:
    print(f"    {lab:>32}  sd(ln th) = {np.sqrt(filt(lv,dd,mode,kap_prior=kp,ke_prior=0.8)[0,0]):.4f}")
print()
t7=np.arange(0.0,D.DOSE_INTERVAL,2.0); lv7=D.quasi_steady_level(t7); dd7=np.asarray(D.delta(t7),float)
print("ONE DOSE CYCLE, n=7")
for mode,kp,lab in [('pub',None,'published'),('del',None,'Delta known'),
                    ('kap',0.20,'Delta scale free, 20% prior'),('kap',None,'Delta scale free, flat')]:
    print(f"    {lab:>32}  sd(ln th) = {np.sqrt(filt(lv7,dd7,mode,kap_prior=kp)[0,0]):.4f}")
