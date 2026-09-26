"""Monte Carlo check of the combination result at 5% of nominal: is the CRB honest there?"""
import numpy as np
from core import *
from fit import fit, record
TH = 150.0
cases = [("2-daily, k_e 2%, b known", np.arange(1.0, 730.0, 2.0), 0.02, 0.0),
         ("weekly+, k_e 2%, b known", np.arange(0.5, 730.0, 7.3), 0.02, 0.0),
         ("2-daily, k_e 2%, b prior 5%", np.arange(1.0, 730.0, 2.0), 0.02, None),
         ("monthly, k_e 0.8 (paper)", np.arange(0.0, 730.0, 30.4), 0.8, None)]
for name, t, kp, bsd in cases:
    crb = med_sd(TH, t, range(12), ke_prior=kp, b_prior=(np.inf if bsd == 0 else None))[0]
    est = []
    for s in range(60):
        m, y, sig, Ls = record(TH, t, s)
        # truth for ke prior centre: draw the prior centre from the prior (honest nuisance)
        rng = np.random.default_rng(5000 + s)
        ke_centre = KE * np.exp(rng.normal(0, kp))
        r = fit(m, train(t), y, sig, L_A, ke_prior=kp, b_sd=bsd, theta0=TH, ke0=ke_centre)
        est.append(np.log(r["theta"]) - np.log(TH))
    est = np.array(est)
    q = np.percentile(est, [2.5, 50, 97.5])
    print(f"{name:32s} CRB {crb:6.3f} | MC sd {est.std():6.3f}  RMSE {np.sqrt(np.mean(est**2)):6.3f}"
          f"  median bias {q[1]:+.3f}  95% range exp -> [{np.exp(q[0])*TH:.0f}, {np.exp(q[2])*TH:.0f}]")
