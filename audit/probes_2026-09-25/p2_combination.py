"""At 5% of nominal: do routes that each fail alone reach 0.20 together?"""
import numpy as np
from core import *

sched = {"monthly (25)": np.arange(0.0, 730.0, 30.4),
         "weekly+1d (104)": np.arange(0.0, 730.0, 7.0) + 0.5 + (np.arange(0,730,7.0)//14)*0.0,
         "every 2 d (365)": np.arange(1.0, 730.0, 2.0)}
# weekly draws land on only 2 phases of a 14 d cycle; use 7.3 d so phases spread
sched["weekly+1d (104)"] = np.arange(0.5, 730.0, 7.3)
TH = 150.0
seeds = range(12)
print("theta = 150 (5% of nominal). cells: median sd(ln theta) over 12 load records")
print("k_e prior sd (log) | " + " | ".join(f"{k:>18s}" for k in sched))
for kp in (0.8, 0.2, 0.1, 0.05, 0.02, 0.0):
    for bmode in ("b prior 5% L", "b known"):
        bp = None if bmode.startswith("b prior") else np.inf
        cells = [med_sd(TH, t, seeds, ke_prior=kp, b_prior=bp)[0] for t in sched.values()]
        lab = "exact" if kp == 0 else f"{kp:.2f} (x{np.exp(kp):.2f})"
        print(f"{lab:>16s} {bmode:>12s} | " + " | ".join(f"{c:18.3f}" for c in cells))
print()
print("Same, with the assay CV halved to 2.5% (sigma_scale approximated by recomputing)")
import core
core.ASSAY_CV = 0.025
for kp in (0.1, 0.05, 0.0):
    cells = [med_sd(TH, t, seeds, ke_prior=kp)[0] for t in sched.values()]
    print(f"k_e prior {kp:5.2f}                | " + " | ".join(f"{c:18.3f}" for c in cells))
core.ASSAY_CV = 0.05
print()
print("What k_e precision reaches 0.20 at 5%, per schedule (b prior as published):")
from scipy.optimize import brentq
for name, t in sched.items():
    f = lambda lk: med_sd(TH, t, seeds, ke_prior=np.exp(lk))[0] - 0.20
    try:
        r = np.exp(brentq(f, np.log(1e-4), np.log(0.8)))
        print(f"  {name:>18s}: ln k_e sd <= {r:.4f}  (k_e known to about {100*r:.1f}%)")
    except ValueError:
        print(f"  {name:>18s}: not reached even with k_e exact")
print()
print("Same question at 10% of nominal (theta=300):")
for name, t in sched.items():
    f = lambda lk: med_sd(300.0, t, seeds, ke_prior=np.exp(lk))[0] - 0.20
    try:
        r = np.exp(brentq(f, np.log(1e-4), np.log(0.8)))
        print(f"  {name:>18s}: ln k_e sd <= {r:.4f}  (k_e known to about {100*r:.1f}%)")
    except ValueError:
        print(f"  {name:>18s}: not reached even with k_e exact")
