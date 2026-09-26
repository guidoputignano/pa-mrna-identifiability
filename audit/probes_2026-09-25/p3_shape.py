"""
The increment is treated as a KNOWN regressor, but its time course is modelled
(single exponential, 3.5 d). What does the estimate do when the true time
course has a different shape and the analyst's regressor does not?
Noiseless records -> the numbers are structural bias only.
"""
import numpy as np
from core import *
from fit import fit, record

TH = TH_NOM
sched = {"monthly (25)": np.arange(0.0, 730.0, 30.4),
         "14.5 d offset (51)": np.arange(0.5, 730.0, 14.5),
         "2-daily (365)": np.arange(1.0, 730.0, 2.0)}
fine = np.arange(0.0, 730.0, 0.01)
ss = (fine > 600) & (fine < 614)
a_pk = train(fine)[ss].max()
a_tr = train(fine)[ss].min()
a_mean = train(fine)[ss].mean()
bi = train_biexp(fine)
b_pk, b_tr, b_mean = bi[ss].max(), bi[ss].min(), bi[ss].mean()
e7 = train(fine, half=7.0)
print(f"analyst's train (3.5 d): steady-state peak {a_pk:.0f}, trough {a_tr:.0f} ({a_tr/a_pk:.3f} of peak), cycle mean {a_mean:.0f}")
print(f"biexponential truth, trough/peak {b_tr/b_pk:.3f}; 7 d single exp trough/peak {e7[ss].min()/e7[ss].max():.3f}")
print()

truths = {}
# (i) biexponential, peak matched; (ii) biexp, cycle-mean matched; (iii) 7 d, peak matched
truths["biexp, peak-matched"] = lambda t: train_biexp(t) * (a_pk / b_pk)
truths["biexp, mean-matched"] = lambda t: train_biexp(t) * (a_mean / b_mean)
truths["7 d half-life, peak-matched"] = lambda t: train(t, half=7.0) * (a_pk / e7[ss].max())
truths["5 d half-life, peak-matched"] = lambda t: train(t, half=5.0) * (a_pk / train(fine, half=5.0)[ss].max())
truths["control: truth = analyst"] = lambda t: train(t)

print("theta_hat/theta on NOISELESS records, paper's priors (k_e 0.8, b 5%), load CV 10%;")
print("median over 8 load records; [rms residual / sigma]. A correct model gives 1.00 [0.00].")
for kcv, klab in ((None, "increment amplitude fixed (the paper's headline footing)"),
                  (0.12, "amplitude free with a 12% prior (Tier 2 characterisation)")):
    print(f"\n  {klab}")
    print(f"  {'true time course':30s} " + " ".join(f"{k:>24s}" for k in sched))
    for name, f in truths.items():
        cells = []
        for t in sched.values():
            r_th, r_rms = [], []
            for s in range(8):
                m, y, sig, Ls = record(TH, t, s, inc_true=f(t), noisy=False)
                r = fit(m, train(t), y, sig, L_A, kappa_cv=kcv)
                r_th.append(r["theta"] / TH); r_rms.append(r["rms"])
            cells.append(f"{np.median(r_th):7.2f} [{np.median(r_rms):4.2f}]  ")
        print(f"  {name:30s} " + " ".join(f"{c:>24s}" for c in cells))
