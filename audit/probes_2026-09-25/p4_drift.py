"""
Things the model holds constant for two years that a real patient may not.
Noiseless records at nominal capacity, paper's priors; the analyst's model is
the paper's (constant theta, constant b, constant increment amplitude).
"""
import numpy as np
from core import *
from fit import fit, record

sched = {"monthly (25)": np.arange(0.0, 730.0, 30.4),
         "14.5 d offset (51)": np.arange(0.5, 730.0, 14.5)}
sd_nom = {k: med_sd(TH_NOM, t, range(12))[0] for k, t in sched.items()}
print("headline CRB sd(ln theta):", {k: round(v, 3) for k, v in sd_nom.items()})

def run(label, make, kcv=None, ref=TH_NOM):
    cells = []
    for name, t in sched.items():
        th, rms = [], []
        for s in range(8):
            kw = make(t, s)
            m, y, sig, Ls = record(TH_NOM, t, s, noisy=False, **kw)
            r = fit(m, train(t), y, sig, L_A, kappa_cv=kcv)
            th.append(r["theta"] / (ref if np.isscalar(ref) else ref(t)))
            rms.append(r["rms"])
        b = np.log(np.median(th)) / sd_nom[name]
        cells.append(f"{np.median(th):5.2f} [{np.median(rms):4.2f}] ({b:+5.1f} sd)")
    print(f"  {label:44s} " + "  ".join(f"{c:>26s}" for c in cells))

print("\ntheta_hat / reference theta   [rms residual / sigma]   (bias in units of the headline sd)")
print(f"  {'':44s} " + "  ".join(f"{k:>26s}" for k in sched))
run("control", lambda t, s: {})
for r in (0.2, 0.4):
    run(f"increment amplitude falls {int(r*100)}% over 2 y (ADA)",
        lambda t, s, r=r: {"inc_true": train(t) * (1 - r * t / 730.0)})
    run(f"  same, amplitude free with 12% prior",
        lambda t, s, r=r: {"inc_true": train(t) * (1 - r * t / 730.0)}, kcv=0.12)
for r in (0.3, -0.3):
    run(f"theta drifts {r:+.0%} over 2 y (vs its mean)",
        lambda t, s, r=r: {"theta_t": TH_NOM * (1 + r * t / 730.0)},
        ref=lambda t, r=r: np.mean(TH_NOM * (1 + r * t / 730.0)))
def rw(t, s, step):
    rng = np.random.default_rng(99 + s)
    return np.cumsum(rng.normal(0, step * L_A, t.size))
for step in (0.02, 0.04):
    run(f"unrecorded load: random walk, {int(step*100)}% of L per visit",
        lambda t, s, step=step: {"b_true": rw(t, s, step)})
run("unrecorded load: +20% of L step at month 12",
    lambda t, s: {"b_true": np.where(t > 365, 0.20 * L_A, 0.0)})
run("unrecorded load: -40% of a 50% non-dietary share, alternate months",
    lambda t, s: {"b_true": np.where((t // 30.4) % 2 == 1, -0.2 * L_A, 0.0)})
