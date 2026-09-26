"""Re-check P3 through the repository's own likelihood (n6_profile_vs_crb.nll) and its own
dose train (design.run.delta), with only the GENERATOR changed. Shares no fitting code with fit.py."""
import sys
from pathlib import Path
#  the repository root, found from this file rather than hard-coded: the probe
#  lives at audit/probes_2026-09-25/, two levels below it.
R = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(R / "audit")); sys.path.insert(0, str(R / "src"))
import numpy as np
from scipy.optimize import minimize
import design.run as DR
import n6_profile_vs_crb as N6
from n4_theta_scope import footing, KM0, KE, TH0, MONTHLY
from common.models import ParamsA
import model_A.run as LOC

def gen(theta, seed, inc_true):
    L, m_bar, _, _ = footing(theta, convention="A")
    sigma = float(LOC.sigma_total(DR.LOAD_FRAC, DR.ASSAY_CV, p=ParamsA(theta=theta, Km=KM0, ke=KE), m_bar=m_bar))
    rng = np.random.default_rng(seed)
    s2 = np.log1p(0.01)
    loads = L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), MONTHLY.size))
    B = KE * KM0 + theta + inc_true - loads
    m = (-B + np.sqrt(B ** 2 + 4.0 * KE * loads * KM0)) / (2.0 * KE)
    return m, loads, sigma, L

inc_model = np.asarray(DR.delta(MONTHLY), float)
# 7-d truth, peak-matched: rebinding LAMBDA in the repo's own module, then restoring it
lam0, d0 = DR.LAMBDA, DR.D_AMP
DR.LAMBDA = np.log(2) / 7.0
DR.D_AMP = 6400.0 * (1 - np.exp(-DR.LAMBDA * 14.0))
inc_true = np.asarray(DR.delta(MONTHLY), float)
DR.LAMBDA, DR.D_AMP = lam0, d0
out = []
for seed in range(8):
    m, y, sigma, L = gen(TH0, seed, inc_true)
    best = None
    for st in ([np.log(TH0), np.log(KM0), 0.0, np.log(KE)], [np.log(TH0) + 1, np.log(KM0) - 1, 0, np.log(KE)],
               [np.log(TH0) - 1, np.log(KM0) + 1, 0, np.log(KE)]):
        r = minimize(N6.nll, st, args=(m, inc_model, y, sigma, L), method="Nelder-Mead",
                     options=dict(maxiter=40000, xatol=1e-10, fatol=1e-12))
        if best is None or r.fun < best.fun: best = r
    out.append(np.exp(best.x[0]) / TH0)
print("repo likelihood, 7 d peak-matched truth, analyst 3.5 d, monthly, noiseless: theta_hat/theta median",
      round(float(np.median(out)), 3), "range", round(min(out), 3), round(max(out), 3))
