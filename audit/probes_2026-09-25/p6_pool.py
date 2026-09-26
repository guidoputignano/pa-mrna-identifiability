"""
A slow peripheral pool behind the plasma marker (Phase 6 section 3.4 names it the
weakest untested assumption). The steady map is unchanged by construction
(k12 m = k21 p at equilibrium), so every steady-state check passes; only the
dynamics change. Daily loads (10% CV), the single-exponential dose train, and a
draw at a fixed hour on each visit day. The analyst fits the paper's model.
"""
import numpy as np
from scipy.integrate import solve_ivp
from core import *
from fit import fit

def simulate(times, seed, k12, k21, theta=TH_NOM):
    days = int(np.ceil(times.max())) + 2
    rng = np.random.default_rng(10_000 + seed)
    s2 = np.log1p(LOAD_CV ** 2)
    Ld = L_A * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), days))
    m0 = level(theta, KM, KE, Ld[0], 0.0)
    y = np.array([m0, k12 / k21 * m0 if k21 > 0 else 0.0])
    out_m, out_L = [], []
    tq = sorted(times)
    j = 0
    def rhs(t, z, L):
        m, p = z
        E = theta + train(np.array([t]))[0]
        return [L - E * m / (KM + m) - KE * m - k12 * m + k21 * p, k12 * m - k21 * p]
    for d in range(days - 1):
        # break at dose times inside the day as well
        edges = [d, d + 1]
        dose = [k * TAU for k in range(int(d // TAU), int((d + 1) // TAU) + 1) if d < k * TAU < d + 1]
        pts = sorted(set(edges + dose))
        for a, b in zip(pts[:-1], pts[1:]):
            evals = [t for t in tq[j:] if a <= t < b]
            sol = solve_ivp(rhs, (a, b), y, args=(Ld[d],), method="LSODA", rtol=1e-9, atol=1e-11,
                            dense_output=bool(evals))
            for t in evals:
                out_m.append(sol.sol(t)[0]); out_L.append(Ld[d]); j += 1
            y = sol.y[:, -1]   # state at the END of the interval, not at the last draw
        if j >= len(tq):
            break
    return np.array(out_m), np.array(out_L)

times = np.round(np.arange(0.0, 730.0, 30.4)) + 0.4   # monthly, whole days, drawn ~10 am
sig, _ = sigma_for(TH_NOM)
print("monthly, noiseless, paper's priors. Pool half-life ln2/k21 = 3 d.")
print(f"{'k12 (/d)':>9} {'pool/plasma':>12} {'theta_hat/theta':>16} {'rms/sigma':>10}")
for k12 in (0.0, 1.0, 3.0, 10.0):
    k21 = np.log(2) / 3.0
    th, rr = [], []
    for s in range(4):
        m, L = simulate(times, s, k12, k21)
        r = fit(m, train(times), L, sig, L_A)
        th.append(r["theta"] / TH_NOM); rr.append(r["rms"])
    print(f"{k12:9.1f} {k12/k21:12.1f} {np.median(th):16.2f} {np.median(rr):10.2f}")
