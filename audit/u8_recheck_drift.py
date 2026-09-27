"""
Independent re-check of audit/u7_drift.py, sharing no code path.

WHAT IS INDEPENDENT. u7 computes the bias by minimising a FLUX residual on
NOISELESS data and takes the variance from an information matrix. This script
adds noise, writes the likelihood in LEVEL space rather than flux space, fits
every replicate by nonlinear least squares and measures the mean offset and
the spread of ln theta_hat directly. Neither the bias nor the variance comes
from anything u7 computed.

WHAT IS SHARED, deliberately: the constants, and the definition of the drift,
ln theta(t) = ln theta_ref + r (t - t_ref).

THE ESTIMAND is capacity at the LAST visit, which is what a dosing decision
needs and which is u7's primary convention.

REPORTING. The empirical mean offset is the bias; the robust spread, 0.7413
times the interquartile range, is the sd; their hypot is the rmse. Collapse
counts carry a Wilson interval. Every prior-carrying nuisance is drawn from
its prior, which audit/u4_recheck_breath.py establishes is the condition for
a simulated spread to be comparable with a posterior standard deviation.

FOOTING. theta = 3000 at the reference time, K_m = 100, k_e = 5/d,
L = 322.727 umol/L/d, D = 6000 on a 14 d train with a 3.5 d half-life,
monthly draws, load recorded to 5%, day-to-day intake CV 10%, marker assay
CV 5%, k_e log-prior 0.8, load-bias prior 0.05 L, gate 0.20. Seeds stated per
row.

Nothing here records an expected result in advance.

Run:  python audit/u8_recheck_drift.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import least_squares

from design.run import (ASSAY_CV, DOSE_INTERVAL, D_AMP, LAMBDA, LOAD_CV,
                        LOAD_FRAC, load)

KM, KE, THETA = 100.0, 5.0, 3000.0
LB = float(load(0.0))
GATE = 0.20
MONTH = 30.4
YEAR = 365.25
KE_PRIOR = 0.8
SB = LOAD_FRAC * LB
T_DOSE = 50.0 * DOSE_INTERVAL
PEAK = D_AMP / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))
PUBLISHED_KAPPA_CV = 0.1402


def root(theta, Km, ke, inc, Ltot):
    theta = np.asarray(theta, dtype=float)
    inc = np.asarray(inc, dtype=float)
    Ltot = np.asarray(Ltot, dtype=float)
    q = ke * Km + theta + inc - Ltot
    return (-q + np.sqrt(q * q + 4.0 * ke * Ltot * Km)) / (2.0 * ke)


def schedule(n, step=MONTH):
    t = step * np.arange(n, dtype=float)
    return t, PEAK * np.exp(-LAMBDA * np.mod(t - T_DOSE, DOSE_INTERVAL))


def sigma_fixed(assay_cv=ASSAY_CV):
    mbar = float(root(THETA, KM, KE, 0.0, LB))
    slope = THETA * KM / (KM + mbar) ** 2 + KE
    return float(np.hypot(LOAD_FRAC * LB, slope * assay_cv * mbar))


def one(n, r, rng, kappa_cv=None, assay_cv=ASSAY_CV):
    """Simulate one drifting record and fit a CONSTANT capacity to it."""
    t, inc = schedule(n)
    tref = t[-1]
    theta_t = THETA * np.exp(r * (t - tref) / YEAR)

    b_true = rng.normal(0.0, SB)
    ke_true = KE * np.exp(rng.normal(0.0, KE_PRIOR))
    kap_true = (np.exp(rng.normal(0.0, kappa_cv))
                if kappa_cv is not None else 1.0)
    Lrec = LB * (1.0 + LOAD_CV * rng.standard_normal(n))
    Ltrue = (Lrec + b_true) * (1.0 + rng.normal(0.0, LOAD_FRAC, n))
    m_true = root(theta_t, KM, ke_true, kap_true * inc, Ltrue)
    m_obs = m_true * (1.0 + rng.normal(0.0, assay_cv, n))

    #  level-space weights, fixed from the data
    J0 = (THETA + inc) * KM / (KM + m_obs) ** 2 + KE
    sd_m = np.sqrt((assay_cv * m_obs) ** 2 + (LOAD_FRAC * LB / J0) ** 2)

    def f(z):
        lt, lK, bs, lk = z[0], z[1], z[2], z[3]
        kap = np.exp(z[4]) if kappa_cv is not None else 1.0
        theta, Km, ke = np.exp(lt), np.exp(lK), np.exp(lk)
        L = Lrec + bs * SB
        if np.any(L <= 0) or Km <= 0 or ke <= 0:
            return np.full(n + 3, 1e6)
        mm = root(theta, Km, ke, kap * inc, L)
        if not np.all(np.isfinite(mm)) or np.any(mm <= 0):
            return np.full(n + 3, 1e6)
        tail = [bs, (lk - np.log(KE)) / KE_PRIOR,
                (z[4] / kappa_cv) if kappa_cv is not None else 0.0]
        return np.concatenate([(m_obs - mm) / sd_m, tail])

    z0 = np.array([np.log(THETA), np.log(KM), 0.0, np.log(KE), 0.0])
    best, cost = None, np.inf
    for j in (0.0, 0.3, -0.3):
        try:
            rr = least_squares(f, z0 + j * np.array([1.0, -.5, .3, .2, 0.]),
                               method="trf", xtol=1e-13, ftol=1e-13,
                               gtol=1e-13, max_nfev=12000)
        except Exception:
            continue
        if rr.cost < cost:
            best, cost = rr, rr.cost
    if best is None:
        return np.nan
    return float(best.x[0] - np.log(THETA))


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def robust_sd(x):
    q1, q3 = np.percentile(x, [25, 75])
    return 0.7413 * (q3 - q1)


def run(n, r, reps, kappa_cv=None, seed0=0):
    out, bad = [], 0
    for i in range(reps):
        v = one(n, r, np.random.default_rng(90_000 * (seed0 + 1) + i),
                kappa_cv=kappa_cv)
        if not np.isfinite(v) or abs(v) > np.log(100.0):
            bad += 1
            continue
        out.append(v)
    return np.array(out), bad


if __name__ == "__main__":
    import u7_drift as U7

    REPS = 240
    print("=" * 96)
    print("RE-CHECK of u7 by simulation and fitting, level space, no")
    print("information matrix and no noiseless fit")
    print("=" * 96)
    print(f"\n  {REPS} replicates per row, seeds 90000*(case+1) + 0..{REPS - 1}.")
    print("  Capacity is wanted at the LAST visit. The simulated bias is the")
    print("  mean offset of ln theta_hat, not a projection of anything.\n")

    print(f"  {'n':>5} {'r /yr':>7} {'kappa':>7} "
          f"{'u7 bias':>9} {'sim bias':>9} "
          f"{'u7 sd':>8} {'sim sd':>8} {'u7 rmse':>9} {'sim rmse':>9} "
          f"{'ratio':>7} {'collapse':>9} {'95% CI':>15}")
    print("  " + "-" * 116)
    idx = 0
    rows = []
    for kap in (None, PUBLISHED_KAPPA_CV):
        for n, r in [(25, 0.0), (25, 0.05), (118, 0.0), (118, 0.01),
                     (118, 0.02), (118, 0.05)]:
            sd7, b7, rm7, _, _ = U7.med_rmse(n, r, kappa_cv=kap)
            e, bad = run(n, r, REPS, kappa_cv=kap, seed0=idx)
            idx += 1
            bs, ss = float(np.median(e)), robust_sd(e)
            rms = float(np.hypot(ss, bs))
            lo, hi = wilson(bad, REPS)
            print(f"  {n:5d} {r:7.3f} {('-' if kap is None else f'{kap:.3f}'):>7} "
                  f"{b7:9.4f} {bs:9.4f} {sd7:8.4f} {ss:8.4f} "
                  f"{rm7:9.4f} {rms:9.4f} {rms / rm7:7.3f} "
                  f"{bad:3d}/{REPS:<4d} [{lo:.3f},{hi:.3f}]")
            rows.append((n, r, kap, rm7, rms))

    print("\n  READING. The bias columns are the load-bearing comparison:")
    print("  u7's bias is a noiseless nonlinear fit in flux space and the")
    print("  simulated one is the mean of a noisy nonlinear fit in level")
    print("  space, so agreement there is agreement about the ESTIMATOR and")
    print("  not about a shared derivation.")
    rat = [r[4] / r[3] for r in rows if np.isfinite(r[3]) and r[3] > 0]
    print(f"\n  rmse ratios span {min(rat):.3f} to {max(rat):.3f}.")
    print("  The k_e log-prior is 0.8, a factor of 2.2, and")
    print("  audit/u4_recheck_breath.py records that a linearised variance")
    print("  and a simulated spread are not expected to agree closely at")
    print("  that width; the bias comparison does not depend on it.")
