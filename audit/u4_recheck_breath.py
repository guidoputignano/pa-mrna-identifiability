"""
Independent re-check of audit/u3_joint_fisher.py, sharing no code path.

WHAT IS INDEPENDENT. u3 is a Fisher calculation: it linearises the flux
balance, builds an information matrix in (ln theta, ln K_m, b, ln k_e, ln c)
and inverts it. This script never forms an information matrix and never
linearises. It simulates records in the OBSERVATION space, with assay error on
the marker in concentration units and recording error on the load in flux
units, fits the five parameters by maximising the same posterior numerically,
and reports the spread of ln theta_hat across replicates. The two agree only
if both the linearisation and the algebra are right.

WHAT IS DELIBERATELY SHARED. The constants, and nothing else: theta, K_m, k_e,
L, the dose train, the prior widths and the gate, all read from
src/common/params.py and src/design/run.py. Re-typing those would be checking
a different model.

THE SIMULATION. For visit i at time t_i:
    recorded load    L_i, constant at L
    true load        (L + b_true) (1 + e_i),   e_i ~ N(0, load_frac^2)
    true marker      the positive root of the balance at (theta, K_m, k_e)
    observed marker  m_i (1 + a_i),            a_i ~ N(0, assay_cv^2)
The load error enters twice, once per visit through e_i and once as the
constant bias b_true, which is the same double entrance q6_sigma_budget.py
records in sigma and in the prior on b.

For breath reading j at the same visit:
    observed y_j = c_true F(m_j) exp(eps_j),   eps_j ~ N(0, breath_cv^2)
with F the oxidised fraction derived in audit/u2_breath_functional.py part A
and c_true an unknown recovery factor, fitted freely with no prior.

REPORTING. The robust spread, 0.7413 times the interquartile range of
ln theta_hat, alongside the plain standard deviation, because a heavy tail
from the nonlinearity inflates the second and this repository has been wrong
that way before. Collapse counts carry a Wilson interval. A replicate counts
as collapsed if the optimiser does not converge or if ln theta_hat leaves
+/- ln 100 of the truth.

FOOTING. As u3: recorded load held at L = 322.727 umol/L/d with the operating
level, sigma and the priors re-solved at each capacity; 25 monthly draws;
k breath readings alternating trough and peak of the dose train; k_e
log-prior 0.8; load-bias prior 0.05 L. Seeds are stated with every number.

Nothing here records an expected result in advance.

Run:  python audit/u4_recheck_breath.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import least_squares

from design.run import (P, delta, load, DOSE_INTERVAL, LAMBDA, D_AMP,
                        LOAD_FRAC, ASSAY_CV)

TH, KM, KE = P.theta, P.Km, P.ke
LB = float(load(0.0))
T = DOSE_INTERVAL
DEC = float(np.exp(-LAMBDA * T))
DPEAK = D_AMP / (1.0 - DEC)
DTROUGH = DPEAK * DEC
KE_PRIOR = 0.8
SB = LOAD_FRAC * LB


def root(theta, Km, ke, inc, Ltot):
    """Positive root of (theta+inc) m/(Km+m) + ke m = Ltot."""
    inc = np.asarray(inc, dtype=float)
    Ltot = np.asarray(Ltot, dtype=float)
    q = ke * Km + theta + inc - Ltot
    return (-q + np.sqrt(q * q + 4.0 * ke * Ltot * Km)) / (2.0 * ke)


def oxidised(theta, Km, ke, m):
    return theta / (theta + ke * (Km + m))


def simulate(theta, n, k, breath_cv, rng, load_frac=LOAD_FRAC,
             assay_cv=ASSAY_CV, c_true=0.62, draw_priors=True,
             ke_prior=KE_PRIOR):
    """
    One synthetic record: observed markers, observed breath, the design.

    EVERY PARAMETER CARRYING A PROPER PRIOR IS DRAWN FROM THAT PRIOR, and
    this is the whole correctness condition of the comparison rather than a
    refinement of it. The Fisher variance u3 reports is a POSTERIOR variance,
    and the spread of a MAP estimate over replicates equals it only when the
    generating distribution of the nuisances is the prior the fit assumes.
    Holding k_e at its prior mean instead, which is what this script did
    first, silently simulates the k_e-exactly-known design and makes the
    Fisher number look pessimistic by a factor of two at nominal capacity.
    `draw_priors=False` reproduces that error on purpose, for the diagnostic.

    theta itself is NOT drawn: it carries a flat prior and is the estimand,
    so the spread reported is conditional on it, which is what the gate is
    stated against.
    """
    t = 30.4 * np.arange(n, dtype=float)
    inc = np.asarray(delta(t), dtype=float)
    if draw_priors:
        b_true = rng.normal(0.0, SB)
        ke_true = KE * np.exp(rng.normal(0.0, ke_prior))
    else:
        b_true, ke_true = rng.normal(0.0, SB), KE
    e = rng.normal(0.0, load_frac, n)
    Ltrue = (LB + b_true) * (1.0 + e)
    m_true = root(theta, KM, ke_true, inc, Ltrue)
    m_obs = m_true * (1.0 + rng.normal(0.0, assay_cv, n))

    if k:
        binc = np.array([DTROUGH if j % 2 == 0 else DPEAK for j in range(k)],
                        dtype=float)
        be = rng.normal(0.0, load_frac, k)
        bL = (LB + b_true) * (1.0 + be)
        bm_true = root(theta, KM, ke_true, binc, bL)
        y = (c_true * oxidised(theta, KM, ke_true, bm_true)
             * np.exp(rng.normal(0.0, breath_cv, k)))
    else:
        binc, y = np.array([]), np.array([])
    return t, inc, m_obs, binc, y


def weights(theta_guess, inc, m_obs, assay_cv=ASSAY_CV, load_frac=LOAD_FRAC):
    """
    Per-visit sd of the observed level, fixed from the DATA and not re-solved
    inside the fit.

    Holding the weights fixed is deliberate and it is what makes the
    comparison with u3 an apples-to-apples one: u3's sigma is a constant of
    the design, not a function of the parameters, so a fit that re-solved its
    own weights would be testing a different estimator as well as a different
    code path. The tangent slope is taken at the observed level with the
    increment carried, which is the local gain of the flux to the level.
    """
    J = (theta_guess + inc) * KM / (KM + m_obs) ** 2 + KE
    return np.sqrt((assay_cv * m_obs) ** 2 + (load_frac * LB / J) ** 2)


def resid(z, inc, m_obs, sd_m, binc, y, breath_cv,
          ke_prior=KE_PRIOR):
    """
    Residual vector whose half-sum-of-squares is the negative log posterior.

    z = (ln theta, ln K_m, b/SB, ln k_e, ln c), every entry O(1), which is
    what lets a scaled trust-region solver work on an objective whose valley
    runs along the fibre.
    """
    lt, lK, bs, lk, lc = z
    theta, Km, ke, c = np.exp(lt), np.exp(lK), np.exp(lk), np.exp(lc)
    Ltot = LB + bs * SB
    if Ltot <= 0 or not np.isfinite(Ltot):
        return np.full(len(m_obs) + len(y) + 2, 1e6)
    m_mod = root(theta, Km, ke, inc, Ltot)
    if not np.all(np.isfinite(m_mod)) or np.any(m_mod <= 0):
        return np.full(len(m_obs) + len(y) + 2, 1e6)
    parts = [(m_obs - m_mod) / sd_m]
    if len(y):
        bm = root(theta, Km, ke, binc, Ltot)
        if not np.all(np.isfinite(bm)) or np.any(bm <= 0):
            return np.full(len(m_obs) + len(y) + 2, 1e6)
        pred = c * oxidised(theta, Km, ke, bm)
        parts.append((np.log(y) - np.log(pred)) / breath_cv)
    parts.append(np.array([bs]))                        # load-bias prior
    parts.append(np.array([(lk - np.log(KE)) / ke_prior]))   # k_e log-prior
    return np.concatenate(parts)


def fit(theta, inc, m_obs, binc, y, breath_cv, c_true=0.62,
        ke_prior=KE_PRIOR):
    """
    Maximum a posteriori by scaled trust-region least squares.

    Nelder-Mead was tried first and is recorded here because it is the
    failure this re-check would otherwise have reported as a finding: on a
    record whose objective has an exactly flat direction, the simplex stalls
    near its start, and starting it at the truth then UNDER-disperses the
    estimates and makes the Fisher number look pessimistic by a factor of
    two. The convergence check below is what catches that.
    """
    sd_m = weights(theta, inc, m_obs)
    z0 = np.array([np.log(theta), np.log(KM), 0.0, np.log(KE),
                   np.log(c_true)])
    if not len(y):
        z0 = z0[:4]

        def f(z):
            return resid(np.append(z, 0.0), inc, m_obs, sd_m, binc, y,
                         breath_cv, ke_prior)
    else:
        def f(z):
            return resid(z, inc, m_obs, sd_m, binc, y, breath_cv,
                         ke_prior)

    best, bestcost = None, np.inf
    for jitter in (0.0, 0.3, -0.3):
        zs = z0.copy()
        if jitter:
            zs = zs + jitter * np.array([1.0, -1.0, 0.5, 0.2, 0.1][:len(zs)])
        try:
            r = least_squares(f, zs, method="trf", x_scale=1.0,
                              xtol=1e-14, ftol=1e-14, gtol=1e-14,
                              max_nfev=20000)
        except Exception:
            continue
        if r.cost < bestcost:
            best, bestcost = r, r.cost
    if best is None:
        return np.nan, False
    #  converged means the scaled gradient is small at the optimum, which is
    #  the property a stalled simplex does not have
    grad = best.jac.T @ best.fun
    ok = bool(np.linalg.norm(grad) < 1e-4 * max(1.0, np.linalg.norm(best.fun)))
    return best.x[0], ok


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


def run(theta, k, breath_cv, reps, n=25, seed0=0, draw_priors=True,
        ke_prior=KE_PRIOR):
    ests, bad = [], 0
    for s in range(reps):
        rng = np.random.default_rng(10_000 * (seed0 + 1) + s)
        t, inc, m_obs, binc, y = simulate(theta, n, k, breath_cv, rng,
                                          draw_priors=draw_priors,
                                          ke_prior=ke_prior)
        lt, ok = fit(theta, inc, m_obs, binc, y, breath_cv,
                     ke_prior=ke_prior)
        if (not ok) or abs(lt - np.log(theta)) > np.log(100.0):
            bad += 1
            continue
        ests.append(lt - np.log(theta))
    e = np.array(ests, float)
    lo, hi = wilson(bad, reps)
    return e, bad, (lo, hi)


if __name__ == "__main__":
    import u3_joint_fisher as U3

    REPS = 300
    print("=" * 92)
    print("RE-CHECK of u3 by simulation and fitting, no information matrix")
    print("=" * 92)
    print()
    print("  25 monthly draws, constant recorded load, load recorded to 5%,")
    print("  marker assay CV 5%, breath readings alternating the dose train's")
    print("  trough and peak. Seeds 10000*(case+1) + 0..%d." % (REPS - 1))
    print()
    print("  THE CONTROL COLUMN IS THE POINT. The k_e log-prior the paper")
    print("  carries is 0.8, a factor of 2.2, which is far outside the range")
    print("  a linearisation is good over, so the two methods are not")
    print("  expected to agree closely at that width. Narrowing the prior to")
    print("  0.2 changes nothing about the algebra and everything about the")
    print("  linearity, so if the two agree there and not at 0.8, the gap is")
    print("  curvature and not an error.")
    print()

    def line(theta, k, cv, ke_prior, draw_priors, idx):
        mbar = float(U3.level(theta, 0.0, LB))
        sig = U3.sigma_at(theta, mbar)
        t = 30.4 * np.arange(25, dtype=float)
        m, inc = U3.record(theta, t, 0.0)
        incs, _ = U3.breath_design(k)
        mb = U3.level(theta, incs, LB) if k else np.array([])
        pred = U3.sd_joint(theta, m, inc, sig, mb, incs, cv,
                           ke_prior=ke_prior)
        e, bad, (lo, hi) = run(theta, k, cv, REPS, seed0=idx,
                               draw_priors=draw_priors, ke_prior=ke_prior)
        rs = robust_sd(e)
        print(f"  {theta:7.0f} {k:4d} {cv:5.0%} {ke_prior:8.2f} "
              f"{str(draw_priors):>6} {pred:10.4f} {rs:10.4f} "
              f"{e.std(ddof=1):9.4f} {e.mean():10.4f} {rs / pred:7.3f} "
              f"{bad:3d}/{REPS:<4d} [{lo:.3f},{hi:.3f}]")
        return rs / pred

    print(f"  {'theta':>7} {'k':>4} {'CV':>5} {'ke prior':>8} {'drawn':>6} "
          f"{'u3 Fisher':>10} {'robust sd':>10} {'plain sd':>9} "
          f"{'mean bias':>10} {'ratio':>7} {'collapse':>8} {'95% CI':>15}")
    print("  " + "-" * 116)

    idx = 0
    print("  -- the paper's own k_e prior, 0.8 --")
    for theta, k, cv in [(3000.0, 0, 0.10), (3000.0, 32, 0.05),
                         (600.0, 34, 0.05), (600.0, 64, 0.05),
                         (300.0, 32, 0.05), (300.0, 70, 0.05)]:
        line(theta, k, cv, 0.8, True, idx); idx += 1

    print("  -- the linearity control: k_e prior 0.2, everything else equal --")
    ctrl = []
    for theta, k, cv in [(3000.0, 0, 0.10), (3000.0, 32, 0.05),
                         (600.0, 34, 0.05), (600.0, 64, 0.05),
                         (300.0, 32, 0.05), (300.0, 70, 0.05)]:
        ctrl.append(line(theta, k, cv, 0.2, True, idx)); idx += 1

    print("  -- the diagnostic: k_e held at its prior mean instead of drawn --")
    for theta, k, cv in [(3000.0, 0, 0.10), (300.0, 70, 0.05)]:
        line(theta, k, cv, 0.8, False, idx); idx += 1

    print()
    print("  READING.")
    print("  1. The 'drawn = False' rows are the error this script made first")
    print("     and are kept as a row rather than deleted. Holding a")
    print("     prior-carrying nuisance at its prior mean simulates the")
    print("     nuisance-known design, and at nominal capacity that")
    print("     reproduces audit/u1_fibre_bypass.py's k_e-exact figure")
    print("     instead of the paper's.")
    print("  2. The control rows carry the verdict. Their ratios are")
    print(f"     {min(ctrl):.3f} to {max(ctrl):.3f}.")
    print("  3. The gate-relevant row is theta = 300 at k = 70, which is the")
    print("     smallest k u3 reports for that capacity at a 5% breath CV.")
