"""
Shared machinery for the October 2026 brief on practical identifiability.

Everything here is built on the repository's own footing, imported from
audit/n4_theta_scope.py rather than restated, so that the numbers this
directory produces are comparable with Section 11 of the paper:

    convention A   the recorded load L stays at the theta = 3000 patient's
                   322.73 and the marker rises as capacity falls
    record         25 monthly draws over two years (t = 0, 30.4, ... 729.6 d)
    intake         day-to-day load variation, lognormal CV 10%, RECORDED
    therapy        the booked dose train, D = 6000 every 14 d, 3.5 d half-life,
                   carried as a KNOWN increment Delta(t) on the enzyme channel
    noise          sigma_total(): 5% assay CV propagated to the flux ordinate at
                   the untreated level, in quadrature with a 5% load-record
                   scatter
    priors         ln k_e ~ N(ln 5, 0.8^2); load bias b ~ N(0, (0.05 L)^2);
                   nothing on ln K_m or ln theta
    "24 records"   SEEDS = range(24): 24 draws of the intake sequence for ONE
                   simulated patient. Not 24 patients.

Coordinates are the repository's, (ln theta, ln K_m, b, ln k_e), so that an
information matrix here can be compared entry by entry with n4's.

Two noise models are carried:

    "repo"     n4's: one flux-space sigma for every draw, evaluated at the
               untreated operating level m_bar
    "perdraw"  the same two error sources, propagated at each draw's own level
               m_i and slope C'(m_i), which includes the increment. This is the
               brief's "Gaussian noise on ln m" read literally: a 5% CV on
               every marker measurement, whatever the marker's level.

The brief's coordinates (theta, k_e, K_m, L) and the repository's differ only
in the load: the brief's L is the repository's recorded load less the bias b.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "audit"))
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

import n4_theta_scope as N4
from n4_theta_scope import (KM0, KE, TH0, D0, GATE, SEEDS, MONTHLY, DENSE,
                            GRID, level, delta_at)
from design.run import LOAD_FRAC, ASSAY_CV, LAMBDA, DOSE_INTERVAL
from common.models import ParamsA
import model_A.run as LOC

COORDS = ("ln theta", "ln K_m", "b", "ln k_e")
I_TH, I_KM, I_B, I_KE = 0, 1, 2, 3
KE_PRIOR = 0.8
TH_5 = 0.05 * TH0
FIG = HERE / "fig"
FIG.mkdir(exist_ok=True)


# --------------------------------------------------------------- the record
def record(theta, D=D0, Km=KM0, ke=KE, times=MONTHLY, seed=0, load_cv=0.10):
    """One simulated intake record on n4's footing, convention A."""
    L, m_bar = N4.anchor(theta, Km, ke, "A")
    times = np.asarray(times, float)
    if load_cv:
        rng = np.random.default_rng(seed)
        s2 = np.log1p(load_cv ** 2)
        loads = L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
    else:
        loads = np.full(times.size, L)
    inc = delta_at(times, D) if D else np.zeros(times.size)
    m = level(theta, Km, ke, loads, inc)
    return dict(theta=theta, Km=Km, ke=ke, L=L, m_bar=m_bar, loads=loads,
                inc=inc, m=m, times=times)


def share(theta, Km=KM0, ke=KE, inc=0.0, m=None, L=None):
    """The enzyme's own share s = theta g(m) / L at the untreated balance."""
    if m is None:
        L, m = N4.anchor(theta, Km, ke, "A")
    return theta * m / (Km + m) / L


def eps_of(m, Km=KM0):
    """eps = K_m/(K_m + m) = d ln g / d ln m."""
    return Km / (Km + m)


# ----------------------------------------------------------- sensitivities
def flux_rows(r):
    """d(clearance)/d(ln theta, ln K_m, b, ln k_e) at each draw, n4's columns."""
    th, K, ke, m, inc = r["theta"], r["Km"], r["ke"], r["m"], r["inc"]
    return np.stack([th * m / (K + m),
                     -(th + inc) * K * m / (K + m) ** 2,
                     np.ones_like(m),
                     ke * m], axis=1)


def slope(r):
    """C'(m_i) = (theta + Delta_i) g'(m_i) + k_e, the clearance slope per draw."""
    K = r["Km"]
    return (r["theta"] + r["inc"]) * K / (K + r["m"]) ** 2 + r["ke"]


def sigmas(r, noise="repo", load_frac=LOAD_FRAC, assay_cv=ASSAY_CV):
    """Flux-space noise per draw."""
    if noise == "repo":
        p = ParamsA(theta=r["theta"], Km=r["Km"], ke=r["ke"])
        s = float(LOC.sigma_total(load_frac, assay_cv, p=p, m_bar=r["m_bar"]))
        return np.full(r["m"].size, s)
    if noise == "perdraw":
        return np.hypot(load_frac * r["L"], slope(r) * assay_cv * r["m"])
    raise ValueError(noise)


def fim(r, noise="repo", ke_prior=KE_PRIOR, b_prior=None, km_prior=None,
        th_prior=None, load_frac=LOAD_FRAC, assay_cv=ASSAY_CV, extra=()):
    """
    Information matrix in (ln theta, ln K_m, b, ln k_e).

    b_prior is the sd of the load-bias prior in flux units; None means n4's
    LOAD_FRAC * L. A prior of None on ln K_m or ln theta means none at all,
    which is the repository's footing. `extra` is a sequence of (row, sd)
    pairs, each an additional scalar observation with sensitivity `row`.
    """
    Q = flux_rows(r)
    sg = sigmas(r, noise, load_frac, assay_cv)
    F = (Q / sg[:, None] ** 2).T @ Q
    sb = load_frac * r["L"] if b_prior is None else b_prior
    pri = np.zeros(4)
    if th_prior:
        pri[I_TH] = th_prior ** -2
    if km_prior:
        pri[I_KM] = km_prior ** -2
    if sb is not np.inf and sb:
        pri[I_B] = sb ** -2
    if ke_prior:
        pri[I_KE] = ke_prior ** -2
    F = F + np.diag(pri)
    for row, sd in extra:
        row = np.asarray(row, float)
        F = F + np.outer(row, row) / sd ** 2
    return F


def cov(F):
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return None
    return np.linalg.inv(F)


def sd_th(F):
    P = cov(F)
    return np.inf if P is None else float(np.sqrt(P[0, 0]))


def corr(P, i, j):
    return float(P[i, j] / np.sqrt(P[i, i] * P[j, j]))


def over_seeds(fn, seeds=SEEDS):
    """(median, min, max) of fn(seed) over the seed sweep, infinities dropped."""
    v = np.array([fn(s) for s in seeds], float)
    v = v[np.isfinite(v)]
    if not v.size:
        return np.inf, np.inf, np.inf
    return float(np.median(v)), float(v.min()), float(v.max())


def sd_theta(theta, seeds=SEEDS, rec_kw=None, **fim_kw):
    rec_kw = rec_kw or {}
    return over_seeds(lambda s: sd_th(fim(record(theta, seed=s, **rec_kw),
                                          **fim_kw)), seeds)


# ------------------------------------------------------ the steady state
def m_ss(theta, Km, ke, L, inc=0.0):
    """Positive root of (theta + inc) g(m) + k_e m = L, any of them arrays."""
    return level(theta, Km, ke, L, inc)


def lnm_rows_fd(r, h=1e-6):
    """
    d ln m_i / d(ln theta, ln K_m, b, ln k_e) by central differences of the
    steady-state solve itself. Shares no algebra with flux_rows: it never
    writes a sensitivity down, it re-solves the balance.
    """
    base = np.array([np.log(r["theta"]), np.log(r["Km"]), 0.0, np.log(r["ke"])])
    out = np.zeros((r["m"].size, 4))
    for j in range(4):
        hj = h * (r["L"] if j == I_B else 1.0)

        def lnm(x):
            th, K, b, ke = np.exp(x[0]), np.exp(x[1]), x[2], np.exp(x[3])
            return np.log(m_ss(th, K, ke, r["loads"] - b, r["inc"]))
        xp, xm = base.copy(), base.copy()
        xp[j] += hj
        xm[j] -= hj
        out[:, j] = (lnm(xp) - lnm(xm)) / (2 * hj)
    return out


def fim_lnm(r, noise="perdraw", ke_prior=KE_PRIOR, load_frac=LOAD_FRAC,
            assay_cv=ASSAY_CV):
    """
    The brief's form: observations of ln m with Gaussian noise. The ln-m noise
    is the flux noise divided by m C'(m), which is what a 5% assay CV and a 5%
    load-record scatter are on the log-marker scale.
    """
    R = lnm_rows_fd(r)
    s_lnm = sigmas(r, noise, load_frac, assay_cv) / (r["m"] * slope(r))
    F = (R / s_lnm[:, None] ** 2).T @ R
    F[I_B, I_B] += (load_frac * r["L"]) ** -2
    F[I_KE, I_KE] += ke_prior ** -2
    return F


# --------------------------------------------------------------- profiling
#  n6_profile_vs_crb.py's likelihood, re-implemented with an analytic-Jacobian
#  least-squares inner solve so that a wide grid and many records are cheap.
#  Same record convention as n6: the marker levels are the regressors, and the
#  ordinate (the recorded load) carries the whole noise sigma.
from scipy.optimize import least_squares


def noisy(r, seed, noise="repo"):
    """n6's observation: the clearance ordinate plus Gaussian sigma."""
    sg = sigmas(r, noise)
    y = (r["theta"] + r["inc"]) * r["m"] / (r["Km"] + r["m"]) + r["ke"] * r["m"]
    rng = np.random.default_rng(10_000 + seed)
    return y + rng.normal(0.0, sg), sg


def _resid(x, lt, r, y, sg, sb, ke_prior, ke0, km_prior=None, km0=KM0):
    lk, b, lke = x
    th, K, ke = np.exp(lt), np.exp(lk), np.exp(lke)
    m, inc = r["m"], r["inc"]
    pred = (th + inc) * m / (K + m) + ke * m + b
    res = [(y - pred) / sg, [b / sb]]
    if ke_prior:
        res.append([(lke - np.log(ke0)) / ke_prior])
    if km_prior:
        res.append([(lk - np.log(km0)) / km_prior])
    return np.concatenate(res)


def _jac(x, lt, r, y, sg, sb, ke_prior, ke0, km_prior=None, km0=KM0):
    lk, b, lke = x
    th, K, ke = np.exp(lt), np.exp(lk), np.exp(lke)
    m, inc = r["m"], r["inc"]
    J = np.zeros((m.size + 1 + bool(ke_prior) + bool(km_prior), 3))
    J[:m.size, 0] = (th + inc) * K * m / (K + m) ** 2 / sg
    J[:m.size, 1] = -1.0 / sg
    J[:m.size, 2] = -ke * m / sg
    J[m.size, 1] = 1.0 / sb
    k = m.size + 1
    if ke_prior:
        J[k, 2] = 1.0 / ke_prior
        k += 1
    if km_prior:
        J[k, 0] = 1.0 / km_prior
    return J


def profile(r, y, sg, lt_grid, ke_prior=KE_PRIOR, sb=None, km_prior=None,
            starts=(-1.5, -0.5, 0.0, 0.5, 1.5, 3.0)):
    """
    -log L profiled over (ln K_m, b, ln k_e) at each fixed ln theta, with
    n6's two priors. Multi-start in ln K_m; each start then warm-starts the
    next grid point, so a branch found anywhere propagates.
    """
    sb = LOAD_FRAC * r["L"] if sb is None else sb
    args = (r, y, sg, sb, ke_prior, KE, km_prior)
    out, best_x = [], []
    prev = None
    for lt in lt_grid:
        cands = [np.array([np.log(KM0) + d, 0.0, np.log(KE)]) for d in starts]
        if prev is not None:
            cands.append(prev)
        best = None
        for x0 in cands:
            s = least_squares(_resid, x0, jac=_jac, args=(lt, *args),
                              method="lm", xtol=1e-12, ftol=1e-12, max_nfev=2000)
            if best is None or s.cost < best.cost:
                best = s
        prev = best.x
        out.append(best.cost)
        best_x.append(best.x)
    return np.array(out), np.array(best_x)
