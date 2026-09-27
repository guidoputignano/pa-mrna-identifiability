"""
Does the profile-marginal disagreement attach to every Fisher figure the
paper computes at a 0.8 log-prior, or only to the ranking it was found on?

THE FINDING IT COMES FROM. audit/u6_recheck_sampling.py profiles the
posterior over ln theta instead of inverting an information matrix and
disagrees with audit/u5_sigma_split.py about which nuisance binds at a low
assay CV. Two optimiser faults were found and fixed and neither changed it.
What did change it was narrowing the k_e log-prior from the paper's 0.8 to
0.2: the two methods then agree to 3% and rank the same.

THE INFERENCE THAT HAS TO BE TESTED. A profile and a marginal coincide
exactly for a quadratic log-posterior and differ otherwise. If the
disagreement vanishes as the prior narrows, the cause is that the posterior
is not quadratic in the ln k_e direction at a width of 0.8, which is a
factor of 2.2. That cause is not a property of the ranking. It is a property
of the PRIOR WIDTH, and every number in this paper computed by inverting an
information matrix at that width inherits it, the 118 included.

This script tests that inference rather than asserting it, three ways.

  A  the gap between profile and marginal against the prior width, at the
     paper's own headline designs rather than at the one the disagreement
     was found on. If the gap closes as the width narrows, at designs that
     have nothing to do with the ranking, the cause is the width.

  B  a direct non-quadraticity measure that needs no second method at all:
     the profile's left and right one-sigma half-widths. They are equal for
     a quadratic and unequal otherwise, so their asymmetry measures the
     departure without reference to any marginal.

  C  what it costs the paper's own headline numbers, stated as the interval
     a reader should carry beside each.

FOOTING. As audit/PHASE11_SAMPLING.md section 0. The two headline designs
are the recommended record, 25 monthly draws at nominal capacity with the
increment exact, which the paper reports as 0.170; and the Tier 1 record,
118 monthly draws with the increment at 14.02%, which the paper reports as
landing on the 0.20 gate. Medians over 8 intake records on
audit/t2_tier1_margin.py's convention, seeds stated.

Nothing here records an expected result in advance.

Run:  python audit/u12_prior_quadratic.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import minimize

from design.run import (ASSAY_CV, DOSE_INTERVAL, D_AMP, LAMBDA, LOAD_FRAC,
                        load)

KM, KE = 100.0, 5.0
LB = float(load(0.0))
GATE = 0.20
MONTH = 30.4
SB = LOAD_FRAC * LB
T_DOSE = 50.0 * DOSE_INTERVAL
PEAK = D_AMP / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))
LOAD_CV = 0.10
PUBLISHED_KAPPA_CV = 0.1402

#: (label, theta, n visits, increment CV or None)
DESIGNS = [
    ("recommended: 25 monthly, nominal", 3000.0, 25, None),
    ("Tier 1: 118 monthly, increment 14.02%", 3000.0, 118,
     PUBLISHED_KAPPA_CV),
    ("disease range: 25 monthly, theta = 300", 300.0, 25, None),
]


def op_level(theta):
    B = KE * KM + theta - LB
    return (-B + np.sqrt(B ** 2 + 4.0 * KE * LB * KM)) / (2.0 * KE)


def sigma_of(theta):
    mbar = op_level(theta)
    slope = theta * KM / (KM + mbar) ** 2 + KE
    return float(np.hypot(LOAD_FRAC * LB, slope * ASSAY_CV * mbar))


def record(theta, n, seed):
    t = MONTH * np.arange(n, dtype=float)
    tau = np.mod(t - T_DOSE, DOSE_INTERVAL)
    inc = PEAK * np.exp(-LAMBDA * tau)
    Ls = LB * (1.0 + LOAD_CV *
               np.random.default_rng(seed).standard_normal(t.size))
    B = KE * KM + theta + inc - Ls
    m = (-B + np.sqrt(B ** 2 + 4.0 * KE * Ls * KM)) / (2.0 * KE)
    return m, inc, Ls


def marginal(theta, m, inc, sigma, ke_prior, kappa_cv=None):
    """sd(ln theta) by inverting the information matrix."""
    gm = m / (KM + m)
    cols = [theta * gm, -(theta + inc) * KM * m / (KM + m) ** 2,
            np.ones_like(m), KE * m]
    pri = [np.inf, np.inf, SB ** 2, ke_prior ** 2]
    if kappa_cv is not None:
        cols.append(inc * gm)
        pri.append(kappa_cv ** 2)
    Q = np.stack(cols, axis=1)
    F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri]) \
        + Q.T @ Q / sigma ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    return float(np.sqrt(np.linalg.inv(F)[0, 0]))


def nlp(lt, theta0, m_dat, inc, Ls, sigma, ke_prior, kappa_cv):
    """
    Profile of the negative log posterior at a fixed ln theta.

    The EXACT flux residual at the observed level, not a first-order
    conversion of a level difference through the local gain, which is the
    error audit/u6_recheck_sampling.py made and REVIEW_LOG R101 records.
    """
    free = ["Km", "b", "ke"] + (["kappa"] if kappa_cv is not None else [])

    def obj(v):
        d = dict(zip(free, v))
        Km = KM * np.exp(d.get("Km", 0.0))
        bs = d.get("b", 0.0)
        ke = KE * np.exp(d.get("ke", 0.0))
        kap = np.exp(d.get("kappa", 0.0))
        if Km <= 0 or ke <= 0:
            return 1e9
        L = Ls + bs * SB
        fl = (np.exp(lt) + kap * inc) * m_dat / (Km + m_dat) + ke * m_dat
        r = 0.5 * np.sum(((fl - L) / sigma) ** 2)
        r += 0.5 * bs ** 2
        r += 0.5 * (d.get("ke", 0.0) / ke_prior) ** 2
        if kappa_cv is not None:
            r += 0.5 * (d.get("kappa", 0.0) / kappa_cv) ** 2
        return float(r)

    r = minimize(obj, np.zeros(len(free)), method="Powell",
                 options=dict(xtol=1e-10, ftol=1e-12, maxiter=4000,
                              maxfev=20000))
    return float(r.fun)


def profile_halfwidths(theta, m, inc, Ls, sigma, ke_prior, kappa_cv=None):
    """(left, right) one-sigma half-widths of the profile in ln theta."""
    lt0 = np.log(theta)
    base = nlp(lt0, theta, m, inc, Ls, sigma, ke_prior, kappa_cv)
    out = []
    for sign in (-1.0, +1.0):
        lo, hi = 1e-4, 6.0
        for _ in range(18):
            mid = 0.5 * (lo + hi)
            if nlp(lt0 + sign * mid, theta, m, inc, Ls, sigma, ke_prior,
                   kappa_cv) - base < 0.5:
                lo = mid
            else:
                hi = mid
        out.append(0.5 * (lo + hi))
    return out[0], out[1]


def over(fn, theta, n, kappa_cv, seeds, *a):
    vals = []
    for s in seeds:
        m, inc, Ls = record(theta, n, s)
        vals.append(fn(theta, m, inc, Ls, sigma_of(theta), *a,
                       kappa_cv=kappa_cv))
    return vals


# ======================================================================
def part_a(seeds=range(8)):
    print("=" * 92)
    print("PART A.  profile against marginal, against the prior width, at the")
    print("         paper's OWN headline designs")
    print("=" * 92)
    print("\n  If the gap closes as the width narrows at designs that have")
    print("  nothing to do with the ranking the disagreement was found on,")
    print("  the cause is the width and not the ranking.\n")
    for lab, theta, n, kap in DESIGNS:
        print(f"  {lab}")
        print(f"    {'sd(ln k_e)':>11} {'marginal':>10} {'profile':>10} "
              f"{'ratio':>8} {'gap':>9}")
        print("    " + "-" * 52)
        for w in (0.8, 0.4, 0.2, 0.1, 0.05):
            mars, profs = [], []
            for s in seeds:
                m, inc, Ls = record(theta, n, s)
                sig = sigma_of(theta)
                mars.append(marginal(theta, m, inc, sig, w, kap))
                l, r = profile_halfwidths(theta, m, inc, Ls, sig, w, kap)
                profs.append(0.5 * (l + r))
            mm, pp = float(np.median(mars)), float(np.median(profs))
            print(f"    {w:11.2f} {mm:10.4f} {pp:10.4f} {pp / mm:8.3f} "
                  f"{abs(pp - mm):9.4f}")
        print()


def part_b(seeds=range(8)):
    print("=" * 92)
    print("PART B.  non-quadraticity measured directly, with no second method")
    print("=" * 92)
    print("\n  The profile's left and right one-sigma half-widths are equal")
    print("  for a quadratic log-posterior. Their ratio measures the")
    print("  departure without reference to any marginal, so it cannot be")
    print("  an artefact of disagreement between two implementations.\n")
    for lab, theta, n, kap in DESIGNS:
        print(f"  {lab}")
        print(f"    {'sd(ln k_e)':>11} {'left':>9} {'right':>9} "
              f"{'right/left':>11} {'marginal':>10}")
        print("    " + "-" * 52)
        for w in (0.8, 0.4, 0.2, 0.1, 0.05):
            ls, rs, ms = [], [], []
            for s in seeds:
                m, inc, Ls = record(theta, n, s)
                sig = sigma_of(theta)
                l, r = profile_halfwidths(theta, m, inc, Ls, sig, w, kap)
                ls.append(l)
                rs.append(r)
                ms.append(marginal(theta, m, inc, sig, w, kap))
            L, R = float(np.median(ls)), float(np.median(rs))
            print(f"    {w:11.2f} {L:9.4f} {R:9.4f} {R / L:11.3f} "
                  f"{float(np.median(ms)):10.4f}")
        print()
    print("  A ratio of 1.000 is a quadratic posterior and the Fisher")
    print("  standard deviation is the whole story there. A ratio far from")
    print("  1 says the Fisher number is a curvature at the optimum and not")
    print("  a summary of the posterior, whatever method computed it.")


def part_c(seeds=range(8)):
    print("=" * 92)
    print("PART C.  what a reader should carry beside each headline number")
    print("=" * 92)
    print("\n  At the paper's own 0.8 prior, the marginal it quotes against")
    print("  the asymmetric profile interval the same posterior actually")
    print("  has. Reported as the interval on theta itself, which is what a")
    print("  clinician reads.\n")
    for lab, theta, n, kap in DESIGNS:
        ls, rs, ms = [], [], []
        for s in seeds:
            m, inc, Ls = record(theta, n, s)
            sig = sigma_of(theta)
            l, r = profile_halfwidths(theta, m, inc, Ls, sig, 0.8, kap)
            ls.append(l)
            rs.append(r)
            ms.append(marginal(theta, m, inc, sig, 0.8, kap))
        L, R, M = (float(np.median(x)) for x in (ls, rs, ms))
        print(f"  {lab}")
        print(f"    quoted marginal sd(ln theta)   {M:.4f}   "
              f"-> symmetric interval x/{np.exp(M):.3f}")
        print(f"    profile half-widths            {L:.4f} down, "
              f"{R:.4f} up")
        print(f"    interval on theta              "
              f"[{theta * np.exp(-L):.0f}, {theta * np.exp(R):.0f}] "
              f"against theta = {theta:.0f}")
        print(f"    against the marginal, the profile is "
              f"{L / M:.2f}x its width DOWNWARD and {R / M:.2f}x UPWARD")
        print(f"    so the symmetric interval UNDERSTATES how far capacity")
        print(f"    may fall below the estimate, which is the unsafe")
        print(f"    direction for a decision about whether capacity clears")
        print(f"    a threshold")
        print()


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
