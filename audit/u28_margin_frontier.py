"""
The frontier again, with the DECISION MARGIN as the target instead of
`ln theta`. Which nuisance binds the decision, as against which binds the
estimate?

THE AUTHOR'S ARGUMENT, AND IT CHECKS OUT. With `theta_crit` at truth,
`K_m = 100` and `m_tol = 40`, the margin's derivatives in log parameters are

    d mu / d ln theta = +3000.00
    d mu / d ln K_m   = -2388.41      (at the UNSAFE amplitude, 3.58x baseline)
    d mu / d ln k_e   =  +700.00

confirmed analytically, by central differences, and against
`design.run.margin_and_gradient`'s own `grad`. **A one per cent relative error in
`K_m` moves the margin 79.6% as much as a one per cent error in theta.**

Set beside what five rounds established: an exactly known `K_m` buys about one
per cent on `sd(ln theta)`. **Both are right and they are about different
quantities.** `mu` is a different function of the parameters than `theta` is, so
a sensitivity ranking does not transfer between them. This repository has been
carrying "K_m does not bind" as though it settled both; it settles one.

THE RATIO IS NOT A CONSTANT, WHICH THE ARGUMENT NEEDS STATED. It carries
`s_peak`: 0.796 at the unsafe amplitude, 0.635 at the safe one, and 0.237 at the
default `1.5 x clearance_A(M_BAR)`. The unsafe side is the costly error, so the
large figure is the relevant one, but it is a figure about one amplitude.

AND SENSITIVITY IS NOT BINDING. What binds is sensitivity times uncertainty, and
the two enter `mu` with OPPOSITE signs in log coordinates, which is the
cancellation `margin_and_gradient`'s own docstring calls the margin theorem. So
the question cannot be settled from the gradient and is settled here by
propagating the full covariance.

WHAT THE GATE IS, stated because `mu` has no inherited one. Two are reported and
neither is smuggled in:
  COMPARABLE   sd(mu)/theta <= 0.20, the same 20% the estimate's gate uses, in
               the units `mu` is already in.
  DECIDABLE    sd(mu) <= |mu| at the unsafe amplitude, one standard deviation
               between the estimate and the boundary it is compared against.

Nothing here records an expected result in advance.

Run:  python audit/u28_margin_frontier.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

import design.run as D
from common.params import _A as P

import u9_ke_frontier as U9

KM, KE, LB = U9.KM, U9.KE, U9.LB
LOAD_FRAC = U9.LOAD_FRAC
BASE = D.clearance_A(D.M_BAR, P)
MT = D.M_TOL
#: the unsafe-side amplitude, which is the costly error and the author's case
S_PEAK = 3.58 * BASE


def mu_at(theta, s_peak=S_PEAK):
    """mu at the true parameters for this capacity. theta_crit does NOT move."""
    return theta - (s_peak - KE * MT) * (KM + MT) / MT


def gates(theta):
    """
    The two gates, BOTH OF WHICH DEPEND ON theta, which the first version got
    wrong.

    MY DEFECT, FIXED HERE. The first run hardcoded the decidable gate as
    |mu| at theta = 3000, which is 343.8, and applied it at 600 and 300 as
    well. But mu = theta - theta_crit and theta_crit is fixed, so |mu| GROWS as
    capacity falls: 343.8 at nominal, 2743.8 at 600, 3043.8 at 300. Using the
    nominal value everywhere made the disease-range rows read "the existing
    prior already suffices" against a gate eight times tighter than the one
    that belongs there. The comparable gate had the same flaw, being pinned at
    0.20 x 3000 rather than 0.20 x theta. **A gate that does not move with the
    quantity it judges is the footing defect wearing a different hat.**
    """
    return {"comparable": 0.20 * theta, "decidable": abs(mu_at(theta))}


def grad_mu(theta, s_peak=S_PEAK):
    """d mu / d(ln theta, ln K_m, b, ln k_e). b does not enter mu."""
    return np.array([theta,
                     -KM * (s_peak - KE * MT) / MT,
                     0.0,
                     KE * (KM + MT)])


def cov(theta, m, inc, sigma, ke_prior=U9.KE_PRIOR, km_prior=None):
    """
    u9's Fisher, with an OPTIONAL prior on ln K_m.

    u9 puts no prior on ln K_m at all, because the record determines it. That is
    right for the estimate and it is the thing under test for the decision, so a
    width can be supplied here and swept.
    """
    gm = m / (KM + m)
    Q = np.stack([theta * gm,
                  -(theta + inc) * KM * m / (KM + m) ** 2,
                  np.ones_like(m),
                  KE * m], axis=1)
    pri = [np.inf,
           np.inf if km_prior is None else km_prior ** 2,
           (LOAD_FRAC * LB) ** 2,
           np.inf if ke_prior is None else ke_prior ** 2]
    F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri]) \
        + Q.T @ Q / sigma ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return None
    return np.linalg.inv(F)


def sd_mu(theta, m, inc, sigma, **kw):
    C = cov(theta, m, inc, sigma, **kw)
    if C is None:
        return np.inf
    g = grad_mu(theta)
    v = float(g @ C @ g)
    return float(np.sqrt(v)) if v > 0 else np.inf


def med_mu(theta, n, **kw):
    sig = U9.sigma_at(theta)
    vals = np.array([sd_mu(theta, *U9.levels(theta, n, s), sig, **kw)
                     for s in U9.SEEDS], float)
    f = vals[np.isfinite(vals)]
    if not f.size:
        return np.inf
    return float(np.median(f))


def frontier_mu(theta, n, gate, which="ke", lo=1e-4, hi=5.0, tol=1e-6):
    """The prior width on `which` at which sd(mu) lands exactly on `gate`."""
    def f(w):
        kw = {"ke_prior": w} if which == "ke" else {"km_prior": w}
        if which == "km":
            kw["ke_prior"] = U9.KE_PRIOR
        return med_mu(theta, n, **kw)
    if f(hi) <= gate:
        return None, "the existing prior already suffices"
    if f(lo) > gate:
        return np.nan, "unreachable: even an exact value misses the gate"
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if f(mid) <= gate:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), ""


if __name__ == "__main__":
    g = grad_mu(P.theta)
    mu_true = P.theta - (S_PEAK - KE * MT) * (KM + MT) / MT
    print("=" * 98)
    print("THE FRONTIER WITH mu AS THE TARGET: which nuisance binds the DECISION?")
    print("=" * 98)
    print(f"\n  s_peak = {S_PEAK:.1f} (unsafe side), mu at truth = {mu_true:+.1f}")
    print(f"  d mu/d ln theta = {g[0]:+.2f}   d mu/d ln K_m = {g[1]:+.2f}   "
          f"d mu/d ln k_e = {g[3]:+.2f}")
    print(f"  |K_m| / |theta| = {abs(g[1])/abs(g[0]):.3f}      "
          f"|k_e| / |theta| = {abs(g[3])/abs(g[0]):.3f}")

    print(f"\n  gates, both of which move with theta:")
    for th in (3000.0, 600.0, 300.0):
        gg = gates(th)
        print(f"    theta = {th:6.0f}   mu = {mu_at(th):+9.1f}   "
              f"comparable sd(mu) <= {gg['comparable']:7.1f}   "
              f"decidable sd(mu) <= {gg['decidable']:7.1f}")

    print()
    print("-" * 98)
    print("PART A. Where the variance of mu comes from, at the recommended "
          "design and the published prior")
    print("-" * 98)
    n = 25
    sig = U9.sigma_at(P.theta)
    m, inc = U9.levels(P.theta, n, U9.SEEDS[0] if hasattr(U9, "SEEDS") else 0)
    C = cov(P.theta, m, inc, sig)
    names = ["ln theta", "ln K_m", "b", "ln k_e"]
    tot = float(g @ C @ g)
    print(f"  sd(mu) = {np.sqrt(tot):.1f}, and sd(ln theta) at the same design "
          f"= {U9.sd_ln_theta(P.theta, m, inc, sig):.4f}")
    print(f"\n  {'term':>22} {'contribution to var(mu)':>24} {'share':>9}")
    for i in range(4):
        c = g[i] * C[i, i] * g[i]
        print(f"  {'diagonal ' + names[i]:>22} {c:24.1f} {100*c/tot:8.1f}%")
    off = 0.0
    for i in range(4):
        for j in range(4):
            if i != j:
                off += g[i] * C[i, j] * g[j]
    print(f"  {'all cross terms':>22} {off:24.1f} {100*off/tot:8.1f}%")
    diag = sum(g[i] * C[i, i] * g[i] for i in range(4))
    print(f"\n  sd from the diagonals alone would be {np.sqrt(diag):.1f}; "
          f"the actual sd is {np.sqrt(tot):.1f}, a factor of "
          f"{np.sqrt(diag/tot):.2f} of cancellation.")
    print(f"\n  The cross terms are the cancellation margin_and_gradient's "
          f"docstring names:\n  theta and K_m enter mu with opposite signs and "
          f"are correlated in the Fisher.")

    print()
    print("-" * 98)
    print("PART B. The frontier on each nuisance, both gates, three capacities")
    print("-" * 98)
    print(f"  {'theta':>7} {'gate':>12} {'sd(mu) now':>11} {'|mu|/sd':>8} "
          f"{'k_e width':>22} {'K_m width':>22}")
    for theta in (3000.0, 600.0, 300.0):
        now = med_mu(theta, U9.TIER1_N)
        for glab, gate in sorted(gates(theta).items()):
            ke_w, ke_n = frontier_mu(theta, U9.TIER1_N, gate, "ke")
            km_w, km_n = frontier_mu(theta, U9.TIER1_N, gate, "km")
            def fmt(w, note):
                if note:
                    return note[:20]
                return f"{w:.4f} ({np.exp(w):.3f}x)"
            z = abs(mu_at(theta)) / now if np.isfinite(now) and now > 0 else float("nan")
            print(f"  {theta:7.0f} {glab:>12} {now:11.1f} {z:8.2f} "
                  f"{fmt(ke_w, ke_n):>22} {fmt(km_w, km_n):>22}")

    print()
    print("-" * 98)
    print("PART C. Against the ESTIMATE's frontier, from u9, at the same designs")
    print("-" * 98)
    print(f"  {'theta':>7} {'k_e width for sd(ln theta) <= 0.20':>38}")
    for theta in (3000.0, 600.0, 300.0):
        w, note = U9.frontier(theta, U9.TIER1_N, U9.PUBLISHED_KAPPA_CV)
        print(f"  {theta:7.0f} {(note[:34] if note else f'{w:.4f} ({np.exp(w):.3f}x)'):>38}")
    print("""
  The comparison is the point. If the k_e width binds the estimate while the
  K_m width binds the margin, the specification is asking for the wrong
  measurement for the decision it is written against, and that is a statement
  about Tier 2 and Tier 3 rather than about either number.""")
