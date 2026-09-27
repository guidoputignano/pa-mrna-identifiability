"""
The profile, the simulated spread and the Tier 1 gap, recomputed on the PAPER'S
committed post-dose footing, so that no ratio is formed across two footings.

THE DEFECT THIS FIXES IS MINE. R142 found that `sd_monthly` is `0.1699` in the
repository's generated macro while `audit/u12_prior_quadratic.py` computes
`0.1695`, and `PHASE8_SIMULATIONS.md:623` documents the difference exactly: the
same quantity with the `t = 0` draw read POST-dose inside the settling window
(`0.1699`, as committed), PRE-dose (`0.1695`), or dropped (`0.1721`). I then
quoted an object spread of `0.2002 / 0.1699 = 1.178`, **which divides a
post-dose marginal into a pre-dose profile and is exactly the error this audit
exists to catch.** The author caught it.

THE RULE THAT FOLLOWS. The object spread is a comparison of two objects at ONE
design. Both objects must be computed on one footing, and since the paper is the
artefact, the paper's footing wins and the audit moves to it.

HOW THE FOOTING IS OBTAINED, and it is not by reimplementing it. The record is
built by `src/design/run.py`'s own `design_record`, which uses the paper's own
`delta()` for the increment, rather than by this file rebuilding the dose train.
**PART A is the anchor: the marginal computed here must reproduce the committed
`sd_monthly` of 0.1699.** If it does not, the footing is not yet matched and
nothing below it means anything.

Nothing here records an expected result in advance.

Run:  python audit/u22_postdose_footing.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from scipy.optimize import minimize

import design.run as D
from common.params import _A as P

KM, KE = float(P.Km), float(P.ke)
LB = float(D.load(0.0))
SB = D.LOAD_FRAC * LB
KE_PRIOR = 0.8

#: the recommended design and Tier 1, in the paper's own time grids
MONTHLY = np.arange(0.0, 730.0, 30.4)
TIER1 = 30.4 * np.arange(118, dtype=float)

#: the committed value this file must reproduce before anything else counts
SD_MONTHLY_COMMITTED = 0.1699


def sigma_of(theta, m_bar=None):
    """Flux-space noise, on the paper's own convention."""
    mb = float(D.M_BAR) if m_bar is None else m_bar
    slope = theta * KM / (KM + mb) ** 2 + KE
    return float(np.hypot(SB, slope * D.ASSAY_CV * mb))


def rows(theta, lv, inc):
    """d(flux)/d(ln theta, ln Km, b, ln ke) at the paper's record."""
    g = lv / (KM + lv)
    return np.stack([theta * g,
                     -(theta + inc) * KM * lv / (KM + lv) ** 2,
                     np.ones_like(lv),
                     KE * lv], axis=1)


def marginal(theta, lv, inc, sigma, ke_prior=KE_PRIOR, kappa_cv=None):
    Q = rows(theta, lv, inc)
    pri = [0.0, 0.0, 1.0 / SB ** 2, 1.0 / ke_prior ** 2]
    if kappa_cv is not None:
        #  Tier 1 as PUBLISHED carries the dose-scale nuisance. Omitting it
        #  gives a different design, and the first version of PART D did,
        #  which is why its marginal was 0.1542 against u12's 0.2000.
        Q = np.column_stack([Q, inc * lv / (KM + lv)])
        pri = pri + [1.0 / kappa_cv ** 2]
    F = np.diag(pri) + Q.T @ Q / sigma ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    return float(np.sqrt(np.linalg.inv(F)[0, 0]))


def nlp(lt, theta, lv, inc, Ls, sigma, ke_prior, kappa_cv=None):
    """
    Negative log posterior profiled over the nuisances at a fixed ln theta.

    The EXACT flux residual at the observed level, which is the form R101
    records as the one that does not collapse.
    """
    #  THE PROFILE MUST FREE EVERY NUISANCE THE MARGINAL CARRIES. The first
    #  version of PART D added the dose-scale column to the marginal and not
    #  here, so a five-parameter marginal was compared with a four-parameter
    #  profile and the gap came out NEGATIVE, -9.1%, which is impossible: a
    #  profile that accounts for non-quadraticity cannot be narrower than the
    #  marginal at the same design. Had that been wired, the guard would have
    #  certified a negative gap.
    n = 3 if kappa_cv is None else 4

    def obj(v):
        Km = KM * np.exp(v[0])
        bs = v[1]
        ke = KE * np.exp(v[2])
        kap = np.exp(v[3]) if kappa_cv is not None else 1.0
        if Km <= 0 or ke <= 0:
            return 1e9
        fl = (np.exp(lt) + kap * inc) * lv / (Km + lv) + ke * lv
        r = 0.5 * np.sum(((fl - (Ls + bs * SB)) / sigma) ** 2)
        r += 0.5 * bs ** 2 + 0.5 * (v[2] / ke_prior) ** 2
        if kappa_cv is not None:
            r += 0.5 * (v[3] / kappa_cv) ** 2
        return float(r)

    return float(minimize(obj, np.zeros(n), method="Powell",
                          options=dict(xtol=1e-10, ftol=1e-12,
                                       maxiter=4000, maxfev=20000)).fun)


def profile(theta, lv, inc, Ls, sigma, ke_prior=KE_PRIOR, kappa_cv=None):
    """(left, right) one-sigma half-widths, and their mean as the sd."""
    lt0 = np.log(theta)
    base = nlp(lt0, theta, lv, inc, Ls, sigma, ke_prior, kappa_cv)
    out = []
    for sign in (-1.0, +1.0):
        lo, hi = 1e-4, 6.0
        for _ in range(18):
            mid = 0.5 * (lo + hi)
            if nlp(lt0 + sign * mid, theta, lv, inc, Ls, sigma,
                   ke_prior, kappa_cv) - base < 0.5:
                lo = mid
            else:
                hi = mid
        out.append(0.5 * (lo + hi))
    return out[0], out[1]


def build(times, seed, load_cv=None):
    """The paper's own record, which is what carries the footing."""
    lv, inc = D.design_record(times, load_cv=load_cv, seed=seed)
    cv = D.LOAD_CV if load_cv is None else load_cv
    Ls = (LB * (1.0 + cv * np.random.default_rng(seed).standard_normal(times.size))
          if cv else np.full(times.size, LB))
    return lv, inc, Ls


def part_a():
    print("=" * 92)
    print("PART A.  the anchor: does the marginal on this footing reproduce the")
    print("         committed sd_monthly of 0.1699?")
    print("=" * 92)
    print("\n  Record built by src/design/run.py's own design_record, so the")
    print("  footing is the paper's rather than a reimplementation of it.\n")
    sig = sigma_of(P.theta)
    vals = []
    for s in D.SEEDS:
        lv, inc, Ls = build(MONTHLY, s)
        vals.append(marginal(P.theta, lv, inc, sig))
    v = np.array(vals, float)
    v = v[np.isfinite(v)]
    med = float(np.median(v))
    print(f"    seeds           {len(v)}")
    print(f"    median          {med:.4f}   [{v.min():.4f}, {v.max():.4f}]")
    print(f"    committed       {SD_MONTHLY_COMMITTED:.4f}")
    print(f"    difference      {abs(med - SD_MONTHLY_COMMITTED):.4f}"
          f"   ({100 * abs(med - SD_MONTHLY_COMMITTED) / SD_MONTHLY_COMMITTED:.2f}%)")
    ok = abs(med - SD_MONTHLY_COMMITTED) < 5e-4
    print(f"\n  FOOTING MATCHED: {'YES' if ok else 'NO'}")
    if not ok:
        print("  Nothing below this line is on the committed footing, and the")
        print("  numbers must not be wired until it matches.")
    return med, ok


def part_b(marg):
    print()
    print("=" * 92)
    print("PART B.  the three objects, all on the SAME footing")
    print("=" * 92)
    sig = sigma_of(P.theta)
    los, his = [], []
    for s in D.SEEDS:
        lv, inc, Ls = build(MONTHLY, s)
        a, b = profile(P.theta, lv, inc, Ls, sig)
        los.append(a)
        his.append(b)
    lo, hi = float(np.median(los)), float(np.median(his))
    prof = 0.5 * (lo + hi)
    print(f"\n  recommended design, {len(MONTHLY)} monthly, nominal capacity\n")
    print(f"    marginal (Fisher)          {marg:.4f}")
    print(f"    profile, left half-width   {lo:.4f}")
    print(f"    profile, right half-width  {hi:.4f}")
    print(f"    profile, mean half-width   {prof:.4f}")
    print(f"    right/left                 {hi / lo:.3f}")
    print(f"\n    OBJECT SPREAD, profile/marginal, both post-dose:  "
          f"{prof / marg:.4f}  ({100 * (prof / marg - 1):.1f}%)")
    print("""
  Both numbers now come from one record built by the paper's own
  design_record, so the ratio is a comparison of two objects at one design
  rather than a comparison of two footings.""")
    return lo, hi, prof


def part_c(marg, prof):
    print()
    print("=" * 92)
    print("PART C.  the simulated sampling spread, on the same footing")
    print("=" * 92)
    print("\n  The third object. Same record, same seeds; the spread is the")
    print("  robust one AND the standard deviation, per section 0 as amended")
    print("  at R132, because the comparison is against a Fisher width.\n")
    sig = sigma_of(P.theta)
    #  A BARE lstsq ON THE FLUX RESIDUAL IS NOT THE ESTIMATOR, and the first
    #  version of this part used one: it omits both priors, so the near-null
    #  direction is unconstrained and the spread came out at 5.37 for a
    #  quantity whose Fisher width is 0.17. Refit all four parameters with the
    #  priors the design actually carries, and draw the prior-carrying
    #  nuisances from their priors, which is R98.
    from scipy.optimize import least_squares
    est = []
    for s_ in D.SEEDS:
        lv, inc, Ls = build(MONTHLY, s_)
        rng = np.random.default_rng(10000 + s_)
        for _ in range(40):
            b_true = SB * rng.standard_normal()
            ke_true = KE * float(np.exp(KE_PRIOR * rng.standard_normal()))
            fl = ((P.theta + inc) * lv / (KM + lv) + ke_true * lv)
            y = fl + sig * rng.standard_normal(lv.size)

            def res(q):
                th, Km, bs, ke = (np.exp(q[0]), KM * np.exp(q[1]),
                                  q[2], KE * np.exp(q[3]))
                pred = (th + inc) * lv / (Km + lv) + ke * lv
                return np.concatenate([(pred - y - bs * SB) / sig,
                                       [q[2], q[3] / KE_PRIOR]])

            try:
                r = least_squares(res,
                                  np.array([np.log(P.theta), 0.0,
                                            b_true / SB,
                                            np.log(ke_true / KE)]),
                                  method="trf", xtol=1e-12, ftol=1e-12,
                                  max_nfev=400)
                est.append(float(r.x[0]) - float(np.log(P.theta)))
            except Exception:
                continue
    e = np.array(est, float)
    e = e[np.isfinite(e)]
    rob = 0.7413 * float(np.subtract(*np.percentile(e, [75, 25])))
    sd = float(np.std(e, ddof=1))
    print(f"    draws                      {e.size}")
    print(f"    robust spread              {rob:.4f}")
    print(f"    standard deviation         {sd:.4f}")
    print(f"    marginal                   {marg:.4f}")
    print(f"    profile                    {prof:.4f}")
    print(f"\n  ordering, all post-dose:  marginal {marg:.4f}, "
          f"simulated {sd:.4f}, profile {prof:.4f}")
    return rob, sd


def part_d():
    print()
    print("=" * 92)
    print("PART D.  Tier 1's gap, on the same footing")
    print("=" * 92)
    print(f"\n  {len(TIER1)} monthly draws, WITH the dose-scale nuisance the")
    print("  published Tier 1 carries. A DIFFERENT design from PART B, and the")
    print("  macro names must say so.\n")
    sig = sigma_of(P.theta)
    KAPPA_CV = 0.1402
    margs, profs = [], []
    for s in D.SEEDS:
        lv, inc, Ls = build(TIER1, s)
        margs.append(marginal(P.theta, lv, inc, sig, kappa_cv=KAPPA_CV))
        a, b = profile(P.theta, lv, inc, Ls, sig, kappa_cv=KAPPA_CV)
        profs.append(0.5 * (a + b))
    m = float(np.median([x for x in margs if np.isfinite(x)]))
    p = float(np.median(profs))
    print(f"    marginal                   {m:.4f}")
    print(f"    profile                    {p:.4f}")
    print(f"    gap                        {100 * (p / m - 1):.1f}%")
    return m, p


if __name__ == "__main__":
    marg, ok = part_a()
    lo, hi, prof = part_b(marg)
    rob, sd = part_c(marg, prof)
    t1m, t1p = part_d()
    print()
    print("=" * 92)
    print("MACRO VALUES, all on the committed post-dose footing")
    print("=" * 92)
    print(f"""
    sd_monthly_recommended      {marg:.4f}    (committed: {SD_MONTHLY_COMMITTED})
    profile_recommended         {prof:.4f}
    profile_lo_recommended      {lo:.4f}
    profile_hi_recommended      {hi:.4f}
    sim_spread_recommended      {sd:.4f}
    object_spread_recommended_pct  {100 * (prof / marg - 1):.1f}
    sd_tier1                    {t1m:.4f}
    profile_tier1               {t1p:.4f}
    profile_gap_tier1_pct       {100 * (t1p / t1m - 1):.1f}

  Every name carries its design, per the author's ruling that wiring freezes
  names as much as values. FOOTING MATCHED: {'YES' if ok else 'NO'}.""")
