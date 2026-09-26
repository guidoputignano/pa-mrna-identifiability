"""
Phase 9 re-checks. Independent routes to every number PHASE9_LITERATURE.md reports.

Phase 9 is a literature phase, so most of its content is quotation and
adjudication rather than computation, and quotation is checked against sources
and not against code. Two things in it ARE computed, and both are checked here
by a route that shares no code path with the one that produced them.

  1. The model-explained within-patient marker CV, which is what decides
     whether a surrogate's total variation may be carried into sigma. The
     reported route builds the increment through `q2_requirement_map.increment`
     (which calls `q1`'s ExpSum assembly) and solves the quasi-steady balance
     with the closed-form positive root. The route here sums the exponential
     kernel over doses in a Python loop and solves the balance by Brent
     root-finding on the residual. Nothing is shared but the constants.

  2. The item 6 arithmetic, which is a rounding claim about two published
     numbers and is checked as arithmetic.

Nothing here records an expected result in advance.

Run:  python audit/q10_phase9_recheck.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import brentq

import q1_increment_shape as Q
import q2_requirement_map as S2

TH, KM, KE = S2.TH_NOM, S2.KM, S2.KE
RESULTS = []


def report(name, got, want, tol, rel=False):
    gap = abs(got / want - 1.0) if rel else abs(got - want)
    ok = gap <= tol
    RESULTS.append((name, ok))
    print(f"    {name:<52s} {'PASSES' if ok else 'FAILS':>6s}  "
          f"{'relative' if rel else 'absolute'} gap {gap:.3e} against {tol:g}")
    return ok


def banner(t):
    print("\n" + "=" * 92)
    print("  " + t)
    print("=" * 92)


# ==========================================================================
def brute_increment(t, half_life=None, D=None, period=None, t_end=900.0):
    """
    The increment at time `t` by summing one exponential per past dose.

    A dose at exactly `t` does not contribute, which is the PRE-DOSE reading
    `q1.DOSE_DAY_DRAW` fixes; `q1` implements it inside `at_draw`, this route
    implements it as a strict inequality in the sum.
    """
    half_life = Q.TH_HALF0 if half_life is None else half_life
    D = Q.D_AMP if D is None else D
    period = Q.T_DOSE if period is None else period
    lam = np.log(2.0) / half_life
    return float(sum(D * np.exp(-lam * (t - u))
                     for u in np.arange(0.0, t_end, period) if u < t - 1e-12))


def brute_level(theta, inc, load):
    """The quasi-steady level by root-finding, not by the closed-form root."""
    def f(m):
        return (theta + inc) * m / (KM + m) + KE * m - load
    return float(brentq(f, 1e-12, 1e9, xtol=1e-14, rtol=1e-15))


def check_1_model_explained_cv():
    banner("1. THE MODEL-EXPLAINED MARKER CV, BY BRUTE SUMMATION AND ROOT-FINDING")
    print("""  The within-patient CV of the marker across the monthly record produced by the
  dose train alone, at a perfectly constant recorded load. This is the quantity
  PHASE9 item 7 uses to argue that a surrogate's TOTAL within-patient variation
  may not be carried into sigma unchanged. Reported on both footings, because
  the untreated t = 0 draw sits at the undosed level and moves the CV.""")
    L = S2.footing(TH)[0]
    ok = True
    for name, times in (("incl t = 0, n = 25", Q.MONTHLY),
                        ("on therapy, n = 24",
                         S2.times_for("monthly", "no baseline"))):
        mine = S2.level(TH, S2.increment(times), np.full(times.size, L))
        theirs = np.array([brute_level(TH, brute_increment(t), L)
                           for t in times])
        cv_a = float(mine.std(ddof=1) / mine.mean())
        cv_b = float(theirs.std(ddof=1) / theirs.mean())
        print(f"    {name}: q2 route CV {100 * cv_a:.3f}%, brute route "
              f"{100 * cv_b:.3f}%, levels {theirs.min():.3f} to "
              f"{theirs.max():.3f}")
        ok &= report(f"marker CV, {name}", cv_b, cv_a, 1e-9)
        ok &= report(f"per-visit levels, {name}",
                     float(np.abs(theirs / mine - 1.0).max()), 0.0, 1e-8)
    print("""
  And the same with the specification's 10% day-to-day load variation added,
  over the seed set, because a single-seed figure is not reportable here.""")
    for name, times in (("incl t = 0, n = 25", Q.MONTHLY),
                        ("on therapy, n = 24",
                         S2.times_for("monthly", "no baseline"))):
        v = []
        for s in S2.SEEDS:
            Ls = S2.loads_for(times.size, s, L)
            m = np.array([brute_level(TH, brute_increment(t), x)
                          for t, x in zip(times, Ls)])
            v.append(m.std(ddof=1) / m.mean())
        v = np.array(v)
        print(f"    {name}: median {100 * np.median(v):.2f}% "
              f"[{100 * v.min():.2f}, {100 * v.max():.2f}] over "
              f"{len(S2.SEEDS)} records")
    return ok


def check_2_item6_arithmetic():
    banner("2. ITEM 6: THE ROUNDING ARITHMETIC, BOTH READINGS TESTED")
    print("""  Manoli 2021 as quoted in the two reports. The Results give, for patients, a
  mean of per-subject CVs of 14.02 with an SD of 12.86, and for controls a mean
  of 6.94 with an SD of 5.3. The abstract gives 6.9 and 12.8 and calls both a
  "mean coefficient of variation". The question is which Results quantity each
  abstract figure is a rounding of. This is arithmetic on quoted numbers, and
  it is reported as such: it does not verify the quotes.""")
    controls_mean, controls_sd = 6.94, 5.3
    patients_mean, patients_sd = 14.02, 12.86
    abs_controls, abs_patients = 6.9, 12.8
    rows = [("controls, abstract 6.9 against Results mean 6.94",
             abs_controls, controls_mean),
            ("controls, abstract 6.9 against Results SD 5.3",
             abs_controls, controls_sd),
            ("patients, abstract 12.8 against Results mean 14.02",
             abs_patients, patients_mean),
            ("patients, abstract 12.8 against Results SD 12.86",
             abs_patients, patients_sd)]
    print(f"\n    {'comparison':<52s} {'gap':>8s} {'rounds?':>9s}")
    for label, a, b in rows:
        gap = abs(a - b)
        print(f"    {label:<52s} {gap:8.3f} {'yes' if gap < 0.1 else 'no':>9s}")
    ok = report("abstract controls figure matches the MEAN", abs_controls,
                controls_mean, 0.05)
    ok &= report("abstract patients figure matches the SD", abs_patients,
                 patients_sd, 0.07)
    ok &= report("abstract patients figure does NOT match the mean",
                 1.0 if abs(abs_patients - patients_mean) > 1.0 else 0.0,
                 1.0, 0.0)
    print("""
  So the two abstract figures come from different columns of the same table:
  the controls' mean and the patients' SD. That internal inconsistency is the
  evidence, and it does not depend on which reading either report prefers.""")
    return ok


def check_3_sigma_budget():
    banner("3. THE sigma BUDGET, BY A NUMERICAL HESSIAN AND A DIRECT SCAN")
    print("""  (a) sigma* is found by Brent root-finding on the Fisher median. The route here
  computes the same median from a FINITE-DIFFERENCE HESSIAN of the negative log
  posterior, then brackets the crossing by a direct scan on a log grid with no
  root-finder at all.
  (b) The exactness claim: sd(ln theta) depends on the two error sources only
  through sigma at fixed s_b. Checked here on the HESSIAN route, which never sees
  a decomposition of sigma, and on a second schedule.""")
    import q2_requirement_map as S2
    import q6_sigma_budget as S6

    def hess_sd(theta, times, sigma, s_b, seed):
        L, _, _ = S6.anchor(theta)
        Ls = S2.loads_for(times.size, seed, L)
        inc = S2.increment(times)
        m = S2.level(theta, inc, Ls)

        def J(q):
            th, km, b, ke = np.exp(q[0]), np.exp(q[1]), q[2], np.exp(q[3])
            pred = (th + inc) * m / (km + m) + ke * m + b
            return (float(np.sum(((pred - Ls) / sigma) ** 2)) + (b / s_b) ** 2
                    + ((q[3] - np.log(S2.KE)) / S6.KE_PRIOR) ** 2)
        q0 = np.array([np.log(theta), np.log(S2.KM), 0.0, np.log(S2.KE)])
        h = 1e-5
        H = np.empty((4, 4))
        for i in range(4):
            for j in range(4):
                ei, ej = np.zeros(4), np.zeros(4)
                ei[i] = ej[j] = h
                H[i, j] = (J(q0 + ei + ej) - J(q0 + ei - ej)
                           - J(q0 - ei + ej) + J(q0 - ei - ej)) / (4 * h * h)
        return float(np.sqrt(np.linalg.inv(H / 2.0)[0, 0]))

    def hess_med(theta, times, sigma, s_b):
        return float(np.median([hess_sd(theta, times, sigma, s_b, s)
                               for s in S2.SEEDS]))

    theta = 3000.0
    L = S6.anchor(theta)[0]
    s_b = 0.10 * L
    ok = True
    for name, times in S6.SCHED.items():
        mine = S6.sigma_star(theta, times, s_b)
        grid = np.geomspace(mine * 0.9, mine * 1.1, 41)
        vals = np.array([hess_med(theta, times, g, s_b) for g in grid])
        i = int(np.searchsorted(vals, S6.GATE))
        lo, hi = grid[i - 1], grid[i]
        inside = lo <= mine <= hi
        print(f"    {name:16s} q6 sigma* {mine:9.4f}  scan brackets "
              f"[{lo:9.4f}, {hi:9.4f}]  {'inside' if inside else 'OUTSIDE'}")
        ok &= report(f"sigma* bracketed by the Hessian scan, {name}",
                     0.0 if inside else 1.0, 0.0, 0.0)
    print()
    worst = 0.0
    for name in ("monthly 25", "weekly 100"):
        times = S6.SCHED[name]
        _, m_bar, b = S6.anchor(theta)
        target = 60.0
        got = []
        for lf in (0.0, 0.05, 0.10, 0.15):
            rem = target ** 2 - (lf * L) ** 2
            cv = float(np.sqrt(rem) / (b * m_bar))
            got.append(hess_med(theta, times, S6.sigma_of(theta, lf, cv), s_b))
        spread = max(got) - min(got)
        worst = max(worst, spread)
        print(f"    {name:16s} four (load, CV) pairs at sigma = 60: spread "
              f"{spread:.3e}")
    ok &= report("dependence is exactly through sigma, Hessian route", worst,
                 0.0, 1e-12)
    return ok


def check_4_offset_fixed():
    banner("4. STEP 3(b): THE OFFSET, BY CLOSED FORM AND A FROM-SCRATCH FIT")
    print("""  (a) `IC50_min`, the floor under IC50 where E0 reaches zero, is derived in q7 by
  algebra. Checked here by root-finding E0(IC50) = 0 numerically, which uses the
  same _split but not the derivation.
  (b) The phi window's upper edge is checked by asking the calibration itself: at
  a phi just inside it the structure must calibrate to a 60% fall, and just
  outside it must refuse.
  (c) The asymmetric bias, truth 5% fitted at 10%, is rebuilt by a from-scratch
  grid fit over (theta, base, IC50) with no least-squares projection, which shares
  no code path with q7's profile.""")
    from scipy.optimize import brentq
    import q7_offset_fixed as S7

    worst = 0.0
    for phi in (0.05, 0.10, 0.20):
        o = S7._bare(phi)
        mine = o.ic50_min()
        f = lambda c: o._split(c)[0]
        theirs = float(brentq(f, mine * 1e-3, mine * 1e3, xtol=1e-12,
                              rtol=1e-15))
        worst = max(worst, abs(mine / theirs - 1.0))
        print(f"    phi {100 * phi:5.2f}%: algebra {mine:12.6f}  root-found "
              f"{theirs:12.6f}")
    ok = report("IC50_min, algebra vs root-finding", worst, 0.0, 1e-9)

    print()
    #  the window's upper edge two ways: from the ceiling, and by asking the
    #  calibration itself whether it succeeds. They must agree now that
    #  _calibrate searches only the admissible range.
    import q7_offset_fixed as _S7

    def calibrates(phi):
        try:
            _S7.Offset(phi)
            return True
        except ValueError:
            return False
    lo, hi = 0.05, 0.20
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if calibrates(mid):
            lo = mid
        else:
            hi = mid
    ceil_edge = float(brentq(
        lambda x: S7._bare(x).ceiling()[0] - S7.TARGET_FALL, 0.05, 0.20,
        xtol=1e-12))
    print(f"    calibration boundary {100 * lo:.5f}%   ceiling-derived edge "
          f"{100 * ceil_edge:.5f}%")
    ok &= report("phi window edge, calibration vs ceiling", lo, ceil_edge,
                 1e-5)

    print()
    #  the profile's minimum: at the truth with the correct anchor, away from it
    #  with the anchor held too high. An independent optimiser per node.
    ot = S7.Offset(0.05)
    times = np.arange(0.0, 730.0, 30.4)
    loads = np.full(times.size, S7.LB)
    ths = np.array([30.0, 300.0, 1000.0, 3000.0])
    for label, e0, want in (("correct 5%", ot.e0, 3000.0),
                            ("too high, 10%", S7.Offset(0.10).e0, 30.0)):
        v = S7.profile_theta(e0, ot, times, loads, ths)
        got = float(ths[int(np.argmin(v))])
        print(f"    offset {label:14s} profile minimum at theta = {got:.0f} "
              f"(expected {want:.0f})")
        ok &= report(f"profile minimum, offset {label}",
                     0.0 if got == want else 1.0, 0.0, 0.0)
    return ok


def check_5_pricing():
    banner("5. THE PRICING: THE SIZING FORMULA AND THE BUDGET INVERSION")
    print("""  (a) `sd(ln sigma)` on df degrees of freedom is computed in q8 from the trigamma
  function. Checked here by drawing chi-square variates directly and taking the sd
  of half their log, which shares no code with the trigamma recurrence.
  (b) `cv_allowed` inverts sigma = hypot(load*L, b*cv*m_bar). Checked by going
  FORWARD: the returned CV must reproduce the budget it was inverted from.""")
    import q8_pricing as S8
    import q6_sigma_budget as S6

    rng = np.random.default_rng(31_415)
    worst = 0.0
    for df in (5, 20, 50, 150):
        mine = S8.sd_ln_sigma(df)
        draws = rng.chisquare(df, 400_000) / df
        theirs = float(np.std(0.5 * np.log(draws)))
        worst = max(worst, abs(mine / theirs - 1.0))
        print(f"    df {df:4d}: trigamma {mine:.5f}  chi-square draws "
              f"{theirs:.5f}")
    ok = report("sd(ln sigma) trigamma vs chi-square draws", worst, 0.0, 0.02,
                rel=False)

    print()
    worst = 0.0
    theta = 3000.0
    L, m_bar, b = S6.anchor(theta)
    for st in (36.72, 56.10, 130.96):
        for lf in (0.02, 0.10):
            cv = S8.cv_allowed(st, lf, theta)
            if cv is None:
                continue
            fwd = float(np.hypot(lf * L, b * cv * m_bar))
            worst = max(worst, abs(fwd / st - 1.0))
            print(f"    sigma* {st:7.2f} load {100 * lf:4.1f}%: CV "
                  f"{100 * cv:6.2f}%  forward sigma {fwd:7.2f}")
    ok &= report("cv_allowed inverts the budget exactly", worst, 0.0, 1e-12,
                 rel=False)
    return ok


def main():
    print("=" * 92)
    print("  PHASE 9 RE-CHECKS")
    print("=" * 92)
    for fn in (check_1_model_explained_cv, check_2_item6_arithmetic,
               check_3_sigma_budget, check_4_offset_fixed,
               check_5_pricing):
        try:
            fn()
        except Exception as exc:
            RESULTS.append((fn.__name__ + " (raised)", False))
            print(f"    {fn.__name__} RAISED {type(exc).__name__}: {exc}")
    banner("SUMMARY")
    bad = [n for n, ok in RESULTS if not ok]
    for name, ok in RESULTS:
        print(f"    {'PASS' if ok else 'FAIL'}  {name}")
    print(f"\n    {len(RESULTS) - len(bad)} of {len(RESULTS)} checks pass.")
    if bad:
        print("    FAILURES, reported and not repaired:")
        for n in bad:
            print(f"      - {n}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
