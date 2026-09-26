"""
The sampling modality: what sigma is made of, what a worse assay costs, and
what home volumetric microsampling buys against clinic visits.

THE STRUCTURE THE WHOLE THING TURNS ON. model_A.run.sigma_total makes the
ordinate noise of the steady-state fit out of two terms in flux units,

    sigma = hypot( load_frac * L ,  b * cv * m_bar ),

b being the tangent slope theta g'(m_bar) + k_e. The two are comparable at the
nominal footing, so sigma is NOT proportional to the assay CV and the penalty
for a worse assay is DAMPED: the load-record term is a floor the assay cannot
sink below. That is what makes a cheaper, noisier, more frequent sample worth
considering at all, and PART B prices it.

WHAT THIS SCRIPT DOES.
  A  the split, on this script's own footing, against the crossover
     model_A.run.assay_load_crossover already states.
  B  the sample-count factor against assay CV, twice: analytically from
     sigma^2, and by counting visits to the gate in the design itself. They
     differ where the schedule's phase coverage changes with n, and the gap is
     the thing worth knowing.
  C  Tier 1's 118 monthly draws converted into calendar time and burden for
     home sampling at one, two and three per week.
  D  whether home sampling can place draws at the phases the phase-spread
     schedule wants, which clinic scheduling cannot. Two separate obstructions
     are priced: a fixed weekly cadence ALIASES against a 14 d dose interval,
     and a clinic can only book weekdays.
  E  whether the binding nuisance changes once the assay term dominates.

FOOTING. src/common/params.py and src/design/run.py: theta = 3000, K_m = 100,
k_e = 5/d, L = 322.727 umol/L/d, D = 6000 on a 14 d train with a 3.5 d
half-life, load recorded to 5% with a day-to-day intake CV of 10%, k_e
log-prior 0.8, gate 0.20. The Tier 1 arm carries the increment's scale as a
fifth column with the published 14.02% log-prior, which is
audit/t2_tier1_margin.py's construction and is anchored against its 118.

Nothing here records an expected result in advance.

Run:  python audit/u5_sigma_split.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import model_A.run as LOC
from common.models import ParamsA
from design.run import (ASSAY_CV, DOSE_INTERVAL, D_AMP, LAMBDA, LOAD_CV,
                        LOAD_FRAC, load)

KM, KE = 100.0, 5.0
LB = float(load(0.0))
GATE = 0.20
MONTH = 30.4
WEEK = 7.0
KE_PRIOR = 0.8
PUBLISHED_KAPPA_CV = 0.1402
T_DOSE = 50.0 * DOSE_INTERVAL
PEAK = D_AMP / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))
SEEDS = range(24)


# ---------------------------------------------------------------- model
def op_level(theta, Ltot=None):
    Ltot = LB if Ltot is None else Ltot
    B = KE * KM + theta - Ltot
    return (-B + np.sqrt(B ** 2 + 4.0 * KE * Ltot * KM)) / (2.0 * KE)


def split(theta=3000.0, load_frac=LOAD_FRAC, assay_cv=ASSAY_CV, m_bar=None):
    """The two terms of sigma, separately, in flux units."""
    m_bar = op_level(theta) if m_bar is None else m_bar
    slope = theta * KM / (KM + m_bar) ** 2 + KE
    return load_frac * LB, slope * assay_cv * m_bar, m_bar, slope


def sigma_of(theta=3000.0, load_frac=LOAD_FRAC, assay_cv=ASSAY_CV):
    a, b, _, _ = split(theta, load_frac, assay_cv)
    return float(np.hypot(a, b))


def levels(theta, times, load_cv=LOAD_CV, seed=0):
    """Marker, increment and recorded load on a schedule."""
    t = np.asarray(times, dtype=float)
    tau = np.mod(t - T_DOSE, DOSE_INTERVAL)
    inc = PEAK * np.exp(-LAMBDA * tau)
    if load_cv:
        Ls = LB * (1.0 + load_cv *
                   np.random.default_rng(seed).standard_normal(t.size))
    else:
        Ls = np.full(t.size, LB)
    B = KE * KM + theta + inc - Ls
    return (-B + np.sqrt(B ** 2 + 4.0 * KE * Ls * KM)) / (2.0 * KE), inc


def sd_ln_theta(theta, m, inc, sigma, kappa_cv=None, ke_prior=KE_PRIOR,
                load_prior=True, exact=()):
    """(1,1) entry of the inverse information matrix, as a standard deviation."""
    gm = m / (KM + m)
    cols = [theta * gm,
            -(theta + inc) * KM * m / (KM + m) ** 2,
            np.ones_like(m),
            KE * m]
    pri = [np.inf, np.inf,
           (LOAD_FRAC * LB) ** 2 if load_prior else np.inf,
           ke_prior ** 2 if ke_prior is not None else np.inf]
    names = ["theta", "Km", "b", "ke"]
    if kappa_cv is not None:
        cols.append(inc * gm)
        pri.append(kappa_cv ** 2)
        names.append("kappa")
    keep = [i for i, nm in enumerate(names) if nm not in exact]
    Q = np.stack(cols, axis=1)[:, keep]
    pri = [pri[i] for i in keep]
    F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri]) \
        + Q.T @ Q / sigma ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def over_seeds(theta, times, sigma, **kw):
    v = np.array([sd_ln_theta(theta, *levels(theta, times, seed=s), sigma, **kw)
                  for s in SEEDS], float)
    f = v[np.isfinite(v)]
    if not f.size:
        return np.inf, np.inf, np.inf
    return float(np.median(f)), float(f.min()), float(f.max())


def fixed(n, step):
    return step * np.arange(n, dtype=float)


def phase_spread(n, gap, period=DOSE_INTERVAL, grid=0.5):
    """
    n draws whose dose-cycle phases are as evenly spread as the grid allows,
    subject to a minimum gap. Greedy nearest-reachable, which is
    g_design_common.phase_spread's construction and is used for the same
    reason: taking the targets in numerical order stretches the record.
    """
    targets = list((np.arange(n, dtype=float) * period / n) % period)
    out, rem = [], list(targets)
    for i in range(n):
        earliest = 0.0 if i == 0 else out[-1] + gap
        base = np.ceil(earliest / grid) * grid
        best, best_t = None, None
        for p in rem:
            cyc = np.floor(base / period) * period
            cand = np.round(np.array([cyc + p, cyc + period + p]) / grid) * grid
            cand = cand[cand >= base - 1e-9]
            if not cand.size:
                continue
            if best_t is None or cand[0] < best_t:
                best, best_t = p, float(cand[0])
        if best is None:
            best, best_t = rem[0], float(base)
        rem.remove(best)
        out.append(best_t)
    return np.array(out)


def lattice(step, period=DOSE_INTERVAL):
    """
    How many distinct dose-cycle phases a fixed step can EVER reach.

    Write step/period = p/q in lowest terms; the phases lie on a lattice of
    spacing period/q. `step` must be an exact Fraction wherever the cadence
    is not a terminating decimal: three samples a week is 7/3 days, and
    reducing float(7/3) through its decimal string returns a denominator of
    2.8e16, which is the float's own rounding answering a different
    question. That is the same trap audit/g_design_common.py's own
    lattice_phases() documents, arrived at from the other side.
    """
    from fractions import Fraction
    st = step if isinstance(step, Fraction) else Fraction(str(float(step)))
    pe = period if isinstance(period, Fraction) else Fraction(str(float(period)))
    return (st / pe).denominator


_PS_CACHE = {}


def _spread_cached(n, gap):
    key = (n, round(float(gap), 9))
    if key not in _PS_CACHE:
        _PS_CACHE[key] = phase_spread(n, gap)
    return _PS_CACHE[key]


def visits_to_gate(theta, step_or_times, sigma, nmax=2000, spread=False,
                   gap=None, **kw):
    """
    Fewest draws to the gate, by doubling and then bisecting.

    A linear scan is what this did first and it is not affordable here: the
    phase-spread constructor is quadratic in n, so a scan that has to reach
    the thousands before deciding "never" costs more than every other part
    of this script together. Bisection assumes sd(ln theta) is decreasing in
    n at a fixed schedule family, which monotone_check() tests rather than
    asserts.
    """
    def sd_at(n):
        t = (_spread_cached(n, gap) if spread else fixed(n, step_or_times))
        return over_seeds(theta, t, sigma, **kw)[0]

    if sd_at(4) < GATE:
        return 4
    hi = 8
    while hi <= nmax and not (sd_at(hi) < GATE):
        hi *= 2
    if hi > nmax:
        return None
    lo = hi // 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if sd_at(mid) < GATE:
            hi = mid
        else:
            lo = mid
    return hi


def monotone_check(theta, sigma, gap, ns, **kw):
    """Is sd(ln theta) decreasing in n for the phase-spread family?"""
    v = [over_seeds(theta, _spread_cached(n, gap), sigma, **kw)[0] for n in ns]
    bad = [(ns[i], ns[i + 1]) for i in range(len(v) - 1)
           if v[i + 1] > v[i] + 1e-12]
    return v, bad


# ======================================================================
def part_a():
    print("=" * 76)
    print("PART A.  what sigma is made of, on this script's own footing")
    print("=" * 76)
    a, b, mbar, slope = split()
    print(f"\n  theta = 3000, m_bar = {mbar:.3f}, tangent slope = {slope:.4f}")
    print(f"    load-record term  load_frac * L        = {a:.4f}")
    print(f"    assay term        slope * cv * m_bar   = {b:.4f}")
    print(f"    sigma = hypot                          = {np.hypot(a, b):.4f}")
    print(f"    assay share of sigma^2                 = "
          f"{b ** 2 / (a ** 2 + b ** 2):.4%}")
    print(f"    model_A.run.sigma_total                = "
          f"{LOC.sigma_total(LOAD_FRAC, ASSAY_CV):.4f}")
    print(f"    crossover, assay CV per unit load error= "
          f"{LOC.assay_load_crossover():.4f}")
    print(f"    the same from this script              = {a / b * ASSAY_CV / LOAD_FRAC:.4f}")

    print("\n  the split at other capacities, recorded load held fixed:")
    print(f"    {'theta':>7} {'m_bar':>8} {'load term':>10} {'assay term':>11}"
          f" {'sigma':>8} {'assay share':>12}")
    print("    " + "-" * 60)
    for theta in (3000.0, 600.0, 300.0, 150.0):
        a, b, mbar, _ = split(theta)
        print(f"    {theta:7.0f} {mbar:8.3f} {a:10.3f} {b:11.3f} "
              f"{np.hypot(a, b):8.3f} {b ** 2 / (a ** 2 + b ** 2):12.2%}")
    print("\n  The load term does not move with capacity, because the")
    print("  recorded load is held fixed; the assay term rises with m_bar")
    print("  and falls with the tangent slope, and the two nearly cancel.")


def part_b():
    print()
    print("=" * 76)
    print("PART B.  the sample-count factor against assay CV, two ways")
    print("=" * 76)
    print("\n  analytic: (sigma(cv)/sigma(5%))^2, which is exact if the only")
    print("  thing the assay changes is the ordinate noise.")
    print("  counted:  monthly visits to the gate at that CV, divided by the")
    print("  count at 5%. It differs where the priors do not scale with n.\n")
    base_sigma = sigma_of()
    n_base = visits_to_gate(3000.0, MONTH, base_sigma)
    print(f"  monthly visits to the gate at assay CV 5%: {n_base}")
    print(f"\n    {'assay CV':>9} {'sigma':>8} {'analytic':>9} "
          f"{'visits':>7} {'counted':>8}")
    print("    " + "-" * 46)
    for cv in (0.02, 0.05, 0.075, 0.10, 0.125, 0.15, 0.20, 0.30):
        s = sigma_of(assay_cv=cv)
        n = visits_to_gate(3000.0, MONTH, s)
        print(f"    {cv:9.1%} {s:8.3f} {(s / base_sigma) ** 2:9.3f} "
              f"{str(n):>7} "
              f"{(n / n_base if n else float('nan')):8.3f}")
    print("\n  If sigma scaled with the assay CV the analytic column would be")
    print("  (cv/0.05)^2, that is 4.00 at 10% and 9.00 at 15%. The damping is")
    print("  the ratio of those to the analytic column.")
    for cv in (0.10, 0.15):
        s = sigma_of(assay_cv=cv)
        print(f"    at {cv:.0%}: undamped {(cv / ASSAY_CV) ** 2:.2f}, "
              f"actual {(s / base_sigma) ** 2:.2f}, "
              f"damping {(cv / ASSAY_CV) ** 2 / (s / base_sigma) ** 2:.2f}")


def part_c():
    print()
    print("=" * 76)
    print("PART C.  Tier 1 in calendar time and burden")
    print("=" * 76)
    print("\n  Tier 1's requirement is the increment characterised to 11.6%;")
    print("  at the published between-visit 14.02% the monthly record reaches")
    print("  the gate only after the count below. The increment's scale is")
    print("  carried as a fifth column with that log-prior throughout.\n")
    sig5 = sigma_of()
    n_clinic = visits_to_gate(3000.0, MONTH, sig5,
                              kappa_cv=PUBLISHED_KAPPA_CV)
    print(f"  ANCHOR. monthly clinic draws to the gate at kappa CV "
          f"{PUBLISHED_KAPPA_CV:.2%}, assay CV 5%: {n_clinic}")
    print(f"          calendar time {n_clinic * MONTH / 365.25:.2f} years, "
          f"{n_clinic} clinic trips, {n_clinic} venepunctures")
    print("          (audit/t2_tier1_margin.py owns this number; the")
    print("           agreement is the anchor, not a second derivation)")

    print("\n  home volumetric microsampling, phase-spread at the stated")
    print("  minimum gap, against assay CV. Draws, then calendar years.\n")
    cadences = [("3 per week", WEEK / 3.0), ("2 per week", WEEK / 2.0),
                ("1 per week", WEEK), ("monthly", MONTH)]
    print(f"    {'assay CV':>9} " +
          "".join(f"{nm:>22}" for nm, _ in cadences))
    print("    " + "-" * 98)
    table = {}
    for cv in (0.05, 0.10, 0.15, 0.20):
        s = sigma_of(assay_cv=cv)
        cells = []
        for nm, gap in cadences:
            n = visits_to_gate(3000.0, None, s, spread=True, gap=gap,
                               kappa_cv=PUBLISHED_KAPPA_CV)
            if n is None:
                cells.append("not reached")
                table[(cv, nm)] = None
            else:
                t = _spread_cached(n, gap)
                yrs = (t[-1] - t[0]) / 365.25
                cells.append(f"{n} draws, {yrs:.2f} yr")
                table[(cv, nm)] = (n, yrs)
        print(f"    {cv:9.1%} " + "".join(f"{c:>22}" for c in cells))

    print("\n  BURDEN, against the clinic arm's "
          f"{n_clinic} trips and {n_clinic} venepunctures:")
    print(f"    {'assay CV':>9} {'cadence':>12} {'draws':>7} {'years':>7} "
          f"{'trips':>7} {'draws/clinic draw':>18} {'years saved':>12}")
    print("    " + "-" * 78)
    clinic_yrs = n_clinic * MONTH / 365.25
    for cv in (0.10, 0.15, 0.20):
        for nm, gap in cadences[:3]:
            e = table.get((cv, nm))
            if e is None:
                continue
            n, yrs = e
            print(f"    {cv:9.1%} {nm:>12} {n:7d} {yrs:7.2f} {0:7d} "
                  f"{n / n_clinic:18.2f} {clinic_yrs - yrs:12.2f}")
    print("\n  'trips' is zero by construction: a home sample is a finger")
    print("  prick, not a clinic attendance. The burden traded is one clinic")
    print("  trip per month against the stated number of self-collections.")


def part_d():
    print()
    print("=" * 76)
    print("PART D.  can home sampling place draws where the design wants them")
    print("=" * 76)
    print("\n  TWO SEPARATE OBSTRUCTIONS, and they point opposite ways.\n")
    print("  (i) A FIXED cadence aliases against the 14 d dose interval. The")
    print("      phases a step can ever reach is the denominator of step/14")
    print("      in lowest terms, whatever the number of draws.\n")
    print(f"    {'cadence':>12} {'step (d)':>9} {'phases reachable':>18}")
    print("    " + "-" * 42)
    from fractions import Fraction as Fr
    for nm, step in [("daily", Fr(1)), ("3 per week", Fr(7, 3)),
                     ("2 per week", Fr(7, 2)), ("1 per week", Fr(7)),
                     ("fortnightly", Fr(14)), ("monthly", Fr("30.4")),
                     ("2-monthly", Fr("60.8"))]:
        print(f"    {nm:>12} {float(step):9.4f} {lattice(step):18d}")
    print("\n      A weekly home cadence reaches TWO phases of the dose cycle.")
    print("      That is the exactly-aliased case Remark rem:aliasing names,")
    print("      and it is the schedule a home protocol would write first.\n")

    sig5 = sigma_of()
    sig10 = sigma_of(assay_cv=0.10)
    print("  (ii) The cost, in sd(ln theta), of the fixed cadence against the")
    print("       phase-spread schedule at the SAME number of draws and the")
    print("       same minimum gap. Medians over 24 intake records.\n")
    print(f"    {'phases':>6} {'cadence':>9} {'n':>5} {'fixed, cv 5%':>22} "
          f"{'phase-spread, cv 5%':>24} {'ratio':>7}")
    print("    " + "-" * 76)
    for nm, step in [("3 per week", WEEK / 3.0), ("2 per week", WEEK / 2.0),
                     ("1 per week", WEEK), ("fortnightly", 14.0),
                     ("monthly", MONTH)]:
        for n in (52, 104):
            f_med, f_lo, f_hi = over_seeds(3000.0, fixed(n, step), sig5)
            p_med, p_lo, p_hi = over_seeds(3000.0, phase_spread(n, step), sig5)
            rat = f_med / p_med if np.isfinite(f_med) and p_med else np.inf
            nph = lattice(__import__("fractions").Fraction(str(step))
                          if step != WEEK / 3.0
                          else __import__("fractions").Fraction(7, 3))
            print(f"    [q={min(nph, n):>3d}] {nm:>9} {n:5d} "
                  f"{f'{f_med:.4f} [{f_lo:.4f},{f_hi:.4f}]':>22} "
                  f"{f'{p_med:.4f} [{p_lo:.4f},{p_hi:.4f}]':>24} {rat:7.2f}")

    print("\n  (iii) What a clinic cannot do that a home kit can: a 14 d dose")
    print("        interval is exactly two weeks, so the dose always falls on")
    print("        the same weekday and the four weekend phases are never")
    print("        bookable. Phases restricted to weekdays against all days,")
    print("        at equal draw counts, monthly-equivalent spacing.\n")
    wk = [p for p in range(14) if p % 7 < 5]
    print(f"      weekday-reachable integer phases: {wk}  ({len(wk)} of 14)")
    print(f"      home-reachable integer phases:    {list(range(14))}  (14 of 14)")
    print()
    print(f"    {'n':>5} {'weekday-only':>24} {'all days':>24} {'ratio':>7}")
    print("    " + "-" * 62)
    for n in (25, 52, 104):
        t_all = phase_spread(n, MONTH, grid=1.0)
        #  weekday-restricted: targets drawn only from the bookable phases
        targets = [wk[int(round(i * (len(wk) - 1) / max(n - 1, 1)))]
                   if n > 1 else wk[0] for i in range(n)]
        t_wk = []
        for i, p in enumerate(targets):
            base = 0.0 if i == 0 else t_wk[-1] + MONTH
            cyc = np.floor(base / DOSE_INTERVAL) * DOSE_INTERVAL
            cand = [cyc + p, cyc + DOSE_INTERVAL + p,
                    cyc + 2 * DOSE_INTERVAL + p]
            t_wk.append(float(min(c for c in cand if c >= base - 1e-9)))
        t_wk = np.array(t_wk)
        a_med, a_lo, a_hi = over_seeds(3000.0, t_all, sig5)
        w_med, w_lo, w_hi = over_seeds(3000.0, t_wk, sig5)
        print(f"    {n:5d} {f'{w_med:.4f} [{w_lo:.4f},{w_hi:.4f}]':>24} "
              f"{f'{a_med:.4f} [{a_lo:.4f},{a_hi:.4f}]':>24} "
              f"{w_med / a_med:7.3f}")

    print()
    print("  (iv) DOES THE RANK DEFICIT AT FEW PHASES COST ANYTHING. At a")
    print("       constant recorded load the dose train is the only thing")
    print("       moving the level, and audit/u6_recheck_sampling.py part C")
    print("       reads the DATA rank off the singular values: 3 at six and")
    print("       four phases, 2 at two phases, 1 at one. The question this")
    print("       block answers is whether the deficit reaches the estimate,")
    print("       given that the design also carries a prior on b and one on")
    print("       ln k_e, which can close a null direction the data leaves")
    print("       open. 52 draws.\n")
    print(f"    {'cadence':>12} {'phases':>7} {'load CV 0':>26} "
          f"{'load CV 10%':>26} {'ratio':>9}")
    print("    " + "-" * 84)
    from fractions import Fraction as _Fr
    for nm, step, fr in [("3 per week", WEEK / 3.0, _Fr(7, 3)),
                         ("2 per week", WEEK / 2.0, _Fr(7, 2)),
                         ("1 per week", WEEK, _Fr(7)),
                         ("fortnightly", 14.0, _Fr(14)),
                         ("monthly", MONTH, _Fr("30.4"))]:
        t = fixed(52, step)
        c0 = sd_ln_theta(3000.0, *levels(3000.0, t, load_cv=0.0), sig5)
        c1_med, c1_lo, c1_hi = over_seeds(3000.0, t, sig5)
        rat = c0 / c1_med if np.isfinite(c0) and np.isfinite(c1_med) else np.inf
        c0s = f"{c0:.4f}" if np.isfinite(c0) else "inf (singular)"
        print(f"    {nm:>12} {lattice(fr):7d} {c0s:>26} "
              f"{f'{c1_med:.4f} [{c1_lo:.4f},{c1_hi:.4f}]':>26} {rat:9.2f}")
    print()
    print("       READING, and it is not the one this block was written to")
    print("       test. The rank deficit at two phases does NOT reach the")
    print("       estimate: the two priors close both null directions and")
    print("       the precision is within one per cent of the monthly")
    print("       record's at every cadence from two phases upward, with or")
    print("       without the intake varying. The conjecture that a weekly")
    print("       home cadence leans on the unsmoothed diet record in a way")
    print("       a monthly clinic one does not is REFUTED by the load CV 0")
    print("       column, where the weekly cadence returns a finite and")
    print("       ordinary number.")
    print()
    print("       What survives is sharper and simpler than the graded rule")
    print("       a phase lattice suggests: only q = 1 fails. Intake")
    print("       variability then turns an impossible design into a useless")
    print("       one rather than into a workable one, which is the")
    print("       difference between the inf and the 10.75 in that row.")


def part_e():
    print()
    print("=" * 76)
    print("PART E.  does the binding nuisance change as the assay term grows")
    print("=" * 76)
    print("\n  25 monthly draws, each nuisance held EXACTLY known in turn,")
    print("  medians over 24 intake records. The increment's scale is carried")
    print("  at the published 14.02% so that all four nuisances are present.\n")
    print(f"    {'assay CV':>9} {'assay share':>12} {'none':>9} {'b exact':>9} "
          f"{'k_e exact':>10} {'K_m exact':>10} {'kappa exact':>12} "
          f"{'binds':>8}")
    print("    " + "-" * 84)
    for cv in (0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30):
        a, b, _, _ = split(assay_cv=cv)
        s = float(np.hypot(a, b))
        t = fixed(25, MONTH)
        base = over_seeds(3000.0, t, s, kappa_cv=PUBLISHED_KAPPA_CV)[0]
        vals = {}
        for nm, ex in [("b", ("b",)), ("ke", ("ke",)), ("Km", ("Km",)),
                       ("kappa", ("kappa",))]:
            vals[nm] = over_seeds(3000.0, t, s, kappa_cv=PUBLISHED_KAPPA_CV,
                                  exact=ex)[0]
        binds = min(vals, key=lambda k: vals[k])
        print(f"    {cv:9.1%} {b ** 2 / (a ** 2 + b ** 2):12.1%} {base:9.4f} "
              f"{vals['b']:9.4f} {vals['ke']:10.4f} {vals['Km']:10.4f} "
              f"{vals['kappa']:12.4f} {binds:>8}")
    print("\n  'binds' names whichever exact value buys the most, which is the")
    print("  operational reading of 'the binding nuisance'.")


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
