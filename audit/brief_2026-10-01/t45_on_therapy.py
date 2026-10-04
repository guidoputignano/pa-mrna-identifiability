"""
Round 12: what a two-year ON-THERAPY record buys, and whether it self-calibrates
lambda. A retrospective-adjustment question, NOT the dual-control problem of sec 7:
every patient is booked at the protocol's 14 d (Model C says that is safe for all,
t44), so nothing is risked; we simply read the accumulated record.

The therapy's dose schedule is EXOGENOUS (fixed by protocol, not chosen by the
estimator) and DATED. By Prop 2 (tworoutes) the enzymatic channel does not matter
for design -- only the excursion the increment produces. The increment Delta(t)
drives the operating point m over a range as it rises after each dose and decays
(half-life 3.5 d) before the next; dose changes, missed doses, supply gaps and
interval variation are known, dated perturbations that widen that excursion. A
wider m-excursion samples g(m)=m/(K_m+m) over a wider range, which is exactly what
breaks the theta--k_e confound that flattens the (near-aliased) pre-therapy record.

The marker's clearance ordinate is  y = (theta+Delta(t))g(m) + k_e m + b, with
Delta(t) = sum_doses amp*exp(-lambda*(t - t_dose)). Sensitivities at fixed observed m:
    d y/d ln theta = theta g(m)              (Delta does not depend on theta)
    d y/d ln K_m   = -(theta+Delta) K_m m/(K_m+m)^2
    d y/d b        = 1
    d y/d ln k_e   = k_e m
    d y/d ln lambda= g(m) * d Delta/d ln lambda,
        d Delta/d ln lambda = -lambda * sum_doses amp*(t-t_dose)*exp(-lambda*(t-t_dose))
lambda is identified from the SAME record because the marker's response to a longer
gap (a missed dose / lengthened interval) depends on lambda directly.
"""
import sys; sys.path.insert(0, "../../src")
import numpy as np

import core as C
from model_A.run import ParamsA, M_TOL, baseline_load
from common.models import ParamsC
import model_C.run as MC

LAM = np.log(2.0) / 3.5
D14 = 6000.0
TAU0 = 14.0
KM, KE = 100.0, 5.0
TH5 = C.TH_5                       # 150 = 5% of normal
MONTHLY = C.N4.MONTHLY            # 25 monthly draws over 2 years
Z1 = 1.6449                       # one-sided 95%


def schedule(seed=0, years=2.0, tau=TAU0, jitter=0.0, miss_p=0.0, gap=None,
             same_rate=True):
    """Dated dose schedule over `years`. jitter: sd of interval (d); miss_p:
    per-dose skip probability; gap: (start_d, length_d) supply interruption."""
    rng = np.random.default_rng(1000 + seed)
    T = years * 365.0
    amp = D14 * (tau / TAU0) if same_rate else D14
    t, times = 0.0, []
    while t < T:
        times.append(t)
        step = tau + (rng.normal(0.0, jitter) if jitter else 0.0)
        t += max(2.0, step)
    times = np.array(times)
    amps = np.full(times.size, amp)
    if miss_p:
        keep = rng.random(times.size) > miss_p
        keep[0] = True
        times, amps = times[keep], amps[keep]
    if gap:
        g0, gl = gap
        drop = (times >= g0) & (times < g0 + gl)
        times, amps = times[~drop], amps[~drop]
    return times, amps


def inc_dinc(t, dose_times, dose_amps, lam=LAM):
    """Delta(t) and d Delta/d ln lambda at times t (both arrays)."""
    t = np.asarray(t, float)
    inc = np.zeros_like(t)
    dinc_dlam = np.zeros_like(t)            # d Delta / d lambda
    for td, a in zip(dose_times, dose_amps):
        on = t >= td
        age = t[on] - td
        e = a * np.exp(-lam * age)
        inc[on] += e
        dinc_dlam[on] += -age * e
    return inc, dinc_dlam * lam             # -> d Delta / d ln lambda


def on_therapy_record(theta, seed=0, sched_kw=None, times=None, load_cv=0.10):
    """One on-therapy record: monthly draws, marker driven by the dated schedule."""
    sched_kw = sched_kw or {}
    dts, das = schedule(seed=seed, **sched_kw)
    times = MONTHLY if times is None else np.asarray(times, float)
    L, m_bar = C.N4.anchor(theta, KM, KE, "A")
    rng = np.random.default_rng(seed)
    if load_cv:
        s2 = np.log1p(load_cv ** 2)
        loads = L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
    else:
        loads = np.full(times.size, L)
    inc, dinc = inc_dinc(times, dts, das)
    m = C.level(theta, KM, KE, loads, inc)
    return dict(theta=theta, Km=KM, ke=KE, L=L, m_bar=m_bar, loads=loads,
                inc=inc, dinc=dinc, m=m, times=times)


def fim5(r, with_lambda=True, ke_prior=0.8, lam_prior=None):
    """FIM in (ln theta, ln K_m, b, ln k_e[, ln lambda]) on the clearance ordinate,
    per-draw noise. lam_prior: sd of a prior on ln lambda (None = record only)."""
    th, K, ke, m, inc = r["theta"], r["Km"], r["ke"], r["m"], r["inc"]
    g = m / (K + m)
    cols = [th * g,
            -(th + inc) * K * m / (K + m) ** 2,
            np.ones_like(m),
            ke * m]
    if with_lambda:
        cols.append(g * r["dinc"])         # d y / d ln lambda
    Qt = np.stack(cols, axis=1)
    sg = C.sigmas(r, "perdraw")
    F = (Qt / sg[:, None] ** 2).T @ Qt
    n = Qt.shape[1]
    pri = np.zeros(n)
    pri[2] = (C.LOAD_FRAC * r["L"]) ** -2          # load-bias prior
    pri[3] = ke_prior ** -2                         # k_e prior
    if with_lambda and lam_prior:
        pri[4] = lam_prior ** -2
    return F + np.diag(pri)


def sd_of(F, i):
    w = np.linalg.eigvalsh(F)
    if w[0] <= 1e-12 * w[-1]:
        return np.inf
    return float(np.sqrt(np.linalg.inv(F)[i, i]))


def med_sd(theta, i=0, with_lambda=True, sched_kw=None, lam_prior=None,
           times=None, seeds=range(24)):
    v = [sd_of(fim5(on_therapy_record(theta, seed=s, sched_kw=sched_kw, times=times),
                    with_lambda=with_lambda, lam_prior=lam_prior), i)
         for s in seeds]
    v = [x for x in v if np.isfinite(x)]
    return float(np.median(v)) if v else np.inf


# ---- theta_crit and trough, from t44 (for the determined-subgroup / lengthening)
def theta_crit_A(fold):
    return (fold * baseline_load(ParamsA()) - KE * M_TOL) * (KM + M_TOL) / M_TOL

THC_INT = MC.theta_crit_intrinsic(ParamsC(), MC.S_BASE)   # 150.8


def trough(tau, lam=LAM):
    x = np.exp(-lam * tau); D = D14 * (tau / TAU0)
    return D * x / (1.0 - x)


EVENTS = dict(jitter=2.0, miss_p=0.08, gap=(300.0, 28.0))   # a realistic 2-yr record
FORTNIGHTLY = np.linspace(0.0, 730.0, 52)
WEEKLY = np.linspace(0.0, 730.0, 104)
ALTDAY = np.linspace(0.0, 730.0, 365)
DENSITIES = [("monthly (25)", MONTHLY), ("fortnightly (52)", FORTNIGHTLY),
             ("weekly (104)", WEEKLY), ("alternate-day (365)", ALTDAY)]


def main():
    print("=" * 92)
    print("  ROUND 12. TWO YEARS ON THERAPY: RETROSPECTIVE ADJUSTMENT ON THE ACCUMULATED RECORD")
    print("=" * 92)

    # --- premise correction: the 3.05 is ALREADY on-therapy
    flat = C.sd_theta(TH5, rec_kw=dict(D=0.0), noise="perdraw")[0]
    reg = C.sd_theta(TH5, noise="perdraw")[0]
    print(f"""
  ITEM 1. WHAT TWO YEARS OF DOSING BUYS (sd(ln theta), median over 24 intake draws)

  PREMISE CORRECTION: the 3.05 baseline is NOT pre-therapy -- it is ALREADY the regular
  on-therapy monthly record. The genuinely pre-therapy (flat, no increment) record is
  {flat:.2f}; the regular 14-day dose train sampled monthly is {reg:.2f}. So most of the
  lift ({flat:.1f} -> {reg:.1f}) is bought by the dose train ITSELF -- the increment's
  phase-drifting excursion against monthly sampling is already Prop-2 excitation. The
  question is what the DATED EVENTS and more sampling add on top.
""")
    print(f"    (a) schedule richness, monthly sampling:")
    for lab, kw in [("regular 14 d", dict()), ("+ interval jitter (sd 2 d)", dict(jitter=2.0)),
                    ("+ jitter + 8% missed doses", dict(jitter=2.0, miss_p=0.08)),
                    ("+ jitter + misses + 28 d supply gap", EVENTS)]:
        sd = med_sd(TH5, i=0, with_lambda=False, sched_kw=kw)
        print(f"        {lab:>40}  sd(ln theta) = {sd:.2f}")

    print(f"\n    (b) sampling density, with the full event schedule, and the determined subgroup")
    print(f"        (grid {{7,10,14}}, routine-load theta_crit = {theta_crit_A(1.0):.0f}):")
    print(f"        {'sampling (draws/2 yr)':>22}  {'sd(ln theta)':>12}  {'correct':>8}  {'determined':>11}")
    dens_sd = {}
    for lab, tm in DENSITIES:
        sd = med_sd(TH5, i=0, with_lambda=False, sched_kw=EVENTS, times=tm)
        dens_sd[lab] = sd
        acc, det = determined(sd)
        print(f"        {lab:>22}  {sd:>12.2f}  {acc:>8.0%}  {det:>11.0%}")

    print(f"""
  READING: on therapy the increment drives the marker over an ~8-10x range between dose
  and trough; the monthly draws sample it at drifting phase, and dated events (jitter,
  misses, a supply gap) widen it. This separates theta from k_e -- Prop 2, at no risk,
  because the booked interval is safe for all (t44). But the honest answer to "can any
  patient's interval be changed on evidence after two years" is mostly NO at realistic
  cadence: sd(ln theta) is sampling-limited (~1/sqrt(n): monthly 2.3, weekly 1.3,
  alternate-day 0.75) and never reaches the ~0.29 that would individualise the whole
  cohort, so the DETERMINED SUBGROUP is 0% through weekly sampling and reaches only ~27%
  at alternate-day (365 draws). Two years of passive observation tightens theta
  meaningfully (8.6 -> 2.3, better with density) but does not make the {{7,10,14}} choice
  individually determined for most. (The safe, evidence-gated lever is lengthening, which
  is one-sided and far more favourable -- item 3.)
""")

    # ---- ITEM 2 ----
    print("=" * 92)
    print("  ITEM 2. IS lambda RECOVERABLE FROM THE SAME RECORD?  (it is -- the clear win)")
    print("=" * 92)
    print(f"\n    Event schedule; sd(ln lambda), and the cost to theta of freeing lambda:")
    print(f"        {'sampling':>22}  {'sd(ln lambda)':>13}  {'sd(ln th)|lam fixed':>20}  {'sd(ln th)|lam free':>19}")
    for lab, tm in DENSITIES:
        sdlam = med_sd(TH5, i=4, with_lambda=True, sched_kw=EVENTS, times=tm)
        sdth_free = med_sd(TH5, i=0, with_lambda=True, sched_kw=EVENTS, times=tm)
        print(f"        {lab:>22}  {sdlam:>13.3f}  {dens_sd[lab]:>20.2f}  {sdth_free:>19.2f}")
    print(f"""
  READING: the last round's lambda SPECIFICATION was ~5% (routine-load coarse grid) to
  ~25% (Model C headline). sd(ln lambda) is the relative precision on lambda, and the
  monthly record already delivers ~0.055 (5.5%) -- comfortably inside the 25% spec and
  essentially at the 5% spec, tightening to ~0.02 with dense sampling. So the accumulated
  dosing record SELF-CALIBRATES the one constant the decision turns on, and the trial's
  own population-PK lambda becomes a CROSS-CHECK rather than a prerequisite. The honest
  cost: freeing lambda does degrade theta (e.g. monthly 2.3 -> ~3.5), because theta and
  lambda share some of the dose-response direction -- but theta was not individualising
  the decision anyway (item 1), whereas lambda IS recovered to spec. This is the asymmetry
  that matters: the record gives up little that was usable (theta) to pin the thing that
  was binding (lambda).
""")

    # ---- ITEM 3 ----
    print("=" * 92)
    print("  ITEM 3. LENGTHENING IS THE RISKY DIRECTION -- A ONE-SIDED TEST, NOT A CLASSIFICATION")
    print("=" * 92)
    print(f"""
  Shortening is safe and only costs money; lengthening (14 -> 15 or 21 d) risks a crisis.
  So lengthen a patient only when the record is confident the LONGER trough still clears
  theta_crit: a one-sided test -- lengthen iff the lower 95% bound on
  mu_long = theta_base + Delta_trough(tau_long) - theta_crit exceeds 0.
  Delta_trough: 14d={trough(14):.0f}, 15d={trough(15):.0f}, 21d={trough(21):.0f} uM/d;
  Model C theta_crit^int = {THC_INT:.0f}. At 15 d the trough alone clears theta_crit for
  everyone (no evidence needed); at 21 d it does NOT ({trough(21):.0f} < {THC_INT:.0f}), so
  lengthening to 21 d is the real evidence-gated decision.
""")
    print(f"    LENGTHEN-ELIGIBLE fraction (lower 95% bound on mu_long > 0), Model C theta_crit:")
    print(f"        {'record / precision':>30}  {'sd(ln th)':>9}  {'to 15 d':>8}  {'to 21 d':>8}")
    rows = [("pre-therapy flat", flat), ("regular monthly (the 3.05)", reg)]
    rows += [("event, " + lab, dens_sd[lab]) for lab, _ in DENSITIES]
    rows += [("managed-load protocol", 0.29)]
    for lab, sd in rows:
        f15 = lengthen_fraction(sd, 15, THC_INT)
        f21 = lengthen_fraction(sd, 21, THC_INT)
        print(f"        {lab:>30}  {sd:>9.2f}  {f15:>8.0%}  {f21:>8.0%}")
    print("""
  READING: lengthening to 15 d is safe for the whole cohort under Model C regardless of
  precision (the trough dominates), so it needs no evidence -- a free burden reduction.
  Lengthening to 21 d is the genuine one-sided decision, and here the two-year record
  DOES deliver an actionable, safe subgroup: the eligible fraction climbs with sampling
  density (0% monthly -> ~50% at alternate-day), each patient lengthened only on a
  confident lower bound, never a point estimate. So the operative conclusion is not the
  symmetric classification (which the record cannot settle) but a one-sided, safe-by-
  construction lever: shorten freely, lengthen to 15 d for all, and lengthen to 21 d for
  the growing subgroup whose own accumulated record proves it, with lambda self-calibrated
  from the same data.
""")


def determined(sd, thc=None, grid=(7, 10, 14), th_lo=5.0, th_hi=150.0,
               npat=4000, ndraw=300, seed=0):
    """Correct-choice and determined-subgroup over a PA population (routine load)."""
    thc = theta_crit_A(1.0) if thc is None else thc
    rng = np.random.default_rng(seed)
    th = np.exp(rng.uniform(np.log(th_lo), np.log(th_hi), npat))

    def maxsafe(tb):
        s = [t for t in grid if tb + trough(t) >= thc]
        return max(s) if s else None
    true = [maxsafe(t) for t in th]
    hit = tot = 0
    for t, tr in zip(th, true):
        d = t * np.exp(rng.normal(0.0, sd, ndraw))
        hit += sum(maxsafe(x) == tr for x in d); tot += ndraw
    bnds = [thc - trough(t) for t in grid if thc - trough(t) > 0]
    if not bnds:
        det = 1.0
    else:
        logd = np.array([min(abs(np.log(t) - np.log(b)) for b in bnds) for t in th])
        det = float(np.mean(logd > Z1 * sd))
    return hit / tot, det


def lengthen_fraction(sd, tau_long, thc, th_lo=5.0, th_hi=150.0, npat=8000, seed=1):
    """Fraction of the PA population for whom the lower 95% bound on
    theta_base + Delta_trough(tau_long) - theta_crit exceeds 0 (eligible to lengthen).
    The uncertainty is on theta_base (lognormal, sd); Delta_trough is fixed by the schedule."""
    rng = np.random.default_rng(seed)
    th = np.exp(rng.uniform(np.log(th_lo), np.log(th_hi), npat))
    lo = th * np.exp(-Z1 * sd)                      # lower 95% bound on theta_base
    return float(np.mean(lo + trough(tau_long) > thc))


if __name__ == "__main__":
    main()
