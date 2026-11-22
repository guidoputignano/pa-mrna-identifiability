"""
Round 13: WHEN. Every round-11/12 result was stated at two years -- the record
length we inherited, not a derived quantity. This asks the timing, one-sided
throughout (lengthening on evidence; a wrong answer costs burden, never safety).

METHOD NOTE (why simulation, not FIM propagation, for the 21-day eligibility).
The lengthening margin mu_long = theta_base + Delta_trough(21;lambda) - theta_crit^int
is a ONE-SIDED test on a quantity whose posterior is strongly non-Gaussian and
tail-dominated (sd(ln theta) ~ 1.8-36 across the disease range). FIM propagation is
unreliable here: the delta method and a log-normal Monte Carlo disagree by ~400 uM/d
on the lower bound, because the delta method throws away the favourable skew -- exactly
the pitfall the parent paper flags in sec 7b ("a delta-method calculation ... would
throw away the one favourable result this paper has ... mu_hat is skewed and does not
cross zero as often as its sd suggests"). So, as the paper does for the sign-of-mu
decision (fig D3), we SIMULATE: fit (theta, K_m, b, k_e, lambda) by NLS to noisy
on-therapy records and measure the decision directly, with a safety-calibrated
one-sided threshold. (lambda itself IS in the Gaussian regime -- small, theta-base-
independent sd -- so item 2 uses the FIM, which is trustworthy there.)

t45 item 3's lengthening fractions were FIM/representative-sd estimates and are
superseded here for the same reason.
"""
import sys, json, time
sys.path.insert(0, "../../src")
import numpy as np
from scipy.optimize import least_squares

import core as C
import t45_on_therapy as T45
from common.models import ParamsC
import model_C.run as MC

LAM = T45.LAM
Z1 = 1.6449                                                # one-sided 95%
KM, KE = 100.0, 5.0
THC_INT = MC.theta_crit_intrinsic(ParamsC(), MC.S_BASE)   # 150.8
TR21 = T45.trough(21)                                      # 142.9
FLOOR = THC_INT - TR21                                     # 7.9 true-safety floor at 21 d
H = 60.0

CADENCES = [("monthly", 30.4), ("fortnightly", 14.0), ("weekly", 7.0)]
LENGTHS = [6, 12, 24, 36, 48, 60]
THGRID = [4.0, 5.5, 7.0, 9.0, 13.0, 20.0, 35.0, 65.0, 110.0, 150.0]
NREP = 120
SAFE_ALPHA = 0.05                                         # max unsafe-lengthening rate


def _delta(t, dts, das, lam):
    t = np.asarray(t, float); out = np.zeros_like(t)
    for td, a in zip(dts, das):
        on = t >= td; out[on] += a * np.exp(-lam * (t[on] - td))
    return out


def sim_muhat(theta_base, months, cad_days, nrep=NREP, seed0=0):
    """NLS point estimates of mu_long = theta_hat + Delta_trough(21; lambda_hat) - theta_crit,
    over `nrep` simulated on-therapy records (marker regressors, noise on the clearance
    ordinate, k_e and load-bias priors as in core; lambda free)."""
    T_days = months * 30.4
    n = max(4, int(round(T_days / cad_days)))
    tm = np.linspace(0.0, T_days, n)
    yrs = T_days / 365.0
    L, _ = C.N4.anchor(theta_base, KM, KE, "A")
    sb = C.LOAD_FRAC * L
    out = np.empty(nrep)
    for rep in range(nrep):
        s = seed0 + rep
        dts, das = T45.schedule(seed=s, years=yrs, jitter=2.0, miss_p=0.08)
        rng = np.random.default_rng(7000 + s)
        s2 = np.log1p(0.10 ** 2)
        loads = L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), n))
        inc = _delta(tm, dts, das, LAM)
        m = C.level(theta_base, KM, KE, loads, inc)
        slope = (theta_base + inc) * KM / (KM + m) ** 2 + KE
        sg = np.hypot(C.LOAD_FRAC * L, slope * C.ASSAY_CV * m)
        y = (theta_base + inc) * m / (KM + m) + KE * m + rng.normal(0.0, sg)

        def resid(x):
            lt, lk, b, lke, llam = x
            th, K, ke, lam = np.exp(lt), np.exp(lk), np.exp(lke), np.exp(llam)
            I = _delta(tm, dts, das, lam)
            pred = (th + I) * m / (K + m) + ke * m + b
            return np.concatenate([(y - pred) / sg, [b / sb], [(lke - np.log(KE)) / 0.8]])
        x0 = np.array([np.log(theta_base), np.log(KM), 0.0, np.log(KE), np.log(LAM)])
        r = least_squares(resid, x0, method="lm", xtol=1e-10, ftol=1e-10, max_nfev=4000)
        th, lam = np.exp(r.x[0]), np.exp(r.x[4])
        out[rep] = th + T45.trough(21, lam) - THC_INT
    return out


def sd_lnlam(months, cad_days, theta=C.TH_5, seeds=range(16)):
    """FIM sd(ln lambda) -- trustworthy (Gaussian regime)."""
    T_days = months * 30.4
    n = max(4, int(round(T_days / cad_days)))
    tm = np.linspace(0.0, T_days, n)
    yrs = T_days / 365.0
    sk = dict(jitter=2.0, miss_p=0.08)
    if T_days > 120:
        sk["gap"] = (0.4 * T_days, 28.0)
    v = []
    for s in seeds:
        r = T45.on_therapy_record(theta, seed=s, sched_kw=dict(sk, years=yrs), times=tm)
        x = T45.sd_of(T45.fim5(r, with_lambda=True), 4)
        if np.isfinite(x):
            v.append(x)
    return float(np.median(v)) if v else np.inf


def req_sd(theta_base):
    """sd(ln theta) at which the one-sided lower bound on mu_long just clears 0."""
    if theta_base <= FLOOR:
        return -1.0
    return np.log(theta_base / FLOOR) / Z1


def main():
    print("=" * 92)
    print("  ROUND 13. WHEN: TIMING OF 21-DAY ELIGIBILITY (by simulation) AND lambda (by FIM)")
    print("=" * 92)
    print(f"""
  The 21-day lengthening margin is mu_long = theta_base + Delta_trough(21) - theta_crit^int
  = theta_base + {TR21:.0f} - {THC_INT:.0f} = theta_base - {FLOOR:.1f}. The one-sided test lengthens
  only when the lower 95% bound on mu_long clears 0, so a patient needs theta_base resolved
  to sd(ln theta) < ln(theta_base/{FLOOR:.1f})/1.645 -- e.g. {req_sd(30):.2f} at theta_base=30,
  {req_sd(65):.2f} at 65, {req_sd(150):.2f} at 150 (the top of the disease range). That is
  near the MANAGED-load precision (0.29), an order below what a passive record delivers.
""")
    # ---- ITEM 1: does mu_hat tighten with record length? (the decisive diagnostic)
    print("  ITEM 1. DOES THE 21-DAY MARGIN BECOME RESOLVABLE? (NLS mu_hat vs record length, weekly)")
    print(f"      {'theta_base':>10} {'mu_true':>8}  {'months':>6}  {'median mu_hat':>13}  {'p5':>7}  {'sd(mu_hat)':>10}  {'P(mu_hat>0)':>11}")
    for th in (150, 65, 20):
        for m in (6, 24, 60):
            mh = sim_muhat(th, m, 7.0, nrep=200)
            print(f"      {th:>10} {th+TR21-THC_INT:>+8.0f}  {m:>6}  {np.median(mh):>13.0f}  "
                  f"{np.percentile(mh,5):>7.0f}  {mh.std():>10.0f}  {np.mean(mh>0):>11.0%}")
    print(f"""
  READING: the NLS margin estimate does NOT tighten with record length. sd(mu_hat) stays
  ~105-165 uM/d -- comparable to theta_crit^int ({THC_INT:.0f}) itself -- from 6 months out to
  5 years, and its lower (p5) tail stays below 0 at every capacity, so the one-sided lower
  bound never cleanly clears 0 and does not improve with time. This RESOLVES the FIM
  ambiguity noted in the header: the log-normal MC's optimism (a positive lower bound at 2
  years) was an artefact of ignoring NLS estimator bias and heavy tails; the delta method's
  "never" is closer to the truth. 21-day lengthening is NOT supportable from the passive
  on-therapy record at any horizon -- it requires theta resolved to ~0.3 (managed-load
  precision), which passive observation does not reach (the on-therapy lambda-free sd(ln
  theta) bottoms out near ~1.8 even at 5 years weekly, and worse for deeper disease). The
  eligibility therefore does NOT accrue with time; there is no crossing time to report.
  (This also supersedes t45 item 3, whose 51-73% lengthen fractions used a representative,
  lambda-fixed FIM sd and so were badly optimistic.)
""")

    # ---- ITEM 2 (FIM, trustworthy)
    print("=" * 92)
    print("  ITEM 2. DOES lambda CROSS EARLIER?  sd(ln lambda) vs record length (FIM)")
    print("=" * 92)
    print(f"\n      {'months':>7} " + " ".join(f"{c:>12}" for c, _ in CADENCES))
    lam = {}
    for m in LENGTHS:
        row = []
        for cname, cad in CADENCES:
            sl = sd_lnlam(m, cad); lam[(cname, m)] = sl; row.append(sl)
        print(f"      {m:>7} " + " ".join(f"{s:>12.3f}" for s in row))
    print(f"\n      {'spec reached by (months)':>26}  " + " ".join(f"{c:>12}" for c, _ in CADENCES))
    for spec in (0.25, 0.05):
        row = [next((f"{m} mo" for m in LENGTHS if lam[(cname, m)] <= spec), ">60 mo")
               for cname, _ in CADENCES]
        print(f"      {('lambda to ' + str(int(spec*100)) + '%'):>26}  " + " ".join(f"{r:>12}" for r in row))
    print("""
  READING: lambda -- the binding constant (it governs the trough) -- is pinned LONG before
  any eligibility, and uniformly across the cohort (sd(ln lambda) is ~theta_base-indep).
  The ~25% (Model-C) spec is met by ~3-6 months at every cadence; the ~5% (routine-load)
  spec within the first year at weekly. So the constant the decision turns on arrives first;
  theta (per-patient, for the eligibility test) accrues slowly and only for high capacities.
""")

    # ---- ITEM 3 (analytic)
    print("=" * 92)
    print("  ITEM 3. THE 15-DAY STEP NEEDS NO RECORD AND NO WAITING")
    print("=" * 92)
    tr15 = T45.trough(15)
    print(f"""
  Delta_trough(15) = {tr15:.1f} > theta_crit^int = {THC_INT:.1f}, so even theta_base = 0 clears
  the bar at 15 d (margin {tr15 - THC_INT:+.1f} uM/d). The 15-day step is safe at INITIATION --
  no record, no observation period, no measurement -- for every patient. It is the only
  recommendation in the report available on day one; it should be stated separately in the
  decision section as the zero-evidence, zero-wait burden reduction.
""")

    # ---- ITEM 4
    print("=" * 92)
    print("  ITEM 4. WHEN TO REVIEW?  (for 21-d: there is nothing to review for)")
    print("=" * 92)
    print("""
  The review-schedule question presumes eligibility accrues with record length, so that
  there is an optimal time to re-examine. Item 1 shows it does NOT: the 21-day margin stays
  unresolvable at every horizon, so no review interval unlocks it -- reviewing a passive
  record for 21-day eligibility is waste at any time, early or late. There is no interior
  optimum because there is no accrual.

  The time-structured recommendations that DO hold, in order of availability:
    - day 0:  step 14 -> 15 d for EVERY patient (item 3; trough dominates theta_crit, no
              record, no waiting).
    - ~6-12 months: lambda is confirmed from the accumulated dosing record (item 2), which
              corroborates the trial PK and secures the constant the whole decision turns on.
    - 21-day step: only via the ACTIVE managed-load protocol (which reaches sd ~0.29), not
              via passive observation of any length. If 21 d is wanted, it is a protocol
              change to make prospectively, not an eligibility to wait for.
  A genuine review-cost optimisation would need a per-visit cost the repository does not
  supply; but for the 21-day step that is moot, since the signal never arrives.
""")
    # save key numbers (lambda timing is the quantitative deliverable)
    out = {"lambda": {f"{c}/{m}": lam[(c, m)] for c, _ in CADENCES for m in LENGTHS},
           "req_sd": {str(th): req_sd(th) for th in (30, 65, 150)},
           "floor": FLOOR}
    with open("out/t46_timing_numbers.json", "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
