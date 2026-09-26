"""
The pricing: what each measurement buys, and what the repeatability study costs.

Three of the four quantities the specification would have to ask for move the
SAME sigma budget, and step 3(a) proved that the design sees them only through
sigma once the load-bias prior is named. So they are not three sweeps, they are
one arithmetic constraint and this script states it as one. The fourth, the
marker's unreachable floor, acts on identifiability instead and is priced apart.

    priced jointly, on the sigma budget
      k_e precision            sets WHICH budget applies, through sigma*
      marker residual scatter  spends the budget
      load-record accuracy     spends the budget, and also moves s_b

    priced separately, on identifiability
      the floor anchor         step 3(b) found it does not buy a usable theta at
                               any marker CV the literature supports, so what is
                               priced here is not the precision it delivers but
                               the marker CV at which it would begin to be worth
                               having at all

AND THE REPEATABILITY STUDY, which is the only study that measures the quantity
sigma actually needs. PHASE9 item 7 found the residual scatter undecomposable
from the published surrogates and left sigma as a bound, 5 to 8% below and 25 to
35% above. Repeat marker draws in treated patients at ONE FIXED dose-cycle phase
with the diet recorded measure the residual directly, because holding the phase
removes the dose cycle and recording the diet removes the load term. This script
asks how many patients and how many repeats that takes.

WHAT THE k_e MEASUREMENT BUYS IS TWO THINGS, not one. It narrows the spread, and
it removes the collapse tail that PHASE9 3a.4 found at the population prior,
about 2% of records failing by orders of magnitude with no warning in the fit.
Both are priced.

Nothing here records an expected result in advance.

Run:  python audit/q8_pricing.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.stats import chi2, beta

import q2_requirement_map as S2
import q6_sigma_budget as S6

SCHED = S6.SCHED
KE_CASES = ((0.8, "population prior, no study"),
            (0.054, "study at 2% with 5% drift"),
            (0.02, "study at 2% with no drift"))
LOADS = (0.02, 0.05, 0.10, 0.15)
SEEDS = S2.SEEDS


def banner(t):
    print("\n" + "=" * 104)
    print("  " + t)
    print("=" * 104)


def cv_allowed(sigma_star, load_frac, theta=3000.0):
    """
    The marker CV a budget allows once the load record has taken its share.

    sigma = hypot(load_frac*L, b*cv*m_bar), so this inverts that exactly. None
    where the load record alone already exceeds the budget.
    """
    L, m_bar, b = S6.anchor(theta)
    rem = sigma_star ** 2 - (load_frac * L) ** 2
    return None if rem <= 0 else float(np.sqrt(rem) / (b * m_bar))


# ==========================================================================
def part_a():
    banner("A. THE THREE sigma AXES, PRICED JOINTLY")
    print("""  One table, not three sweeps. `k_e` picks the budget; the load record and the
  marker CV spend it. Each cell is the largest marker residual scatter admissible
  at that load record, on that schedule, with that `k_e`. `none` means the load
  record alone exceeds the budget, so no assay reaches the gate.""")
    theta = 3000.0
    L = S6.anchor(theta)[0]
    for ke, label in KE_CASES:
        print(f"\n  k_e: {label} ({ke})")
        print(f"    {'schedule':16s} {'sigma*':>8s}" +
              "".join(f"{f'load {int(100 * x)}%':>13s}" for x in LOADS))
        for name, times in SCHED.items():
            st = S6.sigma_star(theta, times, 0.10 * L, ke_prior=ke)
            cells = []
            for lf in LOADS:
                v = cv_allowed(st, lf)
                cells.append(f"{'none':>13s}" if v is None
                             else f"{100 * v:12.2f}%")
            print(f"    {name:16s} {st:8.2f}" + "".join(cells))
    print("""
  Read a row to see the trade: the two axes are interchangeable at equal
  percentages, which is S3's crossover, and the whole row moves when `k_e` moves,
  which is 3a.6's two tables. That is the entire joint pricing, and it is
  arithmetic rather than a simulation, because step 3(a) proved the budget exact.""")


def part_b():
    banner("B. WHAT THE k_e MEASUREMENT BUYS, BOTH THINGS")
    print("""  The study buys a wider budget AND the removal of a failure mode. The second is
  not a refinement of the first: at the population prior about 2% of records
  collapse by orders of magnitude with no warning in the fit, and PHASE9 3a.4
  named the cause, a prior centre drawn several times above the truth. Both are
  reported, with the collapse rate carrying its binomial interval.""")
    theta = 3000.0
    L = S6.anchor(theta)[0]
    s_b = 0.10 * L
    print(f"\n    {'schedule':16s} {'sigma* at 0.8':>13s} {'sigma* at 0.054':>16s} "
          f"{'budget gain':>12s} {'CV gain at load 10%':>20s}")
    for name, times in SCHED.items():
        a = S6.sigma_star(theta, times, s_b, ke_prior=0.8)
        c = S6.sigma_star(theta, times, s_b, ke_prior=0.054)
        va, vc = cv_allowed(a, 0.10), cv_allowed(c, 0.10)
        g = ("n/a" if va is None or vc is None
             else f"{100 * va:.1f}% to {100 * vc:.1f}%")
        print(f"    {name:16s} {a:13.2f} {c:16.2f} {c / a:11.2f}x {g:>20s}")
    print("""
  And the tail, at sigma* on the centre convention, which is what a per-patient
  measurement supplies. 96 records per cell; the interval is Clopper-Pearson.""")
    print(f"\n    {'schedule':16s} {'k_e prior':>10s} {'robust':>8s} "
          f"{'collapses':>10s} {'rate':>8s} {'95% interval':>18s}")
    for name, times in SCHED.items():
        st = S6.sigma_star(theta, times, s_b, ke_prior=0.8)
        for ke in (0.8, 0.054):
            r = S6.mc_at_sigma(theta, times, st, s_b, ke_sd=ke,
                               nuisance="centre")
            k, n = r["collapse"], 96
            lo = 0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1))
            hi = 1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k))
            print(f"    {name:16s} {ke:10.3f} {r['robust']:8.3f} "
                  f"{k:10d} {100 * k / n:7.2f}% "
                  f"{f'[{100 * lo:.2f}%, {100 * hi:.2f}%]':>18s}")


# ==========================================================================
#  The repeatability study.
# ==========================================================================
def sd_ln_sigma(df):
    """
    The relative precision of an estimated standard deviation on `df` degrees of
    freedom, exactly rather than asymptotically.

    s^2/sigma^2 ~ chi2(df)/df, so ln s has the variance of half a log-chi2. The
    asymptotic 1/sqrt(2 df) is printed beside it so the approximation's error is
    visible rather than assumed.
    """
    #  Var(ln s^2) = trigamma(df/2), so sd(ln s) = 0.5 sqrt(trigamma(df/2)).
    return float(0.5 * np.sqrt(_polygamma1(df / 2.0)))


def _polygamma1(x):
    """Trigamma by the recurrence plus an asymptotic tail, so no scipy.special
    import is needed for a one-line quantity."""
    s = 0.0
    while x < 20.0:
        s += 1.0 / (x * x)
        x += 1.0
    inv = 1.0 / x
    return s + inv * (1.0 + inv * (0.5 + inv * (1.0 / 6.0
                                                - inv * inv / 30.0)))


def ci_halfwidth(df, conf=0.95):
    """The 95% interval for sigma as a multiplicative factor, from chi2."""
    a = (1.0 - conf) / 2.0
    lo = np.sqrt(df / chi2.ppf(1.0 - a, df))
    hi = np.sqrt(df / chi2.ppf(a, df))
    return float(lo), float(hi)


def part_c():
    banner("C. THE REPEATABILITY STUDY, SIZED")
    print("""  The design: `n` treated patients, `r` repeat marker draws each, EVERY DRAW AT
  THE SAME DOSE-CYCLE PHASE, with the diet recorded. Holding the phase removes the
  dose cycle, which PHASE9 item 7 measured at 27 to 30% and which is a regressor
  rather than error. Recording the diet removes the load term. What is left in the
  within-patient scatter is the residual plus the analytical component, which is
  what sigma needs and what no published study reports.

  The pooled within-patient variance carries df = n(r - 1). Its precision is
  exact from the chi-square, and the asymptotic 1/sqrt(2 df) is shown beside it so
  the approximation's error is visible.""")
    print(f"\n    {'n':>5s} {'r':>3s} {'df':>5s} {'sd(ln sigma)':>13s} "
          f"{'asymptotic':>11s} {'95% factor on sigma':>21s}")
    for n, r in ((5, 2), (10, 2), (10, 3), (20, 3), (10, 6), (25, 3), (20, 6),
                 (50, 3), (30, 6)):
        df = n * (r - 1)
        lo, hi = ci_halfwidth(df)
        print(f"    {n:5d} {r:3d} {df:5d} {sd_ln_sigma(df):13.4f} "
              f"{1.0 / np.sqrt(2 * df):11.4f} "
              f"{f'x{lo:.2f} to x{hi:.2f}':>21s}")
    print("""
  What precision is ENOUGH is set by the decision it has to support. The schedule
  boundaries are the marker CVs at which one schedule gives way to the next, from
  part A at a 10% load record, and a study resolves a boundary when its interval
  for the true CV excludes it.""")
    theta = 3000.0
    L = S6.anchor(theta)[0]
    bounds = []
    for name, times in SCHED.items():
        st = S6.sigma_star(theta, times, 0.10 * L, ke_prior=0.8)
        v = cv_allowed(st, 0.10)
        if v is not None:
            bounds.append((name, 100 * v))
    print(f"\n    boundaries at a 10% load record, population k_e: " +
          ", ".join(f"{nm} {b:.2f}%" for nm, b in bounds))
    print(f"\n    {'true CV':>8s} {'nearest boundary':>18s} {'gap':>8s} "
          f"{'df needed':>10s} {'n at r=3':>9s} {'n at r=6':>9s}")
    for true_cv in (6.0, 8.0, 12.0, 15.0, 25.0, 30.0):
        near = min(bounds, key=lambda x: abs(x[1] - true_cv))
        gap = abs(near[1] - true_cv) / true_cv
        if gap <= 0:
            continue
        #  the smallest df whose 95% factor interval excludes the boundary
        need = None
        for df in range(2, 20001):
            lo, hi = ci_halfwidth(df)
            if (near[1] > true_cv and hi < near[1] / true_cv) or \
               (near[1] < true_cv and lo > near[1] / true_cv):
                need = df
                break
        if need is None:
            print(f"    {true_cv:7.1f}% {near[0]:>18s} {100 * gap:7.1f}% "
                  f"{'>20000':>10s}")
            continue
        print(f"    {true_cv:7.1f}% {near[0]:>18s} {100 * gap:7.1f}% "
              f"{need:10d} {int(np.ceil(need / 2)):9d} "
              f"{int(np.ceil(need / 5)):9d}")
    print("""
  `n at r = 3` is two repeat draws beyond the first per patient, `n at r = 6` is
  five. Read the column that matches what a clinic would agree to.""")


def part_c2(reps=4000):
    banner("C2. THE SIZING CONFIRMED BY SIMULATION")
    print("""  The chi-square precision assumes the within-patient scatter is normal on the log
  scale and homogeneous across patients. Simulated here instead: `n` patients, `r`
  draws each, log-normal residual with the stated CV, pooled estimate, and the
  empirical coverage of the nominal 95% interval.""")
    rng = np.random.default_rng(8_2026)
    print(f"\n    {'n':>4s} {'r':>3s} {'df':>5s} {'true CV':>8s} "
          f"{'mean CV_hat':>12s} {'sd(ln) exact':>13s} {'sd(ln) MC':>11s} "
          f"{'coverage':>9s}")
    for n, r, cv in ((10, 3, 0.10), (25, 3, 0.10), (10, 6, 0.10),
                     (25, 3, 0.25), (10, 6, 0.25)):
        df = n * (r - 1)
        est, cov = [], 0
        s2 = np.log1p(cv ** 2)
        for _ in range(reps):
            #  each patient's r draws share a level; the scatter is log-normal
            y = rng.normal(0.0, np.sqrt(s2), (n, r))
            within = np.sum((y - y.mean(axis=1, keepdims=True)) ** 2) / df
            est.append(np.sqrt(np.expm1(within)))
            lo, hi = ci_halfwidth(df)
            s_hat = np.sqrt(within)
            if s_hat * lo <= np.sqrt(s2) <= s_hat * hi:
                cov += 1
        e = np.array(est)
        print(f"    {n:4d} {r:3d} {df:5d} {100 * cv:7.1f}% "
              f"{100 * e.mean():11.2f}% {sd_ln_sigma(df):13.4f} "
              f"{np.std(np.log(e)):11.4f} {100 * cov / reps:8.1f}%")


def part_d():
    banner("D. THE FLOOR ANCHOR, PRICED ON ITS OWN TERMS")
    print("""  Step 3(b) removed the question this was going to answer. Fixing the marker's
  unreachable floor restores STRUCTURAL identifiability, rank 3 of 3, and does not
  deliver a usable theta: at phi = 10% and the committed 5% marker CV the fitted
  robust spread is 6.1 in ln theta against a Fisher bound of 0.135, because the
  profile is a needle on a flat floor. So what can be priced is not the precision
  the anchor delivers but the marker CV at which it would begin to be worth
  having, which is a different and harsher question.""")
    print("""
  From 3b.2, fitted, offset held at its own correct value, load bias drawn and
  fitted, 80 records:

      phi = 10%:  CV 5%   robust 6.12      CV 0.5%  robust 3.15
                  CV 1.5% robust 5.40      CV 0.1%  robust 0.121  <- gate met
      phi =  5%:  CV 5%   robust 6.34      CV 0.05% robust 3.44    never met

  So the anchor becomes worth having at a marker CV of about 0.1%, and only at the
  upper admissible anchor. Priced against the other axis:""")
    print(f"\n    {'requirement':44s} {'value':>12s} {'source':>28s}")
    rows = [("marker CV at which the floor anchor pays", "about 0.1%",
             "PHASE9 3b.2, fitted"),
            ("analytical floor on the marker CV", "5 to 8%",
             "PHASE9 item 7, Thompson 2012"),
            ("ratio between them", "50 to 80x", "arithmetic"),
            ("marker CV the surrogates suggest", "25 to 35%",
             "PHASE9 item 7"),
            ("ratio between them", "250 to 350x", "arithmetic")]
    for a, b, c in rows:
        print(f"    {a:44s} {b:>12s} {c:>28s}")
    print("""
  THE PRICE IS THEREFORE NOT PAYABLE. The floor anchor needs an assay two orders of
  magnitude better than the best the literature supports, so verifying R5's
  transplant ratios cannot on its own make theta recoverable under the published
  algebraic structure. What the anchor is still worth is the structural statement,
  which is what the paper's identifiability claim rests on, and that is not
  nothing: it is the difference between a claim that is false and one that is true
  but not yet operational.""")


def part_e(reps=3000):
    banner("E. THE BETWEEN-PATIENT SPREAD OF THE RESIDUAL SCATTER")
    print("""  Part C's design pools the within-patient variance and estimates ONE residual
  scatter. If the scatter differs between patients, no single requirement can be
  written for every patient and the requirement framing needs a per-patient form.
  A pooled design cannot recover that afterwards, so it is a decision to take
  before the study runs rather than a question to file.

  The model: patient i has ln sigma_i ~ N(ln sigma, tau), and r draws give that
  patient's own estimate with sampling sd 0.5 sqrt(trigamma((r-1)/2)). The
  observed spread of ln sigma_hat across patients is tau and that sampling noise
  together, so tau is recovered by subtracting the second from the first.

  THE SAMPLING NOISE PER PATIENT IS SET BY r ALONE, and it is large at small r:""")
    print(f"\n    {'r':>4s} {'df per patient':>15s} "
          f"{'sampling sd(ln sigma_i)':>24s}")
    for r in (2, 3, 4, 6, 8, 12):
        print(f"    {r:4d} {r - 1:15d} {sd_ln_sigma(r - 1):24.3f}")
    print("""
  So more patients cannot buy what more draws per patient buy: n averages the
  noise down in the MEAN but leaves it inside the VARIANCE that tau has to be
  separated from. Simulated below: `detect` is the fraction of studies whose
  method-of-moments tau-hat exceeds zero by more than its own 95% interval, which
  is the question "can this study show the scatter varies between patients at
  all".""")
    rng = np.random.default_rng(6_2026)

    def one(n, r, tau, cv=0.15):
        s2 = np.log1p(cv ** 2)
        lnsig = np.log(np.sqrt(s2)) + rng.normal(0.0, tau, n)
        within = np.empty(n)
        for i in range(n):
            y = rng.normal(0.0, np.exp(lnsig[i]), r)
            within[i] = np.sum((y - y.mean()) ** 2) / (r - 1)
        lh = 0.5 * np.log(within)
        obs = float(np.var(lh, ddof=1))
        samp = sd_ln_sigma(r - 1) ** 2
        return max(obs - samp, 0.0), obs, samp

    print(f"\n    {'n':>4s} {'r':>3s} {'draws':>6s} {'true tau':>9s} "
          f"{'median tau_hat':>15s} {'detect':>8s} {'pooled df':>10s}")
    for n, r in ((10, 3), (25, 3), (100, 3), (10, 6), (25, 6), (16, 6),
                 (10, 12), (25, 12)):
        for tau in (0.0, 0.30):
            ests, det = [], 0
            samp = sd_ln_sigma(r - 1) ** 2
            for _ in range(reps):
                t2, obs, _ = one(n, r, tau)
                ests.append(np.sqrt(t2))
                #  the observed variance exceeds the null by more than the
                #  chi-square upper tail of the null spread
                crit = samp * chi2.ppf(0.95, n - 1) / (n - 1)
                det += obs > crit
            print(f"    {n:4d} {r:3d} {n * r:6d} {tau:9.2f} "
                  f"{np.median(ests):15.3f} {100 * det / reps:7.1f}% "
                  f"{n * (r - 1):10d}")
    print("""
  Read the `detect` column at true tau = 0 as the test's false-positive rate,
  which should sit near 5%, and at tau = 0.30 as its power.""")


def main():
    print("=" * 104)
    print("  THE PRICING")
    print("=" * 104)
    part_a()
    part_b()
    part_c()
    part_c2()
    part_d()
    part_e()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
