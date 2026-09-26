"""
The hierarchical cohort fit, priced.

`k_e` and `K_m` are POPULATION parameters in this model and `theta` is
per-patient. A cohort fit therefore estimates the shared ones FROM THE COHORT
while each `theta_i` stays individual, which obtains them without measuring them
in anyone.

PHASE 10'S REDUCTION MADE THIS URGENT. The paper's second-channel constant is
`k_e/alpha`, not the marker's clearance (PHASE10 R.2), so the tracer route now
needs one IV micro-dose study PER PANEL MARKER plus a closure assumption on the
overflow routes (PHASE9 P.8). A cohort fit estimates `k_e/alpha` directly from
marker records: it is the coefficient multiplying `m` in the balance the records
obey, whatever that coefficient is made of, so the fit never needs the conversion
and never meets `alpha` at all. That is a large advantage arriving from a
direction neither strand was looking, and it is measured here before anything
else.

THE QUESTION THAT DECIDES IT IS NOT THE SAMPLE SIZE. Brief R1 states that the
clearance is patient-specific: 3-HP is renally secreted by OAT1, half of adults
with PA have an eGFR under 60, and the between-patient CV is "probably a third or
more". A cohort fit with a SHARED constant is then misspecified, and a
hierarchical one with a random `k_e_i/alpha_i` pins the POPULATION value while
each patient's own stays uncertain by the between-patient spread. So the sweep
that decides the route is over that spread, and the headline is not "how many
patients" but "how many patients, and does it help at all when the spread is a
third".

WHAT IS COMPUTED
  A  the population constant against cohort size, each patient's own
     (ln theta_i, b_i) profiled out by the exact Schur complement.
  B  what a PATIENT's own constant is then known to, swept over the
     between-patient CV. This is the number their theta_i actually consumes.
  C  the headline cells confirmed by FITTING a simulated cohort jointly, per the
     standing house rule.
  D  fed back into the requirement table beside the micro-dose route.

Throughout, `ke` names the model's second-channel constant `k_e/alpha`, which is
the only thing any marker record can see. Nothing here records an expected result
in advance.

Run:  python audit/q14_hierarchical.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import least_squares

import q2_requirement_map as S2

KM, KE, TH_NOM = S2.KM, S2.KE, S2.TH_NOM
#: the trial's own cohort size is 16; the rest is what a larger study might reach
NS = (4, 8, 16, 24, 32, 48, 64, 96, 128)
#: R1's between-patient CV of the disposal constant, swept. 0 = genuinely shared.
CVS = (0.0, 0.10, 0.20, 0.33, 0.50)
#: the capacities a real cohort would span, as fractions of nominal
COHORT_CAPS = (0.05, 0.08, 0.10, 0.15, 0.20)
#: S2's two targets for the disease range
TARGETS = (0.02, 0.06)
#: the paper's own population prior on ln k_e, used as the population-level prior
POP_PRIOR = 0.8
SCHEDS = ("monthly", "weekly 7.3 d", "every 2 d")


def banner(t):
    print("\n" + "=" * 104)
    print("  " + t)
    print("=" * 104)


def blocks(theta, times, seed):
    """
    One patient's information, split shared / individual.

    Order: shared (ln K_m, ln k_e), individual (ln theta, b). No prior is
    included, so the caller adds each at the level it belongs to.
    """
    L, _, sigma, s_b = S2.footing(theta)
    Ls = S2.loads_for(times.size, seed, L)
    inc = S2.increment(times)
    m = S2.level(theta, inc, Ls)
    g = m / (KM + m)
    Sh = np.stack([-(theta + inc) * KM * m / (KM + m) ** 2, KE * m],
                  axis=1) / sigma
    In = np.stack([theta * g, np.ones_like(m)], axis=1) / sigma
    return Sh.T @ Sh, Sh.T @ In, In.T @ In, s_b


def shared_info(theta, times, seed):
    """What one patient contributes to (ln K_m, ln k_e) once their own
    (ln theta, b) are profiled out. The exact Schur complement."""
    F_ss, F_si, F_ii, s_b = blocks(theta, times, seed)
    F_ii = F_ii + np.diag([0.0, 1.0 / s_b ** 2])
    return F_ss - F_si @ np.linalg.solve(F_ii, F_si.T)


def cohort_info(n, sched, seed0=0):
    """(total shared information, the first patient's own share)."""
    times = S2.times_for(sched, "no baseline")
    tot, first = np.zeros((2, 2)), None
    for i in range(n):
        th = TH_NOM * COHORT_CAPS[i % len(COHORT_CAPS)]
        Fi = shared_info(th, times, seed0 + i)
        tot = tot + Fi
        if first is None:
            first = Fi
    return tot, first


def sd_population(n, sched):
    """sd(ln k_e) for the POPULATION value after an n-patient cohort fit."""
    tot, _ = cohort_info(n, sched)
    F = tot + np.diag([0.0, 1.0 / POP_PRIOR ** 2])
    return float(np.sqrt(np.linalg.inv(F)[1, 1]))


def sd_patient(n, sched, cv):
    """
    sd(ln k_e_i) for ONE patient once the cohort has been fitted.

    Their prior is the population value, known to `sd_population`, widened by the
    between-patient spread (S2.per_patient_prior does exactly this quadrature);
    their own record then sharpens it. With cv = 0 there is no deviation to learn
    and the answer is the population value itself.
    """
    sd_pop = sd_population(n, sched)
    if cv <= 0.0:
        return sd_pop
    _, first = cohort_info(n, sched)
    pri = S2.per_patient_prior(cv, sd_pop)
    F = first + np.diag([0.0, 1.0 / pri ** 2])
    return float(np.sqrt(np.linalg.inv(F)[1, 1]))


# ==========================================================================
def part_a():
    banner("A. THE POPULATION CONSTANT AGAINST COHORT SIZE")
    print("""  sd(ln k_e) for the population value, each patient's own (ln theta_i, b_i)
  profiled out exactly. Capacities are spread over 5 to 20% of nominal, which is
  where a PA cohort would sit. The trial's own n is 16.""")
    print(f"\n    {'n':>5s}" + "".join(f"{s:>16s}" for s in SCHEDS))
    for n in NS:
        row = []
        for sched in SCHEDS:
            v = sd_population(n, sched)
            row.append(f"{v:14.4f}{'*' if v <= 0.06 else ' '}")
        print(f"    {n:5d}" + "".join(row))
    print("""
  `*` at or below 6%, the looser of S2's two targets for the disease range. The
  tighter is 2%. The cohort size each target needs:""")
    print(f"\n    {'schedule':16s}" +
          "".join(f"{f'n for {100 * t:.0f}%':>16s}" for t in TARGETS))
    for sched in SCHEDS:
        cells = []
        for t in TARGETS:
            need = next((n for n in range(2, 20001)
                         if sd_population(n, sched) <= t), None)
            cells.append(f"{need if need else '>20000':>16}")
        print(f"    {sched:16s}" + "".join(cells))


def part_b():
    banner("B. THE QUESTION THAT DECIDES THE ROUTE: IS THE CONSTANT SHARED")
    print("""  R1 puts the between-patient CV of the disposal constant at "probably a third
  or more", with renal OAT1 secretion and a 50% prevalence of eGFR < 60 in adult
  PA as the reasons. A cohort pins the POPULATION value; a patient's theta_i
  consumes THEIR OWN. Monthly, the schedule a cohort study would plausibly run.""")
    print(f"\n    {'CV':>6s} {'meaning':26s}" +
          "".join(f"{f'n={n}':>12s}" for n in (16, 32, 64, 128)))
    for cv in CVS:
        lab = "shared exactly" if cv == 0 else f"between-patient CV {100*cv:.0f}%"
        cells = "".join(f"{sd_patient(n, 'monthly', cv):12.4f}"
                        for n in (16, 32, 64, 128))
        print(f"    {100*cv:5.0f}% {lab:26s}" + cells)
    print("""
  The first row is what a shared-parameter fit returns. Every other row is what
  ONE patient's constant is known to after the cohort is fitted, which is the
  number their theta needs.""")
    print(f"\n    {'CV':>6s} {'floor at n = 20000':>20s} {'meets 6%':>10s} "
          f"{'meets 2%':>10s}   {'what the floor is':s}")
    for cv in CVS:
        f = sd_patient(20000, "monthly", cv)
        own = (sd_patient(20000, "monthly", cv) if cv == 0 else
               float(np.sqrt(np.linalg.inv(
                   cohort_info(1, "monthly")[1]
                   + np.diag([0.0, 1.0 / np.sqrt(np.log1p(cv ** 2)) ** 2])
               )[1, 1])))
        what = ("the cohort itself" if cv == 0 else
                "their own record against a CV-wide prior")
        print(f"    {100*cv:5.0f}% {f:20.4f} "
              f"{'yes' if f <= 0.06 else 'NO':>10s} "
              f"{'yes' if f <= 0.02 else 'NO':>10s}   {what}")


def simulate_and_fit(n, sched, reps, seed0=0):
    """A cohort simulated and fitted jointly: one shared ln K_m, one shared
    ln k_e, a per-patient (ln theta_i, b_i). Returns the shared ln k_e errors."""
    times = S2.times_for(sched, "no baseline")
    errs = []
    for r in range(reps):
        rng = np.random.default_rng(7717 * (seed0 + r) + 13)
        ths = np.array([TH_NOM * COHORT_CAPS[i % len(COHORT_CAPS)]
                        for i in range(n)])
        pat = []
        for i in range(n):
            L, _, sigma, s_b = S2.footing(ths[i])
            Ls = S2.loads_for(times.size, 1000 * r + i, L)
            inc = S2.increment(times)
            b_t = float(rng.normal(0.0, s_b))
            m = S2.level(ths[i], inc, Ls + b_t)
            obs = ((ths[i] + inc) * m / (KM + m) + KE * m + b_t
                   + rng.normal(0.0, sigma, times.size))
            pat.append((m, obs, inc, sigma, s_b))

        def res(q):
            km, ke = np.exp(q[0]), np.exp(q[1])
            out = []
            for i, (m, obs, inc, sigma, s_b) in enumerate(pat):
                th, b = np.exp(q[2 + 2 * i]), q[3 + 2 * i]
                out.append(((th + inc) * m / (km + m) + ke * m + b - obs)
                           / sigma)
                out.append(np.array([b / s_b]))
            out.append(np.array([(q[1] - np.log(KE)) / POP_PRIOR]))
            return np.concatenate(out)

        q0 = [np.log(KM), np.log(KE)]
        for i in range(n):
            q0 += [np.log(ths[i]), 0.0]
        o = least_squares(res, np.array(q0), method="lm", xtol=1e-12,
                          ftol=1e-12, max_nfev=400000)
        errs.append(float(o.x[1] - np.log(KE)))
    return np.array(errs)


def part_c():
    banner("C. THE HEADLINE CONFIRMED BY FITTING")
    print("""  The standing house rule: a Fisher figure that decides a verdict is confirmed by
  fitting before it is reported. The verdict here is a cohort size, so cohorts are
  simulated and fitted jointly on the stacked residual. The marker level `m` is
  held at its simulated value, as everywhere else in this audit, so the fit sees
  the same design the Fisher block does.""")
    print(f"\n    {'n':>5s} {'schedule':13s} {'Fisher':>9s} {'fitted sd':>10s} "
          f"{'fitted/Fisher':>14s} {'median err':>11s} {'reps':>5s}")
    for n, sched, reps in ((16, "monthly", 24), (16, "weekly 7.3 d", 16),
                           (64, "monthly", 12), (128, "monthly", 8)):
        fis = sd_population(n, sched)
        e = simulate_and_fit(n, sched, reps)
        print(f"    {n:5d} {sched:13s} {fis:9.4f} {e.std(ddof=1):10.4f} "
              f"{e.std(ddof=1) / fis:14.3f} {np.median(e):+11.4f} {reps:5d}",
              flush=True)
    print("""
  A fitted spread near the Fisher value means the cohort figure is a budget and
  not a curvature, which is the check 3a.9 ran on the single-patient design.
  Spread is a sample standard deviation, ddof = 1.""")


def part_d():
    banner("D. WHAT IT DOES TO THE REQUIREMENT TABLE")
    print("""  S2's map read across the k_e prior; this reads down the ROUTES that supply it.
  sd(ln theta) at a recorded load, Fisher route, median over S2's 24 load seeds.
  The micro-dose row is PHASE9 P.8's four studies plus closure; the cohort rows
  need neither.""")
    print(f"\n    {'cap':>5s} {'route for the constant':34s} {'k_e sd':>8s}"
          + "".join(f"{s:>13s}" for s in SCHEDS))
    for cap in (0.05, 0.10, 0.20):
        th = TH_NOM * cap
        rows = [("population prior only", POP_PRIOR),
                ("cohort n=16 monthly, shared", sd_population(16, "monthly")),
                ("cohort n=64 monthly, shared", sd_population(64, "monthly")),
                ("cohort n=64 monthly, CV 33%",
                 sd_patient(64, "monthly", 0.33)),
                ("four micro-dose studies + closure", 0.054)]
        for lab, ks in rows:
            cells = "".join(
                f"{S2.fisher_cell(th, S2.times_for(s, 'no baseline'), ks, 'recorded')[0]:13.3f}"
                for s in SCHEDS)
            print(f"    {100*cap:4.0f}% {lab:34s} {ks:8.4f}" + cells)
    print("""
  Gate is 0.20. The comparison that matters is the last two rows of each block:
  what the cohort delivers against what four tracer studies would.""")


def main():
    print("=" * 104)
    print("  THE HIERARCHICAL COHORT FIT, PRICED")
    print("=" * 104)
    part_a()
    part_b()
    part_c()
    part_d()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
