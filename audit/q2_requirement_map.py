"""
S2. The disease-range requirement map.

What precision on ln theta does a record deliver, across the disease range, as
a function of how well k_e is known, how often blood is drawn, and whether the
substrate load is recorded or measured? And if k_e is known only as a
population value, what does the per-patient prior become?

THE MODEL IS CORRECTLY SPECIFIED HERE. The generator and the analyst share the
committed 3.5 d single-exponential increment, so the amplitude convention does
not arise: `dose`, `peak` and `mean` coincide for the committed kernel. The
misspecification questions belong to S1 and S4.

THE DOSE-DAY DRAW, stated once. A draw at a dose instant is read pre-dose
(`q1_increment_shape.DOSE_DAY_DRAW`). No schedule below puts a draw at a dose
instant except the baseline draw at t = 0, which is therefore untreated.

FOOTING. The review's own, so that S2 extends Gap A rather than restating it:

    capacity   theta in {150, 300, 600, 1500, 3000}, 5 to 100% of normal, in
               n4_theta_scope's convention A: the recorded load stays at its
               theta = 3000 value and the marker rises; sigma, the operating
               level and the load-bias prior are re-solved at each theta.
    dose       the booked D = 6000 at every capacity. NOT admissible at 20%
               of normal and below on Section 7's readings R2 and R3; kept so
               that the map isolates identification from dosing.
    schedules  monthly from 30.4 d, n = 24; two-monthly from 60.8 d, n = 12;
               weekly np.arange(0.5, 730, 7.3), n = 100; every 2 days
               np.arange(1, 730, 2), n = 365. All ON THERAPY: no draw before
               the first infusion. Weekly at 7.3 d rather than 7 so the phases
               spread (a 7 d step reaches two phases of a 14 d cycle).
    baseline   every schedule on two footings: `no baseline`, as above, and
               `baseline`, the same draws plus one untreated draw at t = 0
               taken before the first infusion. The review's monthly record
               had its t = 0 draw read POST-dose, inside the settling window,
               which is neither; Part A0 shows its 3.08 is the no-baseline
               figure. The baseline draw is a design choice, and a patient
               already on therapy cannot have one.
    k_e        a log-normal prior of sd 0.8 (the paper's), 0.20, 0.10, 0.05,
               0.02, or exactly known. Brief R1 (k_e and its between-patient CV)
               finds no human 3-HP clearance of any kind, so this is a sweep.
    load       recorded, with the paper's prior b ~ N(0, 0.05 L), or measured,
               in which case b leaves the parameter vector.
    loads      lognormal day-to-day CV 10%, 24 load records per cell.

Every Fisher figure is a median with its min and max over the 24 load records.
The Fisher bound is not trustworthy near the gate at a disease-range capacity
(the review's own finding), so Part D confirms every cell near 0.20 with a
Monte Carlo estimator: the global K_m scan with an exact polish, with the prior
CENTRES drawn from the priors so that the nuisances are honest.

Nothing here records an expected result in advance.

Run:  python audit/q2_requirement_map.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import brentq, least_squares

import n4_theta_scope as N4
import q1_increment_shape as Q

TH_NOM, KM, KE = 3000.0, 100.0, 5.0
CAPS = (0.05, 0.10, 0.20, 0.50, 1.00)
KE_SDS = (0.8, 0.20, 0.10, 0.05, 0.02, 0.0)          # 0.0 = exactly known
#: the schedules ON THERAPY: no draw before the first infusion. The review's
#: weekly and 2-daily records already had none; its monthly and two-monthly
#: records had a t = 0 draw, which is removed here and restored, for every
#: schedule alike, on the `baseline` footing.
SCHED = {"monthly": np.arange(30.4, 730.0, 30.4),
         "two-monthly": np.arange(60.8, 730.0, 60.8),
         "weekly 7.3 d": np.arange(0.5, 730.0, 7.3),
         "every 2 d": np.arange(1.0, 730.0, 2.0)}
FOOTINGS = ("no baseline", "baseline")


def times_for(sched, foot):
    """The draw times. `baseline` adds one untreated draw at t = 0, before the
    first infusion (it is read pre-dose)."""
    base = SCHED[sched]
    return base if foot == "no baseline" else np.concatenate([[0.0], base])
LOADS = ("recorded", "measured")
SEEDS = tuple(range(24))
GATE = 0.20
LOAD_CV = 0.10

#: the committed increment, assembled from the committed kernel and train, read
#: per DOSE_DAY_DRAW. Q.make_analysts' own route, checked against
#: design.run.delta in q1 Part 0.
_TN, _WN = np.arange(0.0, 900.0, 14.0), np.ones(int(np.ceil(900.0 / 14.0)))
_K = Q.kernel_exp(Q.TH_HALF0)


def increment(times):
    return Q.D_AMP * Q.assemble(_K, _TN, _WN, Q.at_draw(times))


def banner(t):
    print("\n" + "=" * 100)
    print("  " + t)
    print("=" * 100)


def footing(theta):
    """(L, m_bar, sigma, s_b) at this capacity, convention A."""
    L, m_bar, sigma, s_b = N4.footing(theta, KM, KE, convention="A")
    return float(L), float(m_bar), float(sigma), float(s_b)


def loads_for(n, seed, L):
    rng = np.random.default_rng(seed)
    s2 = np.log1p(LOAD_CV ** 2)
    return L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), n))


def level(theta, inc, L, Km=KM, ke=KE):
    B = ke * Km + theta + inc - L
    return (-B + np.sqrt(B ** 2 + 4.0 * ke * L * Km)) / (2.0 * ke)


# ==========================================================================
#  The information matrix in (ln theta, ln K_m, b, ln k_e).
# ==========================================================================
def fisher_sd(theta, times, seed, ke_sd, load):
    L, _, sigma, s_b = footing(theta)
    Ls = loads_for(times.size, seed, L)
    inc = increment(times)
    m = level(theta, inc, Ls)
    g = m / (KM + m)
    cols = [theta * g, -(theta + inc) * KM * m / (KM + m) ** 2,
            np.ones_like(m), KE * m]
    pri = [0.0, 0.0, 1.0 / s_b ** 2, None if ke_sd == 0 else 1.0 / ke_sd ** 2]
    keep = [0, 1]
    if load == "recorded":
        keep.append(2)
    if ke_sd != 0:
        keep.append(3)
    Qm = np.stack([cols[i] for i in keep], axis=1)
    F = np.diag([pri[i] for i in keep]) + Qm.T @ Qm / sigma ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 1e-14 * w[-1]:
        return np.inf
    return float(np.sqrt(np.linalg.inv(F)[0, 0]))


def fisher_cell(theta, times, ke_sd, load, seeds=SEEDS):
    v = np.array([fisher_sd(theta, times, s, ke_sd, load) for s in seeds])
    v = v[np.isfinite(v)]
    if not v.size:
        return np.inf, np.inf, np.inf
    return float(np.median(v)), float(v.min()), float(v.max())


def fmt(c, w=20):
    med, lo, hi = c
    if not np.isfinite(med):
        return f"{'inf':>{w}s}"
    star = "*" if med <= GATE else " "
    return f"{med:6.3f}{star}[{lo:5.3f},{hi:6.3f}]".rjust(w)


# ==========================================================================
def part_a0():
    banner("A0. THE FOOTING, REPRODUCED BEFORE IT IS EXTENDED")
    print("  The review's figures, and what its monthly t = 0 draw was.")
    mo = np.arange(0.0, 730.0, 30.4)
    with Q.rebind(Q, DOSE_DAY_DRAW="post-dose"):
        rows_post = [
            ("paper headline, theta 3000, monthly, k_e 0.8, b prior "
             "[t=0 post-dose]", TH_NOM, mo, 0.8, "recorded"),
            ("review Gap A, theta 150, monthly, k_e 0.8, b prior "
             "[t=0 post-dose]", 150.0, mo, 0.8, "recorded")]
        for lab, th, tt, ks, ld in rows_post:
            print(f"    {lab:70s} {fmt(fisher_cell(th, tt, ks, ld))}")
    rows = [("review Gap A, theta 150, monthly, k_e 0.8, b prior, NO baseline",
             150.0, times_for("monthly", "no baseline"), 0.8, "recorded"),
            ("the same WITH a pre-treatment baseline draw",
             150.0, times_for("monthly", "baseline"), 0.8, "recorded"),
            ("review Gap A, theta 150, weekly, k_e 2%, load measured",
             150.0, SCHED["weekly 7.3 d"], 0.02, "measured"),
            ("review Gap A, theta 150, every 2 d, k_e 2%, load measured",
             150.0, SCHED["every 2 d"], 0.02, "measured")]
    for lab, th, tt, ks, ld in rows:
        print(f"    {lab:70s} {fmt(fisher_cell(th, tt, ks, ld))}")


def part_a():
    banner("A. sd(ln theta_hat), FISHER, BY CAPACITY x k_e x SCHEDULE x LOAD "
           "x BASELINE")
    print("  Median [min, max] over 24 load records. `*` marks a median at or "
          "below the 0.20 gate.")
    print("  k_e precision is the sd of the log-normal prior; `exact` is known "
          "exactly.")
    out = {}
    for foot in FOOTINGS:
        print(f"\n  ===== FOOTING: {foot.upper()} =====")
        for cap in CAPS:
            th = TH_NOM * cap
            L, m_bar, sigma, s_b = footing(th)
            print(f"\n  CAPACITY {100 * cap:.0f}% of normal, theta = {th:.0f}"
                  f"   (m_bar {m_bar:.2f}, sigma {sigma:.3f}, b prior "
                  f"{s_b:.3f})")
            print(f"    {'k_e':>6s} {'load':9s}" +
                  "".join(f"{k:>21s}" for k in SCHED))
            for ks in KE_SDS:
                for ld in LOADS:
                    cells = []
                    for sn in SCHED:
                        c = fisher_cell(th, times_for(sn, foot), ks, ld)
                        out[(foot, cap, ks, ld, sn)] = c
                        cells.append(fmt(c, 21))
                    lab = "exact" if ks == 0 else f"{ks:.2f}"
                    print(f"    {lab:>6s} {ld:9s}" + "".join(cells))
    return out


def part_b(out):
    banner("B. THE k_e PRECISION EACH DESIGN NEEDS TO REACH 0.20")
    print("  The prior sd on ln k_e at which the median Fisher bound crosses "
          "0.20, by root-finding")
    print("  on the same 24 load records. `never`: not even with k_e exactly "
          "known. `any`: the")
    print("  paper's 0.8 already suffices.")
    req = {}
    for foot in FOOTINGS:
        print(f"\n  FOOTING: {foot}")
        print(f"    {'capacity':>9s} {'load':9s}" +
              "".join(f"{k:>15s}" for k in SCHED))
        for cap in CAPS:
            th = TH_NOM * cap
            for ld in LOADS:
                cells = []
                for sn in SCHED:
                    tt = times_for(sn, foot)
                    exact = out[(foot, cap, 0.0, ld, sn)][0]
                    loose = out[(foot, cap, 0.8, ld, sn)][0]
                    if exact > GATE:
                        r, txt = 0.0, "never"
                    elif loose <= GATE:
                        r, txt = 0.8, "any"
                    else:
                        f = lambda lk: fisher_cell(th, tt, float(np.exp(lk)),
                                                   ld)[0] - GATE
                        r = float(np.exp(brentq(f, np.log(1e-4), np.log(0.8),
                                                xtol=1e-4)))
                        txt = f"{100 * r:.1f}%"
                    req[(foot, cap, ld, sn)] = r
                    cells.append(f"{txt:>15s}")
                print(f"    {100 * cap:8.0f}% {ld:9s}" + "".join(cells))
    return req


def per_patient_prior(cv_between, sd_mean=0.0):
    """
    The per-patient log-sd when k_e is known only as a population value.

    A patient's ln k_e is the population's ln mean plus a between-patient
    deviation of log-sd sqrt(ln(1 + CV^2)), plus the uncertainty in the
    population mean itself, `sd_mean`, which is independent of it. Nothing
    about the individual is known, so the two add in quadrature.
    """
    return float(np.sqrt(np.log1p(cv_between ** 2) + sd_mean ** 2))


def part_c():
    banner("C. k_e KNOWN ONLY AS A POPULATION VALUE")
    print("""  If k_e is carried as the population's value, the individual's
  deviation from it is the prior. Between-patient CV 0, 10, 20 and 40%, with the
  population mean itself known exactly or to 5% or 10%. Brief R1 puts the
  between-patient CV of a comparable disposal rate at about a third in healthy
  adults (propionate turnover, SD 5.9 on a mean of 17.6) and gives reasons,
  renal OAT1 secretion and a 50% prevalence of eGFR < 60 in adult PA, to expect
  more in patients. It gives no 3-HP figure at all.""")
    print(f"\n    {'between-patient CV':>19s} {'pop. mean sd':>13s} "
          f"{'per-patient prior':>18s}")
    for cv in (0.0, 0.10, 0.20, 0.40):
        for sm in (0.0, 0.05, 0.10):
            print(f"    {100 * cv:18.0f}% {sm:13.2f} "
                  f"{per_patient_prior(cv, sm):18.4f}")
    out = {}
    for foot in FOOTINGS:
        print(f"\n  sd(ln theta_hat) at the per-patient prior, population mean "
              f"known exactly, load measured, {foot}:")
        print(f"    {'capacity':>9s} {'CV':>5s}" +
              "".join(f"{k:>21s}" for k in SCHED))
        for cap in CAPS:
            for cv in (0.0, 0.10, 0.20, 0.40):
                ks = per_patient_prior(cv)
                cells = []
                for sn in SCHED:
                    c = fisher_cell(TH_NOM * cap, times_for(sn, foot), ks,
                                    "measured")
                    out[(foot, cap, cv, sn)] = c
                    cells.append(fmt(c, 21))
                print(f"    {100 * cap:8.0f}% {100 * cv:4.0f}%" +
                      "".join(cells))
    return out


# ==========================================================================
#  The Monte Carlo estimator: global K_m scan, exact polish, honest nuisances.
# ==========================================================================
KM_GRID = np.geomspace(1.0, 1000.0, 401)


def mc_fit(m, inc, s, sigma, s_b, ke_sd, ke_centre, load):
    """
    MAP fit of (theta, K_m[, b][, k_e]) with the stated priors.

    Scan K_m on a grid, solve the linear block exactly at each node with the
    delta-method k_e prior, then polish the EXACT objective (prior on ln k_e)
    inside the scan's box. The structure of q1's estimator, with b and k_e
    optionally known.
    """
    free_b = load == "recorded"
    free_ke = ke_sd != 0
    g = m[None, :] / (KM_GRID[:, None] + m[None, :])
    cols, sd, mu = [g], [None], [0.0]
    if free_ke:
        cols.append(np.broadcast_to(m, g.shape))
        sd.append(ke_sd * ke_centre)
        mu.append(ke_centre)
    if free_b:
        cols.append(np.ones_like(g))
        sd.append(s_b)
        mu.append(0.0)
    y = s[None, :] - inc[None, :] * g
    if not free_ke:
        y = y - ke_centre * m[None, :]
    X = np.stack(cols, axis=-1)
    p, J = Q._linear_solve(X, y, sigma, sd, mu)
    j = int(np.argmin(J))
    th0 = max(float(p[j, 0]), 1e-6)
    ke0 = float(p[j, 1]) if free_ke else ke_centre
    ke0 = ke0 if ke0 > 1e-3 * ke_centre else ke_centre
    b0 = float(p[j, -1]) if free_b else 0.0
    q0 = [np.log(th0), np.log(KM_GRID[j])]
    lb, ub = [-np.inf, 0.0], [np.inf, np.log(1000.0)]
    if free_ke:
        q0.append(np.log(ke0)); lb.append(-np.inf); ub.append(np.inf)
    if free_b:
        q0.append(b0); lb.append(-np.inf); ub.append(np.inf)

    def res(q):
        th, K = np.exp(q[0]), np.exp(q[1])
        i = 2
        ke = np.exp(q[i]) if free_ke else ke_centre
        i += free_ke
        b = q[i] if free_b else 0.0
        r = [((th + inc) * m / (K + m) + ke * m + b - s) / sigma]
        if free_ke:
            r.append([(q[2] - np.log(ke_centre)) / ke_sd])
        if free_b:
            r.append([b / s_b])
        return np.concatenate(r)

    q0 = np.clip(q0, np.array(lb) + 1e-9, np.array(ub) - 1e-9)
    out = least_squares(res, q0, bounds=(lb, ub), method="trf",
                        max_nfev=6000, xtol=1e-13, ftol=1e-13, gtol=1e-13)
    return float(np.exp(out.x[0]))


#: a fit is a COLLAPSE when it misses capacity by more than this factor. Not a
#: bound on anything and nothing is clipped to it: a label, because the RMSE
#: of an estimator with a heavy tail is dominated by its rare collapses and
#: says little about the rest.
COLLAPSE = np.log(4.5)


def mc_cell(theta, times, ke_sd, load, nuisance, reps=96, seed0=0):
    """
    Error of ln theta_hat over `reps` records, with honest nuisances.

    Each record has its own load record and noise, and the nuisances are drawn
    from the priors the estimator carries, in one of two ways that coincide
    for a narrow prior and part company for a wide one:

      "truth"   the patient's TRUE k_e is drawn from the prior, whose centre
                stays at the reference: a POPULATION prior, the patient being
                a draw from the population. i1_estimator's convention.
      "centre"  the true k_e stays at the reference and the prior's CENTRE is
                drawn about it: a per-patient MEASUREMENT of k_e with that
                error, which is what S2b's micro-dose study would supply.

    The load bias b, when recorded, is drawn from its prior about zero in both.
    """
    L, m_bar, sigma, s_b = footing(theta)
    inc = increment(times)
    err = []
    for r in range(reps):
        sd_ = seed0 + r
        rng = np.random.default_rng(10_007 * sd_ + 3)
        Ls = loads_for(times.size, sd_, L)
        z = rng.normal(0.0, ke_sd) if ke_sd else 0.0
        if nuisance == "truth":
            ke_t, ke_c = KE * np.exp(z), KE
        else:
            ke_t, ke_c = KE, KE * np.exp(z)
        b_t = rng.normal(0.0, s_b) if load == "recorded" else 0.0
        m = level(theta, inc, Ls + b_t, ke=ke_t)
        s = Ls + rng.normal(0.0, sigma, times.size)
        th = mc_fit(m, inc, s, sigma, s_b, ke_sd, ke_c, load)
        err.append(np.log(th) - np.log(theta))
    e = np.array(err)
    boot = np.random.default_rng(99).choice(e, (400, e.size))
    rb = np.sqrt(np.mean(boot ** 2, axis=1))
    return dict(sd=float(e.std()), rmse=float(np.sqrt(np.mean(e ** 2))),
                bias=float(np.median(e)), rmse_lo=float(np.percentile(rb, 5)),
                rmse_hi=float(np.percentile(rb, 95)),
                robust=float((np.percentile(e, 84) - np.percentile(e, 16)) / 2),
                collapse=int(np.sum(np.abs(e) > COLLAPSE)), n=e.size)


def primary_nuisance(ke_sd):
    """The paper's 0.8 is a population prior; anything tighter is taken to be
    a per-patient measurement, which is where S2b's precisions come from."""
    return "truth" if ke_sd >= 0.8 else "centre"


def verdict(r):
    return ("MEETS" if r["rmse_hi"] <= GATE else
            "FAILS" if r["rmse_lo"] > GATE else "edge")


def part_d(out, band=(0.10, 0.40)):
    banner("D. MONTE CARLO AT EVERY CELL NEAR THE GATE")
    print(f"  Every cell of Part A whose Fisher median lies in [{band[0]}, "
          f"{band[1]}], re-estimated by the")
    print("  MAP estimator on 96 records with honest nuisances, drawn BOTH "
          "ways (mc_cell). The")
    print("  primary way is marked `>`: `truth` for the paper's 0.8, a "
          "population prior; `centre`")
    print("  for anything tighter, a per-patient measurement. RMSE carries a "
          "bootstrap 90% interval;")
    print("  `robust` is half the 16-84 percentile range; `coll` counts fits "
          "off by more than 4.5x.")
    rows = [(k, v) for k, v in out.items() if band[0] <= v[0] <= band[1]]
    print(f"  {len(rows)} cells qualify.\n")
    print(f"    {'footing':11s} {'cap':>4s} {'k_e':>5s} {'load':9s} "
          f"{'schedule':12s} {'Fisher':>6s} | {'truth-drawn RMSE [90%]':>24s} "
          f"{'rob':>5s} {'cl':>2s} | {'centre-drawn RMSE [90%]':>24s} "
          f"{'rob':>5s} {'cl':>2s} | {'verdict':>7s}")
    res = {}
    for (foot, cap, ks, ld, sn), c in sorted(rows):
        tt = times_for(sn, foot)
        pair = {nu: mc_cell(TH_NOM * cap, tt, ks, ld, nu)
                for nu in ("truth", "centre")}
        res[(foot, cap, ks, ld, sn)] = pair
        prim = primary_nuisance(ks)
        cells = []
        for nu in ("truth", "centre"):
            r = pair[nu]
            mark = ">" if nu == prim else " "
            cells.append(f"{mark}{r['rmse']:6.3f} [{r['rmse_lo']:5.3f},"
                         f"{r['rmse_hi']:6.3f}] {r['robust']:5.3f} "
                         f"{r['collapse']:2d}")
        lab = "exact" if ks == 0 else f"{ks:.2f}"
        print(f"    {foot:11s} {100 * cap:3.0f}% {lab:>5s} {ld:9s} {sn:12s} "
              f"{c[0]:6.3f} | {cells[0]} | {cells[1]} | "
              f"{verdict(pair[prim]):>7s}")
    return res


def part_e(out_c, band=(0.10, 0.40)):
    banner("E. MONTE CARLO FOR THE POPULATION-VALUE PRIORS OF PART C")
    print("  A population value is a population prior by construction, so the "
          "nuisance is drawn")
    print("  `truth`-wise: each patient's k_e is a draw about the population "
          "value. Cells whose")
    print(f"  Fisher median lies in [{band[0]}, {band[1]}].\n")
    rows = [(k, v) for k, v in out_c.items() if band[0] <= v[0] <= band[1]]
    print(f"  {len(rows)} cells qualify.\n")
    print(f"    {'footing':11s} {'cap':>4s} {'CV':>4s} {'schedule':12s} "
          f"{'Fisher':>6s} {'RMSE [90%]':>22s} {'robust':>7s} {'coll':>5s} "
          f"{'verdict':>8s}")
    res = {}
    for (foot, cap, cv, sn), c in sorted(rows):
        r = mc_cell(TH_NOM * cap, times_for(sn, foot), per_patient_prior(cv),
                    "measured", "truth")
        res[(foot, cap, cv, sn)] = r
        print(f"    {foot:11s} {100 * cap:3.0f}% {100 * cv:3.0f}% {sn:12s} "
              f"{c[0]:6.3f} {r['rmse']:6.3f} [{r['rmse_lo']:5.3f},"
              f"{r['rmse_hi']:6.3f}] {r['robust']:7.3f} {r['collapse']:5d} "
              f"{verdict(r):>8s}")
    return res


def main(argv=()):
    print("=" * 100)
    print("  S2. THE DISEASE-RANGE REQUIREMENT MAP")
    print("=" * 100)
    print(f"  Dose-day draws read {Q.DOSE_DAY_DRAW}. The model is correctly "
          "specified throughout.")
    part_a0()
    out = part_a()
    part_b(out)
    out_c = part_c()
    part_d(out)
    part_e(out_c)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(tuple(sys.argv[1:])))
