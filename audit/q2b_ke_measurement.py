"""
S2b. The k_e measurement itself.

`paper/sections/08_specification.tex` proposes, and prices through S2, a
micro-dose intravenous bolus of 13C3-3-hydroxypropionate with serial plasma
sampling and timed urine collection. This script asks what per-patient
precision on k_e that study delivers, and feeds the answer back into S2.

NOTHING PINS THE TRACER'S KINETICS. Brief R1 (k_e and its between-patient CV)
finds no human 3-HP half-life, clearance, volume, or renal fraction in any
population. What it does find are bounds by analogy: blood propionate has a
half-life of 4 to 5 min in controls and 3 to 7 times that in PA; methylmalonic
acid distributes in about total body water; hepatic 3-HP uptake in the perfused
rat liver is a fifth of propionate's. So the half-life is SWEPT, from 5 min to
6 h, and so are everything the design controls.

THE STUDY, as simulated.

    plasma     a bolus at t = 0; N draws, log-spaced from 2 min to the last
               draw the clinic day allows (2 h, 6 h or 24 h); tracer measured
               with a proportional assay error of CV c; draws whose true
               concentration has fallen below the assay's LOQ are lost. The
               bolus is characterised by C0/LOQ, the multiple of the LOQ the
               first concentration sits at, which is what a micro-dose buys.
    model      one compartment, parameters (ln CL, ln V); a two-compartment
               truth and analysis as a sensitivity.
    urine      a timed collection to T_u; the tracer recovered estimates the
               renal clearance as recovered amount over plasma AUC, with a
               collection error of CV u (completeness and assay together).

WHAT THE STUDY MEASURES DEPENDS ON A READING THE PAPER DOES NOT MAKE. In the
paper's balance a tracer at trace level is removed by both routes, at a total
clearance CL = k_e + theta/(K_m + m): the alternative channel and, through
the enzyme's own saturation, the enzymatic one. So:

    reading (ii), biological   3-HP is an end product and not a PCC
                               substrate, so the tracer measures k_e itself
                               (brief R1 argues this: "measured total
                               clearance approximates non-enzymatic k_e").
                               The study's sd(ln CL) becomes the prior on
                               ln k_e.
    reading (i), model-literal the tracer measures k_e + theta/(K_m + m),
                               which is a joint observation on theta and k_e
                               and enters S2's information matrix as a row of
                               its own, not as a prior.

Both are fed back into S2. A third term is added in quadrature because it is
not a property of the study at all: k_e drifting WITHIN a patient between the
measurement and the two-year record (renal decline is the obvious mechanism;
brief R1 reports eGFR under 60 in half of adults with PA). Swept at 0, 5 and
10%.

Nothing here records an expected result in advance.

Run:  python audit/q2b_ke_measurement.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import least_squares

import q2_requirement_map as S2

HALF_LIVES_MIN = (5.0, 15.0, 30.0, 60.0, 120.0, 360.0)
WINDOWS_H = (2.0, 6.0, 24.0)
NS = (4, 6, 8, 12, 16)
CVS = (0.05, 0.10, 0.15)
C0_LOQ = (30.0, 100.0, 1000.0)
FIRST_MIN = 2.0


def banner(t):
    print("\n" + "=" * 100)
    print("  " + t)
    print("=" * 100)


def draw_times(n, window_h, first_min=FIRST_MIN):
    """n draws log-spaced from `first_min` minutes to the end of the window,
    in HOURS."""
    return np.geomspace(first_min / 60.0, window_h, n)


# ==========================================================================
#  The plasma arm.
# ==========================================================================
def one_comp_sd(t_half_h, times, cv, c0_loq):
    """
    sd(ln CL) from a one-compartment fit, parameters (ln CL, ln V).

    ln c(t) = ln(D/V) - (CL/V) t. With lambda = CL/V the sensitivities are
    d ln c / d ln CL = -lambda t and d ln c / d ln V = -1 + lambda t, and the
    error is proportional, so each draw carries 1/cv^2. Draws below the LOQ
    are lost. Returns (sd, draws kept).
    """
    lam = np.log(2.0) / t_half_h
    keep = np.exp(-lam * times) * c0_loq >= 1.0
    t = times[keep]
    if t.size < 2:
        return np.inf, int(t.size)
    J = np.stack([-lam * t, -1.0 + lam * t], axis=1)
    F = J.T @ J / cv ** 2
    if np.linalg.cond(F) > 1e12:
        return np.inf, int(t.size)
    return float(np.sqrt(np.linalg.inv(F)[0, 0])), int(t.size)


def two_comp_conc(t, CL, V1, Q, V2):
    """A unit IV bolus into the central compartment of a two-compartment
    model; the central concentration."""
    k10, k12, k21 = CL / V1, Q / V1, Q / V2
    a = k10 + k12 + k21
    disc = np.sqrt(a * a - 4.0 * k10 * k21)
    alpha, beta = 0.5 * (a + disc), 0.5 * (a - disc)
    A = (alpha - k21) / (V1 * (alpha - beta))
    B = (k21 - beta) / (V1 * (alpha - beta))
    return A * np.exp(-alpha * t) + B * np.exp(-beta * t)


def two_comp_sd(t_half_h, times, cv, c0_loq, v2_over_v1=1.0, q_over_cl=1.0):
    """
    sd(ln CL) from a two-compartment fit, parameters (ln CL, ln V1, ln Q,
    ln V2), with the truth a two-compartment system whose TERMINAL half-life
    is `t_half_h`. V1 = 1 fixes the scale; CL is solved so the terminal rate
    matches.
    """
    V1, V2 = 1.0, v2_over_v1

    def terminal(CL):
        Q = q_over_cl * CL
        k10, k12, k21 = CL / V1, Q / V1, Q / V2
        a = k10 + k12 + k21
        return 0.5 * (a - np.sqrt(a * a - 4.0 * k10 * k21))
    lo, hi = 1e-6, 1e6
    target = np.log(2.0) / t_half_h
    for _ in range(200):
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if terminal(mid) < target else (lo, mid)
    CL = np.sqrt(lo * hi)
    p0 = np.log([CL, V1, q_over_cl * CL, V2])
    c_true = two_comp_conc(times, *np.exp(p0))
    keep = c_true / c_true.max() * c0_loq >= 1.0
    t = times[keep]
    if t.size < 4:
        return np.inf, int(t.size)
    h = 1e-6
    J = np.empty((t.size, 4))
    for j in range(4):
        dp = np.zeros(4)
        dp[j] = h
        J[:, j] = (np.log(two_comp_conc(t, *np.exp(p0 + dp)))
                   - np.log(two_comp_conc(t, *np.exp(p0 - dp)))) / (2 * h)
    F = J.T @ J / cv ** 2
    if np.linalg.cond(F) > 1e12:
        return np.inf, int(t.size)
    return float(np.sqrt(np.linalg.inv(F)[0, 0])), int(t.size)


def one_comp_mc(t_half_h, times, cv, c0_loq, reps=400, seed=0):
    """The Fisher figure confirmed by fitting: sd, RMSE and bias of ln CL_hat."""
    lam = np.log(2.0) / t_half_h
    V = 1.0
    CL = lam * V
    keep = np.exp(-lam * times) * c0_loq >= 1.0
    t = times[keep]
    rng = np.random.default_rng(seed)
    err = []
    for _ in range(reps):
        y = np.log(np.exp(-lam * t) / V) + rng.normal(0.0, cv, t.size)

        def res(p):
            return (np.log(1.0 / np.exp(p[1])) - np.exp(p[0] - p[1]) * t
                    - y) / cv
        out = least_squares(res, [np.log(CL) + 0.2, 0.1], method="lm")
        err.append(out.x[0] - np.log(CL))
    e = np.array(err)
    return float(e.std()), float(np.sqrt(np.mean(e ** 2))), float(np.mean(e))


def part_a():
    banner("A. THE PLASMA ARM: sd(ln CL) FROM THE BOLUS")
    print("  One compartment. Draws log-spaced from 2 min to the end of the "
          "window. Each cell is")
    print("  sd(ln CL) and, in brackets, how many draws survive the LOQ. The "
          "sd is scale-free in")
    print("  the dose; what the dose buys is how long the tracer stays above "
          "the LOQ.")
    out = {}
    for c0 in C0_LOQ:
        for cv in CVS:
            print(f"\n  C0/LOQ = {c0:g}, assay CV {100 * cv:.0f}%")
            print(f"    {'t_half':>7s} {'window':>7s}" +
                  "".join(f"{f'N={n}':>14s}" for n in NS))
            for th in HALF_LIVES_MIN:
                for w in WINDOWS_H:
                    cells = []
                    for n in NS:
                        sd, k = one_comp_sd(th / 60.0, draw_times(n, w), cv, c0)
                        out[(c0, cv, th, w, n)] = (sd, k)
                        cells.append(f"{sd:8.4f} ({k:2d})" if np.isfinite(sd)
                                     else f"{'inf':>8s} ({k:2d})")
                    print(f"    {th:5.0f}mn {w:5.0f} h" +
                          "".join(f"{c:>14s}" for c in cells))
    return out


def part_a2():
    banner("A2. THE FISHER FIGURE CONFIRMED BY FITTING, AND A TWO-COMPARTMENT "
           "TRUTH")
    print("  (a) One-compartment fits to 400 simulated studies against the "
          "Fisher sd, a few cells.")
    for th, w, n, cv, c0 in ((5.0, 2.0, 8, 0.10, 100.0),
                             (30.0, 6.0, 8, 0.10, 100.0),
                             (120.0, 24.0, 12, 0.05, 1000.0),
                             (360.0, 6.0, 6, 0.10, 100.0)):
        t = draw_times(n, w)
        f, _ = one_comp_sd(th / 60.0, t, cv, c0)
        sd, rmse, bias = one_comp_mc(th / 60.0, t, cv, c0)
        print(f"    t1/2 {th:4.0f} min, {w:4.0f} h, N {n:2d}, CV {100 * cv:2.0f}%,"
              f" C0/LOQ {c0:5.0f}: Fisher {f:.4f}  MC sd {sd:.4f} "
              f"RMSE {rmse:.4f} bias {bias:+.4f}")
    print("\n  (b) A two-compartment truth, analysed two-compartment. V2 = V1, "
          "Q = CL, terminal")
    print("  half-life as stated. Assay CV 10%, C0/LOQ 100.")
    print(f"    {'t_half':>7s} {'window':>7s}" +
          "".join(f"{f'N={n}':>14s}" for n in NS) + "   (one-comp, same cell)")
    for th in HALF_LIVES_MIN:
        for w in WINDOWS_H:
            cells = []
            for n in NS:
                sd, k = two_comp_sd(th / 60.0, draw_times(n, w), 0.10, 100.0)
                cells.append(f"{sd:8.4f} ({k:2d})" if np.isfinite(sd)
                             else f"{'inf':>8s} ({k:2d})")
            one, _ = one_comp_sd(th / 60.0, draw_times(12, w), 0.10, 100.0)
            print(f"    {th:5.0f}mn {w:5.0f} h" +
                  "".join(f"{c:>14s}" for c in cells) + f"   {one:.4f}")


# ==========================================================================
#  The urine arm.
# ==========================================================================
def part_b(out):
    banner("B. THE URINE ARM: RENAL CLEARANCE, AND WHAT IT ADDS")
    print("""  Renal clearance is the tracer recovered in urine over the plasma AUC.
  Collected to at least five half-lives, the amount is complete to 3%, so the
  error is the collection's own, CV u, and the plasma AUC's, which is
  sd(ln CL) from Part A because the dose is known. The two are independent:
      sd(ln CL_renal) = sqrt(u^2 + sd(ln CL)^2).
  If k_e is entirely renal, which brief R1 cannot settle (the renal fraction
  is NOT FOUND), the urine arm measures k_e and the plasma arm measures total
  clearance, and their ratio is the renal fraction. If k_e has a non-renal
  part, the urine arm is a lower bound on k_e and not a measurement of it.""")
    print(f"\n    {'plasma sd(ln CL)':>17s}" +
          "".join(f"{f'u={int(100 * u)}%':>10s}" for u in (0.03, 0.05, 0.10,
                                                          0.20)))
    for sd in (0.01, 0.02, 0.03, 0.05, 0.10):
        print(f"    {sd:17.2f}" +
              "".join(f"{np.hypot(u, sd):10.4f}" for u in (0.03, 0.05, 0.10,
                                                           0.20)))


# ==========================================================================
#  What it buys: S2 again, with the study's precision fed in.
# ==========================================================================
def fisher_with_tracer(theta, times, seed, load, sd_tracer, drift,
                       ke_pop=0.8):
    """
    S2's information matrix with the tracer study as a ROW, reading (i).

    The study is run before the first infusion, at the untreated level m0, and
    observes ln CL_total = ln(k_e + theta/(K_m + m0)) with sd `sd_tracer`. The
    within-patient drift enters as extra variance on that row, since what the
    record needs is k_e over two years and not on the day of the study.

    The POPULATION prior on ln k_e is kept, at `ke_pop`, S2's widest: it is a
    property of the literature and not of the study, so reading (i) does not
    abolish it, it merely declines to replace it. Without it the comparison
    would charge reading (i) for information reading (ii) never created.
    """
    L, _, sigma, s_b = S2.footing(theta)
    Ls = S2.loads_for(times.size, seed, L)
    inc = S2.increment(times)
    m = S2.level(theta, inc, Ls)
    g = m / (S2.KM + m)
    cols = [theta * g, -(theta + inc) * S2.KM * m / (S2.KM + m) ** 2,
            np.ones_like(m), S2.KE * m]
    keep = [0, 1, 3] + ([2] if load == "recorded" else [])
    keep.sort()
    Q = np.stack([cols[i] for i in keep], axis=1)
    F = Q.T @ Q / sigma ** 2
    if load == "recorded":
        F[keep.index(2), keep.index(2)] += 1.0 / s_b ** 2
    F[keep.index(3), keep.index(3)] += 1.0 / ke_pop ** 2
    m0 = S2.level(theta, 0.0, L)
    CLe = theta / (S2.KM + m0)
    CL = S2.KE + CLe
    row = np.zeros(len(keep))
    row[keep.index(0)] = CLe / CL
    row[keep.index(1)] = -theta * S2.KM / (S2.KM + m0) ** 2 / CL
    row[keep.index(3)] = S2.KE / CL
    s2 = sd_tracer ** 2 + drift ** 2
    F = F + np.outer(row, row) / s2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 1e-14 * w[-1]:
        return np.inf
    return float(np.sqrt(np.linalg.inv(F)[0, 0]))


def part_c():
    banner("C. WHAT THE STUDY BUYS, FED BACK INTO S2")
    print("""  sd(ln theta_hat), Fisher, median over 24 load records, with the tracer
  study's precision fed in two ways. Reading (ii): the study measures k_e, so
  its sd, with the within-patient drift in quadrature, is the prior on ln k_e.
  Reading (i): it measures k_e + theta/(K_m + m0) at the untreated level, a
  row of its own. Without a baseline draw and with one; `*` at or below 0.20.

  The share of total clearance the enzyme carries at the untreated level,
  which is what separates the two readings:""")
    for cap in S2.CAPS:
        th = S2.TH_NOM * cap
        L = S2.footing(th)[0]
        m0 = S2.level(th, 0.0, L)
        CLe = th / (S2.KM + m0)
        print(f"    {100 * cap:4.0f}% of normal: m0 {m0:6.2f}, enzymatic "
              f"{CLe:7.3f}, k_e {S2.KE:.1f}, enzymatic share "
              f"{100 * CLe / (CLe + S2.KE):5.1f}%")
    out = {}
    for foot in S2.FOOTINGS:
        print(f"\n  FOOTING: {foot}")
        print(f"    {'cap':>4s} {'load':9s} {'study':>6s} {'drift':>6s} "
              f"{'reading':>8s}" +
              "".join(f"{k:>21s}" for k in ("monthly", "weekly 7.3 d",
                                             "every 2 d")))
        for cap in (0.05, 0.10, 0.20):
            th = S2.TH_NOM * cap
            for ld in S2.LOADS:
                for sd_t in (0.02, 0.05, 0.10):
                    for drift in (0.0, 0.05, 0.10):
                        for rd in ("(ii)", "(i)"):
                            cells = []
                            for sn in ("monthly", "weekly 7.3 d", "every 2 d"):
                                tt = S2.times_for(sn, foot)
                                if rd == "(ii)":
                                    c = S2.fisher_cell(th, tt,
                                                       float(np.hypot(sd_t,
                                                                      drift)),
                                                       ld)
                                else:
                                    v = np.array([fisher_with_tracer(
                                        th, tt, s, ld, sd_t, drift)
                                        for s in S2.SEEDS])
                                    c = (float(np.median(v)), float(v.min()),
                                         float(v.max()))
                                out[(foot, cap, ld, sd_t, drift, rd, sn)] = c
                                cells.append(S2.fmt(c, 21))
                            print(f"    {100 * cap:3.0f}% {ld:9s} "
                                  f"{100 * sd_t:5.0f}% {100 * drift:5.0f}% "
                                  f"{rd:>8s}" + "".join(cells))
    return out


def part_d(out, band=(0.10, 0.40)):
    banner("D. MONTE CARLO FOR READING (ii) NEAR THE GATE")
    print("  Reading (ii) makes the study a prior on ln k_e, which S2's MAP "
          "estimator carries")
    print("  directly; the prior's centre is drawn about the truth with the "
          "study's own sd and")
    print("  the drift, i.e. `centre`-wise, because it is a measurement. "
          "Cells in the band only.")
    rows = [(k, v) for k, v in out.items()
            if k[5] == "(ii)" and band[0] <= v[0] <= band[1]]
    print(f"  {len(rows)} cells qualify.\n")
    print(f"    {'footing':11s} {'cap':>4s} {'load':9s} {'study':>6s} "
          f"{'drift':>6s} {'schedule':13s} {'Fisher':>7s} "
          f"{'RMSE [90%]':>22s} {'robust':>7s} {'coll':>5s} {'verdict':>8s}")
    for (foot, cap, ld, sd_t, drift, rd, sn), c in sorted(rows):
        ks = float(np.hypot(sd_t, drift))
        r = S2.mc_cell(S2.TH_NOM * cap, S2.times_for(sn, foot), ks, ld,
                       "centre")
        print(f"    {foot:11s} {100 * cap:3.0f}% {ld:9s} {100 * sd_t:5.0f}% "
              f"{100 * drift:5.0f}% {sn:13s} {c[0]:7.3f} {r['rmse']:6.3f} "
              f"[{r['rmse_lo']:5.3f},{r['rmse_hi']:6.3f}] {r['robust']:7.3f} "
              f"{r['collapse']:5d} {S2.verdict(r):>8s}")


def main():
    print("=" * 100)
    print("  S2b. THE k_e MEASUREMENT ITSELF")
    print("=" * 100)
    out_a = part_a()
    part_a2()
    part_b(out_a)
    out_c = part_c()
    part_d(out_c)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
