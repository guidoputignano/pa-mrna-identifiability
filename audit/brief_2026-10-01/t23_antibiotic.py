"""
Round-3 add-on. Antibiotics as a load term, and as a (downward, safe) excursion.

The disease's propionate load is partly gut-bacterial; antibiotics (metronidazole)
suppress it. The repository KNOWS this in its audit layer but the CORE
identifiability model does not carry antibiotic status:
  - params.py: gut flora >= 22.2 +/- 8.5% of production (fall under metronidazole).
  - m2_antibiotic_cycle.py / Thompson 1990: metabolites fall ~41% (12-76%) on
    metronidazole, a ~21-35 day on/off cycle; "a diet prescription does not record
    antibiotic weeks", so the load record does NOT see it.
  - u15_gut_suppression.py: explores suppression as a k_e relaxation instrument.
Here: (3) the bias an UNMODELLED antibiotic swing puts on theta, and (4) a known
antibiotic STEP (downward, hence safe) as an excursion for theta and k_e.
"""
import numpy as np
from scipy.optimize import least_squares

import core as C
import dyn as D

GUT_SHARE = 0.222          # params.py
THOMPSON_FALL = 0.41       # m2 / Thompson 1990


# --------- (3) unmodelled antibiotic variation: bias on theta ----------
def fit_theta_bfree(loads, inc, m, within_on=None):
    """theta,b free; k_e,K_m FIXED at truth (the identified regime, so the fit is
    stable). Returns (theta_hat, hit_floor)."""
    sb = C.LOAD_FRAC * C.record(C.TH_5, seed=0)["L"]
    sg = np.full(m.size, 0.05 * C.record(C.TH_5, seed=0)["L"])

    def res(x):
        t, b = np.exp(x[0]), x[1]
        return np.r_[(loads - ((t + inc) * m / (C.KM0 + m) + C.KE * m + b)) / sg, [b / sb]]
    s = least_squares(res, [np.log(C.TH_5), 0.0], method="lm", xtol=1e-12, ftol=1e-12)
    th = np.exp(s.x[0])
    return th, th < 0.02 * C.TH_5


def abx_bias(phi, delta, within=True, nmc=80):
    """median ln(theta_hat/theta_true) and fraction of fits driven to theta~0,
    with k_e,K_m fixed (identified regime). within=True: phi of DRAWS suppressed."""
    th = C.TH_5
    out, floor = [], 0
    for s in range(nmc):
        r = C.record(th, seed=s)
        rng = np.random.default_rng(5000 + s)
        if within:
            on = rng.random(r["times"].size) < phi
        else:
            on = (np.ones if rng.random() < phi else np.zeros)(r["times"].size, bool)
        m = C.m_ss(th, C.KM0, C.KE, r["loads"] * np.where(on, 1 - delta, 1.0), r["inc"])
        est, hit = fit_theta_bfree(r["loads"], r["inc"], m)
        out.append(np.log(est / th)); floor += hit
    return float(np.median(out)), floor / nmc


# --------------- (4) a known antibiotic step as an excursion ------------
def step_sds(delta, times, l0frac=0.1):
    """sd(ln theta) [k_e prior] and sd(ln k_e) [no k_e prior] from a KNOWN load
    step of fraction delta (negative = antibiotic suppression, load DOWN)."""
    L = D.LTOT
    ex = D.Experiment([(-1, (1 - l0frac) * L), (0, (1 - l0frac) * L * (1 + delta))],
                      L0=l0frac * L)
    F = D.fim(D.params(ex=ex), ex, times)
    Km_rec, l0p = 0.062, 0.05 * L

    def marg(idx, ke_pr):
        G = F.copy()
        G[1, 1] += (ke_pr ** -2 if ke_pr else 0.0)
        G[2, 2] += Km_rec ** -2
        G[3, 3] += l0p ** -2
        w = np.linalg.eigvalsh(G)
        return np.inf if w[0] <= 1e-12 * w[-1] else float(np.sqrt(np.linalg.inv(G)[idx, idx]))
    return marg(0, 0.8), marg(1, 0.0)


def peak_g(delta, l0frac=0.1):
    L = D.LTOT
    ex = D.Experiment([(-1, (1 - l0frac) * L), (0, (1 - l0frac) * L * (1 + delta))], L0=l0frac * L)
    tt = np.arange(0, 2 + 1e-9, 5 / 60 / 24)
    m, _ = D.simulate(D.params(ex=ex), ex, tt, sens=False)
    mm = m[-1]
    return mm / (C.KM0 + mm)


def main():
    print("=" * 92)
    print("  ADD-ON 3&4. ANTIBIOTICS AS A LOAD TERM, AND AS A DOWNWARD EXCURSION")
    print("=" * 92)

    print("\n  (3a) IS ANTIBIOTIC STATUS IN THE MODEL?")
    print("      CORE identifiability model (core.record / n4, every number in this")
    print("      report): NO. The load is dietary intake, recorded with 10% CV plus a")
    print("      constant bias b; there is no antibiotic term. A diet prescription does")
    print("      not record antibiotic weeks (m2).")
    print(f"      REPO AUDIT layer: YES, characterized but not folded in -- gut share")
    print(f"      {100*GUT_SHARE:.1f}+/-8.5% of production (params.py); metabolites fall "
          f"~{100*THOMPSON_FALL:.0f}% (12-76%)")
    print("      under metronidazole on a 21-35 day cycle (m2, Thompson 1990); u15 treats")
    print("      suppression as a k_e instrument. So L is NOT a patient constant, and the")
    print("      core model's constant-b load is a known idealisation.")

    print("\n  (3b) BIAS ON theta FROM AN UNMODELLED ANTIBIOTIC SWING.")
    print("      Fit with k_e,K_m FIXED at truth (the identified regime -- otherwise the")
    print("      fit is ill-posed and the bias is unbounded, which is itself the point).")
    print("      WITHIN-record cycle (phi of draws suppressed): ln(theta_hat/theta_true)")
    print("      and the fraction of fits driven to theta ~ 0:")
    print(f"      {'depth':>7}" + "".join(f"  phi={p:.0%} (floor)" for p in (0.20, 0.35, 0.50)))
    for delta in (0.30, 0.41, 0.50):
        cells = []
        for p in (0.20, 0.35, 0.50):
            v, fl = abx_bias(p, delta, within=True)
            cells.append(f"{v:+6.2f} ({fl:.0%})")
        print(f"      {delta:>6.0%} " + "  ".join(f"{c:>13}" for c in cells))
    wv, wf = abx_bias(0.5, 0.41, within=False)
    print(f"      WHOLE-record uniform suppression (delta 41%): bias {wv:+.2f}, "
          f"floor {wf:.0%} -- absorbed by the per-record bias b.")
    print("      READING: the bias is DOWNWARD (suppressed-load troughs read as lower")
    print("      capacity) and NOT a modest number. Even a mild 30% cycle pushes ln theta")
    print("      down by ~1-2 and collapses a fifth-to-a-half of fits to theta ~ 0; at the")
    print("      clinically realistic ~41% depth (Thompson) most fits COLLAPSE (69-99%).")
    print("      The cause is that theta is weakly identified, so the near-singular")
    print("      direction amplifies the unmodelled load variation without bound -- a")
    print("      cross-sectional theta estimate is simply not robust to antibiotic status.")
    print("      A UNIFORM whole-record suppression is absorbed by the per-record bias b")
    print("      (so the harm is from VARIATION, as with the load drift in R10). And the")
    print("      longitudinal route (add-on 2) is corrupted too if antibiotic status")
    print("      DIFFERS between the two visits: a load change of delta fakes a capacity")
    print("      change of order (L/theta g)*delta ~ 6*delta, far above the 0.26 contrast")
    print("      precision. So antibiotic status MUST be recorded for either route.")

    print("\n  (4) A KNOWN ANTIBIOTIC STEP AS AN EXCURSION (load DOWN = safe).")
    print("      sd from the transient; daily for two weeks (15 pts) and hourly for two")
    print("      days (49 pts, to catch the ~4 h relaxation). K_m from record, L0 to 5%.")
    daily = np.arange(0, 14 + 1e-9, 1.0)
    hourly = np.arange(0, 2 + 1e-9, 1 / 24)
    print(f"      {'load step':>10}  {'new g':>6} | {'daily-2wk sd(lnth)':>18}  {'sd(lnke)':>9}"
          f" | {'hourly-2d sd(lnth)':>18}  {'sd(lnke)':>9}")
    for delta in (-0.20, -0.30, -0.50):
        dth, dke = step_sds(delta, daily)
        hth, hke = step_sds(delta, hourly)
        print(f"      {delta:>+9.0%}  {peak_g(delta):>6.3f} | {dth:>18.3f}  {dke:>9.3f}"
              f" | {hth:>18.3f}  {hke:>9.3f}")
    print("      UP vs DOWN (the g-direction question), hourly-2d, matched |step| = 30%,")
    print("      and a large +100% up for contrast (g-span = how far g moves):")
    for d in (-0.50, -0.30, +0.30, +1.00):
        th_, ke_ = step_sds(d, hourly)
        L = D.LTOT
        ex = D.Experiment([(-1, 0.9 * L), (0, 0.9 * L * (1 + d))], L0=0.1 * L)
        m, _ = D.simulate(D.params(ex=ex), ex, hourly, sens=False)
        g = m / (C.KM0 + m)
        print(f"        step {d:>+5.0%}: g-span {g.min():.3f}-{g.max():.3f} "
              f"({g.max()-g.min():.3f})  sd(lnth) {th_:.2f}  sd(lnke) {ke_:.3f}")
    print("      READING:")
    print("      - k_e: the antibiotic step is a SAFE, MODERATE k_e instrument, but only")
    print("        if sampled fast enough to see the ~4 h relaxation: hourly gives sd(lnke)")
    print("        ~0.26-0.42 (deeper step = better), daily MISSES it (0.4-0.9, two steady")
    print("        states only). A deliberate larger suppression does better (u15, R15).")
    print("        It is valuable because k_e is the binding nuisance and this adds NO")
    print("        substrate -- it is the safe stage 1 of the two-stage design (R15).")
    print("      - theta: NOT recovered, in EITHER direction (sd(lnth) ~1.5-2.9). The limit")
    print("        is AMPLITUDE, not direction: a downward step is capped at removing the")
    print("        gut share (~20-50% of L), so g moves only ~0.35 -> 0.22-0.28 (span")
    print("        <=0.13). theta's information needs g driven UP to ~0.57 (the 2.3x")
    print("        excursion, g-span ~0.22), which only a large substrate LOAD can do.")
    print("        At matched 30% the two directions barely differ; the low-m collinearity")
    print("        (theta*g ~ (theta/Km)*m, aliased with k_e*m) is a secondary effect under")
    print("        the amplitude cap. So lowering L cannot probe theta -- but it safely")
    print("        probes k_e, which is the half of the problem worth solving.")


if __name__ == "__main__":
    main()
