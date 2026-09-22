"""
Phase 3, the independent re-check. Five key numbers, each recomputed by a
route that shares no code with the script that produced it.

The point is not to run g1-g3 again. It is to get the same answers from
different machinery, so that a wrong column, a wrong prior or a wrong
convention shows up as a disagreement rather than a consistent mistake.

    (1) the phase-to-weekday lock      by calendar arithmetic, not mod 7
    (2) the dose train                 by the closed-form geometric sum
    (3) the settling window            by integrating the balance directly
    (4) sd(ln theta)                   by the curvature of a profiled
                                       NONLINEAR objective, using no analytic
                                       derivative column at all
    (5) g2's floor and c               by an exact two-point solve, using no
                                       curve fit

Nothing here records an expected result in advance: every check prints the two
numbers and their ratio, and the reader decides.

Run:  python audit/g4_recheck.py
"""
import sys
from pathlib import Path
from datetime import date, timedelta

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

from design.run import (P, delta, load, D_AMP, LAMBDA, DOSE_INTERVAL,
                        LOAD_FRAC, quasi_steady_level)
import g_design_common as G
import g2_visits_vs_sd as g2
import g3_calendar as g3


def head(t):
    print("\n" + "=" * 92); print("  " + t); print("=" * 92)


# ------------------------------------------------------------------- check 1
def check_phase_weekday():
    head("(1) DOES A PHASE REALLY FIX ITS WEEKDAY? -- by calendar arithmetic")
    print("  g3 asserts it from 14 = 2 x 7. Here the claim is checked against")
    print("  datetime, over 60 cycles, for every half-day phase and both dose")
    print("  timings. A single drifting phase would refute it.\n")
    bad = 0
    for off in (g3.MORNING, g3.AFTERNOON):
        for p in np.arange(0.0, DOSE_INTERVAL, G.GRID):
            wds = set()
            for k in range(60):
                t = k * DOSE_INTERVAL + p
                d = g3.DOSE_DAY + timedelta(days=int(np.floor(t + off + 1e-9)))
                wds.add(d.weekday())
            if len(wds) != 1:
                bad += 1
                print(f"    DRIFTS: dose offset {off}, phase {p:g} -> {sorted(wds)}")
    print(f"  phases checked 2 x {int(DOSE_INTERVAL/G.GRID)}, "
          f"phases whose weekday drifts over 60 cycles: {bad}")

    print("\n  And the count of bookable slots, recounted from datetime:")
    settle, _ = G.settling_days(tol=0.01)
    for off, lab in ((g3.MORNING, "morning"), (g3.AFTERNOON, "afternoon")):
        n = 0
        for p in np.arange(0.0, DOSE_INTERVAL, G.GRID):
            d = g3.DOSE_DAY + timedelta(days=int(np.floor(p + off + 1e-9)))
            if d.weekday() < 5 and p >= settle - 1e-9:
                n += 1
        print(f"    {lab:>10} dose: datetime says {n:3d}, "
              f"g3 says {len(g3.reachable(off, settle=settle)):3d}")


# ------------------------------------------------------------------- check 2
def check_dose_train():
    head("(2) THE DOSE TRAIN -- by the closed-form geometric sum")
    print("  design.run.delta sums 200 decaying doses. At steady state that sum")
    print("  is D exp(-lambda p)/(1 - exp(-lambda T)), which is computed here")
    print("  and compared. The peak-to-trough range is what f4 found binds, so")
    print("  it is the one number in Phase 3 that must not be a summation slip.\n")
    t0 = 100 * DOSE_INTERVAL
    ph = np.arange(0.0, DOSE_INTERVAL, 0.5)
    closed = D_AMP * np.exp(-LAMBDA * ph) / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))
    code = np.asarray(delta(t0 + ph), float)
    print(f"    worst relative difference over the cycle: "
          f"{np.max(np.abs(code - closed) / closed):.3e}")
    print(f"    peak   Delta/theta  closed {closed[0]/P.theta:.4f}   "
          f"code {code[0]/P.theta:.4f}")
    print(f"    trough Delta/theta  closed {closed[-1]/P.theta:.4f}   "
          f"code {code[-1]/P.theta:.4f}")
    print(f"    full-dial range     closed "
          f"{(closed.max()-closed.min())/P.theta:.4f}   "
          f"code {(code.max()-code.min())/P.theta:.4f}")
    hi = float(np.atleast_1d(delta(np.array([t0 + G.GRID])))[0])
    lo = float(np.atleast_1d(delta(np.array([t0 + DOSE_INTERVAL - G.GRID])))[0])
    print(f"\n    usable range, i.e. dropping only the on-the-dose slot: "
          f"{(hi - lo)/P.theta:.4f}")
    st = G.settling_days(0.01)[0]
    for off, lab in ((g3.MORNING, "morning"), (g3.AFTERNOON, "afternoon")):
        v = [g3.inc_at_phase(x) for x in g3.reachable(off, settle=st)]
        print(f"    g3's {lab:>9}-dose clinic reaches {max(v) - min(v):.4f}")
    print("    A clinic range equal to the usable range means the weekday")
    print("    constraint costs nothing in the quantity that binds.")


# ------------------------------------------------------------------- check 3
def check_settling():
    head("(3) THE SETTLING WINDOW -- by integrating the balance here")
    print("  G.settling_days calls design.run.quasi_steady_check. This integrates")
    print("  the same balance from scratch, with its own quasi-steady root, and")
    print("  reads off where the relative departure first falls below the")
    print("  tolerance and stays there.\n")
    Lb = load(0.0)

    def qs(t):
        inc = np.asarray(delta(np.atleast_1d(t)), float)
        B = P.ke * P.Km + P.theta + inc - Lb
        return (-B + np.sqrt(B ** 2 + 4.0 * P.ke * Lb * P.Km)) / (2.0 * P.ke)

    def rhs(t, y):
        inc = float(np.atleast_1d(delta(np.array([t])))[0])
        return [Lb - (P.theta + inc) * y[0] / (P.Km + y[0]) - P.ke * y[0]]

    # run five cycles in to forget the initial condition, then watch one
    y = float(np.atleast_1d(qs(0.0))[0])
    for c in range(5):
        y = float(solve_ivp(rhs, (c * DOSE_INTERVAL + 1e-9,
                                  (c + 1) * DOSE_INTERVAL - 1e-9), [y],
                            rtol=1e-10, atol=1e-13, method="Radau").y[0, -1])
    t0 = 5 * DOSE_INTERVAL
    sol = solve_ivp(rhs, (t0, t0 + DOSE_INTERVAL - 1e-3), [y], dense_output=True,
                    rtol=1e-11, atol=1e-14, method="Radau")
    tt = t0 + np.geomspace(1e-5, DOSE_INTERVAL - 1e-3, 20000)
    rel = np.abs(sol.sol(tt)[0] - qs(tt)) / qs(tt)
    for tol in (0.01, 0.05):
        bad = np.where(rel >= tol)[0]
        s = (tt[bad[-1]] - t0) if bad.size else 0.0
        after = rel[tt - t0 > s].max()
        ref = G.settling_days(tol=tol)
        print(f"    tol {tol:4.0%}: here {s*24:6.3f} h (worst after {100*after:6.3f}%)"
              f"   repo {ref[0]*24:6.3f} h (worst after {100*ref[1]:6.3f}%)"
              f"   ratio {s/ref[0] if ref[0] else float('nan'):.4f}")


# ------------------------------------------------------------------- check 4
def check_sd_by_profile():
    head("(4) sd(ln theta) -- by the curvature of a profiled NONLINEAR objective")
    print("  G.sd_ln_theta builds analytic derivative columns and inverts the")
    print("  information matrix. Here the same statistical model is written out")
    print("  as an objective in the parameters themselves, three of them are")
    print("  profiled out numerically at each of several ln theta, and sd is read")
    print("  from the curvature of the profile. No derivative column is used, so")
    print("  a wrong column or a misplaced prior cannot cancel.\n")
    settle, _ = G.settling_days(tol=0.01)
    designs = {
        "g3 clinic design, n=12 (pm dose)":
            np.asarray(G.phase_spread(
                12, gap=30.4,
                phases=g3.spread_over(g3.reachable(g3.AFTERNOON, settle=settle), 12)),
                float),
        "monthly, n=12": G.SCHEDULES["monthly (30.4 d)"](12),
        "monthly, n=25": G.SCHEDULES["monthly (30.4 d)"](25),
        "14.5 d, n=25": G.SCHEDULES["14.5 d"](25),
    }
    sd_b = LOAD_FRAC * G.LB
    for name, times in designs.items():
        m, inc, Ls = G.record(times, 0.0)

        def S(x, lt):
            th, km, b, ke = np.exp(lt), np.exp(x[0]), x[1], np.exp(x[2])
            r = (th + inc) * m / (km + m) + ke * m - Ls - b
            return (np.sum(r ** 2) / G.SIG ** 2 + b ** 2 / sd_b ** 2
                    + (x[2] - np.log(G.KE)) ** 2 / 0.8 ** 2)

        def prof(lt):
            best = np.inf
            for s in (1.0, 0.7, 1.4):
                r = minimize(S, [np.log(G.KM * s), 0.0, np.log(G.KE)], args=(lt,),
                             method="Nelder-Mead",
                             options={"xatol": 1e-12, "fatol": 1e-14,
                                      "maxiter": 200000, "maxfev": 200000})
                best = min(best, float(r.fun))
            return best

        lt0, h = np.log(G.TH), 0.02
        curv = (prof(lt0 + h) - 2 * prof(lt0) + prof(lt0 - h)) / h ** 2
        sd_prof = np.sqrt(2.0 / curv) if curv > 0 else np.inf
        sd_filt = G.sd_ln_theta(m, inc)
        print(f"    {name:>34}  filter {sd_filt:.4f}   profile {sd_prof:.4f}   "
              f"ratio {sd_prof/sd_filt:.4f}")
    print("\n  The objective is 2 x the negative log posterior, so the Fisher")
    print("  information is half its second derivative and sd = sqrt(2/S'').")


# ------------------------------------------------------------------- check 5
def check_floor_c():
    head("(5) g2's floor AND c -- by an exact two-point solve, no curve fit")
    print("  sd(n)^2 = floor^2 + c^2/n is two equations in two unknowns, so any")
    print("  two n determine it exactly. If the form is right the answer barely")
    print("  moves with the pair chosen; if it is wrong the pairs disagree and")
    print("  the least-squares fit was averaging a misfit.\n")
    for sched in ("monthly (30.4 d)", "phase-spread"):
        fn = G.SCHEDULES[sched]
        nstar = g2.saturation_n(fn)
        print(f"  {sched}, k_e prior 0.8, load CV 0%, n* = {nstar}")
        _, md, _ = g2.curve(fn, 0.0, 0.8)
        fl, c, rms, _ = g2.fit_floor_c(g2.NS, md)
        sel = g2.NS >= nstar
        fl2, c2, rms2, _ = g2.fit_floor_c(g2.NS[sel], md[sel])
        print(f"    curve_fit, all n      floor {fl:.4f}  c {c:.4f}  rms {rms:.2e}")
        print(f"    curve_fit, n >= {nstar:<3d}   floor {fl2:.4f}  c {c2:.4f}  rms {rms2:.2e}")
        for n1, n2 in ((20, 60), (25, 52), (30, 44), (36, 60)):
            s1 = G.sd_over_seeds(fn(n1), 0.0)[1]
            s2 = G.sd_over_seeds(fn(n2), 0.0)[1]
            cc2 = (s1 ** 2 - s2 ** 2) / (1.0 / n1 - 1.0 / n2)
            ff2 = s1 ** 2 - cc2 / n1
            print(f"    two-point ({n1:2d},{n2:2d})    floor "
                  f"{np.sqrt(max(ff2,0.0)):.4f}  c {np.sqrt(max(cc2,0.0)):.4f}")
        print()


if __name__ == "__main__":
    print(__doc__)
    check_phase_weekday()
    check_dose_train()
    check_settling()
    check_sd_by_profile()
    check_floor_c()
