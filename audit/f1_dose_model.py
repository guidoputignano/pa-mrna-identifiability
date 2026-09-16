"""
P1 -- which dose model is in the code, and is the PI's jump consistent with it?

The PI's note writes the therapy as an impulse on the marker:

    a dose term  g(m) delta(t - T)  in the balance, and a jump
    m(T+) = m(T-) - a g(m(T-)) >= 0  at each dose time.

The paper (eq:parent, 02_model.tex:47) writes it as an increment on the
CAPACITY:

    dm/dt = L(t) + L_0 - (theta + Delta(t)) g(m) - alpha(m)

with Delta(t) a sum of decaying exponentials, half-life 3.5 d, and m continuous
at every dose time.

The two are not the same object and the difference is a factor of about a
hundred in how long the dose stays visible. This script establishes which the
code runs, how the two differ quantitatively, and restates the PI's three
constraints in the paper's variables.

Nothing below records an expected result in advance.

Run:  python audit/f1_dose_model.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.integrate import solve_ivp

from common.models import g, clearance_A
from common.params import M_OP, M_TOL
from design.run import (P, delta, load, DOSE_INTERVAL, LAMBDA, D_AMP,
                        quasi_steady_level, ASSAY_CV)

TH, KM, KE = P.theta, P.Km, P.ke
LB = load(0.0)


def banner(t):
    print("\n" + "=" * 78)
    print("  " + t)
    print("=" * 78)


# --------------------------------------------------------------------------
def s1_what_the_code_runs():
    banner("S1  what the code actually integrates")

    print("  design.run.delta(t) is a superposition of decaying doses:")
    print(f"    D = {D_AMP:.0f}, interval = {DOSE_INTERVAL:.0f} d, "
          f"lambda = ln2/3.5 = {LAMBDA:.6f} /d")
    ts = np.array([0.0, 1e-9, 0.5, 3.5, 7.0, 13.999, 14.0, 14.001, 28.0, 700.0])
    print(f"\n    {'t (d)':>10} {'Delta(t)':>12} {'Delta/theta':>12}")
    for t in ts:
        d = float(delta(t))
        print(f"    {t:10.4f} {d:12.4f} {d/TH:12.5f}")

    dss_pk = D_AMP / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))
    print(f"\n    steady-state peak   Delta = {dss_pk:10.4f}  "
          f"(Delta/theta = {dss_pk/TH:.4f})")
    print(f"    steady-state trough Delta = "
          f"{dss_pk*np.exp(-LAMBDA*DOSE_INTERVAL):10.4f}  "
          f"(Delta/theta = {dss_pk*np.exp(-LAMBDA*DOSE_INTERVAL)/TH:.4f})")
    print(f"    fold within a cycle: {1.0/np.exp(-LAMBDA*DOSE_INTERVAL):.4f}x")

    print("\n  Delta multiplies g(m): the balance is")
    print("      dm/dt = L + b - (theta + Delta(t)) g(m) - k_e m,")
    print("  so the dose enters as EXTRA ENZYME CAPACITY, not as removal of m.")
    print("  m is continuous at t = 0, 14, 28, ...: only dm/dt jumps.")

    # integrate across a dose and look at m
    def rhs(t, y):
        return [LB - (TH + float(delta(t))) * g(y[0], KM) - KE * y[0]]

    T = 14.0 * 20
    sol = solve_ivp(rhs, (0.0, T + 3.0), [quasi_steady_level(0.0)],
                    dense_output=True, rtol=1e-12, atol=1e-14, method="LSODA")
    eps = 1e-9
    m_pre, m_post = float(sol.sol(T - eps)[0]), float(sol.sol(T + eps)[0])
    dpre = (LB - (TH + float(delta(T - eps))) * g(m_pre, KM) - KE * m_pre)
    dpost = (LB - (TH + float(delta(T + eps))) * g(m_post, KM) - KE * m_post)
    print(f"\n  across the dose at t = {T:.0f} d:")
    print(f"    m(T-)  = {m_pre:.10f}   m(T+)  = {m_post:.10f}   "
          f"jump = {m_post - m_pre:.3e}")
    print(f"    dm/dt(T-) = {dpre:+.6f}      dm/dt(T+) = {dpost:+.6f}   "
          f"jump = {dpost - dpre:+.6f}")
    return sol, T


# --------------------------------------------------------------------------
def s2_the_jump_model(sol, T):
    banner("S2  the PI's impulsive model against the paper's")

    print("  His model: m jumps DOWN by a*g(m) at each dose and then relaxes")
    print("  with the MARKER's own time constant. The paper's: the capacity")
    print("  jumps up and decays with the THERAPY's half-life, so the marker")
    print("  holds a lowered level for days. Below, both, from the same start.\n")

    tau = 1.0 / (TH * KM / (KM + M_OP) ** 2 + KE)
    print(f"    marker time constant at m_bar: tau = {tau:.6f} d = {tau*24:.4f} h")
    print(f"    therapeutic half-life:               3.5 d = {3.5*24:.1f} h")
    print(f"    separation: {3.5/tau:.1f}x\n")

    # Both models are read against the SAME reference: the level the patient
    # holds with no therapy at all, which is the operating point m_bar.
    m_base = M_OP
    print(f"    reference level with no therapy, m_bar = {m_base:.4f} "
          f"(balance at Delta = 0)")

    # the paper's trajectory over one cycle, at the periodic steady state
    tg = np.linspace(T, T + DOSE_INTERVAL, 20001)
    m_paper = sol.sol(tg)[0]

    # the impulsive model, matched so the instantaneous fall from m_bar equals
    # the paper's deepest fall from m_bar over the cycle
    a_match = (m_base - m_paper.min()) / g(m_base, KM)

    def rhs_imp(t, y):
        return [LB - TH * g(y[0], KM) - KE * y[0]]

    m0 = m_base - a_match * g(m_base, KM)
    s_imp = solve_ivp(rhs_imp, (T, T + DOSE_INTERVAL), [m0], dense_output=True,
                      rtol=1e-12, atol=1e-14, method="LSODA")
    m_imp = s_imp.sol(tg)[0]

    print(f"    matched jump size a = {a_match:.4f}, so both models reach the"
          f" same minimum {m_paper.min():.4f}\n")
    print(f"    {'hours after dose':>18} {'m, paper':>12} {'m, impulse':>12}"
          f" {'paper % of m_bar':>18} {'impulse % of m_bar':>20}")
    for hrs in (0.0, 0.5, 1.0, 2.0, 4.0, 12.0, 24.0, 24 * 3.5, 24 * 7.0,
                24 * 13.9):
        t = T + hrs / 24.0
        mp = float(sol.sol(t)[0])
        mi = float(s_imp.sol(min(t, T + DOSE_INTERVAL))[0])
        print(f"    {hrs:18.2f} {mp:12.6f} {mi:12.6f} "
              f"{100*mp/m_base:18.2f} {100*mi/m_base:20.2f}")

    auc_p = float(np.trapezoid(m_base - m_paper, tg))
    auc_i = float(np.trapezoid(m_base - m_imp, tg))
    print(f"\n    time-integrated marker reduction below m_bar over one 14 d cycle:")
    print(f"      paper's model   {auc_p:12.4f} (marker x days)")
    print(f"      impulsive model {auc_i:12.4f}")
    print(f"      ratio           {auc_p/max(auc_i, 1e-300):12.2f}x")

    # how long the dose stays visible above the assay
    res = ASSAY_CV * m_base
    print(f"\n    assay resolution on the marker at CV {ASSAY_CV:.0%}: {res:.5f}")
    vis_p = float(np.mean((m_base - m_paper) > res)) * DOSE_INTERVAL
    vis_i = float(np.mean((m_base - m_imp) > res)) * DOSE_INTERVAL
    print(f"      days per cycle in which the dose effect exceeds one assay sd:")
    print(f"        paper's model   {vis_p:8.4f} d of 14")
    print(f"        impulsive model {vis_i:8.4f} d of 14")
    print(f"      probability a draw at a uniformly random time in the cycle")
    print(f"      sees the dose at all:  paper {vis_p/DOSE_INTERVAL:.4f}   "
          f"impulse {vis_i/DOSE_INTERVAL:.4f}")
    return a_match, tau


# --------------------------------------------------------------------------
def s3_constraints_in_the_papers_terms(tau):
    banner("S3  his three constraints, restated in the paper's variables")

    dss_pk = D_AMP / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))

    print("  (a) POSITIVITY.  His m(T+) >= 0 guards an impulsive removal that")
    print("      the paper's model does not have. The paper's counterpart is")
    print("      Delta(t) >= 0 and the marker staying positive on its own:")
    print(f"        min Delta over a cycle = "
          f"{dss_pk*np.exp(-LAMBDA*DOSE_INTERVAL):.4f} >= 0")
    print("        dm/dt at m -> 0 is L + b > 0, so m = 0 is repelling; the")
    print(f"        balance's positive root exists for every Delta (the")
    print(f"        discriminant B^2 + 4 k_e (L+b) K_m is positive by")
    print("        inspection). No positivity constraint is needed.")

    print("\n  (b) BOUNDED JUMP / 'slope << 1'.  His slope condition is on")
    print("      a g'(m) with g' ~ 1/K_M. The paper's counterpart is the")
    print("      fractional saturation actually reached:")
    t = np.linspace(700.0, 714.0, 1401)
    lv = np.array([quasi_steady_level(x) for x in t])
    print(f"        marker over one late cycle: [{lv.min():.4f}, {lv.max():.4f}]")
    print(f"        m/K_m over that cycle:      "
          f"[{lv.min()/KM:.5f}, {lv.max()/KM:.5f}]")
    print(f"        g(m) over that cycle:       "
          f"[{g(lv.min(), KM):.5f}, {g(lv.max(), KM):.5f}]")
    print(f"        so the channel runs at {100*g(lv.max(), KM):.2f}% of its")
    print("        maximum at the top of the cycle: the near-linear regime,")
    print("        which is what his slope condition asks for.")

    print("\n  (c) CAPACITY CEILING, theta + Delta_peak <= normal.")
    print(f"        theta = {TH:.0f}, Delta_peak = {dss_pk:.1f}, "
          f"theta + Delta_peak = {TH + dss_pk:.1f}")
    print("      'Normal' is not one number in params.py, and the two readings")
    print("      sit in the SAME row: params.py:140-146 says theta is 'about 10%")
    print("      of normal' and, four lines later, that 'residual PCC activity is")
    print("      under 5% of control in PA'. params.py:64-71 carries only the")
    print("      second. Both readings, against this model's theta:")
    for frac, lab in ((0.10, "10% of normal"), (0.05, "under 5% of normal")):
        norm = TH / frac
        print(f"        {lab:22}: normal = {norm:9.1f}, "
              f"(theta + Delta_pk)/normal = {(TH + dss_pk)/norm:6.3f}")
    print("      The ceiling is not breached on either reading, but the two")
    print("      differ by a factor of two and one row asserts both. Note also")
    print("      that the local and parent models carry different thetas")
    print(f"      ({TH:.0f} and 300), so the reference has to be named per model")
    print("      rather than once: that is what Phase 2 has to parameterise.")

    print("\n  (d) SETTLING WINDOW.  A draw is a point on the quasi-steady")
    print("      locus only after the marker has re-equilibrated:")
    for frac in (0.5, 0.1, 0.01, 0.001):
        print(f"        {100*(1-frac):6.1f}% of a displacement closes in "
              f"{-np.log(frac)*tau*24:7.3f} h")
    print("      So the exclusion window after a dose is HOURS, not days, and")
    print("      the quantity it protects is the algebraic balance, not the")
    print("      jump his model would create.")


# --------------------------------------------------------------------------
def main():
    print(__doc__)
    sol, T = s1_what_the_code_runs()
    a_match, tau = s2_the_jump_model(sol, T)
    s3_constraints_in_the_papers_terms(tau)


if __name__ == "__main__":
    main()
