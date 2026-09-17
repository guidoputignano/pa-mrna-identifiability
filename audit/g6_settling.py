"""
Phase 3, reconciling three settling times that are all correct.

Three numbers are in circulation for "how long after a dose is a draw no good":

    4.2  h   06_excitation.tex:58, from e5_structural.py:57-61
    3.71 h   PHASE4_PI.md P1f, from f1
    1.30 h   Phase 3 (g1, g3, g_design_common.settling_days)

They differ by more than a factor of three and none of them is wrong. They are
answers to three different questions, and the differences are exactly the two
choices buried in the phrase "settling time": WHICH operating point the rate
is evaluated at, and WHAT displacement is being closed.

    (A) the LINEARISED time constant, and its 99% multiple, at a stated m
    (B) the NONLINEAR time to close 99% of a large gap to a FIXED target
    (C) the time for the ACTUAL post-dose trajectory to come within a stated
        relative tolerance of the MOVING quasi-steady level

This script computes all three from the repository's own balance and shows
which choice accounts for which gap. It then states which one Phase 3 uses and
why, so that one definition is used consistently from here.

Nothing here records an expected result in advance.

Run:  python audit/g6_settling.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from scipy.integrate import solve_ivp

from design.run import P, delta, load, M_BAR, DOSE_INTERVAL
import g_design_common as G

M_STAR = 20.0          # e5_structural.py:41's operating point
LEVELS = [(0.50, "50%"), (0.90, "90%"), (0.99, "99%"), (0.999, "99.9%")]


def head(t):
    print("\n" + "=" * 92); print("  " + t); print("=" * 92)


def g(m):
    return m / (P.Km + m)


def tau_at(m, inc=0.0):
    """The linearised time constant of the balance at level m, in days."""
    return 1.0 / ((P.theta + inc) * P.Km / (P.Km + m) ** 2 + P.ke)


def qs(t, Lb):
    """The quasi-steady level at time t, the positive root of the balance."""
    inc = np.asarray(delta(np.atleast_1d(t)), float)
    B = P.ke * P.Km + P.theta + inc - Lb
    return (-B + np.sqrt(B ** 2 + 4.0 * P.ke * Lb * P.Km)) / (2.0 * P.ke)


# ------------------------------------------------------------------- (A)
def linearised():
    head("(A) THE LINEARISED TIME CONSTANT, AT TWO OPERATING POINTS")
    print("  tau = 1 / [ theta K_m/(K_m+m)^2 + k_e ], which is design/run.py:223.")
    print("  A displacement decays as exp(-t/tau), so closing a fraction f takes")
    print("  -ln(1-f) tau. Nothing here knows about the dose; it is the rate at")
    print("  which the marker returns to WHATEVER level the balance implies.\n")
    print(f"    {'evaluated at':>26} {'tau (h)':>9} " +
          " ".join(f"{lab:>9}" for _, lab in LEVELS))
    print("  " + "-" * 82)
    for m, lab in ((M_BAR, f"m_bar = {M_BAR:g} (the patient)"),
                   (M_STAR, f"m* = {M_STAR:g} (e5's point)")):
        tau = tau_at(m)
        print(f"    {lab:>26} {tau*24:9.3f} " +
              " ".join(f"{-np.log(1-f)*tau*24:9.3f}" for f, _ in LEVELS))
    print("\n  The paper's sentence pairs `about 0.8 h` with `99% in 4.2 h`. The")
    print("  first is tau at m_bar; the second is NOT its 99% multiple, which is")
    print(f"  {-np.log(0.01)*tau_at(M_BAR)*24:.2f} h. P1f's 3.71 h is that multiple. The 4.2 h comes")
    print("  from (B) below, at the other operating point. Two operating points")
    print("  in one sentence, which is the whole of that discrepancy.")


# ------------------------------------------------------------------- (B)
def nonlinear_to_fixed_target():
    head("(B) THE NONLINEAR TIME TO CLOSE A LARGE GAP TO A FIXED TARGET")
    print("  e5_structural.py holds the load constant at the value that puts the")
    print("  steady state at m* = 20, starts the marker away from it, and times")
    print("  the closure of 99% of the gap. The target does not move and the")
    print("  displacement is large, so the rate is the LOCAL one along the whole")
    print("  path, not the one at either endpoint.\n")
    L = P.theta * g(M_STAR) + P.ke * M_STAR
    print(f"    constant load L = {L:.1f} holds m* = {M_STAR:g}\n")
    print(f"    {'start m0':>10} {'gap':>8} " +
          " ".join(f"{lab:>9}" for _, lab in LEVELS))
    print("  " + "-" * 62)
    for m0 in (5.0, 20.0, 60.0, 10.0):
        if m0 == M_STAR:
            continue
        sol = solve_ivp(lambda t, y: [L - P.theta * g(y[0]) - P.ke * y[0]],
                        (0.0, 5.0), [m0], dense_output=True,
                        rtol=1e-11, atol=1e-14, method="Radau")
        tt = np.linspace(0.0, 5.0, 200001)
        gap = np.abs(sol.sol(tt)[0] - M_STAR)
        row = []
        for f, _ in LEVELS:
            i = np.argmax(gap < (1 - f) * gap[0])
            row.append(f"{tt[i]*24:9.3f}")
        print(f"    {m0:10g} {gap[0]:8.2f} " + " ".join(row))
    print("\n  The 4.2 h of 06_excitation.tex:58 is the 99% column of this table.")
    print("  It is a fair answer to `how long does a large displacement take`,")
    print("  and it is slower than (A) at m_bar because the return rate falls as")
    print(f"  m rises: tau is {tau_at(M_BAR)*24:.3f} h at m_bar and {tau_at(M_STAR)*24:.3f} h at m*.")


# ------------------------------------------------------------------- (C)
def post_dose_departure():
    head("(C) THE ACTUAL POST-DOSE DEPARTURE FROM THE MOVING LOCUS")
    print("  This is the design question, and it is not either of the above. A")
    print("  draw is usable when the marker is ON the quasi-steady locus the fit")
    print("  assumes, to within a stated tolerance. Two things make that happen")
    print("  sooner than (A) or (B) suggest:\n")
    Lb = load(0.0)

    def rhs(t, y):
        inc = float(np.atleast_1d(delta(np.array([t])))[0])
        return [Lb - (P.theta + inc) * y[0] / (P.Km + y[0]) - P.ke * y[0]]

    y = float(np.atleast_1d(qs(0.0, Lb))[0])
    for c in range(5):
        y = float(solve_ivp(rhs, (c * DOSE_INTERVAL + 1e-9,
                                  (c + 1) * DOSE_INTERVAL - 1e-9), [y],
                            rtol=1e-10, atol=1e-13, method="Radau").y[0, -1])
    t0 = 5 * DOSE_INTERVAL
    pre = y
    post = float(np.atleast_1d(qs(t0, Lb))[0])
    print(f"    (i) the displacement is SMALL. The pre-dose level is {pre:.4f} and")
    print(f"        the level the new increment implies is {post:.4f}, a gap of")
    print(f"        {abs(pre-post):.4f}, against the {abs(10.0-M_STAR):.0f} and {abs(60.0-M_STAR):.0f} of (B).")
    tau = tau_at(post, inc=float(np.atleast_1d(delta(np.array([t0])))[0]))
    print(f"   (ii) the rate is FASTER there, because the dose raises the")
    print(f"        capacity: tau at the post-dose level with the increment on")
    print(f"        is {tau*24:.3f} h, against {tau_at(M_BAR)*24:.3f} h at m_bar with no dose.")
    print("  (iii) the TARGET MOVES TOWARDS THE TRAJECTORY. The increment decays")
    print("        from the instant it lands, so the locus falls back to meet the")
    print("        marker rather than waiting for it.\n")
    sol = solve_ivp(rhs, (t0, t0 + DOSE_INTERVAL - 1e-3), [y],
                    dense_output=True, rtol=1e-11, atol=1e-14, method="Radau")
    tt = t0 + np.geomspace(1e-5, DOSE_INTERVAL - 1e-3, 20000)
    rel = np.abs(sol.sol(tt)[0] - qs(tt, Lb)) / qs(tt, Lb)
    print(f"    {'tolerance':>12} {'window (h)':>11} {'worst after':>12} "
          f"{'repo (h)':>10}")
    print("  " + "-" * 52)
    for tol in (0.001, 0.01, 0.05, 0.10):
        bad = np.where(rel >= tol)[0]
        s = (tt[bad[-1]] - t0) if bad.size else 0.0
        after = rel[tt - t0 > s].max()
        ref = G.settling_days(tol=tol)[0] * 24 if tol in (0.01, 0.05) else None
        print(f"    {tol:12.1%} {s*24:11.3f} {100*after:11.3f}% "
              f"{('%.3f' % ref) if ref else '--':>10}")
    print("\n  The 0.1% row is not a settling window. It is the whole cycle,")
    print("  because the quasi-steady approximation itself is only good to")
    print("  about a tenth of a per cent ANYWHERE between doses. That bounds")
    print("  how tight a tolerance it is meaningful to ask for, and 1% sits an")
    print("  order of magnitude inside it.")
    return rel, tt, t0


# ------------------------------------------------------------------- verdict
def verdict():
    head("WHICH ONE PHASE 3 USES, AND WHY")
    print("  Phase 3 uses (C) at a 1% tolerance, 1.30 h, everywhere: g1's shaded")
    print("  wedge, g3's reachable-phase filter and g5's slip counting all call")
    print("  g_design_common.settling_days, so they cannot disagree with each")
    print("  other. The reason is that (C) is the only one of the three that")
    print("  answers the question the exclusion rule is for: a sample is")
    print("  excluded because a STEADY-STATE FIT does not apply to it, and (C)")
    print("  measures exactly the departure from that fit's own locus.")
    print("\n  (A) and (B) answer a different and also useful question, how fast")
    print("  the marker returns after being pushed, and they are what licenses")
    print("  the quasi-steady treatment in the first place by separating hours")
    print("  from the 3.5-day therapeutic half-life. Neither is the exclusion")
    print("  window, and using one as the window overstates it by 3x.")
    print("\n  Nothing in the paper is edited by this. 06_excitation.tex:58 is")
    print("  correct as (B); its `0.8 h` and its `4.2 h` are just evaluated at")
    print("  different operating points, which P1f already recorded. The data")
    print("  specification's `no sample within one hour of a dose` sits between")
    print("  the 5% and 1% rows of (C), which is the right order of magnitude")
    print("  for an instruction a clinic has to follow, and is the paper's own")
    print("  choice to keep or revisit.")


if __name__ == "__main__":
    print(__doc__)
    linearised()
    nonlinear_to_fixed_target()
    post_dose_departure()
    verdict()
