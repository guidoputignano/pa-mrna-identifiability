"""
What a 1-13C-propionate breath test measures under THIS model, and what a
two-level ratio therefore determines.

THE QUESTION. The author's algebra assumes the breath reading is
y = c * theta * g(m), with c an unknown subject-specific recovery factor, and
concludes that a ratio at two known marker levels returns K_m,

    K_m = m1 m2 (1 - r) / (r m2 - m1),     r = y1/y2,

because c and theta both cancel. Section sec:why_not_breath of the paper
describes a different observable: the test measures the FRACTION of a labelled
bolus that is oxidised, recovered as 13CO2 in breath. This script derives that
observable from eq:local rather than assuming a form for it, and redoes the
ratio for whatever comes out.

THE DERIVATION, in one line. Label and unlabelled substrate are
indistinguishable, so tracer leaves the pool at the pool's own fractional
disposal rate C(m)/m, and the enzymatic share of each departure is
theta g(m) / C(m). Only the enzymatic route reaches CO2: the alternative
channel is a conjugation-and-excretion flux and carries no carbon to breath.
So for an infinitesimal bolus at a pool held at m,

    F(m) = theta g(m) / (theta g(m) + k_e m) = theta / (theta + k_e (K_m + m)),

and the reading is y = c F(m) with c lumping the label's onward yield to CO2
and bicarbonate recovery. c still cancels in a ratio, which is the property
that matters. But F is NOT proportional to theta g(m), and the combination a
ratio determines is NOT K_m.

WHAT IT IS INSTEAD. F is a Moebius function of m with a single pole, so
r = F(m1)/F(m2) determines exactly one number,

    u = theta / k_e + K_m,

and K_m and theta follow only after the marker record supplies A, B and C:
K_m = A / (B + C u) and theta = C (u - K_m), with A, B, C as in
audit/u1_fibre_bypass.py.

WHAT THIS SCRIPT DOES.
  A  derives F symbolically from the two-channel split and tests the author's
     proportionality assumption against it.
  B  gives both ratio inversions and both amplification factors in closed
     form, each checked by finite difference.
  C  prices the two caveats the author names: a FINITE bolus perturbs the pool
     during the measurement, so F is a residence-weighted average of the
     instantaneous split rather than its value at m_bar. Integrated exactly.
  D  reproduces the author's own figures on the author's own functional, then
     states them on the derived one.
  E  the negative result: the ratio across a load step against the ratio
     across a dose cycle.

FOOTING. As audit/u1_fibre_bypass.py: theta = 3000 nominal, K_m = 100,
k_e = 5/d, recorded load held at L = 322.727 umol/L/d with the operating level
re-solved at each capacity, D = 6000 on a 14 d train with a 3.5 d half-life,
so Delta runs 400 to 6400 at the train's steady state. Breath readings are
taken at the dose cycle's peak and trough unless stated. Breath noise is a
relative CV on each reading, independent between readings, so
sd(ln r) = sqrt(2) CV.

CAVEAT INHERITED, NOT INTRODUCED. This calculation administers the tracer into
the pool the model calls m, which is the plasma MARKER. The paper states as an
assumption that the marker stands in for the enzyme's actual substrate. The
breath test's bolus is propionate, which is on the substrate side of that map,
so every number below carries the observation-model assumption the paper's
discussion names as its first limit, and carries it more heavily than the
marker record does.

Nothing here records an expected result in advance.

Run:  python audit/u2_breath_functional.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

from design.run import P, delta, load, DOSE_INTERVAL, LAMBDA, D_AMP

TH, KM, KE = P.theta, P.Km, P.ke
LB = float(load(0.0))
T = DOSE_INTERVAL

#: the dose train's steady-state extremes, in closed form rather than by
#: sampling the schedule, so the peak and trough are exact.
DEC = float(np.exp(-LAMBDA * T))
DELTA_PEAK = D_AMP / (1.0 - DEC)
DELTA_TROUGH = DELTA_PEAK * DEC


def level(theta, inc, Ltot, Km=KM, ke=KE):
    """Quasi-steady marker, the positive root of the cleared balance."""
    inc = np.asarray(inc, dtype=float)
    Bq = ke * Km + theta + inc - Ltot
    return (-Bq + np.sqrt(Bq ** 2 + 4.0 * ke * Ltot * Km)) / (2.0 * ke)


# ======================================================================
def part_a():
    print("=" * 74)
    print("PART A.  the observable, derived rather than assumed")
    print("=" * 74)
    m, th, Km, ke = sp.symbols("m theta K_m k_e", positive=True)
    g = m / (Km + m)
    enz = th * g
    alt = ke * m
    F = sp.simplify(enz / (enz + alt))
    print("\n  enzymatic flux            :", enz)
    print("  alternative flux          :", alt)
    print("  oxidised fraction F(m)    :", F)
    print("  F, denominator cleared    :", sp.simplify(sp.cancel(F)))

    #  is F proportional to theta*g(m)? The test is whether F / (theta g) is
    #  free of m. If it is not, the author's cancellation argument still holds
    #  for c but the inversion for K_m does not.
    ratio = sp.simplify(F / enz)
    print("\n  F / (theta g(m))          :", ratio)
    dep = sp.simplify(sp.diff(ratio, m))
    print("  d/dm of that              :", sp.simplify(dep))
    prop = sp.simplify(dep) == 0
    print("  F is proportional to theta*g(m):", prop)

    #  what the ratio at two levels determines
    m1, m2, r = sp.symbols("m_1 m_2 r", positive=True)
    u = sp.Symbol("u", positive=True)
    Fr = sp.simplify((F.subs(m, m1)) / (F.subs(m, m2)))
    print("\n  r = F(m1)/F(m2)           :", sp.simplify(Fr))
    Fr_u = sp.simplify(Fr.subs(th, ke * (u - Km)))
    print("  with theta = k_e (u - K_m):", sp.simplify(sp.cancel(Fr_u)))
    print("  -> the ratio sees theta, K_m and k_e only through")
    print("     u = theta/k_e + K_m")

    sol_u = sp.solve(sp.Eq(r, sp.cancel(Fr_u)), u)
    print("\n  inverting for u           :", [sp.simplify(s) for s in sol_u])

    #  and the author's own functional, for comparison
    ya = th * g
    ra = sp.simplify((ya.subs(m, m1)) / (ya.subs(m, m2)))
    sol_k = sp.solve(sp.Eq(r, ra), Km)
    print("\n  on the AUTHOR's functional y = c theta g(m):")
    print("    r                       :", sp.simplify(ra))
    print("    inverting for K_m       :", [sp.simplify(s) for s in sol_k])
    authors = m1 * m2 * (1 - r) / (r * m2 - m1)
    print("    author's stated form    :", authors)
    print("    difference              :",
          sp.simplify(sp.simplify(sol_k[0]) - sp.simplify(authors)))

    #  the derived chain back to theta, using A, B, C from u1
    A, B, C = sp.symbols("A B C", positive=True)
    Km_of_u = sp.simplify(A / (B + C * u))
    th_of_u = sp.simplify(C * (u - Km_of_u))
    print("\n  closing with the record's (A, B, C):")
    print("    K_m   =", Km_of_u)
    print("    theta =", sp.simplify(th_of_u))
    #  check against the fibre of u1: theta(K) = A/K - B - C K
    K = sp.Symbol("K", positive=True)
    chk = sp.simplify(th_of_u.subs(u, (A / K - B - C * K) / C + K) - (A / K - B - C * K))
    print("    substituting the fibre's own u(K) returns theta(K):",
          sp.simplify(chk) == 0)
    return prop


# ======================================================================
def F_inst(m, theta, Km=KM, ke=KE):
    """Instantaneous oxidised fraction at pool level m."""
    return theta / (theta + ke * (Km + m))


def dln_du_closed(u, m1, m2):
    """|dln u / dln r| for the derived functional."""
    return (u + m1) * (u + m2) / (u * (m2 - m1))


def dln_dK_closed(Km, m1, m2):
    """|dln K_m / dln r| for the author's functional."""
    return (Km + m1) * (Km + m2) / (Km * (m2 - m1))


def part_b():
    print()
    print("=" * 74)
    print("PART B.  both amplifications in closed form, checked numerically")
    print("=" * 74)
    print("\n  |dln u  / dln r| = (u + m1)(u + m2) / (u (m2 - m1))")
    print("  |dln K_m/ dln r| = (K + m1)(K + m2) / (K (m2 - m1))")
    print("  the same shape with u in place of K_m, and u = theta/k_e + K_m")
    print("  is much LARGER than K_m at high capacity, which is the whole")
    print("  difference between the two readings.\n")
    print(f"  {'theta':>7} {'u':>9} {'m1':>8} {'m2':>8} {'|dlnu/dlnr|':>12}"
          f" {'fin.diff':>10} {'|dlnK/dlnr|':>12} {'fin.diff':>10}")
    print("  " + "-" * 82)
    for theta in (3000.0, 600.0, 300.0, 150.0):
        m1 = float(level(theta, DELTA_PEAK, LB))
        m2 = float(level(theta, DELTA_TROUGH, LB))
        u = theta / KE + KM
        a_u = dln_du_closed(u, m1, m2)
        a_K = dln_dK_closed(KM, m1, m2)

        h = 1e-6
        def r_of_u(uu):
            return (uu + m2) / (uu + m1)
        fd_u = ((np.log(r_of_u(u * np.exp(h))) - np.log(r_of_u(u * np.exp(-h))))
                / (2 * h))
        def r_of_K(KK):
            return (m1 / (KK + m1)) / (m2 / (KK + m2))
        fd_K = ((np.log(r_of_K(KM * np.exp(h))) - np.log(r_of_K(KM * np.exp(-h))))
                / (2 * h))
        print(f"  {theta:7.0f} {u:9.1f} {m1:8.3f} {m2:8.3f} {a_u:12.3f}"
              f" {1/abs(fd_u):10.3f} {a_K:12.3f} {1/abs(fd_K):10.3f}")
    print("\n  the finite-difference columns are 1/|dln r / dln .|, so they")
    print("  agree with the closed forms when the closed forms are right.")


# ======================================================================
def F_bolus(mbar, theta, d, Km=KM, ke=KE, Ltot=None, tmax_days=2.0):
    """
    The oxidised fraction of a FINITE bolus, integrated exactly.

    The pool is displaced to m_bar + d at t = 0 and relaxes under the same
    balance that set m_bar. Tracer leaves at the pool's fractional disposal
    rate C(m)/m and the enzymatic share of each departure is theta g(m)/C(m),
    so the cumulative oxidised fraction is

        F = (1/d) INT (theta g(m(t)) / m(t)) T(t) dt.

    The background load is held at the value that makes m_bar a steady state,
    which is the quasi-static reading: the marker settles in about 0.8 h
    against a 14 d dose cycle, a separation this repository establishes
    elsewhere (design.run.quasi_steady_check), so the dose train is frozen
    over the tracer's residence.
    """
    if Ltot is None:
        Ltot = theta * mbar / (Km + mbar) + ke * mbar

    def rhs(t, y):
        m, Tr, ox = y
        m = max(m, 1e-12)
        enz = theta * m / (Km + m)
        tot = enz + ke * m
        return [Ltot - tot, -(tot / m) * Tr, (enz / m) * Tr]

    sol = solve_ivp(rhs, (0.0, tmax_days), [mbar + d, d, 0.0],
                    rtol=1e-10, atol=1e-14, method="Radau", dense_output=True)
    return float(sol.y[2, -1] / d), float(sol.y[1, -1] / d)


def part_c():
    print()
    print("=" * 74)
    print("PART C.  the finite bolus: how far F departs from the instantaneous")
    print("=" * 74)
    print("\n  bolus size as a fraction of the operating level. The split")
    print("  fraction falls with m, so a finite bolus biases F DOWN.")
    print("  What matters for the ratio is not the bias but whether it is the")
    print("  SAME at the two levels, because a common factor cancels with c.\n")
    for theta in (3000.0, 300.0):
        m1 = float(level(theta, DELTA_PEAK, LB))
        m2 = float(level(theta, DELTA_TROUGH, LB))
        print(f"  theta = {theta:.0f}   peak level {m1:.3f}, trough {m2:.3f}")
        print(f"    {'d/m_bar':>9} {'F(peak)':>10} {'F(trough)':>10} "
              f"{'r':>9} {'r_inst':>9} {'u_hat':>10} {'u true':>9} "
              f"{'rel err':>9}")
        print("    " + "-" * 80)
        u_true = theta / KE + KM
        for frac in (0.0, 0.01, 0.05, 0.2, 0.5, 1.0, 2.0):
            if frac == 0.0:
                F1, F2 = F_inst(m1, theta), F_inst(m2, theta)
                left = 0.0
            else:
                F1, left1 = F_bolus(m1, theta, frac * m1)
                F2, left2 = F_bolus(m2, theta, frac * m2)
                left = max(left1, left2)
            r = F1 / F2
            r_i = F_inst(m1, theta) / F_inst(m2, theta)
            #  invert r on the instantaneous model, which is what an analyst
            #  who does not model the bolus would do
            u_hat = (m2 - r * m1) / (r - 1.0)
            print(f"    {frac:9.2f} {F1:10.5f} {F2:10.5f} {r:9.5f} "
                  f"{r_i:9.5f} {u_hat:10.2f} {u_true:9.2f} "
                  f"{(u_hat - u_true) / u_true:9.2%}")
        print(f"    (largest tracer fraction still in the pool at 2 d: {left:.2e})")
        print()
    print("  READING. The bolus is administered at the SAME size at both")
    print("  levels, so the perturbation is not a common factor: it displaces")
    print("  the two readings by different relative amounts and the bias in")
    print("  u does NOT cancel. The size of that bias is the number above.")


# ======================================================================
def part_d():
    print()
    print("=" * 74)
    print("PART D.  the author's figures, then the same on the derived form")
    print("=" * 74)
    for cv in (0.05, 0.10):
        sd_lnr = np.sqrt(2.0) * cv
        print(f"\n  breath CV {cv:.0%} per reading, sd(ln r) = {sd_lnr:.4f}, "
              f"peak/trough of one dose cycle")
        print(f"    {'theta':>7} {'m_peak':>8} {'m_trough':>9} {'fold':>6} "
              f"{'sd(lnK) A':>10} {'sd(lnT) A':>10} "
              f"{'sd(ln u) D':>11} {'sd(lnT) D':>10}")
        print("    " + "-" * 80)
        for theta in (3000.0, 600.0, 300.0, 150.0):
            m1 = float(level(theta, DELTA_PEAK, LB))
            m2 = float(level(theta, DELTA_TROUGH, LB))
            fold = m2 / m1
            #  author's reading
            sdK = dln_dK_closed(KM, m1, m2) * sd_lnr
            slope_K = (LB + KE * KM) / theta          # |dln theta / dln K_m|
            sdT_a = slope_K * sdK
            #  derived reading
            u = theta / KE + KM
            sdu = dln_du_closed(u, m1, m2) * sd_lnr
            A = LB * KM
            slope_u = (u / theta) * KE * (A + KE * KM ** 2) / A
            sdT_d = slope_u * sdu
            print(f"    {theta:7.0f} {m1:8.3f} {m2:9.3f} {fold:6.2f} "
                  f"{sdK:10.3f} {sdT_a:10.3f} {sdu:11.3f} {sdT_d:10.3f}")
    print("\n  columns marked A are on the author's assumed y = c theta g(m);")
    print("  columns marked D are on the functional derived in part A.")
    print("  Both propagate through the fibre with A, B and C held EXACT,")
    print("  which is the author's own stated optimism and is corrected in")
    print("  audit/u3_joint_fisher.py by doing the joint Fisher instead.")


# ======================================================================
def part_e():
    print()
    print("=" * 74)
    print("PART E.  the negative result: load step against dose cycle")
    print("=" * 74)
    print("\n  A load step compares two readings at the SAME dose-cycle phase")
    print("  and two load levels; the dose cycle compares two phases at the")
    print("  same load. The phase is held at the train's trough for the load")
    print("  step, so the dose is not doing the work in either arm.\n")
    for cv in (0.05, 0.10):
        sd_lnr = np.sqrt(2.0) * cv
        print(f"  breath CV {cv:.0%}, sd(ln r) = {sd_lnr:.4f}")
        print(f"    {'theta':>7} {'arm':>12} {'m1':>8} {'m2':>8} {'fold':>6} "
              f"{'sd(lnK) A':>10} {'sd(ln u) D':>11} {'sd(lnT) D':>10}")
        print("    " + "-" * 78)
        for theta in (3000.0, 300.0):
            u = theta / KE + KM
            A = LB * KM
            slope_u = (u / theta) * KE * (A + KE * KM ** 2) / A
            arms = {
                "dose cycle": (float(level(theta, DELTA_PEAK, LB)),
                               float(level(theta, DELTA_TROUGH, LB))),
                "load -50%": (float(level(theta, DELTA_TROUGH, 0.5 * LB)),
                              float(level(theta, DELTA_TROUGH, LB))),
                "load -25%": (float(level(theta, DELTA_TROUGH, 0.75 * LB)),
                              float(level(theta, DELTA_TROUGH, LB))),
            }
            for name, (m1, m2) in arms.items():
                fold = m2 / m1
                sdK = dln_dK_closed(KM, m1, m2) * sd_lnr
                sdu = dln_du_closed(u, m1, m2) * sd_lnr
                print(f"    {theta:7.0f} {name:>12} {m1:8.3f} {m2:8.3f} "
                      f"{fold:6.2f} {sdK:10.3f} {sdu:11.3f} "
                      f"{slope_u * sdu:10.3f}")
        print()
    print("  The comparison is a lever-arm comparison and nothing else: both")
    print("  amplifications go as 1/(m2 - m1), so whichever arm moves the")
    print("  marker further wins, on either functional.")


if __name__ == "__main__":
    prop = part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    print()
    print("=" * 74)
    print("PART A verdict: F is proportional to theta*g(m):", prop)
    print("=" * 74)
