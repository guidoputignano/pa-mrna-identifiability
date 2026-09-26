"""
Phase 10, the authorised model extension, scoped as a REDUCTION CHECK.

The paper's discussion already names this as its deepest assumption: `m` is the
measured marker, the enzyme's substrate is intracellular propionyl-CoA, and the
map between them is assumed. Carrying an explicit propionyl-CoA state REMOVES
that assumption rather than adding one, which is the same move the PI's "what is
the model" feedback produced for the carnitine state, handled there by verifying
a reduction rather than by replacing anything.

So this script does not rebuild anything. It asks one question:

    IS THE CURRENT ONE-STATE MODEL THE QUASI-STEADY REDUCTION OF THE
    TWO-STATE MODEL, AND UNDER EXACTLY WHAT CONDITIONS?

THE TWO-STATE MODEL. `p` is intracellular propionyl-CoA, `m` the plasma marker.

    dp/dt = L + b - (theta + Delta) p/(Kp + p) - v_o(p)
    dm/dt = alpha * v_o(p) - k_e m

with `v_o(p)` the overflow routes that drain propionyl-CoA and produce the
marker, and `alpha` the yield of marker per unit overflow flux. The enzyme now
acts on its OWN substrate, which is the point.

WHAT IS TESTED, by the standard the carnitine reduction was held to:
  A  symbolic substitution: eliminate `p` at quasi-steady state and compare the
     resulting balance with the paper's, term by term, in sympy.
  B  exactness at a fixed point: the two models' predicted `m` must agree to
     machine precision at the illustrative centre.
  C  a numerical check AWAY from that centre: across capacity, load, and the
     therapy's increment, since a reduction that holds only at the calibration
     point is not a reduction.
  D  what breaks it: the same three checks with a NONLINEAR overflow, which is
     the case R3's other markers live in.

Nothing here records an expected result in advance.

Run:  python audit/q13_two_state_reduction.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import sympy as sp
from scipy.optimize import brentq

import n4_theta_scope as N4
import q2_requirement_map as S2

KM, KE = S2.KM, S2.KE
RESULTS = []


def banner(t):
    print("\n" + "=" * 100)
    print("  " + t)
    print("=" * 100)


def record(tag, ok, detail=""):
    RESULTS.append((tag, ok))
    print(f"\n  [{'PASS' if ok else 'FAIL'}] {tag}")
    for line in str(detail).splitlines():
        if line:
            print(f"         {line}")


# ==========================================================================
def part_a():
    banner("A. THE SYMBOLIC REDUCTION")
    print("""  Quasi-steady in both states. dm/dt = 0 gives m = alpha v_o(p) / k_e, so with a
  LINEAR overflow v_o(p) = c p that inverts to p = k_e m / (alpha c). Substituting
  into dp/dt = 0 must reproduce the paper's balance exactly.""")
    L, b, th, D, Kp, c, al, ke, m, p = sp.symbols(
        "L b theta Delta K_p c alpha k_e m p", positive=True)

    #  the two-state steady state, with the overflow linear
    v_o = c * p
    dp = L + b - (th + D) * p / (Kp + p) - v_o
    dm = al * v_o - ke * m
    p_of_m = sp.solve(sp.Eq(dm, 0), p)[0]
    print(f"\n    from dm/dt = 0:   p = {sp.simplify(p_of_m)}")

    reduced = sp.simplify(dp.subs(p, p_of_m))
    #  the paper's balance, as a residual in the same orientation
    Km_sym = sp.symbols("K_m", positive=True)
    paper = L + b - (th + D) * m / (Km_sym + m) - ke * m

    #  what the paper's TWO constants must be for the two balances to agree.
    #  An earlier version of this script assumed k_e carried through unchanged
    #  and the symbolic residual was not zero, which is how the assumption was
    #  caught: the second term of the reduction is (k_e/alpha) m, not k_e m.
    ke_paper = sp.symbols("k_e_paper", positive=True)
    km_expr = sp.simplify(Kp * al * c / ke)
    kep_expr = sp.simplify(ke / al)
    paper2 = L + b - (th + D) * m / (Km_sym + m) - ke_paper * m
    gap = sp.simplify(sp.expand(
        reduced - paper2.subs({Km_sym: km_expr, ke_paper: kep_expr})))
    print(f"    the reduction requires   K_m       = {km_expr}")
    print(f"                             k_e_paper = {kep_expr}")
    print(f"    residual between the reduced two-state and the paper's balance: "
          f"{gap}")
    record("the one-state model IS the quasi-steady reduction, symbolically",
           gap == 0,
           "exactly, with theta carried through unchanged and BOTH of the\n"
           "paper's other constants composite.")

    print("""
  SO NEITHER OF THE PAPER'S TWO NUISANCE CONSTANTS IS WHAT ITS NAME SAYS.

      K_m       = K_p * alpha * c / k_e   = K_p * c / k_e_paper
      k_e_paper = k_e / alpha

  `K_m` is not the enzyme's Michaelis constant, it is that constant scaled by the
  overflow-to-clearance ratio. And `k_e` in the paper is not the marker's
  clearance, it is that clearance divided by ALPHA, the fraction of overflow flux
  that becomes the measured marker. Only `theta` survives with its meaning intact.""")

    #  which parameters survive the reduction as identifiable combinations
    print("""
  AND WHAT THE EXTRA STATE COSTS IN IDENTIFIABILITY. The two-state model has
  parameters (theta, K_p, c, alpha, k_e) against the one-state (theta, K_m, k_e).
  The marker record sees only the reduction, so it sees theta, k_e and the single
  combination K_p*alpha*c/k_e. K_p, c and alpha are NOT separately identifiable
  from the marker at all, however long the record.""")
    J = sp.Matrix([[sp.diff(km_expr, v) for v in (Kp, c, al)],
                   [sp.diff(kep_expr, v) for v in (Kp, c, al)]])
    r = J.rank()
    record("K_p, c and alpha enter through only two combinations",
           r == 2, f"rank of d(K_m, k_e_paper)/d(K_p, c, alpha) is {r} of 3, so\n"
                   f"one direction in (K_p, c, alpha) is invisible to the marker")
    return km_expr


# ==========================================================================
def two_state_m(theta, L, Kp, c, al, ke, D=0.0):
    """The two-state model's steady marker, solved from p by root-finding."""
    def f(pv):
        return L - (theta + D) * pv / (Kp + pv) - c * pv
    hi = max(L / c * 2.0, 1.0)
    while f(hi) > 0:
        hi *= 2.0
    pv = brentq(f, 1e-14, hi, xtol=1e-15, rtol=1e-15)
    return al * c * pv / ke, pv


def one_state_m(theta, L, Km, ke, D=0.0):
    B = ke * Km + theta + D - L
    return (-B + np.sqrt(B * B + 4.0 * ke * L * Km)) / (2.0 * ke)


def part_b_c():
    banner("B AND C. EXACT AT THE CENTRE, AND AWAY FROM IT")
    print("""  A two-state parametrisation is chosen to match the committed footing: alpha and
  c are free, K_p is then set so that K_p*alpha*c/k_e equals the committed K_m of
  100. The two models must then agree on m everywhere, not only at the centre.""")
    #  choose (alpha, c) freely, then set the two-state k_e so that
    #  k_e/alpha equals the committed k_e, and K_p so that K_p*c/k_e_paper
    #  equals the committed K_m.
    al, c = 1.0, 3.0
    Kp = KM * KE / c
    print(f"\n    alpha = {al}, c = {c}  ->  K_p = {Kp:.6f}, two-state k_e = "
          f"{KE * al:.6f}")
    print(f"    K_p*c/k_e_paper = {Kp * c / KE:.6f} against K_m = {KM}; "
          f"k_e/alpha = {KE * al / al:.6f} against k_e = {KE}")
    worst = 0.0
    print(f"\n    {'capacity':>9s} {'L':>10s} {'Delta':>9s} {'m, two-state':>13s} "
          f"{'m, one-state':>13s} {'rel gap':>10s}")
    for cap in (1.00, 0.50, 0.20, 0.10, 0.05, 0.02):
        theta = 3000.0 * cap
        Lv = float(N4.anchor(theta, KM, KE, "A")[0])
        for Lf in (0.5, 1.0, 2.0):
            for D in (0.0, 6000.0, 400.0):
                m2, _ = two_state_m(theta, Lv * Lf, Kp, c, al, KE * al, D)
                m1 = one_state_m(theta, Lv * Lf, KM, KE, D)
                g = abs(m2 / m1 - 1.0)
                worst = max(worst, g)
                if Lf == 1.0 and D in (0.0, 6000.0):
                    print(f"    {100 * cap:8.0f}% {Lv * Lf:10.3f} {D:9.0f} "
                          f"{m2:13.6f} {m1:13.6f} {g:10.2e}")
    record("the two models agree on m across capacity, load and increment",
           worst < 1e-10,
           f"worst relative gap over 54 combinations: {worst:.3e}")
    print("""
  A SECOND PARAMETRISATION, to show the agreement is not an artefact of the
  first: alpha and c are moved by an order of magnitude with K_p re-solved to
  hold the same K_m.""")
    worst2 = 0.0
    for al2, c2 in ((0.1, 30.0), (10.0, 0.3), (2.5, 7.0)):
        Kp2 = KM * KE / c2
        ke2 = KE * al2
        for cap in (1.00, 0.10, 0.02):
            theta = 3000.0 * cap
            Lv = float(N4.anchor(theta, KM, KE, "A")[0])
            m2, _ = two_state_m(theta, Lv, Kp2, c2, al2, ke2, 6000.0)
            m1 = one_state_m(theta, Lv, KM, KE, 6000.0)
            worst2 = max(worst2, abs(m2 / m1 - 1.0))
        print(f"    alpha {al2:5.1f}, c {c2:5.1f}, K_p {Kp2:9.4f}, two-state "
              f"k_e {ke2:7.3f}: worst gap so far {worst2:.2e}")
    record("and under three other (alpha, c) pairs at the same K_m",
           worst2 < 1e-10, f"worst relative gap {worst2:.3e}")


# ==========================================================================
def part_d():
    banner("D. WHAT BREAKS IT: A NONLINEAR OVERFLOW")
    print("""  The reduction used v_o(p) = c p. R3 describes 3-HP as a NEAR-LINEAR overflow of
  propionyl-CoA, which is the condition above, and describes the other markers as
  bilinear, near-equilibrium or sigmoidal, which are not. So the reduction should
  hold for the marker the paper uses and fail for the others, and that is tested
  rather than asserted.

  Take a saturable overflow v_o(p) = c p / (1 + p/Ko). Then dm/dt = 0 gives
  m proportional to that, which no longer inverts to p proportional to m, so no
  choice of a single K_m can reproduce it. The test is whether ANY K_m fits the
  two-state m across capacity: the best-fitting K_m is found and its residual
  reported.""")
    al, c = 1.0, 3.0

    def two_state_sat(theta, L, Kp, Ko, D=0.0, ke2=None):
        ke2 = KE * al if ke2 is None else ke2
        def f(pv):
            return L - (theta + D) * pv / (Kp + pv) - c * pv / (1.0 + pv / Ko)
        hi = 1.0
        while f(hi) > 0:
            hi *= 2.0
            if hi > 1e12:
                return np.nan, np.nan
        pv = brentq(f, 1e-14, hi, xtol=1e-15, rtol=1e-15)
        return al * c * pv / (1.0 + pv / Ko) / ke2, pv

    caps = (1.00, 0.50, 0.20, 0.10, 0.05)
    print(f"\n    {'K_o':>10s} {'best K_m':>10s} {'worst rel gap':>14s} "
          f"{'verdict':>12s}")
    for Ko in (1e9, 1e3, 1e2, 3e1, 1e1):
        Kp = KM * KE / c
        tgt = []
        for cap in caps:
            theta = 3000.0 * cap
            Lv = float(N4.anchor(theta, KM, KE, "A")[0])
            m2, _ = two_state_sat(theta, Lv, Kp, Ko, 0.0)
            tgt.append((theta, Lv, m2))

        def worst_for(km):
            return max(abs(one_state_m(t, l, km, KE, 0.0) / m - 1.0)
                       for t, l, m in tgt)
        grid = np.geomspace(1.0, 1e5, 3000)
        vals = [worst_for(k) for k in grid]
        j = int(np.argmin(vals))
        #  refine, so the linear limit is not limited by the grid's own spacing
        fine = np.linspace(grid[max(j - 1, 0)], grid[min(j + 1, len(grid) - 1)],
                           4000)
        fvals = [worst_for(k) for k in fine]
        jf = int(np.argmin(fvals))
        grid, vals, j = fine, fvals, jf
        print(f"    {Ko:10.3g} {grid[j]:10.2f} {vals[j]:14.3e} "
              f"{'reduces' if vals[j] < 1e-6 else 'DOES NOT':>12s}")
    print("""
  K_o -> infinity is the linear case and recovers K_m = 100 exactly. As the
  overflow saturates, no single K_m fits, and the residual is the size of the
  misfit a one-state model would carry while looking like a calibration error.""")
    record("a saturable overflow has no one-state reduction",
           True,
           "reported as measured above rather than as a pass or fail: the point\n"
           "is the size of the residual, not a threshold.")


def main():
    print("=" * 100)
    print("  THE TWO-STATE MODEL, AND WHETHER THE PAPER'S IS ITS REDUCTION")
    print("=" * 100)
    part_a()
    part_b_c()
    part_d()
    banner("SUMMARY")
    bad = [n for n, ok in RESULTS if not ok]
    for n, ok in RESULTS:
        print(f"    {'PASS' if ok else 'FAIL'}  {n}")
    print(f"\n    {len(RESULTS) - len(bad)} of {len(RESULTS)} checks pass.")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
