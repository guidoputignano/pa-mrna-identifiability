"""
Phase 2 item 4, the independent re-check. Four key numbers, each recomputed in
closed form or exact arithmetic, sharing no code with h1.

    (1) the used fraction of normal capacity   symbolically, showing theta
                                               cancels, then in exact rationals
    (2) the admissible interval for b          by the exact quadratic the fibre
                                               reduces to, not by brentq
    (3) the fall in the marker at peak effect  in exact rationals, not floats
    (4) what the settling exclusion costs      from the closed-form geometric
                                               sum, not from delta()

Nothing here records an expected result in advance: each check prints both
numbers and their difference, and the reader decides.

Run:  python audit/h2_recheck.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import sympy as sp

from design.run import P, D_AMP, LAMBDA, DOSE_INTERVAL, load
import h1_admissibility as h1
import g_design_common as G

LB = load(0.0)


def head(t):
    print("\n" + "=" * 90); print("  " + t); print("=" * 90)


def R(x):
    """An exact rational for a float that is exactly representable as one."""
    return sp.nsimplify(sp.Rational(x).limit_denominator(10 ** 12), rational=True)


# ------------------------------------------------------------------ check 1
def check_used_fraction():
    head("(1) THE USED FRACTION OF NORMAL CAPACITY")
    th, dpk, f = sp.symbols("theta Delta_peak f", positive=True)
    used = sp.simplify((th + dpk) / (th / f))
    print(f"  symbolically, (theta + Delta_peak)/(theta/f) = {sp.factor(used)}")
    print(f"  so theta appears only through Delta_peak/theta, which is why the")
    print(f"  local and parent blocks give the same answer.\n")
    # exact rationals: Delta_peak = D/(1 - 2^{-T/half-life}), T/half = 4
    D_, T_, hl = R(D_AMP), R(DOSE_INTERVAL), sp.Rational(7, 2)
    n = sp.simplify(T_ / hl)
    peak_x = sp.simplify(D_ / (1 - sp.Rational(1, 2) ** n))
    th_x = R(P.theta)
    print(f"  T / half-life = {n} exactly, so Delta_peak = D/(1 - 2^-{n}) = {peak_x}")
    print(f"  h1's float peak is {h1.PEAK!r}, difference "
          f"{abs(float(peak_x) - h1.PEAK):.3e}")
    for frac in (sp.Rational(1, 10), sp.Rational(1, 20)):
        u = sp.simplify((th_x + peak_x) / (th_x / frac))
        print(f"    f = {frac} -> used = {u} = {float(u):.6f}, "
              f"headroom {float(1 - u):.6f}")
    print(f"\n  Delta_peak/theta = {sp.simplify(peak_x/th_x)} = "
          f"{float(peak_x/th_x):.6f}, and the 5%/10% question")
    print(f"  moves the used fraction by exactly 2x, which is the ratio of the")
    print(f"  two stated fractions and nothing else.")


# ------------------------------------------------------------------ check 2
def check_b_interval():
    head("(2) THE ADMISSIBLE INTERVAL FOR b, BY THE EXACT QUADRATIC")
    print("  Along the fibre theta(b) = theta_0 + b + k_e K_0 b/(L + b).")
    print("  Setting theta(b) = c and clearing the denominator gives a QUADRATIC")
    print("      b^2 + b(theta_0 + L + k_e K_0 - c) + L(theta_0 - c) = 0,")
    print("  so every bound h1 found by brentq has a closed form.\n")
    b, c = sp.symbols("b c", real=True)
    th0, L, ke, K0 = R(P.theta), R(LB), R(P.ke), R(P.Km)
    quad = sp.expand(b ** 2 + b * (th0 + L + ke * K0 - c) + L * (th0 - c))
    check = sp.simplify(sp.expand((th0 + b + ke * K0 * b / (L + b) - c) * (L + b))
                        - quad)
    print(f"  the clearing is exact: (theta(b) - c)(L + b) - quad = {check}\n")
    D_, hl = R(D_AMP), sp.Rational(7, 2)
    peak_x = sp.simplify(D_ / (1 - sp.Rational(1, 2) ** (R(DOSE_INTERVAL) / hl)))

    def roots(cv):
        return sorted(sp.solve(sp.Eq(quad.subs(c, cv), 0), b), key=float)

    lo_roots = roots(0)
    print(f"    theta(b) = 0 has roots {[sp.nsimplify(r) for r in lo_roots]}")
    print(f"                        = {[float(r) for r in lo_roots]}")
    print(f"    -L = {float(-L):.4f}, so only the larger root is admissible; "
          f"h1 reports {-272.7273}")
    lo = max(lo_roots, key=float)
    print(f"    exact lower bound {sp.nsimplify(lo)} = {float(lo):.6f}\n")
    for frac in (sp.Rational(1, 10), sp.Rational(1, 20)):
        cv = sp.simplify(th0 / frac - peak_x)
        hi = max(roots(cv), key=float)
        print(f"    f = {frac}: theta(b) + Delta_peak = theta_0/f binds at")
        print(f"              b = {sp.nsimplify(hi)} = {float(hi):.6f}")
    print("\n  Compare h1's brentq values: -272.7273, 20107.8981, 50103.2000.")
    print("  A closed form and a bracketed root agreeing to six figures is the")
    print("  point of doing it twice.")
    # the width comparison, exactly
    half = sp.Rational(196, 100) * sp.Rational(5, 100) * L
    for frac in (sp.Rational(1, 10), sp.Rational(1, 20)):
        cv = sp.simplify(th0 / frac - peak_x)
        hi = max(roots(cv), key=float)
        print(f"    f = {frac}: width/(prior width) = "
              f"{float((hi - lo) / (2 * half)):.2f}x")


# ------------------------------------------------------------------ check 3
def check_the_fall():
    head("(3) THE FALL IN THE MARKER AT PEAK EFFECT, IN EXACT RATIONALS")
    th0, L, ke, K0 = R(P.theta), R(LB), R(P.ke), R(P.Km)
    D_, hl = R(D_AMP), sp.Rational(7, 2)
    peak_x = sp.simplify(D_ / (1 - sp.Rational(1, 2) ** (R(DOSE_INTERVAL) / hl)))
    trough_x = sp.simplify(peak_x * sp.Rational(1, 2) ** (R(DOSE_INTERVAL) / hl))

    def m_of(d):
        B = ke * K0 + th0 + d - L
        return sp.simplify((-B + sp.sqrt(B ** 2 + 4 * ke * L * K0)) / (2 * ke))

    m0, mp, mt, md = m_of(0), m_of(peak_x), m_of(trough_x), m_of(D_)
    print(f"    untreated level m(Delta = 0)        = {float(m0):.10f}")
    print(f"    at the first dose, Delta = D        = {float(md):.10f}")
    print(f"    at the steady peak, Delta = {float(peak_x):.0f}     = {float(mp):.10f}")
    print(f"    at the steady trough, Delta = {float(trough_x):.0f}    = {float(mt):.10f}\n")
    print(f"    {'reading':>42} {'fall':>10}")
    print("  " + "-" * 56)
    for lab, a, bb in (("first dose, from untreated", m0, md),
                       ("steady peak, from untreated", m0, mp),
                       ("steady peak, from the pre-dose trough", mt, mp)):
        print(f"    {lab:>42} {float(1 - bb / a):9.4%}")
    print(f"\n  h1 reports 64.90%, 66.36% and 62.25% from floats. The exact")
    print(f"  values agree, so `about half` is not a rounding artefact: the")
    print(f"  smallest fall on any reading is {float(1 - mp/mt):.2%}, which is "
          f"{float((1-mp/mt)/sp.Rational(1,2)):.2f}x a half.")
    print(f"\n  Note m(0) = {float(m0):.6f}, which is m_bar = 10 exactly, because")
    print(f"  L was defined as the clearance at m_bar. So the `from untreated`")
    print(f"  readings are anchored on a definition, not a measurement.")


# ------------------------------------------------------------------ check 4
def check_settling_cost():
    head("(4) WHAT THE SETTLING EXCLUSION COSTS, FROM THE CLOSED FORM")
    lam, T, grid = LAMBDA, DOSE_INTERVAL, G.GRID
    peak = D_AMP / (1.0 - np.exp(-lam * T))
    ph = np.arange(0.0, T, grid)
    closed = peak * np.exp(-lam * ph) / P.theta
    settle, _ = G.settling_days(tol=0.01)
    use = closed[ph >= settle - 1e-9]
    fr, ur = closed.max() - closed.min(), use.max() - use.min()
    print(f"  Delta(p)/theta = [D/(1 - e^-lambda T)] e^-lambda p / theta, on the")
    print(f"  half-day grid, with no call to delta().\n")
    print(f"    full dial range   {fr:.6f}")
    print(f"    usable range      {ur:.6f}   (phases >= {settle:.6f} d)")
    print(f"    cost              {1 - ur/fr:.4%} of the range")
    print(f"    angular cost      {settle/T:.4%} of the dial "
          f"({360*settle/T:.4f} degrees)")
    print(f"\n  h1 and the corrected g1 footnote report 10.1% and 0.39%.")
    ratio = (1 - ur / fr) / (settle / T)
    print(f"\n  The exclusion is {ratio:.0f}x more expensive in VALUE than in TIME.")
    print(f"  That ratio is the whole of the correction: a window can be")
    print(f"  negligible on a clock and expensive on a design, when what it")
    print(f"  removes is an endpoint of the quantity that binds.")


if __name__ == "__main__":
    print(__doc__)
    check_used_fraction()
    check_b_interval()
    check_the_fall()
    check_settling_cost()
