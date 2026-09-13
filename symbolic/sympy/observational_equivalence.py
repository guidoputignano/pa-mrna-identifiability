"""
3.7 -- structural identifiability of the paper's configuration, symbolically.

Configuration claimed by the paper (eq:local):

    m' = L(t) + L_0 - (theta + Delta(t)) g(m; K_m) - k_e m ,    y = m ,

with Delta KNOWN (the dosing record), L KNOWN (the recorded load), L_0 an
UNKNOWN constant, and parameters p = (theta, K_m, k_e, L_0).

Because y = m is observed directly, m' is known too, so at every sampling
instant the model supplies one exact equation

    F(p; m, Delta) = (theta + Delta) g(m; K_m) + k_e m + L_0 = L - m' .

Two parameter vectors are observationally equivalent exactly when F(p) = F(p')
at every (m, Delta) the record visits. This script asks sympy when that can
happen, and does so WITHOUT reference to any parameter value, which is what
makes the answer structural.
"""
import sympy as sp

th, Km, ke, L0 = sp.symbols('theta K_m k_e L_0', positive=True)
tp, Kp, kp, Lp = sp.symbols("theta' K_m' k_e' L_0'", positive=True)
m, d = sp.symbols('m Delta', positive=True)

g = lambda M, K: M / (K + M)
D = ((th + d) * g(m, Km) + ke * m + L0
     - ((tp + d) * g(m, Kp) + kp * m + Lp))


def main():
    print(__doc__)
    print("=" * 78)
    print("THE OBSERVATIONAL-EQUIVALENCE DIFFERENCE")
    print("=" * 78)
    num = sp.simplify(sp.together(D))
    print("  D(m, Delta) = F(p) - F(p') =")
    sp.pprint(num)

    # ---- clear denominators and read off the coefficient of Delta -------
    poly = sp.Poly(sp.numer(sp.cancel(sp.together(D))), d)
    print()
    print("=" * 78)
    print("STEP 1  --  the coefficient of Delta, which only the THERAPY carries")
    print("=" * 78)
    c1 = sp.factor(sp.simplify(poly.coeff_monomial(d)))
    print("  coefficient of Delta^1 (denominators cleared):")
    sp.pprint(c1)
    sol = sp.solve(sp.Eq(c1, 0), Kp, dict=True)
    print(f"\n  vanishing for all m requires: {sol}")
    print("  => K_m is GLOBALLY IDENTIFIABLE as soon as Delta takes two distinct")
    print("     values at two distinct levels. Nothing else in the model can")
    print("     produce a term proportional to Delta.")

    # ---- with K_m pinned, separate the rest ----------------------------
    print()
    print("=" * 78)
    print("STEP 2  --  with K_m pinned, the remaining three separate")
    print("=" * 78)
    D2 = sp.simplify(D.subs(Kp, Km))
    print("  D with K_m' = K_m:")
    sp.pprint(sp.simplify(sp.expand(D2)))
    print()
    print("  This is a linear combination of the three functions of m")
    print("      g(m) = m/(K_m+m),      m,      1")
    print("  with coefficients (theta-theta'), (k_e-k_e'), (L_0-L_0').")
    W = sp.Matrix([[g(mi, Km), mi, 1] for mi in sp.symbols('m1 m2 m3', positive=True)])
    det = sp.factor(sp.simplify(W.det()))
    print("\n  Their 3x3 Wronskian-style determinant at three distinct levels:")
    sp.pprint(det)
    print("\n  Non-zero whenever m1, m2, m3 are distinct and K_m > 0, so the three")
    print("  functions are linearly independent and all three differences vanish.")

    print()
    print("=" * 78)
    print("RESULT")
    print("=" * 78)
    print("  (theta, K_m, k_e, L_0) are ALL GLOBALLY STRUCTURALLY IDENTIFIABLE,")
    print("  for any parameter values, provided the record satisfies")
    print("      (i)  Delta takes at least two distinct values, and")
    print("      (ii) m  takes at least three distinct values.")

    # ---- the degenerate configurations ---------------------------------
    print()
    print("=" * 78)
    print("THE DEGENERATE CONFIGURATIONS, which are the paper's own findings")
    print("=" * 78)
    print()
    print("(a) Delta CONSTANT -- no therapy, or every sample at the same phase")
    Dc = D.subs(d, sp.Symbol('Delta_0', positive=True))
    print("    The Delta term is then absorbable: theta and Delta_0 enter only")
    print("    as the sum (theta + Delta_0), so theta is identifiable only up to")
    print("    that sum, and the coefficient argument of STEP 1 is unavailable.")
    lump = sp.symbols('Theta', positive=True)
    print("    Identifiable instead: the lumped capacity Theta = theta + Delta_0,")
    print("    together with K_m, k_e and L_0.")
    print()
    print("    This is the symbolic form of the aliasing result: sampling at the")
    print("    dose interval holds Delta nearly constant and costs a factor of 22.")
    print()
    print("(b) m CONSTANT -- a single level, however many increments")
    print("    g(m), m and 1 are then three constants, hence linearly DEPENDENT,")
    print("    so STEP 2 fails: (theta, k_e, L_0) are identifiable only through")
    print("    two combinations. K_m still follows from STEP 1 if Delta varies.")
    print("    This is the symbolic form of the single-level result: two")
    print("    increments at one level give rank two of four.")


if __name__ == "__main__":
    main()
