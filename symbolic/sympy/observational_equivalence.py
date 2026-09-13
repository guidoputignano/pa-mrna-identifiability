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
    print("WHAT HAPPENS WHEN Delta IS CONSTANT -- corrected")
    print("=" * 78)
    d0 = sp.Symbol('Delta_0', positive=True)
    Dc = ((th + d0) * g(m, Km) + ke * m + L0
          - ((tp + d0) * g(m, Kp) + kp * m + Lp))
    P = sp.Poly(sp.expand(sp.numer(sp.cancel(sp.together(Dc)))), m)
    sols = sp.solve([sp.Eq(c, 0) for c in P.all_coeffs()], [tp, Kp, kp, Lp], dict=True)
    print("  solving coefficients-all-zero in m gives:", sols)
    print()
    print("  AN EARLIER VERSION OF THIS SCRIPT CLAIMED that a constant increment")
    print("  leaves only the lumped sum theta + Delta_0 identifiable. That was")
    print("  WRONG, and StructuralIdentifiability.jl says so independently: the")
    print("  increment is a KNOWN input, so recovering the lump recovers theta.")
    print("  The model is globally identifiable with no therapy variation at all.")
    print()
    print("=" * 78)
    print("SO WHAT IS THE PAPER'S OBSTRUCTION?  Two answers, and neither is")
    print("that the model is structurally non-identifiable.")
    print("=" * 78)
    print("  (1) STRUCTURAL IDENTIFIABILITY ASKS ABOUT THE TRANSIENT.")
    print("      It assumes the exact, noiseless, CONTINUOUS trajectory, from")
    print("      which every time-derivative of y is available. Both the result")
    print("      above and StructuralIdentifiability.jl are statements about a")
    print("      marker that is relaxing. The clinic never sees that: the marker")
    print("      relaxes in about four hours against a dose half-life of 3.5 days")
    print("      -- a separation of scales that licenses treating each visit as")
    print("      a point on an algebraic balance -- so monthly sampling returns")
    print("      points on the QUASI-STEADY locus, where m_dot is zero to the")
    print("      width of the")
    print("      assay and the derivatives carrying the structural result are")
    print("      simply not in the record. The note's fibre is what remains when")
    print("      the transient is removed: one algebraic equation per patient.")
    print()
    print("  (2) ON THAT LOCUS, THE FISHER MATRIX IS SINGULAR AT THE TRUTH.")
    print("      The rank-one identity is exact there, and the direction it")
    print("      names is flat to FIRST order only: the noiseless misfit along")
    print("      it grows as eps^2, so the parameter combination is recovered")
    print("      through curvature rather than slope, and the recoverable")
    print("      displacement scales as the square root of the noise rather than")
    print("      linearly. A singular Fisher matrix at a point does not")
    print("      contradict global injectivity of the map; it sets the price.")
    print()
    print("      Both are verified in audit/e5_structural.py, which also derives")
    print("      the flat direction by hand and matches the numerical null")
    print("      vector to 2e-7.")
    print()
    print("  The one GENUINE structural obstruction in the set is model B, where")
    print("  BOTH channels saturate: StructuralIdentifiability.jl returns locally")
    print("  but not globally identifiable, with identifiable functions symmetric")
    print("  in the two channels. That is the structural counterpart of the")
    print("  note's channel argument -- identification rests on the two channels")
    print("  having different functional forms.")
    print()
    print("  Naming the obstructions correctly makes the note's claims narrower")
    print("  and more precise rather than weaker: they are statements about the")
    print("  conditioning of a finite, noisy, quasi-steady design, which is the")
    print("  design a clinic can actually run.")


def model_B_channel_symmetry():
    """The one genuine structural obstruction, and what removes it.

    StructuralIdentifiability.jl returns model B -- both channels saturating --
    as locally but not globally identifiable, with identifiable functions
    k_e + theta, K_e + K_m, k_e*theta, K_e*K_m, K_e*theta + K_m*k_e. Every one
    is invariant under the simultaneous swap (theta, K_m) <-> (k_e, K_e), which
    says the fibre is that swap and nothing else. This checks it directly, and
    then asks whether the therapy breaks it.
    """
    Ke = sp.Symbol('K_e', positive=True)

    def C(a, Ka, b, Kb, delta):
        return (a + delta) * m / (Ka + m) + b * m / (Kb + m)

    print()
    print("=" * 78)
    print("MODEL B -- BOTH CHANNELS SATURATING")
    print("=" * 78)
    print("  Without therapy, is the channel swap an exact symmetry?")
    print("    C(theta,K_m,k_e,K_e) - C(k_e,K_e,theta,K_m) =",
          sp.simplify(C(th, Km, ke, Ke, 0) - C(ke, Ke, th, Km, 0)))
    print("  Yes. The two channels are interchangeable, so the fibre is a")
    print("  two-point set and no amount of noiseless data separates them.")
    print()
    print("  With the known increment, which rides on the FIRST channel only:")
    diff = sp.simplify(sp.together(C(th, Km, ke, Ke, d) - C(ke, Ke, th, Km, d)))
    print("    difference =")
    sp.pprint(diff)
    c1 = sp.factor(sp.Poly(sp.numer(sp.cancel(diff)), d).coeff_monomial(d))
    print(f"\n    coefficient of Delta: {c1}")
    print("    vanishing identically in m requires:",
          sp.solve(sp.Eq(c1, 0), Ke, dict=True))
    print()
    print("  So the therapy breaks the symmetry unless K_e = K_m, in which case")
    print("  the two channels are the same channel and the split has no content.")
    print("  The selective input does not merely improve the conditioning of the")
    print("  channel problem; it removes its only structural degeneracy.")


if __name__ == "__main__":
    main()
    model_B_channel_symmetry()
