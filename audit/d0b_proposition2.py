"""
Does Proposition 2 survive the corrected K_m column?

Its proof writes the four columns as
    ( theta g,  -theta K_m m/(K_m+m)^2,  k_e m,  1 )
and argues rank two at both ends. The middle column is the uncorrected one.
With the increment carried it is -(theta + Delta) K_m m/(K_m+m)^2, and at the
LOW end that changes the argument: the column is then proportional to
(theta + Delta_i) m_i rather than to m_i, so it is no longer parallel to the
other two whenever Delta varies across the record.

At the HIGH end nothing changes: g -> 1 makes the theta column constant and
therefore parallel to the load-bias column, whatever Delta does.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
import design.run as D

P = D.P


def cols(m, dlt, corrected):
    eff = P.theta + (dlt if corrected else 0.0)
    return np.stack([P.theta * m / (P.Km + m),
                     -eff * P.Km * m / (P.Km + m) ** 2,
                     P.ke * m,
                     np.ones_like(m)], axis=1)


def sd_theta(m, dlt, corrected, sigma, priors=True):
    """
    The quantity Proposition 2 is about: how well theta is pinned. Rank is the
    wrong diagnostic here -- both ends are LIMITS, not exact degeneracies, so a
    numerical rank at any finite band is four and says nothing.
    """
    Q = cols(m, dlt, corrected)
    F = Q.T @ Q / sigma ** 2
    if priors:
        F = F + np.linalg.inv(np.diag([0.8 ** 2, 0.8 ** 2, 0.8 ** 2, (0.05 * 322.7) ** 2]))
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    C = np.linalg.inv(F)
    return float(np.sqrt(C[0, 0])), float(np.sqrt(C[1, 1]))


def main():
    print(__doc__)
    n = 40
    print("=" * 88)
    print("THE LOW END.  A band at m/K_m = 0.01, with the increment varying over a")
    print("dose cycle (Delta/theta from 0 to 2) as the trial schedule delivers it.")
    print("=" * 88)
    import model_A.run as LOC
    sig = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)
    m = np.geomspace(0.005, 0.02, n) * P.Km
    print(f"  {'':>40} | {'sd(ln theta)':>13} {'sd(ln K_m)':>12}")
    for lab, dlt in (("Delta absent (a fixed operating point)", np.zeros(n)),
                     ("Delta varying over the dose cycle", np.linspace(0.0, 2.0, n) * P.theta)):
        for corrected in (False, True):
            a, b = sd_theta(m, dlt, corrected, sig)
            tag = "corrected" if corrected else "as printed"
            print(f"  {lab:>40} | {tag:>10} | {a:11.4f} {b:12.4f}")
    print()
    print("  Proposition 2(1) as printed says rank two. With the increment varying and")
    print("  the column corrected it is rank THREE: the therapy partly breaks the")
    print("  sub-saturation degeneracy, because it moves the K_m column without moving")
    print("  the theta column.")

    print()
    print("=" * 88)
    print("THE HIGH END.  A band at m/K_m = 100.")
    print("=" * 88)
    m = np.geomspace(50.0, 200.0, n) * P.Km
    print(f"  {'':>40} | {'sd(ln theta)':>13} {'sd(ln K_m)':>12}")
    for lab, dlt in (("Delta absent", np.zeros(n)),
                     ("Delta varying", np.linspace(0.0, 2.0, n) * P.theta)):
        for corrected in (False, True):
            a, b = sd_theta(m, dlt, corrected, sig)
            tag = "corrected" if corrected else "as printed"
            print(f"  {lab:>40} | {tag:>10} | {a:11.4f} {b:12.4f}")
    print()
    print("  Proposition 2(2) is untouched: theta enters only through L_0 - theta and")
    print("  no increment schedule separates them, because g -> 1 kills the K_m column")
    print("  entirely and makes the theta column constant.")

    print()
    print("=" * 88)
    print("WHAT IT IS WORTH, at the operating point the paper actually uses")
    print("=" * 88)
    t = np.arange(0.0, 730.0, 2.0)
    lv, dd = D.quasi_steady_level(t), np.asarray(D.delta(t), dtype=float)
    import model_A.run as LOC
    sig = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)
    print(f"    m/K_m over the record: {lv.min()/P.Km:.4f} to {lv.max()/P.Km:.4f}"
          f"  (a decade below K_m, so Proposition 2(1) is the relevant end)")
    for corrected in (False, True):
        f = float(np.sqrt(D._filter(lv, sig, dd if corrected else None)[-1][0, 0]))
        print(f"    {'increment carried' if corrected else 'increment dropped':>20}:"
              f"  sd(ln theta) = {f:.4f}")


if __name__ == "__main__":
    main()
