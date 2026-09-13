"""
The excitation functional, corrected.

Section 6 builds everything on the ratio of the two sensitivities,

    (dC/d ln K_m) / (dC/d ln theta)  =  -phi(m),     phi = K_m/(K_m+m),

and defines the excitation as the DROP IN phi across the visited band. That
ratio is computed from the uncorrected K_m column. The correct one is

    (dC/d ln K_m) / (dC/d ln theta)  =  -(1 + Delta(t)/theta) * phi(m) ,

so what a design has to vary is the PRODUCT, and there are two factors to move.
Moving m is the excursion route and needs a decompensation. Moving Delta is the
therapy route, and the dosing schedule already does it.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
import design.run as D
from common.excitation import phi, excitation, excitation_max

P = D.P


def ratio_range(m, dlt):
    r = (1.0 + dlt / P.theta) * phi(m, P.Km)
    return float(r.min()), float(r.max()), float(r.max() - r.min())


def main():
    print(__doc__)
    t = np.arange(0.0, 730.0, 2.0)
    lv = D.quasi_steady_level(t)
    dd = np.asarray(D.delta(t), dtype=float)

    print("=" * 78)
    print("THE TWO FACTORS, over the record the trial schedule already produces")
    print("=" * 78)
    ph = phi(lv, P.Km)
    fa = 1.0 + dd / P.theta
    print(f"  marker levels            {lv.min():8.3f} to {lv.max():8.3f}"
          f"   (m/K_m {lv.min()/P.Km:.4f} to {lv.max()/P.Km:.4f})")
    print(f"  phi(m)                   {ph.min():8.4f} to {ph.max():8.4f}"
          f"   range {ph.max()-ph.min():.4f}   ratio {ph.max()/ph.min():.3f}x")
    print(f"  1 + Delta/theta          {fa.min():8.4f} to {fa.max():8.4f}"
          f"   range {fa.max()-fa.min():.4f}   ratio {fa.max()/fa.min():.3f}x")
    lo, hi, rng = ratio_range(lv, dd)
    print(f"  the product              {lo:8.4f} to {hi:8.4f}   range {rng:.4f}")
    print()
    print(f"  E as Section 6 defines it (the drop in phi alone): "
          f"{excitation(lv.min(), lv.max(), P.Km):.4f}")
    print(f"  the drop in the product, same record:              {rng:.4f}"
          f"   -- a factor of {rng/(ph.max()-ph.min()):.1f}")
    print(f"  against the E ~ 0.50 the decision was said to need.")

    print()
    print("=" * 78)
    print("WHERE THE TWO ROUTES REACH E = 0.50")
    print("=" * 78)
    print("  excursion route: carry the marker up from the operating point")
    for k in (2, 5, 8, 20, 40):
        print(f"    to {k:3d} m_bar:  drop in phi = "
              f"{excitation(D.M_BAR, k*D.M_BAR, P.Km):.4f}"
              f"   ({100*excitation(D.M_BAR, k*D.M_BAR, P.Km)/0.50:.0f}% of the requirement)")
    print()
    print("  therapy route: hold the marker where it is and vary the increment")
    for a in (0.25, 0.5, 1.0, 2.0, 2.13, 4.0):
        r = (1.0 + a) * phi(D.M_BAR, P.Km) - (1.0 + 0.0) * phi(D.M_BAR, P.Km)
        print(f"    Delta/theta 0 to {a:4.2f}:  drop in the product = {r:.4f}"
              f"   ({100*r/0.50:.0f}% of the requirement)")
    print()
    print(f"  The trial schedule already delivers Delta/theta from "
          f"{dd.min()/P.theta:.2f} to {dd.max()/P.theta:.2f}.")
    print("  It is administering the excitation and nobody was reading it as such.")


if __name__ == "__main__":
    main()
