"""
Item 11. The flux observable, as a specification (conditional analysis only).

No claim is made about which biological species is measurable or what any assay
reads. The observable is treated as a specification and the question is what it
would have to satisfy.

(1) A generic observable O = c * theta * g(m) with an unknown calibration c.
    sd(ln theta_hat) at 5% as a function of the reading precision Y and the
    precision Z on c; the 0.20 contour; the one-sentence requirement.

(2) The reparameterisation: if the reading is the fraction oxidised
    s = theta g(m) / L rather than the flux itself, then (s to Y%, L to X%)
    gives the flux as s*L. Does the 0.20 crossing survive, and at what X, Y?

(3) What the repository's breath test, which reads a fraction oxidised with an
    unknown recovery, can and cannot deliver -- stated plainly, no speculation
    about other assays.
"""
import numpy as np

import core as C
import t5_value as T5
import t12_flux_mc as T12

GATE = 0.20


# -------- (2) the (s, L) pair: fraction-oxidised reading + a load measurement
def share_at(th, K, ke, b, Lrec):
    m0 = C.m_ss(th, K, ke, Lrec - b)
    return th * m0 / (K + m0) / (Lrec - b)


def share_row(r, h=1e-6):
    """s and d s / d(ln theta, ln K_m, b, ln k_e) at the untreated baseline."""
    th, K, ke, Lrec = r["theta"], r["Km"], r["ke"], r["L"]
    base = np.array([np.log(th), np.log(K), 0.0, np.log(ke)])
    s0 = share_at(th, K, ke, 0.0, Lrec)
    g = np.zeros(4)
    for j in range(4):
        hj = h * (Lrec if j == C.I_B else 1.0)
        xp, xm = base.copy(), base.copy()
        xp[j] += hj
        xm[j] -= hj

        def sh(x):
            return share_at(np.exp(x[0]), np.exp(x[1]), np.exp(x[3]), x[2], Lrec)
        g[j] = (sh(xp) - sh(xm)) / (2 * hj)
    return s0, g


def sd_sL(Y, X):
    """median sd(ln theta) from the record + one s-reading to Y% + L to X%."""
    def one(seed):
        r = T5.base_record(seed, baseline=True)
        s0, g = share_row(r)
        return C.sd_th(C.fim(r, b_prior=max(X, 1e-7) * r["L"],
                             extra=[(g, Y * s0)]))
    return C.over_seeds(one)[0]


def crossing_in_Y(f, lo=0.01, hi=0.60):
    if f(lo) > GATE:
        return None
    if f(hi) <= GATE:
        return hi
    for _ in range(46):
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if f(mid) <= GATE else (lo, mid)
    return lo


def main():
    print("=" * 96)
    print("  ITEM 11.  THE FLUX OBSERVABLE AS A SPECIFICATION (CONDITIONAL)")
    print("=" * 96)

    # ---- (1) O = c theta g(m): sd over (Y, Z), the 0.20 contour ----------
    print("\n  (1) GENERIC OBSERVABLE O = c * theta * g(m), calibration c unknown to Z.")
    Ys = [0.05, 0.10, 0.15, 0.175, 0.20, 0.25]
    Zs = [1e-7, 0.05, 0.10, 0.15, 0.20]
    print("      sd(ln theta), median over 24 seeds (rows Z = sd(ln c), cols Y):")
    print("        Z \\ Y   " + "".join(f"{y*100:6.1f}%" for y in Ys))
    worst = 0.0
    for Z in Zs:
        row = []
        for Y in Ys:
            v = T12.fim_seeds(Y=Y, cal_sd=(None if Z < 1e-6 else Z))
            row.append(v)
            # hypot collapse check against a single calibrated reading at Y_eff
            yeff = float(np.hypot(Y, 0.0 if Z < 1e-6 else Z))
            ref = T12.fim_seeds(Y=yeff, cal_sd=None)
            worst = max(worst, abs(v - ref))
        lab = "exact" if Z < 1e-6 else f"{Z:.2f}"
        print(f"      {lab:>7}   " + "".join(f"{v:6.3f}" for v in row))
    print(f"      the 5x5 FIM equals a single calibrated reading at "
          f"Y_eff = hypot(Y, Z): worst gap over the grid {worst:.1e}")
    #  the contour: hypot(Y, Z) <= Y*, Y* the single-reading crossing
    ystar = T5.crossing(lambda y: T5.sd5(flux_Y=y), 0.05, 0.40)
    print(f"      single-reading FIM crossing Y* = {100*ystar:.2f}% "
          f"(Monte Carlo Y* = 18.4%, t12).")
    print(f"      => 0.20 CONTOUR:  hypot(Y, Z) <= {100*ystar:.1f}% (FIM) "
          f"/ 18.4% (Monte Carlo).")
    print(f"      REQUIREMENT SENTENCE: an observable proportional to the enzymatic")
    print(f"      flux, read to Y% with calibration known to Z%, reaches sd(ln theta)")
    print(f"      = 0.20 at 5% of normal if and only if sqrt(Y^2 + Z^2) <= ~18-19%.")
    print(f"      (Averaging n readings divides Y by sqrt(n) but NOT Z: a single")
    print(f"       shared calibration c is the floor. Three readings at Y need")
    print(f"       Z <= {100*ystar:.1f}% on its own to pass.)")

    # ---- (2) the (s, L) reparameterisation -------------------------------
    print("\n  (2) REPARAMETERISATION: fraction oxidised s = theta g(m)/L, with L.")
    print("      (s to Y%, L to X%); flux = s*L. sd(ln theta), median over seeds:")
    Xs = [1e-7, 0.01, 0.02, 0.05]
    print("        X \\ Y   " + "".join(f"{y*100:6.1f}%" for y in Ys))
    for X in Xs:
        row = [sd_sL(Y, X) for Y in Ys]
        lab = "exact" if X < 1e-6 else f"{X*100:.0f}%"
        print(f"      L {lab:>5}  " + "".join(f"{v:6.3f}" for v in row))
    print("      0.20 crossing in Y, per load precision X:")
    for X, lab in ((1e-7, "L exact"), (0.01, "L to 1%"), (0.02, "L to 2%"),
                   (0.05, "L to 5%")):
        yc = crossing_in_Y(lambda Y: sd_sL(Y, X))
        if yc is None:
            print(f"        {lab:>8}: never reaches 0.20 for Y in [1%, 60%]")
        else:
            print(f"        {lab:>8}: s to Y <= {100*yc:.2f}%   "
                  f"(hypot(Y,X) = {100*np.hypot(yc, X):.1f}%)")
    print("      So the (s, L) pair reaches target on the SAME locus: s and L")
    print("      enter as ln(flux) = ln s + ln L, so a load measurement to X% plays")
    print("      exactly the role of calibrating c to Z%. There is no free lunch in")
    print("      switching to the fraction oxidised; it only moves the calibration")
    print("      burden onto an independent measurement of L.")

    # ---- (3) the breath test, stated plainly -----------------------------
    print("\n  (3) WHAT THE REPOSITORY'S BREATH TEST CAN AND CANNOT DELIVER.")
    print("      Finding (round 2): u2_breath_functional.py models it as y = c*F(m),")
    print("      F = theta g/(theta g + k_e m) = the fraction oxidised = s, with the")
    print("      recovery c UNKNOWN. On the specification above:")
    print("      - It can deliver the fraction oxidised s (a relative reading), i.e.")
    print("        a Y. It CANNOT deliver the absolute flux: c is its unknown Z.")
    print("      - With c entirely unknown (Z -> inf), one reading is worth nothing")
    print("        (t12): sqrt(Y^2 + Z^2) -> inf, off the contour.")
    print("      - It reaches target ONLY if BOTH (i) its recovery c is calibrated to")
    print("        Z <= ~18% AND (ii) the load L is measured to X, with hypot of the")
    print("        reading's Y, the calibration Z and the load X <= ~18-19%.")
    print("      Whether such a calibration exists is a domain question, not answered")
    print("      here. No other assay is assumed.")


if __name__ == "__main__":
    main()
