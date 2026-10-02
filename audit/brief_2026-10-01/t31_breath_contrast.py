"""
Item 3: the breath test as the MONITORING instrument.

The one reproducibility claim in Shchelochkov 2021 is for the 1-13C-propionate
breath test, NOT plasma C3; Thompson 1990 called in-vivo oxidation "both
reproducible and prognostically useful". The breath test reads the fraction of
propionate routed through PCC -- i.e. in-vivo oxidation capacity DIRECTLY -- so
(unlike the ex vivo assay, t29) it does NOT inherit the in-vitro/in-vivo mapping
scatter. Its recovery factor c is subject-specific (u3): fatal for cross-
sectional classification, but it CANCELS in a within-patient contrast, exactly
like k_e, K_m, b and kappa.

So compute the longitudinal contrast with the breath test as the observable and
compare to the marker's 0.13 (R-D).

The breath observable is the enzyme's share S = theta g(m) / L (core.share),
which is what the repo already uses for the fraction oxidised. Its sensitivity to
theta is computed from the repo model (convention A: L fixed, m re-solved). The
per-reading repeatability CV_S is NOT a repo figure; the literature gives only
between-subject spreads (Manoli controls 47.9 +/- 4.26%, CV ~9%; Thompson
controls mean 81.9%, range 69.4-101.0, n=5) and a qualitative 'good intrasubject
reproducibility', so CV_S is SWEPT (5-10%) and flagged.
"""
import numpy as np

import core as C


def share_and_elasticity(theta=C.TH_5):
    """S = theta g(m)/L at the untreated balance, and d ln S / d ln theta under
    convention A (L fixed, m re-solved), by central finite difference on theta."""
    import n4_theta_scope as N4

    def S_of(th):
        L, m = N4.anchor(th, C.KM0, C.KE, "A")   # convention A: same diet
        return th * m / (C.KM0 + m) / L, L, m
    S0, L0, m0 = S_of(theta)
    h = 1e-5
    Sp = S_of(theta * (1 + h))[0]
    Sm = S_of(theta * (1 - h))[0]
    elas = (np.log(Sp) - np.log(Sm)) / (2 * h)
    return S0, elas, L0, m0


def main():
    print("=" * 92)
    print("  ITEM 3. THE BREATH TEST AS THE MONITORING INSTRUMENT")
    print("=" * 92)

    S0, elas, L0, m0 = share_and_elasticity()
    print(f"\n  At 5% capacity (repo model, convention A): fraction oxidised S = {S0:.3f}")
    print(f"    operating point L = {L0:.1f}, m = {m0:.1f}")
    print(f"    sensitivity  d ln S / d ln theta = {elas:.3f}")
    print("    (so a 1% change in whole-body capacity moves the breath signal by"
          f" ~{abs(elas):.2f}%.")
    print("     S is read DIRECTLY in vivo -- no in-vitro/in-vivo mapping, unlike t29.)")

    print("\n  WITHIN-PATIENT CONTRAST sd(ln theta_2 - ln theta_1), breath test:")
    print("    per window: sd = CV_S / (|elasticity| * sqrt(n_reads));"
          "  contrast = sqrt(2) * that")
    print("    (c, k_e, K_m all cancel between windows -- patient constants)\n")
    reads = (1, 3, 6)
    print(f"      {'CV_S per reading':>18}  " + "".join(f"n={n:<6}" for n in reads))
    for cv in (0.05, 0.08, 0.10):
        cells = []
        for n in reads:
            contrast = np.sqrt(2) * cv / (abs(elas) * np.sqrt(n))
            cells.append(f"{contrast:7.3f}" + ("* " if contrast <= 0.20 else "  "))
        print(f"      {int(cv*100):>14}%     " + "".join(cells))

    print(f"\n  marker-based contrast for comparison (R-D): 0.26 monthly / 0.13 weekly"
          " / 0.07 daily")
    print("""
  READING: at a plausible ~8% within-subject repeatability, a SINGLE breath
  reading per window gives a contrast of ~0.13 -- matching 104 weekly blood
  draws -- and 3-6 readings beat it (~0.05-0.07). The breath test reaches the
  monitoring target with a handful of readings rather than a dense draw schedule,
  AND it reads in-vivo capacity directly, so it is immune to BOTH the calibration
  problem (c cancels in the contrast) and the in-vitro mapping (t29) that sink the
  ex vivo route cross-sectionally.

  DETECTABLE CHANGE at 80% power (min |Delta ln theta| = 2.80 * contrast):""")
    for cv, n in ((0.08, 1), (0.08, 3), (0.05, 3)):
        contrast = np.sqrt(2) * cv / (abs(elas) * np.sqrt(n))
        mdd = 2.80 * contrast
        print(f"    CV_S {int(cv*100)}%, {n} reads/window: contrast {contrast:.3f} -> "
              f"detect +{100*(np.exp(mdd)-1):.0f}% / -{100*(1-np.exp(-mdd)):.0f}%")

    print("""
  CAVEAT (flagged): the per-reading CV_S is a LITERATURE-informed sweep, not a
  repo-measured or even a published within-subject figure -- Shchelochkov's
  'good intrasubject reproducibility' and Thompson's 'reproducible' are
  qualitative. The elasticity 0.89 IS from the repo model at 5%. The headline --
  breath beats the marker for monitoring and dodges both calibration problems --
  holds unless the true within-subject CV_S is much worse than ~10% (at 15% a
  single reading gives ~0.24, still a few reads reach target). This is a material
  revision to R-D: IF a within-subject CV_S is measured, the monitoring
  instrument should probably be the breath test, not plasma 3-HP.
""")


if __name__ == "__main__":
    main()
