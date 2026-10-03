"""
Round 10, item 1: R-D (the monitoring contrast) against the only empirical anchor
for serial within-patient variation.

No PA serial CV exists. The closest anchor (round-10 literature) is healthy adults:
within-individual same-day C3 CV = 26.4% (range 8.2-61.1), against 7.9% analytical
and 39.2% between-subject (Thompson 2012, Daily Variation of Serum Acylcarnitines).
R-D's 0.13 assumed a ~5% per-draw assay. Recompute the contrast with the per-draw
marker scatter set to the SERIAL CV, swept over 8.2-61.1%, and report the draws per
window needed for 0.20 and the detectable change.

The serial CV is the per-draw INDEPENDENT marker scatter (biological + preanalytical
+ analytical); it enters as the assay term in the per-draw flux noise. The load-
record bias b stays a shared constant that cancels in the contrast (R-D).
"""
import numpy as np

import core as C
import contrast as K

ZPOW = 1.959963985 + 0.841621234   # 80% power, two-sided 5%


def contrast_at(cv, n, win=730, seeds=range(24)):
    """Weekly-style contrast with per-draw marker scatter = cv (load bias shared)."""
    def one(s):
        r1 = C.record(C.TH_5, seed=2 * s, times=np.linspace(0, win, n))
        r2 = C.record(C.TH_5, seed=2 * s + 1, times=np.linspace(0, win, n))
        L = r1["L"]
        F1 = C.fim(r1, noise="perdraw", ke_prior=0.0, b_prior=np.inf, assay_cv=cv)
        F2 = C.fim(r2, noise="perdraw", ke_prior=0.0, b_prior=np.inf, assay_cv=cv)
        return K.window_contrast([F1, F2], ["theta", "Km", "b", "ke"],
                                 priors={"ke": 0.8, "b": 0.05 * L})
    return float(np.median([one(s) for s in seeds]))


def main():
    print("=" * 92)
    print("  ROUND 10 ITEM 1. R-D AGAINST THE SERIAL-CV ANCHOR (Thompson: C3 CVi 26.4%, 8.2-61.1)")
    print("=" * 92)
    print("\n  Weekly contrast (n=104/window) as the per-draw marker CV rises from the")
    print("  repo's ~5% assay to the empirical serial CV:")
    print(f"      {'per-draw marker CV':>20}  {'contrast sd':>12}  {'draws for 0.20':>15}  {'detect +/-':>18}")
    base_n = 104
    for cv, lab in [(0.05, "repo assay 5%"), (0.082, "serial low 8.2%"),
                    (0.264, "serial mean 26.4%"), (0.611, "serial high 61.1%")]:
        sd = contrast_at(cv, base_n)
        n20 = base_n * (sd / 0.20) ** 2        # contrast ~ 1/sqrt(n)
        mdd = ZPOW * 0.20                       # detectable change AT the 0.20 gate
        # detectable change at the achieved sd:
        mdd_here = ZPOW * sd
        up = 100 * (np.exp(mdd_here) - 1); dn = 100 * (1 - np.exp(-mdd_here))
        print(f"      {lab:>20}  {sd:>12.3f}  {n20:>15.0f}  {('+%.0f%% / -%.0f%%' % (up, dn)):>18}")

    print("""
  READING: at the repo's 5% assay the weekly contrast is ~0.13 (R-D). At the only
  empirical serial CV (26.4%, healthy adults), the SAME 104-draw window gives a
  contrast of ~0.48, and reaching 0.20 needs of order ~590 draws per window --
  not deliverable in a monitoring programme. At the low end of the range (8.2%) it
  is ~0.18 at 104 draws (borderline); at the high end (61.1%) it is hopeless.
  So R-D's 0.13 is a best-case floor at an assay CV far below the only measured
  serial variability.
""")
    print("  Draws per window for a 0.20 contrast, across the serial-CV range:")
    print(f"      {'serial CV':>12}  {'draws/window':>14}")
    for cv in (0.082, 0.15, 0.264, 0.40, 0.611):
        sd104 = contrast_at(cv, 104)
        n20 = 104 * (sd104 / 0.20) ** 2
        print(f"      {int(cv*1000)/10:>11}%  {n20:>14.0f}")

    print("""
  TWO-SIDED UNCERTAINTY (stated honestly, as the anchor is weak in both directions):
  - The 26.4% is HEALTHY ADULTS, same-day, UNSTANDARDISED, on a small labile plasma
    C3 pool -- it captures feeding/activity swing and likely OVERSTATES the week-to-
    week stable-state CV a standardised PA clinic would see.
  - But PA patients have grossly elevated C3 (thousands of nmol/L vs ~450), and a
    larger, less labile pool may vary LESS in relative (CV) terms -- so the true PA
    serial CV could be below 26.4%.
  - No PA serial CV exists in any published dataset (the two candidate papers hold
    only summary statistics or control-scaled Z-scores; the one cohort with the raw
    data, Moderna n=50, withdrew follow-up biomarkers "due to high variability").
  So R-D should be given as a RANGE over serial CV, not a point at ~0.13: the
  monitoring result is deliverable (<~0.20 at ~100-200 draws) only if the true PA
  serial CV is near the LOW end (<~10%); at the measured mean it is not. The
  decisive measurement -- a standardised within-patient CV in stable PA -- has
  never been made, and it is the single number the monitoring claim turns on.
""")


if __name__ == "__main__":
    main()
