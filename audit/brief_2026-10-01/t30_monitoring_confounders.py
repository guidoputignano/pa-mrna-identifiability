"""
Item 2: the monitoring result (R-D) against the hostile empirical literature.

R-D's 0.13 weekly contrast is a WITHIN-DRAW-NOISE-ONLY floor: it assumes the only
non-shared variance between the two windows is the per-draw assay+load noise
(5%+5%), and that every patient-level nuisance (k_e, K_m, b, kappa, c) is
perfectly shared. The literature says that last assumption is the weak point:

  - No published within-patient CV exists for 3-HP, C3 or C3/C2 in stable PA.
  - The Moderna n=50 natural-history cohort collected these every 6 months for
    36 months and restricted analysis to BASELINE "due to small sample sizes and
    high variability during follow-up" -- the empirical analogue of the
    monitoring claim, and it did not deliver a usable follow-up signal.
  - C3 rises with levocarnitine dose (r=0.41) and with declining eGFR
    (Shchelochkov 2021). Both drift over months.
  - Preanalytical standardisation (time of day, fed state, timing vs dosing) is
    weak or unreported in every PA cohort.

These confounders are exactly the NON-SHARED drift the contrast is vulnerable to.
They map onto the two drift axes already stress-tested in t24:
  - declining eGFR  ->  the marker's renal clearance k_e drifts  (k_e-drift axis)
  - levocarnitine dose change / preanalytical  ->  the C3 calibration and the
    load-record bias drift                                        (b-drift axis)
This script re-expresses t24's breakpoints in those terms and states precisely
what raw data would test the model's per-draw-noise prediction.
"""
import numpy as np

import t24_contrast_stress as T24


def main():
    print("=" * 92)
    print("  ITEM 2. THE MONITORING CONTRAST vs THE CONFOUNDERS THE LITERATURE NAMES")
    print("=" * 92)

    base_w = T24.med_sd(104, 730)
    base_m = T24.med_sd(25, 730)
    print(f"\n  Within-draw-noise-only floor (R-D): weekly {base_w:.3f}, monthly {base_m:.3f}.")
    print("  This is the BEST case: it assumes k_e, b and the calibrations are shared")
    print("  exactly between windows. The literature confounders break that sharing.\n")

    print("  (a) eGFR decline -> k_e drift (renal clearance of the marker):")
    print(f"      {'k_e drift (log units)':>22}  {'weekly contrast':>15}")
    for dk in (0.0, 0.03, 0.06, 0.10):
        v = T24.med_sd(104, 730, drift_ke=dk)
        print(f"      {dk:>22.2f}  {v:15.3f}{'   <= gate' if v <= 0.20 else '   > gate'}")
    print("      Shchelochkov 2021: plasma C3 and 2-MC rise significantly as eGFR falls.")
    print("      In a growing child renal clearance tracks size and eGFR can decline over")
    print("      months, so a k_e drift of order 0.05-0.10 between windows is plausible --")
    print("      which is already past the weekly breakpoint (~0.035).")

    print("\n  (b) levocarnitine dose change / preanalytical -> calibration (b) drift:")
    print(f"      {'b drift (fraction L)':>22}  {'weekly contrast':>15}")
    for db in (0.0, 0.015, 0.03, 0.05):
        v = T24.med_sd(104, 730, drift_b=db)
        print(f"      {db:>22.3f}  {v:15.3f}{'   <= gate' if v <= 0.20 else '   > gate'}")
    print("      Shchelochkov 2021: plasma C3 rises with levocarnitine dose (r=0.41).")
    print("      A carnitine dose change or a shift in sampling conditions (fed/fasted,")
    print("      time since dose) between windows moves the calibration; past ~0.015-0.03")
    print("      L the weekly contrast exceeds the gate.")

    print("""
  REFRAME OF R-D: the 0.13 is a model-predicted FLOOR under perfect nuisance
  sharing, and it is EMPIRICALLY UNTESTED. The one cohort that collected the right
  data longitudinally (Moderna n=50) found the markers too variable to use past
  baseline -- a direct, if indirect, challenge. The contrast's own stress tests
  show it survives only small between-window drift, and the literature names real
  month-scale confounders (eGFR, carnitine, preanalytical) that plausibly exceed
  those breakpoints. So the monitoring claim must be stated as CONDITIONAL on
  short intervals, a stable regimen (fixed carnitine, stable renal function) and
  standardised sampling -- and as awaiting an empirical within-patient CV.
""")
    print("  WHAT WOULD TEST IT (the most valuable outstanding item):")
    print("""  From Haijes et al. 2020 (11 PA, longitudinal deep phenotyping) and Cappuccio
  et al. (serial C3 per patient), extract, for each patient, the serial 3-HP (or
  C3) values from CLINICALLY STABLE periods only (no decompensation), with:
    - sampling date/time, fed vs fasted state, time since last feed and last dose,
    - levocarnitine dose and renal function at each draw, same lab/assay throughout.
  Then (i) compute the raw within-patient CV of the marker in stable periods and
  compare it to the model's per-draw 5% assay + 5% load assumption; (ii) partition
  it into within-window (days) vs across-window (months) variance -- the first is
  the per-draw noise the 0.13 rests on, the second is the non-shared drift that
  inflates the contrast; (iii) only the across-window, confounder-adjusted
  residual tests whether sd(ln theta_2 - ln theta_1) really approaches ~0.13.
  No published study has done this; the raw material exists in those two papers.""")


if __name__ == "__main__":
    main()
