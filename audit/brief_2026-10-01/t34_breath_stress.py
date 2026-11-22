"""
Round 7, item 3: stress-test the breath test the way R-D stress-tested the marker.

t31 found a single breath reading at ~8% repeatability matches 104 weekly draws.
That is the project's strongest result and had one round. Give it R-D's treatment:
provenance of the 8%, c-drift breakpoint, control-spread uncertainty, confounder
inheritance, and the cross-sectional behaviour with c uncalibrated vs calibrated.

The elasticity d ln S/d ln theta = 0.888 at 5% is from t31 (repo model).
"""
import numpy as np

import core as C
import t31_breath_contrast as T31

ELAS = abs(T31.share_and_elasticity()[1])   # 0.888


def contrast(cv_s, n):
    return np.sqrt(2) * cv_s / (ELAS * np.sqrt(n))


def main():
    print("=" * 92)
    print("  ROUND 7 ITEM 3. BREATH TEST STRESS TESTS (elasticity = %.3f)" % ELAS)
    print("=" * 92)

    print("""
  (a) PROVENANCE OF THE 8%. It is NOT a measured within-subject figure. The
      literature gives only BETWEEN-subject spreads: Manoli 2021 controls
      47.9 +/- 4.26% (CV 8.9%), heterozygotes 41.6 +/- 5.05% (CV 12.1%); Thompson
      1990 controls mean 81.9%, range 69.4-101.0 (n=5). Within-subject
      reproducibility is documented only QUALITATIVELY ('good intrasubject
      reproducibility', Shchelochkov 2021; 'reproducible', Thompson 1990). So the
      8% is a between-subject CV used as a stand-in and SWEPT 5-10%; the headline
      needs a measured within-subject CV_S to stand, and that number does not yet
      exist. This is the same evidence gap as the marker's missing within-patient
      CV (R-D) -- but here at least the qualitative claim is on record for the
      right quantity.""")

    print("  (b) c-DRIFT BREAKPOINT. c cancels only if it is STABLE between windows.")
    print("      Subject-specific is not temporally stable: c is the 1-13C recovery,")
    print("      set by gut absorption + oxidation, so it can drift with diet, gut flora")
    print("      or an antibiotic course. A between-window drift d(ln c) biases the")
    print("      contrast by d(ln c)/elasticity; rmse = hypot(base, d(ln c)/elasticity).")
    print(f"      {'base (reads @CV_S)':>24}  {'c-drift that breaks 0.20':>26}")
    for cv, n, lab in ((0.08, 1, "1 read @ 8%"), (0.08, 3, "3 reads @ 8%"),
                       (0.05, 3, "3 reads @ 5%")):
        base = contrast(cv, n)
        dc = ELAS * np.sqrt(max(0.20 ** 2 - base ** 2, 0.0))
        print(f"      {lab:>24}  d(ln c) = {dc:.3f}  (c drift {100*(np.exp(dc)-1):+.0f}%)")
    print("      So the breath contrast tolerates a ~14-18% between-window drift in the")
    print("      recovery c before the bias alone reaches target -- wider than the marker's")
    print("      k_e/b tolerances (R-D), but NOT unlimited, and c-drift is plausible across")
    print("      an antibiotic course or a gut-flora change (the same events that corrupt")
    print("      the load record, R-E). c must be held as fixed as k_e and the diet.")

    print("\n  (c) CONTROL-SPREAD UNCERTAINTY (Thompson n=5). For MONITORING it is")
    print("      irrelevant -- c and the control mean cancel in the within-patient")
    print("      contrast. It matters only for a CROSS-SECTIONAL absolute reading, where")
    print("      the control mean carries SE = spread/sqrt(n_controls):")
    sd_ctrl = (np.log(101.0) - np.log(69.4)) / 4.0   # Thompson range as +/-2SD, log
    print(f"      Thompson controls: ln-spread ~ {sd_ctrl:.3f} (range 69.4-101.0 as +/-2SD),")
    for nc in (5, 10, 30):
        print(f"        n={nc:>2} controls: control-mean SE = {sd_ctrl/np.sqrt(nc):.3f}")
    print("      At Thompson's n=5 the control-mean SE alone is ~0.04 -- small next to the")
    print("      cross-sectional c-calibration problem below, but it says a 5-control")
    print("      reference is thin; ~15-30 controls make it negligible.")

    print("\n  (d) CONFOUNDER INHERITANCE. The breath test reads the OXIDISED FRACTION of")
    print("      an administered dose, normalised to measured CO2 production (indirect")
    print("      calorimetry, Manoli). So it does NOT inherit the marker-pool confounders")
    print("      that move C3: renal clearance (eGFR) and carnitine dose act on the")
    print("      accumulated marker, not on the oxidation fraction. It IS sensitive to")
    print("      fed/fasted state and activity (substrate competition, REE) -- but those")
    print("      are controllable WITHIN a test protocol (fasting, rest, REE normalisation),")
    print("      unlike the uncontrolled month-scale drift of eGFR/carnitine. Net: fewer")
    print("      month-scale confounders than C3, at the cost of a standardised test visit.")

    print("\n  (e) CROSS-SECTIONAL behaviour, c uncalibrated vs calibrated (the R11 contour).")
    print("      A single breath reading gives S = c * fraction(theta,...); cross-")
    print("      sectionally c is the SAME uncalibrated multiplier as the flux route.")
    print("      sd(ln theta) = hypot(CV_S, Z)/elasticity, Z = sd(ln c):")
    print(f"      {'Z = sd(ln c)':>16}  {'CV_S=5%':>9}{'8%':>9}{'10%':>9}")
    for Z, zl in ((0.0, "c known"), (0.10, "0.10"), (0.20, "0.20"), (np.inf, "c free")):
        cells = ""
        for cv in (0.05, 0.08, 0.10):
            if not np.isfinite(Z):
                cells += f"{'inf':>9}"
            else:
                cells += f"{np.hypot(cv, Z)/ELAS:9.3f}"
        print(f"      {zl:>16}  {cells}")
    print("      So cross-sectionally the breath test has the SAME calibration problem as")
    print("      every direct route (R11/R-G): c known -> reaches ~0.06-0.11; c uncalibrated")
    print("      -> worthless. Its advantage is entirely longitudinal, where c cancels. The")
    print("      cross-sectional verdict (only a CALIBRATED direct in-vivo read reaches")
    print("      target) is unchanged; the breath test is that read, and its calibration is")
    print("      the open requirement -- but UNLIKE the ex vivo route it has no in-vitro/")
    print("      in-vivo mapping gap (t29): it reads capacity in vivo by construction.")


if __name__ == "__main__":
    main()
