"""
Round 16: can human HEPATIC PCC turnover be bridged from what we have?

We have the human FIBROBLAST PCC half-life (Chandler 1985, 3.6 d) and mouse values
in both tissues. The bridge: human_liver = human_fibroblast * (liver/fibroblast
tissue ratio), with the ratio taken from mouse. The decision turns on which side of
the round-15 crossing (t_half* = 3.5452 d) the bridged estimate falls.

THE DATA, AND WHY THE RATIO IS NOT CLEAN (this is itself the answer).
  human fibroblast PCC : 3.6 d   Chandler 1985, [3H]biotin pulse-chase (DEGRADATION)
  mouse whole-liver med: 3.6 d   Rolfs 2021 SILAC; 3.5 d Price 2010 15N; 4.26 d Kim 2012 mito
  mouse fibroblast med : 2.5 d   Swovick 2018 dynamic-SILAC, CULTURED primary dermal
  mouse liver PCC      : 3.8 d intrinsic / ~7 d effective (Attarwala 2023); 5.3-6.7 d (Jiang 2020)
No single study/method has BOTH mouse fibroblast and mouse liver, and no mouse
FIBROBLAST PCC exists -- so a PCC-specific tissue ratio is unavailable, and the only
ratio is a cross-study ratio of tissue MEDIANS. Worse, it is confounded: the mouse
fibroblast median (2.5 d) is from CULTURED, dividing cells, where division dilutes
protein and inflates the apparent turnover rate, whereas the liver median is in-vivo.
So the face-value median ratio over-states how much faster fibroblasts turn over.

Three defensible readings of the liver/fibroblast ratio r:
  A  cross-study medians at face value : r = 3.6/2.5 = 1.44 (liver SLOWER)
  B  like-for-like degradation measures: human fibroblast PCC 3.6 d ~= mouse liver
     median 3.6 d ~= mouse liver PCC intrinsic 3.8 d  ->  r ~= 1.0
  C  "liver is the fastest-turnover tissue in vivo" (Rolfs), so hepatic <= fibroblast
     ->  r ~= 0.8 (liver FASTER)
The readings disagree on the DIRECTION. That disagreement is the result.

Uncertainty carried (all multiplicative / log-normal):
  - human fibroblast anchor: sd(ln) = 0.20 (single 1985 value, reproduced in 2 lines)
  - the tissue ratio r: sd(ln) = 0.35 (cross-study, cross-method, culture-vs-vivo)
  - per-PROTEIN deviation from the tissue-median ratio (a median applied to ONE
    protein): sd(ln) = 0.35 (within-tissue protein half-lives span ~1 order; Swovick
    cross-tissue conservation is only moderate)
"""
import numpy as np
from scipy.stats import norm

CROSS = 3.5452                    # round-15 crossing half-life (21-day safe-for-all floor)
FIB = 3.6                         # human fibroblast PCC, Chandler 1985
SD_ANCHOR, SD_RATIO, SD_PROTEIN = 0.20, 0.35, 0.35

READINGS = {
    "A  cross-study medians (r=1.44, liver slower)": 1.44,
    "B  like-for-like degradation (r=1.0)": 1.00,
    "C  liver fastest in vivo (r=0.8, liver faster)": 0.80,
}


def p_below(r_center, n=400000, seed=0, extra_ratio_sd=0.0):
    """P(human hepatic PCC half-life < crossing), log-normal MC."""
    rng = np.random.default_rng(seed)
    ln_hep = (np.log(FIB * r_center)
              + rng.normal(0, SD_ANCHOR, n)
              + rng.normal(0, np.hypot(SD_RATIO, extra_ratio_sd), n)
              + rng.normal(0, SD_PROTEIN, n))
    hep = np.exp(ln_hep)
    return float(np.mean(hep < CROSS)), float(np.median(hep)), np.percentile(hep, [5, 95])


def main():
    print("=" * 92)
    print("  ROUND 16. BRIDGING HUMAN HEPATIC PCC TURNOVER, AND P(below the 3.5452 d crossing)")
    print("=" * 92)
    print(f"""
  Bridge: human_liver t_half = human_fibroblast (3.6 d) * (liver/fibroblast ratio r),
  r taken from mouse. Crossing = {CROSS} d (round 15). Uncertainty: anchor sd(ln)=0.20,
  ratio sd(ln)=0.35, per-protein sd(ln)=0.35 -> combined sd(ln) ~ {np.hypot(np.hypot(SD_ANCHOR,SD_RATIO),SD_PROTEIN):.2f}
  (a factor ~{np.exp(np.hypot(np.hypot(SD_ANCHOR,SD_RATIO),SD_PROTEIN)):.1f} one-sigma band -- deliberately wide).
""")
    print(f"      {'ratio reading':>48}  {'bridged median':>14}  {'5-95% band':>18}  {'P(hepatic<cross)':>16}")
    for lab, r in READINGS.items():
        p, med, band = p_below(r)
        print(f"      {lab:>48}  {med:>12.1f} d  {f'[{band[0]:.1f}, {band[1]:.1f}]':>18}  {p:>15.0%}")

    # combined honest reading: center at like-for-like 1.0 but widen the ratio sd to
    # span the A/C direction disagreement (ln 1.44 and ln 0.8 are +-0.36 and -0.22)
    p, med, band = p_below(1.0, extra_ratio_sd=0.30)
    print(f"\n      {'COMBINED (center r=1, ratio sd widened to span A-C)':>48}  "
          f"{med:>12.1f} d  {f'[{band[0]:.1f}, {band[1]:.1f}]':>18}  {p:>15.0%}")

    print(f"""
  READING. There is NO clean paired mouse fibroblast-vs-liver dataset for PCC: the
  fibroblast numbers are cultured (Swovick) or [3H]biotin (Chandler), the liver numbers
  in-vivo (Rolfs/Price/Kim) or exogenous-mRNA (Jiang/Attarwala); no study has both, and
  no mouse fibroblast PCC exists. The only available ratio is a cross-study ratio of
  tissue MEDIANS, confounded by culture-vs-in-vivo (division dilutes cultured protein and
  inflates its apparent turnover). The three defensible readings DISAGREE ON DIRECTION --
  liver slower (A), equal (B), or faster (C) -- so the bridged estimate spans ~2.9-5.2 d
  and P(hepatic < crossing) spans ~{p_below(0.8)[0]:.0%} to ~{p_below(1.44)[0]:.0%}, centred near 50%. The bridge
  therefore CANNOT determine which side of the 3.5452 d crossing human hepatic PCC falls
  on; it is consistent with either. The striking like-for-like coincidence -- human
  fibroblast PCC (3.6 d) = mouse liver median (3.6 d) = mouse liver PCC intrinsic (3.8 d)
  -- argues the tissue effect is small (r ~= 1), which lands the estimate ~3.6 d, almost
  exactly ON the crossing, where the probability is ~50/50 by construction.

  CONCLUSION: the bridge is uninformative about the side of the crossing -- which is a
  real answer: the trial record is the only route. Round 12 shows the on-therapy record
  reads lambda at 7-18 sigma (separating 3.6 from 7 d decisively), and it separates the
  two sides of a 3.5452 d crossing a fortiori. The single decision-relevant measurement --
  human HEPATIC PCC turnover -- cannot be bridged from fibroblast data and must be READ,
  either from a dedicated hepatic measurement or from NCT04159103's accumulated record.
""")


if __name__ == "__main__":
    main()
