"""
Round 17: the hepatic correction runs the OTHER way (a correction to round 15/16).

Rounds 15 (R-N) and 16 (R-O) carried the user's earlier steer that human HEPATIC PCC
"likely turns over FASTER than fibroblast" because liver is the fastest-turnover TISSUE
(Rolfs). That used the WHOLE-LIVER median (2.4-3.6 d). But PCC is a mitochondrial-MATRIX
enzyme, and the apt comparator is the liver-MITOCHONDRIAL median, which is SLOWER:

  fibroblast proteome median        ~1.7 d   (human dermal, quiescent dSILAC)
  human fibroblast PCC (Chandler)    3.6 d   -> already in the SLOW tail of its own dist.
  whole mouse liver median         2.4-3.6 d (Rolfs/Price/Hasper)
  mouse liver MITOCHONDRIAL median   4.26 d  (Kim 2012) <- PCC's compartment, SLOWER
  human bulk-liver FSR anchor        ~4 d    (turnover time, in vivo biopsy)

So Chandler's 3.6 d sits in the slow tail of the fibroblast distribution and matches the
liver-mito median; the compartment-anchored best estimate for human hepatic PCC is
~3-6 d, central ~4.3 d -- i.e. AT or ABOVE the round-15 crossing (3.5452 d), not below.
The correction runs AWAY from the thresholds: the 21-day step is MORE likely free-for-all
than round 15 implied. But the band still straddles, so it does not settle the side.

Same wide log-normal band as round 16 (combined sd(ln)~0.53) so only the CENTRE moves
(from the r=1 bridge's 3.6 d to the mito-compartment anchor 4.3 d); also a tighter band
(sd 0.35) treating the 3-6 d range as ~+/-1 sigma.
"""
import numpy as np
from scipy.stats import norm

CROSS = 3.5452
D14, TAU0 = 6000.0, 14.0
THC_INT = 150.775            # Model C intrinsic theta_crit (round 15, verified)
CENTRAL = 4.3               # compartment-anchored central estimate (liver-mito context)
BAND = (3.0, 6.0)


def trough(tau, thalf):
    lam = np.log(2.0) / thalf
    x = np.exp(-lam * tau)
    return D14 * (tau / TAU0) * x / (1.0 - x)


def p_below(center, sd, n=400000, seed=0):
    rng = np.random.default_rng(seed)
    hep = np.exp(np.log(center) + rng.normal(0, sd, n))
    return float(np.mean(hep < CROSS))


def p_below_analytic(center, sd):
    return float(norm.cdf((np.log(CROSS) - np.log(center)) / sd))


def main():
    print("=" * 92)
    print("  ROUND 17. THE HEPATIC CORRECTION RUNS THE OTHER WAY (PCC is mitochondrial-matrix)")
    print("=" * 92)
    print(f"""
  Compartment anchors (half-lives): fibroblast proteome median ~1.7 d; human fibroblast
  PCC 3.6 d (slow tail); whole-liver median 2.4-3.6 d; liver-MITOCHONDRIAL median 4.26 d
  (Kim 2012, PCC's compartment); human bulk-liver FSR ~4 d. Best estimate for human
  hepatic PCC: {BAND[0]}-{BAND[1]} d, central ~{CENTRAL} d -- ABOVE the {CROSS} d crossing.
""")
    print(f"  P(human hepatic PCC half-life < {CROSS} d crossing):\n")
    print(f"      {'estimate':>34}  {'P(below) MC':>12}  {'P(below) analytic':>18}")
    for lab, c, sd in [("central 4.3 d, wide band sd(ln)=0.53", CENTRAL, 0.53),
                       ("central 4.3 d, tight band sd(ln)=0.35", CENTRAL, 0.35),
                       ("round-16 bridge centre r=1 -> 3.6 d (for contrast)", 3.6, 0.53)]:
        print(f"      {lab:>34}  {p_below(c, sd):>11.0%}  {p_below_analytic(c, sd):>17.0%}")
    print(f"\n  P(below) across the 3-6 d band endpoints (wide band sd=0.53):")
    print(f"      {'central':>10}  {'P(below)':>9}  {'21d step':>10}")
    for c in (3.0, 3.6, 4.3, 5.0, 6.0):
        t21 = trough(21, c)
        step = "free-for-all" if t21 > THC_INT else "needs theta"
        print(f"      {c:>8.1f}d  {p_below(c,0.53):>8.0%}  {step:>12} (Dtr21={t21:.0f} vs {THC_INT:.0f})")
    print(f"""
  READING: with the mitochondrial-compartment anchor the central estimate is ~4.3 d, above
  the crossing, so P(hepatic below the crossing) is ~30% (29-36% over the band widths), not
  ~50-65%. The 21-day safe-for-all step is therefore MORE likely to hold (~65-70%) than
  round 15's "toward every threshold" framing implied -- the directional claim there was
  wrong because it used the whole-liver median rather than the mitochondrial-matrix one, and
  PCC is a matrix enzyme. BUT the band still straddles: at the fast end of the plausible
  range (3-3.6 d) the 21-day step still needs theta, and P(below) is still ~36-49%. So the
  correction reverses the direction without settling the side -- the decision is now more
  likely free than round 15 said, but is still not resolved by the literature. The bridge
  narrows the question (compartment context + four converging anchors put it in 3-6 d) but
  does not settle it (no PCC-specific number in any supplementary table; within-tissue
  spread 2-4 orders; cross-tissue R^2 only 0.3-0.65), so a MEDIAN-based ratio is only weakly
  transferable to one protein. As in round 16, the only route that settles the side is to
  READ lambda from the trial record (round 12: 7-18 sigma).
""")


if __name__ == "__main__":
    main()
