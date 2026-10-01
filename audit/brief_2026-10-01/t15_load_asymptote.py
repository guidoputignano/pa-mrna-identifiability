"""
Item 13. Does the fall asymptote? The load-record term, bias against scatter.

The 'load record' part of the 0.39 is two different things:
    - a per-draw RECORD SCATTER (LOAD_FRAC = 5% of L in the flux-ordinate noise,
      core.sigmas), which averages down with n, and
    - a single systematic BIAS b (one unknown per patient, prior 5% of L),
      which the user expected not to average down.

This separates them, reports the n -> infinity limit with k_e exact, and costs
an independent measurement of L (which tightens the bias prior).
"""
import numpy as np

import core as C

L0 = C.record(C.TH_5, seed=0)["L"]
NS = (25, 100, 365, 730, 1460, 5000, 20000, 100000)


def sd(n, ske=1e-4, bX=None, lf=C.LOAD_FRAC, ac=C.ASSAY_CV, seeds=range(8)):
    kw = dict(noise="perdraw", ke_prior=ske, load_frac=lf, assay_cv=ac)
    if bX is not None:
        kw["b_prior"] = max(bX * L0, 1e-9)

    def one(s):
        r = C.record(C.TH_5, seed=s, times=np.linspace(0.0, 730.0, int(n)))
        return C.sd_th(C.fim(r, **kw))
    return C.over_seeds(one, seeds=seeds)[0]


def main():
    print("=" * 96)
    print("  ITEM 13.  THE LOAD-RECORD TERM: BIAS vs SCATTER, AND THE n -> inf LIMIT (k_e EXACT)")
    print("=" * 96)

    base = sd(25)
    print(f"\n  (1) DECOMPOSITION OF THE {base:.3f} AT n = 25 (k_e exact, K_m free):")
    only_assay = sd(25, lf=1e-6, bX=1e-7)       # no load at all
    no_scatter = sd(25, lf=1e-6)                # bias kept, scatter off
    no_bias = sd(25, bX=1e-7)                   # scatter kept, bias pinned (L exact)
    km_pin = sd(25, bX=1e-7, lf=1e-6)           # as only_assay
    print(f"      full (5% assay, 5% load scatter, 5% bias prior):  {base:.3f}")
    print(f"      remove the per-draw load SCATTER (LOAD_FRAC->0):   {no_scatter:.3f}"
          f"   (scatter worth {base - no_scatter:+.3f})")
    print(f"      remove the systematic BIAS b (measure L exactly):  {no_bias:.3f}"
          f"   (bias worth    {base - no_bias:+.3f})")
    print(f"      remove BOTH load terms (assay + K_m only):         {only_assay:.3f}")
    print(f"      -> at n = 25 the per-draw SCATTER costs more than the BIAS.")
    print(f"      -> the user's premise (bias is the binding term) does not hold here.")

    print(f"\n  (2) n -> inf WITH k_e EXACT. Does sd asymptote above 0.20?")
    print(f"      {'n':>8}  {'b 5% prior':>11}  {'b pinned (L exact)':>18}  {'sd*sqrt(n), b 5%':>18}")
    for n in NS:
        a = sd(n)
        b = sd(n, bX=1e-7)
        print(f"      {n:>8}  {a:>11.4f}  {b:>18.4f}  {a*np.sqrt(n):>18.3f}")
    print(f"      sd*sqrt(n) settles to a constant: the fall is ~1/sqrt(n), limit 0.")
    print(f"      NEITHER the bias nor the scatter creates a positive floor: both are")
    print(f"      learned from the dose-train modulation of the operating point m.")
    print(f"      => 0.151 at 730 draws is NOT flattening toward a floor. Item 2 is")
    print(f"         not answered by an asymptote; it is answered by feasibility (item 12).")

    print(f"\n  (3) THE SAME WITH k_e NOT EXACT (the realistic case).")
    print(f"      {'n':>8}  {'s_ke=0.10':>10}  {'s_ke=0.20':>10}  {'s_ke=0.40':>10}")
    for n in (25, 365, 2000, 20000, 100000):
        print(f"      {n:>8}  " + "  ".join(f"{sd(n, ske=s):>10.4f}" for s in (0.1, 0.2, 0.4)))
    print(f"      Still -> 0, but only once the dose train has identified k_e, which")
    print(f"      needs ~1e4-1e5 draws. At any clinical n the k_e prior sets the level")
    print(f"      (item 12). So there is no floor, but no feasible escape either.")

    print(f"\n  (4) COSTING AN INDEPENDENT MEASUREMENT OF L (tighten the bias prior to X%).")
    print(f"      {'X% on L':>10}  {'n=25, k_e exact':>15}  {'n=25, k_e to 0.2':>16}  {'n=365, k_e exact':>16}")
    for X, lab in ((0.05, "5% (prior)"), (0.02, "2%"), (0.01, "1%"), (1e-7, "exact")):
        print(f"      {lab:>10}  {sd(25, bX=X):>15.3f}  {sd(25, ske=0.2, bX=X):>16.3f}"
              f"  {sd(365, bX=X):>16.3f}")
    print(f"      Measuring L helps the finite-n load term (0.38 -> 0.28 at n=25, k_e")
    print(f"      exact) but does not move the n -> inf limit, and does almost nothing")
    print(f"      when k_e is the binding nuisance (k_e to 0.2 column).")

    print(f"\n  (5) WHERE L BINDS, AND WHAT MEASURING IT TO ~1-2% BUYS (three places):")
    print(f"      (i) the 0.39 partition  : 0.38 -> 0.30 at n=25 (this script, part 4).")
    print(f"      (ii) the washout intercept L0 (t7): the single biggest rung on the")
    print(f"           knowledge ladder -- withdraw-half goes 9.4 -> 1.5 -> 0.63 as L0")
    print(f"           is pinned to 5% then 1% (t7_design.txt). But the washout route")
    print(f"           is ruled out on safety (item 15), so this lever is unavailable.")
    print(f"      (iii) the flux route (item 11): the flux reading needs ABSOLUTE")
    print(f"           calibration; measuring L to X% is one way to supply it via the")
    print(f"           (s, L) pair. See t16.")


if __name__ == "__main__":
    main()
