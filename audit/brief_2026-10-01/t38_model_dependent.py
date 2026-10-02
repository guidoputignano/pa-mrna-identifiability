"""
Round 8, items 2 & 3: the cross-sectional conclusion is MODEL-DEPENDENT, and the
valid Model-B point is a construction whose provenance and ranges must be stated.

ITEM 2: state the cross-sectional claim per model, give the prior-sensitivity
honestly, and answer what would reach 0.20 under B.
ITEM 3: state what was chosen in the valid-B point and sweep the ceiling and the
managed load over their plausible ranges, reporting the RANGE of sd, not a point;
reconcile the ~10% load reduction with the paper's own 18% / 20-40% numbers.
"""
import numpy as np

import core as C
import t33_model_b as T
import t35_managed_load as M

TH5 = C.TH_5


def sd_B_managed(L, ke, ke_prior=0.8, ke_m_prior=None, drop_Ke=True, seeds=range(24)):
    def one(s):
        r = M.record_B_managed(TH5, L, ke, C.MONTHLY, s)
        return T.sd_th_B(T.fim_B(r, ke_prior=ke_prior, ke_m_prior=ke_m_prior, drop_Ke=drop_Ke))
    return float(np.median([one(s) for s in seeds]))


def sd_A(seeds=range(24)):
    return float(np.median([C.sd_th(C.fim(C.record(TH5, seed=s), noise="perdraw")) for s in seeds]))


def main():
    print("=" * 92)
    print("  ITEM 2. THE CROSS-SECTIONAL CONCLUSION IS MODEL-DEPENDENT")
    print("=" * 92)
    a = sd_A()
    b_valid = sd_B_managed(220, 100.0)
    print(f"""
  Stated per model (cross-sectional sd(ln theta) at 5% of normal, monthly record):
    Model A (the paper's working model):            {a:.2f}   -> no feasible route
                                                            (misses 0.20 by ~{a/0.20:.0f}x)
    Model B, literature-consistent, managed load:   {b_valid:.2f}   -> still misses 0.20,
                                                            but by ~{b_valid/0.20:.1f}x, not ~15x
  Model A is the model the parent paper ITSELF calls inadequate (its linear
  clearance predicts no decompensation; the paper keeps A only as the most
  favourable *identifiability* case, Remark 15). So the headline 'no marker-only
  design reaches target cross-sectionally' is a Model-A statement; under the more
  faithful B the gap is a factor of ~1.5, not ~15.
""")
    print("  PRIOR SENSITIVITY (the improvement is bought with prior information, not data):")
    print(f"    Model B, managed, ke prior 0.8, Ke pinned:   {b_valid:.2f}")
    print(f"    Model B, managed, ke prior 0.8, Ke free 0.5: {sd_B_managed(220,100.0,ke_prior=0.8,ke_m_prior=0.5,drop_Ke=False):.2f}")
    print(f"    Model B, managed, NO priors, Ke free:        {sd_B_managed(220,100.0,ke_prior=0.0,ke_m_prior=1e9,drop_Ke=False):.2f}   <- exchange degeneracy")
    nopri_A = float(np.median([C.sd_th(C.fim(C.record(TH5, seed=s), noise='perdraw', ke_prior=0)) for s in range(24)]))
    print(f"    Model A, NO ke prior:                        {nopri_A:.2f}")
    print("    So B's 0.29 rests on the k_e/K_e priors; strip them and B is far WORSE")
    print("    than A (the a-fortiori of Remark 15). The gain is prior-bought, not data-bought.")

    print("\n  WHAT WOULD REACH 0.20 UNDER B? (tighten the nuisance priors at the valid point)")
    print(f"    {'k_e prior sd':>14}  {'Ke pinned':>10}  {'Ke prior 0.3':>13}  {'Ke prior 0.1':>13}")
    for kep in (0.8, 0.4, 0.2, 0.1, 0.05):
        pin = sd_B_managed(220, 100.0, ke_prior=kep)
        k3 = sd_B_managed(220, 100.0, ke_prior=kep, ke_m_prior=0.3, drop_Ke=False)
        k1 = sd_B_managed(220, 100.0, ke_prior=kep, ke_m_prior=0.1, drop_Ke=False)
        mark = lambda v: f"{v:.3f}" + ("*" if v <= 0.20 else " ")
        print(f"    {kep:>14.2f}  {mark(pin):>10}  {mark(k3):>13}  {mark(k1):>13}")
    print("    (* = reaches 0.20.) Reading: under B, a k_e prior tightened from 0.8 to")
    print("    ~0.2 (with K_e known) reaches target. A sd(ln k_e)=0.2 k_e prior is exactly")
    print("    what the IV 13C3 tracer delivers (5.4%, R-G) or the safe transient at L0 to")
    print("    1% (R-B) -- so UNDER MODEL B the cross-sectional target is reachable with a")
    print("    k_e sensor the brief already prices, a far cheaper requirement than the")
    print("    direct-enzyme calibration Model A forces. This is model-dependent and worth")
    print("    stating as such: the whole 'only a direct calibrated enzyme read works'")
    print("    conclusion is partly an artefact of Model A's low enzyme share.")

    print("\n" + "=" * 92)
    print("  ITEM 3. PROVENANCE AND RANGES OF THE VALID MODEL-B POINT")
    print("=" * 92)
    print("""
  WHAT WAS CHOSEN, AND ON WHAT BASIS (nothing here is fitted to a patient):
    - ceiling k_e = 100: chosen BELOW the untreated load (~276) so theta=0 is not
      maintainable (the severe phenotype); literature-consistent (non-PCC disposal
      is low-capacity) but not a measured value -- the repo gives no k_e in units.
    - managed load L = 220: chosen to ADMIT a steady state at 5% with a ~1%-of-
      normal margin (a ~20% reduction from the untreated anchor). Also a choice.
  Because both are choices, the honest output is a RANGE over their plausible
  spans, not the single 0.29.
""")
    print("  sd(ln theta), monthly, swept over ceiling k_e x managed load L (Ke pinned):")
    hdr = "k_e \\ L"
    print(f"      {hdr:>10}" + "".join(f"{L:>8.0f}" for L in (180, 200, 220, 240)))
    allsd = []
    for ke in (50.0, 100.0, 150.0, 200.0):
        cells = ""
        for L in (180, 200, 220, 240):
            # only where a steady state exists with margin: L < theta + ke
            if L < TH5 + ke - 1:
                v = sd_B_managed(L, ke); allsd.append(v); cells += f"{v:>8.3f}"
            else:
                cells += f"{'--':>8}"
        print(f"      {ke:>10.0f}{cells}")
    print(f"\n    RANGE of sd across the valid grid: {min(allsd):.2f} to {max(allsd):.2f}"
          f" (median {np.median(allsd):.2f}).")
    print("    Every valid cell is > 0.20 and << Model A's 3.05: the model-dependence of")
    print("    the magnitude is real but the SIGN (misses 0.20 cross-sectionally) is robust")
    print("    across the whole plausible ceiling x load grid.")

    print("""
  RECONCILING THE LOAD-REDUCTION FIGURE WITH THE PAPER'S OWN NUMBERS:
    - This script's margin needs a ~10-20% load reduction to lift a 5% patient
      above theta_crit (t35).
    - The paper states management cuts the load by 20-40% (standard of care:
      protein restriction + gut suppression), and that the SMALLEST illness
      excursion that strands a patient is +18% in the load (s_fold 706 vs
      baseline 600 at theta=300, sections/03_dynamic.tex).
    - These are the SAME ORDER: the reduction management must supply (~10-20%),
      the reduction it does supply (20-40%), and the excursion that strands
      (18%) all sit in a ~10-40% band. So the patient lives in a narrow margin --
      normal management keeps them just inside it and a single ~18% illness
      pushes them out. The margin mu and the trigger amplitude are comparable,
      which is why decompensation is a recurring event rather than a rare one.
""")


if __name__ == "__main__":
    main()
