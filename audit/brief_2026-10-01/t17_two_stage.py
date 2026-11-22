"""
Item 15. The safety threshold as a hard constraint, and the two-stage design.

(1) The crisis threshold, stated and sourced, and flagged as a MODELLING
    construct rather than a clinical measurement.
(2) The excursion sweep for theta with the admissible region marked (figure data):
    every load-up excursion that reaches 0.20 sits past the threshold.
(3) The smallest SAFE excursion that pins k_e to ~0.2 log units. lambda ~ k_e
    once the enzyme is working hard, so a mild transient fixes k_e where it
    cannot fix theta -- but only as well as the unrecorded load L0 is known.
(4) The two-stage cost: a safe k_e transient, then n steady-state records for
    theta, at n = 24, 52, 104, 365. Does it reach 0.20?
"""
import numpy as np

import core as C
import dyn as D
import t7_design as T7
import t14_ke_sampling as T14

L = D.LTOT
KM_REC = T7.KM_REC
FOLD = 1.177                      # results/numbers/boundary.tex, printed by t0
FOLD_LOAD = FOLD * L
FOLD_M = 276.0                    # boundary.tex fold_m, the marker at the saddle


def marg(ex, tt, idx, priors):
    F = D.fim(D.params(ex=ex), ex, tt, priors=priors)
    w = np.linalg.eigvalsh(F)
    return np.inf if w[0] <= 1e-12 * w[-1] else float(np.sqrt(np.linalg.inv(F)[idx, idx]))


def two_stage_theta(ske, n):
    """n draws over two years for theta, k_e known to ske from stage 1."""
    def one(s):
        return C.sd_th(C.fim(T14.record_n(n, s), noise="perdraw", ke_prior=ske))
    return C.over_seeds(one)[0]


def main():
    print("=" * 96)
    print("  ITEM 15.  SAFETY AS A HARD CONSTRAINT, AND THE TWO-STAGE DESIGN")
    print("=" * 96)

    print("\n  (1) THE CRISIS THRESHOLD.")
    print(f"      results/numbers/boundary.tex: fold_load_ratio = {FOLD}, i.e. the")
    print(f"      parent model (Model C / ParamsC) loses its steady state when the")
    print(f"      load exceeds {FOLD:.3f} x baseline = {FOLD_LOAD:.0f}. (That is a LOAD")
    print(f"      RATIO, which is model-independent; the fold MARKER {FOLD_M:.0f} is Model")
    print(f"      C's own, not Model A's -- at the fold load Model A's marker is only")
    print(f"      ~64. So the admissibility bound used below is f <= {FOLD:.3f} x load.)")
    print(f"      PROVENANCE AND STATUS: src/model_A/run.py:35-40 states the ILLNESS")
    print(f"      ratio (~1.5x baseline load) as 'the physiological claim', and that")
    print(f"      Model C's FOLD 'sits at 1.18x' -- a saddle-node of ParamsC's own R0")
    print(f"      setting. So 1.177 is a MODELLING construct (where this model has no")
    print(f"      steady state), not a measured clinical crisis level. The repository's")
    print(f"      one clinical-labelled figure is M_TOL = 40 (a 'tolerable metabolite")
    print(f"      level', run.py:34); note the untreated 5% marker is already 54 > 40,")
    print(f"      so M_TOL is a target, not a ceiling. The hard bound used here is the")
    print(f"      substrate one: an ADDED-substrate excursion must keep the load below")
    print(f"      {FOLD:.2f}x. Substrate-reducing excursions (withdrawal) and the drug")
    print(f"      channel (no substrate) are admissible by construction.")

    print("\n  (2) THE EXCURSION SWEEP FOR THETA, WITH ADMISSIBILITY (figure data).")
    print(f"      load multiple f -> untreated marker, admissible?, sd(ln theta)")
    print(f"      (dense 5-min transient to the new level for 48 h, joined to the")
    print(f"      monthly record, with the k_e prior + K_m from the record + L0 to 1%):")
    print(f"      {'f x load':>9}  {'marker':>7}  {'admissible':>10}  {'sd(ln theta)':>12}")
    sweep = []
    tt = np.arange(0, 2.0 + 1e-9, 5 / 60 / 24)
    for f in (0.5, 0.75, 1.0, 1.177, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0):
        m_f = float(C.m_ss(C.TH_5, C.KM0, C.KE, f * L))
        adm = f <= FOLD + 1e-9
        if f <= 1.0:   # withdrawal / baseline: step DOWN, m drops
            ex = D.Experiment([(-1, 0.5 * L), (0, (f - 0.5) * L)], L0=0.5 * L)
        else:          # load up, held
            ex = D.Experiment([(-1, 0.5 * L), (0, 0.5 * L + (f - 1) * L)], L0=0.5 * L)
        F = D.fim(D.params(ex=ex), ex, tt)
        sdth = T7.joint_record(F, {"ke": 0.8, "Km": KM_REC, "L0": 0.01 * L})
        sweep.append((f, m_f, adm, sdth))
        print(f"      {f:>9.3f}  {m_f:>7.1f}  {'yes' if adm else 'NO':>10}  {sdth:>12.3f}")
    reach = [f for f, m, a, s in sweep if s <= 0.20]
    adm_reach = [f for f, m, a, s in sweep if s <= 0.20 and a]
    print(f"      load multiples reaching 0.20: {reach}; of these ADMISSIBLE: "
          f"{adm_reach if adm_reach else 'NONE'}")
    print(f"      -> every load-up excursion that reaches 0.20 needs the marker far")
    print(f"         past the fold {FOLD_M:.0f}; the admissible band never reaches 0.20.")

    print("\n  (3) THE SMALLEST SAFE EXCURSION THAT PINS k_e.")
    print(f"      A withdrawal (marker DROPS, always safe), dense 5-min for 48 h,")
    print(f"      K_m from the record, L0 known to the stated precision:")
    exW = D.Experiment([(-1, 0.5 * L), (0, 0.0)], L0=0.5 * L)   # full diet withdrawal
    mW, _ = D.simulate(D.params(ex=exW), exW, tt, sens=False)
    print(f"      withdrawal: marker {mW.min():.1f}-{mW.max():.1f} (drops, safe),")
    for l0, lab in ((0.05 * L, "L0 to 5%"), (0.01 * L, "L0 to 1%"), (1e-6, "L0 exact")):
        print(f"        {lab:>9}: sd(ln k_e) = {marg(exW, tt, 1, {2: KM_REC, 3: l0}):.4f}, "
              f"sd(ln theta) from the same transient = "
              f"{marg(exW, tt, 0, {2: KM_REC, 3: l0}):.3f}")
    print(f"      Even a MILD withdrawal suffices for k_e: 25% of the diet (marker")
    mild = D.Experiment([(-1, 0.5 * L), (0, 0.375 * L)], L0=0.5 * L)
    mm, _ = D.simulate(D.params(ex=mild), mild, tt, sens=False)
    print(f"        {mm.min():.1f}-{mm.max():.1f}) gives sd(ln k_e) = "
          f"{marg(mild, tt, 1, {2: KM_REC, 3: 0.05 * L}):.3f} at L0 5%.")
    print(f"      So the user's hypothesis holds: a safe transient pins k_e (~0.2 log")
    print(f"      units) where it is nowhere near pinning theta (~1.2). But k_e's")
    print(f"      precision is capped by the UNRECORDED load L0: ~0.2 at L0 5%, ~0.06")
    print(f"      at L0 1%. Tightening k_e past 0.06 needs L0 measured -- item 13.")

    print("\n  (4) THE TWO-STAGE COST: safe k_e transient, then n records for theta.")
    print(f"      sd(ln theta) at 5%, median over seeds, k_e fixed at the stage-1 value:")
    print(f"      {'draws (2 yr)':>12}  {'k_e~0.20 (L0 5%)':>16}  {'k_e~0.06 (L0 1%)':>16}"
          f"  {'k_e exact (ref)':>15}")
    for n in (24, 52, 104, 365):
        a = two_stage_theta(0.20, n)
        b = two_stage_theta(0.06, n)
        c = two_stage_theta(1e-4, n)
        print(f"      {n:>12}  {a:>16.3f}  {b:>16.3f}  {c:>15.3f}")
    print(f"      VERDICT: the two-stage design is SAFE and feasible, but it does NOT")
    print(f"      reach 0.20. With a realistic L0 (5%) it lands near 0.9 at every n;")
    print(f"      with L0 measured to 1% it reaches ~0.25-0.30 only with hundreds of")
    print(f"      draws, still short of 0.20. k_e to ~0.2 is simply not tight enough")
    print(f"      (item 12), and the safe transient cannot do better without also")
    print(f"      solving the load-measurement problem.")


if __name__ == "__main__":
    main()
