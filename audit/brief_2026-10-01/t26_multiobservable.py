"""
The joint multi-observable FIM, built from the observables the REPOSITORY defines
(audit/PHASE10_SENSORS.md, audit/q15_sensors.py, audit/u3_joint_fisher.py). Repo
values only; gaps flagged.

THE TEN CANDIDATES (PHASE10_SENSORS.md step 1). Five are DROPPED because the model
carries only the marker state and has no propionyl-CoA pool, so every observable
that reads that pool is inexpressible:
  1  1-13C-propionate breath test   DROPPED (no propionyl-CoA state; S not a flux)
  2  ex vivo PCC activity           KEPT  A = kappa*theta*exp(eps)  -> theta DIRECTLY
  3  plasma 3-HP                     KEPT  = the marker m (the record itself)
  4  2-methylcitrate (2-MC)         DROPPED (needs propionyl-CoA x oxaloacetate)
  5  C3 (propionylcarnitine)        DROPPED (needs propionyl-CoA; free carnitine obs.)
  6  C3/C2 ratio                    DROPPED (needs propionyl-CoA vs acetyl-CoA)
  7  propionylglycine               DROPPED (needs propionyl-CoA; glycine obs.)
  8  IV 13C3-3-HP micro-dose        KEPT  CL_total = k_e + theta/(K_m+m0) -> k_e
  9  whole-body propionate turnover KEPT  R_a = L + b                 -> the LOAD
  10 timed urine + paired plasma    KEPT  CL_renal = f_renal*k_e      -> k_e (bound only)

So C3/propionylcarnitine (cand. 5) is NOT a re-read of the 3-HP marker -- it reads
the UPSTREAM propionyl-CoA pool, which would be INDEPENDENT information, but the
model cannot express it. The three expressible independent sensors each buy
precision on ONE existing parameter: ex vivo PCC -> theta, IV tracer -> k_e,
turnover -> L+b. REPO PRECISIONS (q15_sensors.py / PHASE10): ex vivo PCC 20%
(assay repeatability; nuisances kappa, tau_map NOT reported in any human study),
IV 13C3 tracer 5.4% (P.8's four micro-dose studies, alpha closure assumed free),
turnover 10% (a turnover study at the load record's CV). The breath test is
modelled anyway in u3 as a functional ratio S = theta/(theta+k_e(K_m+m)) with an
unknown recovery c carried FREE, swept at CV 5-10% (not a repo-measured figure).
"""
import numpy as np

import core as C
import dyn as D
import t16_flux_observable as T16

TH = C.TH_5


def sd_theta(sd_th=np.inf, sd_ke=C.KE_PRIOR, turnover=None, seeds=C.SEEDS):
    """sd(ln theta) at 5% monthly, with sensors as diagonal-precision buys
    (q15's framing). sd_th: ex vivo PCC on ln theta; sd_ke: external ln k_e;
    turnover: fractional precision on L+b, COMBINED with the record's 5% prior."""
    def one(s):
        r = C.record(TH, seed=s)
        L = r["L"]
        b_pr = C.LOAD_FRAC * L
        if turnover is not None:
            b_pr = 1.0 / np.sqrt(1.0 / (C.LOAD_FRAC * L) ** 2 + 1.0 / (turnover * L) ** 2)
        th_pr = None if not np.isfinite(sd_th) else sd_th
        ke_pr = 0.0 if sd_ke == 0 else sd_ke
        return C.sd_th(C.fim(r, noise="perdraw", ke_prior=ke_pr, b_prior=b_pr,
                             th_prior=th_pr))
    return C.over_seeds(one, seeds)[0]


def main():
    print("=" * 96)
    print("  THE JOINT MULTI-OBSERVABLE FIM (repo observables, repo precisions), 5% MONTHLY")
    print("=" * 96)

    base = sd_theta()
    print(f"\n  record (plasma 3-HP) alone: sd(ln theta) = {base:.3f}")
    print("\n  EACH INDEPENDENT SENSOR, AT ITS REPO PRECISION:")
    rows = [("ex vivo PCC activity (-> theta, 20%)", dict(sd_th=0.20)),
            ("ex vivo PCC at 10%", dict(sd_th=0.10)),
            ("IV 13C3 tracer (-> k_e, 5.4%)", dict(sd_ke=0.054)),
            ("IV 13C3 tracer, k_e exact", dict(sd_ke=1e-4)),
            ("turnover (-> L+b, 10%)", dict(turnover=0.10)),
            ("turnover, L exact", dict(turnover=1e-6))]
    for lab, kw in rows:
        print(f"    {lab:38s} {sd_theta(**kw):.3f}")

    print("\n  COMBINATIONS (does any sensor set reach 0.20 at 5% monthly?):")
    combos = [("k_e 5.4% + turnover 10% (nuisance sensors only)", dict(sd_ke=0.054, turnover=0.10)),
              ("k_e exact + turnover exact (nuisances, best case)", dict(sd_ke=1e-4, turnover=1e-6)),
              ("ex vivo PCC 20% + k_e 5.4%", dict(sd_th=0.20, sd_ke=0.054)),
              ("ex vivo PCC 20% + k_e 5.4% + turnover 10% (all three)", dict(sd_th=0.20, sd_ke=0.054, turnover=0.10))]
    for lab, kw in combos:
        v = sd_theta(**kw)
        print(f"    {lab:52s} {v:.3f}{'  <= GATE' if v <= 0.20 else ''}")

    print("\n  THE FLUX PAIR (breath fraction-oxidised S + turnover L+b), repo precisions:")
    for X, lab in ((1e-7, "L exact"), (0.10, "L via turnover (10%)")):
        yc = T16.crossing_in_Y(lambda Y: T16.sd_sL(Y, X))
        print(f"    if S were absolutely known: pair reaches 0.20 at S-precision <= "
              f"{100*yc:.1f}% ({lab})")
    print("    BUT the breath S has an UNKNOWN recovery c (carried free, u3), so a single")
    print("    reading is worth nothing and the pair needs many readings: u3 finds ~70")
    print("    breath readings at 10% capacity (CV 5%), more at 5%. Turnover cleanly")
    print("    supplies the L half (10% is ample); the breath recovery is the unresolved")
    print("    half -- exactly the calibration requirement of R11.")

    print("\n  ITEM 3 RECOMPUTED WITH TURNOVER'S OWN 10% PRECISION (not a hypothetical):")
    print(f"    cross-sectional theta, record + turnover 10%: {sd_theta(turnover=0.10):.3f}")
    print("      (turnover 10% is LOOSER than the record's own 5% load-bias prior, so it")
    print("       barely tightens b: combined ~4.5%. Measuring the load does ~nothing for")
    print("       cross-sectional theta, because k_e -- not the load -- is the binding")
    print("       nuisance. This is the q15 candidate-9 result.)")
    L = D.LTOT
    tt = np.arange(0, 2.0 + 1e-9, 5 / 60 / 24)
    exW = D.Experiment([(-1, 0.5 * L), (0, 0.0)], L0=0.5 * L)
    for l0f, lab in ((0.10, "turnover 10%"), (0.05, "5%"), (0.01, "1%")):
        F = D.fim(D.params(ex=exW), exW, tt, priors={2: 0.062, 3: l0f * L})
        print(f"    k_e from the safe transient, load known to {lab:>12}: "
              f"sd(ln k_e) = {np.sqrt(np.linalg.inv(F)[1, 1]):.3f}")
    print("      So turnover at its own 10% gives k_e only to ~0.29 from the transient --")
    print("      the ~0.06 of R-B needed the load to 1%, which turnover does not deliver.")

    print("\n  LONGITUDINAL CONTRAST UNDER SENSORS:")
    print("    The contrast already CANCELS the shared k_e, b and K_m, so the k_e tracer")
    print("    and the turnover add nothing to it (they tighten what already cancels). A")
    print("    per-visit ex vivo PCC at 20% gives a contrast of sqrt(2)*0.20 ~ 0.28, WORSE")
    print("    than the record's own 0.13. So no sensor improves monitoring: the record")
    print("    alone is the best instrument for CHANGE; sensors only help CLASSIFICATION.")

    print("\n  VERDICT (item 6):")
    print("    - Only the direct enzyme assay (ex vivo PCC) reaches 0.20 cross-sectionally")
    print("      at 5% monthly, at its repo 20% (0.200), alone. The nuisance sensors do")
    print("      NOT: k_e 5.4% + turnover 10% = 0.45; even k_e exact + L exact stays ~0.4")
    print("      monthly (it needs dense sampling too, q15).")
    print("    - So the report's cross-sectional claim needs ONE addition: besides the")
    print("      calibrated flux reading (~19%), a direct ex vivo PCC enzyme assay (repo-")
    print("      priced at 20%) also reaches target -- it measures theta itself. Both are")
    print("      DIRECT measurements, not inferences from the accumulation marker; the")
    print("      marker-only conclusion stands. The nuisance sensors do not overturn it.")
    print("    - The biggest untapped lever is NOT a new assay but a MODEL EXTENSION: a")
    print("      propionyl-CoA state would make the breath test and the C3/2-MC/")
    print("      propionylglycine markers expressible (5 of 10 candidates), the only")
    print("      observables that read the enzyme's substrate pool. That is outside the")
    print("      current model (PHASE10 1.1(b)).")


if __name__ == "__main__":
    main()
