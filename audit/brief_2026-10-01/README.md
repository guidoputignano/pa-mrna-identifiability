# Brief of 1 October 2026: practical identifiability of enzyme capacity

The findings are in **[`REPORT.md`](REPORT.md)**. This directory holds the code behind them.

It builds on the repository rather than beside it. The footing (convention A, the monthly record, σ, the priors, the dose train and the 24 intake seeds) is imported from `audit/n4_theta_scope.py`, and `core.fim` reproduces `n4.sd_ln_theta` to 1e-13. Nothing here writes into `results/` or `paper/`, so the artefact guard is untouched: it globs `audit/*.py`, not subdirectories.

| script | brief task | runtime |
|---|---|---|
| `t0_footing.py` | the repository constants and files the report quotes, printed from source | 1 s |
| `core.py` | the repository's FIM, the ln-m form, profiling | — |
| `t1_identifiability.py` | 1: U-curve, correlations, elasticity, profile | 15 s |
| `t2_cohort.py` | 2: joint and hierarchical cohort fit, with Monte Carlo | 2 min |
| `t3_variance.py` | 3: the 1/n versus n-independent split, and 15× sampling | 5 s |
| `t4_set.py` | 4: set identification and the one-sided bound | 1 min |
| `t5_value.py` | 5: what each extra measurement is worth | 5 s |
| `dyn.py` | the ODE with forward sensitivities, for tasks 6 and 7 | — |
| `t6_estimators.py` | 6: FD, integral, NLS, EKF and iterated EKF smoother | ~15 min |
| `t6b_ma1.py` | 6: the MA(1) claim, on the equation error | 30 s |
| `t7_design.py` | 7: washout and drug-forcing design | 1 min |
| `t8_floor.py` | follow-up 1: decomposing the 0.39, the n-dependence, the k_e share curve | 15 s |
| `t9_omega.py` | follow-up 2: Task 2 at N = 24, and whether 24 records could identify ω | 8 min |
| `t10_factorial.py` | follow-up 3: the TB-A to TB-B factorial, and the 5% design frontier | 8 min |
| `t6c_iv.py` | follow-up 4: IV, bias-compensated LS and generalised TLS for the FD equation | 6 min |
| `t6d_ekf.py` | follow-up 5: the one-pass EKF diagnosed | 6 min |
| `t11_fim_check.py` | follow-up 6: FIM against Monte Carlo, and the Task 5 spot-checks | 9 min |
| `t12_flux_mc.py` | follow-up 6: the Monte Carlo crossing for the direct flux reading (no recheck) | 8 min |
| `t13_conditioning.py` | round 3, item 9: what the likelihood conditions on (per-draw L and dose) | 2 s |
| `t19_residuals.py` | round 3, item 10: residual diagnostics as a power analysis of misspecification | 6 min |
| `t16_flux_observable.py` | round 3, item 11: the flux observable O=c·θ·g(m) as a specification; (s,L) reparam | 2 min |
| `t14_ke_sampling.py` | round 3, item 12: sd(ln θ) over (sd(ln k_e), draws/year) | 2 min |
| `t15_load_asymptote.py` | round 3, item 13: the load-bias asymptote and the value of measuring L | 2 min |
| `t17_two_stage.py` | round 3, item 15: the safety threshold and the two-stage design | 3 min |
| `t18_set_id.py` | round 3, item 16: set identification elevated — what is recoverable | 10 s |
| `t20_saturation.py` | round-3 close-out A: saturation regime, integral collapse, excursion = saturation | 1 min |
| `t21_sensitivity.py` | round-3 close-out B: sensitivity of the conclusions, flux contour across capacity | 3 min |
| `t22_longitudinal.py` | round-3 add-on: sampling at realistic k_e, and longitudinal (monitoring) precision | 2 min |
| `t23_antibiotic.py` | round-3 add-on: antibiotics as an unmodelled load term, and as a safe downward excursion | 3 min |
| `t24_contrast_stress.py` | round-3 close: stress-testing the longitudinal contrast (drift, antibiotics, power, minimum programme) | 2 min |
| `t25_drug_cv_marker.py` | round-3 final: drug response, assay-CV x draws surface, and the falsifiable marker-variability prediction | 2 min |
| `t26_multiobservable.py` | multi-observable: joint FIM over the repo's sensor panel (ex vivo PCC, IV tracer, turnover, breath) | 2 min |
| `t27_exvivo_kappa.py` | the ex vivo PCC calibration (kappa) contour, and scoping the propionyl-CoA model extension | 1 min |
| `t28_ratio_route.py` | the ratio route (capacity as a fraction of normal cancels kappa; sd over healthy-control N), and the kappa-difference tolerance | 1 min |
| `t29_mapping_scatter.py` | round 6, item 1: the ratio route fails for whole-body capacity — per-patient in-vitro/in-vivo mapping scatter (Manoli R²) does not cancel | 1 min |
| `t30_monitoring_confounders.py` | round 6, item 2: the monitoring contrast vs the literature confounders (eGFR→k_e, carnitine→calibration); the empirical challenge and falsification spec | 1 min |
| `t31_breath_contrast.py` | round 6, item 3: the breath test as the monitoring instrument (reads in-vivo capacity directly, c cancels; one reading ≈ 104 draws) | 1 min |
| `t32_massbalance.py` | round 6, item 4: the clearance split (84% non-PCC at 5%) vs the PA literature; Model A's unbounded k_e·m and Model B's saturating fix | 1 min |
| `t33_model_b.py` | round 7, item 1: Model A vs Model B side by side — steady-state existence (θ_crit), enzyme share, cross-sectional sd, longitudinal contrast | 1 min |
| `t34_breath_stress.py` | round 7, item 3: breath-test stress tests — 8% provenance, c-drift breakpoint, control spread, confounder inheritance, cross-sectional c-contour | 1 min |
| `r3_figs.py` | round 3 figures (flux contour, k_e×sampling, asymptote, excursion/safety) | 2 min |
| `*_recheck3.py` | independent verifier re-derivations of the round-3 scripts | — |
| `*_recheck.py` | an independent verifier's re-derivation of each follow-up script, sharing no code with it | 3–25 min |
| `plotstyle.py` | the repository's rcParams with a validated palette | — |

To run them all, from this directory:

```bash
pip install -r ../../requirements.txt
for s in t0_footing t1_identifiability t2_cohort t3_variance t4_set t5_value t7_design t6_estimators t6b_ma1 \
         t8_floor t9_omega t10_factorial t6c_iv t6d_ekf t11_fim_check t12_flux_mc \
         t13_conditioning t14_ke_sampling t15_load_asymptote t16_flux_observable t17_two_stage \
         t18_set_id t19_residuals t20_saturation t21_sensitivity t22_longitudinal t23_antibiotic t24_contrast_stress t25_drug_cv_marker t26_multiobservable t27_exvivo_kappa t28_ratio_route \
         t29_mapping_scatter t30_monitoring_confounders t31_breath_contrast t32_massbalance \
         t33_model_b t34_breath_stress; do
  python -W error $s.py > out/$s.txt
done
python -W error r3_figs.py   # round-3 figures
```

Figures go to `fig/` and printed output to `out/`. Every number in `REPORT.md` appears in one of the `out/*.txt` files. `out/q15_sensors.txt` and `out/u3_joint_fisher.txt` are the repository's own Phase-10 sensor scripts (`audit/q15_sensors.py`, `audit/u3_joint_fisher.py`), run unchanged. `out/n6_profile_vs_crb.txt` is the repository's own `audit/n6_profile_vs_crb.py`, run unchanged (`python -W error audit/n6_profile_vs_crb.py` from the repository root), because the report quotes its 708 = 24%.
