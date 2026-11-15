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
| `t35_managed_load.py` | round 8, items 2–3: θ_crit moves with the managed load (margin μ); the only valid Model B column (low ceiling + managed load); ParamsB miscalibration | 1 min |
| `t36_reconcile_remark15.py` | round 8, items 1,4: reconcile "B better for θ" with the paper's Remark 15 (marginal CRB vs joint conditioning + exchange degeneracy); share mechanism; corr(θ,k_e) vs share | 1 min |
| `contrast.py` | round 8: the sharing contract in one place — `window_contrast` ties every patient constant between windows unless given a drift; `check_contract` fails if any is left independent | — |
| `t37_contrast_audit.py` | round 8, item 1: the contract test + an audit of every FIM contrast number in the report through the helper (R-D clean; only the Model-B 0.167→0.088 ever changed) | 1 min |
| `t38_model_dependent.py` | round 8, items 2–3: cross-sectional sd per model, prior-sensitivity, what reaches 0.20 under B; valid-point provenance and ceiling×load sweep | 2 min |
| `t39_affine_check.py` | round 8, item 4: the paper's own affine condition (eq:affine) run on the repo's data — passes in the dosed record, fails over the excursion to the resting balance | 1 min |
| `t40_kepoint_robustness.py` | round 8 close: the ~0.19 under B swept across the ceiling×load grid under tightened priors, and sd(ln θ) vs sd(ln K_e) — the gate is a low-ceiling point needing K_e known to ~30% | 1 min |
| `t41_model_c.py` | round 9, item 2: the two-observable (marker + free carnitine) FIM under Model C — does an observed competing channel identify θ? Only with its shape constants (K_c, R₀, M_i, ρ₀) known | 2 min |
| `t42_serial_cv.py` | round 10, item 1: R-D against the only empirical serial-CV anchor (C3 CVi 26.4%) — contrast and draws/window for 0.20 across 8.2-61.1% | 1 min |
| `t43_plasma_pool.py` | round 10, item 2: the plasma-vs-pool compartment defect in Model C — free carnitine never earns its keep; the partition is a fifth unmeasured constant | 3 min |
| `t44_interval_decision.py` | round 11: the dosing-INTERVAL decision as a trough condition θ_base+Δ_trough>θ_crit — the interval tolerance on θ, compared to the record; the coarse grid (and all of Model C) is determined with no new measurement; the decision turns on λ, which is unmeasured | 1 min |
| `t45_on_therapy.py` | round 12: retrospective interval adjustment on a two-year on-therapy record — what the dated dose schedule buys for θ (sampling-limited, premise correction: 3.05 is already on-therapy), whether λ self-calibrates from the same record (it does, ~5.5%), and the one-sided lengthening test (15 d free, 21 d evidence-gated) | 1 min |
| `t46_timing.py` | round 13: WHEN — timing of each step. 21-day eligibility by NLS simulation (not FIM: the margin is tail-dominated, §7b skew) — it does NOT accrue with record length, superseding t45's optimism; λ crosses its spec early (FIM); the 15-day step needs no record; no review unlocks 21 d (needs the active protocol) | 3 min (diagnostic) |
| `t47_fim_audit.py` | round 14: the FIM-vs-simulation audit. States the method lesson (log-scale sd reliable, natural-scale/one-sided decision tail-dominated); spot-checks by NLS the two at-risk claims (round-11 determined-subgroup 61%→54%, cross-sectional at 5%) — both survive; scans the posterior correlations for hidden favourable cancellations (θ–k_e, θ–λ both already simulation-handled; no third) | 2 min |
| `t48_lambda_literature.py` | round 15: λ re-grounded on literature (human 3.6 d Chandler 1985; mouse ~7 d Attarwala). Which λ the trough sees (2-comp → effective 3–5 d); the decision at 3.6/3.8/7 d and the fail half-lives (14d 2.62, 15d 2.75, coarse 3.32, 21d 3.55); every headline flips inside the 3-fold band (λ the dominant uncertainty); the record discriminates 3.6 vs 7 d at 7–18 σ (a test, not self-calibration) | 1 min |
| `t49_hepatic_bridge.py` | round 16: can human hepatic PCC turnover be bridged from the fibroblast value via the mouse tissue ratio? No clean paired dataset exists; the three defensible ratio readings disagree on direction, so P(hepatic < 3.5452 d crossing) spans ~24–65 %, centred ~50 % — the bridge cannot settle the side, so the trial record is the only route | 1 min |
| `t50_hepatic_mito.py` | round 17: the hepatic correction runs the OTHER way — PCC is mitochondrial-matrix, so the liver-mito median (4.26 d) is the comparator, not whole-liver; central estimate ~4.3 d, ABOVE the crossing. P(hepatic < crossing) ~30 % (not ~50 %), so the 21-day step is ~65–70 % free-for-all, though the band straddles. Also records that the human λ (Gao model) is not public — the on-therapy record is the only accessible route | 1 min |
| `t51_pencil_geometry.py` | round 18: the PI's line-pencil geometry (each steady-state draw is a line g(m)θ+m·k_e=L in the (θ,k_e) plane) and what it adds — the near-concurrence region IS the fixed-K_m FIM ellipse, the −0.99 IS near-parallelism, Prop 1's affine degeneracy IS the narrow-band slope collapse, the K_m family IS Remark 2's fibre; plus (2) a model-free uncertainty from pairwise-intersection scatter (heavy-tailed, a third route to the round-13/14 tail), (3) a concurrency χ² adequacy test (calibrated, power vs a biased load, blind to a second channel), (4) Remark 6's TLS weighting reproduces the FIM/MLE | 1 min |
| `t52_knife_resolution.py` | round 19: can the on-therapy record settle which side of the 3.5452 d crossing a patient is on? Per patient, only if the true half-life is clear of it (a 4.3 d truth: an 80 % chance of ≥2σ within a year of monthly or 6 months of weekly draws, against the K_c = 100 line); a 3.6 d truth is never placed (≤0.6σ even after 5 years weekly). Pooling: 23–77 patients at 50 % power, and only if λ is shared. Corrects R-O/R-P's "a fortiori". Re-evaluates t46's FIM sd(ln λ) at each true half-life (reproduces t46 exactly at 3.5 d) | 1 s |
| `t53_crossing_shape.py` | round 19: the 21-day crossing under Model C's unmeasured shape constants (K_c, M_i, ρ₀; R0 re-derived as the repo does). Over the box corners the crossing moves about 3.29–4.08 d (3.30–3.97 d over 300 sampled sets); 3.6 d clears it in 80 % of sampled sets and at 6 of 8 corners (one at a time it fails at K_c ≤ 50 or M_i = 150); 4.3 d clears everywhere; 14/15-day crossings ≤ 3.06/3.21 d. Also shows the repo's 12-draw collapse_check understates the joint θ_crit spread (factor 1.4 vs 2.4 at the corners) | 1 min |
| `t54_dose_timing.py` | round 20: how late a fortnightly dose can be (trough condition, and the two-state parent model simulated through a late or skipped dose), and whether an estimate of theta decides a longer interval per patient; precision recomputed at each capacity (the uncertainty is about +-64 uM/d at every capacity with both tracers, so 0.44 holds at 5 % only) | 2 min |
| `pi_note/pi_note.tex`, `pi_note/pi_figs.py` | technical note for the PI (LaTeX, 8 pages): the four results at 5 % of normal and the dosing result as a specification on the PCC half-life and the delivered increment; `pi_figs.py` draws its four figures at print size (Figure 1 recomputed with t1's own functions). Build: `cd pi_note && python -W error pi_figs.py && pdflatex pi_note.tex` (twice) | 1 min |
| `story_figs.py` | round 19: the two plain-language charts for `story_so_far.html` (every route to capacity on one scale; how long to wait) | 10 s |
| `r18_figs.py` | round 18 figures: (1) the pencil in (θ,k_e) — narrow (near-parallel smear) vs wide excursion (concurs), on shared axes, lines coloured by m_i, with the slope-identity inset; (2) the K_m fibre; (3) the knife-edge — safety margin μ=Δ_trough(21)−θ_crit vs PCC protein half-life, crossing zero at 3.5452 d, assumed 3.5 d (−8 µM/d) / measured 3.6 d (+10 µM/d) / 3–6 d band; 1.5 % is a half-life gap (the report's lead figure) | 1 min |
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
         t33_model_b t34_breath_stress t35_managed_load t36_reconcile_remark15 \
         t37_contrast_audit t38_model_dependent t39_affine_check t40_kepoint_robustness t41_model_c t42_serial_cv t43_plasma_pool t44_interval_decision t45_on_therapy t46_timing t47_fim_audit t48_lambda_literature t49_hepatic_bridge t50_hepatic_mito t51_pencil_geometry t52_knife_resolution t53_crossing_shape t54_dose_timing; do
  python -W error $s.py > out/$s.txt
done
python -W error r3_figs.py    # round-3 figures
python -W error r18_figs.py   # round-18 pencil figure
python -W error story_figs.py # round-19 story charts
```

Figures go to `fig/` and printed output to `out/`. Every number in `REPORT.md` appears in one of the `out/*.txt` files. `out/q15_sensors.txt` and `out/u3_joint_fisher.txt` are the repository's own Phase-10 sensor scripts (`audit/q15_sensors.py`, `audit/u3_joint_fisher.py`), run unchanged. `out/n6_profile_vs_crb.txt` is the repository's own `audit/n6_profile_vs_crb.py`, run unchanged (`python -W error audit/n6_profile_vs_crb.py` from the repository root), because the report quotes its 708 = 24%.
