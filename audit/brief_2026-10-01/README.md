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
| `*_recheck.py` | an independent verifier's re-derivation of each follow-up script, sharing no code with it | 3–25 min |
| `plotstyle.py` | the repository's rcParams with a validated palette | — |

To run them all, from this directory:

```bash
pip install -r ../../requirements.txt
for s in t0_footing t1_identifiability t2_cohort t3_variance t4_set t5_value t7_design t6_estimators t6b_ma1 \
         t8_floor t9_omega t10_factorial t6c_iv t6d_ekf t11_fim_check t12_flux_mc; do
  python -W error $s.py > out/$s.txt
done
```

Figures go to `fig/` and printed output to `out/`. Every number in `REPORT.md` appears in one of the `out/*.txt` files. `out/n6_profile_vs_crb.txt` is the repository's own `audit/n6_profile_vs_crb.py`, run unchanged (`python -W error audit/n6_profile_vs_crb.py` from the repository root), because the report quotes its 708 = 24%.
