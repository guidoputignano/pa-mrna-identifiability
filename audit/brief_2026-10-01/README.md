# Brief of 1 October 2026: practical identifiability of enzyme capacity

The findings are in **[`REPORT.md`](REPORT.md)**. This directory holds the code behind them.

It builds on the repository rather than beside it. The footing (convention A, the monthly record, σ, the priors, the dose train and the 24 intake seeds) is imported from `audit/n4_theta_scope.py`, and `core.fim` reproduces `n4.sd_ln_theta` to 1e-13. Nothing here writes into `results/` or `paper/`, so the artefact guard is untouched: it globs `audit/*.py`, not subdirectories.

| script | brief task | runtime |
|---|---|---|
| `core.py` | the repository's FIM, the ln-m form, profiling | — |
| `t1_identifiability.py` | 1: U-curve, correlations, elasticity, profile | 15 s |
| `t2_cohort.py` | 2: joint and hierarchical cohort fit, with Monte Carlo | 2 min |
| `t3_variance.py` | 3: the 1/n versus n-independent split, and 15× sampling | 5 s |
| `t4_set.py` | 4: set identification and the one-sided bound | 1 min |
| `t5_value.py` | 5: what each extra measurement is worth | 5 s |
| `dyn.py` | the ODE with forward sensitivities, for tasks 6 and 7 | — |
| `t6_estimators.py` | 6: FD, integral, NLS, EKF and iterated EKF smoother | ~40 min |
| `t7_design.py` | 7: washout and drug-forcing design | 1 min |
| `plotstyle.py` | the repository's rcParams with a validated palette | — |

To run them all, from this directory:

```bash
pip install -r ../../requirements.txt
for s in t1_identifiability t2_cohort t3_variance t4_set t5_value t7_design t6_estimators; do
  python -W error $s.py > out/$s.txt
done
```

Figures go to `fig/` and printed output to `out/`. Every number in `REPORT.md` appears in one of the `out/*.txt` files.
