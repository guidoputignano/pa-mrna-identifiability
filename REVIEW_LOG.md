# Review log

Where the mathematics of Section 3 stands, so that the next session does not
have to re-derive it. Every number below is produced by `src/` and asserted in
`tests/test_models.py`; none is typed into the paper by hand.

## What each commit closed

| commit | closed |
|---|---|
| `11c86ba` | The sign error in §3.6 item 2, which reversed its own conclusion: both parameters enter the margin positively, the quadratic terms cancel along the fibre, and the margin is linear where `theta` is quadratic. Also: `a = theta g(m_bar)^2` carries `Km`, so "the intercept carries `theta` alone" was false and the admissible set is a fibre spanning 169x; `det M = g(m)g(m')(m'-m)` exactly, with `g(m_bar)^2` its midpoint approximation, 8% out at a 55% excursion; Figure 8's counting assumes one observation per patient; Figure 5's pair was matched in value only; the 1.4 load-error factor is the flux ratio `C/theta g`; (16) is Model A's `theta_crit` while Figure 10 plots Model B's; and `m_bar` denoted three quantities, the tolerable level now `\mtol`. |
| `7731e8d` | `s_peak` stated absolutely disagreed across the note (4.65x baseline in Model A against 1.18-1.50x in Model C); it is now a fold of the patient's own baseline. That reframed the margin result: the sign of `mu` is determined outside a computable window of `s_peak`, `(3.25, 3.77)x` baseline at `m_tol = 40`, 14.8% wide, closing to 1.2% at `m_tol = 15` and opening past 100% at `m_tol = 320`. Also fixed the abstract's surviving "carried entirely by the intercept", the warrant for "given `Km`" (Assumption 2 grants sharedness, not knowledge), and the band convention (`Delta` is full width). |
| `98c65e8` | Evaluated (7), which the note derived and never computed. `sd(theta_hat)/theta = sd(a_hat)/a` exactly, so the relative precision is identified although `theta` is not. `sigma` and `n` defined for the steady-state design, `n` counting instrument events. Curvature resolvability derived with design constant `sqrt(180)`. Errors-in-variables simulated. The null anchor given a variance. Model C shown to have no fibre given `alpha`. |
| `f981cad` | Five corrections to the above, all mine. `sigma` had been quoted in flux units from a constant documented as a percentage, implying a load recorded to 1.55%; both criteria are now functions of the load error. The precision result does not close the location question. The anchor cost is a ratio scaling as `1/theta`, crossing one at `theta = 605` and `1100`, so it amplifies at Model C's `theta = 300`. The errors-in-variables sign flip was an inadmissible fit (`ke < 0` in 96% of them), detectable from the fit alone. "No fibre" is given `alpha`, and `Kc` costs a factor of 1.52 in recovered `theta` over its own arbitrary range. |

## Still open

1. **A defensible bound on `sigma`.** Every conclusion in §3.7 is conditional on how well the substrate load is recorded, and the note never bounds it. Prescribed protein intake is documented; gut flora and unrecorded triggers are not. Until this is bounded from the records the identifiability answer is a family of curves rather than a number.
2. **The hazard reformulation of §7.** `s_peak` and `m_tol` are now the binding arbitrariness rather than identifiability, and `lambda_MDE(theta) = nu P(A > A_crit(theta))` is what removes both. The natural history study already fits a hazard of that shape on a measured biomarker; what is proposed is substituting a mechanistic `theta`. Not attempted.
3. **Whether Model C's four constants are a net gain.** `R0`, `Kc`, `Mi` and `rho0` replace `ke`. They remove Model A's fibre given `alpha` but `Kc` alone costs 1.52x in recovered `theta` and 1.56x in `theta_crit`. Model C's identifiability is not analysed anywhere in the note.
4. **The marker choice.** §3.2 defers it and §7 poses it as a design question answerable on `m_bar/Delta` at a fixed clinical excursion. What is missing is the paired marker-and-management data to run the comparison, and §5.6 adds that no marker is currently known to track the state through a decompensation at all.

## Two failure modes worth watching

Both recurred within one review cycle:

- **A ratio reported as a constant.** The 1.4 load-error amplification and the 0.36 anchor attenuation were both quoted as properties of the model. Both are ratios in `theta` and both invert for the severe patient.
- **A unit carried across scales.** `ASSAY_CV` is a percentage and was used as absolute flux. The definition site now says so.
