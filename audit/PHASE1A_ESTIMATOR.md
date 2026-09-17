# Phase 1, step 1. The estimator

Two scripts. **No paper edits, no merge, no changes to existing constants or printed
numbers.** `PHASE4_PROPOSED_TEXT.md` stays unapplied.

| | what it does | output |
|---|---|---|
| `audit/i1_estimator.py` | the variable-projection fit, against two others, in four parts | tables |
| `audit/i2_recheck.py` | four of its claims, by routes sharing no code path | tables |

**Footing, once, for every number below.** Local model constants as
`src/common/params.py` carries them: `theta = 3000`, `K_m = 100`, `k_e = 5`, `L = 322.7273`,
`D = 6000`, `T = 14 d`. Design: the paper's own monthly record `np.arange(0, 730, 30.4)`,
`n = 25`, load CV 10%. Noise: the repository's, ordinate sd
`sigma_total(0.05, 0.05) = 21.9612`. Parameters in the paper's coordinates, `ln theta` and
`ln K_m`. Every Monte Carlo figure is over 8 independent seeds with the load record held
fixed within a seed, reported as the median with the min-max range beside it.

---

## The estimator

Hold `K_m` fixed and the balance is **linear in the other three**:

```
s_k  =  theta g_k  +  k_e m_k  +  b ,        g_k = m_k/(K_m + m_k),
```

with `s_k` the observed ordinate less `Delta_k g_k`. So the problem has exactly **one**
nonlinear parameter. Scanning it on a grid and solving the 3x3 normal equations at each node
gives the **global** minimiser of the same objective the repository's `least_squares`
minimises locally, with no start, no convergence test and no local minimum to fall into.

Two others are fitted for comparison: the `beta` form of P3 (linear in four parameters, no
grid, at the price of reweighting every residual by `(K_m + m_k)`), and the repository's own
nonlinear fit reproduced from `design/run.py:628-641`.

---

## (A) Noiseless data: all three are exact

| estimator | theta | K_m | k_e | b |
|---|---|---|---|---|
| the truth | 3000.000000 | 100.000000 | 5.000000 | 0.000000 |
| varpro (global grid) | 3000.000239 | 99.999999 | 4.999998 | -0.000003 |
| beta (P3 linear) | 3000.000000 | 100.000000 | 5.000000 | 0.000000 |
| nonlinear (the repo's) | 3000.000000 | 100.000000 | 5.000000 | -0.000000 |

Worst relative error 8.0e-8 in `theta`, 8.0e-9 in `K_m`. This rules out a coding error in any
of the three before any is asked a harder question; it makes no other claim.

---

## (B) Noisy data, and the finding that reframes the rest

**The four-parameter fit on this design is not identified without priors.** The information
filter's bound with no priors at all is `sd(ln theta) = 5.75` (8 seeds 4.50 to 6.95), a
factor of `e^5.75`. On one of those records the unpenalised information matrix has
eigenvalues `[1.5e-4, 1.9e-3, 1.1e2, 5.5e3]` and condition number `3.7e7` (`i2` check 1).

| footing | filter's bound | varpro | beta | nonlinear |
|---|---|---|---|---|
| **unpenalised**, nuisances at prior means | 5.75 (4.50-6.95) | bias 1.11, sd 1.13 | bias 0.95, sd 0.83 | bias **-6.11**, sd **8.35** |
| **penalised** with the paper's priors, nuisances at prior means | 0.1698 (0.1688-0.1709) | bias 0.0003, sd **0.0872** | cannot carry the priors | bias 0.0003, sd **0.0872** |
| **penalised**, nuisances **drawn from** those priors | 0.1704 (0.1678-0.1853) | bias 0.0619, sd **0.1856** | cannot carry the priors | bias 0.0619, sd **0.1856** |

*bias is `mean(ln theta_hat) - ln theta`; sd is of `ln theta_hat`; medians over 8 seeds,
60 replicates each.*

**Three results.**

**(a) The paper's `0.1697` is a prior-averaged posterior sd, and the realised precision
depends on where the patient sits.** At the prior means the same design and the same
estimator give **0.0872**, half the reported figure. With the nuisances drawn from their own
priors it gives **0.1856**, 9% above it, with a bias of +0.062 in `ln theta`. The reported
number is neither a floor nor a ceiling: it is an average over a population of patients the
prior describes.

**(b) The `beta` route cannot carry the paper's priors.** They are Gaussian on `b` and on
`ln k_e`, and neither is a linear function of `beta`, so the linear form is unpenalised
whatever the caller wants. Combined with its reweighting of every residual by `(K_m + m_k)`,
which means it is not minimising the paper's objective either, this is the price of the
linear route. It is a good diagnostic and a good way to see the structure. It is not the
estimator to report from.

**(c) With the priors on, the grid buys nothing, and without them it buys everything.**
Penalised, varpro and the repository's nonlinear fit agree to seven digits (bias 0.000319556
against 0.000319558): the penalised objective is well behaved and the local search finds the
global optimum. Unpenalised, the nonlinear fit is catastrophic, sd 8.35 against a bound of
5.75, because it wanders in a nearly flat four-dimensional objective. varpro's own 1.13 there
is not a triumph either; it is the `K_m` grid's hard bounds acting as a crude prior.

### Start dependence

| jitter (log sd of the start) | worse than the grid | median excess objective | bias in `ln theta` |
|---|---|---|---|
| 0 | 0.0% | 8.4e-14 | -0.0003 |
| 0.25 | 0.0% | 4.5e-14 | -0.0003 |
| 0.5 (the repository's) | 0.0% | 5.1e-14 | -0.0003 |
| 1.0 | 0.0% | 5.9e-14 | -0.0003 |
| 2.0 | **15.0%** | 5.7e-14 | 0.0274 (8 seeds -0.584 to 0.040) |

**At the repository's own jitter of 0.5 the local fit never loses to the grid**, on 320
replicates. It starts to at a jitter of 2.0, where one replicate in seven stops at a strictly
larger objective and the seed range on the bias opens to -0.58. So the repository's fit is
sound as configured, and the margin is about a factor of four in start dispersion.

---

## (C) The errors-in-variables bias at the paper's assay CV

`design/run.py:598-604` states this gap and declines to measure it: the assay error is placed
on the **ordinate**, where the Fisher matrix puts it, rather than on the abscissa where it
physically sits, and the docstring says fitting against a perturbed `m` *"would make the fit
look better than it is"*. Both halves are testable.

To first order the two models agree **by construction**: a level error `dm` enters the
balance as `b_loc dm`, and `sigma_total` is exactly `hypot(load_frac L, b_loc cv m_bar)` with
`b_loc = 29.7934`, which `i2` check 4 verifies symbolically against `R'(m_bar) = 29.793388`.
So any gap is second order.

Measured with **common random numbers**, so the two models differ only in where the assay
error is placed, 12000 paired replicates per row:

| assay CV | paired gap in `ln theta` | resolved at 2 se? | predicted (`i2` check 4) | sd ratio, EIV/ordinate |
|---|---|---|---|---|
| 2% | +0.000635 +/- 0.000404 | no | +0.000030 | 0.958 |
| **5% (the paper's)** | **+0.001521 +/- 0.001009** | **no** | **+0.000185** | **0.917** |
| 10% | +0.001681 +/- 0.002007 | no | +0.000739 | 0.901 |
| 20% | **-0.014673 +/- 0.003967** | **yes, 3.7 se** | +0.002954 | 0.904 |

**The answer is that the bias is negligible and the repository's choice is vindicated.** At
the paper's assay CV the gap is **+0.0015 in `ln theta`, 1.5 standard errors from zero after
12000 paired replicates**, and the analytic second-order prediction is **+0.000185**, which is
**0.11% of the 0.17 the paper reports**. The prediction scales exactly as `cv^2`
(`pred/cv^2 = 0.0739` at every CV, to four figures), which is what a second-order effect must
do.

**Where the expansion stops working, and it is not where the paper operates.** At 2%, 5% and
10% the measured gap is positive, of the predicted order, and never resolved from zero. At
**20%** it is **-0.0147 +/- 0.0040**, resolved at 3.7 standard errors and of the **opposite
sign** to the second-order prediction. So the expansion describes the effect up to about a
10% assay CV and fails beyond it, where terms past `cv^2` take over. Even there the gap is 6%
of that row's own sd of 0.24, so it is a limit of the *explanation*, not a new problem for the
design. It is reported because a prediction that works at 5% and breaks at 20% should be said
to break.

**The docstring's other claim is confirmed, and it matters more.** The EIV sd is **0.92x** the
ordinate model's at the paper's CV, and 0.90x at 10% and 20%. Fitting against the level the
assay reported does make the fit look better than it is, by about 8%, because the model is
evaluated where the noise put it. That is a reason to keep the repository's construction, not
to change it.

**A sign correction of my own.** An earlier draft of `i2` argued from `R'' < 0` that the EIV
fit would return a **low** `theta`. That is backwards: the fit drives the residual to zero, so
a residual that sits low must be met by **raising** the model, and the model increases in
`theta`. The estimating equation settles it, not the curvature alone. The measured gaps are
positive, as the corrected derivation requires.

---

## (D) At a constant load the estimator returns the whole fibre

Noiseless, load CV 0%. The variable-projection profile over `K_m` is the **exact** profile of
the objective, because the other three parameters are solved rather than searched.

**It is flat.** Over 601 nodes spanning three decades, the objective runs from `5.3e-26` to
`1.5e-19`, where one unit is one sigma of one residual.

| `K_m` node | theta_hat | k_e_hat | b_hat | fibre b | fibre theta | match | admissible? |
|---|---|---|---|---|---|---|---|
| 10 | 6354.54545 | 5.00000 | -2904.54545 | -2904.54545 | 6354.54545 | exact | yes |
| 30 | 4103.03030 | 5.00000 | -753.03030 | -753.03030 | 4103.03030 | exact | yes |
| 100 | 3000.00000 | 5.00000 | -0.00000 | -0.00000 | 3000.00000 | exact | yes |
| 300 | 1784.84849 | 5.00000 | 215.15152 | 215.15152 | 1784.84848 | exact | yes |
| 1000 | -1790.45449 | 5.00000 | 290.45455 | 290.45455 | -1790.45455 | exact | **no** |

*`b` in the estimator's sign convention, an offset on the model flux, which is the negative
of the fibre's `b`. Worst relative departure from the fibre 1.8e-8.*

**So the estimator does not fail at a constant load and does not return a wrong answer. It
returns the correct answer, which is a curve.** Software that reports a point here is
reporting its start, its grid's first node, or its stopping rule.

`i2` check 2 proves the flatness is an **identity** rather than a property of 601 sampled
nodes. Symbolically, the residual along the fibre is the true residual times a factor that
does not depend on the data:

```
resid_fibre / resid_true  =  (K_0 + m)/(K + m)
```

At a constant load the true residual is zero at every sample, so the fibre's is too, at every
`K` and every `m`. This agrees with the independently derived identity in `PHASE4_STRESS.md`,
`R_fibre = R_true [1 + K_0 b/(K_0 L + (L+b) m)]`; the two forms are algebraically the same.

**Where the fibre stops describing a patient** is Phase 2 item 4's bound, read in the other
coordinate. Admissible `b` in `[-272.7273, 20107.8981]` at the 10% capacity reference maps to
admissible `K_m` in `[1.5796, 645.4545]`; at 5% it is `[0.6400, 645.4545]`. The profile is
flat over the whole grid and only that sub-interval is a patient.

**With the priors on, the profile acquires a minimum, and the curvature is entirely the
prior's.** The penalised objective is minimised at `K_m = 100.0000` with `b = -1.8e-12`,
rising by 1 over `K_m` in `[93.3254, 107.1519]`. Since the data's own profile is at most
`1.5e-19` across the entire grid, `J_penalised` and `J_penalty` differ by at most `5.0e-23`:
**the data contributed nothing to that curvature and the `b` prior supplied all of it.**

That is the estimator-side statement of what `PHASE3_DESIGN.md` section 4(a) found on the
design side, and it is why a constant-load record must not be reported as an estimate with an
interval, whatever a fitting routine returns.

---

## The re-check (`i2`)

| # | claim | route | result |
|---|---|---|---|
| 1 | varpro's minimiser is global | brute force on a **different** separation: profile `(theta, b)`, scan `K_m` **and** `k_e` on a 400x400 grid | brute force cannot beat varpro; difference **+1.4e-4**, which is the 2-D grid paying for its own 1.3% resolution |
| 2 | the constant-load profile is flat | symbolically, not sampled | `resid_fibre/resid_true = (K_0+m)/(K+m)`, an identity; and `1.8e-22` at four `K_m` values that are not grid nodes |
| 3 | the calibrated sd exceeds the filter's because of nonlinearity | shrink sigma and the priors together and watch the ratio | 1.247 -> 1.087 -> 1.040 -> 1.024 -> 1.021 as the scale goes 1 -> 0.03 |
| 4 | the EIV gap's sign and size | second-order expansion projected through the same information matrix | `R'(m_bar) = 29.793388` against `b_loc = 29.793388`; predicted bias `+0.000185` at 5%, scaling as `cv^2` to four figures |

Check 3 is the one that matters for what the paper reports: **the filter is right to first
order and the 9% excess is the model's nonlinearity**, not a specification mismatch between
the filter and the fit.

---

## What this step establishes, for step 2

1. **Report the penalised estimator.** Unpenalised, this design does not identify `theta` at
   all (bound 5.75). Every number the paper quotes is conditional on the priors, and that
   should be said where the number is.
2. **Use the grid, or check that the local fit did not need it.** At the repository's own
   start dispersion the local fit is safe, on 320 replicates. The grid is free insurance and
   makes the flat case visible instead of hiding it behind a returned point.
3. **`0.1697` is a posterior sd averaged over the prior.** The same design gives 0.087 at the
   prior mean and 0.186 with the nuisances drawn. Step 2's coverage study has to say which of
   these it is calibrating against, and the answer is the third row.
4. **The errors-in-variables gap is not worth carrying.** It is 0.11% of the reported sd at
   the paper's assay CV by prediction, and not resolved from zero at 12000 paired replicates.
   It becomes visible only at a 20% assay CV, where it reverses sign and is still only 6% of
   that row's own width.
5. **The constant-load case must not be reported as a point with an interval.** The estimator
   returns a curve, and any interval around a point on it is the prior's width, not the
   data's.

## Status

`i1` and `i2` exit 0. No paper file touched, no constant changed, no printed number in the
paper altered.

**Stopping here for go-ahead before step 2 (items 2 and 2b, uncertainty).**
