# Phase 8. The simulation programme

Marks, as `PHASE4_SUMMARY.md` uses them: **[proven]** mathematics, exact or symbolic.
**[measured]** computed, footing stated. **[judgement]** a call for the author.

**The paper has not been edited and no constant in `src/common/params.py` has moved.**
Everything that had to vary is assembled in the audit script from an explicit kernel and an
explicit dose train, and the committed shape is reproduced as one member of that family.

Scripts: `audit/q1_increment_shape.py` (S1). Independent re-check: `audit/q9_recheck.py`.

---

## 0. The footing every number below is on

Unless a line says otherwise:

| what | value | where from |
|---|---|---|
| capacity | `theta = 3000`, the **normal-capacity** illustration | `src/common/params.py` through `i1_estimator` |
| `K_m`, `k_e`, load | `100`, `5`, `L = 322.7273` | same |
| noise | `sigma = 21.9612`, the ordinate model of `design/run.py:598-604` | same |
| priors | `b ~ N(0, 16.1364)` and `ln k_e ~ N(ln 5, 0.8)`, the paper's own two | same |
| design | monthly, `n = 25` over 729.6 d, load CV `10%` | the paper's own |
| dose train | uniform 14 d, `D_AMP = 6000`, steady peak `6400` | `design/run.py` |
| estimator | global scan over `K_m` in `[1, 1000]`, exact polish inside that box | `i1_estimator.varpro_exact`, extended |
| seeds | **24**, and every median below carries its min and max | the repository's `SEEDS` |

**The amplitude convention is `peak`-matched** except where section S1.6 varies it: each shape
is scaled so the steady peak increment of the nominal train equals the committed train's
`6400`. That holds the capacity ceiling of Section 7 fixed across shapes, which is why it is
the headline. It is a **[judgement]**, and S1.6 prices it.

### 0a. Where the shape constants come from, and that they are not human

Every constant that defines a shape here is a **mouse** parameter or an assumption.
Brief **R2** (the increment time course) finds **no human hepatic PCC protein or activity time
course at any time point**, so the half-life is carried as a sweep and named against R2
throughout, and so is the measurement error in S1.8. The two-compartment rates
`kdeg = 0.0075 /h`, `kq = 0.00474 /h`, `ke0 = 0.242 /h` and `Slope = 55.6` are Attarwala et
al. 2023 Table 1, which `audit/k1_sources.py` finding 2 already established are mouse.

---

## S1. The increment time course: analyst against truth

### S1.0 What the re-check caught before anything was reported

Three defects, all found by `audit/q9_recheck.py` against routes sharing no code path with
`q1`, and all recorded here rather than quietly repaired.

1. **The profile was left unpolished**, so it carried the delta-method form of the `k_e`
   prior (Gaussian on `k_e` with sd `0.8 k_e`) rather than the paper's prior on `ln k_e`.
   Down the profile's low tail the compensating direction is `k_e` **rising**, past five times
   its reference, and there the two are not the same prior: the delta-method form charges
   about seven times what the log prior does. Intervals came out about a third too narrow.
   Found by check 5, against `n6_profile_vs_crb.nll`.
2. **The `1.92` cut was applied to an objective that is twice the negative log likelihood.**
   `i1_estimator.obj_exact`'s scale is `2 x` `n6`'s `nll`, so the 95% cut on that scale is
   `3.84`. Same check.
3. **The polish walked `K_m` past the top of its own scan grid** on the badly misspecified
   cells, where `theta g(m)` collapses into `(theta/K_m) m` and stops being separable from the
   `k_e m` column, returning capacities of order `1e7`. Found by check 4 once that check was
   made domain-consistent. The polish is now bounded to the scan's box, so the scan and the
   polish are one estimator; the unbounded fit reaches `K_m = 2607` and
   `theta_hat/theta = 20.6` on an Attarwala record, and check 4 prints that beside the bounded
   answer rather than hiding it.

**All 21 re-checks now pass.** The ones that matter most are not (1) to (3) but these:

- the two-compartment kernel by hand from the Laplace transform against `q1`'s
  eigen-decomposition, worst relative difference `< 1e-9`. **[proven]**
- the committed kernel and train against `design.run.delta`, worst relative difference
  `1.4e-16` over 0 to 600 d, which is what makes "the committed shape" one member of this
  family rather than a second spelling of it. **[measured]**
- `q1`'s fit against `i1_estimator.varpro_exact` on the well-specified case, worst relative
  difference in `theta` `4.7e-6` over 6 records, and in `K_m` `3.1e-7`. **[measured]**
- the same two-compartment kernel by eigen-decomposition against `solve_ivp`, the route
  `audit/k1_sources.py` takes, `9.3e-12` over 0 to 300 d. **[measured]**
- the published dose arithmetic in exact rationals: `346` doses, dose-weighted mean
  `1419/3460 = 0.4101156069` mg/kg. **[proven]**

### S1.1 The headline grid **[measured]**

Nominal 14 d train, peak-matched, `theta = 3000`, monthly `n = 25`, 24 seeds.
`theta_hat/theta` above 1 is capacity read **too high**. `cover` is the fraction of seeds whose
95% profile interval contains the truth, against a nominal 95%; with 24 seeds its own standard
error at 95% is 4.4 points, so it separates 95% from 50% and not 95% from 88%.

| truth | analyst | `theta_hat/theta` median [min, max] | RMSE `ln theta_hat` | RMS resid / sigma | cover |
|---|---|---|---|---|---|
| exp 2.5 d | fixed 3.5 d | **0.842** [0.686, 1.033] | 0.188 | 1.137 [0.88, 1.51] | 91.7% * |
| exp 2.5 d | half-life free | 0.977 [0.841, 1.164] | 0.082 | 0.968 | 100% |
| exp 2.5 d | two-comp free | 1.178 [1.007, 1.533] | 0.219 | 1.879 | 100% |
| **exp 3.5 d (null)** | **fixed 3.5 d** | **0.975** [0.797, 1.184] | 0.085 | **0.980** [0.72, 1.32] | **100%** |
| exp 3.5 d | half-life free | 1.006 [0.704, 1.204] | 0.143 | 0.964 | 100% |
| exp 3.5 d | two-comp free | 1.326 [1.035, 1.699] | 0.310 | 1.809 | 100% |
| exp 5.0 d | fixed 3.5 d | **1.310** [1.053, 1.618] | 0.302 | 1.080 [0.80, 1.33] | **100%** |
| exp 5.0 d | half-life free | 1.047 [0.479, 1.339] | 0.243 | 0.969 | 100% * |
| exp 5.0 d | two-comp free | 1.682 [1.313, 2.146] | 0.532 | 1.676 | 58.3% |
| exp 7.0 d | fixed 3.5 d | **1.934** [1.475, 2.469] | 0.686 | 1.184 [0.83, 1.43] | 16.7% |
| exp 7.0 d | half-life free | 1.071 [0.440, 1.602] | 0.304 | 0.978 | 100% * |
| exp 7.0 d | two-comp free | 2.247 [1.761, 2.954] | 0.820 | 1.463 | 0% |
| Attarwala 2-comp | fixed 3.5 d | **5.637** [2.440, 8.609] | 1.730 | 2.347 [2.12, 2.83] | 4.2% |
| Attarwala 2-comp | half-life free | **0** on all 24 seeds | n/a | 1.593 | 0% |
| Attarwala 2-comp | two-comp free | 1.051 [0.730, 1.210] | 0.145 | 0.931 [0.68, 1.33] | 100% |

`*` marks a cell where the interval failed to close below on some seeds, so the truth is
contained for a reason that is not precision: 10 of 24 on the exp 2.5 d row and 1 of 24 on
each half-life-free row. Every other 100% in the table is a closed interval on every seed.

**The review's gap B is confirmed and its range is a statement about the amplitude
convention, not only about the shape.** At the committed 3.5 d assumption, a truth at 5.0 d
reads capacity `1.31` high and a truth at 7.0 d reads it `1.93` high, which brackets the
review's `1.3 to 2.0` at the same peak-matched convention. A truth at 2.5 d reads it `0.84`
**low**: the bias is monotone in the truth's half-life and changes sign at the assumption.

**The new number is the Attarwala truth: `5.64` [2.44, 8.61], not `1.3 to 2.0`.** The review
reported the two-compartment truth inside the same `1.3 to 2.0` band as the 5 to 7 day
exponentials. Read at the peak-matched convention with the published rates it is about three
times worse than that, and it stays between `5.1` and `6.0` under all three conventions
(S1.6). **The two-compartment truth is not a slightly slower exponential and should stop being
grouped with one.**

### S1.2 Freeing the amplitude does not help, and on one truth it is worse **[measured]**

The review reports that freeing the amplitude with a 12% prior does not help. Checked rather
than inherited, on the same records:

| truth | fixed 3.5 d | fixed 3.5 d + free amplitude (12% prior) |
|---|---|---|
| exp 2.5 d | 0.842 [0.686, 1.033] | 0.834 [0.680, 1.020] |
| exp 3.5 d | 0.975 [0.797, 1.184] | 0.971 [0.797, 1.178] |
| exp 5.0 d | 1.310 [1.053, 1.618] | 1.308 [1.053, 1.612] |
| exp 7.0 d | 1.934 [1.475, 2.469] | 1.928 [1.471, 2.459] |
| Attarwala | 5.637 [2.440, 8.609] | **0.111** [0.096, 0.138] |

**Confirmed on the exponential truths**, where the two agree to under half a per cent and the
residual is identical to three decimals: a free amplitude buys nothing at all, because the
misfit is in the shape and a scale cannot reach it.

**Not confirmed on the Attarwala truth, where it is a different failure rather than the same
one.** The fit pays `8.3` prior standard deviations to drive the amplitude to
`kappa = 0.0093` [-0.0002, 0.0169],
switches the increment off, and returns a capacity a factor of nine **too low** instead of a
factor of five too high. The residual falls from `2.35` to `1.54 sigma`, so the fit looks
**better** while the answer is wrong in the opposite direction. **[measured]** That a 12%
prior is soft enough to permit this is worth saying plainly: the freedom does not fix the
bias, it relocates it, and it makes the residual less useful as a warning.

### S1.3 Freeing the half-life works, and only inside the right family **[measured]**

Carried as a free parameter, the half-life is recovered and the bias goes away, on every
single-exponential truth:

| truth | half-life recovered, median [min, max] | `theta_hat/theta` |
|---|---|---|
| exp 2.5 d | 2.41 d [1.93, 3.19] | 0.977 [0.841, 1.164] |
| exp 3.5 d | 3.33 d [2.53, 4.92] | 1.006 [0.704, 1.204] |
| exp 5.0 d | 4.61 d [3.26, 7.35] | 1.047 [0.479, 1.339] |
| exp 7.0 d | 6.63 d [4.23, 9.56] | 1.071 [0.440, 1.602] |

It costs precision, as the review said: the min-max spread on `theta_hat/theta` widens from
`[1.48, 2.47]` to `[0.44, 1.60]` at the 7.0 d truth, and the RMSE at the null rises from
`0.085` to `0.143`. That is the price of the fifth parameter and it is worth paying.

**On the Attarwala truth it fails completely, and it fails quietly.** The best single
exponential is `14.28 d` [12.55, 15.86], and with the shape that slow the increment explains
the whole record: `theta_hat` is driven to the floor on **all 24 seeds**, and the RMS residual
is `1.59 sigma`, which is **lower** than the `2.35` the fixed analyst leaves. So the wrong
family, given a free half-life, returns "this patient has no residual capacity" and shows a
better-looking fit while doing it.

**The symmetric statement, which is the one to carry.** The two-compartment analyst recovers
the Attarwala truth almost exactly, `1.051` [0.730, 1.210] at a residual of `0.931`
[0.675, 1.328], which is the null cell's own range. It is then **biased on every exponential
truth**, from `1.18` at
2.5 d to `2.25` at 7.0 d. **No analyst family in this grid is safe. The one that matches the
truth recovers capacity and every other one is biased**, and the direction of the bias depends
on which pair you are in, not on any property of the estimator.

### S1.4 How loud is the warning **[measured]**

The residual reference is stated once, because the obvious one is wrong. `sqrt((n-p)/n)` is
the **unpenalised** value, `0.917` at `n = 25, p = 4`, and the paper's fit is penalised, so a
well-specified fit sits above it. Re-check 8 confirms both halves: the same records fitted
with the priors removed give `0.947`, against `sqrt((n-p)/n) = 0.917`, and the penalised fit
gives `0.972`, which is the larger. **The reference to read against is the measured null
cell, `0.980` [0.723, 1.324].**

| truth, fixed 3.5 d analyst | RMS resid / sigma | against the null's 0.980 |
|---|---|---|
| exp 3.5 d (null) | 0.980 [0.723, 1.324] | . |
| exp 5.0 d | 1.080 [0.799, 1.332] | **+10%**, inside the null's own seed range |
| exp 7.0 d | 1.184 [0.827, 1.427] | +21%, inside the null's own seed range |
| exp 2.5 d | 1.137 [0.882, 1.505] | +16%, inside the null's own seed range |
| Attarwala | 2.347 [2.117, 2.834] | **+140%**, clear of it on every seed |

**This is the finding the paper's `sec:structuredep` claim needs qualifying against, and it
cuts both ways.** A 5 to 7 day truth read as 3.5 d raises the residual by 10 to 21%, and the
null cell's own seed-to-seed range is `0.723 to 1.324`, a spread of `+/-30%`. **So a single
record cannot tell the two apart: the misspecification is genuinely inaudible in exactly the
diagnostic a reader would use.** But the two-compartment truth is at `2.35 sigma` on every one
of 24 seeds, more than three times the null's worst. **That one is audible, loudly, and the
paper's "the residual of that fit can be under a fifth of one sigma" does not describe it.**

The coverage column carries the same split and adds one thing the residual does not:

- at a **5.0 d** truth the bias is `1.31` and the coverage is **100%**. The interval is wide
  enough to cover a 31% error, so neither the residual nor the interval reports anything.
  **That is the worst cell in the table, because nothing at all goes wrong visibly.**
- at a **7.0 d** truth coverage collapses to `16.7%` and at the Attarwala truth to `4.2%`. By
  then the residual is also moving, so the two diagnostics fail together rather than
  independently.

### S1.5 The trial's real, irregular dose train **[measured]**

Built from what the Nature 2024 interim analysis publishes and
`audit/PHASE6_REVIEWS_SEP2026.md` section 8 records: cohort 1 Q3W and cohorts 2 to 5 Q2W,
0.30 to 0.90 mg/kg, a delayed dose, a skipped dose, a de-escalation, and the one slip the
paper quantifies (dose 3 at day 52 against a 42 d nominal). mg/kg is put on the committed
scale by the dose-weighted mean of the published 346 doses, `0.41012` mg/kg, so `D_AMP = 6000`
stands for `0.410` mg/kg. The **size** of the dose-2 delay is **ASSUMED** at 7 d, being the
published 10 d slip on a 21 d interval scaled to a 14 d one, and the composite train is a
**construction** carrying every published feature in one record rather than a transcription of
any one participant. **[judgement]**

**(a) An irregular train does not rescue the shape problem.** With the analyst holding the
real train, so that only the shape is wrong:

| truth | nominal Q2W | C5 de-escalated | C1 Q3W then Q2W | composite |
|---|---|---|---|---|
| exp 7.0 d | 1.934 | 2.293 | 1.642 | 1.720 |
| Attarwala | 5.637 | 6.270 | 6.614 | 2.167 |

Better excitation moves it by a few tens of per cent in either direction and closes nothing.
**The shape error is not a conditioning problem, so a better-conditioned design does not
touch it.**

**(b) Getting the dose train wrong is a larger lever than getting the shape wrong, and unlike
the shape it is audible.** With the analyst assuming the nominal Q2W train while the truth ran
the real one, at the **correctly specified** 3.5 d shape:

| train the truth ran | `theta_hat/theta` | RMS resid / sigma |
|---|---|---|
| C5 de-escalated | 0.779 [0.629, 1.075] | 1.934 [1.52, 2.22] |
| C1 Q3W then Q2W | **6.426** [5.219, 8.174] | 2.308 [2.10, 3.09] |
| composite | **5.285** [4.054, 6.153] | 2.582 [2.19, 3.20] |

A Q3W patient analysed as Q2W reads capacity **six times** too high with the shape assumption
correct, against `1.93` for a 7 d shape error on the right train. **[measured]** But the
residual is `2.3` to `2.6 sigma`, clear of the null's worst seed, so this failure announces
itself. **That is the useful asymmetry: the record error is bigger and detectable, the model
error is smaller and hidden.** It also means a design that only records "biweekly" is worse
than one that records the shape wrongly but the dates correctly.

### S1.6 The amplitude convention, which is a footing and not a result **[judgement]**

Comparing two shapes requires saying what is held equal and there is no right answer, so all
three are computed. `dose` holds the single-dose peak, so each administration deposits the
same and a slower kernel accumulates more. `peak` holds the steady peak increment, so the
capacity ceiling does not move. `mean` holds the cycle-mean exposure.

`theta_hat/theta`, fixed 3.5 d analyst:

| truth | dose | peak | mean |
|---|---|---|---|
| exp 2.5 d | 0.881 | 0.842 | 0.619 |
| exp 3.5 d | 0.975 | 0.975 | 0.975 |
| exp 5.0 d | 1.226 | 1.310 | 1.605 |
| exp 7.0 d | 1.704 | 1.934 | 2.619 |
| Attarwala | 6.024 | 5.637 | 5.134 |

**The direction is a finding; the size is a statement about the convention as much as about
the shape.** At a 7.0 d truth the bias runs `1.70` to `2.62` across the three, a spread of
`1.5x`. The review's `1.3 to 2.0` sits inside the `dose` and `peak` columns for the 5 to 7 day
band. The Attarwala row is the stable one, `5.1` to `6.0`, which is a second reason to report
it separately from the exponentials.

**Any number in this programme that compares two shapes must name its convention.** S2 onward
should inherit `peak`, on the ground that it is the only one of the three under which the
admissibility check of Section 7 does not move when the shape does.

### S1.7 The disease range, where a bias figure is not the operative statement **[measured]**

At `theta = 150` on `n4_theta_scope`'s convention A, footing re-solved:
`L = 322.7273`, `m_bar = 54.0231`, `sigma = 22.1774`, `b_loc = 5.6323`.

**The 95% profile interval fails to close below on 24 of 24 seeds in every cell but one**,
which is `sec:scope`'s own result reproduced here from a different direction. A `cover`
percentage on such a row is not a coverage statement: an interval unbounded below contains the
truth for a reason that has nothing to do with the estimate. The point estimate is likewise
unusable: `theta_hat` is at the floor on 24 of 24 seeds against the exp 2.5 d truth and on 10
of 24 against the well-specified one.

The one thing worth reading is the residual, and it says the misspecification is **not** more
hidden here:

| truth, fixed 3.5 d analyst | RMS resid at `theta = 3000` | at `theta = 150` |
|---|---|---|
| exp 3.5 d (null) | 0.980 | 0.957 |
| exp 5.0 d | 1.080 | 1.175 |
| exp 7.0 d | 1.184 | 1.421 |
| Attarwala | 2.347 | 3.185 |

**The warning is if anything louder in the disease range, and it is warning about a quantity
that was not recoverable there in the first place.** Note that the booked `D = 6000` is not
admissible at `theta = 150` on readings R2 or R3 of Section 7, so this block is a
misspecification measurement at an inadmissible dose and not a design recommendation.

### S1.8 What measurement would pin the shape

Part 1 says the shape decides the answer and the fit does not reliably announce a wrong one,
so the shape has to come from outside the marker record. Three things are fixed first because
each changes the answer: the **readout** (`delta`, the increment itself, against `capacity`,
`theta + Delta`, which is what a functional assay returns and is damped by
`Delta/(theta + Delta)`); the **accessible window**, because phase 0 is the instant of the
infusion and is both where the families differ most and the phase a clinic cannot book; and
the **measurement CV**, swept because **brief R2 supplies none** and `REVIEW_LOG.md` **R4**
already retracts the one in-vivo repeatability figure the September reviews offered.

**One phase identifies nothing, at any number of measurements. [proven]** The amplitude has to
be carried and profiled out, because nothing independently fixes how much increment one dose
delivers. Two phases is the minimum and, the model being linear in the two log parameters,
also the optimum, so `N` splits evenly between two phases.

**The closed form, confirmed by re-check 7 to 5e-11. [proven]**
`sd(ln t_half) = sqrt(2) / (lambda |tau_1 - tau_2|)`, minimised at the two **ends** of the
accessible window. The phase-independent part of the sensitivity is absorbed by the free
amplitude, which is why only the separation of the two phases matters.

**The gate.** A truth of 5.0 d read as 3.5 d is already the `1.31` bias of S1.1, so a useful
measurement separates them at two standard deviations: `sd(ln t_half) <= 0.178`.

**(i) Estimating the half-life, within the single-exponential family. [measured]**
`N` is the total number of increment measurements, half at each phase.

| readout | earliest phase | best pair | sd per unit CV | N at 10% CV | at 20% | at 30% | at 50% |
|---|---|---|---|---|---|---|---|
| `delta` | 1 h | 0.04 and 13.96 d | 0.513 | 2 | 2 | 2 | 6 |
| `delta` | 6 h | 0.25 and 13.75 d | 0.529 | 2 | 2 | 2 | 6 |
| `delta` | 1 d | 1.00 and 13.00 d | 0.595 | 2 | 2 | 4 | 6 |
| `capacity` | 1 h | 0.04 and 8.40 d | 2.279 | 4 | 14 | 30 | 82 |
| `capacity` | 6 h | 0.25 and 8.53 d | 2.342 | 4 | 14 | 32 | 88 |
| `capacity` | 1 d | 1.00 and 9.02 d | 2.587 | 6 | 18 | 38 | 106 |

**(ii) Telling the two families apart, which is the harder question and the one S1.1 says
actually decides the answer. [measured]** Both shapes are given their best free amplitude
first, so what is left is shape alone; the statistic is the difference of two log readings and
the target is a 3-sigma call.

| readout | earliest phase | best pair | separation per unit CV | N at 10% CV | at 20% | at 30% | at 50% |
|---|---|---|---|---|---|---|---|
| `delta` | 1 h | 13.96 and 0.04 d | 1.359 | 2 | 2 | 2 | 4 |
| `delta` | 6 h | 13.75 and 0.25 d | 0.679 | 2 | 2 | 4 | 10 |
| `delta` | 1 d | 13.00 and 3.48 d | 0.515 | 2 | 4 | 8 | 18 |
| `capacity` | 1 h | 13.96 and 0.04 d | 0.515 | 2 | 4 | 8 | 18 |
| `capacity` | 6 h | 13.75 and 0.25 d | 0.149 | 10 | 34 | 74 | 202 |
| `capacity` | 1 d | 13.00 and 3.78 d | 0.087 | 24 | 96 | 214 | 592 |

**The answer, stated for a protocol.**

- **Two measurements, at the two ends of the dose cycle**, one as soon after an infusion as
  the clinic can book and one immediately before the next dose. Not four spread evenly: the
  amplitude absorbs everything a middle phase carries.
- **If the readout is the increment itself, two are enough** to pin the half-life to the gate
  at any measurement CV up to 30%, and to separate the two families at 3 sigma at any CV up to
  30%. Brief **R2** reports that this assay has never been run in a human, so this row prices
  an assay that does not exist.
- **If the readout is functional**, `theta + Delta`, the same two measurements still pin the
  half-life at 10% CV but the family question needs **10 measurements at 10% CV and 34 at
  20%**, and it needs the early draw: moving the earliest accessible phase from one hour to
  six costs a factor of **3.4** in the separation, and to one day a factor of **5.9**.
- **The two questions want different draws and only one of them is flexible.** The half-life
  wants the ends of the cycle. The family wants the hours immediately after an infusion, where
  one shape jumps and the other rises, and the whole discriminating signal is there: at one
  hour the two peak-matched shapes differ by `ln(2comp/exp) = -1.01`, at six hours by `-0.08`,
  and the departure has all but vanished.

**The practical consequence.** The trial already brings a patient in for the infusion, so the
early draw is the one the schedule gives away, and it is the one that carries the family
question. This is the same draw Section 8's Tier 1 asks for, for a different reason, and the
two asks should be booked together rather than separately. **[judgement]**

### S1.9 What S1 changes for S2 onward

1. **Carry the two-compartment truth separately from the 5 to 7 day exponentials.** Its bias
   is `5.6`, not `1.3 to 2.0`, and it is the one truth the residual actually reports.
2. **Name the amplitude convention on every cross-shape number.** `peak` is proposed as the
   standard, since it is the only one under which the admissibility check does not move.
3. **The dose train is a bigger lever than the shape and belongs in the specification.**
   A Q3W patient analysed as Q2W reads capacity six times high with the shape assumption
   correct. `src/design/run.py`'s `DOSE_INTERVAL = 14.0` annotation, "the biweekly schedule
   selected for mRNA-3927", is one of two published schedules and neither is adhered to.
4. **The residual is a usable warning for a wrong family and not for a wrong half-life**, so
   S4's joint stress scenario should report both and not one.
5. **Freeing the half-life is not a safe default.** It removes the bias inside the right family
   and, in the wrong one, returns capacity zero at a *better*-looking residual.

---

## What is open after S1

- **The shape itself.** Nothing here decides it. Brief **R2** finds no human hepatic PCC
  protein or activity time course, so the whole grid is a sensitivity analysis and not a
  correction. S1.8 prices the measurement that would close it.
- **The measurement CV on an increment reading**, which is swept for the same reason and which
  **R2** would pin.
- **The size of the dose-2 delay** in the C5 train, which the Nature paper does not give. Swept
  in `q1`'s `irregular_train(delay=...)`, held at 7 d for the tables above.
- **Whether the paper's `sec:structuredep` sentence should be qualified.** The passage says
  the residual of a misspecified fit "can be under a fifth of one sigma". S1.4 measures `+10%`
  to `+21%` above the null for a wrong half-life, which is consistent with that sentence's
  spirit, and `2.35 sigma` on every seed for the published two-compartment truth, which is
  not. **No text is proposed here and `paper/` is untouched.** This is for the author.
- **S2 to S5 have not been run.** S1 stops here as instructed, because S1 may change the
  footing the rest uses, and points 1, 2 and 4 above are exactly such changes.
