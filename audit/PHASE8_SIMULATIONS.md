# Phase 8. The simulation programme

Marks, as `PHASE4_SUMMARY.md` uses them: **[proven]** mathematics, exact or symbolic.
**[measured]** computed, footing stated. **[judgement]** a call for the author.

**The paper has not been edited and no constant in `src/common/params.py` has moved.**
Everything that had to vary is assembled in the audit script from an explicit kernel and an
explicit dose train, and the committed shape is reproduced as one member of that family.

Scripts: `audit/q1_increment_shape.py` (S1), `audit/q1b_reconcile.py` (S1 against the
25 September probes). Independent re-check: `audit/q9_recheck.py`. The probes themselves:
`audit/probes_2026-09-25/`.

**This file was revised on 26 September 2026.** The first S1 report (commit `043a3ec`)
carried numbers and four conclusions that the reconciliation against the probes withdrew.
What changed and why is section S1.R, first; `REVIEW_LOG.md` R54 records it. Everything
after S1.R is the corrected S1.

---

## 0. The footing every number below is on

Unless a line says otherwise:

| what | value | where from |
|---|---|---|
| capacity | `theta = 3000`, the **normal-capacity** illustration | `src/common/params.py` through `i1_estimator` |
| `K_m`, `k_e`, load | `100`, `5`, `L = 322.7273` | same |
| noise | `sigma = 21.9612`, the ordinate model of `design/run.py:598-604` | same |
| priors | `b ~ N(0, 16.1364)`, `ln k_e ~ N(ln 5, 0.8)`, and where the dose scale is free, `ln kappa ~ N(0, 0.12)` | the paper's own, and `g_design_common.sd_ln_theta`'s coordinate for `kappa` |
| design | monthly, `np.arange(0, 730, 30.4)`, `n = 25`, load CV `10%` | the paper's own |
| dose train | uniform 14 d from `t = 0`, `D_AMP = 6000`, steady peak `6400` | `design/run.py` |
| **dose-day draw** | **read PRE-dose**: a draw at a dose instant is taken before the infusion, for the truth and the analyst alike. Only `t = 0` is affected on this schedule | new, S1.R |
| estimator | global scan over `K_m` in `[1, 1000]`, exact polish inside that box | `i1_estimator.varpro_exact`, extended |
| seeds | **24**, and every median below carries its min and max | the repository's `SEEDS` |

**The amplitude convention of the tables below is `peak`-matched**, as the first report's
was, so that the correction and the convention change are not confounded. S1.6 gives all
three; S1.9 records the author's decision that `dose` becomes primary from S2 on.

### 0a. Where the shape constants come from, and that they are not human

Every constant that defines a shape here is a **mouse** parameter or an assumption. Brief
**R2** (the increment time course) finds **no human hepatic PCC protein or activity time
course at any time point**, so the half-life is carried as a sweep and named against R2
throughout, and so is the measurement error in S1.8. The two-compartment rates
`kdeg = 0.0075 /h`, `kq = 0.00474 /h`, `ke0 = 0.242 /h` and `Slope = 55.6` are Attarwala et
al. 2023 Table 1, which `audit/k1_sources.py` finding 2 established are mouse.

---

## S1.R Reconciliation against the probes, and what it withdrew

### The question

For the Attarwala two-compartment truth `q1` reported `theta_hat/theta = 5.64` [2.44, 8.61]
at a `2.35 sigma` residual. The probe `p3_shape.py`, on noiseless records, gives `1.37` to
`1.40` at `0.33` to `0.37 sigma`. On the exponential truths the two agreed. So either the two
scripts model different objects or one of them is wrong. `audit/q1b_reconcile.py` puts both
on identical noiseless records (the probe's own 8 load records) and fits each record with
**both** estimators, so that a modelling difference and an estimator error cannot be
confused.

### 1. The two kernels **[proven]**

| | `q1` | probe |
|---|---|---|
| protein block | `[[-(kq+kdeg), kq], [kq, -kq]]` | the same |
| rates | `kdeg = 0.0075 /h`, `kq = 0.00474 /h`, times 24, so **per day** | the same |
| protein half-lives | **1.99 d and 11.81 d** | **1.99 d and 11.81 d** |
| where the dose enters | **the effect compartment**, `Ce(0) = 1`, `ke0 = 0.242 /h`, `Slope = 55.6` | **PCC directly**, a unit impulse |
| third mode | **0.119 d**, the effect compartment | none |
| shape after an infusion | **rises** from zero over hours | **jumps** |
| mRNA PK | none beyond the effect compartment | none |

Both share the two half-lives `PHASE5_ROBUSTNESS.md` finding 3 gives. They differ in one
thing: `q1` reads Attarwala's Suppl. Table S1(B) with its effect compartment, and the probe
without it.

**Where each amplitude convention is matched.** In `q1`, `peak` is the steady-state peak of
the infinite nominal train and `mean` the steady-state cycle mean (the impulse area over the
interval), both closed form; `dose` is the **single-dose** peak of one impulse, which is a
first-cycle quantity. The probe matches `peak` and `mean` on days 600 to 614 of a 730 d
train, which is steady state to within `exp(-600 ln2/11.81)`. So both scripts match `peak`
and `mean` at steady state, and `q1`'s `dose` is the only first-cycle convention in either.

### 2. The same kernel and convention: the fitters agree **[measured]**

On identical records the two estimators agree to **5e-7** in the median, on every row of
`q1b` R2: both kernels, three schedules, amplitude fixed and freed. Reproduced here:

| kernel | draw | schedule | probe fitter | `q1` fitter |
|---|---|---|---|---|
| probe | as run | monthly | 1.37 [1.36, 1.38] at 0.36 sigma | 1.37 [1.36, 1.38] at 0.36 sigma |
| probe | as run | 14.5 d | 1.40 [1.39, 1.41] | 1.40 [1.39, 1.41] |
| probe | as run | 2-daily | 1.36 [1.35, 1.36] | 1.36 [1.35, 1.36] |
| `q1` | as run | monthly | **5.81** [4.02, 9.03] at 2.20 sigma | **5.81** [4.02, 7.33] at 2.20 sigma |

**So there is no estimator error in the fixed-amplitude cells. The 5.8 is a property of the
record, and specifically of one draw in it.**

### 3. What carries the 5.8: one draw at the instant of an infusion **[measured]**

The monthly record `np.arange(0, 730, 30.4)` puts exactly one draw at a dose instant in two
years, `t = 0`, because `30.4 k` is a multiple of 14 only at `k = 35`. `design.run.delta`
reads it on the **post**-dose side (`t >= t_k`), and so did `q1`. At that draw the analyst's
increment is `6000`, the probe's jump kernel is `5541`, and `q1`'s rising kernel is **zero**.
Everything else in the record is within the kernels' ordinary disagreement.

`theta_hat/theta`, fixed 3.5 d analyst, peak-matched, 8 noiseless records:

| kernel | as run (post-dose at t = 0) | **pre-dose at t = 0** | t = 0 dropped |
|---|---|---|---|
| probe (jump) | 1.37 | 1.31 | 1.39 |
| `q1` (effect compartment) | **5.81** | **1.14** | 1.23 |

The probe's own 14.5 d schedule has the same kind of draw, at `t = 392`
(`0.5 + 14.5 x 27`); it moves `q1`'s kernel from `2.23` as run to `1.25` pre-dose.

**Why the pre-dose reading is the right one, and not a convenience.** The draw is inside the
marker's own settling window, `0.054 d` (1.3 h) by `g_design_common.settling_days`, so the
quasi-steady model under test does not describe it on either side. A clinic draws before it
infuses. And it is exactly the instant at which two defensible readings of the published
kernel disagree by the whole dose, a jump against a rise, so a result carried by it is a
result about that instant and not about the shape. **[judgement]**, adopted as the S1 footing
(`DOSE_DAY_DRAW = "pre-dose"` in `q1`).

**The same draw moves the exponential truths too, and that was not visible before.**
`q1b` R4, fixed analyst, peak-matched: 7 d is `2.02` as run, `1.38` pre-dose and `1.78` with
the draw dropped; 5 d is `1.36`, `1.17` and `1.29`. A pre-dose draw at `t = 0` is an
**untreated baseline** measurement, which pins `theta` through the untreated balance with no
increment in it at all, and that is why it shrinks the shape bias by more than dropping the
post-dose draw does. **A pre-treatment baseline draw is worth about half the 7 d shape bias.
[measured]** That is a design finding, and the trial already takes pre-treatment draws.

### 4. The free-amplitude cell: an error in `q1` **[measured]**

With the 12% prior on the dose scale, on identical records, the probe's fitter returned
`5.74` and `q1`'s returned `0.11`. The objectives differed: the probe carries the prior on
`ln kappa`, and `q1` carried it on `kappa`, the delta-method form. **The repository's own
convention is `ln kappa`**: the dose-scale column in `g_design_common.sd_ln_theta` is
`inc g`, the derivative in `ln kappa`, and its prior is stated on that coordinate. At
`kappa = 0.01` the linear form charges 68 and the log form 1470, which is the whole of why the
first run could "switch the increment off". **Same class of error as the `k_e` prior defect
caught before the first report.** Fixed; the fitters now agree on every interior seed to
`1.2e-7` (`q9` check 9). The two seeds where they still differ, by `0.041` and `0.223`, are
the two where `q1`'s `K_m` sits on its box bound of 1000 and the probe's runs past it (1041,
1271), which is the domain difference `q9` check 4 already measures and not the prior.

**The first version of `q9` check 9 compared every seed and FAILED at 0.223.** It is now made
on interior seeds and prints the box seeds; this paragraph is the record of the failure.

### 5. A second error, found on the way **[measured]**

Re-running `q9` on the pre-dose footing, check 5 (the profile interval against `n6`'s grid)
failed at `2.0e-2` on one upper edge, `4026` against `4108.6`. With `theta` held high on a
profile, the linear inner solve can return `k_e <= 0`; the polish then started from
`ln(1e-12)`, where `exp(ln k_e)` has no gradient, and was accepted because the start's own
objective is infinite. **The profile jumped from 1.53 to 553 between two adjacent `theta`**,
and the upper interval edge came out short. Fixed by restarting a non-positive `k_e` from the
prior's centre; the two profiles now agree to four decimals out to `theta = 8000`, and check 5
passes at `7.1e-4`. The first report's upper interval edges were affected wherever this
triggered, which makes its coverage figures conservative, not optimistic.

### 6. The wrong-family cells, reconciled **[measured]**

Monthly, 8 noiseless records, single-exponential analyst with the half-life free:

| analyst | kernel | as run | pre-dose |
|---|---|---|---|
| half-life free, dose amplitude fixed (probe p3b; `q1`'s `dose`) | probe | **1.45** (probe 1.46) | 1.31 |
| same | `q1` | 0.00, both fitters | **1.13**, both fitters |
| half-life free, amplitude re-normalised to the steady **peak** (what `q1` Part 1 ran) | probe | 1.51 | 1.30 |
| same | `q1` | 0.00 | **1.11** |
| half-life and amplitude both free, flat (probe p3c) | probe | **1.01** [0.75, 1.83] (probe 0.90 on 6 seeds) | 0.30 [0.18, 0.43] |
| same | `q1` | 0.04 | 0.29 [0.17, 0.42] |

Three things follow. **The "capacity zero on all 24 seeds" of the first report was the `t = 0`
draw**, on both fitters: it goes to `1.13` pre-dose. **`q1`'s Part 1 did not run the probe's
analyst**: it re-normalised each candidate half-life to the steady peak where p3b holds the
per-dose amplitude, and the two differ by a few per cent. And **the both-free analyst is
poorly determined**: its 8-seed range as run is `0.75` to `1.83`, the probe's `0.90` sits
inside it, and the pre-dose baseline draw moves it to `0.30`. It should not be quoted as a
point.

### 7. What the reconciliation withdrew from the first report

| first report said | corrected | status |
|---|---|---|
| Attarwala truth `5.637` [2.440, 8.609], "not in the 1.3 to 2.0 band" | **`1.131` [1.036, 1.366]**, between the 3.5 d and 5 d exponentials | **withdrawn** |
| "the two-compartment truth should stop being grouped with a slightly slower exponential" | it behaves like one, and is grouped | **withdrawn** |
| "audible, `2.35 sigma` on every seed; the paper's under-a-fifth-of-sigma sentence does not describe it" | `1.079` [0.743, 1.524] against the null's `0.981` [0.706, 1.324]: **inaudible** | **withdrawn** |
| free amplitude "switches the increment off, capacity nine times low" | `1.126`, identical to the fixed analyst | **withdrawn**; `q1`'s kappa prior was in the wrong coordinate |
| half-life free on the two-compartment truth, "capacity zero on all 24 seeds" | **`1.129`** [0.886, 1.286] | **withdrawn** |
| two-compartment analyst "biased on every exponential truth, 1.18 to 2.25"; "no analyst family is safe" | **1.05 to 1.14 on every truth** | **withdrawn**; see S1.3 |
| 7 d truth `1.934`, 5 d `1.310`, "confirming the review's 1.3 to 2.0" | **`1.358` and `1.161`** at the pre-dose footing | numbers change; direction stands |
| "the family question wants the hours after an infusion" (S1.8) | true under the effect-compartment reading **only** | **qualified**, S1.8 |
| getting the dose train wrong costs more than the shape, and is audible | still so: Q3W as Q2W reads **6.52** at `2.20 sigma` | **stands** |
| freeing the amplitude does not help on the exponential truths | still so, and now on every truth | **stands, widened** |

---

## S1. The increment time course: analyst against truth (corrected)

### S1.0 The checks **[measured]**

`audit/q9_recheck.py`, **24 of 24 pass**, each by a route sharing no code with what it
checks. The ones that bear on this revision:

- two-compartment kernel by hand from the Laplace transform against `q1`'s eigen-decomposition,
  `< 1e-9` **[proven]**; and against `solve_ivp`, `9.3e-12`.
- the committed kernel and train against `design.run.delta`, `1.4e-16`.
- `q1`'s fit against a from-scratch least squares on the same objective and domain, objective
  gap `3.7e-15`.
- `q1`'s profile interval against `n6`'s objective on a fixed grid, `7.1e-4`.
- check 9, the free-amplitude fit against the probe's fitter, interior seeds, `1.2e-7`.
- check 10, the dose-day draw: pre-dose every kernel is zero at `t = 0`; post-dose the
  committed analyst is `design.run.delta` there, `1.5e-16`.

And the three defects the checks caught **before** the first report, still fixed: the
unpolished profile, the `1.92` cut on a doubled objective, and the polish leaving the `K_m`
box. Three were caught **after** it, and are S1.R 4, S1.R 5 and the dose-day draw.

### S1.1 The headline grid **[measured]**

Nominal 14 d train, peak-matched, pre-dose at `t = 0`, `theta = 3000`, monthly `n = 25`,
24 seeds. `theta_hat/theta` above 1 is capacity read **too high**. `cover` is the fraction of
95% profile intervals containing the truth; with 24 seeds its standard error at 95% is 4.4
points.

| truth | analyst | `theta_hat/theta` median [min, max] | RMSE | RMS resid / sigma | cover |
|---|---|---|---|---|---|
| exp 2.5 d | fixed 3.5 d | **0.926** [0.838, 1.134] | 0.101 | 1.139 [0.86, 1.50] | 100% * |
| exp 2.5 d | half-life free | 0.977 [0.900, 1.168] | 0.067 | 0.974 | 100% |
| exp 2.5 d | two-comp free | 1.140 [0.795, 1.296] | 0.156 | 0.951 | 100% |
| **exp 3.5 d (null)** | **fixed 3.5 d** | **0.986** [0.900, 1.183] | 0.070 | **0.981** [0.71, 1.32] | **100%** |
| exp 3.5 d | half-life free | 0.996 [0.806, 1.162] | 0.081 | 0.978 | 100% |
| exp 3.5 d | two-comp free | 1.072 [0.738, 1.209] | 0.133 | 0.934 | 100% |
| exp 5.0 d | fixed 3.5 d | **1.161** [1.066, 1.377] | 0.171 | 1.139 [0.80, 1.49] | **100%** |
| exp 5.0 d | half-life free | 1.029 [0.768, 1.229] | 0.103 | 0.981 | 100% |
| exp 5.0 d | two-comp free | 1.068 [0.714, 1.239] | 0.141 | 0.937 | 100% |
| exp 7.0 d | fixed 3.5 d | **1.358** [1.267, 1.634] | 0.336 | 1.471 [1.03, 1.86] | 83.3% |
| exp 7.0 d | half-life free | 1.030 [0.765, 1.286] | 0.118 | 0.980 | 100% |
| exp 7.0 d | two-comp free | 1.060 [0.706, 1.314] | 0.144 | 0.942 | 95.8% |
| Attarwala 2-comp | fixed 3.5 d | **1.131** [1.036, 1.366] | 0.155 | 1.079 [0.74, 1.52] | 100% |
| Attarwala 2-comp | half-life free | 1.129 [0.886, 1.286] | 0.143 | 1.047 | 100% |
| Attarwala 2-comp | two-comp free | 1.051 [0.730, 1.210] | 0.145 | 0.931 | 100% |

`*` 5 of 24 intervals open below, so that 100% is partly a failure to bound rather than
precision. Every other 100% is closed intervals throughout.

**Gap B stands in direction and is smaller than the review reported, on this footing.** A
truth slower than the assumption reads capacity high and a faster one low, monotonically:
`0.926` at 2.5 d, `1.161` at 5 d, `1.358` at 7 d. The review's `1.3 to 2.0` for 5 to 7 d
becomes `1.16 to 1.36` once the `t = 0` draw is a baseline rather than a post-dose instant.
S1.6 shows the size is also convention-dependent, and S1.R 3 that it depends on whether the
record has a pre-treatment draw at all.

**The two-compartment truth is not different in kind.** `1.131` sits between the 3.5 d null
and the 5 d truth, at a residual indistinguishable from the null. **By the author's rule (a),
it is grouped with the exponentials from here on**: it behaves like a single exponential of
about 4 d, which is what a single-exponential analyst with its half-life free recovers from it
(`3.61 d` [2.30, 4.89]).

### S1.2 Freeing the amplitude does not help, on any truth **[measured]**

| truth | fixed 3.5 d | + free amplitude, 12% prior on `ln kappa` |
|---|---|---|
| exp 2.5 d | 0.926 | 0.917 |
| exp 3.5 d | 0.986 | 0.985 |
| exp 5.0 d | 1.161 | 1.154 |
| exp 7.0 d | 1.358 | 1.334 |
| Attarwala | 1.131 | 1.126 |

The review's report is confirmed on every truth, and the first report's exception on the
Attarwala truth is withdrawn (S1.R 4). A scale cannot reach a misfit that is in the shape.

### S1.3 Carrying the shape: the richer family is the robust one **[measured]**

**Half-life free** removes the bias inside the single-exponential family, recovering the
half-life each time (`2.41`, `3.45`, `4.84`, `6.73 d` against 2.5, 3.5, 5, 7), and costs
precision modestly: RMSE at the null `0.070` fixed against `0.081` free. On the
two-compartment truth it does **not** help, `1.129` against the fixed analyst's `1.131`,
because no single exponential is that shape.

**The two-compartment family with free rates is within `1.05` and `1.14` on every truth
tested**, the exponentials included, at a residual at or below the null's and coverage of at
least 95.8%. It pays for that robustness: its RMSE is `0.13` to `0.16` everywhere, about
twice the fixed analyst's `0.070` when the fixed analyst happens to be right. **So the
trade is a factor of about two in precision for bias that stays under 14% across a factor of
2.8 in half-life and both family members.** The first report's "no analyst family is safe"
is withdrawn; this is the replacement. **[judgement]** on whether the trade is worth it, which
is S2's question and not S1's.

### S1.4 How loud is the warning **[measured]**

The reference is the measured null, `0.981` [0.706, 1.324], not `sqrt((n-p)/n) = 0.917`,
because the paper's fit is penalised (`q9` check 8: unpenalised `0.95`, penalised `0.97`).

| truth, fixed 3.5 d analyst | RMS resid / sigma | against the null |
|---|---|---|
| exp 5.0 d | 1.139 [0.80, 1.49] | +16%, inside the null's seed range |
| Attarwala | 1.079 [0.74, 1.52] | +10%, inside it |
| exp 2.5 d | 1.139 [0.86, 1.50] | +16%, inside it |
| exp 7.0 d | 1.471 [1.03, 1.86] | +50%; the median clears the null's worst seed |

**A wrong shape is inaudible for every truth within a factor of 1.5 of the assumption, and
the two-compartment truth among them.** Only the 7 d truth is marginally audible. The paper's
`sec:structuredep` sentence, that a misspecified fit's residual "can be under a fifth of one
sigma", **needs no qualification from S1**; the first report's claim that it did is
withdrawn. The worst cell remains the 5 d truth: `1.161` at 100% coverage and a residual
inside the null's range, so nothing visible goes wrong.

### S1.5 The trial's real, irregular dose train **[measured]**

Construction as before, `audit/PHASE6_REVIEWS_SEP2026.md` section 8, mg/kg on the committed
scale through the dose-weighted mean of the 346 published doses, `0.41012` mg/kg; the size of
the dose-2 delay ASSUMED at 7 d and the composite train a construction. **[judgement]**

**(a) An irregular train does not rescue the shape problem**, analyst holding the real train:

| truth | nominal | C5 de-escalated | C1 Q3W then Q2W | composite |
|---|---|---|---|---|
| exp 7.0 d | 1.358 | 1.418 | 1.211 | 1.314 |
| Attarwala | 1.131 | 1.187 | 1.078 | 1.280 |

**(b) Getting the dose train wrong is the larger lever, and it is audible.** Analyst assumes
the nominal Q2W train, correct 3.5 d shape:

| train the truth ran | `theta_hat/theta` | RMS resid / sigma |
|---|---|---|
| C5 de-escalated | 0.855 [0.749, 1.121] | 1.710 [1.36, 2.12] |
| C1 Q3W then Q2W | **6.517** [4.176, 8.033] | 2.203 [1.96, 2.94] |
| composite | **5.109** [4.119, 5.911] | 2.584 [2.17, 3.19] |

Unchanged in substance by the reconciliation. By the author's decision (c) this enters the
specification as a record requirement, S1.9.

### S1.6 The amplitude convention **[measured]**

`theta_hat/theta`, fixed 3.5 d analyst:

| truth | dose | peak | mean |
|---|---|---|---|
| exp 2.5 d | 0.967 | 0.926 | 0.690 |
| exp 3.5 d | 0.986 | 0.986 | 0.986 |
| exp 5.0 d | 1.080 | 1.161 | 1.442 |
| exp 7.0 d | **1.145** | **1.358** | **2.000** |
| Attarwala | 1.046 | 1.131 | 1.282 |

The 7 d bias runs `1.145` to `2.000` across conventions, a factor of 1.75, so any
cross-shape number must name its convention. **`dose` is the convention the Tier 2
calibration measurement pins** (S1.9), and it gives the smallest biases in the table.

### S1.7 The disease range **[measured]**

At `theta = 150`, footing re-solved in `n4_theta_scope`'s convention A:
`L = 322.7273`, `m_bar = 54.0231`, `sigma = 22.1774`, `b_loc = 5.6323`.

| truth, fixed 3.5 d analyst | `theta_hat/theta` | resid | cover |
|---|---|---|---|
| exp 2.5 d | 0 on 24 of 24 | 1.752 | 0%, 24 open below |
| exp 3.5 d (null) | 0.214 [0, 3.37], 11 at floor | 0.973 | 54.2%, 24 open below |
| exp 5.0 d | **6.24** [5.62, 7.23] | 1.582 | **0%**, closed |
| exp 7.0 d | **6.71** [5.75, 7.80] | 2.344 | **0%**, closed |
| Attarwala | **6.09** [5.28, 7.27] | 1.471 | **0%**, closed |

The null reproduces `sec:scope`: no lower bound and a point estimate that is not usable. The
new thing is the slow truths: **the interval closes, excludes the truth on every seed, and
the estimate is six times high.** A slower increment than assumed leaves more of it
unexplained where the enzyme's own flux is a twentieth of the load, and the fit hands it to
`theta`. The residual is audible there, `1.5` to `2.3 sigma`. The booked `D = 6000` is not
admissible at this capacity on readings R2 or R3 of Section 7, so this is a
misspecification measurement at an inadmissible dose and not a design recommendation.

### S1.8 What measurement would pin the shape **[measured]**, with one qualification

Unchanged numerically, since it is computed on steady-state phases, which the dose-day draw
does not touch.

**Estimating the half-life within the single-exponential family** needs two measurements at
the two ends of the accessible window, `sd(ln t_half) = sqrt(2)/(lambda |tau_1 - tau_2|)`
**[proven]**, re-checked to 5e-11. Against the gate that separates 3.5 d from 5 d at two sd
(`<= 0.178`), reading the increment directly needs **2** measurements at any CV up to 30%; a
functional readout (`theta + Delta`) needs **4** at 10% CV and **14** at 20%.

**Telling the families apart depends on which reading of the published kernel is right.**
`q1b` R6 re-ran it on both kernels:

| readout | kernel | earliest phase 1 h | 6 h | 1 d | N at 20% CV, 1 h / 6 h |
|---|---|---|---|---|---|
| increment | effect compartment | sep 1.36 | 0.68 | 0.52 | 2 / 2 |
| increment | jump | sep 0.62 | 0.60 | 0.53 | 2 / 4 |
| functional | effect compartment | sep 0.51 | **0.15** | 0.09 | 4 / **34** |
| functional | jump | sep **0.13** | 0.13 | 0.12 | **44** / 44 |

Under the effect-compartment reading the discriminating signal is the rise in the first hours
after an infusion, and the first report's "the family question wants the hours after an
infusion" holds. **Under the jump reading it does not**: the window barely matters, the best
pair moves to about 3.3 d and the end of the cycle, and a functional readout needs about 44
measurements at 20% CV whatever the window. **The first report's S1.8 conclusion is qualified
to the effect-compartment reading.** The half-life result is not affected, since it is within
one family.

### S1.9 The author's decisions, and what S1 now hands S2

The author answered three questions on 26 September, conditionally on this reconciliation.

**(a) Group the two-compartment truth with the exponentials.** Resolved by S1.1: it lands at
`1.131`, between 3.5 d and 5 d, and is inaudible. Grouped.

**(b) The amplitude convention is the one the calibration measurement pins.** Tier 2 is a
same-day pair of breath tests across one dose, so it measures the **per-dose jump** in
capacity. For an instantaneous-rise kernel that is the delivered dose `D`. **That is `q1`'s
`dose` convention**, which holds the single-dose peak equal to `D_AMP`; for every exponential
kernel the single-dose peak **is** the jump at the infusion. It exists and needs no adding.
One qualification the author should see: under the effect-compartment reading the protein
rises over a few hours, so a breath-test pair spanning one dose measures the jump only if its
second test falls after the rise (`ke0` gives a 2.9 h half-life, so about 12 h after the
infusion is within 6% of it). **[judgement]** `dose` becomes primary from S2 on, with `peak`
and `mean` as the sensitivity column in S4.

**(c) The dose train enters the specification as a record requirement**: administered dates,
times and amounts, including delays, skipped doses and dose changes, used as the regressor,
not a new design variable. S1.5 (b) prices ignoring it: up to `6.5x` on capacity.

**Not started.** The author's instruction was to proceed to S2 only if the reconciliation left
S1's conclusions unchanged. It did not (S1.R 7), so S2 has not been started.

### A flag for the author, not a finding

The paper's own headline design numbers are computed on the same record, with the `t = 0` draw
read post-dose inside the settling window. `q1b` R5 prices it through
`g_design_common.sd_over_seeds`: **`0.1699`** [0.1684, 0.1709] as committed, **`0.1695`** with
that draw read pre-dose, `0.1721` with it dropped. **The headline does not depend on it.** It
matters for misspecification, not for the well-specified bound. No text proposed; `paper/` is
untouched.

---

## What is open after S1

- **The shape itself**, and **which reading of Attarwala's kernel is right**. Nothing here
  decides either. Brief **R2** finds no human hepatic time course. S1.8 says the second
  question matters for how the first would be measured.
- **The measurement CV on an increment reading**, swept, and pinned only by **R2**.
- **The size of the dose-2 delay** in the C5 train, not published; held at 7 d.
- **Whether to carry the two-compartment family as the analyst.** S1.3 says it buys bias
  under 14% on every truth for about twice the RMSE. That is S2's to price against the
  gate.
- **S2 to S5** have not been run.
