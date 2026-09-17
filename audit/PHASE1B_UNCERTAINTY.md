# Phase 1, step 2. Uncertainty

Two scripts. **No paper edits, no merge, no changes to existing constants or printed
numbers.** `PHASE4_PROPOSED_TEXT.md` stays unapplied.

| | what it does | output |
|---|---|---|
| `audit/j1_uncertainty.py` | four interval methods on a small grid of conditions | tables |
| `audit/j2_recheck.py` | four of its claims, by routes sharing no code path | tables |

**Footing, once.** Local model constants as `src/common/params.py` carries them.
**Design: the paper's monthly record `np.arange(0, 730, 30.4)`, `n = 25`, and nothing else.**
Ordinate noise sd `sigma_total(0.05, 0.05) = 21.9612`. Parameters in the paper's coordinates.
`b` is the **balance's** `b`, an unrecorded addition to the load, so `b >= 0` says the diet
record cannot over-state intake. **The main setting draws the nuisances from their priors**;
`at prior means` is a reference row only. **"Load prior off" means the load-bias prior only**;
the `k_e` log-prior of 0.8 always stays, because `i1`'s ladder shows flattening it costs a
factor of 23 and the fully prior-free design is not identified. 30 replicates x 8 seeds per
row, 199 bootstrap resamples per replicate, load record fixed within a seed. Every figure is
the median over seeds with the min-max range.

**Four methods.** `filter` = Wald, `theta_hat +/- 1.96 sd` from the inverse information
matrix at the fit. `profile` = `{theta : J(theta) - J_min <= 3.841}`. `calibrated` = the same
set with the threshold replaced by the 95th percentile of a parametric-bootstrap LR
distribution (test inversion). `percentile` = the 2.5/97.5 percentiles of `ln theta_hat*`
over the same bootstrap.

**A diagnostic this study needed and did not start with.** An interval that runs off the end
of the `ln theta` grid is not an interval of that width, it is an **unbounded** one truncated
by the grid, and reporting its width as the grid's reads a failure as a result. Every table
below flags the fraction of censored intervals, and a coverage of 1.000 next to a censoring
flag means nothing.

---

## (A) The two parametric bootstraps

Load CV 10%, load prior off, dose exact, no bound, nuisances drawn. **The only difference
between the blocks is what the bootstrap holds fixed when it resimulates.**

| method | nuisances FIXED at estimates | | nuisances REDRAWN from priors | |
|---|---|---|---|---|
| | coverage | width | coverage | width |
| filter | 0.850 | 0.6471 | 0.850 | 0.6471 |
| profile | 0.900 | 0.7137 | 0.900 | 0.7137 |
| calibrated | **0.683** | 0.3440 | **0.967** | 1.2681 (23% censored) |
| percentile | **0.683** | 0.3563 | 0.817 | 0.7892 |
| bootstrap LR 95% | **1.000** (0.857-1.224) | | **8.301** (8.058-10.402) | |

*`filter` and `profile` do not use the bootstrap, so they are identical between the blocks.*

**The fixed-nuisance parametric bootstrap is not usable here.** It halves the width and
delivers 68% coverage for a nominal 95%. The reason is visible in the LR row: resimulating at
`b_hat` and `k_e_hat` throws away exactly the variation the nuisances contribute, so the LR
distribution collapses to a 95th percentile of 1.0 against a `chi^2_1` of 3.841.

**Redrawing from the priors over-corrects**, to 8.30, and buys nominal coverage at 1.8x the
width, with a quarter of the intervals unbounded. It is the less wrong of the two and it is
not right.

---

## (B) The main grid

All rows: nuisances drawn, redrawn bootstrap, `n = 25`.

| load CV | load prior | dose | filter | profile | calibrated | percentile | empirical LR 95% |
|---|---|---|---|---|---|---|---|
| 2% | off | exact | 0.850 / 0.995 | **1.000** / 1.93 **(78% cens.)** | 1.000 / 2.25 (93% cens.) | 0.949 / 1.82 | 2.04 (1.75-5.42) |
| 2% | off | CV 11.7% | 0.917 / 1.075 | **1.000** / 1.96 **(78% cens.)** | 0.983 / 2.12 (77% cens.) | 0.966 / 1.81 | 1.90 (1.26-4.83) |
| 2% | **on** | exact | 0.883 / 0.654 | **0.900** / 0.679 | 0.950 / 1.051 | 0.767 / 0.757 | 6.92 (3.17-29.8) |
| 2% | **on** | CV 11.7% | 0.900 / 0.776 | **0.917** / 0.791 | 0.917 / 0.939 | 0.767 / 0.758 | 4.63 (2.13-15.6) |
| 10% | off | exact | 0.850 / 0.647 | **0.900** / 0.714 | 0.967 / 1.268 (23% cens.) | 0.817 / 0.789 | 5.36 (1.43-16.0) |
| 10% | off | CV 11.7% | 0.917 / 0.770 | **0.933** / 0.817 | 0.933 / 1.025 (12% cens.) | 0.817 / 0.789 | 4.07 (1.09-14.4) |
| 10% | **on** | exact | 0.883 / 0.654 | **0.900** / 0.676 | 0.950 / 1.059 | 0.783 / 0.754 | 6.92 (3.01-31.8) |
| 10% | **on** | CV 11.7% | 0.917 / 0.776 | **0.917** / 0.788 | 0.933 / 0.936 | 0.783 / 0.751 | 4.63 (2.03-16.2) |

*cells are `coverage / median width`; nominal coverage 0.95. Bias +0.069 and RMSE 0.194-0.201
in every row with the load prior on; bias -0.044 and RMSE 0.425 in the CV 2% prior-off rows.*

**Five results.**

**(a) `chi^2_1` is the wrong threshold, and by a factor that is not constant.** The empirical
95th percentile of the LR at the **true** `theta` runs **4.07 to 6.92** across the usable
rows, against 3.841, and its 8-seed range reaches 31.8. So the profile interval is
systematically a little too narrow, which is exactly the 0.90-0.93 coverage it delivers. **No
fixed threshold is right**, because the right one varies by a factor of ten across load
records.

**(b) The profile is nevertheless the best of the four**, in every usable row. It is never
worse than 0.900, never censored except where the design itself has failed, and needs no
bootstrap.

**(c) The Wald interval is 5 to 7 points worse than the profile at the same footing**, and
the gap is asymmetry: the profile is not symmetric about `theta_hat`, so a `+/- 1.96 sd`
summary of it loses coverage. The gap closes when the dose-scale prior is added (0.917 in
both), because that widens the profile and makes it more nearly quadratic.

**(d) The percentile bootstrap should not be used.** 0.767 to 0.817 in every row with the
load prior on, at a width close to the profile's. It is the only method that is both narrow
and wrong.

**(e) `load CV 2%` with the load prior off is a broken design, not a hard one.** 78% of
profile intervals are unbounded. Its coverage of 1.000 is that censoring, and its RMSE of
0.425 is twice any other row's. This is the same design that `i1`'s ladder put at 0.1723 on
the paper's footing; removing the load prior at a load CV of 2% is what breaks it.

---

## (C) The reference row: nuisances at their prior means

Load CV 10%, load prior off, dose exact, `b = 0` and `k_e = k_e0` exactly.

| filter | profile | calibrated | percentile | bias | RMSE | empirical LR 95% |
|---|---|---|---|---|---|---|
| 1.000 / 0.689 | 1.000 / 0.744 | 1.000 / 1.42 (20% cens.) | 1.000 / 0.828 | 0.0048 | **0.0934** | **0.879** (0.553-1.258) |

**Everything over-covers, at every method.** The empirical LR 95% is 0.879 against 3.841, so
the `chi^2_1` threshold is more than four times too large here. A patient sitting at the
prior mean gets an interval that is roughly twice as wide as it needs to be, and an RMSE of
0.093 against 0.194 in the drawn case. This is the other half of `i1`'s finding that 0.1697
is a prior-average: **the same interval recipe over-covers the typical patient and
under-covers the population.**

---

## (D) The `b >= 0` bound

Load CV 10%, load prior off, dose exact, `k_e` still drawn.

| `b_true` | bound | filter | profile | calibrated | percentile | bias | RMSE | on the bound | empirical LR 95% |
|---|---|---|---|---|---|---|---|---|---|
| 0 | off | 0.900 / 0.686 | **0.950** / 0.759 | 0.967 / 1.356 | 0.867 / 0.824 | 0.029 | 0.188 | 0.0% | 4.67 (2.29-8.54) |
| 0 | **on** | 0.883 / 0.638 | 0.900 / 0.661 | 0.967 / 1.065 | 0.850 / 0.720 | 0.051 | 0.179 | **38.3%** (13-47%) | 5.38 (2.12-12.8) |
| 16.14 | off | 0.900 / 0.684 | **0.933** / 0.755 | 0.967 / 1.378 | 0.867 / 0.827 | 0.029 | 0.187 | 0.0% | 4.67 (2.28-8.53) |
| 16.14 | **on** | 0.883 / 0.656 | 0.933 / 0.690 | 0.967 / 1.139 | 0.833 / 0.726 | 0.046 | 0.178 | **26.7%** (7-43%) | 5.09 (2.09-11.5) |

**The bound buys RMSE and costs coverage.** It narrows the profile by 9% at `b_true = 0` and
takes coverage from 0.950 to 0.900; RMSE improves from 0.188 to 0.179. That is the classic
boundary trade and it is the classic boundary problem: **with `b_true = 0`, 38% of fits land
exactly on the bound**, where `chi^2_1` is not the reference distribution for any of these
methods. Even at `b_true` one prior sd away, 27% still land on it.

**It is not a large effect and it is not free.** If the bound is imposed, the coverage
shortfall has to be stated with it.

---

## (E) The constant-load case

| load CV 0% | filter | profile | calibrated | percentile | bias | RMSE | `K_m` on a grid edge | empirical LR 95% |
|---|---|---|---|---|---|---|---|---|
| load prior **on** | 0.883 / 0.653 | 0.900 / 0.678 | 0.950 / 1.047 | 0.767 / 0.752 | 0.069 | 0.194 | 0.0% | 6.94 (3.23-29.4) |
| load prior **off** | 1.000 / **928.85** (382-4059) | 1.000 / 3.0000 **(100% censored)** | 1.000 / 3.0000 (98% cens.) | 1.000 / 3.52 | 0.507 | **1.331** | **89.3%** | 0.629 (0.291-1.112) |

**With the load prior on, a constant load behaves like any other row**, because the prior is
supplying all the curvature (`i1` part D proved the data's own profile is flat to `1.5e-19`).

**With it off, every method fails and three of them fail while reporting 100% coverage.** The
Wald width is 929 in `ln theta` and reaches 4059. The profile and calibrated widths of exactly
3.0000 are the width of my own `ln theta` grid: **the intervals are unbounded**. 89% of fits
sit on a `K_m` grid edge. And the empirical LR 95% is **0.629**, comfortably *below* 3.841, so
a reader checking only that diagnostic would conclude the threshold was conservative.

**That is the single most important row in this study.** A constant-load record with no load
prior produces an estimate, an interval, a coverage of 1.000 and an LR diagnostic that all
look acceptable, and every one of them is an artefact of where the computation was truncated.

---

## The re-check (`j2`)

| # | claim | route | result |
|---|---|---|---|
| 1 | the Wald sd and the profile describe one object | finite-difference the profile's curvature instead of inverting a matrix | ratios 0.925 to 0.983 across six replicates |
| 2 | the LR statistic is correctly formed | test it where its law is **exact**: `K_m` held, all priors flat, model linear | mean **1.0138**, 95% **3.8844** against `chi^2_1`'s 1 and 3.8415, quantile ratios within 5% at 20000 replicates; and `j1`'s solver agrees with a from-scratch projection computation to **9.4e-13** |
| 3 | the `b >= 0` active set is the constrained optimum | SLSQP, which knows nothing about active sets | worst objective difference **1.4e-14**, and the on-bound cases agree exactly |
| 4 | the profile intervals are the set they claim | direct evaluation on a grid ten times finer, no interpolation, no walking | endpoints agree to within the fine grid's own spacing of 1.9e-3; every accepted set connected |

**Check 2 failed twice before it worked, and both failures were mine.** The first draft
asserted `chi^2_1` for the full nonlinear problem, which has a condition number of `3.7e7`
and no reason to obey it. The second said "priors off" while the `k_e` prior was hard-wired
on in `build`, which is what shrank the statistic to 0.29 of `chi^2_1` and looked like a bug
in the statistic. The prior was the bug, in the check. With the `k_e` prior genuinely flat and
`K_m` held, the statistic sits on `chi^2_1`.

That sequence is also the answer to the question the check was asked: **the departures `j1`
reports are properties of the penalised, nuisance-carrying problem, not of the arithmetic.**
`j2` check 2 shows the same code giving mean LR **1.0138** (linear, flat priors), **1.0133**
(`K_m` free, flat priors) and **0.3212** (`K_m` free, `k_e` prior on). The prior is what
moves it, and one nonlinear direction costs nothing at all.

*(Reported faithfully: the Kolmogorov-Smirnov p against `chi^2_1` is 0.0066 in the linear
block. At 20000 replicates that test resolves a 1% departure, and the 1.4% excess in the mean
is 1.4 standard errors. The check is that the statistic sits on `chi^2_1` to within a few per
cent, not that it is indistinguishable from it at any sample size.)*

---

## The recommendation

**Report the profile interval with the `chi^2_1` threshold, and report it as approximate.**
It is the best of the four in every usable row, it needs no bootstrap, it degrades visibly
rather than silently, and its failure mode is a 2 to 5 point coverage shortfall rather than a
27 point one.

**And state these four conditions with it, because each one is load-bearing:**

1. **It is conditional on the priors, and the paper should say which.** The `k_e` log-prior of
   0.8 is doing the identification. The load-bias prior matters at `n = 12` and at a constant
   load, and costs 1.5% at `n = 25` (`i1`'s ladder).
2. **It covers 90-93%, not 95%.** The empirical LR 95% is 4.1 to 6.9 against 3.841. If a
   stated 95% matters, widen the threshold to the empirical value on the design in hand; the
   paper's own `sd_of` machinery can compute it.
3. **It is a statement about the population the prior describes, not about the patient in
   front of you.** At the prior mean the same recipe over-covers (1.000) with an RMSE of
   0.093; with the nuisances drawn it under-covers (0.900) with an RMSE of 0.194.
4. **It is void at a constant recorded load without the load prior.** There the method returns
   an unbounded interval that looks finite, and no diagnostic in this study catches it except
   the censoring flag.

**Do not report the percentile bootstrap** (0.77-0.82 everywhere). **Do not report a
parametric bootstrap resimulated at the fitted nuisances** (0.683). **Do not report the Wald
interval** where the profile is available: it is the same object linearised, and the
linearisation costs 5 to 7 points of coverage for a 5% saving in width.

**If a bootstrap must be used**, redraw the nuisances from their priors. It reaches nominal
coverage, at 1.8x the width and with up to a quarter of the intervals unbounded, which is a
high price for a number that the profile already supplies.

## Addendum: test inversion, the boundary cutoff, seeds and Monte Carlo error

`audit/j3_testinversion.py`, `audit/j4_recheck.py`.

### What `calibrated` was, and what test inversion is

**They are not the same thing, and the first pass ran only the first.** `j1`'s `calibrated`
bootstraps at `theta_hat`, forms `LR* = J*(theta_hat) - J*_min` once, and uses that single 95th
percentile at **every** `theta_0`. Proper test inversion simulates under the **constrained**
fit at each `theta_0` separately, so the cutoff `q(theta_0)` varies along the grid. On one
replicate it runs from **1.30 to 13.00**, which is why one threshold cannot serve. It is also
orthogonal to `fixed vs redrawn`, which is about the nuisances.

`n = 25` monthly, nuisances drawn, dose exact, 21 `theta_0` nodes x 99 datasets each,
160 replicates = 8 diet records x 20. Nominal 0.95.

| row | profile (3.841) | `j1` calibrated | **test inversion** |
|---|---|---|---|
| CV 10%, prior off | 0.881 +/- 0.026, w 0.708, 2% cens | 0.950 +/- 0.017, w 1.268, 19% | **0.931 +/- 0.020, w 0.976, 5%** |
| CV 10%, prior on | 0.863 +/- 0.027, w 0.677, 0% | 0.931 +/- 0.020, w 1.076, 2% | **0.944 +/- 0.018, w 1.471, 51%** |
| CV 2%, prior off | 0.975 +/- 0.012, w 1.948, 75% | 0.975 +/- 0.012, w 2.256, 89% | **0.930 +/- 0.020, w 1.755, 68%** |

**Test inversion does better on the headline number and the re-check shows the number is not
earned.** `j4` check 1 puts its simulated threshold at `theta_true` against the empirical 95th
percentile of the LR there, measured over 200 replicates per diet record: the ratio is **0.255
to 0.427**. The textbook fixed-nuisance form calibrates against a threshold two and a half to
four times too small, exactly as `j1`'s fixed bootstrap did. Its coverage comes from `q`
growing steeply in the tails, and from intervals that are unbounded 5% to 68% of the time.

**Redrawing the nuisances inside the inversion largely fixes the calibration**, moving the
ratio to **1.11 to 1.77**. So the nuisance handling governs calibration in test inversion just
as it did in the simple bootstrap, and the constrained fit is not the right thing to hold
fixed here.

**Recommendation, revised.** Still report the profile with `chi^2_1`, and now state its
measured coverage. Test inversion buys 5 to 8 points of coverage and is more *consistent*
across rows (0.930 to 0.950 against the profile's 0.863 to 0.975), but it pays with intervals
that have no upper bound in a large minority of cases and a threshold whose own calibration is
off by a factor of three in the form usually written down. **Use it as a check on the profile,
not as the reported statement.** If a stated 95% is required, the honest route is to widen the
profile's threshold to the empirical LR 95% on the design in hand, which the paper's own
machinery can compute.

### The boundary cutoff

`b >= 0` imposed, `k_e` drawn, load CV 10%, prior off. 160 replicates.

| `b_true` | on the bound | profile 3.841 | **mixture 2.706** | test inversion |
|---|---|---|---|---|
| 0 (on the bound) | 40.0% (10-55%) | **0.938 +/- 0.019**, w 0.656 | 0.887 +/- 0.025, w 0.545 | 0.950, w 1.478, 58% cens |
| 16.14 (one prior sd) | 25.0% (5-50%) | **0.938 +/- 0.019**, w 0.683 | 0.900 +/- 0.024, w 0.569 | 0.950, w 1.525, 70% cens |

**The boundary correction makes coverage worse, and it should.** The 50:50 mixture of `chi^2_0`
and `chi^2_1` is the right law when the parameter **being tested** sits on a boundary. Here the
parameter on the boundary is a **nuisance** and `theta` is interior, so narrowing the cutoff
from 3.841 to 2.706 just narrows the interval and loses five points of coverage. `j4` check 2
confirms the constant itself is right (2.7055 is the `chi^2_1` 90th percentile, and 2 million
simulated draws give `P = 0.94998`). **Do not apply it here.**

### What a seed varies, and how many diet patterns reach the target

**A seed is a diet record.** The load realisation is drawn once per seed and held; the
replicates inside it vary the observation noise and the nuisances. `j4` check 3 confirms it on
the arrays. So every seed range quoted in Phase 1 is variation **across diet patterns**, a
design property, and Monte Carlo error is reported separately as a binomial standard error.

*(A "load CV 10%" record delivers a realised CV of 0.078 to 0.084 on 25 points, which is part
of why diet patterns differ.)*

**32 diet records, 60 replicates each. Fraction reaching RMSE <= 0.20:**

| setting | median RMSE | range | fraction |
|---|---|---|---|
| CV 10%, prior off, dose exact | 0.2001 | 0.1381-0.2692 | **47% of 32** |
| CV 10%, prior **on**, dose exact | 0.1922 | 0.1362-0.2557 | **59% of 32** |
| CV 2%, prior off, dose exact | 0.5015 | 0.2882-0.9679 | **0% of 32** |
| CV 2%, prior **on**, dose exact | 0.1916 | 0.1336-0.2536 | **59% of 32** |
| CV 10%, prior off, dose CV 11.7% | 0.2000 | 0.1385-0.2692 | **50% of 32** |

**On the paper's own footing the design reaches its own target for 59% of diet patterns, not
for a typical one.** The 8-seed ranges quoted elsewhere in Phase 1 are the same quantity at a
quarter of the resolution and understate how often it misses.

### Monte Carlo error

Every coverage figure above rests on **160 replicates** (8 diet records x 20), binomial
standard error **+/- 0.017 to 0.027** at these coverages. The main `j1` grid rests on 240
(8 x 30), `+/- 0.014` to `0.019`. Those standard errors treat the replicates as independent;
they are not independent across diet records, so the per-seed range beside each figure is the
honest measure of how much a different patient's diet would move it.

---

## Status

`j1` to `j4` exit 0. No paper file touched, no constant changed, no printed number in the
paper altered.
