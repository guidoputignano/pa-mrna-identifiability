# Summary of the audit

Marks: **[proven]** mathematics, exact or symbolic. **[measured]** computed, footing stated.
**[judgement]** a call for the author. Detail in `PHASE4_PI.md`, `PHASE4_STRESS.md`,
`PHASE3_DESIGN.md`, `PHASE2F_ADMISSIBILITY.md`, `PHASE1A_ESTIMATOR.md`,
`PHASE1B_UNCERTAINTY.md`.

**The paper has not been edited.** Proposed wording is held in `PHASE4_PROPOSED_TEXT.md`,
eight passages, none written into `paper/`.

## 1. What is proven

**The constant-load degeneracy is an exact curve, not curvature.** [proven] Solving the
balance for the increment that produces a level `m` gives, with no remainder,
`Delta(m) = A/m + B + C m`, `A = K_m(L+b)`, `B = L+b-theta-k_e K_m`, `C = -k_e`. Three numbers
are all a constant-load record contains, so four unknowns leave a curve:
`K_m(b) = K_0 L/(L+b)`, `theta(b) = theta_0 + b + k_e(K_0-K_m(b))`, `k_e` fixed. Checked
symbolically, in exact rationals, and by the data-free identity
`residual_on_curve/residual_true = (K_0+m)/(K+m)`. A mutation-checked test is in the suite.

**The estimator returns that curve.** [proven] A variable-projection fit has one nonlinear
parameter; its objective is flat to `1.5e-19` across three decades of `K_m`, the fitted
parameters landing on the curve to `1.8e-8`.

**The linear reparametrisation is right, and is the PI's idea.** [proven] Multiplying by
`(K_m+m)` gives a form linear in `beta = (K_m, K_m b, b-theta-k_e K_m, -k_e)`. Symbolic
remainder zero.

**A visit's dose-cycle phase fixes its weekday, permanently.** [proven] 14 days is two whole
weeks. Checked against a calendar over 60 cycles: none of 56 half-day phases drifts.

**The capacity ceiling holds on both readings.** [proven] The used fraction is
`f(1+Delta_peak/theta)`, so `theta` cancels: exactly `47/150` at 10%, `47/300` at 5%.

## 2. What changes in the paper

| # | where | what |
|---|---|---|
| 1-2 | `06_excitation.tex:71-75`, `:497-516` | the `sqrt(noise)` price goes; the degeneracy is an exact curve [proven] |
| 3-4 | `07_design.tex:141`, `:417-419` | the mechanism clause; `rem:no_delta` [proven] |
| 5 | `06_excitation.tex:371` | the elimination drops the bias, so `eq:nullvec` is a `b = 0` value [proven] |
| 6 | inside `rem:rankvsprec` | a sentence false at two of its own four rows [measured] |
| 7 | `07_design.tex:25`, `A1_supplementary.tex:516` | "a `K_m` prior cannot help" has an inverted mechanism [proven] |
| 8 | `params.py`, the `D` row | "falls by about half" is 62-66% on every reading [measured] |

**New from Phase 1: how an interval on `theta` is stated.** [measured] Four conditions must sit
beside it, and each is load-bearing. See section 6.

## 3. What is conditional

**The diet record must carry real day-to-day variation.** [measured] It is the *recorded* `L_k`
that enters, so a record imputed from a prescription carries none. Without the load prior a
constant load is *exactly singular*, and a 2% recorded CV cannot reach the target in 60 visits.

**The priors do the identification, and `k_e` does most of it.** [measured] Monthly `n = 25`,
CV 10%, dose exact: paper's footing **0.1698**; load prior off **0.1723**; `k_e` flat
**3.9570**; both off **5.7432**. This *confirms* the stated footing rather than undermining it.

**The dose scale must be measured, and has not been.** [measured] The required CV, recomputed
over 8 diet records, is **0.1168** (0.1141-0.1187); printed 0.1167. At that CV monthly `n = 25`
reaches 0.20 **with no margin**, and Phase 3's 12-visit calendar needs **36 visits**. The
repository's own between-visit row is *mean 14.0%, range 0.4-50.8%*, an upper bound on what is
needed.

**The fortnightly result needs therapy to start at `t = 0`.** [measured] Already on therapy,
its `Delta/theta` range is **exactly zero** and it is exactly singular at a constant load.

**The weekday calendar needs exact 14-day dosing at a fixed hour.** [measured] Slips move
phases by up to a tenth of a cycle; precision barely follows, the map does not survive, so the
calendar must be re-derived from the actual infusion dates.

## 4. What stays unchanged

**`0.1697` stands on its stated footing.** [measured] A posterior sd: 0.0875 at the prior
means, 0.1796 with the nuisances drawn from those priors, RMSE **0.1922** across 32 diet
records (0.1362-0.2557).

**`+/- 36.9` is right as computed.** [measured] My earlier contrary claim is **withdrawn**:
`chi^2/eps^4` is parametrisation-dependent and the paper uses log parameters. Redrawn in them,
two-daily gives `2.0900e-6` and `36.82` against `2.08e-6` and `36.9`.

**Putting the assay error on the ordinate is vindicated.** [measured] The errors-in-variables
bias at the paper's assay CV is `+0.000185` predicted, 0.11% of the reported sd, and
`+0.0015 +/- 0.0010` over 12000 paired replicates, not resolved from zero.

**The three settling times are all correct**, differing by operating point (4.2 vs 3.71 h) and
by displacement (3.71 vs 1.30 h); the design needs 1.30 h. The exclusion costs 0.39% of the
dial and **10.1% of the reachable `Delta/theta` range**, because it removes the peak.
**Every admissibility constraint holds**, at both capacity readings. [measured]

## 5. Open questions for the author

**(a) Does the exact curve go in, and where?** [judgement] Passages 1, 2, 5 are drafted.

**(b) 5% or 10% of normal capacity?** [judgement] Changes no pass or fail; changes the
admissible range of the unidentified direction by 2.5x, `[-272.7, 20107.9]` to
`[-272.7, 50103.2]`. Needed only if that range is quoted as a number.

**(c) Should the phase clock be the main design figure?** [judgement] 14 d and 14.5 d are
identical on a calendar axis and differ twentyfold in the answer. Folded onto the dose cycle,
one is a single spoke and the others fill the dial. `fig_E1` is drawn and unused.

**(d) Reword the `D` row, or leave it?** [judgement] Nothing downstream breaks if it stays.

## 6. The one recommendation on uncertainty

**Report the profile interval with the `chi^2_1` threshold, and state its measured coverage
rather than calling it 95%.** [measured] Over 160 replicates per row (8 diet records x 20,
binomial se `+/- 0.02`) it covers **0.881** at CV 10% with the load prior off and **0.863**
with it on, against a nominal 0.95. It is never censored there and needs no bootstrap.

Four things beside it: [measured] conditional on the stated priors; **covers about 86 to 88%
at 10% recorded load variation, not 95%**; describes the prior's population, not one patient
(RMSE 0.093 at the prior means against 0.192 drawn); **void at a constant load without the load
prior**, where every method returns an unbounded interval while reporting 100% coverage.

**Read no coverage figure without its censoring column.** The profile's 0.975 at a 2% recorded
CV with the load prior off is not good coverage: about three-quarters of those intervals are
unbounded (78% on the main grid, 75% in the addendum's smaller run), so the hits are intervals
with no upper end. That row is a failed design, not a well-covered one.

**Do not use** the percentile bootstrap (0.77-0.82) or a parametric bootstrap at the fitted
nuisances (0.683). **Do not apply the boundary correction** [proven]: the 50:50
`chi^2_0`/`chi^2_1` mixture is for a parameter *being tested* on a boundary, and here that
parameter is a nuisance; applying it takes coverage from 0.938 to 0.887.

**On test inversion**, which the brief asked for and the first pass did not run: it covers
0.930-0.950 against the profile's 0.863-0.975, better and more consistent. **I still would not
report it** [judgement]: its intervals have no upper bound in 5-68% of cases, and its threshold
at the true `theta` is **0.26-0.43** of the empirical one in the textbook fixed-nuisance form,
so the coverage is not earned by calibration. Redrawing the nuisances inside it moves that
ratio to 1.11-1.77 and is the better variant. Use it to check the profile, not to replace it.

**If a stated 95% is required**, widen the threshold to the empirical LR 95% on the design in
hand: 4.1 to 6.9 here, not 3.84. [measured]

## 7. The design answer, in one line

**On the paper's own footing the monthly `n = 25` design reaches its 0.20 target for 59% of 32
SIMULATED diet records drawn at a 10% load CV** (median RMSE 0.1922, range 0.1362-0.2557).
[measured] At a 2% recorded CV without the load prior, **none of 32**.

**That 59% is conditional on the assumed 10% variation and is not a property of real
patients.** No real diet record enters anywhere in this audit; the load realisations are
lognormal draws at a CV the repository assumes. What the figure establishes is that even under
that assumption the design misses its own target for four records in ten, which is a weaker
claim than "it reaches 0.20" and a stronger one than any single median. The binding
requirements, in order: the `k_e` prior, real recorded diet variation, a measured dose scale.
