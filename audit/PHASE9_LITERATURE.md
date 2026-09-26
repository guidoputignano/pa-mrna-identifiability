# Phase 9. The 25 September literature reports, checked

Marks, as `PHASE4_SUMMARY.md` uses them: **[proven]** mathematics, exact or symbolic.
**[measured]** computed, footing stated. **[judgement]** a call for the author.

The reports themselves are in `refs/review-20260925/`, verbatim, with their mapping and three
process warnings in that directory's `README.md`. **Nothing in them is a finding of this
repository and nothing here moves a constant.** `src/common/params.py` is untouched, no file under
`paper/` is edited, and the only computation this file reports is rebound in a context manager or
passed as an argument.

Independent re-check: `audit/q10_phase9_recheck.py`, 7 checks, all passing.

---

## 0. Egress is blocked, so items 1 to 5 have not been attempted

The verification of the reports' quantitative claims against primary sources needs outbound access
to publisher and index records. It is not available in this container:

```
api.crossref.org:443        CONNECT rejected, gateway answered 403
eutils.ncbi.nlm.nih.gov:443 CONNECT rejected, gateway answered 403
```

Both are policy denials by the environment's network settings, confirmed in the proxy's own
failure log, not transport errors and not rate limits. The person running the session changes this
under Network access in the environment's settings, either by raising the access level or by adding
those two hosts to the allowed domains.

**So step 2 items 1 to 5 are not started, and this file records no verdicts for them.** There is no
table of `NOT CHECKED` rows here, because a step that has not happened is not a step with null
results. `NOT CHECKED` and `NOT FOUND` are different claims and neither is asserted for any of:
R3's two variability surrogates, R5's transplant ratios, R1's propionate turnover and
kidney-disease prevalence, R1's OAT1 substrate claim, and R4's Spearman r.

**What is below needs no network.** Item 6 is arithmetic on figures both reports quote, and item 7
is a decomposition whose only computed input is this repository's own model.

---

## 6. The breath-test conflict: the repository's retraction stands **[measured]**

### The conflict

`REVIEW_LOG.md` R4 retracted the 12.8% breath-test CV on the grounds that the abstract reported a
standard deviation as a mean. The two reports disagree about that.

**R3** (`wf-f8c2cfcf`, §1) quotes the abstract verbatim, "Test repeatability was robust in controls
and in MMA subjects (mean coefficient of variation 6.9% and 12.8%, respectively)", and reads it at
face value: the 12.8% is "a mean of per-subject CVs; the ±12.86% quoted elsewhere in the paper is
the SD of that distribution". Its pending-reference table repeats this as "breath CV is mean of
per-subject CVs", marked VERIFIED.

**R4** (`wf-71022447`, §§ on repeatability) quotes the Results instead: "the CV for
1-13C-propionate dose metabolized in 15 mut patients, who had repeated breath testing 2-3 times
each, ranged from 0.39% to 50.84% with an average of 14.02 ± 12.86%", and separately gives controls
as "6.94 ± 5.3% ... abstract rounds these to 6.9% and 12.8%".

### The arithmetic, which settles it

| comparison | gap | is it a rounding |
|---|---|---|
| controls: abstract 6.9 against Results mean 6.94 | 0.04 | **yes** |
| controls: abstract 6.9 against Results SD 5.3 | 1.60 | no |
| patients: abstract 12.8 against Results mean 14.02 | 1.22 | no |
| patients: abstract 12.8 against Results SD 12.86 | 0.06 | **yes** |

1. **R3's reading is arithmetically impossible.** If the abstract's 12.8% were the mean of
   per-subject CVs, it would equal the Results' 14.02, and it does not. It equals 12.86, the
   standard deviation. **[measured]**
2. **The decisive evidence is stronger than either report noticed: the abstract's two figures come
   from different columns of the same table.** The controls' 6.9 is their mean, 6.94, and the
   patients' 12.8 is their SD, 12.86. A reader who assumes the abstract is internally consistent
   will take both as means, which is exactly the error the retraction identified, and the
   inconsistency is visible without leaving the paper. **[measured]**
3. **So `REVIEW_LOG.md` R4 stands, and there is no retraction of a retraction.** R3 is wrong on
   this point, and wrong in a specific and instructive way: it trusted a label in an abstract over
   numbers in a Results section. That is the same failure mode the September pack produced and the
   reason the retraction exists. **[measured]**
4. **R4 is not right either, it is merely less wrong.** Its own sentence "abstract rounds these to
   6.9% and 12.8%" asserts that both abstract figures are roundings of the two means it just
   quoted, when 14.02 does not round to 12.8. It reproduces the correct numbers and draws no
   conclusion from their mismatch. **[measured]**

### What remains open on item 6

Everything above rests on two quotations neither of which this repository has verified: R4's
Results sentence and R4's controls figure. **The adjudication is therefore conditional on those
quotes being accurate, and it must be confirmed at the source when egress allows.** If R4
misquoted the Results, the whole argument changes. What does not depend on the quotes is the
direction of the standing retraction, which no reading here overturns. **[judgement]** on how much
weight the conditional carries; the arithmetic itself is **[measured]**.

Neither figure may be used as a marker CV in any case, for the reason R3 gives and R4 confirms
independently: it is a breath-test repeatability CV across visits months to years apart, in MMA,
not a within-patient CV of any plasma marker.

---

## 7. The surrogates cannot be carried into sigma unchanged **[measured]**, and cannot be decomposed **[judgement]**

### What the brief got wrong, and why it matters

The earlier framing took R3's 25 to 35% as variation the model does not explain, and concluded the
monthly design fails. **The first half stands and the second does not follow.** In this model a
large part of a treated patient's marker variation is produced by the dose cycle and by the
recorded diet, and both are regressors. Variation a regressor explains is not error, and putting it
in sigma counts it twice.

### How much the model already explains **[measured]**

The marker's within-patient CV across the monthly record, produced by the dose train alone at a
perfectly constant recorded load, θ = 3000. **It is one quantity on three conventions for the
t = 0 draw, and it is recorded that way because a number without its convention is this
repository's own named failure** (`REVIEW_LOG.md` R54, S1.R 3a):

| t = 0 draw | n | levels | CV, sample sd | CV, population sd |
|---|---|---|---|---|
| read post-dose, `design.run.delta`'s own convention, so every draw is treated | 25 | 3.510 to 8.747 | 28.63% | **28.05%** |
| read pre-dose, the convention S1 adopted, so the t = 0 draw is untreated at m = 10.000 | 25 | 3.544 to 10.000 | 28.35% | **27.78%** |
| dropped, the on-therapy footing | 24 | 3.544 to 8.747 | 27.22% | **26.65%** |

With the specification's 10% day-to-day load variation added, over 24 records, the on-therapy
footing gives a median 30.0% [25.6, 33.5] and the 25-draw pre-dose footing 30.6% [27.0, 35.2],
sample sd.

The per-visit levels were checked by summing the kernel over doses and solving the balance by root
finding rather than by the closed-form root, agreeing to 1.3e-10 (`q10` check 1). **[measured]**

**Two reconciliations, both resolved, and neither a disagreement about the model.** The author's
independent figures are 28.0%, 27.8% and 26.6% on those three conventions, with the same level
ranges to three decimals. The 0.6 point gap against the first column is entirely the sample-versus-
population standard deviation: multiplying each sample-sd figure by `sqrt((n-1)/n)` gives 28.05,
27.78 and 26.65 against the author's 28.0, 27.8 and 26.6. **So the three conventions agree exactly
and the remaining choice is a definition, not a result.** The population sd is quoted in bold above
because it is the author's; the sample sd is the more conservative reading for a requirement and
both are given. An earlier draft of this section described the author's figures as mixing two
footings. That was wrong, it is withdrawn here, and the cause was that this file compared a
post-dose level range against pre-dose CVs.

**The point survives on every convention: roughly 27 to 30% of within-patient marker variation in a
treated patient on this design is dose cycle and recorded diet, which the model fits rather than
absorbs into sigma.** **[measured]**

### The decomposition, surrogate by surrogate

| | Thompson 2012, C3 | Manoli 2021, serum MMA |
|---|---|---|
| population | healthy adults | MMA, not PA |
| marker | propionylcarnitine | methylmalonic acid |
| timescale | daily, within-subject | between-visit, 2 months to 4 years |
| total within-patient CV | 26.4% (range 8.2 to 61.1) | 31.38 ± 24.55% (range 1.5 to 75.3) |
| (a) analytical | **7.9%, reported separately** | not reported in either report |
| biological remainder | 25.2% by quadrature | not separable |
| (b) what a model like this would explain | **not reported.** R3 attributes the variation to "feeding and activity"; feeding drives both C3 and the substrate load this model records, so an unknown and probably large part of the 25.2% is a regressor here | **not reported, and less separable.** No diet standardisation over an interval of years, and treatment changes, growth and intercurrent illness are inside the figure |
| (c) residual, the only part belonging in sigma | **cannot be separated** | **cannot be separated** |
| dose train | none. These subjects received no dose train, so none of this variation is the dose-cycle component above | none |
| treated at all | no | yes, dietary and pharmacological management, but no dose train |

1. **Neither surrogate can be decomposed into explained and unexplained parts.** Thompson separates
   analytical from biological and stops there. Manoli separates nothing. The quantity sigma needs,
   the scatter left after the dose cycle and the recorded load are fitted, is not reported by
   either and cannot be recovered from what they do report. **[judgement]**, and it is the central
   result of item 7.
2. **So the honest statement is a bound and not a value.** Sigma is bounded below by the analytical
   component and above by the total within-patient variation, with **no measured value in
   between**:
   - **lower bound 5 to 8%.** The design's 5% is analytical only. Thompson's analytical figure for
     a comparable assay is 7.9%, so 5% is the optimistic end of the analytical range rather than
     its centre. **[measured]**
   - **upper bound 25 to 35%**, the surrogates' totals, which are cross-disease or cross-population
     and contain components this model fits. **[measured]**
   This is weaker than either end taken alone, and it is what the evidence supports.
3. **The surrogates do not license subtracting 27% from 26.4%.** The dose-cycle CV is a property of
   the treated design, not a component of an untreated subject's variation, so the two numbers are
   not on the same axis and the arithmetic of differencing them is meaningless. What the comparison
   does establish is that the surrogates' totals are the same order as what the model explains, so
   carrying them into sigma unchanged would double-count on any reading. **[measured]** on the
   sizes, **[judgement]** on the conclusion.
4. **The design's 5% is not defensible as a central case, and neither is 25 to 35%.** There is no
   central case. That is a statement about the evidence, and the decision it forces on the author
   is whether the specification quotes a bound rather than a number. **[judgement]**

### Shared across markers, or marker-specific

The paper's panel remark rests on the residual being marker-specific, because only then can a panel
average it away. R3's evidence bears on this and does not settle it.

1. **Part of the residual must be shared, and that part is not averageable.** Variation in the
   propionyl-CoA pool itself propagates to every marker derived from it. R3's own functional forms
   make this explicit: 3-HP is a near-linear overflow of the pool, and C3 tracks it through a
   near-equilibrium exchange. **[judgement]**, from the reports' mechanism and not from data.
2. **Part must be marker-specific, and R3 names the confounders.** C3 is rescaled by free carnitine
   status, routinely low in PA and therapeutically supplemented. 2-MC needs oxaloacetate, so it
   carries the TCA state. Propionylglycine is made by GLYAT with sigmoidal kinetics and depends on
   glycine. These are different quantities, not scaled copies. **[measured]** as a quotation of
   R3, which marks them VERIFIED.
3. **The split cannot be estimated from anything R3 found.** Separating shared from specific
   requires repeated simultaneous measurements of two or more markers in the same patients at the
   same times, and R3 reports no within-patient variability figure for any marker in PA in any
   matrix at any timescale. **So the assumption the panel remark rests on is neither supported nor
   refuted, and the honest description of it is untested.** **[judgement]**
4. A consequence worth stating before any panel is costed: the fasting result R3 quotes, that
   fasting drops plasma C3/C2 by 62% driven mostly by C3, means nutritional state moves several
   markers together. **A shared driver of that size would defeat panel averaging**, and it is the
   kind of thing a standardised sampling protocol addresses more cheaply than a second assay.
   **[judgement]**

---

## Step 3(a). The requirement as a budget on sigma

`audit/q6_sigma_budget.py`. The two-axis surface S3 measured, and the tolerances reproduced above,
both dissolve into one number per schedule. Re-check: `q10` checks 3, by a numerical Hessian and a
direct scan with no root-finder. 12 of 12 overall.

### 3a.0 The budget is exact, not tight **[proven]**, **[measured]**

The design sees the marker's scatter and the load record's error **only** through the combined
ordinate noise, once the load-bias prior `s_b` is named. The reason is structural: `load_frac`
enters the information matrix twice, through `sigma = hypot(load_frac·L, b·CV·m_bar)` and through
`s_b = load_frac·L`, and the second entrance is not part of sigma. Holding `s_b` therefore leaves
`load_frac` nowhere else to act.

Five `(load_frac, marker CV)` pairs constructed to give exactly sigma = 60, `s_b` held at 0.10 L,
monthly, θ = 3000, from a marker-only budget (CV 20.14%, load 0) to an almost load-only one (CV
0.71%, load 18.58%):

| load record | marker CV | sigma | sd(ln θ̂) |
|---|---|---|---|
| 0.00% | 20.14% | 60.000000 | 0.2610866262 |
| 5.00% | 19.40% | 60.000000 | 0.2610866262 |
| 10.00% | 16.98% | 60.000000 | 0.2610866262 |
| 15.00% | 11.90% | 60.000000 | 0.2610866262 |
| 18.58% | 0.71% | 60.000000 | 0.2610866262 |

**Spread 3.3e-16, machine precision.** **[measured]**, and **[proven]** in the sense that the
Schur-complement argument makes it an identity rather than a coincidence: the noise enters only as
`QᵗQ/sigma²`.

**And the author's residual 4% is exactly the `s_b` channel, with nothing else in it.** With `s_b`
moving as `load_frac` does, sigma* spreads 2.57% (monthly) to 3.86% (every 2 days) across load
records of 2 to 15%. With `s_b` held at 0.10 L, the spread is **0.00e+00 at every schedule**.
**[measured]**

So the budget is exact once `s_b` is stated, and the only reason it looked tight is that a load
record does two jobs at once.

### 3a.1 The budgets **[measured]**

sigma*, the ordinate noise at which each schedule reaches `sd(ln θ) = 0.20`, `s_b` at 0.10 L:

| schedule | n | band low, 2472 | normal, 3000 | band high, 18814 |
|---|---|---|---|---|
| monthly | 25 | 20.80 | **36.72** | 21.75 |
| fortnightly, off-period | 51 | 26.75 | **45.99** | 27.97 |
| weekly, phase-spread | 100 | 40.25 | **69.69** | 43.14 |
| every 2 days | 365 | 74.68 | **130.96** | 79.57 |

The author's own figures at normal capacity are 35.5 to 37.0, 47.4 to 48.3, 68.7 to 71.3 and 129.3
to 133.9, against this session's `s_b`-moving ranges of 36.4 to 37.3, 45.3 to 46.9, 68.5 to 71.0 and
128.4 to 133.4. **Three of four overlap; fortnightly does not**, and the likely cause is the t = 0
convention, since this session reads dose-day draws pre-dose and the 14.5 d schedule has a
dose-instant draw at t = 392. The ordering and every conclusion below are unaffected. **[measured]**

**The budget roughly halves at both band edges**, and by almost the same factor at each, which is
worth noting because S3 found the band closing from the top: the budget does not distinguish the two
edges the way the band's width does. **[measured]**

### 3a.2 Where the committed footing and the evidence bound sit **[measured]**

At θ = 3000, converting each reading to sigma. The committed footing reproduces from
`design/run.py`'s own constants at **21.96**, against its stated 22.0.

| reading | load 5% | load 10% | load 15% |
|---|---|---|---|
| analytical floor, CV 5% | 22.0 | 35.5 | 50.6 |
| analytical floor, CV 8% | 28.8 | 40.1 | 54.0 |
| surrogate lower, CV 25% | 76.2 | 81.2 | 88.8 |
| surrogate upper, CV 35% | 105.5 | 109.2 | 115.0 |

**The whole finding in one line, at a 10% load record:**

| schedule | sigma* | committed 22 | analytical floor, 5% then 8% | surrogates, 25% then 35% |
|---|---|---|---|---|
| monthly | 36.72 | meets | **meets / fails** | fails / fails |
| fortnightly | 45.99 | meets | meets / meets | fails / fails |
| weekly | 69.69 | meets | meets / meets | fails / fails |
| every 2 days | 130.96 | meets | meets / meets | **meets / meets** |

1. **The committed sigma of 22 clears every schedule's budget, monthly included.** The design as
   booked is not in difficulty; the question is entirely whether 22 is the right sigma. **[measured]**
2. **Monthly is viable only at the optimistic end of the analytical floor.** At CV 5% the floor is
   35.5 against monthly's 36.72, a margin of 3%. At CV 8%, Thompson's own analytical figure for a
   comparable assay, the floor is 40.1 and monthly fails. **So monthly does not merely require the
   residual scatter to sit at its analytical floor, it requires the floor itself to be at the
   optimistic end of its range.** **[measured]**
3. **The surrogates' upper bound clears only the densest schedule**, as the author read it: 81.2 to
   109.2 against 130.96 for every-2-days and 69.69 for weekly. Weekly fails even at the surrogates'
   lower end. **[measured]**
4. **Nothing in the evidence distinguishes the three middle schedules**, because the bound is wide
   enough to contain all of them. That is the bound's fault, not the design's. **[judgement]**

### 3a.3 The load record alone can close the gate **[measured]**

**The band edges throughout this section are θ = 2472 and θ = 18814**, Section 6.2's own working band,
82% and 627% of the normal-capacity illustration, as `n5_scope_contribution.py` records them.

sigma has a floor the marker cannot get below: at marker CV zero it is still `load_frac·L`. Where
that floor exceeds a schedule's budget, **no assay however good reaches the gate**.

| capacity | load 2% (6.5) | load 5% (16.1) | load 10% (32.3) | load 15% (48.4) |
|---|---|---|---|---|
| normal, 3000 | all ok | all ok | all ok | **monthly and fortnightly OVER** |
| band low, 2472 | all ok | all ok | **monthly and fortnightly OVER** | weekly also OVER |
| band high, 18814 | all ok | all ok | **monthly and fortnightly OVER** | weekly also OVER |

1. **The author's 15% claim is confirmed at normal capacity**: the floor is 48.4 against monthly's
   36.72 and fortnightly's 45.99, so neither reaches 0.20 at any marker scatter including zero.
   **[measured]**
2. **At both band edges it happens at a 10% load record, not 15%.** The floor of 32.3 exceeds
   monthly's 20.80 and fortnightly's 26.75 there. **So at the edges of the working band the paper's
   own stated load requirement already closes the gate for both sparse schedules, whatever the
   assay does.** That is a stronger statement than the 15% one and it applies at the load accuracy
   the specification actually asks for. **[measured]**

### 3a.4 The Monte Carlo, on both nuisance conventions **[measured]**

An earlier draft of this section reported RMSE 0.33 to 0.42 at sigma* and attributed the miss to the
`k_e` prior's collapse tail as a disease-range mechanism. **Both halves needed correcting, and the
author caught both.** S2 established that the two honest ways of drawing the `k_e` nuisance answer
different questions and carried both; this section now does the same, and every figure is labelled
with its convention AND its seed count, because the second turned out to matter as much as the
first.

- **truth**: the patient's true `k_e` is drawn from the prior and the prior's centre stays at the
  reference. A population prior, with the patient a draw from the population.
- **centre**: the true `k_e` stays at the reference and the prior's **centre** is drawn. A
  per-patient measurement of `k_e` carrying that error, which is what the micro-dose study supplies.

At sigma*, where the Fisher bound reads exactly 0.200, `s_b` at 0.10 L, 96 records, `k_e` prior 0.8:

| schedule | convention | Fisher | RMSE | robust sd | collapses |
|---|---|---|---|---|---|
| monthly | truth | 0.200 | 0.220 | 0.200 | 0 |
| monthly | centre | 0.200 | 0.351 | 0.200 | 2 |
| fortnightly | truth | 0.200 | 0.224 | 0.223 | 0 |
| fortnightly | centre | 0.200 | 0.382 | 0.268 | 1 |
| weekly | truth | 0.200 | 0.241 | 0.229 | 0 |
| weekly | centre | 0.200 | 0.327 | 0.241 | 1 |
| every 2 days | truth | 0.200 | 0.214 | 0.201 | 0 |
| every 2 days | centre | 0.200 | 0.419 | 0.251 | 1 |

1. **On the population convention the budget is honoured to within 7 to 21%**, with zero collapses on
   every schedule. **[measured]**
2. **The disagreement with the author's independent figure was seed count, not convention.** The
   author's code fixes the truth at `k_e` and draws the prior centre, which is the **centre**
   convention, so this file's earlier relabelling of their 0.220 as the truth convention was wrong
   and is withdrawn. Their own code gives RMSE 0.220 at 80 draws, 0.307 at 96 with one collapse,
   4.125 at 160 with three, and 3.372 at 240, while the robust spread stays between 0.186 and 0.213
   throughout. The 0.220 was 80 draws that happened to miss the tail. Their 1 in 96 and 3 in 160
   agree with this session's 2 in 96, and their collapsing records are the ones handed a centre 3.7 to
   9.8 times the truth, which confirms the mechanism below by an independent route. **[measured]**
3. **RMSE is not a usable summary at the population prior, and that is the finding rather than a
   caveat.** It moves from 0.220 to 4.125 on seed count alone while the robust spread does not move at
   all. So the result is **two facts, not one**:
   - **at sigma* the budget is honoured in the bulk**, robust spread 0.200 against the Fisher 0.200 at
     monthly, and 0.24 to 0.27 on the denser schedules; and
   - **about 2% of records collapse by orders of magnitude, with no warning in the fit.** At this
     session's 2 in 96 that is 2.1%, 95% Clopper-Pearson [0.25%, 7.3%]; at the author's 3 in 160,
     1.9% [0.39%, 5.4%]. **The interval is wide and the rate is not pinned by any of these seed
     counts.** **[measured]**
4. **The earlier claim that the collapses are the disease-range mechanism was wrong, and so is the
   expectation that there are none at normal capacity.** They occur at normal capacity, and the
   records are identifiable: on the centre convention, monthly, the collapsing records are exactly
   those where the drawn `k_e` prior centre sits several times **above** the truth, 3.7 and 5.0 in this
   session's seed set. The prior then pulls `k_e·m` up until it absorbs load that θ should carry, and
   θ̂ falls, by 1.5 and 1.9 in ln units. At the band's lower edge, θ = 2472, the same two records
   collapse to the floor, 52 and 36 in ln units. At the band's **upper** edge, θ = 18814, and at
   θ = 300, there are none. **[measured]**
5. **So the mechanism is the opposite of a disease-range one.** It bites where sigma* is large
   relative to the enzymatic share of the balance, which is normal capacity and the lower band edge,
   and not where capacity is lowest: at θ = 300 sigma* is only 2.01, so the data pin θ whatever the
   prior does. **[judgement]** on the mechanism, **[measured]** on where it occurs.
6. **This is a second and independent reason to measure `k_e`, and a different kind of reason.**
   Measuring it does not merely narrow the spread, **it removes the tail**: at the realistic 0.054 the
   collapse count is zero on every schedule. A population prior leaves a per-patient failure mode
   whose cause is now named and which no goodness-of-fit check would catch. **[measured]**
7. **Quoting rule that follows, and it is S2's own.** Lead with the robust spread and the collapse
   count with its binomial interval. No RMSE in this phase is reportable without its nuisance
   convention AND its seed count, and a budget derived from a Fisher bound should be read against the
   robust spread, since the bound has no tail to match. **[judgement]**

### 3a.6 Two budget tables, not one conservative table **[measured]**

sigma* was solved at the 0.8 population prior throughout. **Solving it at S2b's realistic
per-patient 0.054 gives a different table, not a tighter version of the same one**, and the author
is right that the difference exceeds the spacing between adjacent schedules. `CV` is at a 10% load
record; `none` means the load record alone exceeds the budget.

| capacity | schedule | sigma* at k_e 0.8 | sigma* at k_e 0.054 | ratio | CV at 0.8 | CV at 0.054 |
|---|---|---|---|---|---|---|
| band low, 2472 | monthly | 20.80 | 57.73 | **2.78** | none | 16.24% |
| band low, 2472 | fortnightly | 26.75 | 71.23 | 2.66 | none | 21.54% |
| band low, 2472 | weekly | 40.25 | 108.07 | 2.69 | 8.16% | 35.00% |
| band low, 2472 | every 2 days | 74.68 | 203.87 | 2.73 | 22.85% | 68.30% |
| normal, 3000 | monthly | 36.72 | 56.10 | **1.53** | 5.88% | 15.40% |
| normal, 3000 | fortnightly | 45.99 | 69.83 | 1.52 | 11.00% | 20.79% |
| normal, 3000 | weekly | 69.69 | 105.95 | 1.52 | 20.73% | 33.87% |
| normal, 3000 | every 2 days | 130.96 | 199.36 | 1.52 | 42.60% | 66.03% |
| band high, 18814 | monthly | 21.75 | 21.88 | **1.006** | none | none |
| band high, 18814 | every 2 days | 79.57 | 80.05 | 1.006 | 22.91% | 23.07% |

1. **The author's 1.5 at normal capacity reproduces**: 36.72 to 56.10, against their 36.66 to 55.00.
   **[measured]**
2. **The ratio is not a constant, and that is the finding.** Measuring `k_e` is worth a factor of 2.7
   at the band's lower edge, 1.5 at normal capacity, and **nothing at all, 1.006, at the upper
   edge**, where the enzymatic term dominates the balance so tightly that `k_e`'s precision stops
   mattering. **[measured]**
3. **The cell the author names changes verdict as predicted.** At a 10% load record with marker
   scatter at 15%, sigma is about 55: that is above monthly's budget at the population prior, 36.72,
   and just below it at a measured `k_e`, 56.10. **So monthly's verdict in that cell is decided
   entirely by which prior the budget was solved at, and a table quoted without naming the prior is
   not interpretable.** **[measured]**
4. Since a measured `k_e` is what S2b's study is for, the second column is the one a specification
   written around that study should quote, and the first is the one to quote if the study is not
   done. **They are not upper and lower bounds on one number.** **[judgement]**

### 3a.7 The fortnightly gap, settled **[measured]**

This session's fortnightly sigma* is 45.99 and the author's is 47.4 to 48.3. The first hypothesis was
the dose-instant draw at t = 392, where `0.5 + 14.5k` is a multiple of 14. **Tested, and it is not
the cause:**

| variant | sigma* |
|---|---|
| as is, pre-dose | 45.99 |
| t = 392 moved to 392.25 | 46.47 |
| t = 392 dropped, n = 50 | 45.02 |
| as is, read post-dose | 46.58 |

None reaches 47.4, so the convention accounts for at most 0.6 of a gap of 1.4 to 2.3. **The
schedules are not the same object**, and the difference is the first draw:

| schedule | n | sigma* |
|---|---|---|
| `arange(0.5, 730, 14.5)`, this session | 51 | 45.99 |
| `arange(0.0, 730, 14.5)`, t = 0 included | 51 | **47.56** |
| `arange(7.0, 730, 14.5)`, trough phase | 50 | 50.35 |
| `linspace(0.5, 730, 51)` | 51 | 51.13 |

**`arange(0.0, 730, 14.5)` lands inside the author's range.** So the author's fortnightly record
starts at t = 0 and, on the pre-dose reading, carries an untreated baseline draw; this session's
starts half a day in and carries none. That is the baseline-draw effect S1.R identified and S4
priced, worth 1.6 in sigma* here, and not a dose-instant artifact. **Both are legitimate records and
the table must say which.** The budgets quoted in 3a.1 are the no-baseline ones. **[measured]**

**One warning found on the way, and not part of this question.** A 14.0 d draw schedule against the
14 d dose train puts every draw at one cycle phase, the increment is then constant across the record,
and sigma* falls to 0.64, a factor of 72. That is `OPEN_QUESTIONS` item (c), 14 d and 14.5 d
identical on a calendar axis and twentyfold apart in the answer, appearing in the budget.
**[measured]**

### 3a.8 The practical form, for a protocol writer

A protocol writer spends on diet records and assays, not on sigma, so the budget's practical form is
the marker CV each budget allows at each load record. At θ = 3000, `s_b` at 0.10 L, **sigma* solved
at the 0.8 population prior**, no baseline draw; 3a.6 gives the measured-`k_e` column.

| schedule | sigma* | load 2% | load 5% | load 10% | load 15% |
|---|---|---|---|---|---|
| monthly | 36.72 | 12.13% | 11.07% | 5.88% | none |
| fortnightly | 45.99 | 15.28% | 14.45% | 11.00% | none |
| weekly | 69.69 | 23.29% | 22.75% | 20.73% | 16.83% |
| every 2 days | 130.96 | 43.90% | 43.62% | 42.60% | 40.84% |

The pairing the author noticed is now arithmetic rather than coincidence: monthly's 11.07% at a 5%
load record and fortnightly's 11.00% at a 10% one are the same budget spent differently, and
monthly's 5.88% at a 10% load record is the same budget again. **[measured]**


---

## Step 3(b). Identifiability with the marker's unreachable offset fixed

`audit/q7_offset_fixed.py`, 6 of 6 internal structural checks passing. Phase 5 proved that carrying
the offset `E0` free is what makes θ unidentifiable under the published algebraic structure, and that
fixing it restores identifiability. What Phase 5 could not do was say what to fix it at. R5 supplies
a human anchor from liver transplantation, 5 to 10% of the pre-treatment level with a range to about
30%, against a published value that is a mouse estimate. **Nothing in R5 is verified here**: its two
transplant ratios are step 2 item 2, which has not run.

**The offset's definition, stated because the structure admits three.** The observable is
`m(E) = [E0 + base·f(E)]·(L+b)/L0` with `f(E) = (IC50 + qE)/(IC50 + E)` and `q = 1 − Imax`. A
transplant replaces the whole hepatic pool, which is `E → ∞`, so what a transplant ratio measures is
`(E0 + base·q) / m(untreated)`, **not `E0` alone**: at `Imax = 0.999` the incomplete-inhibition term
`base·q` rides along inside it. This section fixes that ratio, called φ, and reports the `E0` it
implies. The structure reproduces `k2_structures.PubB` to 8e-13 at `PubB`'s own ratio, which is
11.277% on this definition.

### 3b.0 A constraint that had to be found first: φ is bounded both ways **[measured]**

**The offset and the paper's calibrated fall are not independent.** Therapy can only remove what is
reachable, so a large unreachable fraction caps the fall the structure can produce. And `E0` cannot
be negative, which puts a floor under IC50 at `IC50_min = qθ(1−φ)/(φ−q)`; since the fall falls
monotonically in IC50 above that floor, the ceiling on the fall is attained exactly at `IC50_min`.

| φ | IC50 at least | admissible ceiling on the fall | 60% reachable |
|---|---|---|---|
| 0.26% | 1870.1 | 56.64% | **no** |
| 1.00% | 330.0 | 65.12% | yes |
| 2.00% | 154.7 | 65.64% | yes |
| 5.00% | 58.2 | 64.28% | yes |
| 10.00% | 27.3 | 61.10% | yes |
| 12.77% | 20.7 | 59.26% | **no** |
| 20.00% | 12.1 | 54.40% | **no** |
| 30.00% | 7.0 | 47.62% | **no** |

1. **The φ window compatible with the paper's own 60% fall is 0.343% to 11.660%.** **[measured]**
2. **R5's 5 to 10% anchor fits. R5's 20% sensitivity value does not**, and neither does 30%: at 20%
   the structure's best achievable fall is 54.4%, so it cannot be calibrated to the paper's 60% at
   any IC50. **The author's sweep therefore has two admissible points and one inadmissible one, and
   that is a property of the published structure rather than of the evidence.** **[measured]**
3. **The constraint is a window and not a ceiling, and its lower edge is set by `Imax`, not by
   biology.** As φ falls toward `q = 1 − Imax` the unreachable part stops being `E0` and becomes the
   incomplete-inhibition term, `IC50_min` diverges, and a large IC50 caps the fall just as a large φ
   does. **[measured]**, **[judgement]** on the reading.
4. **The published structure sits 0.383 points inside the upper edge**, at 11.277% against 11.660%.
   That is a second knife edge of the same kind `REVIEW_LOG.md` R15 recorded for `Imax = 0.999`, and
   it is worth stating because it means the published calibration is close to being unable to produce
   its own headline fall. **[measured]**

**And the window redone at `Imax = 1`, as R15's own correction requires.** The edges rest on
`Imax = 0.999`, a constant the source fixed rather than estimated, and R15 retracted a verdict over
exactly that. At `Imax = 1` the residual `q` vanishes, so the unreachable ratio is `E0/m(untreated)`
with no incomplete-inhibition term inside it, and `IC50_min` falls to zero because `E0` can no longer
be driven negative:

| | φ window for a 60% fall | φ = 5% | φ = 10% | φ = 20% |
|---|---|---|---|---|
| `Imax = 0.999` | 0.343% to **11.660%** | ceiling 64.28%, calibrates | 61.10%, calibrates | 54.40%, **inadmissible** |
| `Imax = 1` | 0% to **11.875%** | 64.68%, calibrates | 61.28%, calibrates | 54.47%, **inadmissible** |

**The substantive conclusion survives R15's correction**: 5% and 10% calibrate on both readings and
20% does not, its best achievable fall being 54.4 to 54.5%. What moves is the lower edge, which
**exists only because `q > 0`** and is therefore an artefact of the fixed `Imax` rather than a
constraint on biology. A coincidence worth flagging so nobody reconciles the two wrongly: the figure
11.875% appeared in conversation as a wrong `Imax = 0.999` edge, from the admissibility bug below, and
is separately the correct `Imax = 1` edge. **[measured]**

**An error of mine, corrected before any of this was reported.** The first version of the ceiling
maximised the fall over all IC50 and returned 64.7% at φ = 5% and 54.5% at φ = 20%, attained at
`IC50 = 0.001` where `E0` is about −28,500. A negative unreachable term is not a structure. The
admissible ceilings above are lower, the conclusion about 20% survives, and the non-monotone lower
edge only appeared once admissibility was enforced.

### 3b.1 Fixing the offset restores STRUCTURAL identifiability **[proven]**, **[measured]**

What this section establishes is a rank condition and an exact fibre argument. **It does not establish
that θ can be recovered from a record**, which is 3b.2's question and comes out differently.

A constant-load record fixes only `P = E0 + base·q`, `U = P·θ + (E0 + base)·IC50` and
`V = IC50 + θ`. Symbolically, `rank ∂(P,U,V)/∂(E0, base, IC50, θ) = 3` against 4 unknowns, a
one-parameter fibre; with `E0` held, `rank ∂(P,U,V)/∂(base, IC50, θ) = 3` against 3, identifiable.
**[proven]**

Walking that fibre at φ = 10%, with θ re-solved from the invariants at each `IC50 = R`:

| IC50 = R | θ = V − R | base | E0 | max relative change in the record, E0 free | same, E0 held |
|---|---|---|---|---|---|
| 160 | 3040 | 180.2 | 0.820 | 2.0e-16 | 4.1e-03 |
| 640 | 2560 | 45.05 | 0.955 | 2.0e-16 | 1.1e-02 |
| 1600 | 1600 | 18.02 | 0.982 | 2.0e-16 | 1.4e-02 |
| 3040 | 160 | 9.48 | 0.991 | 2.0e-16 | 1.5e-02 |

**With the offset free, θ slides from 160 to 3040, a factor of 19, at a bit-identical record. With it
held, the same walk moves the record by a factor of 3.8e13 more.** **[measured]**

### 3b.2 What the anchor buys: locally identifiable is not identified **[measured]**

**This section replaces a second wrong one, and the correction removes its headline.** The earlier
draft read the information matrix and reported that fixing the offset at φ = 10% takes `sd(ln θ)` from
infinite to 0.135, inside the gate, calling it the strongest consequence R5 has for this repository.
The author asked for the factor of 36 between the two admissible anchors to be confirmed by the
profile route that had already caught one artefact here. **It does not survive, and neither does the
0.135.**

The information matrix, as before:

| φ | marker CV | E0 free | E0 fixed, Fisher |
|---|---|---|---|
| 5% | 5% | infinite | 4.817 |
| 10% | 5% | infinite | 0.1347 |

And the same cells fitted, records simulated with the load bias drawn from its prior and fitted with
it, the offset held at its own correct value, 80 records:

| φ | marker CV | Fisher | fitted robust spread | gate met |
|---|---|---|---|---|
| 5% | 5% | 4.817 | 6.34 | no |
| 5% | 1.5% | 1.461 | 5.75 | no |
| 5% | 0.5% | 0.532 | 4.94 | no |
| 5% | 0.05% | 0.233 | 3.44 | no |
| 10% | 5% | **0.135** | **6.12** | **no** |
| 10% | 1.5% | 0.055 | 5.40 | no |
| 10% | 0.5% | 0.042 | 3.15 | no |
| 10% | 0.1% | 0.040 | 0.121 | **yes** |

1. **The Fisher bound at φ = 10% is not achievable at any marker CV the literature supports.** At the
   committed 5% the estimator's robust spread is 6.1 in ln θ, a factor of 45 above the bound, so θ is
   unknown to within orders of magnitude. **[measured]**
2. **It is not an optimiser failure.** Starting the fit exactly at the truth gives the same answer,
   robust spread 5.47 at φ = 10% and a 5% marker CV. The likelihood genuinely has no useful curvature
   there: the optimiser starts on the truth and walks away because a better fit exists far from it.
   **[measured]**
3. **The mechanism is a needle.** The noiseless profile with the correct offset is 1.6e-16 at the
   truth and only 1.3e-03 at a quarter of it. So fixing the offset produces a minimum that is sharp
   in an infinitesimal neighbourhood and sits on a nearly flat floor, and the information matrix
   measures the needle while a real record sees the floor. A 5% marker CV produces noise of 5e-02,
   forty times the floor's depth. **[measured]**
4. **So the gate is met only at a marker CV near 0.1%**, which is **two orders of magnitude below the
   5 to 8% analytical floor** item 7 established. There is no assay in prospect at that level.
   **[measured]**
5. **And the factor of 36 between the anchors is not real at realistic noise: both fail.** It appears
   only below about a 0.5% marker CV, where φ = 10% works and φ = 5% still does not. **[measured]**
6. **Phase 5's "fixing E0 restores identifiability" is true in the rank sense and misleading in the
   practical one.** Rank 3 of 3 is what 3b.1 proves and it is correct. **What it licenses is a claim
   about structure, not a claim that θ can be recovered**, and this file's earlier draft made the
   second claim from the first. **[judgement]**

### 3b.2b The pattern this makes, and it is the second instance **[judgement]**

The author's observation on the recommendation is decisive and it generalises. At φ = 5%, the safe end
of R5's range, the anchor is safe against the bias of 3b.3 and useless for precision. At φ = 10% it
was thought usable and is unsafe, and 3b.2 now shows it is not usable either. **There is no value in
R5's stated span at which the anchor is both safe and useful.**

**That is exactly the shape S5 found for the ex vivo prior on θ**: a wide prior is safe and useless, a
tight one useful and unsafe, with no width at which it is both. Two independent external anchors, on
two different parameters, in two different structures, have produced the same trade. **It is worth
naming as a pattern rather than recorded twice as a coincidence: an anchor borrowed through an
unverified mapping buys precision only in proportion to how much of its own error it injects, so the
mapping's verification is not a refinement of the anchor, it is the anchor.** **[judgement]**

**And it changes what step 2 item 2 is for.** Verifying R5's transplant ratios stops being a check
that a range is roughly right and becomes **a requirement to pin a value**, because a 5-to-10 span
spans useless to unsafe. The specification should say that: the anchor is admissible only if the
mapping is pinned, not merely bounded. **[judgement]**

### 3b.3 Getting the anchor wrong does not bias θ, it un-identifies it **[measured]**

**This section replaces a wrong one.** An earlier draft fitted a point estimate with the offset held
at the wrong value and reported a bias, θ̂/θ = 0.141, a sevenfold underestimate. That number was an
artifact of the fit's own θ grid: the objective has no interior minimum there, so the fit returned a
value near its lower bound and the "bias" was the bound. An independent `least_squares` from many
starts drove θ̂ to zero at the same residual, which is how the artifact was caught. The profile is
what answers the question, and it says something worse than a bias.

The best achievable relative residual at each θ with the offset **held**, `(base, IC50)` optimised
exactly at every node. Truth is φ = 5%, `E0 = 0.45159`, θ = 3000:

| offset held at | θ=30 | θ=300 | θ=1000 | θ=2000 | **θ=3000, the truth** | θ=6000 |
|---|---|---|---|---|---|---|
| its own 5%, **correct** | 3.20e-04 | 3.14e-04 | 2.92e-04 | 2.32e-04 | **1.6e-16** | 6.9e-02 |
| 10%, **too high** | 3.21e-03 | 3.22e-03 | 3.24e-03 | 3.33e-03 | **4.09e-03** | 7.9e-02 |
| 2%, **too low** | 2.59e-03 | 2.59e-03 | 2.57e-03 | 2.52e-03 | **2.38e-03** | 6.1e-02 |

1. **With the correct anchor the profile has a sharp minimum exactly at the truth**, at machine
   precision, and within 2% of that minimum only at the truth itself. That is 3b.1's identifiability
   made numerical. **[measured]**
2. **With the anchor held too high, θ is not identified at all.** The minimum sits at θ = 30, the
   profile is flat within 2% from θ = 30 to θ = 1000, a factor of 33, and **the true value fits worse
   than every point below it**. So a wrong anchor in the unsafe direction gives back exactly the flat
   profile Phase 5 found with the offset free, and points away from the truth while doing it.
   **[measured]**
3. **With the anchor held too low the minimum stays at the truth.** So the asymmetry is real and it is
   sharper than a bias: **too low costs precision, too high costs identifiability.** **[measured]**
4. **And none of it is visible in the fit.** The too-high row's residual is 3.2e-03 against the
   correct row's 1.6e-16, but a 5% marker CV produces noise of 5e-02, more than ten times larger, so
   no goodness-of-fit test on a real record would separate them. **[measured]**
5. **So a high anchor is not an anchor at all**, and the reason is not bias but the loss of the very
   identifiability the anchor was meant to supply. This was drafted as a recommendation to fix at or
   below the bottom of R5's range; 3b.2 removes that recommendation, because the bottom of the range
   does not deliver a usable θ either. **[judgement]**, on **[measured]** sizes.
6. **That sharpens what step 2 item 2 is for.** It is not asking whether R5's transplant ratios are
   roughly right. **It is asking for a defensible LOWER bound on the unreachable fraction**, because
   any value above the truth returns the problem the anchor was meant to solve. **[judgement]**

### 3b.4 The factor of 5.8 does not move, and why that is the answer **[judgement]**

The author asked what fixing the offset does to the factor of 5.8 S2b priced for the observation-model
assumption. **It does not move it, and the reason is structural rather than numerical.**

S2b's 5.8 is the ratio of what a per-patient `k_e` prior buys (7.6 on `sd(ln θ̂)`) to what a joint
`k_e + θ/(K_m+m₀)` observation buys (1.3), **in the PAPER structure**, the quasi-steady balance
`(θ+Δ)m/(K_m+m) + k_e·m = L + b`. The offset `E0` is a parameter of **PUB_B**, the published
algebraic structure, which has no `k_e` term at all: its marker is an algebraic function of the
enzyme with no alternative clearance route to be uncertain about. So there is no cell of one
calculation in which both quantities appear, and no trade-off between them to compute. **[judgement]**

What the two do have in common is that both are observation-model questions, and they act through
different channels:

| | what is uncertain | channel | size |
|---|---|---|---|
| S2b's 5.8, PAPER structure | whether the marker is itself cleared enzymatically | **variance** on `sd(ln θ̂)` | factor 5.8 |
| 3b's anchor, PUB_B structure | how much of the marker is unreachable | **bias** on θ̂, undetectable | factor 7 if set too high |

**So a specification has to choose a structure before either number means anything, and the honest
statement is that the two largest observation-model risks the audit has found are not commensurable
and are not additive.** Computing a cross-structure analogue, a PAPER-structure balance carrying an
unreachable floor, would make them comparable; that is a new model and it is filed rather than run.
**[judgement]**


---

## The pricing: what each measurement buys, and what the repeatability study costs

`audit/q8_pricing.py`. Three of the four quantities a specification would have to ask for move the
same sigma budget, and step 3(a) proved the design sees them only through sigma once `s_b` is named,
so they are priced as one arithmetic constraint rather than three sweeps. The fourth, the marker's
unreachable floor, acts on identifiability and is priced apart.

### P.1 The three sigma axes, jointly **[measured]**

`k_e` picks the budget; the load record and the marker CV spend it. Each cell is the largest marker
residual scatter admissible, θ = 3000, `s_b` at 0.10 L:

| k_e | schedule | sigma* | load 2% | load 5% | load 10% | load 15% |
|---|---|---|---|---|---|---|
| population, 0.8 | monthly | 36.72 | 12.13% | 11.07% | 5.88% | none |
| population, 0.8 | fortnightly | 45.99 | 15.28% | 14.45% | 11.00% | none |
| population, 0.8 | weekly | 69.69 | 23.29% | 22.75% | 20.73% | 16.83% |
| population, 0.8 | every 2 days | 130.96 | 43.90% | 43.62% | 42.60% | 40.84% |
| measured, 0.054 | monthly | 56.10 | 18.70% | 18.03% | 15.40% | 9.52% |
| measured, 0.054 | fortnightly | 69.83 | 23.34% | 22.80% | 20.79% | 16.89% |
| measured, 0.054 | weekly | 105.95 | 35.50% | 35.15% | 33.87% | 31.63% |
| measured, 0.054 | every 2 days | 199.36 | 66.88% | 66.70% | 66.03% | 64.91% |

1. **The whole table is one identity**, `sigma = hypot(load_frac·L, b·CV·m̄)` inverted, so it needs no
   simulation and can be checked by hand. **[measured]**
2. **The drift term stops mattering at normal capacity.** A study at 2% with 5% drift gives 56.10 and
   one with no drift 56.16, a difference of 0.1%, which is S4.3's negative result arriving in the
   pricing. So **the specification does not need to ask for a repeated `k_e` measurement at normal
   capacity**; S2b.5's requirement for one is disease-range only. **[measured]**

### P.2 What the `k_e` measurement buys, both things **[measured]**

| schedule | sigma* at 0.8 | at 0.054 | budget gain | marker CV at load 10% | collapses at 0.8 | at 0.054 |
|---|---|---|---|---|---|---|
| monthly | 36.72 | 56.10 | 1.53x | 5.9% to 15.4% | 2 in 96, 2.1% [0.25, 7.3] | **0** [0, 3.8] |
| fortnightly | 45.99 | 69.83 | 1.52x | 11.0% to 20.8% | 1 in 96, 1.0% [0.03, 5.7] | **0** [0, 3.8] |
| weekly | 69.69 | 105.95 | 1.52x | 20.7% to 33.9% | 1 in 96 | **0** |
| every 2 days | 130.96 | 199.36 | 1.52x | 42.6% to 66.0% | 1 in 96 | **0** |

1. **Measuring `k_e` buys a factor of 1.5 on the budget, which nearly triples the admissible marker
   CV at a 10% load record**, from 5.9% to 15.4% at monthly. **[measured]**
2. **And it removes the tail, which is the part a variance figure does not show.** Zero collapses on
   every schedule at 0.054 against 1 or 2 in 96 at the population prior. The robust spread also
   halves, 0.200 to 0.104 at monthly. **[measured]**
3. **The intervals are wide and honest**: zero in 96 bounds the rate at 3.8%, not at zero. **A study
   that wanted to demonstrate the tail is gone would need several hundred records.** **[measured]**

### P.3 The repeatability study, sized **[measured]**

The design: `n` treated patients, `r` repeat marker draws each, **every draw at the same dose-cycle
phase**, with the diet recorded. Holding the phase removes the dose cycle, which item 7 measured at 27
to 30% and which is a regressor rather than error; recording the diet removes the load term. What is
left is the residual plus the analytical component, which is what sigma needs and what no published
study reports. The pooled within-patient variance carries `df = n(r−1)`.

| n | r | df | sd(ln σ̂) | 95% factor on σ |
|---|---|---|---|---|
| 10 | 2 | 10 | 0.235 | ×0.70 to ×1.75 |
| 10 | 3 | 20 | 0.162 | ×0.77 to ×1.44 |
| 25 | 3 | 50 | 0.101 | ×0.84 to ×1.24 |
| 10 | 6 | 50 | 0.101 | ×0.84 to ×1.24 |
| 50 | 3 | 100 | 0.071 | ×0.88 to ×1.16 |

Confirmed by simulating log-normal residuals: at df 20 and 50 the empirical `sd(ln σ̂)` matches the
trigamma value to 2% and the nominal 95% interval covers at 94.5 to 95.6%. **[measured]**

**What precision is enough is set by the decision.** The schedule boundaries at a 10% load record and
the population prior are monthly 5.88%, fortnightly 11.00%, weekly 20.73%, every 2 days 42.60%. A
study resolves a boundary when its interval excludes it:

| if the true CV is | nearest boundary | relative gap | df needed | n at r=3 | n at r=6 |
|---|---|---|---|---|---|
| 6% | monthly, 5.88% | 2.0% | 4429 | 2215 | 886 |
| 8% | monthly, 5.88% | 26.5% | 15 | **8** | **3** |
| 12% | fortnightly, 11.00% | 8.4% | 230 | 115 | 46 |
| 15% | fortnightly, 11.00% | 26.7% | 15 | **8** | **3** |
| 25% | weekly, 20.73% | 17.1% | 45 | 23 | 9 |
| 30% | weekly, 20.73% | 30.9% | 10 | **5** | **2** |

1. **The study is cheap unless the truth sits near a boundary.** Ten patients with three draws each,
   thirty blood draws in total, resolves the CV to ×0.77 to ×1.44 and settles every case above that
   is more than about 25% clear of a boundary. **This is by far the cheapest item in the
   specification.** **[measured]**
2. **It is unboundedly expensive when the truth sits on a boundary**, 2215 patients at a true CV of
   6% against monthly's 5.88%. That is not a reason to skip it: **the same study tells you which case
   you are in**, since a df-20 interval that straddles a boundary is itself the answer that the
   schedule cannot be chosen on this evidence. **[judgement]**
3. **Run it at df 20 first and decide then.** Thirty draws buys either a decision or the knowledge
   that a decision needs more, and no other item in the specification has that property.
   **[judgement]**

### P.4 The floor anchor, priced on its own terms **[measured]**, **[judgement]**

Step 3(b) removed the question this was going to answer. Fixing the floor restores structural
identifiability and does not deliver a usable θ, so what is priced is **the marker CV at which the
anchor would begin to be worth having**:

| requirement | value | source |
|---|---|---|
| marker CV at which the floor anchor pays | about **0.1%** | 3b.2, fitted |
| analytical floor on the marker CV | 5 to 8% | item 7, Thompson 2012 |
| ratio | **50 to 80x** | arithmetic |
| marker CV the surrogates suggest | 25 to 35% | item 7 |
| ratio | **250 to 350x** | arithmetic |

**The price is not payable.** The floor anchor needs an assay two orders of magnitude better than the
best the literature supports, so verifying R5's transplant ratios cannot on its own make θ
recoverable under the published algebraic structure. What the anchor is still worth is the structural
statement the paper's identifiability claim rests on, and that is not nothing: **it is the difference
between a claim that is false and one that is true but not yet operational.** **[judgement]**

### P.5 What the pricing says, in order of cost

1. **The repeatability study, about 30 blood draws.** It measures the one quantity sigma needs, it is
   the cheapest item by a wide margin, and it returns either a schedule decision or the knowledge
   that none is available. **Run it first.** **[judgement]**
2. **The `k_e` micro-dose study**, one bolus with a 24 h window and 8 or more draws per patient
   (S2b.1). It buys a factor of 1.5 on the budget and removes a 2% per-patient failure mode. At
   normal capacity it need not be repeated. **[measured]**
3. **The load record to better than 5%**, which inside the booked two years is the binding
   requirement S3 found, and which at 15% closes the gate for the two sparse schedules whatever the
   assay does (3a.3). **[measured]**
4. **The label-fate measurement** inside the `k_e` study, which settles the observation-model
   assumption S2b priced at a factor of 5.8 and costs no extra study (S2b.4). **[judgement]**
5. **The transplant-ratio verification**, which is worth doing for the structural claim and will not
   on its own recover θ. **[judgement]**


---

## What Phase 9 has and has not done

**Done.** The reports are recorded with their mapping and process warnings. Item 6 is adjudicated
and the repository's standing retraction is confirmed, conditionally on two quotations. Item 7's
decomposition is attempted and its central result is that the decomposition is not possible from the
published record, which converts the specification's sigma from a number into a bound.

**Not done, and not to be inferred.** Items 1 to 5 of step 2 are not started, for want of egress:
the nine approved hosts were re-tested on 26 September and all nine still return a 403 on CONNECT.
No constant has moved. Step 3(b) is done, with the qualification that its 20% case turns out to be
inadmissible under the published structure rather than merely unlikely.

**Reproduced, and it raises one question before step 3.** The author's four schedule tolerances at
normal capacity over two years, the largest total residual scatter at which `sd(ln theta)` reaches
0.20, were reproduced here to locate their footing rather than taken on trust:

| schedule | n | author | this footing, load recorded to 5% | this footing, load recorded to 10% |
|---|---|---|---|---|
| monthly | 25 | 11.0% | 10.9% | **5.9%** |
| fortnightly, off-period | 51 | 15.2% | 14.7% | 11.0% |
| weekly, phase-spread | 100 | 23.2% | 23.1% | 20.7% |
| every 2 days | 365 | 44.4% | 44.2% | 42.6% |

The author's figures reproduce on a footing that reads the t = 0 draw post-dose, carries both the
load-recording and marker terms in sigma, and holds **the load record to 5%**. **At the paper's own
stated load requirement of about 10%, monthly's tolerance halves, from 11.0% to 5.9%.**
**[measured]**

That matters for the shape of the requirement table, and it sharpens the author's reading rather
than softening it. At a 10% load record, monthly demands a residual scatter of 5.9%, which is
**below the 5 to 8% analytical floor established above**, so monthly does not merely require the
scatter to sit near its floor, it requires the floor itself to be at the optimistic end. And the
coincidence the author noticed, monthly's tolerance landing on the load requirement's own figure,
holds one step along: it is 11.0% at a 5% load record and 11.0% is the fortnightly tolerance at a
10% one. **The two axes trade against each other exactly as S3 measured, so the requirement table
has to fix one axis to state the other, and which value it fixes changes the monthly verdict from
"demanding" to "impossible".** That is a decision for the author and step 3 waits on it.
**[judgement]**

**Filed rather than pursued.**

- Whether a standardised sampling protocol, fasting state and time of day fixed, would narrow item
  7's upper bound enough to matter. It bears on both the sigma bound and the panel question.
- **A PAPER-structure balance carrying an unreachable floor**, which would make the offset anchor and
  the `k_e` observation-model assumption commensurable (3b.4). It is a new model.
- **Re-solving sigma\* at the measured `k_e` prior throughout**, rather than carrying two tables. P.1
  shows the two differ by a factor of 1.5, so the choice of which to quote is a specification
  decision and not an audit one.
- **The between-patient spread of the marker's residual scatter**, which P.3's pooled design does not
  estimate and which decides whether one requirement can be written for every patient.
