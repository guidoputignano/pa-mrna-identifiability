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

### 3a.4 The Monte Carlo, and whose miss it is **[measured]**

A budget the estimator does not honour is not a budget, so each schedule was re-run through the MAP
estimator at its own sigma*, 96 records, `s_b` at 0.10 L. **At the paper's own 0.8 population prior
on `ln k_e` the estimator misses badly**: RMSE 0.33 to 0.42 where the Fisher bound reads exactly
0.200, with one or two collapses per cell.

**That miss is the `k_e` prior's, not the budget's.** Re-running the identical cells with the `k_e`
prior at 0.054, the realistic per-patient precision S2b says the micro-dose study delivers, and
changing nothing else, sigma* left where the 0.8 prior put it:

| schedule | sigma* | Fisher | RMSE at k_e 0.8 | robust | RMSE at k_e 0.054 | robust | collapses |
|---|---|---|---|---|---|---|---|
| monthly | 36.72 | 0.200 | 0.351 | 0.200 | **0.110** | 0.104 | 0 |
| fortnightly | 45.99 | 0.200 | 0.382 | 0.268 | **0.150** | 0.144 | 0 |
| weekly | 69.69 | 0.200 | 0.327 | 0.241 | **0.142** | 0.159 | 0 |
| every 2 days | 130.96 | 0.200 | 0.419 | 0.251 | **0.132** | 0.125 | 0 |

1. **Once `k_e` is measured per patient, the budget is honoured and is conservative.** RMSE 0.110 to
   0.150 against a gate of 0.20, no collapses anywhere. **[measured]**
2. **At the population prior the estimator misses sigma*'s gate by a factor of about 1.8, and the
   robust spread shows why:** it is 0.200 at monthly, exactly the Fisher figure, while the RMSE is
   0.351. The difference is a collapse tail, which is the wide `k_e` prior's signature that S2
   identified, not a property of the noise budget. **[measured]**
3. **So the budget is conditional on the same thing everything else in Phase 8 was conditional on.**
   sigma* is a valid requirement to state provided `k_e` is measured; quoted against a population
   prior it is optimistic by about 1.8 and the specification should not quote it that way.
   **[judgement]**
4. A direction not pursued, and filed: sigma* was solved at the 0.8 prior throughout, so **the
   budget for a measured-`k_e` patient is larger than the numbers above**, and every verdict in
   3a.2 is therefore conservative. Re-solving sigma* at 0.054 would raise all four budgets and could
   move monthly's verdict against the analytical floor. That is a new computation and is filed, not
   run. **[judgement]**

### 3a.5 The practical form, for a protocol writer

A protocol writer spends on diet records and assays, not on sigma, so the budget's practical form is
the marker CV each budget allows at each load record. At θ = 3000, `s_b` at 0.10 L:

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

## What Phase 9 has and has not done

**Done.** The reports are recorded with their mapping and process warnings. Item 6 is adjudicated
and the repository's standing retraction is confirmed, conditionally on two quotations. Item 7's
decomposition is attempted and its central result is that the decomposition is not possible from the
published record, which converts the specification's sigma from a number into a bound.

**Not done, and not to be inferred.** Items 1 to 5 of step 2 are not started, for want of egress:
the nine approved hosts were re-tested on 26 September and all nine still return a 403 on CONNECT.
No constant has moved. Step 3(b), identifiability with the marker offset fixed at 5, 10 and 20% of
baseline, is not begun.

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

**Filed rather than pursued.** Whether a standardised sampling protocol, fasting state and time of
day fixed, would narrow the upper bound enough to matter. It bears on both the sigma bound and the
panel question, and it is a new sub-question.
