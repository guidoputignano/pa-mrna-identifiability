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

## What Phase 9 has and has not done

**Done.** The reports are recorded with their mapping and process warnings. Item 6 is adjudicated
and the repository's standing retraction is confirmed, conditionally on two quotations. Item 7's
decomposition is attempted and its central result is that the decomposition is not possible from the
published record, which converts the specification's sigma from a number into a bound.

**Not done, and not to be inferred.** Items 1 to 5 of step 2 are not started, for want of egress.
No constant has moved. The two reruns are not begun: item 7 was to supply the central case for the
marker-variability rerun and instead establishes that there is no central case, so what the rerun
should hold sigma at is now a decision for the author and not a reading from the literature.

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
