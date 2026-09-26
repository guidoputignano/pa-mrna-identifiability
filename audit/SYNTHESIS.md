# Synthesis, by claim

**What this is.** One entry per claim the paper makes or the audit added, in the order a reader meets
them, with everything needed to decide whether the claim can be written today. It replaces reading
eight phase files in sequence. It does not reproduce them: each entry points at the file that holds
the working.

**Marks**, as `PHASE4_SUMMARY.md` uses them. **[proven]** mathematics, exact or symbolic.
**[measured]** computed, footing stated. **[judgement]** a call for the author.

**Footing**, stated in the same breath as every measured claim, because in this project it has three
times been the whole of a disagreement: the design (schedule, draws, capacity), the priors (`k_e`,
the load offset), and the conventions (spread divisor, `t = 0` draw, nuisance draw).

**The standing footing**, assumed where an entry does not say otherwise: 25 monthly draws over two
years on `arange(0, 730, 30.4)` **including a t = 0 draw**, medians over 24 simulated intake records,
a `k_e` log-prior of 0.8, a 5% assay CV, a 5% load-record accuracy, and `sd(ln theta-hat)` against a
gate of 0.20.

**"Unverified" throughout means one thing**: the six literature reports in `refs/review-20260925/`
were never checked against their primary sources, because this container has no egress. Phase 9
checked what needed no network and stopped. Any claim resting on a number from those reports is
marked, and the specific number is named.

---

## I. The structure. What the model says before any data arrive

### 1. The steady record determines only three invariants, and the fourth parameter is free

> At a constant recorded load the increment depends on the parameters only through `(A, B, C)`, three
> invariants that three distinct marker levels determine exactly. Four unknowns against three
> invariants leave a one-parameter family of exact observational twins.

**[proven]**, `06_excitation.tex` eq. `invariants` and `fibrecl`. Footing: none needed, it is algebra.
**Unchanged by phases 8 to 10.** Evidence: `paper/sections/06_excitation.tex:28-47`.
**Depends on no unverified number.**

**Phase 10 adds a reading, not a correction** (`audit/PHASE10_PROPOSED_TEXT.md`): `C` is minus the
non-enzymatic rate constant, `A/m` the sub-saturation penalty, `B` the shortfall. The fibre then
follows in two lines. Drafted, unapplied.

### 2. The therapy's dose train is the excitation, and a constant increment removes it

> The identification rests on the variation of the increment, not its presence. A switch to
> continuous administration would silently remove the only excitation the design has.

**[proven]** on the mechanism, **[measured]** on the magnitudes. Footing: the committed 14 d train,
3.5 d half-life. **Unchanged.** Evidence: `paper/sections/07_design.tex:320-378`.
**Depends on no unverified number.**

### 3. The one-state model is the exact reduction of a two-state model, and two of its three constants are composites

> Eliminating intracellular propionyl-CoA from a two-state model returns the paper's balance
> exactly. Only `theta` survives with its meaning intact: the paper's `k_e` is `k_e/alpha`, the
> overflow disposal expressed in marker units, and its `K_m` is `K_p c / k_e_paper`, not the enzyme's
> Michaelis constant.

**[proven]**, symbolic residual 0 and 6.7e-15 over 54 numeric combinations.
**NEW IN PHASE 10.** Before: the paper treats `k_e` as the marker's clearance and `K_m` as the
enzyme's constant throughout. Evidence: `audit/q13_two_state_reduction.py` (5 of 5),
`audit/PHASE10_SENSORS.md` R.0 to R.5.
**Depends on no unverified number for the algebra.** `alpha`'s closure assumption depends on R3's
marker panel being complete, which is unverified.

**What it costs when wrong**: a saturable overflow at `K_o = 100` fits `K_m = 69.48` against a true
100, a 16% residual that **reads as a calibration error and not as a missing state**. Fourth instance
of this project's own pattern.

---

## II. The design. What a record can deliver

### 4. The gate passes at normal capacity on ordinary outpatient care

> Monthly sampling across dose cycles for two years, 25 draws, reaches `sd(ln theta-hat)` = 0.170,
> inside the 20% target. Two-monthly, 13 draws, reaches 0.188 and also passes. Six-monthly, 4 draws,
> reaches 0.457 and does not.

**[measured]**, standing footing. **Unchanged, and re-verified twice in Phase 10**: 0.1699 on the
paper's own grid, and all three six-monthly figures are generated in `tab_D_designs.tex`.
Evidence: `results/tables/tab_D_designs.tex`, `audit/q18_section8.py`, `audit/q19_loadbearing.py`.
**Depends on no unverified number.**

### 5. In the disease range the clearance must be measured in the patient

> A population prior of any plausible width fails below 20% of normal capacity. `k_e` must be known
> to about 2% at 5% of nominal, 5% at 10% and 10% at 20%.

**[measured]**, footing as standing but sweeping the `k_e` prior. **CONFIRMED AND SHARPENED IN
PHASE 9**: at the population prior there is **no statable budget on either nuisance convention**,
the per-cell envelope running 0.040 to 1.148 of the Fisher figure, so measuring `k_e` is the
precondition for the requirement framework existing at all, not merely a tightening of it.
Evidence: `audit/q2_requirement_map.py`, `audit/PHASE9_LITERATURE.md` 3a.
**Depends on R1** for the claim that no human `k_e` exists in any population, which is unverified.

### 6. The load record's requirement is per-visit scatter, not the mean offset

> Intake must be recorded to better than about 10%. Separating the channels at the 20% row, the
> per-visit scatter alone gives 0.283 and the systematic offset alone 0.172 against a 0.170 baseline,
> so 98% of the degradation is scatter and 2% is offset. A weighed-intake calibration of the mean
> buys almost nothing.

**[measured]**, standing footing at nominal capacity. **Unchanged in value and now guarded**: every
one of those numbers was a prose literal owned by a generator the build never ran, and all four are
now `\num{}` macros (R87). Evidence: `audit/t4_bias_vs_scatter.py`, `audit/q18_section8.py`.
**Depends on no unverified number.**

**One number in this cluster is unreconciled**: `06_excitation.tex:286`'s "0.288 at 20%" against the
generator's 0.286 (R88). The 1% end, 0.159, agrees.

### 7. Tier 2: the increment's repeatability is the binding requirement, and the published bound fails it

> At a 14.0% coefficient of variation the recommended record returns 0.212 against a gate of 0.20,
> and reaches the gate only at 118 monthly draws. A same-day pair must be a factor of 1.21 tighter
> than a between-visit one, a 17% reduction. Beyond an increment CV of 18.7% no schedule reaches
> 0.20 at all.

**[measured]**, standing footing at nominal. **Unchanged in value and now guarded** (R86, R87): this
whole cluster was owned by `t2_tier1_margin.py`, which declared that nothing of it reached the paper
and which the build never ran. Evidence: `audit/t2_tier1_margin.py`, `results/numbers/tier1.tex`.
**Depends on R3's between-visit repeatability figure**, the 14.02%, which is unverified; the 12.8%
widely quoted is the standard deviation of per-subject coefficients and not their mean.

### 8. Six proposed routes to more precision all fail, and they fail the same way

> Sampling harder, sampling the transient, measuring the increment, modulating the diet, recording
> the diet better, and measuring the diet instead: each attacks the parameter count, and what binds
> is `sigma` at nominal and the marker's sensitivity below it.

**[measured]** per route, **[judgement]** on the common cause. **Route 6 was retracted in Phase 10
and the retraction was then withdrawn** (R82): it is correct in all three clauses, 13 measured visits
matching 25 recorded at nominal, 6 reaching the gate, and 18 at 2% of nominal. Evidence:
`audit/t9_measured_load.py`, `audit/t8_load_by_capacity.py`, `audit/PROVENANCE_SWEEP.md` 11.3.
**Route 4's two percentages are unreconciled** (R88): the paper's 7% and 32% against the generator's
6.3% and 24.1%, on a footing its own 2.32 confirms.

---

## III. What Phase 8 added, and what it withdrew

### 9. A wrong increment shape biases the estimate but does not break it

> Peak-matched with an untreated baseline draw, the committed analyst reads capacity 7% low against a
> 2.5 d truth and 36% high against 7 d, monotonic in the truth's half-life, and the residual rises
> with the bias, so a real record would announce the misspecification.

**[measured]**. **NEW, AND IT REPLACED A LARGER CLAIM**: the review's "1.3 to 2.0 for 5 to 7 d"
becomes 1.16 to 1.36. Evidence: `audit/PHASE8_SIMULATIONS.md` S1, `audit/q1_increment_shape.py`.
**Depends on R2** for the published two-compartment kernel, unverified.

### 10. The observation model's weakest assumption is worth a factor of 5.8

> That the marker is not itself cleared enzymatically. Under the literal reading no schedule reaches
> the gate at any capacity in the disease range.

**[measured]**. **NEW IN PHASE 8**, and the first quantification of an assumption the paper states
without pricing. Evidence: `audit/PHASE8_SIMULATIONS.md` S5.
**Depends on R3** for whether the end-product reading is the right one, unverified.

### 11. The stresses do not compound

> A wrong shape, a 15% marker CV, the irregular trial train and `k_e` drift all at once give a 2 to
> 13% bias at normal capacity, and two of the stresses partly cancel.

**[measured]**. **NEW IN PHASE 8.** Evidence: `audit/q4_joint_stress.py`.
**Depends on no unverified number.**

---

## IV. What Phase 10 added about sensors and cohorts

### 12. The cohort route escapes `alpha` entirely and still cannot supply the constant

> A cohort fit estimates `k_e/alpha` directly and never meets `alpha`. On the trial's monthly record
> the disease range's looser 6% target needs 723 patients and its 2% target 6491. And at R1's
> between-patient CV of about a third, a patient's own constant floors at 0.317 with an infinite
> cohort.

**[measured]**, confirmed by fitting, fitted/Fisher 0.617 to 1.168. **NEW IN PHASE 10.** Evidence:
`audit/q14_hierarchical.py`, `audit/PHASE10_SENSORS.md` C.
**Depends on R1's "probably a third or more"**, unverified, and carries the stated assumption that
`alpha` is constant across patients or perfectly anticorrelated with `k_e`, which would make the
floor worse.

### 13. Of ten candidate sensors, five are not expressible and four of the rest are dominated

> Every observable that reads the propionyl-CoA pool has no expression in a model whose only state is
> the marker. Of what remains, propionate turnover is worth 7 to 15%, timed urine is dominated by the
> tracer study because `f_renal` sits inside it, and only ex vivo activity plus a `k_e` measurement is
> a combination worth buying.

**[measured]**. **NEW IN PHASE 10**, and one figure inside it was corrected twice (R80, R84): the
b-removal channel alone is worth under 1%, which is what the paper says in two places, and the full
load measurement is worth 7 to 15%. Both are right; they are different channels. Evidence:
`audit/q15_sensors.py`, `audit/q16_corrections.py`.

### 14. The ex vivo assay does not dominate the record at the width the literature supports

> At an assay log-sd of 1.0, which is what R4's Spearman 0.682 and `R^2` 0.327 support, the two-year
> record supplies 83% of the joint precision at nominal and 37% at 20% capacity, and the pair still
> misses the gate everywhere in the disease range.

**[measured]**. **NEW IN PHASE 10, AND IT WITHDREW ITS OWN FIRST VERSION** (R74), which was computed
at an assay width of 0.20 and concluded the opposite. Evidence: `audit/q16_corrections.py` part A.
**Depends on R4** for the width, unverified, and this is the single largest dependency in the file:
the conclusion inverts between 0.20 and 1.0.

---

## V. The machinery, which is a claim about the repository and not about the biology

### 15. 343 of the paper's numbers are checked and 2472 are not

> A number written `\num{}` is owned by a generator and compared on every build. A number typed into
> the prose is owned by nobody. 206 literals in 130 claims assert a computed result that nothing
> checks, and five generators the build never ran owned 21 of them.

**[measured]**. **NEW IN PHASE 10.** Evidence: `audit/PROVENANCE_SWEEP.md`, `audit/q17_*.py`,
`audit/q21_unreachable.py`. 22 of those literals are now macros and the five generators run (R87).

### 16. Checking found nothing wrong with the paper; connecting it to its own code found two things

> Five "does not reproduce" verdicts were the checker's, not the paper's. Wiring the generators into
> the build then produced two genuine disagreements immediately.

**[measured]** on both halves. Evidence: `audit/PROVENANCE_SWEEP.md` head, R88.

---

## A. What can be written today

Depends on nothing unverified.

| claim | where |
|:--|:--|
| 1, the three invariants and the fibre | `06_excitation.tex`, proven |
| 2, the dose train as excitation | `07_design.tex`, proven |
| 3, the two-state reduction's algebra | `audit/q13`, proven |
| 4, the gate at normal capacity | `tab_D_designs.tex`, guarded |
| 6, the load record's per-visit requirement | now guarded, one number unreconciled |
| 8, routes 1, 2, 3, 5 and 6 | now guarded |
| 11, the stresses do not compound | `audit/q4` |
| 15 and 16, the provenance findings | `PROVENANCE_SWEEP.md` |

**The Phase 10 proposed passage on `A`, `B` and `C`** is writable today except for its closing block,
which depends on the author's decision about whether the reduction enters the paper at all.

## B. What waits for the literature verification, and on which number

| claim | the number it waits on | report |
|:--|:--|:--|
| 5, `k_e` must be measured per patient | that no human `k_e` exists in any population | **R1** |
| 12, the cohort route's floor | the between-patient CV, "probably a third or more" | **R1** |
| 7, Tier 2's requirement | the between-visit repeatability, 14.02% | **R3** |
| 10, the observation model's 5.8 | whether the end-product reading is right | **R3** |
| 3, `alpha`'s closure | whether the marker panel exhausts the overflow routes | **R3** |
| 9, the shape misspecification | the published two-compartment kernel | **R2** |
| **14, the ex vivo comparison** | **the assay's width, 1.0 against 0.20** | **R4** |

**Item 14 is the one to verify first.** Its conclusion inverts between the two widths, and no other
entry in this table changes a verdict rather than a magnitude.

`refs/review-20260925/README.md` records three process warnings about these reports, including that
three of them claim branches and pull requests that do not exist.

## C. What was withdrawn, and why

| withdrawn | cause | record |
|:--|:--|:--|
| The review's 1.3 to 2.0 shape bias | recomputed at 1.16 to 1.36 against the probes | Phase 8 S1.R, R54 |
| The 0.135 Fisher figure for the offset anchor | refuted by fitting, robust spread 6.1 | Phase 9 3b.2 |
| "The record adds 0.2 to 2.9% to an ex vivo assay" | computed at an assay width R4 rules out | R74 |
| "Candidate 9 is worth 1 to 2%" | priced the offset channel only, not `sigma` | R80 |
| "S2's two load columns are not two designs" | a naming collision; the reverse is true | R80 |
| **"Route 6's eighteen visits does not reproduce"** | **the generator prints exactly 18; the sweep stopped at 5%** | **R82** |
| "Route 6's prose describes a different computation" | it names one of three simultaneous changes | R82 |
| "`2.659` is unowned and unprinted" | both `t8` and `t9` print it; the pass never ran them | R83 |
| "`07_design:113` is unowned" | generated in `tab_D_designs.tex`, merely rounded | R79 |
| "`t2` is genuine scoping work" | it owns Tier 2's whole requirement cluster | R86 |
| "Nothing pending can change the four-script finding" | a prediction about work not yet done | R86 |

**Nine of the eleven are the audit correcting itself, not the paper.** That is the shape of this
phase and the reason the provenance sweep's head says what it says.
