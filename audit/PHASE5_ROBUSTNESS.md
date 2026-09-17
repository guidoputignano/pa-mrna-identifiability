# Phase 5. Robustness to model structure

**No paper edits, no merge, no constant moved.** Step 0 was a gate. It is recorded below as
it stood when it was raised, unedited, so that what was missing and what the sources then
supplied can be read against each other. Everything from *Phase 5 continued* onward is the
work done after you supplied the five PDFs.

---

## Step 0 (GATE). What the repository holds for the two published models

Everything below was found in the repository. Nothing was fetched and nothing is inferred
from memory.

### `pkpd_translational`: Attarwala, Lumley, Liang, Ivaturi, Senn (2023)

*Translational Pharmacokinetic/Pharmacodynamic Model for mRNA-3927, an Investigational
Therapeutic for the Treatment of Propionic Acidemia*, Nucleic Acid Therapeutics 33(2):141-147,
`doi:10.1089/nat.2022.0036`. Full CrossRef metadata in `paper/refs.bib:71`.

**Text held: two verbatim sentences**, in `paper/claims.csv` C46:

> "PCC protein expression rate was modeled as a linear function of plasma PCCA/B mRNA
> concentration using a 2-compartment indirect response PD model with an empirical effect
> compartment"
>
> "PCC protein-mediated reduction of 2-MC, 3-HP, and C3/C2 ratios was modeled using a direct
> sigmoidal maximum inhibition (Imax) model"

**Numbers held**, all in `paper/refs.sources:104-112`:

| quantity | value | unit as recorded |
|---|---|---|
| PCC-protein degradation rate `kdeg` | 0.0075 | /h (half-life 3.85 d) |
| plasma 3-HP disease baseline | 60.1 (14.4% RSE) | umol/L |
| between-subject variability on that baseline | 51.1 | % |
| fully corrected floor `E0` | ~0.0987 | umol/L |
| `IC50` for 3-HP | 37.5 | ug/g PCC |
| "slope" | 55.6 | **no unit recorded** |

**One structural statement about the parameter table**, C46's note: Table 1 "carries `E0`,
`base`, `IC50` and a fixed `Imax` and no biomarker turnover constant at all."

**PROVENANCE WARNING, the repository's own** (`paper/refs.sources:113-116`):

> "these figures are transcribed from the author's parameter-grounding review of September
> 2026, **not read off the paper in this repository's own pass**. The citekey is confirmed to
> be the right article; **the numbers await a check against the primary source**."

**What is missing.**

1. **The equations.** Not one differential or algebraic equation is held. The two sentences
   name a structure; they do not write it down.
2. **The PK half.** No mRNA/LNP disposition model, no two-compartment parameters, and no
   definition of the "empirical effect compartment" or its transfer to PCC protein.
3. **The `Imax` value.** C46 says it is fixed; the value is not recorded anywhere.
4. **What "slope 55.6" is.** No unit, and three readings are open: a Hill exponent on the
   sigmoidal `Imax`; the linear slope from plasma mRNA to PCC expression rate; or a
   steepness parameter on the effect compartment. A Hill exponent of 55.6 would be a step
   function, so the reading changes the model qualitatively. **This is the single most
   blocking gap.**
5. **Units on `IC50`.** "37.5 ug/g PCC" does not say ug of what per g of what, nor in which
   compartment, so it cannot be put into the paper's units.
6. **Whether the three markers share parameters** or each carries its own set, and whether
   `base` is per-marker.
7. **The residual-error model** and whether the 51.1% is on the baseline only.

### `dose_selection_poster`: Gao, He, Li, Liang, Luo, Natenshon, Patel, Real, Van, Ward, Zhang

*Model-Informed Dose Selection for mRNA-3927 in Participants with Propionic Acidemia (PA)*,
International Society of Pharmacometrics, `doi:10.70534/JVQU8233`. Full CrossRef metadata in
`paper/refs.bib:73`. **It is a conference abstract.**

**Text held: one verbatim sentence**, in `paper/claims.csv` C46:

> "An indirect inhibition Kin model was used to describe the PD response of 3-HP and 2-MC."

**Facts held** (`paper/refs.sources:119-121`, and claims C17 / C18): 0.6 mg/kg biweekly,
selected by simulation in 1000 virtual participants aged 0-36 y over body-weight
distributions, with 3-HP and 2-MC as the response; body weight is a covariate; **genotype is
not**.

**What is missing: everything structural beyond the model-class name.** The repository says so
itself, in C46's note:

> "it carries no PD equations and no biomarker parameters"

So there is no `kin`, no `kout`, no inhibition function, no `Imax` or `IC50`, no baseline, no
units, and no random-effect structure.

---

## What this means for Steps 1 to 4

**Steps 3 and 4 do not need the missing material.** Identifiability of a structure, and the
theorem about when a quasi-steady rational balance becomes linear in a reparametrisation, are
questions about the **form** of the equations. The forms are named unambiguously enough to
write down: (A) turnover with inhibited production and first-order loss; (B) indirect response
on an enzyme state with an algebraic saturating map to the marker.

**Step 2 does need it**, or an explicit substitute. Generating data from A and B and fitting
our model requires numbers. With what is held I can only supply them myself, and that would
make the misspecification test a test of **my** choice of parameters rather than of the
published models.

**One further gap, not about the references.** `StructuralIdentifiability.jl` is **not
installed in this container** (`julia` is absent; `symbolic/julia/` holds the two scripts from
earlier phases but they cannot be run here). Step 3's transient case can be done in `sympy` by
resultants and Groebner elimination, which is what Phase 0 used and what `f5`/`f6` already do,
but it will not be the same tool the paper cites. Say if you want Julia installed.

---

## The ask

**Stopping here, as the gate requires.** *(Answered. You supplied all of it; see finding 1
and finding 5 below. The alternative offered at the end of this section is what was run:
A and B calibrated to this paper's operating point, flagged in every table.)* To proceed
properly I need:

1. **`pkpd_translational` full text**, and above all **the model-equations section and Table 1**
   with units. The blocking item is what "slope 55.6" denotes.
2. **`dose_selection_poster`**: the poster itself if anything exists beyond the abstract. If
   the abstract is all there is, say so and I will treat (A) as *a named structure with no
   published parameters*.

**An alternative, if you would rather not wait.** I can run Steps 1, 3 and 4 now on the two
structures as named, and run Step 2 with A and B **calibrated to this paper's own operating
point** rather than to published values, flagged in every table as my parameterisation and not
theirs. That still answers the question Step 2 asks, which is whether the design rules and the
`theta` estimate survive a change of structure. It does not answer whether they survive *the
published models*, and the report would say so in those words.

Tell me which, and I will proceed.

---

# Phase 5 continued. The five source checks, and steps 1 to 2

**Status of the gate.** The gate is cleared. You supplied the Attarwala full text with
Supplementary Table S1 and Figures S1 to S3, the Nature trial paper with its extended data,
and Baek et al. 2024. Every blocking item listed above is now answered from the primary
sources, and the answers are below. All five PDFs were read in this pass with `pdftotext
-layout`, and Extended Data Table 3 was read off a 200 dpi render of page 11 because it is a
graphic in the text layer.

**Unchanged rules.** No paper edits, no merge, no constant moved. Every correction below is
*proposed and unapplied*. `PHASE4_PROPOSED_TEXT.md` stays unapplied. No scheduled check-ins
were armed and no agents were spawned.

**Scripts.** `audit/k1_sources.py` (the five checks and the reruns), `audit/k2_structures.py`
(steps 1 and 2), `audit/k3_recheck.py` (the independent re-check, 6 of 6 pass).

Findings 2, 3 and 4 come first, as you asked, because they bear on the paper's physical
relevance whatever the robustness test says.

---

## Finding 2 (PROVEN from the source). The Table 1 PD estimates are mouse, and four repository rows read them as human

**What the source says.** The PD half of Attarwala's model is fitted in
`Pcca-/- (A138T)` mice, `N = 111`. The section is headed *"PK/PD model for plasma 2-MC, 3-HP,
and C3/C2 ratio in PCC-deficient mice"*. The parameter definitions say the baselines are
*"plasma levels of corresponding biomarkers that are susceptible to suppression by PCC protein
expressed through mRNA-3927 in mice"*. Supplementary Table S1 scales on `(BW/0.025)`, a 25 g
reference animal.

**And the scaling sentence covers PK only.** The paper's own words:

> "The preclinical PK/PD model described earlier was extrapolated to humans based on
> allometric scaling of model parameters. Based on the PCCA/B mRNA PK model described earlier,
> estimated allometric exponent was 0.63 for CL12 and CL32 and 1.1 for CL23 and CL20.
> Allometric exponent for the volume of distribution parameters (ie, V and V2) was fixed to 1."

Every parameter named there is a clearance or a volume. Baek's companion paper says the same in
one sentence: *"Each preclinical PK/PD model was extrapolated to humans based on allometric
scaling of PK parameters."* So `E0`, `base`, `IC50`, `Imax`, `Slope`, `kdeg`, `kq` and `ke0`
are mouse estimates, and the `51.1%` and `73.8%` inter-individual variabilities are between
mice of one inbred strain.

**Where the repository reads them as human.** `results/tables/tab_params.tex` is generated
from `src/common/params.py` by `make figures`, so each row below is one edit and not two.

| place | what it says | why it is wrong |
|---|---|---|
| `params.py:60-63` → `tab_params.tex:19` (the `$\bar m$` row) | "the translational model puts the untreated plasma 3-HP baseline at $60.1$ $\mu$M and the fully corrected floor near $0.1$ $\mu$M, so a managed **patient** sits between the two" | both numbers are mouse. The floor is `E0 = 0.0987`, also mouse. |
| `params.py:107-113` → `tab_params.tex:26` (the `$\lambda$` row) | "**Corroborated**: the translational PK/PD model estimates a PCC-protein degradation rate of $0.0075$ h$^{-1}$, a half-life of $3.85$ d, so this is within $10\%$ of a published value" | mouse, and see finding 3: `3.85` d is one rate constant inside the model, not the model's protein half-life. |
| `params.py:135-136` → `tab_params.tex:30` (the `CV` row) | "Biological variability is far larger, the translational model carrying $51\%$ between subjects on the 3-HP baseline" | `51.1%` is between **mice**, so if anything it is a lower bound on between-patient variability, and the sentence deploys it as the large figure. |
| `params.py:82-83` (the `$k_e$` row) | "The nearest route is the first-order 3-HP elimination implicit in the published indirect-response model~\cite{pkpd_translational}" | not a species error but a **structure** error. Supplementary Table S1(C) is algebraic in the marker, `3-HP = E0 + base*(1 - Imax*PCC/(IC50+PCC))`. There is no marker state and no elimination constant anywhere in that model, so nothing can be extracted from it. `claims.csv` C46 already records this ("no biomarker turnover constant at all"). The model that *does* carry a first-order marker loss is the other one, and it publishes no value for it. |
| `refs.sources:102-116` | "fitted largely to preclinical data and **allometrically scaled to human**", plus the provenance warning | true of the PK, false of the PD, and the six figures beneath it are all PD. |

**Proposed corrections, unapplied.** Four annotation edits and one source note. No value moves:
`60.1`, `0.0075` and `51%` are all prose, and `M_BAR` stays `10`.

1. `$\bar m$` row: after "translational model", insert "of PCC-deficient **mice**"; and change
   "so a managed patient sits between the two" to "so a managed patient is read on a scale of
   the same order, with the caveat that the two anchors are murine. The nearest human anchor is
   the newborn dried-blood-spot 99.9th centile of $38.13$ $\mu$mol/L~\cite{dbs_3hp_cutoff},
   which is within a factor of $1.6$ of the mouse baseline."
2. `$\lambda$` row: replace "**Corroborated**: ... within $10\%$ of a published value rather
   than free" with "**Partly corroborated, and only early.** The translational model of
   PCC-deficient mice estimates `kdeg` $= 0.0075$ h$^{-1}$~\cite{pkpd_translational}, but that
   is one rate constant inside a two-compartment protein model whose own stated protein
   half-life is $\sim 7$ d. The model's local half-life passes $3.5$ d at about five days after
   a dose and reaches $9.0$ d by the end of a fortnightly interval, so $\lambda$ tracks the
   published decline over the first week and decays faster than it thereafter."
3. `CV` row: change "the translational model carrying $51\%$ between subjects" to "the
   translational model carrying $51\%$ between **mice** of one inbred strain, which bounds
   between-patient variability from below rather than above".
4. `$k_e$` row: replace "the first-order 3-HP elimination implicit in the published
   indirect-response model~\cite{pkpd_translational}" with "the first-order marker loss of the
   turnover model~\cite{dose_selection_poster,mrna_moa}. The translational
   model~\cite{pkpd_translational} does **not** supply one: its marker is an algebraic function
   of the enzyme with no state and no elimination constant."
5. `refs.sources`: change "allometrically scaled to human" to "whose **PK** parameters were
   allometrically scaled to human; the PD parameters below were **not** scaled and are mouse
   estimates", and **lift the provenance warning**: all six figures are confirmed against
   Table 1 in this pass, with two unit misprints noted under finding 1.

---

## Finding 3 (MEASURED). "3.85 d" is not the published protein half-life, and the corroboration is early-only

Supplementary Table S1(B) gives the protein two compartments:

```
d(PCC)/dt = Ksyn + Kq*(PCP - PCC) - Kdeg*PCC
d(PCP)/dt = Kq*(PCC - PCP)
```

with `Kdeg = 0.0075 /h` and `Kq = 0.00474 /h`, both from Table 1. (S1 as printed writes
`Ksyn*Ce` with `Ksyn = Ce*Slope` on the next line, which is quadratic in `Ce` and
dimensionally wrong. The text says synthesis rises *linearly* in mRNA, so the intended term is
`Ksyn` alone. The eigenvalues do not depend on which reading is taken, because the input enters
the forcing and not the matrix.)

| reading | value |
|---|---|
| `ln2/Kdeg` alone, which is what the repository quotes | **3.85 d** |
| fast eigenvalue | 1.99 d |
| slow (terminal) eigenvalue | 11.81 d |
| peak to half-peak of the simulated single-dose profile | 2.79 d |
| local half-life 14 d after the peak | **9.04 d** |
| the paper's own statement | **~7 d** |

No single constant reproduces the `~7 d`. The protein declines biexponentially, which is
exactly what the text says the second compartment was added for, so its half-life is a function
of time: `2.54 d` at one day after the peak, `3.53 d` at five days, `6.34 d` at ten, `9.04 d` at
fourteen, tending to `11.81 d`. It passes `7 d` at about `10.9 d` after the peak. **The `~7 d`
is consistent with the constants, read from the middle of that decline. There is no arithmetic
disagreement with the source.**

What this does to `$\lambda$`: the local half-life passes `3.5 d` at about **4.9 d** after the
peak, so over roughly the first week `$\lambda$` tracks the published model closely. Over the
window `$\lambda$` has to govern, a whole 14 d dose interval, the model's local half-life is
`9.04 d` (a ratio of `2.6`) and the terminal figure is `11.81 d` (a ratio of `3.4`). A
single-exponential increment with a 3.5 d half-life therefore **decays faster late in the
interval than the published model's protein does**, which is the direction that makes the
inter-dose excursion look larger than it is. What the annotation cannot support as written is
that `3.85 d` is "a half-life of PCC protein".

---

## Finding 4 (PROVEN from the source). The trial reports a 20 to 29% fall in two of five cohorts, not "about half"

Extended Data Table 3 of the Nature paper, read off a render of page 11. 3-HP, mean percentage
change from baseline after the third dose:

| cohort | regimen | n | 3-HP mean % change (SD) |
|---|---|---|---|
| 1 | 0.30 mg/kg Q3W | 4 | **-29 (25)** |
| 2 | 0.30 mg/kg Q2W | 3 | **-3 (26)** |
| 3 | 0.45 mg/kg Q2W | 3 | **-20 (24)** |
| 4 | 0.60 mg/kg Q2W | 3 | **+5 (57)** |
| 5 | 0.90 mg/kg Q2W | 3 | **-20 (41)** |

These match your figures exactly. The mean of the cohort means is `-13.4%`; the range runs from
`-29%` to `+5%`; there is no dose trend, the largest fall being in the lowest-dose cohort.

Against this the repository's `D_AMP = 6000` is calibrated to a fall of `64.90%` at the
first-cycle peak, `66.36%` at the steady-state peak, and `62.25%` from the patient's own
pre-dose trough to the following peak, which is the `62-66%` already recorded unapplied in
`PHASE4_PROPOSED_TEXT.md`. **That is two to three times the reported falls, and the reported
falls are not "about half".**

**The assumption about when "after the third dose" was sampled, stated explicitly.** *The
Nature paper does not say.* Extended Data Table 3's caption carries no time point; the Methods
define the endpoint only as "changes in blood 2-MC and 3-HP levels from baseline (predose
levels) to post-dose levels"; the Reporting section names no biomarker sampling schedule. The
one sampling time the paper does name anywhere is a **predose** draw ("the predose sample for
dose 3 (52 days after the first dose)"), which is a trough. I therefore carry **both**
readings throughout, and the choice turns out to matter more than the fall itself.

**What the between-patient SD implies.** The SDs are 24 to 57 percentage points on cohorts of
three or four, against mean changes of 3 to 29 points. Two consequences, both for the paper's
noise model rather than its conclusions:

- The `5%` analytical assay CV in `params.py` is nowhere near the whole story. Even the
  smallest cohort SD, 24 points, is five times it. The repository already says so in the `CV`
  row's annotation, and finding 2 corrects the figure it says it with.
- The SD is **between patients**, over three or four people, so it confounds assay repeatability,
  within-patient day-to-day biological variability, and genuine between-patient response
  heterogeneity, and nothing in the table separates them. It is an upper bound on each. With
  `n = 3` the SD itself has about `50%` relative uncertainty, so `24` and `57` are not
  distinguishable from one another. What can be said is that a single 3-HP measurement carries a
  variability of *at least* the same order as the treatment effect it is meant to show, which is
  the paper's own argument for repeated visits and is consistent with `LOAD_CV = 0.10` being an
  optimistic figure rather than a pessimistic one.

---

## The reruns

Every number comes from `g_design_common`, unchanged, with `design.run.LAMBDA` and `D_AMP`
rebound in a context manager that restores them. Nothing in `src/` is edited and no committed
constant moves. `n` is the number of monthly visits to reach `sd(ln theta) <= 0.20`, searched to
`n = 60` (five years). The `sd` column is at `n = 12` with the 8-seed spread.

First, the dose that produces each target fall, at the committed half-life:

| reading of "after dose 3" | target fall | `D` | `Delta` at the sample |
|---|---|---|---|
| peak | 20% | 766.0 | 816.8 |
| peak | 30% | 1311.1 | 1398.1 |
| trough (next predose) | 20% | 12255.3 | 816.8 |
| trough (next predose) | 30% | 20976.9 | 1398.1 |

Then the design numbers:

| configuration | peak `Delta/theta` | `Delta/theta` range, monthly n=12 | `sd(ln theta)` at n=12, CV 10% | visits to 0.20 |
|---|---|---|---|---|
| **as committed** (3.5 d, `D = 6000`) | 2.13 | 1.8170 | 0.1888 [0.1857, 0.1919] | **11** |
| finding 3 alone: 7 d half-life | 2.67 | 1.7820 | 0.2657 [0.2597, 0.2780] | **30** |
| finding 4: 20% fall at the **peak** | 0.27 | 0.2320 | 0.3682 [0.3591, 0.3821] | **never by 60** |
| finding 4: 30% fall at the **peak** | 0.47 | 0.3970 | 0.2744 [0.2669, 0.2805] | **42** |
| finding 4: 20% fall at the **trough** | 4.36 | 3.7112 | 0.1835 [0.1808, 0.1865] | **10** |
| finding 4: 30% fall at the **trough** | 7.46 | 6.3523 | 0.1879 [0.1853, 0.1921] | **12** |
| findings 3+4: 20% at the peak, 7 d | 0.28 | 0.1848 | 0.5154 [0.4912, 0.5429] | **never by 60** |
| findings 3+4: 30% at the trough, 7 d | 1.89 | 1.2655 | 0.2639 [0.2578, 0.2756] | **30** |

**Three things this says.**

1. **The sampling-time assumption dominates, and we do not have it.** Under the trough reading
   the trial's 20 to 30% fall needs a *larger* dose than the repository assumes, not a smaller
   one, because the trough is `6.7%` of the peak increment. The design then gets slightly
   **better** than committed, 10 to 12 visits against 11. Under the peak reading it gets much
   worse, 42 visits at 30% and unreachable within five years at 20%. The two readings of one
   unstated sampling time move the answer from "better than we said" to "not feasible". This is
   the single largest uncertainty the reruns expose, and it is a gap in the trial report, not in
   the repository.
2. **The 7 d half-life alone roughly triples the record needed**, from 11 monthly visits to 30,
   with the `Delta/theta` range almost unmoved (1.817 to 1.782). So it is not the excursion size
   that changes; it is where in the cycle the excursion sits relative to a monthly visit grid.
3. **Load CV still barely moves any of it**, which is the prior-conditional finding already
   recorded in Phase 3: the load prior is carrying that comparison, and the seed spread is
   narrower than the gap between rows in every configuration.

---

## Finding 1 (PROVEN from the source). The transfer function: a linear slope, a fixed `Imax`, and two misprinted units

**(a) `Slope 55.6` is a linear coefficient, not a Hill exponent.** This was the gate's single
most blocking item and it is settled. Supplementary Table S1(B) gives `Ksyn = Ce * Slope`, and
the text says *"indirect response PD model where the synthesis rate of PCC protein increased
linearly with an increase in plasma PCC mRNA concentration through an effect compartment"*.
There is no exponent in the PCC equation at all. The only saturating term in the entire model is
the `Imax` expression on the marker, in part (C). The RSE is `1.93%`.

**(b) `Imax` is fixed at `0.999`.** Table 1, footnote m: *"Maximum inhibition (Imax) was fixed at
0.999."* So the marker's floor is `E0`, not a fitted maximum effect: as PCC grows without bound,
`3-HP -> E0 + base*(1-Imax) = 0.0987 + 0.0601 = 0.1588`.

**(c) The slope's unit as printed is `mg`, not `ug`.** Table 1 prints
`Slope [(mg/g/h)/(mg/mL)]`. The repository's `refs.sources` gives the `IC50` in `ug/g`. The ratio
inside the slope is scale-free if both masses carry the same prefix, so `55.6` is unchanged
either way, but the annotation should quote the unit as printed. The repository currently records
**no unit at all** for the slope, which was gap 4 at the gate.

**(d) `"mg/g"` in the `IC50` rows is a misprint for `ug/g`, on the paper's own evidence.** The
Fig. 2B axis label, read in this pass, is `"PCC protein, ug/g"`, with plotted liver PCC running
to about 700 on that axis. The three `IC50`s are printed `21.0 / 37.5 / 32.1` "mg/g". On the
figure's own axis those sit inside the plotted range; read as `mg/g` they would be a
thousandfold above every point the model was fitted to and no dose could reach half effect.
**The repository already has this right.**

**(e) `"mM/L"` in the marker rows is a misprint for `umol/L`.** Three independent reasons. The
3-HP baseline is printed `60.1` "mM/L", which as millimolar is a plasma organic acid at 60 mM.
The same text gives the assay LOQ as `25` "mM/L", and an LOQ at `0.42` of the disease baseline is
an ordinary assay whereas an LOQ at 25 mM is not an assay at all. The repository's own population
anchor, the newborn dried-blood-spot 99.9th centile, is `38.13 umol/L`, a factor `1.58` from
`60.1`. **The repository already has this right.** "mM/L" is in any case not a unit, since mM is
already per litre.

**Proposed addition to `refs.sources`, unapplied**, folding (a) to (e) into the
`pkpd_translational` annotation and replacing the provenance warning: record the slope as
`55.6 (mg/g/h)/(mg/mL) as printed`, record `Imax fixed at 0.999`, record that the two unit
strings are misprints with the evidence above, and record that the marker equation is algebraic.

---

## Finding 5 (PROVEN from the source, with one judgement). What Models A and B are

**Model A is a named structure with no published parameters.** It is named twice and written
down nowhere:

- Gao (ACoP2025, `dose_selection_poster`), quoted in `claims.csv` C46: *"An indirect inhibition
  Kin model was used to describe the PD response of 3-HP and 2-MC."*
- Baek et al. 2024, Methods, read in this pass: *"For PD biomarkers, classic indirect response
  models with inhibition of kin were developed for 3-HP and 2-MC/MMA, incorporating PCC/hMUT
  enzyme-induced decreases in organic acid accumulation."*

Baek's **Table 1 carries no PD parameter values at all**. It is a dose-projection table: LNP
formulation, starting dose, FIH dose species, NOAEL, safety-relevant species, safety margin. So
it cannot supply numbers either. **Every number in our Model A is therefore ours, calibrated to
this repository's operating point, and is flagged as such in every table produced by
`k2_structures.py`.**

**Registration: Baek is already a source here.** It is `mrna_moa`,
`doi:10.1038/s41467-024-47460-9`, at `refs.sources:200-203`. No new citekey is needed, and
adding one would create a duplicate. What is missing is the structure sentence above: the
repository currently attributes Model A's structure to the conference abstract alone, when the
same structure is stated in a peer-reviewed Methods section it already cites.

**Proposed addition to the `mrna_moa` annotation, unapplied:** add the quoted Methods sentence
and the note that Table 1 is a dose-projection table carrying no PD parameters; and add
`mrna_moa` alongside `dose_selection_poster` at the two places in the paper that attribute the
turnover structure (`sections/03_dynamic.tex:176-179`, `sections/A1_supplementary.tex:41-44`).
This strengthens a citation rather than correcting one: C46 records the Gao attribution as
`matches-with-drift` precisely because an abstract was carrying it alone.

**Model B is Attarwala's structure with Table 1's values, unscaled**, because the paper rescales
nothing in this block to human (finding 2). Supplementary Table S1(B) and (C):

```
d(Ce)/dt  = ke0*(C13 - Ce)
d(PCC)/dt = Ksyn + kq*(PCP - PCC) - kdeg*PCC ,   Ksyn = Slope*Ce
d(PCP)/dt = kq*(PCC - PCP)
3-HP      = E0 + base*(1 - Imax*PCC/(IC50 + PCC))
```

**The structural point that step 2 turns on: Model B's marker is algebraic in the enzyme.** It
has no state, so it cannot express a substrate load at all, and neither can Model A as named.
Where the load enters is our choice in both, and that is stated in every table.

---

# Step 1. The two published structures as simulators

`audit/k2_structures.py`. A naming trap is stated at the top of that file and is repeated here:
this repository already has `src/model_A/` (the paper's local model) and `src/model_C/` (the
parent), and they have nothing to do with the two published structures. The file writes `PUB_A`,
`PUB_B` and `PAPER` throughout.

**Where the load enters, flagged as our assumption.** Neither published structure has a
substrate load in it. Both were fitted to a fixed-diet animal colony, so the load never varied
and nothing in either paper says where it would go. We choose:

- `PUB_A`: the load is the marker's **production rate**, `k_in = L + b`. This is the
  least-committal choice available, because `k_in` already is a production rate and the load
  already is a rate of substrate arrival.
- `PUB_B`: the marker is an instantaneous function of the enzyme, so a load can only be imposed.
  We scale the whole marker by `(L+b)/L0`, that is, at fixed enzyme the marker is proportional to
  substrate supply. This is the weakest assumption that gives the structure a load at all.

**What is published and what is ours, per constant.**

| constant | source |
|---|---|
| `Imax = 0.999` | **published**, Attarwala Table 1 note m; dimensionless, transfers unchanged |
| `E0/base = 0.0987/60.1` | **published**, Attarwala Table 1; a ratio, so it transfers unchanged |
| `IC50` on our enzyme axis | **ours**. `37.5 ug/g` is a liver protein concentration and `theta` is a flux capacity, and the untreated PA-mouse PCC level that would fix the ratio is not published |
| `k_out` (PUB_A only) | **ours**, and it follows once `IC50` is fixed by the operating point |

**The calibration has two roots, and both are carried.** Matching one fall at one enzyme level is
a quadratic in `IC50`, and both roots reproduce the target exactly, so no conclusion below rests
on a branch choice. The two branches differ enormously in what they imply for the marker's own
memory: at a 30% calibration `PUB_A` is either `k_out = 0.058/d` (marker half-life 12 d) or
`k_out = 25.8/d` (39 minutes).

**A structural result that falls out of the calibration, before any fitting.** The published
inhibition shape has a **hard ceiling** on the fall it can produce between the untreated enzyme
and the treated steady-state peak, maximised over all `IC50`:

| structure | ceiling on the fall | at `IC50` |
|---|---|---|
| `PUB_A` (no floor term) | **65.65%** | 167.9 |
| `PUB_B` (with the published `E0/base` floor) | **64.13%** | 272.7 |

The repository's committed calibration is `64.90%` at the first-cycle peak and `66.36%` at the
steady-state peak. **The steady-state figure is above both ceilings and the first-cycle figure is
above `PUB_B`'s.** So with Attarwala's own `Imax` and floor ratio, neither published structure can
reproduce the fall the repository assumes, on any `IC50`, at a `3.13`-fold enzyme increase. This
is an independent argument, from the published shape rather than from the trial, that `D_AMP` is
set too high. The reruns above are the quantitative version.

---

# Step 2. Fit the paper's model to records the published structures generate

Five generators (`PAPER` as the control, `PUB_A` and `PUB_B` on both calibration branches), five
schedules at `n = 12`, load CV `0`, `2%` and `10%`, three calibration targets (`60%`, `30%`,
`20%`), 8 daily diet records per cell. The fit is `i1_estimator`'s variable-projection scan
polished on the exact objective, unchanged. Full tables in `k2_structures.py`'s output; the
`60%` calibration is quoted here because it is nearest the repository's committed dose.

**The control row works.** `PAPER` returns `rms/sigma = 0.000` and `theta_hat/theta = 1` on every
schedule and every load CV, except `14 d, already on therapy` at constant load, where the fit is
still exact but `theta_hat` collapses to the estimator's lower clip. That is P5's exactly
singular case, and it is reproduced here from a completely separate code path.

## (a) Is the misspecification detectable? Structure only, no observation noise

RMS residual of the paper's best fit, in units of the repository's own `sigma = 21.9612`, with
the noise switched off, so the number is structure and nothing else. `theta_hat/theta` beside it.
Load CV `0`. In the full tables the `2%` and `10%` columns keep every residual inside the
8-seed range of the CV `0` cell. `theta_hat` holds too, with one exception worth its own line:
on `14 d, already on therapy` **any** load variation breaks the exact singularity, and what the
fit then returns depends on the structure. Under `PAPER` it recovers the truth
(`theta_hat/theta = 1`). Under `PUB_B` it jumps from the lower clip to about `1e6` times the
truth. Under `PUB_A` it stays at the clip. So load variation rescues that schedule only when the
model is right, which is a caveat the paper does not currently carry.

| generator | 14 d | 14 d, on therapy | 14.5 d | monthly | phase-spread |
|---|---|---|---|---|---|
| `PAPER` | 0.000, `1` | 0.000, clip | 0.000, `1` | 0.000, `1` | 0.000, `1` |
| `PUB_B` IC50 63.75 | **0.000**, `2.03` | 0.000, clip | **0.068**, `1.53` | **0.180**, `1.34` | **0.170**, `1.32` |
| `PUB_B` IC50 1167 | **0.000**, `1.34` | 0.000, clip | **0.006**, `1.31` | **0.016**, `1.29` | **0.015**, `1.29` |
| `PUB_A` IC50 22.93 | 0.624, clip | 0.000, clip | **2.051**, `2.3e7` | **1.961**, `2.34e7` | **1.960**, `2.31e7` |
| `PUB_A` IC50 1230 | 0.134, clip | 0.000, clip | **2.591**, `1.69e7` | **2.469**, `12.1` | **2.423**, `16.7` |

"clip" is `theta_hat = 1e-9`, the estimator's own lower bound: the fit wants `theta <= 0` and
cannot go there.

**Two different answers, and the difference is the marker's memory.**

- **`PUB_B` is invisible.** Attarwala's algebraic structure produces records the paper's model
  fits to between `0.000` and `0.19` of one `sigma` at this calibration. At the repository's
  noise level nothing in the residuals could tell the two structures apart, on any schedule.
  **But `theta_hat` comes out `1.29` to `2.03` times the truth**, and `K_m` moves with it: at the
  `30%` calibration the fitted `K_m` runs `462` to `979` against a true `100`. So the paper's
  model absorbs `PUB_B` completely and returns a wrong answer with no warning in the fit.
- **`PUB_A` is loud, and loudest on exactly the schedules the paper recommends.** The turnover
  structure gives an RMS of `2.0` to `2.6` sigma on `14.5 d`, monthly and phase-spread, and
  `theta_hat` runs away by six to seven orders of magnitude into the linear limit where
  `theta g(m)` merges with `k_e m`. On the two `14 d` schedules the residual is small (`0.000` to
  `0.62`) and `theta_hat` collapses to zero instead.

**And the `PUB_B` bias grows sharply as the dose falls, which ties this straight to finding 4.**
`theta_hat/theta` under `PUB_B`, both branches, at load CV `0`, with the monthly residual beside
it:

| calibration | `14 d` | `14.5 d` | monthly | phase-spread | monthly `RMS/sigma` |
|---|---|---|---|---|---|
| 60% fall (near the committed `D`) | 1.34 to 2.03 | 1.31 to 1.53 | **1.29 to 1.34** | 1.29 to 1.32 | 0.016 to 0.180 |
| 30% fall | 4.30 to 11.4 | 4.29 to 6.67 | **4.29 to 4.81** | 4.29 to 4.67 | 0.001 to 0.248 |
| 20% fall | 7.30 to 20.5 | 7.29 to 11.7 | **7.29 to 8.30** | 7.29 to 8.04 | 0.000 to 0.179 |

The residual does not grow with the bias. It shrinks. So **if the trial's 20 to 30% fall is the
right scale, a `PUB_B` truth would make the paper's estimator wrong by a factor of four to
twenty, and the fit would look better than at the dose the repository currently assumes.** That
is the single most consequential number in this report, and it depends on the same unstated
sampling time as everything else in the reruns.

**The `PUB_A`/`PUB_B` contrast is the other useful result.** The schedules the paper recommends are the ones that
*expose* a turnover marker; the schedule locked to the dose period hides it and quietly returns
`theta = 0`. So the paper's design advice does two jobs at once, and the second was not claimed:
a phase-spread or off-period record is not only better conditioned, it is the only one of these
five that would have told the analyst the structure was wrong.

## (b) With the repository's own observation noise

`sd(ln theta)` as the fitter would report it, median over 8 seeds (which carry both the diet
record and the noise draw). `sing` means every seed returned a singular information matrix at its
own fit, so no precision can be reported at all.

| generator | 14 d | 14.5 d | monthly | phase-spread |
|---|---|---|---|---|
| `PAPER` | 7.198 | 0.2885 | **0.1893** | **0.1894** |
| `PUB_B` IC50 63.75 | 4.990 | 0.2861 | **0.1894** | **0.1900** |
| `PUB_B` IC50 1167 | 8.964 | 0.2912 | **0.1952** | **0.1955** |
| `PUB_A` IC50 22.93 | `sing` | `sing` | `sing` | `sing` |
| `PUB_A` IC50 1230 | 355.1 | `sing` | 0.6936 | 1.0527 |

Load CV `0`. Across `2%` and `10%` these move by under `2%` in the `14.5 d`, monthly and
phase-spread columns. The `14 d` column is not stable in any structure (`PAPER` runs
`7.198 / 7.063 / 4.454` across the three CVs), which is the near-singular schedule behaving as
the paper says it does.

**The schedule ranking survives both structures.** Monthly and phase-spread beat `14.5 d` beats
`14 d`, in that order, under `PAPER` and under both `PUB_B` branches, with the monthly and
phase-spread figures agreeing to within about `3%` across structures (`0.1893`, `0.1894`,
`0.1952` on the monthly record). Under `PUB_A` the
reported precision is either unavailable or four times worse, and the ranking within the usable
cells is unchanged. **So the design advice transfers; the parameter estimate does not.**

A calibration note on the residual column of the full tables: a well-specified fit of `n = 12`
points for 4 parameters has an expected `RMS/sigma` of `sqrt(8/12) = 0.816`, and the
well-fitting rows sit at `0.78` to `0.87`. That is the figure to read the residual column
against, and it confirms the fit is absorbing the noise rather than the structure.

## (c) Does the fibre survive a change of structure?

The width of the `K_m` profile within 1 unit of the minimum of the same penalised objective the
estimator minimises, at a constant recorded load and no noise. The scan runs `0.3` to `3e4`;
"edge" means the flat region runs off the end of it. This is the paper's central obstruction
made measurable, and `A1_supplementary.tex:51-52` asserts exactly what is being tested here:
*"The obstruction survives a change of structure; this parametrisation of it does not."*

| generator | 14 d | 14 d, on therapy | 14.5 d | monthly | phase-spread |
|---|---|---|---|---|---|
| `PAPER` | 35.6 to edge | **0.3 to edge** | 88.5 to 114.7 | 91.7 to 109.9 | 91.7 to 110.2 |
| `PUB_B` IC50 63.75 | 46.6 to edge | **0.3 to edge** | 120.8 to 164.9 | 118.4 to 144.4 | 117.0 to 143.6 |
| `PUB_B` IC50 1167 | 44.1 to edge | **0.3 to edge** | 115.7 to 154.3 | 119.4 to 145.7 | 119.1 to 146.1 |
| `PUB_A` IC50 22.93 | 36.5 to 77.8 | **0.3 to edge** | 1.13e4 to edge | 1.07e4 to edge | 9877 to edge |
| `PUB_A` IC50 1230 | 68.7 to 538 | **0.3 to edge** | 5243 to edge | 858 to 3680 | 1014 to 1.18e4 |

**Both halves of the supplement's sentence are confirmed, and one of them more strongly than the
paper claims.**

- **The obstruction survives.** `14 d, already on therapy` is flat across the **whole scan** for
  **every** generator, on both branches of both structures. The exactly singular case is not an
  artefact of the paper's rational parametrisation; it is a property of a record whose samples all
  sit at the same dose-cycle phase, and it holds when the marker is algebraic in the enzyme and
  when it is a turnover state. `14 d` from `t = 0` runs to the edge for `PAPER` and both `PUB_B`
  branches, and breaking the period lock bounds `K_m` for `PAPER` and both `PUB_B` branches, in
  the same direction and by a comparable factor.
- **The parametrisation does not survive.** The bounded intervals sit in the wrong place. Under
  `PUB_B` the fitted `K_m` is `116` to `165` against a true `100` at the `60%` calibration, and
  `354` to `1086` at the `30%` calibration. Under `PUB_A` the flat region is three to four orders
  of magnitude out. A bounded interval is not a correct interval.

---

# The independent re-check

`audit/k3_recheck.py`. Every check reaches the same number by a route that shares no code with
the script it checks. **6 of 6 pass**, and two of them failed first and were the reason two real
errors were found.

| check | independent route | result |
|---|---|---|
| 1. PCC eigenvalues | exact `Fraction` characteristic polynomial, `Decimal` at 60 digits, against numpy's eigensolver | agree to `1e-9` |
| 2. `IC50` calibration roots | the closed form of the same quadratic, against the grid bisection `k2` uses | agree to `1e-6` on both structures at three targets |
| 3. Structure ceilings | stationarity of the ratio by derivative bisection, against `k2`'s grid maximisation | agree to `1e-7` |
| 4. `D` calibrated to a target fall | bisection on the balance itself, against `k1`'s quadratic root formula | falls reproduce to `1e-8`; the untreated level returns `10.0000000000` |
| 5. `PUB_A`'s integration | variation of constants with Simpson quadrature on a dose-aligned grid, against `solve_ivp` | now agree to `3e-10` (slow branch) and `4e-8` (fast branch) |
| 6. Control row | the paper's model fitted to its own data on every schedule | `rms/sigma < 3e-13`, `theta_hat/theta = 1`, `K_m = 100`, except the two constant-load fibre cases |

**Check 5 failed twice, and both failures were real.**

1. The first failure (errors of `2e-4` to `5e-2` relative) was in **`k2`**: `PUB_A` was being
   integrated with one `solve_ivp` call across the whole record, and the right-hand side jumps at
   every dose time and every day boundary. An adaptive Runge-Kutta run straight through a jump
   loses its order there. `PUB_A.marker` now restarts the integration at each breakpoint, which
   fixed the accuracy **and** cut the runtime of the whole step-2 matrix by more than an order of
   magnitude.
2. The residual failure was in **`k3` itself**: the last quadrature node of an interval ending at
   a dose time was evaluating the enzyme *after* the dose, because `enzyme(t)` is
   right-continuous. Taking the left limit there brought the two routes to `1e-10`. The check was
   wrong, and the fix is recorded in the file rather than silently applied.

Neither error changed a reported conclusion, because both were found before any table was
written. They are recorded because the point of the re-check is that it catches things.

---

---

# The follow-up: admissibility of the trough reading, and the dose sensitivity table

`audit/k4_admissibility.py`, re-checked in `audit/k6_recheck.py`.

## Does anything constrain when "after the third dose" was sampled?

**The trial paper: no, and this pass checked exhaustively.** Extended Data Table 3's caption
carries no time point. The Methods define the endpoint only as "changes in blood 2-MC and
3-HP levels from baseline (predose levels) to post-dose levels". There is **no biomarker
time-course figure**: Fig. 3 is PK exposure to PCCA mRNA and Fig. 4 is MDE profiles, so 3-HP
appears in the article exactly once as a number. The only sampling time named anywhere is a
**predose** draw, in the antibody paragraph, "the predose sample for dose 3 (52 days after
the first dose)". That is a trough, it is about anti-PEG titres rather than 3-HP, and it is
the only evidence there is. (In passing: cohort 1 is Q3W, so a third dose on schedule falls
at day 42, not 52. Nominal and actual timing are not the same thing in this trial.)

**The registry (NCT04159103): not checked, and this is a real open item rather than a null
result.** Outbound network access is blocked by this environment's egress policy. Both
`curl https://clinicaltrials.gov/api/v2/studies/NCT04159103` and a `WebFetch` of the same
record return a 403 CONNECT denial, logged by the proxy itself as *"gateway answered 403 to
CONNECT (policy denial or upstream failure)"*. The registry's secondary outcome measures
carry a `timeFrame` field, and if it names a visit for the 3-HP endpoint it settles the
question outright. **It needs a session with egress, or the record pasted in.**

**Baek's Source Data: not available**, for the same reason. Fig. 2's caption says "Source
data are provided as a Source Data file"; that file was not among the five PDFs.

**What Baek's Fig. 2 caption does say**, read in this pass:

> "Blood was collected from mice prior to treatment and on Days 2, 8, 14, 22, and 28 after
> the first IV dose, and 2 days after the second IV dose"

on a regimen of "two IV bolus doses ... on Days 0 and 28". So the preclinical PD data the
translational model was fitted to is a **full time course across a 28-day cycle**, not a
trough-only design: Day 2 is near that cycle's PD peak and Day 28 is its trough. The
programme's own preclinical convention was to sample the whole cycle. It does **not** tell us
which point of the human cycle Extended Data Table 3 reports, because the trial's biomarker
sampling schedule is not published.

**One thing the same figure does settle, for finding 1(e).** Fig. 2e's y-axis is
`3-HP (uM)`, plotted 0 to 200, in the **same** `Pcca-/- (A138T)` mouse model Attarwala
fitted. Same company, same strain, same analyte, plotted in micromolar. Attarwala's Table 1
prints the 3-HP baseline as `60.1` "mM/L". That reading is now confirmed as `umol/L` by an
**external** source and not only by the internal LOQ argument.

## The dose sensitivity table

`D_AMP` is **not moved**. Every row rebinds `design.run.D_AMP` and `LAMBDA` inside a context
manager that restores them. `n` is monthly visits to `sd(ln theta) <= 0.20` at load CV 10%,
median over 8 seeds, searched to `n = 60`.

| half-life | row | `D` | `Δ_peak/θ` | fall at peak | `m` at peak | Δ/θ range, n=12 | `sd@12` | `n` |
|---|---|---|---|---|---|---|---|---|
| 3.5 d | **committed `D = 6000`** | 6000.0 | 2.133 | 66.36% | 3.364 | 1.817 | 0.1888 | **11** |
| 3.5 d | 20% fall at the peak | 766.0 | 0.272 | 20.00% | 8.000 | 0.232 | 0.3682 | **>60** |
| 3.5 d | 30% fall at the peak | 1311.1 | 0.466 | 30.01% | 6.999 | 0.397 | 0.2744 | **42** |
| 3.5 d | 20% fall at the trough | 12255.3 | 4.357 | 80.15% | 1.985 | 3.711 | 0.1835 | **10** |
| 3.5 d | 30% fall at the trough | 20976.9 | 7.458 | 87.37% | 1.263 | 6.352 | 0.1879 | **12** |
| 7 d | committed `D = 6000` | 6000.0 | 2.667 | 71.16% | 2.884 | 1.782 | 0.2657 | **30** |
| 7 d | 20% fall at the peak | 622.3 | 0.277 | 20.25% | 7.975 | 0.185 | 0.5154 | **>60** |
| 7 d | 30% fall at the peak | 1065.2 | 0.473 | 30.33% | 6.967 | 0.316 | 0.3798 | **>60** |
| 7 d | 20% fall at the trough | 2489.4 | 1.106 | 50.51% | 4.949 | 0.739 | 0.2814 | **37** |
| 7 d | 30% fall at the trough | 4260.9 | 1.894 | 63.64% | 3.636 | 1.265 | 0.2639 | **30** |

And admissibility on every row. `capacity` is `h1_admissibility.py`'s constraint 5 at both
readings of the residual-capacity fraction; it binds at `Δ_peak/θ = 9.0` (10% reading) and
`19.0` (5%). `B@1` and `A@1` are the two published inhibition shapes' ceilings at that row's
own enzyme ratio, with `Imax` set to 1 rather than the fixed 0.999.

| half-life | row | used@10% | used@5% | capacity | fall | `B@1` | `A@1` | bracketed |
|---|---|---|---|---|---|---|---|---|
| 3.5 d | committed | 0.3133 | 0.1567 | OK | 66.36% | 64.97% | 68.09% | yes |
| 3.5 d | 20% peak | 0.1272 | 0.0636 | OK | 20.00% | 19.91% | 21.40% | yes |
| 3.5 d | 30% peak | 0.1466 | 0.0733 | OK | 30.01% | 29.72% | 31.79% | yes |
| 3.5 d | **20% trough** | 0.5357 | 0.2679 | **OK** | 80.15% | 78.43% | 81.33% | yes |
| 3.5 d | **30% trough** | 0.8458 | 0.4229 | **OK** | 87.37% | 85.63% | 88.18% | yes |
| 7 d | committed | 0.3667 | 0.1833 | OK | 71.16% | 69.63% | 72.73% | yes |
| 7 d | 20% peak | 0.1277 | 0.0638 | OK | 20.25% | 20.16% | 21.67% | yes |
| 7 d | 30% peak | 0.1473 | 0.0737 | OK | 30.33% | 29.82% | 31.68% | yes |
| 7 d | **20% trough** | 0.2106 | 0.1053 | **OK** | 50.51% | 49.63% | 52.53% | yes |
| 7 d | **30% trough** | 0.2894 | 0.1447 | **OK** | 63.64% | 61.50% | 65.02% | yes |

## Are the trough rows physically admissible? Yes, on both tests

**Saying so plainly, because the first version of this check said the opposite and was
wrong.** Neither the capacity ceiling nor the published shape excludes any trough row.

**(a) The capacity ceiling: every trough row passes, on both references.** The four sit at
`Δ_peak/θ` of `4.357`, `7.458`, `1.106` and `1.894` against a bound of `9.0` at the 10%
reading. The worst row, a 30% trough fall at the committed half-life, uses `84.58%` of normal
capacity. That is inside the bound, and it is also a strong clinical claim: one fortnightly
dose lifting a patient from 10% of normal PCC capacity to 85% of it at every peak. The
trial's own outcomes do not support that reading, since the relative risk of a metabolic
decompensation event fell by 70% rather than to near-normal, and the biomarkers moved by
`-29%` to `+5%`.

**(b) The shape ceiling: it does not exclude them either, and this is the correction.** At
`Imax` fixed to its published `0.999` the ceilings sit a fraction of a percentage point
*below* the paper's fall on most rows, which reads as an exclusion. But Attarwala **fixed**
`Imax` at 0.999 rather than estimating it (Table 1, footnote m), so a verdict that turns on
it is a verdict about a modelling convenience. Rerunning each ceiling at `Imax = 1` separates
the two questions:

- **`PUB_A`'s ceiling at `Imax = 1` is exactly `Δ/(θ+Δ)`**, the `IC50 → 0` limit of the
  inhibition ratio, verified symbolically in `k6` check 8. It **clears** the paper's fall on
  every row by 1.2 to 2.0 percentage points. `PUB_A` is not excluded at any dose. The 0.999
  cap costs it 1.0 to 2.4 points, and that is the entire reason it appeared to fail.
- **`PUB_B`'s ceiling falls short on every tabulated row at both values of `Imax`**, by 0.09
  to 1.74 points. That shortfall is structural: it is the published floor
  `E0/base = 0.0987/60.1` riding along and blunting every fall, not the `Imax` cap.

**The result that matters is the bracket.** On every row, at `Imax = 1`, the paper's fall
lies strictly between the two published shapes' ceilings: `PUB_B` a little too shallow,
`PUB_A` with headroom. The paper's dose-response is bracketed by the published structures
rather than outside them, and a shortfall of under two percentage points on one side cannot
exclude a dose.

**One scope limit on that bracket**, found by `k6` check 9 and not by `k4`. `PUB_A`'s ceiling
is above the paper's fall at every dose from `Δ_peak = 100` to `40000`. `PUB_B`'s **crosses**
it at `Δ_peak = 443.6`, that is `Δ_peak/θ = 0.148`: below that dose `PUB_B`'s shape can do
*more* than the paper's marker, not less, and the bracket reverses. Every row in the table
sits above the crossover, the smallest at `Δ_peak = 816.8`, so the bracket holds throughout
the table. It is not a property of the three forms at every dose, and `k4`'s first draft said
it was.

**What is actually against the trough reading**, then, is two things, neither of them
admissibility:

1. At the committed half-life it needs `D` to be **2.0x and 3.5x** the committed value
   (12255 and 20977 against 6000), producing peak marker falls of 80% and 87%, while the
   trial reports falls two to three times *smaller* than the committed calibration already
   is. That is coherent only if the trial sampled at the very bottom of a very deep
   excursion.
2. The clinical claim in (a).

**And the 7 d half-life mostly dissolves objection 1.** A longer half-life makes the trough a
much larger fraction of the peak (25.0% at 7 d against 6.25% at 3.5 d), so the same trough
target needs far less dose: `D = 2489` and `4261`, peak falls of 51% and 64%, capacity at 21%
and 29% of normal, both bracketed. **So 7 d plus the trough reading is the most internally
coherent combination on the table**, and it is also the one that costs the design the most,
30 to 37 monthly visits against the committed 11. Both of its ingredients are the reading the
primary sources lean towards, and neither is established.

---

# Step 3. Is the enzyme capacity identifiable under each structure?

`audit/k5_identifiability.py`, re-checked in `audit/k6_recheck.py`. Three routes that share
no code: exact reparametrisation in sympy, the numerical rank of the parameter-to-observation
Jacobian, and constructing the unidentifiable family explicitly and walking it.
`StructuralIdentifiability.jl` is still absent from this container and is not needed: every
model here reduces to a rational map of low degree, and the reduction is **verified** rather
than asserted.

## 3a. Each structure has an exact reparametrisation, and each carries three numbers

**`PUB_B`** (Attarwala, algebraic). With `q = 1 - Imax` and `E = θ + Δ`, sympy returns the
marker as a ratio of two polynomials of degree 1 in `Δ`, that is a **Möbius function**:

```
m(Δ) = (P Δ + U)/(Δ + V) ,   P = E0 + base·q
                             U = P·θ + (E0 + base)·IC50
                             V = IC50 + θ
```

residual after substitution exactly `0`. Re-derived in `k6` check 1 by **exact rational
interpolation** with no floating point anywhere: three points determine `(P, U, V)` as exact
fractions, and four more agree exactly.

**`PUB_A`** (turnover). The forcing term is the same shape, `f(θ+Δ) = (qΔ + W)/(Δ + V)` with
`W = qθ + IC50` and `V = θ + IC50`. There is **no free multiplicative constant** in front of
it: the equation is `dm/dt = (L+b)·f - k_out·m` and `(L+b)` is the recorded load, so `f` is
identified as a function rather than up to a scale. Inverting, with `Imax` fixed,

```
IC50 = (W - qV)/Imax ,   θ = V - IC50 = (V - W)/Imax
```

both recovered symbolically. `k6` check 3 reaches the same closed form by a Gröbner
elimination of `IC50`, which leaves the generator `V - W + qθ - θ`.

**`PAPER`**. Re-derived from scratch by the same method so all three are counted one way:
`Δ(m) = A/m + B + Cm` with `A = K_m(L+b)`, `B = L+b-θ-k_e K_m`, `C = -k_e`. Three invariants,
four parameters. This is P5's result and it is unchanged.

## 3b. The numerical rank of the sensitivity map agrees exactly

24 phase-spread visits. `PAPER` and `PUB_B` are differentiated by the **complex step**, exact
to machine precision; `PUB_A` by integrating its **forward sensitivity system**. Columns are
scaled before the SVD.

| load | structure | params | rank | deficit | smallest singular value |
|---|---|---|---|---|---|
| constant | `PAPER` | 4 | **3** | 1 | 3.17e-16 |
| constant | `PUB_A` | 4 | **4** | 0 | 1.58e-04 |
| constant | `PUB_B` | 5 | **3** | 2 | 5.58e-17, 4.38e-17 |
| varying | `PAPER` | 4 | **4** | 0 | 1.02e-03 |
| varying | `PUB_A` | 4 | **4** | 0 | 7.96e-03 |
| varying | `PUB_B` | 5 | **4** | 1 | 3.15e-17 |

And the null directions read off the algebra exactly:

- `PAPER`, constant load: `θ=+0.9250, K_m=-0.1124, k_e=+0.0000, b=+0.3629`. **`k_e` is
  exactly zero**, which is P5's fibre ("`k_e` unchanged") recovered by a completely different
  route. Varying the load: full rank, the fibre closes.
- `PUB_B`, varying load: `E0=-0.0010, base=+0.9784, IC50=-0.1461, theta=+0.1461, b=-0.0000`.
  **`IC50` and `θ` exactly equal and opposite**, which is `V = IC50 + θ` held fixed, and
  **`b` exactly zero**, which says load variation bought `b` and nothing else.

`k6` check 4 confirms every rank by **Richardson-extrapolated central differences** instead
of the complex step, with spectral gaps of 5.6e4 to 1.9e10 at the predicted rank.

**Two bugs this leg found, both mine, both fixed.** The first version used central differences
at `h = 1e-6` and reported `PAPER` as full rank at a constant load, contradicting P5's exact
fibre; the fibre is exact and the difference quotient was not. The second used
`k2_structures.PubA.marker`, which reads `θ` from a module-level function and therefore has
**no `θ` dependence at all**, so the `θ` column was identically zero and the reported rank
deficit was an artefact of the harness rather than of the model.

## 3c. `PUB_B`'s fibre, constructed and walked

Parametrising by `IC50 = R` and solving the rest in closed form,
`θ = V - R`, `base = (U - PV)/(R·Imax)`, `E0 = P - q·base`:

| `IC50 = R` | `θ = V - R` | `base` | `E0` | max relative change in the record |
|---|---|---|---|---|
| 30.6 | 3033.1 | 888.1 | 0.2396 | 1.44e-14 |
| 306.4 | 2757.4 | 88.81 | 1.0389 | 7.60e-15 |
| 1531.9 | 1531.9 | 17.76 | 1.1100 | 7.50e-15 |
| 63.75 (truth) | 3000.0 | 426.81 | 0.7009 | 1.08e-14 |
| 2757.4 | 306.4 | 9.868 | 1.1179 | 7.60e-15 |
| 3060.7 | 3.064 | 8.890 | 1.1188 | 7.60e-15 |

**`θ` runs from 3.06 to 3033, a factor of 990, with the record unchanged to 1.4e-14.** Run
again at a 10% varying load, the figures are the same to a digit: load variation does **not**
close this fibre, because `L` enters only as the common factor `(L+b)/L0`.

`k6` check 2 reaches the same conclusion by Gröbner elimination: eliminating
`(E0, base, IC50)` from the three invariant equations leaves **no generator** involving only
`θ`. An empty elimination ideal is the definition of `θ` being unconstrained.

## 3d. `PUB_A`'s profile in `θ` is not flat

The counterpart. Moving `θ`, sliding `IC50 = V - θ` to hold the one invariant a null direction
would have to preserve, and **re-optimising `k_out`** by a golden-section search at each trial
`θ`:

| `θ` | `IC50 = V - θ` | best `k_out` | mismatch |
|---|---|---|---|
| 750 | 2273 | 18.98 | 3.196e-01 |
| 2700 | 322.9 | 2.887 | 2.476e-01 |
| 2970 | 52.93 | 0.5681 | 7.994e-02 |
| **3000 (truth)** | **22.93** | **0.2768** | **2.893e-11** |
| 3030 and above | negative | | inadmissible |

A gap of `2.8e9` between the truth and the nearest admissible alternative. `θ` is decisively
identifiable under `PUB_A`, and its admissible range is `(0, V)` with the record picking one
point out of it.

---

# Step 4. The general statement

What the three reductions have in common is not the algebra, which differs, but the **count**.

## The three-number bound

Fix the recorded load. The entire content of a record is then the **function** that sends the
therapeutic increment to the marker level. In each structure that function belongs to a
family of dimension exactly **three**:

| structure | the family | three numbers |
|---|---|---|
| `PAPER` | `Δ(m) = A/m + B + Cm` | `A = K_m(L+b)`, `B = L+b-θ-k_e K_m`, `C = -k_e` |
| `PUB_B` | `m(Δ) = (PΔ+U)/(Δ+V)` | `P = E0+base·q`, `U = Pθ+(E0+base)IC50`, `V = IC50+θ` |
| `PUB_A` | `f(Δ) = (qΔ+W)/(Δ+V)` | `q` fixed, `W = qθ+IC50`, `V = θ+IC50` |

In every case the map is a rational function of degree (1,1) in one of the two variables, and
such a function carries three numbers and no more. **That is a property of the form and not
of the biology, which is why the obstruction survives a change of structure: it is the same
count each time.**

## What differs is what each structure has to pay for out of those three

| structure | free parameters (`Imax` fixed) | what a record carries | deficit |
|---|---|---|---|
| `PAPER` | 4: `θ, K_m, k_e, b` | 3 at constant load | **1** |
| `PAPER` | 4 | 4 once the load varies (`b` enters the invariants) | **0** |
| `PUB_A` | 4: `θ, IC50, k_out, b` | 2 from `(W,V)` + `k_out` from the transient + `b` | **0** |
| `PUB_B` | 5: `θ, IC50, E0, base, b` | 3 from `(P,U,V)` + `b` from the load | **1** |

## Three consequences

**1. The paper's fibre is generic; its cure is specific.** Every structure here is one number
short at a constant load. What differs is where the fourth comes from. For the paper it is
the load, because `b` enters the invariants. That is P5, and it is why a diet record smoothed
to a weekly mean destroys the identification.

**2. A marker state buys a number; a marker scale spends one.** `PUB_A` is identifiable and
`PUB_B` is not, and they differ in exactly two places. `PUB_A`'s marker is a **state**, so the
transient pins `k_out` separately from the steady map: one number bought. `PUB_B`'s marker
carries **two free scale constants**, `E0` and `base`, which enter the invariants only
through `P` and `E0+base`: one number spent. Net, `PUB_A` is exactly determined and `PUB_B`
is one short. That is the thing to look for in any third structure.

**3. Where `θ` sits in the deficit is not a detail.** Under `PUB_B`, `θ` appears only inside
`V = IC50 + θ` and inside `U`, and 3c walks it over its whole positive range at a
bit-identical record. So under `PUB_B`, `θ` is not merely imprecise, it is **not
identifiable**, at any design, at any number of visits, with or without load variation.

**That is the exact explanation of what step 2 measured.** Fitting the paper's model to
`PUB_B` data returned `θ̂` biased by 1.3x at the committed dose and up to 20x at the trial's
own, with a residual under a fifth of one `σ`. The estimator was not failing. It was
reporting a point on a fibre along which the data are silent, and the paper's parametrisation
chose which point.

**The statement for the paper, in one sentence:** `θ` is identified relative to a structure,
not absolutely, and the structure has to supply a fourth number from somewhere (a varying
load, a marker transient, or an external anchor) before any record can determine it.

---

# The independent re-check for the follow-up and steps 3 to 4

`audit/k6_recheck.py`. **9 of 9 pass.** Three of them failed first, and all three failures
were real.

| check | independent route | result |
|---|---|---|
| 1. `PUB_B` is Möbius | exact rational interpolation, no floating point | `(P,U,V)` exact, matches sympy |
| 2. `PUB_B`'s fibre | Gröbner elimination of `(E0, base, IC50)` | no generator constrains `θ` |
| 3. `PUB_A` pins `θ` | the same elimination, `IC50` only | `θ = (V-W)/Imax` in closed form |
| 4. the ranks | Richardson extrapolation instead of the complex step | spectral gaps 5.6e4 to 1.9e10 |
| 5. `PUB_A` sensitivities | finite differences at each parameter's own scale | agree to 1.5e-9 |
| 6. two `PUB_A` integrators | `k2`'s simulator against `k5`'s sensitivity solver | agree to 2.7e-10 |
| 7. the capacity cancellation | symbolic, `used = frac·(1 + Δ_peak/θ)` | `θ` cancels; binds at 9.0 and 19.0 |
| 8. the `Imax=1` ceiling | symbolic `IC50 → 0` limit | exactly `Δ/(θ+Δ)`, matches the grid |
| 9. the bracket | swept 61 doses instead of tabulating 10 | holds above `Δ_peak = 443.6`, reverses below |

**The three failures.**

1. **Check 5 found a real error in `k5`.** The forward sensitivity system's **initial
   conditions** set `∂m0/∂IC50` and `∂m0/∂θ` to zero. They are not zero: `m0` depends on both
   through `f(E(t0))`. The omission left a transient worth `1e-3` relative, which is exactly
   what the check measured. Now carried; both columns agree to `1e-9`. The rank conclusion
   was unaffected, since a `1e-3` error in two columns cannot create rank, but the
   sensitivity values were wrong.
2. **Check 9 refuted a claim in `k4`.** `k4`'s first draft called the bracket "a property of
   the three forms". Sweeping 61 doses instead of tabulating 10 shows `PUB_B`'s ceiling
   crosses the paper's fall at `Δ_peak = 443.6`, below which the bracket reverses. `k4` now
   carries the crossover and the scope limit.
3. **Check 5 also failed a second time on its own step size**, and that one was the check's
   fault, not the code's. With the integrator at `rtol = 1e-12` and `b` entering on the load's
   scale of `322.7`, a step of `1e-5` made the difference quotient integrator noise. The step
   is now taken relative to each parameter's own scale and Richardson-extrapolated. Recorded
   rather than silently fixed, because a re-check that quietly retunes until it passes is not
   a re-check.


# Proven, measured, judgement

| claim | status |
|---|---|
| The slope is a linear coefficient and not a Hill exponent; `Imax` is fixed at `0.999` | **proven** from Suppl Table S1(B) and Table 1 footnote m |
| `"mg/g"` and `"mM/L"` are misprints for `ug/g` and `umol/L` | **proven** from the paper's own Fig. 2B axis and its own LOQ |
| The Table 1 PD estimates are mouse and were not scaled to human | **proven** from the paper's own scaling sentence and Baek's |
| Five repository places present them as human, or cite a non-existent extraction route | **proven**, located by line |
| `3.85 d` is one rate constant, not the model's protein half-life; the `~7 d` is consistent with the constants | **measured**, `k1_sources.py`, re-checked exactly |
| The trial's 3-HP fall after dose 3 is `-29` to `+5%` across five cohorts | **proven**, Extended Data Table 3 |
| `D_AMP = 6000` implies a `62-66%` fall, two to three times the reported figures | **measured** |
| The design numbers under a 7 d half-life and under a 20-30% fall | **measured**, table above, 8-seed spread quoted |
| Model A has no published parameters; Baek is already registered as `mrna_moa` | **proven** from Baek's Methods and Table 1, and from `refs.sources:200-203` |
| ~~The published inhibition shape cannot reach the repository's calibrated fall~~ | **RETRACTED**, `REVIEW_LOG.md` R15. It was a knife edge on `Imax = 0.999`, a constant the source FIXED rather than estimated. At `Imax = 1` the turnover shape clears the fall at every dose |
| The repository's fall is BRACKETED by the two published shapes' ceilings on every tabulated dose | **measured**, and swept independently over 61 doses |
| That bracket reverses below `Delta_peak = 443.6`, so it is not a property of the forms at every dose | **measured**, `k6` check 9, which refuted `k4`'s own wording |
| `PUB_B` misspecification is invisible in the residuals and biases `theta` by `1.29` to `2.03` at the `60%` calibration, and by `4.3` to `20.5` at the trial's own `20-30%` | **measured** |
| `PUB_A` misspecification is loud on the recommended schedules and silent on the dose-period one | **measured** |
| The schedule ranking survives both structures; the parameter estimate does not | **measured** |
| The exactly singular case survives both structures on every branch | **measured**, which is the first direct test of `A1_supplementary.tex:51-52` |
| Where the load enters `PUB_A` and `PUB_B` | **judgement**, stated in the file header and in every table |
| Which calibration branch is "the" one | **not judged**. Both are carried everywhere |
| That the trial's samples were troughs | **not judged**. The paper does not say, and both readings are carried, and they disagree about whether the design is feasible |
| Every trough row passes the capacity ceiling at both references, the worst at `84.58%` of normal | **measured**, `h1`'s constraint 5 recomputed per row |
| The trial paper names no sampling time for Extended Data Table 3, and has no biomarker time-course figure | **proven**, checked exhaustively in this pass |
| The registry's `timeFrame` field would settle it | **not checked**. Outbound egress is blocked here; both `curl` and `WebFetch` get a 403 CONNECT denial |
| Baek's mouse PD was sampled at Days 2, 8, 14, 22, 28 of a 28 d cycle, a full time course | **proven**, Fig. 2 caption |
| Baek Fig. 2e plots 3-HP in `uM` in the same mouse strain, confirming Attarwala's `"mM/L"` is `umol/L` | **proven**, and this is an EXTERNAL confirmation of what was an internal argument |
| `PUB_B` is exactly Moebius in `Delta`, so a record carries three numbers | **proven** symbolically, and by exact rational interpolation |
| Under `PUB_B`, `theta` is NOT IDENTIFIABLE at any design, with or without load variation | **proven**. Groebner elimination leaves no relation on `theta`; the fibre was walked over a `990`-fold range at a record unchanged to `1.4e-14` |
| Under `PUB_A`, `theta` IS identifiable, at a constant load, with no external information | **proven**, `theta = (V-W)/Imax` in closed form; profile gap `2.8e9` |
| The rank counts (`PAPER` 3 of 4 then 4 of 4; `PUB_A` 4 of 4; `PUB_B` 3 of 5 then 4 of 5) | **measured** by complex step, re-checked by Richardson |
| `PAPER`'s null direction has `k_e` exactly zero, which is P5's fibre by a new route | **measured** |
| `PUB_B`'s null direction has `IC50` and `theta` exactly opposite and `b` exactly zero | **measured** |
| The three-number bound, and that a marker state buys a number while a marker scale spends one | **proven** for these three structures. Whether it generalises beyond degree-(1,1) maps is **judgement** |
| That `PUB_B`'s non-identifiability EXPLAINS step 2's bias | **judgement**, though a well-supported one: the mechanism predicts both the size of the bias and the absence of a residual signal |

## What I did not do

- **No paper edit, no constant change, no merge.** Every correction above is proposed text.
  `PHASE4_PROPOSED_TEXT.md` is still unapplied.
- **Steps 3 and 4 are now done** (they were not, when the section above was first
  written; the entry is kept so the order of work is legible).
- **No test was added to `tests/test_models.py`.** The natural candidates are the two structure
  ceilings and the `PUB_A`/`PUB_B` calibration identities, since both are exact arithmetic on
  published constants. They belong there if and when the simulators become part of the paper
  rather than of the audit, and that is a decision for after you have read this.
- `StructuralIdentifiability.jl` is still not installed in this container. It turned out not
  to matter: every structure here reduces to a rational map of degree (1,1), the reduction
  is verified three ways, and a differential-algebra package would answer the same
  question more slowly. It would still be worth running as a fourth opinion.

## The questions, and where they now stand

1. **The sampling time: carry both readings.** Decided. It is carried throughout, and the
   two readings still disagree about feasibility: `10` to `12` monthly visits under the
   trough reading at the committed half-life, against `42` or unreachable under the peak
   reading. **The one thing that would settle it and has not been checked is the registry's
   `timeFrame` field for the 3-HP secondary outcome**, which this container cannot reach.
   That is the single highest-value open item in this report and it costs one fetch.
2. **`D_AMP` does not move.** Decided. The sensitivity table above is what stands in its
   place, with admissibility on every row. One of the two arguments for moving it is
   withdrawn (`REVIEW_LOG.md` R15); the argument from the trial's reported falls stands.
3. **The `theta` bias goes in the paper.** Decided. Drafted, unapplied, as passage 9 of
   `audit/PHASE4_PROPOSED_TEXT.md`, together with a pointer in the discussion's "What would
   change the answer" list. Steps 3 and 4 turned it from a measurement into a proof: under
   that structure `theta` is not identifiable at all, so the bias is not estimator failure
   but a parametrisation choosing a point on a silent fibre.
4. **The five mouse-versus-human corrections** are drafted, unapplied, as passage 10.

## What remains open

- **The registry `timeFrame`.** See 1.
- **Baek's Source Data file**, for the same reason. It would give the mouse 3-HP time course
  as numbers rather than as a figure, which would let the translational model be simulated
  end to end rather than through our own calibration.
- **Whether `lambda` moves to `ln 2/7`.** Not proposed here. Passage 10b states the
  qualification and keeps the constant. Moving it costs `11 -> 30` monthly visits, and every
  Phase 3 figure is computed at `ln 2/3.5`.
- **Whether the `PUB_A`/`PUB_B` simulators become part of `src/` and get tests.** The natural
  test candidates are the two calibration identities and the closed form
  `theta = (V-W)/Imax`, all exact arithmetic. Not added, because they are audit code until
  you say otherwise.
