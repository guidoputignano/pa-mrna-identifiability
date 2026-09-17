# The September 2026 review pack, checked

Four reviews arrived (3-HP response timescales; prior art on periodic dosing, sampling design
and diet-input models; non-dietary propionate sources; dietary intake variability; breath-test
repeatability). The 3-HP review is present **twice, verbatim**, so the pack is five documents
and four reviews.

Three of the five say their retrieval tools failed and that they are written from memory with
confidence labels. One (breath tests) plainly had working retrieval: it carries verbatim
quotes and section pointers.

**What changed since `OPEN_QUESTIONS.md` was written: network egress works in this container.**
`doi.org`, `api.crossref.org` and NCBI E-utilities all answer. The verification rule in
`refs/literature_pending.md` can now be executed rather than deferred. Nine identifiers were
resolved for this note. None of the rest was touched, and nothing has been moved into
`paper/refs.sources`.

---

## 1. Verified against the registrar

| identifier | resolves to | verdict |
|---|---|---|
| PMID 2233273 | Thompson GN et al., *Metabolism* 1990;39(11):1133-7 | **Confirms this repository, contradicts the review.** See 2.1 |
| PMID 481943 | Hsia YE, Scully KJ, Rosenberg LE, *Pediatr Res* 1979;13(6):746-51 | Confirms "<5%", with a caveat. See 3.2 |
| PMID 29021961 | Daly A et al., *Mol Genet Metab Rep* 2017;13:83-89 | Confirms the feared reading. See 3.1 |
| PMID 1909779 | Barshop BA et al., *Pediatr Res* 1991;30(1):15-22 | **Contradicts a sentence in `08_specification.tex`.** See 3.3 |
| `10.1088/1752-7155/5/4/046001` | Wagner D et al., *J Breath Res* 2011;5(4):046001 | Real. The closest published thing to Tier 1. See 3.3 |
| `10.3390/nu15234903` | Pinto A et al., *Nutrients* 2023;15(23):4903 | Real. First author is Pinto, not Evans. See 3.1 |
| `10.1007/s10928-009-9114-z` | **Nyberg J, Karlsson M, Hooker A**, *JPKPD* 2009;36(2):125-145 | This repository's attribution is right, the review's is wrong. See 2.2 |
| `10.1038/s41467-020-19156-3` | Jiang L et al., *Nat Commun* 2020;11 | Real |
| `10.1038/s41467-024-47460-9` | Baek R et al., *Nat Commun* 2024;15 | Real, and already at `refs.sources:201` |

---

## 2. Errors in the pack

### 2.1 The odd-chain fraction. The review is wrong and this repository is right

The non-dietary review puts odd-chain fatty acids at "minor (single-digit %)" and gut bacteria
at "~20-25%". Thompson 1990's abstract, fetched:

> Amino acid catabolism accounted for a mean of 51.7% (range, 24.5% to 66.4%) of total
> propionate production. The mean decrease in propionate production after metronidazole was
> 22.2% +/- 8.5 (P less than .02); this percentage is likely to represent the **minimum**
> propionate production attributable to gut bacteria. **Approximately 30% of total propionate
> production was unaccounted for, and is likely to arise primarily from odd-chain fatty acid
> catabolism in the fasting state.**

So `refs/literature_pending.md`'s "gut bacteria >=22%, odd-chain fatty acids ~30%" is exactly
right, including the `>=`, and the review's single-digit figure is off by a factor of about
five. Two further points the review drops: the study is **six patients in the fasting state**,
and the 30% is an **unaccounted residual attributed by elimination**, not a measurement. The
review also cites Bain 1988's metronidazole effect (urinary metabolites, ~40-60%) alongside
Thompson's (whole-body production, 22.2%) without noting these are different quantities.

**This row needed no change and still needs none.** It is recorded because the review would
have moved it in the wrong direction.

### 2.2 A DOI that resolves to the wrong paper, and two misattributions

- `10.1016/j.celrep.2017.11.049`, given by the prior-art review as An et al., "Systemic
  Messenger RNA Therapy as a Treatment for Methylmalonic Acidemia", **resolves to a barrel-cortex
  neuroscience paper** (Crandall et al., *Cell Rep* 21:3065-3078). The correct identifier is
  `10.1016/j.celrep.2017.11.081`, *Cell Rep* 2017;21(12):3548-3558, PMID 29262333. Three digits
  wrong, and it resolves to something real, which is the worst case for a hand check.
- The clinical review's "An D et al., *Cell Reports Medicine* 2019" does not exist as described.
  The 2019 follow-up is *EBioMedicine* 2019;45:519-528 (PMID 31303505) with a different author
  line; the Schneller/Frassetto line belongs to the 2017 *Cell Reports* paper.
- `10.1007/s10928-009-9114-z` is **Nyberg, Karlsson and Hooker**, as `literature_pending.md`
  already has it. The prior-art review attributes it to "Duffull SB, et al." The registrar
  settles it.

### 2.3 The breath-test CV, again

The clinical review repeats "CV 12.8% in MMA patients, 6.9% in controls". `REVIEW_LOG.md` **R4**
already retracted that: Manoli's Results give `14.02 +/- 12.86%`, and `12.86` is the standard
deviation quoted in the abstract as though it were the mean. The review has inherited the
source's own error a second time. **Do not import it.** Note that the breath-test review, the
one with working retrieval, independently found the same discrepancy and flagged it, which is
corroboration of R4 from outside this repository.

### 2.4 A timescale inference that does not follow

The 3-HP review concludes that capacity-driven changes equilibrate over days, citing
post-transplant decay curves. That does not distinguish a **slow marker state** from a **fast
marker state tracking a capacity that itself arrives over days**. Graft function ramps over the
first post-operative week (the review's own Di Campli 2003 entry shows aminopyrine recovery
climbing to day 5), so the transplant curve is consistent with either. The review's dialysis
rebound observation is the one that does bear on the question, and it points at slow peripheral
pools rather than a slow plasma compartment. See 3.4.

---

## 3. What actually changes

### 3.1 The load recording error, and the one number the specification does not sweep. **LARGEST**

Daly 2017's abstract, fetched, says patients were **"prescribed"** natural protein intakes; it
is a cross-sectional survey of centre practice, not a record of intake as eaten. The dietary
review's stronger claim, that **no PA study reports repeated as-eaten daily intake paired with
biomarkers**, is consistent with everything found here. So the feared reading in
`OPEN_QUESTIONS.md` section 2 holds: a PA dietary record, as the literature has it, is a
prescription.

This turns out **not** to break the design the way that row anticipated, and the reason is worth
recording because it is the opposite of what was expected. `audit/m1_loadbias_prior.py`, two
independent code paths agreeing to `1e-14`, 24 seeds:

| load-bias prior | record varies at 10% CV | record **constant** |
|---|---|---|
| 5% (as published) | `0.1660` | `0.1659` |
| 40% | `0.1709` | `0.1911` |
| flat | `0.1717` | **`0.2659`** |

At a constant recorded load with the prior removed entirely the answer is `0.2659`, not
infinity. **The dose train breaks Proposition rank's singularity**: the increment varies, so the
columns are not constant even when the load is. `08_specification.tex` currently says the null
direction "is filled by the 5% prior on the load bias, and only real day-to-day variation lets
the data fill it instead", and that overstates it for `theta`: without the prior and without any
variation the data still give `0.266`, which is a degraded result, not a prior-supported one.
The sentence should be scoped. It is a softening, not a retraction.

**The real exposure is elsewhere, and it is not currently swept.** `spec_numbers()` computes
`sigma = LOC.sigma_total(LOAD_FRAC, ASSAY_CV)` **once, at 5%,** and reuses it on every row. The
table's "load record kept to 1% or to 20%" rows vary the day-to-day **variation** in the record,
not the recording **error**. The error enters only through `sigma`:

| recording error | `sigma` | `sd(ln theta)` median, 24 seeds |
|---|---|---|
| 1% | `15.24` | `0.1566` |
| 5% | `21.96` | `0.1660` (the published headline) |
| 10% | `35.54` | `0.1916` |
| 20% | `66.24` | **`0.2637`, fails the target** |
| 40% | `129.95` | **`0.4023`, twice the target** |

So the design meets `0.20` only if the substrate load is recorded to better than roughly 11 to
12%. A reader will take "the load record gives `0.170` whether it is kept to `1%` or to `20%`"
as robustness to the recording error, and it is not that.

And the review supplies the number to put against it. *Nutrients* 2023;15(23):4903 (verified;
Pinto first author, not Evans) is the only study of its kind: 16 children with PKU, three
one-day diaries per month for six months, weekly blood spots. Prescribed against actual natural
protein was a **median +40%**, actual running at 167% of prescribed. If a PA record is a
prescription and sits a comparable distance from what is eaten, the recording error is roughly
40%, not 5%, and the headline becomes `0.40`.

**This is the change to make.** Add the `sigma` sweep to the specification, state the
requirement as a recording error below about 12%, and say where the 40% analogue comes from.
Nothing else in the paper moves. No constant is touched.

### 3.2 Open question (b), 5% or 10%, can be closed but not on this evidence alone

Hsia 1979 verified: "Each of seven lines from patients with propionic acidemia had clearly
detectable PCC activity which was **less than 5%** of that in control lines." It favours 5%, as
`literature_pending.md` says. The caveat the review does not state: these are **cultured
fibroblasts**, seven lines, partially purified enzyme. That is not in vivo hepatic residual
capacity, which is what `theta` is. It is support, not settlement. Still the PI's call.

### 3.3 Tier 1. Two things the pack gets right, and one sentence of the paper it breaks

**`08_specification.tex` says "the cohort is methylmalonic acidemia, not propionic acidemia".
That is true of Manoli 2021 and false of the literature.** Barshop 1991, verified: five PA
patients, `1-13C`-propionate, and

> Recovery of `13CO2` and half-time of `1-13C`-propionate in PA were significantly less than
> normal... Recovery of `13CO2` was **well correlated with clinical severity in PA**, but did
> not correlate in MMA.

That is stronger for the paper than what it currently claims, not weaker: the assay tracks
severity in PA specifically, and does not in MMA. The correct statement is that no **repeatability**
figure exists for PA, which the review also says and which is what Tier 1 actually needs.
Barshop reports no CV.

**Wagner 2011 is closer to Tier 1 than anything the paper cites.** 33 of 58 subjects received a
second `1-13C`-propionate test **1 to 3 days** later, against Manoli's 1 to 3 weeks and 2 months
to 4 years, and the two values were "on average within 10% of each other". Same tracer. Against
it: a B12-deficiency population, a `%dose/hr` readout at 10 or 20 min rather than the 120-min
cumulative, and the review itself notes the printed confidence interval is internally
inconsistent (an upper bound of 9.8% under a point estimate of 10%). It is a lead to read in the
primary text, not a number to import. If it holds, it sits just inside the paper's 11.7%
requirement, which would make Tier 1 a confirmation rather than a gamble.

**Tier 1 gains a sample size.** The review's arithmetic, re-derived here and correct: at the
control within-subject CV of 6.9%, the repeatability coefficient for a single pre/post pair is
`2.77 x 6.9 = 19%`; `k` replicates per timepoint give `2.77 x 6.9 / sqrt(k)`, so `13.5%` at
`k = 2` and `11.0%` at `k = 3`; a paired group study at 80% power for a 12% change needs about
6 subjects at the control CV and about 21 at the patient CV. Tier 1 currently specifies "one
cohort, two breath tests bracketing a single dose, on one day" and states no `n`.

One caveat to carry with those numbers. They size the detection of a 12% **change**, whereas the
paper's 11.7% is a requirement on the precision with which the **increment** is characterised.
The binding quantity is the increment as a fraction of baseline oxidation, which is the thing
Tier 1 exists to measure, so the sizing cannot be fixed without assuming it. The review's `n`
is conservative for this purpose, and should be quoted as a bracket rather than a design.

### 3.4 The quasi-steady separation is now the weakest assumption in the design, and the review says why

`07_design.tex` rests on a separation of 104 between the marker's `0.81` h time constant and the
`3.5` d therapeutic half-life, and `rem:transient` already names the monotone-tracking assumption
as the weakest link. The 3-HP review's finding is that **no human plasma 3-HP half-life,
clearance or renal excretion figure exists**, which confirms `PHASE4_SUMMARY.md` section 5a and
adds nothing new. What is new is the dialysis observation: the pool refills after intermittent
haemodialysis stops, from ongoing extrahepatic production and redistribution from tissue
compartments. **A slow peripheral pool feeding the plasma compartment is exactly the structure
that would break the 104-fold separation**, and it is a different objection from the one
`rem:transient` already records.

This does not change a number. It changes what `rem:transient` should say: the untested
assumption is not only that the marker tracks the substrate monotonically, it is that the plasma
pool has no slow compartment behind it, and the dialysis rebound literature is the one human
observation that bears on it. Barshop's "differences between MMA and PA may indicate different
distribution of propionate pools" is the same suggestion from the other direction.

### 3.5 Prior art. No change, one confirmation

The prior-art review's conclusion that "stroboscopic sampling" and "phase-folded sampling" have
no pharmacological precedent agrees with `literature_pending.md`'s existing instruction: no
single source states the dose-synchronous claim verbatim, so the paper should say so rather than
imply one exists. The van Noort exclusion stands. The dialysis residual-clearance literature
(Daugirdas, Depner, Obi) and C-peptide deconvolution (Van Cauter, Sparacino and Cobelli,
Mortensen) are genuine adjacent prior art that this repository has not previously named, and are
worth citing in the discussion as the nearest precedents for estimating a residual function from
a periodic schedule. None of them uses the identification framing, so none of them costs the
paper its contribution.

---

## 4. What was not done

- **No prose was edited.** The paper and supplement remain frozen.
- **No constant was moved.** `LOAD_FRAC` is still `0.05`.
- **Nothing moved into `paper/refs.sources`.** The nine identifiers above are verified and could
  move, but that is a registration decision, not an audit finding.
- **The `sigma` sweep is not in `spec_numbers()`.** It lives in `audit/m1_loadbias_prior.py`
  only, pending a decision on 3.1.
- Wagner 2011 was verified as a record, not read. The 10% figure is quoted from the review.

---
---

# Part 2. The full verification pass, and what reading the primary sources found

Part 1 was written after resolving nine identifiers. This part completes
`refs/literature_pending.md` and reports what came out of reading the primary sources rather
than the reviews' descriptions of them. **Several of the most consequential findings below are
in papers no review in the pack cited.**

## 5. `refs/literature_pending.md` is now fully verified

Every row resolved against CrossRef or PubMed. **All 16 rows are correct as written.**

**Methods, 11 rows, all confirmed:**

| row | registrar says | note |
|---|---|---|
| Raue 2009 | *Bioinformatics* 25(15):1923-1929 | exact |
| Bates & Watts 1980 | *JRSS-B* 42(1):1-16, "Relative Curvature Measures of Nonlinearity" | exact |
| Self & Liang 1987 | *JASA* 82(398):605-610 | exact |
| Andrews 2000 | *Econometrica* 68(2):399-405 | exact |
| Lavielle & Aarons 2016 | *JPKPD* 43(1):111-122 | exact |
| Janzén 2016 | *Front Physiol* 7, "Parameter Identifiability of Fundamental Pharmacodynamic Models" | family name carries an accent, **Janzén** |
| Golub & Pereyra 2003 | *Inverse Problems* 19(2):R1-R26 | exact; 1973 founding paper confirmed at *SIAM J Numer Anal* 10(2):413-432 |
| Nyberg & Karlsson 2009 | **Nyberg, Karlsson, Hooker**, *JPKPD* 36(2):125-145 | third author is Hooker |
| **Walter & Pronzato 1990** | *Automatica* 26(2):195-213 | **the "Verify (review)" row. Confirmed** |
| **Mehra 1974** | *IEEE TAC* **19(6):753-768** | **the citation trap. Confirmed exactly as this file has it** |
| **Segel & Slemrod 1989** | *SIAM Rev* 31(3):446-477 | **the second "Verify (review)" row. Confirmed** |

**Held references and the exclusion:**

- `natural_history` `10.1186/s13023-026-04251-3` resolves to Schwahn, Berry, Vernon et al.,
  *Orphanet J Rare Dis* 21(1), 2026. Already at `refs.sources:36`. Correct.
- `breath_test_mma` confirmed at *Genet Med* 23(8):1522-1533.
- **van Noort 2024** resolves to "Navigating the landscape of parameter identifiability
  methods: a workflow recommendation for model development", *CPT PSP* 13(7):1170-1179. It is
  a methods-workflow review. **The exclusion stands and is now checked at the source.**

**Clinical, 4 rows, all confirmed**, and three of them upgraded from "would support" to a
verbatim quotation. See Part 1 section 2.1 (Thompson), 3.1 (Daly) and 3.2 (Hsia).

**Nothing in `literature_pending.md` was wrong.** Every row can move into `paper/refs.sources`
on a registration decision.

## 6. Four further misattributions in the pack, none previously reported

Added to the three in Part 1 section 2.2, these make **seven** citation errors across the pack.

| the pack says | the registrar says |
|---|---|
| "Ogungbenro K, Aarons L, et al., *JPKPD* **2007**", `10.1007/s10928-006-9009-1` | **Gueorguieva**, Aarons, Ogungbenro, Jorga, Rodgers, Rowland, *JPKPD* **2006;33(2):97-124**. Title correct, first author, year and volume wrong |
| "Sharma/Kumar et al. 2025", `10.1038/s41598-025-09110-y` | **Devi**, Pasanna, Ayoob, Sachdev, Thomas, Fiehn, **Kurpad**, *Sci Rep* 15(1). Real and topical, attribution invented |
| "Al-Dirbashi/Chakrapani et al.", PMID 25943030 | **Turki**, Murthy, Ueda, Cheng et al., *Mol Genet Metab* 2015;115(2-3):78-83 |
| "Kasicka-Jonderko et al. 2018", PMID 29409350 | **Bożek**, Kamińska, Kasicka-Jonderko, Krusiec-Świdergoł. Kasicka-Jonderko is third author |

Smaller: "Perri" for PMID 15797668 (actually **Barbosa**, *Nutrition* 2005;21(3):289-94);
"Dijkman SC" for **van Dijkman SA**; "Evans/MacDonald group" for *Nutrients* 2023;15:4903,
whose first author is **Pinto**. Verified as correct and unremarkable: Baumgartner 2014,
Cummings 1987, Sasaki 2023, Kasicka-Jonderko 2007 and 2011, Afolabi 2011 and 2013, Choi 1997,
Ziegler 1996, Gharaibeh 1998, Di Campli 2003, Fuller 2000, Daugirdas 2010, Obi 2018
(*NDT* 33(3):530-539, the article number the review flagged), Daugirdas & Bolasco 2025,
Shafi 2016, Van Cauter 1992, Sparacino & Cobelli 1996, Mortensen 2007, Ruan 2019,
Qi 2021 (*Clin Transl Sci* 14(5):1894-1905, the volume the review flagged), Feillet 2008,
Kaufman 1999, Langenbeck 2001, Griffin & Bradshaw 2019, Rossato & Fuchs 2021, van Rijn 2011,
Evans 2022, the mFood paper, and **Haitjema 2026** (*Mol Genet Metab* 148(4):110183, real).

## 7. A correction to Part 1. Wagner 2011 cannot carry the weight I gave it

Part 1 section 3.3 said that if Wagner's 10% figure held, Tier 1 would become a confirmation
rather than a gamble. **Having now read the full text (PMC3204151), it does not hold.** The
reproducibility paragraph in full:

> "A subset of all subjects (33 of 58) received a second breath test 1 to 3 days after initial
> testing to evaluate the reliability of the BBT." ... "The Kappa value for the agreement
> between tests was 0.7415 (95% CI: 0.4664 - **1.0166**). The BBT value of the two tests was on
> average within 10% of each other (95% CI: **7.5 - 9.8%**)."

A kappa cannot exceed 1, and a confidence interval cannot exclude its own point estimate.
**Two impossible statistics in one paragraph.** The kappa is also categorical agreement on a
deficient/not-deficient call, not a continuous CV, and the readout is `%dose/hr` at 10 or 20
minutes in a B12 population. Wagner remains the **shortest-interval repeat of a
`1-13C`-propionate breath test in the literature**, and that fact is worth one sentence, but
its numbers are not usable. **Tier 1 stays where the paper had it.**

## 8. The trial's dose train is not a fixed biweekly schedule

Read from the Nature 2024 paper (PMC11156579), which the repository holds as `trial_phase12`.

- **Cohort 1 is Q3W, not Q2W.** Cohort 1: 0.30 mg/kg **every three weeks**, n = 4. Cohorts 2-5:
  0.30, 0.45, 0.60, 0.90 mg/kg **every two weeks**, n = 3 each. Dose amplitude spans **threefold**.
- **Five of 16 participants changed dosing regimen.** One had a dose decrease (0.90 to 0.60
  after grade 3 pancreatitis), **skipped dose 3 entirely**, and had dose 2 delayed by
  hospitalisation. Two cohort-1 participants moved to Q2W, two cohort-2 participants to 0.60
  Q2W, one cohort-5 participant to 0.60 Q2W.
- **The schedule slips.** For a cohort-1 (Q3W) participant the paper records "the predose sample
  for dose 3 (**52 days after the first dose**)". Nominal Q3W puts dose 3 at day 42. That is
  about **ten days of slippage on a 21-day nominal interval**.
- Per-regimen dose counts are published: 73 doses at 0.30 Q3W, 136 at 0.30 Q2W, 68 at 0.45 Q2W,
  45 at 0.60 Q2W, 24 at 0.90 Q2W, total 346.

`src/design/run.py` annotates `DOSE_INTERVAL = 14.0` as "the biweekly schedule selected for
mRNA-3927". That is one of two published schedules, and neither is adhered to exactly.

**Which way this cuts is worth stating carefully, because it is mostly good.** A ragged,
amplitude-varying, occasionally-skipped dose train is **better excitation** than a uniform one,
and slippage makes `rem:aliasing`'s failure mode (every draw at the same phase) much harder to
fall into by accident. What it costs is the paper's framing: the claim that "the **fixed**
biweekly schedule the trial already runs turns out to be the identification experiment" should
become a claim about the dose train as administered, which is better for the argument and
requires the sentence to change. `D_AMP`, a single assumed amplitude, is also standing in for
five published dose levels.

## 9. Extended Data Table 3. Narrowed, with a test the repository can run on numbers it already has

`OPEN_QUESTIONS.md` section 2 calls this "the single largest swing in the design" and says the
public record is exhausted. It is narrower than that:

- The table's own title is "**Biomarker percentage change from baseline after the third dose**".
- The Methods define the comparison: "changes in blood 2-MC and 3-HP levels from baseline
  (**predose levels**) to post-dose levels".
- The Extended Data tables are published as **images** (`41586_2024_7266_Tab3_ESM.jpg`), which is
  why no timepoint appears in any machine-readable version.

The phase of the post-dose reading is still unstated. **But the table reports the Q2W and Q3W
cohorts separately**, and that gives a discriminator the repository can apply to the numbers it
already has: if the post-dose sample is at a fixed clock time after dose 3, the Q2W and Q3W
cohorts sit at the **same** phase; if it is the following predose, they sit at **14 and 21 days**,
different phases of a decay with a 3.5-day half-life. A systematic Q2W-versus-Q3W offset in the
tabulated falls is evidence for the trough reading; its absence is evidence for a fixed clock
time. This does not need the sponsor.

## 10. Three published corrections, none cited here. One of them settles R4

- *Nature* **Author** Correction, 629(8011):E10, `10.1038/s41586-024-07413-0`.
- *Nature* **Publisher** Correction, 630(8017):E13, `10.1038/s41586-024-07646-z`.
- *Genet Med* Correction, 23(11):2233, `10.1038/s41436-021-01188-9`, to `breath_test_mma`.

The Manoli correction was read in full (PMC8553604). In its entirety: "Due to a processing
error, the online graphical abstract was not given. The original article has been corrected."

**It does not touch the coefficient of variation.** So the abstract's `12.8%`, which is the
Results' `12.86` standard deviation misreported as a mean, **stands uncorrected in the published
record**. `REVIEW_LOG.md` **R4 is not superseded**, and the reviews will keep inheriting the
error because the source has never fixed it.

## 11. Tier 1 is cheaper and better supported than the paper assumes

All from Manoli's own Methods and Figure 1 (PMC8354855), read here for the first time.

1. **A metabolic cart is not needed for a short test.** Figure 1a: at the 0.5 mg/kg dose, "when
   predicted CO2 rather than the VCO2 measured by the metabolic cart was used for calculations,
   **results were identical up to the 60-minute timepoint**". Indirect calorimetry is the main
   logistical obstacle to running Tier 1 outside a research unit, and a pre/post pair read at
   60 minutes does not need it.
2. **The signal peaks early.** Time of maximum enrichment **22.5 ± 12.9 minutes**. Tier 1 does
   not need the 2-hour protocol.
3. **The dose is small and the result is dose-invariant.** 0.5 mg/kg, twenty-fold below
   Barshop's 9.7 mg/kg, with 48.8%, 56.03% and 55.91% oxidised at 9.7, 5 and 0.5 mg/kg.
4. **The paper's central Tier 1 argument can now be cited instead of asserted.**
   `08_specification.tex` says between-visit variation is an upper bound on same-day error
   "because it contains drift a same-day pair removes". Manoli's Discussion says exactly that in
   the source's own words: oral bolus results are affected by "fasting state, gut motility, and
   microbiome differences, but also, importantly, by the **dilution of the tracer in a large
   endogenous metabolite pool** and the **effects of renal function on its clearance**,
   parameters that **can vary greatly between follow-up clinic visits**."
5. **Manoli states the sensitivity claim Tier 1 rests on**: the test "can detect smaller
   enzymatic activity changes after an isolated kidney transplant or B12 therapy when tested
   **pre- and postintervention in individual subjects**".
6. **A protocol conflict to resolve.** Wagner found an **8-hour fast essential** (a 1-hour fast
   produced 4/10 false positives). Manoli used "a minimum of **3-hour** fast... **as suggested by
   the study of Wagner et al.**", which is shorter than the minimum Wagner established. Tier 1
   should fix the fast at 8 hours or longer and say why.

## 12. The biological-variation gap is real for 3-HP and 2-MC, and **not** real for serum MMA

Both reviews call within-subject variability for these markers a genuine, unfilled gap. For 3-HP
and 2-MC that is correct. For serum MMA it is not, and the data are in the paper this repository
already cites:

- **Between-visit** serum MMA CV in the repeat-tested subjects: range **1.5% to 75.3%, average
  31.38 ± 24.55%** (against the breath test's 14.02 ± 12.86% in the same subjects).
- **Within-week**: Manoli Figure 1h plots serum MMA means and standard deviations for **57
  subjects with more than one measurement during a week-long CRC stay**.

That is a within-subject biological-variation dataset for the serum marker in this disease,
n = 57, which both reviews missed. It does not close the 3-HP gap, but "no data exist for any of
these markers" is too strong and should not be written that way.

## 13. The finding no review produced: an unrecorded antibiotic cycle costs a factor of 1.8

Two primary sources, neither cited anywhere in the pack.

**Baumgartner 2014 (PMC4180313), the regimen, verbatim:** "Metronidazole, 10-20 mg/kg/day
divided in 2-3 times... In order to avoid the development of drug-resistant colonies, **1-2
weeks of therapy alternating with 2-3 weeks off or alternating every month is advisable**."
This closes the "verify exact dose/schedule" item.

**Thompson GN, Chalmers RA, Walter JH, Bresson JL, Lyonnet SL, Reed PJ, Saudubray JM, Leonard
JV, Halliday D, "The use of metronidazole in management of methylmalonic and propionic
acidaemias", *Eur J Pediatr* 1990;149(11):792-6, PMID 2226555, `10.1007/BF01957284`.** Verified.
**Correction to this section as first written: this paper was already registered here as
`metronidazole_pa`, with these figures already in its annotation, and Baumgartner was already
held as `pa_guidelines`. Neither is new to the repository. What is new is the conjunction, and
the heading "two primary sources, neither cited anywhere in the pack" is true only of the pack.**
**Four PA and five MMA children.** "Urinary excretion of propionate metabolites fell
significantly during the treatment in all subjects, the mean decrease being **41% (range 12-76,
P less than 0.01)**, while mean plasma propionate was reduced from **45.0 to 25.1 umol/l**
(P less than 0.05)." This is the source of `literature_pending.md`'s "metronidazole ~-40%", now
verified with an n and a range, and it is a **different paper** from the Thompson 1990
*Metabolism* partitioning study. It is also, incidentally, a human plasma propionate
concentration in PA/MMA before and after a load intervention, which the 3-HP review lists as
absent.

**Put together: in a patient on the standard regimen the non-dietary load `b` is not constant.
It is close to a square wave of order 40% amplitude on a 21 to 35 day period**, sitting
alongside a 14 or 21 day dose train and a 30.4 day sampling grid.

The paper knows about cyclical antibiotics and treats them as one of the **pre-therapy
instruments** (`A1_supplementary.tex`, and the "three readings" in `07_design.tex`, which
correctly conclude they act below `K_m` and move the patient the wrong way). **It never asks
what a concurrent antibiotic cycle does to the dose-train design.** `audit/m2_antibiotic_cycle.py`,
24 seeds, monthly draws:

| antibiotic period | cycle **in** the load record | cycle **not** in the record |
|---|---|---|
| no cycle at all | `0.1660` | `0.1660` |
| 21 d | `0.1709` | **`0.2974`** |
| 28 d | `0.1716` | **`0.3044`** |
| **30.4 d** (the sampling interval) | `0.1665` | `0.1665` |
| 35 d | `0.1712` | **`0.3024`** |
| 60 d | `0.1731` | **`0.3128`** |

Three things follow.

1. **Recording the antibiotic weeks costs nothing and removes a factor of 1.8.** Recorded, the
   cycle neither helps nor hurts (`0.166` to `0.173`), which is exactly what the paper's own
   "three readings" predicts for a perturbation an order of magnitude below `K_m`. Unrecorded,
   it takes the design from `0.166` to about `0.30` and **fails the `0.20` target by half**.
   The instruction is one line and it is free: **record antibiotic weeks alongside intake**.
   It is worth as much as the entire load-precision discussion in Part 1 section 3.1.
2. **The aliasing case is the benign one here, which is the opposite of the dose-train result.**
   At a period equal to the sampling interval every draw lands at the same phase, the swing
   becomes a constant offset, and the load-bias column absorbs it. A reader who has just been
   told that aliasing destroys the dose-train result will assume aliasing is always the enemy.
   On an unrecorded nuisance input it is the rescue.
3. **The two audit scripts corroborate.** `m2`'s unrecorded case runs at an effective
   `sigma` near `68`, and `m1`'s independent `sigma` sweep gives `0.264` at `sigma = 66.2`.
   Different code paths, consistent answers.

## 14. Two smaller things worth carrying

- **The trial's efficacy figure needs its interval.** The abstract's "70% reduction in the risk
  of metabolic decompensation events" is `1 - 0.30`, and the paper gives the relative risk as
  **0.30 (95% CI 0.066 to 1.315), P = 0.0927** overall, and **0.18 (95% CI 0.011 to 2.907),
  P = 0.1832** among the biweekly cohorts. **Neither is significant and both intervals cross 1
  substantially.** Every review in the pack quoted the 70% without the interval.
- **Dense post-dose sampling already happens in this trial, on the wrong analyte.** "Blood
  samples were collected before dosing and at the indicated time points for up 21 days after
  dose 1" for mRNA PK. `rem:transient` asks for exactly this schedule on the marker, and says
  "whether that is clinically feasible is a different question". The sampling infrastructure and
  consent for it demonstrably exist. That strengthens the remark.

## 15. What Part 2 did not do

- **No prose edited, no constant moved, nothing registered.** Unchanged from Part 1.
- **The Extended Data Table 3 discriminator in section 9 was not run.** It needs the tabulated
  per-cohort falls, which are in an image the repository already holds and this container does
  not.
- **`audit/m2_antibiotic_cycle.py` models the antibiotic cycle as a square wave on the load.**
  That is the shape the regimen implies, not a shape anyone has measured.

---
---

# Part 3. The discriminator, run; and what is now discharged

## 16. Section 9's discriminator fails, and section 9 is withdrawn

The per-cohort values of Extended Data Table 3 were supplied by the author from the Nature PDF.
The test proposed in section 9 was run with **its power printed before the table was read**,
which `audit/m3_edt3_discriminator.py` enforces structurally: Part A executes before Part B.

**Power, established first.** The only contrast free of the dose confound is cohort 1
(0.30 mg/kg Q3W, `n = 4`) against cohort 2 (0.30 mg/kg Q2W, `n = 3`). Cohorts 3 to 5 differ from
cohort 1 in dose, interval **and** age together.

| marker | SE of the difference | df | smallest difference detectable at 80% power |
|---|---|---|---|
| 3-HP | 19.5 pp | 4.3 | **70.8 pp** |
| 2-MC | 21.2 pp | 4.9 | **74.2 pp** |
| C3 | 17.2 pp | 4.3 | **62.6 pp** |
| n-PG | 17.0 pp | 5.0 | **59.5 pp** |

Against that, the two conventions differ by `8.1 pp` at the assumed `D_AMP = 6000` and by
`22.5 pp` at `D_AMP = 24000`, four times the assumed value. The fixed-clock convention predicts
`1.2 pp`. **The test is underpowered by a factor of three to thirty and cannot separate the two
sampling conventions at any plausible amplitude.**

**The data, read second, agree.** `p = 0.249` (3-HP), `0.964` (2-MC), `0.215` (C3), `0.911`
(n-PG). The two markers that move at all move `-26` and `-25` pp, which is the **direction
opposite** to the trough hypothesis: a Q3W arm sampled at its next predose should show a
*smaller* fall, not a larger one. There is no dose trend to lean on either: over the four Q2W
dose levels the dose-versus-fall correlation is `-0.37`, `-0.50`, `-0.16` and **`+0.67`**, on
four points.

**Part D of the script adds nothing beyond what `PHASE5_ROBUSTNESS.md` already records.** All
twenty cohort-by-marker means sit above the `-61.8%` a one-day-post-dose reading predicts at
`D_AMP = 6000`, but only three of twenty reach `p < 0.05`, the cohorts within a marker are not
independent, and the comparison is conditional on an assumed amplitude that a smaller value
would move up to meet the data. Phase and amplitude cannot be separated by it.

**Recorded as `REVIEW_LOG.md` R18.** `OPEN_QUESTIONS.md` section 2 stands: the sampling time
must come from the sponsor. What is now pinned is only that the table is a change from a
**predose** baseline measured **after the third dose**.

Stating the power first is what produced the answer before the data were read. Had the table
been read first, `-26` and `-25` pp on two of four markers would have looked like a finding.

## 17. `refs/literature_pending.md` is discharged

All 16 rows resolved and all 16 were correct. Registered into `paper/refs.sources`, which now
holds **46 identifiers, no duplicate keys and no duplicate identifiers**, and
`tools/fetch_refs.py --check` resolves every one with an empty `paper/refs.todo`.

Two rows needed no new entry because the repository already held them, **which corrects section
13 above**: Thompson 1990 *Metabolism* is `propionate_sources` and Manoli 2021 is
`breath_test_mma`. Two further sources this pass turned out to be already held with correct
annotations, `pa_guidelines` and `metronidazole_pa`, the latter already carrying the `41%`,
the `12-76` range and the `45.0` to `25.1` umol/L plasma fall. **The metronidazole finding is a
conjunction of two sources the repository had, not a discovery of sources it lacked.**

Twelve method identifiers, four clinical ones and three corrections are new. Annotations were
extended on three existing entries:

- `trial_phase12` now carries the MDE relative risk **with its interval**, `0.30 (95% CI
  0.066-1.315, P = 0.0927)` overall and `0.18 (95% CI 0.011-2.907, P = 0.1832)` biweekly, and
  records that the abstract's "70% reduction" is `1 - 0.30` and is **not significant**, with an
  instruction not to quote it bare. It also now records the non-uniform dose train and what
  Extended Data Table 3's title does and does not fix.
- `breath_test_mma` now carries its published correction and the finding that the correction
  restores a graphical abstract and **leaves the CV error standing**, so `R4` is not superseded.
- `pa_guidelines` now carries the metronidazole regimen verbatim, which was an open
  verification item, and the dialysis rebound that passage 14 rests on.

**`van_noort` was deliberately not registered.** It resolved and was read; it is a workflow
review of identifiability methods and is not support for the dose-phase result. Keeping it out
of `refs.bib` means it cannot be cited here by accident. The reason is recorded in
`refs.sources` itself.

## 18. State of the gate

- **No `.tex` file has been touched.** `git status` confirms it. `paper/refs.bib` regenerated
  from `refs.sources`, which is its defined behaviour; no `.tex` cites any new key, so the built
  PDF's bibliography is unchanged. The PDF was not rebuilt.
- **No constant moved.** `LOAD_FRAC` is `0.05`, `D_AMP` is `6000`, `DOSE_INTERVAL` is `14.0`.
- **Five proposed passages (11 to 15) are held unapplied** in `PHASE4_PROPOSED_TEXT.md`,
  alongside the ten that were already there. None is in the paper.
- **Two retractions recorded**, `R17` (the `spec_numbers()` sweep) and `R18` (my own
  discriminator proposal).
- Three audit scripts added: `m1_loadbias_prior.py`, `m2_antibiotic_cycle.py`,
  `m3_edt3_discriminator.py`. None is in `src/`, none has a test, and that is deliberate per
  `OPEN_QUESTIONS.md` section 5.
