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
