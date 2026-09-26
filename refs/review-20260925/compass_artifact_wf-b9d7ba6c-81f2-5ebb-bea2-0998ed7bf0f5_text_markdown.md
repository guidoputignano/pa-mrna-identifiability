# k_e (non-enzymatic 3-HP clearance): can it be pinned down well enough to identify θ? — Literature review for `01_ke_clearance.md`

**Scope.** Whether the first-order clearance term `k_e` in the quasi-steady balance `s = (θ + Δ(t))·m/(K_m+m) + k_e·m + b` — the disposal of plasma 3-hydroxypropionate (3-HP) by every route that does not pass through propionyl-CoA carboxylase (PCC) — can be constrained to the ~2% (at 5% residual capacity) or ~6% (at 10%) per-patient precision the 25 Sep 2026 robustness review says θ-recovery requires. Every quantitative claim below was checked against a primary source that was actually opened; each is labelled VERIFIED (read in the primary source), SECONDARY (seen only in a review/citing paper), or NOT FOUND. Species and population are stated because the repository has previously mistaken a mouse estimate for a human one and an SD for a mean.

---

## TL;DR

- **No human elimination kinetics for 3-HP exist.** There is no measured plasma half-life, systemic clearance, volume of distribution, fraction excreted unchanged, or turnover rate for 3-hydroxypropionate in any human population — PA, MMA, or healthy. `k_e` therefore cannot currently be sourced from the literature at any precision, let alone the 2–6% the design needs, and the paper's factor-of-2.2-per-SD prior is the honest state of knowledge, not a placeholder that a citation can replace.
- **`k_e` is patient-specific, not a population constant.** 3-HP is renally secreted by the organic anion transporter OAT1 (its plasma level rises in Oat1-knockout mice),\[1\] 4/8 (50%) of adults with PA had eGFR <60 mL/min/1.73 m² per Shchelochkov OA et al., Genet Med 2019;21(12):2830–2835 (PMID 31249402): "Among adults with PA, 4/8 (50%) had eGFR <60 mL/min/1.73 m²"; and by cystatin-C 25/30 (82%) had eGFR <90 (CKD stages 2–4), and carnitine status diverts propiogenic flux into a competing excretory route — three independent reasons the non-enzymatic disposal rate varies between patients and within a patient over time. A single population `k_e` prior cannot be correct for every patient.
- **The one feasible way to measure `k_e` per patient is an IV micro-dose of ¹³C₃-3-HP with serial plasma sampling**, mirroring the bolus [1-¹³C]propionate studies already done safely in PA/MMA children; a timed urine collection alone measures only the renal fraction of `k_e·m` and would underestimate total non-enzymatic disposal.

---

## 1. Summary table

| # | Quantity the model needs | Best available value | Population / species | Mean or SD? | Label |
|---|---|---|---|---|---|
| 1 | 3-HP plasma half-life | none | human (any) | — | **NOT FOUND** |
| 2 | 3-HP systemic clearance (mL/min) | none | human (any) | — | **NOT FOUND** |
| 3 | 3-HP volume of distribution | none | human (any) | — | **NOT FOUND** |
| 4 | 3-HP fraction excreted unchanged (renal) | none quantified | human (any) | — | **NOT FOUND** |
| 5 | 3-HP turnover/disposal rate | none | human (any) | — | **NOT FOUND** |
| 6 | Rat perfused-liver 3-HP uptake | 806 ± 149 nmol·min⁻¹·(g dry wt)⁻¹ (2 mM 3-HP alone); ≈20% of propionate uptake\[2\] | rat liver, perfused | mean ± SD | VERIFIED |
| 7 | Surrogate: blood propionate half-life | 4–5 min controls; 3–7× longer in PA | human PA + controls | point values | VERIFIED |
| 8 | Surrogate: propionate turnover, healthy adults | 17.6 μmol·kg⁻¹·h⁻¹ (SD 5.9)\[3\] | healthy adults | mean (SD) | VERIFIED |
| 9 | Surrogate: propionate turnover, patients | 0.41–4.98 μmol·min⁻¹·kg⁻¹ | human PA + MMA | range | VERIFIED |
| 10 | 3-HP is an OAT1 (renal) substrate | plasma 3-HP elevated in Oat1-knockout\[1\] | mouse | qualitative | VERIFIED (mouse) |
| 11 | CKD prevalence, adult PA | eGFR <60 in 4/8 (50%) adults; <90 in 25/30 (82%) | human PA | proportions | VERIFIED |
| 12 | Total measured urinary metabolite excretion, PA | fed 1.45, fasting 2.98 μmol·kg⁻¹·h⁻¹\[4\] | human PA (n=3) | means | VERIFIED |
| 13 | Total measured urinary metabolite excretion, MMA | fed 38.1, fasting 54.6 μmol·kg⁻¹·h⁻¹\[4\] | human MMA (n=5) | means | VERIFIED |
| 14 | Urinary 3-HP/methylcitrate pattern day-to-day | "wide variations" (qualitative)\[5\] | human PA (n=8) | qualitative | VERIFIED |
| 15 | ¹³C₃-3-HP labelled compound availability | not offered as a stocked catalogue item; synthesis routes published for related ¹³C-hydracrylates | reagent | — | SECONDARY |
| 16 | Precedent: IV bolus ¹³C-propionate in children | 30 μmol/kg IV bolus, plasma at 5–90 min\[6\] | human PA/MMA children | protocol | VERIFIED |

---

## 2. Evidence

### 2.1 Direct human 3-HP elimination kinetics — none exist

The most authoritative statement comes from the only group to have studied 3-HP disposal experimentally. In the rat-liver perfusion paper (Wilson KA et al., *Am J Physiol Endocrinol Metab* 2017;313:E413–E428, DOI 10.1152/ajpendo.00105.2017, PMID 28634175), the authors write that "There is much less information on mammalian metabolism of 3HP, which is present at very low concentrations in normal human plasma and urine."\[2\] They cite no human clearance, half-life or turnover value, and neither does any other source located. **VERIFIED that the datum is absent** (i.e., a genuine "no human data exist" finding), not merely un-found.

The classic human tracer work — Ando T, Rasmussen K, Nyhan WL, Hull D, "3-Hydroxypropionate: significance of β-oxidation of propionate in patients with propionic acidemia and methylmalonic acidemia," *Proc Natl Acad Sci USA* 1972;69(10):2807–2811 (DOI 10.1073/pnas.69.10.2807) — gave IV [1-¹⁴C]propionate and detected labelled 3-HP in urine,\[7\] establishing the β-oxidation shunt qualitatively,\[7\] but reports no half-life, clearance, Vd or fractional urinary excretion of 3-HP. Population: human PA and MMA patients. **VERIFIED (qualitative only).**

### 2.2 The closest quantitative surrogate: rat perfused-liver 3-HP uptake

Wilson et al. 2017 measured net 3-HP uptake by isolated recirculating-perfused livers of overnight-fasted male Sprague-Dawley rats (values per gram dry weight, mean ± SD, n = 5–6). Basal 3-HP release by control liver was 0.86 nmol·min⁻¹·(g dry wt)⁻¹ with accumulation of only 3 μM after 1 h;\[2\] with 2 mM [¹³C₃]3-HP supplied, 3-HP uptake was 806 ± 149 nmol·min⁻¹·(g dry wt)⁻¹, rising to 947 ± 476 and 1,517 ± 411 when 5 mM propionate was co-perfused.\[2\] Critically, "The basal rate of 3HP uptake was only 20% of that of propionate uptake,"\[2\] and "the net rate of 3HP release by the liver was <3% of the rate of propionate uptake."\[2\] The paper also shows the two metabolites compete: "The rate of propionate uptake decreased by one-half in the presence of 3HP."\[2\] Species: rat, isolated liver. **VERIFIED.** Interpretation: hepatic 3-HP disposal is real but Wilson KA et al., Am J Physiol Endocrinol Metab 2017;313:E413–E428 (PMID 28634175): "The basal rate of 3HP uptake was only 20% of that of propionate uptake," with 3HP-uptake variability "coefficient of variation = 42%" in the presence of propionate, and the liver is not a single-pass sink for 3-HP the way it is for propionate — consistent with a slower, distributed, partly renal `k_e`.

### 2.3 Propionate kinetics as an upper-bound surrogate for `k_e`

Because 3-HP shares the propionyl-CoA pool, propionate kinetics bound the timescale. Barshop BA et al., "Metabolism of 1-¹³C-Propionate in Vivo in Patients with Disorders of Propionate Metabolism," *Pediatr Res* 1991;30(1):15–22, gave a 30 μmol/kg IV bolus (or 100 μmol/kg orally) of sodium 1-¹³C-propionate\[6\] and sampled heparinised blood at 5, 10, 15, 20, 30, 45, 60 and 90 min.\[6\] "The half-life of ¹³C-propionate was between 4 and 5 min in each of the three control experiments"; "In the patients with PA, the values ranged from about 3 to 7 times greater than control."\[6\] Whole-body propionate turnover in patients was 0.41–4.98 μmol·min⁻¹·kg⁻¹.\[6\] Population: human PA/MMA children + controls; blood propionate measured by GC-MS. **VERIFIED (full text read).**

For a between-subject spread in a related first-order disposal parameter, the nearest primary datum is Thompson GN et al.'s [²H₅]/[¹³C]propionate infusion in 4 healthy adults: mean propionate turnover 17.6 μmol·kg⁻¹·h⁻¹ (SD 5.9),\[3\] i.e. a between-subject CV of roughly one-third in healthy people with intact kidneys. Population: healthy adults. **VERIFIED (abstract).** This is the best-supported analogue for a `k_e` CV; no CV for 3-HP disposal itself can be given.

### 2.4 3-HP is renally handled (OAT1) — so renal function drives `k_e`

3-HP is a substrate of the renal organic anion transporter OAT1/SLC22A6. In Nigam SK et al., "The Organic Anion Transporter (OAT) Family: A Systems Biology Perspective," *Physiol Rev* 2015;95(1):83–123 (DOI 10.1152/physrev.00025.2013), Table 2 lists "3-Hydroxypropionate — Elevated plasma concentration" among metabolites disturbed in the Oat1 knockout,\[1\] citing reference 74\[1\] = Eraly SA et al., "Decreased renal organic anion secretion and plasma accumulation of endogenous organic anions in OAT1 knock-out mice," *J Biol Chem* 2006;281(9):5072–5083 (PMID 16354673). Species: mouse knockout; the evidence is a knockout metabolomic phenotype, not a direct Km/Vmax transport assay. **VERIFIED (mouse, indirect).** Implication for the model: the renal component of `k_e` scales with functional OAT1-mediated secretion and glomerular filtration, both of which fall with the CKD that is common in PA.

### 2.5 CKD is common in adult PA — a first-order driver of between-patient `k_e` variation

Shchelochkov OA et al., "Chronic kidney disease in propionic acidemia," *Genet Med* 2019;21(12):2830–2835 (PMID 31249402): "Among adults with PA, 4/8 (50%) had eGFR <60 mL/min/1.73 m²,"\[8\] and by cystatin-C eGFR "25/30 subjects (82%) had eGFR <90 mL/min/1.73 m² (chronic kidney disease stages 2–4)."\[8\] The paper also warns that creatinine-based eGFR performs poorly in PA because of low muscle mass and low-protein diet.\[8\] Population: human PA. **VERIFIED (abstract + figure legend).** Because renal secretion is part of `k_e`, a 50% CKD prevalence with unreliable eGFR means the renal component of `k_e` is both variable and hard to measure indirectly.

### 2.6 Carnitine status opens a competing excretory route

L-carnitine conjugates accumulating propionyl-CoA to propionylcarnitine, which is renally excreted (Roe CR et al., "L-carnitine enhances excretion of propionyl coenzyme A as propionylcarnitine in propionic acidemia," *J Clin Invest* 1984, DOI 10.1172/JCI111387: "Treatment with L-carnitine greatly enhanced the formation and excretion of short-chain acylcarnitines").\[9\] Population: human PA (n=3) + controls. **VERIFIED (abstract).** This diverts propiogenic flux upstream of 3-HP formation, so a patient's carnitine dose and free-carnitine status perturb the mapping between the true propiogenic load and the recorded 3-HP — effectively a time-varying contribution to both `b` and `k_e`.

### 2.7 Urinary excretion of the propionate metabolites — can a timed urine measure `k_e·m` directly?

Quantitative 24-h data are sparse and are reported as summed metabolite excretion, not per-analyte plasma-referenced clearances. Thompson GN, Chalmers RA, "Increased urinary metabolite excretion during fasting in disorders of propionate metabolism," *Pediatr Res* 1990;27(4):413–416, measured methylmalonate, methylcitrate, propionylglycine and 3-HP in 5 MMA and 3 PA children:\[4\] "the mean total measured metabolite excretion during fasting (fed 38.1 μmol/kg/h, fasting 54.6 in methylmalonic acidemia, fed 1.45, fasting 2.98 in propionic acidemia)."\[4\] Population: human PA + MMA. **VERIFIED (abstract).** Note the order-of-magnitude lower total excretion in PA than MMA, and the large fed-vs-fasting swing — evidence that urinary output tracks metabolic state, not a stable `k_e`.

Day-to-day within-patient variability is qualitatively large: Duran M et al., "The variability of metabolite excretion in propionicacidaemia," *Clin Chim Acta* 1978;82:93–99, analysed random urines from 8 PA patients\[5\] and found "wide variations" in the metabolite pattern.\[5\] **VERIFIED (abstract).**

**Could a timed urine collection measure `k_e·m` directly?** Only partially. Urinary 3-HP excretion rate divided by plasma 3-HP gives the *renal* clearance component of `k_e`, but (i) 3-HP disposal is also hepatic/peripheral (§2.2), so renal excretion is a lower bound on total non-enzymatic disposal; (ii) OAT1-mediated secretion means renal handling itself is patient- and CKD-dependent; and (iii) no study has published paired timed-urine-plus-plasma 3-HP to establish what fraction of non-enzymatic disposal is renal. That renal fraction is **NOT FOUND** in the literature and would have to be measured.

### 2.8 Where 3-HP sits in the mRNA-3927 programme

3-HP and 2-methylcitrate are the pharmacodynamic biomarkers of the mRNA-3927 (PCCA+PCCB dual mRNA) phase 1/2 trial (NCT04159103;\[10\] Interim analyses, *Nature* 2024, PMC11156579). As of the 31 May 2023 cut-off, 16 participants received mRNA-3927 across 5 dose levels (cohort 1 n=4; cohorts 2–5 n=3), totalling 346 IV doses over 15.69 person-years, per Koeberl D et al., Nature 2024;628:872–877 (PMC11156579). Moderna's Phase 1/2 interim slide deck (mRNA-3927, 24 May 2024) states verbatim: "Baseline 3-HP biomarker levels highly variable across patients… Pre-treatment 3-HP values are highly variable between patients and within patients" (each patient had four pre-treatment 3-HP blood draws), which is exactly the between-patient heterogeneity that an unshared, patient-specific `k_e` would produce. Population PK/PD modelling of the interim trial data (ScienceOpen, DOI 10.70534/JVQU8233) reports PCCA mRNA "typical clearance of 0.266 L/hr and volume of distribution of 13.9 L," with a 0.6 mg/kg biweekly dose giving meaningful 3-HP/2-MC reductions. **VERIFIED (abstract/slide deck).**

### 2.9 Haemodialysis / haemofiltration clearance of 3-HP — not reported for 3-HP; MMA analogue available

No dialysis clearance for 3-HP itself was found (**NOT FOUND**). The closest kinetic analogue is methylmalonic acid: a single haemodialysis session gave "a plasma clearance of 103 ml/min and Vd 0.48 L/kg, which approximates the total body free water space… followed by rapid reaccumulation of methylmalonic acid over 24 h to the predialysis concentration"\[11\] (Sloan JL / detailed kinetics paper, PMID 24961826, mut⁰ MMA adult). Population: human MMA (n=1). **VERIFIED (abstract).** The rapid post-dialysis rebound to baseline is the key qualitative lesson: for a small, water-soluble, rapidly-produced organic acid, extracorporeal removal does not durably change the steady state — consistent with a disposal term dominated by continuous metabolic/renal `k_e` rather than by any removable pool. The Vd ≈ total body water for MMA is a reasonable structural prior for 3-HP's own Vd (both are small polar monocarboxylic acids), but it is a cross-metabolite inference, not a 3-HP measurement.

---

## 3. Mapping to the model

The balance term is `k_e·m`, a flux (in the load units of `s`, i.e. μmol·kg⁻¹·h⁻¹) obtained by multiplying plasma 3-HP concentration `m` (μmol·L⁻¹) by `k_e`. Dimensionally `k_e` is therefore a **body-mass-normalised clearance** (L·kg⁻¹·h⁻¹) — the sum of renal secretion/filtration and non-PCC metabolic disposal of 3-HP. Each measurable quantity would translate as follows.

- **IV micro-dose of ¹³C₃-3-HP + serial plasma (preferred).** The terminal plasma decay constant λ gives total 3-HP clearance CL = λ·Vd; dividing by body mass yields `k_e` directly. In a disease-range patient (5–10% residual capacity) the PCC route carries little flux, so measured total clearance ≈ non-enzymatic `k_e` — the cleanest available estimator. Precedent: the Barshop 30 μmol/kg IV bolus with 5–90 min plasma sampling (§2.3) achieved a resolvable half-life for propionate; an analogous 3-HP schedule (denser early points, out to ~2–3 h given the slower hepatic uptake in §2.2) should resolve λ. This is the only route with a realistic path to the ~2% per-patient precision the design demands, and even then precision will be limited by Vd estimation.
- **Timed urine + paired plasma.** Urinary 3-HP excretion rate ÷ plasma 3-HP = renal clearance CL_renal, a *component* of `k_e`. Usable as a lower bound and to track within-patient change, but it omits hepatic/peripheral disposal (§2.2) and is itself CKD-sensitive (§2.4–2.5). Cannot on its own deliver total `k_e`.
- **Population prior from surrogates.** Propionate turnover CV ≈ one-third in healthy adults (§2.3) is the best-supported analogue for between-patient spread; combined with 50% adult CKD (§2.5) the true `k_e` CV in a PA cohort is plausibly larger. This supports keeping `k_e` patient-specific (a per-patient measured value or a hierarchical prior), not a shared constant, and confirms the current factor-2.2 prior is not narrow enough for θ-identifiability at 5% capacity.
- **Vd structural prior.** In the absence of a 3-HP Vd, MMA's dialysis-derived Vd ≈ 0.48 L/kg ≈ total body water (§2.9) is a defensible placeholder for converting a measured λ to CL, flagged explicitly as a cross-metabolite assumption.

**Bottom line for the design:** `k_e` cannot be fixed from the literature to 2–6%; it must be measured per patient, and the IV ¹³C₃-3-HP decay study is the only method that estimates the whole of `k_e` rather than a fraction of it.

---

## 4. Pending-reference table (format of `refs/literature_pending.md`)

| Citation | DOI / PMID | What it supplies | Status | Registered in `paper/refs.sources`? |
|---|---|---|---|---|
| Wilson KA et al., *Am J Physiol Endocrinol Metab* 2017;313:E413–E428 | 10.1152/ajpendo.00105.2017 / 28634175 | rat perfused-liver 3-HP uptake; "no human data" statement | to register | check before adding |
| Ando T, Rasmussen K, Nyhan WL, Hull D, *PNAS* 1972;69:2807–2811 | 10.1073/pnas.69.10.2807 / 4507604 | qualitative β-oxidation shunt to urinary 3-HP | to register | check before adding |
| Barshop BA et al., *Pediatr Res* 1991;30:15–22 | — / (Pediatr Res 30:15) | blood propionate t½; IV bolus schedule precedent | to register | check before adding |
| Thompson GN, Chalmers RA, *Pediatr Res* 1990;27:413–416 | — / 2352453 | 24-h urinary metabolite excretion, PA vs MMA | to register | check before adding |
| Duran M et al., *Clin Chim Acta* 1978;82:93–99 | 10.1016/0009-8981(78)90030-X | within-patient urinary variability | to register | check before adding |
| Nigam SK et al., *Physiol Rev* 2015;95:83–123 | 10.1152/physrev.00025.2013 | 3-HP as OAT1 substrate (via ref 74) | to register | check before adding |
| Eraly SA et al., *J Biol Chem* 2006;281:5072–5083 | — / 16354673 | primary Oat1-KO 3-HP elevation | to register | check before adding |
| Shchelochkov OA et al., *Genet Med* 2019;21:2830–2835 | — / 31249402 | CKD prevalence in PA | to register | check before adding |
| Roe CR et al., *J Clin Invest* 1984;73:1785–1788 | 10.1172/JCI111387 | carnitine → propionylcarnitine excretion | to register | check before adding |
| MMA HD kinetics paper, *Mol Genet Metab* 2014 | — / 24961826 | MMA dialysis clearance/Vd; rebound | to register | check before adding |
| Interim analyses mRNA-3927, *Nature* 2024 | — / PMC11156579 | 3-HP/2-MC as trial PD biomarkers; baseline variability | check | likely already registered |

---

## 5. Three findings that most change the design

1. **`k_e` has never been measured in humans, for 3-HP or by any surrogate that pins it to 2–6%.** The design's most load-bearing prior is unsupported by literature. This is a genuine gap, not a search failure — the one group that studied 3-HP disposal states outright that mammalian 3-HP metabolism is poorly characterised. The paper should stop treating `k_e` as knowable from references and treat it as a quantity to be measured.
2. **`k_e` is intrinsically patient-specific**, because 3-HP is renally secreted via OAT1\[1\] and CKD affects ~50% of adult PA patients,\[2\] with carnitine status adding a competing excretory route. A single shared `k_e` prior is therefore mis-specified; the model needs either a per-patient measured `k_e` or a hierarchical patient-level parameter, and the between-patient CV is plausibly ≥⅓ (the healthy-adult propionate-turnover analogue) and larger in a CKD-laden cohort.
3. **Only an IV ¹³C₃-3-HP plasma-decay study can measure the whole of `k_e`**; a timed urine gives only its renal fraction. The bolus IV [¹³C]propionate studies already performed safely in PA/MMA children establish the ethical and technical precedent, so the study is feasible — but the labelled 3-HP is not an off-the-shelf reagent and would need synthesis, and Vd must still be assumed (MMA's ≈ total-body-water value is the best placeholder).

## 6. What remains unknown

- 3-HP plasma half-life, systemic clearance, Vd, and fraction excreted unchanged in humans — all **NOT FOUND**.
- The renal fraction of total non-enzymatic 3-HP disposal (i.e. how much of `k_e` a timed urine would capture) — **NOT FOUND**.
- Between-patient and within-patient CV of `k_e` for 3-HP specifically — **NOT FOUND**; only propionate-turnover analogues exist.
- Haemodialysis/haemofiltration clearance of 3-HP and its post-dialysis rebound — **NOT FOUND** (only the MMA analogue is available).
- Commercial availability and radiochemical/isotopic purity of ¹³C₃-3-HP suitable for IV human microdosing — not confirmed; synthesis appears necessary.
- Whether the mRNA-3927 trial captured any paired plasma+urine or any tracer 3-HP kinetics that could retrospectively bound `k_e` — not established from public sources.

*Prepared for `refs/reviews_2026-10/01_ke_clearance.md`. Do not re-register sources already in `paper/refs.sources`; verify each row of §4 against that file before adding. `SUMMARY.md` and `README.md` are known to be out of date and were not relied upon.*

## Sources

1. [The Organic Anion Transporter (OAT) Family: A Systems Biology Perspective](https://journals.physiology.org/doi/full/10.1152/physrev.00025.2013)
2. [Inter-relations between 3-hydroxypropionate and propionate metabolism in rat liver: relevance to disorders of propionyl-CoA metabolism](https://journals.physiology.org/doi/full/10.1152/ajpendo.00105.2017)
3. [Measurement of propionate turnover in vivo using sodium \[2H5\]propionate and sodium \[13C\]propionate - ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/0009898189900739)
4. [Increased Urinary Metabolite Excretion during Fasting in Disorders of Propionate Metabolism](https://www.nature.com/articles/pr199094)
5. [The variability of metabolite excretion in propionicacidaemia - ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/000989817890030X)
6. [Metabolism of 1-13C-Propionate in Vivo in Patients with ...](https://escholarship.org/content/qt3jk0n822/qt3jk0n822_noSplash_948821ffb98b55153dd4bed342a4ffa3.pdf)
7. [3-hydroxypropionate: significance of -oxidation of propionate in patients with propionic acidemia and methylmalonic acidemia - PubMed](https://pubmed.ncbi.nlm.nih.gov/4507604/)
8. [Chronic kidney disease in propionic acidemia - PubMed - NIH](https://pubmed.ncbi.nlm.nih.gov/31249402/)
9. [JCI - L-carnitine enhances excretion of propionyl coenzyme A as propionylcarnitine in propionic acidemia.](https://www.jci.org/articles/view/111387)
10. [Interim analyses of a first-in-human phase 1/2 mRNA trial for propionic acidaemia](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11156579/)
11. [A detailed analysis of methylmalonic acid kinetics during hemodialysis and after combined liver/kidney transplantation in a patient with mut (0) methylmalonic acidemia - PubMed](https://pubmed.ncbi.nlm.nih.gov/24961826/)
