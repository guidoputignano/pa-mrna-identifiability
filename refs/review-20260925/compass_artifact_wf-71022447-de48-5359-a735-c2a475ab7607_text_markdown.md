# Ex vivo → in vivo anchor for residual PCC capacity (θ) in propionic acidemia

*Literature review 04 for `guidoputignano/system-identification-pa`; branch `reviews/2026-10-ex-vivo-anchor`; prepared 25 September 2026. Scope: can an ex vivo PCC activity serve as a usable prior on the in vivo maximal velocity θ of propionyl-CoA carboxylase, and how uncertain is the ex-vivo→in-vivo mapping? Also: does a same-day repeatability figure exist for the 13C-propionate breath test in PA?*

## TL;DR

- **The ex-vivo→in-vivo mapping uncertainty cannot be estimated as a defensible log-scale SD from the published record.** No study measures hepatic PCC Vmax against the same patient's fibroblast or lymphocyte PCC activity, so the quantity the therapy targets (θ = hepatic PCC capacity) has never been calibrated against the ex vivo assay. The single quantitative anchor that exists is indirect and from MMA, not PA: whole-body 1-13C-propionate oxidation versus fibroblast 14C-propionate incorporation correlated at Spearman r = 0.682 (R² = 0.327, n = 11),\[1\] i.e. roughly two-thirds of the variance is unexplained — an informal lower bound implying at least a several-fold (≈0.3–0.5 log10, i.e. ~0.7–1.1 natural-log SD) residual scatter that conflates assay, tissue and disease.
- **No same-day (within-day, test–retest) repeatability figure for the 13C-propionate breath test in PA exists.** The best available repeatability data are from the MMA cohort (Manoli 2021): intrasubject CV 6.9% in controls and 12.8% in patients\[1\] across visits weeks-to-years apart, with near-perfect analytic agreement between IRMS and the BreathID device (r = 0.996). The PA severity study (Shchelochkov 2021) reports only "good intrasubject reproducibility"\[2\] in a supplemental figure without a numeric same-day CV.
- **An ex vivo activity is nonetheless a legitimate, and probably the only practical, external anchor on θ**, because in the disease range the plasma-marker record cannot bound θ from below and residual activity does track severity, though the difference between mild and severe phenotypes may be only a few percentage points, with very mild forms retaining approximately 3% of WT at the lower end (Totsune et al. 2026, Mol Ther Nucleic Acids); but any prior built from it must carry a wide, largely unquantified variance and should be treated as order-of-magnitude, not calibrated.

## Key findings

1. **Residual PCC activity in PA cells is low, assay-dependent, and only weakly tied to genotype.** Across the classic and modern assays, patient fibroblast/lymphocyte PCC activity clusters below ~5% of control for classic disease and in the ~1–8% band for milder/late-onset disease, but the boundary between "mild" and "severe" can be a few percentage points, and no genotype–activity law is established.
2. **The record contains a quantitative ex-vivo↔in-vivo link only in MMA, and only for intact-cell flux vs whole-body oxidation** (Spearman r = 0.682). There is no equivalent PA figure, and nothing that maps the *enzyme Vmax* (θ) specifically.
3. **Hepatic-vs-cell proportionality for PCC is unproven.** No paired liver-and-cell PCC activity dataset exists; the only paired same-patient tissue/cell enzyme data in the pathway are for MMA mutase (Morrow 1975) and are qualitative (cofactor-responsive yes/no), not a proportional correlation.
4. **The breath test is reproducible enough to be a pharmacodynamic readout, but its PA repeatability is characterised only indirectly** through the MMA method-development data and a qualitative PA supplemental figure.

## Details — sourced evidence

### Question 1 — Distribution of residual PCC activity in PA cells (assay, reproducibility, genotype, biotin, temperature)

**Hsia, Scully & Rosenberg 1979 — VERIFIED.** *Human propionyl CoA carboxylase: some properties of the partially purified enzyme in fibroblasts from controls and patients with propionic acidemia.* Pediatr Res 13:746–751. PMID 481943.
- Species/tissue/population: human cultured skin fibroblasts; 7 PA patient lines + controls. Assay: partial purification + carboxylase assay.
- Verbatim: "Each of seven lines from patients with propionic acidemia had clearly detectable PCC activity which was less than 5% of that in control lines."\[3\]\[4\]
- Verbatim (temperature lability): "Mutant enzyme had a pH optimum, ionic requirements, and substrate Km's similar to those of control PCC, but was distinctly more labile to both cold and heat."\[3\]\[4\]
- Relevance: sets the classic "<5% of control" figure and shows mutant PCC is thermolabile — a warning that ex vivo handling temperature can move the measured activity.

**Wolf & Rosenberg 1978 — VERIFIED.** *Heterozygote expression in propionyl coenzyme A carboxylase deficiency. Differences between major complementation groups.* J Clin Invest 62:931–936. PMID 711858; PMCID PMC371850.
- Species/tissue/population: human skin fibroblasts and peripheral blood leukocytes; obligate heterozygotes (6 pccA, 12 pccC) + controls.
- Verbatim: "Mean PCC activity in fibroblast extracts from pcc A heterozygotes was 52% of that in controls, whereas mean PCC activity in pcc C heterozygotes was indistinguishable from that of controls. Similar results were obtained with extracts of peripheral blood leukocytes."\[5\]\[6\]
- Relevance: the one paper showing fibroblast and leukocyte PCC activity track each other (cell-to-cell consistency), though it involves no liver and gives no correlation coefficient.

**Wolf & Rosenberg 1979 (biotin) — VERIFIED.** *Stimulation of Propionyl CoA and β-Methylcrotonyl CoA Carboxylase Activities in Human Leukocytes and Cultured Fibroblasts by Biotin.* Pediatr Res 13:1275–1279.
- Species/tissue/population: human leukocytes (incl. a biotin-supplemented PCC-deficient infant) and control fibroblasts.
- Verbatim: "PCC and βMCC activities were increased 70–75% and 120–140%, respectively, when confluent normal fibroblasts were incubated in medium containing 10% human serum and 0.1–1.0 μg/ml biotin."\[7\]
- Verbatim: "whereas the propionyl CoA carboxylase (PCC) activity in peripheral blood leukocytes of a 13-month-old girl with isolated PCC deficiency who had been given biotin supplements (4.5 mg/day) for 6 months was only 3% of control, her βMCC activity was greatly elevated at 20 times control activity."\[7\]
- Relevance: culture biotin can shift measured PCC by tens of percent even in normal cells — a second, controllable source of ex vivo variance.

**Kobayashi et al. 2017 (UPLC-MS/MS lymphocyte assay) — VERIFIED.** *Determination of methylmalonyl coenzyme A by ultra high-performance liquid chromatography tandem mass spectrometry for measuring propionyl coenzyme A carboxylase activity in patients with propionic acidemia.* J Chromatogr B. (ScienceDirect S1570023216306456; PMID 28189105.)
- Species/tissue/population: human peripheral lymphocytes; 20 healthy individuals + 6 PA patients. Assay: enzyme reaction → UPLC-MS/MS quantitation of methylmalonyl-CoA.
- Verbatim: "Methylmalonyl-CoA was not detected in one PA patient with a severe form of the disease, but the remaining PA patients with mild disease showed residual activities (3.3–7.8%)."\[8\]
- Verbatim (assay reproducibility): intra-assay CVs "8.2%, 8.7%, and 5.1%, respectively, and the inter-assay CVs were 13.6%, 10.5%, and 5.9%, respectively."\[9\]
- Relevance: modern PA-specific assay; gives both the mild-PA residual band (3.3–7.8%) and the analytic CV of the assay itself (~5–14%).

**Yorifuji et al. / Current Perspectives (Japan, PCCB mild variant) — VERIFIED (secondary quoting Yorifuji primary).** *Current Perspectives on Neonatal Screening for Propionic Acidemia in Japan.* Int J Neonatal Screen 2021;7(3):35. PMCID PMC8293189.
- Population: Japanese newborns with the common mild PCCB variant (Y435C).
- Verbatim: "it has been difficult to conclude them as false positive subjects, for PCC activities in their lymphocytes were reported to be as low as 3–8% of normal control value."\[10\]
- Relevance: even asymptomatic/very mild PA sits at 3–8% lymphocyte activity — i.e. the informative dynamic range of the ex vivo assay across the whole clinical spectrum is compressed into single-digit percentages.

**Totsune et al. 2026 (PCCA pseudoexon ASO) — VERIFIED.** *From N-of-1 to versatility in propionic acidemia: Antisense oligonucleotide-mediated skipping of a constitutive PCCA pseudoexon.* Totsune E, Maeda Y, Kikuchi A, et al. Mol Ther Nucleic Acids 2026;37(2). PMCID PMC13101689.
- Species/tissue/population: human PA fibroblasts (index case A1 + 7 additional PA lines). Assay: UHPLC-MS/MS PCC activity.
- Verbatim: "PCC activity in A1 fibroblasts also profoundly decreased to 0.5% of that in WT fibroblasts."\[11\]
- Verbatim (genotype): "No other relationship was found between PCC activity and genotype."\[11\]
- Verbatim (severity granularity): "the difference between mild and severe phenotypes may be only a few percentage points… very mild forms of PA… retain PCC activity of approximately 3% of WT at the lower end of the\[11\] range… Conversely, a late-onset PA case with PCC activity of 2.24% of WT… has been reported to develop severe cardiomyopathy."
- Relevance: directly quantifies how compressed and noisy the activity–phenotype relation is (a 2.24%-WT case being severe while ~3%-WT cases are near-asymptomatic), and confirms no usable genotype→activity law.

**GeneReviews (Shchelochkov/Carrillo/Venditti) — VERIFIED (secondary/authoritative synthesis).** *Propionic Acidemia.* NCBI Bookshelf NBK92946.
- Verbatim: "no genotype-phenotype correlation is known [Pena et al 2012, Forny et al 2021]."\[12\]
- Verbatim: "The homozygous PCCB c.1606A>G (p.Asn536Asp) pathogenic variant in Amish and Mennonite populations is associated with high residual propionyl-coenzyme A carboxylase (PCC) activity."\[12\]
- Relevance: the field's consensus that genotype does not pin down activity, and the one named high-residual-activity variant.

### Question 2 — Quantitative ex vivo ↔ in vivo relation

**Manoli et al. 2021 (Part 1, MMA) — VERIFIED.** *1-13C-propionate breath testing as a surrogate endpoint to assess efficacy of liver-directed therapies in methylmalonic acidemia (MMA).* Genet Med 23:1522–1533. DOI 10.1038/s41436-021-01143-8; PMID 33820958; PMCID PMC8354855.
- Species/tissue/population: 57 MMA subjects (incl. 19 transplant recipients) + 16 controls;\[13\] separately, diagnostic skin-fibroblast 14C-propionate incorporation available for 11.
- Verbatim (the key ex-vivo↔in-vivo number): "The percent of control 14C-propionate incorporation correlated with the 60-minute 1-13C-propionate propionate oxidation (Spearman r = 0.682, R2 = 0.327, P =\[1\] 0.025)."\[14\]
- Relevance: the ONLY published quantitative cell-vs-in-vivo correlation in a propionate disorder. It is MMA, and it compares intact-cell flux (14C incorporation) — not enzyme Vmax — to whole-body oxidation. R² = 0.327 ⇒ ~67% of variance unexplained.

**Barshop et al. 1991 (PA, in vivo vs lymphocyte enzyme) — VERIFIED.** *Metabolism of 1-13C-propionate in vivo in patients with disorders of propionate metabolism.* Pediatr Res 30:15–22. DOI 10.1203/00006450-199107000-00004; PMID 1909779.
- Verbatim (in vivo vs clinical severity, PA): "For PA, slope m = −2.57, regression coefficient r = 0.849, and significance p = 0.00047; for MMA, m = 0.033, r = 0.015, and p = 0.89."
- Verbatim (lymphocyte PCC in PA, % control of PCC/MCC ratio): patient P.S. 0.08%, E.A. 0.28%, L.S. 0.86%, E.N. 0.22% (Table 2, "% Control").
- Relevance: in PA, in vivo oxidation tracks clinical severity strongly (r = 0.849) but the paper does NOT regress in vivo oxidation on the lymphocyte enzyme activity — so it gives no ex-vivo→in-vivo slope for PA.

**Shchelochkov et al. 2021 (Part 2, PA severity modeling) — VERIFIED.** *Severity modeling of propionic acidemia using clinical and laboratory biomarkers.* Genet Med 23:1534–1542. DOI 10.1038/s41436-021-01173-2; PMID 34007002; PMCID PMC8354856.
- Species/tissue/population: 40 PA patients (5 transplanted). In vivo 1-13C-propionate oxidation.
- Verbatim: "in vivo 1-13C-propionate oxidation, ranked 20/301, clearly separated liver-transplanted participants from severely affected PA patients, and did not correlate with cystatin C–based eGFR."\[2\]
- Verbatim: "This suggests that in vivo propionate oxidation may be linked to the residual activity of propionyl-CoA carboxylase activity."\[2\]
- Relevance: establishes in vivo oxidation as the best-separating PA biomarker and explicitly hypothesises (does not measure) its link to residual PCC activity.

### Question 3 — In vivo 13C-propionate oxidation in PA specifically; dose, sampling, repeatability

**Barshop 1991 — VERIFIED (PA-specific values).**
- Dose/route: "administered i.v. at 30 pmol/kg [30 µmol/kg] or orally at 100 Frnol/kg [100 µmol/kg] body weight";\[15\] breath sampled to 180 min; recovery normalised to predicted resting energy expenditure (no indirect calorimeter).
- Verbatim (PA phenotype of the signal): "There was very little generation of 13CO2 in patients with PA, and the curves were virtually flat."
- Verbatim (biotin within-patient change, milder PA patient L.S.): "the conversion of 1-13C-propionate to 13CO2 was increased in L.S. during biotin treatment (from 21.7 to 30.5% at 3 h)."
- Repeatability: patients were re-tested on/off biotin or carnitine but NOT in a same-day duplicate design; no CV reported. Newer than Barshop: Manoli 2021 and Shchelochkov 2021 (above) are the current PA/MMA references.

**Manoli 2021 — VERIFIED (method, dose, repeatability).**
- Dose/route: "A dose of 0.5 mg/kg body weight… administered orally or through a gastrostomy tube… over 2 minutes"\[1\] (a 20-fold reduction from Barshop's 9.7 mg/kg); VCO2 measured by metabolic cart; δ13C by IRMS and by Exalenz BreathID.
- Verbatim (repeatability): "At 120 minutes, the coefficient of variation (CV) for dose metabolized within each subject ranged from 1.83% to 13.33%, with an average of 6.94 ± 5.3%"\[1\] (7 healthy volunteers, repeated 3× at **1- to 3-week** intervals).
- Verbatim (patients): "the CV for 1-13C-propionate dose metabolized in 15 mut patients, who had repeated breath testing 2–3 times each, ranged from 0.39% to 50.84% with an average of 14.02 ± 12.86%"\[1\] (intervals **2 months to 4 years**).
- Verbatim (analytic agreement): IRMS vs BreathID "correlation coefficient was r = 0.996 (P < 0.0001)… 95% limits of agreement ranging between −2.386% and 3.266%."
- Relevance: this is the source of the "6.9% / 12.8%" repeatability figures. They are across-visit, not same-day, and in MMA.

### Question 4 — Is hepatic PCC activity proportional to fibroblast/lymphocyte activity?

**Answer: NOT FOUND — no numeric paired liver-vs-cell PCC dataset or correlation exists.**

**Kalousek, Darigo & Rosenberg 1980 — VERIFIED (liver enzyme characterised, but not paired).** *Isolation and characterization of propionyl-CoA carboxylase from normal human liver.* J Biol Chem 255:60–65. PMID 6765947.
- Species/tissue: normal postmortem human liver; ~3000-fold purification, 26% yield.\[16\] Km ATP 0.08 mM, propionyl-CoA 0.29 mM, bicarbonate 3.0 mM.\[16\]
- Relevance: defines the hepatic enzyme's kinetics but provides no fibroblast/leukocyte comparison and no tissue-distribution activity table.

**Kalousek, Orsulak & Rosenberg 1983 — VERIFIED (paired material, but immunochemical not quantitative).** *Absence of cross-reacting material in isolated propionyl CoA carboxylase deficiency: nature of residual carboxylating activity.* Am J Hum Genet 35:409–420. PMID 6859037; PMCID PMC1685634.
- Species/tissue/population: PA-patient fibroblast extracts and fetal liver homogenates (pccA, pccC, one bio) + controls. Endpoint: immunotitration + protein-concentration dependence.
- Verbatim: "When propionyl CoA carboxylation was assayed at high protein concentration in a fetal liver homogenate from a pcc C patient, the apparent PCC activity was comparable to that found in normal human fetal liver."\[17\]
- Verbatim: "the residual propionyl CoA carboxylating activity found in isolated PCC deficiency represents another carboxylase that can utilize propionyl CoA as a substrate rather than a mutant form of PCC."\[18\]
- Relevance: a caution — a large part of measured "residual PCC activity" in deficient extracts is not PCC at all, and is protein-concentration-dependent, further degrading the ex vivo assay as an estimator of true θ.

**McKeon, Eanes, Fall, Tasset & Wolf 1980 — SECONDARY (full text not opened).** *Immunological studies of propionyl CoA carboxylase in livers and fibroblasts of patients with propionic acidemia.* Clin Chim Acta 101:217–233. PMID 6766827.
- Studied both liver and fibroblasts from PA patients but for CRM/Km, not for a quantitative liver-vs-cell activity comparison; no correlation coefficient reported in any citing source. Flagged as the single most likely place a numeric paired PCC table could exist; not verifiable here.

**Morrow, Mahoney, Mathews & Lebowitz 1975 — VERIFIED (paired, but MMA mutase and qualitative).** *Studies of Methylmalonyl Coenzyme A Carbonylmutase Activity in Methylmalonic Acidemia. I.* Pediatr Res 9:641–644. DOI 10.1203/00006450-197508000-00006.
- Species/tissue/population: human; fibroblasts from 15 MMA patients, postmortem liver/kidney from 6 of the same patients.\[19\]
- Verbatim: "In the liver or kidney extracts of postmortem tissues, mutase activity responded to AdoCbl supplementation if fibroblast mutase activity from that patient had responded, and failed to respond if fibroblast activity failed to respond."\[19\]
- Relevance: the closest existing same-patient tissue-vs-cell enzyme dataset in the pathway; concordance is qualitative (cofactor-responsive yes/no), for MMA mutase, not PCC, and yields no proportionality constant.

## Synthesis — best estimate of the ex-vivo→in-vivo mapping uncertainty

**A defensible numeric log-scale SD for the fibroblast/lymphocyte-PCC → hepatic-PCC-Vmax (θ) mapping cannot be extracted from the literature, because the mapping has never been measured.** The reasoning:

1. **The target quantity is unmeasured against the anchor.** θ is hepatic PCC maximal velocity. No study reports hepatic PCC activity paired with the same patient's cell PCC activity (Question 4, NOT FOUND). The only paired tissue/cell enzyme data (Morrow 1975, MMA mutase) are qualitative.
2. **The one quantitative cell-vs-in-vivo correlation is a proxy, not the mapping.** Manoli 2021 gives Spearman r = 0.682 (R² = 0.327)\[1\] between fibroblast 14C-propionate *incorporation* and whole-body oxidation, in MMA. Interpreting R² = 0.327 as the fraction of variance explained leaves ~0.67 unexplained. If (and only if) one treats the two measures as linearly related on the log scale with that R², the residual SD is √0.67 ≈ 0.82 of the total spread — heuristically a **~0.3–0.5 log10 (≈0.7–1.1 natural-log) residual SD**, i.e. a several-fold uncertainty. This is a SECONDARY, order-of-magnitude inference: it conflates (a) intact-cell flux vs enzyme Vmax, (b) fibroblast vs hepatocyte, and (c) MMA vs PA, each of which adds variance. It should therefore be read as a **lower bound** on the true θ-mapping uncertainty, not an estimate of it.
3. **Assay and pre-analytic variance alone are already large.** The PCC assay's own analytic CV is ~5–14% (Kobayashi 2017); measured activity moves by tens of percent with culture biotin (Wolf & Rosenberg 1979) and is thermolabile in mutant cells (Hsia 1979); and a substantial fraction of low-level "residual PCC activity" is a different carboxylase entirely (Kalousek 1983). These stack on top of the tissue-mapping uncertainty.
4. **The informative range is compressed.** Across the whole PA clinical spectrum, cell activity spans roughly 0–8% of control, with severe vs near-asymptomatic separated by "only a few percentage points" and at least one late-onset case at 2.24% of WT developing severe cardiomyopathy (Totsune et al. 2026, Mol Ther Nucleic Acids; PMC13101689). A prior on θ built from an ex vivo value near the assay floor is therefore intrinsically wide.

**Practical recommendation for the model.** Use an ex vivo PCC activity as a weakly-informative, order-of-magnitude prior on θ — it is the only external anchor available besides k_e and it does bound θ where the marker record cannot — but assign it a deliberately wide log-scale prior. Given the evidence above, a **natural-log SD of order ≥1 (i.e. roughly a 3-fold 1-σ band, wider at the tails)** is the most honest choice, and the report should state plainly that this figure is an assumption informed by an MMA proxy rather than a measured PA cell→liver calibration. Do not treat the ex vivo value as a tight prior; sensitivity analysis over the prior SD is warranted.

## Same-day repeatability of the 13C-propionate breath test in PA

**No same-day (within-day) repeatability figure for the PA breath test exists in the literature (NOT FOUND).** What exists:
- **MMA, across-visit (Manoli 2021, VERIFIED):** intrasubject CV 6.94 ± 5.3% in controls (repeats 1–3 weeks apart) and 14.02 ± 12.86% in patients (repeats 2 months–4 years apart); abstract rounds these to 6.9% and 12.8%. Analytic method agreement IRMS-vs-BreathID r = 0.996, limits of agreement −2.4% to +3.3%.
- **PA (Shchelochkov 2021, VERIFIED):** "good intrasubject reproducibility shown in Supplemental Fig. 3a–c"\[2\] — qualitative only; no numeric same-day CV in the main text.
- **PA (Barshop 1991, VERIFIED):** repeat tests were done on/off biotin/carnitine, not as same-day duplicates; no CV.

So the best available surrogate for expected PA test–retest noise is the MMA figure of ~7% (controls) to ~13% (patients), acknowledging it is across-visit and in a different disease. A dedicated same-day PA repeatability study is a genuine gap.

## What remains unknown (caveats)

- **No cell→liver PCC calibration in PA (or any human).** The proportionality constant and its scatter between hepatic θ and ex vivo PCC activity are simply unmeasured. McKeon 1980 (Clin Chim Acta 101:217–233) is the single most likely place a numeric paired PCC table might exist and should be retrieved in full before concluding definitively.
- **The Manoli r = 0.682 anchor is MMA, intact-cell flux, and n = 11** — three reasons it may understate the true PA enzyme-Vmax mapping uncertainty.
- **No same-day PA breath-test CV**, so pharmacodynamic change thresholds under mRNA-3927 — ModernaTX's dual LNP-encapsulated mRNA encoding the PCCA and PCCB subunits, evaluated in the first-in-human Phase 1/2 open-label dose-optimization study NCT04159103 (Nature 2024, DOI 10.1038/s41586-024-07266-7; 346 intravenous doses over 15.69 person-years, with a 70% reduction in the risk of metabolic decompensation events among 8 participants who reported them in the 12-month pretreatment period) — currently borrow MMA across-visit repeatability.
- **Genotype cannot substitute for the assay** (no genotype–activity law; GeneReviews, Totsune et al. 2026).
- **The ex vivo assay is itself noisy and partly non-specific** (biotin, temperature, protein concentration, non-PCC carboxylase), so even a perfect tissue mapping would not make it a precise anchor.

## Pending-reference table (format of `refs/literature_pending.md`)

| bibkey | citation | DOI / PMID | maps to repo ref | relevance to θ anchor | label |
|---|---|---|---|---|---|
| hsia_1979_fibroblast_pcc | Hsia YE, Scully KJ, Rosenberg LE. Human propionyl CoA carboxylase: properties of the partially purified enzyme in fibroblasts from controls and PA patients. Pediatr Res 1979;13:746–751. | PMID 481943 | pcc_residual_activity | PA fibroblast PCC <5% control; mutant enzyme cold/heat labile | VERIFIED |
| wolf_1978_heterozygote_pcc | Wolf B, Rosenberg LE. Heterozygote expression in PCC deficiency. J Clin Invest 1978;62:931–936. | PMID 711858 / PMC371850 | pcc_residual_activity | fibroblast≈leukocyte PCC (52% in pccA heterozygotes) | VERIFIED |
| wolf_1979_biotin_pcc | Wolf B, Rosenberg LE. Stimulation of PCC and βMCC by biotin. Pediatr Res 1979;13:1275–1279. | (Nature pr1979460) | pcc_residual_activity | biotin raises measured PCC 70–75% in normal fibroblasts | VERIFIED |
| kobayashi_2017_uplc_pcc | Determination of methylmalonyl-CoA by UPLC-MS/MS for measuring PCC activity in PA. J Chromatogr B 2017. | PMID 28189105 | pcc_residual_activity | mild-PA lymphocyte residual 3.3–7.8%; assay CV 5–14% | VERIFIED |
| yorifuji_japan_mild_pccb | Current Perspectives on Neonatal Screening for PA in Japan. Int J Neonatal Screen 2021;7:35. | PMC8293189 | pcc_residual_activity | very-mild PCCB variant lymphocyte activity 3–8% control | VERIFIED |
| totsune_2026_aso_pcca | Totsune E, Maeda Y, Kikuchi A, et al. From N-of-1 to versatility in PA: ASO skipping of a PCCA pseudoexon. Mol Ther Nucleic Acids 2026;37(2). | PMC13101689 | pcc_residual_activity | activity 0.5–~3% WT; no genotype–activity law; 2.24%-WT severe case | VERIFIED |
| barshop_1991_13c_pa | Barshop BA, et al. Metabolism of 1-13C-propionate in vivo in disorders of propionate metabolism. Pediatr Res 1991;30:15–22. | PMID 1909779 | breath_test_mma | PA in vivo oxidation vs severity r=0.849; dose 30 µmol/kg IV / 100 µmol/kg PO; lymphocyte PCC 0.08–0.86% control | VERIFIED |
| manoli_2021_breath_mma | Manoli I, et al. 1-13C-propionate breath testing… MMA. Genet Med 2021;23:1522–1533. | 10.1038/s41436-021-01143-8 / PMID 33820958 | breath_test_mma | fibroblast 14C incorporation vs breath r=0.682; repeatability CV 6.9%/12.8%; dose 0.5 mg/kg | VERIFIED |
| shchelochkov_2021_pa_severity | Shchelochkov OA, et al. Severity modeling of PA. Genet Med 2021;23:1534–1542. | 10.1038/s41436-021-01173-2 / PMID 34007002 | breath_test_mma | PA in vivo oxidation best-separating biomarker; hypothesised link to residual PCC | VERIFIED |
| kalousek_1980_liver_pcc | Kalousek F, Darigo MD, Rosenberg LE. Isolation and characterization of PCC from normal human liver. J Biol Chem 1980;255:60–65. | PMID 6765947 | pcc_liver_kinetics | hepatic PCC kinetics (Km propionyl-CoA 0.29 mM); no cell pairing | VERIFIED |
| kalousek_1983_crm | Kalousek F, Orsulak MD, Rosenberg LE. Absence of cross-reacting material in isolated PCC deficiency. Am J Hum Genet 1983;35:409–420. | PMID 6859037 / PMC1685634 | pcc_liver_kinetics | residual "PCC" activity is largely a different carboxylase; protein-conc dependent | VERIFIED |
| mckeon_1980_liver_fibroblast | McKeon C, et al. Immunological studies of PCC in livers and fibroblasts of PA patients. Clin Chim Acta 1980;101:217–233. | PMID 6766827 | pcc_liver_kinetics | paired liver+fibroblast material but CRM/Km only; possible paired activity table — needs full text | SECONDARY |
| morrow_1975_mutase_paired | Morrow G, et al. Methylmalonyl-CoA mutase activity in MMA: clinical, hepatic, fibroblast correlation. Pediatr Res 1975;9:641–644. | 10.1203/00006450-197508000-00006 | (new) | only paired same-patient tissue/cell enzyme data in pathway; MMA mutase; qualitative concordance | VERIFIED |

## Process note

Work performed on branch `reviews/2026-10-ex-vivo-anchor`; this document is the intended contents of `refs/reviews_2026-10/04_ex_vivo_anchor.md`. No files under `paper/`, `src/`, or `params.py` were edited and no constant was moved, per the task rules. New pending references above (notably `morrow_1975_mutase_paired`, `kalousek_1983_crm`, `manoli_2021_breath_mma`, `shchelochkov_2021_pa_severity`) are proposed for registration alongside the already-registered `pcc_residual_activity`, `pcc_liver_kinetics`, and `breath_test_mma`. Open the PR and stop.

## Sources

1. [1-13C-propionate breath testing as a surrogate endpoint to assess efficacy of liver-directed therapies in methylmalonic acidemia (MMA)](https://www.nature.com/articles/s41436-021-01143-8)
2. [Severity modeling of propionic acidemia using clinical and laboratory biomarkers - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8354856/)
3. [Human propionyl CoA carboxylase: some properties of the partially purified enzyme in fibroblasts from controls and patients with propionic acidemia - PubMed](https://pubmed.ncbi.nlm.nih.gov/481943/)
4. [Human Propionyl CoA Carboxylase: Some Properties of the Partially Purified Enzyme in Fibroblasts from Controls and Patients with Propionic Acidemia](https://www.nature.com/articles/pr1979206)
5. [Heterozygote expression in propionyl coenzyme A carboxylase deficiency. Differences between major complementation groups - PubMed](https://pubmed.ncbi.nlm.nih.gov/711858)
6. [pmc.ncbi.nlm.nih.gov](https://pmc.ncbi.nlm.nih.gov/articles/PMC371850)
7. [Stimulation of Propionyl CoA and β-Methylcrotonyl CoA Carboxylase Activities in Human Leukocytes and Cultured Fibroblasts by Biotin](https://www.nature.com/articles/pr1979460)
8. [Determination of methylmalonyl coenzyme A by ultra high-performance liquid chromatography tandem mass spectrometry for measuring propionyl coenzyme A carboxylase activity in patients with propionic acidemia - ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S1570023216306456)
9. [Determination of methylmalonyl coenzyme A by ultra high-performance liquid chromatography tandem mass spectrometry for measuring propionyl coenzyme A carboxylase activity in patients with propionic acidemia - PubMed](https://pubmed.ncbi.nlm.nih.gov/28189105)
10. [Current Perspectives on Neonatal Screening for Propionic Acidemia in Japan: An Unexpectedly High Incidence of Patients with Mild Disease Caused by a Common PCCB Variant - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8293189/)
11. [From N-of-1 to versatility in propionic acidemia: Antisense oligonucleotide-mediated skipping of a constitutive PCCA pseudoexon - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC13101689/)
12. [Propionic Acidemia - GeneReviews® - NCBI Bookshelf - NIH](https://www.ncbi.nlm.nih.gov/books/NBK92946/)
13. [1-13C-propionate breath testing as a surrogate endpoint to assess efficacy of liver-directed therapies in methylmalonic acidemia (MMA)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8354855/)
14. [1-13C-propionate breath testing as a surrogate endpoint to assess efficacy of liver-directed therapies in methylmalonic acidemia (MMA) - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8354855/)
15. <https://escholarship.org/content/qt3jk0n822/qt3jk0n822_noSplash_948821ffb98b55153dd4bed342a4ffa3.pdf>
16. [Isolation and characterization of propionyl-CoA carboxylase from normal human liver. Evidence for a protomeric tetramer of nonidentical subunits - PubMed](https://pubmed.ncbi.nlm.nih.gov/6765947)
17. [Absence of cross-reacting material in isolated propionyl CoA carboxylase deficiency: nature of residual carboxylating activity - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC1685634/)
18. [pmc.ncbi.nlm.nih.gov](https://pmc.ncbi.nlm.nih.gov/articles/PMC1685634/?page=12)
19. [Studies of Methylmalonyl Coenzyme A Carbonylmutase Activity in Methylmalonic Acidemia. I. Correlation of Clinical, Hepatic, and Fibroblast Data](https://www.nature.com/articles/pr197585)
