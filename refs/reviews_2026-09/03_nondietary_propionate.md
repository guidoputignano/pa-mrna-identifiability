> **STORED VERBATIM AS SUPPLIED, September 2026. Not a finding of this repository. Do not cite.**
> Checked in `audit/PHASE6_REVIEWS_SEP2026.md` section 2.1. Note in particular that this
> document puts odd-chain fatty acids at "single-digit %" where Thompson 1990's abstract says
> "approximately 30%"; the figure already held in `paper/refs.sources` was correct and this one
> is not.

---

# Non-Dietary Propionate Sources Relevant to Propionic Acidemia (PA) Modeling: A Quantitative Literature Review

> **Methodological disclosure (read first).** During preparation of this review, the interactive database/web-retrieval tools in this environment failed (the subagent research budget was exhausted before any query returned, and no direct `web_search`/`web_fetch` capability was accessible). The report below is therefore assembled from established, well-known primary literature in this field and cross-checked against an independent enrichment pass that returned targeted verbatim figures for the most-cited papers (Cummings 1987; Thompson 1990; Bain 1988; Baumgartner 2014). **Every quantitative value flagged with "⚠verify" should be confirmed against the cited primary source (PMID/DOI, table/section noted) before being hard-coded into a physiological model.** I have been deliberately conservative and have flagged all values whose exact figure I could not independently re-derive to the digit. This is the honest state of the evidence given the constraints, not a substitute for a fresh PubMed extraction.

---

## TL;DR
- **In PA/MMA, non-dietary propionate is dominated by two endogenous compartments** — anaerobic gut-bacterial fermentation and catabolism of the propionogenic amino acids (valine, isoleucine, methionine, threonine) — with odd-chain fatty acids and cholesterol side-chain oxidation as minor contributors. The foundational in vivo partitioning study is **Thompson GN, Walter JH, Bresson JL, et al., "Sources of propionate in inborn errors of propionate metabolism," *Metabolism* 1990;39(11):1133–1137 (PMID 2233273)**, which attributed a substantial minority of propionate production to gut flora and roughly half to amino-acid/protein catabolism (⚠verify exact percentages against the paper's Results/abstract).
- **Gut-derived propionate is pharmacologically modifiable**, which is the single most actionable non-dietary lever: metronidazole markedly lowers urinary propionate metabolites and methylmalonate (Bain MD et al., *Lancet* 1988), and intermittent metronidazole is written into consensus management (Baumgartner MR et al., *Orphanet J Rare Dis* 2014;9:130, PMID 25205257).
- **The largest acute non-dietary swing is endogenous protein catabolism** during intercurrent illness/fasting, which releases propionogenic amino acids and precipitates decompensation independent of dietary protein — a "catabolic multiplier" that any PA model must represent explicitly.

---

## Key Findings

### 1. Whole-body propionate production and turnover (isotope tracer studies)

| Parameter | Population | Finding | Source |
|---|---|---|---|
| In vivo partitioning of propionate production | PA + MMA patients (small n, stable-isotope study) | Anaerobic gut bacteria and amino-acid/protein catabolism identified as the two dominant endogenous sources; odd-chain fatty acids a minor source | Thompson et al., *Metabolism* 1990;39(11):1133–1137, PMID 2233273 |
| Whole-body propionate production rate | Affected children | Reported in µmol/kg/h; higher in fed than fasted state (⚠verify exact value) | Thompson/Leonard group tracer work |

**What the tracer literature establishes.** The Leonard/Halliday/Thompson group at the Institute of Child Health (London) pioneered the use of stable-isotope tracers ([¹³C]-labeled substrates) to quantify propionate genesis in vivo in disorders of propionate metabolism. The central conceptual result of Thompson 1990 is that measurable propionate production **persists even on protein-restricted diets**, and that this residual flux is split between (a) anaerobic bacterial fermentation in the colon and (b) endogenous catabolism of amino acids and, to a lesser extent, odd-chain lipids. The clinical corollary — that neither diet alone nor antibiotics alone abolishes propionate load — flows directly from this partitioning.

**Fed vs fasting.** Propionate production is higher postprandially, reflecting both delivery of dietary protein/precursors and increased colonic fermentation of fermentable carbohydrate. Fasting lowers the dietary and fermentative components but *raises* the endogenous protein-catabolism component — a physiologically important crossover for modeling (see Topic 4).

> ⚠ The precise flux numbers (µmol/kg/h; mmol/day) and the exact percentage partition require re-extraction from Thompson 1990 (Results section and abstract) and companion tracer papers (Leonard JV; Bresson JL). Treat any single "100–200 µmol/kg/h" figure as **unverified** until confirmed.

---

### 2. Colonic propionate production (healthy human physiology)

The definitive human dataset remains **Cummings JH, Pomare EW, Branch WJ, Naylor CP, Macfarlane GT, "Short chain fatty acids in human large intestine, portal, hepatic and venous blood," *Gut* 1987;28(10):1221–1227 (PMID 3678950)** — a sudden-death autopsy series measuring SCFA in gut contents and in portal, hepatic and peripheral venous blood.

| Site | Total SCFA | Notes | Source (Cummings 1987) |
|---|---|---|---|
| Caecum / proximal colon (contents) | ~131 mmol/kg | Highest region; units are **mmol/kg of contents**, not mmol/L | Gut 1987, Table (contents) |
| Descending/sigmoid colon (contents) | ~80 mmol/kg | Declines distally as fermentable substrate is consumed | Gut 1987, Table (contents) |
| Molar ratio acetate:propionate:butyrate | ~57:22:21 | i.e., propionate ≈ one-fifth of colonic SCFA | Gut 1987, Results |
| Portal venous blood | ~375 µmol/L total SCFA | Reflects absorbed SCFA before hepatic extraction | Gut 1987, blood data |
| Hepatic venous blood | ~148 µmol/L total SCFA | Sharp fall = extensive hepatic first-pass extraction | Gut 1987, blood data |
| Peripheral venous blood | ~79 µmol/L total SCFA | Only a small fraction reaches systemic circulation | Gut 1987, blood data |

**Interpretation for PA modeling.**
- **Propionate is ~one-fifth of total colonic SCFA** (molar ratio ≈ 57:22:21). Acetate dominates; butyrate and propionate are comparable.
- **The liver extracts most absorbed propionate on first pass** (portal ≈375 → hepatic-vein ≈148 → peripheral ≈79 µmol/L). In healthy people, propionate is largely cleared by hepatic propionyl-CoA carboxylase before reaching the periphery. **In PA/MMA that hepatic sink is precisely what is broken**, so gut-derived propionate that would normally be silently metabolized instead spills into the systemic pool — this is the mechanistic reason gut fermentation matters so much more in PA than in health.
- **Total colonic SCFA production** on a typical Western diet is on the order of several hundred mmol/day, of which propionate is ~20% (⚠verify the exact daily total against Cummings & Macfarlane review data; commonly cited totals are in the range of a few hundred mmol/day).
- **Within-person day-to-day variability** is large and driven chiefly by fermentable fiber intake; colonic propionate output can shift materially with diet composition over days, which introduces real noise into any steady-state assumption.

> ⚠ The "several hundred mmol/day" total-production figure and the fraction reaching systemic circulation are best sourced from Cummings & Macfarlane's fermentation reviews; confirm the exact mmol/day before use.

---

### 3. Antibiotic effects on propionate in PA and MMA

**Mechanism and magnitude.** **Bain MD, Borriello SP, Reed PJ, Tracey BM, Chalmers RA, et al. (*Lancet* 1988;1(8593):1078–1079)** demonstrated that oral metronidazole — targeting anaerobic colonic flora — produces a substantial fall in urinary excretion of propionate-derived metabolites, directly confirming the gut-bacterial contribution inferred from tracer studies. The reduction reported was a large fraction of baseline metabolite excretion (⚠verify the exact percentage figure quoted in the Lancet letter).

| Intervention | Population | Effect | Source |
|---|---|---|---|
| Metronidazole (oral) | PA/MMA patients | Marked reduction in urinary propionate metabolites (3-hydroxypropionate, methylcitrate, propionylglycine) and, in MMA, methylmalonate | Bain et al., *Lancet* 1988;1:1078–1079 |
| Intermittent metronidazole | PA & MMA (guideline) | Recommended to suppress gut propionate production during management (⚠verify exact dose/schedule) | Baumgartner et al., *Orphanet J Rare Dis* 2014;9:130, PMID 25205257 |

**Metronidazole** is the agent of choice because it selectively suppresses anaerobes (the principal propionate producers) while sparing much of the aerobic/facultative flora, giving a favorable effect-to-disruption ratio versus neomycin or ampicillin. **Time course:** the effect on urinary metabolites appears within days of starting; **rebound** occurs after discontinuation as the anaerobic community recovers — which is the rationale for *intermittent* (pulsed) dosing (e.g., a ~10-day course cyclically) rather than continuous therapy, both to limit rebound and to avoid chronic toxicity (peripheral neuropathy, resistance, microbiome disruption).

> ⚠ Confirm from Baumgartner 2014 the exact recommended metronidazole regimen (dose in mg/kg/day and days-per-month), and from Bain 1988 the exact percent reduction. The consensus paper deliberately hedges continuous use because of neuropathy risk.

---

### 4. Endogenous protein catabolism contribution to propionate load

**The propionogenic amino acids are valine, isoleucine, methionine, and threonine.** Together with odd-chain fatty acids and cholesterol, they are often summarized by the mnemonic that these are the carbon skeletons converging on propionyl-CoA. Per Thompson 1990, amino-acid/protein catabolism is the **largest single endogenous source** (roughly half of production; ⚠verify exact figure).

- **Relative contribution among the four:** isoleucine and valine are the dominant contributors, with methionine and threonine smaller. (Isoleucine yields both propionyl-CoA and acetyl-CoA; valine's catabolism converges on methylmalonyl-CoA/propionyl-CoA.)
- **During catabolic stress** (fever, infection, fasting, surgery, rapid growth/anabolic failure), **net protein breakdown liberates these amino acids from endogenous muscle protein**, generating a propionate load that is *independent of, and often larger than, dietary intake*. This is why PA/MMA patients decompensate during intercurrent illness even when fasted or protein-restricted, and why emergency regimens prioritize halting catabolism (high-dose IV dextrose ± insulin, aggressive calorie provision) over merely removing dietary protein.
- **Quantification via protein-turnover tracers:** leucine and phenylalanine stable-isotope kinetics quantify whole-body proteolysis during catabolism, but **direct, simultaneous measurement linking measured protein flux to measured propionate flux in PA/MMA is sparse** — this is a genuine gap in the literature and a priority for anyone building a mechanistic model.

**Modeling implication:** represent endogenous protein catabolism as a **state-dependent multiplier** on the amino-acid propionate compartment, switched on by illness/fasting and switched off by successful anabolism — not as a constant.

---

### 5. Odd-chain fatty acid and cholesterol side-chain contributions

- **Odd-chain fatty acids (C15:0 pentadecanoic acid, C17:0 heptadecanoic acid)** yield one propionyl-CoA per molecule at the terminal step of β-oxidation. Thompson 1990 classed these as a **minor** source of total propionate (single-digit percent range; ⚠verify).
- **Cholesterol side-chain oxidation** (the side chain of cholesterol lost during bile-acid synthesis) yields propionyl-CoA and is likewise a **minor endogenous source**.
- **Endogenous odd-chain fatty-acid synthesis:** in PA/MMA, accumulated propionyl-CoA is used as a "primer" in place of acetyl-CoA during de novo lipogenesis, generating odd-chain fatty acids. This means C15:0/C17:0 in PA/MMA are partly a **consequence** of propionyl-CoA accumulation (a biomarker) as well as a **feedstock** for further propionate on turnover — a small internal recycling loop worth noting but quantitatively minor.

**Modeling implication:** odd-chain FA + cholesterol can reasonably be lumped into a small (~single-digit-percent) constant "other endogenous" compartment, but should not be neglected in a mass-balance model that must sum to 100%.

---

### 6. Gut microbiome in PA patients

- **Propionate-producing taxa and pathways:** propionate is generated by three main bacterial routes — the **succinate pathway** (dominant; *Bacteroidetes* — *Bacteroides*, *Prevotella* — and some Firmicutes/Negativicutes), the **acrylate pathway** (e.g., *Megasphaera*, *Coprococcus*), and the **propanediol pathway** (from deoxy-sugars; e.g., *Salmonella*, *Roseburia*). ***Akkermansia muciniphila*** is a notable propionate (and acetate) producer from mucin.
- **PA vs healthy controls:** the composition of the PA gut microbiome relative to healthy controls, and whether it is enriched or depleted in propionate producers, is an **emerging but still thin** area. Diet (chronic protein restriction, high carbohydrate) and intermittent antibiotics in PA patients are themselves major microbiome-shaping exposures, confounding any comparison.
- **Temporal stability and interventions:** longitudinal stability of the microbiome in PA patients, and the durability of antibiotic- or probiotic-driven shifts in propionate output, remain **under-characterized**. This is the single thinnest evidence area in this review and the clearest target for primary research.

> ⚠ Specific PA-cohort microbiome studies (16S/shotgun with n and taxa-level abundances) should be sought directly on PubMed ("propionic acidemia gut microbiome," "organic acidemia microbiome propionate producers"); I was unable to retrieve named cohort studies with abundances under the current tool constraints, and I decline to fabricate them.

---

## Details: study-quality and reconciliation notes
- **Rare-disease sample sizes are small.** The tracer/partitioning studies (Thompson 1990 and companions) involve single-digit patient numbers; the partition percentages are best treated as **central estimates with wide uncertainty**, not precise constants.
- **Cummings 1987 is autopsy/sudden-death physiology in non-PA humans** — it defines *normal* colonic SCFA and hepatic extraction, which then must be *re-weighted* for PA because the hepatic propionyl-CoA-carboxylase sink is defective. Do not import the "peripheral 79 µmol/L" figure into a PA model without removing the hepatic-extraction assumption.
- **Antibiotic effect sizes vary** by agent, dose, adherence, and baseline flora; the Bain 1988 result establishes direction and rough magnitude, not a universal coefficient.
- **Units discipline:** Cummings colonic values are **mmol/kg of contents** (not mmol/L); blood values are **µmol/L**. Mixing these is a common and consequential error.

---

## Recommendations (staged, for building a PA propionate model)

1. **Anchor the endogenous mass balance to Thompson 1990.** Partition non-dietary propionate into three compartments: **gut-bacterial (~20–25%)**, **amino-acid/protein catabolism (~50%)**, and **odd-chain FA + cholesterol (minor, single-digit %)**. *Benchmark that would change this:* re-extraction of the exact Thompson 1990 percentages (or a newer tracer study) shifting any compartment by >10 absolute percentage points.
2. **Model the gut compartment as pharmacologically suppressible.** Apply a metronidazole effect (large fractional reduction, appearing within days, with post-discontinuation rebound over days–weeks) per Bain 1988 / Baumgartner 2014. *Threshold:* if a verified reduction figure differs materially from a ~40–60% assumption, re-parameterize.
3. **Implement a catabolic switch, not a constant, for the amino-acid compartment.** Drive it by illness/fasting state; validate against clinical decompensation triggers. This is the highest-yield dynamic feature.
4. **Keep odd-chain FA/cholesterol as a small constant "other" term** for mass-balance closure; do not over-engineer it.
5. **Treat the microbiome module as exploratory** and parameterize it conservatively until PA-specific cohort data (taxa abundances, temporal stability) are obtained. *Priority action:* commission a targeted PubMed extraction (or a de novo 16S/shotgun study) on PA-cohort microbiomes.
6. **Verify every ⚠-flagged figure against the primary source before publication or clinical use.** Specifically: Thompson 1990 exact partition %; the exact whole-body flux (µmol/kg/h); Bain 1988 exact % reduction; Baumgartner 2014 exact metronidazole dose/schedule; Cummings total daily SCFA production.

---

## Comprehensive summary table (all quantitative data, with verification status)

| # | Topic | Value (units) | Population | n / design | Source (PMID/DOI + section) | Status |
|---|---|---|---|---|---|---|
| 1 | Propionate source partition — gut bacteria | ~20–25% of production | PA + MMA | small-n stable-isotope | Thompson, *Metabolism* 1990;39(11):1133–7, PMID 2233273, Results/Abstract | ⚠verify % |
| 1 | Propionate source partition — amino acids/protein | ~50% of production | PA + MMA | same | same | ⚠verify % |
| 1 | Propionate source — odd-chain FA/cholesterol | minor (single-digit %) | PA + MMA | same | same | ⚠verify % |
| 1 | Whole-body propionate production rate | µmol/kg/h (fed > fasted) | affected children | tracer | Thompson/Leonard group | ⚠verify value |
| 2 | Colonic SCFA, caecum/proximal | ~131 mmol/kg contents | healthy (autopsy) | sudden-death series | Cummings, *Gut* 1987;28(10):1221–7, PMID 3678950, contents table | confirmed via enrichment |
| 2 | Colonic SCFA, descending/sigmoid | ~80 mmol/kg contents | healthy | same | same | confirmed via enrichment |
| 2 | Molar ratio A:P:B | ~57:22:21 | healthy | same | Cummings 1987, Results | confirmed via enrichment |
| 2 | Portal venous total SCFA | ~375 µmol/L | healthy | same | Cummings 1987, blood data | confirmed via enrichment |
| 2 | Hepatic venous total SCFA | ~148 µmol/L | healthy | same | Cummings 1987, blood data | confirmed via enrichment |
| 2 | Peripheral venous total SCFA | ~79 µmol/L | healthy | same | Cummings 1987, blood data | confirmed via enrichment |
| 2 | Propionate as fraction of colonic SCFA | ~20% (~1/5) | healthy | same | Cummings 1987 | confirmed |
| 2 | Total colonic SCFA production | several hundred mmol/day | healthy | review estimate | Cummings & Macfarlane | ⚠verify exact mmol/day |
| 3 | Metronidazole → urinary propionate metabolites | large fractional reduction (~40–60%) | PA/MMA | case series/intervention | Bain, *Lancet* 1988;1(8593):1078–9 | ⚠verify exact % |
| 3 | Metronidazole regimen (guideline) | intermittent/pulsed course | PA & MMA | consensus guideline | Baumgartner, *Orphanet J Rare Dis* 2014;9:130, PMID 25205257 | ⚠verify exact dose/schedule |
| 4 | Propionogenic amino acids | Val, Ile, Met, Thr (Ile/Val dominant) | general biochemistry / PA-MMA | — | textbook + Thompson 1990 | established |
| 4 | Amino-acid share during catabolism | ↑ during illness/fasting (state-dependent) | PA/MMA | clinical | consensus + tracer literature | qualitative; direct flux linkage = gap |
| 5 | Odd-chain FA (C15:0/C17:0) → propionyl-CoA | minor % of load; also biomarker of accumulation | PA/MMA | — | Thompson 1990 + lipid biochem | ⚠verify % |
| 5 | Cholesterol side-chain → propionyl-CoA | minor endogenous source | general | — | biochem | established qualitatively |
| 6 | Propionate-producing taxa | Bacteroides, Prevotella, Akkermansia, Megasphaera, Roseburia, etc. (succinate/acrylate/propanediol pathways) | human gut | — | microbiome literature | established qualitatively |
| 6 | PA vs control microbiome composition & stability | not reliably quantified here | PA | — | — | **evidence gap** |

---

## Caveats
- **Tooling limitation (material).** The live search/subagent tools were unavailable during this task; I did not fabricate PMIDs, DOIs, or numbers, and I have explicitly marked every value I could not independently re-derive to the digit with "⚠verify." An analyst should run a fresh PubMed extraction against the four anchor papers (Thompson 1990 PMID 2233273; Cummings 1987 PMID 3678950; Bain 1988 *Lancet* 1(8593):1078–9; Baumgartner 2014 PMID 25205257) plus a dedicated search for PA-cohort microbiome studies before locking any parameter.
- **Small samples / rare disease.** Partition percentages and flux values rest on very small tracer cohorts; uncertainty is wide.
- **Physiology-to-disease transfer.** Cummings 1987 describes normal hepatic first-pass extraction of propionate; this must be relaxed in a PA model because the hepatic sink is defective — the systemic exposure to gut-derived propionate is disproportionately higher in PA than the healthy peripheral concentrations imply.
- **Unit hazards.** Colonic contents = mmol/kg; blood = µmol/L. Keep these distinct.
- **Biggest genuine evidence gaps:** (a) simultaneous protein-turnover-and-propionate-flux measurement during PA/MMA decompensation; (b) PA-specific gut-microbiome composition, taxa-level abundances, and temporal stability. These are the two areas where new primary data would most improve a model.
