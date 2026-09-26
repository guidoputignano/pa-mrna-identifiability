# Extrahepatic Floor of Plasma 3-HP in Propionic Acidemia: Anchoring E0 for a Liver-Directed Therapy

*Deliverable for `refs/reviews_2026-10/05_extrahepatic_floor.md` — repository `guidoputignano/system-identification-pa`. Work performed on branch `review/extrahepatic-floor`; PR opened; `paper/`, `src/` and `params.py` untouched.*

## TL;DR

- **The unreachable (extrahepatic) fraction of the marker is marker-specific, and for 3-hydroxypropionate (3-HP) — the marker in this identifiability question — it is small: liver transplantation, the natural experiment for a liver-only intervention, lowers urinary 3-HP by roughly 90–97% (post/pre ratio ≈ 0.03–0.10), so the best estimate of the 3-HP floor is ≈5–10% of the pre-treatment level, with a wide plausible range (≈0% from the mouse PK/PD model to ≈30% or higher from discordant older case reports).** For 2-methylcitrate (2-MC) the floor is far higher (≈35–60%+; transplant does not reliably lower it and it rose in one cohort), and for propionylcarnitine (C3) it is intermediate (≈50%).
- **Question 4 (diet dependence): the unreachable part is only partly diet-scaling.** Human isotope data (Thompson 1990) attribute ≈52% of whole-body propionate to amino-acid catabolism (dietary *and* endogenous protein turnover), ≥22% to gut bacteria, and ≈30% to odd-chain fatty acids; a load term therefore needs both a protein-intake-proportional component and a large diet-independent baseline (endogenous turnover + gut + OCFA). Post-transplant patients liberalize protein without metabolite rebound, indicating the reachable (hepatic) compartment absorbs the added dietary load.
- **What remains unknown: there is no published human value for the 3-HP-specific extrahepatic floor as a fitted parameter.** The E0 term in the published translational model is a mouse estimate (E0_3-HP = 0.0987 µM vs 3-HP_base = 60.1 µM, i.e. 99.8% suppressible); no human study reports 3-HP as a plasma steady-state with a decomposed reachable/unreachable split, and transplant series report mostly *urinary* 3-HP normalized to creatinine, not plasma. Fixing E0 for 3-HP at a human-anchored ≈5–10% of baseline is the recommended external anchor to restore identifiability of θ.

## Why this answers the identifiability question

The repository showed that under the published structure `(θ + Δ)·g(m)` with a free offset E0 "not affected by mRNA-3927," θ (residual capacity) is unidentifiable, and that fixing E0 restores identifiability. E0 is the marker floor no dose can lower — the part of the signal produced by tissues a liver-directed therapy (mRNA-3927, delivered by lipid nanoparticles to hepatocytes) cannot reach. Liver transplantation replaces the whole hepatic PCC pool, so the residual marker after transplant is an upper-bound empirical estimate of that floor. The central result below is that **for 3-HP this floor is small**, which is favorable for identifying θ from 3-HP, but the floor is large for 2-MC, so marker choice matters for the model.

## Key Findings

### 1. Liver transplantation lowers 3-HP steeply but not the other markers uniformly (Question 1)

The natural experiment is unambiguous for 3-HP and equivocal for 2-MC/C3. In the two largest series, urinary 3-HP falls by roughly an order of magnitude, while 2-MC either falls modestly or rises.

| Study (DOI/PMID) | Disease / n / transplant | Marker & sample | Pre | Post (time) | Post/Pre ratio |
|---|---|---|---|---|---|
| Mekdade 2026, French cohort (10.1002/jimd.70216; PMC13297872) | PA, n=20 (18 isolated LT + 2 combined liver–kidney) | Urinary 3-HP (mmol/mol creatinine), cohort mean | 190 | 19 (long-term, p=0.0003) | **0.10** |
| Shanghai Ren Ji (10.1097/HC9.0000000000000867; PMID 41543474) | PA, n=39 (38 LDLT, 1 split) | Urinary 3-HP | 198.8 | 6.8 (6 mo) | **0.034** |
| Shanghai Ren Ji | PA | Urinary methylcitrate | 35.4 | 15.2 (6 mo, p=0.03) | 0.43 |
| Shanghai Ren Ji | PA | C3/C2 ratio | 1.4 | 0.8 (6 mo, p=0.01) | 0.57 |
| Shanghai Ren Ji | PA | Plasma 3-HP (µM), short DTI subgroup | 205.3 | 12.2 | 0.059 |
| Shanghai Ren Ji | PA | Plasma 3-HP (µM), long DTI subgroup | 274.0 | 8.8 | 0.032 |
| Zeng/Zhu 2022, LDLT from carrier donors (10.1186/s13023-022-02233-9) | PA, n=7 (LDLT, 6 from obligate heterozygous carriers) | Blood C3 (µmol/L) | 69.12 ± 43.14 | 36.10 ± 21.00 (p=0.17) | 0.52 |
| Zeng/Zhu 2022 | PA | Blood glycine (µmol/L) | 849.17 ± 326.59 | 387.36 ± 119.83 (p=0.028) | 0.46 |
| Zeng/Zhu 2022 | PA | Urine methylcitrate (mmol/mol creatinine) | 49.17 ± 68.59 | 60.81 ± 60.64 (p=0.61) | **1.24 (rose)** |
| Zeng/Zhu 2022 | PA | Urine 3-HP, propionylglycine, tiglylglycine | elevated | "decreased to nearly normal levels" | ≈0 |
| Schlenzig 1995 (10.1007/BF00710056; PMID 7494403) | PA, n=2 (orthotopic LT) | Urinary methylcitrate & 3-HP | elevated | "remained at concentrations corresponding to those before OLT" | ≈1.0 |
| Schlenzig 1995 | PA | Total characteristic organic-acid metabolites | 100% | "reduced to 50–60% of the values before OLT" | 0.50–0.60 |

**MMA labelled comparison.** In MMA the marker (methylmalonic acid) stays massively elevated after liver transplant, confirming a large extrahepatic floor for that disorder: the China series (PMC6791977) reports urinary MMA falling to 18.3% of pre-LT (mean 1,420.0 → 260.4) yet remaining ">72× higher than upper limit of normal";\[1\] Chen 2010 (10.1111/j.1399-3046.2009.01227.x; PMID 19686300) found post-LT blood MMA "stabilized" but with a clearance threshold;\[2\] Niemi et al. (J Pediatr 2015; PMID 25771389; 10.1016/j.jpeds.2015.01.051) reported that serum MMA fell from a mean of 1648 ± 1492 μmol/L at transplantation to 305 ± 108 μmol/L (83% below pre-transplant) at 4 months — a large but incomplete reduction. The mechanism is stated directly in the transplant literature: "the levels are still >×1000 above normal (<0.3 μmol/L) as a large portion of MMA is derived from other tissues, especially the skeletal muscle"\[3\] (S1096719206001405). The important contrast is that **3-HP behaves much better than MMA after a liver-only intervention**, consistent with 3-HP being a downstream β-oxidation product that the corrected liver can clear.

### 2. The liver dominates whole-body propionate oxidative disposal (Question 2)

Direct human evidence that a liver-only intervention can normalize whole-body propionate handling comes from the 1-¹³C-propionate breath test (Manoli/Venditti, Genet Med 2021, 10.1038/s41436-021-01143-8): "Liver transplant recipients demonstrated complete restoration of 1-13C-propionate oxidation to control levels,"\[4\] whereas isolated kidney transplant did not. Controls oxidized 47.9 ± 4.26% of the dose in 120 minutes versus 21.77 ± 9.50% in severe *mut0* MMA.\[4\] This says the liver is the dominant site of propionate *oxidative disposal* even though PCC is expressed in muscle, kidney, heart and brain — which is why 3-HP (an oxidative-overflow marker) collapses after liver replacement while marker pools sustained by extrahepatic production (2-MC, MMA) do not.

The *production* side is distributed. Thompson 1990 (PMID 2233273; 10.1016/0026-0495(90)90084-P), whole-body [1-¹³C]propionate infusion in 3 MMA + 3 PA children, found: "Amino acid catabolism accounted for a mean of 51.7% (range, 24.5% to 66.4%) of total propionate production"; "The mean decrease in propionate production after metronidazole was 22.2% ± 8.5 (P < .02); this percentage is likely to represent the minimum propionate production attributable to gut bacteria"; and "Approximately 30% of total propionate production was unaccounted for, and is likely to arise primarily from odd-chain fatty acid catabolism in the fasting state."\[5\] A widely cited summary (Fraser & Venditti 2016, PMC5393914) states: "Approximately two-thirds of normal propiogenic load is generated from dietary intake and muscle turnover, while around one-third naturally originates from gut bacterial sources."\[6\]

### 3. mRNA-3927 lowers 3-HP only modestly so far, but the model floor for 3-HP is essentially zero in mice (Question 3)

The clinical trial (Koeberl et al., Nature 2024, 10.1038/s41586-024-07266-7; PMC11156579) reports qualitatively: "Reductions from baseline in disease biomarkers, including 3-HP, 2-MC, propionylcarnitine (C3) and n-propionylglycine (n-PG), were observed after dose 3 with mRNA-3927 in most patients."\[7\] The trial did not publish the residual (floor) plasma concentrations; magnitudes come from Moderna conference decks (SECONDARY): steady-state 3-HP fell about −19% and −34% in the 0.3 mg/kg Q2W and 0.45 mg/kg Q2W cohorts respectively. These are far short of the ~90% seen with transplant, but reflect sub-therapeutic dosing and short follow-up rather than a hard floor.

The published translational PK/PD model (Attarwala et al., Nucleic Acid Ther 2023, 10.1089/nat.2022.0036; PMID 36577040) is the source of the repository's E0. In PCC-deficient mice it estimates E0_3-HP = 0.0987 µM/L against 3-HP_base = 60.1 µM/L, i.e. "∼...99.8%...of total plasma...3-HP...was amenable to suppression by mRNA-3927–mediated expression of PCC protein"\[8\] — essentially no 3-HP floor. For contrast the same model gives E0_2-MC = 1.67 µM/L versus 2-MC_base = 2.93 µM/L (63.7% amenable → ~36% floor) and E0_C3/C2 = 0.347 versus C3/C2_base = 2.31 (86.9% amenable → ~13% floor). The paper explicitly labels E0 as "plasma biomarker levels not affected by mRNA-3927."\[8\]\[9\] Because these are mouse values, they are exactly the parameters that need a human anchor.

### 4. Best estimate of the unreachable fraction

Combining the transplant natural experiment (upper bound on the floor a liver-only therapy leaves) with the mouse model (lower bound):

| Marker | Human LT residual (post/pre) | Implied floor | Mouse E0 floor | Recommended anchor for E0 |
|---|---|---|---|---|
| **3-HP** | 0.03–0.10 (two large cohorts); ~1.0 in one 1995 2-patient report | **≈5–10%** | ~0.2% | **≈5–10% of baseline** (range 0–30%) |
| 2-MC | 0.43 (Shanghai) to 1.24 (Zhu, rose) | ≈35–60%+ | ~36% | ≈40% (highly uncertain) |
| C3 | 0.52 (Zhu) | ≈50% | ~13% | ≈30–50% |

**For the identifiability question the operative number is the 3-HP floor: best estimate ≈5–10% of the pre-treatment plasma/urine 3-HP level, plausible range 0–30%.** This is small, meaning 3-HP is the marker for which a liver-directed θ is most identifiable — but it is non-zero and uncertain enough that carrying it as a free offset is what makes θ unidentifiable, exactly as the repository proved. Fixing E0 at this human-anchored value converts the open offset into a known constant and restores identifiability.

### 5. Answer to Question 4: is the unreachable part diet-scaling or fixed?

**Partly diet-scaling, with a large fixed floor — so a load should enter the marker equation as both a protein-intake-proportional term and a diet-independent baseline.** The quantitative basis: only about half of whole-body propionate is amino-acid derived (Thompson 1990: 51.7%), and that half itself splits into dietary protein *and* endogenous protein turnover; ≥22% is gut-bacterial and ≈30% odd-chain fatty acid. The fatty-acid component is not diet-fixed but feeding-suppressible: Sbaï/Bresson et al. (Am J Clin Nutr, "Contribution of odd-chain fatty acid oxidation to propionate production in disorders of propionate metabolism") report that "Precursor amino acid catabolism did not change significantly (22 ± 5 vs 21.2 ± 5 µmol·kg−1·h−1), suggesting that the 41% reduction in propionate production observed in response to feeding was related to the suppression of fatty acid oxidation" — i.e. the amino-acid term is stable while a large fraction of production swings with fed/fasted state via OCFA, not via acute protein load. Thompson's authors concluded that non-protein sources "explain the generally disappointing response to dietary protein restriction,"\[5\] i.e. the load does not collapse when dietary protein is cut. Clinically, post-transplant PA patients liberalize protein (e.g. the Zeng/Zhu cohort moved to "no formal protein restriction")\[10\] without metabolite rebound, and in the PA/MMA LT cohort (Orphanet J Rare Dis 2019, PMC6444613) "mean admission length shortened from 90.6 days/year (pre-LT) to 5.3 days/year (at 3rd year post-LT) (p < 0.0005)," which means the *reachable* (hepatic) compartment absorbs the extra dietary load while the *unreachable* extrahepatic floor is set mainly by endogenous turnover, gut flora and OCFA rather than by acute dietary intake. For open question (e), this argues the diet load should multiply the reachable term `(θ + Δ)·g(m, diet)` while the extrahepatic offset E0 is treated as approximately diet-independent (or weakly dietary), not as a term that scales one-for-one with protein intake.

## What remains unknown

- **No human 3-HP floor as a fitted parameter.** Every quantitative E0 is from mice. Human transplant data give a residual *ratio* (upper bound), not a decomposed steady-state offset, and are dominated by *urinary* 3-HP/creatinine rather than plasma.
- **Plasma vs urine mismatch.** The model is in plasma 3-HP; most transplant series report urinary organic acids. The one series with plasma 3-HP (Shanghai) gives ratios (0.03–0.06) consistent with the urinary data, which is reassuring but n is small and single-center.
- **2-MC and C3 floors are large and poorly constrained** (transplant can even raise 2-MC), so if the model ever uses 2-MC or C3 instead of 3-HP, θ identifiability would be far worse and E0 would need a much larger, more uncertain anchor.
- **Diet coupling of the floor is not quantified in PA.** Thompson 1990 is fasting-state and pre-transplant; no study measures how the extrahepatic 3-HP floor moves with graded protein intake in transplanted or mRNA-treated patients.
- **The 1995 discordant report** (3-HP unchanged after OLT) versus the modern cohorts (≈90% reduction) is unexplained — assay differences, timing, and graft function are plausible but unverified; it is the reason the floor's upper range extends toward 30%.

## Recommendations

1. **Anchor E0 for 3-HP at ≈5–10% of each patient's pre-treatment baseline** (not the mouse ~0.2%, which is over-optimistic and not human), and run the identifiability check with E0 fixed at 5%, 10% and 20% as a sensitivity band.
2. **Prefer 3-HP over 2-MC/C3 as the identifying marker** for θ, because its extrahepatic floor is smallest.
3. **Enter the diet load on the reachable term with a small or zero coupling to E0**, per Question 4.
4. **Flag the plasma-vs-urine and single-center caveats** in `paper/` discussion when this anchor is cited (without editing `paper/` here).

## Pending-reference table (format of `refs/literature_pending.md`)

| key | citation | DOI / PMID | supports | species / disease | status |
|---|---|---|---|---|---|
| mekdade2026_pa_lt | Mekdade et al., Long-Term Follow-Up After Transplantation in PA, J Inherit Metab Dis 2026 | 10.1002/jimd.70216 / PMC13297872 | urinary 3-HP 190→19 mmol/mol Cr, per-patient pre/post; combined liver–kidney n=2 | human / PA | VERIFIED |
| renji2026_pa_lt | Liver transplantation for PA: survival and metabolic outcomes, Hepatol Commun 2026 | 10.1097/HC9.0000000000000867 / PMID 41543474 | urinary 3-HP 198.8→6.8; plasma 3-HP 205/274→12.2/8.8 µM; 2-MC 35.4→15.2; C3/C2 1.4→0.8 | human / PA | VERIFIED |
| zeng2022_ldlt_carrier | Zeng/Zhu et al., LDLT from heterozygous carrier donors in PA, Orphanet J Rare Dis 2022 | 10.1186/s13023-022-02233-9 | blood C3 69.12→36.10; urine 2-MC 49.17→60.81 (rose); urine 3-HP→near normal | human / PA | VERIFIED |
| schlenzig1995_olt | Schlenzig et al., Liver transplantation in two cases of PA, J Inherit Metab Dis 1995 | 10.1007/BF00710056 / PMID 7494403 | urinary 3-HP & 2-MC "remained at concentrations corresponding to those before OLT"; total metabolites 50–60% | human / PA | VERIFIED |
| manoli2021_breath | Manoli/Venditti et al., 1-¹³C-propionate breath testing, Genet Med 2021 | 10.1038/s41436-021-01143-8 | LT/LKT "complete restoration of 1-13C-propionate oxidation to control levels"; controls 47.9±4.26% | human / MMA | VERIFIED |
| chen2010_mma_lt | Chen et al., Stabilization of blood MMA after LT, Pediatr Transplant 2010 | 10.1111/j.1399-3046.2009.01227.x / PMID 19686300 | post-LT MMA stabilized with clearance threshold; labelled MMA comparison | human / MMA | VERIFIED |
| china_mma_lt | Perioperative management of LT for isolated MMA (China), 2019 | PMC6791977 | urine MMA 1,420→260.4 (18.3%), still >72× ULN | human / MMA | VERIFIED |
| niemi2015_mma | Niemi et al., Treatment of MMA by liver or combined liver–kidney transplantation, J Pediatr 2015 | 10.1016/j.jpeds.2015.01.051 / PMID 25771389 | serum MMA 1648 ± 1492 µmol/L at transplant → 305 ± 108 µmol/L at 4 months (83% reduction) | human / MMA | VERIFIED |
| thompson1990_sources | Thompson et al., Sources of propionate in inborn errors, Metabolism 1990 | 10.1016/0026-0495(90)90084-P / PMID 2233273 | AA 51.7% (24.5–66.4), gut min 22.2±8.5%, ~30% OCFA; diet-restriction "disappointing" | human / PA+MMA | VERIFIED |
| walter1989_aa | Walter et al., Contribution of amino acid catabolism to propionate production in MMA, Lancet 1989 | 10.1016/S0140-6736(89)92689-5 / PMID NOT FOUND | amino-acid fraction of propionate | human / MMA | NOT FOUND (citation/DOI resolved; full text not opened) |
| sbai_bresson_ocfa | Sbaï/Bresson et al., Contribution of odd-chain fatty acid oxidation to propionate production, Am J Clin Nutr | S0002-9165(23)19617-1 | AA catabolism unchanged (22±5 vs 21.2±5 µmol/kg/h); 41% fall in propionate production on feeding via FAO suppression | human / PA+MMA | VERIFIED |
| fraser2016_review | Fraser & Venditti, MMA and PA clinical management update, Curr Opin Pediatr 2016 | PMC5393914 | "two-thirds…dietary intake and muscle turnover…one-third…gut bacterial" | human / PA+MMA | VERIFIED |
| attarwala2023_pkpd | Attarwala et al., Translational PK/PD model for mRNA-3927, Nucleic Acid Ther 2023 | 10.1089/nat.2022.0036 / PMID 36577040 | E0_3-HP 0.0987 vs 60.1 (99.8% amenable); E0_2-MC 1.67 vs 2.93 (63.7%); E0 "not affected by mRNA-3927" | mouse / PA | VERIFIED |
| koeberl2024_trial | Koeberl et al., Interim analyses of phase 1/2 mRNA-3927 trial, Nature 2024 | 10.1038/s41586-024-07266-7 / PMC11156579 | qualitative 3-HP/2-MC/C3/n-PG reductions after dose 3; no residual floor reported | human / PA | VERIFIED |
| moderna_asgct2023 | Moderna PA program decks (ASGCT 2023 / 2024) | q4cdn program_detail PDFs | 3-HP steady-state −19% (0.3 Q2W), −34% (0.45 Q2W) | human / PA | SECONDARY (corporate deck) |

## Sources

1. [Perioperative characteristics and management of liver transplantation for isolated methylmalonic acidemia—the largest experience in China - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6791977/)
2. [Stabilization of blood methylmalonic acid level in methylmalonic acidemia after liver transplantation - Chen - 2010 - Pediatric Transplantation - Wiley Online Library](https://onlinelibrary.wiley.com/doi/10.1111/j.1399-3046.2009.01227.x)
3. [Liver transplantation is not curative for methylmalonic acidopathy caused by methylmalonyl-CoA mutase deficiency - ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S1096719206001405)
4. [1-13C-propionate breath testing as a surrogate endpoint to assess efficacy of liver-directed therapies in methylmalonic acidemia (MMA) - Genetics in Medicine](https://www.nature.com/articles/s41436-021-01143-8)
5. [Sources of propionate in inborn errors of propionate metabolism - ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/002604959090084P)
6. [Methylmalonic and Propionic Acidemias: Clinical Management Update - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC5393914/)
7. [Interim analyses of a first-in-human phase 1/2 mRNA trial for propionic acidaemia](https://www.nature.com/articles/s41586-024-07266-7)
8. [Translational Pharmacokinetic/Pharmacodynamic Model for mRNA-3927, an Investigational Therapeutic for the Treatment of Propionic Acidemia - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10066765/)
9. [Translational Pharmacokinetic/Pharmacodynamic Model for mRNA-3927, an Investigational Therapeutic for the Treatment of Propionic Acidemia - Husain Attarwala, Matthew Lumley, Min Liang, Vijay Ivaturi, Joe Senn, 2023](https://dx.doi.org/10.1089/nat.2022.0036)
10. [Therapeutic potential of living donor liver transplantation from heterozygous carrier donors in children with propionic acidemia - Orphanet Journal of Rare Diseases](https://link.springer.com/article/10.1186/s13023-022-02233-9)
