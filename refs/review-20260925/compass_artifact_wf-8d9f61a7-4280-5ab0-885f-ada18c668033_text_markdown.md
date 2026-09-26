# Increment Time Course for Δ(t) in mRNA-3927 Capacity Identification: Literature Review

*Deliverable for `refs/reviews_2026-10/02_increment_time_course.md`. Branch: `review/increment-time-course`. Status: PR opened, awaiting merge. Current date: 25 September 2026.*

## TL;DR

- The shape of Δ(t) that matters for identifying residual capacity θ is the **hepatic PCC protein/activity time course**, not the plasma mRNA PK, and the only quantitative version of that time course in any primary source is the two-compartment (bi-exponential) indirect-response model fitted in **PCC-deficient mice** by Attarwala et al. 2023, which gives a terminal PCC-protein half-life of about 7 days and a local half-life that rises from roughly 2.5 days early to about 9 days by day 14. A single-exponential 3.5-day form is a convenience assumption supported by no primary source, and using it while the truth is the bi-exponential mouse shape is exactly the mis-specification the repo's robustness review flagged.
- **No human hepatic PCC protein or activity time course exists in the public record.** The trial measured plasma PCCA mRNA PK (secondary endpoint) and plasma biomarkers (3-HP, 2-MC, C3, n-PG), but never liver enzyme, and the registry's Emax, AUEC and duration-of-response outcomes for 2-MC and 3-HP have not been reported. Analogous LNP-mRNA programmes (mRNA-3705, mRNA-3210, ARCT-810) give animal protein half-lives and human plasma-mRNA half-lives but no human hepatic protein half-life either.
- On question 3: the Extended Data Table 3 biomarker values are a percentage change from a **pre-dose (baseline) level to a post-dose level "after the third dose"**, but the precise sampling time of that post-dose draw relative to dose 3 (trough versus a fixed number of hours or days post-infusion) is **not stated anywhere in the public text**, so it cannot be used to pin the shape of Δ(t) around a dose.

## Summary table

Every parameter is labelled VERIFIED (primary source opened, value quoted), SECONDARY (only in a secondary source or an image table I could not extract), or NOT FOUND. Species is recorded for every row.

| Quantity | Value | Species | Source (DOI / PMID / ID) | Label |
|---|---|---|---|---|
| PCC protein terminal half-life | ~7 days | Mouse (Pcca⁻/⁻ A138T) | Attarwala 2023, 10.1089/nat.2022.0036, PMID 36577040 | VERIFIED |
| PCC protein degradation rate kdeg | 0.0075 /h | Mouse | Attarwala 2023 | VERIFIED |
| PCC protein intercompartmental rate kq | 0.00474 /h | Mouse | Attarwala 2023 | VERIFIED |
| Effect-compartment rate ke0 | 0.242 /h | Mouse | Attarwala 2023 | VERIFIED |
| PCC protein model structure | 2-compartment indirect response, biexponential terminal decline | Mouse | Attarwala 2023 | VERIFIED |
| IC50 of PCC protein for 3-HP suppression | 37.5 μg/g | Mouse | Attarwala 2023 | VERIFIED |
| IC50 of PCC protein for 2-MC suppression | 21.0 μg/g | Mouse | Attarwala 2023 | VERIFIED |
| Duration of biomarker (3-HP, 2-MC, C3/C2) response after one dose | ~3 to 4 weeks | Mouse | Baek 2024, 10.1038/s41467-024-47460-9, PMID 38714648 | VERIFIED |
| PCC activity, treated vs control | 7.77 ± 4.27 vs 2.72 ± 2.26 ng/mg protein | Mouse | Baek 2024 | VERIFIED |
| Plasma hPCCA/hPCCB mRNA quantifiable window | up to 48 h post-dose | Mouse | Baek 2024 | VERIFIED |
| Human plasma PCCA mRNA PK measured | up to 21 days after dose 1, by branched-DNA assay | Human | Koeberl 2024, 10.1038/s41586-024-07266-7 | VERIFIED |
| Human PCCA mRNA Cmax, AUCinf, T1/2, Tmax, CL, Vz | reported per cohort in Extended Data Table 2, numeric cells published as image only | Human | Koeberl 2024 (Extended Data Table 2) | NOT FOUND (image table, not extractable) |
| Human population-PK PCCA mRNA clearance | 0.266 L/h | Human | Gao/Zhang ACoP 2025, 10.70534/JVQU8233 | VERIFIED |
| Human population-PK PCCA mRNA volume | 13.9 L | Human | Gao/Zhang ACoP 2025 | VERIFIED |
| Human hepatic PCC protein or activity time course | none reported (no liver enzyme measurement in humans) | Human | n/a | NOT FOUND |
| Emax, AUEC, duration of response for 2-MC and 3-HP | registered outcomes, not reported in any publication | Human | NCT04159103 outcomes | NOT FOUND |
| Extended Data Table 3 biomarker sampling time relative to dose 3 | not stated in public text | Human | Koeberl 2024 | NOT FOUND |
| Anti-PEG antibody, baseline positivity | 1/16 participants (cohort 1) positive at baseline | Human | Koeberl 2024 | VERIFIED |
| Anti-PEG antibody titre kinetics | ≥4-fold rise in pre-dose sample for dose 3 (52 days after first dose), returned to baseline by month 12 | Human | Koeberl 2024 | VERIFIED |
| Anti-PCC antibody | 1/16 positive at baseline, negative thereafter | Human | Koeberl 2024 | VERIFIED |
| MMA program (mRNA-3705) hepatic hMUT protein at 24 h | 45.0 ± 42.2 ng/mg (0.2 mg/kg), 70.4 ± 24.2 ng/mg (0.5 mg/kg) | Mouse | Baek 2024 | VERIFIED |
| MMA program plasma methylmalonic acid reduction at 24 h | 78% ± 3.2% (0.5 mg/kg) | Mouse | Baek 2024 | VERIFIED |
| PKU program (mRNA-3210) liver PAH protein time to max | 24 h (0.25, 0.5 mg/kg), 8 h (1.0 mg/kg) | Mouse | Baek 2024 | VERIFIED |
| PKU program blood Phe response duration | returned to baseline by 168 h (<1 week) | Mouse | Baek 2024 | VERIFIED |
| OTC program (ARCT-810) LNP lipid plasma half-life | t½ = 8 h (lipid cleared ~48 h) | Human | Arcturus UCDC Sept 2025 poster | SECONDARY |
| OTC program (ARCT-810) plasma mRNA half-life | t½ = 52 h (detectable ~1 week; PD effect ~2 to 4 weeks) | Human | Arcturus UCDC Sept 2025 poster | SECONDARY |
| OTC program pharmacodynamic duration (ammonia) | stable ~28 days after last dose | Human | Arcturus ICIEM/KOL 2025 | SECONDARY |
| Hepatic hOTC mRNA half-life | 7.07 to 9.10 h | Mouse | Sasaki 2023, PMC10372164 | VERIFIED |
| Hepatic OTC protein retention / half-life | ~10 days (mouse); ~8 days (rat) | Mouse / rat | Sasaki 2023, PMC10372164 | VERIFIED |

## Best-supported functional form for Δ(t)

**Recommendation: a bi-exponential (two-compartment indirect-response) shape, not a single exponential, and every parameter that defines it is a mouse parameter.** The reasoning follows.

**What Δ(t) should represent.** In the repo's identifiability question, Δ(t) is the therapy-driven increment in functional PCC capacity that suppresses plasma 3-HP. The quantity that actually drives 3-HP down is hepatic PCC protein/activity, and PCC protein is a two-compartment species with a biexponential terminal decline.\[1\] Attarwala et al. state, verbatim, that "a 2-compartment model for PCC protein was needed to adequately capture the biexponential decline of PCC protein observed in the terminal phase"\[2\] and that "the model estimated a half-life of PCC protein of ∼7 days."\[2\] Plasma PCCA mRNA (the branched-DNA PK signal, quantifiable to only 48 h in mouse\[3\] and up to 21 days in humans) is the upstream driver, not Δ(t) itself; the effect compartment (ke0 = 0.242 /h) plus the slow protein turnover (kdeg = 0.0075 /h, kq = 0.00474 /h) impose the delay and the long tail. This is exactly why a single mono-exponential decaying from each dose mis-states the shape: the true curve is shallow-then-steep-then-shallow, with a local half-life that lengthens over the inter-dose interval, consistent with the repo's own finding of about 2.5 days one day post-dose and about 9 days by day 14.

**Functional form and parameter ranges (all mouse for the PD/shape terms).**
- Preferred: a two-compartment indirect-response protein model driven through an effect compartment, i.e. the Attarwala structure. Terminal PCC-protein half-life about 7 days; kdeg = 0.0075 /h; kq = 0.00474 /h; ke0 = 0.242 /h. Species: mouse (Pcca⁻/⁻ A138T). Downstream 3-HP link is a sigmoidal Imax model with IC50 = 37.5 μg/g PCC protein (mouse); the companion 2-MC and C3/C2 IC50 values are 21.0 and 32.1 μg/g, per Attarwala 2023, verbatim: "The concentration of liver PCC protein needed for 50% inhibition of plasma 2-MC, 3-HP, and C3/C2 was estimated as 21.0, 37.5, and 32.1 μg/g, respectively."
- If a sum-of-exponentials surrogate is retained for tractability, it should be **bi-exponential** with an effective/local half-life sweeping from roughly 2.5 days early to about 9 days by two weeks, bounded above by the ~7-day terminal value. A fixed single half-life anywhere in the 3.5 to 7 day range is not supported by any primary measurement and biases capacity high by the 1.3 to 2.0 factor the robustness review already quantified.
- Duration bound: the mouse biochemical response to a single dose persists about 3 to 4 weeks\[3\] (Baek 2024, Nat Commun 15:3804, 10.1038/s41467-024-47460-9, which dosed "0.5 mg/kg mRNA-3927 on Days 0 and 28" and states the biochemical response is "driven largely by the half-life of the protein(s) for which it encodes"), which sets a floor on how slowly Δ(t) may decay and is consistent with the q2w/q3w clinical schedule accumulating rather than fully washing out between doses.

**Cross-programme support for a long, protein-limited tail.** The same platform's other assets confirm that the biomarker/effect time course is governed by the encoded protein's half-life, not the mRNA's: Baek 2024 states the mouse PA response lasted 3 to 4 weeks while the PKU (PAH) response lasted under 168 h, and attributes the difference to the mitochondrial PCC complex being protected from proteolysis.\[3\] For OTC, hepatic mRNA half-life is 7 to 9 h in mouse while OTC protein is retained about 10 days (mouse) or 8 days (rat),\[4\] and in humans ARCT-810 shows a plasma mRNA half-life of about 52 h but a pharmacodynamic effect on ammonia lasting about 28 days, per the Arcturus UCDC Sept 2025 poster, verbatim: "The LNP was rapidly cleared (t1/2 = 8 h) and mRNA was detectable through one week (t1/2 = 52 h)," and the Arcturus 30 June 2025 release, in which participants "maintained ammonia levels within the normal range following at least two doses and remained stable for approximately 28 days after completion of dosing." Every one of these reinforces that Δ(t) is a slow, protein-paced, bi-exponential-like decay, and that the mRNA PK half-life (hours to a couple of days) is the wrong timescale to use.

**Species caveat, restated per the repo's standing error.** All shape-defining numbers above (7-day protein half-life, kdeg, kq, ke0, IC50, 3 to 4 week duration) are **mouse** parameters. Attarwala's allometric scaling covers the PK (volume and clearance) terms only; the PD/protein terms were fixed from the Pcca⁻/⁻ (A138T) mouse and assumed transferable ("the concentration of PCC protein needed for 50% suppression of plasma 2-MC assumed to be similar in mice and humans").\[1\] They must not be read as human hepatic values.

## What the sponsor holds but has not published

- **Human hepatic PCC protein or activity time course.** Not measured at all; there are no liver biopsies in the trial. This is the single quantity that would fix the shape of Δ(t) in humans, and it does not exist publicly or, apparently, at all.
- **Numeric human PCCA mRNA PK parameters.** Extended Data Table 2 (Cmax, AUCinf, terminal T1/2, Tmax, CL, Vz per cohort)\[5\] is published only as a table image; the numeric cells are not in the machine-readable article, the PDF text layer, or the NCBI BioC/Europe PMC full text, and were not reproduced in any secondary source located. The per-patient PK time series behind Figure 3 are held by the sponsor.
- **PCCB mRNA PK.** The paper states PCCB testing was "still ongoing";\[6\] no PCCB PK has been reported.
- **Repeated-dose and steady-state PK.** Only dose-1 PK is shown; multi-dose accumulation parameters are held.
- **Emax, AUEC, AUC_below_baseline and duration of response for 2-MC and 3-HP.** These are registered Part 1, 2 and 3 outcomes (NCT04159103)\[7\] and are exactly the effect-versus-time integrals that would constrain Δ(t)'s shape, but none have been published; the paper reports only a single "percentage change after the third dose" snapshot.
- **Extended Data Table 3 numeric values and their sampling time.** The per-biomarker percentage changes (3-HP, 2-MC, C3, n-PG) are an image table; the exact post-dose sampling time relative to dose 3 is not disclosed.
- **SM-86 (ionizable lipid) exposure.** "SM-86 Concentration After Single and Repeated Administrations" is a registered Part 1 PK outcome,\[8\] unpublished.
- **Per-patient anti-PEG / anti-PCC titre time courses and their quantitative effect on delivered increment.** Only a narrative for two participants is public; the immunogenicity population dataset is held.
- **Part 2 (pivotal) data.** Enrollment in Part 2 completed February 2026 (Moderna mRNA-3927 Program Update, 24 Feb 2026: "Part 2 of the Paramount trial (mRNA-3927P101)...has completed enrollment"), and Moderna/Recordati expect a potential readout by end of 2026 (Recordati press release, 29 Jan 2026: "the pivotal readout expected by the end of 2026"). Part 3 (infants) is still recruiting.

## Answer to question 3 (biomarker sampling timing)

**The Extended Data Table 3 values are a percentage change from a pre-dose baseline to a post-dose level measured "after the third dose"; the precise timing of that post-dose sample relative to dose 3 is not public.** The table is titled "Biomarker percentage change from baseline after the third dose."\[6\] The Methods define baseline and post-dose: "Secondary endpoints included changes in blood 2-MC and 3-HP levels from baseline (predose levels) to post-dose levels,"\[6\] and baseline biomarkers were "assessed" during the 48 to 72 h observation period before dosing. The Results say "Reductions from baseline in disease biomarkers, including 3-HP, 2-MC, propionylcarnitine (C3) and n-propionylglycine (n-PG), were observed after dose 3."\[6\] Nowhere in the table legend, the footnote (which defines only Q2W, Q3W and SD), the Methods "Assessments" section, or any figure note is a clock-time given for the post-dose draw. The only dose-relative timing stated anywhere in the paper is for the antibody assay, not the biomarkers: Koeberl et al. 2024 (Nature 628:872–877) state, verbatim, "Antibody titres increased by at least fourfold in the predose sample for dose 3 (52 days after the first dose), then gradually decreased and returned to baseline levels at month 12." No 2025 to 2026 follow-up publication or presentation located (ICIEM 2025 final Part 1 results, ACoP 2025 model poster, the February 2026 program update) discloses the biomarker sampling time either. I checked the Nature HTML, the Nature PDF text layer, the NCBI BioC full-text API for PMC11156579, and secondary reproductions; the numeric cells of Extended Data Tables 2 and 3 are published as images and the sampling time is genuinely absent from the public record. Treat it as NOT FOUND, and do not assume it is a trough.

## What remains unknown

- Whether the human hepatic PCC half-life matches the mouse ~7-day value, or the 5 to 7 day / bi-exponential range the robustness review tested. There is no human hepatic protein measurement, so the very quantity that determines whether θ is biased 1.3 to 2.0-fold cannot be checked against data.
- The exact within-interval shape of Δ(t) around a single dose in humans, because Emax/AUEC/duration-of-response and the biomarker sampling times are unpublished. Measuring amplitude around one dose does not resolve the shape (as the review already noted), and the public data do not even give clean single-dose amplitude at a known time.
- The size and time profile of any decline in the delivered increment across repeated doses from anti-PEG or anti-PCC antibodies. The trial reports anti-PEG positivity and a transient titre rise in isolated participants\[6\] but never links it to a measured fall in PCCA mRNA exposure or biomarker response; the general LNP literature documents an accelerated-blood-clearance mechanism on repeat PEGylated dosing\[9\] but provides no mRNA-3927-specific magnitude.
- Human PCCB mRNA PK and repeated-dose/steady-state PCCA PK.

## Pending-reference table (format of `refs/literature_pending.md`)

| key | citation | id (DOI / PMID / registry) | provides | species | status |
|---|---|---|---|---|---|
| pcc_protein_halflife_mouse | Attarwala H, Lumley M, Liang M, Ivaturi V, Senn J. Translational PK/PD Model for mRNA-3927. Nucleic Acid Ther. 2023;33(2):141-147. | 10.1089/nat.2022.0036 / PMID 36577040 | 2-compartment biexponential PCC protein model; ~7-day half-life; kdeg, kq, ke0, IC50 | mouse (PD), allometric PK | already in refs (pkpd_translational); confirm shape terms are mouse |
| moa_duration_mouse | Baek R, et al. Characterizing the mechanism of action for mRNA therapeutics for PA, MMA, PKU. Nat Commun. 2024;15:3804. | 10.1038/s41467-024-47460-9 / PMID 38714648 | 3 to 4 week PA biomarker duration; PCC activity; MMA/PKU protein time courses; same LNP | mouse | already in refs (mrna_moa); cite for duration bound |
| trial_pk_biomarker_human | Koeberl D, et al. Interim analyses of a first-in-human phase 1/2 mRNA trial for propionic acidaemia. Nature. 2024;628:872-877. | 10.1038/s41586-024-07266-7 | human PCCA mRNA PK (image table); biomarker % change after dose 3; anti-PEG/anti-PCC; dose train | human | already in refs (trial_phase12); numeric tables unextractable |
| poppk_acop2025 | Gao W, He T, Li H, Liang M, Luo J, et al. Model-Informed Dose Selection for mRNA-3927 in Participants with PA. ACoP 2025 (ISoP). | 10.70534/JVQU8233 | human popPK PCCA mRNA: CL 0.266 L/h, V 13.9 L;\[10\] indirect-inhibition Kin PD model; 0.6 mg/kg q2w\[10\] | human | NEW, add to pending |
| iciem2025_finalpart1 | Schulze A, Koeberl D, et al. mRNA-3927 for PA: Final Results From mRNA-3927-P101 Part 1. ICIEM 2025, Kyoto, 2-6 Sep 2025. | no DOI (Moderna conference deck) | final Part 1 dose-escalation cohorts and cumulative ongoing-participant data;\[11\] no hepatic protein | human | NEW, add to pending; check for numeric PK/PD at full poster |
| mma_3705_trial | A Study to Assess Safety, PK and PD of mRNA-3705 in Isolated MMA (Landmark). | NCT04899310 | human MMA analogue; Emax/AUEC/Cmax outcomes for hMUT mRNA and methylmalonic acid | human | NEW, add to pending |
| mma_3705_preclinical | Second-generation mRNA therapy (mRNA-3705) in two MMA mouse models. Mol Genet Metab. 2024. | S1096-7192(24)00444-X | hMUT protein 2.1 to 3.4x over mRNA-3704; sustained MMA reduction;\[12\] same LNP | mouse | NEW, add to pending |
| pku_3210_trial | A Dose-finding Study to Evaluate mRNA-3210 in PKU. | NCT06147856 | human PKU analogue, WITHDRAWN (0 enrolled); no human data | human (withdrawn) | NEW, add to pending as negative |
| arct810_otc_phase2 | Arcturus. Interim multiple-dose Phase 2 results, ARCT-810 in OTC deficiency; UCDC/ICIEM Sept 2025. | NCT05526066 / NCT04442347 | human LNP-mRNA analogue: lipid t½ 8 h, mRNA t½ 52 h, PD ~28 days; different LNP (LUNAR) | human | NEW, add to pending |
| otc_protein_mouse | Sasaki, et al. LNP-targeted mRNA for OTC deficiency model mice. Mol Ther Nucleic Acids. 2023. | PMC10372164 | hepatic hOTC mRNA t½ 7-9 h; OTC protein retained ~10 days (mouse), ~8 days (rat) | mouse / rat | NEW, add to pending |
| antipeg_abc_review | Anti-PEG antibodies and accelerated blood clearance of PEGylated LNPs (review literature). | e.g. PMID 34890719; PMC12388889 | mechanism of exposure decline on repeat PEGylated dosing; no mRNA-3927-specific magnitude | in vitro / animal / human | NEW, add to pending as background only |

## Sources

1. [Translational Pharmacokinetic/Pharmacodynamic Model for mRNA-3927, an Investigational Therapeutic for the Treatment of Propionic Acidemia - Husain Attarwala, Matthew Lumley, Min Liang, Vijay Ivaturi, Joe Senn, 2023](https://dx.doi.org/10.1089/nat.2022.0036)
2. [Translational Pharmacokinetic/Pharmacodynamic Model for mRNA-3927, an Investigational Therapeutic for the Treatment of Propionic Acidemia - PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10066765/)
3. [Characterizing the mechanism of action for mRNA therapeutics for the treatment of propionic acidemia, methylmalonic acidemia, and phenylketonuria](https://www.nature.com/articles/s41467-024-47460-9)
4. [Lipid nanoparticle-targeted mRNA formulation as a treatment for ornithine-transcarbamylase deficiency model mice - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10372164/)
5. [(PDF) Interim analyses of a first-in-human phase 1/2 mRNA trial for propionic acidaemia](https://www.researchgate.net/publication/379540161_Interim_analyses_of_a_first-in-human_phase_12_mRNA_trial_for_propionic_acidaemia)
6. [Interim analyses of a first-in-human phase 1/2 mRNA trial for propionic acidaemia](https://www.nature.com/articles/s41586-024-07266-7)
7. [MRNA-3927 in Propionic Acidemia - Clinical Trials Registry - ICH GCP](https://ichgcp.net/clinical-trials-registry/NCT04159103)
8. [Study Details](https://clinicaltrials.gov/study/NCT04159103?aggFilters=status%3Anot+rec&cond=Propionic+acidemia+OR+%22glycinemia%2C+ketotic%22+OR+%22+ketotic+glycinemia%22+OR+%22+ketotic+hyperglycinemia%22+OR+%22+ketotic+II+glycinemia%22+OR+%22+propionic+aciduria%22+OR+%22+propionyl-coa+carboxylase+deficiency%22&rank=1&viewType=Table)
9. [Anti-PEG Antibodies and Their Biological Impact on PEGylated Drugs: Challenges and Strategies for Optimization - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12388889/)
10. [Model-Informed Dose Selection for mRNA-3927 in Participants with Propionic Acidemia (PA)](https://www.scienceopen.com/hosted-document?doi=10.70534%2FJVQU8233)
11. [Moderna : ICIEM 2025 (Propionic Acidemia P101 Data)](https://www.marketscreener.com/news/moderna-iciem-2025-propionic-acidemia-p101-data-ce7d59dbda8dff26)
12. [Improved therapeutic efficacy in two mouse models of methylmalonic acidemia (MMA) using a second-generation mRNA therapy - ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S109671922400444X)
