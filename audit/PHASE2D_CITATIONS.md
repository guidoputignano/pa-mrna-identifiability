# 2D.7 — citation sweep. Results and Tables only, never abstracts.

Two standing rules, adopted after three instances of the same failure:

1. **Every quantitative claim is checked against the source's Results or Tables**,
   never its abstract.
2. **Every row records which section of the source the value came from.**

Sources swept so far: `trial_phase12`, `natural_history`, `mde_claims_us`,
`breath_test_mma`, `amd_phenotyping`. These carry 16 of the 22 quantitative claims in
`paper/claims.csv`.

---

## Non-reconciling: 4 items

### N1 — `C16` / `trial_phase12`. A missing confidence interval that spans 1.

Register says: *"Phase 1/2 reports a 70 percent reduction in decompensation risk,
single-arm, computed within the eight participants who reported an event in the
pretreatment year."* Sourced to the abstract.

Results say:

> *"Across all cohorts in the combined dose-optimization and extension studies,
> participants had an overall relative risk of MDEs of **0.30 (95% confidence interval
> (CI), 0.066–1.315)** during the treatment period compared to the pretreatment period."*

Two corrections. The relative risk is reported **across all cohorts**, not within the eight.
And **the interval includes 1**: the headline efficacy signal is not statistically
significant, which the register's phrasing conceals. The biweekly subgroup is worse still —
RR 0.18 (95% CI 0.011–2.907; **P = 0.1832**).

This matters to this note's argument rather than merely to its bookkeeping. The paper asks
what a record can deliver about a per-patient decision; the trial's own event-based
efficacy estimate cannot exclude no effect on 16 patients. That is an argument *for* a
mechanistic per-patient estimand, and it should be made from the interval, not from "70%".

### N2 — `C05`, `C06`, `C07` / `mde_claims_us`. Abstract-rounded values, missing CIs.

| register | Results |
|---|---|
| 0.5 PPY overall | **0.44 (0.38–0.50)** |
| 1.4 PPY at 3–6 y | **1.38 (1.09–1.67)**, N = 32 pairs |
| 0.1 PPY at 13–17 y | **0.14 (0.03–0.26)**, N = 24 pairs |
| 0.2 PPY in adults | **0.22 (0.16–0.27)**, N = 106 pairs |

The **0.1 and 0.2 figures appear only in the abstract**; Results give 0.14 and 0.22.

And a substantive consequence the register does not carry: the 13–17 interval
[0.03, 0.26] **overlaps** the adult interval [0.16, 0.27]. C07 asserts *"the age profile is
U-shaped"*. The trough and the adult rise are not separated by these data, on 24 and 106
pairs. The U-shape is supported at the top (3–6 y at 1.38 is clearly above everything else)
and not at the bottom.

### N3 — `C48` / `breath_test_mma`. A standard deviation reported as a mean, by the source.

Already recorded in `REVIEW_LOG.md` R4. Abstract: *"mean coefficient of variation 6.9% and
12.8% respectively."* Results: controls *"average of 6.94 ± 5.3%"*, patients *"ranged from
0.39% to 50.84% with an average of **14.02 ± 12.86%**"*. The abstract's 12.8 is the SD. The
figure is also between-visit over 2 months to 4 years, not test–retest. The published
Correction restores a graphical abstract and changes no numbers.

### N4 — `C37` / `amd_phenotyping`. A column mislabelled by the source.

Haijes 2020 Table 3 is headed `R2` while the Methods specify Spearman correlation tests,
and the column contains negative entries. These are Spearman rho, not $R^2$. The register
inherits the source's mislabel. Plasma 3-HP against ammonia is rho = 0.66, n = 46.

---

## Reconciles cleanly against Results/Tables: 12 items

- `C08`, `C43`, `C49`, `C19` / `natural_history` — 49/50 (98.0%) and mean annualised 1.0
  (1.6); age strata 2.6 (2.7), 2.4 (4.2), 0.2 (0.4); early- vs late-onset 1.2 (1.8) vs 0.3
  (0.4); RR for a 50% lower 3-HP 0.7 (0.6–0.9, p<0.01), C3 0.6 (0.4–0.9, p=0.017), 2-MC 1.0
  (0.7–1.3, p=0.84). All in Results, *"Analysis of MDEs"* and *"Correlation between
  biomarkers and MDEs"*.
- `C36`, `C50` / `trial_phase12` — anti-PEG fourfold at the pre-dose sample for dose 3 (day
  52), back to baseline at month 12, grade 1–2 IRRs on 11/43 (25.6%) doses; the four
  biomarkers reduced after dose 3, quantified in **Extended Data Table 3**. Results.
- `C11`, `C12`, `C26`, `C29` / `saudi84`; `C01` / `pcc_structure`; `C51` /
  `carnitine_excretion_pa` — unchanged, and already sourced to Results in the register.

## A finding about the register itself

**The "three times higher" claim in `natural_history` sits in the Discussion, not the
Results**, and cites reference [20] = Ryu 2012, which contains no quantitative 3-HP value.
A register that recorded the source *section* would have caught this at entry: a Discussion
sentence carrying someone else's number is the highest-risk row type there is, and the
sweep should treat it as such.

## Not yet swept

`pathophys_complications` (C02, C03 — already `matches-with-drift`, full text not
accessible), `dose_selection_poster` / `pkpd_translational` (C46 — already
`matches-with-drift`), the three newborn-screening sources (C10 — counts already recorded
and the statistical point already made), `mde_emr_us` (C31), `genereviews_pa` (C29).
