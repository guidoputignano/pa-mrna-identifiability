# Literature pending verification

**Nothing in this file is registered.** These references arrived from two literature reviews
(clinical and methods) supplied in September 2026. **None has been verified in this
repository**: network egress is blocked in this container, so no DOI here has been resolved,
no landing page opened, and no bibliographic field checked against a publisher record.

**The rule.** A reference moves from this file to `paper/refs.sources` only after its DOI has
been resolved somewhere with egress and the resolved record matches the citation. That is the
same standard `refs.sources` already applies, where the fields are fetched rather than typed.
Until then nothing here may be cited in `paper/`, and nothing here changes a constant or a
conclusion.

**How to read the confidence column.** It records what the *review* claimed, not what this
repository checked. "CONFIRMED (review)" means the review stated it had seen the DOI on a
publisher, PubMed or MathSciNet page. "Verify (review)" means the review itself flagged the
DOI as unchecked. Both still need verification here.

---

## Already held. Do not re-register these

| citekey | DOI | note |
|---|---|---|
| `natural_history` | `10.1186/s13023-026-04251-3` | The methods review lists this as new. **It is already in `refs.sources:36`** with exactly this DOI. No action. |
| `breath_test_mma` | `10.1038/s41436-021-01143-8` | Manoli et al. 2021, *Genet Med* 23(8):1522-1533, already at `refs.sources:97`. **See the conflict below.** |

### The breath-test conflict

The clinical review reports "CV 12.8% in MMA patients, 6.9% in controls". **`REVIEW_LOG.md`
R4 retracted that figure in this repository**, on two grounds:

1. Manoli's Results give **`14.02 +/- 12.86%`**, range `0.39` to `50.84`. The abstract reports
   `12.86`, the **standard deviation**, as though it were the mean. The review has inherited
   the source's own error.
2. The figure is **between-visit, over two months to four years**, with the authors naming
   changing renal status as a cause. It is therefore an **upper bound** on the same-day
   repeatability that Tier 1 of `08_specification.tex` asks for, not a measurement of it.

The review's other claim, that **no PA breath-test cohort exists**, agrees with the paper as
written and needs no change.

---

## Methods references

Each row records what the reference would support if verified. The mapping is this
repository's, not the review's.

| reference | DOI | confidence | what it would support here |
|---|---|---|---|
| Raue et al. 2009, *Bioinformatics* 25(15):1923 | `10.1093/bioinformatics/btp358` | CONFIRMED (review) | The profile-likelihood machinery of `j1_uncertainty.py`, and the structural-vs-practical distinction the flat `K_m` profile turns on |
| Bates & Watts 1980, *JRSS-B* 42(1):1 | `10.1111/j.2517-6161.1980.tb01094.x` | CONFIRMED (review) | The curvature-versus-tangent decomposition. **This is the classical statement of what `REVIEW_LOG.md` R14 retracted**: the paper's `eps^4` is the curvature of an exact fibre against its own tangent |
| Self & Liang 1987, *JASA* 82(398):605 | `10.1080/01621459.1987.10478472` | CONFIRMED (review) | The boundary LRT mixture. **This is the citation for why the boundary correction must NOT be applied here** (section 6): the mixture is for a parameter *being tested* on a boundary, and here that parameter is a nuisance |
| Andrews 2000, *Econometrica* 68(2):399 | `10.1111/1468-0262.00114` | CONFIRMED (review) | Bootstrap inconsistency on a boundary. **The explanation for the measured bootstrap failures** in section 6: percentile `0.77-0.82`, parametric at fitted nuisances `0.683` |
| Lavielle & Aarons 2016, *JPKPD* 43(1):111 | `10.1007/s10928-015-9459-4` | CONFIRMED (review) | Population-versus-individual identifiability. Supports section 6's "describes the prior's population, not one patient" |
| Janzen et al. 2016, *Front Physiol* 7:590 | `10.3389/fphys.2016.00590` | CONFIRMED (review) | Structural identifiability of the standard PD model families. The nearest prior art to Phase 5 steps 3 and 4 |
| Golub & Pereyra 2003, *Inverse Problems* 19(2):R1 | `10.1088/0266-5611/19/2/201` | CONFIRMED (review) | Variable projection, which is exactly `i1_estimator.py`'s estimator. The 1973 founding paper is `10.1137/0710036` |
| Nyberg & Karlsson 2009, *JPKPD* 36(2):125 | `10.1007/s10928-009-9114-z` | CONFIRMED (review) | Joint optimal design on dose and sample times: the design counterpart of the phase-spread schedule |
| Walter & Pronzato 1990, *Automatica* 26(2):195 | `10.1016/0005-1098(90)90116-Y` | Verify (review: "confidence high") | The qualitative design step, choosing an input/output configuration so that parameters are identifiable at all. The general form of the dose-phase result |
| Mehra 1974, *IEEE TAC* 19(6):753 | `10.1109/TAC.1974.1100701` | CONFIRMED (review) | Optimal input design. **Citation trap the review flags: this is 19(6), 753-768, and is often miscited as "Vol. 29" or "AC-21, 55-64"** |
| Segel & Slemrod 1989, *SIAM Rev* 31(3):446 | `10.1137/1031091` | Verify (review: "confidence high") | The validity condition for the quasi-steady-state assumption the whole balance rests on |

### Do not cite for the dose-phase result

**van Noort et al. 2024, *CPT Pharmacometrics Syst Pharmacol*, `10.1002/psp4.13148`.** It uses
"aliasing" for parameters with **similar sensitivity profiles**, not for sampling synchronous
with dosing. It is not support for the dose-phase result and must not be cited as such. The
review's own conclusion is that **no single source states the dose-synchronous claim
verbatim**, so it has to be assembled from the input-design references above, and the paper
should say so rather than imply a single citation exists.

---

## Clinical references

| reference | identifier | what it bears on |
|---|---|---|
| Daly et al. 2017, *Mol Genet Metab Rep* 13:83 | PMID 29021961 | That PA diet records are **prescriptions, not intake as eaten**. If it holds, the `10%` recorded load CV is unsupported and a prescription-based record is a **constant** recorded load, which is the exactly singular case |
| Thompson et al. 1990, *Metabolism* 39:1133 | n/a (6 patients, fasting, mixed PA/MMA) | That about **half of propionate production is non-dietary** (gut bacteria `>=22%`, odd-chain fatty acids `~30%`). Bears on the scale of `b` (prior `S_B = 16.14` against `L = 322.73`, `b_true = 0`) and on its constancy (metronidazole `~-40%`, illness, fasting) |
| Hsia et al. 1979 | PMID 481943 | That residual PCC activity in PA is typically **<5% of control**, which favours the 5% reference in open question (b). `params.py` states the fraction two ways four lines apart and this audit left the choice to the PI |
| Manoli et al. 2021, *Genet Med* 23(8):1522 | PMID 33820958 | Already held as `breath_test_mma`. **See the conflict above before using any number from it** |

**One review finding needs no citation because it is an absence.** No human data exist on
plasma 3-HP half-life, clearance or test-retest CV, and absolute plasma 3-HP in `umol/L` is
unpublished for PA. This **confirms** the `k_e` row's "priority gap" and adds a third gap on
the same quantity. It is recorded in `audit/PHASE4_SUMMARY.md` section 5a.

---

## What verification means for this file

For each row: resolve the DOI, confirm the landing page matches author, year, journal, volume
and pages, and confirm the quoted claim appears where the review says it does. Then move the
row into `paper/refs.sources` with an annotation in that file's own style, which records what
the source says rather than only that it exists. The two "Verify (review)" rows and the Mehra
volume number are the ones most likely to be wrong.
