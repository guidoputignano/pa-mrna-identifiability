# Open questions, and how to pick this up

Written at the end of the Phase 0 to Phase 5 audit, September 2026, as a handoff to a fresh
session. **Everything here is unresolved.** What is resolved lives in `audit/PHASE*.md` and
`REVIEW_LOG.md`.

**Start here:** `audit/PHASE4_SUMMARY.md` is the two-page consolidated summary and its section
0 is the cut to read first (structure-robust results against model-dependent ones). This file
is only the list of what is still open.

---

## 0. Standing rules, which the next session should inherit

These were set by the author and held throughout. They are not negotiable defaults.

- **The paper and supplement are frozen.** No prose edit until a change is reported and
  approved. `audit/PHASE4_PROPOSED_TEXT.md` holds ten proposed passages and **none is applied**.
- **No constant moves** without an explicit decision. `src/common/params.py` is the register.
- **No agent swarms and no scheduled check-ins.** Work solo and in the foreground.
- **Footing and seed spread with every number.** A single-seed figure quoted without its range
  has twice been the error here (`REVIEW_LOG.md` C5, and the stress test's own finding).
- **One independent re-check before each report**, sharing no code path with the thing it
  checks. Every `*_recheck.py` in `audit/` is one, and they have caught real errors in every
  phase; several are recorded in `REVIEW_LOG.md`.
- **Stop at gates.** Report and wait rather than proceeding.
- **Conventions:** audit scripts in `audit/`, findings in `audit/PHASE*.md`, constants
  registered in `src/common/params.py` with tests in `tests/test_models.py`, retractions in
  `REVIEW_LOG.md`, no expected result written into a script in advance.
- **Style:** the author asks for no em-dashes. Use a full stop to start a new sentence or a
  comma to continue one.

---

## 1. Decisions the author owes (from `PHASE4_SUMMARY.md` section 5)

| # | question | what turns on it |
|---|---|---|
| **(a)** | Does the exact constant-load curve go into the paper, and where? | Passages 1, 2 and 5 of `PHASE4_PROPOSED_TEXT.md` are drafted for it. This is the `R14` retraction reaching the paper. |
| **(b)** | **5% or 10% of normal capacity?** | `params.py` states the residual fraction two ways four lines apart. Changes no pass or fail; changes the admissible range of the unidentified direction by `2.5x`. **New, unverified:** Hsia et al. 1979 (PMID 481943) favours **5%**. |
| **(c)** | Should the phase clock be the main design figure? | `fig_E1` is drawn and unused. 14 d and 14.5 d are identical on a calendar axis and differ twentyfold in the answer. |
| **(d)** | Reword the `D` row, or leave it? | "falls by about half" is `62` to `66%` on every reading. Nothing downstream breaks if it stays. |
| **(e)** | **Where does dietary load enter a direct-response model?** | **The largest one open.** See section 2. |

### (e) in full, because it is the one a fresh session needs the context for

A marker that is an algebraic function of the enzyme has **no load input**, so any load must be
imposed. Attarwala's model was fitted in a fixed-diet animal colony and states no coupling.
The choice is not cosmetic. With `u = (L+b)/L0`:

| coupling | `theta` identifiable? | structural residual |
|---|---|---|
| `(E0 + base*f)*u`, the whole marker scales | **no** | `0.18 sigma`, inaudible |
| `E0 + base*u*f`, only the suppressible pool | **YES, rank 5 of 5** | `0.22-0.24 sigma`, inaudible |
| `E0*u + base*f`, only the unaffected floor | no | `0.95-1.15 sigma`, audible |
| `E0 + base*f + gamma*u`, a separate source | no | `0.93-1.13 sigma`, audible |

So the answer decides **both** whether `theta` is identifiable at all **and** whether a misfit
would be noticed. The audit's judgement, flagged as judgement and not a finding, is that
`E0 + base*u*f` is the most natural reading of Attarwala's own definitions (`base` is "amenable
to suppression", `E0` is "not affected by mRNA-3927"), so a larger substrate load should lift
the suppressible pool and not the unaffected floor. **Nothing in either source says this.**
Full working in `audit/k7_insertions.py`, re-checked in `audit/k8_recheck.py`.

---

## 2. Data that must come from outside the public record

| what | why it is needed | who has it |
|---|---|---|
| **When was Extended Data Table 3 sampled?** | The single largest swing in the design. Trough reading: `10` to `12` monthly visits to `sd(ln theta) = 0.20`. Peak reading: `42`, or not within five years. | **The sponsor.** The article names no timepoint and has no biomarker time-course figure; the registry time frame is "Baseline (predose levels) to postdose levels ... up to Week 40", a window. **The public record is exhausted.** The economical ask is **`Emax` and `AUEC` per cohort for 3-HP**, which the registry shows Moderna computes and which would pin the dose calibration without the timing question arising. |
| **A human dose-to-marker calibration** | `D_AMP = 6000` implies a `62-66%` fall against reported per-cohort falls of `-29%` to `+5%`. The constant is deliberately **not moved**; `audit/PHASE5_ROBUSTNESS.md` carries a ten-row sensitivity table with admissibility on every row in its place. | Same request as above. |
| **The enzyme half-life** | `lambda = ln 2/3.5 d` is annotated as corroborated by `kdeg = 0.0075/h`, which is one rate constant inside a two-compartment model whose own stated protein half-life is `~7 d`. Moving `lambda` costs `11 -> 30` monthly visits and moves every Phase 3 figure. Passage 10b states the qualification and keeps the constant. | Published, but the reading needs a decision. |
| **Whether PA diet records are intake or prescription** | **Possibly the most consequential of all.** If records are prescriptions, the recorded load is **constant**, which is the exactly singular case (rank 3 of 4). The whole design rests on the load varying. | Unverified literature: Daly et al. 2017, PMID 29021961. |
| **Baek's Source Data file** | Would let the translational model be simulated end to end instead of through our calibration. | Not yet in `refs/data/`. **See section 3.** |

---

## 3. The Baek Source Data file: a standing instruction

The author said they would download it into `refs/data/`. **It had not arrived when this
session ended.**

**Do not fit it.** `audit/k9_sourcedata_prior.py` was written **before** the file arrives,
deliberately, so it is a yardstick and not a rationalisation. It records, from the design
alone:

- It **can** separate a marker state from none on the **slow** calibration branch (shape gaps
  `80.2%` and `49.7%` against a best-case standard error of `20.9%`).
- It **cannot** on the fast branch (`16.1%` and `2.6%`; a 39-minute half-life on a 2-day grid).
  **So a null result is not evidence for the algebraic structure.**
- It **cannot identify `theta`** under the algebraic structure at all: a fixed-chow colony has a
  constant load, which is a rank deficit of two regardless of `n`.

**The first pass when it arrives is descriptive only:** per-animal series or per-timepoint
summaries (if summaries, most of the value is gone), how many animals and timepoints actually
carry a 3-HP value, whether the pre-treatment sample is present for every animal, and the units
and assay LOQ. **Then report and wait.**

---

## 4. Literature held unverified

`refs/literature_pending.md`. Twelve methods references and four clinical ones arrived in
September 2026 from two literature reviews. **None is registered in `paper/refs.sources` and
none has been verified here**, because network egress is blocked in this container. The file
records what each would support and the rule for moving a row into `refs.sources`.

**Two things a fresh session must not walk into:**

1. **The breath-test figure is already retracted here.** The clinical review gives "CV 12.8% in
   MMA patients, 6.9% in controls" (Manoli et al. 2021). `REVIEW_LOG.md` **R4** records that
   Manoli's Results give `14.02 +/- 12.86%` and the abstract reports the **standard deviation**
   as though it were the mean, and that the figure is **between-visit over two months to four
   years**, so it is an upper bound on the same-day quantity Tier 1 needs. **Do not import it.**
2. **van Noort et al. 2024 (`10.1002/psp4.13148`) is not support for the dose-phase result.** It
   uses "aliasing" for parameters with similar sensitivity profiles, not for sampling
   synchronous with dosing. The reviews' own conclusion is that **no single source states the
   dose-synchronous claim verbatim**, so the paper should say so rather than imply one exists.

---

## 5. Work deliberately not done

- **No test was added for the Phase 5 simulators.** The natural candidates are the two
  calibration identities and the closed form `theta = (V-W)/Imax`, all exact arithmetic. They
  belong in `tests/test_models.py` if and when the simulators become part of `src/` rather than
  of the audit. That is a decision, not an oversight.
- **`StructuralIdentifiability.jl` is not installed** in this container. It turned out not to
  matter, since every structure here reduces to a rational map of degree (1,1) and the
  reduction is verified three ways, but it would be worth running as a fourth opinion.
- **Nothing beyond Phase 5 step 4 was attempted.**

---

## 6. Retractions a fresh session should know about

Read `REVIEW_LOG.md` for the full entries. The four most recent:

| # | what was retracted or scoped |
|---|---|
| **R14** | "the degenerate combination is recovered through curvature, at a square-root-of-noise price". On the quasi-steady locus at a constant load the direction is **not recovered at all**: `Delta(m) = A/m + B + Cm` gives three invariants for four parameters, an exact one-parameter fibre. This one **is** in the paper's prose in three places and passages 1, 2 and 5 are drafted for it. |
| **R15** | "the published inhibition shape cannot reach the repository's calibrated fall". A knife edge on `Imax = 0.999`, a constant the source **fixed** rather than estimated. At `Imax = 1` the turnover shape clears the fall at every dose. Replaced by a better statement: the paper's fall is **bracketed** by the two published shapes' ceilings on every tabulated dose. |
| **R16** | Four Phase 5 statements that were true of the configuration measured and were reported as unconditional: `theta`'s non-identifiability under the algebraic structure (insertion-dependent), the bias figures (prior-dependent, `1.3-1.7` with the paper's priors and `5-9` without), the detection claim (a turnover property, not a design property), and the three-number bound (proved for a stated class, not general). |
| **R4** | The breath-test CV. See section 4. |

**The failure mode all four share** is reporting the scope of a measurement as the scope of a
claim. It is worth watching for.
