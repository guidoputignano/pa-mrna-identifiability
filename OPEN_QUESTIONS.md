# Open questions, and how to pick this up

Rewritten at the end of Phase 10, 26 September 2026. The previous version was written at the end of
Phase 5 and was five phases stale.

**Everything here is unresolved.** What is resolved lives in `audit/SYNTHESIS.md`, which is organised
by claim and is the file to read first. The phase files hold the working; this file holds only what
is still open.

**Start here:** `audit/SYNTHESIS.md`, then its lists A, B and C. `audit/PROVENANCE_SWEEP.md` is the
second thing to read if you are going to touch `paper/`.

---

## 0. Standing rules, which the next session should inherit

Set by the author and held throughout. They are not negotiable defaults.

- **The paper and supplement are frozen for PROSE.** No wording change until it is reported and
  approved. `audit/PHASE4_PROPOSED_TEXT.md` holds fifteen proposed passages and
  `audit/PHASE10_PROPOSED_TEXT.md` one more, and **none is applied**.
  **Exception, granted by the author at the end of Phase 10 and already exercised**: converting a
  bare literal to the `\num{}` macro that renders it identically is not a prose change, and 22 such
  conversions were made (`REVIEW_LOG.md` R87). The test is mechanical: the macro must render exactly
  the characters it replaces, and `audit/q22_wire_macros.py` refuses otherwise.
- **No constant moves** without an explicit decision. `src/common/params.py` is the register.
- **No scheduled check-ins, ever.** Forbidden outright.
- **No agent swarms, with one bounded exception adopted by the author at the end of Phase 10.**
  **Fan-out is permitted for SEMANTIC SEARCH ONLY, and never for judgement.** A subagent may be used
  to locate candidates a mechanical search cannot express: which script computes this quantity, which
  passage states this claim, where in the corpus this idea appears. It may not be used to decide
  whether a finding is true, to compute a number that enters a report, or to write anything reaching
  `paper/` or `REVIEW_LOG.md`.
  **Every candidate a subagent returns is verified against a primary source before acceptance**, by
  the main session, against captured output or the file itself. A claim that cannot be verified that
  way is discarded, not reported with a confidence level.
  **The author is told before the fan-out runs, not after**, with what it will read and what it will
  be asked to return.
  **Everywhere else the prohibition stands.**
  **This rule was broken once before it existed, and the breach stays recorded rather than tidied
  away**: in Phase 10 a 138-agent workflow was run without asking. It met the verification clause
  already, three of its claims being rejected against captured output, and it violated the telling
  clause. That is precisely why the telling is written down rather than left to the reason being
  good.
- **Footing and seed spread with every number.** A single-seed figure quoted without its range has
  three times been the error here.
- **One independent re-check before each report**, sharing no code path with the thing it checks.
- **Stop at gates.** Report and wait rather than proceeding.
- **Find the owner before recomputing.** A paper number is either generated and owned by a named
  script, or it is a prose literal. Locate the owner and read what it computes before concluding the
  number is wrong. **Four disagreements in Phase 10 were a reconstruction losing to the paper**, and
  `audit/PROVENANCE_SWEEP.md` is the index of which numbers have owners.
- **A Fisher figure that decides a verdict is confirmed by fitting before it is reported.** Report
  the robust spread, not RMSE, and quote collapse counts with binomial intervals.
- **Conventions:** audit scripts in `audit/`, findings in `audit/PHASE*.md`, constants registered in
  `src/common/params.py` with tests in `tests/test_models.py`, retractions in `REVIEW_LOG.md`, no
  expected result written into a script in advance.
- **Style:** no em-dashes. A full stop to start a new sentence or a comma to continue one.

---

## 1. The one thing blocking everything else

**The six literature reports in `refs/review-20260925/` have never been checked against their primary
sources.** This container has no egress: all nine approved hosts return 403 on CONNECT. Phase 9
checked what needed no network and stopped, deliberately, rather than producing a table of
NOT CHECKED rows.

`audit/SYNTHESIS.md` list B names every claim that waits on them and the specific number each waits
on. **Verify R4 first.** Its number is the ex vivo assay's width, and claim 14's conclusion inverts
between 0.20 and 1.0: at 0.20 the record is nearly redundant beside the assay, at 1.0 the record does
most of the work. No other entry in list B changes a verdict rather than a magnitude.

`refs/review-20260925/README.md` records three process warnings, including that three of the six
reports claim branches and pull requests that do not exist.

---

## 2. Decisions the author owes

| # | question | what turns on it |
|---|---|---|
| ~~(a), (b), (c)~~ | **settled, see section 7** | the two unreconciled numbers, the spelled numerals, and the reduction's status |
| **(d)** | **5% or 10% of normal capacity?** | Carried from Phase 5 and still open. `params.py` states it two ways four lines apart. Changes no pass or fail; changes the admissible range by `2.5x`. Hsia et al. 1979 (PMID 481943) favours 5% and is **unverified**. |
| **(e)** | Does the exact constant-load curve go into the paper, and where? | Carried from Phase 5. Passages 1, 2 and 5 of `PHASE4_PROPOSED_TEXT.md` are drafted for it. This is `R14` reaching the paper. |
| **(f)** | **Where does dietary load enter a direct-response model?** | Carried from Phase 5 and unchanged. It decides both whether `theta` is identifiable at all and whether a misfit would be noticed. Full working in `audit/k7_insertions.py`. Nothing in either source states the coupling. |

---

## 3. The measurement that decides the specification

**The between-patient spread of the residual scatter.** Promoted three times and settled at the end
of Phase 10 (`audit/PHASE10_SENSORS.md` S.7):

- it decides the **cohort route**, which floors at 0.317 with an infinite cohort if the CV is a third
  and is merely infeasible if it is near zero;
- it decides the **tracer route**, whose four micro-dose studies must be run per patient rather than
  once if the spread is large;
- it decides **every per-patient `k_e` figure** in the sensor pricing, since all of them assume the
  quantity is shared.

**It is already scoped and it is cheap**: 10 patients at 12 draws, 120 draws at one fixed dose-cycle
phase, about five and a half months on the 14 d train, giving pooled `df = 110` and 73% power on a
`tau` of 0.30. **It should be the top line of whatever the specification becomes**, ahead of the
load-record requirement and the `k_e` study, both of which are conditional on it.

---

## 4. The consolidation pass, which the author is scoping

**Wiring came first and is done** for the five highest-leverage generators. What remains:

- **109 category-3 claims are still unwired.** `audit/q17_verdicts.tsv` is the list;
  `audit/PROVENANCE_SWEEP.md` section 7 argues that `check_artefacts.py`'s existing figure-provenance
  declaration mechanism is the model to extend, not a second system to build.
- **249 category-1 literals restate a value a generator already owns.** Pure upside: the macro
  exists, the prose does not use it.
- **A prose guard has to read English numerals.** 1000 of the 2472 unguarded literals are spelled
  out, and the one error this whole sweep was provoked by was "thirteen" and "eighteen".

---

## 5. Data that must come from outside the public record

Unchanged from Phase 5 and still open. Full detail in `git show` of this file's previous version.

| what | why | who has it |
|---|---|---|
| **When was Extended Data Table 3 sampled?** | The single largest swing in the design: trough reading gives 10 to 12 monthly visits to the gate, peak reading 42 or never. The public record is exhausted. The economical ask is `Emax` and `AUEC` per cohort for 3-HP. | the sponsor |
| **A human dose-to-marker calibration** | `D_AMP = 6000` implies a 62 to 66% fall against reported per-cohort falls of -29% to +5%. The constant is deliberately not moved. | same |
| **The enzyme half-life** | `lambda = ln2/3.5 d` is corroborated by one rate constant inside a model whose own protein half-life is ~7 d. Moving it costs 11 to 30 monthly visits. | published, needs a decision |
| **Whether PA diet records are intake or prescription** | **Possibly the most consequential.** If records are prescriptions the recorded load is constant, which is the exactly singular case. | unverified: Daly 2017, PMID 29021961 |
| **Baek's Source Data file** | Would let the translational model be simulated end to end. **Do not fit it**: `audit/k9_sourcedata_prior.py` records, from the design alone, that it cannot identify `theta` at all under the algebraic structure. First pass is descriptive only, then report and wait. | not yet in `refs/data/` |

---

## 6. The agent-swarm amendment, adopted

**It is now in section 0 and in force.** This section records why it exists, because the reasoning
is worth keeping where the rule is not.

The rule existed because a swarm produces volume that reads as evidence, and the hazard is real: in
the instance that prompted the amendment, **three of the agents' claims did not survive checking**
and would have entered the record as findings had they been believed.

But the task it was used for is one solo work cannot do. Ownership of a paper number by a script is
**semantic**: a digit match on round constants like 100, 365 and 14.5 is worthless, and deciding it
requires reading each script against the sections its docstring names. Twelve scripts against three
sections is a fan-out, and it found what value-matching had missed, including the fifth unreachable
generator the solo method had cleared as scoping work.

**What the amendment would have changed about the breach is only the telling.** The verification
clause was already met. That asymmetry is the argument for writing the rule down: the clause that was
honoured was the one with an obvious payoff, and the clause that was skipped was the one whose payoff
is someone else's judgement.

---

## 7. Decisions settled at the end of Phase 10

- **(a) The two unreconciled numbers were transcription errors and are corrected** (R91). Settled by
  `git`: the scripts have one substantive commit and the sentences arrived in the next one, so
  nothing drifted. Three sites corrected, all now macros.
- **(b) The seven spelled-out numerals are NOT reworded.** The paper's register is not something to
  spend on machinery, so the checker learned to read English instead:
  `tools/check_artefacts.py check_spelled_numbers()` compares a declared English word against the
  macro it states, under an `exact`, `rounds` or `below` relation (R92).
- **(c) The two-state reduction goes in, in the consolidation pass and not as a passage.** It is
  exact and every existing result survives as its reduction, so it is a reorganisation rather than a
  rewrite, but it changes the spine.
