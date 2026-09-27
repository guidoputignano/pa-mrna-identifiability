# What it would take to land this

*Asked for by the author before round 6. This is an assessment, not a proposal to apply
anything. Nothing under `paper/` is changed by this file.*

Five rounds of Phase 11 sit on a draft PR that has never merged, alongside unapplied text from
Phase 4 and Phase 10. This file says what is ready, what needs a decision, and what the
consolidation pass would actually contain. **Counted rather than estimated**, by unbounded search
over the three proposed-text files.

---

## 1. The inventory

| set | passages | macros used | status |
|---|---|---|---|
| `audit/PHASE4_PROPOSED_TEXT.md` | **15** | **0** | prose only |
| `audit/PHASE10_PROPOSED_TEXT.md` | **1** | **0** | prose only |
| `audit/PHASE11_PROPOSED_TEXT.md` | **11** proposals (plus 3 "not proposed" sections) | **29 distinct**, of which **25 do not exist** | blocked on wiring |

**So the split is clean and it is not what the volume suggests.** Sixteen of the twenty-seven
passages need **no new numbers at all**. Every passage that needs a number is in Phase 11, and
all of them are blocked by the same mechanical thing.

---

## 2. Ready as they stand: 16 passages, needing only approval of the wording

All of `PHASE4` §1 to §15 and `PHASE10` §1. They change prose, cite existing results, and
introduce no macro. Under section 0 the paper is frozen for prose, so the only gate is the
author's approval, passage by passage.

**One caveat, and it is real.** `PHASE4`'s passages were drafted against the paper **as it was
several phases ago**. Since then Phase 10 corrected three numbers in `paper/` and converted 22
bare literals to `\num{}` macros. **A line-and-context re-check is required before any of them is
applied**, because a proposal quoting a line number or a sentence that has moved will apply
cleanly and wrongly. That re-check is mechanical and is item 1 of the consolidation pass. It is
also exactly the class of error section 0's shallow-clone and bounded-search rules exist to catch:
stale evidence that looks current.

**What they contain, by target file:**

| target | passages |
|---|---|
| `06_excitation.tex` | PHASE4 §1, §2, §5, §6; PHASE10 §1 |
| `07_design.tex` | PHASE4 §3, §4, §7, §14 |
| `08_specification.tex` | PHASE4 §11, §12, §13, §15 |
| `A1_supplementary.tex` | PHASE4 §7, §9 |
| `src/common/params.py` | PHASE4 §8, §10 |

---

## 3. Blocked on a decision from the author

Two decisions gate text that is otherwise drafted. Two more gate nothing textual and should not
hold the pass.

| | decision | what it gates | can it wait? |
|---|---|---|---|
| **(g1)** | Should a decade-long record carry a drift term? | `PHASE11` §3 and §6, which both state Tier 1's missing margin against a drift of order 1% a year. The passages can be written either way but not both. | **No.** It is the cheapest finding of the five rounds and it is one word, yes or no. |
| **(g2)** | Which object does "the binding nuisance" mean? | `PHASE11` §7, the caveat attaching to every Fisher figure. Section 0 now records that the profile, the asymmetry and the robust-versus-sd gap are **one non-Gaussian posterior with three symptoms**, which makes this one decision rather than three. | **No**, and it is now cheaper than when it was raised. |
| **(g7)** | **Which `rho` is "the crisis regime" in the threshold sentence, and does `-0.39` still stand?** | `07b_decision.tex:70-71`. The paper names two regimes, `d6` has four rows, and `-0.39` does not reproduce from any of them at seed 11. Needs the Phase 2D-era code or the author. | **No.** It sits in the same edit as `PHASE11` §14 and the seed spreads. |
| (g5) | Does the model need the catabolic coupling? | Nothing in the proposed text. It gates *future* pricing of prescribed load reduction. | **Yes.** |
| (g6) | What is `q`? | Nothing. The author has already ruled it sizes the hazard's wording, not the protocol, since a metered load gives `rho = 1` either way. | **Yes**, and it should be researched in parallel. |

**(g3) and (g4) are already settled** and their passages, `PHASE11` §10 and §11, are drafted with
the author's own margin target of 0.17 on the profile.

---

## 4. Blocked on the wiring pass: 29 macros, 25 of which do not exist

Every `PHASE11` passage that carries a number is written against a `\num{}` macro that has to be
generated first. This is stated three times in that file and it is the single largest blocker.

**The four that exist:** `disease_factor`, `ke_needed_factor`, `sd_monthly`,
`tier1_visits_at_pub`.

**The 25 that do not**, grouped by the audit script that already computes each value:

| group | macros | owner |
|---|---|---|
| breath elasticity and pair counts | `breath_el_nominal`, `breath_el_disease`, `marker_el_nominal`, `marker_el_disease`, `breath_k_600`, `breath_k_300`, `breath_k_300_pub` | `u2`, `u11` |
| the `k_e` frontier | `ke_frontier_3000_pct`, `ke_frontier_600_pct`, `ke_frontier_300_pct`, `ke_frontier_600_factor`, `ke_frontier_300_factor`, `ke_prior_pct`, `ke_spec_3000_pct` | `u9`, `u13` |
| drift | `drift_cross_pct`, `drift_premium_lo`, `drift_premium_hi` | `u7` |
| the non-Gaussian posterior | `profile_nominal`, `profile_lo_nominal`, `profile_hi_nominal`, `profile_gap_tier1_pct`, `sim_spread_nominal`, `quad_width`, `object_spread_pct` | `u12` |
| the margin | `margin_target` | author's decision, `u11` |

**Plus four new ones** from this round's `PHASE11` §14: `rho_pole_q`, `rho_above_fail_pct`,
`rho_above_fail_lo`, `rho_above_fail_hi`.

**The work is not the arithmetic.** Every value is already computed by a script in `audit/`. The
work is promoting each into `results/numbers/*.tex` through a generator the build runs, in the
repository's own pattern (`%  generated by <script> -- do not edit`, then
`\expandafter\def\csname paper@num@<key>\endcsname{...}`), and adding those generators to
`make figures` so `tools/check_artefacts.py` validates them. Phase 10 already did this for five
generators, taking the guard from 14 reachable generators and 328 keys to 19 and 364, so **the
pattern and its cost are both known.**

---

## 5. New work this round created, which must land with §14

`audit/u19_seed_defect.py`, recorded at R134. `07b_decision.tex` reports **fourteen prose
literals**, lines 58 to 61 and 70 to 71, in a section whose other numbers are macro-wired, and
**every one is a single seed**. Section 0 requires a spread with every number, and this is the
decision section.

The twelve at lines 58 to 61, regenerated at twelve seeds:

| `rho` | paper | mean | min | max | max/min |
|---|---|---|---|---|---|
| 1.00 | **1.030** | 0.884 | **0.289** | 1.079 | **3.73x** |
| 0.50 | **0.900** | 0.686 | **0.000** | 0.968 | **7119x** |
| 0.00 | 0.086 | 0.039 | 0.009 | 0.141 | 15.4x |
| -0.25 | 0.009 | 0.007 | 0.005 | 0.011 | 2.47x |

**All four published values are better than their own twelve-seed mean**, meaning the estimator
looks more accurate in print than it is across seeds. The `rho = 0.5` cell is the worst: the
published 0.900 hides a seed on which the estimate **collapses to 0.000**. The error-rate columns
are by contrast stable, 19.1% [16.8, 24.5] and 41.6% [37.5, 44.5] at `rho = 1`.

**And the two thresholds at lines 70 to 71 are worse than single-seed: one of them does not
reproduce.** The owner is **not** `d6_regimes.py`, which a first attempt assumed. An unbounded
search finds the literal at `audit/PHASE2D_FINDINGS.md:150`, so **Phase 2D is the source and `d6`
is a later re-derivation**. Re-deriving from `d6` at seed 11, cost ratio 1:1:

| regime | **1:1** | 10:1 | 30:1 | 100:1 |
|---|---|---|---|---|
| managed, `rho = 1` | **0.067** | escalate all | escalate all | escalate all |
| partial, `rho = 0.5` | -0.467 | escalate all | escalate all | escalate all |
| crisis pinned, `rho = 0` | *escalate none* | escalate all | escalate all | escalate all |
| crisis stopped, `rho = -0.25` | *escalate none* | **-0.400** | **-0.400** | **-0.400** |

**The number is there and it is in the wrong column for the sentence.** `-0.400` sits at
`rho = -0.25` and at **10:1 and worse**, within one grid step of the published `-0.39`. The
sentence says *"at equal costs ... the crisis regime wants it at $-0.39\theta$; at a cost ratio of
ten to one or worse every regime wants every patient escalated"*, and **both clauses fail in a way
that looks like a crossing**: at equal costs both crisis rows want *escalate none*, and at 10:1 the
protein-stopped row wants `-0.400` rather than "escalate all", which is the very row the `-0.39`
came from.

Still **not asserted as a paper error**, since the Phase 2D-era code may differ, but the question
is now a two-line one for the author: **is "at equal costs" an error for "at ten to one or worse",
and does the second clause need its exception?** Details and the correction to my own first reading
at R135 and R136.

**Consequence for the pass, and it is now larger than a wiring job.** `PHASE11` §14, the `rho`
pole, attaches to this paragraph. Applying it alone would add a fifteenth unwired number to the
one paragraph in the section with this problem. **§14, the twelve seed spreads, and the
resolution of the crisis threshold are one edit or none of them.** The third of those needs
either the Phase 2D-era code or the author, so **this item now has a decision inside it** and
belongs in section 3's list rather than section 4's.

---

## 6. What the consolidation pass would contain

In this order, because each step's output is the next step's input.

1. **Staleness re-check** of all 16 prose passages against the current `paper/`, line by line.
   Mechanical, no judgement, and it must precede any application. Output: each passage marked
   applies-cleanly, needs-rebasing, or superseded.
2. **The two decisions**, (g1) and (g2). One word and one choice of object. Everything in
   `PHASE11` §3, §6 and §7 waits on them.
3. **The wiring pass**: 29 generators' worth of keys into `results/numbers/`, added to
   `make figures`, validated by `check_artefacts.py`. Largest single item, entirely mechanical,
   and it unblocks all 11 `PHASE11` passages at once.
4. **The `07b_decision` seed fix**: wire the fourteen literals and give each its spread. Lands
   with `PHASE11` §14.
5. **Apply**, in dependency order: the 16 prose passages, then the 11 `PHASE11` passages, then
   §14 with the seed fix.
6. **One full build and guard run**: `make figures`, `make test`, `tools/check_artefacts.py`, and
   a document compile. The compile is the one step this container cannot do, having no TeX
   distribution, and the author has done it before for Phase 10.

**What is NOT in the pass, deliberately.** Nothing about the `k_e` measurement routes, the
glucose arm's precision numbers, the sampling modality, or the frontier's threshold onsets. Those
are audit findings that inform a future study design and belong in `OPEN_QUESTIONS.md` and the
`PHASE11_*` documents until a real assay and a real value of `q` exist. **Five rounds produced
eight findings documents and eleven proposed passages; the ratio is the point, not a failure.**

---

## 7. The honest assessment

**The blocker is not judgement and not volume. It is one mechanical pass and two one-line
decisions.** Sixteen passages need only a staleness check and approval. Eleven need 29 macros
that are already computed and merely unpublished. Two decisions gate three of those eleven.

**The risk in landing is concentrated in step 1**, the staleness re-check, because a Phase 4
proposal that quotes a line which has since moved will apply cleanly and wrongly, and nothing in
the output will look malformed. That is the same failure mode as the shallow clone and the
bounded search, and it is the third instance of the pattern, so it should be run as a check with
an explicit pass or fail per passage rather than as a reading.

**The risk in NOT landing** is that `PROPOSED_TEXT` has now grown for five rounds against a
moving paper, and every round makes step 1 more expensive.
