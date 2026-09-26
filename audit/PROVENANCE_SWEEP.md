# The numeric provenance sweep

**What this is.** Every numeric literal in `paper/` that is not inside a `\num{}` macro, classified,
with category 3 listed. Commissioned after R75 found that route 6's `2.659`, `thirteen`, `six` and
`eighteen` are hardcoded prose rather than generated macros, so the artefact guard never compared
them with anything and one of them had been wrong since it was written.

Nothing in `paper/` was touched. The output is the list.

Scripts: `audit/q17_numeric_provenance.py` (scan), `audit/q17_classify.py` (classify),
`audit/q17_owners.py` (does a generator exist), `audit/q17_reproduce.py` (does it reproduce).
Data: `audit/q17_literals.txt`, `audit/q17_category1.txt`, `audit/q17_category3.txt`,
`audit/q17_claims.tsv`, `audit/q17_verdicts.tsv`.

---

## 0. What the guard actually sees

`tools/check_artefacts.py` regenerates the tables, the figures **and** `results/numbers/*.tex`, then
compares each with what is committed. So **a number written `\num{name}` is owned by a named
generator, runs in the build, and is checked.** A number typed into the prose is owned by nobody and
is compared with nothing, ever. That is the whole of the hole, and it is not small:

| | count |
|:--|---:|
| `\num{}` references in `paper/` (guarded) | **343** |
| digit literals outside `\num{}` | 1472 |
| spelled-out numerals (`thirteen`, `eighteen`) | 1000 |
| **unguarded literals in all** | **2472** |

**Spelled-out numerals are in the census because route 6's two wrong-or-unchecked values are words.**
A scan over digits would have missed the finding that prompted this sweep. That is the sharpest
single lesson here: the hole is wider than a regex over digits, and any guard built from this list
has to read English numerals.

---

## 1. The four categories

| | category | count | what it is |
|---|:--|---:|:--|
| 1 | GENERATED-ELSEWHERE, retyped | **218** | the value also exists in `results/`, so a generator owns the quantity and the prose retypes it by hand |
| 2 | STRUCTURAL | **1017** | equation constants, TikZ drawing coordinates, cross-references, years, counts of sections |
| 3 | **ASSERTING-COMPUTED** | **206** | prose-only and stating a computed outcome. **The list the author asked for** |
| 4 | HARMLESS | **1031** | prose-only, asserting nothing computed |

**Category 1 deserves a line of its own even though it was not asked for.** 218 literals restate a
value a generator already owns. Each is a silent copy that can drift from its source with nothing
noticing, which is R23's failure mode again, and unlike category 3 **every one of them has an obvious
fix**: replace the literal with the `\num{}` macro that already exists. It is the cheapest
improvement available to the consolidation pass.

**How categories 3 and 4 were separated**, since the boundary is the judgement in this exercise. A
literal is category 3 only if it is (a) not inside a `tikzpicture`, (b) not a hyphenated compound
(`eight-fold`, `six-monthly`, `one-dimensional`), (c) in a sentence carrying a result verb (buys,
beats, reaches, leaves, moves, costs, is worth, achieves), and (d) shaped like a computed output: a
decimal with two or more places, a percentage, an explicit factor, or an integer attached to a
countable design quantity. The thresholds are deliberately strict, because an over-inclusive list
buries the members that matter. **The cost is that the list is a floor, not a ceiling.**

---

## 2. Category 3: 206 literals in 130 claims

A claim is one line of prose; several literals in one sentence are one claim. Full rows in
`audit/q17_verdicts.tsv`.

### Does a generator exist that could own it?

The criterion is checkable rather than interpretive: **a generator exists if the same paragraph
already cites at least one `\num{}` macro.** That macro is written by a named script that already
runs in `make figures` and already emits numbers into that paragraph, so adding one more is an edit
to an existing generator, not a new one.

| | claims |
|:--|---:|
| a generator already serves the paragraph | **52** |
| no generator reaches the paragraph at all | **78** |

Which generator is already present, where one is:

| generator | claims |
|:--|---:|
| `results/numbers/design.tex` | 32 |
| `results/numbers/transient.tex` | 12 |
| `results/numbers/scope.tex` | 11 |
| `structure`, `fibre`, `lattice`, `phaseclock`, `decision`, `cohort` | 11 |

### Does the value still reproduce?

A category-3 literal is by definition not generated, so there is nothing to diff it against. What
can be done is to run the generators and search their printed output: **these scripts print far more
than they write as macros, so a literal appearing in its own generator's stdout is a value the code
still produces and merely never wired to a macro.** Fourteen generators were run for this.

| verdict | claims | meaning |
|:--|---:|:--|
| REPRODUCES | **26** | printed by a generator serving that section |
| NO OWNER, printed somewhere | **53** | not owned, but the value appears in the generated corpus |
| **NOT IN OUTPUT** | **18** | a generator serves the paragraph, ran, and did not print it |
| **PARTIAL** | **8** | same, for some of the sentence's numbers |
| **NO OWNER, UNPRINTED** | **25** | no generator reaches it and no generator prints it |

**NOT IN OUTPUT is not a verdict of error.** The generator may simply never print that quantity. It
is a verdict that **nothing in the build checks it**, which is the question asked. Route 6's
`eighteen` is exactly this: NOT IN OUTPUT, and then shown wrong by direct recomputation in q16 part
E. So the 51 claims in the bottom three rows are the queue, not the indictment.

---

## 3. The queue: 26 claims whose own generator ran and did not print them

Densest in Section 11's six-route enumeration, which is the section R75 already retracted from.

| location | values | verdict |
|:--|:--|:--|
| `11_scope.tex:201` | 0.3 | NOT IN OUTPUT |
| `11_scope.tex:204` | 7, 32, **2.32** | PARTIAL |
| `11_scope.tex:207` | **1.81, 1.22, 1.46** | NOT IN OUTPUT |
| `11_scope.tex:210` | **2.659** | PARTIAL |
| `11_scope.tex:213` | **thirteen** | NOT IN OUTPUT |
| `11_scope.tex:215` | **2.659** | PARTIAL |
| `11_scope.tex:216` | **thirteen, eighteen** | NOT IN OUTPUT |
| `11_scope.tex:105` | fifteen | NOT IN OUTPUT |
| `11_scope.tex:166` | 2.96 | NOT IN OUTPUT |
| `11_scope.tex:172` | 2.75, 3.07 | NOT IN OUTPUT |
| `07_design.tex:113` | **0.457, 0.412, 0.510** | NOT IN OUTPUT |
| `07_design.tex:181` | 0.999 | NOT IN OUTPUT |
| `07_design.tex:329` | 6.5, 87, thirteen | PARTIAL |
| `07_design.tex:333` | 29, seven | PARTIAL |
| `07_design.tex:613` | eight | NOT IN OUTPUT |
| `08_specification.tex:308` | 41, 12, four, five | PARTIAL |
| `09_validation.tex:103` | 0.170 | NOT IN OUTPUT |
| `06_excitation.tex:192` | 5.5 | NOT IN OUTPUT |
| `06_excitation.tex:214` | 0.55 | NOT IN OUTPUT |
| `06_excitation.tex:301` | 0.5, 2, 0.86 | PARTIAL |
| `A1_open.tex:139` | 1.79, 1.41 | NOT IN OUTPUT |
| `A1_open.tex:141` | 1.91, 2.02 | NOT IN OUTPUT |
| `A1_open.tex:175` | 1.79, 0.20 | PARTIAL |
| `A1_supplementary.tex:725` | six | NOT IN OUTPUT |
| `00_abstract.tex:24` | Fifteen | NOT IN OUTPUT |
| `01_question.tex:46` | 17, 95 | NOT IN OUTPUT |

**`07_design.tex:113` is the one to look at first after the route list.** It is the gate claim, "the
gate passes, and it passes on ordinary outpatient care", and its three numbers `0.457, 0.412, 0.510`
are a median and an interval that no generator printed. A result of that standing should not be a
literal.

## 4. And 25 claims with no owner and no match anywhere

| location | values |
|:--|:--|
| `08_specification.tex:255` | 0.283, 0.172, 0.170, 98 |
| `08_specification.tex:270` | **2.659** |
| `08_specification.tex:126` | 0.212, 0.20 |
| `08_specification.tex:361` | 120, 200 |
| `08_specification.tex:372` | 2.32, 0.20 |
| `08_specification.tex:289` | 167 |
| `04_local.tex:525`, `:526`, `:690`, `:714` | 0.170; 1.91, 2.02; 0.170; 169 |
| `03_dynamic.tex:129`, `:294`, `:352`, `:474` | 0.35; 10, 18; 1.1, 5, 4.8, 7; 1.4 |
| `02_model.tex:165` | 0.8, 0.5, 0.25, 1.9, 4.8 |
| `02_system.tex:126` | 50, 0.6, 0.017 |
| `06_excitation.tex:103`, `:148` | 0.587, 0.209; 0.058 |
| `07b_decision.tex:70` | 0.07 |
| `10_results.tex:216` | 1.1, 10, 5, 4.8, 10, 7 |
| `11_scope.tex:84` | 30, 88 |
| `A1_supplementary.tex:317`, `:390`, `:538`, `:705` | 25, 0.0370; 1.91, 2.02; 81; 0.8, 255, 1225, forty |

**`2.659` appears in three places across two sections** (`08_specification.tex:270`,
`11_scope.tex:210`, `:215`) and is owned by nothing in any of them. One value, hand-copied three
times, underpinning both route 5 and route 6.

---

## 5. Two clusters checked by direct recomputation

The stdout method answers "is it checked". For the two highest-stakes clusters the stronger question
was answered by recomputing the quantity.

### `08_specification.tex:255`, the load-record decomposition: REPRODUCES

The paper's third-tightest requirement rests on it, and it has no generator. Recomputed:

| quantity | paper | recomputed |
|:--|---:|---:|
| baseline, both channels at 5% | 0.170 | 0.172 |
| per-visit scatter alone, 20% | 0.283 | 0.293 |
| systematic offset alone, 20% | 0.172 | 0.174 |
| scatter / offset share of the degradation | 98% / 2% | 96.8% / 1.6% |
| offset's share at nominal / half / a fifth | 1.8 / 9.9 / 26.2% | 1.6 / 10.1 / 27.9% |
| scatter threshold with the offset held at 5% | 10.4% | 9.73% |

**Every figure reproduces to within a few per cent**, consistent with a different seed set or draw
convention rather than an error. **So the hole is real and this cluster is sound**, which is worth
saying as loudly as the retraction: an unguarded number is unchecked, not wrong.

### `11_scope.tex:213-216`, route 6: two clauses hold, one does not

Carried in from R75, and restated here because it is the worked example the rest of the list is
calibrated against.

- **thirteen measured-load visits match twenty-five recorded, at nominal**: reproduces, 12 here.
- **six reach the target**: reproduces exactly.
- **`2.659` at 5% of nominal**: reproduces, 2.655 at 1% accuracy.
- **"the advantage shrinks... from thirteen visits to eighteen"**: **does not reproduce.** Swept from
  nominal to 2% of nominal the matching count runs **12, 11, 10, 12, 12, 12**. The advantage is flat
  and is strongest at 20% of nominal. Route 6's conclusion, that the route opens nothing below
  nominal, is separately true, so a supporting number is wrong and the finding it supports is not.
- **And route 6's prose describes a different computation from the one that produced its numbers.**
  It says measuring the load "removes `b` from the parameter vector rather than bounding it", which
  is the b-only reading and gives **24 visits against 25, not 13**. The numbers come from moving
  `LOAD_FRAC`, which drives both channels. **Sections 7 and 11 contradict each other on the
  mechanism, and Section 7 is the correct one**: "what the load record has to deliver is per-visit
  precision entering through `sigma` and not a calibrated mean offset".

---

## 6. What the sweep says about the machinery, rather than about any number

1. **The guard's coverage is 343 against 2472, and its blind spot is not random.** It sees exactly
   the numbers somebody already decided to generate. Everything a sentence asserts in passing is
   invisible, and sentences are where conclusions live.
2. **206 literals in 130 claims assert a computed result that nothing checks.** 51 of them are not
   printed by any generator that serves them.
3. **78 of the 130 sit in paragraphs no generator reaches at all.** For those, closing the hole is
   not wiring but deciding which script should own the quantity, which is design work.
4. **218 further literals restate values a generator already owns**, and those are pure upside: the
   macro exists, the prose just does not use it.
5. **A guard built from this would have to read English numerals**, which the current one does not
   and no digit-based check ever would.

---

## 7. The repository has already solved this once, for figures

**This is the most useful thing the sweep found, and it is not a number.** `check_artefacts.py`
already carries a provenance mechanism for numbers **drawn on figures**, built after the phase
overlay printed "25 phases" on three panels whose lattices are 35, 28 and 28: the header was
`min(q, N)` with `N` the visit count, so three panels printed their own sample size back at the
reader as a property of the schedule, and nothing caught it because the number was computed inside
the figure module at draw time.

Its design is exactly what prose needs, and it already states why:

> The machine can find every number drawn at a call site and it can verify that an expression reads
> the register. It CANNOT tell a claim from an axis annotation ... So the rule is enforced as a
> DECLARATION.

Each figure module carries a `DRAWN` dict mapping every interpolation to one of `registered:`,
`own:`, `input:`, `geometry:`, `text:` or `UNREGISTERED:`. A number drawn and not declared fails; a
declaration for a number no longer drawn fails, so the list cannot rot into a blanket exemption; and
**`UNREGISTERED:` is a claim with no generator yet, counted and printed on every run, so the debt is
a number that can be driven to zero rather than a silence.** The current run reports 29 drawn
numbers, 13 registered, **0 unregistered**.

**So the answer to "what should own the 130 claims" already exists in this repository, applied to a
different surface.** The four categories in this document map onto that vocabulary directly:
category 1 is `registered:` waiting to be wired, category 2 is `geometry:` and `text:`, category 3 is
`UNREGISTERED:`, and category 4 needs no declaration. **The 130-claim list is this project's prose
debt stated in the form its own figure guard already uses**, and the 51 in the queue are what
`UNREGISTERED:` would print on the first run.

The one thing prose needs that figures did not: **the declaration has to cover English numerals**,
because `thirteen` and `eighteen` carried route 6's claim and no interpolation-based scan sees them.

---

**The limit of this sweep, stated so it is not over-read.** The category 3 / 4 boundary is
heuristic and tuned to be strict, so the list is a floor. Multi-line sentences are split across
claims. And "REPRODUCES" means a generator printed the value, which is weaker than recomputing the
quantity the prose asserts; only the two clusters in section 5 were checked that way.
