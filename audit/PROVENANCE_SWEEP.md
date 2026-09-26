# The numeric provenance sweep

## What this document found, in one paragraph

**Across 130 claims this sweep found no wrong number in the paper by checking, and five wrong numbers
in the checker.** Every "does not reproduce" it produced was mine: a schedule built by hand instead
of read from the table that owns it, a t = 0 convention, a capacity sweep stopped at 5% of nominal
twice, and a reconstruction of route 6 that concluded a correct clause was false. Route 6 is right in
all three clauses. **The defect this document exists to report is entirely in the machinery**, not in
the paper's arithmetic: 343 numbers are generated and checked, 2472 are not, and five generators that
the build never runs own 21 of the paper's sentences.

**One amendment, made after the above was drafted and honest about the order.** The author asked for
that paragraph and then directed the wiring of those five generators into the build. **The wiring
immediately found two genuine paper-against-generator disagreements**, `t3`'s "7% and 32%" against a
computed 6.3% and 24.1%, and `t4`'s "0.288" against a computed 0.286 (R88). So the accurate statement
is narrower and more interesting than the one drafted: **checking the paper against reconstructions
found nothing wrong with it five times over, and connecting it to its own generators found two things
wrong immediately.** That is the argument for the machinery, not against the paper, and it is why
wiring came before rewording.

---

**What this is.** Every numeric literal in `paper/` that is not inside a `\num{}` macro, classified,
with category 3 listed. Commissioned after R75 found that route 6's `2.659`, `thirteen`, `six` and
`eighteen` are hardcoded prose rather than generated macros, so the artefact guard never compared
them with anything and one of them had been wrong since it was written.

Nothing in `paper/` was touched. The output is the list.

Scripts: `audit/q17_numeric_provenance.py` (scan), `audit/q17_classify.py` (classify),
`audit/q17_owners.py` (does a generator exist), `audit/q17_reproduce.py` (does it reproduce),
`audit/q18_section8.py` (the convention, settled), `audit/q19_loadbearing.py` (the load-bearing
subset, verified), `audit/q20_rounding.py` (the rounding hole in the classifier).

**Sections 8 to 10 are later than sections 0 to 7 and correct two things in them.** They are appended
rather than merged, because what the first pass got wrong is part of the result.
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

---

## 8. The rounding hole, and the claim this sweep named wrongly

**`07_design.tex:113` was named in section 3 as "the one to look at first", on the ground that its
`0.457 (0.412 to 0.510)` was "a median and an interval that no generator printed". That was wrong.**
The committed table carries

    six-monthly clinic visits & 4 & ... & 0.4570 [0.4121, 0.5097]

and `results/numbers/design.tex` carries the same three values. The prose rounds them, `0.4570` to
`0.457` and `0.5097` to `0.510`, and the classifier matched literal **strings**, so every rounded
restatement of a generated value fell through into category 3.

`audit/q20_rounding.py` redoes the match with rounding tolerance, claiming a match only at **three or
more significant digits**, below which `0.45` meeting `0.5` is a coincidence and not a restatement.

| | before | after |
|:--|---:|---:|
| category 1, generated and retyped | 218 | **249** |
| category 3, asserting-computed | 206 | **175** |

**31 literals move.** Among them `09_validation.tex:103`'s `0.170`, which `\num{sd_monthly}` already
holds at `0.1699`; two of the three numbers in Section 8's decomposition, `0.170` and `0.172`; and
`11_scope.tex:172`'s `3.07` and `:207`'s `1.46`.

**The correction matters in both directions.** It shrinks the unowned pile by 15% and it grows the
cheap-to-fix pile, and it means the sweep's own headline overstated the danger of the claim it chose
to lead with. That claim reproduces, is guarded, and its defect is duplication.

---

## 9. The convention on Section 8, settled rather than described

Section 5 reported the load decomposition as reproducing "within a few per cent, consistent with a
seed convention", with the requirement threshold 7% apart. **That is not a verdict.** The convention
is identified here.

**It is the t = 0 convention, again.** `n4_theta_scope.MONTHLY` is `arange(0.0, 730.0, 30.4)`: 25
draws **starting at t = 0**, one before the first infusion. The first check used
`arange(30.4, 30.4 * 26, 30.4)`, dropping that draw and running 30 days past the window. On the
paper's own grid `n4.over_seeds(theta=3000)` returns **0.1699**, the paper's 0.170 exactly.

Redone on the paper's grid, seeds and median convention (`audit/q18_section8.py`):

| quantity | paper | first check | on the paper's convention |
|:--|---:|---:|---:|
| baseline, both channels at 5% | 0.170 | 0.172 | **0.170** (−0.0%) |
| per-visit scatter alone, 20% | 0.283 | 0.293 | **0.284** (+0.5%) |
| systematic offset alone, 20% | 0.172 | 0.174 | **0.172** (+0.1%) |
| scatter / offset share | 98 / 2% | 96.8 / 1.6% | **97.3 / 1.9%** |
| offset share, nominal / half / a fifth | 1.8 / 9.9 / 26.2% | 1.6 / 10.1 / 27.9% | **1.9 / 9.9 / 25.1%** |
| **scatter threshold, the stated requirement** | **10.4%** | 9.73% (−7%) | **10.27% (−1.3%)** |
| target met at every offset 1% to 500% | yes | not tested | **yes** |

**The verdict is REPRODUCES, not "consistent".** The threshold closed from 7% apart to 1.3%, which is
the residual of a bisection on a median over 24 seeds. The one cell still moving is the offset's
share at a fifth of nominal, 25.1 against 26.2.

**So the finding stands and is now properly earned: Section 8's decomposition is unowned by any
generator and correct.** Unguarded means unchecked, not wrong, and a sweep that only reported its
failures would misdescribe the repository.

---

## 10. The load-bearing subset, verified

`audit/q19_loadbearing.py`. Scope set by the author: any queued claim stating **a precision, a visit
count, a requirement threshold, or a pass-or-fail verdict**, since those are what the consolidation
pass will restate and lean on. No wiring, because wiring prose immediately before a pass that will
rewrite that prose is work thrown away. Everything on the convention pinned in section 9.

| claim | paper | here | verdict |
|:--|:--|:--|:--|
| `07_design:113` six-monthly median and interval | 0.457 (0.412, 0.510) | generated, all three | **REPRODUCES**, and is category 1 |
| two-monthly / fortnightly visit counts | 13 / 53 | 13 / 53 | **REPRODUCES** |
| `2.659`, the perfect-record limit at 5% of nominal | 2.659 | 2.668 (+0.4%) | **REPRODUCES** |
| route 6: thirteen measured visits match 25, nominal | 13 | 12 | **REPRODUCES** |
| route 6: six reach the target, nominal | 6 | 5 | **REPRODUCES** |
| **route 6: the advantage shrinks to eighteen visits** | **18** | **8 to 12, flat** | **DOES NOT REPRODUCE** |
| route 5: factor at nominal | 1.81 | 1.809 | **REPRODUCES** |
| route 5: disease-range band | 1.22 to 1.46 | 1.22 at ~25%, 1.457 at 2% | **REPRODUCES** |
| `09_validation:103` | 0.170 | 0.1699 | **REPRODUCES**, and is category 1 |

**`2.659` is asserted at three sites across two sections and every one of them is right.** The limit
is approached from above, 3.067 at a 5% record, 2.695 at 1%, 2.668 at 0.01%. The defect is that one
value is hand-copied three times with nothing owning it.

**Route 6's third clause remains the only thing in this entire sweep that does not reproduce**, and
it is now confirmed on the paper's own grid rather than a hand-rolled one: the matching visit count
runs 12, 8, 8, 10, 10 from nominal to 5% of nominal, flat and never approaching 18.

### What cannot be recomputed here, and why

| location | claim | why |
|:--|:--|:--|
| `11_scope:201` | 0.3% / 88%, increment CV 30% to exact | `t6_cohort_with_kappa.py` has the machinery and runs, but prints no such pair; the paper's sweep is not the script's |
| `11_scope:204`, `08_spec:372` | 7% / 32% / 2.32, prescribed diet swing | `t3_diet_modulation.py` exists and is **not in `make figures` at all** |
| `08_spec:126` | 0.212 at a 14.0% CV, 118 draws against 11 | same increment-CV machinery |
| `08_spec:308` | 41 / 12 / four / five, antibiotic weeks | descriptive, so outside this subset by the author's rule |

**The finding under those four is worth more than any of their values.** `t3_diet_modulation.py`
appears nowhere in `make figures`. The guard reports "14 numbers-file generator(s) checked for
reachability from `make figures`", so **a generator that writes no numbers file is not reachable, not
checked, and its results enter the paper only as literals.** A generator existing is not the same as
a generator running, and the reachability check already in this repository is the thing that would
say which.

### The tally

**Ten verdicts recorded: nine reproduce, one does not, and five further claims cannot be recomputed
with a stated reason.** The one that does not is route 6's third clause, already retracted in R75.
(The table above groups `07_design:113`'s three values and the two visit counts into single rows;
`audit/q19_loadbearing.py` records them separately, which is where the count of ten comes from.)

**And twice in this exercise a "does not reproduce" was mine and not the paper's**: the six-monthly
schedule, where I built the grid by hand instead of using the table's, and route 5's band, where I
swept 5% to 20% and the paper's upper endpoint sits at 2%. Both were caught by asking where the
paper's own value comes from before concluding it is wrong. That is the habit this list should be
read with.

---

## 11. The declaration is not a mechanism, and route 6 was right all along

`audit/q21_unreachable.py`. R80 reported `t3_diet_modulation.py` as one instance of a generator the
build never runs whose numbers are in the paper as literals. **It is one of four, and finding the
other three overturns R75's retraction.**

### 11.1 The 12 scripts that declare they do not reach the paper

Of 111 scripts in `audit/`, **12 declare in their docstrings that nothing of theirs reaches the
paper** ("Scoping only. Nothing here reaches the paper."). They split exactly in two:

| | in `make figures` | writes a numbers file | owns a paper claim |
|:--|:--|:--|:--|
| `t1_transient_worth`, `t5_band_vs_ke`, `t6_cohort_with_kappa`, `t7_structure_transient`, `t11_phase_lattice` | **yes** | **yes** | no |
| `t10_vandermonde`, `t12_triples` | no | no | no |
| **`t2_tier1_margin`** | no | no | **YES** |
| **`t3_diet_modulation`, `t4_bias_vs_scatter`, `t8_load_by_capacity`, `t9_measured_load`** | **no** | **no** | **YES** |

**CORRECTION, made after the adversarial workflow finished.** An earlier version of this section put
`t2_tier1_margin` in the middle row as genuine scoping work. **It is not: it owns the whole of Tier
2's requirement cluster**, `08_specification.tex:125-151`, verified here against its captured output:

| the paper says | `t2` prints |
|:--|:--|
| "at a 14.0% CV the record returns 0.212 against a gate of 0.20" | `14.02%   0.2120 [0.2108, 0.2135]` |
| "reaches the gate only at 118 monthly draws, nearly ten years" | `visits needed: 118 ... 9.8 years` |
| "against eleven with the increment known exactly" | `exact  0.1696  ...  11` |
| the floor row "0.151, 0.160, 0.174, 0.207, 0.248" | `0.1510 | 0.1603 | 0.1740 | ... | 0.2482` |
| "the floor crosses the target at an increment CV of 18.7%" | `the floor reaches the gate at an increment CV of 18.7%` |

**So it is five scripts, not four, and the second load-bearing requirement in the specification is
owned the same way the first is.** Section 8's load record comes from `t4` and `t8`; Tier 2's
increment repeatability comes from `t2`. Neither generator is run by the build.

This also closes a row section 10 left open: `08_specification.tex:126` was listed there as "cannot
be recomputed, same increment-CV machinery as `11_scope:201`". The machinery is `t2`, not `t6`, and
`t2` prints the answer.

**The declaration is wrong in both directions.** Five scripts say nothing of theirs reaches the paper
while being wired into `make figures` and writing macros the paper cites, which is harmless because
those numbers are guarded. Four say the same thing while owning prose the paper asserts, which is not.
**A docstring is not a mechanism.** What the build runs is decided by `make figures` and what the
paper cites is decided by `\num{}`, and neither reads English.

### 11.2 The four that own paper claims

None writes a `results/numbers` file, so none is among the 14 the guard checks for reachability. The
guard is not failing; it was never pointed at this class.

**`t4_bias_vs_scatter.py` owns Section 8's load decomposition**, which carries the paper's
third-tightest requirement:

| the paper says | the script prints |
|:--|:--|
| "the per-visit scatter alone gives 0.283" | `at 20% in the SCATTER alone   0.2829` |
| "the systematic offset alone gives 0.172" | `at 20% in the BIAS alone      0.1718` |
| "against a 0.170 baseline" | `bias to 5%:  0.1696` |
| "98% of the degradation is scatter" | `the scatter channel accounts for 0.1133  (98%)` |

**And `t4`'s decomposition is not stated once but restated across sections, all of it unowned.**
`06_excitation.tex:286-290` carries the same result again: "runs 0.159 at a 1% recording error to
0.288 at 20%", "98% of that is the per-visit scatter and 2% is the systematic offset", and "with the
scatter at 5% the target is met at every bias from 1% to 500%". `t4` prints `1%  0.1589`,
`20%  0.2863`, and, verbatim in its own output, `the gate holds at EVERY bias from 1% to 500%
(0.1691 at 1%, 0.1733 at 500%)`. **So one unreachable generator supplies the same finding to two
sections, nine literals in all, and no macro cites any of them.**

**`t8_load_by_capacity.py` owns route 5's factors**: `100% of nom  1.808` against the paper's 1.81,
and `20% of nom 1.215 ... 2% of nom 1.464` against its "1.22 to 1.46 across the disease range". So the
band runs from a fifth of nominal down to two per cent, which section 10 inferred and this confirms
from the owning generator.

**And `t8` owns a second Section 8 sentence, which resolves the last residual in section 9 above.**
Section 9 verified the decomposition against the paper and left one cell moving, the offset's share
at a fifth of nominal, 25.1 here against the paper's 26.2. `t8` prints

```
       theta  5/5 baseline  offset at 20%  scatter at 20%  offset share
 100% of nom        0.1699         0.1720          0.2842          1.8%
  50% of nom        0.3228         0.3308          0.3960          9.9%
  20% of nom        0.8104         0.8408          0.8958         26.2%
```

**1.8, 9.9 and 26.2, exactly the paper's three figures.** So that residual was not seed granularity
either; it was my reconstruction again, and the owner matches the paper to every digit.

**`2.659` is printed by BOTH `t8` and `t9`**, at 2.6595. One value, two unreachable generators, three
sites in the paper, no macro.

**`t9_measured_load.py` owns route 6 entirely**, and this is the finding that matters:

```
        theta  recorded n=25   measured n to match   measured n to gate
  100% of nom         0.1699                    13                    6
   50% of nom         0.3228                    12                > 120
   20% of nom         0.8104                    13                > 120
   10% of nom         1.5828                    15                > 120
    5% of nom         3.0642                    17                > 120
    2% of nom         7.4406                    18                > 120
```
and, at 5% of nominal, `to 0.1%, n=25:  2.6595`.

**A fifth site, and it reconciles R80's two figures.** `06_excitation.tex:309` says "measuring the
load instead of recording it is worth under one per cent". `t9` prints, at nominal,
`recorded 5%  0.1699` against `b gone, sig 5%  0.1693`, which is 0.35%. **So the paper separately and
correctly prices the b-removal channel ALONE at under one per cent**, which is the same channel R80's
first pass measured at 0.3 to 1.5% and then treated as an error. Both figures were right: the paper
states the b-only channel here and the full accuracy move in route 6, and the two differ by a factor
of twenty because they are different channels. A spelled-out numeral again, so no digit scan sees it.

### 11.3 R75's retraction is withdrawn: route 6 is correct in every clause

| route 6's clause | R75 | the owning generator |
|:--|:--|:--|
| thirteen measured visits match twenty-five recorded | reproduces (12) | **13, exactly** |
| six reach the target | reproduces exactly | **6, exactly** |
| `2.659` at 5% of nominal | reproduces | **2.6595, exactly** |
| **"the advantage shrinks ... from thirteen visits to eighteen"** | **RETRACTED as wrong** | **13 at nominal, 18 at 2% of nominal. Correct.** |

**The clause is right and the retraction was wrong.** The count rises monotonically as capacity falls,
13, 12, 13, 15, 17, 18, which is exactly "shrinks rather than reverses". R75 swept only to 5% of
nominal and built the schedule by hand; the paper's range runs to 2%.

**And R75's second charge against route 6 also falls.** R75 said its prose "describes a different
computation from the one that produced its numbers", because "removes `b` from the parameter vector"
is the b-only reading. `t9`'s own docstring is explicit: **"THREE THINGS CHANGE AT ONCE, and they are
not independent. (1) b disappears... (2) sigma stops being a diet record's accuracy and becomes a
tracer measurement's precision. (3) The steady record closes on itself."** The prose names the first
of three and omits the second. **That is an incomplete description, not a wrong one**, and Sections 7
and 11 do not contradict each other: Section 7 says `sigma` carries the requirement, route 6 says `b`
is removed, and measuring the load does both.

**What survives from R75 is the provenance defect and nothing else.** Route 6's numbers are correct,
computed by a real generator, and stated in the paper as literals that no macro cites and no build
step regenerates. That was always the finding worth having.

### 11.4 Four times, the same error

R77 reported `2.659` as "unowned and unprinted", meaning no generator printed it. `t9` prints it at
four significant figures. R77's reproduction pass ran fourteen generators and `t9` was not among them,
because `t9` is not in `make figures` **and that is precisely the defect being reported**. The method
inherited the very gap it was measuring.

That is the fourth time in this phase that a reconstruction of mine disagreed with the paper and the
paper was right:

1. the six-monthly schedule, built by hand instead of taken from `tab_D_designs.tex`
2. Section 8's decomposition, on the wrong t = 0 convention
3. route 5's band, swept to 5% when its endpoint sits at 2%
4. route 6's eighteen visits, swept to 5% when its endpoint sits at 2%, and on a hand-rolled grid

**Every one has the same cause: recomputing a quantity instead of finding the generator that owns it.**
The discipline that would have prevented all four is the one this sweep exists to enable, which is
worth stating plainly: **find the owner before recomputing, and if no owner can be found, say so
before concluding the number is wrong.** The sweep's own list is the tool for finding the owner, and I
built it and then did not use it four times.
