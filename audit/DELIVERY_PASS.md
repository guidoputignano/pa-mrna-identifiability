# The delivery pass: the paper as an outside reader meets it

*Read from `paper/main.pdf` alone, 45 pages, 28,226 words, built at `fd02ea8` with 0 LaTeX errors and
0 unresolved references. `audit/`, the supplement and every `.py` file were closed for this pass, per
the author's instruction. Nothing here is checked against the code, because the point of the pass is
what the document says by itself.*

---

## 1. The single claim an outside reader would say the paper makes

> **The mRNA dosing schedule is already an identification experiment, and here are the rules for
> reading it. The quantity it was meant to identify cannot be recovered in the patients who have the
> disease.**

That is one claim with a hinge in the middle, and **the paper knows it.** §8.1 opens: *"The
contribution is the design rules, with a scope condition attached, and not the capacity estimate."*

**But the title, the running head and the first paragraph of the abstract all sell the first half.**
"Estimating Residual Enzymatic Capacity in Propionic Acidemia" names the thing the paper concludes
cannot be estimated in propionic acidemia. The abstract opens on 42 times the excitation and reaches
the factor of 16 below the recoverable band in its third paragraph, after the design rules.

A reader who stops at the title takes away the opposite of the result. A reader who finishes takes
away that the method is sound and the application is out of range. **The gap between those two
readers is the paper's central delivery problem**, and it is a problem of order rather than of
content: every qualification is present, correct and late.

---

## 2. Three reasons a referee would reject it

### (a) Every precision is computed at a capacity no patient has

The paper reports `sd(ln θ̂) = 0.1699` and builds a three-tier specification around a gate of 0.20.
**Both are at θ = 3000, which §7.1 states is normal capacity.** Residual PCC activity in propionic
acidemia is under 5% of control, so the patient sits at θ ≤ 150.

§7.1 says this plainly: *"the illustration cannot be read as a description of the patient."* **It says
it on page 34 of 45**, after four sections of design work, a sampling schedule, an admissibility
analysis and a decision estimand, all evaluated at the wrong operating point.

A referee reads that as the paper demonstrating its apparatus on a case chosen because it works. The
defence is in the text, and it arrives after the reader has formed the objection.

### (b) The decision that motivates the paper is not decided by the measurement, and the paper concedes it

The stated motivation is a dosing interval. §5 is titled "The decision estimand" and its own summary
is: *"Everything above says the therapy channel makes θ identifiable; it does not follow that it makes
the sign of the margin decidable, and it does not."*

§8.1 goes further: near normal capacity *"the decision is blocked not by identifiability but by its own
boundary, which is stipulated rather than measured and which no better estimate of θ and no larger
cohort can repair; in the disease range both block."*

**So the measurement does not decide the decision at any capacity**: at normal capacity the threshold
is missing, in the disease range the estimate is. A referee will ask what the paper is for, and the
honest answer, which §8.1 gives, is that it is for the design rules. That answer is inconsistent with
the framing everywhere above §8.

### (c) The specification's inputs are assumed, and the paper says so in three separate places

- `k_e`: *"the one quantity this design asks for and cannot cite"* (§6.3). Its prior width of 0.8 is
  listed in Table 2 as *"chosen for illustration."* §6.5: the band's **location** depends on `k_e`'s
  value, not only its width, and quartering it moves the band to 24%–647% while quadrupling it means
  the gate is never met at any capacity.
- Tier 2, the same-day repeatability that everything the therapy channel buys is conditional on,
  **has never been measured**, and the nearest published figure is between-visit, in the wrong
  disease, and on the wrong side of the requirement.
- The cohort sweep of §6.6 runs over *"every range assumed and none measured."*

Each disclosure is honest. **Their sum is a quantitative specification whose every input is a
literature guess or an illustration**, presented in tables with four significant figures. A referee
will say the precision of the presentation is not earned by the precision of the inputs.

---

## 3. What I went looking for and could not find

**Two candidate absences were searched unbounded across the full text with no result limit, and the
first one was mine rather than the paper's.**

**Not a gap, and I would have reported it as one.** I looked for how `k_e` would actually be obtained,
expecting the paper to name it as binding and then leave it. **§6.3 specifies the study**: a micro-dose
intravenous bolus of ¹³C₃-3-hydroxypropionate with serial plasma sampling and timed urine collection,
with the note that the obstacle is GMP re-synthesis and ethics rather than chemistry. The absence
claim was wrong and is not made.

**The real gap: there is not one observation from one patient anywhere in the paper.**
Searched unbounded: `real patient` 0, `observed data` 0, `clinical data` 0, `held-out` 0, `we fit` 0.
Every number is a forward computation on an assumed model with literature constants, or a simulation
over assumed ranges. §6.1 concedes it exactly: *"Exactly one row is a measurement rather than a
computation."*

**This is the thing a referee will ask for first and the paper cannot supply.** It is not a flaw in
the reasoning. It is that a paper claiming a clinical record suffices to identify a parameter has
never been shown a clinical record. The design rules are derived, the scope limit is derived, the
specification is derived, and the one empirical contact point, a between-visit repeatability in
methylmalonic acidemia, is cited to show that a requirement is **not** met.

**The corollary for delivery**: the paper's strongest move is to stop presenting as a measurement
study and present as what it is, a design and identifiability analysis that tells an ongoing trial
what to record. That is order B.

---

## 4. Two abstracts and two section orders, same content, no new claims

### A. As it stands: the negative result on capacity

**Order**

| § | title | what it does |
|---|---|---|
| 1 | Introduction | the question, and the answer previewed |
| 2 | Model and methods | the dynamic model, the equilibrium reduction, the local model |
| 3 | The excitation condition | what a record must do for θ to be recoverable |
| 4 | The sampling design | the schedule, admissibility, sources of excitation |
| 5 | The decision estimand | **the decision is blocked** |
| 6 | The data specification | three tiers ordered by dependency |
| 7 | Limits of recoverability | **the scope limit, and every precision was at normal capacity** |
| 8 | Discussion | what the five results amount to |

**The shape**: method, then design, then two negative results at §5 and §7, then a discussion that
reinterprets everything above in their light. The reader meets the qualifications at pages 24 and 34
of 45.

**Abstract A** is the one in the document. Its logic, made explicit: the trial already runs the
experiment; here is what transfers; those rules hold at every capacity but the estimate does not;
the band is 82% to 627% of normal and the disease is a factor of 16 below it; what survives is `K_m`
and an upper bound; the paper ends in a three-tier specification.

**Its cost**: the concessive turn arrives in the third of four paragraphs, after the reader has been
given a factor of 42 and a factor of 27 to hold. The title contradicts the conclusion.

---

### B. Leading with what is identifiable, the design that gets it, and `k_e` as the quantity that decides it

**Order**

| § | title | drawn from |
|---|---|---|
| 1 | Introduction: what is recoverable from a clinical record, and what is not | current §1 + §8.1 |
| 2 | Where capacity is recoverable, and why the boundary is an identity | current §7 and §7.1, **moved to the front** |
| 3 | What is identifiable at every capacity | `K_m` to 6–7%, the upper bound on θ, from §6 and §7 |
| 4 | The excitation condition | current §3 |
| 5 | The sampling design | current §4 |
| 6 | `k_e`: the constant that places the band, and what would measure it | §6.3's second half, §6.5, Table 11, **gathered** |
| 7 | The data specification | current §6 |
| 8 | What the estimate decides, and what it does not | current §5 |
| 9 | Model and methods | current §2, **moved to the back** |
| 10 | Discussion | current §8 |

**The two structural moves**: the scope limit leads instead of arriving at page 34, so no precision is
ever quoted before the range it is valid in is stated; and `k_e`, currently split across §6.3, §6.5 and
Table 11, becomes one section, because it is the single constant that sets where the band sits, how
precisely θ can be estimated, and which of the two tiers is better value in the population the other
is hardest to run in.

**Abstract B**

> The clinical decision in mRNA enzyme replacement for propionic acidemia is a dosing interval, set by
> residual enzymatic capacity θ, which nobody measures. **Whether a clinical record can supply it is
> settled by an identity rather than by a study design.** What the record carries about `ln θ` is the
> flux through the enzyme, `∂s/∂ln θ = θ g(m) = L − k_e m`, so its sensitivity to capacity is
> capacity's own share of the load: 84.5% at normal, 29.2% at a tenth of normal, 16.3% at the 5% that
> defines the disease. **θ is recoverable between 82% and 627% of normal, and propionic acidemia sits
> a factor of 16 below that band.** Fifteen times the sampling moves the lower edge by 9%, from 82% to
> 75%: the limit is the identity, not the schedule.
>
> **What survives there is worth stating first.** `K_m` is determined to 6–7% at every capacity, the
> record still supplies an upper bound on θ, and the design rules hold throughout: carry the dosing
> record as a regressor rather than a nuisance; never sample at the dose interval or a multiple of it,
> which costs a factor of 27 against monthly draws and is recovered by moving the draw half a day;
> record intake as eaten and to better than about 10%, a requirement on the record's accuracy and not
> its variability; assume nothing about `K_m`. The therapy is a known input acting selectively on θ and
> supplies 42 times the excitation the marker's own band does, so the trial is already running the
> experiment.
>
> **One constant decides where the band sits, and nothing supplies it.** The non-enzymatic clearance
> `k_e` is assumed here at a width chosen for illustration. Quartering it moves the band to 24%–647%
> and quadrupling it means no capacity meets the gate; tightening its prior from the published width
> to 5% is worth a factor of 2.1 at normal capacity and 6.9 in the disease range, though even an exact
> `k_e` leaves 0.388 there against a target of 0.20. A micro-dose ¹³C₃-3-hydroxypropionate bolus with
> serial plasma and timed urine would return it in one study.
>
> **Identification is not what binds the decision in any case.** The threshold the decision turns on is
> stipulated rather than measured, and no better estimate of θ and no larger cohort repairs that. The
> paper therefore delivers a three-tier specification ordered by dependency: one blood draw three
> minutes after an infusion, which separates two published structures by 16.3 standard deviations and
> decides which model the other tiers calibrate; then the same-day repeatability of a ¹³C-propionate
> breath test, which nobody has published; then the record itself. **Every number here is a forward
> computation on literature constants, and exactly one row of the specification is a measurement.**

**What B buys**: the scope limit is the second sentence rather than the third paragraph, so no reader
carries a precision without its range. `k_e` is promoted from a table row to the quantity the result
turns on. The title's promise is replaced by the paper's actual contribution. And the last sentence
makes the absence of patient data the author's statement rather than the referee's.

**What B costs**: it opens on a negative and needs the "what survives" paragraph to carry the reader.
It also concedes the empirical gap in the abstract, which is a choice about how the paper wants to be
attacked and is the author's to make.

---

## 5. The one thing both orders need and neither currently has

**A sentence saying what a reader should do on Monday.** Both versions end in a three-tier
specification ordered by dependency, and neither says which tier the reader is expected to run, at what
cost, or what changes if it reports. Tier 1 is one blood draw and it gates everything else, which is
the strongest practical sentence in the document, and it is currently a subsection heading at §6.2.
