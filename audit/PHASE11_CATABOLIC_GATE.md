# Phase 11, round 4. A result of the audit's own is withdrawn, its substitute is priced, and the re-check found a defect that reversed the verdict

*Proposals only. Nothing under `paper/` is changed by this file.*

Scripts: `audit/u15_gut_suppression.py`, `audit/u16_recheck_gut.py`.
`REVIEW_LOG.md` R116 to R124.

**Read section 5 before quoting any number from sections 3 or 4.** The independent
re-check found a real defect in this round's own code, the correction moved
`sd(ln k_e)` by a factor of 2.6, and the verdict reversed. Every figure below is
post-correction.

---

## 1. The fast is withdrawn, and the paper had already priced it

Round 3 priced a **50% fast** at the trough as the safe lowering step for a relaxation
burst, and reported `sd(ln k_e) = 0.1121` at a tenth of nominal at a 1% assay and sixteen
admissions. **In propionic acidemia a fast is not a load reduction.** Acute restriction
triggers catabolism, which releases the precursor amino acids from muscle and raises the
load.

**The paper says so in its own words.** `01_biology.tex:151-156`: the dangerous events are
triggered by "infection, **fasting**, surgery, trauma", and "the patient stops eating, the
body enters catabolism, and skeletal muscle protein is degraded for energy... the propionate
load rises sharply precisely when the patient is eating nothing". `:160` gives the management
consequence, intravenous glucose to suppress catabolism rather than restricting protein.

**And `07b_decision.tex:55-62` has already priced a prescribed cut.** Writing the recorded
load as `L_base + rho (s_true - L_base)`, `rho <= 0` is "crisis protocol, where prescribed
intake is cut while catabolic muscle breakdown supplies the load". At `rho = -0.25`, protein
stopped, **the sign error is 100% and the median `theta_hat/theta` is 0.009**. The fast arm
assumed a genuine reduction that the record tracks, which is `rho = 1`.

> So the withdrawn numbers do not merely lack support. **They contradict an existing paper
> result about the same intervention**, and the existing result is the pessimistic one. The
> figures are struck through in `PHASE11_KE_MEASUREMENT.md` section 2.6 rather than deleted,
> so that the correction is visible to the next reader.

**This was the author's catch, not the audit's**, and verifying it made it stronger than
stated: the premise was that the model cannot see catabolism, and what the repository
actually shows is that one part of it can and already does.

---

## 2. (g5) Does the identification model need a catabolic coupling? A GATE, not a result

`eq:local` carries `L` as an input dialable in either direction. The catabolic response lives
in the prose, and as `rho` it lives in a **recording asymmetry**: `rho` says how faithfully
the *record* tracks the true load, not that a prescribed cut *causes* the true load to rise.
Those are different objects. A physiological coupling would make `L` a function of the
prescription rather than a free input.

**The audit's reading, which is not the decision:**

| | |
|---|---|
| **Not needed** for anything the paper currently claims | Its load variation is spontaneous day-to-day intake at 10%, and its decompensations are *already* modelled as load spikes, which is what a fast provokes. Nothing in print depends on `L` being dialable downward. |
| **Needed** before any *prescribed* reduction is priced | Round 3 priced one and got it wrong by two orders of magnitude on the paper's own `rho` result. `rho` is the hook a coupling would attach to. |

**Nothing that prices a prescribed load reduction should proceed until this is settled.** The
route in section 3 is a reduction of a different kind, and section 3.1 is the argument that
it is exempt.

---

## 3. The substitute, and what its amplitude rests on

### 3.1 Why gut-derived suppression is exempt from section 2's gate

The property needed is a load reduction that does **not** provoke catabolism, so that the
true load genuinely falls and the record tracks it, which is `rho = 1`. Suppressing
gut-derived propionate removes substrate **at source**: the patient keeps eating, no energy
deficit is created, no muscle is degraded. It is also routine management rather than a
protocol to be invented.

**This is a different use from the one the audit closed earlier, and the earlier negative
does not transfer.** `audit/u2_breath_functional.py` part E found a load step *worse* than the
dose cycle, for a **breath ratio**, because a ratio wants the widest marker-level contrast and
the dose train supplies a larger one. Here the step is not widening a contrast. It is the input
whose **relaxation** carries `k_e`, placed at the trough where `k_e` owns the relaxation rate
(R113: `k_e/J` is 0.552 at the trough against 0.076 just after a dose). Different quantity,
different design, and the earlier verdict says nothing about it.

### 3.2 Two registered numbers, and only one of them is an amplitude

| source | number | what it measures | usable as a load-step amplitude? |
|---|---|---|---|
| `src/common/params.py`, `s_base` row | gut flora **"at least 22.2 +/- 8.5%"** | a share of whole-body propionate **production**, from the fall under metronidazole | **yes.** A share of production is exactly what bounds how much of `L` can be removed. |
| `audit/m2_antibiotic_cycle.py`, Thompson 1990 (PMID 2226555) | **41%**, range 12 to 76 | the fall in urinary propionate **metabolites** on metronidazole | **no.** A fall in excretion is not a fall in production. |

> **The provenance is load-bearing here**, because the verdict differs between the two. The
> tables below sweep across both and mark which is which, rather than picking one.

### 3.3 The amplitude sweep, at an ideal step

24 hourly draws at the trough, Tier 1 record, one admission, 5% assay. Record-only floors
**0.6565** at `theta = 600` and **0.5430** at 300.

| amplitude | source | `theta = 600` | `theta = 300` |
|---|---|---|---|
| 10.0% | | 0.6452 | 0.5302 |
| 12.0% | Thompson range, low | 0.6403 | 0.5239 |
| 13.7% | that share minus one sd | 0.6354 | 0.5176 |
| **22.2%** | **`params.py` gut share of production** | **0.6029** | **0.4762** |
| 30.7% | that share plus one sd | 0.5623 | 0.4277 |
| 41.0% | Thompson metabolite fall, *not a share* | 0.5102 | 0.3711 |
| 50.0% | | 0.4660 | 0.3279 |
| 76.0% | Thompson range, high | 0.3586 | 0.2349 |

Monotone throughout. It was **not** monotone before the section 5 correction, which turned
over at 50%, and that non-monotonicity was an artefact of the defect rather than physics.

### 3.4 The onset sweep, at the registered amplitude

Metronidazole is dosed over one to two weeks and flora suppression develops over days, while
the marker relaxes in 2.65 h at a tenth of nominal. The burst spans four onset time constants
or 24 h, whichever is longer, at 24 draws.

| onset `tau` | burst span | `theta = 600` | `theta = 300` | |
|---|---|---|---|---|
| ideal step | 1.0 d | 0.6029 | 0.4762 | |
| 1 h | 1.0 d | 0.5959 | 0.4696 | |
| **3 h** | 1.0 d | **0.5924** | **0.4677** | interior optimum |
| 6 h | 1.0 d | 0.5990 | 0.4774 | |
| 12 h | 2.0 d | 0.6100 | 0.4928 | |
| 1 day | 4.0 d | 0.6275 | 0.5188 | metronidazole |
| 3 days | 12.0 d | 0.6509 | 0.5326 | metronidazole |
| 1 week | 28.0 d | 0.6541 | 0.5295 | metronidazole |

**There is an interior optimum at about 3 h, and an ideal step is not the best input.** Over a
fixed 24 h window an instantaneous step has finished relaxing in about 3 relaxation times and
the remaining draws sit on the quasi-steady locus, where R112 established there is no `k_e`
information at any density. An input that keeps moving for a few hours keeps the marker off
that locus for more of the window.

**And the slow end does not converge to the record-only floor**, which an earlier draft of
this work asserted. At a one-week onset the 28-day window with 24 draws gives 0.5295 against
the floor's 0.5430. The residual 2.5% is not the suppression's relaxation at all: it is 24
extra samples spanning two dose cycles, which is a small record extension. At that onset the
design has stopped being a relaxation measurement.

---

## 4. What would have to be true, and the actionable number

Best case: a 1% burst assay and sixteen admissions. Frontiers **0.1567** at `theta = 600` and
**0.0633** at 300 (`audit/u9_ke_frontier.py`).

| amplitude | onset | `theta = 600` | `theta = 300` | what it rests on |
|---|---|---|---|---|
| **22.2%** | **6 h** | **0.0806** ✓ | **0.0588** ✓ *(but see below)* | production share, defensible |
| 22.2% | 1 day | 0.1833 | 0.1269 | *neither frontier* |
| 30.7% | 6 h | 0.0569 ✓ | 0.0385 ✓ | that share plus one sd |
| 30.7% | 1 day | 0.1222 ✓ | 0.0858 | |
| 41.0% | 6 h | 0.0538 ✓ | 0.0327 ✓ | metabolite fall, *not a share* |
| 41.0% | 1 day | 0.0808 ✓ | 0.0579 ✓ | |

> **At the defensible amplitude the route reaches the frontier at a fifth of nominal**, if
> the suppression establishes within hours. At metronidazole's actual time constant it misses
> by 2.0x. So the amplitude is not what blocks this route. **The onset is.**

**The cell at a tenth of nominal is NOT established, and the reason is this audit's own
(g4) lesson.** 0.0588 against 0.0633 is a **7% margin**, and section 5.6 puts the spread
between the two independent methods at **9% to 33%** on these cells, with the independent one
giving **0.0781 for this exact cell, which misses**. A number that lands within 7% of a gate
carries no margin and fails under any object but the one it was computed against, which is
exactly what (g4) says about the recommended design's 118 draws. **So: reaches at
`theta = 600` with a 1.9x margin that survives the method spread comfortably; not established
either way at `theta = 300`.**

### The slowest onset that still reaches the frontier

This is the number to hand to whoever chooses the agent, because the amplitude is capped by
biology while the onset is a property of the agent.

| amplitude | `theta` | frontier | **slowest onset that reaches it** |
|---|---|---|---|
| **22.2%** | 600 | 0.1567 | **16.4 h** |
| **22.2%** | **300** | **0.0633** | **6.6 h** |
| 30.7% | 600 | 0.1567 | 38.8 h |
| 30.7% | 300 | 0.0633 | 13.6 h |
| 41.0% | 600 | 0.1567 | 60.5 h |
| 41.0% | 300 | 0.0633 | 28.2 h |

**The requirement and the research item meet here.** At the registered share the agent must
act within about 6.6 h at a tenth of nominal. At one standard deviation above it, 13.6 h. At
the top of Thompson's range, 28.2 h, which a one-day agent satisfies. **So whether the
registered route works is a question about the gut-derived share, which is already on the
research list**, and not a question requiring a new intervention. That is a materially
different answer from "the route needs an agent nobody has".

**These thresholds inherit the same caveat.** They are computed from `u15`'s Fisher, and the
independent method runs 9% to 33% wider, which shortens every threshold. The column at
`theta = 600` has room for that; the column at `theta = 300` does not, and its entries should
be read as upper bounds rather than as figures.

### Which input binds depends on the regime, and two of this audit's own readings were wrong

Each input swept over its registered range at a tenth of nominal, as a ratio of largest to
smallest so the two regimes are comparable.

| regime | amplitude span | onset span | which binds |
|---|---|---|---|
| 5% assay, 1 admission | x1.41 | x1.14 | amplitude, by 1.2x |
| **1% assay, 16 admissions** | x4.34 | **x5.94** | **onset**, by 1.4x |

With a clinic assay and one admission the record dominates the posterior, both inputs move it
little, and the amplitude moves it more. With a 1% assay and sixteen admissions the burst
dominates and the onset binds.

> **`u15`'s docstring asserted that the onset decides the route. A later revision of the same
> file asserted the opposite, that the amplitude does, and computed a span ratio to support
> it. Both are recorded in the file and both are wrong as stated**, because each holds in one
> regime only. The first was right in the regime that decides the verdict and wrong in the
> cheap one; the second was right in the cheap one and wrong in the regime that decides the
> verdict. The lesson is narrow and worth keeping: **a "which input dominates" claim is not
> well posed without naming the regime**, and the span ratio that looks like evidence for one
> is evidence only at the operating point where it was computed.

---

## 5. The independent re-check found a defect, and it reversed the verdict

### 5.1 What the re-check is

`audit/u16_recheck_gut.py` re-derives the same quantity sharing no code path with `u15`:

| | `u15` | `u16` |
|---|---|---|
| residual space | flux | **level** |
| `sd(ln k_e)` from | inverse Fisher | **Monte Carlo sampling spread** |
| sensitivities | integrated alongside the state | **none; the fitter differences** |
| integrator | `scipy` LSODA | **own fixed-step RK4** |
| quasi-steady level | `u14.qss` | **own quadratic solve** |
| estimator | implicit, linearised | **`least_squares`, trf** |

It is built to avoid two errors this audit made before: **R98**, holding a prior-carrying
nuisance at its prior mean, so `b` and `k_e` are drawn from their priors every replicate; and
**R101**, evaluating the clearance at the model's own root, so the record residual compares
the predicted **level** with the **observed** level.

### 5.2 The discrepancy, and the three candidates

The first cells came out 1.06 to 1.23 times `u15`'s Fisher, which curvature could explain. The
strongest arm came out **0.500 times it**, 0.0781 against 0.1561, and **a sampling spread
cannot beat a correct Fisher**. That direction is a defect, not curvature.

Three candidates were tested rather than argued:

1. **The declared sigma convention** (`u16` weights the record with the pointwise slope,
   `u15` with the slope at the mean level). **Refuted twice.** Substituting `u15`'s own
   convention made the gap *worse*, 1.324 against 1.232 on the record-only cell and 0.473
   against 0.501 on the strong one.
2. **A prior mixture**, that an inverse Fisher is conditional at one `k_e` while a Monte
   Carlo draws over the prior. **Refuted by algebra.** For a linear Gaussian model the
   marginal spread of (estimate minus truth) over a truth drawn from the prior equals the
   posterior sd exactly: `s²τ²/(s²+τ²) = σ²`. So drawing from the prior should give a ratio
   of one, and 1.23 is not explained. **A diagnostic run in support of this candidate
   reintroduced R98's own error** and its conclusion is withdrawn; see R120.
3. **Dropping the record entirely**, which isolated it: burst-only gave **0.1304 against
   `u15`'s 0.3538, a ratio of 0.368**, the worst of all. The defect is in the burst.

### 5.3 The defect

Every sensitivity column of `u15`'s burst against a central finite difference of the level
trajectory alone:

| column | `u15` ODE | finite difference | max relative disagreement |
|---|---|---|---|
| `dln m/dln theta` | -0.26676 | -0.26676 | 3.9e-09 |
| `dln m/dln Km` | 0.46997 | 0.46997 | 2.9e-09 |
| **`dln m/db`** | **0.0042062** | **0.0034897** | **2.2e-01** |
| `dln m/dln ke` | -0.56398 | -0.56398 | 3.1e-09 |

**The load-bias forcing was wrong.** `u15` inherited a forcing of exactly `1.0` from
`u14.burst_rows`, where it is **correct** for what `u14` models and says it models: a
prescribed change of *known absolute size*, with the unknown bias added after it, so
`dL/db = 1`. A gut-derived share is not that. It removes a **proportion** of the load actually
being produced, so if true production exceeds the recorded value then the gut share of it is
proportionally larger too:

```
    L(t) = (L_tot + b)(1 - amp f(t))      =>      dL/db = 1 - amp f(t)
```

Three columns agreed to 3e-09 and that one to 2.2e-01. **This is the third time in this audit
that a sensitivity disagreed with its own forward model** (R101 twice, now R119), so the
finite-difference check is now a permanent part of the script as `part_f` rather than a
one-off diagnostic. After the fix every column agrees to better than 3e-08.

### 5.4 What it changed

| | before | after | independent MC |
|---|---|---|---|
| burst-only `sd(ln k_e)`, 22.2%, 6 h, 1%, 16 adm | 0.3538 | **0.1365** | **0.1304** |
| with the record, `theta = 300` | 0.1561 | **0.0588** | *see below* |

The correction makes the burst **2.6x more informative**, because a smaller `b` column means
less of the `k_e` signal can be explained away along the near-collinear `b` direction.

> **The verdict reversed.** Before the correction the route reached neither frontier at the
> defensible amplitude and the conclusion was that the amplitude blocks it. After, it reaches
> both given a fast enough onset, and the conclusion is that the onset blocks it at
> metronidazole's time constant. **An earlier version of this document, and the paragraph
> first written into `OPEN_QUESTIONS.md` section 11, stated the pre-correction verdict.** Both
> are corrected; the record of the change is here and at R119.

**The re-check earned its keep.** Its purpose was to confirm a number and instead it found the
number was wrong, in the direction that mattered, on the one arm that decided the verdict. The
weak arms' ratios of 1.06 to 1.23 would not have raised the question on their own.

### 5.6 What the re-check licenses after the fix, and what it does not

| cell, `theta = 300` | `u15` Fisher | `u16` Monte Carlo | ratio |
|---|---|---|---|
| record only, no burst | 0.5430 | 0.6690 | 1.232 |
| 22.2%, ideal step | 0.4762 | 0.5688 | 1.194 |
| 22.2%, 1 day onset | 0.5188 | 0.5645 | 1.088 |
| 41%, ideal step | 0.3711 | 0.4144 | 1.117 |
| **22.2%, 6 h, 1%, 16 adm** | **0.0588** | **0.0781** | **1.327** |

**Every ratio now exceeds one, which is the only direction available.** A sampling spread can
be wider than a local quadratic and cannot be narrower. The impossible direction is gone, and
it was the thing that located the defect.

**What that licenses:** `u15`'s ranking, its structure, and any claim whose margin exceeds
the 9% to 33% method spread. The `theta = 600` verdict, with a 1.9x margin, is licensed.

**What it does not:** the `theta = 300` headline. 0.0588 against 0.0633 is a 7% margin and
`u16` gives 0.0781 for the same cell, which misses. **That cell is not established either
way.** The honest statement is that the route reaches the frontier at a fifth of nominal and
that a tenth is undecided at this precision.

**Neither file's absolute level is the reference.** They differ in residual space and in the
record's sigma convention, and both are declared. Deciding between them is a larger question
than this round, and is the same question (g2) already asks about which object a verdict is
stated against.

### 5.5 Does the same convention question touch `u14`?

`u14`'s forcing of `1.0` is self-consistent with its stated model, so `u14` is not wrong the
way `u15` was. Its fast arms are withdrawn on other grounds. What survives and could be
affected is section 2.5's safe **raising** step. For a raising step the factor is
`(1 + step) > 1`, which makes the `b` column larger, which should make `k_e` worse, so the
reported conclusion that it "buys almost nothing" should be unchanged or reinforced. That is
tested in `REVIEW_LOG.md` R121 rather than asserted.

---

## 6. The two assay gates are one question in two settings

The home-microsampling CV and the 1% inpatient assay were recorded as separate gates. They
are the same question, **how precisely can this marker be measured**, and one search answers
both. Merged as **(G-ASSAY)** in `OPEN_QUESTIONS.md` section 11.

| setting | the number needed | what turns on it |
|---|---|---|
| **home volumetric microsampling** | any CV at all; none exists | Whether three samples a week finishing Tier 1 in 2.11 years at 10%, or 4.25 at 15%, beats 9.82 in clinic. One a week does not pay at 15% or worse. |
| **an inpatient series**, one line, one run | whether 1% is reachable; asserted, unsourced | The burst route entirely. At a 5% assay and one admission the route is nowhere near any frontier; the 1% assay with sixteen admissions is what puts it in the regime where it reaches. |

The two settings differ in what makes a tight CV plausible, so one source will not answer
both, but one search will find both.

---

## 7. What is not settled

- **(g5) is open** and everything about prescribed load reduction waits on it. Section 3.1
  argues gut-derived suppression is exempt; that argument is not the same as the coupling
  being specified.
- **The suppression's delivered amplitude is treated as exactly known.** A delivered step
  with its own error degrades every figure in sections 3 and 4, and nothing here bounds it.
  R119 is a reminder that how the unknown load bias interacts with a prescribed change is not
  a detail: getting it wrong moved the answer 2.6x.
- **The gut-derived share is a lower bound**, "at least 22.2 +/- 8.5%". Section 4 shows the
  verdict turns on where in that range the truth sits, so the bound is not sufficient.
- **The onset time constant of any specific agent is not sourced here.** Metronidazole acting
  over days comes from `m2_antibiotic_cycle.py`'s docstring, not from a pharmacodynamic study
  of flora suppression, and 6.6 h is a requirement this work states rather than a threshold
  any agent has been measured against.
- **A 1% marker assay remains asserted**, now as (G-ASSAY) setting two.
- **Safety is assumed, not established.** Suppression is routine management, which is not the
  same as safe at the amplitude and repetition sixteen admissions imply.
- **The burst inherits the paper's observation-model assumption**, and inherits it through the
  relaxation rather than the level.
- Nothing here is checked against a primary source; the container has no egress.
