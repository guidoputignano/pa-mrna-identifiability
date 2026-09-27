# Phase 11, round 5. The manoeuvre that has all three properties, the prize that is not where it looked, and a failure mode the paper's own sweep never samples

*Proposals only. Nothing under `paper/` is changed by this file.*

Scripts: `audit/u17_glucose_burst.py`, `audit/u18_rho_glucose.py`, re-check in
`audit/u16_recheck_gut.py`. `REVIEW_LOG.md` R125 to R131.

---

## 0. A process defect first, because it caused five wrong claims

Five claims of mine across this effort were wrong for one reason: **reasoning from a bounded
search as though it were complete.** Twice it was a shallow clone answering a history question
confidently and wrongly (R93). Three more times it was a grep truncated at its default head
limit, where a high-frequency match filled the quota and hid the low-frequency one that
answered the question. `01_biology.tex:151-160` names fasting and catabolism explicitly and
`params.py` registers the gut share with the metronidazole fall attached; `decompensation`
matches crowded both out.

**A truncated search and an exhaustive one look identical in the output**, which is what makes
this class survive re-reading. The rule added to `OPEN_QUESTIONS.md` section 0:

> **Never assert that something is ABSENT from a bounded search.** Before claiming the
> repository lacks a term, a concept or a treatment, run the search **unbounded** and **state
> what was searched**: the patterns, the paths, and the fact that no result limit was applied.

Every search in this round's scripts states its scope and its match count. Doing so
immediately paid: the unbounded read of `params.py` for the partition returned the sentence
**"a protein-restricted diet lowers the first two deliberately"**, which is the additivity
section 4 prices, in the register's own words rather than as an assumption of mine.

---

## 1. The manoeuvre, and why it is the first candidate with all three properties

`01_biology.tex:160` gives the standard emergency regimen as intravenous glucose to suppress
catabolism together with protein withdrawal. Run **electively** rather than in response to
illness, it is a load reduction with three properties, and no earlier candidate had more than
two:

| | fast onset | no catabolic penalty | measured load |
|---|---|---|---|
| a 50% fast (R116, **withdrawn**) | yes | **no**, it provokes the load it removes | no |
| metronidazole (R118) | **no**, bacterial turnover over days | yes | no |
| **elective glucose + withdrawal** | **yes**, hours | **yes**, by design | **yes**, infused |

The third property is the one this round was sent to price, and section 3 finds it worth
almost nothing. The first two are what carry the result.

---

## 2. The amplitude, and the arm that reaches the frontier

`src/common/params.py`, `s_base` row, read in full rather than paraphrased:

> "amino-acid catabolism a mean of **51.7%** (range 24.5 to 66.4%), gut flora at least
> **22.2 +/- 8.5%** by the fall under metronidazole, and roughly 30% unaccounted for ...;
> **a protein-restricted diet lowers the first two deliberately.**"

Glucose with protein withdrawal targets the first share, an osmotic bowel preparation targets
the second. Each anchor is an **upper bound** on its arm, since obligate protein turnover and
residual flora are not zero, so every table sweeps below it.

**Best case, 1% burst assay, sixteen admissions, 6 h onset, infused load.** R124 put the
spread between this Fisher and `u16`'s independent Monte Carlo at 9% to 33%, so a margin under
about **1.33x is not a result**.

| arm | `theta` | `sd(ln k_e)` | frontier | margin |
|---|---|---|---|---|
| bowel prep alone, 22.2% | 600 | 0.0842 | 0.1567 | **1.86x** ✓ |
| bowel prep alone, 22.2% | 300 | 0.0642 | 0.0633 | 0.99x, *undecided* |
| **glucose + withdrawal, 51.7%** | **600** | **0.0548** | 0.1567 | **2.86x** ✓ |
| **glucose + withdrawal, 51.7%** | **300** | **0.0281** | 0.0633 | **2.25x** ✓ |
| both, 73.9% | 600 | 0.0698 | 0.1567 | 2.24x ✓ |
| both, 73.9% | 300 | 0.0367 | 0.0633 | 1.72x ✓ |

> **The glucose arm settles the cell R124 had to withhold.** Round 4's tenth-of-nominal result
> was 0.0588 against 0.0633, a 7% margin inside the method spread, and was reported as
> undecided. At the catabolism share the margin is **2.25x**, which clears the spread with room
> to spare. **That is the round's result, SUBJECT TO the re-check flagged immediately below.**

> ### RESOLVED: the re-check was comparing the wrong statistic, and both cells stand
>
> The independent Monte Carlo initially returned **0.828** at the `theta = 600` headline, the
> direction a sampling spread cannot go against a correct Fisher, which is the signature that
> located the R119 defect. **It was not a defect this time. It was me comparing a robust
> spread with a standard deviation.** Full diagnosis in section 5.5 and at R132.
>
> Re-checked with the standard deviation, which is what a Fisher width is:
>
> | arm | `theta` | Fisher | MC sd | frontier | **margin on the MC sd** |
> |---|---|---|---|---|---|
> | glucose + withdrawal, 51.7% | 600 | 0.0553 | 0.0578 | 0.1567 | **2.71x** ✓ |
> | glucose + withdrawal, 51.7% | 300 | 0.0316 | 0.0387 | 0.0633 | **1.64x** ✓ |
>
> **Both clear the 1.33x method spread on the independent method's own number.** The
> conservative figures to quote are therefore **2.71x and 1.64x**, not the Fisher's 2.86x and
> 2.25x, and the cell R124 withheld is settled at 1.64x rather than 2.25x.

---

## 3. The load-record prize is real, and it is not in the burst

**Sweeping the burst's own load-record term from the diet record's 5% to an exactly known
infused load is worth 1.00x to 1.02x.** Not a rounding artefact: the same sweep at three
amplitudes and two capacities never exceeds 1.02x.

**The diagnosis is structural, not numerical.** `LOAD_FRAC` enters the burst in exactly one
place, the prior on its load bias, **because a burst observation is a marker level with assay
noise and carries no load measurement at all**. It is the *record* whose per-visit sigma is
`hypot(load term, assay term) = hypot(16.14, 14.90)`, and the record is still the ordinary
outpatient diet record whatever happens during an admission. So the premise is right, the 5% is
half of sigma, and it is half of the **record's** sigma, which an infused burst does not touch.

**Where the prize does live**, with the glucose burst held fixed and the record's load term
swept:

| the record's load-record error | `theta = 600` | `theta = 300` |
|---|---|---|
| as published, 5% | 0.0548 | 0.0281 |
| halved, 2.5% | 0.0481 | 0.0254 |
| a fifth, 1% | 0.0416 | 0.0233 |
| a hundredth, 0.05% | 0.0389 | 0.0226 |
| known to 1e-4, and 1e-6 (converged) | **0.0389** | **0.0226** |

**So it is worth 1.41x at a fifth of nominal and 1.24x at a tenth**, which is real and modest,
and it buys that only if **every one of the 118 routine outpatient visits** has its load
measured rather than recalled. The glucose protocol changes the admission, not the outpatient
record, so it does not deliver this. It is a separate and much larger proposal.

*An error made and caught in the course of this: the first version set the load-bias prior
precision to zero at the exact-load limit, which makes `b` a free parameter with no prior,
less information rather than more. It showed as a non-monotone column, the exact-load row
landing worse than the 1% row. Recorded at R128.*

---

## 4. An interior optimum, so the two arms do not add as design

`u17` PART H, fine sweep at a 1% assay, sixteen admissions, 6 h onset:

| amplitude | `theta = 600` | `theta = 300` |
|---|---|---|
| 30% | 0.0546 | 0.0383 |
| 35% | 0.0489 | 0.0307 |
| **40%** | **0.0484** | 0.0274 |
| **45%** | 0.0505 | **0.0267** |
| 50% | 0.0536 | 0.0276 |
| 51.7%, catabolism mean | 0.0548 | 0.0281 |
| 65% | 0.0641 | 0.0331 |
| 73.9%, both arms | 0.0698 | 0.0367 |
| 85% | 0.0764 | 0.0409 |

**The optimum is at 40% at a fifth of nominal and 45% at a tenth**, and the catabolism share
at 51.7% sits just past it, close enough to be near-optimal. **73.9% is clearly worse.**

**Mechanism, and it is a quantity already in this audit.** Suppressing more load drives the
marker further off enzyme saturation, which raises the relaxation rate `J` and lowers `k_e`'s
share of it. That is the same `k_e/J` that R113 used to reject the post-dose placement, acting
here through the amplitude instead of through the increment.

> **So the register's additivity is real as biology and is the wrong thing to want as design.**
> "A protein-restricted diet lowers the first two deliberately" is true, and a design that
> lowers both lands at 73.9%, past the optimum, and measures `k_e` **worse** than glucose alone.
> **Adding a bowel preparation to the glucose arm degrades the measurement.**

> ### WITHDRAWN as stated: the two methods disagree on the SIGN of a small effect
>
> The optimum's depth depends on how the burst's load bias is treated, and the independent
> method does not reproduce it even after the statistic is corrected.
>
> | amplitude | Fisher, b known | Fisher, shared b | **MC sd**, shared b |
> |---|---|---|---|
> | `theta = 300`, 35.0% | 0.0307 | 0.0349 | 0.0549 |
> | `theta = 300`, 51.7% | 0.0281 | **0.0316** | **0.0387** |
> | `theta = 300`, 73.9% | 0.0367 | **0.0324** | **0.0344** |
> | `theta = 600`, 51.7% | 0.0548 | **0.0553** | **0.0578** |
> | `theta = 600`, 73.9% | 0.0698 | **0.0589** | **0.0533** |
>
> With the bias known the 51.7% to 73.9% gap is **31%** at a tenth of nominal, which is the
> table above. With one bias shared it is **2.5% the other way**. And the Monte Carlo, on the
> standard deviation, puts 73.9% **better** at both capacities, by 11% and 8%.
>
> **So "adding a bowel preparation degrades the measurement" is withdrawn.** The two methods
> disagree on the sign of an effect of a few per cent, which is inside their disagreement
> elsewhere, so **the amplitude beyond about 50% is indeterminate at this precision**. The
> sensitivities are not the cause: they agree with finite differences to 3e-09 at every
> amplitude tested, including 85%.
>
> **The decision does not need this resolved.** Glucose with withdrawal alone already clears
> both frontiers at 2.71x and 1.64x, so there is **no precision argument for adding a second
> intervention**, in either direction. A bowel preparation should be decided on safety, burden
> and its own onset, not on `k_e` precision. R133.

An osmotic preparation remains worth an arm **on its own**, where its 22.2% needs the tightest
onset of any candidate, and section 5 gives that number.

---

## 5. The onset requirement, and a convention correction to `u15`

**`u15` let the burst window grow as four onset time constants.** That silently hands a slow
input a free record extension: at a one-week onset the window is 28 days and the 24 draws span
two dose cycles, which is a record extension at a better assay and not a relaxation measurement
at all. **R118 noticed the symptom**, that the slow end sat *under* the record-only floor
instead of converging to it, and diagnosed it correctly without fixing it. An admission is an
admission, so the window is now fixed at one inpatient day.

**The slowest onset that still reaches the frontier**, 1% assay, sixteen admissions:

| arm | `theta` | fixed 1 d window | `u15`'s growing window |
|---|---|---|---|
| bowel prep alone, 22.2% | 600 | **14.9 h** | 15.5 h |
| bowel prep alone, 22.2% | 300 | **5.8 h** | 5.8 h |
| glucose + withdrawal, 51.7% | 600 | **50.0 h** | beyond 13 d |
| glucose + withdrawal, 51.7% | 300 | **30.5 h** | 45.0 h |
| both, 73.9% | 600 | **76.0 h** | beyond 13 d |
| both, 73.9% | 300 | **48.2 h** | beyond 13 d |

The two columns agree where the threshold is near one day and diverge by more than an order of
magnitude where the growing window runs to weeks. **The right column is the artefact.** Round
4's reported 6.6 h and 16.4 h at 22.2% become **5.8 h and 14.9 h**, a small change because
those thresholds sat near a day already, so **round 4's conclusions are unaffected** and only
its long-onset rows were ever contaminated.

**Against the manoeuvres:** glucose acts over hours, so **30.5 h is cleared comfortably**. An
osmotic preparation acts over a gut transit, which is the same order as its 5.8 h requirement
and is therefore the one arm whose onset is genuinely in question. Metronidazole acts over
days and clears none of them, and R118 established that this is a property of bacterial
turnover rather than of the drug, so no change of formulation or route reaches it.

### 5.5 The re-check's own statistic was wrong, and it touches a standing rule

`u16` reported **0.7413 x IQR**, the robust spread, because `OPEN_QUESTIONS.md` section 0 says
*"report the robust spread, not RMSE"*. That rule exists so that collapsed fits do not dominate
a reported spread, and for that purpose it is right. **A Fisher width is a standard deviation,
and 0.7413 x IQR equals a standard deviation only for a Gaussian.** Where the sampling
distribution is heavy tailed the robust spread sits below the sd, and comparing it with a
Fisher manufactures ratios below one, **the direction that is supposed to be impossible and
that correctly located the real defect at R119.**

Both statistics, at six cells, with the skew of the pooled draws:

| amplitude | `theta` | Fisher | robust | **sd** | sd of mid 98% | skew | **sd/Fisher** |
|---|---|---|---|---|---|---|---|
| 22.2% | 300 | 0.0588 | 0.0781 | 0.0869 | 0.0740 | -0.46 | 1.477 |
| 35.0% | 300 | 0.0349 | 0.0476 | 0.0549 | 0.0455 | -0.72 | 1.571 |
| 51.7% | 600 | 0.0553 | 0.0466 | 0.0578 | 0.0452 | +0.06 | **1.046** |
| 51.7% | 300 | 0.0316 | 0.0298 | 0.0387 | 0.0315 | -0.61 | **1.224** |
| 73.9% | 300 | 0.0324 | 0.0241 | 0.0344 | 0.0270 | +0.35 | **1.062** |
| 73.9% | 600 | 0.0589 | 0.0351 | 0.0533 | 0.0418 | +1.01 | 0.904 |

**The robust spread is below the Fisher at four of six cells; the sd is above it at five of
six.** The sixth, at the most extreme amplitude, is 0.904 with a skew of 1.01, where no single
number summarises the distribution well and which is an amplitude section 4 does not recommend
anyway. **The apparent impossible direction was the statistic throughout.**

**A second candidate was tested and refuted, in the direction opposite to the one needed.** An
inverse Fisher is conditional on `k_e = 5` while the Monte Carlo averages over a 0.8-wide
prior, and the R119 algebra showing those coincide assumes the information is the same for
every truth, which here it is not. But the prior-weighted spread comes out **above** the
at-the-mean value, by 1.19x to 2.34x, so this effect predicts a Monte Carlo **wider** than the
Fisher and cannot produce one narrower.

**This reaches backwards, and in the safe direction.** R124's ratios of 1.088 to 1.327 were
robust against Fisher, so they **understate**. Round 4's withheld cell, 0.0588 against a
frontier of 0.0633, re-checks at 0.0869 on the sd rather than 0.0781, a ratio of 1.477 instead
of 1.327. **It misses by more than reported, so round 4's decision to withhold it was right and
is reinforced.**

**Proposed refinement to a section 0 rule, which is the author's to accept or decline.** The
rule as written, *"report the robust spread, not RMSE"*, is sound for reporting and unsound for
comparison against a Fisher. Suggested wording: *report both, use the standard deviation when
comparing a sampling spread against a Fisher width, and use the robust spread to detect
collapse.* Stated as a proposal and **not applied to section 0**, unlike the unbounded-search
rule, which the author instructed.

---

## 6. `rho` is not a property of the chart, and the paper never asks what its own glucose does to it

### 6.1 The object, and two registered numbers that agree without having been made to

`src/design/run.py:1071` generates the record as `s_rec = L_base + rho (s_true - L_base)`, and
`07b_decision.tex:53-62` sweeps `rho` as an **exogenous property of the protocol**. The sweep
is reproduced exactly from its owner: 1.027 against the paper's 1.030, 0.903 against 0.900,
0.086 against 0.086, 0.009 against 0.009.

But the paper's own words say `rho <= 0` is "crisis protocol, where prescribed intake is cut
while **catabolic muscle breakdown** supplies the load". **So `rho` is set by how much
catabolism the protocol allows**, and `01_biology.tex:160` names the intervention that
suppresses catabolism. The paper states both and connects neither.

With `d` the dietary share withdrawal removes, `F` the unsuppressed catabolic load fold, and
`q` the fraction of the catabolic response the glucose suppresses:

```
    rho(q) = -d / ((F - 1)(1 - q) - d),      pole at  q* = 1 - d/(F-1)
```

**A consistency test the repository passes.** `params.py`'s illness row anchors the
conservative catabolic fold at **F = 1.50**, from a 227% to 325% rise in 3-methylhistidine in
sepsis and trauma. The paper's `rho = -0.25` then implies **d = 0.100**, a withdrawal removing
10% of the load. For a patient **already on a protein-restricted diet**, which is the therapy
under study, a dietary share that small is exactly right. **Two registered numbers, chosen
independently in different sections, agree.**

That puts the pole at **q\* = 0.800**.

### 6.2 The verdict neither survives nor reverses. It changes which error it makes

Each `q` mapped to its `rho` and put through the paper's own estimator:

| `q` | `rho` | treat unneeded | fail to treat | median `theta_hat/theta` |
|---|---|---|---|---|
| 0.00 | -0.25 | 100.0% | 0.0% | 0.009 |
| 0.50 | -0.67 | 100.0% | 0.0% | 0.009 |
| 0.70 | -2.00 | 100.0% | 0.0% | 0.010 |
| 0.78 | -10.00 | 100.0% | 0.0% | 0.010 |
| | | | | *pole at q\* = 0.80* |
| 0.82 | +10.00 | 0.0% | **100.0%** | 1.531 |
| 0.90 | +2.00 | 0.0% | **100.0%** | 0.776 |
| 0.95 | +1.33 | 0.0% | **93.8%** | 0.938 |
| 1.00 | +1.00 | 17.5% | 41.8% | 1.027 |

> **Partial glucose protection is a failure mode the paper's sweep never samples.** That sweep
> stops at `rho = 1` and never enters `rho > 1`. Below the pole the estimator treats everyone,
> which is the paper's published catastrophe. Above it, the estimator **fails to treat
> everyone**, which is the dangerous direction, and it does so totally rather than at the
> benign cell's 41.8%.

**Seed spread on the three cells that carry the claim**, eight seeds at 400 replicates, because
a reversal on one seed is not a reversal:

| `q` | `rho` | treat unneeded | fail to treat | median `theta_hat/theta` |
|---|---|---|---|---|
| 0.70 | -2.00 | 100.0% [100.0, 100.0] | 0.0% [0.0, 0.0] | 0.010 [0.009, 0.011] |
| 0.90 | +2.00 | 0.0% [0.0, 0.0] | 100.0% [100.0, 100.0] | 0.812 [0.686, 0.933] |
| 1.00 | +1.00 | 18.8% [16.8, 21.5] | 40.7% [37.5, 44.0] | 0.926 [0.484, 1.079] |

**The reversal is clean**: zero spread on both sides, 100% to 0% and 0% to 100% across every
seed. The finding does not rest on a seed.

**And the spread caught something about the paper's own number.** `07b_decision` reports the
benign cell's median `theta_hat/theta` as **1.03**, which is its single seed 11. The eight-seed
range is **0.484 to 1.079** with a mean of 0.926. The published figure sits near the top of its
own range, and the bottom of that range is a two-fold underestimate of capacity rather than the
3% error the text implies. This is not a wrong number, it is a single-seed number quoted without
its spread, in the row the section leads with. R131.

### 6.3 What that means for the protocol, and it is actionable

`rho > 1` means the **record overstates** the load change: the chart records the full
prescribed cut while the true load moved less, because the glucose absorbed most of the
catabolic response. **Under an infused, measured load that cannot happen**, because the record
*is* the infusion and `rho ≈ 1` by construction.

> **So the hazard is specific and avoidable: running glucose protection while still recording
> intake from the diet chart.** The same manoeuvre with the load taken from the pump is benign.
> This is the one lever in the repository that moves `rho` at all. Every other entry in
> `07b_decision`'s sweep is a fact about the chart; this one is a fact about the protocol.

### 6.4 The limit of the paper's parameterisation, measured rather than asserted

The generator holds the marker levels **fixed** and scales only the record, so section 6.2
prices the record-corruption channel and **not** the excitation loss. A physiological `rho`
also moves the true load, which moves the marker:

| `q` | `rho` | true load deviation | excitation retained |
|---|---|---|---|
| 0.00 | -0.25 | 0.4000 | 1.00x |
| 0.50 | -0.67 | 0.1500 | 0.37x |
| 0.78 | -10.00 | 0.0100 | 0.02x |
| 0.82 | +10.00 | -0.0100 | -0.02x |
| 1.00 | +1.00 | -0.1000 | -0.25x |

**The two channels pull in opposite directions**: as `q` rises the record becomes more
faithful while the excitation shrinks, and at the pole the excitation is zero while `rho` is
infinite. So section 6.2's numbers near and beyond the pole are **optimistic** and are not a
verdict there. The sign of the effect, a fail-to-treat regime at `rho > 1`, does not depend on
the pole and stands.

**And the deviation changing sign past `q*` is the connection between the two halves of this
round.** Beyond the pole the true load genuinely *falls*, which is precisely the elective load
reduction section 2 prices as a `k_e` burst. **One manoeuvre, two estimands, opposite
requirements on the record**: the burst wants the largest true load change it can safely get,
and the decision wants the record to track it.

---

## 7. What is not settled

- **(g5) is still open**, and this round does not close it. Sections 1 and 6 argue the glucose
  arm is exempt from the catabolic objection *because suppressing catabolism is what the
  glucose is for*, which is an argument and not a specified coupling.
- **`q`, the fraction of the catabolic response glucose actually suppresses, is nowhere in this
  repository.** Section 6 sweeps it because nothing here can supply it, and section 6.2 shows
  the verdict turns on which side of 0.8 it falls. **This is the round's new gate.**
- **`F = 1.50` is a conservative illness fold, not a fast-specific one**, and the register's
  own row says the 3-methylhistidine bound "is looser than it looks".
- **The catabolism share is an upper bound** on what withdrawal plus glucose removes. Section 4
  shows the optimum is at 40 to 45%, so the arm is near-optimal at its anchor and this matters
  less than it would have.
- **(G-ASSAY) still gates everything**, and the 1% inpatient assay is what puts the arm in the
  regime where it reaches at all.
- **Safety of an elective admission** with glucose and total protein withdrawal, repeated
  sixteen times, is assumed and not established. It is routine as *emergency* management; doing
  it electively and repeatedly to a stable patient is a different proposition.
- Nothing here is checked against a primary source; the container has no egress.
