# Phase 11, round 2. The `k_e` requirement frontier, the carnitine lever, and what could measure it

*Proposals only. Nothing under `paper/` is changed by this file.*

Scripts: `audit/u9_ke_frontier.py`, `audit/u10_carnitine_modulation.py`.

---

## 0. Why this replaces the `K_m` bypass

Round 1 found that closing the fibre buys about five per cent and that **`k_e` is the
load-bearing external quantity**. That turns the scope limit from a verdict into a
requirement, and the requirement is the subject of this file.

**Footing** as `audit/PHASE11_KM_BYPASS.md` section 0 throughout. The Tier 1 arm is 118
monthly draws with the increment's scale at the published 14.02%, on
`audit/t2_tier1_margin.py`'s own record convention, which this work reproduces exactly:
**118 draws, `sd(ln theta) = 0.2000`** (`audit/u11_trade_curve.py` anchor).

**Independent confirmation of round 1, on the author's footing.** At 118 visits the author
reports an exact `K_m` worth 1.007 / 1.010 / 1.011 and an exact `k_e` worth 3.96 / 9.36 /
8.70 at `theta = 3000 / 600 / 300`, against this work's 1.056 / 1.037 / 1.058 and 2.09 /
6.43 / 7.43 at 25 visits. The two differ by footing and agree in direction and in which
quantity binds.

---

## 1. The frontier

**The `k_e` log-prior width at which the record lands exactly on the 0.20 gate.** Solved by
bisection on the width, medians over 24 intake records.

| capacity | % of nominal | Tier 1 record (118 draws) | recommended record (25 draws) |
|---|---|---|---|
| 3000 | 100% | **0.8001**, i.e. the existing prior exactly | 0.9815 |
| 600 | 20% | **0.1567** (15.7%) | 0.1508 |
| 300 | 10% | **0.0633** (6.3%) | unreachable at any width |
| 150 | 5% | unreachable at any width | unreachable at any width |

The author reports 0.8 suffices at nominal, 0.177 at `theta = 600` and 0.0665 at
`theta = 300`. Reproduced.

**The nominal row is 0.8 to four figures and that is not a coincidence.** 118 is *defined*
as the draw count landing on 0.20 at a 0.8 `k_e` prior, so the frontier at nominal capacity
must return the prior it was constructed from. It is an internal consistency check that the
bisection passes, not a finding.

**The curve behind it**, Tier 1 record, medians over 24 records:

| `sd(ln k_e)` | theta=3000 | theta=600 | theta=300 |
|---|---|---|---|
| 0.80 | 0.2000 | 0.7048 | 1.2062 |
| 0.40 | 0.1521 | 0.4165 | 0.7890 |
| 0.20 | 0.1375 | 0.2377 | 0.4485 |
| 0.10 | 0.1336 | 0.1553 | 0.2614 |
| 0.05 | 0.1326 | 0.1258 | 0.1812 |
| 0.02 | 0.1323 | 0.1160 | 0.1510 |
| 1e-4 | 0.1323 | 0.1141 | 0.1445 |

Each column floors: the increment's prior does not shrink with `k_e`, so even an exact
`k_e` leaves 0.132, 0.114 and 0.145. **`theta = 150` is unreachable because its floor sits
above the gate**, which is why that row has no frontier rather than a large one.

> **The scope limit, restated as a requirement.** Not *"capacity is unrecoverable in the
> severe patient"*, which is what `11_scope` says now, but *"capacity is recoverable at a
> tenth of nominal if and only if the non-PCC clearance is known to about six per cent, and
> is unrecoverable at a twentieth at any precision"*. The first is a verdict a reader can
> only accept; the second is a specification a reader can cost.

---

## 2. Prescribed carnitine modulation: a real gain, and an order of magnitude too small

### 2.1 Why the paper has not already closed it

`eq:input_carn` writes the carnitine channel in exactly the form the design wants,
`k_e(t) = k_e + kappa u_c(t)`, and `sec:selective` calls it one of the two **selective**
perturbations. Section 1 has just made `k_e` the binding quantity. `sec:whoto` does close
carnitine, but on two arguments neither of which reaches a **prescribed** swing: that
carnitine cannot touch the `theta` sensitivity's first factor, and that the dose is
adjusted in response to the patient's state and so is endogenous. A dose a protocol
prescribes is exogenous by construction, exactly as the dose train is.

### 2.2 The parameterisation, and the one that was wrong

Writing `f` for the carnitine-mediated share of the alternative route,
`k_e(t) = k_e (1 + f·amp·s(t))` with `s` the square wave of
`audit/t3_diet_modulation.py` at its 21 d period.

**The first attempt carried the two unknowns as the non-carnitine baseline and the scale,
each at the paper's 0.8 prior. That is wrong and it is wrong silently**: at a share of 0.8
it puts a 0.8 prior on a term carrying a fifth of the clearance and another on the term
carrying four fifths, which constrains the *total* better than the paper's single 0.8 prior
on it. It reported gains up to **2.97x**, and the gain was the reparameterisation.

The honest coordinates are the **total** and the **share**, with `ln K` carrying the paper's
own 0.8 prior because it is a prior on the same object. The check this buys is decisive and
the first version failed it: **at zero amplitude the arm must reproduce the paper's design
exactly.** It now does, over twelve cells, to `0.00e+00` — 0.169646, 0.809698, 1.580037,
3.058346 at `theta = 3000 / 600 / 300 / 150`.

### 2.3 What it buys

Gain is `sd(ln theta)` at amplitude zero over `sd` at that amplitude.

| reading of "kappa" | gain range |
|---|---|
| carnitine share only | 1.010 to 1.338 |
| plus the mRNA dose scale at 0.8 | **1.001 to 1.309** |
| plus the mRNA dose scale at 14.02% | 1.006 to 1.334 |

**The author's 1.00 to 1.04 reproduces, and it is a nominal-capacity result.** With the
mRNA dose scale carried at 0.8, `theta = 3000` gives 1.007 at a share of 0.2 and 1.018 at
0.8, because that prior then dominates everything the swing could buy. In the disease range
the same reading gives 1.21 to 1.31, because there the binding term is the alternative
clearance and not the increment. Both numbers are right and they are about different
capacities. *"kappa" is ambiguous in this repository's own vocabulary* —
`t2_tier1_margin.py` uses it for the mRNA dose scale — and the ambiguity is worth exactly
this factor of ten in the gain.

### 2.4 The gap, which is the whole answer

Best case throughout: share 0.8, amplitude 50%, the good period.

| theta | no swing | best swing | gain | factor still short of 0.20 |
|---|---|---|---|---|
| 3000 | 0.1696 | 0.1358 | 1.249 | meets it |
| 600 | 0.8097 | 0.6052 | 1.338 | **3.0x** |
| 300 | 1.5800 | 1.2052 | 1.311 | **6.0x** |
| 150 | 3.0583 | 2.3574 | 1.297 | **11.8x** |

**So the negative result is not "the swing buys nothing".** It buys a real factor of up to
1.34. It is that the gap is a factor of six, and a third is not a sixth. Stating it the
first way invites a reader to stack the swing with something else; stating it this way says
how much that something else would have to be worth.

### 2.5 The proposed mechanism is refuted, and the correction is more useful

The mechanism offered was that modulating an input known up to scale **pins the scale and
not the baseline**, so `theta` gains nothing. The columns say otherwise. At `theta = 300`,
`f = 0.5`, as the amplitude goes from 0 to 50%:

| amplitude | `sd(ln theta)` | `sd(ln K)`, the total | `sd(ln f)`, the share |
|---|---|---|---|
| 0.00 | 1.5800 | 0.7228 | 0.8000 |
| 0.15 | 1.3420 | 0.6088 | 0.6117 |
| 0.30 | 1.2586 | 0.5713 | 0.5464 |
| 0.50 | 1.2253 | 0.5594 | 0.5241 |

**Both improve, by almost exactly the same amount.** Neither is preferentially pinned. The
real mechanism is quantitative rather than structural: the swing does buy information about
the total alternative clearance, by about a third, and `sd(ln theta)` in the disease range
is governed by `sd(ln K)`, so it improves by about a third too. Nothing is being cancelled.
The gain is real and it is small.

### 2.6 Judged on the frontier, which is the fair test

Carnitine modulation is a way of measuring `k_e` *in situ* with no new instrument, so the
right question is how close it brings `sd(ln K)` to section 1's frontier.

| theta | frontier | share | amp | clearance swing | `sd(ln K)` | factor short |
|---|---|---|---|---|---|---|
| 600 | 0.1567 | 0.8 | 0.50 | 40% | 0.5670 | 3.6x |
| 600 | 0.1567 | 0.8 | 1.00 | 80% | 0.5428 | 3.5x |
| 300 | 0.0633 | 0.8 | 0.50 | 40% | 0.5452 | 8.6x |
| 300 | 0.0633 | 0.8 | 1.00 | 80% | 0.4687 | **7.4x** |

An amplitude of 1.00 is a carnitine dose swinging between zero and twice the prescribed
mean, the physical ceiling on this lever and almost certainly not admissible. Even there the
total clearance is short of the frontier by 7.4x.

> **Record this before anyone tries it.** Prescribed carnitine modulation is a real but weak
> `k_e` measurement. The requirement it has to meet is an order of magnitude tighter than the
> largest physically possible swing delivers.

**Clinical safety is assumed, not established**, exactly as `t3_diet_modulation.py` assumes
it for diet. Whether a prescribed carnitine swing of any size is admissible is not a
question these scripts can answer, and the answer does not matter if the swing cannot close
the gap.

**The aliasing trap applies here as it does to diet.** A modulation whose period is the dose
interval lands every draw at the same phase of it.

---

## 3. What could measure `k_e`, and whether a timed collection is it

**The modelling question was put first and it decides the costing.** *Is the model's `k_e`
the same object as the urinary clearance a timed collection measures?* The answer is **no**,
on three separate counts, and they do not overlap.

### (i) Units

`k_e = CL / V_d`, so the conversion needs a **volume of distribution for 3-HP**, which is
not measured either. This repository already performs exactly this conversion for its own
analogue: the `k_e` row of `src/common/params.py` takes methylmalonic acid at 103 mL/min and
divides by a `V_d` of 0.48 L/kg to reach a rate. **The `V_d` does as much work as the
clearance**, and it is the term nobody has for this marker.

### (ii) Channel coverage

`src/common/params.py` registers `k_e` with `kind=LUMPED`. That file's own taxonomy says a
lumped constant *"absorbs several distinct mechanisms into one number"* and therefore
**"cannot be checked against a direct measurement even in principle"**. A urinary collection
sees the renal channel; any non-renal disposal of 3-HP sits inside `k_e` and outside the
measurement. *The repository had already answered this question about its own constant, in
the register, and the answer was not carried into the design discussion.*

### (iii) Species

The parent model's alternative channel is carnitine conjugation, which esterifies
propionyl-CoA to **propionylcarnitine**. That is a different molecule from 3-HP. The channel
the parent says dominates the alternative route does not remove the marker at all, and a
3-HP collection cannot see it even in principle.

### What that leaves, quantitatively

`CL_r / V_d` is a **lower bound** on `k_e`, not an estimate of it. Writing
`k_e = CL_r / (V_d · f)` with `f` the renal share of the non-PCC route, the three
uncertainties add in quadrature, so **each term alone** must sit below the frontier:

| theta | frontier | each term alone under | all three equal, each under |
|---|---|---|---|
| 600 | 0.1567 | 15.7% | 9.05% |
| 300 | 0.0633 | **6.3%** | **3.66%** |

**A `V_d` for 3-HP known to a few per cent, and a renal share known to a few per cent, are
both further from the published record than the clearance itself.** So the route does not
reach the disease-range frontier, and the reason is not the measurement's own precision.

### What the route is still worth, and it is not nothing

- It bounds `k_e` **from below**, where the methylmalonic-acid analogue in `params.py` bounds
  it from above. `sec:scopeprec` shows the working band moving across a factor of eight in
  the **value** of `k_e`, which is a different question from its precision and one a bound
  does address.
- The renal share `f` is itself measurable, as the fraction of an administered tracer
  recovered in urine, from the same collection.
- Nothing in it needs new technology, which is the property that made it worth asking about.

### What goes on the gate list

**Not "the achievable precision on a timed collection", which is only the first of three
terms.** What goes on the list is the pair:

1. **a volume of distribution for 3-HP**, and
2. **the renal share of its non-PCC disposal**.

Neither was located here and this container cannot check. They join the home-microsampling
assay CV from round 1.

---

## 4. What is not settled

- **Nothing here is checked against a primary source**, `OPEN_QUESTIONS.md` section 1 still
  holding.
- **The carnitine-mediated share `f` is swept, not known.** Every figure in section 2 is
  stated across `f` from 0.2 to 0.8 for that reason.
- **The dose-to-`k_e` relationship is taken as linear**, which is `eq:input_carn`'s own form.
  If it saturates, every gain in section 2 is an overestimate and the negative result gets
  stronger.
- **Safety of any prescribed carnitine swing is assumed**, as above.
