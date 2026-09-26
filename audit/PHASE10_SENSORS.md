# Phase 10. Which observable identifies theta

Marks, as `PHASE4_SUMMARY.md` uses them: **[proven]** mathematics, exact or symbolic.
**[measured]** computed, footing stated. **[judgement]** a call for the author.

**Nothing under `paper/` is edited and no constant in `src/common/params.py` has moved.**

Section 7 of the paper says where the disease-range limit has to be attacked: "on sigma, or on an
observable that still responds to theta where the marker has stopped". Phase 9 spent itself on sigma.
This phase is the second clause.

**Every NUMBER taken from the 25 September reports is an unverified bound**, named as such, pending
the verification session; their MECHANISMS and functional forms are used as structure, which is what
they are good for. `refs/review-20260925/README.md` records why.

---

## Step 1. The measurement equations, before any precision is computed

This step is a gate. Each candidate is written down in this model's parameters and units, or dropped
with a reason. **A dropped row with a stated reason is worth more than a modelled row resting on an
invented equation**, and three rows are dropped below.

### 1.0 The model, and the two things a sensor can attach to

The balance is

    (theta + Delta(t)) * g(m) + k_e * m = L + b,        g(m) = m / (K_m + m)

with `m` the marker, `L` the recorded substrate load, `b` its bias, and `Delta(t)` the therapy's
increment. Write the **enzymatic flux** `F = (theta + Delta) g(m)` and the **non-enzymatic flux**
`k_e m`. At a draw, `F + k_e m = L + b` identically.

**The whole of Phase 10 turns on one line of algebra.** The record's sensitivity to `ln theta` is

    d/d ln theta [ theta g(m) ]  where  theta g(m) = L + b - k_e m

so **theta is recovered from the marker record as a difference of two large numbers**, one of which,
`k_e m`, is known only through a prior. Concretely, inverting the balance for theta from a marker
reading requires

    theta = (L + b - k_e m) / g(m)            the marker route, needs k_e
    theta = S * (L + b) / g(m)                the flux route, needs S and NOT k_e

where `S = F / (L + b)` is the enzymatic share of disposal. **[proven]**

**What that costs in the disease range, arithmetically.** Taking `k_e` 50% higher than the truth:

| capacity | theta from the marker route with k_e 1.5x | truth |
|---|---|---|
| 100% of nominal | 2725.0 | 3000.0 |
| 5% of nominal | **-235.1** | 150.0 |

At normal capacity a 50% error in `k_e` is a 9% error in theta. **At 5% of nominal the same error
drives the estimate negative**, which is the difference-of-large-numbers made concrete and is why an
observable that never forms that difference is worth looking for. **[measured]**

**The flux route's identity is exact**, `theta = S (L+b) / g(m)`, reproducing theta to machine
precision at every capacity tested, with `k_e` appearing nowhere. **[proven]** It trades a dependence
on `k_e` for a dependence on `L`, which is candidate 9's business.

**And the sensitivity orders the other way round.** `S = 1 - k_e m / (L+b)`, so
`dS/d ln theta = -(k_e/(L+b)) dm/d ln theta`:

| capacity | S | dS/d ln theta | theta*g(m), the marker's sensitivity |
|---|---|---|---|
| 100% | 0.8451 | 0.1418 | 272.7 |
| 50% | 0.7173 | 0.2280 | 231.5 |
| 20% | 0.4723 | **0.2832** | 152.4 |
| 10% | 0.2917 | 0.2274 | 94.1 |
| 5% | 0.1630 | 0.1447 | 52.6 |
| 2% | 0.0697 | 0.0666 | 22.5 |
| 1% | 0.0357 | 0.0349 | 11.5 |

**The share's sensitivity peaks at about 20% of nominal and is non-monotone; the marker's falls
monotonically by a factor of 24 over the same range.** These reproduce the author's motivating figures
exactly (0.142, 0.283, 0.227, 0.145, 0.035 at 100, 20, 10, 5 and 1%). **[measured]**

**This is a steady-state share and not yet a measurement.** Whether any instrument reads it is what
the rest of step 1 decides, and it is where two of the three drops happen.

---

### 1.1 Candidate 1: the 1-13C-propionate breath test **[judgement]**, and it does not derive cleanly

**What it is a function of.** A dose of 1-13C-propionate enters the propionyl-CoA pool. Carboxylation
by PCC sends the labelled carbon to methylmalonyl-CoA, succinyl-CoA and the TCA cycle, where it
appears as 13CO2 and is exhaled. The cumulative percentage of the dose recovered as breath 13CO2 over
a window W is the reading.

**The honest measurement equation has three factors and the model supplies only the first:**

    PDR_cum(W)  =  S_tracer  x  R(W)  x  (1 - ret)

- `S_tracer`, the fraction of the tracer committed to the carboxylation route rather than to any
  other disposal. This is the quantity the model can speak to.
- `R(W)`, the fraction of carbon committed to oxidation that has actually appeared in breath by W.
  This is TCA-cycle and bicarbonate-pool kinetics, not in this model.
- `(1 - ret)`, bicarbonate retention, the label fixed in body carbon pools and never exhaled.

**What has to be assumed to write it down at all, and this is the gate.**

1. **That `S_tracer` equals the model's steady-state `S`.** It does not, in general. `S` is a ratio of
   fluxes at a steady state; `S_tracer` is the fate of a bolus perturbing that steady state, read over
   a finite window. They coincide only if the tracer equilibrates with the endogenous pool before
   being disposed of, and if the bolus is small enough not to move `g(m)`. The second is what
   "micro-dose" buys and is checkable; **the first requires a propionyl-CoA pool turnover time that
   this model does not have, because the model has no propionyl-CoA state at all.**
2. **That the reading is whole-body.** Breath 13CO2 integrates PCC flux in every tissue. The model's
   `theta` is a whole-body maximal velocity (`n4_theta_scope`), so that matches. But the therapy's
   `Delta` is hepatic, and R5's whole point is that the extrahepatic fraction is not negligible. **So
   the breath test reads `theta_total` while `Delta` moves `theta_hepatic`**, and the two enter the
   same equation only if their ratio is known.
3. **That `R(W)` and `ret` are known or cancel.** They do not cancel in a single absolute reading. The
   conventional fix is an assumed bicarbonate recovery factor, which is a borrowed constant of exactly
   the kind `refs/review-20260925/README.md` says not to import unchecked.

**And the fatal one, which is structural rather than a missing constant.** The model's enzymatic term
is `(theta + Delta) g(m)`, a Michaelis-Menten clearance **of the marker**. Under the end-product
reading that Phase 9 3b made primary, the marker 3-HP is NOT the enzyme's substrate; propionyl-CoA is.
**So in this model there is no expression for the flux of propionyl-CoA through PCC**, which is
exactly what the breath test measures. The quantity `S = 1 - k_e m/(L+b)` is the share of the
MARKER's disposal that runs through the enzymatic term, which under the primary reading is not a
physical flux at all.

**Verdict: the breath test cannot be written down honestly in this model's parameters.** Not for want
of a constant, but because the model lacks the state the instrument reads. **[judgement]**

**What that does NOT mean.** It does not mean the breath test is uninformative; it means this model
cannot say how informative it is. Two routes would make it expressible, and both are new models
rather than new sensors:

- **(a)** Adopt the model-literal reading, in which the marker is the substrate and `S` is a real
  flux share. Then everything above derives, and the price is the reading Phase 9 3b.4 priced at a
  factor of 5.8 on `sd(ln theta)` and declined to adopt.
- **(b)** Extend the model with an explicit propionyl-CoA state, `P`, with the enzyme acting on `P`
  and the marker formed from `P` by overflow. Then `S` is well defined, the breath test attaches to
  it, and R3's functional forms for the other markers attach to it too.

**Route (b) is what this phase should recommend and cannot itself run.** It is a structural change to
the model, it is the thing that would make five of the ten candidates expressible at once, and it
should be a decision before Phase 10 steps 2 to 4 are attempted. **[judgement]**

---

### 1.2 Candidate 2: ex vivo PCC activity **[measured]**, expressible, already priced

**Measurement equation.** `A_exvivo = kappa * theta * exp(eps)`, with `kappa` the cell-to-whole-body
mapping and `eps ~ N(0, tau_map)` its scatter. It touches `theta` and nothing else in the balance.

**Nuisances of its own:** `kappa` and `tau_map`, neither of which any human study reports (R4).

**Status: already priced, as a prior, in S5**, and Phase 9 3b.2b named the resulting pattern. Carried
here as a sensor rather than re-derived. Its single-parameter attachment is why S5's quadrature
identity held. **[measured]**

---

### 1.3 Candidates 3 to 7: the marker panel

**Candidate 3, plasma 3-HP.** The current design and the baseline every other row is measured
against. Measurement equation: the balance itself, `m` observed with the noise Phase 9 priced.
Expressible by construction. **[proven]**

**Candidates 4 to 7 share one obstruction and it is the same one as the breath test.** R3 gives their
mechanisms, and every one of them is a function of the propionyl-CoA pool, which this model does not
carry:

| candidate | R3's functional form | new unmeasured state | is that state observable |
|---|---|---|---|
| 4, 2-MC | bilinear: propionyl-CoA x oxaloacetate, via citrate synthase | oxaloacetate | no, and it carries the whole TCA state |
| 5, C3 | near-equilibrium, scaled by free carnitine | free carnitine | **yes**, C0 is on the same assay |
| 6, C3/C2 | near-equilibrium ratio, propionyl-CoA against acetyl-CoA | acetyl-CoA | no, and it moves with feeding |
| 7, propionylglycine | sigmoidal through GLYAT, glycine-dependent | glycine | **yes**, on an amino-acid panel |

1. **All four are functions of propionyl-CoA, and the model's only state is the marker.** To write any
   of them down, the model needs the `P` state of route (b). **So candidates 4 to 7 are dropped at
   step 1, with the same reason as candidate 1**, and they are recoverable together by the same
   structural change. **[judgement]**
2. **Two of the four new states are themselves observable**, free carnitine and glycine, which matters
   because it decides whether the new nuisance is free or measured. C3 and propionylglycine are
   therefore the more promising pair, and 2-MC and C3/C2 the less, which **inverts R3's own
   recommendation** to add C3/C2 as the second marker. R3 recommended it for having a distinct
   functional form; on identifiability grounds the distinct form is a liability unless the state it
   depends on is measured, and acetyl-CoA is not. **[judgement]**
3. **The paper's panel remark is about noise averaging and is not what any of this is about.** Phase 9
   item 7 already found the shared-versus-specific split unestimable. **Averaging and structural
   information are different benefits; only the first was ever priced, and the second cannot be priced
   until the model carries `P`.** **[judgement]**

---

### 1.4 Candidates 8 to 10: the nuisance sensors

**Candidate 8, IV 13C3-3-HP micro-dose.** Measurement equation: `CL_total = k_e + theta/(K_m + m_0)`
under the model-literal reading, or `CL_total = k_e` under the end-product reading. **Expressible,
already derived and priced in S2b**, including the factor of 5.8 that separates the two readings.
Carried here unchanged. **[measured]**

**Candidate 9, whole-body propionate turnover by isotope dilution.** Measurement equation:
`R_a = L + b`, the total appearance rate of propionate, which is the balance's right-hand side.
**Expressible, and it is the sensor the flux route needs**, since `theta = S (L+b)/g(m)` requires
`L + b`. Nuisance of its own: the tracer's own steady state and the usual dilution assumptions.
`audit/t9_measured_load.py` already treats a measured load and is the place to attach it. **[proven]**

**Candidate 10, timed urine with paired plasma.** Measurement equation:
`CL_renal = U_rate / m = f_renal * k_e`, with `f_renal` the renal fraction of non-enzymatic disposal.
**Expressible but it measures a product**, and R1 records that `f_renal` is NOT FOUND in any
population. So it bounds `k_e` from below and does not measure it. Already stated in S2b.3 and carried
unchanged. **[measured]**

---

### 1.5 What step 1 returns

| # | candidate | expressible in this model | why not |
|---|---|---|---|
| 1 | 1-13C-propionate breath test | **no** | the model has no propionyl-CoA state, and under the primary reading `S` is not a physical flux |
| 2 | ex vivo PCC activity | yes | attaches to `theta` alone |
| 3 | plasma 3-HP | yes | the balance itself |
| 4 | 2-MC | **no** | needs propionyl-CoA and oxaloacetate |
| 5 | C3 | **no** | needs propionyl-CoA; free carnitine is observable |
| 6 | C3/C2 | **no** | needs propionyl-CoA and acetyl-CoA |
| 7 | propionylglycine | **no** | needs propionyl-CoA; glycine is observable |
| 8 | IV 13C3-3-HP | yes | derived in S2b |
| 9 | propionate turnover | yes | it is `L + b` |
| 10 | timed urine | yes, as a bound | measures `f_renal * k_e` |

**Five of ten candidates are dropped, and all five for the same reason.** The model's only state is
the marker; every observable that reads the propionyl-CoA pool, which is every flux observable and
every alternative marker, has no expression in it.

**So the interesting half of this phase does not exist in the current model, and the author asked to
be told that in one report rather than after four.** **[judgement]**

**What does exist, and is worth running as steps 2 to 4 without any structural change:** candidates 2,
3, 8, 9 and 10, which is the marker record plus the three nuisance sensors plus the ex vivo anchor.
That set can answer "marker alone against marker with a measured `k_e` against marker with a measured
load", and the load row is new: **no phase has yet priced measuring `L` instead of `k_e`**, and the
algebra above says it is the substitution the flux route would have made.

**What would change the answer** is route (b), a propionyl-CoA state. It makes candidates 1, 4, 5, 6
and 7 expressible at once, it is the only way to price structural information rather than noise
averaging, and it is a model change and therefore the author's to authorise.

**And the hierarchical family is untouched by all of this.** `k_e` and `K_m` are population parameters
and `theta` is per-patient, so a cohort fit estimates the shared ones without measuring them in
anyone. It needs no new observable and no new state, and it is the first thing a pharmacometrics
reviewer asks. It is worth running whatever is decided about route (b).

---

## The reduction verdict

`audit/q13_two_state_reduction.py`, 5 of 5 checks passing. Authorised as a reduction check rather
than a rebuild, and scoped to that: **no sensor is priced below.**

### R.0 The two-state model

    dp/dt = L + b - (theta + Delta) p/(K_p + p) - v_o(p)
    dm/dt = alpha * v_o(p) - k_e m

`p` is intracellular propionyl-CoA, `m` the plasma marker, `v_o(p)` the overflow routes that drain
`p`, and `alpha` the fraction of overflow flux appearing as the **measured** marker. **The enzyme now
acts on its own substrate**, which is the point: the paper's discussion names the `p`-to-`m` map as
its deepest assumption, and carrying `p` explicitly removes that assumption rather than adding one.

**New parameters:** `K_p`, `c` (the overflow rate constant, where `v_o = c p`), `alpha`. **What is now
assumed that was not before:** that both states are quasi-steady on the sampling timescale, and that
the overflow is linear in `p`. The first is a timescale claim and checkable; the second is what R3
describes for 3-HP, "a near-linear overflow of propionyl-CoA", and is exactly the condition R.3 tests.

### R.1 The reduction holds, exactly **[proven]**

`dm/dt = 0` gives `p = k_e m / (alpha c)`. Substituting into `dp/dt = 0`:

    (theta + Delta) * m/(K_m + m) + (k_e/alpha) * m = L + b

**which is the paper's balance exactly**, sympy residual **0**, provided

    K_m       = K_p * alpha * c / k_e   =   K_p * c / k_e_paper
    k_e_paper = k_e / alpha

Numerically the two models agree on `m` to **6.7e-15** across 54 combinations of capacity (100% down
to 2%), load (0.5x to 2x) and increment (0, 400, 6000), and under four independent `(alpha, c)`
parametrisations holding the same reduced constants. **[proven]**, **[measured]**

**So every existing result in this project survives as the reduction of a better-founded model.** That
is a strictly better position than the paper is in today, and it is the outcome the author's scoping
hoped for.

### R.2 But neither of the paper's two nuisance constants means what its name says **[proven]**

This is the finding, and it was not visible before the state was made explicit.

| the paper's constant | what it actually is |
|---|---|
| `theta` | **unchanged.** The enzyme's maximal velocity on its own substrate. |
| `K_m` | **not** the enzyme's Michaelis constant. It is `K_p * c / k_e_paper`, that constant scaled by the overflow-to-clearance ratio. |
| `k_e` | **not** the marker's clearance. It is `k_e / alpha`, the clearance divided by the fraction of overflow flux that becomes the measured marker. |

1. **Only `theta` survives with its meaning intact**, which is fortunate, since it is the estimand.
   Every statement this project has made about `K_m` or `k_e` is a statement about a composite.
   **[proven]**
2. **`K_p`, `c` and `alpha` enter through only two combinations**, rank 2 of 3, so one direction in
   them is invisible to the marker record however long it runs. Adding the state adds no parameter
   the marker can estimate. **[proven]**
3. **And this lands directly on Phase 9's precondition.** Phase 9 concluded that measuring `k_e` is
   not one improvement among several but the precondition for a statable requirement. **The
   micro-dose tracer study of S2b measures the marker's own clearance, `k_e`. The model's parameter is
   `k_e / alpha`.** So the study delivers the numerator and `alpha` is a new unmeasured nuisance
   sitting on exactly the quantity Phase 9 wanted to 2%. If `alpha` were 0.3, the model's constant is
   3.3 times the measured one. **[proven]** on the algebra; the size of `alpha` is unmeasured.
4. **And the physical reading, which makes `alpha` intuitive rather than surprising.** Eliminating
   `p` gives `c p` identically equal to `(k_e/alpha) m`. **So the paper's second channel is not a
   clearance of the marker at all.** It is the disposal of propionyl-CoA down the overflow route,
   expressed in marker units, and `alpha` appears because the model carries a flux out of a pool it
   does not represent, divided by the fraction of that flux that happens to be measured. Read that
   way the `alpha` dependence is not a technicality: it is the exchange rate between a pool the model
   omits and the one quantity it observes. **[proven]**
5. **`alpha` is measurable, and what measures it is the marker panel.** Each overflow marker's
   formation flux is its clearance times its level, so
   `alpha = k_e,3HP * m_3HP / sum_j (k_e,j * m_j)`: measuring the panel's levels **and** each
   marker's clearance closes the overflow budget and returns `alpha`. **[proven]**
6. **That estimator carries a CLOSURE ASSUMPTION, and it must be named.** It is right only if the
   measured markers **exhaust** the overflow routes. The panel R3 describes covers four: 3-HP, 2-MC,
   C3 and propionylglycine. Routes known to exist outside it include urinary free propionate,
   propionyl-CoA incorporated into odd-chain fatty acids, and any conjugate not on the assay.
   **The error is directional, not symmetric**: an unmeasured route makes the denominator short, so
   `alpha` comes out **too high** and the model's `k_e/alpha` **too low**. If the panel captures a
   fraction `f` of the overflow, the estimate of `alpha` is inflated by `1/f` and `k_e/alpha`
   deflated by `f`, so a panel capturing 70% of the overflow understates the model's constant by
   30%. **This is the same failure mode as R5's floor**: a ratio whose denominator is assumed
   complete. **[proven]** on the direction and the factor, **[judgement]** on which routes matter.
7. **And it is not one measurement.** That estimator needs the clearance of **every** marker in the
   panel, not only 3-HP, and R1 establishes that no human clearance exists for even one of them. So
   converting the micro-dose study's result into the model's parameter needs **as many tracer studies
   as there are panel markers**, each with its own synthesis, its own window and its own assay. The
   pricing in section P treats the `k_e` study as a single item and that is now wrong; the corrected
   row is in P.8. **[judgement]**
8. **Which answers Phase 10's original question from an unexpected direction.** The panel's value is
   not noise averaging, which is what the paper priced, and not the different functional forms, which
   is what R3 recommended it for. **It is that the sum of the overflow routes is what converts a
   measurable tracer clearance into the model's parameter.** That is structural information, it is the
   thing that was never priced, and it is now the panel's strongest justification. **[judgement]**

### R.3 What breaks the reduction, and it is the other markers **[measured]**

The reduction used `v_o = c p`. With a saturable overflow `v_o = c p/(1 + p/K_o)`, no single `K_m`
reproduces the two-state marker. The best-fitting `K_m` across capacity, and its residual:

| K_o | best-fitting K_m | worst relative misfit | reduces |
|---|---|---|---|
| 1e9, the linear limit | 100.00 | 7.5e-07 | **yes** |
| 1000 | 97.05 | 1.1e-02 | no |
| 100 | 69.48 | 1.6e-01 | no |
| 30 | 19.22 | 6.7e-01 | no |
| 10 | 1.58 | 9.5e-01 | no |

1. **The linear limit recovers `K_m = 100` exactly**, which is the reduction reappearing as a limit
   and is the check on this test rather than a result of it. **[measured]**
2. **As the overflow saturates, the one-state model absorbs the misfit into `K_m`** and does so
   badly: at `K_o = 100` the best `K_m` is 69 instead of 100, a 31% error in a constant, and the
   residual is still 16%. **A reader would see that as a calibration problem rather than as a missing
   state.** **[measured]**
3. **So the reduction holds for the marker the paper uses and fails for the ones R3 describes
   differently**, which is the convergence that makes the whole picture consistent: the other markers
   are informative precisely because they are not reductions of the same one-state form.
   **[judgement]**

### R.4 Where this sits in the reduction ladder **[judgement]**

The repository already maintains reductions rather than replacements: `REVIEW_LOG.md` **R22** records
the carnitine one, where `k` cancels from the reduced curve and survives only as the relaxation rate,
verified to 2.45e-16 over four orders of magnitude. **This belongs beside it, not in a category of its
own**, and gets the same three-part treatment.

| | the carnitine reduction (R22) | the propionyl-CoA reduction (here) |
|---|---|---|
| what it assumes | the conjugation pool is quasi-steady | both `p` and `m` quasi-steady, **and the overflow linear in `p`** |
| where it is exact | `k` cancels exactly; worst spread 2.45e-16 | residual 0 symbolically, 6.7e-15 numerically over 54 cells |
| what it costs when wrong | `k` still sets the relaxation rate, so a fast-pool assumption is checkable from dynamics | **a saturable overflow has no one-state form at all**: at `K_o = 100` the best `K_m` is 69 against 100, with a 16% residual |

**And the 16% residual is invisible as a missing state.** It presents as a `K_m` that will not
calibrate, which a reader fixes by moving `K_m`, which is exactly what the fit does. **That is the
fourth instance of this project's own pattern**: R14's curvature read as recovery, R15's knife edge
read as a shape limit, Phase 9 3b's needle read as identifiability, and now a missing state read as a
calibration error. In every case the arithmetic was right and the object it described was not the one
named. **[judgement]**


### R.5 What the verdict licenses, and what it does not

**Licensed.** The five dropped candidates are now expressible, because `p` exists: the breath test
attaches to `theta p/(K_p + p)`, and 2-MC, C3, C3/C2 and propionylglycine attach to `p` through R3's
functional forms. Steps 2 to 4 can be written against the two-state model.

**Not licensed, and this is the qualification the author should see before authorising more.** Adding
the state does not add identifiability by itself: R.2 point 2 shows the marker record sees no new
parameter. **The gain is entirely that new SENSORS become attachable**, and each of them brings its
own new nuisance, `alpha` being the first and the one that now sits on Phase 9's precondition.

**And one thing that must be re-examined rather than assumed to carry over.** Phase 9's entire σ
budget, and the S2b factor of 5.8, were derived with `k_e` meaning the marker's clearance. Under the
reduction that parameter is `k_e/alpha`. **The numbers do not change, since the reduction is exact and
the parameter is whatever the fits estimated, but their INTERPRETATION does**: "measure `k_e` to 2%"
becomes "know `k_e/alpha` to 2%", which is a harder requirement than the micro-dose study alone
satisfies. That re-reading is the first thing steps 2 to 4 should carry. **[judgement]**


---

## The hierarchical cohort fit, priced

`audit/q14_hierarchical.py`, re-checked by `audit/q11_phase10_recheck.py` check 3, 6 of 6.

Run first, ahead of steps 2 to 4, because the reduction gave the cohort route a large new advantage
from a direction neither strand was looking. The cohort route estimates the coefficient multiplying
`m` in the balance the records obey, which IS `k_e/alpha`, whatever that coefficient turns out to be
made of. So it never performs the conversion and never meets `alpha` at all, while the tracer route
now needs four micro-dose studies plus a closure assumption to reach the same quantity (P.8).

**The advantage is real and it does not rescue the route.** The cohort route cannot deliver the
number at any cohort size a PA study could reach, and if R1 is right about the between-patient spread
it cannot deliver it at any cohort size whatever.

### C.0 What is computed

Each patient contributes information about the shared `(ln K_m, ln k_e)` and about their own
`(ln theta_i, b_i)`. The individual pair is profiled out by the exact Schur complement,
`F_ss - F_si F_ii^-1 F_is`, with the patient's own `b` prior inside `F_ii`; the population prior of
0.8 on `ln k_e` is added once, at the population level, not once per patient. Capacities are spread
over 5 to 20% of nominal, which is where a PA cohort would sit. The trial's own cohort is 16.

The audit has not settled on one recommended schedule for the main record, so all three are priced,
with monthly the trial's own.

### C.1 The population value against cohort size **[measured]**

| n | monthly | weekly 7.3 d | every 2 d |
|---:|---:|---:|---:|
| 16 (the trial) | 0.3546 | 0.1845 | 0.1028 |
| 32 | 0.2627 | 0.1332 | 0.0733 |
| 64 | 0.2003 | 0.0955 | 0.0514 |
| 128 | 0.1432 | 0.0675 | 0.0365 |

And the cohort each target needs, by a doubling search and bisection on a quantity proven monotone
in `n`:

| schedule | n for 6% | n for 2% |
|:--|---:|---:|
| monthly | **723** | **6491** |
| weekly 7.3 d | 163 | 1461 |
| every 2 d | 48 | 429 |

**On the trial's own monthly record, the disease range's looser target needs 723 patients and its
tighter one needs 6491.** The largest cohort appearing anywhere in the six literature reports is 60
participants. So on the schedule the trial actually runs, the cohort route is not expensive, it is
**out of reach by one to two orders of magnitude**. **[measured]**, **[judgement]** on the reach.

Only the every-2-day record brings the looser target inside a plausible study, 48 patients. That is
48 patients giving a draw every two days for two years, 365 draws each, 17,520 draws in total. It is
not obviously cheaper than four micro-dose studies.

### C.2 The question that actually decides it, and the cohort loses on it **[measured]**

R1 states the disposal constant is patient-specific: 3-HP is renally secreted by OAT1, half of adults
with PA have an eGFR under 60, and the between-patient CV is "probably a third or more". A cohort fit
pins the POPULATION value. A patient's `theta_i` consumes THEIR OWN. The gap between the two is the
between-patient spread, and no cohort size closes it.

`sd(ln k_e_i)` for one patient once the cohort is fitted, monthly:

| between-patient CV | n=16 | n=64 | n=128 | floor at n=4096 | meets 6% |
|:--|---:|---:|---:|---:|:--|
| 0, shared exactly | 0.3546 | 0.2003 | 0.1432 | 0.0251 | yes |
| 10% | 0.3596 | 0.2218 | 0.1735 | **0.1027** | NO |
| 20% | 0.3945 | 0.2777 | 0.2418 | **0.1982** | NO |
| 33%, R1's figure | 0.4598 | 0.3693 | 0.3443 | **0.3166** | NO |
| 50% | 0.5564 | 0.4901 | 0.4731 | 0.4549 | NO |

**At R1's "a third or more" the floor is 0.3166 and the looser target is 0.06, so the cohort route
misses it by a factor of five with an infinite cohort.** Even a CV of 10%, well below what R1 argues
for, floors at 0.1027 and still misses. Only the exactly-shared row, which R1 gives reasons to
disbelieve, ever reaches a target, and then only the looser one.

**What the floor is made of, and why more patients cannot remove it.** Once the cohort is large the
population value is known essentially exactly, and the patient's own constant is known only to what
their own record says about it against a prior of the between-patient width. The cohort contributes
nothing further. **A larger cohort measures the population better and the patient no better.**

The cost of ignoring this is not a wider interval but a wrong one. A shared-`k_e` fit run on a
population whose CV is a third is misspecified, and 3a's whole apparatus prices variance, not that.

### C.3 The headline confirmed by fitting **[measured]**

Per the standing house rule. Cohorts simulated and fitted jointly, one shared `ln K_m`, one shared
`ln k_e`, a `(ln theta_i, b_i)` per patient, on the stacked residual. Spread is a sample standard
deviation, ddof = 1.

| n | schedule | Fisher | fitted sd | fitted/Fisher | median err |
|---:|:--|---:|---:|---:|---:|
| 16 | monthly | 0.3546 | 0.2189 | 0.617 | -0.0607 |
| 16 | monthly, start jittered | 0.3546 | 0.2211 | 0.624 | -0.0607 |
| 16 | weekly 7.3 d | 0.1845 | 0.1944 | 1.053 | +0.0141 |
| 64 | monthly | 0.2003 | 0.1769 | 0.883 | +0.0416 |
| 128 | monthly | 0.1432 | 0.1672 | 1.168 | +0.0029 |

**fitted/Fisher runs 0.617 to 1.168, so the cohort figures are a budget and not a curvature**, the
same verdict 3a.9 reached for the single-patient design.

**The first cell came out BELOW its own Fisher figure, which the house rule obliges naming rather
than reporting past.** That is either a biased estimator or an optimiser that never left its start,
and the two are distinguished by starting every shared and individual parameter a log-sd of 0.5 from
the truth. The jittered row reproduces it to 0.002 with an identical median error, so the optimum is
global and the gap is bias: at `n = 16` monthly the estimator trades a median error of -0.061 for a
narrower spread. **The Fisher figure is conservative there, not wrong, which means C.1's cohort sizes
are if anything generous to the cohort route and it fails anyway.** **[measured]**

### C.4 What it does to the requirement table **[measured]**

S2's map read ACROSS the `k_e` prior. This reads DOWN the routes that could supply it.
`sd(ln theta)`, recorded load, Fisher route, median over S2's 24 load seeds. Gate 0.20.

| capacity | route for the constant | `k_e` sd | monthly | weekly 7.3 d | every 2 d |
|---:|:--|---:|---:|---:|---:|
| 5% | population prior only | 0.800 | 3.068 | 2.194 | 1.424 |
| 5% | cohort n=16 monthly, shared | 0.355 | 1.578 | 1.411 | 1.132 |
| 5% | cohort n=64 monthly, shared | 0.200 | 0.969 | 0.904 | 0.816 |
| 5% | cohort n=64 monthly, CV 33% | 0.369 | 1.635 | 1.453 | 1.153 |
| 5% | four micro-dose studies + closure | 0.054 | 0.463 | 0.367 | 0.313 |
| 10% | population prior only | 0.800 | 1.581 | 1.224 | 0.842 |
| 10% | cohort n=16 monthly, shared | 0.355 | 0.787 | 0.725 | 0.617 |
| 10% | cohort n=64 monthly, shared | 0.200 | 0.485 | 0.451 | 0.419 |
| 10% | cohort n=64 monthly, CV 33% | 0.369 | 0.815 | 0.748 | 0.631 |
| 10% | four micro-dose studies + closure | 0.054 | 0.243 | 0.185 | **0.155** |
| 20% | population prior only | 0.800 | 0.812 | 0.707 | 0.549 |
| 20% | cohort n=16 monthly, shared | 0.355 | 0.391 | 0.369 | 0.340 |
| 20% | cohort n=64 monthly, shared | 0.200 | 0.245 | 0.224 | 0.213 |
| 20% | cohort n=64 monthly, CV 33% | 0.369 | 0.405 | 0.382 | 0.350 |
| 20% | four micro-dose studies + closure | 0.054 | **0.138** | **0.097** | **0.077** |

**Not one cohort row meets the gate anywhere in the table.** The best of them, n=64 shared at 20%
capacity on a 2-daily record, reaches 0.213 against a gate of 0.20 and misses. The micro-dose row
meets it at 20% capacity on every schedule and at 10% capacity on the 2-daily record.

**And the two rows to compare are the last two of each block**, because those are the two routes as
they would actually be run: a cohort of 64 with R1's spread, against four micro-dose studies. At 10%
capacity that is **0.815 against 0.243, a factor of 3.4**. At 20% it is 0.405 against 0.138, a factor
of 2.9. The cohort route is not close.

### C.5 The verdict, and the shape it has **[judgement]**

**The cohort route's advantage over the tracer route is real, large, and irrelevant.** It genuinely
never meets `alpha`: `k_e/alpha` is exactly the coefficient a marker record identifies, so the
closure assumption that P.8 showed the tracer route cannot escape simply does not arise. That is a
clean win on the axis the reduction exposed. It loses on a different axis entirely, arithmetic, and
loses there by a factor of three to five with the largest cohort the disease could supply.

So the answer to "how many patients at the recommended schedule to pin `k_e/alpha` to 2 to 6 per
cent" is: **723 and 6491 on the trial's monthly record if the constant is shared, and no number at
all if it is not.** The requirement table is unchanged by the cohort route because the cohort route
cannot enter it.

**What this does NOT license.** It does not restore the tracer route to where P.8 left it. Four
micro-dose studies plus a closure assumption is still four studies and an assumption, and C.4's
micro-dose rows are computed at `k_e` sd 0.054 as though the closure assumption were free, which it
is not. **Both routes to the model's second constant are now in difficulty, and the cohort route's
difficulty is the harder kind**: the tracer route is expensive and leans on an assumption that could
be tested by widening the panel, while the cohort route runs into a floor that no amount of money
moves.

**Where it leaves Phase 9's conclusion.** Strengthened, and in the direction Phase 9 already pointed.
Phase 9 concluded that measuring the second constant is the precondition for the requirement
framework existing at all. The obvious cheap alternative to measuring it, borrowing it from a cohort,
has now been priced and does not work. **[judgement]**

### C.6 The item this raises, filed rather than run

If the between-patient CV is what kills the cohort route, then the CV is the quantity worth measuring
first, and P.6's repeatability design, 10 patients at 12 draws, is already scoped to estimate exactly
that. **The order is now: measure the between-patient spread before committing to either route**,
since a CV near zero makes a 723-patient cohort merely infeasible rather than pointless, and a CV
near a third rules the cohort route out permanently and makes the four tracer studies the only route
there is. That is a cheap decision-relevant measurement standing in front of an expensive one.

---

## Steps 2 to 5: the five survivors, priced

`audit/q15_sensors.py`. The hierarchical strand is done, above.

### S.0 The framing that makes them comparable

The record's information about `(ln theta, ln K_m, b, ln k_e)` is a fixed 4x4 block, and every
surviving sensor does the same kind of thing to it: it supplies an independent measurement of ONE of
those four parameters. So the phase is not "which sensor is best" but **which parameter is worth
buying precision in, and at what price**, and the four sensors buy three different parameters.

| candidate | what it measures | the parameter it buys |
|---|---|---|
| 2, ex vivo PCC activity | the enzyme | `ln theta`, directly |
| 8, IV 13C3-3-HP | the marker's disposal | `ln k_e`, which is `k_e/alpha` (R.2) |
| 9, propionate turnover | `L + b` | `b`, the offset |
| 10, timed urine | `f_renal * k_e` | `ln k_e`, partially |

**And the framing shows why candidate 9 was never priced.** S2's requirement map has a `recorded`
load column and a `measured` load column and treats them as two different studies. **They are the two
endpoints of candidate 9**: `recorded` is a turnover measurement of infinite width, `measured` one of
zero width. Everything between was never computed, and the sensor lives there.

### S.1 Candidate 9 buys nothing, and that is the phase's clearest result **[measured]**

`sd(ln theta)`, monthly, median over 24 load seeds:

| capacity | no turnover study | at 10% | at 2% | **exact** |
|---:|---:|---:|---:|---:|
| 5% | 3.068 | 3.060 | 3.029 | **3.022** |
| 10% | 1.581 | 1.577 | 1.563 | **1.560** |
| 20% | 0.812 | 0.811 | 0.804 | **0.804** |

**Knowing `L + b` EXACTLY improves `sd(ln theta)` by 1 to 2%.** Not at a plausible precision, at
perfect precision, on every schedule and at every capacity. There is no precision at which this
sensor is worth buying, because there is no precision at which it does anything.

**Why, and it is worth stating because it also reprices S2's own two columns.** The offset is already
pinned by its own prior `s_b = 0.10 L`, and what limits `theta` is `k_e`, not the load. The load
channel was never the binding constraint. **So S2's `recorded` and `measured` load columns differ by
1 to 2% and are not two designs at all**, which no phase had noticed.

**The row the brief most wanted priced comes back flat, and that is an answer.** It removes a study
from the candidate set rather than adding one. **[measured]**

### S.2 The k_e sensor has a floor, and the monthly record sits under it **[measured]**

Candidates 8 and 10 both buy `ln k_e`. Sweeping the precision they deliver:

| capacity | schedule | 0.80 prior | 10% | 5% | 2% | **exact** |
|---:|:--|---:|---:|---:|---:|---:|
| 5% | monthly | 3.068 | 0.597 | 0.454 | 0.405 | **0.395** |
| 5% | every 2 d | 1.424 | 0.483 | 0.299 | 0.214 | **0.193** |
| 10% | monthly | 1.581 | 0.305 | 0.239 | 0.217 | **0.213** |
| 10% | every 2 d | 0.842 | 0.240 | 0.148 | 0.107 | **0.097** |
| 20% | monthly | 0.812 | 0.164 | 0.136 | 0.127 | **0.126** |

**At 5% and 10% capacity the monthly record misses the 0.20 gate even with `k_e` known exactly**,
0.395 and 0.213. No tracer study of any quality rescues it, because the limit is the record, not the
sensor. The 2-daily record passes at both, barely at 5% (0.193).

**So the tracer route's value is schedule-dependent in a way P.8's pricing did not show.** Four
micro-dose studies at 5.4% buy the 2-daily record a pass at 10% capacity and buy the monthly record
nothing at all below 20%. **The schedule decision comes before the sensor decision.** **[measured]**

**And candidate 10 cannot reach even that**, because `f_renal` is inside it:

| `f_renal` known to | floor on `ln k_e` | 10% cap, monthly | 20% cap, monthly |
|:--|---:|---:|---:|
| a factor of 3 (95% interval) | 0.275 | 0.631 | 0.315 |
| a factor of 2 (95% interval) | 0.174 | 0.435 | 0.223 |
| +/- 5% | 0.054 | 0.243 | 0.138 |

R1 gives no `f_renal` for 3-HP at all. To match candidate 8's 5.4%, the renal fraction would have to
be known to about 5%, which is a harder measurement than the clearance it is standing in for.
**Candidate 10 is dominated by candidate 8 on its own terms** and survives only as the bound step 1
called it. **[measured]**

### S.3 The result the project should be most uncomfortable with **[measured]**

An ex vivo assay measures `ln theta` directly, so it ALONE gives `sd(ln theta)` equal to its own
precision, with no record at all. What does the two-year marker record add on top of it?

| capacity | assay alone | + monthly record | + weekly | + 2-daily |
|---:|---:|---:|---:|---:|
| 5%, 20% assay | 0.2000 | 0.1996 (**0.2%**) | 0.1992 (0.4%) | 0.1981 (1.0%) |
| 10%, 20% assay | 0.2000 | 0.1984 (**0.8%**) | 0.1974 (1.3%) | 0.1946 (2.7%) |
| 20%, 20% assay | 0.2000 | 0.1942 (**2.9%**) | 0.1924 (3.8%) | 0.1879 (6.0%) |
| 20%, 30% assay | 0.3000 | 0.2814 (6.2%) | 0.2762 (7.9%) | 0.2632 (**12.3%**) |

**Two years of monthly plasma 3-HP, added to a single ex vivo assay of ordinary quality, improves the
estimate of residual capacity by between 0.2 and 2.9 per cent.** The best cell anywhere in the sweep
is 12.3%, and it needs a poor assay, the densest schedule and the top of the disease range together.

**This is not an argument for the assay and it is not a refutation of the record.** The ex vivo route
was priced in q5 and carries its own problem, that a fibroblast's activity is not the patient's in
vivo capacity, and nothing here measures that transfer. **What it is, is the marker record's
information content stated on an absolute scale instead of a relative one.** Every earlier phase
compared schedules against schedules and sensors against sensors. This compares the whole record
against one blood draw sent to a reference laboratory, and the record is worth a few per cent.

**Where that leaves the project's question.** "Can residual capacity be identified from plasma 3-HP
under mRNA-3927" now has a sharper form: **it can, to within the ex vivo assay's own precision, and
almost entirely because of the assay.** If no assay is available the record cannot do it alone at 5
or 10% capacity under any sensor (S.2). If an assay is available the record is nearly redundant. The
band where the record is the decisive instrument is narrow and this phase has not found it.
**[measured]**, **[judgement]** on the reading.

### S.4 The combinations, and which are redundant **[measured]**

Candidate 2 at 20%, candidate 8 at 5.4%, candidate 9 at 10%:

| combination | 5%, monthly | 5%, 2-daily | 10%, monthly | 10%, 2-daily |
|:--|---:|---:|---:|---:|
| record alone | 3.068 | 1.424 | 1.581 | 0.842 |
| c2 only | 0.200 | 0.198 | 0.198 | 0.195 |
| c8 only | 0.463 | 0.313 | 0.243 | 0.155 |
| c9 only | 3.060 | 1.414 | 1.577 | 0.836 |
| **c2 + c8** | **0.184** | **0.168** | **0.154** | **0.122** |
| c2 + c9 | 0.200 | 0.198 | 0.198 | 0.195 |
| c8 + c9 | 0.449 | 0.307 | 0.237 | 0.152 |
| c2 + c8 + c9 | 0.183 | 0.168 | 0.153 | 0.121 |

**Only one pair is worth anything, `c2 + c8`**, and even it improves on `c2` alone by 8 to 22%.
**Every combination containing c9 equals the same combination without it to three figures.** Adding
candidate 9 to candidate 2 changes nothing at all, 0.200 to 0.200.

### S.5 Both verdicts confirmed by fitting **[measured]**

Per the standing house rule, since both decide a verdict. Robust spread (1.4826 x MAD) with the
collapse count, 96 records per cell, `centre` nuisance convention, monthly.

| capacity | configuration | Fisher | fitted | fitted/Fisher | collapses |
|---:|:--|---:|---:|---:|---:|
| 10% | record alone | 1.581 | 1.340 | 0.847 | 29/96 |
| 10% | + c9 exact turnover | 1.560 | 1.196 | 0.767 | 28/96 |
| 10% | + c2 at 20% | 0.198 | 0.203 | 1.022 | 0/96 |
| 10% | + c2 at 20% + c9 exact | 0.198 | 0.209 | 1.056 | 0/96 |
| 20% | record alone | 0.812 | 0.713 | 0.878 | 20/96 |
| 20% | + c9 exact turnover | 0.804 | 0.736 | 0.915 | 21/96 |
| 20% | + c2 at 20% | 0.194 | 0.188 | 0.967 | 0/96 |
| 20% | + c2 at 20% + c9 exact | 0.194 | 0.193 | 0.993 | 0/96 |

**fitted/Fisher runs 0.767 to 1.056, so both verdicts stand on the fitted footing.**

**On candidate 9, the fitted rows do not match the Fisher rows as tidily as the Fisher rows match each
other, and that is worth stating rather than smoothing.** At 10% the c9 row is 11% better fitted where
Fisher said 1.3%; at 20% it is 3% WORSE fitted, which cannot be a real improvement. At 96 reps with
20 to 29 collapses the robust spread carries its own noise of about that size. **So the honest reading
is that the c9 rows are indistinguishable from the rows above them on both footings**, which is the
verdict, just not by the clean margin the Fisher table suggested.

**On candidate 2 the two footings agree closely**, 0.203 against 0.198 and 0.188 against 0.194, with
**zero collapses in every c2 row against 20 to 29 without it**. That is a second thing the assay buys
which the spread alone does not show: it removes the collapse mode entirely. **[measured]**

### S.6 What steps 2 to 5 return

1. **Candidate 9, propionate turnover, is dropped on measurement.** Exact knowledge of `L + b` is
   worth 1 to 2%, and S2's two load columns collapse into one.
2. **Candidate 10, timed urine, is dominated by candidate 8** unless `f_renal` is known to about 5%,
   which R1 does not supply and which is harder than the measurement it substitutes for.
3. **Candidate 8's value is schedule-dependent**: it buys the 2-daily record a pass at 10% capacity
   and buys the monthly record nothing below 20% capacity, because at 5 and 10% the monthly record
   misses the gate with `k_e` known exactly.
4. **Candidate 2, ex vivo activity, dominates everything**, and the uncomfortable corollary is that
   the two-year marker record adds 0.2 to 2.9% on top of it. The record's remaining distinct value is
   removing nothing the assay does not already remove, except that the assay removes the collapse
   mode the record alone has at 20 to 30%.
5. **Only `c2 + c8` is a combination worth buying**, and the phase's practical recommendation is
   therefore an ex vivo assay plus the densest affordable record, with the tracer study bought only
   if the record is dense.

**[measured]** on 1 to 4, **[judgement]** on 5.
