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

## Steps 2 to 4, and the hierarchical family

Not started. Step 1 is a gate and this file is the report at it.
