# Identifiability of a therapeutic threshold in propionic acidemia

A control-theoretic study of individualised dosing for an mRNA
enzyme-replacement therapy in a monogenic metabolic disease.

**The clinical decision** is a dosing interval. **The quantity that determines
it** is the patient's residual enzymatic capacity `theta`, which is never
measured. **The question** is how much of it can be known, from data it is
ethical to collect, before the decision has to be made.

---

## The model, in one line

```
dm/dt = s(t) - theta * g(m) - [ alternative clearance ] ,    g(m) = m/(Km + m)
```

`g` is the **Michaelis-Menten saturation function**. In the standard biochemical form the
velocity of a single-substrate enzymatic reaction is `v = V_max * [S] / (Km + [S])`, so
`theta` is the patient's residual `V_max`, `Km` is the Michaelis constant, and `g(m)` is
the fractional saturation. The same form is the Monod function in microbial kinetics, and
the Hill function with `n = 1`.

Three hypotheses about the bracketed term, each answering a different question:

| | alternative clearance | equilibria | question |
|---|---|---|---|
| **A** | `ke * m`, linear and unbounded | one, never lost | can `theta` be identified? |
| **B** | `ke * h(m)`, Michaelis-Menten too | one, lost above the ceiling | under what condition is a tolerable steady state lost? |
| **C** | `k * c(t) * m`, carnitine depleted | two, colliding at a fold | why does recovery require the load driven *below* baseline, not just the trigger removed? |

They are alternative hypotheses, not stages of a development, and they are not
ranked by realism. Each is here for a different and specific reason.

**C is the one that makes `theta_crit` well posed.** This is the point that
matters most and it is easy to miss, because C is usually introduced as "the
most biologically faithful of the three", which undersells it. In A and B the
critical capacity has to be *defined* by a convention:

```
theta_crit = (s_peak - ke*m_bar) / g(m_bar)
```

with `m_bar` a clinical judgement about what metabolite level is tolerable and
`s_peak` an unmeasured illness amplitude. In C it is *derived*: it is the theta
at which the attractor and the separatrix merge at a saddle-node, and no
tolerability convention appears anywhere in the construction. A cannot offer
this because its equilibrium is never lost; B cannot because its equilibrium
escapes to infinity instead of colliding with a second one. Only C has a second
equilibrium for the first to collide with.

**A is the most favourable case, not the most accurate.** It is the working
model for identifiability precisely because its two clearance channels have
maximally different functional forms. If identification fails there it fails a
fortiori elsewhere. Its known defect is that it predicts no crisis at all.

**B buys the load ceiling and pays for it.** It is what is needed to say
anything about a steady state being lost, and the price is that both channels
acquire the same functional form and become far harder to separate.

---

## Results so far

**Model A is globally structurally identifiable.** The obstruction is
conditioning, not structure. Knowing `m(t)` and `s(t)` determines the rational
function `F(m) = theta*g(m) + ke*m`, whose coefficients recover all three
parameters.

**Three distinct metabolite levels are required.** A patient held at a single
steady state gives rank 1.

**Conditioning degrades quadratically** as the patient is stabilised. A patient
varying by a few per cent between visits presents a problem conditioned at 1e5
to 1e6, beyond assay precision.

**Cohort size alone never helps.** N patients at steady state sharing the
population parameters give N equations in N+2 unknowns. Deficiency 2 for every N.

**Two anchors close it exactly.** A null genotype supplies `ke` directly, since
its equation loses the enzymatic term. One patient observed at two distinct
loads supplies `Km`. Both are present in the Saudi cohort.

**Biological fidelity costs identifiability.** At a 10% excursion the condition
number rises from 1.1e5 in Model A to 4.8e7 in Model B, a factor of ~450, because
saturating both channels gives them the same functional form.

**Models B and C differ in a checkable prediction.** Under B, cessation of the
precipitating illness suffices for recovery. Under C, the load must be driven
below baseline. This is answerable from clinical notes with no new measurement.

**`theta_crit` can be defined intrinsically, and only in Model C.** It is the
theta at which, at a given load, the maximum of the clearance curve falls to
that load: attractor and separatrix merge at the saddle-node. Above it two
equilibria exist, below it none. `m_bar` does not appear. The construction
removes the clinical convention from the definition of the threshold; it does
not remove `s_peak`, which survives in the separate question of whether an
episode is large enough to cross the separatrix. Only the hazard reformulation
in open item 2 removes that.

**The attractor-to-separatrix distance is the safety margin**, and it was
computed nowhere in this repository before. It says how far a patient can be
pushed before the crisis becomes self-sustaining, closes continuously to zero
at `theta_crit`, and diverges as theta approaches the baseline load, where the
patient becomes unconditionally safe at that load. At the illustrative
`theta = 300` it is 398, against `theta_crit = 167`. Neither A nor B admits
this quantity at all.

**Model B gives a condition, not a frequency.** The answer it returns is the
inequality `s > theta + ke`, written in two quantities (`s_peak`, `m_bar`) that
are never measured, and it is silent on how often the threshold is crossed.
Frequency is the one thing the records hold reliably. Closing that gap needs
the hazard reformulation in open item 2, not a better threshold.

**The identification determinant does not depend on which channel is perturbed.**
Perturb `k_e` with carnitine or `theta` with the therapy: the 2x2 system for
`(theta, ke)` has the same coefficient matrix either way, and

```
det = g(m)*m' - g(m')*m = m*m'*(m'-m) / ((Km+m)*(Km+m'))
```

depends only on the two metabolite levels. The perturbation enters the
right-hand side, not the matrix. So the analysis can be done now, on the input
that is in every chart, and it transfers unchanged to the input that is in no
chart. What the therapy would add is statistical, not structural: an exogenous
schedule instead of one confounded by indication.

---

## Layout

```
paper/          LaTeX source, split by section
  main.tex        build target
  sections/       one file per section; 06_results pulls in the figures
  note_structural_identifiability.tex   standalone note with the propositions
src/
  common/         model definitions and plotting style, shared
  model_A/run.py  identifiability, conditioning, population anchors
  model_B/run.py  load ceiling, threshold, cost of the correction
  model_C/run.py  bistability, separatrix, intrinsic theta_crit, safety margin
symbolic/       structural identifiability, three independent tools
  julia/          StructuralIdentifiability.jl, 8 cases
  maple/          SIAN
  daisy/          DAISY
results/
  figures/        every figure in the paper, regenerated by `make figures`
  tables/         LaTeX tables generated from code
tests/          consistency checks asserting the claims in the paper
```

No figure in the paper is drawn by hand. Everything under `results/` is
regenerated from `src/`.

---

## Running

```bash
pip install -r requirements.txt

make figures     # regenerate all figures and tables (seconds)
make paper       # build paper/main.pdf
make test        # assert the structural claims still hold
make all
```

Structural identifiability needs Julia and is slower:

```bash
julia -e 'using Pkg; Pkg.add("StructuralIdentifiability")'
julia symbolic/julia/run_identifiability.jl
```

The two cases that matter are **A2** (constant input) and **A4** (null-genotype
anchor). A2 should return the degenerate combination `theta/Km + ke` from
`find_identifiable_functions` rather than the three parameters separately; A4
should show `ke` becoming identifiable. Together they are the mechanical
confirmation of the hand analysis.

---

## Open items

Roughly in order of how much they would change the conclusions.

1. **Uncertainty in `s`.** Both anchors assume the substrate load is known. It is
   only partially recorded. The current amplification factor is ~1.4; bounding it
   properly is the first extension.

2. **Reformulating the crisis as a first-passage problem.** *Probably the single
   most valuable next step, and not started.* Model C's intrinsic threshold has
   already removed `m_bar` from the definition of `theta_crit`. `s_peak` is still
   there, and no sharper deterministic threshold will remove it, because the
   quantity being asked for is not a threshold at all. The concrete structure:

   ```
   lambda_MDE(theta) = nu * P( episode amplitude > critical amplitude | theta )
   ```

   where `nu` is the arrival rate of precipitating events (infection, fasting,
   surgery, dietary error) and the amplitude distribution is **estimated from the
   data rather than postulated**. The critical amplitude is the smallest excursion
   in `s` that carries the patient across the separatrix from the attractor, which
   Model C already computes: it is the margin plotted in `fig_C5_margin`, expressed
   in load rather than in metabolite.

   What this buys, point by point:

   - **It removes `s_peak` from the definition of the threshold.** A single
     worst-case amplitude is replaced by a distribution, so nothing has to be
     assumed about the severity of the worst illness a patient will meet.
   - **It uses the reliable observable.** The endpoint becomes dated hospital
     admissions, which are coded consistently, instead of sparse metabolite
     measurements whose relation to the enzyme's substrate is itself an
     assumption.
   - **It treats the ~70% of patients with no events correctly.** They are
     right-censored observations carrying real information about `theta`, not
     missing data to be dropped or zeros to be explained away.
   - **It makes the count model derived rather than assumed.** At present a
     Poisson or negative-binomial count model would be imposed on the admission
     counts by fiat. Under this construction the counting process falls out of the
     arrival process and the crossing probability, and its overdispersion has a
     mechanism behind it.

   It would also complete the answer to Model B's question. Model B currently
   returns a condition and the records contain a frequency; `lambda_MDE(theta)` is
   a frequency, and it is the quantity the two can be compared in.

3. **Optimal experimental design under a state constraint.** Conditioning scales
   as the inverse square of the excursion, and the excursion is bounded above by
   safety. How much information can be bought without approaching `theta_crit`?

4. **Identifiability of the decision rather than the parameters.** The clinical
   quantity is the sign of `theta - theta_crit`. The two parameters enter that
   margin with *opposite* signs, whereas the identifiable combination carries them
   with the same sign, so errors accumulate rather than cancel. The margin is
   plausibly *worse* conditioned than either parameter. To be established.

5. **Invariance to monotone rescaling.** If every variant's activity is distorted
   by the same increasing map, and `theta_crit` is estimated empirically on that
   same scale, the decision is preserved. This would mean an assay need only be
   *consistent*, not absolutely calibrated. Fails if the distortion differs by
   mechanism class.

6. **Multi-parameter estimator design.** The scalar Immersion and Invariance
   construction needs a sign condition. With three parameters entering at once it
   becomes a matrix condition and existence is not guaranteed.

7. **Defending carnitine as an instrument.** The selective input in the records
   is endogenous: the dose responds to the patient's state. The quasi-exogenous
   variation the argument leans on (weight-based recalculation after growth,
   supply interruptions, changes of formulation) has to be argued from the
   records patient by patient. This is a question about the charts, not about
   the model.

---

## Caveats a reader should hold

Parameter values are illustrative, not fitted. Their role is to fix the
qualitative structure and the orders of magnitude.

`m` in the model is the measured plasma biomarker; the enzyme's actual substrate
is intracellular propionyl-CoA, which is never measured. The map between them is
an assumption.

`k_e` is treated as a population parameter on the grounds that it does not involve
PCC, but carnitine is also administered therapeutically at patient-specific doses,
so in real data it is partly clinician-manipulated. That manipulation is the
working perturbation of the analysis and its weakness at the same time: the dose
is adjusted in response to the patient's state, so it is confounded by indication
and only its quasi-exogenous component can be used.

Model C's depletion mechanism is a hypothesis. A plain two-state depletion model
gives saturation, not decline; decline requires that replenishment itself be
impaired as the metabolite accumulates. The check is whether free carnitine
collapses during decompensation, which is measured routinely.
