# Identifiability of a therapeutic threshold in propionic acidemia

A control-theoretic study of individualised dosing for an mRNA
enzyme-replacement therapy in a monogenic metabolic disease.

**The clinical decision** is a dosing interval. **The quantity that determines
it** is the patient's residual enzymatic capacity `theta`, which is never
measured. **The question** is how much of it can be known, from data it is
ethical to collect, before the decision has to be made.

---

## The model, at three resolutions

One model. The letters A, B and C are gone, because two of the three were the same model
read at different resolutions and the third was a hypothesis about one term.

```
dm/dt = L(t) + L0 - (theta + Delta(t)) * g(m) - alpha(m),    g(m) = m/(Km + m)
```

`g` is the **Michaelis-Menten saturation function** and `theta` **absorbs** the maximal
velocity: it *is* the patient's residual `V_max`, not a fraction of a population one.
`Delta(t)` is the therapeutic increment, known from the dosing schedule. `L0` is the
**unrecorded** part of the propiogenic load, carried as an unknown — which turns out to be
what makes the problem well posed at all.

| resolution | obtained by | what it is for |
|---|---|---|
| **the dynamic model** | as written | the trajectory, the Gramian, the excitation condition |
| **the equilibrium reduction** | `dm/dt = 0` | the fold, `theta_crit`, the safety margin |
| **the local identification model** | `alpha(m) -> ke*m`, over the visited band | the fibre, the margin theorem, the design |

**The third is NOT a reduction of the first**, and the check is one line. The tangent
intercept of `theta*g + alpha` at the operating point is
`a = theta*g^2 + [alpha - m*alpha']`, so the local model's identity `a = theta*g(m_bar)^2`
holds only when the alternative channel is a straight line through the origin. At the
illustrative constants the second bracket is 374 against the enzymatic term's 61: the
alternative channel owns six-sevenths of the intercept. The local model is therefore the
*parent with one term replaced*, and it is retained because that replacement is the **most
favourable** case: a curved `alpha` contributes to both fitted coefficients and closes the
fibre, a straight one contributes to neither and opens it.

A fourth structure, `alpha(m) = ke*h(m)` with a ceiling instead of a fold, survives as one
remark. It is the null against which the fold is testable, and it prices biological
fidelity in the identification (a factor of ~450 in condition number).

---

## Results

### The chain

**The decision is a sign, not a value.** What matters clinically is the sign of
`mu = theta - theta_crit` for a patient under that patient's own load.

**A band determines the affine fit and nothing more.** Over the range a managed patient
visits, total clearance is indistinguishable from a line at assay resolution.

**The fibre's width is observable, in closed form.**
`span = (1 + b*m_bar/a)^2 = (C(m_bar)/a)^2` — the recorded load over the fitted intercept,
squared. 169 at the illustrative values, and **computable per patient without knowing a
single parameter**.

**The sign survives what the value does not — asymmetrically.** The `Km^2` terms cancel
along the fibre, so `mu` is linear where `theta` is quadratic. Simulating the accessible
design, `theta_hat` spans three orders of magnitude while the sign of `mu` is right 99.8%
of the time for a patient *with* margin and wrong **53%** of the time for one without, at
the same distance from the boundary. **The estimator fails towards declaring patients
safe**, because the abscissa noise inflates `theta_hat` and `mu_hat` inherits it.

**Both ends of the saturation curve are degenerate, and they fail differently.** Below
`Km` only the lumped ratio `theta/Km + ke` is identifiable. Above it, `theta` enters only
through `L0 - theta`, so with the load partly unrecorded the target parameter is confounded
with the load. Losing a nuisance parameter is survivable; losing the target is not.

**The excitation is maximised straddling `Km`, with a closed form.** With
`phi(m) = Km/(Km+m)`, the band `[lo, hi]` delivers `E = phi(lo) - phi(hi)`, maximised at
`sqrt(lo*hi) = Km` exactly, where `E_max(r) = (sqrt(r)-1)/(sqrt(r)+1)` for fold range `r`.
**The requirement is a band around `Km`, not a threshold above it.**

**Time does not buy it. Placement does.** An augmented-state estimator with **no process
noise on `theta`**, under the trial's biweekly schedule over two years:

| design | n | sd(ln theta) |
|---|---|---|
| six-monthly clinic visits | 4 | 0.545 |
| one dose cycle, every 2 days | 7 | 0.540 |
| every 2 days for two years | 365 | 0.489 |
| *the same 365 samples straddling `Km`* | 365 | **0.004** |

Ninety times the sampling buys a tenth of a decade. Recording the load ten times better
(0.5%) buys another 5%. Placement buys a factor of 122. The dense design was expected to
collapse the covariance and **did not** — that is the result, not a tuning failure.

**Only a decompensation approaches the excitation required.** Against a requirement of
`E ~ 0.50`: inter-dose recovery 0.049, carnitine dose steps 0.062, antibiotics 0.034,
supply interruptions 0.047, labelled-propionate bolus 0.076 — and a decompensation 0.354.
The routine sources deliver a fifth to a quarter of what their own fold range *could*,
because they sit an order of magnitude below `Km`. **The excitation the decision needs is
available only from the event the therapy exists to prevent.**

**No cohort rescues it.** At the accessible design the sign error near the boundary is 41%,
and correcting a fraction for misclassification inflates its variance by `(1-2p)^-2`.
Placing the fraction with `mu < 0` within ten points needs **about 3,100 patients**. The
phase 1/2 trial enrolled 16.

### Supporting results

**The determinant is blind to the channel — given `Km`, and only given `Km`.** At a matched
trajectory the `(theta, ke)` Gramians differ by under 5%, always favouring the therapy,
because the sensitivity forcings are functions of the state alone and only the shared decay
rate differs: `J_therapy - J_carnitine = -kappa*u*g(m) <= 0`. With `Km` unknown — which is
the paper's own situation — it fails by six orders of magnitude in `det G` and 46x to 807x
in `sd(ln theta)`: the `Km` sensitivity carries `theta + Delta(t)`, a prefactor only the
therapy modulates. **Carnitine buys the excursion; it does not buy the curvature.**

**The intercept variance has a floor no band width beats.** `Var(a_hat)` has a `1/n` term
the published derivation dropped, and `m >= 0` caps a band at twice the operating point, so
`sd(theta)/theta >= c(n)*sigma/(a*sqrt(n))` with `c(n) <= 2`. The design question splits
into a floor on the number of perturbation events and, above it, an excursion requirement.

**Assay variability is the larger noise term whenever the load is well recorded.** The two
balance at `CV = 1.08 x (load error)`. A marker panel averages the assay term as `J^-1/2`
and cannot touch the load term, so **buy the load record before the panel**.

**`theta_crit` is nearly a function of one clinical claim.** Holding the calibration that
fixes `R0`, sweeping `Kc`, `Mi` and `rho_0` over their whole ranges moves `theta_crit` by a
factor of 1.39. **The margin is not so protected**: it runs 509 to 3275 over an 8x range in
`Mi`, so Figure 16's 954 is an order of magnitude, not a value.

**The marker is 3-HP.** RR 0.7 (95% CI 0.6–0.9, p < 0.01) per 50% *lower* concentration — a
risk marker in the expected direction. C3 is excluded structurally (carnitine raises it
while lowering the organic acids, reversing the perturbation's sign); 2-MC on the endpoint
(RR 1.0, p = 0.84), its strength being as a correlate of ammonia *during* an event, which
is the arm the model excludes.

**The conflation is measured, not argued.** Propionate oxidation by breath test is
near-identical in patients whose serum markers differ widely, repeating at CV 6.9%
(controls) and 12.8% (patients). The cohort is MMA, not PA, and no PA breath-test series
has been located.

## Layout

```
paper/          LaTeX source, split by section
  main.tex        build target
  sections/       one file per section; 06_results pulls in the figures
  note_structural_identifiability.tex   standalone note with the propositions
src/
  common/models.py      the kinetics and the three clearance structures
  common/params.py      every constant, with a source column; generates tab_params
  common/dynamics.py    the sensitivity system and the observability Gramian
  common/excitation.py  phi, E, E_max, and the design quantities
  model_C/run.py  THE DYNAMIC MODEL and its equilibrium reduction
  model_A/run.py  THE LOCAL IDENTIFICATION MODEL: the fibre, the margin, the floor
  model_B/run.py  the ceiling hypothesis (one remark in the paper)
  design/run.py   the excitation figure, the covariance collapse, the decision
                  characteristic, the excursion table, the data specification
symbolic/       structural identifiability, three independent tools
  julia/          StructuralIdentifiability.jl, 8 cases
  maple/          SIAN
  daisy/          DAISY
tools/
  fetch_refs.py   builds paper/refs.bib from CrossRef and PubMed
  claims_to_tex.py  renders paper/claims.csv as the claim register appendix
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
make refs        # fetch paper/refs.bib, render the claim register
make paper       # build paper/main.pdf (runs refs first)
make test        # assert the structural claims still hold
make audit       # non-zero while any claim is contradicted or unsourced,
                 # any citekey is missing from refs.bib, or any \cite in the
                 # paper has no row in the register
make all
```

`make audit` currently **fails**, on 14 unresolved claims. That is the
state of the audit, not a broken build. See `REFERENCES.md`: no citation
in this repository is typed by hand, and no factual claim goes into the
paper without a row in `paper/claims.csv` saying whether the source
actually says it. Seven rows are contradicted by the work they cite, and
those sentences are left standing until each is deliberately revised.

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
   - **It uses the recurrent-event structure instead of flattening it.**
     Patients do not have one event each: 49 of 50 participants in the natural
     history cohort reported at least one, at a mean annualised frequency of
     1.0, and the rate varies enough between them that the study carries a
     participant-level random effect. A hazard formulation uses every event and
     every event-free interval a patient contributes, and makes the spread
     between patients a parameter rather than noise.
   - **It makes the count model derived rather than assumed.** At present a
     Poisson or negative-binomial count model would be imposed on the admission
     counts by fiat. Under this construction the counting process falls out of the
     arrival process and the crossing probability, and its overdispersion has a
     mechanism behind it.

   The construction is partly pre-empted and the precedent belongs in the
   statement of it: the natural history study already fits a hazard of this
   shape, a GAMM on 3-HP with a smooth term in days since the previous event and
   a participant random effect. What is proposed is not that machinery but the
   substitution of a mechanistic `theta`, with the critical amplitude supplied by
   Model C, for a measured biomarker as the covariate carrying the risk.

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
