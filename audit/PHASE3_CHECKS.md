# Phase 3.0.5 and 3.0.6 — reported before any prose.

**3.0.5's gate fires, and against Option 2.** At the patient's own operating point the
model-based rate is not merely no better than reading the rate off the marker — it is
**13.5× worse**. The model's only exclusive capability is counterfactual prediction.

Reproduced by `audit/e3_baseline_and_scope.py`.

---

## 3.0.5 (i) — at the operating point, the direct read-off wins outright

Monthly design, n = 25, assay CV 5%.

| route | sd(ln m) | × β | 95% width |
|---|---|---|---|
| **direct read-off, no model, no θ** | 0.0100 | **0.0051** | **1.02×** |
| model-based, via θ̂ | 0.1354 | **0.0697** | **1.31×** |

The model reconstructs a level that was **measured directly**, and pays parameter
uncertainty to do it. Averaging 25 monthly assays gives the level to 1%; the fitted
parameters give it to 13.5%.

This is not a close call and it should not be argued away. The natural history study
already publishes the association; applying it to a measured level requires nothing this
paper contributes.

## 3.0.5 (ii) — the model's one exclusive capability is counterfactual

| counterfactual | model sd(ln rate) | direct |
|---|---|---|
| load 1.5× baseline | 0.0711 | **no answer possible** |
| load 3.0× baseline | 0.0906 | **no answer possible** |

"What would this patient's rate be if the dose were halved, or if their load rose by
half?" has no observation behind it, so the direct read-off is silent. **This is the only
column where the model has something the data do not**, and whether it is worth a paper
depends on whether that counterfactual is the clinical question. It is a different
question from the one the paper opened with.

## 3.0.5 (iii) — the map's own uncertainty dominates both routes

RR 0.70 (95% CI 0.60–0.90) per halving gives β = −log₂(RR) = **0.5146, 95% CI 0.152 to
0.737 — the exponent is known only to a factor of 4.8.** Neither route can reduce it:

| m / m_ref | contribution from the map alone |
|---|---|
| 1.25 | 0.0333 |
| 1.50 | 0.0605 |
| **2.00** | **0.1034** |
| 3.00 | 0.1639 |
| 5.00 | 0.2402 |

At a doubling of the marker the map alone contributes 0.103 — **larger than the model's
entire parameter uncertainty of 0.070.** Any argument that the model's precision matters
is an argument about the second-largest term.

## 3.0.5 (iv) — and the map has no anchor, which blocks both routes equally

`rate(m) = rate_ref · (m/m_ref)^β` needs both constants. The study reports a cohort mean
annualised MDE frequency of 1.0 (SD 1.6), so `rate_ref` exists. **It reports no
concentration range and no reference 3-HP value, so `m_ref` does not.** Without it neither
route yields a rate — only a ratio between two levels. The model adds nothing to that.

---

## 3.0.6 — where the estimand fails, in observable units

The marker rises towards a ceiling as capacity falls: as θ → 0, m → L₀/k_e = **6.45× the
managed level**. The elasticity dies with it.

| θ (% of nominal) | marker, × managed | \|elasticity\| |
|---|---|---|
| 100% | 1.00 | 0.915 |
| 50% | 1.82 | 0.807 |
| 20% | 3.41 | 0.537 |
| 10% | 4.57 | 0.321 |
| 5% | 5.40 | 0.173 |
| 2% | 6.00 | 0.072 |
| 1% | 6.22 | 0.036 |

| elasticity falls below | at θ | marker reads |
|---|---|---|
| 0.50 | 17.9% of nominal | **3.6× managed** |
| 0.25 | 7.5% | **5.0× managed** |
| 0.10 | 2.8% | **5.8× managed** |
| 0.05 | 1.4% | **6.1× managed** |

**The clinician-facing statement, and it belongs wherever the estimand is claimed:** once a
patient's managed marker sits above about three to four times the reference level, a
further fall in residual capacity barely moves it. Two patients differing tenfold in
capacity are then indistinguishable on the marker.

**But the failure is asymmetric, and that asymmetry is the one point favourable to
Option 2.** What dies is the estimate of *capacity*. The *rate* prediction does not,
because the rate is read off the marker and the marker is what is measured. So the scope
limitation is on reading θ out of the marker for severe patients — which is an argument
for the rate estimand and against the capacity one, in exactly the population the paper
opened by caring about.

---

## The rate map's provenance, to be recorded in the text as-is

Verified against Results, not the abstract:

- **Unadjusted GAMM.** Three terms: 3-HP concentration, a smooth in days since the previous
  MDE, a per-participant random-effect smooth. No adjustment for age, disease onset, or
  treatment.
- **Unresolved whether 3-HP entered as baseline-only or time-varying.** Repeated measures
  exist (6-monthly to month 36, plus during MDEs), and a figure caption says "only baseline
  biomarkers are presented", but the source does not say which governs the GAMM.
- **No concentration range**, no explicit log-linear specification, so the per-halving
  coefficient cannot be extrapolated beyond a range that is itself unstated.
- **No N** for the regression models.
- **The authors' own limitation:** the biomarkers "may not represent their overall severity",
  and "further investigation is required to determine whether they are sufficient to
  reflect the risk of MDE in future clinical trials."

Direction is as expected for a toxic metabolite, so this is not a paradoxical protective
RR. But it is unadjusted, and severity and treatment intensity both act in the same
direction, so the association is an upper bound on any causal reading.

**This is not an established dose–response and must not be presented as one.**

---

## Recommendation

**Return to Option 1**, on the evidence above, and narrow the paper's claim.

Option 2's warrant was that a rate uses the precision on θ̂ that a sign discards. 3.0.5
shows it does not use it: the rate is better estimated **without** θ̂ at the operating
point, and the map's own exponent dominates anyway. The reframing argument against the
threshold survives intact — the evidence is a smooth dose-response with no threshold — but
that argument does not require the paper to adopt the rate as its estimand.

What the paper can defensibly claim is narrower than either option promised:

1. **The identification result**, which stands on its own: the rank-one identity and its
   closure condition, the dose train supplying 42× the marker band, and sampling on dosing
   days costing a factor of 22 while half a day off recovers it.
2. **That identification is not sufficient for the decision**, demonstrated rather than
   asserted — with the 24× threshold sensitivity and the sign-flip contour as the evidence,
   and `m_tol` attached to every decision number.
3. **That the threshold framing is itself questionable**, warranted by the cited smooth
   dose-response, stated as a finding rather than adopted as a replacement.

That is a smaller paper than Phase 3 envisaged and it is the one the evidence supports.
The decision is yours; nothing is written until you take it.

---

# 3.0.7 — can controls supply the population constants? **Route closed.**

## (a) Controls stay in the linear regime, by two independent arguments

**From the model.** A control has roughly ten times a PA patient's residual capacity
(`tab_params`: the patient's θ is "~10% of normal"). The load required to carry a control
to the neighbourhood of K_m:

| target m | m/K_m | load required | × a PA patient's baseline |
|---|---|---|---|
| 10 | 0.10 | 2 777 | 8.6 |
| 50 | 0.50 | 10 250 | 31.8 |
| **100** | **1.00** | **15 500** | **48.0** |

And where a control actually sits at ordinary loads: **m/K_m = 0.011** at the patient's own
baseline, 0.022 at twice it, 0.118 at ten times. Deep in the linear regime, where
g(m) ≈ m/K_m and only the lumped ratio θ/K_m is identifiable — Proposition 2(1).

**From the published protocol, which is decisive on its own.** The dose in the only
1-¹³C-propionate series located is **0.5 mg/kg body weight, administered orally**, and the
authors state they optimised the method using *"a 20-fold reduced dose"* to *"minimize
propionate exposure in MMA subjects"*.

**A dose chosen to minimise exposure is, by construction, not a dose that drives the system
towards saturation.** It is a tracer. It is also a bolus rather than a sustained load, so it
does not create the steady state the identification is built on. The route is closed on
magnitude and on kind.

Record it as closed.

## (b) Moot, and it would not have helped

Not required, since (a) fails. The check is cheap and confirms 2D.0 from the other side —
monthly design, minimal footing:

| K_m prior | sd(ln θ) |
|---|---|
| free (what the paper does) | 0.1697 |
| a factor of 2.2 | 0.1696 |
| 10% | 0.1663 |
| 1% — what a control study might buy | **0.1597** |

**A control-derived K_m tight to 1% buys six per cent.** k_e remains the binding direction.
So even had (a) passed, the constant it would have supplied is not the one that binds.

*(A fifth row, K_m fixed exactly, returned a singular matrix under the audit script's own
conditioning test. That is the test firing on a badly scaled matrix rather than a finding,
and it is omitted rather than reported as one.)*

---

# 3.7 — structural identifiability, symbolically. **The machine disagrees with two of my own comments, and it is right.**

## What was run

Julia 1.10.4 with `StructuralIdentifiability.jl` (40 dependencies, 193 s precompile),
on `symbolic/julia/run_identifiability.jl`; independently, `sympy` on
`symbolic/sympy/observational_equivalence.py`; and `audit/e5_structural.py` to settle what
the disagreement actually was.

## The verdicts

| model | what it encodes | global |
|---|---|---|
| A1 | generic time-varying input | θ, K_m, k_e **globally** |
| A2 | constant input, θ free | θ, K_m, k_e, s₀ **globally** |
| A3 | two patients, shared (K_m, k_e), constant inputs | all **globally** |
| A4 | as A3 with a null-genotype anchor | all **globally** |
| A5 | carnitine as a known input on the alternative channel | +κ, all **globally** |
| A6 | therapy as a known input on the enzymatic channel | all **globally** |
| **B** | **both channels saturating** | **locally, not globally** |
| C | two-state carnitine depletion, c(t) unobserved | all **globally** |

## Three disagreements, reported rather than reconciled

### (a) My sympy commentary was wrong, and the machine was right

The script had asserted that a **constant** therapeutic increment leaves only the lumped
sum θ + Δ₀ identifiable. It does not. Δ is a **known** input — it is the dosing record —
so recovering the lump recovers θ. Solving coefficients-all-zero in sympy returns the
single solution `{K_m′: K_m, L_0′: L_0, k_e′: k_e, θ′: θ}`. The model is globally
identifiable with no therapy variation at all. The script's own pre-written "EXPECTED"
comments were wrong on the same point.

This is the third correction of my own reasoning this session, and the first where an
external tool caught it rather than a re-run of my own code.

### (b) The Julia encodings A2–A4 do not encode what their comments claimed

The comments said A2 "encodes a patient held at steady state" and predicted deficiency 2.
That was a category error on my part, and `audit/e5_structural.py` Q1 settles it:

> **A constant input is not a constant state.** The framework assumes generic initial
> conditions, so m(t) relaxes *towards* its steady state and visits an interval on the way.
> From m(0) = 5 the trajectory spans [5.000, 20.000]; from m(0) = 60 it spans
> [20.000, 60.000]. Time constant ≈ 0.8 h; 99% of the gap closed in **4.19 hours**.

Proposition `prop:global` then applies directly and the verdict is unsurprising. The note's
steady-state claim is about a patient *already* at steady state, which supplies one
algebraic equation in three unknowns. Q2 exhibits the fibre explicitly: nine
(K_m, k_e, θ) triples spanning K_m over 9× and k_e over 4× that all reproduce the observed
(s, m) = (600, 20) with residual exactly 0. **Deficiency 2, as the note states.** The two
results are answers to different questions, and the ODE framework cannot express the
second one.

### (c) Model B is the one genuine structural obstruction — and it confirms the paper's channel argument

B is **locally but not globally** identifiable, with identifiable functions
`ke + theta`, `Ke + Km`, `ke*theta`, `Ke*Km`, `Ke*theta + Km*ke` — all symmetric under
exchanging the two channels. When both channels saturate, capacity and alternative
clearance are interchangeable. This is the structural counterpart of the paper's channel
argument: identification rests on the two channels having *different functional forms*.
Nothing in the paper depended on this; it is a finding, not a repair.

## What the paper's obstructions therefore are

Every obstruction the paper reports survives, but its **type** must be named correctly.
Two mechanisms, and neither is structural non-identifiability:

**(1) Structural identifiability asks about the transient, and the clinic never sees it.**
The framework assumes the exact, noiseless, continuous trajectory, from which every
time-derivative of y is available. The marker's time constant is 0.8 h against a dose
half-life of 3.5 days — a 104× separation, which is what licenses treating each visit as a
point on an algebraic balance rather than a state of an ODE — so monthly sampling returns
points on the **quasi-steady locus**, where ṁ is zero to the width of the assay. The derivatives carrying the structural result are not in the record.
The fibre is what remains once the transient is removed.

**(2) On that locus the Fisher matrix is singular at the truth, and flat only to first order.**
`e5_structural.py` Q3, on a 28-day biweekly record sampled daily with four free parameters
(θ, K_m, k_e, load bias):

| | |
|---|---|
| singular values | 8.083e−4, 1.694e−5, 1.222e−8, **3.911e−15** |
| s_min/s_max | **4.84e−12** |
| flat direction (column-scaled) | (+0.1904, −0.6942, −0.0000, +0.6942) |

Noiseless misfit along that direction, in ε:

| ε | rms misfit | misfit/ε² |
|---|---|---|
| 0.010 | 4.53e−14 | 4.53e−10 |
| 0.030 | 4.03e−13 | 4.48e−10 |
| 0.100 | 4.46e−12 | 4.46e−10 |

Exactly quadratic over a decade. The direction is flat to **first order** and curved at
second. A singular Fisher matrix at a point does not contradict global injectivity of the
map — it sets the price: the combination is recovered through curvature rather than slope,
so the recoverable displacement scales as the **square root** of the noise rather than
linearly.

## The rank-one identity, derived by hand and matched to the numerics

This was the check reserved for the PI. It is done, by a route sharing no code with the
design module. On the quasi-steady locus (θ+Δ)u + k_e m = L + b with u = m/(K_m+m):

```
dG/dtheta = u,   dG/dK_m = -(theta+Delta) u(1-u)/K_m,   dG/dk_e = m,   dG/db = -1
```

Eliminating the increment with (θ+Δ)u = L − k_e m and using m(1−u) = K_m u, a combination
(a, b, c, e) vanishes for **every** Δ iff

```
c = 0                    (m = K_m u/(1-u) is not polynomial in u)
a + bL/K_m + b k_e = 0   (coefficient of u)
e + bL/K_m = 0           (constant)
```

so, taking b = −K_m,  **n = (L + k_e K_m, −K_m, 0, L)**.

| | θ | K_m | k_e | b |
|---|---|---|---|---|
| predicted, column-scaled | +0.1904 | −0.6942 | +0.0000 | +0.6942 |
| computed from the SVD | +0.1904 | −0.6942 | −0.0000 | +0.6942 |

Max absolute difference **1.79e−07**; ‖Jn‖/(‖J‖‖n‖) = **1.53e−11**. Dividing through by θL
gives c_θ = (1 + K_m k_e/L)/θ and c_Km = −1/L, which is the identity the paper states.
**Confirmed.** Note that k_e does not appear in the null direction — the degeneracy is
between capacity, affinity and the load bias only.

## Consequence for the paper

The claims are unchanged in substance and narrower in type. `prop:rank`, the factor-of-22
aliasing rule and the fibre are results about the **conditioning of a finite, noisy,
quasi-steady design** — the design a clinic can actually run — not about the structure of
the model. Stating that makes them more precise, not weaker, and it removes a claim the
paper never made but a reader could have inferred.

The standalone companion `paper/note_structural_identifiability.tex` already says exactly
this in its abstract ("the obstruction is one of conditioning, not of structure"), so no
correction to it is required. The wrong comments were in the two symbolic scripts, and both
are now corrected in place.

---

# 3.8 — is the dose train itself a design variable? **No. About 13:1 against.**

`audit/e4_dose_design.py`. Figure of merit sd(ln θ) on the minimal footing (θ and K_m free,
priors only on ln k_e and the load record). Reported alongside it: mean Δ delivered, and
the **trough** — the minimum increment between doses, which is the quantity a dosing
interval exists to protect.

### Amplitude modulation, at constant mean dose

| schedule | sd(ln θ) | mean Δ/θ | trough/θ |
|---|---|---|---|
| fixed, 6000 every 14 d (the trial schedule) | 0.1697 | 0.723 | **0.128** |
| alternating ±25% | 0.1676 | 0.724 | 0.104 |
| alternating ±50% | 0.1647 | 0.726 | 0.075 |
| alternating ±75% | 0.1615 | 0.727 | **0.045** |
| random amplitude, sd 50% | 0.1683 | 0.764 | 0.008 |

### Dose interval, at constant delivered dose rate

| schedule | sd(ln θ) | mean Δ/θ | trough/θ |
|---|---|---|---|
| every 7 d | 0.2190 | 0.720 | 0.253 |
| every 10 d | 0.1863 | 0.720 | 0.203 |
| **every 14 d (the trial schedule)** | **0.1697** | 0.723 | **0.128** |
| every 21 d | 0.1617 | 0.724 | 0.047 |
| every 28 d | 0.1587 | 0.729 | **0.016** |

### The exchange rate

Every change that improves identifiability lowers the trough, and the exchange is lopsided:

- ±75% amplitude modulation: **4.8%** better on sd(ln θ) for a **65%** cut in the trough.
- 28-day interval: **6.5%** better on sd(ln θ) for an **87%** cut in the trough.

That is roughly **13:1 against** in both cases. The mean increment is essentially constant
across every row (0.72–0.76 θ), so the loss is concentrated exactly where it matters — the
gap a patient sits in before the next dose.

### The recommendation, stated explicitly as the handoff requires

**The identifiability optimum conflicts with the therapeutic one, and the therapeutic one
wins.** No dose modulation is recommended. The sampling schedule is free to optimise — that
is where the factor-of-22 aliasing result applies and where half a day off the dosing day
recovers it — and the dose schedule is not. The direction of the shorter-interval rows is
worth naming for the same reason: a 7-day interval is *worse* for identification and
*better* for the trough, so a clinician who shortens the interval on therapeutic grounds
pays an identification cost rather than gaining one.
