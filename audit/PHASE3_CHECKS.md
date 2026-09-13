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

---

# 12.1 — does the Model B result generalise? **Yes, completely — and the two-channel case was hiding a limit.**

Derived in `audit/e6_channel_exchange.py`, checked independently in
`symbolic/julia/run_channel_exchange.jl`. Nothing below is asserted; the derivation is in
the script and the tool never saw it.

## The setting, stated once

Let the observed ordinate be a sum of *n* channels sharing one functional form `f`, each
with its own amplitude `a_i`, its own shape constant `K_i`, and its own **known** exogenous
input `Δ_i(t)` riding on the amplitude:

```
C(m, t) = Σ_i ( a_i + Δ_i(t) ) · f(m; K_i)
```

The unknown parameter is the **list of pairs** `p = ((a_1,K_1), …, (a_n,K_n))`. The inputs
are indexed by **channel position** and do not move when the parameter is permuted. That
is the entire content of "an input on one channel".

## The general result

Permuting the parameter pairs and reindexing the sum gives, for **any** `f` and any `n`:

```
C(p^π) − C(p)  =  Σ_j ( Δ_{π⁻¹(j)} − Δ_j ) · f(m; K_j)
```

The amplitudes cancel identically; the difference collects entirely on the **shape**
functions. Verified symbolically for every permutation at n = 3 with `f` left as an
uninterpreted function — so this step does not depend on the channel's functional form at
all.

If the family `{f(·;K_i)}` is linearly independent — which for Michaelis–Menten with
distinct positive `K_i` follows by multiplying through by `∏(K_j+m)` and evaluating at
`m = −K_i` — then each coefficient must vanish separately.

> **Theorem (channel exchange), first form — incomplete, see below.** A permutation π of the
> parameter pairs is an observational symmetry **iff** `Δ_{π(i)} ≡ Δ_i` for every *i*. The
> residual symmetry group is the Young subgroup `∏_B S_{|B|}`, with *B* running over the
> blocks of channels sharing an input, and the fibre has **`∏_B |B|!` points**.

This is correct for **pure permutations of the pairs**, and it is what the exhaustive check
below tests. It is **not** the operative statement, and an adversarial pass caught the gap
before it reached the paper — see *The correction* two subsections down.

Checked exhaustively against the symbolic difference over all 3! permutations:

| input pattern | predicted fibre | found |
|---|---|---|
| no input anywhere | 6 | 6 |
| one input, channel 1 | 2 | 2 |
| the **same** input on channels 1 and 2 | 2 | 2 |
| **different** inputs on channels 1 and 2 | 1 | 1 |
| a different input on each | 1 | 1 |

## The correction: the input factor is about *variation*, not presence

The first form asks whether permuting the pairs leaves `C` unchanged. The operative question
is whether **some other parameter vector** gives the same `C` — and the amplitudes are free.
Since `a_i` is itself an additive constant on channel *i*, and `Δ_i` adds to it, any
**constant** part of `Δ_i` is indistinguishable from `a_i`, so a permutation can be repaired
by shifting amplitudes.

With a known **constant** input `c` on channel 1, sympy solves the swapped system exactly:

```
b_1 = a_2 − c,   b_2 = a_1 + c        residual: 0
```

**A constant known input breaks nothing.** With an input taking two distinct values, the
`d`-linear coefficient is `K_2 − K_1`, so with distinct shapes there is no solution and the
exchange is broken. Demanding one `(b_1, b_2)` at `d = 100` and `d = 400`: no solution.

> **Theorem (channel exchange), corrected.** The exchange of channels *i* and *j* is an
> observational symmetry unless **both** `Δ_i − Δ_j` is **non-constant in t** **and**
> `f(·;K_i) ≢ f(·;K_j)`. The residual group is the permutations preserving the inputs **up
> to an additive constant**. Hypotheses: `K_i` pairwise distinct, `f` injective in `K`,
> `{f(·;K_i)}` linearly independent, and the pairs `(a_i,K_i)` pairwise distinct (so the
> action is free — otherwise the fibre is *smaller* than the count).

**The tool confirms the correction, and this is the sharper test.** Two saturating channels
with a known **constant** input on channel 1 return **locally**, not globally, identifiable,
with identifiable functions

```
a1 + a2,  K1 + K2,  K1·K2,  a1·a2 + 100·a2,  K1·a2 + K2·a1 + 100·K2
```

— precisely the invariants of the shifted swap `(a1,K1) ↔ (a2−100, K2)`, each verified by
substitution. With a **time-varying** input the same model returns globally identifiable.

## The mechanism, in the paper's own idiom

```
break of the exchange (i j)  requires   Δ_i − Δ_j non-constant   AND   f(·;K_i) ≢ f(·;K_j)
                                        the input factor               the shape factor
```

**Both must hold, and the input factor is about variation.** This is the third degeneracy in
the note to close on *variation* rather than on magnitude or precision: the sensitivity
ratio `−(1 + Δ/θ)·φ(m)` is a product of a therapy factor and a marker factor;
`prop:rank` closes on load **variability**, not load precision; `prop:exchange` closes on
input **variation**, not input presence. Three different degeneracies, and in each the
separator is a difference that moves.

## Two consequences the two-channel case conceals

**(1) A single selective input resolves the exchange iff n = 2.** At n ≥ 3 it leaves
`S_{n−1}` untouched: the fibre falls from `n!` to `(n−1)!`, not to 1. Exhibited at n = 3:
`[(300,100),(50,20),(8,700)]` and `[(300,100),(8,700),(50,20)]` are identical at every m and
every Δ.

**(2) Identical inputs on two channels separate neither.** A selective input must be
selective *pairwise*.

## Independent confirmation — the tool was not told the answer

`StructuralIdentifiability.jl`, on five models written with no expected results in the
comments (per the R6 rule):

| model | returned | fibre | predicted |
|---|---|---|---|
| 2 saturating channels, no input | all **locally** | 2 | 2! = 2 ✓ |
| 2 saturating, input on channel 1 | all **globally** | 1 | 1! = 1 ✓ |
| **3 saturating, input on channel 1 only** | `a1, K1` **globally**; `a2, a3, K2, K3` **locally** | 2 | 2! = 2 ✓ |
| 3 saturating, distinct inputs on 1 and 2 | all **globally** | 1 | 1 ✓ |
| 1 saturating + 1 **linear**, no input | all **globally** | 1 | 1 ✓ |

The third row is the discriminating one and it confirms the claim in a sharper sense than
asked: the tool identifies *which* parameters stay ambiguous, and they are exactly those of
the two channels that share an input. Its identifiable functions for that model are

```
a1,  K1,          a2+a3,  K2+K3,  a2·a3,  K2·K3,  K2·a3 + K3·a2
```

— channel 1 pinned by its input, and then **precisely the five S₂ invariants** of channels
2 and 3. Those five are structurally the same list Model B returned
(`ke+theta`, `Ke+Km`, `ke·theta`, `Ke·Km`, `Ke·theta+Km·ke`). Model B's invariant list was
the S₂ invariant ring all along, and it reappears here embedded in a larger model, sitting
on exactly the channels the input fails to separate.

## Scope: this is a limit on the note, not a result in its favour

The note's model escapes the exchange for a reason close to bookkeeping: its two channels are
not two members of one parameterised family. The enzymatic term carries `(θ, K_m)` and the
alternative term carries `k_e` alone, so **there is no permutation of `(θ, K_m, k_e)` that
exchanges them** and nothing to break. The tool agrees — saturating plus linear, no input,
returns globally identifiable.

**An earlier draft of this record argued that from the limit `lim_{K_e→∞} k_e K_e m/(K_e+m)
= k_e m`, i.e. that "the shapes differ". That reasoning is poor and is withdrawn.** The
limit is singular — it holds `k_e K_e` fixed while `K_e` grows, and `K_e = ∞` is not a value
of a positive parameter — so it does not put the two channels in one family, and "the shapes
differ" is not what the argument turns on.

**And the limit cuts the other way, which is the point worth keeping.** Model B contains the
note's model only as a boundary case, and at any *finite* `K_e ≠ K_m` with no selective input
the exchange is present and exact. So the note's freedom from it is a property **of the
linearisation**, not evidence that the linearisation was harmless. That is a scope limit the
note owes a reader, and it is now stated as one.

It is not fatal, and the reason is itself a result: the note's therapy is a decaying dose
train, so Δ **varies**, and a varying selective input breaks the exchange whenever the shapes
differ. A note built on a saturating alternative channel would still identify — but it would
depend on the therapy for something this note gets for free, and **a constant infusion would
not do it**.

## Verdict

**It generalises, and it is now stated as the general mechanism** (`prop:exchange`, with the
affinity/capacity result as one instance) — with three qualifications stated in the same
breath, because each is false if dropped:

1. **"A selective input resolves channel exchange" is false as a general claim.** True only
   at n = 2.
2. **"An input is present" is not the condition.** The input must *vary*; a constant one is
   absorbed by the amplitudes.
3. **Model B is a scope limit on this note, not corroboration of it.** An earlier draft had
   it the other way round, and the paper is corrected.

Four of the five claims in this record's first draft were refuted by an adversarial pass
before reaching the paper; the fifth (that the Michaelis–Menten difference is the general
formula's instance) survived. Two further findings from that pass were checked and
**rejected**: that `rem:ceiling`'s condition number is contradicted by the exact degeneracy
(it is not — local identifiability *means* a full-rank Jacobian, so a finite condition number
is the correct local description and the discrete fibre is an additional, compatible fact),
and that the note's constants sit on the unbreakable `K_e = K_m` case (they do not; `K_e =
300` against `K_m = 100` — the finding confused `K_e` with `K_c`, which is a different
constant in a different model).

---

# 14.1 — the excitation budget. **GATE ANSWERED IN THE NEGATIVE.**

`audit/e7_excitation_budget.py`. One footing throughout: θ and K_m free, log-prior 0.8 on
ln k_e, load recorded to 5%, 25 monthly draws over two years. Medians over 24 seeds with
the full range, because a single draw is not a result (R5).

## The three scenarios

| scenario | sd(ln θ̂) | vs target 0.20 |
|---|---|---|
| **(1)** load variability only, no therapy, at the note's own 10% CV | **1.77** [1.30, 2.61] | **9× worse** |
| (1) at 20% CV | 1.05 | 5.2× worse |
| (1) at an implausible 40% CV | 0.62 | 3.1× worse |
| **(2)** as (1) plus carnitine steps, 25–100% in k_e, gain **known** | **1.91 – 2.06** | **worse than (1)** |
| **(3)** the dose train | **0.170** [0.168, 0.171] | **passes** |

Scenario (3) is essentially insensitive to load CV (0.1699 at 10%, 0.1701 at 20%): once the
dose train is present it dominates the excitation budget entirely. Without it, load CV is
the only lever and is not enough at any value.

## Scenario 2 goes *backwards*, and the mechanism is checked rather than asserted

Carnitine acts on k_e, so it cannot touch the first factor of
`−(1 + Δ/θ)·φ(m)` at all. It can only move `φ(m)`, and only indirectly, by shifting the
level. Measured:

| k_e step | m* | φ moves | q_θ = θm/(K_m+m) loses |
|---|---|---|---|
| +25% | 9.60 | +0.37% | −3.66% |
| +50% | 9.23 | +0.70% | −7.05% |
| +100% | 8.58 | +1.31% | **−13.1%** |
| +200% | 7.52 | +2.30% | −23.0% |

**About ten to one against.** Raising k_e lowers the quasi-steady marker, and the θ
sensitivity falls with it — so carnitine buys a small move in the weak factor by pushing the
patient *deeper into the linear regime*, which is where `prop:excitation` says only the
lumped ratio survives. It does not merely supply weak excitation; it moves the patient the
wrong way along the axis that carries the information.

For contrast, `1 + Δ/θ` ranges **1.126 to 3.125** over a continuous dose cycle, a factor of
2.78, while φ moves only from 0.9091 to 0.9177. *(That 2.78 is the continuous-cycle
extremum and is not the same statistic as the paper's published 2.64, which is the range
over its own sampled record — 1.186–3.133, drop 1.948. The paper quotes the record-based
figure throughout and this script's number is not substituted for it.)*

## What scenario (1) would need

**(a) More variable intake:** a day-to-day CV of **229%** of mean intake. That is not a diet.

**(b) More draws, at 10% CV:**

| draws over two years | sd(ln θ̂) |
|---|---|
| 24 (monthly) | 1.89 |
| 104 (weekly) | 1.01 |
| 730 (**daily**) | **0.67** |

**The scaling is worse than 1/√n and the obvious extrapolation is wrong.** Independent
samples would put 730 draws at ≈0.33; the computed value is 0.67, twice that. Information in
the degenerate direction saturates, because consecutive days resample one load distribution
rather than adding independent directions. So daily venepuncture for two years lands at 3.3×
the target, and a 1/√n estimate of what would be required *understates* it.

## The conclusion, in its stronger form

Neither route exists. **The therapy is not a design improvement — it is the design's
condition of existence.** There is no second specification for untreated patients, and the
note should not imply one is available. Every design result here is conditional on a patient
receiving a therapy that varies over the sampling window, which today means a trial
participant.

---

# 14.3 — the cohort for the question that is not blocked

`audit/e8_cohort_unblocked.py`. The ~2,200 figure is for the **decision** question, which is
blocked for reasons unrelated to cohort size. The design specification makes a different
claim — that a monthly record identifies a patient's capacity to the stated precision — and
that claim is a *within-patient* one. Patients are needed only to cover heterogeneity and to
confirm the prediction.

## Achieved precision across a heterogeneous cohort

4,000 simulated patients, log-uniform over θ (9-fold), K_m (9-fold), k_e (the prior's own
factor of 2.2) and managed load (4-fold). **Every range is assumed, not measured.**

| percentile | sd(ln θ̂) |
|---|---|
| 5th | 0.080 |
| 25th | 0.132 |
| **50th** | **0.192** |
| 75th | 0.329 |
| 95th | 0.767 |

**Just over half — 52.5% — meet the 0.20 target**, and the median sits almost exactly on it.
The headline 0.17 is the illustrative parameter set, which sits near the middle of the
distribution rather than at its optimum.

**What drives the spread**, as log–log correlations with the achieved sd: `Δ/θ` **+0.56**
(and θ itself −0.56, the same fact), `K_m` **+0.59**, `k_e` +0.44, load −0.19.

That first one is a design statement, not a nuisance: **the dose is fixed in absolute terms
while residual capacity is not**, so a patient with more residual capacity receives less
*relative* excitation and is identified less well. The design is self-weakening in exactly
the patients who need the therapy least.

## Patients required

| purpose | patients |
|---|---|
| confirm the predicted precision to ±40% relative | 5 |
| **confirm the predicted precision to ±20% relative** | **14** |
| confirm to ±15% | 24 |
| bound the fraction meeting the target to ±20 points | 24 |
| **bound that fraction to ±10 points** | **96** |
| bound it to ±5 points | 384 |

Confirming the design works needs **~14 patients** — inside the phase 1/2 trial of sixteen.
Bounding *for how many* patients it works to ±10 points needs ~96, which is not. The
decision question needs ~2,200.

**Two orders of magnitude separate the blocked question from the unblocked one**, and the
specification should ask for the affordable one.

---

# 15.1 — reconciling scenario 1 against `prop:rank`. **The gate answer stands, and the contradiction was in the paper.**

`audit/e9_rank_vs_precision.py`.

## (a) The variation reaches the balance, not the noise — but my 40% row was unphysical

Levels solve the quasi-steady balance at the per-sample load, so the marker moves as it
physically would: at 10% CV, L spans [240, 430] and m spans [7.29, 13.70]. **The mechanism
was exercised.**

But the sweep drew the load from a **Gaussian**, and at 40% CV that returns `L = −7.19` and
`m = −0.205`. A negative load and a negative marker are not a record. Every scenario-1 number
is now on a **lognormal** load of the same CV, which moves the 40% row from 0.70 to **0.62**
and leaves 10% essentially unchanged (1.79 → 1.77). Corrected in the table above and in the
paper.

## (b) The load-bias column *was* prior-held — and removing it makes scenario 1 worse

| load CV | dose train | prior ON | prior OFF |
|---|---|---|---|
| 0.5% | no | 34.8 | 2726 |
| 2% | no | 8.57 | 336 |
| **10%** | **no** | **1.77** | **13.8** |
| 40% | no | 0.62 | 0.94 |
| 0.5% | yes | 0.1696 | 0.860 |
| 10% | yes | 0.1699 | 0.174 |

The prior was helping. **So 14.1's negative answer was the conservative one and it stands** —
prior-free, the untreated record is an order of magnitude further from usable than reported.

## (c) Rank restored, conditioning not. **This is the reconciliation.**

Singular spectrum of the column-normalised design matrix, seed 3:

| load CV | **no dose train** s_min/s_max | sd | **with dose train** s_min/s_max | sd |
|---|---|---|---|---|
| 0% | 5e−48 | ∞ | 1.8e−15 | 0.1696 |
| 0.5% | 9.7e−10 | **28.1** | 1.2e−3 | 0.1696 |
| 2% | 6.4e−8 | 6.92 | 1.2e−3 | 0.1696 |
| 10% | 8.9e−6 | **1.30** | 1.6e−3 | 0.1697 |
| 40% | 7.7e−4 | 0.38 | 1.6e−2 | 0.1690 |

Variation lifts the smallest singular value off zero exactly as `prop:rank` says — **forty-four
orders of magnitude** across the left column — and the precision is still 28 at half a per
cent and 1.3 at ten. **Restoring rank and achieving usable precision are different things,
and `prop:rank` delivers only the first.** Both statements were true; the paper ran them
together. Now stated as `rem:rankvsprec`.

One entry cuts the other way and is named in the paper: with the dose train at constant load
the matrix is *still* rank deficient yet reads 0.1696 — the load prior filling the deficient
direction, which is exactly what `rem:second_order` exists to diagnose. At any real variation
the data fill it instead, so the headline does not depend on that prior but **does** depend on
intake actually varying.

## The published 0.054 and 0.023 do not reproduce. **Withdrawn.**

| design | k_e | load prior | CV 2% | CV 10% |
|---|---|---|---|---|
| monthly (24) | prior 0.8 | off | 0.2343 | 0.1721 |
| monthly (24) | prior 0.8 | on | 0.1696 | 0.1697 |
| **two-daily (365)** | **known** | **off** | **0.0539** | **0.0227** |

The published figures match only the **two-daily design with k_e held exactly known** — a
footing `sec:whypriors` disowns and which `_filter`'s own docstring reserves "for the audit,
not for a headline". On the claimed footing they are **0.234 and 0.172**. Both were also
computed **with the dose train present**, so neither measured what load variation buys alone.

**That is the whole of the apparent contradiction**, and it was a defect in the paper rather
than in scenario 1.

---

# 15.2 — the saturation claim, derived. **My stated mechanism was wrong and is withdrawn.**

Scenario 1, 10% CV, n from 24 to 300,000:

| n | sd(ln θ̂) | local exponent | data info, weak dir | prior info, same dir |
|---|---|---|---|---|
| 24 | 1.304 | — | 1.80e−7 | 1.011 |
| 730 | 0.643 | 0.126 | 2.03e−6 | 1.055 |
| 10,000 | 0.478 | 0.125 | 2.66e−5 | 1.049 |
| 300,000 | 0.321 | 0.064 | 7.59e−4 | 1.052 |

The exponent is 0.06–0.24 and never approaches the independent-sample 0.5 across five orders
of magnitude in n.

**But the reason is not saturation.** I wrote that "consecutive days resample one load
distribution rather than adding independent directions". The data's information in the weak
direction **grows**, roughly in proportion to n — 1.8e−7 to 7.6e−4. It never saturates. What
it never does is **arrive**: the load prior contributes ~1.05 in that direction and does not
move, so even at n = 3×10⁵ the data supply three orders of magnitude less than the prior
already holds.

**The correct statement: the weak direction is prior-held at every n a study could reach, and
the precision falls only through the directions the data do reach.** Sub-√n here is a
statement about a fixed prior, not about correlated samples. The conclusion — that no
achievable n closes it — is unchanged; the mechanism is different, and the paper now says the
right one.
