# The binding nuisance for the decision is not the one that binds the estimate

*`audit/u28_margin_frontier.py`. The author's derivation, checked first, then the frontier re-run with
the decision margin as the target instead of `ln theta`.*

---

## 1. The derivatives, confirmed three ways

| quantity | analytic | central differences | `margin_and_gradient`'s own `grad` |
|---|---|---|---|
| `d mu / d ln theta` | **+3000.00** | +3000.00 | +3000.00 |
| `d mu / d ln K_m` | **−2388.41** | −2388.41 | −2388.41 |
| `d mu / d ln k_e` | +700.00 | +700.00 | +700.00 |

**A one per cent relative error in `K_m` moves the margin 79.6% as much as a one per cent error in
capacity itself.** The author's figures are right.

**One refinement the argument needs carried with it: the ratio is not a constant.** It depends on
`s_peak`, which sets `theta_crit`:

| amplitude | `theta_crit` | `|dmu/dlnKm| / |dmu/dlntheta|` |
|---|---|---|
| unsafe side, 3.58× baseline | 3343.8 | **0.796** |
| safe side, 2.98× baseline | 2666.0 | 0.635 |
| the default, 1.5× `clearance_A` | 994.3 | 0.237 |

The unsafe side is the costly error, so 0.796 is the relevant figure, but it is a figure about one
amplitude and should be quoted with it.

## 2. Why the gradient cannot settle it, and what does

`theta` and `K_m` enter `mu` with **opposite signs** in log coordinates and are correlated in the
Fisher, which is the cancellation `margin_and_gradient`'s docstring calls the margin theorem. At the
recommended design, published priors:

| term | contribution to `var(mu)` | share |
|---|---|---|
| diagonal `ln theta` | 258,572 | 455% |
| diagonal `ln K_m` | 27,274 | 48% |
| diagonal `ln k_e` | 313,233 | 552% |
| **all cross terms** | **−542,280** | **−955%** |

**`sd(mu)` from the diagonals alone would be 774.0. It is 238.3: a factor of 3.25 of cancellation.**
So sensitivity times uncertainty overstates by more than three times and the question has to be
answered by propagating the full covariance, which is what the frontier below does.

## 3. The frontier with `mu` as the target

Two gates, both of which **move with `theta`**, because `mu = theta − theta_crit` and `theta_crit` is
fixed, so `|mu|` grows as capacity falls:

| `theta` | `mu` | comparable gate, `0.20 theta` | decidable gate, `|mu|` |
|---|---|---|---|
| 3000 | −343.8 | 600.0 | 343.8 |
| 600 | −2743.8 | 120.0 | 2743.8 |
| 300 | −3043.8 | 60.0 | 3043.8 |

| `theta` | `sd(mu)` | `|mu|/sd` | `k_e` width reaching `0.20 theta` | `K_m` width reaching `0.20 theta` |
|---|---|---|---|---|
| 3000 | 183.4 | 1.87 | 4.408 (82×), **looser than published** | already suffices |
| **600** | 145.8 | **18.81** | **unreachable: even exact misses** | **0.0566 (1.058×)** |
| **300** | 140.8 | **21.62** | **unreachable: even exact misses** | **0.0169 (1.017×)** |

Against the estimate's frontier from `u9`, same designs: `k_e` widths of **0.8001, 0.1567, 0.0633**.

---

## 4. Two findings, and they point in opposite directions

### (a) The inversion is real, and it is the author's

**In the disease range an exactly known `k_e` cannot bring the margin to the comparable gate, and
`K_m` can, at 1.058× and 1.017×.** The estimate's frontier is the mirror image: `k_e` is what has to
be tightened, to 1.170× and 1.065×, and `K_m` is worth about one per cent.

**So the binding nuisance for the decision is `K_m` where the binding nuisance for the estimate is
`k_e`.** Both of this repository's earlier statements were right and neither settled the other,
because `mu` is a different function of the parameters than `theta` is and **a sensitivity ranking
does not transfer between them.** "K_m does not bind" has been carried as though it settled both.

**What that says about the specification**, which is the point rather than either number: Tier 2 and
the `k_e` companion measurement are sized for identifying capacity. **If the deliverable is the
decision, the measurement to buy is `K_m`** — and `K_m` is the one quantity the paper says the record
already determines, to 6–7% at every capacity, which is 1.065× in log terms and **tighter than either
requirement above.**

### (b) A positive result the paper does not state anywhere, ON A CONDITION STATED FIRST

**THE RESULT IN THE FORM IT SHOULD BE PUBLISHED IN, which is the author's:**

> **Moving either convention moves all three ratios in the same direction, so the ORDERING
> survives any single convention and the MAGNITUDES do not.**

Everything below is that sentence's working.

**THE CONDITION, BEFORE THE RESULT, BECAUSE THE RESULT IS ONLY AS GOOD AS IT.** The whole of
finding (b) is distance to a fixed boundary. The author verified that it reproduces from
`theta_crit` alone, and it does, exactly:

| `theta` | `mu = theta - 3343.8` | `sd(mu)` | `|mu|/sd` |
|---|---|---|---|
| 3000 | −343.8 | 183.4 | **1.87** |
| 600 | −2743.8 | 145.8 | **18.82** |
| 300 | −3043.8 | 140.8 | **21.62** |

**So there is no subtle mechanism here and the finding should not be dressed as one.** A patient
two thousand seven hundred units below the critical capacity does not need a precise estimate to
be classified. That is the entire content, and it is worth stating plainly rather than as a
property of the information matrix.

**Which makes `theta_crit`'s provenance the condition on the result.**
`theta_crit = (s_peak − k_e·m_tol)(K_m + m_tol)/m_tol`, and its three inputs are:

| input | what it is | status |
|---|---|---|
| `m_tol = 40` | *"tolerable metabolite level, a clinical convention"*, `src/model_A/run.py:35` | **a convention, and the paper calls it the more troubling of the two** |
| `s_peak` | the amplitude of an illness | **a convention, never measured** |
| `K_m` | the half-saturation | **the parameter that diverges by eight orders of magnitude under crisis protocol** (R157) |

**The paper already owns this and says it more sharply than a caveat would.**
`03_dynamic.tex:376-378`: the local reduction *"imports two quantities that are never measured:
`m_tol`, a clinical convention about what metabolite level is tolerable, and `s_peak`, the
amplitude of an illness. **The convention is the more troubling of the two, because it is not a
fact about the patient at all.**"* It then notes that the equilibrium reduction does not need it.

**One attribution corrected, and the paper's REASON is the part worth carrying**, not the
attribution. "More troubling" applies to **`m_tol`**, not to `s_peak`, and the paper says why:
*"because it is not a fact about the patient at all."* `s_peak` is at least a claim about how ill
someone gets. `m_tol` is a claim about what a clinic will tolerate.

**So the boundary that decidability is measured from is partly set by a quantity that describes
nobody.** That is a sharper statement than "both are conventions", and it is the one the proposed
text should carry: the distance `|mu|` is large in the disease range, and part of what it is a
distance *from* is an administrative choice.

**So the result reads:** *conditional on a boundary whose position is set by two clinical
conventions and one parameter this audit has shown to be unstable, the sign of `mu` is robust in
the disease range and fragile near nominal capacity.* **The sign is robust; the position of the
line it is measured from is not.** Moving `m_tol` or `s_peak` moves `theta_crit` and moves every
ratio in the table with it, in the same direction for all three capacities, so the ORDERING
survives any single convention and the MAGNITUDES do not.

### Tested against the limit where the fibre is exactly null

**The author's challenge, and it is the right one to put to this result.** Along the fibre
`d mu/db = 9.9500` against `d theta/db = 2.5493`, confirmed analytically and by central differences
at every capacity: **the margin is 3.903 times more exposed to the unidentified direction than
capacity is.** Yet `sd(mu)` comes out small, and small only through a 3.25x cancellation. A
quantity more exposed to the fibre emerging better determined, by a cancellation, is the shape of a
result that might depend on something incidental.

`audit/u29_fibre_limit.py` sweeps the load variation to zero, where Proposition 4 is exactly
singular, and watches both together.

| `theta` | `sd(ln theta)` 0.10 to 0 | `sd(mu)` 0.10 to 0 | `|mu|/sd` 0.10 to 0 | cancellation 0.10 to 0 |
|---|---|---|---|---|
| 3000 | 0.1543 to 0.1542, **x1.00** | 183.4 to 210.4, **x1.15** | 1.87 to 1.63 | 4.00 to 3.51 |
| 600 | 0.6988 to 0.7331, **x1.05** | 145.9 to 185.6, **x1.27** | 18.81 to **14.79** | 4.31 to 3.61 |
| 300 | 1.2062 to 1.2930, **x1.07** | 140.8 to 178.9, **x1.27** | 21.62 to **17.01** | 3.81 to 3.26 |

**THE RESULT HOLDS, DEGRADED AND BOUNDED.** `sd(mu)` is the more exposed of the two, growing by
1.15 to 1.27 against 1.00 to 1.07, so the dependence is real and must be stated. But **neither
diverges**, both converge to finite limits by `load_cv ~ 0.002`, `|mu|/sd` falls by at most **21%**
and then plateaus, and the cancellation weakens by 12 to 16% rather than by orders. **The ordering
across capacities is untouched.**

**AND THE SWEEP FOUND WHAT ACTUALLY HOLDS THE DESIGN TOGETHER, which is neither candidate.**

| `theta` | `load_cv` | load-bias prior | min/max eigenvalue |
|---|---|---|---|
| 3000 | 0.10 | present | 2.13e-07 |
| 3000 | 0.10 | **removed** | 7.08e-08 |
| 3000 | 0 | present | 1.45e-07 |
| 3000 | **0** | **removed** | **1.41e-19** |
| 300 | **0** | **removed** | **-1.20e-18** |

**At a constant load with the load-bias prior removed, the matrix is singular to machine
precision**, confirming Proposition 4 directly. The dose train varies the levels by **162%** at
nominal and **550%** at a twentieth, and **does not break the fibre**, because the fibre is a
direction in parameter space that level variation alone does not separate.

**So the dependence to state is on `LOAD_FRAC`, a registered constant, and not on `LOAD_CV`.** What
keeps the constant-load design finite is the prior on the load bias; the load variation is what
makes the design **prior-free** rather than what makes it possible. That is a different sentence
from the one the challenge anticipated, and a more useful one, because `LOAD_FRAC = 0.05` is a
number the specification already asks a study to achieve.

### The result itself

**`|mu|/sd(mu)` rises from 1.87 at nominal capacity to 18.81 at `theta = 600` and 21.62 at 300.** The
margin becomes **more** decidable as capacity falls, because `|mu|` grows while `sd(mu)` slightly
shrinks, 183.4 to 140.8.

That bears directly on Section 5's negative result, which is stated as though decidability failed
everywhere. **On this calculation it fails near normal capacity, at 1.87 standard deviations, and not
in the disease range.** Section 7's scope limit is about recovering `theta`; it does not carry over to
the sign of `mu`, for the same reason (a) does not.

---

## 5. What is conditional on what, stated so neither finding is over-read

**Finding (b) does not depend on a gate at all.** `|mu|/sd(mu)` is a ratio of computed quantities.

**Finding (a) depends on the comparable gate `0.20 theta`, which is MY construction and not
inherited.** `mu` has no gate in the literature or in this paper. On the decidable gate, `|mu|`, the
published priors already suffice at all three capacities, so no frontier exists to report. The
inversion is therefore a statement about how hard it is to pin `mu` to within a fifth of capacity,
which is the closest analogue to the estimate's own gate, and it is not a statement that the decision
currently fails.

**Both are Fisher calculations at the true parameters**, so they inherit every caveat the paper
already attaches to its Fisher figures, and `sd(mu)` here is a marginal width rather than a profile.
The three symptoms of non-Gaussianity recorded in `OPEN_QUESTIONS.md` section 0 apply to it
unchanged.

**A defect of mine, caught and fixed before any of this was read.** The first run hardcoded the
decidable gate at `|mu|` for `theta = 3000` and applied it at 600 and 300, where the true gate is
eight times looser. Both disease-range rows then read "the existing prior already suffices" against
the wrong number. **A gate that does not move with the quantity it judges is the footing defect
wearing a different hat**, and it is the eighth defect of this kind in my own checking code this
sitting.
