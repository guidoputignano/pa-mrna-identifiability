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

### (b) A positive result the paper does not state anywhere

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
