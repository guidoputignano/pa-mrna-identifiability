# Phase 11, round 2. The (n, k) trade curve, and a caveat that attaches to every Fisher figure

*Proposals only. Nothing under `paper/` is changed by this file.*

Scripts: `audit/u11_trade_curve.py`, `audit/u12_prior_quadratic.py`.

---

## 1. The trade curve, which is the object a protocol chooses on

`audit/u3_joint_fisher.py` fixes the record at 25 monthly draws and asks how many breath
tests reach the gate. That reads as *"how many breath tests instead of a design"*. A
protocol chooses both, so the decision-relevant object is the curve of `(draws, breath
tests)` pairs that reach the gate: two points on it are exchangeable, and a point off it is
not a design.

**Anchor.** On `audit/t2_tier1_margin.py`'s own record convention the `k = 0` column
returns **118 monthly draws** with `sd(ln theta) = 0.2000` at `n = 118` and 0.2001 at
`n = 117`, reproducing `08_specification.tex:126` exactly. Matching that convention
mattered: `u3`'s record starts the dose train at the first dose and draws a lognormal
intake, where `t2` takes the train at steady state and draws a Gaussian one, and the two
give 141 draws against 118. **A footing difference, not a disagreement**, and it is stated
here because it is the size of the effect the whole Tier 1 row is about.

### Tier 1 footing, increment at 14.02%, 5% breath CV

| theta | k=0 | k=1 | k=2 | k=4 | k=8 | k=16 | k=32 | k=64 | k=128 |
|---|---|---|---|---|---|---|---|---|---|
| 3000 | **118** (9.8y) | **118** (9.8y) | 106 (8.8y) | 92 (7.7y) | 72 (6.0y) | 52 (4.3y) | 35 (2.9y) | 19 (1.6y) | 12 (1.0y) |
| 600 | — | — | — | — | — | — | 59 (4.9y) | 4 (0.3y) | 4 (0.3y) |
| 300 | — | — | — | — | — | — | — | 450 (37.5y) | 63 (5.2y) |

### Increment known exactly, 5% breath CV

| theta | k=0 | k=2 | k=8 | k=32 | k=64 | k=128 |
|---|---|---|---|---|---|---|
| 3000 | 11 (0.9y) | 11 | 10 | 7 | 7 | 6 (0.5y) |
| 600 | — | — | — | 95 (7.9y) | 4 (0.3y) | 4 |
| 300 | — | — | — | — | 94 (7.8y) | 6 (0.5y) |

*A dash is "not reached at 4000 monthly draws", which is 333 years.*

### What each pair of breath tests is worth, Tier 1 at nominal

| k | draws | years | draws saved | per breath test |
|---|---|---|---|---|
| 0 | 118 | 9.8 | | |
| 2 | 106 | 8.8 | 12 | 6.0 |
| 4 | 92 | 7.7 | 14 | 7.0 |
| 8 | 72 | 6.0 | 20 | 5.0 |
| 16 | 52 | 4.3 | 20 | 2.5 |
| 32 | 35 | 2.9 | 17 | 1.1 |
| 64 | 19 | 1.6 | 16 | 0.5 |

**The first sixteen breath tests are worth five to seven monthly draws each and the next
forty-eight are worth about one.** That is the shape a protocol needs, and it is invisible
in a table at fixed `n`.

### The dashes are the frontier showing through

In the disease range the curve does not decline: it is a dash at every `n` up to a
threshold in `k`, then drops abruptly. At `theta = 600`, nothing at `k = 16` and 59 draws
at `k = 32`. **No number of draws helps below that threshold**, because `sd(ln k_e)` is
held at 0.8 by its prior and the record cannot move it. Above the threshold enough breath
contrasts have tightened the `k_e`-containing combination `u = theta/k_e + K_m` to cross
`audit/PHASE11_KE_ROUTE.md`'s frontier, and then the draws do the rest.

**So the trade curve and the frontier are the same statement in two coordinates**, and
breath tests are a substitute for the external `k_e` measurement rather than an addition to
the record.

---

## 2. `k = 1` buys exactly nothing, and that is why breath tests come in pairs

Tested rather than asserted, at three capacities, three record lengths and two breath CVs:

| theta | n | CV | k=0 | k=1 | difference | k=2 |
|---|---|---|---|---|---|---|
| 3000 | 25 | 5% | 0.169636 | 0.169636 | 5.0e-16 | 0.169010 |
| 600 | 25 | 5% | 0.817143 | 0.817143 | 1.6e-14 | 0.576289 |
| 300 | 25 | 5% | 1.612199 | 1.612199 | 3.7e-13 | 0.823406 |
| 300 | 118 | 5% | 0.690126 | 0.690126 | 1.4e-14 | 0.549323 |

**Largest difference over eighteen cells: `3.673e-13`.** Not close — identical.

The unknown subject-specific recovery factor is a free parameter with no prior, so a single
reading determines it and nothing else. Only the contrast between two readings carries
information about capacity.

> **Breath tests are bought two at a time.** A protocol that books an odd number books one
> it cannot use. The `k = 2` column shows what the pair is worth: at `theta = 300` and 25
> draws it halves `sd(ln theta)`, from 1.612 to 0.823.

---

## 3. The profile-marginal caveat attaches to every Fisher figure at the 0.8 prior

### 3.1 The inference, and why it had to be tested rather than asserted

Round 1 reported a disagreement about which nuisance binds, and reported it unresolved. The
narrow-prior control pointed at the resolution: the disagreement vanishes as the `k_e`
prior narrows, which is the signature of a log-posterior that is **not quadratic** at a
width of 0.8, a factor of 2.2. If that is the cause, then it is a property of the **prior
width** and not of the ranking, and every number computed by inverting an information
matrix at that width inherits it.

**That is a larger claim than the one it came from, so it is tested three ways, at the
paper's own headline designs rather than at the one the disagreement was found on.**

### 3.2 Profile against marginal, against the prior width

| design | `sd(ln k_e)` | marginal | profile | ratio |
|---|---|---|---|---|
| **recommended: 25 monthly, nominal** | **0.80** | **0.1695** | **0.2002** | **1.182** |
| | 0.40 | 0.1093 | 0.1124 | 1.029 |
| | 0.20 | 0.0880 | 0.0884 | 1.004 |
| | 0.10 | 0.0818 | 0.0820 | 1.002 |
| **Tier 1: 118 monthly, increment 14.02%** | **0.80** | **0.2000** | **0.2270** | **1.135** |
| | 0.40 | 0.1522 | 0.1536 | 1.009 |
| | 0.20 | 0.1376 | 0.1376 | 1.000 |
| **disease range: 25 monthly, theta = 300** | 0.80 | 1.5820 | 3.3936 | 2.145 |
| | 0.20 | 0.4801 | 0.5499 | 1.145 |
| | 0.05 | 0.2337 | 0.2388 | 1.022 |

**The gap closes to nothing by a width of 0.2, at designs that have nothing to do with the
ranking.** The cause is the width.

### 3.3 Non-quadraticity measured directly, with no second method

The profile's left and right one-sigma half-widths are equal for a quadratic log-posterior.
Their ratio measures the departure **without reference to any marginal**, so it cannot be an
artefact of two implementations disagreeing.

| design | `sd(ln k_e)` | left | right | right/left |
|---|---|---|---|---|
| recommended, 25 monthly, nominal | **0.80** | **0.2679** | **0.1327** | **0.495** |
| | 0.40 | 0.1227 | 0.1022 | 0.833 |
| | 0.20 | 0.0896 | 0.0872 | 0.973 |
| | 0.10 | 0.0821 | 0.0820 | **1.000** |
| Tier 1, 118 monthly | 0.80 | 0.2787 | 0.1755 | 0.630 |
| | 0.20 | 0.1372 | 0.1380 | 1.006 |
| disease range, theta = 300 | 0.80 | **does not close** | 0.7871 | — |

At the paper's own prior the recommended design's posterior is **twice as wide downward as
upward**, and the asymmetry is gone by a width of 0.1. In the disease range the downward
half does not close at all, which **independently reproduces `11_scope`'s own sentence**
that *"the profile interval does not close below at all"*, from a different construction.

### 3.4 Three methods, and the truth is between them

At the recommended design, `theta = 3000`:

| object | value | source |
|---|---|---|
| Fisher marginal, what the paper quotes | **0.1695** | `u5`, `u12` |
| simulated sampling spread of the MAP estimator | **0.1777** | `u4`, 300 replicates, robust spread |
| profile half-width, mean of the two sides | **0.2002** | `u12` |

**The Fisher number is optimistic by 5% against a simulation and by 18% against a profile,
and these are three different objects rather than three estimates of one.** The Fisher
inverse is a curvature at the optimum; the simulation is a sampling spread; the profile is
a likelihood interval. They coincide for a quadratic posterior and this one is not quadratic
at a 0.8 prior.

### 3.5 What a reader should carry beside each headline number

| design | quoted `sd(ln theta)` | profile interval on theta | direction of the error |
|---|---|---|---|
| recommended, 25 monthly, nominal | 0.1695 | **[2295, 3426]** against 3000 | 1.58x wider down, 0.78x narrower up |
| Tier 1, 118 monthly | 0.2000 | [2270, 3576] against 3000 | 1.39x down, 0.88x up |
| disease range, theta = 300 | 1.5820 | [1, 659] against 300 | unbounded down |

**The symmetric interval errs in the unsafe direction.** It understates how far capacity may
fall *below* the point estimate and overstates how far it may rise above, and a decision
about whether capacity clears a threshold is a decision about the lower tail.

### 3.6 What this does and does not overturn

**It does not overturn any number.** Every figure stands as what it is: a Fisher standard
deviation at a stated prior width. `rem:curvature` already tells a reader that the marker
side's interval is about twice as wide as its spread, so the paper is not silent on the
distinction.

**What it adds is that the gap is a function of the `k_e` prior width, and that at the
width the paper carries it is large enough to matter at the gate.** The recommended
design's headline 0.170 becomes 0.200 under a profile of the same posterior, which is the
target exactly. That is the sentence worth adding, and `audit/PHASE11_PROPOSED_TEXT.md`
drafts it.

---

## 4. What is not settled

- **Which object the 0.20 gate is stated against** is not settled anywhere in the paper, and
  section 3 makes the choice consequential. That is a question for the author, not a defect
  to fix.
- The profile computations use Powell on the exact flux residual; the two optimiser faults
  found in round 1 (`REVIEW_LOG` R101) are fixed in both, and the bisection depth is 18
  steps on a bracket of `[1e-4, 6]`.
- The disease-range left half-width hits the search ceiling of 6.0 rather than diverging
  provably; it is reported as "does not close" for that reason.
