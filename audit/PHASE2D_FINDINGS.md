# Phase 2D — what the computation found before the prose was written

Four things changed the shape of 2D. All are reproducible from `audit/`.

## 1. The excitation functional had the wrong number of factors

Section 6 built everything on `(dC/dlnK_m)/(dC/dlntheta) = -phi(m)`. That ratio is
computed from the uncorrected K_m column. The correct one is

    -(1 + Delta(t)/theta) * phi(m)

so the quantity a design must move is a **product**, and there are two factors. On the
record the trial schedule already produces (`audit/d0d_generalised_E.py`):

| factor | range | as a ratio | drop |
|---|---|---|---|
| phi(m), the marker band | 0.9213–0.9675 | 1.05x | 0.046 |
| 1 + Delta/theta, the dose train | 1.186–3.133 | 2.64x | **1.948** |

**The dose train supplies 42x the excitation the marker band does, on the same record**,
and 3.9x the E ≈ 0.50 the decision was said to need. An increment reaching Δ/θ = 0.55
meets it alone. The trial has been administering the excitation and reading it as a
nuisance.

Proposition 2 itself survives but its proof does not: with a varying Δ the low-end limit
is rank three rather than two, and what it gives back is K_m (sd 0.587 → 0.209), not
theta. The high end is untouched, because g → 1 makes the theta column constant and
parallel to the load-bias column whatever Δ does.

## 2. The gate rests on ONE prior, and it is not the forbidden one

`audit/d11_priorfree.py`. Two-daily × 2 years, increment carried, load CV 10%:

| priors | sd(ln theta) |
|---|---|
| all four (theta, K_m, k_e, load) | 0.1474 |
| **without the ln K_m prior** | **0.1474** — unchanged |
| without the ln theta prior | 0.1499 |
| without the load-bias prior | 0.1483 |
| **without the ln k_e prior** | **0.6451** |
| none | 1.3153 |

A ln K_m prior of sd 3.0 — a factor of twenty — still gives 0.1474. **`rem:no_km` is not
violated: the result does not depend on knowing K_m at all.** The whole gate rests on
knowing k_e within a factor of about three (log-sd 1.08 exactly).

So the paper claims **one** prior: k_e to within a factor of ~2, which the literature
supports. theta and K_m float. Every number below is on that minimal footing.

## 3. On the minimal footing, ordinary outpatient sampling passes

`audit/d12_minimal_footing.py`, load CV 10%, priors only on ln k_e (0.8) and the load (5%):

| design | Δ dropped | Δ carried | ratio | gate |
|---|---|---|---|---|
| six-monthly (n=4) | 2.9906 | 0.4177 | 7.2 | fail |
| one dose cycle (n=7) | 4.7378 | 0.2187 | 21.7 | fail |
| **monthly × 2y (n=24)** | 2.0920 | **0.1697** | 12.3 | **PASS** |
| weekly × 2y (n=105) | 2.1196 | 0.1528 | 13.9 | PASS |
| two-daily × 2y (n=365) | 1.0295 | 0.1499 | 6.9 | PASS |

Required precision on the delivered increment: **11.7%** monthly, 14.3% weekly, 14.6%
two-daily. Close to the 15.9% reported in 2C, now on a footing the paper is entitled to.

## 4. The sharpest design rule in the paper: do not sample at the dose interval

Sampling lands at the same phase of every cycle when the interval divides or equals the
dose interval, the increment stops varying, and the therapy channel closes outright:

| sampling interval | n | range of Δ/θ | sd(ln theta) | |
|---|---|---|---|---|
| 2 d | 365 | 1.9476 | 0.1499 | PASS |
| 7 d | 105 | 1.6333 | 0.1528 | PASS |
| **14 d = the dose interval** | 53 | **0.1333** | **3.5050** | **fail** |
| **14.5 d** | 51 | 1.9861 | **0.1612** | PASS |
| **28 d = 2x the dose interval** | 27 | 0.1333 | 3.5873 | fail |
| 30.4 d (monthly) | 25 | 1.8438 | 0.1697 | PASS |
| 60 d | 13 | 1.9352 | 0.1868 | PASS |

**Half a day off the dose interval is the difference between 3.51 and 0.16, a factor of
22.** And sampling as sparse as every 60 days — thirteen draws over two years — passes.
This is what makes the specification a small study rather than a large one.

---

## Two estimator facts that must be stated, not smoothed

**The prior-free Monte Carlo does not reproduce the Fisher, and neither is simply right.**
On the monthly design the filter reports 0.1697 and a prior-free nonlinear least squares
started half a decade away reports sd(ln theta) of 5–8 (`audit/d9_loadvariation.py`). The
filter is a local object at the truth; the global likelihood has a long flat valley that a
started-away optimiser does not climb out of, because along the near-null direction the
gradient is third order. Both numbers are reported. The design is defensible on the
Fisher/MAP reading and is not defensible as a prior-free global fit, and the paper says so.

**The published decision table's 52.5% is a seed outlier.** Twelve seeds at 400 replicates
give the managed-regime errors as **6.0% ± 1.1%** (treat unneeded) and **44.5% ± 3.5%**
(fail to treat) at |mu| = 339. The table prints 6.0% and 52.5%; the first is the mean, the
second is 2.3 sd above it. Quote the spread.

## The decision asymmetry reverses between regimes

`audit/d6_regimes.py`, recording fraction rho: `s_rec = L_base + rho (s_true - L_base)`.

| regime | treat unneeded | fail to treat | median theta_hat / theta |
|---|---|---|---|
| managed (rho = 1) | 10.5% | 42.5% | 1.079 |
| partial recording (rho = 0.5) | 81.8% | 17.5% | 0.605 |
| crisis, intake pinned (rho = 0) | 78.0% | 21.8% | 0.030 |
| crisis, protein stopped (rho = -0.25) | 100.0% | 0.0% | 0.0002 |

The mechanism is visible in the median: crisis protocol drives theta_hat to zero because
recorded intake falls while the marker climbs, so the fit must attribute the marker to lost
capacity. "Fails towards declaring patients safe" holds only while the record tracks the
load.

---

## 5. The decision does not improve, and that is the note's ending

`audit/d13_decision_therapy.py`, 12 seeds × 400 replicates, prior-free (`fit_ensemble`):

| design | κ | treat unneeded | fail to treat | cohort N |
|---|---|---|---|---|
| fixed operating point, n=9 (old) | 0 | 6.9% ± 1.7 | 44.0% ± 3.2 | ~4900 |
| **monthly across dose cycles, n=25** | 0 | 19.1% ± 2.1 | 41.6% ± 2.5 | **~2200** |
| the same, increment known to 11.7% | 0.117 | 19.1% ± 1.7 | 42.3% ± 2.9 | — |

**The therapy channel identifies θ and does not make μ decidable.** Near the boundary the
new design errs on about two in five, the same as the old one. The cohort falls from ~4900
to ~2200 — a real gain, still impossible against a trial of sixteen. And the dose-precision
requirement does not bind here at all: degrading Δ from exact to 11.7% moves the sign error
by less than its seed-to-seed spread. **The 11.7% is a requirement on identifying θ, not on
deciding μ.**

This is not a contradiction with finding 3. Finding 3 is a posterior sd under the k_e
prior; this is a prior-free global fit. Both are reported, and the gap between them is the
honest shape of the result.

### The crisis reversal, on the therapy-channel design

| recording fraction ρ | treat unneeded | fail to treat | median θ̂/θ |
|---|---|---|---|
| 1 (record tracks the load) | 17.5% | 41.8% | 1.027 |
| 0.5 (partly recorded) | 99.0% | 0.0% | 0.903 |
| 0 (intake pinned at baseline) | 100% | 0.0% | 0.086 |
| −0.25 (protein stopped) | 100% | 0.0% | 0.009 |

Sharper than on the old design: the reversal is complete by ρ = 0.5. Regret-equalised
thresholds: managed +0.07θ at equal costs, crisis −0.39θ; at a cost ratio of 10:1 or worse
every regime wants every patient escalated, which is the decision-theoretic statement that
the test has stopped contributing.

### The published table's 52.5% is a seed outlier

Twelve seeds give the managed-regime errors as **6.9% ± 1.7%** and **44.0% ± 3.2%**, range
[37.2%, 49.2%]. The table printed 6.0% and 52.5%; the second is outside the range of twelve
seeds. Quote the spread.
