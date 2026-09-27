# Phase 11, item 2. The sampling modality: assay CV, home microsampling, and what binds

*Proposals only. Nothing under `paper/` is changed by this file.*

Scripts: `audit/u5_sigma_split.py`. Independent re-check:
`audit/u6_recheck_sampling.py`.

---

## 0. Footing

As `audit/PHASE11_KM_BYPASS.md` section 0: `theta = 3000`, `K_m = 100`, `k_e = 5` /d,
`L = 322.727` uM/d, `D = 6000` on a 14 d train with a 3.5 d half-life, load recorded to 5%,
day-to-day intake CV 10%, `k_e` log-prior 0.8, load-bias prior `0.05 L`, gate 0.20, medians
over 24 intake records with the seed range printed. The Tier 1 arm carries the increment's
scale as a fifth column at the published 14.02%, which is
`audit/t2_tier1_margin.py`'s construction.

**Anchor.** `u5` reaches the gate at **118** monthly clinic draws at that increment CV,
which is the number `t2_tier1_margin.py` owns and `08_specification.tex:126` states. The
agreement is the anchor, not a second derivation.

---

## 1. The split, confirmed, and the damping it creates

```
sigma = hypot( load_frac * L , slope * cv * m_bar ),   slope = theta g'(m_bar) + k_e
```

At the nominal footing, `m_bar = 10.000`, `slope = 29.7934`:

| term | value |
|---|---|
| load-record term, `0.05 x 322.727` | **16.1364** |
| assay term, `29.7934 x 0.05 x 10` | **14.8967** |
| `sigma` | **21.9612** |
| assay share of `sigma^2` | **46.01%** |
| crossover, assay CV per unit load error | **1.0832** |

The author reports 16.14 and 14.90. Confirmed. The crossover reproduces
`model_A.run.assay_load_crossover()` and the paper's own 1.08.

The split is remarkably stable across capacity, because the recorded load is held fixed
while the assay term's two factors move oppositely: assay share 46.01% at `theta = 3000`,
43.64% at 600, 45.22% at 300, 47.06% at 150.

**Measured rather than evaluated** (`u6` part A, 20000 simulated draws, seed 20260926,
each error source switched on alone): load term **16.0214**, assay term **15.0436**, both
together **22.1434** against a hypot of 21.9772. Ratios to the formula 0.9929, 1.0099 and
1.0083. The 1% on the assay term is the formula evaluating it at the untreated operating
point while the dose train visits levels from 3.4 to 8.9.

## 2. The sample-count factor, and how damped it is

| assay CV | sigma | analytic `(sigma/sigma_5%)^2` | monthly visits to the gate | counted factor |
|---|---|---|---|---|
| 2.0% | 17.201 | 0.614 | 6 | 0.545 |
| **5.0%** | **21.961** | **1.000** | **11** | **1.000** |
| 7.5% | 27.562 | 1.575 | 16 | 1.455 |
| **10.0%** | **33.883** | **2.380** | **23** | **2.091** |
| 12.5% | 40.587 | 3.416 | 34 | 3.091 |
| **15.0%** | **47.514** | **4.681** | **43** | **3.909** |
| 20.0% | 61.733 | 7.902 | 76 | 6.909 |
| 30.0% | 90.825 | 17.104 | 165 | 15.000 |

The author reports 2.4 at 0.10 and 4.7 at 0.15. Confirmed analytically (2.380, 4.681), and
the count in the actual design is **lower still** (2.091, 3.909).

If `sigma` scaled with the assay CV the factor would be `(cv/0.05)^2`, that is 4.00 and
9.00. **The damping is 1.68 at a 10% assay and 1.92 at 15%**, and the counted damping is
1.91 and 2.30. The load-record term is a floor the assay cannot sink below, and that is
what makes a cheaper, noisier, more frequent sample worth considering at all.

**Re-check by simulation**, and it validates the surface rather than the count. Two
attempts to locate the gate crossing by simulation failed, and both are recorded in `u6`'s
`part_b` docstring because the second failed quietly.

- Interpolating between whichever grid points bracketed the gate reported Monte Carlo noise
  as a count: at 160 replicates the simulated curve is **not monotone in `n`**, so it
  crosses the gate, comes back and crosses again.
- Fitting `sd = c n^p` across the grid looked principled and was not. With a per-point noise
  of about 8% in the spread, the fitted exponents came out **-0.29, -0.10 and -0.22**
  against the -0.5 a count factor assumes, and the implied counts 8.8, 14.0 and 64.3, for
  "factors" of 1.60 and 7.32. Those are not estimates of anything.

Resolving a 2.4x count factor by simulation needs the spread resolved to a few per cent at
each `n`, which is thousands of replicates per cell. What is affordable, and is the stronger
check, is to compare the simulated spread with the Fisher **at the same `(n, assay CV)`**:
the count factor is a deterministic consequence of that surface. 400 replicates per cell,
nonlinear fitting, nine cells:

| assay CV | n=12 | n=24 | n=48 |
|---|---|---|---|
| 5% | 0.1884 / 0.1976 (1.049) | 0.1709 / 0.1631 (0.954) | 0.1606 / 0.1546 (0.963) |
| 10% | 0.2316 / 0.2569 (1.109) | 0.1963 / 0.2051 (1.045) | 0.1744 / 0.1519 (0.871) |
| 15% | 0.2896 / 0.3307 (1.142) | 0.2324 / 0.2572 (1.107) | 0.1952 / 0.2127 (1.089) |

*(Fisher / simulated (ratio); bootstrap standard errors on the simulated column 0.011 to
0.018; 400/400 kept in every cell.)*

**Ratios span 0.871 to 1.142** over a fourfold range of `n` and a threefold range of assay
CV, so the surface the factor is read off is confirmed to about 14%.

## 3. Tier 1 in calendar time and burden

**The clinic arm:** 118 monthly draws, **9.82 years**, 118 clinic trips, 118 venepunctures.

**Home volumetric microsampling**, phase-spread at the stated minimum gap:

| assay CV | 3 per week | 2 per week | 1 per week | monthly |
|---|---|---|---|---|
| 5% | 119 draws, 0.88 yr | 120, 1.23 yr | 119, 2.34 yr | 119, 9.93 yr |
| **10%** | **285 draws, 2.11 yr** | 289, 2.95 yr | 284, 5.63 yr | 285, 23.88 yr |
| **15%** | **567 draws, 4.25 yr** | 561, 5.75 yr | 560, 11.12 yr | 567, 47.64 yr |
| 20% | 953 draws, 7.17 yr | 952, 9.74 yr | 952, 18.90 yr | 953, 80.15 yr |

Burden traded, against 118 trips and 118 venepunctures:

| assay CV | cadence | draws | years | clinic trips | draws per clinic draw | years saved |
|---|---|---|---|---|---|---|
| 10% | 3 per week | 285 | 2.11 | 0 | 2.42 | **7.71** |
| 10% | 1 per week | 284 | 5.63 | 0 | 2.41 | 4.19 |
| 15% | 3 per week | 567 | 4.25 | 0 | 4.81 | **5.57** |
| 15% | 1 per week | 560 | 11.12 | 0 | 4.75 | -1.29 |
| 20% | 3 per week | 953 | 7.17 | 0 | 8.08 | 2.65 |
| 20% | 1 per week | 952 | 18.90 | 0 | 8.07 | -9.08 |

**The trade is favourable and it is favourable because of the damping.** Home sampling
three times a week at a 15% assay, three times worse than the clinic's, finishes Tier 1 in
4.25 years against 9.82, with no clinic attendances at all. At a 20% assay it still
finishes 2.65 years sooner. One per week does not: at 15% and worse it is slower than the
clinic arm, because 2.4 to 8 times the draws at one seventh the cadence is a longer record,
not a shorter one.

The monthly home column is the control: 119 draws in 9.93 years at a 5% assay, against the
clinic's 118 in 9.82, which is the same design measured twice.

## 4. Phase placement: both effects computed, both smaller than expected

### 4.1 The aliasing lattice

A fixed cadence can only ever reach `q` dose-cycle phases, `q` being the denominator of
`step/14` in lowest terms:

| cadence | step (d) | phases reachable |
|---|---|---|
| daily | 1 | 14 |
| 3 per week | 7/3 | **6** |
| 2 per week | 7/2 | **4** |
| 1 per week | 7 | **2** |
| fortnightly | 14 | **1** |
| monthly | 30.4 | 35 |

*(Three a week is 7/3 days. Reducing `float(7/3)` through its decimal string returns a
denominator of 2.8e16, which is the float's rounding answering a different question. The
same trap `g_design_common.lattice_phases()` documents, met from the other side.)*

### 4.2 The cliff is between one phase and two, and nowhere else

52 draws, assay CV 5%, medians over 24 intake records:

| phases | cadence | n | fixed cadence | phase-spread | ratio |
|---|---|---|---|---|---|
| 6 | 3 per week | 52 | 0.1593 [0.1590,0.1597] | 0.1607 [0.1602,0.1612] | 0.99 |
| 4 | 2 per week | 52 | 0.1587 [0.1584,0.1591] | 0.1606 [0.1603,0.1610] | 0.99 |
| **2** | **1 per week** | 52 | **0.1584** [0.1581,0.1588] | 0.1606 [0.1601,0.1611] | **0.99** |
| **1** | **fortnightly** | 52 | **10.7530** [9.2099,13.0194] | 0.1606 [0.1603,0.1611] | **66.94** |
| 35 | monthly | 52 | 0.1602 [0.1599,0.1605] | 0.1607 [0.1602,0.1612] | 1.00 |

**Phase-spreading buys nothing at any cadence reaching two phases or more**, and it is
marginally worse, because enforcing a minimum gap stretches the record. The catastrophe is
confined to `q = 1`.

### 4.3 A conjecture tested and refuted

`u6` part C reads the **data** rank off the singular values at a constant recorded load: 3
at six and four phases, **2 at two phases**, 1 at one. That invited the conjecture that a
weekly home cadence leans on the unsmoothed diet record in a way a monthly clinic one does
not.

**It does not.** At load CV 0 the weekly cadence returns **0.1583**, a finite and ordinary
number: the two priors close both null directions at no measurable cost.

| cadence | phases | load CV 0 | load CV 10% | ratio |
|---|---|---|---|---|
| 3 per week | 6 | 0.1591 | 0.1593 [0.1590,0.1597] | 1.00 |
| 2 per week | 4 | 0.1585 | 0.1587 [0.1584,0.1591] | 1.00 |
| 1 per week | 2 | **0.1583** | 0.1584 [0.1581,0.1588] | 1.00 |
| fortnightly | 1 | **inf (singular)** | 10.7530 [9.2099,13.0194] | inf |
| monthly | 35 | 0.1600 | 0.1602 [0.1599,0.1605] | 1.00 |

What survives is simpler than a graded lattice rule: **only `q = 1` fails**, and intake
variability there turns an impossible design into a useless one rather than a workable one,
which is the difference between the `inf` and the 10.75.

### 4.4 The weekday constraint, priced

A 14 d dose interval is exactly two weeks, so the dose always falls on the same weekday and
four of the fourteen integer phases are never bookable: weekdays reach `[0,1,2,3,4,7,8,9,10,11]`
against home sampling's all fourteen.

| n | weekday-only | all days | ratio |
|---|---|---|---|
| 25 | 0.1720 [0.1707,0.1739] | 0.1711 [0.1695,0.1725] | 1.005 |
| 52 | 0.1608 [0.1605,0.1613] | 0.1603 [0.1598,0.1606] | 1.004 |
| 104 | 0.1550 [0.1549,0.1551] | 0.1548 [0.1547,0.1549] | 1.001 |

**Half a per cent.** The phase freedom home sampling has and clinic scheduling does not is
worth 0.1% to 0.5%. The case for home sampling is the calendar time in section 3, and it is
not the phases.

## 5. The binding nuisance does change, and the load-record prior is never it

25 monthly draws, increment at 14.02% so all four nuisances are present, each held exactly
known in turn, medians over 24 intake records:

| assay CV | assay share of `sigma^2` | none | `b` exact | `k_e` exact | `K_m` exact | `kappa` exact | binds |
|---|---|---|---|---|---|---|---|
| 1.0% | 3.3% | 0.2052 | 0.2050 | **0.1403** | 0.1574 | 0.1612 | `k_e` |
| 2.0% | 12.0% | 0.2061 | 0.2059 | **0.1416** | 0.1578 | 0.1622 | `k_e` |
| **5.0%** | **46.0%** | 0.2120 | 0.2118 | **0.1503** | 0.1598 | 0.1696 | **`k_e`** |
| **10.0%** | **77.3%** | 0.2319 | 0.2316 | 0.1774 | **0.1658** | 0.1937 | **`K_m`** |
| 15.0% | 88.5% | 0.2615 | 0.2612 | 0.2147 | **0.1736** | 0.2282 | `K_m` |
| 20.0% | 93.2% | 0.2981 | 0.2977 | 0.2580 | **0.1821** | 0.2692 | `K_m` |
| 30.0% | 96.8% | 0.3837 | 0.3833 | 0.3537 | **0.2003** | 0.3616 | `K_m` |

**Yes, it changes.** `k_e` binds while the load-record term still carries `sigma`, `K_m`
once the assay term does, and **the crossover sits between a 5% and a 10% assay CV, which
is exactly where the assay share of `sigma^2` goes from 46% to 77%.**

**And the load bias is never the binding nuisance at any assay CV.** Holding `b` exactly
known buys between 0.0002 and 0.0004 throughout. The load-record prior that
`08_specification.tex:367` reports stuck near 1.05 stays stuck, at every assay CV, and
nothing about the sampling modality moves it. It is not what a protocol should attack.

### 5.1 The re-check disagrees on the ranking, and the disagreement is informative

`u6` part D profiles the posterior over `ln theta` by Powell minimisation instead of
inverting an information matrix. At the paper's `k_e` log-prior of 0.8 it puts `K_m` first
at **every** assay CV, including the two where `u5` puts `k_e` first.

Two optimiser faults were found and fixed before this was treated as a disagreement: a
least-squares solver applied to the square root of the objective, whose Jacobian is
singular at the optimum; and a flux residual formed as a first-order conversion of a level
difference through the local gain rather than evaluated exactly. Neither changed the
verdict.

The control did. Re-running **both** methods at a `k_e` log-prior of 0.2:

| assay CV | profile: none / `k_e` ex. / `K_m` ex. | marginal: none / `k_e` ex. / `K_m` ex. |
|---|---|---|
| 2.0% | 0.1472 / 0.1428 / 0.0672 | 0.1465 / 0.1416 / 0.0689 |
| 5.0% | 0.1556 / 0.1515 / 0.0723 | 0.1549 / 0.1503 / 0.0741 |
| 10.0% | 0.1827 / 0.1793 / 0.0845 | 0.1813 / 0.1774 / 0.0874 |
| 20.0% | 0.2670 / 0.2649 / 0.1110 | 0.2607 / 0.2580 / 0.1165 |

**They agree to 3% everywhere, and both rank `K_m` first.** A profile is not a marginal and
the two coincide only for a Gaussian posterior; a log-prior of 0.8 is a factor of 2.2 and
that direction is not Gaussian. The same width is what `audit/u4_recheck_breath.py`
identifies as the source of its own 0.9-to-1.2 spread.

**So the honest reading is split.** The `K_m`-binds verdict at a high assay CV is robust to
the choice of method. The `k_e`-binds verdict at a low assay CV is **not**: it holds for the
marginal, which is the object the design's own gate is stated in and which `u4` validates
against simulated sampling spread, and it inverts under a profile at the paper's own prior
width. Anything built on the ranking at low assay CV should say which object it means.

## 6. What is not settled

- **No assay CV for home volumetric microsampling is taken from any source.** Section 3
  sweeps it precisely because nothing here can supply it. Which column applies is an
  external measurement, and it is the single input that decides whether the trade in
  section 3 is worth making.
- The 118-draw anchor inherits everything `t2_tier1_margin.py` inherits, including the
  14.02% between-visit figure that `08_specification.tex` already flags as an upper bound
  on the wrong side of the requirement.
- Section 5's split verdict is a real ambiguity, not a rounding one.
