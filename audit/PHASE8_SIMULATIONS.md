# Phase 8. The simulation programme

Marks, as `PHASE4_SUMMARY.md` uses them: **[proven]** mathematics, exact or symbolic.
**[measured]** computed, footing stated. **[judgement]** a call for the author.

**The paper has not been edited and no constant in `src/common/params.py` has moved.**
Everything that had to vary is assembled in the audit script from an explicit kernel and an
explicit dose train, and the committed shape is reproduced as one member of that family.

Scripts: `audit/q1_increment_shape.py` (S1), `audit/q1b_reconcile.py` (S1 against the
25 September probes), `audit/q1c_baseline_footings.py` (Gap B on two baseline footings, and
the provenance of `ke0`), `audit/q2_requirement_map.py` (S2). Independent re-check:
`audit/q9_recheck.py`. The probes themselves: `audit/probes_2026-09-25/`.

**This file was revised twice on 26 September 2026.** The first S1 report (commit `043a3ec`)
carried numbers and conclusions that the reconciliation against the probes withdrew. The
second revision, at the author's request, separates the three error corrections from a
design change the first revision had run together with them: an untreated baseline draw.
What changed, sorted by kind, is section S1.R; `REVIEW_LOG.md` R54 records it. The author has
independently reproduced the corrected numbers with the probe, the `t = 0` draw read
pre-dose.

---

## 0. The footing every number below is on

Unless a line says otherwise:

| what | value | where from |
|---|---|---|
| capacity | `theta = 3000`, the **normal-capacity** illustration | `src/common/params.py` through `i1_estimator` |
| `K_m`, `k_e`, load | `100`, `5`, `L = 322.7273` | same |
| noise | `sigma = 21.9612`, the ordinate model of `design/run.py:598-604` | same |
| priors | `b ~ N(0, 16.1364)`, `ln k_e ~ N(ln 5, 0.8)`, and where the dose scale is free, `ln kappa ~ N(0, 0.12)` | the paper's own, and `g_design_common.sd_ln_theta`'s coordinate for `kappa` |
| design | monthly, `np.arange(0, 730, 30.4)`, `n = 25`, load CV `10%` | the paper's own |
| **baseline** | **two footings.** `baseline`: the `t = 0` draw taken before the first infusion, an untreated measurement. `on therapy`: no untreated draw, the `t = 0` draw removed, `n = 24`. S1.1b gives Gap B on both; the other tables are on `baseline` unless they say otherwise | the author, 26 Sep |
| dose train | uniform 14 d from `t = 0`, `D_AMP = 6000`, steady peak `6400` | `design/run.py` |
| **dose-day draw** | a draw at a dose instant is **never read post-dose**, which would put it inside the marker's settling window. It is read pre-dose, for truth and analyst alike; on this schedule only `t = 0` is affected, and it then becomes the baseline draw | S1.R 3a |
| estimator | global scan over `K_m` in `[1, 1000]`, exact polish inside that box | `i1_estimator.varpro_exact`, extended |
| seeds | **24**, and every median below carries its min and max | the repository's `SEEDS` |

**The amplitude convention of S1.1 to S1.5 is `peak`-matched**, as the first report's was,
so that the corrections and the convention change are not confounded. **S1.1b, the Gap B
statement the author asked for, leads with `dose`**, which becomes the primary convention from
S2 on (S1.9); S1.6 gives all three.

### 0a. Where the shape constants come from, and that they are not human

Every constant that defines a shape here is a **mouse** parameter or an assumption. Brief
**R2** (the increment time course) finds **no human hepatic PCC protein or activity time
course at any time point**, so the half-life is carried as a sweep and named against R2
throughout, and so is the measurement error in S1.8. The two-compartment rates
`kdeg = 0.0075 /h`, `kq = 0.00474 /h`, `ke0 = 0.242 /h` and `Slope = 55.6` are Attarwala et
al. 2023 Table 1, which `audit/k1_sources.py` finding 2 established are mouse.

---

## S1.R Reconciliation against the probes, and what it withdrew

### The question

For the Attarwala two-compartment truth `q1` reported `theta_hat/theta = 5.64` [2.44, 8.61]
at a `2.35 sigma` residual. The probe `p3_shape.py`, on noiseless records, gives `1.37` to
`1.40` at `0.33` to `0.37 sigma`. On the exponential truths the two agreed. So either the two
scripts model different objects or one of them is wrong. `audit/q1b_reconcile.py` puts both
on identical noiseless records (the probe's own 8 load records) and fits each record with
**both** estimators, so that a modelling difference and an estimator error cannot be
confused.

### 1. The two kernels **[proven]**

| | `q1` | probe |
|---|---|---|
| protein block | `[[-(kq+kdeg), kq], [kq, -kq]]` | the same |
| rates | `kdeg = 0.0075 /h`, `kq = 0.00474 /h`, times 24, so **per day** | the same |
| protein half-lives | **1.99 d and 11.81 d** | **1.99 d and 11.81 d** |
| where the dose enters | **the effect compartment**, `Ce(0) = 1`, `ke0 = 0.242 /h`, `Slope = 55.6` | **PCC directly**, a unit impulse |
| third mode | **0.119 d**, the effect compartment | none |
| shape after an infusion | **rises** from zero over hours | **jumps** |
| mRNA PK | none beyond the effect compartment | none |

Both share the two half-lives `PHASE5_ROBUSTNESS.md` finding 3 gives. They differ in one
thing: `q1` reads Attarwala's Suppl. Table S1(B) with its effect compartment, and the probe
without it.

**Where each amplitude convention is matched.** In `q1`, `peak` is the steady-state peak of
the infinite nominal train and `mean` the steady-state cycle mean (the impulse area over the
interval), both closed form; `dose` is the **single-dose** peak of one impulse, which is a
first-cycle quantity. The probe matches `peak` and `mean` on days 600 to 614 of a 730 d
train, which is steady state to within `exp(-600 ln2/11.81)`. So both scripts match `peak`
and `mean` at steady state, and `q1`'s `dose` is the only first-cycle convention in either.

### 2. The same kernel and convention: the fitters agree **[measured]**

On identical records the two estimators agree to **5e-7** in the median, on every row of
`q1b` R2: both kernels, three schedules, amplitude fixed and freed. Reproduced here:

| kernel | draw | schedule | probe fitter | `q1` fitter |
|---|---|---|---|---|
| probe | as run | monthly | 1.37 [1.36, 1.38] at 0.36 sigma | 1.37 [1.36, 1.38] at 0.36 sigma |
| probe | as run | 14.5 d | 1.40 [1.39, 1.41] | 1.40 [1.39, 1.41] |
| probe | as run | 2-daily | 1.36 [1.35, 1.36] | 1.36 [1.35, 1.36] |
| `q1` | as run | monthly | **5.81** [4.02, 9.03] at 2.20 sigma | **5.81** [4.02, 7.33] at 2.20 sigma |

**So there is no estimator error in the fixed-amplitude cells. The 5.8 is a property of the
record, and specifically of one draw in it.**

### 3. The t = 0 draw was two things, and only one of them an error **[measured]**

The monthly record `np.arange(0, 730, 30.4)` puts exactly one draw at a dose instant in two
years, `t = 0`, because `30.4 k` is a multiple of 14 only at `k = 35`. `design.run.delta`
reads it on the **post**-dose side (`t >= t_k`), and so did `q1`. At that draw the analyst's
increment is `6000`, the probe's jump kernel is `5541`, and `q1`'s rising kernel is **zero**.

`theta_hat/theta`, fixed 3.5 d analyst, peak-matched, 8 noiseless records (`q1b` R2):

| kernel | post-dose at t = 0 | t = 0 dropped | pre-dose at t = 0 |
|---|---|---|---|
| probe (jump) | 1.37 | 1.39 | 1.31 |
| `q1` (rising) | **5.81** | **1.23** | 1.14 |

Reading those columns in order separates two changes that the first revision of this file ran
together, and which the author has asked to be kept apart.

**(a) The fix: not reading a draw at an infusion instant as post-dose.** Read post-dose, the
draw sits inside the marker's own settling window, `0.054 d` (1.3 h) by
`g_design_common.settling_days`, where the quasi-steady model under test does not describe
the marker at all. **That is wrong for any kernel, not only the rising one**; it is merely
catastrophic for the rising one, because there the two readings of the published kernel
disagree by the whole dose. Removing it moves the rising kernel from `5.81` to `1.23` and the
jump kernel from `1.37` to `1.39`. **This is an error correction.** The probe's own 14.5 d
schedule has the same kind of draw, at `t = 392` (`0.5 + 14.5 x 27`).

**(b) The design change: an untreated baseline draw.** Read pre-dose rather than dropped,
the `t = 0` draw becomes an **untreated baseline**: a measurement of the balance with no
increment in it, which pins `theta` directly. That moves the rising kernel from `1.23` to
`1.14` and, much more, the slow exponentials (7 d, peak: `1.78` dropped, `1.38` pre-dose).
**This is not an error correction. It is a different record**, and a patient already on
therapy when the record starts cannot have it.

So the first report's `1.934` (7 d) and `1.310` (5 d) were **not wrong arithmetic**. They
belong to a record with no valid baseline, carrying one draw that should not have been there;
with that draw removed and no baseline added, the same cells read `1.744` and `1.267` (S1.1b,
noisy, 24 seeds). The first report's `5.637` for the rising kernel **was** the invalid draw:
without it, on the same no-baseline footing, `1.196`.

The corrected S1 tables below are on the `baseline` footing unless they say otherwise, and
S1.1b gives Gap B on both.

### 4. The free-amplitude cell: an error in `q1` **[measured]**

With the 12% prior on the dose scale, on identical records, the probe's fitter returned
`5.74` and `q1`'s returned `0.11`. The objectives differed: the probe carries the prior on
`ln kappa`, and `q1` carried it on `kappa`, the delta-method form. **The repository's own
convention is `ln kappa`**: the dose-scale column in `g_design_common.sd_ln_theta` is
`inc g`, the derivative in `ln kappa`, and its prior is stated on that coordinate. At
`kappa = 0.01` the linear form charges 68 and the log form 1470, which is the whole of why the
first run could "switch the increment off". **Same class of error as the `k_e` prior defect
caught before the first report.** Fixed; the fitters now agree on every interior seed to
`1.2e-7` (`q9` check 9). The two seeds where they still differ, by `0.041` and `0.223`, are
the two where `q1`'s `K_m` sits on its box bound of 1000 and the probe's runs past it (1041,
1271), which is the domain difference `q9` check 4 already measures and not the prior.

**The first version of `q9` check 9 compared every seed and FAILED at 0.223.** It is now made
on interior seeds and prints the box seeds; this paragraph is the record of the failure.

### 5. A second error, found on the way **[measured]**

Re-running `q9` on the pre-dose footing, check 5 (the profile interval against `n6`'s grid)
failed at `2.0e-2` on one upper edge, `4026` against `4108.6`. With `theta` held high on a
profile, the linear inner solve can return `k_e <= 0`; the polish then started from
`ln(1e-12)`, where `exp(ln k_e)` has no gradient, and was accepted because the start's own
objective is infinite. **The profile jumped from 1.53 to 553 between two adjacent `theta`**,
and the upper interval edge came out short. Fixed by restarting a non-positive `k_e` from the
prior's centre; the two profiles now agree to four decimals out to `theta = 8000`, and check 5
passes at `7.1e-4`. The first report's upper interval edges were affected wherever this
triggered, which makes its coverage figures conservative, not optimistic.

### 6. The wrong-family cells, reconciled **[measured]**

Monthly, 8 noiseless records, single-exponential analyst with the half-life free:

| analyst | kernel | as run | pre-dose |
|---|---|---|---|
| half-life free, dose amplitude fixed (probe p3b; `q1`'s `dose`) | probe | **1.45** (probe 1.46) | 1.31 |
| same | `q1` | 0.00, both fitters | **1.13**, both fitters |
| half-life free, amplitude re-normalised to the steady **peak** (what `q1` Part 1 ran) | probe | 1.51 | 1.30 |
| same | `q1` | 0.00 | **1.11** |
| half-life and amplitude both free, flat (probe p3c) | probe | **1.01** [0.75, 1.83] (probe 0.90 on 6 seeds) | 0.30 [0.18, 0.43] |
| same | `q1` | 0.04 | 0.29 [0.17, 0.42] |

Three things follow. **The "capacity zero on all 24 seeds" of the first report was the `t = 0`
draw**, on both fitters: it goes to `1.13` pre-dose. **`q1`'s Part 1 did not run the probe's
analyst**: it re-normalised each candidate half-life to the steady peak where p3b holds the
per-dose amplitude, and the two differ by a few per cent. And **the both-free analyst is
poorly determined**: its 8-seed range as run is `0.75` to `1.83`, the probe's `0.90` sits
inside it, and the pre-dose baseline draw moves it to `0.30`. It should not be quoted as a
point.

### 7. What changed, sorted by kind

**Error corrections**, which change the first report's numbers because the first report was
wrong:

| error | what it did | fixed by |
|---|---|---|
| a draw at an infusion instant read post-dose, inside the settling window | carried the rising two-compartment cell to `5.637` at `2.35 sigma`; drove the half-life-free analyst to capacity zero on that truth; biased the two-compartment analyst on every exponential truth | not reading such a draw post-dose (S1.R 3a) |
| the dose-scale prior on `kappa` rather than `ln kappa` | let the free-amplitude fit "switch the increment off", `0.111` | the prior on `ln kappa` (S1.R 4) |
| a profile stranded near `k_e = 0` | shortened some upper interval edges | restarting `k_e` from the prior's centre (S1.R 5) |

**A design change**, which moves numbers because it is a different record, not because
anything was wrong:

| change | what it does |
|---|---|
| an untreated baseline draw at `t = 0`, taken before the first infusion | roughly halves the slow-truth bias: 7 d, dose-matched, `1.529` without it and `1.145` with it (S1.1b) |

**What the corrections withdrew**, on either footing:

| first report said | now | status |
|---|---|---|
| rising two-compartment truth `5.637`, "not in the 1.3 to 2.0 band", "should stop being grouped" | `1.196` without baseline, `1.131` with, peak-matched; between the 3.5 d and 5 d exponentials | **withdrawn**, grouped |
| "audible, `2.35 sigma` on every seed; the paper's under-a-fifth-of-sigma sentence does not describe it" | `1.01` to `1.08 sigma` against the null's `0.97` to `0.98` | **withdrawn** |
| free amplitude "switches the increment off, capacity nine times low" | identical to the fixed analyst | **withdrawn** |
| half-life free on the two-compartment truth, "capacity zero on all 24 seeds" | `1.129` (baseline footing) | **withdrawn** |
| two-compartment analyst "biased 1.18 to 2.25"; "no analyst family is safe" | `1.05` to `1.14` on every truth (baseline footing) | **withdrawn**, S1.3 |
| "the family question wants the hours after an infusion" (S1.8) | holds only on the reading that drops plasma-mRNA PK | **qualified**, S1.8 |
| 7 d `1.934`, 5 d `1.310` | **not withdrawn**: a record without a valid baseline, `1.744` and `1.267` with the invalid draw removed | footing named |
| an unrecorded dose train costs more than a wrong shape, and is audible | Q3W read as Q2W, `6.52` at `2.20 sigma` | **stands** |
| freeing the amplitude does not help | on every truth now | **stands, widened** |

---

## S1. The increment time course: analyst against truth (corrected)

### S1.0 The checks **[measured]**

`audit/q9_recheck.py`, **24 of 24 pass**, each by a route sharing no code with what it
checks. The ones that bear on this revision:

- two-compartment kernel by hand from the Laplace transform against `q1`'s eigen-decomposition,
  `< 1e-9` **[proven]**; and against `solve_ivp`, `9.3e-12`.
- the committed kernel and train against `design.run.delta`, `1.4e-16`.
- `q1`'s fit against a from-scratch least squares on the same objective and domain, objective
  gap `3.7e-15`.
- `q1`'s profile interval against `n6`'s objective on a fixed grid, `7.1e-4`.
- check 9, the free-amplitude fit against the probe's fitter, interior seeds, `1.2e-7`.
- check 10, the dose-day draw: pre-dose every kernel is zero at `t = 0`; post-dose the
  committed analyst is `design.run.delta` there, `1.5e-16`.

And the three defects the checks caught **before** the first report, still fixed: the
unpolished profile, the `1.92` cut on a doubled objective, and the polish leaving the `K_m`
box. Three were caught **after** it, and are S1.R 4, S1.R 5 and the dose-day draw.

### S1.1 The headline grid **[measured]**

Nominal 14 d train, peak-matched, pre-dose at `t = 0`, `theta = 3000`, monthly `n = 25`,
24 seeds. `theta_hat/theta` above 1 is capacity read **too high**. `cover` is the fraction of
95% profile intervals containing the truth; with 24 seeds its standard error at 95% is 4.4
points.

| truth | analyst | `theta_hat/theta` median [min, max] | RMSE | RMS resid / sigma | cover |
|---|---|---|---|---|---|
| exp 2.5 d | fixed 3.5 d | **0.926** [0.838, 1.134] | 0.101 | 1.139 [0.86, 1.50] | 100% * |
| exp 2.5 d | half-life free | 0.977 [0.900, 1.168] | 0.067 | 0.974 | 100% |
| exp 2.5 d | two-comp free | 1.140 [0.795, 1.296] | 0.156 | 0.951 | 100% |
| **exp 3.5 d (null)** | **fixed 3.5 d** | **0.986** [0.900, 1.183] | 0.070 | **0.981** [0.71, 1.32] | **100%** |
| exp 3.5 d | half-life free | 0.996 [0.806, 1.162] | 0.081 | 0.978 | 100% |
| exp 3.5 d | two-comp free | 1.072 [0.738, 1.209] | 0.133 | 0.934 | 100% |
| exp 5.0 d | fixed 3.5 d | **1.161** [1.066, 1.377] | 0.171 | 1.139 [0.80, 1.49] | **100%** |
| exp 5.0 d | half-life free | 1.029 [0.768, 1.229] | 0.103 | 0.981 | 100% |
| exp 5.0 d | two-comp free | 1.068 [0.714, 1.239] | 0.141 | 0.937 | 100% |
| exp 7.0 d | fixed 3.5 d | **1.358** [1.267, 1.634] | 0.336 | 1.471 [1.03, 1.86] | 83.3% |
| exp 7.0 d | half-life free | 1.030 [0.765, 1.286] | 0.118 | 0.980 | 100% |
| exp 7.0 d | two-comp free | 1.060 [0.706, 1.314] | 0.144 | 0.942 | 95.8% |
| Attarwala 2-comp | fixed 3.5 d | **1.131** [1.036, 1.366] | 0.155 | 1.079 [0.74, 1.52] | 100% |
| Attarwala 2-comp | half-life free | 1.129 [0.886, 1.286] | 0.143 | 1.047 | 100% |
| Attarwala 2-comp | two-comp free | 1.051 [0.730, 1.210] | 0.145 | 0.931 | 100% |

`*` 5 of 24 intervals open below, so that 100% is partly a failure to bound rather than
precision. Every other 100% is closed intervals throughout.

**Gap B stands in direction and is smaller than the review reported, on this footing.** A
truth slower than the assumption reads capacity high and a faster one low, monotonically:
`0.926` at 2.5 d, `1.161` at 5 d, `1.358` at 7 d. The review's `1.3 to 2.0` for 5 to 7 d
becomes `1.16 to 1.36` once the `t = 0` draw is a baseline rather than a post-dose instant.
S1.6 shows the size is also convention-dependent, and S1.R 3 that it depends on whether the
record has a pre-treatment draw at all.

**The two-compartment truth is not different in kind.** `1.131` sits between the 3.5 d null
and the 5 d truth, at a residual indistinguishable from the null. **By the author's rule (a),
it is grouped with the exponentials from here on**: it behaves like a single exponential of
about 4 d, which is what a single-exponential analyst with its half-life free recovers from it
(`3.61 d` [2.30, 4.89]).

### S1.1b Gap B on two baseline footings, all three conventions **[measured]**

`audit/q1c_baseline_footings.py`. Fixed 3.5 d analyst, amplitude fixed and matched to each
convention, `theta = 3000`, monthly, noisy, 24 seeds: median `theta_hat/theta` [min, max],
the RMS residual in sigma, and profile-interval coverage. The rising and jump rows are the two
readings of the published two-compartment kernel (S1.R 1).

**Therapy starts within the record, with a pre-dose baseline draw** (`n = 25`):

| truth | **dose** | peak | mean |
|---|---|---|---|
| exp 2.5 d | **0.967** [0.875, 1.187] 1.13s 100% | 0.926 [0.838, 1.134] 1.14s 100% | 0.690 [0.625, 0.836] 1.21s 17% |
| exp 3.5 d (null) | **0.986** [0.900, 1.183] 0.98s 100% | 0.986 | 0.986 |
| exp 5.0 d | **1.080** [0.993, 1.273] 1.17s 100% | 1.161 [1.066, 1.377] 1.14s 100% | 1.442 [1.319, 1.754] 1.10s 100% |
| exp 7.0 d | **1.145** [1.065, 1.342] 1.60s 100% | 1.358 [1.267, 1.634] 1.47s 83% | 2.000 [1.829, 2.525] 1.26s 38% |
| 2-comp, rising | **1.046** [0.960, 1.255] 1.09s 100% | 1.131 [1.036, 1.366] 1.08s 100% | 1.282 [1.171, 1.567] 1.05s 100% |
| 2-comp, jump | **1.211** [1.105, 1.471] 1.09s 100% | 1.292 [1.177, 1.579] 1.07s 100% | 1.309 [1.191, 1.602] 1.07s 100% |

**Already on therapy, no untreated draw** (`n = 24`):

| truth | **dose** | peak | mean |
|---|---|---|---|
| exp 2.5 d | **0.911** [0.760, 1.180] 1.06s 96% | 0.871 [0.727, 1.124] 1.07s 92% | 0.638 [0.535, 0.811] 1.12s 4% |
| exp 3.5 d (null) | **0.978** [0.822, 1.235] 0.97s 100% | 0.978 | 0.978 |
| exp 5.0 d | **1.185** [0.990, 1.488] 1.04s 100% | 1.267 [1.054, 1.602] 1.04s 100% | 1.555 [1.277, 2.017] 1.03s 100% |
| exp 7.0 d | **1.529** [1.261, 1.943] 1.09s 71% | 1.744 [1.424, 2.260] 1.08s 67% | 2.384 [1.888, 3.288] 1.05s 8% |
| 2-comp, rising | **1.113** [0.920, 1.436] 1.01s 100% | 1.196 [0.988, 1.555] 1.01s 100% | 1.345 [1.106, 1.772] 1.01s 100% |
| 2-comp, jump | **1.273** [1.052, 1.685] 1.02s 100% | 1.353 [1.115, 1.803] 1.02s 100% | 1.370 [1.128, 1.828] 1.02s 100% |

A sensitivity row for "on therapy" places the 25 monthly draws 100 dose intervals into the
train, so no draw is first-cycle: every cell is within `0.05` of the `on therapy` table
(7 d, dose: `1.489`). So removing the `t = 0` draw is a fair stand-in for a patient already on
therapy.

**The noiseless probe agrees.** `q9` check 11 rebuilds the 7 d and jump rows from the probe's
own truths and fitter, noiseless, on both footings, and matches `q1`'s route to `6e-4`: 7 d
`1.144 / 1.378 / 2.05` with baseline and `1.557 / 1.780 / 2.44` without; jump `1.226 / 1.310
/ 1.327` and `1.311 / 1.394 / 1.411`. Those are the figures the author reported from the
probe.

**What the table says.** On the primary `dose` convention, a truth half-life of 5 to 7 d
read as 3.5 d costs `1.08` to `1.15` with a baseline draw and `1.19` to `1.53` without one.
The two readings of the two-compartment kernel sit inside that band on both footings. The
baseline draw is worth most where the truth is slowest: it takes 7 d from `1.529` to `1.145`.
**Every one of these is at a median residual inside the null's own seed range, `0.71` to
`1.32`, except 7 d with a baseline**, at `1.60 sigma` dose-matched and `1.47` peak-matched: the
misfit a baseline draw exposes at 7 d is also the one that shows. **[measured]**

### S1.2 Freeing the amplitude does not help, on any truth **[measured]**

| truth | fixed 3.5 d | + free amplitude, 12% prior on `ln kappa` |
|---|---|---|
| exp 2.5 d | 0.926 | 0.917 |
| exp 3.5 d | 0.986 | 0.985 |
| exp 5.0 d | 1.161 | 1.154 |
| exp 7.0 d | 1.358 | 1.334 |
| Attarwala | 1.131 | 1.126 |

The review's report is confirmed on every truth, and the first report's exception on the
Attarwala truth is withdrawn (S1.R 4). A scale cannot reach a misfit that is in the shape.

### S1.3 Carrying the shape: the richer family is the robust one **[measured]**

**Half-life free** removes the bias inside the single-exponential family, recovering the
half-life each time (`2.41`, `3.45`, `4.84`, `6.73 d` against 2.5, 3.5, 5, 7), and costs
precision modestly: RMSE at the null `0.070` fixed against `0.081` free. On the
two-compartment truth it does **not** help, `1.129` against the fixed analyst's `1.131`,
because no single exponential is that shape.

**The two-compartment family with free rates is within `1.05` and `1.14` on every truth
tested**, the exponentials included, at a residual at or below the null's and coverage of at
least 95.8%. It pays for that robustness: its RMSE is `0.13` to `0.16` everywhere, about
twice the fixed analyst's `0.070` when the fixed analyst happens to be right. **So the
trade is a factor of about two in precision for bias that stays under 14% across a factor of
2.8 in half-life and both family members.** The first report's "no analyst family is safe"
is withdrawn; this is the replacement. **[judgement]** on whether the trade is worth it, which
is S2's question and not S1's.

### S1.4 How loud is the warning **[measured]**

The reference is the measured null, `0.981` [0.706, 1.324], not `sqrt((n-p)/n) = 0.917`,
because the paper's fit is penalised (`q9` check 8: unpenalised `0.95`, penalised `0.97`).

| truth, fixed 3.5 d analyst | RMS resid / sigma | against the null |
|---|---|---|
| exp 5.0 d | 1.139 [0.80, 1.49] | +16%, inside the null's seed range |
| Attarwala | 1.079 [0.74, 1.52] | +10%, inside it |
| exp 2.5 d | 1.139 [0.86, 1.50] | +16%, inside it |
| exp 7.0 d | 1.471 [1.03, 1.86] | +50%; the median clears the null's worst seed |

**A wrong shape is inaudible for every truth within a factor of 1.5 of the assumption, and
the two-compartment truth among them.** Only the 7 d truth is marginally audible. The paper's
`sec:structuredep` sentence, that a misspecified fit's residual "can be under a fifth of one
sigma", **needs no qualification from S1**; the first report's claim that it did is
withdrawn. The worst cell remains the 5 d truth: `1.161` at 100% coverage and a residual
inside the null's range, so nothing visible goes wrong.

### S1.5 The trial's real, irregular dose train **[measured]**

Construction as before, `audit/PHASE6_REVIEWS_SEP2026.md` section 8, mg/kg on the committed
scale through the dose-weighted mean of the 346 published doses, `0.41012` mg/kg; the size of
the dose-2 delay ASSUMED at 7 d and the composite train a construction. **[judgement]**

**(a) An irregular train does not rescue the shape problem**, analyst holding the real train:

| truth | nominal | C5 de-escalated | C1 Q3W then Q2W | composite |
|---|---|---|---|---|
| exp 7.0 d | 1.358 | 1.418 | 1.211 | 1.314 |
| Attarwala | 1.131 | 1.187 | 1.078 | 1.280 |

**(b) Getting the dose train wrong is the larger lever, and it is audible.** Analyst assumes
the nominal Q2W train, correct 3.5 d shape:

| train the truth ran | `theta_hat/theta` | RMS resid / sigma |
|---|---|---|
| C5 de-escalated | 0.855 [0.749, 1.121] | 1.710 [1.36, 2.12] |
| C1 Q3W then Q2W | **6.517** [4.176, 8.033] | 2.203 [1.96, 2.94] |
| composite | **5.109** [4.119, 5.911] | 2.584 [2.17, 3.19] |

Unchanged in substance by the reconciliation. By the author's decision (c) this enters the
specification as a record requirement, S1.9.

### S1.6 The amplitude convention **[measured]**

`theta_hat/theta`, fixed 3.5 d analyst:

| truth | dose | peak | mean |
|---|---|---|---|
| exp 2.5 d | 0.967 | 0.926 | 0.690 |
| exp 3.5 d | 0.986 | 0.986 | 0.986 |
| exp 5.0 d | 1.080 | 1.161 | 1.442 |
| exp 7.0 d | **1.145** | **1.358** | **2.000** |
| Attarwala | 1.046 | 1.131 | 1.282 |

The 7 d bias runs `1.145` to `2.000` across conventions, a factor of 1.75, so any
cross-shape number must name its convention. **`dose` is the convention the Tier 2
calibration measurement pins** (S1.9), and it gives the smallest biases in the table.

### S1.7 The disease range **[measured]**

At `theta = 150`, footing re-solved in `n4_theta_scope`'s convention A:
`L = 322.7273`, `m_bar = 54.0231`, `sigma = 22.1774`, `b_loc = 5.6323`.

| truth, fixed 3.5 d analyst | `theta_hat/theta` | resid | cover |
|---|---|---|---|
| exp 2.5 d | 0 on 24 of 24 | 1.752 | 0%, 24 open below |
| exp 3.5 d (null) | 0.214 [0, 3.37], 11 at floor | 0.973 | 54.2%, 24 open below |
| exp 5.0 d | **6.24** [5.62, 7.23] | 1.582 | **0%**, closed |
| exp 7.0 d | **6.71** [5.75, 7.80] | 2.344 | **0%**, closed |
| Attarwala | **6.09** [5.28, 7.27] | 1.471 | **0%**, closed |

The null reproduces `sec:scope`: no lower bound and a point estimate that is not usable. The
new thing is the slow truths: **the interval closes, excludes the truth on every seed, and
the estimate is six times high.** A slower increment than assumed leaves more of it
unexplained where the enzyme's own flux is a twentieth of the load, and the fit hands it to
`theta`. The residual is audible there, `1.5` to `2.3 sigma`. The booked `D = 6000` is not
admissible at this capacity on readings R2 or R3 of Section 7, so this is a
misspecification measurement at an inadmissible dose and not a design recommendation.

### S1.8 What measurement would pin the shape **[measured]**, with one qualification

Unchanged numerically, since it is computed on steady-state phases, which the dose-day draw
does not touch.

**Estimating the half-life within the single-exponential family** needs two measurements at
the two ends of the accessible window, `sd(ln t_half) = sqrt(2)/(lambda |tau_1 - tau_2|)`
**[proven]**, re-checked to 5e-11. Against the gate that separates 3.5 d from 5 d at two sd
(`<= 0.178`), reading the increment directly needs **2** measurements at any CV up to 30%; a
functional readout (`theta + Delta`) needs **4** at 10% CV and **14** at 20%.

**Telling the families apart depends on which reading of the published kernel is right.**
`q1b` R6 re-ran it on both kernels:

| readout | kernel | earliest phase 1 h | 6 h | 1 d | N at 20% CV, 1 h / 6 h |
|---|---|---|---|---|---|
| increment | effect compartment | sep 1.36 | 0.68 | 0.52 | 2 / 2 |
| increment | jump | sep 0.62 | 0.60 | 0.53 | 2 / 4 |
| functional | effect compartment | sep 0.51 | **0.15** | 0.09 | 4 / **34** |
| functional | jump | sep **0.13** | 0.13 | 0.12 | **44** / 44 |

Under the effect-compartment reading the discriminating signal is the rise in the first hours
after an infusion, and the first report's "the family question wants the hours after an
infusion" holds. **Under the jump reading it does not**: the window barely matters, the best
pair moves to about 3.3 d and the end of the cycle, and a functional readout needs about 44
measurements at 20% CV whatever the window. **The first report's S1.8 conclusion is qualified
to the effect-compartment reading.** The half-life result is not affected, since it is within
one family.

**Where `ke0` comes from, and what the rising reading leaves out.** The timing conclusion rests
on how fast the protein rises, and in the rising reading that is set by `ke0`.

- **Value and source.** `ke0 = 0.242 /h`, Attarwala et al. 2023, *Nucleic Acid Ther*
  33(2):141-147, Table 1, RSE 3.35%; transcribed in `audit/k1_sources.py` as `KE0`. **[proven
  from the source, as transcribed]**
- **Species: mouse.** The paper scales only the PK parameters to human by allometry; `ke0`,
  `Slope`, `kdeg` and `kq` are the Pcca-/- (A138T) mouse fit, unscaled
  (`k1_sources.py` finding 2, brief R2).
- **Role: not the decay of a bolus.** Suppl. Table S1(B) writes `d(Ce)/dt = ke0 (C13 - Ce)`,
  an effect compartment **equilibrating with plasma mRNA `C13`**, which has its own
  two-compartment PK. `q1`'s rising kernel starts `Ce` at 1 and lets it decay at `ke0`: that
  is the published structure **with the plasma-mRNA stage removed**, and it is why the rise is
  over in hours (`ke0` half-life 2.86 h, single-dose peak at 13 h). The probe's jump kernel
  removes the effect compartment too. **Both readings drop a stage the published model has.**

`q1c` Part B restores the stage as a sweep, since no human plasma-mRNA half-life is
established in this repository. Brief R2 quotes a human population-PK PCCA mRNA clearance of
0.266 L/h and volume 13.9 L (Gao/Zhang, ACoP 2025), a one-compartment half-life of 36.2 h,
**not verified here** and carried as one point of the sweep. Functional readout, peak-matched:

| kernel | single-dose peak at | family separation per CV, window from 1 h | 6 h | 1 d |
|---|---|---|---|---|
| jump (probe) | 0 h | 0.131 | 0.130 | 0.125 |
| rising, no mRNA stage (`q1`) | 13.1 h | 0.515 | 0.149 | 0.087 |
| rising, mRNA half-life 2 h | 17.4 h | 0.680 | 0.281 | 0.078 |
| rising, mRNA half-life 12 h | 39.6 h | 0.719 | 0.546 | 0.130 |
| rising, mRNA half-life 36.2 h (R2's human popPK) | **72.5 h** | 0.772 | **0.689** | 0.379 |
| rising, mRNA half-life 72 h | 104.5 h | 0.756 | 0.707 | 0.497 |

**With the plasma-mRNA stage restored, the protein takes days to peak, not hours, and the
families become easier to tell apart, not harder, and less dependent on the window.** At R2's
human mRNA half-life the separation from a six-hour window is `0.69` per unit CV against `0.15`
without the stage: roughly a 3-sigma call from two functional readings at 20% CV rather than 34.
**[measured]** So the timing conclusion depends on two things this repository cannot settle:
whether the protein jumps or rises, and, if it rises, on the plasma-mRNA half-life, which is a
human PK quantity brief R2 reports and this repository has not verified. The earlier
statement "the family question wants the hours after an infusion" holds only for a fast rise;
for a slow one, any draw in the first day or two is informative.

### S1.9 The author's decisions, and what S1 now hands S2

The author answered three questions on 26 September, conditionally on this reconciliation.

**(a) Group the two-compartment truth with the exponentials.** Resolved by S1.1: it lands at
`1.131`, between 3.5 d and 5 d, and is inaudible. Grouped.

**(b) The amplitude convention is the one the calibration measurement pins.** Tier 2 is a
same-day pair of breath tests across one dose, so it measures the **per-dose jump** in
capacity. For an instantaneous-rise kernel that is the delivered dose `D`. **That is `q1`'s
`dose` convention**, which holds the single-dose peak equal to `D_AMP`; for every exponential
kernel the single-dose peak **is** the jump at the infusion. It exists and needs no adding.
One qualification the author should see: under a rising reading, a same-day breath-test pair
spanning one dose measures the per-dose jump only if its second test falls after the rise. On
`q1`'s reading, with no plasma-mRNA stage, the single-dose peak is at 13 h. With the stage
restored at R2's human mRNA half-life it is at **72 h**, and a same-day pair would see a
fraction of the jump (S1.8). So "same-day" and "per-dose jump" are the same thing only for a
fast rise. **[judgement]** `dose` becomes primary from S2 on, with `peak` and `mean` as the
sensitivity column in S4.

**(c) The dose train enters the specification as a record requirement**: administered dates,
times and amounts, including delays, skipped doses and dose changes, used as the regressor,
not a new design variable. S1.5 (b) prices ignoring it: up to `6.5x` on capacity.

**The author reviewed the reconciliation, reproduced its numbers independently, asked for the
fixes and the design change to be separated, and released S2** with `dose` primary. S4 carries
both baseline footings and both kernels.

### A flag for the author, not a finding

The paper's own headline design numbers are computed on the same record, with the `t = 0` draw
read post-dose inside the settling window. `q1b` R5 prices it through
`g_design_common.sd_over_seeds`: **`0.1699`** [0.1684, 0.1709] as committed, **`0.1695`** with
that draw read pre-dose, `0.1721` with it dropped. **The normal-capacity headline does not
depend on it.** At 5% of normal it does, and in the correctly specified model: S2 Part A0
gives the review's monthly `3.07` as the figure with **no** baseline, and `2.38` with an
untreated baseline draw. No text proposed; `paper/` is untouched.

---

## S2. The disease-range requirement map

`audit/q2_requirement_map.py`. Re-checked by `q9` check 12: the Fisher map against the probe's
own `core.fisher_sd` to 2.3% on the medians, and the Monte Carlo estimator against the probe's
fitter to `1.5e-7` in `ln theta_hat` on every record.

### S2.0 Footing **[measured]**

**The model is correctly specified**: generator and analyst share the committed 3.5 d
increment, so the amplitude convention does not arise. **The `t = 0` convention**: a draw at a
dose instant is read pre-dose; no schedule here meets a dose instant except the optional
baseline draw at `t = 0`. The review's own footing otherwise: convention A (the recorded load
stays put and the marker rises), the booked `D = 6000` at every capacity (not admissible at 20%
of normal and below on Section 7's R2 and R3, kept so the map isolates identification from
dosing), load CV 10%, 24 load records per cell.

| schedule | draws | first draw |
|---|---|---|
| monthly | 30.4 d, n = 24 | 30.4 d |
| two-monthly | 60.8 d, n = 12 | 60.8 d |
| weekly | 7.3 d from 0.5 d, n = 100 | 0.5 d |
| every 2 d | odd days, n = 365 | 1 d |

**Every schedule is on two footings**: `no baseline`, as listed, and `baseline`, the same
draws plus one untreated draw at `t = 0`. That had to be done for the comparison to be fair,
and it is itself a finding: the review's monthly record carried a `t = 0` draw and its weekly
and 2-daily records did not, so its schedule comparison at 5% of normal was partly a baseline
comparison. Its footing reproduces: `0.170` at the paper's headline cell, **`3.07` at 5% of
normal monthly with no baseline** (the review's `3.08`), `0.171` weekly and `0.117` every 2 d
with `k_e` to 2% and the load measured. **With a baseline draw the monthly 5% figure is
`2.38`.**

### S2.1 The k_e precision each design needs, Fisher **[measured]**

The log-sd of the `k_e` prior at which the median Fisher bound reaches 0.20. `never`: not even
with `k_e` exactly known. `any`: the paper's 0.8 suffices.

| capacity | load | monthly | two-monthly | weekly | every 2 d |
|---|---|---|---|---|---|
| 5% | recorded | never | never | never | 1.1% (1.2%) |
| 5% | measured | never | never | 3.1% (3.2%) | 4.2% |
| 10% | recorded | never (**2.5%**) | never | 6.4% | 7.9% |
| 10% | measured | 4.7% (6.2%) | never (**2.7%**) | 8.4% | 9.2% |
| 20% | recorded | 14.8% (15.3%) | 7.8% (13.1%) | 17.6% | 18.6% |
| 20% | measured | 16.0% (16.5%) | 9.7% (14.5%) | 18.6% | 19.9% |
| 50% | either | 46% | 40% (44%) | 50% | 52 to 53% |
| 100% | either | any | any | any | any |

No baseline, with the baseline figure in brackets where it differs. **The baseline draw
matters at 10% of normal on the sparse schedules** and hardly at all on the dense ones, which
already have hundreds of draws.

### S2.2 The Fisher bound near the gate is not the estimator: the Monte Carlo **[measured]**

Every cell of the map whose Fisher median lies in `[0.10, 0.40]`, 214 of them, re-estimated
by the MAP estimator on 96 records each. **The nuisances have to be drawn, and there are two
honest ways to draw them, which part company at a wide prior** (`q2` `mc_cell`):

- **truth-drawn**: the patient's true `k_e` is a draw about a fixed prior centre, a
  **population** prior. `i1_estimator`'s convention; primary for the paper's 0.8.
- **centre-drawn**: the true `k_e` is fixed and the prior's centre is a noisy **measurement** of
  it. Primary for every tighter prior, which is what S2b's measurement supplies.

At the paper's headline cell (normal capacity, monthly, `k_e` 0.8, load measured, 96
records) the two give RMSE **`0.192`** and **`0.295`** against a Fisher bound of `0.169`,
identically by both fitters; the first reproduces `i1`'s `0.1796`. The second carries a lower
tail: when the prior's centre is drawn well above the truth, the prior hands the load to `k_e`
and `theta_hat` collapses towards zero. **That tail is the MAP, not a search failure**: the
probe's fitter finds the same collapse. RMSE is therefore reported beside a robust spread
(half the 16 to 84 percentile range) and a count of fits off by more than 4.5x.

**What the Monte Carlo says, verdict on the upper end of the RMSE's bootstrap 90% interval:**

| capacity | what reaches 0.20 in the Monte Carlo | what the Fisher map promised that the Monte Carlo does not deliver |
|---|---|---|
| **5%** | **only a dense record with `k_e` to about 2% and the load measured**: every 2 d, RMSE `0.121` to `0.125`; weekly with `k_e` exact, `0.150` to `0.173` | weekly at `k_e` 2% with the load measured, Fisher `0.167` to `0.171`, is **edge**, `0.186` to `0.198`. Every 2 d with the load recorded, Fisher `0.19` to `0.21`, is edge at `0.22` to `0.23`. Monthly never, even with `k_e` exact and the load measured (`0.25` to `0.33`) |
| **10%** | every 2 d with `k_e` to 5%, either load; weekly with `k_e` to 5% and the load measured (`0.14` to `0.15`); monthly with `k_e` to 2% and a baseline draw (`0.17`) | weekly with the load recorded at `k_e` 5% is edge (`0.19` to `0.20`); monthly at `k_e` 2% without a baseline is edge (`0.18` to `0.21`); two-monthly is edge at best, with a baseline, and fails without |
| **20%** | every schedule with `k_e` to 10%; two-monthly edge without a baseline | `k_e` to 20% is edge on every 2 d and fails elsewhere, `0.19` to `0.28` |
| **50%** | `k_e` to 20% on monthly and two-monthly | the paper's 0.8 fails, with collapses on several cells |
| **100%** | every tighter prior on every schedule | the paper's 0.8 is **edge** everywhere, truth-drawn `0.18` to `0.22` against Fisher `0.15` to `0.19`; with no baseline draw the dense schedules show an occasional collapse even truth-drawn (1 in 96 on weekly and every 2 d) |

**The Fisher bound is optimistic by up to about 30% at 5 and 10% of normal** on the cells that
meet or nearly meet, and close to honest from 20% up, where the ratios straddle 1. The median bias of `ln theta_hat` is negative at low
capacity, `-0.02` to `-0.08`: the estimator reads capacity slightly low.

### S2.3 k_e known only as a population value **[measured]**

If `k_e` is carried as a population value, a patient's deviation from it is the prior: a
between-patient CV of 10, 20 and 40% is a per-patient log-sd of `0.100`, `0.198` and `0.385`,
and uncertainty in the population mean adds in quadrature (`0.141` for 10% CV with the mean to
10%). Brief **R1** puts the between-patient CV of the nearest comparable disposal rate at about a
third in healthy adults and gives reasons to expect more in patients (renal OAT1 secretion;
eGFR under 60 in half of adults with PA); it has no 3-HP figure. Monte Carlo, truth-drawn, load
measured:

| capacity | between-patient CV that still reaches 0.20 |
|---|---|
| 5% | **none**, not even 10% (Fisher `0.43` to `0.54`) |
| 10% | **none reliably**: 10% CV is edge on every 2 d (`0.207` to `0.209`) and fails elsewhere |
| 20% | **10%** on every schedule but two-monthly without a baseline (edge); 20% is edge on weekly and every 2 d and edge or failing on monthly |
| 50% | **40%** on every schedule but two-monthly without a baseline (edge) |
| 100% | **40%** on every schedule |

**So below about 20% of normal a population value for `k_e` does not identify capacity at any
CV brief R1 would call plausible, and `k_e` has to be measured in the patient.** The per-patient
measurement has to reach about 2% at 5% of normal and about 5% at 10% (S2.1, S2.2). That is
S2b's question.

### S2.4 What S2 hands on

1. **The requirement at 5% of normal is `k_e` measured in the patient to about 2%, the
   substrate load measured, and a dense record**, every 2 days, over two years. Weekly draws are
   at the edge with `k_e` exact. **[measured]**
2. **At 10% the requirement relaxes to `k_e` to about 5% and weekly draws**, either load; at 20%
   to `k_e` to 10% and monthly draws. **[measured]**
3. **A population prior on `k_e` is not enough below 20% of normal.** **[measured]**
4. **An untreated baseline draw is worth having on sparse schedules** at 10% of normal: it moves
   the monthly requirement from never to `k_e` 2.5% with the load recorded. **[measured]**
5. **The Fisher bound should not be quoted near the gate below 20% without the Monte Carlo
   beside it**, and the nuisance convention has to be named, because at the paper's own 0.8 the
   two honest ones differ by half. **[measured]**

---

## S2b. The `k_e` measurement itself

`audit/q2b_ke_measurement.py`. S2 priced the whole requirement on a per-patient measurement
of `k_e`: about 2% at 5% of normal, about 5% at 10%, about 10% at 20%. This section asks
whether the study `paper/sections/08_specification.tex` proposes, a micro-dose intravenous
bolus of ¹³C₃-3-HP with serial plasma sampling and timed urine collection, can deliver those
numbers. Re-checks: `q9_recheck.py` check 13, three routes, all passing; 30 of 30 overall.

The `dose` amplitude convention, the committed 3.5 d increment and the pre-dose dose-day draw
throughout, as in S2.

### S2b.0 What is not known, and is therefore swept **[judgement]**

Brief **R1** finds **no human 3-HP half-life, clearance, volume or renal fraction**, in any
population. So the tracer's kinetics are swept rather than assumed: terminal half-life 5, 15,
30, 60, 120 and 360 min, bracketing blood propionate's 4 to 5 min in controls and its 3 to 7
fold lengthening in PA. The design is swept alongside it: 4 to 16 draws log-spaced from 2 min
to the end of a 2 h, 6 h or 24 h window, assay CV 5, 10 or 15%, and the bolus expressed as
C0/LOQ, the multiple of the assay's limit the first concentration sits at, which is what a
micro-dose buys.

### S2b.1 The plasma arm can reach 2%, but only when the window matches the half-life **[measured]**

`sd(ln CL)`, one compartment, parameters (ln CL, ln V), draws below the LOQ lost. The sd is
scale-free in the dose; what the dose buys is how long the tracer stays readable.

| truth t½ | window | N = 4 | N = 8 | N = 16 |
|---|---|---|---|---|
| 5 min | 2 h | 0.045 | 0.023 | 0.016 |
| 30 min | 2 h | 0.025 | 0.018 | 0.013 |
| 30 min | 24 h | 0.152 | 0.025 | 0.017 |
| 120 min | 2 h | 0.075 | 0.066 | 0.053 |
| 360 min | 2 h | 0.254 | 0.224 | 0.178 |
| 360 min | 24 h | 0.025 | 0.020 | 0.015 |

(C0/LOQ 30, assay CV 5%; the full sweep is in the script's output.)

1. **The 2% S2 asks for at 5% of normal is reachable**, at 5% assay CV, with 8 or more draws,
   *provided the window covers the half-life*. **[measured]**
2. **A window shorter than the terminal half-life is the failure mode, and more draws do not
   repair it**: at t½ 360 min in a 2 h clinic window, 16 draws still give 0.18. Each extra
   draw buys only the assay error's √N, while the truncation bias on the slope does not move.
   **[measured]**
3. **Since the half-life is unknown by a factor of 70, the window cannot be chosen in advance.**
   A 24 h collection covers every half-life in the sweep at 8 draws or more, and is the only
   design that does. **[judgement]**
4. The Fisher figure is the fitted one: 400 simulated studies per cell reproduce it to within
   3%, with bias under 0.003 in ln CL. **[measured]**

### S2b.2 A two-compartment tracer is the real risk **[measured]**

With a two-compartment truth analysed two-compartment (V₂ = V₁, Q = CL), the same designs give:

| terminal t½ | window | N = 6 | N = 12 | N = 16 | one-comp, N = 12 |
|---|---|---|---|---|---|
| 5 min | 2 h | 0.293 | 0.202 | 0.186 | 0.038 |
| 30 min | 2 h | 0.042 | 0.029 | 0.025 | 0.030 |
| 120 min | 2 h | 0.825 | 0.479 | 0.442 | 0.118 |
| 360 min | 2 h | 55.6 | 23.6 | 20.4 | 0.395 |
| 360 min | 24 h | 0.051 | 0.034 | 0.029 | 0.034 |

1. **A truncated window does not merely degrade a two-compartment estimate, it destroys it.**
   At t½ 360 min in 2 h the sd is 20 and over: the distribution and elimination phases are not
   separated, and CL is unidentified. **[measured]**
2. **A very short half-life is the other end of the same problem**: at t½ 5 min a 2 min first
   draw is already past the distribution phase, and sd stays near 0.2 however many draws follow.
   **[measured]**
3. **Both failures are cured by matching the window, not by adding draws**: 24 h and N ≥ 8
   returns 0.03 to 0.04 at every half-life in the sweep, within a third of the one-compartment
   figure. **[measured]**
4. **So the study has to be designed as if the tracer were two-compartment**, since nothing
   establishes that it is not, and a one-compartment design that happens to be wrong is not
   conservative: it is 0.44 where it believed it was 0.12. **[judgement]**

### S2b.3 The urine arm measures the renal arm, and cannot close the gap **[measured]**

Renal clearance is the tracer recovered over the plasma AUC, so its error is the collection's
own CV and the plasma arm's in quadrature. With the plasma arm at 0.02, a 5% collection gives
0.054 and a 10% collection 0.102: **the collection, not the assay, sets the renal figure.**
And **R1 does not establish that `k_e` is entirely renal**, so the urine arm is a lower bound
on `k_e`, not a measurement of it. Its role is the renal fraction, which is a different
quantity from the one S2 needs. **[measured]**, **[judgement]** on the second point.

### S2b.4 The two readings price the observation-model assumption **[measured]**, **[judgement]** on which is primary

At trace level the paper's own balance removes the tracer by both routes, at
`CL = k_e + θ/(K_m + m₀)`. Two readings follow:

- **(ii), end product. PRIMARY.** The enzyme acts on propionyl-CoA, not on the marker, so 3-HP
  is not a PCC substrate, the tracer measures `k_e` itself, and the study's sd becomes the prior
  on `ln k_e`. This is the reading S2's pricing assumes and the one to carry. **[judgement]**
- **(i), model-literal.** The model *writes* the enzymatic term as a clearance of the marker, so
  taken literally the tracer measures `k_e + θ/(K_m + m₀)`: one joint observation on θ and
  `k_e`, entering as a row rather than a prior. The population prior 0.8 is kept, since it
  belongs to the literature and not to the study.

Reading (i) is not an arbitrary alternative. It is the observation-model assumption that the
paper's own discussion already names as its weakest, written out and priced. So the gap between
the two readings is not a choice between a premise and a sensitivity: **it is the first
quantification of what that assumption is worth.**

The enzyme's share of total clearance at the untreated level is 16% at 5% of normal, 29% at 10%
and 47% at 20%, so under (i) the study is mostly measuring `k_e` anyway, but it never separates
the two. `sd(ln θ̂)`, Fisher, median over 24 records, no baseline draw, load recorded, study 2%,
no drift:

| capacity | reading (ii) | reading (i) | S2 with no study (0.8) |
|---|---|---|---|
| 5%, monthly | 0.405 | 2.33 | 3.07 |
| 5%, every 2 d | 0.214 | 1.29 | 1.42 |
| 10%, monthly | 0.217 | 1.35 | 1.58 |
| 20%, monthly | 0.127 | 0.78 | 0.81 |

**The assumption is worth a factor of 5.8 on `sd(ln θ̂)`**: reading (ii) buys 7.6 against no
study, reading (i) buys 1.3, and under (i) no schedule at any capacity reaches the 0.20 gate.
**[measured]** That is the price of the weakest assumption in the observation model, and it is
the largest single number S1 or S2 has produced.

**And the same study can settle the reading, rather than asserting it.** A labelled bolus goes
somewhere: measuring **the fraction of label recovered unchanged against the fraction appearing
in downstream products** distinguishes a marker that is only excreted from one that re-enters
the pathway. The specification should say so, which turns a stated premise into a measured one
at no extra study. **[judgement]**

### S2b.5 Within-patient drift in `k_e` is an unknown, and it binds harder than the study

**Drift is a sweep, not a finding.** Nothing measures the within-patient stability of this
clearance over months to a year: **R1** reports no such study, and the 0, 5 and 10% range here
is chosen because renal decline is the obvious mechanism and **R1** reports eGFR under 60 in
about half of adults with PA. **[judgement]** on the range. What follows from a given drift is
arithmetic. **[measured]**

Reading (ii), Fisher, load measured, no baseline draw, study and drift in quadrature:

| capacity | schedule | study 2%, drift 0 | 2%, 5% | 2%, 10% | 5%, 10% | 10%, 10% |
|---|---|---|---|---|---|---|
| 5% | every 2 d | 0.117 | 0.249 | 0.439 | 0.475 | 0.579 |
| 5% | weekly | 0.171 | 0.280 | 0.471 | 0.510 | 0.629 |
| 10% | every 2 d | 0.061 | 0.124 | 0.219 | 0.238 | 0.293 |
| 10% | monthly | 0.177 | 0.208 | 0.280 | 0.297 | 0.351 |

**A 10% drift alone puts 5% of normal beyond the gate on every schedule, however good the
study.** At 5% drift only the densest record survives. **Once drift reaches 10%, tightening
the study from 5% to 2% buys about 8%** (0.475 to 0.439 at 5% of normal, every 2 days): the two
enter in quadrature, so the larger term governs, and below the drift the study's own precision
stops mattering. **[measured]**

### S2b.6 The Monte Carlo, at the gate **[measured]**

218 Fisher cells fall in the 0.10 to 0.40 band and were re-run through S2's MAP estimator, 96
records each, with the prior's centre drawn about the truth because the study is a measurement
(S2's `centre` convention). At 5% of normal, **exactly one cell MEETS the gate**: an untreated
baseline draw, the load measured, the study at 2%, no drift, every 2 days (RMSE 0.125
[0.109, 0.142]). Weekly at the same settings is `edge` (0.186). Every 5% cell with any drift
fails. At 10% the meeting cells need the study at 2 to 5% and drift at 0 to 5%; at 20% the
gate is met monthly across the study and drift sweep except the 10%/10% corner, which is
`edge`. Collapses appear only in the failing cells. **[measured]**

### S2b.7 What S2b hands back to S2, and to the specification

1. **The study can deliver 2 to 5% on `ln k_e`, but only with a 24 h window and at least 8
   draws**, designed for a two-compartment tracer. A clinic-day window is a design failure
   whose size depends on a half-life nobody has measured. **[measured]**
2. **The first thing to measure is the half-life itself**, in a handful of patients, because
   it is what the main study's window must be chosen against. **[judgement]**
3. **Reading (i) versus (ii) decides whether the study is worth running at all.** It is a
   biochemical question, not a design one. **[measured]** on the sizes, **[judgement]** on the
   premise.
4. **Within-patient drift in `k_e` binds harder than the study below 10% of normal**, on the
   swept range. The design question is not yet "how often to repeat the study": it is **whether
   the drift is predictable**. If it tracks kidney function, which is measured routinely in this
   disease and declines in it, a covariate replaces a repeat, and the covariate is already in
   every clinic record. Only if it does not track anything measured does repetition become the
   answer. Neither is in the specification. **[judgement]**
5. **Drift must be carried at normal capacity too**, not only in the disease range: the same
   `k_e` prior carries the headline precision, so S4 sweeps it there. **[judgement]**
6. **S2's requirement at 5% of normal survives only in its most favourable corner** and only if
   `k_e` is stable to better than 5% over two years. S2's "about 2%" was the right number; what
   S2b adds is that it is not sufficient on its own. **[measured]**

---

## S3. The marker CV against the load record's accuracy

`audit/q3_marker_cv_load.py`. `tab_D_loadprec` prices the load record's accuracy along one axis,
at the committed 5% marker CV, and the paper concludes from it that the design reaches 0.20 "only
if intake is recorded to better than about 10%". That sentence is conditional on a marker CV
nobody has measured in this population: brief **R3** puts it between an assay's repeatability and
an outpatient's biological variation, and the repository's own 12% question (`REVIEW_LOG` R4) is
unresolved. So the requirement is a surface and S3 measures it, on `tab_D_loadprec`'s own load
grid (1, 5, 10, 20, 40%) against marker CV 5 to 30%.

**Conventions, stated once.** The `dose` amplitude convention. Dose-day draws read pre-dose, and
every schedule is **on therapy**: the first draw is one dose interval in, so no draw is an
untreated baseline. The one exception is the reproduction row in S3.1, which uses the published
table's own record starting at t = 0. Capacities are normal, θ = 3000, and the working band's
edges, 2472 and 18814; the disease range is S2's and is not repeated. `k_e` carries the paper's
own 0.8 population prior throughout, so S3 measures the marker and load axes and not S2's.

Re-checks: `q9` check 14, three routes (a numerical Hessian of the log posterior against the
analytic information matrix, agreeing to 2e-8; σ against the probe's independently coded footing
formula, to 4e-14; the visits-to-gate crossings confirmed by the Hessian route). 33 of 33 overall.

### S3.1 The published table, reproduced, and what the marker CV does to σ **[measured]**

| load recorded to | σ | S3, published record | `tab_D_loadprec` | S3, on therapy |
|---|---|---|---|---|
| 1% | 15.2 | 0.1591 | 0.1590 | 0.1603 |
| 5% | 22.0 | 0.1695 | 0.1699 | 0.1719 |
| 10% | 35.5 | 0.1973 | 0.1991 | 0.2034 |
| 20% | 66.2 | 0.2815 | 0.2876 | 0.2970 |
| 40% | 129.9 | 0.4904 | 0.5049 | 0.5248 |

**The reproduction is to 0.9% up to the 10% row and drifts to 2.9% at the 40% row**, where the
published route's Kalman filter and this information matrix part company on a record whose noise
is four times the signal. That is a difference in route, not a correction, and it does not move
any conclusion below. **[measured]**

**Dropping the untreated t = 0 draw costs about 3% at the 10% row** (0.1973 to 0.2034), which is
S2b's baseline-draw effect at normal capacity, where it is small. **[measured]**

σ has two terms, `hypot(load_frac × L, b × CV × m̄)`. At normal capacity, load recorded to 10%,
the load term is 32.3 and the assay term runs 14.9 at CV 5% to 89.4 at CV 30%: **the marker CV
overtakes the load record between CV 10% and 15%.** **[measured]**

### S3.2 The surface: visits to 0.20 **[measured]**

The cell is the smallest number of draws whose median over 24 intake records meets the gate;
`never` means not by 240 draws, which is 20 years monthly. Normal capacity, monthly:

| load to | CV 5% | CV 10% | CV 15% | CV 20% | CV 25% | CV 30% |
|---|---|---|---|---|---|---|
| 1% | 6 | 18 | 33 | 46 | 86 | never |
| 5% | 11 | 24 | 35 | 50 | 98 | never |
| 10% | 29 | 33 | 42 | 64 | never | never |
| 20% | 68 | 91 | never | never | never | never |
| 40% | never | never | never | never | never | never |

Two-monthly, normal capacity, is the same surface in fewer visits and more calendar time: 6, 16,
23, 55, never, never on the 1% row, and the 10% row reaches 42 visits, which is 7 years.

1. **The paper's "better than about 10%" is right only at a 5% marker CV, and only just.** At
   CV 5% the 10% load row needs 29 monthly draws, so 2.4 years, against the booked 24. At CV 15%
   the same row needs 42 draws, 3.5 years, and at CV 25% it never arrives. **[measured]**
2. **A 20% load record is not recoverable by sampling longer** at any marker CV above 10%, and a
   40% record never is, on any schedule or capacity in the sweep. **[measured]**
3. **The band edges cost about a factor of two in visits and the upper edge costs more than the
   lower.** At CV 5%, load 5%: 29 draws at θ = 2472, 11 at normal, 30 at 18814; at CV 10% the
   upper edge needs 94 against normal's 24. **[measured]**
4. **Two-monthly sampling is not a saving.** It reaches the gate in about the same number of
   draws and takes twice as long, and it loses the CV 20% and 25% cells that monthly keeps. The
   draw count is what buys precision here, not the calendar. **[measured]**

### S3.3 Which axis binds, and where they trade places **[measured]**

**The two-year frontier**, the loosest setting on each axis that still meets the gate inside the
booked 24 monthly draws, against the loosest that meets it at all by 240:

| marker CV | loosest load, ≤ 24 draws | loosest load, ≤ 240 draws |
|---|---|---|
| 5% | 5% | 20% |
| 10% | 5% | 20% |
| 15% | — | 10% |
| 20% | — | 10% |
| 25% | — | 5% |
| 30% | — | — |

(Normal capacity. At both band edges the ≤ 24 column is 1% at CV 5% and empty above it.)

1. **Inside two years the load record has to be better than 5%, not 10%**, at any marker CV, and
   at the band edges better than 1%. The paper's 10% is a 240-draw figure, not a 24-draw one.
   **[measured]**
2. **At CV 30% nothing on the load axis rescues the design**, at any capacity, on any schedule,
   by 240 draws. **If the marker CV is really near the unresolved 12% figure the design survives;
   if it is 30% it does not, and no amount of dietary bookkeeping changes that.** **[measured]**
3. **Neither axis dominates: they cross, and the crossing is a simple rule.** σ's two terms are
   equal at `CV = load_frac × L / (b × m̄)`, which comes out at **CV ≈ the load fraction**, to
   within 10% across all three capacities (load 1% ↔ CV 1.1%, 5% ↔ 5.5%, 10% ↔ 11.0%, 20% ↔
   21.9%, 40% ↔ 43.8% at the lower edge, and slightly tighter at the upper). **So whichever of
   the two numbers is larger, in percent, is the one to spend money on.** That is a design rule
   the specification does not currently state. **[measured]**
4. Across their own sweeps at a fixed 24-draw record, the load axis multiplies sd by 1.7 to 3.3
   at normal capacity and the marker axis by 1.2 to 2.3, so **the load axis has the wider range,
   but only because its sweep runs to 40%**: at equal percentages they are interchangeable, which
   is what point 3 says. **[measured]**

### S3.4 The working band moves with both axes **[measured]**

Section 6.2's band, θ ∈ [2472, 18814], is quoted at one marker CV and one load accuracy. Re-solved
on the booked monthly two-year record, on therapy, by scanning θ on a 120-point log grid (so each
edge is located to about 6%):

| load to | CV 5% | CV 10% | CV 15% and above |
|---|---|---|---|
| 1% | [2485, 25676] | [2810, 10214] | — |
| 5% | [2642, 16699] | [2988, 7988] | — |
| 10% | [3177, 7512] | — | — |
| 20% and above | — | — | — |

1. **The band is not a property of the system, it is a property of the measurement.** At CV 5%
   and load 1% it spans a factor of 10.3; at CV 10% and load 5% a factor of 2.7; at CV 15% it is
   empty at every load accuracy in the sweep. **[measured]**
2. **The upper edge moves far more than the lower.** From CV 5% to CV 10% at load 1% the lower
   edge moves 13% and the upper edge falls by 60%. The band closes from the top. **[measured]**
3. **S3's own lower edge at CV 5% and load 1%, 2485, reproduces the quoted 2472 to 0.5%**, which
   is the consistency check that matters here. S3's upper edge, 25676, sits above the quoted
   18814 because Section 6.2's band is solved at a different load accuracy and on a record that
   includes the t = 0 draw. **[measured]**

### S3.5 What S3 hands on

1. **"Record intake to better than about 10%" needs its condition attached.** It holds at a 5%
   marker CV over about 2.4 years of monthly draws. Inside the booked two years the requirement
   is better than 5%, and at the band edges better than 1%. **[measured]**
2. **The marker CV deserves the same status in the specification as the load record**, since the
   two are interchangeable percent for percent, and it is currently the less constrained of the
   two. Brief **R3** should be asked for the within-subject biological variation of the marker in
   stable outpatients, not only the assay's repeatability. **[measured]**, **[judgement]** on the
   brief.
3. **A single design rule covers both axes:** spend on whichever of the marker CV and the load
   recording error is the larger percentage, until they are equal. **[measured]**
4. **Above a 25% marker CV the design does not close at normal capacity**, and above 15% the
   working band is empty. If the 12% question resolves upward, that is the finding that matters
   most. **[measured]**
5. **S4 inherits the marker CV as its noise axis** and should carry 10 to 15% as the brief says,
   which S3 places just past the crossover with a 10% load record. **[judgement]**

---

## What is open after S1

- **The shape itself**, and **which reading of Attarwala's kernel is right**. Nothing here
  decides either. Brief **R2** finds no human hepatic time course. S1.8 says the second
  question matters for how the first would be measured.
- **The measurement CV on an increment reading**, swept, and pinned only by **R2**.
- **The size of the dose-2 delay** in the C5 train, not published; held at 7 d.
- **Whether to carry the two-compartment family as the analyst.** S1.3 says it buys bias
  under 14% on every truth for about twice the RMSE. That is S2's to price against the
  gate.
- **Whether the marker is itself cleared enzymatically**, the observation-model assumption the
  discussion calls weakest. S2b.4 prices it at a factor of 5.8 on `sd(ln θ̂)`. The label-fate
  measurement in S2b.4 would settle it inside the study already proposed.
- **Within-patient stability of marker clearance**, and **whether its drift has a covariate**,
  kidney function the obvious candidate. Both are in `refs/brief_questions.md` for **R1**.
- **The human 3-HP half-life**, unknown by a factor of 70, which is what the tracer study's
  sampling window must be chosen against (S2b.1, S2b.2).
- **Within-patient stability of `k_e` over two years.** S2b.5 makes it the binding constraint
  below 10% of normal, and nothing in the specification addresses it.
- **The marker's within-subject biological variation in stable outpatients**, which S3 makes as
  important as the load record and which brief **R3** has not been asked for in that form. Above
  25% the design does not close at normal capacity; above 15% the working band is empty.
- **S4 and S5** have not been run.

### A boundary on this phase

S2b crossed from checking the paper's robustness into designing a new study, which is where
scope creep starts. **S3, S4 and S5 are finished as the brief writes them**, with the two
additions S2b earned (the label-fate measurement, and drift carried at normal capacity in S4),
and then this phase is written up and stops. **There is no S2c.** A new sub-question that
appears on the way is recorded here, or in `refs/brief_questions.md` if it is a question for the
literature, and is not run.
