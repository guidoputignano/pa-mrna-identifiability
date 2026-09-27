# Phase 11, item 1. The `K_m` bypass, the breath test's actual functional, and the joint Fisher

*Proposals only. Nothing under `paper/` is changed by this file. The wording this work
suggests for `11_scope.tex` is drafted in `audit/PHASE11_PROPOSED_TEXT.md` and is
**not applied**.*

Scripts: `audit/u1_fibre_bypass.py`, `audit/u2_breath_functional.py`,
`audit/u3_joint_fisher.py`. Independent re-check: `audit/u4_recheck_breath.py`.

---

## 0. The footing every number below is on

`src/common/params.py` and `src/design/run.py` throughout:
`theta = 3000` at nominal, `K_m = 100` uM, `k_e = 5` /d, recorded load held fixed at
`L = clearance_A(10) = 322.727` uM/d, `D = 6000` on a 14 d train with a 3.5 d half-life,
so `Delta` runs 400 to 6400 at the train's steady state. Capacity sweeps hold the
**recorded load fixed** and re-solve the operating level, `sigma` and the load-bias prior at
each capacity, which is Table `tab:scope`'s own convention. `k_e` log-prior 0.8, load-bias
prior `0.05 L`, `sigma = model_A.run.sigma_total(0.05, 0.05) = 21.9612`, gate 0.20.
Marker record 25 monthly draws unless stated. Load CV 0 for the headline tables and 10%
with 24 seeds wherever a spread is quoted.

**Anchor.** `u3`'s marker block reproduces `audit/g_design_common.py`'s `sd_ln_theta` to
`0.000e+00` over nine comparisons (constant load, and 10% load CV over eight seeds). The
two are the same design, not two designs that happen to agree.

---

## 1. The fibre bypass is exactly right, and it is worth about five per cent

### 1.0 CORRECTION: the algebra is the paper's own, and this work re-derived it

**`eq:invariants` at `06_excitation.tex:31` already publishes `A`, `B` and `C` with these
definitions and with `C = -k_e`**, together with the sentence that at a constant recorded
load `(A,B,C)` is exactly what the data contain; `eq:fibrecl` already publishes the fibre,
parametrised by the load bias. `audit/PHASE4_PI.md` section 3.2 owns the derivation and
`REVIEW_LOG` R26 and R27 extend it. Section 1.1 below is therefore a **re-derivation by a
second route**, not a discovery, and `OPEN_QUESTIONS.md` section 0's "find the owner before
recomputing" was not followed. The two parametrisations are verified to be the same curve:
substituting `eq:fibrecl`'s `K_m(b)` into this work's `theta(K_m)` returns `eq:fibrecl`'s
`theta(b)` with a `sympy` residual of **`0`**.

**What is new is sections 1.2 and 1.3**, which are numerical: that the record's null space
IS that fibre to machine precision, and that closing it buys five per cent.

### 1.1 The algebra holds, symbolically

Clearing the denominator in the quasi-steady balance gives a quadratic in `m` whose
coefficients are the only things a constant-load record can determine:

```
C m^2 + (Delta - B) m - A = 0
A = (L + b) K_m      B = (L + b) - theta - k_e K_m      C = k_e
```

`sympy` confirms three identities (`u1` part A, all `True`):

- the triple reproduces the quadratic exactly (all three coefficient residuals `0`);
- inverting `(A, B, C)` for `(theta, b, k_e)` at a **given** `K_m` returns
  `theta = A/K_m - B - C K_m`, `b = A/K_m - L`, `k_e = C`, so the author's inversion is
  exact;
- the fibre's slope has the closed form
  `dln theta / dln K_m = -(L + b + k_e K_m) / theta`, confirmed by finite difference to
  `3e-10` or better at every capacity.

That slope is `0.2742 / 1.3712 / 2.7424 / 5.4849` at `theta = 3000 / 600 / 300 / 150`. It
scales as `1/theta`, so the fibre is ten times steeper at a tenth of nominal.

### 1.2 The record's null space IS the fibre, exactly

The data's information matrix (no priors), in scaled coordinates
`(ln theta, ln K_m, b/L, ln k_e)`:

| theta | n | eigenvalues | `lam_min/lam_max` | `lam_2/lam_max` | `abs(cos)` with the fibre tangent |
|---|---|---|---|---|---|
| 3000 | 25 | 4.96e-13, 6.35e-04, 1.19e+02, 1.10e+04 | 4.52e-17 | 5.80e-08 | **1.00000000** |
| 3000 | 400 | 2.62e-10, 1.20e-02, 1.62e+03, 1.94e+05 | 1.35e-15 | 6.18e-08 | **1.00000000** |
| 600 | 25 | -1.34e-13, 6.33e-02, 1.77e+02, 9.00e+03 | -1.49e-17 | 7.04e-06 | **1.00000000** |
| 300 | 25 | 1.61e-12, 1.29e-01, 2.14e+02, 8.38e+03 | 1.92e-16 | 1.53e-05 | **1.00000000** |

A relative smallest eigenvalue at machine level says the null space is exactly
one-dimensional; the cosine says which direction it is. **The author's premise is
confirmed and it is confirmed exactly.** One external `K_m` restores rank.

### 1.3 But restoring rank is not buying precision

`sd(ln theta)` at 25 monthly draws with one quantity held **exactly** known:

| held exactly known | theta=3000 | theta=600 | theta=300 |
|---|---|---|---|
| nothing (the paper's footing) | 0.1696 | 0.8171 | 1.6122 |
| `K_m` exact | 0.1606 | 0.7879 | 1.5239 |
| `k_e` exact | 0.0813 | 0.1270 | 0.2169 |
| `b` exact | 0.1691 | 0.8143 | 1.6064 |
| `K_m` and `k_e` exact | 0.0669 | 0.1268 | 0.2121 |
| `k_e` and `b` exact | 0.0801 | 0.1069 | 0.1681 |
| `K_m` and `b` exact | 0.1444 | 0.6775 | 1.2275 |
| all three exact | 0.0242 | 0.0588 | 0.1037 |

**An exact `K_m` is worth a factor of 1.056, 1.037 and 1.058.** The reason is in the table
above it: the second eigenvalue is 5.8e-8 to 1.5e-5 of the largest, and `K_m` does not
touch it. Closing the fibre by any route costs the same: `b` exact, which also closes it,
buys 0.3%.

**`k_e` is the load-bearing external quantity, not `K_m`**: an exact `k_e` is worth a
factor of 2.09, 6.43 and 7.43. That reproduces, from an independent construction, the
paper's own Table `tab:keprior` figure of "2.1 at nominal capacity and 6.9 in the disease
range".

> **The first finding is therefore a qualified yes.** One external `K_m` does close the
> fibre, exactly, as the author says. It buys five per cent, because the fibre is not what
> binds.

---

## 2. The breath test does not measure what the ratio assumed, and the correction matters

### 2.1 The functional, derived rather than assumed

Label and unlabelled substrate are indistinguishable, so tracer leaves the pool at the
pool's own fractional disposal rate and the enzymatic share of each departure is
`theta g(m) / C(m)`. Only the enzymatic route carries carbon to CO2. So for an
infinitesimal bolus at a pool held at `m`,

```
F(m) = theta g(m) / (theta g(m) + k_e m) = theta / (theta + k_e (K_m + m))
```

and the reading is `y = c F(m)`, with `c` lumping the label's onward yield and bicarbonate
recovery.

**`F` is not proportional to `theta g(m)`** (`u2` part A, verdict `False`; the ratio
`F/(theta g)` has nonzero `d/dm`). `c` still cancels in a ratio, which is the property that
matters and which survives intact. But `F` is a Moebius function of `m` with a single pole,
so a two-level ratio determines exactly one number, and it is not `K_m`:

```
r = F(m1)/F(m2) = (u + m2)/(u + m1),      u = theta/k_e + K_m
u = (m2 - r m1)/(r - 1)
```

The chain back to capacity needs the record's own `(A, B, C)`:
`K_m = A/(B + C u)` and `theta = C (u - K_m)`, and substituting the fibre's own `u(K)`
returns `theta(K)` identically (`sympy`, `True`).

The author's stated form for their own functional is confirmed exactly: `sympy`'s
inversion minus `m1 m2 (1-r)/(r m2 - m1)` is `0`.

### 2.2 Both amplifications, in closed form

```
|dln u   / dln r| = (u + m1)(u + m2) / (u (m2 - m1))
|dln K_m / dln r| = (K + m1)(K + m2) / (K (m2 - m1))
```

the same shape with `u` in place of `K_m`. Both confirmed by finite difference to the
digits printed. Since `u = theta/k_e + K_m` is 700 at nominal against `K_m = 100`, the
amplification is far larger at high capacity: 128.4 against 20.3.

### 2.3 The author's own figures, reproduced exactly, then corrected

At a 5% breath CV per reading, peak and trough of one dose cycle
(`sd(ln r) = sqrt(2) x CV`):

| theta | m_peak | m_trough | fold | sd(ln K_m), author's functional | sd(ln theta), author's chain | sd(ln u), derived | sd(ln theta), derived chain |
|---|---|---|---|---|---|---|---|
| 3000 | 3.364 | 8.911 | 2.65 | **1.435** | **0.394** | 9.081 | 27.007 |
| 600 | 4.483 | 24.801 | 5.53 | **0.454** | 0.622 | 0.869 | 4.063 |
| 300 | 4.677 | 31.229 | **6.68** | **0.366** | **1.003** | 0.524 | 3.563 |
| 150 | 4.780 | 35.642 | 7.46 | 0.326 | 1.786 | 0.393 | 4.347 |

The author reports 1.44, 0.46 and 0.37 for the first of those columns and 0.40 and 1.01 for
the second, and a dose-train fold of 6.6. **Every figure reproduces.** The propagation the
author used is exactly `|dln theta/dln K_m| = (L + b + k_e K_m)/theta`, which is 1.2 above.

**The correction.** On the functional the model actually implies, the same propagation
gives 27.0, 4.06, 3.56 and 4.35. The direction the author noted survives for the measured
combination, `sd(ln u)` falling monotonically as capacity falls (9.08, 0.87, 0.52, 0.39),
but **not for `theta`**, which is non-monotone with a minimum near `theta = 300`.

### 2.4 The finite bolus, integrated exactly

Integrating the pool and the tracer together (`u2` part C; residual tracer below 5.4e-07 at
2 d), the inferred `u` is biased **down**, and the bias does **not** cancel with `c`,
because the same absolute bolus displaces the two levels by different relative amounts:

| bolus as a fraction of the level | 0.01 | 0.05 | 0.20 | 0.50 | 1.00 | 2.00 |
|---|---|---|---|---|---|---|
| relative error in `u_hat`, theta=3000 | -0.52% | -2.56% | -9.49% | -20.76% | -34.34% | -51.05% |
| relative error in `u_hat`, theta=300 | -0.52% | -2.54% | -9.41% | -20.45% | -33.60% | -49.57% |

The error is almost independent of capacity and is a pure function of the bolus fraction,
which makes it correctable if the bolus size relative to the pool is known and a first-order
problem if it is not.

**Caveat inherited, not introduced.** This calculation administers tracer into the pool the
model calls `m`, the plasma marker. The paper states as an assumption that the marker stands
in for the enzyme's actual substrate; the breath test's bolus is propionate, on the
substrate side of that map. Every number in this section carries that assumption **more
heavily than the marker record does**, and that is a reason to treat the section as scoping
rather than as a specification.

---

## 3. The number that was asked for: the joint Fisher, and the smallest `k`

Parameter vector `(ln theta, ln K_m, b, ln k_e, ln c)`, the recovery factor carried **free
with no prior**. Breath readings are attached to visits the marker record already contains,
so they cost an assay and not a blood draw, and they alternate the dose train's trough and
peak, which `u2` part E establishes is the widest lever available.

**`k = 1` contributes exactly nothing** at every capacity and every CV (0.170 against 0.170,
0.817 against 0.817, 1.612 against 1.612). The cancellation of `c` is not assumed anywhere
in `u3`; it falls out of carrying `ln c` as a free parameter.

`sd(ln theta)`, 25 monthly draws, constant recorded load, full likelihood:

| breath CV | theta | k=0 | k=2 | k=4 | k=8 | k=16 | k=32 | k=64 | **smallest k** |
|---|---|---|---|---|---|---|---|---|---|
| 5% | 3000 | 0.170 | 0.169 | 0.168 | 0.167 | 0.165 | 0.160 | 0.153 | **0** |
| 5% | 600 | 0.817 | 0.576 | 0.471 | 0.366 | 0.275 | 0.203 | 0.151 | **34** |
| 5% | 300 | 1.612 | 0.823 | 0.629 | 0.469 | 0.349 | 0.263 | 0.206 | **70** |
| 10% | 600 | 0.817 | 0.730 | 0.666 | 0.576 | 0.471 | 0.366 | 0.275 | **133** |
| 10% | 300 | 1.612 | 1.230 | 1.035 | 0.823 | 0.629 | 0.469 | 0.349 | **279** |
| 14% | 600 | 0.817 | 0.769 | 0.728 | 0.664 | 0.573 | 0.468 | 0.363 | **261** |
| 14% | 300 | 1.612 | 1.379 | 1.225 | 1.028 | 0.817 | 0.623 | 0.465 | **547** |

The `measured`-level convention (treating the breath visit's level as a known covariate,
which is the author's own "two known marker levels") gives 36 / 62, 141 / 246 and 277 / 482
for the same cells, so the two conventions bracket the answer within about 12%.

**Lengthening the marker record barely moves it.** At `theta = 300` and a 5% breath CV the
smallest `k` is 78 at `n = 13`, 70 at 25, 67 at 50, 54 at 100 and 47 at 200. Sixteen years
of monthly draws buys a reduction from 70 breath tests to 47.

**Breath alone is not a route**: with the marker block dropped entirely, `sd(ln theta)` is
infinite at every capacity, because the contrasts determine `u` and nothing else. That is
the sense in which the question "how few breath tests does a marker record need" is the right
one and "does breath alone suffice" is not.

### 3.1 Why so many: where the breath row points

In the scaled coordinates, the breath contrast at `theta = 300` is
`[-0.4806, -0.1982, 0.8309, -0.1982]` with `abs(cos)` of **0.632** against the null (fibre)
direction and **0.254** against the second-weakest. At `theta = 600` it is 0.849 and 0.082.
The contrast points mostly along the direction that 1.3 shows is worth five per cent when
closed, and only weakly along the one that binds. That is the whole explanation for the
counts above.

### 3.2 Independent re-check

`u4` simulates records in observation space, fits five parameters by scaled trust-region
least squares and reports the robust spread over 300 replicates. No information matrix, no
linearisation.

| theta | k | breath CV | u3 Fisher | robust sd | ratio | collapse (Wilson 95%) |
|---|---|---|---|---|---|---|
| 3000 | 0 | 10% | 0.1696 | 0.1777 | 1.047 | 0/300 [0.000,0.013] |
| 3000 | 32 | 5% | 0.1605 | 0.1638 | 1.021 | 0/300 [0.000,0.013] |
| 600 | 34 | 5% | 0.1979 | 0.2093 | 1.058 | 1/300 [0.001,0.019] |
| 600 | 64 | 5% | 0.1506 | 0.1755 | 1.165 | 0/300 [0.000,0.013] |
| 300 | 32 | 5% | 0.2632 | 0.2655 | 1.009 | 7/300 [0.011,0.047] |
| **300** | **70** | **5%** | **0.1998** | **0.1819** | **0.911** | 6/300 [0.009,0.043] |

Ratios span 0.911 to 1.165 at the paper's own `k_e` log-prior of 0.8, and 0.879 to 1.175 at
a 0.2 control. The plain standard deviation runs to 0.48 where the robust spread is 0.27,
which is why the robust spread is the one quoted.

**One error is kept as a row rather than deleted.** The first version of `u4` held `k_e` at
its prior mean instead of drawing it from the prior. A Fisher variance is a **posterior**
variance, and a simulated spread equals it only when the generating distribution of the
nuisances is the prior the fit assumes. Holding `k_e` fixed silently simulates the
`k_e`-exactly-known design: it returned 0.0782 against `u1`'s `k_e`-exact figure of 0.0813,
and would have been reported as the Fisher being pessimistic by a factor of two.

---

## 4. The negative result, recorded so nobody retries it

A ratio taken across a **load step** is worse than one taken across the **dose cycle**, and
it is a pure lever-arm effect: both amplifications go as `1/(m2 - m1)`.

At `theta = 300` and a 10% breath CV per reading, with the load-step arm read at the dose
train's trough so the dose is not doing the work in either arm:

| arm | m1 | m2 | fold | sd(ln K_m), author's functional | sd(ln u), derived |
|---|---|---|---|---|---|
| dose cycle | 4.677 | 31.229 | **6.68** | **0.732** | 1.048 |
| load cut 50% | 14.521 | 31.229 | **2.15** | **1.272** | 1.766 |
| load cut 25% | 22.601 | 31.229 | 1.38 | 2.637 | 3.577 |

The author reports 1.275 against 0.735, and folds of 6.6 and 2.1. Reproduced to three
figures. The conclusion holds on the derived functional too, and at `theta = 3000` as well
(32.5 against 27.0 in `sd(ln theta)`).

---

## 5. The regime complementarity, and what it does to `11_scope`

`11_scope` closes with: *"no schedule, no excursion and no precision on the dose recovers a
quantity the observable has stopped responding to"*, and, earlier, names where one would
have to act: *"on `sigma`, or on an observable that still responds to `theta` where the
marker has stopped."*

**The breath observable is exactly that observable, and it does respond.** At the untreated
operating point, recorded load held fixed:

| theta | marker elasticity `dln m / dln theta` | breath elasticity `dln F / dln theta = 1 - F` |
|---|---|---|
| 3000 | -0.9154 | 0.1549 |
| 600 | -0.5367 | 0.5277 |
| 300 | -0.3210 | 0.7083 |
| 150 | -0.1729 | 0.8370 |

The two run in opposite directions and **cross at about 20% of nominal**. Over the range
3000 to 150 the marker's elasticity falls by 5.3 and the breath's rises by 5.4.

The usable signal is smaller than that, because `c` must cancel, so only the **contrast**
across the dose cycle is available:

| theta | F(peak) | F(trough) | contrast | signal-to-noise per pair at 5% CV |
|---|---|---|---|---|
| 3000 | 0.8530 | 0.8464 | -0.0067 | 0.111 |
| 600 | 0.5346 | 0.4902 | -0.0444 | 1.225 |
| 300 | 0.3644 | 0.3138 | -0.0506 | 2.114 |
| 150 | 0.2226 | 0.1811 | -0.0415 | 2.916 |

**So the complementarity is real, and it is not enough.** The breath contrast's
signal-to-noise rises by a factor of 26 from nominal to a twentieth of nominal, exactly
against the marker's trend, and the joint Fisher still needs 70 readings at `theta = 300`
and a 5% breath CV, 279 at 10% and 547 at the only CV anybody has published.

The claim in `11_scope` that needs qualifying is therefore **not** the verdict, which
survives. It is the stated **mechanism**, "a quantity the observable has stopped responding
to", which is true of the marker and false as a general statement, and the section's own
"an observable that still responds to `theta`" now has a named candidate with a number
against it. Proposed wording is in `audit/PHASE11_PROPOSED_TEXT.md`, unapplied.

---

## 6. What is not settled

- **Nothing here is checked against a primary source.** `OPEN_QUESTIONS.md` section 1 still
  holds: the container has no egress. The breath test's functional is derived from *this
  repository's own model*, not from the assay literature, and the one published
  repeatability figure (14.0% between-visit, `breath_test_mma`) is used only as a CV, on
  trust, exactly as `08_specification.tex` already uses it.
- **The bolus size relative to the marker pool is unknown** and the model's units are
  illustrative, so section 2.4 is a sensitivity and not a correction.
- **The observation-model assumption is heavier here than elsewhere**, as section 2.4 says.
- The `measured`/`modelled` bracket in section 3 is 12% and has not been narrowed.
