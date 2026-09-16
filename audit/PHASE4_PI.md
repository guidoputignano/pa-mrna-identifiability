# Phase 4. the PI's notes, and a second reader's hypothesis, audited against the code

Six points were put as hypotheses and each was tested against `src/`, the paper and the
earlier audit records. Nothing was built on any of them first. Every number below comes
from a committed script and is reproducible:

```
python audit/f1_dose_model.py     # P1  the dose model
python audit/f2_euler.py          # P2  Euler, and what a between-visit difference contains
python audit/f3_linearisation.py  # P3  his multiplied-out form and his constrained solve
python audit/f4_sampling.py       # P4  the sampling rule
python audit/f5_fibre.py          # P5  curvature, or an exact fibre
```

No script records an expected result in advance. The four numbers quoted from the paper in
`f5_fibre.py` are labelled as quotations, at the file and line they are quoted from, and
the script says whether it reproduces them.

**The headline.** Of seventeen testable claims below, **six hold**, **three are already in
the paper**, and **eight fail**. The eight are not evenly weighted. Three of them
(P1a, P1b, P3b) share one root cause, which is that his balance is written for an impulsive
dose model and carries no `Delta_k m_k` term. Two more (P3d, P3f) are properties of his
recipe rather than of the biology, and both are repairable. P2b and P4b are self-contained.
P5 is the one that matters.

The largest result is P5, and it goes against the paper rather than against him. The
premise he restates is the paper's own (`PI_NOTE.md` p.3, from `06_excitation.tex:66-78`):
the degenerate combination is "recovered through curvature" at a `sqrt(noise)` price. The
second reader's hypothesis is **confirmed**, symbolically and numerically: the combination
is not recovered at all. The fibre is exact, the true profile is flat to machine precision,
and the `eps^4` the paper reports is the curvature of that fibre relative to its own tangent
line. This retires three sentences in Section 6 and one cross-reference in Section 7. It
does **not** touch the paper's headline precision, which is earned at a varying load where
the fibre is closed.

---

## Verdict table

| # | The claim | Verdict | Evidence |
|---|---|---|---|
| P1a | The dose is an impulse `g(m) delta(t-T)` on the balance | **fails as a description of this model** | `src/design/run.py:68-77`: `delta(t)` is a sum of decaying exponentials, half-life 3.5 d, multiplying `g(m)`. `f1` S1 |
| P1b | The marker jumps, `m(T+) = m(T-) - a g(m(T-))` | **fails, and is not consistent with a jump in `Delta`** | across the dose `m` moves by `-4.9e-7` (integrator noise) and `dm/dt` by `-490.6`. `f1` S1 |
| P1c | The jump needs a positivity guard `m(T+) >= 0` | **already in the paper, in a different form** | `Delta >= 0` (min 400 per cycle), and `m=0` is repelling since `dm/dt -> L+b > 0`. `f1` S3a |
| P1d | Bounded jump, slope `<< 1 ~ 1/K_M` | **holds, and is already satisfied by a wide margin** | over a late cycle `m/K_m` in `[0.034, 0.089]`, `g(m) <= 0.082`. `f1` S3b |
| P1e | `theta + Delta_peak <= normal capacity` | **holds, and exposes a real inconsistency in `params.py`** | `theta + Delta_pk = 9400`; `0.313` of normal read one way, `0.157` the other, and `src/common/params.py:140-146` asserts both. `f1` S3c |
| P1f | Samples inside the post-dose settling window must be excluded | **holds; the window is hours, not days** | 99% closes in 3.71 h, 99.9% in 5.57 h. `f1` S3d |
| P2a | Explicit Euler is unsafe across the jump | **holds, and understates it** | it is unsafe everywhere on the record: `h_max = 0.516 h` against a two-daily gap of 48 h. `f2` S1, S2 |
| P2b | A finite difference between visits carries information | **fails** | the visits already sit on the locus to 3.1% of one assay sd; the difference is biased by 28 noise sds as an estimate of `dm/dt`. `f2` S3 |
| P2c | Euler is used in the computation | **already excluded** | no explicit-Euler path exists; `tests/test_models.py:1064-1082` is the test that forbids it. `f2` S4 |
| P3a | Multiplying by `(K_m+m)` linearises the balance | **holds, and is the right move** | it gives exactly the second reader's `beta`. `f3` S1 |
| P3b | His five-term form matches the paper's balance | **fails in three places** | no `Delta m` term; a spurious `L_0` column; and the `psi` terms carry the wrong sign. `f3` S1, S2 |
| P3c | The constraint `first * last = fourth` | **holds** | `c1 c5 = K_M a = c4` identically, and it is what splits the duplicated `m` column. `f3` S3 |
| P3d | "Keep the positive root" | **fails as a rule** | the quadratic is degenerate (`A = 9.7e-32`), both roots are positive, and the spurious one is `1.0e31`. `f3` S4 |
| P3e | His solve, with the right root and `Delta` restored | **holds, and equals the true constrained minimiser** | at a varying load it returns `(100, 0, 3000, 500, 5)` exactly. `f3` S4 |
| P3f | The linear route runs with `psi` unknown | **fails** | `K_M psi` is a product of unknowns. Linear only after the `beta` reparametrisation. `f3` S3 |
| P4a | Synchronous sampling is bad | **already in the paper** | 14, 28, 42, 56 d give 3.51, 3.59, 3.63, 3.68. `f4` S1, and `PHASE2D_FINDINGS.md` §4 |
| P4b | Sampling at 1/2 of the interval is also bad | **fails** | 7 d gives 0.1528 and passes. `f4` S2 |
| P5 | Recovered through curvature, `sqrt(noise)` price | **fails** | the second reader is right: the fibre is exact, the true profile is flat to machine precision, and the `eps^4` is the fibre's curvature relative to its own tangent. `f5` S1-S4 |
| P6 | (what the notes do not cover) | see below | |

---

## P1. The dose model

**What the code runs.** `src/design/run.py:68-77`. `Delta(t)` is `sum_k D exp(-lambda(t-kT))`
over past dose times, `D = 6000`, `T = 14` d, `lambda = ln2/3.5 = 0.198042` /d. At the
periodic steady state the increment peaks at `6400` (`Delta/theta = 2.1333`) and troughs at
`400` (`0.1333`), a 16-fold swing inside one cycle. `design_record` (`:92-109`) never
integrates: it solves the quasi-steady quadratic `k_e m^2 + (k_e K_m + theta + Delta - L) m
- L K_m = 0` for its positive root. There is no jump on `m` anywhere in `src/`.

**Is his jump consistent with a jump in `Delta`?** No, and the difference is the whole
design. A jump in `Delta` leaves `m` continuous and only kinks `dm/dt`: integrating across
the dose at `t = 280` d, `m` moves by `-4.9e-7` (the integrator's own tolerance) while
`dm/dt` moves by `-490.6`. `m` then relaxes to a *new, lower* quasi-steady level and stays
near it for days, because what decays is the increment, on a 3.5 d half-life.

Under his jump the marker is removed instantaneously and then relaxes back with the
*marker's* time constant, `0.806 h`. Matching the two so that both reach the same minimum
(`a = 72.62`), over one 14 d cycle:

| | paper's model | impulsive model |
|---|---|---|
| time-integrated marker reduction below `m_bar` | 48.754 | 0.216 |
| days per cycle the effect exceeds one assay sd | 14.000 of 14 | 0.085 of 14 |
| probability a draw at a random time sees the dose | 1.0000 | 0.0061 |

A factor of **226** in the integrated effect, and a monthly clinic draw would see the dose
about once in 160 visits. The separation that makes this work is `3.5 d / 0.806 h = 104x`.
His model is not a simplification of the paper's; it removes the therapy channel.

**His constraints, restated in the paper's terms.**

- *Positivity.* `Delta(t) >= 0`, satisfied with a minimum of 400 per cycle. The marker needs
  no guard: `dm/dt -> L + b > 0` as `m -> 0`, and the discriminant `B^2 + 4 k_e (L+b) K_m`
  is positive for every `Delta`, so a positive root always exists.
- *Bounded jump, slope `<< 1`.* The paper's counterpart is the fractional saturation
  actually reached. Over a late cycle `m` runs `[3.364, 8.909]`, so `m/K_m` runs
  `[0.0336, 0.0891]` and `g(m)` runs `[0.0325, 0.0818]`. The channel runs at 8.2% of
  `V_max` at the top of the cycle. This is the near-linear regime his condition asks for,
  and it holds with about an order of magnitude to spare.
- *Capacity ceiling.* `theta + Delta_peak = 9400`. **`params.py` does not define "normal"
  once, and the two readings are in the same row.** `src/common/params.py:140-146` says
  `theta` is "about 10% of normal" and, four lines later, that "residual PCC activity is
  under 5% of control in PA". `src/common/params.py:64-71` carries only the second. Read
  the first way, normal is 30000 and the peak sits at 0.313 of it; read the second, normal
  is 60000 and the peak sits at 0.157. The ceiling is not breached either way, but the two
  differ by a factor of two, one row asserts both, and the local and parent models carry
  different `theta` (3000 and 300) so the reference has to be named per model rather than
  once. That is Phase 2's item and it is a genuine find of his.
- *Settling window.* 50% of a displacement closes in 0.558 h, 90% in 1.855 h, 99% in
  3.710 h, 99.9% in 5.565 h. So the exclusion window is **hours**, and what it protects is
  the algebraic balance, not a jump.

*Minor, and a question rather than a finding.* `06_excitation.tex:56-58` pairs "the marker's
own time constant is about 0.8 h" with "99% of a displacement closes in 4.2 h". The 0.8 h is
computed at `m_bar = 10`; the 4.2 h comes from `e5_structural.py` relaxing towards
`m* = 20`. At `m_bar` the 99% time is 3.71 h. Two operating points in one sentence.

---

## P2. Euler

**Time constant and stability limit.** The decay rate is
`lambda(m, Delta) = (theta+Delta) K_m/(K_m+m)^2 + k_e`:

| where | `m` | `Delta` | `lambda` (/d) | `tau` (h) | `h_max = 2/lambda` (h) |
|---|---|---|---|---|---|
| no therapy, operating point | 10.000 | 0 | 29.793 | 0.806 | 1.611 |
| cycle trough in `m`, peak dose | 3.364 | 6400 | **92.981** | 0.258 | **0.516** |
| cycle peak in `m`, trough dose | 8.910 | 400 | 33.666 | 0.713 | 1.426 |
| the instant after a dose | 8.910 | 6400 | 84.249 | 0.285 | 0.570 |

The stiffest point is the *trough*, not the dose instant. Explicit Euler needs
`h < 0.0215 d`. Against that:

| visit gap | `h/h_max` | `|1 - h lambda|` per step |
|---|---|---|
| two-daily | 93.0 | 185.0 |
| weekly | 325.4 | 649.9 |
| the dose interval | 650.9 | 1300.7 |
| monthly | 1413.3 | 2825.6 |
| two-monthly | 2789.4 | 5577.9 |

Run at a monthly step, his recursion reaches `-4.6e5` in three steps and puts the marker
negative on the first. At a weekly and a two-daily step it overflows in 8 and 12 steps. At
a stable step it needs 22315 steps to cross 120 days and then agrees with the balance
(7.1206 against 7.1302). So his warning is right and understated: the step is not
"unsafe across the jump", it is unstable at every clinic spacing and it leaves the physical
region immediately.

**Does a between-visit difference carry information about the locus?** No, on two counts.

1. It is not measuring `dm/dt`. On the monthly design the true `|dm/dt|` at the visits has
   median 0.4327 /d, while `(m_{k+1}-m_k)/tau` has median 0.0362 /d, 8% of it. The bias is
   0.4100, which is **28 times** the noise sd on a difference
   (`sqrt(2) sigma_m / tau = 0.01452`) and 93% of `dm/dt` itself. A gap longer than the
   dose cycle averages the locus's own drift over whole cycles.
2. There is nothing left over for it to measure. At the visits the ODE sits on its own
   algebraic locus to `max |m_ODE - m_balance| = 9.66e-3`, which is **3.1% of one assay sd**
   (0.312). The derivative information the structural result lives on is not merely noisy
   here; it is two orders of magnitude below the assay.

Worth stating plainly, because the paper's phrasing at `06_excitation.tex:62-63` is loose:
`dm/dt` at a visit is **not** zero (median 0.43 /d). What is true is that the *level* is on
the locus to 3% of an assay sd, because the locus itself drifts as the increment decays.

**Is there an explicit-Euler path?** No. A case-insensitive grep over `src/`, `audit/`,
`tests/` and `tools/` for the word or for a fixed-step update returns exactly one hit
outside `f2_euler.py` itself: `tests/test_models.py:1066`, the docstring of the test that
*forbids* it. Integrators actually requested: LSODA 8 times, Radau 3 times.

---

## P3. His linearisation

**The paper's own multiplied-out balance.** Multiplying
`(theta+Delta) m/(K_m+m) + k_e m = L + b` by `(K_m+m)`:

```
0 = -K_m L - K_m b + k_e m^2 + m (Delta + K_m k_e - L - b + theta)
```

and putting the terms that are known sample by sample on one side,

```
Delta_k m_k - L_k m_k = K_m L_k  +  K_m b  +  m_k (b - theta - k_e K_m)  -  k_e m_k^2
```

which is **linear in `beta = (K_m, K_m b, b - theta - k_e K_m, -k_e)` against the regressors
`(L_k, 1, m_k, m_k^2)`**. That is the second reader's `beta` exactly, verified symbolically
(`f3` S1, remainder 0).

**Mapping his symbols.** Against `K_M psi + psi m - L_0 + theta m + a K_M m + a m^2 = 0`:

| his term | paper's term | |
|---|---|---|
| `K_M psi` | `K_m (L+b)` | matches if `psi` is the **total** load `L+b` |
| `psi m` | `(L+b) m` | matches |
| `- L_0` | (nothing) | the balance has **no free constant column**; the constant is `K_m(L+b)`, already present |
| `+ theta m` | `- theta m` | his `theta` is `-theta` |
| `+ a K_M m` | `- k_e K_m m` | matches with `a = -k_e` |
| `+ a m^2` | `- k_e m^2` | matches with `a = -k_e` |
| (absent) | `- Delta m` | **his form has no dose term** |

The last row is the consequential one. `Delta_k m_k` is the therapy channel, which supplies
42x the excitation the marker band does (`PHASE2D_FINDINGS.md` §1). Dropped, his design is
the no-dose-train column of `rem:rankvsprec`. It belongs on the **known** side: `Delta_k` is
read off the dosing record and `m_k` off the assay.

There is also a **sign slip**. His `psi` terms carry the same sign as the clearance terms;
multiplying `theta g(m) + a m = psi` through gives them opposite signs. Run exactly as
written, every coefficient returns negative and his "positive root" rule selects nothing.
`f3` corrects the sign and says so, because the alternative is to report a transcription
error rather than test the method.

**The five questions asked.**

- *Still linear with `psi` unknown?* **No.** `K_M psi = K_M L_k + K_M b`, and `K_M b` is
  unknown times unknown. `psi m = L_k m + b m` is linear but lands on the same column as
  `theta m` and `a K_M m`. The form is linear only after the `beta` reparametrisation, in
  which `K_M b` is carried as one coefficient and never split.
- *Are `psi K_M` and `L_0` separately identifiable?* **No**, when `psi` is a single unknown
  constant: both multiply the constant column, so only `K_M psi - L_0` is identified.
  Adding `L_0` to the model adds nothing the data can see.
- *Where does `b` appear?* Never in a column of its own. It reaches the fit as `K_m b` (the
  constant column) and inside `(b - theta - k_e K_m)` (the `m` column).
- *What is the `m` coefficient really?* `(b - theta - k_e K_m)`. Not `theta`, not `-theta`.
  The load bias and the alternative channel are both inside it, so reading that coefficient
  as a capacity is a three-way confusion.
- *Are columns 3 and 4 the same regressor?* **Yes**, literally. Both are `m_k`. Only
  `theta + a K_M` is identified by the columns, and his constraint is what splits them.
  The constraint is correct: `c1 c5 = K_M a = c4`.

**Rank of his design `[psi, -1, m, m, m^2]`:**

| reading | load CV | rank | `s_min/s_max` |
|---|---|---|---|
| `psi = L_k`, recorded and varying | 0% | 3 of 5 | 1.76e-17 |
| `psi` one unknown constant | 0% | 3 of 5 | 1.97e-17 |
| `psi = L_k`, recorded and varying | 10% | **4 of 5** | 5.29e-18 |
| `psi` one unknown constant | 10% | 3 of 5 | 2.08e-17 |

**Yes, the unknown-`psi` degeneracy is the same object as `prop:rank` and as the
constant-load fibre of P5.** At a constant `psi` the deficiency is two: the duplicated `m`
column (closed by his constraint) and the collapse of the `psi` column onto the constant
column, which is `prop:rank` and which no constraint on the coefficients repairs.

**His solve, implemented as written.** Conditioning on a scalar `w` is only possible for two
of the five coefficients, because the duplicated `m` column makes the `v`-block singular
otherwise:

| `w` held | `v`-block `s_min/s_max` | |
|---|---|---|
| `c1 = K_M` | 7.75e-18 | unformable |
| `c2 = L_0` | 7.25e-18 | unformable |
| `c3 = theta` | 1.53e-2 | formable |
| `c4 = a K_M` | 1.53e-2 | formable |
| `c5 = a` | 4.96e-17 | unformable |

Taking `w = c4`, the constraint becomes `A w^2 + B w + C = 0` with
`(A, B, C) = (9.70e-32, -1.000000, 500.0000)`. **`A` vanishes identically**, because `v(w)`
moves only in the coefficient of the column `w` came from. So the equation is linear with
the single root `w = 500`, and the "second root" `1.03e31` is an artefact of dividing by
`A`. Both are positive, so "keep the positive root" does not disambiguate, and taken at
face value it selects the artefact.

With the right root, on a monthly design at load CV 10%, noiseless:

| case | `K_M` | `L_0` | `theta` | `a K_M` | `a` | residual SS |
|---|---|---|---|---|---|---|
| truth | 100 | 0 | 3000 | 500 | 5 | |
| **`Delta` restored** | **100.000000** | **-0.000000** | **3000.000000** | **500.000000** | **5.000000** | 3.3e-18 |
| `Delta` dropped | -6.144816 | 1093.785 | 167.4956 | -157.5039 | 25.632 | 1.26e5 |
| `Delta` restored, `b_true = 32.27` | 100.000000 | 3227.000 | 2967.730 | 500.000000 | 5.000000 | 2.4e-18 |
| `Delta` restored, marker CV 5% | 107.5683 | -2139.178 | 5042.861 | -1355.104 | -12.598 | 7.1e7 |
| `Delta` restored, load CONSTANT | 42.7026 | 18491.42 | 3286.487 | 213.5132 | 5.000000 | **7.1e-19** |

Reading the rows:

1. **His method works**, on two conditions: the `Delta m` term restored and the recorded load
   varying. It then returns the truth exactly and coincides with the true constrained
   least-squares minimiser to 7 figures. That is not an accident: the residual depends only
   on `(c1, c2, c3+c4, c5)`, so the unconstrained fit in those four is unique and the
   constraint splits `c3` from `c4` uniquely. His solve **is** the constrained minimiser.
2. Dropping `Delta` puts `K_M` negative and `theta` at 5.6% of truth. The fit leaves the
   physical region, which is at least detectable from the fit alone.
3. With `b_true = 32.27` the parametrisation recovers everything: `L_0 = K_m b = 3227` and
   `c3 = theta - b = 2967.73`, so `b = c2/c1` and `theta = c3 + c2/c1`. A genuinely useful
   property of his form that the paper's parametrisation does not display as directly.
4. **At the paper's marker CV the linear route collapses**: `a = -12.60` against a truth of
   5. `m` appears in the regressors and in `m^2`, so this is an errors-in-variables problem
   and the naive linear fit is badly biased. This is the first item of Phase 1 and it is the
   reason a variable-projection fit is worth having.
5. At a constant load the residual is **7.1e-19**, an exact fit, at
   `K_m = 42.70`, `b = 433.2` (from `L_0/K_M = 18491/42.70`), `k_e = 5.000000`. That point
   is on the P5 fibre: `K_m(433.2) = 100 * 322.727/755.9 = 42.70` and
   `theta - b = 3286.5`. Which point it lands on is decided by round-off, not by data.

---

## P4. Sampling

The table reproduces `PHASE2D_FINDINGS.md` §4 on d12's minimal footing (`theta` and `K_m`
free, priors only on `ln k_e` at 0.8 and the load record at 5%), load CV 10%, 730 d:

| step (d) | step/`T` | `n` | distinct phases | `Delta/theta` range | `sd(ln theta)` | |
|---|---|---|---|---|---|---|
| 1.0 | 1/14 | 730 | 14 | 1.9810 | 0.1487 | PASS |
| 2.0 | 1/7 | 365 | 7 | 1.9476 | 0.1499 | PASS |
| 3.5 | 1/4 | 209 | 4 | 1.8833 | 0.1508 | PASS |
| **7.0** | **1/2** | 105 | 2 | 1.6333 | **0.1528** | **PASS** |
| 10.0 | 5/7 | 73 | 7 | 1.9352 | 0.1565 | PASS |
| **14.0** | **1** | 53 | **1** | **0.1333** | **3.5050** | **fail** |
| 14.5 | 29/28 | 51 | 28 | 1.9861 | 0.1612 | PASS |
| 21.0 | 3/2 | 35 | 2 | 1.6021 | 0.1637 | PASS |
| **28.0** | **2** | 27 | **1** | **0.1333** | **3.5873** | **fail** |
| 30.4 | 2.171 | 25 | 25 | 1.8438 | 0.1697 | PASS |
| **42.0** | **3** | 18 | **1** | **0.1333** | **3.6267** | **fail** |
| **56.0** | **4** | 14 | **1** | **0.1333** | **3.6834** | **fail** |
| 60.0 | 30/7 | 13 | 7 | 1.9352 | 0.1868 | PASS |

**His first claim holds.** It is the paper's own rule and it extends cleanly to every
multiple of `T`.

**His second claim fails.** `T/2 = 7 d` passes at 0.1528. So do `T/3`, `T/4`, `T/7` and
`T/14`. The half-interval is the *worst of the passing* schedules, not a failing one.

**The correct rule, in the two quantities asked for.** A step `h` visits
`q = ` denominator of `h/T` in lowest terms distinct phases of the dose cycle, densely when
`h/T` is irrational. But the phase **count** is not what binds. Holding `n = 25` and the
count at exactly 2, and moving only the separation of the two phases:

| two phases, separated by | `Delta/theta` range | `sd(ln theta)` |
|---|---|---|
| 0.5 d | 0.3219 | 1.2449 |
| 1.0 d | 0.4927 | 0.6443 |
| 2.0 d | 0.7874 | 0.3389 |
| 3.5 d | 1.1333 | 0.2239 |
| 7.0 d | 1.6333 | 0.1700 |

A factor of 7.3 in the precision at a constant phase count. Conversely, 2, 3, 4, 7 and 14
evenly spread phases at fixed `n` give 0.1709, 0.1698, 0.1701, 0.1746, 0.1855: adding
phases buys nothing once the range is covered, and eventually costs, because the visits get
spent on phases that are close together in `Delta`.

So: **the binding quantity is the range of `Delta/theta` the visited phases span, and the
phase count matters only through `q >= 2`.** Every failing schedule has range 0.1333 (which
is not zero only because the dose train is still building up at the start of the record);
every passing schedule has range `>= 1.6021`. The widest available is 2.0000.

**One thing the paper does not say, and should.** "Move the schedule half a day" is not
symmetric:

| step (d) | offset from `T` | `Delta/theta` range | `sd(ln theta)` |
|---|---|---|---|
| 13.50 | -0.50 | 1.9953 | 0.1604 |
| 13.95 | -0.05 | 1.8738 | 0.1858 |
| 13.99 | -0.01 | 1.8748 | 0.1839 |
| 14.00 | 0 | 0.1333 | 3.5050 |
| 14.01 | +0.01 | 0.1998 | **1.4583** |
| 14.05 | +0.05 | 0.8166 | **0.3178** |
| 14.50 | +0.50 | 1.9861 | 0.1612 |

Stepping just **short** of the dose interval walks each visit backwards into the low-`Delta`
tail and keeps most of the range. Stepping just **long** keeps every visit on the post-dose
plateau where `Delta` has barely moved, and 14.05 d already fails the 0.20 gate while
13.95 d passes. If a protocol has to be near the dose interval, it should be *under* it.

---

## P5. The premise: curvature, or an exact fibre

**The second reader is right.** Both symbolically and numerically, and the consequence is
larger than the hypothesis as stated.

**Symbolic (`f5` S1, sympy).** Solving the balance for the increment that produces a level
`m`,

```
Delta(m) = A/m + B + C m ,   A = K_m (L+b),  B = L + b - theta - k_e K_m,  C = -k_e
```

with remainder exactly 0. `{1/m, 1, m}` are independent as functions of `m`, so
**`(A, B, C)` is everything a constant-load record of `(m_k, Delta_k)` pairs can contain**:
four unknowns, three invariants, a one-parameter fibre. Fixing `(A, B, C)` at the truth and
solving gives

```
K_m(b) = K_0 L/(L+b) ,    theta(b) = theta_0 + b + k_e (K_0 - K_m(b)) ,    k_e unchanged
```

which is the second reader's curve, and its tangent at `b = 0` is `prop:rank`'s null vector
exactly (`difference = [0, 0, 0, 0]`, asserted in the script).

**Numeric (`f5` S2, monthly design, ordinate sd 21.9612).**

| `b` displacement | chi2 along the straight line `n` | chi2 along the exact curve |
|---|---|---|
| 0.01 | 4.41e-15 | 1.86e-25 |
| 1.0 | 4.44e-07 | 1.85e-25 |
| 30.0 | 4.29e-01 | 2.19e-25 |
| 36.9 | 1.03e+00 | 2.25e-25 |
| 200.0 | 4.08e+03 | 4.55e-25 |

Along the line the local slope `d ln chi2 / d ln eps` is 4.0000 over four decades, exactly
reproducing the paper's power law. **Along the curve chi2 is machine zero at every
displacement, out to `b = 200`, which is 62% of the recorded load.**

**The true profile (`f5` S3).** Re-optimising `(theta, K_m, k_e)` at each fixed `b`:

| `b` | chi2 along the line | chi2, TRUE profile | `theta_hat` | `K_m_hat` | `k_e_hat` |
|---|---|---|---|---|---|
| 0 | 1.8e-25 | 1.2e-25 | 3000.0000 | 100.0000 | 5.00000 |
| 30 | 4.29e-01 | 1.5e-25 | 3072.5258 | 91.4948 | 5.00000 |
| 100 | 8.80e+01 | 2.1e-25 | 3218.2796 | 76.3441 | 5.00000 |
| 200 | 4.08e+03 | 3.1e-25 | 3391.3043 | 61.7391 | 5.00000 |
| 300 | 2.41e+05 | 4.3e-25 | 3540.8759 | 51.8248 | 5.00000 |
| -200 | 2.81e+02 | 1.9e-26 | 1985.1852 | 262.9630 | 5.00000 |

The optimiser lands on the closed-form fibre point to every digit printed. **The true
profile is flat. The `eps^4` is the curvature of an exact fibre relative to its own tangent
line, and it carries no information.**

**Was the `+/-36.9` a line or a profile? A line.** Three checks:

1. The printed pair is internally consistent with a quartic:
   `(3.84 / 2.08e-6)^(1/4) = 36.86`. A quartic in the displacement is what a straight line
   through a curved fibre gives, and is not what a profile gives.
2. `+/-31.6` reproduces exactly as `1.96 * 0.05 * 322.7273 = 31.6273`, the load prior's own
   95% half-width. Design-independent.
3. **The two numbers are not on the same design.** On the accessible monthly design the
   same line construction gives `chi2/eps^4 = 4.44e-7` and a half-width of **`+/-54.24`**,
   not 36.9. `2.08e-6` matches a design of roughly 100 visits. The closest in the
   sweep is weekly x 2 y at `1.91e-6` and `+/-37.67`, which is 8% out in the constant and
   2% out in the half-width, so the exact design is not pinned. What is pinned is that it is
   **not the monthly one**: that would need `4.44e-7`, a factor of 4.7 away. So
   `rem:second_order` compares a line profile from a design of roughly weekly size against a
   design-independent prior width, in defence of a covariance quoted on the **monthly**
   design. On the monthly design the same two numbers are 54.2 against 31.6, a factor of
   1.72, not "within a fifth".

This is the `R8` failure mode again (`REVIEW_LOG.md`): a number that is right about its own
computation and attached to a sentence with a different footing.

**With noise (`f5` S4, 200 replicates, constant load, no priors, free 4-parameter fit).**
23% of replicates return `theta_hat <= 0`; `k_e_hat` has median `-0.879` against a truth of
5. Relative inter-quartile ranges:

| quantity | relative IQR |
|---|---|
| `A = K_m(L+b)` | **0.084** |
| `B` | 0.307 |
| `C = -k_e` | 17.39 |
| `b` | 13.90 |
| `theta` | 3.56 |
| `K_m` | 0.96 |

Only the invariant `A` is pinned. This is what an exactly flat direction looks like in
practice, and it is not what "recovered through curvature" predicts.

**A second-order point that cuts the other way (`f5` S3c).** `C = -k_e` is an *exact*
invariant, so `k_e` is structurally identified at a constant load with no fibre at all. It
is still barely estimable: over the marker band the design visits, `[3.51, 8.75]`, a
2.49-fold range, `cond{1/m, 1, m} = 59.9`, and the Monte Carlo IQR above is 17.4. So
"exactly identified" and "usefully estimated" are as far apart here as `rem:rankvsprec`
already says, and this time the gap runs in the opposite direction. The paper's `k_e` prior
is not redundant.

**Where the fibre ends (`f5` S3).** Nowhere, inside the model: it runs to `b -> +infinity`
with `theta`, `K_m` and `k_e` all positive. `theta(b)` reaches zero at `b = -215.15`,
`-66.67%` of the recorded load. What closes it is a bound on `b`:

| bound on `|b|` | `theta` interval | fold |
|---|---|---|
| 5% of `L` (the paper's load prior) | [2957.55, 3039.95] | 1.028 |
| 10% | [2912.17, 3077.73] | 1.057 |
| 20% | [2810.45, 3147.88] | 1.120 |
| 50% | [2338.64, 3328.03] | 1.423 |

and with `b >= 0` imposed, `[3000.00, 3039.95]`, `[3000, 3147.88]`, `[3000, 3328.03]`. So
the practical damage at the paper's own 5% load prior is a 2.8% fold in `theta`, which is
small. **The prior is doing that work, and the paper should say so rather than say the
curvature is.**

**Does a varying load close it (`f5` S5)?** Yes, and here are the rates:

| load CV | rms residual on the fibre at `b = 30` | 95% TRUE profile half-width in `b` |
|---|---|---|
| 0% | 2.06e-12 | `> 1e5` |
| 0.5% | 0.1642 | 2108.5 |
| 2% | 0.5694 | 508.8 |
| 10% | 2.9314 | 96.0 |

Even at CV 10% the data's own interval in `b` is `+/-96`, three times looser than the
prior's `+/-31.6`. It does not matter much for `theta` (`PHASE2D_FINDINGS.md` section 2, two-daily x 2 y at
load CV 10%: 0.1474 with the load prior, 0.1483 without), but "at any real day-to-day variation the
data fill it instead" (`rem:rankvsprec`) is generous: the data close the *rank* at 0.5% and
only overtake the prior's *width* somewhere well past 10%.

**Is the fibre an artefact of the algebraic locus (`f5` S6)?** No. Integrating the real ODE
at the truth and at fibre points and sampling at the visits:

| `b` | max `|m_fibre - m_true|` | as a fraction of one assay sd |
|---|---|---|
| 30 | 7.70e-4 | 0.25% |
| 100 | 2.14e-3 | 0.69% |
| 200 | 3.46e-3 | 1.11% |

The fibre is exact for the balance and survives the ODE two orders of magnitude below the
assay.

**Reconciling with `StructuralIdentifiability.jl` (`f5` S7).** No contradiction, and the
paper already has the right reconciliation at `06_excitation.tex:55-65`. The Julia verdict
is about the transient from a generic initial condition, where every time derivative of the
output is available. On the locus `dm/dt = 0`, so every higher derivative is `0 = 0` and
those equations are simply absent. `tau = 0.0336 d` against a 3.5 d half-life, a separation
of 104x, is what puts the clinic on the locus. What changes is only the sentence that
follows: the obstruction is not "real but priced at `sqrt(noise)`". At a constant load it is
an exact fibre priced at whatever bound you can put on `b`; at a varying load it is closed
by the data.

**Three fibres, and they are different objects.** Worth writing down because the paper now
has three and calls two of them the same word:

| | unknowns | observables | dimension | `k_e` along it |
|---|---|---|---|---|
| `rem:fibre` (`04_local.tex:176`) | `(theta, K_m, k_e)` | `(a, b)` of an affine tangent | 1 | **varies** |
| `prop:rank` / this one | `(theta, K_m, k_e, b)` | `(A, B, C)` on the locus | 1 | **fixed** |
| `rem:modelB` (`06_excitation.tex:306`) | `(theta, K_m, k_e, K_e)` | saturating-channel exchange | 0 (two points) | exchanged |

---

## P6. What the notes do not address

His note is entirely about the **local identification model**. Nothing in it touches:

1. **His own second question in `PI_NOTE.md`, whether the reduction error can be bounded a
   priori.** The 374-against-61 gap is `eq:parent_intercept` (`03_dynamic.tex:242`): the
   tangent's intercept is `theta g(m_bar)^2 + [alpha(m_bar) - m_bar alpha'(m_bar)]`, and the
   second bracket vanishes only when `alpha` is a straight line through the origin. None of
   P1 to P4 bears on it. Worse, **P3 makes it sharper in the wrong direction**: multiplying
   by `(K_m+m)` is exact only because the alternative channel is `k_e m`. Under the parent's
   `alpha(m) = R_0 u/(u + rho_0(1+m/M_i))` the multiplied-out form is not polynomial in `m`
   at all, so his whole linear route is available only in the model the paper says is *not*
   a reduction of the parent. That is a real observation to take back to him: the estimator
   he proposes is structurally tied to the hypothesis whose reduction error he wants
   bounded.
2. **The cohort result.** Median achieved precision 0.192 with quartiles 0.132 and 0.329
   over 4000 simulated patients, just over half meeting the 0.20 target, and the shortfall
   systematic rather than random (`PI_NOTE.md` p.2).
3. **The decision problem.** `theta_crit` moving by a factor of 24 over the defensible range
   of its two chosen constants, the sign flipping for the same patient with perfect data,
   and the sign-error rates staying at about two in five on the new design
   (`PHASE2D_FINDINGS.md` §5).
4. **The missing measurement.** `Delta` characterised to about 12%, never measured, which is
   the single requirement the paper cannot discharge.
5. **`kappa`**, never assigned, and the carnitine route's two equations in four unknowns.
6. **`rem:modelB`**, the one configuration `StructuralIdentifiability.jl` does not return as
   globally identifiable.

His three requests (physical acceptability, a design output in visits against sd, and the
phase clock) are all Phase 2 and Phase 3 work and none of them is blocked by anything above.

---

## What P5 costs the paper, if you accept it

Nothing in the headline. `0.1697` on the monthly design is computed at load CV 10% with the
dose train and a `k_e` prior, where the fibre is closed and the load prior is not
load-bearing for `theta` (0.1474 with it against 0.1483 without, on the two-daily x 2 y
design at load CV 10%, `PHASE2D_FINDINGS.md` section 2). The design recommendation, the sampling
rule, the 42x excitation and the specification are all untouched.

What fails is a mechanism sentence and a diagnostic, in four places. Listed here, **not
edited**, pending your go-ahead:

| file:line | what it says | why it fails |
|---|---|---|
| `06_excitation.tex:66-78` | "flat only to first order ... recovered through curvature rather than through slope ... the recoverable displacement scales as the square root of the noise level ... Halving the assay CV buys a factor of `sqrt2` along this direction" | the true profile is flat, so nothing is recovered and halving the CV buys nothing along this direction |
| `06_excitation.tex:497-516` (`rem:second_order`) | the `+/-31.6` against `+/-36.9` comparison, and "`chi2 ~ eps^4` means the data do pin it, at second order, and the filter is understating what the record contains" | the diagnostic gives the wrong answer here: it reads a curved exact fibre as data-pinned. Separately, the two numbers are on different designs |
| `07_design.tex:141` | "along the near-null direction of `prop:rank` the gradient is third order" | true of the line, vacuous on the fibre, where every derivative is zero |
| `07_design.tex:417-419` (`rem:no_delta`) | "the second is `rem:second_order`, where a first-order filter reports a curvature-held direction as prior-held" | exactly inverted. The direction *is* prior-held; it is the `eps^4` diagnostic that calls it curvature-held |

`rem:rankvsprec` (`06_excitation.tex:439-490`) is **correct as it stands** ("that is the
load-record prior filling the deficient direction"). Only its forward reference to
`rem:second_order` needs to go.

Two audit records carry the same claim and will need retraction notes:
`PHASE3_CHECKS.md:301-302` and `PHASE2C.md:141-150`. So does the closing commentary of
`audit/e5_structural.py` `q3`, whose *computation* is correct and whose *interpretation* is
the one being retracted.

---

## Questions to take back to him

1. **Is the impulsive dose model a proposal or a reading of ours?** If he believes the
   therapy removes metabolite rather than adding capacity, that is a different biology and
   worth settling first: mRNA-3927 delivers PCC mRNA, which is a capacity increment, and
   the 226x difference in integrated effect means the two models do not share a design.
2. **In his `K_M psi + psi m - L_0 + ...`, is `psi` the total load or the recorded one, and
   is `L_0` intended to be our unrecorded load `L_0`?** The two cannot both be free: if
   `psi` is the total, `L_0` has no column; if `psi` is the recorded load, `L_0` is
   `K_m b` and not a load at all.
3. **Was the `Delta_k m_k` term dropped deliberately?** If he was deliberately analysing the
   no-therapy case, his conclusions are right for that case and his solve returns
   `K_M = -6.14`, which is the correct diagnosis of it.
4. **Does he agree the sign on the `psi` terms is a transcription slip?** As written the
   solve has no positive root.
5. **"1 or 1/2 of the interval": which did he mean?** At `T` the design fails by a factor of
   22; at `T/2` it passes. If he meant something like "any rational ratio with a small
   denominator", the data say the opposite: small `q` is fine, and what matters is the
   `Delta` range.
6. **One `params.py` row says `theta` is "about 10% of normal" and, four lines later, "under 5% of control".** Which
   reference does he want the acceptability check written against? The two differ by a
   factor of two and they are the denominator of the ceiling constraint he asked for.
7. **Does he accept the exact fibre, and with it that the `sqrt(noise)` sentence goes?** The
   honest replacement is that at a constant load the deficiency is closed by a bound on the
   load bias, and at any real day-to-day variation it is closed by the data. Both are
   defensible and neither needs curvature.
8. **His second question, the a priori reduction bound, is untouched by this round**, and
   P3 tightens it: the linear route he proposes exists only for the straight-`alpha`
   hypothesis. Does he want that pursued next, or the design outputs first?

---

## Status

Phase 0 is complete and nothing has been built on top of it. No prose has been edited, no
figure regenerated, no constant changed, and no retraction written into `REVIEW_LOG.md`
yet. Phases 1, 2 and 3 are waiting at their gates.

Two things are drafted and deliberately not committed to the paper's apparatus, because
both would embed a Phase 0 conclusion before the gate:

- a test in `tests/test_models.py` asserting that the constant-load fibre is exact (that
  `chi2` along `K_m(b) = K_0 L/(L+b)`, `theta(b) = theta_0 + b + k_e(K_0 - K_m(b))` stays at
  machine zero out to `b = 200`, while `chi2` along the tangent line is quartic). It locks a
  fact about the model rather than a sentence about it, and it is the natural guard against
  this finding regressing;
- a retraction entry, `R14`, for `REVIEW_LOG.md`, in the form the existing entries use: what
  was claimed, what is true, how it was found, and what class of failure it belongs to. On
  the repository's own taxonomy it is a new class and worth naming: **a diagnostic that
  answers a question the model cannot pose.** The `eps^2`-against-`eps^4` test is a correct
  computation of the wrong object, because it profiles along a straight line through a
  curved exact fibre, and every numerical check in this repository passes it. What caught
  it was a reader proposing a different parametrisation, which is the same detector that
  caught `R8`: two descriptions of the same object that do not agree.
