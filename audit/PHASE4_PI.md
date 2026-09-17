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
python audit/f6_fibre_adversarial.py   # P5  the adversarial checks, reproducible
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
second reader's hypothesis is **confirmed**: on the quasi-steady locus at a constant
recorded load the degeneracy is an **exact one-parameter fibre**, the true profile is flat,
and the `eps^4` the paper measures is the curvature of that fibre relative to its own
tangent line. Two clauses in Section 6 are false, one cross-reference in Section 7 is
inverted, and `eq:nullvec`'s derivation drops the load bias. It does **not** touch the
paper's headline precision, which is earned at a varying load where the fibre is closed.

**This finding was then adversarially stress-tested before anything was built on it**
(`audit/PHASE4_STRESS.md`): eight independent lenses, each prompted to refute it, each
re-run by a second agent. Claims (1) to (4) survived all eight. **The audit's own claim (5),
about where the paper's `+/-36.9` came from, was refuted and is withdrawn in its strongest
form**, along with seven smaller errors of this document's, including a bug in
`audit/f5_fibre.py`, two single-seed numbers quoted without their spread, and a
justification for claim (2) that was a non-sequitur. All of them are recorded as C1 to C11
at the end of the P5 section. The core survived; this audit did not.

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
| P5 | Recovered through curvature, `sqrt(noise)` price | **fails** | the second reader is right: on the locus at constant load the fibre is exact, the true profile is flat, and the `eps^4` is its curvature against its own tangent. Confirmed by eight adversarial lenses. `f5` S1-S4, `PHASE4_STRESS.md` |
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

> **Clarification received, and what it retires.** The PI's equation was an illustrative
> sketch of the method, not a formal model. On his instruction it is no longer analysed as
> a model of its own: the questions about `psi`, the missing `Delta_k m_k` term and the
> signs are withdrawn, and the symbol-mapping table below is kept only as the record of how
> the sketch was read, not as a list of defects in it. **What survives, and is his, is the
> idea: multiply the balance out to get a linear form.** That is where the `beta`
> reparametrisation comes from, and it is right. The rest of P3 stands on its own terms --
> the multiplied-out balance, the exact `beta`, and the constrained solve -- and none of it
> depended on the sketch being formal.

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

**Verdict: the core holds and one sub-claim was mine and was wrong.** The finding was put
through eight independent adversarial lenses, each prompted to refute it and each re-run by a
second agent (`audit/PHASE4_STRESS.md`). Claims (1) to (4) survived every one. Claim (5), the
forensic claim about where the paper's `+/-36.9` came from, was refuted in its strongest form
and survives only in a weaker one. The corrections are set out at the end of this section and
in the stress-test record; the numbers below are the ones that survived.

### What is claimed

**(1)** Solving the quasi-steady balance for the increment that produces a level `m` gives,
with no remainder,

```
Delta(m) = A/m + B + C m ,   A = K_m(L+b),  B = L+b-theta-k_e K_m,  C = -k_e
```

**(2)** The record determines `(A, B, C)` and nothing else, and determines them
whenever it visits three distinct levels. Both halves need saying, because the first
draft of this section gave a reason that carries neither. *Nothing more*: `Delta(m)`
depends on the parameters only through `(A, B, C)`, and since `A = K_m(L+b) > 0` and
`C = -k_e < 0` the map is strictly decreasing from `(0, infinity)` onto the reals, hence a
bijection, so the level a known `Delta` produces is a function of `(A, B, C)` alone, with no
branch, ordering or positivity freedom left over. *Nothing less*: three distinct levels pin
`(A, B, C)`, because `det[[1/m_i, 1, m_i]] = -(m_1-m_2)(m_1-m_3)(m_2-m_3)/(m_1 m_2 m_3)`, a
Vandermonde, is non-zero exactly when the `m_i` differ. So four unknowns against three
invariants leave a one-parameter fibre,

```
K_m(b) = K_0 L/(L+b) ,  theta(b) = theta_0 + b + k_e (K_0 - K_m(b)) ,  k_e unchanged.
```

**(3)** The fibre is exact on the locus, not second order.

**(4)** So the `eps^4` the paper measures is the curvature of that fibre relative to its own
tangent line, and "recovered through curvature" and the `sqrt(noise)` price are false.

**(5)** The paper's `+/-36.9` was computed along a ray, not as a profile.

### The evidence, by lens

**Symbolic, two independent eliminations.** Solving for `Delta` and collecting in
`{1/m, 1, m}` leaves remainder `0`. Eliminating `m` instead, by a resultant between the
truth's quadratic and a candidate's and requiring it to vanish identically in the known
`Delta`, sympy returns the single solution `{A2: A1, B2: B1, k_e2: k_e}`, and the Groebner
basis of the coefficient ideal contains `(B1-B2)^3`, `(k_e-k_e2)^4` and `k_e^2 (A1-A2)^2`.
**There is no fourth invariant.** The Jacobian of `(theta, K_m, k_e, b) -> (A, B, C)` has rank
3 with minors `{-(L+b), 0, K_m, K_m k_e + L + b}`, and the analytic Jacobian of the level map
on the real records has fourth singular value `9.0e-17` (monthly), `9.4e-16` (weekly),
`5.1e-16` (two-daily) against a first of `~1.99`.

**An identity, not a coincidence.** `R_fibre(m) = R_true(m) * [1 + K_0 b/(K_0 L + (L+b) m)]`,
so the fibre's residual vanishes wherever the truth's does, for every `b`, every `Delta` and
every `m`. It follows immediately that **any** true profile leaving `theta` and `K_m` free is
exactly flat.

**Exact arithmetic settles the artefact question.** The worry is that `chi^2 ~ 1e-25` is
cancellation in the closed-form quadratic root. It is not. sympy returns the integer `0`, not
something small. In exact rationals (`clearance_A(10) = 3550/11`) the residual is the rational
`0` in eighteen cases, at fibre points `K_m = 8875/97, theta = 298035/97` (`b = 30`) and
`K_m = 1420/23, theta = 78000/23` (`b = 200`). Decisively, re-running the pipeline in mpmath
at 15 / 30 / 60 / 120 digits gives `chi^2 = 1.86e-25 / 2.94e-55 / 1.49e-115 / 2.08e-235`:
**a fall of 210 orders of magnitude as precision goes from 16 to 120 digits.** A zero
manufactured by cancellation cannot track working precision like that.

**The observable itself is invariant.** Stronger than the load-scale argument the audit first
made: the predicted marker is bit-identical along the fibre, `max|m(b) - m(0)| = 0.0` in
floating point at `b = 1, 10, 30, 100, 200, 1000, -100, -200`. Every likelihood whose noise
acts on `m` with parameter-free noise constants is therefore exactly flat, whatever its shape.
Confirmed directly: a log-normal (CV 5%) `m`-scale profile moves `chi^2` by at most `4.7e-4`
out to `b = +5000`, against `chi^2(0) = 19.82`; a consistently heteroscedastic load-scale fit
carrying its Jacobian term returns `-26.77721254` at `b = 0, 30, 100, 200, -100`, identical to
eight decimals.

**The `sqrt(noise)` claim, refuted quantitatively.** At a constant load, halving the assay CV
buys nothing at all: the true profile half-width is infinite at `sigma`, `sigma/2` and
`sigma/4`. Where the data do work, at load CV 0.5 / 2 / 10%, the half-width scales as
`sigma` **linearly**, not as its square root: the ratio per halving is `0.43` to `0.54`
against `0.7071` for a square-root law and `0.500` for a first-order one. (An earlier draft
quoted an exponent of `0.98`. That was one seed and two digits, and the robust statement is
the qualitative one.)

**The ODE does not rescue it.** The departure of the ODE record from the locus is perfectly
systematic (the pattern at `b = 1` and at `b = 1000` is the same unit vector to
`cos = 0.99999926`), so it could in principle be estimable. It is not: re-optimising the other
three parameters absorbs `99.9934%` of it, and the profile Fisher information about `b` from
the transient alone gives `sd(b_hat) = 3.90e7`, which is 2.4 million times the load prior's
`+/-31.6`. Matching the prior would take `1.4e14` monthly draws.

**No profile construction reproduces `36.9`.** A sweep of 88 designs x 7 sigma conventions x 3
residual scales, plus restricted profiles (holding each of `K_m`, `theta`, `k_e` and every
subset), profiles with each prior inside `chi^2`, profiles at a varying load, profiles on the
ODE record, and profiles with the assay error moved onto the abscissa. Three land near 36.9
and all three are disqualified on inspection; the closest, a `theta`-only profile on an offset
monthly grid giving 36.84, is exactly **quadratic** in `b` (so it cannot also yield a constant
`chi^2/eps^4`) and holds `K_m` exactly known, which `07_design.tex:35-39` forbids.

**The exact fibre is not already in the paper.** `K_m(b) = K_0 L/(L+b)` appears nowhere
outside the audit. And there is a reason it is not derivable from the paper as written: the
elimination at `06_excitation.tex:371` drops the bias, writing `(theta+Delta) u = L - k_e m`
instead of `L + b - k_e m`, so `eq:nullvec` is the `b = 0` value of a field that in fact
depends on `b`. Carrying `b` through gives `n = (L + b + k_e K_m, -K_m, 0, L + b)`, and
integrating **that** field reproduces the fibre with symbolic difference exactly `0`;
integrating the paper's frozen-`L` vector instead gives `K_m = K_0 exp(-b/L)`, which leaves the
fibre at second order. That is a real defect in the proof and it should be fixed whatever is
decided about the rest.

### What the fibre costs, and what closes it

Nothing in the model closes it: the fibre runs to `b -> +infinity` with `theta`, `K_m` and
`k_e` all positive. `theta(b)` reaches zero at **`b = -272.7273`**, `-84.51%` of the recorded
load. What closes it is a bound on `b`:

| bound on `\|b\|` | `theta` interval | fold |
|---|---|---|
| 5% of `L` (the paper's load prior) | [2957.55, 3039.95] | 1.028 |
| 10% | [2912.17, 3077.73] | 1.057 |
| 20% | [2810.45, 3147.88] | 1.120 |
| 50% | [2338.64, 3328.03] | 1.423 |

So at the paper's own 5% load record the practical damage is a 2.8% fold in `theta`. **The
prior is doing that work, and the paper should say so rather than say the curvature is.**

A varying load closes it from the data, and here are the rates (monthly design, true profile
in `b`, `Delta chi^2 = 3.84`):

| load CV | rms residual on the fibre at `b = 30` | 95% profile half-width in `b`, 8 seeds |
|---|---|---|
| 0% | 2.06e-12 | unbounded |
| 0.5% | 0.123 | 2108 to 4155, median 3500 |
| 2% | 0.491 | 424 to 690, median 593 |
| 10% | 2.430 | 80 to 134, median 112 |

Even at CV 10% the data's own interval is three to four times looser than the prior's
`+/-31.6`. Those are eight-seed ranges, not single realisations: an earlier draft of this
table quoted one seed, and that seed is the **minimum** at every CV, so it understated the
width by a third to a half. That is correction C5 again, in the one place it would have
reached the paper, and it was the completeness critic rather than any of the eight lenses
that caught it.

### With noise, and why `k_e` is the interesting case

200 replicates at each of 8 seeds, monthly design, no priors, constant load. Quoted with the
spread across seeds, because an earlier version of this section quoted one:

| quantity | truth | min | median | max |
|---|---|---|---|---|
| median `b_hat` | 0 | 147.3 | 399.8 | 774.4 |
| median `theta_hat` | 3000 | 6905 | 7889 | 8646 |
| median `K_m_hat` | 100 | 27.8 | 42.8 | 69.6 |
| median `k_e_hat` | 5 | **-0.88** | 4.67 | **14.34** |
| fraction with `k_e_hat < 0` | | 0.415 | 0.490 | 0.505 |
| fraction with `theta_hat <= 0` | | 0.220 | 0.245 | 0.290 |

The sign of the fitted `k_e` is a coin toss and its median is not quotable from one seed. What
is stable is the collapse. Relative inter-quartile ranges: `A = K_m(L+b)` **0.078**, `B` 0.303,
`C = -k_e` 16.9, `b` 12.6, `theta` 4.36, `K_m` 0.857. Only the invariant `A` is pinned.

**`k_e` is exactly identified and still barely estimable.** `C = -k_e` is an exact invariant,
so at a constant load `k_e` is structurally identified with no fibre at all. Over the marker
band the design visits, `[3.51, 8.75]`, a 2.49-fold range, `cond{1/m, 1, m} = 59.9` and the
Monte Carlo IQR above is 16.9. So "exactly identified" and "usefully estimated" are as far
apart here as `rem:rankvsprec` says, and this time the gap runs the other way. The paper's
`k_e` prior is not redundant.

### Reconciling with `StructuralIdentifiability.jl`

No contradiction, and the paper already has the right reconciliation at
`06_excitation.tex:55-65`. The Julia verdict is about the transient from a generic initial
condition, where every time derivative of the output is available. On the locus `dm/dt = 0`, so
every higher derivative is `0 = 0`. `tau = 0.0336 d` against a 3.5 d half-life, a separation of
104x, is what puts the clinic on the locus. What changes is only the sentence that follows: the
obstruction is not "real but priced at `sqrt(noise)`". At a constant load it is an exact fibre
priced at whatever bound you can put on `b`; at a varying load it is closed by the data.

### Three fibres, and they are different objects

| | unknowns | observables | dimension | `k_e` along it |
|---|---|---|---|---|
| `rem:fibre` (`04_local.tex:176`) | `(theta, K_m, k_e)` | `(a, b)` of an affine tangent | 1 | **varies** |
| `prop:rank` / this one | `(theta, K_m, k_e, b)` | `(A, B, C)` on the locus | 1 | **fixed** |
| `rem:modelB` (`06_excitation.tex:306`) | `(theta, K_m, k_e, K_e)` | saturating-channel exchange | 0 (two points) | exchanged |

### Corrections the stress test forced, on this audit's own work

**C1. Claim (5)'s design and coordinates were wrong.** `chi^2/eps^4` is a second-order
quantity and is therefore **parametrisation dependent**. `prop:rank` is stated in log
parameters (`06_excitation.tex:332`, "the three design columns ... in `(ln theta, ln K_m,
L_0)`") and the information filter's columns are derivatives in those coordinates. The audit
drew the line in natural `(theta, K_m)` coordinates. Redrawn in the paper's own:

| design | n | natural `chi2/eps^4` | half-width | **log `chi2/eps^4`** | **half-width** |
|---|---|---|---|---|---|
| **two-daily x 2 y** | **365** | 6.499e-06 | 26.66 | **2.0900e-06** | **36.15** |
| weekly x 2 y | 105 | 1.908e-06 | 35.68 | 5.879e-07 | 49.30 |
| monthly x 2 y | 25 | 4.438e-07 | 50.18 | 1.436e-07 | 69.41 |

The paper prints `2.08e-6` and `+/-36.9`. The log / two-daily entry matches to 0.5%, and the
provenance is documented rather than inferred: `audit/c2b_degeneracy.py:95` runs
`np.arange(0.0, 730.0, 2.0)`, `audit/PHASE2C.md:92-94` attaches the `2.08e-6` to
`sd = 0.0251`, and `audit/PHASE2A.md:35` shows `0.0251` is the two-daily row. **So
"a design of roughly weekly size", "not the monthly one", "`+/-54.2`" and "a factor of 1.72,
not within a fifth" are all withdrawn.** On the design the number was computed for, the
paper's "the two agree to within a fifth" is right as computed: `36.86/31.63 = 1.17`.

**C2. A weaker design charge survives.** `rem:second_order` claims "every covariance this note
reports survives", and the paper points the remark at the `0.1696` entry of `rem:rankvsprec`
(`06_excitation.tex:482-485`), which is the **monthly** design. On monthly the same log
construction gives `69.41` against the prior's `31.63`, a ratio of `2.19`. The half-width runs
`36.2 / 49.3 / 69.4` across two-daily / weekly / monthly, so one number cannot defend every
covariance the note reports.

**C3. The `+/-54.2` was a quartic extrapolation, and so is the paper's `+/-36.9`.** `chi^2`
on the line is not exactly quartic: `chi^2/eps^4` drifts about 6% over four decades and about
40% by the time `eps` reaches the quoted half-width. The exact `Delta chi^2 = 3.84` crossings
are in the table above. This cuts **for** the surviving half of claim (5): the printed pair is
internally the quartic extrapolation, `(3.84/2.08e-6)^(1/4) = 36.86` against a printed `36.9`,
and a true profile does not produce a constant `chi^2/eps^4` to read a half-width off.

**C4. A bug in `audit/f5_fibre.py`.** The bisection for `theta(b) = 0` inverted its bracket
update and printed `-215.1515`, which is where `K_m` reaches 300 and where `theta` is still
1785. The correct root is `-272.7273`. Fixed with `brentq`.

**C5. One S4 number was a seed artefact.** The `k_e` median of `-0.879` is the minimum over
eight seeds, and the median across them is `+4.67`. The quotable statistic is that the sign is
a coin toss, negative in 42-51% of replicates. The `theta_hat <= 0` rate does survive the seed
sweep (0.22 to 0.29). `S4` now quotes spreads. This repository's own rule, from `PHASE2D`, is
to quote the spread, and the first version of this section broke it.

**C6. Scope, in both directions.** Too wide: "flat only to first order"
(`06_excitation.tex:66`) is a **true** and reproducible statement about `prop:rank`'s null
*vector* (local power 2.0000 to 2.0002 over a decade), and only becomes false read as a claim
about the likelihood's degeneracy, which is a curve. The same applies to `07_design.tex:141`.
Only "recovered through curvature rather than through slope" and "the recoverable displacement
scales as the square root of the noise level" are unambiguously false. Too narrow:
`rem:rankvsprec` is **not** wholly correct, as this document elsewhere concedes. Its sentence
"At any real day-to-day variation the data fill it instead and the same `0.1696` is earned
rather than assumed" is false at the 0.5% and 2% rows the same remark prints, where prior-free
gives `0.8604` and `0.2630` against `0.1696`. Separately, "Halving the assay CV buys ... a
factor of 2 along every other direction" (`06:74-75`) fails for a reason independent of P5,
because the load-record term floors `sigma_total`.

**C7. Two qualifiers the claims need.** Claim (2)'s "everything" is everything **the
quasi-steady locus** contains, and the count is generic: it needs at least three distinct
increments and `L + b != 0`. (A claimed further rank drop on the synchronous 14 d and 28 d
designs was a finite-difference artefact; built analytically at 60 digits those designs are
rank 3 with a 24-order gap.) Claim (3)'s "machine zero" is a locus statement: on the ODE record
the same quantity is `2.28e-3` along the fibre and `8.71e-5` as a true profile at `b = 200`,
four orders below any inferential threshold but not machine zero.

**C8. "Carries no information" overstates, by two channels.** Neither rescues the paper.
Under errors-in-variables with an unrecorded per-visit load fluctuation, the lever arm
`J = d(load)/dm` varies along the fibre (33.7 to 51.9 at the high-`m` end between `b = 0` and
`b = 200`) even though `m` does not, so the variance transmitted from load to marker differs
between twins. Marginalised by Gauss-Hermite quadrature this gives a genuine 95% interval in
`b` of about `[-120, +383]` (400-replicate median, 95.8% coverage), four to twelve times looper
than the prior, a variance channel rather than a mean-curvature one, and with a noise scaling
that **saturates** rather than following `sqrt(noise)`. Second, a record sampling inside the
post-dose settling window does break the fibre: nine draws at 5 to 120 minutes cut the
half-width to `+/-108` and 225 such draws reach `+/-24.2`. That design is excluded by P1f and
by `06_excitation.tex:56-65`, but replacement text should not say flatly that the combination
is never recovered.

**C9. Two more passages belong in the scope, found by the judges.**
`07_design.tex:25` and `A1_supplementary.tex:516` say that `K_m` lies inside the exact
degeneracy, so a prior on it cannot help, and that this is "a reason, not a convenience".
The stated mechanism is inverted, and the fibre is what proves it: along the fibre
`b = L(K_0/K_m - 1)`, so a `ln K_m` prior of sd `s` bounds `b` to
`[L(e^{-1.96s} - 1), L(e^{+1.96s} - 1)]`. At `s = 0.8` that is `[-255, +1225]`, and at
`s = 0.0477` it matches the load record's `+/-31.6` exactly. Numerically, with the load
prior removed and `k_e` free, a `K_m` log-prior alone gives `sd(ln theta)` of
`0.2820` (monthly) at `s = 0.8` and `0.1854` at `s = 0.2`, against a singular filter with no
`K_m` prior at all. So a `K_m` prior does help, and the empirical conclusion survives for a
different reason: one tight enough to matter is far outside anything the literature
supports. Treat it as `07_design.tex:141` was treated, replacing the mechanism clause and
keeping the conclusion.

**C10. C8 is too generous to itself.** The errors-in-variables channel it records exists
only while the load record's error is held exactly known at 5% of `L`. Profiling that scale
as a free nuisance collapses the discrimination by three orders of magnitude, and there is
no 95% crossing at any admissible `b`. So that channel is the load prior again rather than
the data, and C8's "[-120, +383]" must carry the condition.

**C11. The design pinning in C1 is good to about 2%, not 0.5%.** Walked in log coordinates
on the two-daily design, `chi^2/eps^4` runs `2.0399e-6` to `2.1172e-6` across the prior
conventions the repository uses (pure data, `c2b_degeneracy.py`'s own `0.8` priors, the
minimal footing, `k_e` free). All of them sit on the printed `2.08e-6` and none of them sits
near the natural-coordinate `6.5e-6`, so the conclusion is unaffected and the quoted
agreement should be `2%`.

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

**Nothing in the headline, and this was attacked directly rather than assumed.** On the
monthly design at load CV 10% the load-bias prior is worth **1.4%** of `sd(ln theta)`
(0.1697 with it, 0.1721 without; across seven seeds 0.1688-0.1709 with, 0.1710-0.1810
without), and 0.6% to 3.5% on every other design in `tab:designs`. The 42x excitation
(1.9390 against 0.0462), the factor-of-22 sampling rule (21.7x, and 25.3x prior-free), the
11.7% dose-precision requirement and the cohort numbers all reproduce. **No printed number
in the paper changes.**

What fails is two clauses, a diagnostic, a cross-reference and a step in one proof. Listed
here, **not edited**, pending your go-ahead. The scope is narrower than this document's
first draft had it, and wider in one place.

| file:line | what to change | why |
|---|---|---|
| `06_excitation.tex:71-75` | "**recovered through curvature rather than through slope**" and "**the recoverable displacement scales as the square root of the noise level**" | both false. The profile is flat, so nothing is recovered; and where the data do work the half-width scales linearly in `sigma`, not as its square root |
| `06_excitation.tex:74-75` | "Halving the assay CV buys a factor of `sqrt2` along this direction and a factor of `2` along every other" | the first half is false for the reason above. The second half fails **independently of P5**, because the load-record term floors `sigma_total` |
| `06_excitation.tex:497-516` (`rem:second_order`) | the diagnostic, and "`chi2 ~ eps^4` means the data do pin it, at second order" | it reads a curved exact fibre as data-pinned, so it returns the wrong answer on the problem it was written for. The `+/-36.9` is a ray, not a profile |
| `07_design.tex:141` | "because ... the gradient is third order" | the *mechanism clause* only. The statement is true of the line; the operative fact is start dependence on a flat fibre |
| `07_design.tex:417-419` (`rem:no_delta`) | "where a first-order filter reports a curvature-held direction as prior-held" | exactly inverted. The direction *is* prior-held; the `eps^4` diagnostic is what calls it curvature-held |
| `06_excitation.tex:371` | `(theta+Delta) u = L - k_e m` | drops the bias. It should be `L + b - k_e m`, so `eq:nullvec` as printed is the `b = 0` value of a `b`-dependent field. Independent of the rest, and worth fixing either way |

**Not to be changed, against this document's first draft.** `06_excitation.tex:66`, "flat
only to first order", is a true and independently reproducible statement about
`prop:rank`'s null *vector*. It is incomplete, not false, and calling it wrong was an
overreach.

**To be changed, also against this document's first draft.** `rem:rankvsprec`
(`06_excitation.tex:439-490`) was called "correct as it stands" here. It is not. Its
sentence "At any real day-to-day variation the data fill it instead and the same `0.1696`
is earned rather than assumed" is false at two of the four rows the same remark prints: at
load CV 0.5% and 2% the prior-free figures are `0.8604` and `0.2630` against `0.1696`. The
reading of the zero-CV row as the prior filling the direction *is* correct and should stay.

Two audit records carry the retracted claim and will need notes: `PHASE3_CHECKS.md:301-302`
and `PHASE2C.md:141-150`. So does the closing commentary of `audit/e5_structural.py` `q3`,
whose *computation* is correct and whose *interpretation* is the one being retracted.

---

## Questions to take back to him

1. **Is the impulsive dose model a proposal or a reading of ours?** If he believes the
   therapy removes metabolite rather than adding capacity, that is a different biology and
   worth settling first: mRNA-3927 delivers PCC mRNA, which is a capacity increment, and
   the 226x difference in integrated effect means the two models do not share a design.
   *(Three questions stood here -- on `psi`, on the missing `Delta_k m_k` term, and on the
   signs. They are withdrawn: the equation was a sketch of the method, and a sketch does not
   owe a complete term list. The `beta` reparametrisation it pointed at is kept and
   credited.)*
2. **"1 or 1/2 of the interval": which did he mean?** At `T` the design fails by a factor of
   22; at `T/2` it passes. If he meant something like "any rational ratio with a small
   denominator", the data say the opposite: small `q` is fine, and what matters is the
   `Delta` range.
3. **One `params.py` row says `theta` is "about 10% of normal" and, four lines later, "under 5% of control".** Which
   reference does he want the acceptability check written against? The two differ by a
   factor of two and they are the denominator of the ceiling constraint he asked for.
4. **Does he accept the exact fibre, and with it that the `sqrt(noise)` sentence goes?** The
   honest replacement is that at a constant load the deficiency is closed by a bound on the
   load bias, and at any real day-to-day variation it is closed by the data. Both are
   defensible and neither needs curvature.
5. **His second question, the a priori reduction bound, is untouched by this round**, and
   P3 tightens it: the linear route he proposes exists only for the straight-`alpha`
   hypothesis. Does he want that pursued next, or the design outputs first?

---

## Status

Phase 0 is complete, and it has now been adversarially stress-tested
(`audit/PHASE4_STRESS.md`). No prose has been edited, no figure regenerated and no constant
changed. Phases 1, 2 and 3 are waiting at their gates.

Committed alongside this report, because both lock facts rather than sentences:

- `tests/test_models.py::test_the_constant_load_deficiency_is_an_exact_fibre_not_a_curvature`,
  which asserts that the fibre reproduces the same levels from the same increments, that the
  misfit is machine zero along it, that its tangent at `b = 0` is `prop:rank`'s null vector,
  and that the straight line gives a quartic. Six mutations of the fibre formula were tried
  against it and all six fail it;
- `REVIEW_LOG.md` `R14`, the retraction, in the form the existing entries use.

Held, unapplied, pending your reading: `audit/PHASE4_PROPOSED_TEXT.md`, replacement wording
for the six passages listed above. The paper stays frozen.

On the repository's taxonomy this is a new failure class and worth naming: **a diagnostic
that answers a question the model cannot pose.** `audit/e5_structural.py` `q3` computes
`eps^2` correctly, from committed code; the paper quotes it accurately; the power law is real
to four decades. What is wrong is the inference, and it is wrong because the test profiles
along a straight line in a space whose relevant object is curved. No numerical check here
could have caught it, because there is no number in it to disagree with.

The audit's own failure has the same shape and is recorded as C1: `chi^2/eps^4` is
second-order and therefore parametrisation dependent, and the audit measured it in one
coordinate system while the paper had computed it in another. Both errors are a quantity
being right about its own computation and wrong about what it is a computation of, which is
`R8`, and both were caught the same way `R8` was, by two descriptions of one object failing
to agree.
