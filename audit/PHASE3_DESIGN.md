# Phase 3. The design outputs

Four scripts, four figures-or-tables, one independent re-check. **No paper edits, no merge,
no changes to existing constants or printed numbers.** `PHASE4_PROPOSED_TEXT.md` stays
unapplied.

| | what it answers | output |
|---|---|---|
| `audit/g1_phase_clock.py` | where in the dose cycle a schedule actually lands | `results/figures/fig_E1_phaseclock.{pdf,png}` |
| `audit/g2_visits_vs_sd.py` | how many visits, and what precision they buy | `results/figures/fig_E2_visits.{pdf,png}` + table |
| `audit/g3_calendar.py` | the inverse map: target phase to calendar date | `results/figures/fig_E3_calendar.{pdf,png}` + tables |
| `audit/g4_recheck.py` | the five key numbers, by routes sharing no code | table |

Shared machinery is in `audit/g_design_common.py` so the three cannot drift.

## Footing, once, for every number below

Information filter in **the paper's coordinates** `(ln theta, ln K_m, b, ln k_e)`, the
columns `src/design/run.py:324-329` differentiates in and the ones `prop:rank` is stated in
(`06_excitation.tex:332`). Load record 5%, `sigma = sigma_total(0.05, 0.05) = 21.9612`,
`L = 322.7273`, load prior half-width `1.96 x 0.05 x L`. **`k_e known` DROPS the column**
rather than giving it a flat prior, which is `_filter`'s convention and `e9`'s
(`e9_rank_vs_precision.py:49-58`). Visits lie on a half-day grid. **Every figure that
depends on a random draw is the median over 8 seeds and the min-max spread is quoted beside
it.** Nothing in any script records an expected result in advance; the one outside number,
a second reader's estimate of the monthly floor and `c`, is labelled where it is used.

---

## 1. The phase clock (`g1`, fig_E1)

Visits folded onto the dose cycle: angle `2 pi (t mod 14)/14`, dose at 12 o'clock, radius
increasing with visit order, colour `Delta/theta`, settling window shaded. Four panels at
`n = 25`.

| schedule | phases | `Delta/theta` range | horizon | sd(ln theta), median of 8 seeds (min-max) |
|---|---|---|---|---|
| 14 d (the dose interval) | **1** | 0.1333 | 336.0 d | **4.7760** (3.5141-5.4718) |
| 14.5 d | 25 | 1.8019 | 348.0 d | 0.1739 (0.1732-0.1757) |
| monthly (30.4 d) | 25 | 1.8438 | 729.6 d | 0.1698 (0.1688-0.1709) |
| phase-spread | 25 | 1.8528 | 737.5 d | 0.1707 (0.1698-0.1723) |

*load CV 10%, `k_e` log-prior 0.8.*

The figure makes the design rule visible in one look: the 14 d schedule is a single spoke,
the other three fill the dial. Two things worth stating that the dial shows and prose does
not.

**The settling window is not a design constraint.** It is **1.30 h** at a 1% tolerance and
**0.91 h** at the repository's own 5%, which is `1.39` degrees of the dial. None of the four
schedules has to be moved to avoid it. That is a finding, not a formality: the specification
carries the exclusion, and it now has a price, which is nothing.

**The 14 d schedule's only excitation is the dose train still filling up.** Its
`Delta/theta` range is 0.1333, which is exactly `(6400 - 6000)/3000`. That is the approach of
the geometric sum to its steady state, not any spread of phase. That is why its sd is
around 5 rather than exactly singular, and why the sd falls so slowly with `n`: 5.6347 at
`n = 4` to 4.4503 at `n = 60` (median of 8 seeds, prior 0.8, CV 10%), a factor of 1.27 for
fifteen times the visits. A design resting on that is resting on a transient.

---

## 2. Visits against precision (`g2`, fig_E2)

`sd(ln theta_hat)` against `n` for each schedule, at `k_e` prior 0.8, 0.3, 0.05 and known,
at load CV 0, 2 and 10%, with the 0.20 target marked, fitted to
`sd(n) = sqrt(floor^2 + c^2/n)`.

The fit only describes a schedule **once that schedule has covered the dose cycle**. Below
that, sd falls for two reasons at once, more visits *and* more increment range, and the
fit absorbs both. So the table reports `n*`, the first `n` at which the `Delta/theta` range
reaches 95% of its own maximum, and fits twice.

### At the paper's own footing, `k_e` prior 0.8

| schedule | `n*` | floor | `c` | rms | visits to 0.20, **observed** | implied by the fit |
|---|---|---|---|---|---|---|
| 14 d | 4 | 4.51 | 7.06 | 6.1e-2 | **never** (>60) | never |
| 14.5 d | 30 | **0.1518** | **0.3763** | 9.3e-4 | **20** | 8 |
| monthly (30.4 d) | 36 | **0.1504** | **0.3947** | 1.0e-4 | **11** | 9 |
| phase-spread | 60 | **0.1496** | **0.4128** | 9.7e-4 | **10** | 10 |

*load CV 10%; the CV 0% and 2% rows differ in the fourth decimal. floor and `c` are the
saturated fit where it exists; phase-spread does not saturate inside the sweep, and its
all-n fit has rms 9.7e-4, so the form holds there anyway.*

**Three results.**

**(a) The `k_e` prior sets the floor, and nothing else does.** Monthly's floor goes
0.1504 (prior 0.8) -> 0.0590 (0.3) -> 0.0204 (0.05) -> 0.0181 (known), while `c` stays at
0.394-0.396 throughout. Visits buy `c/sqrt(n)`; they cannot buy the floor. Tightening the
`k_e` prior from 0.8 to 0.3 takes monthly from 11 visits to 6.

**(b) Load variability moves neither floor nor `c`.** At every schedule and prior, CV 0%,
2% and 10% agree to the fourth decimal. This is the paper's own rank-versus-precision
distinction landing on the design: variation is what restores *rank*, and once rank is
restored more of it buys almost no *precision*.

**(c) Observed and fit-implied visit counts disagree, and observed is the honest one.**
The saturated fit answers "how many visits if the first `n` had the coverage the large-`n`
design has", and for a fixed interval they do not. 14.5 d needs 20 visits, not the 8 its
fit implies. The one schedule where the two agree is **phase-spread**, which is the whole
point of it: coverage is built in at every `n` rather than accumulated.

### Against the second reader's estimate

Their monthly floor 0.145 and `c` 0.45, at prior 0.8: **floor confirmed** (0.1503, 1.04x),
**`c` 12% high** (0.3941, 0.88x). Their `c` matches the all-n fit (0.482, 1.07x), but that
fit has rms 1.3e-2 against the saturated fit's 1.2e-4, so it is describing cycle coverage,
not sample size.

### The caveat the load-variation benefit rests on

Every gain attributed to load CV assumes **the diet record captures real day-to-day
variation in intake**. It is the recorded `L_k` that enters the design columns, so a record
imputed from a prescription, or averaged to a weekly mean, carries no variation into the
fit however much the patient's intake actually varies, and turns the CV 10% rows into the
CV 0% rows. The paper already asks for this as a row of the specification
(`06_excitation.tex`, after `rem:rankvsprec`); these figures price it.

---

## 3. The inverse map, target phase to calendar date (`g3`, fig_E3)

g1 and g2 argue in phase. A clinic books in dates. Running the map the other way turns up
the sharpest practical finding of Phase 3, and it is arithmetic:

> **The dose interval is 14 days, which is exactly two weeks. So a visit's phase fixes its
> weekday, permanently.** A target phase is either always bookable or never bookable. There
> is no rescheduling it into the working week.

Verified against `datetime` over 60 cycles for all 28 half-day phases and both dose timings:
**zero phases drift** (`g4` check 1).

### What a weekday clinic can reach

| dose given | bookable slots | `Delta/theta` reachable | range | of the usable dial | reaches the trough? |
|---|---|---|---|---|---|
| morning | 19 of 27 | 0.2188 - 1.9322 | 1.7135 | 96.0% | no |
| afternoon | 19 of 27 | **0.1472** - 1.9322 | **1.7850** | **100%** | **yes** |

*at steady state; the usable dial is the half-day grid less the on-the-dose slot.*

**The time of day the dose is given is a design variable.** Both timings leave the same
number of slots. They are not the same slots: an afternoon dose puts phase 13.5, the
trough and one end of the quantity `f4` found binds, on a Monday morning, and a morning dose
loses it to the weekend. It costs nothing and it is worth asking for.

### One booked design, 12 visits

Targets chosen from the 19 reachable phases, taken nearest-reachable first, minimum gap
30.4 d, first dose Monday (an arbitrary one; only the weekday matters):

| # | phase | date | day | when | `Delta/theta` | | # | phase | date | day | when | `Delta/theta` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.5 | 2026-01-06 | Tue | am | 1.8114 | | 7 | 7.5 | 2026-07-14 | Tue | am | 0.4831 |
| 2 | 3 | 2026-02-05 | Thu | pm | 1.1774 | | 8 | 10 | 2026-08-13 | Thu | pm | 0.2944 |
| 3 | 6.5 | 2026-03-09 | Mon | am | 0.5888 | | 9 | 1.5 | 2026-09-16 | Wed | am | 1.5851 |
| 4 | 9 | 2026-04-08 | Wed | pm | 0.3589 | | 10 | 4 | 2026-10-16 | Fri | pm | 0.9661 |
| 5 | 13.5 | 2026-05-11 | Mon | am | 0.1472 | | 11 | 8 | 2026-11-17 | Tue | pm | 0.4375 |
| 6 | 2 | 2026-06-10 | Wed | pm | 1.4356 | | 12 | 10.5 | 2026-12-18 | Fri | am | 0.2667 |

Horizon 346.5 d (11.4 months), 12 distinct phases, `Delta/theta` range 1.6642, every visit
on a clinic day outside the settling window. `sd(ln theta)` **0.1929** at CV 0%, 0.1932
(8 seeds 0.1890-0.1958) at CV 10%.

### Phase-spread against fixed intervals, at equal `n`

The fair comparator is not the fixed schedule as written but the fixed schedule **as a
clinic would actually run it**: each off-clinic visit moved to the nearest bookable slot.
Moving one moves its phase, so it is a cost, not a formality.

At `n = 12`, load CV 10%, `k_e` prior 0.8, median of 8 seeds:

| schedule | horizon | phases | range | off-clinic | sd (CV 0%) | sd (CV 10%, min-max) |
|---|---|---|---|---|---|---|
| 14 d | 154.0 d | 1 | 0.1333 | **12/12** | 5.1569 | 5.0155 (3.67-5.84) |
| ... snapped | 154.5 d | 1 | 0.1208 | 0/12 | 5.0027 | 4.8443 (3.49-5.67) |
| 14.5 d | 159.5 d | 12 | 1.2822 | 3/12 | 0.2899 | 0.2879 (0.280-0.296) |
| ... snapped | 158.5 d | 9 | 1.0496 | 0/12 | 0.3239 | 0.3205 (0.311-0.332) |
| monthly | 334.4 d | 12 | 1.8170 | **4/12** | 0.1890 | 0.1888 (0.186-0.192) |
| ... snapped | 333.5 d | 9 | 1.7135 | 0/12 | 0.1927 | 0.1927 (0.189-0.195) |
| phase-spread, whole dial | 347.5 d | 12 | 1.8375 | 3/12 | 0.1909 | 0.1913 (0.187-0.194) |
| ... snapped | 347.5 d | 10 | 1.5927 | 0/12 | 0.1938 | 0.1943 (0.190-0.196) |
| **phase-spread, clinic (pm dose)** | 346.5 d | **12** | 1.6642 | **0/12** | **0.1929** | **0.1932** (0.189-0.196) |
| phase-spread, clinic (am dose) | 347.5 d | 12 | 1.5927 | 0/12 | 0.1952 | 0.1951 (0.193-0.198) |

**Five things this settles.**

1. **A third of a monthly schedule falls on a weekend** (4 of 12, and 8 of 25 at
   `n = 25`). The paper's comparator schedules were never checked against a calendar.
2. **Snapping is cheap, about 2%.** Monthly 0.1890 to 0.1927. The weekday constraint does
   not threaten the design.
3. **It is cheap for a reason worth noticing.** Snapping monthly costs it **3 of its 12
   distinct phases** (-25%) but only 5.7% of its `Delta/theta` range, and sd moves 2%,
   tracking the range rather than the phase count. That is an independent confirmation of
   `f4`'s finding, arrived at by a route that was not looking for it.
4. **The clinic phase-spread ties the snapped monthly** (0.1929 against 0.1927) at
   essentially the same horizon (346.5 d against 333.5 d), while keeping all 12 distinct
   phases where the snapped monthly keeps 9. What it gains is not speed. It is that the
   design is *specified in phases*, so it reaches its target by construction rather than by
   luck, and it never asks the clinic for a Saturday.
5. **Nothing rescues the 14 d schedule.** Snapped, all 12 visits become bookable and sd
   stays at 4.84. The settling objection to "fortnightly, on the dose" is real but
   secondary. The rank objection is the one doing the work.

---

## 4. The independent re-check (`g4`)

Five numbers, each recomputed by machinery that shares no code with the script that produced
it. **It found one error and one misreport, both now fixed.**

| # | what | route | result |
|---|---|---|---|
| 1 | phase fixes weekday | `datetime`, 60 cycles, not mod-7 arithmetic | 0 of 56 phases drift; slot counts 19/19 both timings |
| 2 | the dose train | closed-form geometric sum vs the 200-term loop | worst relative difference **1.9e-16** |
| 3 | settling window | the balance integrated here, own quasi-steady root | 1.299 h vs the repo's 1.299 h (ratio 1.0002); 0.909 vs 0.907 at 5% |
| 4 | `sd(ln theta)` | curvature of a **profiled nonlinear** objective, no analytic column anywhere | ratios **0.9984-0.9999** across four designs |
| 5 | floor and `c` | exact two-point solve, no curve fit | four pairs give monthly floor 0.1496-0.1507, `c` 0.390-0.407 |

**Check 4 is the one that matters.** `sd_ln_theta` builds analytic derivative columns and
inverts an information matrix; the profile writes the same posterior out in the parameters
themselves, profiles three of them out numerically, and reads sd from the curvature. A wrong
column or a misplaced prior could not cancel between the two. They agree to within 0.16%.

### What the re-check caught

**(i) An error, now fixed.** `g3` Part A evaluated `Delta/theta` **on the first dose cycle
instead of at steady state**, understating every phase's increment by about 7% (peak read
2.0000 where it is 2.1333; the usable range read 1.6734 where it is 1.7850). Check 2 caught
it against the closed form. Parts B, C and the figures were never affected, because they read
actual visit times, where the transient is real. Fixed in `inc_at_phase`; the *conclusions*
of Part A, including "an afternoon dose reaches 100% of the usable range", survive unchanged
because both endpoints moved together.

**(ii) A misreport, now fixed.** `g2`'s `n at 0.20` column read **"never" for phase-spread**
-- an artefact of the saturated fit returning nan when `n*` lands at the end of the sweep,
not a finding. The column now falls back to the all-n fit and, more importantly, reports the
**observed** crossing read straight off the sweep beside the fit-implied one. The sweep was
also densified to integer resolution through `n = 14`, where the crossings are. A
`saturation_n` docstring claim that monthly "is covered almost at once" was wrong, `n*` is
36, and has been replaced with what the function measures.

---

## What I would report, and what I would not

- **Book the clinic phase-spread with an afternoon dose.** 12 visits, 11.4 months, every
  visit bookable, `sd(ln theta) = 0.193`. The calendar is in `g3` Part B.
- **Quote the saturated floor and `c`, not the all-n fit**, monthly `floor 0.150`,
  `c 0.394` at prior 0.8. And quote the **observed** visit count, not the fit-implied one.
- **Do not quote a precision without the `k_e` prior beside it.** It is the floor, and the
  floor is most of the answer at the visit counts a trial can afford.
- **Do not claim load variability buys precision.** It buys rank. Once rank is there, CV 0%
  and CV 10% are the same number to four decimals.

## Status

`g1`-`g4` all exit 0. Three figures written to `results/figures/`. No paper file touched, no
constant changed, no printed number in the paper altered, `PHASE4_PROPOSED_TEXT.md` still
unapplied. The PI's clarification is recorded at the head of `PHASE4_PI.md` P3 and the three
questions it retires are withdrawn from the questions list.

**Stopping here for go-ahead before Phase 2 item 4 (physical admissibility).**
