# Phase 2, item 4. Physical admissibility

Two scripts and six tests. **No paper edits, no merge, no changes to existing constants or
printed numbers.** `PHASE4_PROPOSED_TEXT.md` stays unapplied.

| | what it does | output |
|---|---|---|
| `audit/h1_admissibility.py` | the eight checks, with a pass/fail table | table |
| `audit/h2_recheck.py` | four key numbers, in closed form or exact rationals | table |
| `tests/test_models.py` | six tests, mutation-checked six ways | `pytest` |

**Footing.** Local model constants as `src/common/params.py` carries them: `theta = 3000`,
`K_m = 100`, `k_e = 5`, `D = 6000`, `T = 14 d`, `lambda = ln2/3.5`, `L = load(0) = 322.7273`,
load record 5%. `Delta_peak = D/(1 - e^{-lambda T}) = 6400` exactly, since `T`/half-life is
exactly 4. The fibre of P5 is used in its own coordinates. Nothing here depends on a random
draw, so no seed spread arises.

---

## The table

**19 pass, 1 fail, 20 checks.** Every admissibility constraint proper holds, at both
readings of the capacity reference, and the recommended Phase 3 design satisfies all of them
as constructed.

| # | check | verdict | evidence |
|---|---|---|---|
| 1 | every rate positive | pass | `theta` 3000, `K_m` 100, `k_e` 5, `D` 6000, `lambda` 0.198042, `L` 322.727 |
| 2 | `Delta(t) >= 0` on 60 cycles | pass | min 375.104, max 6400 |
| 3 | `Delta(t) <=` its steady peak | pass | max 6400.000000 against 6400.000000 |
| 4 | `Delta` returns to its trough each cycle | pass | 400.000000 at the end of cycle 50; peak/trough 16.0000 = `2^(T/half-life)` |
| 5 | the balance has exactly one positive root | pass | product of roots `-L K_m/k_e = -6454.5455 < 0` at every `Delta` |
| 6 | the closed form solves the balance | pass | worst residual 1.24e-9 against `L K_m = 32272.7` |
| 7 | `m(t) > 0` across the cycle | pass | 3.3638 to 8.9715 |
| 8 | `L + b > 0` at `b = 0` | pass | `L = 322.7273` |
| 9 | `L + b > 0` across the load prior | pass | `b` in [-31.6273, 31.6273], `L + b` in [291.1000, 354.3545] |
| 10 | `theta + Delta_peak <= normal`, **both references** | pass | used 0.3133 at 10%, 0.1567 at 5% |
| 11 | physics bounds `b` on both sides | pass | [-272.7, 20107.9] at 10%, [-272.7, 50103.2] at 5% |
| 12 | the load prior is tighter than admissibility | pass | 63.3 wide against 20380.6 and 50375.9 |
| 13 | the booked design is outside the settling window | pass | 0 of 12 inside |
| 14 | every settling collision is removable by a half-day offset | pass | 6 of 6 schedules clean at `t0 = 0` or `t0 = 0.5 d` |
| 15 | `D` makes the marker "fall by about half" | **fail** | 62.25% to 66.36% on every reading |

*(Checks 1 and 8-9 expand to several rows in the script's own table; the counts above are
its rows, not this table's.)*

---

## The capacity ceiling, with the reference as a parameter

**The reference is not chosen here.** `src/common/params.py` states the residual capacity two
ways four lines apart, on the parent `theta` row: *"about 10% of normal"* and *"residual PCC
activity is under 5% of control"*. Both are carried through.

A second ambiguity has to be named first. The fraction is stated on the **parent** row
(`theta = 300`); `D` and `Delta` are defined in the **local** block (`theta = 3000`).
Evaluating the ceiling at all means carrying a fraction across models. That is an assumption,
and it is made explicit rather than silently.

| model | `theta` | `Delta_peak` | `theta + peak` | fraction | normal | used | headroom |
|---|---|---|---|---|---|---|---|
| local | 3000 | 6400.00 | 9400.00 | 10% | 30000 | **0.3133** | 0.6867 |
| local | 3000 | 6400.00 | 9400.00 | 5% | 60000 | **0.1567** | 0.8433 |
| parent | 300 | 640.00 | 940.00 | 10% | 3000 | **0.3133** | 0.6867 |
| parent | 300 | 640.00 | 940.00 | 5% | 6000 | **0.1567** | 0.8433 |

**The ceiling is not breached on either reading, and the cross-model question does not change
the answer.** `h2` check 1 shows why symbolically: the used fraction is
`f (1 + Delta_peak/theta)`, so `theta` cancels, and `Delta_peak/theta = 32/15 = 2.1333` is the
same in both blocks once `Delta` scales with `theta`. In exact rationals the used fraction is
**47/150** at 10% and **47/300** at 5%.

So the 5%/10% question changes the answer by exactly a factor of two, which is the ratio of
the two stated fractions and nothing else. **What it costs is the headroom**, 0.687 against
0.843, and headroom is what a bound on the unidentified direction is spent against.

---

## Bounds on `b`, and what they are worth

Along P5's exact fibre the data cannot separate `(theta, K_m, b)`, so admissibility is the
only thing that bounds `b`. It does bound it, on both sides.

| constraint | binds at |
|---|---|
| `K_m(b) > 0` and `L + b > 0` | `b = -L = -322.7273` |
| `theta(b) > 0` | `b = -3000/11 = -272.7273` (exact; this is the binding lower bound) |
| `theta(b) + Delta_peak <= normal` at 10% | `b = 20107.8981` |
| `theta(b) + Delta_peak <= normal` at 5% | `b = 50103.2000` |

`h2` check 2 re-derives all four in closed form: setting `theta(b) = c` and clearing the
denominator gives the exact quadratic `b^2 + b(theta_0 + L + k_e K_0 - c) + L(theta_0 - c) = 0`,
whose roots agree with `h1`'s `brentq` to six figures.

| reference | admissible `b` | width | against the load record's 95% interval |
|---|---|---|---|
| 10% | [-272.7273, 20107.8981] | 20380.63 | **322x** |
| 5% | [-272.7273, 50103.2000] | 50375.93 | **796x** |

*The load record's own 95% interval is [-31.6273, 31.6273], width 63.2545.*

**Two conclusions, and they matter for Phase 3 section 4(a).**

Admissibility **does** close the fibre, so the constant-load case is not unbounded in the way
"exactly singular" suggests. But the interval it leaves is 322x to 796x the one the 5% load
record assumes, so **on these numbers the record, not physics, is carrying the
identification.** Replacing the prior with the admissibility bound would not rescue a
constant-load design; it would only make its width finite.

And this is the concrete cost of not settling the reference: **the pass/fail answers do not
move, but the admissible range of the unidentified direction does, by 2.5x.** A bound used as
a yes-or-no does not need the reference settled. A bound used quantitatively does.

---

## The settling exclusion, and a correction to Phase 3

The check found two things I had wrong, both in this audit's own work rather than in the
repository.

**(i) The collision is never structural.** Every schedule in `g1` and `g2` starts *on* the
dose, at `t = 0`, so its first draw is at phase 0 by construction. A half-day offset clears
all six, **including the two fortnightly ones**: phase 0 to phase 0.5 is 12 h against a 1.30 h
window.

So the settling exclusion is **not a second objection to "fortnightly, on the dose"**. An
offset start clears the window and leaves the rank deficiency exactly where it was. The rank
objection is the whole of the case against that schedule, which is worth knowing, because a
reader offered two objections may discount both.

**(ii) The window is negligible in time and not in value.** Phase 3's `fig_E1` said excluding
it "costs a schedule almost nothing". That is true of the **angle** and false of the
**value**, because the phase it removes is phase 0, the **peak**, which is one endpoint of the
`Delta/theta` range `f4` found binds.

| | cost |
|---|---|
| of the dial | **0.3865%** (1.3913 degrees) |
| of the reachable `Delta/theta` range | **10.1264%** (1.986121 down to 1.784999) |

The exclusion is **26x more expensive in value than in time**. `h2` check 4 computes both
from the closed-form geometric sum, calling `delta()` nowhere. `fig_E1`'s footnote and
`PHASE3_DESIGN.md` section 1 are corrected.

The Phase 3 recommended design already pays this correctly: its most peak-ward phase is 0.5,
worth 1.9322, not phase 0's 2.1333. It is at the constrained optimum, not below it.

---

## The `D` calibration: checked, not edited

`src/common/params.py`'s `D` row states `D = 6000` was chosen so the marker *"falls by about
half at peak effect"*. **This is reported and not edited.** The claim needs two things named
that the row does not name: fall *from* what, and *which* peak.

| reading | from | to | fall |
|---|---|---|---|
| first dose `D = 6000`, from untreated | 10.0000 | 3.5099 | **64.90%** |
| steady peak `6400`, from untreated | 10.0000 | 3.3638 | **66.36%** |
| steady peak, from the pre-dose trough | 8.9106 | 3.3638 | **62.25%** |

`h2` check 3 recomputes these in exact rationals and they agree, so it is not a rounding
artefact. **On every reading the fall is larger than half**, 62% to 66%, and the smallest is
1.24x a half.

Two things stop this being serious. The row's own escape clause is the honest reading:
*"Only the excursion it produces is used"*. Nothing in the design depends on the fall being
exactly half, and the excursion is what `f4` and Phase 3 measure, as `Delta/theta` running
2.1333 at peak to 0.1333 at trough. And the "from untreated" readings are anchored on a
definition rather than a measurement: `m(Delta = 0) = 10.000000` exactly, because `L` was
defined as the clearance at `m_bar = 10`.

So the row is an order-of-magnitude justification whose order is right and whose word is
wrong by a factor of about 1.3. **It is an existing printed justification and it is the PI's
to keep or reword.**

---

## Tests

Six added to `tests/test_models.py`, one per constraint:

- `test_every_rate_in_every_parameter_set_is_positive`
- `test_the_increment_is_non_negative_and_never_exceeds_its_steady_peak`
- `test_the_balance_has_exactly_one_positive_root_at_every_increment`
- `test_the_total_load_stays_positive_across_the_load_prior`
- `test_the_capacity_ceiling_holds_at_both_readings_of_the_reference` (parametrised over the
  reference, so neither reading is chosen in the test either)
- `test_admissibility_bounds_the_constant_load_fibre_on_both_sides`

**Mutation-checked six ways, all six killed:** `D` raised tenfold (ceiling); `Delta` scaled
1% past its geometric limit; `quasi_steady_level` switched to the negative branch; the load
record widened to 80% so `L + b` can go negative; the therapeutic half-life stretched to 60
days so `Delta` accumulates past the ceiling; and the fibre map changed to
`K_m(b) = K_0 L/(L + 2b)`.

---

## Status

`h1` and `h2` exit 0. **84 tests pass**: 77 before, plus 6 new tests, which run as 7 cases because the ceiling test is parametrised over both references. No paper file touched, no
constant changed, no printed number in the paper altered, `PHASE4_PROPOSED_TEXT.md` still
unapplied.

**Two open questions for the PI**, both already on the `PHASE4_PI.md` list:

1. **Which capacity reference?** It does not change any pass or fail. It changes the
   admissible range of the unidentified direction by 2.5x, so it has to be settled before
   that range is quoted.
2. **`D` "falls by about half"**: keep as an order-of-magnitude statement, or reword to the
   62-66% the constant actually produces?

**Stopping here, as instructed.**
