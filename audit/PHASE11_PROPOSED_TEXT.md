# Phase 11 proposed text. Three passages, NONE APPLIED

`OPEN_QUESTIONS.md` section 0: *the paper and supplement are frozen for prose. No wording
change until it is reported and approved.* This file reports. Nothing in `paper/` is
touched, and `audit/PHASE4_PROPOSED_TEXT.md` and `audit/PHASE10_PROPOSED_TEXT.md` remain
the other two unapplied sets.

Each passage names its file and line, quotes the current text, gives the proposed text, and
says which script owns each number so that a wiring pass can reach it.

---

## 1. `paper/sections/11_scope.tex`, the closing paragraph of `sec:scopeprec`

**This is the one the work asks for.** The verdict survives; the stated mechanism does not.

### Current

> This is the severely affected patient. The population in which the dosing question
> is most urgent is the population in which the record cannot answer it. This is not a
> statement about sampling effort: no schedule, no excursion and no precision on the dose
> recovers a quantity the observable has stopped responding to. The therapy channel of
> Section~\ref{sec:excitation} does not rescue it either: it supplies excitation in
> $\theta$ against $K_m$, not sensitivity of $m$ to $\theta$.

### What is wrong with it

Nothing in the verdict. The difficulty is the clause *"a quantity the observable has
stopped responding to"*, which is true of **the marker** and false as a statement about
observables in general. The section itself already says where one would have to act, *"on
$\sigma$, or on an observable that still responds to $\theta$ where the marker has
stopped"*, and there is now a named candidate with a number against it: the breath test's
oxidised fraction responds to capacity **better** as capacity falls, by almost exactly the
factor the marker loses.

| $\theta$ | marker, $\mathrm{d}\ln m/\mathrm{d}\ln\theta$ | breath, $\mathrm{d}\ln F/\mathrm{d}\ln\theta$ |
|---|---|---|
| 3000 | -0.9154 | 0.1549 |
| 600 | -0.5367 | 0.5277 |
| 300 | -0.3210 | 0.7083 |
| 150 | -0.1729 | 0.8370 |

The two cross at about 20% of nominal. Leaving the sentence as it stands invites a reader
to conclude that the whole class is closed, which this work finds is not the reason it
fails.

### Proposed, replacing the third sentence and adding one

> This is the severely affected patient. The population in which the dosing question
> is most urgent is the population in which the record cannot answer it. This is not a
> statement about sampling effort: no schedule, no excursion and no precision on the dose
> recovers a quantity \emph{the marker} has stopped responding to. The therapy channel of
> Section~\ref{sec:excitation} does not rescue it either: it supplies excitation in
> $\theta$ against $K_m$, not sensitivity of $m$ to $\theta$. \textbf{Nor is it a statement
> about every observable.} The oxidised fraction a breath test reads responds to capacity
> \emph{better} as capacity falls, its elasticity running $\num{breath_el_nominal}$ at
> nominal against the marker's $\num{marker_el_nominal}$ and
> $\num{breath_el_disease}$ against $\num{marker_el_disease}$ at a tenth of nominal, the
> two crossing near a fifth of nominal. What closes that route is not insensitivity but
> arithmetic downstream of it, and Remark~\ref{rem:breathcount} gives the count.

### Numbers, and their owners

| macro | value | owner |
|---|---|---|
| `breath_el_nominal` | 0.1549 | `audit/u2_breath_functional.py` |
| `marker_el_nominal` | -0.9154 | `audit/u2_breath_functional.py` |
| `breath_el_disease` | 0.7083 | `audit/u2_breath_functional.py` |
| `marker_el_disease` | -0.3210 | `audit/u2_breath_functional.py` |

*None of these is currently a generated macro. Wiring them is a `q22`-style pass and is not
part of this proposal.*

---

## 2. `paper/sections/11_scope.tex`, a new remark after the six closed routes

The six routes are all routes to precision on $\theta$ **from the marker record**. A
seventh, from a second observable, is now priced and belongs beside them rather than
inside them.

### Proposed, new

> \begin{remark}[A seventh route, from a second observable, and what it costs]
> \label{rem:breathcount}
> The six routes above are routes to precision on $\theta$ from the marker record. A
> breath test is a route from a different observable, and it is the one the closing
> paragraph of this section points at. Three things have to hold for it and the third is
> what decides it.
>
> It measures a \emph{fraction}, not a flux. Under \eqref{eq:local} the oxidised fraction
> of a labelled bolus is $F = \theta/(\theta + k_e(K_m+m))$, so a reading is
> $y = c\,F(m)$ with $c$ an unknown subject-specific recovery. It is \emph{not}
> proportional to $\theta\,g(m)$.
>
> The unknown recovery cancels in a ratio, which is what makes the test usable at all, and
> a ratio at two levels then determines one number: $u = \theta/k_e + K_m$, not $K_m$.
> Capacity follows only with the record's own coefficients, $K_m = A/(B+Cu)$ and
> $\theta = C(u-K_m)$.
>
> And the count is the answer. Joint Fisher information from the recommended monthly record
> and $k$ breath readings placed at the dose train's trough and peak: at nominal capacity
> the record already meets the gate and $k=0$; at a fifth of nominal $k = \num{breath_k_600}$;
> at a tenth $k = \num{breath_k_300}$, all at a $5\%$ repeatability per reading. At the
> $14\%$ that is the only figure published the last of those becomes
> $\num{breath_k_300_pub}$. \textbf{A single reading contributes exactly nothing}, the free
> recovery factor absorbing it, and lengthening the marker record from thirteen draws to
> two hundred moves the count only from $78$ to $47$.
>
> So the route is not closed by insensitivity, which is the closing paragraph's point, and
> it is closed by arithmetic, which is this one's.
> \end{remark}

### Numbers, and their owners

| macro | value | owner |
|---|---|---|
| `breath_k_600` | 34 | `audit/u3_joint_fisher.py` |
| `breath_k_300` | 70 | `audit/u3_joint_fisher.py` |
| `breath_k_300_pub` | 547 | `audit/u3_joint_fisher.py` |

### Where it interacts with what is already written

`09_validation.tex` `sec:why_not_breath` gives three reasons the breath test is a reference
standard and not a replacement: capacity is not margin, the test is not where the patients
are, and scale transfer is not established. **All three survive and none is quantitative.**
This remark adds a fourth that is, and the fourth is independent of the other three: even
granting the transfer, and granting that the instrument exists, the count is what it is.

Note also that `sec:why_not_breath`'s framing, *"the breath test measures how fast a patient
oxidises a labelled bolus"*, is the one this work had to make precise. A rate and a fraction
are different functionals and they determine different combinations of the parameters. No
change is proposed to that sentence, but a wiring pass should know the two are not
interchangeable.

---

## 3. `paper/sections/08_specification.tex`, a sentence on the span of the Tier 1 record

Tier 1's requirement is 118 monthly draws, which is 9.82 years, and nothing in the
specification says that its own constancy assumption is being asked to hold that long.
`04_local.tex` warrants constancy **over the afebrile range**.

### Current, at `:126`

> and reaches the gate only at $\num{tier1_visits_at_pub}$ monthly draws, which is nearly
> ten years of follow-up, against eleven with the increment known exactly.

### Proposed, appending one sentence

> and reaches the gate only at $\num{tier1_visits_at_pub}$ monthly draws, which is nearly
> ten years of follow-up, against eleven with the increment known exactly.
> \textbf{That span is itself an assumption and the paper has not priced it elsewhere.}
> Section~\ref{sec:local} warrants $\theta$ constant over the afebrile range; this row
> spends a decade of it, and because the count is by construction the one that lands
> exactly on the gate, it carries no margin: a log-linear drift in capacity of
> $\num{drift_cross_pct}\%$ a year is enough for no number of monthly draws to reach
> $0.20$. Carrying a drift term costs between $\num{drift_premium_lo}\%$ and
> $\num{drift_premium_hi}\%$ when there is no drift and removes the exposure, so a record
> of this length should carry one.

### Numbers, and their owners

| macro | value | owner |
|---|---|---|
| `drift_cross_pct` | 0.92 | `audit/u7_drift.py` |
| `drift_premium_lo` | 0.3 | `audit/u7_drift.py` |
| `drift_premium_hi` | 1.6 | `audit/u7_drift.py` |

**Caution on `drift_cross_pct`.** `audit/u8_recheck_drift.py` finds a smaller bias than
`u7` at the smallest drift rates (-0.0303 against -0.0458 at $r = 0.01$), which would move
the crossing up. If this passage is adopted it should say *of order one per cent a year*
rather than quote 0.92, unless the gap is closed first. Recorded here rather than rounded
away.

---

## 4. What is deliberately NOT proposed

- **No change to any of the four claims, the specification table, or any existing number.**
  Nothing in Phase 11 contradicts a figure the paper states; the two places it touches
  existing prose both add a qualification rather than correct an error.
- **No change to `sec:why_not_breath`'s three reasons.** They are unaffected, as section 2
  above says.
- **No change to the `0.20` gate, the `k_e` prior, or any constant in
  `src/common/params.py`.** `OPEN_QUESTIONS.md` section 0: no constant moves without an
  explicit decision, and nothing here needs one to move.
- **No wiring.** Every macro named above would have to be generated before the passages
  could be built, and that is a separate pass with its own guard, exactly as
  `audit/q22_wire_macros.py` did for the twenty-two conversions already made.
