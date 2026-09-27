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

### The breath group's structure, settled by the author and not yet wired

The seven breath macros are **deferred rather than guessed** (R143): their values are
row-selections from tables indexed by breath CV, and naming the row is a choice rather than a
lookup. **The author has now made that choice, and it is the structure the rest of the paper
already uses:**

> **Quote the PUBLISHED repeatability as the headline and $5\%$ as the requirement.** The sourced
> instrument is what a reader must plan against, so its figure leads; the tighter one belongs
> beside it as what a dedicated assay would have to reach.

So `breath_k_*_pub` is the headline in each sentence and `breath_k_*` at $5\%$ is the requirement
stated next to it, rather than the reverse. **The values are still not wired**, because the
published repeatability itself rests on R4, whose own width the sweep flagged, and the rows cannot
be promoted until which published figure is meant is settled.

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

---

## 5. `paper/sections/11_scope.tex`: the scope limit as a requirement, not a verdict

**Round 2's principal proposal.** The section currently states a verdict a reader can only
accept. The same computation states a requirement a reader can cost, and the requirement is
the more useful object because it names what would have to be true.

### Current, the opening of `sec:scopeprec`'s closing argument

> \textbf{The disease sits a factor of $\num{disease_factor}$ below the lower edge of the
> band the design works in.}

### Proposed, adding a sentence after it

> \textbf{The disease sits a factor of $\num{disease_factor}$ below the lower edge of the
> band the design works in.} \textbf{That edge is a function of how well the alternative
> clearance is known, and saying so turns the limit into a requirement.} On the Tier~1
> record the gate is met at a tenth of nominal if and only if $k_e$ is known to
> $\num{ke_frontier_300_pct}\%$, and at a fifth of nominal to
> $\num{ke_frontier_600_pct}\%$, against the $\num{ke_prior_pct}\%$ the literature
> supplies; at a twentieth of nominal no precision on $k_e$ suffices, because the
> increment's own prior already floors the design above the gate there. So the limit is
> two different statements at two capacities: a \emph{measurement} requirement between a
> fifth and a tenth of nominal, and a \emph{structural} one below it.

### Numbers, and their owners

| macro | value | owner |
|---|---|---|
| `ke_frontier_300_pct` | 6.3 | `audit/u9_ke_frontier.py` |
| `ke_frontier_600_pct` | 15.7 | `audit/u9_ke_frontier.py` |
| `ke_prior_pct` | 80 | `src/common/params.py`, `KE_PRIOR_LOG` |

---

## 6. `paper/sections/08_specification.tex`: 118 carries no margin, promoted from an open question

**This belongs in the paper and not only in `OPEN_QUESTIONS.md` section 8.** It is a
criticism of how the specification is CONSTRUCTED, and it will outlive the drift number
that provoked it: drift is merely the first unmodelled thing to exhaust a margin that was
never there.

### Current, at `:126`

> and reaches the gate only at $\num{tier1_visits_at_pub}$ monthly draws, which is nearly
> ten years of follow-up, against eleven with the increment known exactly.

### Proposed, replacing section 3's earlier draft of this passage

> and reaches the gate only at $\num{tier1_visits_at_pub}$ monthly draws, which is nearly
> ten years of follow-up, against eleven with the increment known exactly.
> \textbf{That count is defined as the one that lands on the target, so it carries no
> margin against anything the model leaves out.} A design solved to meet a bound exactly
> meets it only under the model it was solved in, and three things this paper does not
> model each exhaust it on their own: a drift in capacity of
> $\num{drift_cross_pct}\%$ a year over the decade the record spans
> (Appendix~\ref{sec:supp}); the difference between a Fisher standard deviation and a
> likelihood interval of the same posterior, which at the prior widths carried here is
> $\num{profile_gap_tier1_pct}\%$; and the span assumption itself, since
> Section~\ref{sec:local} warrants $\theta$ constant over the afebrile range and this row
> spends ten years of it. \textbf{The requirement should be stated with a margin, or
> stated as the boundary it is.}

### Numbers, and their owners

| macro | value | owner |
|---|---|---|
| `drift_cross_pct` | of order 1 | `audit/u7_drift.py`, `audit/u8_recheck_drift.py` |
| `profile_gap_tier1_pct` | 13.5 | `audit/u12_prior_quadratic.py` |

**This supersedes section 3 above**, which proposed the same passage with the drift finding
alone. Section 3's caution on quoting 0.92 rather than "of order one per cent" still stands.

---

## 7. A caveat that attaches to every Fisher figure, and where it should live

`audit/PHASE11_TRADE_AND_PRIOR.md` section 3 establishes that the gap between a Fisher
standard deviation and a likelihood interval of the same posterior is a function of the
$k_e$ prior width, and that at the 0.8 the paper carries it is 13 to 18 per cent at the
headline designs. The recommended record's $0.170$ becomes $0.200$ under a profile, which
is the target exactly.

**This does not belong in `11_scope`, and it does not belong beside any one number.** It is
a property of the whole method, so the proposal is a single remark in
`04_local.tex` beside `rem:curvature`, which already tells a reader the marker side's
interval is about twice as wide as its spread.

### Proposed, new, beside `rem:curvature`

> \begin{remark}[Which interval, and how much it matters at this prior width]
> \label{rem:whichinterval}
> Every precision in this paper is the square root of a diagonal entry of an inverse
> information matrix: a curvature at the optimum. That equals the posterior's standard
> deviation, and a likelihood interval's half-width, only when the log-posterior is
> quadratic, and at the $k_e$ log-prior of $\num{ke_prior_pct}\%$ carried throughout it is
> not. Profiling the same posterior over $\ln\theta$ gives half-widths of
> $\num{profile_lo_recommended}$ below and $\num{profile_hi_recommended}$ above at the recommended
> design, against a quoted $\num{sd_monthly}$: the interval is asymmetric, it is wider
> downward, and its mean half-width crosses the $0.20$ target. The asymmetry is a function
> of the prior width and not of the design, disappearing by a width of
> $\num{quad_width}$, and it therefore attaches to every figure computed at the width this
> paper uses, the Tier~1 count included. \textbf{The direction matters for a decision.} A
> symmetric interval understates how far capacity may fall below the estimate, and whether
> capacity clears a threshold is a question about the lower tail.
> \end{remark}

### Numbers, and their owners

| macro | value | owner |
|---|---|---|
| `profile_lo_recommended` | 0.268 | `audit/u12_prior_quadratic.py` |
| `profile_hi_recommended` | 0.133 | `audit/u12_prior_quadratic.py` |
| `quad_width` | 0.1 | `audit/u12_prior_quadratic.py` |

---

## 8. `paper/sections/09_validation.tex`: breath tests come in pairs

One sentence, because it is a protocol fact with an exact justification.

### Proposed, appending to the "Data per patient" paragraph

> \textbf{Breath measurements are booked in pairs.} The subject's tracer recovery is an
> unknown of its own and a single reading determines it and nothing else, so the first
> reading of a pair carries exactly zero information about capacity: the joint information
> matrix at $k=1$ equals the one at $k=0$ to machine precision
> (audit/\texttt{u11\_trade\_curve.py}, largest difference $3.7\times10^{-13}$ over
> eighteen designs). An odd number books one measurement that cannot be used.

---

## 9. Still NOT proposed

- **No change to any existing number**, as in round 1.
- **No change to the $0.20$ gate**, though section 7 makes which object it is stated
  against a question the author owes an answer to.
- **No wiring.** Every macro named above has to be generated before these passages build.
- **Nothing about carnitine reaches the paper.** `audit/PHASE11_KE_ROUTE.md` section 2 is a
  negative result about a lever the paper never proposes, so it belongs in the audit record
  and in `OPEN_QUESTIONS.md`, not in `sec:whoto`, whose existing carnitine argument is
  untouched and correct on its own ground.

---

## 10. DECISION (g3), SETTLED BY THE AUTHOR: the gate is stated against the profile

**The author's ruling**, recorded here because it changes what several numbers mean: *the
gate is the profile.* It is the object that does not assume the quadratic approximation
`audit/u12_prior_quadratic.py` part B shows fails at the paper's own prior width, and the
asymmetry runs the unsafe way, so a symmetric interval understates the lower tail. All
three objects are to be recorded at the recommended design; the **gate** is the profile.

**This supersedes section 7 above**, which proposed the caveat as a remark without
resolving which object the gate means. The remark is still wanted; what changes is that it
now states an answer rather than opening a question.

### What it does to the headline claim, stated plainly because it is not small

| design | Fisher marginal | simulated spread | **profile, the gate** | verdict against 0.20 |
|---|---|---|---|---|
| recommended, 25 monthly, nominal | 0.1695 | 0.1777 | **0.2002** | **does not meet it** |
| Tier 1, 118 monthly | 0.2000 | — | **0.2270** | does not meet it |

**The recommended design no longer passes its own gate**, by 0.0002. That is the
consequence of the ruling and it should be stated as such rather than absorbed.

### Proposed, replacing section 7's remark

> \begin{remark}[Which interval the target is stated against]
> \label{rem:whichinterval}
> Every precision in this paper is the square root of a diagonal entry of an inverse
> information matrix: a curvature at the optimum. That equals a likelihood interval's
> half-width only when the log-posterior is quadratic, and at the $k_e$ log-prior of
> $\num{ke_prior_pct}\%$ carried throughout it is not. Profiling the same posterior over
> $\ln\theta$ at the recommended design gives half-widths of $\num{profile_lo_recommended}$
> below and $\num{profile_hi_recommended}$ above, against a Fisher $\num{sd_monthly}$ and a
> simulated sampling spread of $\num{sim_spread_recommended}$. The asymmetry is a function of
> the prior width and not of the design, disappearing by a width of $\num{quad_width}$, so
> it attaches to every figure computed at the width this paper uses.
> \textbf{The target of $0.20$ is stated against the profile}, because it is the object
> that does not assume the approximation just shown to fail, and because a symmetric
> interval understates how far capacity may fall below the estimate, which is the tail a
> threshold decision turns on. On that reading the recommended design returns
> $\num{profile_recommended}$ and does not meet the target; the Fisher figure quoted elsewhere
> in this paper is retained as the design quantity it is, and is not the gate.
>
> \textbf{That failure is not a new defect, and the two facts belong in one sentence.} A
> design whose size is solved from $\operatorname{sd} = 0.20$ lands on the target by
> construction, so it carries no margin and will fail under \emph{any} object but the one
> it was solved against. The recommended record and the Tier~1 count are both solved that
> way. So what the stricter object exposes is the missing margin of
> Section~\ref{sec:spec}, not a separate problem with the estimate, and the repair is the
> one stated there: size the design against a target with room in it.
> \end{remark}

| macro | value | owner |
|---|---|---|
| `profile_recommended` | 0.2002 | `audit/u12_prior_quadratic.py` |
| `sim_spread_recommended` | 0.1777 | `audit/u4_recheck_breath.py` |

---

## 11. DECISION (g4), SETTLED BY THE AUTHOR: Tier 1 restated with a margin

**The author's ruling**: restate Tier 1 with a margin rather than as the boundary, and use
the frontier as the currency, since a count defined to land exactly on the target has no
way to express one.

**Why the frontier is the right currency.** "118 draws" is the solution of an equation, so
any margin added to it is arbitrary. The frontier is a requirement on an input, so a margin
on it is a statement about how well something must be known, which is the form a protocol
can act on and a measurement can be judged against.

### Proposed, replacing section 6 above

> \paragraph{Stated as a requirement rather than as a boundary.} The count above is the
> solution of $\operatorname{sd}(\ln\hat\theta) = 0.20$, so it lands on the target by
> construction and carries no margin against anything this model leaves out. Three such
> things are known and each exhausts it alone: a drift in capacity of order one per cent a
> year over the decade the record spans; the gap between a Fisher standard deviation and
> the profile the target is now stated against, $\num{profile_gap_tier1_pct}\%$ at this
> footing (Remark~\ref{rem:whichinterval}); and the constancy assumption itself, since
> Section~\ref{sec:local} warrants $\theta$ constant over the afebrile range and this row
> spends ten years of it.
>
> \textbf{The requirement is therefore stated on the input rather than on the count.} At
> nominal capacity the design reaches the target when $k_e$ is known to
> $\num{ke_frontier_3000_pct}\%$; specifying it to $\num{ke_spec_3000_pct}\%$ leaves the
> margin the count cannot express. At a fifth and a tenth of nominal the same requirement
> is $\num{ke_frontier_600_pct}\%$ and $\num{ke_frontier_300_pct}\%$, and at a twentieth
> no precision on $k_e$ suffices at all. A protocol should carry the input requirement and
> derive the count, not the reverse.

| macro | value | owner |
|---|---|---|
| `ke_frontier_3000_pct` | 80.0 | `audit/u9_ke_frontier.py` |
| `ke_frontier_600_pct` | 15.7 | `audit/u9_ke_frontier.py` |
| `ke_frontier_300_pct` | 6.3 | `audit/u9_ke_frontier.py` |

### The margin's size, settled by the author with a principle rather than a taste

**Size it to at least cover the spread among the objects the gate could have been stated
against.** At the recommended design those are Fisher **0.1695**, simulated sampling spread
**0.1777** and profile **0.2002**, a spread of **18.1%**. A design should therefore target
about **0.17 on the profile** rather than 0.20 on the marginal.

That is a floor and not a total. **Drift needs its own allowance on top**, because it is a
bias and not an object-choice, so it adds in a different way: `audit/u7_drift.py` reports
`rmse = hypot(sd, bias)`, and a drift of order one per cent a year over a decade contributes
a bias the object spread does not contain.

> \textbf{Sizing the margin.} The target is stated against the profile
> (Remark~\ref{rem:whichinterval}), and a design sized to meet it exactly carries no
> allowance for anything outside the model it was sized in. The floor adopted here is the
> spread among the three objects the same posterior admits, $\num{object_spread_recommended_pct}\%$ at
> the recommended design, so a design is sized to $\num{margin_target}$ on the profile
> rather than to $0.20$. Drift is an allowance on top of that floor and not inside it,
> being a bias rather than a choice of object.

| macro | value | owner |
|---|---|---|
| `object_spread_recommended_pct` | 18.1 | `audit/u12_prior_quadratic.py`, `audit/u4_recheck_breath.py` |
| `margin_target` | 0.17 | the author's ruling, on that spread |

---

## 12. `src/common/params.py`, the `k_e` row: the capacity qualifier

The only change proposed to the register, and it is the one `audit/u13` establishes.

### Current, at `:206`

> Section~\ref{sec:designsec} needs it within a factor of about three, and \textbf{no
> independent measurement in humans supplies it}.

### Proposed

> Section~\ref{sec:designsec} needs it within a factor of $\num{ke_needed_factor}$
> \emph{at nominal capacity}, and within $\num{ke_frontier_600_factor}$ and
> $\num{ke_frontier_300_factor}$ at a fifth and a tenth of it, which is where the clinical
> question is; below that no precision on this constant suffices, because the increment's
> own prior floors the design above the target. \textbf{No independent measurement in
> humans supplies any of them.}

| macro | value | owner |
|---|---|---|
| `ke_needed_factor` | 2.67 | `src/design/run.py`, already wired |
| `ke_frontier_600_factor` | 1.17 | `audit/u9_ke_frontier.py` |
| `ke_frontier_300_factor` | 1.07 | `audit/u9_ke_frontier.py` |

**Note the rounding trap.** At these widths a factor and a log-sd diverge: 1.07 as a factor
is 6.3% as a log-sd, and a reader who reads "a factor of 1.07" as "about one" has lost the
requirement entirely. The register may be better served stating the log-sd directly, which
is the form the frontier is computed in and the form a prior is specified in.

---

## 13. Still NOT proposed

- **No change to any existing number.** §12 adds qualifiers to a correct figure.
- **Nothing about the `k_e` measurement routes reaches the paper.**
  `audit/PHASE11_KE_MEASUREMENT.md` §2 prices three of them and closes none cleanly; that
  belongs in the audit record and in `OPEN_QUESTIONS.md` until one of them is costed
  against a real assay.
- **No wiring**, as before.

---

## 14. `paper/sections/07b_decision.tex`: the other side of the paper's own qualification

**Status: proposed, unapplied. This is the one passage of the round the author has asked to be
drafted for print.**

### 14.1 Why it is an extension and not a new result

The existing paragraph already carries the qualification this text completes. Its own words:

> ``The estimator fails towards declaring patients safe'' holds **only while the record tracks
> the load**; under the protocol a decompensation actually triggers, the same estimator declares
> every patient unsafe.

That sentence names one departure from `rho = 1` and prices it, `rho <= 0`. It does not name the
other. **The record can also OVERSTATE the load change**, `rho > 1`, and then the estimator
returns to declaring patients safe, no longer at the benign cell's rate but at 100%. The paper's
sweep runs `-0.25, 0, 0.5, 1` and stops, so it never enters that regime.

**And the regime is reachable by a protocol, not by a filing error.** It is what an incompletely
effective catabolism-suppressing regimen produces: the chart records the full prescribed cut
while the glucose absorbs most of the catabolic response, so the recorded change exceeds the true
one. That is standard emergency management, named in `01_biology.tex:160`, run electively.

### 14.2 The coincidence that makes it quotable

Two registered numbers pin the crossing without having been chosen with each other in view.
Writing `d` for the dietary share withdrawal removes, `F` for the load fold an unsuppressed
catabolic response would produce, and `q` for the fraction of that response suppressed:

```
    rho(q) = -d / ((F - 1)(1 - q) - d),        pole at  q* = 1 - d/(F - 1)
```

`params.py`'s illness row anchors `F = 1.50`, conservatively, from a 227% to 325% rise in
3-methylhistidine excretion in sepsis and trauma. The paper's own `rho = -0.25` then implies
`d = 0.100`, a withdrawal removing 10% of the load, **which is what a patient already on a
protein-restricted diet should show**. Those two together put the crossing at **`q* = 0.800`**.

### 14.3 Proposed text, to follow the existing sentence

> The qualification cuts both ways, and the second direction is reachable by protocol rather
> than by charting. Write the catabolic response the regimen suppresses as a fraction $q$ of
> what an unsuppressed fast would release; with the load fold of Table~\ref{tab:params} and the
> dietary share that $\rho = -0.25$ implies, $\rho$ passes through a pole at
> $q^\ast = \num{rho_pole_q}$, where suppression exactly cancels the prescribed cut and the true
> load does not move. Below it the record understates the change and the estimator escalates
> every patient, as above. Above it the record \emph{overstates} the change, and the estimator
> returns to declaring patients safe, no longer at the benign rate but at
> $\num{rho_above_fail_pct}\%$, \textbf{with no seed-to-seed variation at all}. A partially
> effective catabolism-suppressing regimen therefore lands in the more dangerous of the two
> failures, and the sweep of this paragraph does not reach it because it stops at $\rho = 1$.
> **The mitigation is a recording choice rather than a physiological one: when the load is
> delivered and metered rather than recalled, the record is the infusion and $\rho = 1$ by
> construction, whichever side of $q^\ast$ the suppression falls on.**

### 14.4 What must be true before this is applied

- **The numbers need generating and wiring.** The passage above is written with `\num{}` macros
  for that reason, and **an earlier draft of it used bare literals**, which is exactly the defect
  section 14.5 reports in the paragraph it attaches to. Caught and corrected here rather than
  carried in. **Two** new macros, `rho_pole_q` and `rho_above_fail_pct`, and not four: the
  seed range is **zero across all eight seeds**, so `\num{lo} to \num{hi}` would render as
  "100.0 to 100.0" and destroy the finding it was meant to protect. That clause stays prose,
  per the section 0 rule added at R143.
- **The excitation channel is not priced.** `src/design/run.py:1071` holds the marker levels fixed
  and scales only the record, so the 100% figure prices record corruption alone. A physiological
  `rho` also shrinks the true load deviation, which shrinks the marker excursions; the retained
  excitation runs 1.00x, 0.37x, 0.02x and then changes sign across the pole. **Near the pole the
  100% is therefore optimistic**, and the proposed text says "above it" rather than "at it" for
  that reason. The sign of the effect does not depend on the pole.
- **`q` itself is unmeasured**, which is gate (g6). It sizes how loudly the hazard is stated and
  not whether the mitigation works, since the pump gives `rho = 1` either way.

### 14.5 The paragraph this attaches to is entirely prose literals, and single-seed

Established in `audit/u19_seed_defect.py` and recorded at R134. An unbounded search of
`07b_decision.tex` shows lines 37, 40 and 43 macro-wired and lines 58 to 61 and 70 to 71 as bare
literals, fourteen numbers in all, every one at `seed=11`. **So this proposal should not be
applied on its own.** It belongs in one edit with the wiring and the spreads, or the section gains
a fifteenth unwired number.
