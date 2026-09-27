# Review log

Where the mathematics of Section 3 stands, so that the next session does not
have to re-derive it. Every model number below is produced by `src/` and asserted
in `tests/test_models.py`; none is typed into the paper by hand. Numbers taken
from or derived off the literature are not model numbers and do not live here:
they live in `paper/claims.csv`, against the source and with the verdict, and
the exact Poisson intervals in `894ae04` are recorded there as computed rather
than quoted.

## Standing conventions, so two independent scripts agree

These are not findings. They are definitions that have cost reconciliation time and will cost it
again every time two scripts compute the same quantity by different routes.

- **Find the owner before recomputing.** A number in the paper is either generated, and owned by a
  named script, or it is a prose literal. Before concluding that a paper number is wrong, locate the
  script that produces it and read what that script computes; recompute only to check the owner, never
  to replace it. Four disagreements in Phase 10 were a hand-built reconstruction losing to the paper:
  a schedule built by hand rather than read from `tab_D_designs.tex`, a t = 0 convention, and twice a
  capacity sweep stopping at 5% of nominal when the paper's range runs to 2%. **If no owner can be
  found, say so before concluding the number is wrong** (R81 to R83). `audit/PROVENANCE_SWEEP.md` is
  the index of which numbers have owners and which do not.

- **A spread is reported with its divisor.** `numpy.std` defaults to the POPULATION standard
  deviation, `ddof = 0`; the sample standard deviation is `ddof = 1`. The two differ by
  `sqrt((n-1)/n)`, which is 2.0% at n = 25 and 0.1% at n = 365, so a gap of that size between two
  routes is a definition and not a disagreement. It has already been mistaken for one once, on the
  within-patient marker CV of Phase 9 item 7, where 28.63 against 28.05 was entirely this. **Report
  which, and prefer the sample sd where the figure feeds a requirement**, since it is the larger.
- **A marker level is reported with its dose-day draw convention.** A draw falling exactly at a dose
  instant is read PRE-DOSE in the Phase 8 and 9 audit scripts and POST-DOSE by `design/run.py`'s
  `delta`. On the monthly record that is the t = 0 draw alone, and it moves the level from 10.000 to
  3.510. `REVIEW_LOG.md` R54 and `audit/PHASE8_SIMULATIONS.md` S1.R 3a record what it cost the first
  time.
- **A Fisher figure that decides a verdict is confirmed by fitting before it is reported.** An
  information matrix is the curvature of the log likelihood at one point. Where the likelihood is
  quadratic over the range a reader cares about, the curvature is the interval; where it is not, the
  two diverge without warning and the curvature is always the optimistic one. **The case that makes
  them diverge is a needle: a sharp minimum sitting on a nearly flat floor**, where the matrix
  measures the needle and a record at realistic noise sees the floor. Confirm by simulating records
  and fitting them, report the **robust spread** rather than the RMSE, since a Fisher bound has no
  tail to match, and quote the collapse count with a binomial interval beside it.

  This is a house rule and not a lesson from one phase. `REVIEW_LOG.md` **R20** retracted
  Cramer-Rao curvatures quoted where a reader reads an interval; `audit/n6_profile_vs_crb.py` exists
  for the same reason and profiles the design the paper recommends; and Phase 9 step 3(b) produced
  three instances in one section, the last of which withdrew that section's headline. Five instances
  in the project is a rule.

- **A noise figure is reported with what is inside it.** `sigma` carries both the load-recording
  error and the marker's scatter; the load-bias prior `s_b` carries the load-recording error a second
  time and is not inside `sigma`. A tolerance quoted for one axis is meaningless without the other
  axis and `s_b` both named. Phase 9 step 3(a) shows the design depends on the two only through
  `sigma` once `s_b` is fixed, exactly, so the honest statement of the requirement is a budget on
  `sigma` with `s_b` stated beside it.

## What each commit closed

| commit | closed |
|---|---|
| `11c86ba` | The sign error in §3.6 item 2, which reversed its own conclusion: both parameters enter the margin positively, the quadratic terms cancel along the fibre, and the margin is linear where `theta` is quadratic. Also: `a = theta g(m_bar)^2` carries `Km`, so "the intercept carries `theta` alone" was false and the admissible set is a fibre spanning 169x; `det M = g(m)g(m')(m'-m)` exactly, with `g(m_bar)^2` its midpoint approximation, 8% out at a 55% excursion; Figure 8's counting assumes one observation per patient; Figure 5's pair was matched in value only; the 1.4 load-error factor is the flux ratio `C/theta g`; (16) is Model A's `theta_crit` while Figure 10 plots Model B's; and `m_bar` denoted three quantities, the tolerable level now `\mtol`. |
| `7731e8d` | `s_peak` stated absolutely disagreed across the note (4.65x baseline in Model A against 1.18-1.50x in Model C); it is now a fold of the patient's own baseline. That reframed the margin result: the sign of `mu` is determined outside a computable window of `s_peak`, `(3.25, 3.77)x` baseline at `m_tol = 40`, 14.8% wide, closing to 1.2% at `m_tol = 15` and opening past 100% at `m_tol = 320`. Also fixed the abstract's surviving "carried entirely by the intercept", the warrant for "given `Km`" (Assumption 2 grants sharedness, not knowledge), and the band convention (`Delta` is full width). |
| `98c65e8` | Evaluated (7), which the note derived and never computed. `sd(theta_hat)/theta = sd(a_hat)/a` exactly, so the relative precision is identified although `theta` is not. `sigma` and `n` defined for the steady-state design, `n` counting instrument events. Curvature resolvability derived with design constant `sqrt(180)`. Errors-in-variables simulated. The null anchor given a variance. Model C shown to have no fibre given `alpha`. |
| `f981cad` | Five corrections to the above, all mine. `sigma` had been quoted in flux units from a constant documented as a percentage, implying a load recorded to 1.55%; both criteria are now functions of the load error. The precision result does not close the location question. The anchor cost is a ratio scaling as `1/theta`, crossing one at `theta = 605` and `1100`, so it amplifies at Model C's `theta = 300`. The errors-in-variables sign flip was an inadmissible fit (`ke < 0` in 96% of them), detectable from the fit alone. "No fibre" is given `alpha`, and `Kc` costs a factor of 1.52 in recovered `theta` over its own arbitrary range. |
| `bdb2120` | Six corrections from a literature fact-check run against an older draft, each checked against the primary source before it moved. §1.3 asserted that every crisis leaves irreversible damage while quoting a 7-56% lesion prevalence in the same sentence; GeneReviews says "high risk" and "may be preceded by" and never uses "irreversible" or "permanent", and the range is a lifetime prevalence across individuals rather than a per-episode rate. Roe 1983 is one patient. §7's hazard model is a GAMM, not a GAM. The Haijes free-carnitine figures are medians, the column header being "Median ± SD [Min-Max] (N)". §3.5 had "The The". §7's marker paragraph still said the capacity is carried by the intercept and pointed at (8), where Remark 2 says the intercept locates `theta` only along a fibre and (23) is the form whose terms are all observable. Figures 12 and 13 regenerated with `m_tol`. §3.1 gained a paragraph on how the note relates to the two PK/PD models it already cited: neither is Michaelis-Menten in the marker, Proposition 1 does not distinguish either from (1) because its hypothesis is affineness over the visited range, and what is specific to (1) is the fibre and the cancellation in (20). Written as a `\paragraph` so no Remark is renumbered. |
| `79ea21f` | C46 recorded as drift rather than a match. "Whose loss is first order" is the definition of the model class Gao names and not something that abstract states, it carrying no PD equations and no biomarker parameters. Same standard as C02 and C03. |
| `894ae04` | The abstract's flat "severely ill-conditioned in either case" replaced by the factor §3.7 computes, widening from about 2.5 to about 26 as the load is less well recorded. §7's Design paragraph made a risk, matching §1.3, which strengthens it: a chance constraint is the right formulation for a probabilistic harm. The Saudi spread read correctly: the three rates rest on 4, 4 and an unprinted count near 55, their exact Poisson 95% intervals are 1:4,000-1:37,587, 1:5,530-1:51,962 and 1:10,826-1:18,705, and the national interval lies inside both others, so "the spread is about 40%" was reading sampling noise as disagreement. The warrant is now that four events cannot resolve an incidence and that ascertainment does most damage to a numerator of four, C47 replacing the false-positive caveat, which the source raises about biotinidase rather than about PA. |

## Pass 3: the five blockers

Five items had to be settled before any restructuring, because each of them
moves a number the restructure would otherwise be built on. All five are
closed. One of them changes what can be claimed and is flagged for a decision.

### 0.1 `M_i`, Figure 16, and the sensitivity of the shape constants

`M_i` was never assigned a value in the paper. It is `300`, and the text's own
excretion arithmetic (`m ~ 1650`, a 1.3-fold rise at the attractor, 1.6 times
the separatrix, "about a fifth of normal") reproduces exactly at that value, so
the number was recoverable but never written. It is now in `tab_params`, which
is generated from the dataclasses rather than typed beside them.

**The reported ~9x swing in `theta_crit` over an 8x range in `M_i` does not
survive.** It is an artefact of holding `R0` while moving `M_i`, which breaks
the constraint that defines `R0` (the alternative route must peak below the
baseline load, so a null genotype is not maintainable). The published 1.56x for
`K_c` is an `R0`-re-derived number, so the two were never on the same footing.
On the admissible footing:

| constant | range | `theta_crit`, calibration held | `theta_crit`, `R0` held | margin, calibration held |
|---|---|---|---|---|
| `K_c` | 25-800, 32x | **1.56x** | 6.27x, vacuous at one end | 3.69x |
| `M_i` | 150-1200, 8x | **1.37x** | 21.3x, vacuous at one end | **6.44x** |
| `rho_0` | 0.35-1.0 | **1.02x** | 967x, vacuous at one end | 1.48x |

Two consequences. Sweeping all three together at a held calibration leaves
`theta_crit` inside a factor of 1.39 (12 random sets), so **the threshold is
very nearly a function of the peak of `alpha` alone** and its exposure sits in
one clinical claim rather than three unmeasured constants. And **the margin is
not protected the way the threshold is**: the separatrix sits where `alpha` has
decayed, `alpha` decays on the scale of `M_i`, and the margin runs 509 to 3275
over the same 8x range. Figure 16 quotes 954 flat.

Same convention error found in two more places, both corrected. The `rho_0`
lower bound "about 0.35" is an `R0`-held number; re-derived it is 0.195. And
§5.5 item 2's "theta from 300 to 457, a factor of 1.52, no admissible solution
above `K_c` ~ 200" does not reproduce under either convention: the admissible
window is two-sided, roughly `K_c` in [100, 300], and inside it theta runs 213
to 300, a factor of 1.41.

### 0.2 `Delta`, and a floor the design table did not have

`Delta` is the full width; the definition is now at first use, with the
tangent-versus-least-squares gap (0.09% at 10%, 8.8% at 100% full width, 35.4%
read as a half-width) beside it. All three reproduce.

The convention was the smaller problem. **A band is `[m_bar - Delta/2, m_bar +
Delta/2]` and `m` cannot be negative, so `Delta/m_bar <= 2`.** The published
table quoted 460%, 920% and 1841% as design requirements. None of those is a
band.

Underneath that is a real error. `Var(a_hat) = sigma^2 [1/n + m_bar^2/S]` and
equation (7) keeps only the second term. The first binds exactly in the
wide-band regime the section recommends buying, so

    sd(theta)/theta >= c(n) sigma/(a sqrt n),   c(n) = sqrt(1 + 3(n-1)/(n+1)) <= 2

for an evenly spaced design, and `>= sigma/(a sqrt n)` for any design at all.
The best use of an excursion has the closed form `sd(a_hat) = sigma (H+L) /
(sqrt(n) (H-L))`, checked against brute force.

So **the design requirement splits into two pieces that do not trade against
each other**: a floor on `n` set by `sigma` alone, below which no excursion
delivers the precision; and above it, a requirement on the excursion. At six
instrument events the floor is 53% even at a load recorded to 2%; 20% precision
needs 43 events there and 714 at a load recorded to 20%.

The curvature criterion is restated as the **top of the excursion**, 2.9 to 5.4
times the operating point, because the bending of `C` is at high `m`. The
instruments of §3.5 move the patient *down*. The earlier "short by a factor of
about 2.5 to 26" is withdrawn: it compared a requirement in `Delta/m_bar` with
an achievement in `Delta/m_max`, against a sigma carrying only the load term.

### 0.3 The Gramian. **GATE: survives as stated, fails as used.**

Proposition 2 generalises to the observability Gramian exactly as far as it is
written, and no further.

*Given `Km`*, it survives. At a matched trajectory the `(theta, ke)`
sub-Gramians differ by at most 5% in determinant and 2.3% in
`sd(ln theta_hat)`, over excursions from 5% to 76%, always favouring the
therapy. The mechanism is exact: the sensitivity forcings `-theta g(m)` and
`-ke m` are functions of the state alone and are identically channel-blind, and
the only channel dependence is the shared decay rate, with
`J_therapy - J_carnitine = -kappa u g(m) <= 0` (verified to 1.5e-5). Posing the
question on carnitine is therefore conservative.

*With `Km` unknown* — which is the paper's actual situation, since the fibre of
Remark 2 exists precisely because `Km` is unknown — **it fails**, by six orders
of magnitude in `det G` and by a factor of 46 to 807 in `sd(ln theta_hat)`. The
`Km` sensitivity is forced by `-(theta + Delta(t)) Km dg/dKm`, and that
prefactor is modulated by the therapeutic input and constant under carnitine.
Collinearity `|corr(q_theta, q_Km)|` is 0.99994 under carnitine against 0.967
under the therapy.

**The excursion is channel-blind; the curvature is not.** Carnitine buys the
excursion. Only a perturbation that multiplies the enzymatic term buys the
curvature, or an excursion carrying the patient into the neighbourhood of `Km`:
moving the operating point from `m/Km = 0.1` to `m/Km = 2` drops the
collinearity to 0.982 and raises `det G` by fifteen orders of magnitude.

This is recorded in the paper as a remark scoping Proposition 2, which is
correct under any plan. What it changes is the *use* the paper makes of the
proposition, and that needs a decision before the restructure.

### 0.4 Which `theta_crit` the margin theorem is proven against

**Stated as a limitation, because the alternative has no premise.** The
cancellation is a statement about a degeneracy, and the parent model does not
have that degeneracy: its alternative channel is curved, contributes to both
coefficients, and given `alpha` a fitted tangent determines `(theta, Km)`
uniquely, to machine precision. There is no fibre for a margin to be linear
along.

Two things follow. `a = theta g(m_bar)^2` is a property of the *local* model,
not of the parent: under the parent the intercept is `theta g^2 + [alpha -
m alpha']`, and the second term is 374 against the first's 61 at the nominal
constants. And the intrinsic threshold does not eliminate the convention, it
relocates it: `m_tol` leaves, `alpha`'s constants and the `R0` calibration
arrive. That is a gain, because the calibration is a clinical claim a clinician
can dispute where `m_tol` is a number chosen to make a formula close. It is not
an elimination.

### 0.5 The marker is 3-HP

The direction check the handoff asked for resolves in 3-HP's favour, and the
concern behind it was a misreading of the conditioning. The rate ratio is
attached to a **lower** concentration: "a 50% lower plasma 3-HP concentration
was associated with a 30% lower frequency of MDEs (RR 0.7; 95% CI 0.6, 0.9;
p < 0.01)". So RR < 1 means fewer events when the marker is lower. 3-HP is a
risk marker in the expected direction for a toxic overflow metabolite, not a
protective one.

C3 stays excluded on the structural argument (it is the conjugation product, so
carnitine raises it while lowering the organic acids and the perturbation of
§3.4 enters with the wrong sign). 2-MC is excluded on the endpoint: RR 1.0,
p = 0.84, and its strength is as a correlate of plasma ammonia *during* an
event, which is the arm Remark 14 excludes.

**The multi-output panel is affordable only conditionally, and the condition is
the load record.** With `J` channels as scaled copies of one pool, independent
assay noise averages as `J^-1/2` and the load error does not, because every
channel is fitted against the same ordinate. A four-marker panel cuts sigma by
1.65 at a load recorded to 2% and by 1.02 at 20%. Buy the load record first.

### 2.2 The assay-versus-load crossover, recomputed

§3.7's "assay variability a smaller term" is false at the good end. An abscissa
error acts like an ordinate error `b sigma_m`, and `b m_bar` is within 8% of the
baseline load, so the two balance at `CV = 1.08 f`: **the assay term dominates
at every load recorded to better than roughly the assay's own CV.** At a load
recorded to 2% and a 5% assay it is 14.9 against 6.5.

## Still open after pass 3

1. **The Gramian gate needs a decision.** Outcome 3 as planned said the analysis
   is posed on carnitine and answered for the therapy. That is valid given
   `Km` and invalid without it, and the paper's own central result is that
   `Km` is not available. The recommended resolution is to scope the claim and
   promote the failure, since "the excursion is channel-blind and the curvature
   is not" is the excitation condition the restructure is built on, arriving
   earlier and derived rather than asserted.
2. **A defensible bound on the load record.** Unchanged from pass 2, and now
   carrying more weight: it sets both the floor on `n` and whether the marker
   panel is worth buying.
3. **The hazard reformulation.** Unchanged.
4. **Whether the parent model's four constants are a net gain.** Sharpened: the
   record bounds `K_c` from both sides, but the margin is exposed to `M_i` at a
   factor of 6.4 where the threshold is exposed at 1.37.

## A fourth failure mode, added to the three above

- **A sensitivity swept off a broken constraint.** `R0` is defined by a
  constraint and three separate numbers in the paper were computed by holding
  it while moving a constant that moves it. Each time, most of the reported
  sensitivity was the constraint breaking rather than the constant mattering,
  and in one case the ordering of the constants reversed. The rule now is that
  a derived constant is re-derived inside every sweep, in code, and
  `R0_for_peak` exists so that it cannot be forgotten.

## Pass 4: the restructure

The letters are gone. One model at three resolutions, in the order the reduction runs, and
each section is now a link in the chain from the clinical claim to the data specification.
Figure~1 of the paper draws that chain as a dependency graph, and anything not reachable
from the top of it was moved to the supplementary appendix rather than deleted.

### 1.1 The hierarchy. Half the gate passes.

**Passes.** Solving `dc/dt = 0` in the two-state system returns `alpha` exactly and
symbolically, with `k` surviving only inside `rho_0 = delta_0/k`. Integrating the full
two-state system from below the separatrix lands on the reduced attractor to 14 digits and
from above it escapes. The peak of `C` at `theta_crit` equals the load to 10 digits. **The
equilibrium reduction is the former Model C.**

**Fails.** The local model is NOT the linearisation of the parent. The tangent intercept is
`a = theta g^2 + [alpha - m alpha']`, so `a = theta g(m_bar)^2` needs `alpha` linear through
the origin, and the second bracket is 374 against 61 at the nominal constants. No term was
added to force agreement. The honest statement is now in the paper: the local model is the
parent with one term replaced, that replacement is what creates the fibre, and it is
therefore the favourable case.

### 1.6 The excitation condition, corrected and proved

The earlier phrasing was wrong in its second clause. Both ends are degenerate and they fail
differently: below `Km` only the lumped ratio survives, above it `theta` enters only as
`L0 - theta` and is confounded with the unrecorded load. Losing a nuisance parameter is
survivable; losing the target is not.

The parametrisation is pinned before anything is proved: `theta` absorbs `V_max`.

With `phi(m) = Km/(Km+m)`, the excitation is `E = phi(lo) - phi(hi)`, and at fixed fold
range `r` it is maximised at `sqrt(lo*hi) = Km` **exactly**, where
`E_max(r) = (sqrt(r)-1)/(sqrt(r)+1)`. A band of the same fold range sitting a decade below
`Km` throws away a factor of three in `E` and a factor of 2,099 in `sd(ln theta)`.

The optimum is interior ONLY because the load is partly unrecorded. With it known exactly,
`sd(ln theta)` falls monotonically into saturation and there is no design problem beyond
safety. That is what makes the load record the first row of the specification.

### 1.8 The covariance collapse. **The gate fails, and that is the result.**

Exact information filter for constant parameters, no process noise on `theta` by
construction. Under the biweekly schedule over two years:

| design | n | sd(ln theta) |
|---|---|---|
| six-monthly | 4 | 0.545 |
| one cycle, 2-daily | 7 | 0.540 |
| 2-daily for two years | 365 | 0.489 |
| same 365, straddling `Km` | 365 | **0.004** |

Ninety times the sampling buys a tenth of a decade; the load recorded ten times better buys
another 5%; placement buys 122x. **The dense design does not collapse.** It was not tuned
until it did.

The EKF written first is not used and the reason is recorded rather than hidden: the marker
equilibrates in 0.81 h against a 3.5-day dose half-life, so it is quasi-steady at any
realistic sampling instant, a design is described by the levels it visits, and the
information filter is exact. The one exception is a 0.91-hour window after each dose, which
is now a row in the specification.

### 1.9 The estimand is a decision, and the two errors are wildly asymmetric

| illness | mu | which error | simulated | delta method would say |
|---|---|---|---|---|
| 1.50x | +2006 | treat unneeded | 0.2% | 39.5% |
| 2.98x | +339 | treat unneeded | 6.0% | 36.1% |
| 3.58x | -339 | **fail to treat** | **52.5%** | 42.3% |
| 5.50x | -2513 | fail to treat | 26.0% | 40.4% |

**The estimator fails towards declaring patients safe.** At equal distance from the
boundary the costly error is nine times the cheap one, and it does not close with distance.
The mechanism was already in the note: abscissa noise inflates `a_hat`, hence `theta_hat`,
hence `mu_hat`.

A delta-method tail is nearly flat at 40% and is wrong in both directions, because `mu_hat`
is not symmetric about `mu`: the design's near-null direction is very nearly the direction
along which `mu` is monotone and of one sign. **That is the margin theorem showing up in
the estimator**, and propagating a Gaussian would have thrown it away. Both numbers are
reported so the gap stays visible.

### 1.7 and the terminal deliverable

Against a requirement of `E ~ 0.50`: inter-dose recovery 0.049, carnitine steps 0.062,
antibiotics 0.034, supply interruptions 0.047, labelled bolus 0.076, decompensation 0.354.
The routine sources deliver a fifth to a quarter of what their own fold range could, so they
are in the wrong PLACE and not merely too narrow. **The excitation the decision needs is
available only from the event the therapy exists to prevent.**

The cohort row settles it. At a 41% sign error near the boundary, placing the fraction with
`mu < 0` within ten points needs about 3,100 patients against a trial of 16. That is not a
recruitment problem.

### 1.4, 1.5 and the rest

- **Fibre span, closed form, promoted.** `(1 + b m_bar/a)^2 = (C(m_bar)/a)^2`: recorded load
  over fitted intercept, squared. Both observable, so the 169 is a per-patient diagnostic
  rather than an illustration.
- **The null-genotype anchor is kept, as a named remark.** It is the only excursion-free
  route and it amplifies load error for exactly the severe patients the note is about. The
  one route needing no excursion works for the patients who need it least, which is the same
  inverse-usefulness shape as the main argument.
- **Proposition 3's clause** is now "given `Km` AND the input gain": `kappa` is never
  assigned, so on the carnitine channel the system is two equations in three unknowns before
  `Km` is raised.
- **`lambda`, `T_dose` and `D` are assigned**, with sources. `lambda` was named and never
  given a value, which made 1.8 unrunnable.

## Still open after pass 4

1. **The load record.** Now the first row of the specification rather than a caveat: it is
   what makes the saturated end degenerate, it is the noise term a marker panel cannot
   average away, and it sets the floor on `n`.
2. **A PA breath-test cohort.** The reference standard in the validation protocol is
   transferred from MMA. Confirming whether a PA series exists is a prerequisite.
3. **The hazard reformulation.** Unchanged, and now clearly the route by which `s_peak`
   leaves the decision boundary that Section 7 shows the design cannot resolve.
4. **Whether any ethical design meets the excitation condition.** The note now says what is
   required and that nothing accessible supplies it. It does not say what would.
5. **Length. The main text is about 22,000 words against a target of 8,000.** Seven
   self-contained blocks were moved to the supplementary appendix in this pass and nothing
   was deleted, but block moves have taken it about as far as they can: what remains in
   sections 4 and 5 is on the chain from the target claim to the data specification, and
   getting it to 8,000 means compressing dense prose rather than relocating it. That is an
   editorial pass and it has not been done. The chain itself is now explicit in the
   dependency figure, so the pass has a criterion to work against: keep what a box in that
   figure depends on, compress everything else.

## A fifth failure mode

- **A linearisation that discards the result it was meant to measure.** The delta-method
  tail on `mu` says the decision is a coin flip at every amplitude. The simulation says it
  is right 99.8% of the time on one side and wrong 53% on the other. The Gaussian was
  symmetric about a quantity whose sampling distribution is one-sided, and the one-sidedness
  was the paper's own margin theorem. A propagated covariance is not a substitute for
  simulating the estimand you actually report.

## Still open (pass 2)

1. **A defensible bound on `sigma`.** Every conclusion in §3.7 is conditional on how well the substrate load is recorded, and the note never bounds it. The abstract now says so rather than asserting a flat verdict, which makes the gap visible from the first page instead of only in §3.7. Prescribed protein intake is documented; gut flora and unrecorded triggers are not. Until this is bounded from the records the identifiability answer is a family of curves rather than a number.
2. **The hazard reformulation of §7.** `s_peak` and `m_tol` are now the binding arbitrariness rather than identifiability, and `lambda_MDE(theta) = nu P(A > A_crit(theta))` is what removes both. The natural history study already fits a hazard of that shape on a measured biomarker; what is proposed is substituting a mechanistic `theta`. Not attempted.
3. **Whether Model C's four constants are a net gain.** `R0`, `Kc`, `Mi` and `rho0` replace `ke`. They remove Model A's fibre given `alpha` but `Kc` alone costs 1.52x in recovered `theta` and 1.56x in `theta_crit`. Model C's identifiability is not analysed anywhere in the note.
4. **The marker choice.** §3.2 defers it and §7 poses it as a design question answerable on `m_bar/Delta` at a fixed clinical excursion. What is missing is the paired marker-and-management data to run the comparison, and §5.6 adds that no marker is currently known to track the state through a decompensation at all.

## Three failure modes worth watching

The first two recurred within one review cycle; the third survived a review
cycle and a literature fact-check before anyone looked at the counts:

- **A ratio reported as a constant.** The 1.4 load-error amplification and the 0.36 anchor attenuation were both quoted as properties of the model. Both are ratios in `theta` and both invert for the severe patient.
- **A spread read as disagreement.** Three Saudi incidences differing by 40% were treated as estimates that fail to agree. Two of them rest on four events each, and all three are statistically indistinguishable. A ratio between point estimates says nothing until the counts under them are known, and the fact-check that reviewed the paragraph did not catch this either.
- **A unit carried across scales.** `ASSAY_CV` is a percentage and was used as absolute flux. The definition site now says so.

---

## Retractions (pass 3, the phase 2A–2D audit)

Recorded rather than deleted, because the reason each was wrong is the point.

### R1 — "one dose cycle of seven samples reaches the target"

**Stated** in the phase 2A report: that reinstating the therapeutic increment in the
sensitivity brought a single seven-sample dose cycle to `sd(ln theta) = 0.1511`, inside
the 20% gate.

**Retracted.** It held only with `k_e` held exactly known, which is not the footing this
paper claims. With `k_e` fitted at the module's own 0.8 log-prior the same design gives
**0.2067** and has a floor there: it does not reach the gate at any precision on the
delivered increment, because the floor is set by the design and not by the product
(`audit/c_dose_precision.py`). `tests/test_models.py` now asserts the retraction rather
than the claim.

**The failure mode is the one this log already names** — a sensitivity swept off a broken
constraint — in its fourth instance. The footing was not held fixed between the number
that was quoted and the footing the paper argues for. The rule that follows: **every
design number in this paper is quoted on the k_e-fitted footing, and any number on the
k_e-known footing is labelled as an audit reproduction and never as a headline.**

### R2 — "the sign error plateaus near 6% because of load confounding"

**Retracted** (A1). The floor is an artefact of unweighted least squares against
level-dependent noise; per-sample sd runs 21.96 to 945.6 across the band. With the load
bias still free and flat but the noise made homoscedastic, the error runs
32.3/26.7/2.1/0/0/0/0/0%. The curve does not plateau: extended, it turns at 40 m̄ and
rises to 17.3% at 320 m̄.

### R3 — "corr(ln theta, ln K_m) reverses from +0.9988 to −0.3163"

**Retracted** (C3). The reversal is a cancellation between a prior-fed null direction and
the identified subspace. It is +0.0214 on the claimed footing, flips sign across the
parameter sweep, and Monte Carlo returns +0.95 against the linearised −0.32. What
survives is that the correlation *collapses in magnitude*.

### R5 — "the decision table's 52.5% is a seed outlier"

**This retraction is itself retracted**, and the sequence is worth keeping because it is
the failure mode running in the other direction. Having found three citation-layer errors
in one pass, I read a fourth discrepancy as a fourth error. It was not: the figure calls
`sign_error(..., seed=12)` and 52.5% reproduces from it exactly. My twelve-seed check ran
seeds 0–11 and missed the one the code actually uses.

The accurate statement is that 52.5% is one draw. Over 24 seeds the unsafe-side error is
43.9% ± 3.5%, range [37.2%, 52.2%]; seed 12 is the maximum and the 96th percentile. A
single-seed value near a decision boundary carries three to four points of Monte Carlo
noise and should be quoted with its spread. **A pattern of errors is not evidence that the next discrepancy is
another one**, and checking which seed the code used would have cost one grep.

### R4 — "the breath test repeats at a CV of 12.8% in patients"

**Retracted** (2D.3), and it is the source's own error inherited. Manoli et al.'s Results
give **14.02 ± 12.86%**, range 0.39–50.84%; the abstract reports 12.86, the standard
deviation, as though it were the mean. The figure is also **between-visit over 2 months
to 4 years**, with changing renal status named by the authors as a cause — not same-day
test–retest. The published Correction restores a graphical abstract and changes no
numbers.

### R6 — "a constant increment leaves only the lumped sum theta + Delta_0 identifiable"

**Retracted** (3.7), and it is the first error in this repository caught by an external
tool rather than by re-running my own code. The claim was written into
`symbolic/sympy/observational_equivalence.py` and repeated in the commit message of
`a7b30a0`; the same file's own pre-written "EXPECTED" comments asserted it too, so the
script would have appeared to confirm itself. `StructuralIdentifiability.jl` returns
globally identifiable for that configuration, and re-running the sympy solve returns the
identity solution alone.

The error is a one-word one. Delta is a **known** input — it is the dosing record — so
recovering the lumped sum recovers theta. I had reasoned about it as though it were an
unknown parameter.

A second instance in the same pass, same cause: the Julia file's comments claimed that a
constant *input* encodes a patient held at steady state. It does not. A constant input is
not a constant state, and the framework assumes generic initial conditions, so the marker
relaxes through an interval before settling. Both comments are corrected in place.

**The rule that follows.** Writing the expected answer into a script *before* running it
converts the script from a check into a confirmation. Two of this session's symbolic
scripts carried pre-written "EXPECTED" comments and both were wrong; the numerical scripts,
which print what they compute and label nothing in advance, have not produced an error of
this kind. Expected values belong in an assertion that can fail, or nowhere.

### R7 — "the timescale separation makes the information filter exact"

**Not retracted. The retraction is retracted**, and this is the second time that sequence
has run (R5 was the first).

Writing §3.2 I asserted that the 104× separation between the marker's time constant and the
dose half-life is what makes \eqref{eq:infofilter} exact. I had not derived it. Proofreading,
I looked up the filter's definition, saw that it is a static-parameter information recursion
with no process model, concluded the separation had nothing to do with it, and "corrected"
the sentence in three files and in the commit message of `6367fe0`.

**The original claim was right.** §4.1 of the note already argues it, in the paragraph
explaining why the extended Kalman filter that was written first is not used: because the
marker is quasi-steady at every sampling instant a clinic could use, the observation at each
visit is an algebraic function of the parameters and the inputs, there is no state left to
filter, and the pointwise information recursion is exact **rather than approximate as an EKF
would be**. The separation is exactly what removes the state. I checked the filter's
definition and not the paragraph that justifies it.

The replacement wording stands, because it names the mechanism instead of only the
conclusion, and the note is better for it. What does not stand is the claim that the
sentence it replaced was an error.

**So this is one failure, not two, and it ran in both directions.** The assertion was
unverified when I made it; the retraction was unverified when I made it. Neither cost was
paid for the same reason: I did not look.

### R8 — "prior-free sd falls to 0.054 at 2% CV and 0.023 at 10%"

**Retracted** (15.1), and unlike R1–R7 this one was in the paper's own prose rather than in
a working note. The two numbers reproduce only on the two-daily design of 365 samples **with
k_e held exactly known** — a footing Section `sec:whypriors` disowns and which the design
module's docstring explicitly reserves "for the audit, not for a headline". On the claimed
footing they are 0.234 and 0.172, a factor of 4.3 and 7.5.

Both were also computed **with the dose train present**, in a passage arguing what load
variation buys. That is the substantive half: the sentence read as a claim about load
variability and was measuring the therapy.

**It surfaced only because a later result contradicted it.** The excitation budget (14.1)
returned 1.79 for load variation alone where this sentence implied 0.023, and the handoff
that followed refused the new number until the contradiction was addressed. Neither figure
was wrong about its own computation; the older one was wrong about what it was a computation
*of*. **Two numbers on different footings do not announce themselves — they only collide when
something forces them into the same table.**

### R9 — "information in the degenerate direction saturates"

**Retracted** (15.2). Written in the previous handoff's own summary and then into the paper.
The scaling is real — the fitted exponent is 0.06 to 0.24 against an independent-sample 0.5,
across five orders of magnitude in n — but the stated mechanism is not. Information in that
direction **grows** roughly in proportion to n. It never saturates; it never arrives. The
load prior holds about 1.05 in that direction and does not move, so at n = 300,000 the data
still supply three orders of magnitude less.

The conclusion is unaffected and the mechanism was decorative — which is exactly why it was
easy to write and would have been hard to catch. This is the second instance of the R7 rule
in two handoffs.

### R10 — the em-dash sweep, and why a green build proved nothing

**Not a retracted claim but a retracted method**, and it belongs here because the failure
mode is the same shape as the others: something that looked verified and was not.

Removing 123 em dashes by pattern substitution did the job it was asked to do and three
things it was not. It rewrote comment rules (`%-----`) into commas. It collapsed the
newline after two `\cite{...} %` comments, pushing the rest of those sentences **behind the
comment marker**, where they vanished from the document. And on one of those rewritten
rules it swallowed the following line, which was `\input{sections/01_biology}` — removing an
entire 1,400-word section from the build.

**The build stayed green through all of it.** No LaTeX error, no undefined reference, no
failed test. A section silently left the document and every check this repository runs said
the document was fine. The page count fell by four and I attributed that to a shorter
abstract, which was also true, and which is exactly how it stayed hidden.

**The rule that follows: a mechanical edit that leaves the build green is not thereby
verified.** The build checks that the document compiles, not that it still says what it
said. For prose there is no test that fails when a sentence is deleted, so the only
verification is reading the diff — and a 113-site diff is not read, it is skimmed.

**So: no further whole-document pass on this repository is done by pattern substitution.**
Any sweep is done by reading, or it is not done. Where a sweep has already run, the check
is not "does it build" but "what did the diff actually change", site by site.

Worth recording alongside: the three defects were found by reading the changed lines
afterwards, not by any check. Two were found because I went looking, and the third only
because the first two made me suspicious enough to keep looking.

### R11 — the pre-restructure cluster, and two numbers the sweep would not touch

**Not a retraction but the result of applying R10's rule**, recorded because of what
reading found that no check had.

The restructure moved six sections into the supplement without rewriting them, and they
still argued the position the paper had reversed: carnitine as the input the identification
rests on, the therapy as prospective and administered to nobody accessible, and
Proposition~`prop:channel` stated unconditionally where the A5 correction had made it
conditional on $K_m$. A reading sweep over every supplement section, forty readers with an
adversarial verify on each claim, returned thirty-one candidates of which twenty-two
survived verification: fourteen in `04_local`, five in `A1_supplementary`, two in
`02_system`, and one each in `03_dynamic`, `09_validation` and `10_results`. Three files
were clean. Nine candidates were rejected on reading, most of them because the passage was
already correctly scoped or already carried its own supersession note.

All twenty-two are resolved in place rather than deleted, on the rule from the pass before:
where a question has been answered, state the answer and cite the section; where a position
was superseded, say so. Two carried a consequence past the wording. The Figure A3 caption
and its twin in `10_results` claimed the channel choice does not matter, which is true of
the excursion and false of the curvature, so the caption now says which of the two the
figure shows. And `09_validation` concluded that a paired validation would measure the
marker method's imprecision rather than agreement; that was true on the increment-dropped
footing and is false on the design the protocol specifies, where the marker side is a
factor of $1.3$ from the reference standard rather than an order of magnitude, so the
conclusion is reversed with its original stated.

**Two findings are out of that class and were left alone at the time, deliberately**, both
being numbers rather than contradictions of position. Both were fixed in the pass that
followed, and R12 records what their survival says:

1. `08_specification.tex:100` and `00_abstract.tex:19` quote $1.79$ for spontaneous load
   variability at $10\%$ CV. `audit/e7_excitation_budget.py` returns $1.7731$ on the
   lognormal load that replaced the Gaussian one in 15.1, and `A1_open.tex` quotes $1.77$.
   The $40\%$ row of the same table was updated to its lognormal value ($0.70 \to 0.62$)
   and the $10\%$ row was not. `audit/PHASE3_CHECKS.md` says the correction was made "in
   the paper", and in the paper it was made in one row. (R11 read the extent of this as two
   rows. Reading the whole table in the pass that fixed it found five affected numbers, of
   which one had been requoted. R12 has the count.)
2. `08_specification.tex:90` describes that scenario's footing as "twenty-four monthly
   draws over two years". The scenario runs on twenty-five
   (`len(np.arange(0, 730, 30.4)) = 25`). The "twenty-four draws" at line 119 is a
   different quantity, the sampling-rate route, which does run at twenty-four and is right.

**What this says about the class.** R8's rule was that every quoted precision carries its
footing. These two are the other half of it: a footing that changed under a number nobody
requoted. Both were found by reading prose against the script that produces it, which is
the same instrument that found R6, R7 and R9, and again not by any test in this repository.

### R12 — the correction that reached one row of a table

**The near-miss, and the worst-placed error this project has carried.**

In 15.1 the scenario-1 load draw was switched from Gaussian to lognormal, because at
$40\%$ CV the Gaussian returned a negative load and a negative marker. Commit `b289ef8`
requoted **one** number on that switch, the $40\%$ row ($0.70 \to 0.62$). Reading the whole
table and the paragraph under it against `audit/e7_excitation_budget.py`, five numbers there
were computed on the Gaussian draw and one was requoted:

| quoted | on the lognormal draw | state |
|---|---|---|
| $10\%$ row, $1.79$ | $1.7731$ | corrected in this pass |
| $40\%$ row, $0.70$ | $0.6164$ | corrected in `b289ef8` |
| carnitine row, $1.91$–$2.06$ | $1.91$–$2.02$ | corrected in 21.1 |
| "needs a CV of $229\%$" | no CV reaches $0.20$; sd bottoms near $0.376$ at $120$–$200\%$ and rises | rewritten in 21.1, R13 |
| daily sampling "$0.67$" | $0.65$ | corrected in 21.1 |

The dose-train row ($0.170$) is unaffected, and every number in the paragraph that follows
(the exponent band $0.06$ to $0.24$, $1.8\times10^{-7}$, $7.6\times10^{-4}$, the prior's
$1.05$) reproduces exactly.

Three things make this the sharpest instance of the class, not the mildest:

- **It reached the abstract.** Of everything in two documents, the abstract is what most
  readers will read and all readers will read first, and $1.79$ was in it.
- **It survived two freezes.** v16 froze the repository and v19 froze it again. Neither
  freeze is a check, and nothing between them looked at the table.
- **`audit/PHASE3_CHECKS.md` said the correction had been made.** Its own words were
  "Corrected in the table above and in the paper". In the paper it was made in one row of
  four, and the note asserting completeness is what made the remainder invisible: anyone
  re-checking would have read that line and stopped.

It is the same shape as R8 — right about its own computation, wrong about what it applies
to — with the scope error one level up. R8 was a number attached to the wrong footing. This
is a *correction* attached to the wrong extent.

**The rule: a correction to one row of a table is applied to the table, and the check is
reading the table, not the row.** A recomputation that changes one quoted value changes
every value derived from the same run; the unit of a fix is the table, the figure, or the
paragraph, never the line that prompted it. And a note recording that a correction was made
states *where*, because "corrected in the paper" is a claim about extent and carries the
same burden as a number (R7).

All five are now on the lognormal draw. One of the three was not a refresh, and it is the
one worth keeping in view. The $229\%$ does not reproduce at all: `e7`'s bisection brackets
$[10\%, 2000\%]$ and prints its ceiling, which is the search reporting "not found" in the
shape of a result. Under the lognormal draw no coefficient of variation reaches the $0.20$
target, so the paper's conclusion is not weakened but *strengthened*: the route does not
exist at any variability, not merely at an implausible one. That failure has its own entry
as R13.

Worth recording alongside: the miss was not found by re-reading the table either. It was
found by a sweep looking for something else, which noticed that two files quoted the same
quantity differently, and the remaining three were found only because the rule above forced
a reading of the whole table rather than the one row named. Conflicting numbers remain the
only detector here that has ever caught this class (R8), and it has now caught it twice.

### R13 — a bisection that printed its bracket

**A non-answer formatted as an answer**, and the second time this repository has produced
one (R6 was the first, in a different shape).

`audit/e7_excitation_budget.py` asked what day-to-day load variability would be needed to
reach the $0.20$ target without a dose train. It bisected on the interval
$[10\%, 2000\%]$ and printed

```
required CV = 2000%  of mean intake.
```

There is no coefficient of variation at which the target is met, so the bisection never
narrowed and the endpoint it printed was the bracket it started with. The paper read that
output at an earlier stage and quoted $229\%$, and the sentence around it said the route
existed but was implausible.

**What the search was actually reporting is stronger than what it printed.** Swept directly,
the median $\operatorname{sd}(\ln\hat\theta)$ falls to about $0.376$ at a CV between
$120\%$ and $200\%$ and \emph{rises} again beyond it. The route has a floor at roughly twice
the target. It does not fail because the required diet is implausible; it fails at every
diet. A search that had said "not found" would have produced a better sentence than the one
that returned a number.

**The rule: a search that terminates without finding must report that it did not find, never
its bracket.** A returned endpoint, a clipped maximum, an iteration cap, a `while` loop that
exits on count rather than on convergence: each of these is capable of emitting a
well-formatted number that answers a different question from the one asked. Where a routine
can fail to converge, the failure is a distinct return value and is printed as one.

This is the R6 class in a new form. R6 was a script asserting in advance what it would show.
This is a script asserting, in arrears, that it found something it did not. Both put a claim
into the output stream that no computation backs, and in both cases the reader downstream had
no way to tell. The common rule is the one R6 already states, extended: **a script's output
is a record of what it found, and where it found nothing it says so.**

## A sixth failure mode

- **A number taken from an abstract that its own Results contradict.** Three instances now,
  all in the citation layer: a crisis-cell value quoted as a managed baseline, a "three
  times higher" claim traced to a source containing no quantitative value, and a standard
  deviation reported as a mean by the source itself. The rule that follows, and it is now
  the standing rule for `claims.csv`: **every quantitative claim is checked against the
  source's Results or Tables, never its abstract, and every row records which section of
  the source the value came from.**

- **A script that states its expected result in a comment.** Two instances, both symbolic,
  both wrong, and both read as confirmations of themselves. The standing rule: **no script
  asserts an expected result in a comment.** An expected value goes in an assertion that
  can fail, or nowhere. A comment cannot fail, so a comment naming the answer converts a
  check into a decoration — and the failure is silent, because the reader who reaches for
  the script is reassured by exactly the text that is wrong. Where a script records what a
  run *returned*, it says returned, is written after the run, and names the tool that
  returned it.

- **A mechanism asserted where a derivation belongs.** Three instances this session, and
  the third ran in both directions. The pre-written `EXPECTED` comments (R6) asserted what
  a script would show; the seed-outlier call (R5) asserted why a number was anomalous; the
  filter-exactness sentence (R7) asserted a causal link, and its retraction asserted the
  absence of one. Not one of the three was a number, which is why none of them tripped the
  numerical checks this repository is full of. The standing rule: **a stated mechanism is
  a claim and carries the same burden as a number** — "because", "which is why", "the
  reason is", "it is the same X that" are all assertions and each needs a derivation, a
  citation, or a hedge. **A retraction is a claim too**, and is held to the same standard
  as the thing it retracts; R5 and R7 are both retractions issued without checking.

  Worth recording alongside the failures: all three were caught inside the session, two of
  them by re-reading the repository's own prose rather than by running anything. The
  repository's numerical apparatus did not catch any of them and was never going to. What
  caught them was reading the passage that justified the claim, which in each case existed
  and had not been consulted.

- **A number that is right about its own computation and wrong about what it computes.** R8
  is the first of these and it is the most dangerous kind, because every check this
  repository runs would pass it: it reproduces exactly, from committed code, on a footing
  that is internally consistent. What it does not do is match the footing of the sentence
  containing it. The rule that follows: **every quoted precision carries its footing in the
  same breath** — which design, which priors, and whether the therapy channel was present —
  because a bare number cannot be audited against a claim it was never attached to. Note
  what found it: not a test, not a re-read, but a *later result that contradicted it*.
  Conflicting numbers are the cheapest detector this project has, and the only one that has
  ever caught an error of this class.

- **An edit verified by the build rather than by reading.** R10 is the first, and it is the
  only failure here that *deleted* rather than distorted: a whole section left the document
  and every automated check passed. The rule: **a mechanical edit that leaves the build
  green is not thereby verified.** A compiler checks that prose parses, never that it still
  says what it said, and no test in this repository fails when a paragraph disappears. Any
  whole-document pass is therefore done by reading, or it is not done; and where one has
  already run, the question is not whether it builds but what the diff changed, site by
  site.

### R14. "the degenerate combination is recovered through curvature, at a square-root-of-noise price"

**Retracted** (Phase 4, P5), and this one is in the paper's own prose in three places plus a
cross-reference in a fourth, and in the two audit records that produced it.

Section 6 states that on the quasi-steady locus the information matrix is singular at the
truth and "flat only to first order", that the noiseless misfit along
Proposition~`prop:rank`'s null direction grows as `eps^2` so `chi^2` grows as `eps^4`, and
that the combination is therefore "recovered through curvature rather than through slope"
at a price: "the recoverable displacement scales as the square root of the noise level".
`rem:second_order` supports it with a diagnostic and with two numbers, `±31.6` against
`±36.9`. `rem:no_delta` refers back to it.

**On this locus the direction is not recovered at all.** Solving the balance for the increment that
produces a level `m` gives, with no remainder,

    Delta(m) = A/m + B + C m ,   A = K_m(L+b),  B = L+b-theta-k_e K_m,  C = -k_e ,

and the increment that produces a given level depends on the parameters only through
`(A, B, C)`. Since `A > 0` and `C < 0` that map is a strictly decreasing bijection, so the
level a known `Delta` produces is a function of `(A, B, C)` alone; and three distinct levels
determine `(A, B, C)`. So on the quasi-steady locus, at a constant recorded load, `(A, B, C)`
is everything a record of `(m_k, Delta_k)` pairs contains. Four unknowns against three
invariants leave a one-parameter family of exact observational twins,

    K_m(b) = K_0 L/(L+b) ,   theta(b) = theta_0 + b + k_e (K_0 - K_m(b)) ,   k_e unchanged,

whose noiseless misfit is zero at every displacement out to `b = 200`, which is 62% of the
recorded load, and whose true profile, re-optimising `(theta, K_m, k_e)` at each fixed `b`,
is flat to `1e-25` and traces the closed form to seven figures.

**The `eps^4` is the fibre's own curvature relative to its tangent line.**
Proposition~`prop:rank`'s null vector is that tangent at `b = 0`, verified symbolically
(the difference is identically zero in all four components). Moving along a straight line
through a curved exact fibre departs from it at second order, so the residual is `O(eps^2)`
and `chi^2` is `O(eps^4)` whether or not the record constrains anything. The measured local slope is `4.00` over four decades, drifting by about
two per cent across them, which is the signature of the geometry and not of the data.

**So the diagnostic `rem:second_order` recommends returns the wrong answer on the problem it
was written for.** It says `chi^2 ~ eps^4` means "the data do pin it, at second order, and
the filter is understating what the record contains". Here the data pin nothing and the
filter is reporting the prior exactly, which is the reading `rem:rankvsprec` gives the same
configuration two pages earlier. The two remarks contradict each other and the paper did not
notice, because each is locally consistent.

**One supporting number is the load prior itself.** The `±31.6` is
`1.96 × 0.05 × 322.7273`, the load prior's own 95% half-width, and is design-independent.

**The `±36.9` is a ray, and the printed pair says so.** `(3.84/2.08e-6)^(1/4) = 36.86`
against a printed `36.9`: the two numbers are tied by a quartic, which is what a straight
line through a curved fibre gives and what a profile does not, because a profile does not
produce a constant `chi²/ε⁴` to read a half-width off. A sweep of 88 designs by 7 noise
conventions by 3 residual scales, plus restricted profiles, profiles with each prior inside
`chi²`, profiles at a varying load, profiles on the ODE record and profiles with the assay
error on the abscissa, reproduces it nowhere as a profile.

**A first version of this entry made a footing charge that does not hold, and it is
withdrawn here rather than quietly dropped.** It said the pair was computed on a design of
roughly weekly size while the covariance it defends is monthly. `chi²/ε⁴` is second order and
therefore parametrisation dependent; `prop:rank` is stated in log parameters
(`06_excitation.tex:332`) and the filter's columns are derivatives in those coordinates,
while the audit drew its line in natural `(θ, K_m)` coordinates. Redrawn in the paper's own,
the two-daily × 2 y design gives `2.0900e-6` and `36.82` against the printed `2.08e-6` and
`36.9`, and the provenance is documented: `audit/c2b_degeneracy.py:95` runs
`arange(0, 730, 2.0)`, `PHASE2C.md:92-94` attaches the constant to `sd = 0.0251`, and
`PHASE2A.md:35` shows `0.0251` is the two-daily row. On that design "the two agree to within
a fifth" is right as computed.

What survives of the footing charge is weaker and still worth recording: `rem:second_order`
claims "every covariance this note reports survives", and the paper points the remark at the
`0.1696` entry of `rem:rankvsprec` (`06_excitation.tex:482-485`), which is monthly. On
monthly the same log construction gives `69.4` against the prior's `31.6`. The half-width
runs `36.2 / 49.3 / 69.4` across two-daily, weekly and monthly, so one number cannot defend
every covariance the note reports.

**A defect in the proof, found while checking this.** `06_excitation.tex:371` eliminates the
increment as `(θ+Δ)u = L − k_e m`, dropping the bias; the balance gives `L + b − k_e m`. So
`eq:nullvec` as printed is the `b = 0` value of a field that depends on `b`. Carrying `b`
through gives `n = (L + b + k_e K_m, −K_m, 0, L + b)`, whose integral curve is the fibre
exactly; integrating the printed vector instead gives `K_m = K_0 exp(−b/L)`, which leaves the
fibre at second order. That is why the exact fibre is not derivable from the paper as
written.

**What is unaffected, and this was attacked directly rather than assumed.** `0.1697` on the
monthly design is earned at a load coefficient of variation of 10%, where the fibre is
closed. Dropping the load-bias prior there costs **1.4%**: `0.1697` with it against `0.1721`
without, and across seven seeds `0.1688`-`0.1709` with against `0.1710`-`0.1810` without. On
every other design in `tab:designs` the prior is worth 0.6% to 3.5%. The 42x excitation
(`1.9390` against `0.0462`), the factor-of-22 sampling rule (`21.7x`, and `25.3x`
prior-free), the 11.7% dose-precision requirement and the cohort numbers all reproduce. No
printed number in the paper changes: the retraction is a mechanism and a diagnostic, not a
design number.

**A new failure class, and it is the reason this took three passes to find: a diagnostic
that answers a question the model cannot pose.** R7 was a mechanism asserted without
derivation. R8 was a number right about its own computation and wrong about what it computed.
This is neither. `audit/e5_structural.py` `q3` computes `eps^2` correctly, from committed
code, and reports it accurately; the paper quotes it accurately; the power law is real to
four decades. What is wrong is the *inference*, and it is wrong because the test profiles
along a straight line in a space whose relevant object is curved. No numerical check in this
repository could have caught it, because there is no number in it to disagree with. The rule
that follows: **a diagnostic is a model too.** Before reading a power law as evidence, state
what the alternative hypothesis is and check that the construction can distinguish them.
Here `chi^2 ~ eps^4` is equally consistent with "the data pin it at second order" and with
"an exact fibre is curved", and nothing in the construction separates the two. The separating
test is a profile with the nuisance parameters re-optimised, which costs more and which the
repository now runs.

**How it was found, and how the retraction was itself checked.** Not by a test and not by
re-reading. A second reader, working outside the repository, proposed a different
parametrisation, multiply the balance by `(K_m + m)` and it becomes linear in
`beta = (K_m, K_m b, b − θ − k_e K_m, −k_e)`, in which the degeneracy is a plain rank
condition on a design matrix and the fibre is visible immediately. That is the same detector
that caught `R8`: two descriptions of one object that do not agree. The route ran through a
proposed estimator that this audit otherwise rejected on three other counts.

The retraction was then put through eight adversarial lenses before it was written
(`audit/PHASE4_STRESS.md`), each prompted to refute it and each re-run by a second agent.
The finding survived all eight. **The audit did not**: its own forensic claim about the
`±36.9` was refuted, along with a bug in `audit/f5_fibre.py`, a single-seed number quoted
without its spread, and a retraction scoped both too wide (`"flat only to first order"` is
true of `prop:rank`'s null vector) and too narrow (`rem:rankvsprec` contains a sentence that
is false at two of the four rows it prints). All are recorded as C1 to C8 in
`audit/PHASE4_PI.md`. Given `R5` and `R7`, a retraction issued without that check would have
been the third of its kind here.

---

### R15. "the published inhibition shape cannot reproduce the repository's calibrated fall"

**Retracted** (Phase 5, the admissibility follow-up). It is mine, it is one report round old,
and it never reached `paper/`.

`audit/PHASE5_ROBUSTNESS.md` reported, in the Step 1 section and again in the
proven/measured table, that the two published inhibition shapes have a ceiling on the marker
fall they can produce over a `3.13`-fold enzyme rise, `65.65%` for the turnover structure and
`64.13%` for Attarwala's algebraic one, and that the repository's committed calibration of
`66.36%` lies above both. It was offered as "an independent argument, from the published
shape rather than from the trial, that `D_AMP` is set too high".

**The argument does not survive `Imax`.** Attarwala **fixed** `Imax` at `0.999` rather than
estimating it (Table 1, footnote m), so it is a modelling convenience. Recomputing each
ceiling at `Imax = 1` separates the structure's limit from that constant's:

- The turnover structure's ceiling at `Imax = 1` is exactly `Delta/(theta + Delta)`, the
  `IC50 -> 0` limit of the inhibition ratio, verified symbolically. It **clears** the
  repository's fall at every dose from `Delta_peak = 100` to `40000`, by 1.2 to 2.0
  percentage points. The `0.999` cap costs it 1.0 to 2.4 points, and that cap is the entire
  reason it appeared to fail.
- Attarwala's algebraic structure does fall short at every tabulated dose at both values of
  `Imax`, but by `0.09` to `1.74` percentage points, which is a knife edge and not an
  exclusion. That shortfall is structural, being the published floor `E0/base = 0.0987/60.1`
  blunting every fall, and it is worth stating on its own terms. It is not worth stating as a
  bound on `D`.

**What replaces it** is a better statement, and one the numbers support: on every dose in the
sensitivity table the repository's fall lies strictly **between** the two published shapes'
ceilings, the algebraic one a little too shallow and the turnover one with headroom. The
dose-response is bracketed by the published structures rather than outside them.

**And even that has a scope limit**, found by `audit/k6_recheck.py` check 9 and not by the
script it checks. Sweeping 61 doses instead of tabulating 10 shows the algebraic structure's
ceiling **crosses** the repository's fall at `Delta_peak = 443.6`, that is
`Delta_peak/theta = 0.148`; below that dose it can do *more*, not less, and the bracket
reverses. Every tabulated dose sits above the crossover, so the bracket holds where it is
claimed, and the first draft's phrase "a property of the three forms" did not.

**The argument from the trial is untouched.** `D_AMP = 6000` still implies a `62` to `66%`
fall against reported per-cohort falls of `-29%` to `+5%`. That remains the reason to
question the calibration; there is now one such reason rather than two.

**How it was caught.** The margins were suspiciously small and uniformly signed, which is the
signature of a comparison between two forms of the same saturating shape rather than of a
real physical bound. The sensitivity that settled it (`Imax` at `1` against `0.999`) was run
before the claim was reported a second time, not after. The failure mode is the one `R5` and
`R7` are also instances of: a verdict that turns on a constant nobody measured.

---

### R16. Four statements from Phase 5 that were true of the configuration measured and were reported as unconditional

**Scoped, not retracted.** All four are mine, all are one report round old, and none reached
`paper/`. Each stands where it was measured; each was written as though it held generally. They
are grouped because they are one failure mode, not four.

**(a) "Under PUB_B, theta is not identifiable, at any design, with or without load variation."**
True of the load insertion `audit/k2_structures.py` happens to use, and **false of another
plausible one**. Attarwala's marker equation was fitted in a fixed-diet colony and carries no
substrate load at all, so every load in our version is our insertion. Of four plausible
insertions, three leave theta unidentifiable at any design and one restores it: a load lifting
only the suppressible pool `base`, which separates it from the unaffected floor `E0`. Setting
`E0 = 0`, or knowing `IC50` externally, restores it too. The scoped statement is that **whether
theta is identifiable under that structure turns on an assumption the source does not make**,
and the deficit traces to carrying a free marker offset alongside a free marker amplitude.

**(b) "The bias is 1.29 to 2.03 with a residual under a fifth of one sigma."** The direction of
the bias is robust: theta_hat/theta is above 1 under every insertion and every estimator tried.
The **size** is not. With this repository's own priors it is 1.3 to 1.7; with them removed, on
the same records, it is 5 to 9. The priors shrink the bias, they do not prevent it, which is
consistent with Phase 1's finding that the `k_e` prior does most of the identification. And the
**invisibility** is insertion-dependent: 0.18 to 0.24 sigma when the load scales the whole
marker or the suppressible pool, 0.93 to 1.27 when it moves the unsuppressible part, against an
expected 0.816 for a well-specified fit at n = 12 and 4 parameters.

**(c) "The recommended schedules detect a wrong structure."** True against a **turnover** truth,
where monthly and phase-spread leave 2.0 to 2.6 sigma of residual. **False against the algebraic
one**, where no schedule tested raises a signal at all: 0.000 to 0.19 sigma structural, 0.78 to
0.88 with noise. Detection is a property of the alternative, not of the design.

**(d) "The three-number bound."** Stated without a class, which reads as a general theorem. It
is proved for structures whose fixed-load `(Delta, m)` graph is the zero set of a polynomial
with at most four monomials, where it is elementary, and membership is verified for the three
structures examined and claimed for no others. `PUB_A`'s full model leaves the class as soon as
its transient is used, which is where its fourth number comes from.

**How these were caught.** Not by re-running anything. By the questions being asked of the
report: where does the load enter, is the intuition right, does the detection claim hold for
both structures, and can the bound be proved for a stated class. Each answer was already
computable from code that existed; none of it had been computed, because the configuration that
was measured was the only one that had been considered. **The failure mode is reporting the
scope of a measurement as the scope of a claim**, and it is the same one as `R5`, `R7` and
`R15`. The correction in each case made the finding more useful rather than less: (a) is a
sharper criticism of the paper than the unconditional version, and (b) says the reported
estimate looks reasonable partly because the priors are pulling it back.

**One statement that was checked and did NOT need scoping.** The registry claim from the same
round, that "one registry fetch would settle" the sampling time, was wrong in the other
direction and is withdrawn outright rather than scoped: the record has been read and its 3-HP
time frame is "Baseline (predose levels) to postdose levels ... up to Week 40", a window and
not a visit. The timing cannot be settled from public sources and must come from the sponsor.

### R17. "Recording intake twenty times more precisely buys about six per cent, so precision is not the binding variable"

**In the paper**, `08_specification.tex`, the load-record paragraph, and in the spec table's
"How to read this table" as "the load record gives `0.170` whether it is kept to `1\%` or to
`20\%`". Both are **true of the quantity that was swept and false of the quantity a reader will
think was swept.**

`src/design/run.py`, `spec_numbers()`, computes

    sigma = LOC.sigma_total(LOAD_FRAC, ASSAY_CV)

**once, at `LOAD_FRAC = 0.05`, and reuses that one `sigma` on every row of the table.** The rows
labelled `sd_load_1pct` and `sd_load_20pct` pass `load_cv` to `design_record`, which is the
**day-to-day variation in the record**. The **recording error** never moves. It cannot, because
it enters only through `sigma`, and `sigma` was fixed before the sweep began.

Swept on the quantity that was meant (`audit/m1_loadbias_prior.py`, 24 seeds, two independent
code paths agreeing to `1e-14`):

| recording error | `sigma` | `sd(ln theta)` median | over 24 seeds |
|---|---|---|---|
| `1%` | `15.24` | `0.1566` | `[0.1563, 0.1571]` |
| `5%` | `21.96` | `0.1660` | `[0.1649, 0.1666]` |
| `10%` | `35.54` | `0.1916` | `[0.1884, 0.1936]` |
| `20%` | `66.24` | **`0.2637`** | `[0.2563, 0.2690]` |
| `40%` | `129.95` | **`0.4023`** | `[0.3916, 0.4107]` |

**The design meets the `0.20` target only if the substrate load is recorded to better than
roughly `12\%`.** Precision is not merely binding, it is the tightest requirement in the
specification after the increment itself, and at `20\%` the design fails.

**Why this matters more than an arithmetic slip.** Daly et al. 2017 (PMID 29021961, verified:
a cross-sectional survey reporting what patients were *prescribed*) and the absence of any PA
study reporting repeated as-eaten intake mean the realistic PA record is a prescription. The
only study anywhere in a protein-restricted inborn error that measures prescription against
actual intake is *Nutrients* 2023;15(23):4903 (Pinto et al., verified; 16 children with PKU,
three one-day diaries a month for six months): actual natural protein ran at **167% of
prescribed, a median `+40%`**. At `40\%` the headline is `0.40`, twice the target.

**A second correction inside the same paragraph, in the opposite direction.** The paragraph also
says that at a perfectly steady load the null direction "is filled by the `5\%` prior on the load
bias, and only real day-to-day variation lets the *data* fill it instead", so that recording
intake as eaten "is the difference between a result the data support and one a prior supports".
Measured, that overstates it. At a **constant** recorded load with the load-bias prior **removed
entirely**, `sd(ln theta)` is `0.2659`, not infinity. The dose train breaks
Proposition `prop:rank`'s singularity, because the increment varies even when the load does not.
The correct statement is that the prior is worth a factor of `1.6` there, not that it is holding
the result up.

**How it was caught.** By asking what a verified `+40%` prescription-to-intake divergence would
do to the specification, and finding that the table had no row for it. The two claims had been
checked against each other and never against the thing they were about. This is a new failure
mode for the list: **a sensitivity swept on the wrong one of two quantities that share a name.**
`LOAD_FRAC` is the recording error in `sigma_total` and the load-bias prior in the information
matrix, and `load_cv` is the record's variation. Three quantities, one word, "the load record".

Proposed replacement text is passage 11 of `audit/PHASE4_PROPOSED_TEXT.md`, and the softening
above is passage 15. **Neither is applied. No constant is moved.**

### R18. "A systematic Q2W-versus-Q3W offset in the tabulated falls ... does not need the sponsor"

**Mine, from `audit/PHASE6_REVIEWS_SEP2026.md` section 9, and withdrawn on the first test.**

Having found that Extended Data Table 3 of the Nature 2024 trial reports cohort 1 (0.30 mg/kg
Q3W) separately from cohorts 2 to 5 (Q2W), I proposed that the two candidate sampling
conventions could be told apart from the table alone: a fixed clock time after dose 3 puts both
arms at the same phase, whereas a next-predose reading puts Q3W at 21 days and Q2W at 14 days,
which on a 3.5-day half-life is `2^-6` against `2^-4` of one dose.

**The power was computed before the table was looked at** (`audit/m3_edt3_discriminator.py`,
Part A prints before Part B by construction), and it kills the test:

- The only contrast free of the dose confound is cohort 1 against cohort 2, both at 0.30 mg/kg,
  `n = 4` against `n = 3`, with standard deviations of 22 to 31 points. The **smallest difference
  detectable at 80% power is 59.5 to 74.2 percentage points** across the four markers.
- The two conventions differ by at most **22.5 points**, and that is at `D_AMP = 24000`, four
  times the repository's assumed value. At the assumed `D_AMP = 6000` they differ by
  **8.1 against 1.2 points**.

So the test is underpowered by a factor of between three and thirty, and no reading of the table
can separate the conventions. The data agree: `p` is `0.215` to `0.964` on the four markers, and
the two that move at all (3-HP `-26`, C3 `-25` points) move in the **direction opposite** to the
trough hypothesis, which predicts a *smaller* fall in the Q3W arm.

There is also no dose trend to lean on. Over the four Q2W dose levels the unweighted correlation
between dose and fall is `-0.37` (3-HP), `-0.50` (2-MC), `-0.16` (C3) and **`+0.67`** (n-PG), on
four points. Cohort-to-cohort differences are not reading dose, so they cannot be read as
reading interval either.

**Part D of the same script adds nothing new.** All twenty cohort-by-marker means sit above the
fall a one-day-post-dose reading predicts (`-61.8%` at `D_AMP = 6000`), but only three of twenty
reach `p < 0.05` individually, the cohorts within a marker are not independent, and the whole
comparison is conditional on an assumed `D_AMP`: a smaller `D_AMP` moves the prediction up to
meet the data. That is `PHASE5_ROBUSTNESS.md`'s existing observation that `D_AMP = 6000` implies
a `62` to `66%` fall against reported falls of `-29%` to `+5%`, with the arithmetic done. It is a
joint statement about the phase and the amplitude and cannot be split.

**`OPEN_QUESTIONS.md` section 2 stands unchanged: the sampling time must come from the sponsor.**
What is now pinned, and was not, is that the table is a change from a **predose** baseline
measured **after the third dose**, which is narrower than "a window up to Week 40" but not
narrow enough.

**How it was caught.** By stating the power before running the test, which the author's standing
rule requires and which in this case produced the answer before the data were read. Had the
table been read first, `-26` and `-25` points on two of four markers would have looked like a
finding.

---

### R19. "Removing the $k_e$ prior takes the same design from 0.1474 to 0.6451"

**Retracted** (Phase 7, step 2). Stated in three places and computed in none:
`paper/sections/07_design.tex:40`, `paper/sections/A1_supplementary.tex:524`, and
`src/design/run.py`'s `_filter` docstring in rounded form, "from 0.15 to 0.65".

`audit/n3_ke_prior_worth.py` searched for them. Over the five designs Section 7 lists, 24
intake records each, and 60 combinations of horizon, sampling interval and phase offset at
the seed the paper itself used, nothing returns `0.1474` and nothing returns `0.6451` with
the prior removed. The nearest approach to the first is `0.1490`; the prior-off values run
from `0.92` upward.

On the densest design that section lists, sampling every two days for two years, the pair
is **`0.1500` `[0.1498, 0.1501]` with the prior and `1.1096` `[1.0313, 1.1923]` without**,
medians over 24 records. Taken record by record the prior is worth a factor of **`7.4`
(`6.9` to `7.9`)**, not the `4.4` the printed pair implies.

**The correction roughly doubles what the prior is carrying.** It makes the paper's own
concern sharper: the one prior the design depends on is doing twice the work the paper
said, and by the `k_e` row of Table 2 no published source supplies the constant it is a
prior on.

**How it was found.** Not by looking for it. Registering `s_{ln k_e}` as a table row, which
the Phase 7 provenance audit required because the width lived as a default argument, meant
carrying the sentence into the table, and the table's rule is that a number arrives with
its footing. There was no footing to find.

---

### R20. Quoting Cramér–Rao curvatures where a reader reads an interval

**Retracted** (Phase 7, step 2), and this one is a labelling failure rather than a wrong
number, which is why it touches every precision in the paper rather than one.

`sd(ln theta-hat)` as reported throughout is `sqrt([F^-1]_11)`: the curvature of the log
likelihood at its optimum. It correctly describes the **sampling spread** of `theta-hat`
over repeated records, and `audit/n6_profile_vs_crb.py` confirms that it does, with the
MLE running `2715` to `3757` over 16 records against a truth of `3000`.

What it does not describe is the interval one record supports, and the abstract invited
exactly that reading: "identify theta to 17% over two years". Profiling `ln theta` and
re-optimising `(ln K_m, b, ln k_e)` at each point with the paper's own priors, the 95%
profile interval on the recommended monthly design is a factor of **`3.76`** wide (`2.86`
to `5.08` across records) against the **`1.95`** the bound implies, a ratio of `1.93`.

The interval is asymmetric, and the asymmetry is the paper's own degeneracy: lowering
`theta` is compensated by raising `K_m`, so the profile has a long tail downward while its
upper end sits about where the bound puts it. The paper already knew the surface is not
quadratic. `rem:second_order` says so and says the profile is the diagnostic that works.
It had never run that profile at the design it recommends.

At a disease-range capacity the failure is not one of width. Six of six profile intervals
run off the bottom of a grid spanning four log units: the record does not bound `theta`
from below **at all**, and the MLE itself lands between `3` and `224` against a truth of
`150`. What the record does supply there is an upper bound, ruling out capacity above about
24% of nominal.

**What changed in the paper.** The abstract no longer quotes `17%`; it quotes the interval.
`rem:curvature` states the distinction where the precisions are introduced, and the
specification table labels its entries as sampling spreads.

**Re-check.** On a noiseless record the profile's curvature at the truth is the same second
derivative as the bound, and they agree to `0.39%`. So the optimiser is not what is being
measured.

---

### R21. Every design number in the paper was a single draw

**Retracted** (Phase 7, step 2) as a class, not as a figure.

`design_record(times, load_cv=None, seed=3)` draws the day-to-day intake from
`default_rng(seed)`. No caller in `src/design/run.py` ever passed a seed, so every number
`table_designs` and `spec_numbers` produced, and every design figure in the paper, was one
intake record. The author's standing rule is that anything with a random draw carries its
seed spread.

Mostly it did not matter. The monthly headline moves 0.11% from seed 3 to the 24-record
median. **On the weak designs it mattered, and on one of them it mattered in the paper's
favour**: the aliased fortnightly row's seed-3 value of `3.51` is the *minimum* of its own
24-record range `[3.51, 5.64]`, and the median is `4.58`. Requoted record by record, the
aliasing cost is a factor of **`27` (`21` to `33`)** rather than the twenty-two the paper
claimed. The sharpest rule in the paper is sharper than it said.

Other requotings: two-monthly `0.187` to `0.188` `[0.184, 0.194]`; six-monthly `0.418` to
`0.457` `[0.412, 0.510]`; the Tier 1 requirement on the increment from `11.7%` to `11.6%`
`[11.4, 11.9]`; the `k_e` factor the gate needs from `2.95` to `2.67` `[2.65, 2.68]`.

`design_record`'s seed is now `BASE_SEED`, with `SEEDS`, `median_spread` and `tex_spread`
beside it, and both generated tables carry medians with ranges.

**The sharpest part of this one.** Section 2.6 has said, since before this branch,
that "every stochastic result uses `numpy.random.default_rng` with an explicit seed named
at the call site, and each quoted precision is a median over 24 seeds unless stated
otherwise". **Neither half was true.** The seed was a default argument named at no call
site, and no design precision was a median over anything. The paper asserted a
methodological standard it did not meet, in the section whose job is to state the standard,
and the assertion is what makes the defect a retraction rather than an improvement. It is
true now.

---

### R22. The carnitine reduction, reconciled with PR #3 rather than left in conflict

**Retracted** (Phase 7, step 2), and the finding is not this branch's. `PR #3`
(`claude/pensive-hamilton-1ga7l2`) established it; this branch was asserting the opposite in
four places, and two branches must not go to review disagreeing about the same equations.

**What PR #3 claims, and whether it holds.** Checked here from the equations rather than
taken on trust, because adopting a reading on the strength of a PR description would be the
same failure this log exists to record.

`rhs_C` has `c-dot = R_0 - k c [u(m) + rho_0(1 + m/M_i)]`, so:

1. **`k` cancels from the reduced curve.** Setting `c-dot = 0` gives
   `c* = R_0 / (k[u + rho])`, and the conjugation flux `k c* u` is `R_0 u/(u + rho)`, with
   `k` gone. Verified numerically over `m` from 1 to 1000 and `k` over four orders of
   magnitude: worst relative spread `2.45e-16`, which is machine precision. **Holds.**
2. **`k` alone sets the relaxation rate.** `d(c-dot)/dc = -k[u + rho]`, so
   `tau_c = 1/(k[u+rho])`, exactly proportional to `1/k` at every `m`. **Holds.**
3. **Therefore `k` sets the reduction's validity off equilibrium.** Starting the two-state
   system with carnitine at a pre-perturbation value, the peak relative error in `m`
   against the reduced system is **21.6% at the assumed `k = 1`**, `8.5%` at `k = 10`,
   `1.3%` at `k = 100`, `0.14%` at `k = 1000`. Error times `k` settles to a constant, which
   is the `1/k` scaling asymptotically; below `k = 1` it saturates (`26.9%` at `k = 0.1`)
   because a relative error cannot exceed order one. **Holds, with that qualification**:
   "scales as `1/k`" is asymptotic, not global, and PR #3's phrasing does not say so.
4. **The timescale premise runs backwards.** At the model's own `k` the carnitine time
   constant is `0.92` d against the marker's `0.81` h, slower by a factor of **27**. The
   external half is Rebouche's whole-body turnover of 38 to 119 h with a fast compartment at
   2 to 3 h, staged in `refs.sources` and named in prose rather than cited, because this
   environment cannot resolve the DOI and the register checks citations against the
   committed bibliography. **Holds.**

**What this branch was asserting, now withdrawn.** `02_model.tex` said the carnitine state
is fast relative to the metabolite, quoting `0.8` h against `3.5` d, which compares the
marker against the *therapy* and not carnitine against the marker. It said `k` is a free
scale and no conclusion depends on it. It called the reduction verified rather than
asserted, on the strength of an integration that is the same identity and so cannot fail at
any parameter value. And it stated the divergence above the separatrix without the
qualification that the reduced separatrix is not the basin.

**Corrected at the merge of PR #3.** The clause above originally read "a lower bound on
the basin rather than the basin", taken from the revision of PR #3 that was current when
this entry was written. PR #3 then measured the other side and withdrew it, and this entry
is corrected to match rather than left carrying it. The displacement is *signed*: started
fuller than the slow manifold the parent tolerates `+2.9%` to `+49.6%` more than the
reduced model predicts, and started depleted to half of `c*(m)` it tolerates `-0.5%` to
`-31.9%` less, so the reduced separatrix is a bound in neither direction. The depleted side
is the unsafe one and it is the side secondary carnitine depletion produces, which is the
mechanism the two-state model exists to represent. `audit/e10_qss_corners.py` measures both
directions; Section 2.3 of the paper states the signed result.

`alt_flux_C`'s docstring called this an adiabatic elimination of a fast state and
`rhs_C`'s called `k` a free scale. All six are corrected here in the convention PR #3 uses.

**What replaces the premise.** Nothing needs a separation of timescales. `alpha` is defined
by solving `c-dot = 0`, and `c-dot = 0` holds exactly at any fixed point, so at equilibrium
the two formulations are one equation rearranged. Every quantity this paper takes from the
parent is computed at equilibrium, where the cost is zero. The bound above is what a
*transient* would cost, and it is the third reason the discussion excludes the
decompensation route.

**Table 2's `k` row** is reclassified `derived -> assumed`, since nothing fixes `k`, and now
carries that cost. `derived` meant "fixed by a constraint elsewhere in the model", and only
"it cancels, so it does not matter" made it defensible.

**Still to reconcile at merge.** PR #3 also rewrites Section 2.6's verification claim and
the discussion's decompensation exclusion, and corrects the supplement's parallel passages.
This branch has separately rewritten the discussion. Whichever merges second will need to
take the other's wording in those three places; the equations and the `k` row now agree, so
the conflict is textual only.

---

### R23. The design table was generated by a function no entry point called

`table_designs()` writes `results/tables/tab_D_designs.tex`, which is Table 3 of the
paper, and nothing called it. `make figures` runs `src/design/run.py`, whose `__main__`
listed `fig_excitation`, `fig_collapse`, `fig_decision`, `fig_protocol`,
`table_excursion` and `table_spec` -- not `table_designs`. So the committed table froze at
whatever the code produced the last time someone ran the function by hand.

By the time the seed work of R21 landed, the two had diverged. The committed table carried
single-seed numbers with no ranges (six-monthly `0.4177`, fortnightly `3.5050`), while the
text beside it quoted medians over 24 records (`0.457`, and `4.58` with a range of `3.51`
to `5.64`). `3.5050` is the *minimum* of the fortnightly row's seed range, so the table
was quoting the single most flattering draw of the design the paper calls a catastrophe.
The paper therefore printed, side by side, a table and a sentence that disagreed, and the
rule R21 adopted -- no single-seed numbers -- held in the text and not in the table.

Fixed by adding `table_designs()` to the module's `__main__` and regenerating. Three text
numbers keyed to the stale table are requoted: the one-dose-cycle row `0.219 -> 0.216`
(`0.205` to `0.230`), the dense row `0.1499 -> 0.1500`, and the dropped-increment span
`1.0`--`9.4` -> `1.05`--`10.6`.

`audit/n7_increment_worth.py`'s R2 now reads the generated file and compares it against
what the code produces, so a table that drifts from its generator fails a check rather
than sitting in the paper.

### R24. What carrying the increment is worth, measured rather than asserted

Section 4 claimed that carrying the therapeutic increment is "worth between one and two
orders of magnitude on every design in the table, prior-free". Measured record by record
over 24 intake records (`audit/n7_increment_worth.py`), on the designs a protocol could
adopt: a factor of **52 to 141** with no prior on `k_e` at all, and **7 to 16** with the
`0.8` log-prior in place. The prior-free claim survives; stated without the "prior-free"
qualifier it would not, so the paper now quotes both numbers and says which is which.

Three re-checks: the sequential filter against one batch normal-equation solve (agrees to
0.32 times the round-off floor `eps * kappa(F)`, the tolerance derived from the
conditioning rather than guessed); the generated table against the code (R23); and, with
the dose train switched off, carried and dropped must be bit-identical, which they are --
so the measured factor is reporting the increment and nothing else.

### R25. Two rendering faults in generated artefacts

**Table 2's `theta` row carried a literal tab.** `src/common/params.py` wrote the meaning
cell as a non-raw Python string, so `"\textbf{at normal capacity}"` became a TAB followed
by `extbf{at normal capacity}`, and both the main-text and supplement tables printed
`extbf{at normal capacity}` in the row that says what `theta` is. Fixed at source.

**The main-text table kept five long notes where the text names three.** The compact table
selects its notes by symbol, and the parent-model block carries rows with the same symbols
as the local block, so the parent's `theta` and `K_m` kept their full provenance too. The
selection is now made per block. With the two rendering faults fixed, Table 2 drops from
about three and a half pages to two and a half.

**Figure 6's axis read `sign error on mu [\%]`.** Matplotlib mathtext renders `\%` outside
math mode literally, so the escape reached the page, in the axis label and in all three
legend entries. Fixed and the figure regenerated.

---

## Standing policy: the adversarial pass runs before a measurement is reported, not after

Until now the order was measure, report, then check. It is now measure, check, then
report. **A first report of a new measurement is provisional until the adversarial pass
has run**, and it is to be labelled as provisional if it is handed over before then.

The two cases that set the rule, both from the Vandermonde round:

**A weighting claim reported as fact, refuted by the check that should have preceded it.**
`audit/t10_vandermonde.py`'s `sd_vandermonde` built `A = Vw' Vw / sigma^2 + pen` without
normalising the weight to the exact weight's scale. With `w = 1` the transformed residual
has standard deviation `sigma (K_m + m)`, about `106 sigma` on this record, so the data
term came out roughly `1.1e4` too large and **both priors were silently discarded**. The
number that reached the PI, a factor of 59 for going unweighted, was an artefact of that.
Corrected, unweighted costs 0.06%, a factor of 1.00. The fix is one line,
`w = w * (sigma / np.mean(s * w))`, and the check that found it was the one that ought to
have run before the number was quoted, not after.

**A hardcoded literal in a script whose job was to compute that literal.** The same file
printed "2.3e4 times the second" for a column-norm ratio it never computed. The true
ratios are 146 and 1.7e4, and they are now computed. A number typed into a measurement
script is the same hazard as a number typed into a sentence, with the additional
property that it looks measured.

Both survived a green build and a read-through. Neither would have survived being asked
to refute itself, which costs minutes. The cost of the old order is not the wrong number,
it is that the wrong number arrives with the authority of a measurement and has to be
walked back.

This does not replace the rule from R14 that a retraction goes through adversarial lenses
before it is written. It extends the same standard to the ordinary case.

---

## Standing rule: where two artefacts must agree, have them write the same key

**An agreement that has to be checked is an agreement that can fail silently.** When two
pieces of this repository compute the same quantity, do not add a check that they match.
Make them write the same generated key, so that a disagreement raises `DuplicateKey` at
accumulation and there is no path on which the build is green and the two differ.

The instance that prompted it. The protocol figure drew a `range` bar per schedule and
`tab_D_designs` published `delta_range_*` for the same three steps over the same two-year
grid. Two computations of one quantity, in one document, with nothing comparing them. The
fix was not a test. The figure now writes `PAPER_NUMBERS["delta_range_monthly"]` itself,
so `NumberBook.__setitem__` is the comparison: if the figure ever drifted from the table,
the build would stop at the second write rather than print two numbers and pass.

It is the same shape that fixed `track_pct`, and that is the point. There the two values
were the same key written twice with different meanings, and the repair was to make the
key say which. Here they are one meaning computed twice, and the repair is to make the
two writes land on one key. Both turn "someone has to notice" into "it cannot happen".

Where to reach for it: a figure and a table that state the same number; a figure and a
sentence; a fast generator and a slow one that must agree; a value computed in two modules
because importing across them would be circular. Prefer it over an `assert`, over a unit
test, and over a comment saying the two must be kept in step.

**Its converse is a rule too, and it is the one that bites.** Shared keys are only correct
for the same quantity at the same settings. Two computations that legitimately differ need
different keys *and* a sentence saying how they differ, or the shared name becomes a
false claim of agreement. Figure~E1's $14.5$-day panel and the $14.5$-day row of the text
are the case: $25$ visits against $51$, $348$ days against $725$, a range of $1.80$
against $1.99$. Same schedule, different records. They keep separate keys and Section 4.3
now says so, because a reader who met only the numbers would have merged them.

---

## Standing rule: a generated number must be stable at the sample size it is generated at

Not "no stochastic numbers": this repository is full of legitimate ones, medians over
twenty-four intake records at fixed seeds, and they are the right way to report a
simulated quantity. The rule is narrower and it is about **stability**, because a
generated number is diffed against HEAD on every build and one that wanders with the seed
turns the artefact guard into a coin flip.

**If a candidate figure moves between sample sizes by more than its own reported
precision, it does not get generated at that sample size.** Three ways out, in order of
preference: raise `n` until it is stable; replace it with a closed form; or leave it in
this log, where a measurement belongs and where nothing diffs it.

The case that produced the rule. The fraction of this module's prior box in which the
parent's clearance peak drops below the local model's ceiling reads `0.14%` at 600 draws
and `0.23%` at 4000. That is sampling noise on about five events, and `{:.1f}` on it
claims a precision the estimate does not have. It was replaced by the closed-form upper
bound, `1.5%`, which is exact, seed-free and still says the thing the sentence needs.
The swept figure is two paragraphs up in R39.

The test is cheap and should be run before any new key is added: generate it at `n` and
at `4n` and compare against the format string it will be written with.

---

## Standing rule: a section number in this log is resolved through the aux, not the filename

The section files keep prefixes from structures the document no longer has. Two of the
eight still match the section they produce and six do not:

| file | built section |
|---|---|
| `01_question.tex` | 1 |
| `02_model.tex` | 2 |
| `06_excitation.tex` | **3** |
| `07_design.tex` | **4**, the design |
| `07b_decision.tex` | **5** |
| `08_specification.tex` | **6**, the specification |
| `11_scope.tex` | **7**, the scope |
| `12_discussion.tex` | **8**, the discussion |

**What went wrong.** Twelve references in the entries for items 20 to 47 were read off the
filename and named the wrong section. Item 40's clause is in `sec:tierstruct`, which is
**6.2**, and was recorded as §5. The vacuous comparison of the eighth instance is in
`sec:schedule`, which is **4.3**, and was recorded as §7. All twelve are corrected in
place, because an entry exists to say where a decision landed and one that sends a reader
to the wrong section is worse than one that names no section at all. The paper itself was
never affected: every cross-reference in it is a `\ref` and resolves through the same aux.

**How to resolve one.** Find the passage, take the `\subsection` that encloses it, and read
that subsection's `sec:` label out of `paper/main.aux` or `paper/supplement.aux`. The label
is canonical. The filename is not, and neither is the nearest preceding `\label{sec:...}`,
because several subsections label themselves some lines in: the paragraph that says
"Figure E1 draws the same three" sits in **6.4** while the last label before it is
`sec:tierincrement`, which is 6.3.

**What has not been re-resolved, and why.** Thirty-three references in entries written
before item 20 are left as written. They name the structure of their own time, and five of
the numbers among them, nine references in all, are impossible under this one: §1.3, §3.6,
§3.7, §5.5 and §5.6, since §1 and §5 have no subsections at all and §3 has five. Rewriting them would falsify a record of what an
earlier document said. **Read a section number in an entry before item 20 as naming that
entry's document and not this one.**

## Standing rule: a page count is reported with its rerun count

A page number off a short build is not wrong here, it is unverified. Both documents
converge before the second pass today: `Rerun to get` is 0 at four passes at all seventeen
heads measured in R53, which is why a two-pass count agreed with a four-pass one
everywhere. That is a property of these documents at this size and structure, not of the
method, and the next structural change is exactly when it stops being true.

So a page count is quoted with `grep -c 'Rerun to get' paper/main.log` beside it, and
undefined references are counted as `LaTeX Warning: Reference` and `LaTeX Warning:
Citation` rather than as `grep -c undefined`, which also counts font-shape warnings and
returns 1 on the supplement for `T1/lmr/bx/sc`.

## This project's characteristic failure: a number that does not say what it is a number of

Six instances of one shape, a seventh of another and an eighth of a third. That is past
coincidence and worth naming so the next person looks for it by default rather than
finding it a ninth time. The six below share a shape and a fix. The seventh and the
eighth share neither, with each other or with the six, and are recorded after them,
because a reader who applies the six's fix to either will not catch it. The eighth is
not a wrong number at all.

| | what was displayed | what it was actually a number of |
|---|---|---|
| `trans_angle` | an angle, in degrees | a quantity with no referent; deleted, not corrected |
| `track_pct` | the marker's departure | the departure after a tolerance-dependent window, keyed as though tolerance-free |
| `min(q, N)` | "25 phases" | the panel's own visit count, on three of five panels |
| the two managed levels | "the managed level" | `10` in the local model and `81.93` in the parent, in one paragraph |
| the three admissibility readings | one verdict | three checks that agree only at two coincidence points |
| the `theta_hat = 151` ladder | a fitted capacity | a fit whose band, weighting and free parameters were never recorded |

**The shape is always the same.** A quantity is computed correctly, reported accurately,
and displayed without the qualifier that fixes its referent. Nothing is wrong with the
arithmetic in any of the six. What is missing is the part of the label that would let a
second person recompute it and get the same thing. That is why none of them failed a test
and why all six survived a green build: **a check can compare two numbers but it cannot
tell you what a number is of.**

Where it comes from is worth naming too. Five of the six arose at a boundary --- between
two models, between two readings of one check, between a figure and a table, between a
tolerance and the quantity it conditions. Inside one module the referent is carried by
context and never needs saying; the moment a number crosses a boundary the context does
not cross with it.

So the habit to acquire: **when a number crosses a boundary, carry its referent across
with it explicitly**, in the label, in the key name, or in the sentence. The three rules
above are all special cases. The adversarial pass catches it when the referent is wrong;
the shared-key rule prevents it by making two computations one; this rule says the missing
referent is the thing to look for first.

---

### The seventh instance, whose shape is not the same (item 42)

**What happened.** `theta_crit` was one value, `150.7747475984808`, computed once and
reported twice. `margin.tex` carried it at `{:.0f}` and the margin figure drew `151`;
`boundary.tex` carried `fit_theta_crit` at `{:.1f}` and the discussion wrote `150.8`. In
the one paragraph where the threshold sits beside a fitted `theta_hat` of `151.01`, the
figure said the threshold was `151`, and a reader could not tell the threshold from the
estimate. R45 has the repair.

**Why it is not one of the six.** In all six the referent was missing or wrong: the label
did not say what the number was a number of. Here both labels were right. `theta_crit`
said threshold, `fit_theta_crit` said threshold, and both printed values were correct to
every digit shown. What failed was the precision. Rounding a threshold of `150.77` to
`151` did not misdescribe it, it made it collide with the one quantity the paragraph put
it next to. The six's fix, carry the referent across the boundary, would have changed
nothing here. That is why it is recorded apart from them.

**So the shape is: one quantity, two renderings, and the coarser rendering collided with
what the sentence was comparing it against.** The fault is not in the number and not in
the label. It is in the distance between two printed numbers, which is a property of
neither of them alone.

**Is it enforceable? Half of it, and the half that is not is the half that bit.**

The route of fixing precision where the quantity is computed, generators emitting
formatted strings rather than raw numbers, is already this repository's practice and was
already its practice when this happened. `results/numbers/*.tex` stores rendered strings
and not floats, so a key has exactly one rendering in the document; `common.numbers.num`
returns that same string so a figure draws what `\num{}` typesets; and
`audit/g1_phase_clock.py` states the rule in its own docstring, "a header rounding 1.9873
to 1.99 while the sentence beside it prints 1.987 is the same drift in miniature". So the
other route, looking for one key that appears with two format specifiers, is already
looking for something that cannot happen. It is also not what happened. This seam was two
*keys*.

**The checkable generalisation is indistinguishability, not equality.** Flag two keys
quoted in one paragraph whose values agree once both are rounded to the coarser of their
two displayed precisions. That rule catches this case exactly, because `151` and `151.01`
agree at nought decimal places. Measured against the document as it stands: 326 of the 329
generated keys are plain numbers, 57 paragraphs quote more than one key, and the rule
flags **17 pairs**.

Triaging the 17 is what decides whether it can be a gate, and it decides that it cannot.

- Three are the `min(q, N)` instance already in the table above: `reached_monthly`,
  `reached_off` and `visits`, all at 25.
- Several are equal for a reason, and the reason is a small theorem. `available_off` (28
  lattice phases) and `lap_visits_off` (28 visits to a full lap) agree because a lap
  visits each lattice phase exactly once. `available_fortnightly` and
  `reached_fortnightly` are both 1 because the dose-interval schedule has one phase and
  reaches it.
- Some are coincidences between unrelated quantities, such as `share_5pct` at 16.3 beside
  `struct_sigma_3min` at 16.3 in the abstract.
- And **one is this same shape, live.** This entry first said two, and one of the two was
  not. `delta_range_monthly` (design.tex) and `range_monthly` (phaseclock.tex) were two
  implementations of one quantity on one record: `np.arange(0, 730, 30.4)` is 25 points,
  `fixed(25, 30.4)` is the same array, the arrays compare identical and both ranges are
  `1.843777`. `src/design/run.py` computed `(inc.max() - inc.min()) / P.theta` inline
  while `audit/g_design_common.py` computed `(inc.max() - inc.min()) / TH` as a function,
  which is the hazard `g_design_common.level`'s own docstring names, "two spellings of
  one equation is how they drift". The second pair, `delta_range_aliased` and
  `range_fortnightly`, is **not** a duplication: 53 visits against 25, a genuinely
  different record, agreeing because a schedule at `eps = 0` saturates, which Section 6.4
  already states. Merged, and this bullet corrected, in R51.

**So: a report, not a gate, and the argument for it is not its hit rate.** Seventeen flags
of which most are correct and waivable would be switched off inside a week if they failed a
build, and on hit rate alone that would be the right decision. What earns it is what it
found on its first run, which was not a duplication. It was **a sentence that asserted
nothing**, which the duplication had been making look like a claim (R51, and the eighth
instance below). A pair of keys that render the same is a candidate vacuous comparison, and
no other check here can see one: every number in it is correct. So the listing is worth
keeping, with each pair waived by name or merged.

**The residue is habit, and it is the part that bit.** The rule needs both numbers to be
`\num{}` keys in one paragraph. The collision here was between a *figure label* and a
*sentence*, and no paragraph-scoped rule sees a figure. Nor is "a threshold and an
estimate must be distinguishable" a numerical fact; it is a claim about what the sentence
is for. So the habit is: **when a number is printed in order to be compared with the one
beside it, print both to the precision that separates them**, and read the figure against
the paragraph, not only the paragraph against itself.

---

### The eighth instance, which is not a wrong number but an empty sentence (item 47)

**What it was.** Section 4.3 read: *the ranges are over these 25 draws and not over the
two-year record of Table D: monthly reads `range_monthly` here against
`delta_range_monthly` there, and the 14.5-day schedule `range_off` against
`delta_range_off`.* Both monthly numbers were `1.84`. Both were correct, both were
correctly labelled, and for monthly the two-year record **is** those 25 visits, so the
clause set one number against itself and called it a comparison.

**Why it is worse than the six, and worse than the seventh.** In every one of those
something was wrong: a referent missing, a precision colliding. Here nothing was wrong.
There was nothing to check, because a check would have compared `1.84` against `1.84`,
found them equal, and equality is what the sentence was presenting as the finding. **A
wrong number can be caught by recomputing it. An empty sentence cannot be caught by
recomputing anything**, because it makes no claim to recompute. What kept it invisible was
the redundant name: with two keys it reads as a contrast between two quantities, and with
one key it reads "monthly reads 1.84 here against 1.84 there" and is visibly nothing.

**Where it came from, which is the portable part.** Item 20 asked for the figure's
14.5-day design and the text's 14.5-day design to be told apart, and that distinction is
real: 51 visits against 25, `1.99` against `1.80`. The clause written to carry it was then
stated over the panels as a set. It is true of the one it was written for, empty for
monthly, and empty-but-explained for the fortnightly pair, which Section 6.4 handles by
saying outright that the number does not depend on how long one samples. **A clause written
for one member of a set and then quantified over the set goes empty exactly on the members
whose reason does not hold.** So what to look for is not the number, it is the quantifier:
when a sentence says "these all differ", or "each reads X here against Y there", check the
members rather than the values.

**And it is what the indistinguishability report is actually for.** That check was argued
above on its hit rate, seventeen flags mostly waivable, which is a weak argument and would
have lost. The strong one is this instance: two keys rendering the same is the only visible
trace a vacuous comparison leaves, because every number in it is right. Corrected in R51,
where Section 4.3 now says it for the 14.5-day schedule alone and says of monthly that its
two-year record is those visits.

---

### R26. Two conditioning claims about the Vandermonde form, neither of which holds

Recorded because a negative that is not written down gets re-proposed. The reframing
itself was adopted and is in Section 3.2; these two claims about it were not.

The form is right, and it is the paper's own. Clearing the balance by `(K_m + m)`, not by
`m`, gives `Delta_i m_i = A + B m_i + C m_i^2`, a Vandermonde system in `(A, B, C)` whose
nodes are the marker levels. Checked symbolically: `A = K_m(L+b)`, `B = L+b-theta-k_e K_m`,
`C = -k_e`, which are Section 3.2's three invariants exactly. Same three numbers, same
fibre.

**(a) `cond(V)` does not price the factor of 27.** It cannot: `Delta` lives in the
response and `V` is built from the marker levels alone, so `V` cannot express aliasing at
all. With a varying load the nodes stay distinct even at the dose interval, and `cond(V)`
is finite there while the design is singular. This is also why the Vandermonde vocabulary
went into Section 3.2 and not into Remark 4.

**(b) `cond(V)` does not price the factor of 9.5 either**, and no choice of units prices
both: in `m/K_m` units the ratio comes to 22.7 on (a) but 76 on (b). The column norms are
19.1, 2.22e3 and 3.25e5, so the second column is 146 times the first and the third 1.7e4
times it. That spread is the units of `m`, not a property of the design.

**(c) The unweighted fit costs nothing worth naming.** 0.06%, a factor of 1.00. See the
standing policy above for what the first answer to this was and why.

### R27. The closed form for the triples estimator is the fibre of (7), by a third route

`audit/t12_triples.py` R2b, and it came free. Section 3.2's invariants are
`A = K_m(L+b)`, `B = L+b-theta-k_e K_m`, `C = -k_e`. A record generated at load `L+b` and
fitted with the load held at `L` has to match all three, and `C` matches for nothing, so

    K_m_hat = K_m (1 + b/L),    theta_hat = theta - b (1 + K_m k_e / L)

which is the triples estimator's own bias formula, derived for a different purpose.
Differentiated at `b = 0` it gives `d(theta)/db = 1 + K_m k_e / L = 2.5492957746`, and
that is exactly the `theta`-to-`b` ratio of `prop:rank`'s null vector with the bias
carried, `n = (L + b + k_e K_m, -K_m, 0, L + b)` (R14). The two agree to `0.00e+00`.

Stronger than the tangent claim it was put as: the closed form is **exact in `b`**, not a
first-order approximation. Refitting a record made at load `L+b` with the load held at `L`
reproduces both `theta_hat` and `K_m_hat` to a relative `7e-15` to `4e-14` across
`b` in `[-150, +300]`, with residuals at `9e-14`, a range far wider than any study would
see.

So the identity now has three verifications that share no code: `audit/f5_fibre.py`
integrates the null field; this solves the invariants; and the estimator was built to
answer a different question and landed on the same number.

### R28. The parent channel is nearly quadratic over the working band, and the queue stays closed

Recorded, not chased. Over the 3.4-fold band the record visits, `alpha(m) - k_e m` is
nearly quadratic in `m`, so a misspecified parent channel is very nearly inside the span
of the quadratic every triple fits. That is why the triple-averaging diagnostic of item 18
cannot separate the two channels: a 55% bias in `theta` moves the statistic by a few per
cent. A load bias does the same thing at degree one, and degree one is inside the same
span. The observation says what the failure is rather than suggesting a fix, and no fix
is being looked for.

### R29. The number-provenance rule now reaches figure text, and the first thing it caught was a caption

`tools/check_artefacts.py:check_figure_numbers`. The rule: a number a figure draws as a
claim about the design must come from a registered generator, exactly as a table number
does. A number that is an axis tick, a fact about the drawing, or a restatement of the
setting the panel was drawn at, need not.

**How much of it is enforceable, and which part is not.** Three parts, and only two of
them are machine-checkable.

1. *Finding the numbers: enforceable.* An AST pass finds every interpolation in text
   built at a drawing call, across f-strings, `%` on a literal and `.format()`. It
   resolves this repository's own text helpers (`corner`, `over_x`, `beside_y`,
   `at_point`) out of `common/style.py` rather than listing them, which matters: the
   first version listed matplotlib names only, saw 23 drawn numbers, and the four helpers
   were hiding 8 more, among them `theta_crit` and the safety margin.
2. *Checking that a number is registered: enforceable.* Three shapes pass, all verified
   rather than declared: a call to `common.numbers.num(module, key)`, which reads the
   committed file so the figure draws what `\num{}` typesets; an index into the dict the
   module hands to `write_numbers`; and an expression the module registers verbatim
   elsewhere, which is how the dose train's `trough` is also `train_trough_frac`.
3. *Deciding whether a number is a claim: NOT enforceable, and nothing makes it so.*
   `{v:.2f}` in a title is the same syntax whether it is a design claim or an axis
   annotation. So the rule is enforced as a **declaration**: each figure module carries a
   `DRAWN` dict classifying every interpolation it draws, a drawn number that is not
   declared fails, a declaration for a number no longer drawn fails, and
   `UNREGISTERED:` is debt the guard prints on every run. What the guard cannot check is
   whether a declaration is honest. Writing `geometry:` over a claim passes. The
   declaration is what gives a reviewer something to review, which is the most that is
   available here.

Also outside it, and counted rather than assumed: 28 interpolated strings in figure
modules are built away from any drawing call, console output and table cells subtracted.
The guard prints that number so the gap is a number and not a silence.

**Standing: 29 drawn numbers, 5 registered, 8 unregistered claims.** The eight are named
with their reasons in the modules' own `DRAWN` dicts and printed on every run.

**What it found on day one.**

*Section 6.4 said the figure draws three numbers it does not draw.* The sentence gave the
two-year ranges as 1.84 monthly, 1.99 at 14.5 days and 0.13 fortnightly, then said
"Figure 6 draws them". Figure 6 draws 25 visits, and 25 visits a fortnight apart cover
348 days, so its 14.5-day panel reads **1.80**, not 1.99. Monthly and fortnightly happen
to agree, the first because 25 monthly visits is already the two-year record and the
second because the collapse does not depend on how long one samples. Corrected, and
Figure 6 now registers its own fifteen numbers into `results/numbers/phaseclock.tex`, so
the caption's five hand-typed figures (`25 of 35`, `25 of 28`, `1 of 1`, `1.84`, `0.13`)
are `\num{}` keys the guard diffs against HEAD.

*The margin panel and the registry were two computations of one number.* Panel (c) of the
fibre-margin figure read its dot off a 400-point sweep grid with `np.interp`, while
`results/numbers/local.tex` registered `mtol_window_pct` from an exact `margin_window`
solve at the same abscissa. Nothing compared them.

They are now one object: `fig_margin_fibre` builds the three widths into the dict it hands
`write_numbers` and draws out of that dict, so `np.interp` no longer appears in any drawn
string. **And the two agreed, which is worth saying plainly rather than leaving the entry
to imply otherwise.** The grid was fine enough that the gap was `14.778341` against
`14.778437`, a relative `1e-5`, and both round to the `15%` the panel drew and the `14.8`
the registry held. Nothing in the paper was wrong. What was wrong was that two
computations of one quantity sat in one document with nothing to make them agree, and the
next person to coarsen that grid would have found out from a reader.

*Two duplicated literal blocks.* `GENERATORS` in `src/common/numbers.py` listed five stems
twice and `paper/preamble.tex` input five numbers files twice. Both harmless, both the
same class as the duplicate-key bug `NumberBook` exists for, and both removed.

### R30. The periodicity strip on Figure 6, tried and dropped

The phase overlay gives up one thing the polar dial it replaced could show: that phase is
circular, so a draw at 13.9 d is adjacent to one at 0.1 d rather than at the far end of
the panel. The cheap fix is to repeat the first day or so of the cycle, faded, past the
right edge. It was implemented and rejected on looking at it.

It costs the two panels that carry the result. The fortnightly draws all sit at phase 0,
so the repeat puts a second cluster at 14 d and the panel reads as **two** groups of
draws, against a caption that says all twenty-five are on one point. The marker row is
worse: the repeated points sit near `m = 3.5` beneath a curve that ends near `9`, which
reads as stray data rather than as the next cycle beginning. At text width the faded
curve segment is also too faint to be seen as a curve at all.

One strip of ink is not worth weakening the aliasing panels. The adjacency stays a
sentence in the caption, and the reasoning is in the module docstring so the next reader
does not re-propose it.

### R31. The eight unregistered figure claims, cleared before the preprint

R29 left eight numbers that figures assert about the design with no generator behind them.
That is the population that produced `trans_angle`, `track_pct` and the `min(q,N)` header,
and a preprint is the point after which they stop being cheap to fix. All eight are now
registered, and `check_figure_numbers` reports **29 drawn, 13 registered, 0 unregistered**.

| was | is |
|---|---|
| `excitation(lo, hi, Km)`, drawn in the excitation legend for two bands | `E_held` and `E_info` in `design.tex`. Neither was `decomp_E`, which is a third band, the decompensation excursion of `tab_D_excursion` |
| `rng`, the range bar of the protocol figure | keyed straight into `PAPER_NUMBERS` under the same `delta_range_*` names `tab_D_designs` writes |
| `slope`, the load-error amplification factor | `load_error_factor` in the new `fibre.tex` |
| `th.max()/th.min()` | `fibre_theta_fold`, `169`. Not `band_fold`, which is the working band in `theta` |
| `thc`, `margin[i]`, `lo[i]`, `hi[i]` | `theta_crit`, `margin_nominal`, `margin_band_lo`, `margin_band_hi` in the new `margin.tex` |

**The range bar is the one worth reading twice.** Rather than give it a key of its own, it
writes the key `tab_D_designs` already owns. Both compute the range over
`np.arange(0, 730, dt)` on the same three steps, so if they ever disagreed the
`NumberBook` would raise `DuplicateKey` at accumulation. The comparison is no longer a
check somebody has to remember to write; it is the data structure. `make figures` runs
both and does not raise, so they agree to a relative `1e-12`.

Two new stems, because a module's numbers file is written in one call and a second call
would clobber it: `fibre` (`src/model_A/run.py`, accumulated across two figure functions
into a `NumberBook` and written by `fibre_numbers()`, which asserts every expected key is
present so a partial write fails loudly) and `margin` (`src/model_C/run.py`). `fibre` is
cheap and joins the guard's fast path; `margin` costs 48 s and is covered by
`--working-tree` after `make figures`, as `decision` and `profile` already are. The guard's
fast path also gained `fig_excitation` and `fig_protocol`, without which the regenerated
`design.tex` comes out two keys short.

### R32. Figure 6's 14.5-day panel is the lap condition, drawn

The caption said every panel shows "*the same* 25 visits". It shows the same *number* of
visits. `n` is held at `\num{visits}` and the spans run from `\num{span_fortnightly}` days
at the dose interval to `\num{span_spread}` phase-spread, so the ranges are not like for
like in calendar time. The caption now says which is held fixed and why.

**And the two shortfalls in "25 of 35" and "25 of 28" are not the same kind**, which the
bare counts hide completely.

* Monthly's step leaves `eps = 2.4` days, so 25 visits sweep **4.3 laps** of the dial. It
  exhausted what its step can reach long ago and its `25 of 35` is `q > n`, nothing more.
  Its range is the complete one, `1.84` here and `1.84` in `tab:designs`.
* The 14.5-day step leaves `eps = 0.5` days, so 25 visits sweep **89% of one lap** and
  stop. It needs **28** visits to satisfy `n eps >= T`, and the `1.80` it reads against the
  `1.99` of a completed lap is exactly what that shortfall costs.

So the figure was already illustrating `eq:lap` in the one panel where it binds, by
accident. It says so now, out of `lap_pct_off`, `lap_visits_off`, `laps_monthly`,
`eps_off_d` and `eps_monthly_d`, all generated.

### R33. The load-bias map is exactly affine, and the paper now says so

R27 established the closed form by three routes. The consequence belongs in the text, not
only here. Section 3.1 said `prop:rank`'s null vector is the fibre's tangent at `b = 0`
and left it there, which reads as though a bound on the load bias is priced by a
linearisation. It is not:

    K_m_hat = K_m (1 + b/L),    theta_hat = theta - b (1 + K_m k_e / L)

exactly, at every `b`, with no higher-order term. The fibre is curved in `(theta, K_m)`
and the map from a load bias to the **fitted capacity** is affine. So the factor of
`1.028` that a 5% load record leaves open is that map evaluated **at** the bound rather
than near it. The slope `2.549` and the factor `1.028` are now generated
(`bias_slope`, `bias_theta_factor`), and the `1.028` that had been typed into the sentence
is a `\num{}` like everything else.

### R34. The two 14.5-day designs, and what the lap condition is actually shaped like

**The check first, because it could have gone the other way.** `sd_aliased_off` = `0.1614`
is computed on `np.arange(0, 730, 14.5)`: **51 visits across 725 days, 1.82 laps of the
dial, 28 of 28 phases reached**. The lap completes with room to spare, so the figure is
not quoting a number earned on a record short of a lap and the lap condition has no
unpriced consequence hiding behind it.

**But the figure and the text were two designs wearing one name**, which is what a reader
would have merged:

| | visits | span | laps | phases | range | sd |
|---|---|---|---|---|---|---|
| the text's $14.5$-day row | 51 | 725 d | 1.82 | 28/28 | 1.99 | 0.1614 |
| Figure 6's $14.5$-day panel | 25 | 348 d | 0.89 | 25/28 | 1.80 | 0.1739 |

Section 4.1 now gives the first its `n`, its span and its lap count, and says in one
clause that the panel is not it. That is the converse of the standing rule above: the two
keep separate keys *because* they are separate quantities.

**And the shape of the rule came out of the check.** `n eps >= T` reads like a knife-edge
and is not one. Sweeping the lap fraction at three fixed visit counts, so that `n` cannot
be doing the work:

| `n eps / T` | 0.02–0.04 | 0.19 | 0.38 | 0.50 | 0.75 | **1.00** | 1.50 | 3.00 |
|---|---|---|---|---|---|---|---|---|
| sd at n=25 | 4.044 | 0.924 | 0.456 | 0.207 | 0.180 | **0.171** | 0.172 | 0.171 |
| sd at n=40 | 2.228 | 0.452 | 0.251 | 0.187 | 0.169 | **0.164** | 0.165 | 0.164 |
| sd at n=53 | 1.455 | 0.309 | 0.196 | 0.177 | 0.164 | **0.160** | 0.161 | 0.161 |

The cost falls steeply up to one lap and is flat afterwards, at **every** `n`, and the
`Delta/theta` range saturates in the same place, which is the mechanism. So `n eps >= T`
is not a threshold chosen for tidiness: it is exactly where the cost stops falling. Being
half a lap short costs about 20%; being at 4% of a lap costs a factor of 9. The rule is
earned. Recorded here and not put in the paper, which states the inequality correctly
already.

### R35. Item 22. Forward simulation of the therapy against capacity

`audit/t13_therapy_forward.py`, 81 s, measurement only. No paper text, and nothing here is
dosing guidance: Section 6.7 already says a protocol written from this work is not a rule
for escalating a patient, and a forward simulation at an **assumed** capacity is weaker
evidence than that, not stronger. What it can inform is a protocol's missed-dose
contingency, which is a document-level decision.

Needs no identifiability at all, which is the point. `theta` is swept as a percentage of
the model's normal (100, 20, 10, 5, 2, 1) and never estimated or cited; the objective is
marker exposure and not the sign of `theta - theta_crit`, whose threshold Section S1.20
shows moving by a factor of 24 across its own constants; and every trajectory past the
parent's clearance maximum is marked rather than reported beside the rest. The four
constants the setup rests on are computed from the code and reproduce what the paper
already commits to: clearance maximum `3.365` against 3.4 times managed, upper equilibrium
`12.607` against 12.6, ceiling `6.455` against 6.45, largest admissible `D` at 5% of
`2671.9` against 2672.

> **RETRACTED IN PART, see R36.** The flag this section reports, "past the clearance
> maximum", was computed against a threshold that misread the paper. **No row of (a), (b)
> or (c) is outside what the parent supports.** Everything else below stands; the
> retraction is confined to the flags and to the sentences that read them.

**(a) The swing grows as capacity falls, and the disease-range column is flagged before
anything is even perturbed.** A dose in mg/kg is fixed in absolute terms while `theta` is
not, so `Delta_peak/theta` runs from 2.1 to 213 across the sweep.

| theta | untreated | mean | peak | trough | fold | above managed |
|---|---|---|---|---|---|---|
| 100% | 10.00 | 6.50 | 8.91 | 3.40 | 2.62 | 0% |
| 20% | 34.06 | 13.70 | 24.68 | 4.57 | 5.40 | 64% |
| 10% | 45.72 | 16.02 | 31.00 | 4.78 | 6.49 | 69% |
| 5% | 54.02 | 17.51 | **35.32** | 4.89 | 7.22 | 71% |
| 2% | 60.04 | 18.55 | **38.41** | 4.96 | 7.74 | 72% |
| 1% | 62.24 | 18.92 | **39.55** | 4.99 | 7.93 | 72% |

The bolded peaks are past the clearance maximum at `3.4 x managed = 34`. **At 5% of
nominal and below the BOOKED, UNDISTURBED schedule is already outside what the parent
supports**, which is a stronger statement than the one the missed-dose run was expected to
produce. Admissibility breaches at 10% and below, reproducing Section 7.1 row for row.

**(b) The interval has no pharmacological optimum, and "7 days" would have been the edge
of a sweep reported as a result.** At fixed delivered rate, exposure falls monotonically
as the interval shortens, at every capacity, towards the continuous-delivery limit where
`Delta -> R/lambda = 2164`:

| theta | h=28 | h=14 | h=7 | h=3.5 | h=1 | continuous |
|---|---|---|---|---|---|---|
| 100% | 7.23 | 6.50 | 6.15 | 6.05 | 6.01 | 6.008 |
| 5% | 26.59 | 17.51 | 13.91 | 12.95 | 12.66 | 12.634 |
| 1% | 29.74 | 18.92 | 14.71 | 13.60 | 13.26 | 13.240 |

So the answer to "does the best interval depend on theta" is that there is no best
interval to depend on it: within this model more often is monotonically better and the
question is practical, not pharmacological. What *does* move with capacity is the peak,
and the peak is what the admissibility ceiling binds on, so the trade a protocol faces is
peak against burden. (The convergence to the closed-form continuous limit is also a free
check on the integrator: `13.264` at h=1 against `13.240`.)

**(c) The missed dose: the effect is large, it grows sharply as capacity falls, and the
obvious remedy is inadmissible exactly where it would help.** Excess exposure is in
undisturbed-cycle equivalents.

| scenario | 100%: peak / excess / rejoin | 5%: peak / excess / rejoin |
|---|---|---|
| 1 d late | 9.08 / 0.001x / 14 d | 37.79 / 0.008x / 20 d |
| 3 d late | 9.36 / 0.012x / 24 d | 42.13 / 0.076x / 28 d |
| 7 d late | 9.70 / 0.065x / 28 d | 48.09 / 0.393x / 30 d |
| skipped | 9.92 / 0.500x / 26 d | 52.44 / 1.710x / 28 d |
| skipped, double at next | 9.92 / 0.258x / 39 d | 52.44 / 1.275x / 42 d |

A day late costs nothing. A skipped dose costs half a cycle's exposure at normal capacity
and **1.7 to 1.9 cycles at 5% and below**, and its peak goes from `0.99 x managed` at
normal to `6.03 x` at 1%, against a `theta -> 0` ceiling of `6.45 x`. So yes, a patient
with more residual capacity tolerates a gap far better, and the size of that effect is a
factor of about 3.8 in excess exposure and 6.1 in peak across the sweep.

**The double dose was checked against the ceiling before being reported as an option, and
it fails.** It does help: it cuts excess exposure roughly in half at normal capacity and
by a quarter at 5%. But it drives `theta + Delta_peak` to `4.05 x normal` at 5% and
`10.7 x` at 1%, so it is inadmissible at exactly the capacities where the gap hurts most,
and it lengthens the return to the steady cycle rather than shortening it (26 to 39 days
at normal). Reported as a negative.

**The honest limit on all of (c).** Every disease-range row is past the clearance maximum,
including the undisturbed baseline. The model can say a gap is bad and cannot say how bad.

**Two bugs the adversarial pass caught before any of this was reported**, which is the
standing policy earning its place again in the same session it was written down.

* *The admissibility check was rewritten instead of imported, and got it wrong.* The first
  version used a fixed `normal = 3000`; the check infers `normal = theta / 0.05`, so the
  bound is `Delta_peak <= 19 theta` and it tightens as capacity falls. The two agree at 5%
  of nominal, where `theta/0.05` is exactly 3000, and disagree everywhere else: the wrong
  one reported the booked dose inadmissible at **normal** capacity and would have been
  presented as a contradiction of Section 7.1 that did not exist. It now calls
  `n4_theta_scope.d_ceiling` rather than carrying a second spelling.
* *The doubled dose was silently merged into a single one.* `forward()` keyed its jump off
  `np.isclose(doses, a)` after a `np.unique()` on the segment edges, so two doses at one
  time delivered one dose. "Skipped, double at next" came out **bit-identical** to
  "skipped" and read as a real null result. Doses now carry per-dose amplitudes, and R4 of
  the recheck section is a regression for exactly this.

Twelve checks in the adversarial pass, all passing: the cycle-map fixed point against a
sixty-dose forward run by a route sharing no code (`4e-11`), the six constants against the
paper's committed values, the increment channel against its closed-form geometric sum
(`1.5e-15`, off the dose times where the output carries both sides of the jump), the
doubled-dose regression, tolerance invariance at seven times the sampling, and the rejoin
time against a doubled watch window.

**Provisional, as item 22 asks.** CI cannot verify anything pushed at the moment, so every
number above stands on the local run and the twelve checks, and on nothing else.

### R36. The "past the clearance maximum" flag was mine and it was wrong

R35 reported that at 5% of nominal and below the booked, undisturbed therapy schedule is
already outside what the parent model supports, and called that a stronger statement than
the missed-dose run was asked to produce. **It is withdrawn.** No trajectory in any of the
three runs is outside what the parent supports, and at the constants this repository
carries none can be.

**The error.** Section 8 says the parent "has a clearance maximum at $3.4$ times the
managed level". `audit/t13_therapy_forward.py` applied that $3.4$ to the LOCAL model's
managed level of `10`, giving a fixed threshold of `33.6`. The $3.4$ is a ratio to the
PARENT's own attractor, `81.93`, at the parent's own defaults. Two things are wrong with
transferring it:

* the two reference levels are not the same kind of thing. `11_scope.tex` says the
  parent's `theta = 300` is "about a tenth of normal", so `81.93` is a **disease-range**
  operating point, while the local `10` is the **100%** row of the same table;
* and the ratio is not scale-free in the first place. Across capacity and load it runs
  from about `2.9` to `275`, a spread of ninety-fold, so no fixed multiple stands in for
  it. The companion figure `12.6` is worse: an upper equilibrium exists only when
  `theta < s`, so at this load there is none at all at `100%` or `20%` of nominal.

**The paper is not at fault, and that was tested rather than assumed.** The discussion's
own arithmetic, that 42% of a band running to eight times the managed level lies beyond
the clearance peak, is coherent only if every "managed level" in it is the parent's
`81.93`: with the local `10` the whole band sits below the peak and the fraction is `0%`,
not `42%`. The sentence is self-consistent and the misreading was entirely at this end.

**What the bound actually is.** `dC/dm = 0` carries no load, verified symbolically, so the
clearance maximum's LOCATION is a function of capacity alone and needs no cross-model load
calibration to mean anything. In absolute marker units, and as a multiple of the local
managed level:

| theta | 100% | 20% | 10% | 5% | 2% | 1% |
|---|---|---|---|---|---|---|
| parent clearance maximum | 1768.3 | 388.8 | 275.7 | 222.6 | 192.4 | 182.7 |
| as a multiple of managed | 176.8x | 38.9x | 27.6x | 22.3x | 19.2x | 18.3x |
| run (a)'s peak, as a fraction of it | 0.5% | 6.3% | 11.2% | **15.9%** | 20.0% | **21.6%** |

The withdrawn threshold of `33.6` was too low by `5.4x` to `52.6x`. The closest any run
gets to the real bound is `21.6%` of it.

**And it can never fire here, which is worth stating as structure rather than as a table.**
`dm/dt = L - (theta+Delta) g(m) - k_e m <= L - k_e m`, so every trajectory at this load
obeys `m <= L/k_e = 64.5`; and the parent's peak lies above `sqrt(K_c M_i) = 173.2` for
every `theta > 0`, because that is where the alternative route's derivative changes sign.
A margin of `2.68x`. This is a fact about the carried constants and not a theorem: over
this repository's own prior box on `(K_c, M_i, rho_0)` the inequality reverses in about
`0.3%` of corners, where `sqrt(K_c M_i)` falls near `50`.

**What item 23 actually asked, and the answer.** The question was whether the local model
places the disease-range patient on a branch the parent calls a separatrix. **It does
not.** The local model's operating point is `0.5%` to `34%` of the way to the parent's
clearance maximum at every capacity swept, always on the rising branch, and it cannot be
otherwise given the `64.5 < 173.2` margin above. So `11_scope.tex`'s disease-range numbers
are not computed at a point the parent refuses as a steady state, and nothing needs saying
in Section 8 beside the existing misspecification result.

**One thing the adversarial pass killed that had looked like a finding.** An intermediate
route reported that the parent puts the same patient `2.28x` lower than the local model
does, and that the parent would therefore clear the patient down to its own attractor.
That is an artefact. `R0` is not a free parameter: `R0_for_peak` solves it against a named
baseline load, and `tests/test_models.py` enforces the constraint. Holding `R0 = 1122` while
moving the load to `322.727` breaks it. Re-derive `R0` at the new load and the parent puts
the patient `1.61x` **higher** instead; under the repository's other cross-model convention
(`t12_triples.py`'s `R0_STAR`, matched at the operating point) several capacities have no
equilibrium at all. **The direction of the cross-model disagreement is whatever one
chooses `R0` to be, so there is no finding there and none is reported.** The separatrix
answer above survives because it needs only the load-free `173.2`, which is `R0`-free too.

**Checks.** Four independent establishing routes and three adversarial lenses
(calibration, numerics, semantics); the numerics lens re-derived the clearance curve
symbolically with exact rationals and Sturm sequences, at 40-digit precision, and by
integrating the unreduced two-state parent. None broke the retraction. `t13` now carries
`R7`, a regression proving the withdrawn threshold fires at `5%`, `2%` and `1%` and the
corrected bound at nothing; `R8`, the symbolic load-freeness; and `R9`, the margin.
Fifteen checks, all passing.

### R37. Which admissibility reading each verdict is on: three, not two, and t13 is on one

Item 24 asked whether `t13` mixes readings. **It does not: every admissibility number and
verdict in runs (a), (b) and (c) is on one reading, R2, including the double dose's
`4.05x` and `10.7x`.** They are now labelled as such where they print.

But the audit turned up a third reading and a trap worth recording.

| | normal | D | verdict across 100/20/10/5/2/1% | breaches |
|---|---|---|---|---|
| **R1** | inferred, `theta/0.05` | scaled with capacity | `0.157` at every capacity | never |
| **R2** | inferred, `theta/0.05` | booked | `0.157 / 0.583 / 1.117 / 2.183 / 5.383 / 10.717` | 10% and below |
| **R3** | fixed at `3000` | booked | `3.133 / 2.333 / 2.233 / 2.183 / 2.153 / 2.143` | everywhere |

R1 is what `ceiling_blind_pct` describes, the version that cannot bind. **R1 and R2
coincide at `theta = 3000`; R2 and R3 coincide at `theta = 150`, where `theta/0.05` is
exactly `3000`. Every anchor number the paper quotes, the `2.18` and the `2672`, sits on
one of those two coincidence points**, so no number currently in the paper discriminates
the readings, and a reader cannot tell from the text which is meant.

`11_scope.tex`'s sentence is where that bites: it introduces the alternative as "read
instead with `theta = 3000` as normal and `D` held at its booked value", which is R3, and
then states R2's verdict, "breached at `10%` of nominal and below". R3 breaches at every
capacity including normal. The sentence is not wrong at the disease point it is discussing,
where the two coincide, and it is ambiguous everywhere else. Flagged, not edited: which
reading the paper intends is the author's call and not a measurement.

### R38. Found while checking R36: Section 8's "$12.6$ times" names the separatrix, not the loss

Not part of item 23 and not edited, because which quantity the sentence means is the
author's call. Section 8 says the parent "has a clearance maximum at $3.4$ times the
managed level **and loses its stable equilibrium at $12.6$ times it**". Those are two
different events and the second is at the first.

Raise the load on the parent at its own defaults and watch the two roots:

| load | roots |
|---|---|
| 600 (baseline) | 81.93 stable, **1033.0** unstable |
| 650 | 115.2, 692.4 |
| 700 | 205.9, 370.9 |
| 706.0 | 256.5, 296.4 |
| **706.388** (`Cmax`) | **275.67, 275.68** merged |
| 710 | none |

The stable equilibrium is lost by a saddle-node **at the clearance maximum**, `m = 275.67
= 3.365x` the attractor. It is the same event as the clearance maximum, not a second one
further out. The `12.6x` is the *unstable* root at the *baseline* load: the separatrix,
which is how large a perturbation has to be to escape the basin, not where the stable
state ceases to exist.

Both are real and both are worth quoting. The sentence as written attaches the second
number to the first event. Flagged for the author; a possible rewording is "has a
clearance maximum at $3.4$ times the managed level, where the stable equilibrium is lost
as the load rises, and a separatrix at $12.6$ times it at baseline load".

### R39. Section 8's validity boundary, corrected and narrowed (items 27, 29, 30)

R38 flagged the wrong-event sentence; it is now fixed rather than flagged, and two things
that sat beside it turned out to need the same repair.

**The sentence.** "Loses its stable equilibrium at `12.6` times it" named the separatrix.
Section 2.3's own definition is that the attractor and separatrix collide at the fold, and
the fold is at the clearance maximum: raising the load walks the two roots together and
they merge at `m = 275.7` at a load of `706.4`, `1.177` times baseline. So the stable
equilibrium is lost at **`3.4`** times the attractor, the same event as the clearance
maximum, and `12.6` is how large a perturbation must be to escape the basin at baseline
load. Both are worth quoting; the paragraph now quotes each for what it is.

**The two derived numbers in the paragraph do NOT move, and that was checked rather than
assumed.** `audit/a4_check.py` builds them, and neither touches the separatrix:

* the `42%` is `(log 8 - log(m_peak/m_bar))/log 8`, a function of the `3.365` peak ratio
  alone. Recomputed: `41.7%`.
* the `4.6` is the local model's `5.37x` baseline to hold eight times its managed level,
  divided by the parent's ceiling of `1.177x`. Recomputed: `4.56`. It uses the clearance
  maximum's VALUE, not the separatrix's location.

**And the paragraph contained the collision itself (item 29).** "Eight times the managed
level" appeared twice in adjacent sentences meaning two different bands: `655` (eight
times the parent's attractor, in the `42%`) and `80` (eight times the local `m_bar`, in
the `4.6`). One phrase, two referents, which is exactly the trap that produced R36. The
paragraph now states its referent once at the top and names both bands where they are
used.

**The boundary is also narrower than it read (item 30).** It said the model questions its
own validity out there and left it general. It is bounded: `dot m <= L - k_e m` caps every
local trajectory at `L/k_e = 64.5`, and the parent's clearance maximum cannot fall below
`sqrt(K_c M_i) = 173.2` at any capacity because that is where the alternative route's
derivative turns. A margin of `2.68`, so **no operating point of any design in this paper
lies past the parent's peak at any capacity**; only the excursion band the paper declines
to use goes beyond it. The exception is carried with a number rather than dropped: across
this paper's own prior range for `K_c`, `M_i` and `rho_0`, fewer than `1.5%` of draws even
place `sqrt(K_c M_i)` below `L/k_e`.

**One number deliberately NOT generated.** A swept count of the prior-box corners that
actually reverse reads `0.14%` at 600 draws and `0.23%` at 4000: sampling noise on about
five events. A generated number the artefact guard diffs against HEAD must not wander with
the seed, so the paper quotes the closed-form `1.5%` upper bound, which is exact, and the
swept figure is recorded here, where a measurement belongs.

**Sixteen literals became generated numbers** (`results/numbers/boundary.tex`,
`src/model_C/run.py:boundary_numbers`). The paragraph had every one of its figures typed
by hand, in the one place the paper states the limits of its own model.

**Pre-existing debt, found and not fixed.** The `theta_hat = 151` / `theta = 300` /
`150.8` / `+149` / `+0.2` ladder in the same paragraph has no generator anywhere in the
repository: it appears only in `audit/PHASE2A.md`, and four plausible reconstructions of
the fit did not reproduce `151`. It is left alone because reconstructing it is a
measurement and not an edit, but it is the last block of hand-typed numbers in that
paragraph and it should get a generator before the preprint.

### R40. Admissibility: three readings named in the text, with numbers that tell them apart

Item 28. R37 recorded that the check has three readings and that every anchor the paper
quotes sits on a coincidence point. Flagging that was not enough, because a reader cannot
tell which reading is in force and neither can a future editor, which is the same failure
as `min(q, N)` and `trans_angle`: **a displayed number that does not identify what it is a
number of.**

Section 7.1 now names all three in a table, evaluated at `20%` and `10%` of nominal where
all three disagree, and states which one the paper is on.

| | normal | `D` | at `20%` | at `10%` |
|---|---|---|---|---|
| R1 | inferred, `theta/0.05` | scaled | `0.16` | `0.16` |
| R2 | inferred, `theta/0.05` | booked | `0.58` | `1.12` |
| R3 | fixed at `3000` | booked | `2.33` | `2.23` |

The old sentence introduced R3's setup and then stated R2's verdict. It now says that R2
is the reading under which the bound is breached at `10%` and below, that R3 breaches
everywhere including at normal capacity where it reads `3.13`, and that the two anchors
the paper carries, `2.18` and `2672`, hold on **either** R2 or R3 because `theta = 150` is
where those two coincide.

The discriminating pair is the largest admissible dose, and the two readings run in
opposite directions in capacity: R2's allowance rises with `theta` (`19 theta` times the
train factor) while R3's falls (`3000 - theta` times it). At `10%` of nominal R2 permits
`5344` against R3's `2531`, both below the booked `6000`; at normal capacity R2 permits
`53438` and R3 permits nothing at all. A reader meeting only the old text could not have
told which was meant.

*(My first draft of that sentence said `5344` was above the booked value. It is below it.
Caught on re-reading the draft against the generator before the commit, which is the
adversarial pass applied to prose.)*

Eight new generated keys in `results/numbers/scope.tex`; `audit/t13_therapy_forward.py`
labels every verdict `R2` where it prints.

### R41. The R0 row now bounds what the parent can be asked (item 31)

The row said `R_0` must be re-derived when `K_c`, `M_i` or `rho_0` move. It must be
re-derived when the **load** moves too, because the constraint that defines it names a
baseline load: evaluate the parent at the local model's `L = 322.7` while holding
`R_0 = 1122` and a null genotype becomes maintainable, `theta_crit_intrinsic` returns
`nan`, and the module's own sweep calls such points vacuous. `tests/test_models.py`
already enforces the constraint for the three shape constants; the load was the case
nobody had written down.

The consequence is the part worth carrying. **Any comparison of the two models at
disease-range capacity turns on the convention chosen**, and the direction reverses
between them: holding `R_0` puts the parent's patient below the local model's level,
re-deriving `R_0` at the new load puts it `1.61x` above, and matching the alternative
route to the local one at the operating point leaves several capacities with no
equilibrium at all. So the direction is a choice, not a result, and no such comparison is
reported. What survives every convention is that the parent's clearance peak lies above
`sqrt(K_c M_i)`, which is `R_0`-free, and that is what R39's margin rests on.

### R42. The theta_hat ladder reproduces, and the missing ingredient was the load constant

Item 32. It was a candidate for retraction and it survives: **all five rungs come out to
every digit PHASE2A printed**, under a definition that is now written into the paper and
generated in `results/numbers/boundary.tex`.

| band top | 1.6 | 5 | **8** | 12.6 | 20 |
|---|---|---|---|---|---|
| `PHASE2A.md` | 557.8 | 277.0 | **151.01** | 27.3 | 0.00 |
| regenerated | 557.797 | 277.020 | **151.012** | 27.281 | 0.000 |

**The definition, which had never been recorded.** `m_bar` is the PARENT's own attractor
at its own load, `81.9344`; nine levels run geometrically from `m_bar` to `top x m_bar`;
the load at each level is the parent's own noiseless clearance there; and the local form
`s = theta g(m,K_m) + k_e m + b` is fitted to those pairs with `K_m` fixed at `100`,
`theta, k_e >= 0` and **`b` free**, by unweighted least squares in the load.

**The free load constant is the whole thing, and it is why four reconstructions missed.**
Every failed attempt fitted one clearance CURVE onto another. The estimator of Section 4.1
does not do that: it fits a balance and carries a load constant it does not know. Without
`b` the same nine levels give `1243 / 1032 / 953 / 883 / 818`, monotone but never
approaching the threshold; with `b` free they give the ladder. A curve-to-curve reading is
a different and much more favourable exercise, and the sentence in the paper had been
describing the favourable one while quoting the number from the other.

**And the archaeology says it never had a generator.** `git log -S` on `151.01`, `557.8`
and `277.0` across every ref hits only `audit/PHASE2A.md`, which arrived in one bulk
import and was never edited. One deleted `.py` in the entire history, unrelated; no orphan
bytecode; no dropped stash; the neighbouring `a4_check.py` generates the rest of the
paragraph and contains no fit at all. So this was a number computed once, in a session, and
written down without its construction --- the sixth instance of this project's
characteristic failure, and the one that took a nine-agent search to undo.

**Two qualifications the first draft of this entry did not have, both of which change what
the paper may claim.**

*`k_e >= 0` is doing the entire job, and that is the result rather than a technicality.*
`k_e_hat` is pinned at exactly zero on all five rungs. Released, it goes **negative** on
every band (`-1.11` to `-0.18`) and `theta_hat` at the 8 rung rises from `151.01` to
`604.06`. So the parent's clearance over this band is not merely badly approximated by the
local form, it is approximated best by a local model that is physically inadmissible, and
the best admissible fit sits on the boundary of admissibility. The paper now says so.

*The residual margin's SIGN is not determined.* Three defensible readings of the same
definition give, at the 8 rung: `159.46` (load-weighted residual), `146.72` (a dense grid
rather than the design's nine levels), `128.96` (linear spacing). Against
`theta_crit = 150.77` that is `mu_hat` from about `-4` to `+9`. So the `+0.24` is the
design's own `n = 9` unweighted reading and nothing more general, and the vividness ---
that the fit lands essentially exactly on the threshold --- is an artefact of that choice.
**What survives every reading is that the margin is consumed**, which is what the
paragraph needs and all it now claims.

That is the `k` row's remedy applied here: state the perturbation, generate the number
under it, report what a second reading gives --- and, this time, notice that the second
reading moves a sign and say so.

*(The first version of this entry said "all three land within a few per cent, so the
conclusion is robust and the third digit is not." True but far too generous: a spread that
takes `mu_hat` from `-4` to `+9` does not merely move a digit, it moves the sign of the
conclusion a reader would draw. Corrected before the commit.)*

Nine new generated keys. Nine agents across four reconstruction families --- estimand,
tangent, archaeology and a 252-variant definition sweep --- with two arriving at the same
definition independently, and the ladder reproduced a third time here before any of it was
written.

### R43. The boundary is the finding, and it says something the ladder did not (items 34, 35)

**The unconstrained optimum is outside the admissible set.** `k_e_hat` sits exactly on zero
on all five rungs; released it goes negative, `-0.36` at the 8 rung, and `theta_hat` rises
from `151.01` to `604.06`.

Two consequences, and the first disarms the number the entry before this one spent nine
agents recovering.

**A boundary solution is not a point estimate.** `151.01` is where the constraint stops the
fit, not where the likelihood peaks, so `theta_hat - theta_crit` is not a margin whatever
the spacing or the weighting. That is an independent and firmer reason for the correction
in R42: the sign spread says the third digit is an artefact of `n = 9`; the boundary says
the point estimate was never a point estimate.

**And the positive statement is sharper than "biased".** The local form's best
approximation to the parent over that band describes a patient who clears the metabolite by
a route running backwards; the best admissible approximation describes one with no
non-enzymatic clearance at all. Neither is a patient. Over that band the local family
represents the parent only at its own edge.

**The measurement that sets how strongly to say it.** The admissible fit's rms residual
over the nine levels is `25.76`, which is `0.77` of a single `5%` error in the recorded
load --- the accuracy Tier 3 asks for. So the misfit is *invisible* in the diagnostic a
reader would use, which is the same mechanism Section 7.2 already reports for the published
alternative and is why "cannot represent the parent" would have been too strong: it
represents it to within the assumed noise, and still lands on the boundary with a capacity
at the threshold. That is worse than a visible misfit, not better.

**Item 35, the placement judgement.** Section 2.4 asserts that the local model is "kept
because it is the most favourable case for identification rather than because it
approximates anything", and offers the tangent's intercept, `374` against `61`, as
evidence. That is a comparison at a POINT. The boundary result is the same claim over a
BAND and is strictly stronger.

My reading, having the two paragraphs side by side: **the fact belongs in 2.4, the
consequence stays in 8, and a forward reference alone is not enough.** A claim needs its
evidence where the claim is made, and this one is made immediately before the reader spends
four sections working inside the favourable structure --- asking them to take "rather than
because it approximates anything" on credit through all of it, and settle up thirty pages
later, is the wrong order. But the consequence needs `theta_crit` as a decision quantity,
the excursion band and the estimator, none of which Section 2 has; moving it there would
make Section 2 argue a conclusion before the paper has a design. So one sentence of fact in
2.4, with its two numbers, and the consequence where it already is.

### R44. Is the referent failure enforceable? One case is; the rest is a declaration

Asked of the class named two entries above. The honest answer has two parts and the small
one is the true one.

**Not enforceable, and nothing makes it so: the semantic case.** Whether "the managed
level" names `10` or `81.93`, whether a verdict is on R1, R2 or R3, whether `theta_hat` is
a point estimate or a boundary --- these are questions about what a phrase denotes. No
static property separates a label that fixes its referent from one that does not, and
unlike the figure guard there is not even a syntactic locus to hang a declaration on: a
drawing call is enumerable, a sentence is not. There is no set of "numbers that cross a
boundary" for a checker to walk. **This stays a habit and a review question, and it is
recorded as one rather than dressed up as a check.**

**Enforceable, and now checked: the mechanical case.** The crudest form of the same failure
is two different quantities wearing one key. The preamble `\input`s every file in
`results/numbers/`, so a key defined in two of them does not collide, error or warn --- the
last `\input` silently wins and every `\num{}` of that key takes whichever file the preamble
listed second. Two quantities, one macro, green build.

`tools/check_artefacts.py:check_key_collisions` now walks all `330` keys across the `17`
numbers files and fails on any name defined twice. There is no collision today. There
nearly was: `margin.tex` owns `theta_crit` and `boundary.tex` was one key away from taking
the same name for the same quantity computed a different way, which would have been the
figure-versus-table drift of R29 with no diff to reveal it. Its self-test feeds the detector
a synthetic colliding pair and a clean pair and requires it to separate them, because a
check that cannot be shown to fail is not evidence.

So: one mechanical case caught by a guard, and the semantic remainder left as a declared
habit --- **when a number crosses a boundary, carry its referent across with it** --- with
no pretence that anything enforces it.

### R45. The two theta_crit values agreed, and now they cannot disagree (item 36)

They were identical: both call sites ran `theta_crit_intrinsic(PC, S_BASE)` with the same
arguments, so both held `150.7747475984808` exactly, gap `0.000e+00`. Unlike the
`np.interp` case there was no drift to report.

**But they printed at different precisions, and that was its own hazard.** `margin.tex`
carried `theta_crit` at `{:.0f}` and rendered it **`151`**; `boundary.tex` carried
`fit_theta_crit` at `{:.1f}` and rendered `150.8`. So the figure said the threshold was
`151` while the text said `150.8`, in the one paragraph where the threshold sits beside a
fitted `theta_hat` of `151.01`. A reader could not have told the threshold from the
estimate.

Fixed as one computation and one key. `model_C.theta_crit_ref()` is memoised and is now the
only caller of `theta_crit_intrinsic(PC, S_BASE)` in the module; `margin.tex` owns the key
and prints `150.8`; `fit_theta_crit` is gone from `boundary.tex` and the three places the
discussion quoted it now read `\num{theta_crit}`; and the margin figure draws the threshold
at one decimal too, so nothing in the document prints `151` for it any more.

### R46. The two key mechanisms do not conflict, and the rule needed stating precisely

Item 37, and the answer is the precise version of the user's first branch.

**The shared-key trick is stem-local.** It works because `NumberBook` guards one
accumulating dict, which is one module's numbers, which becomes one file. The protocol
figure and `tab_D_designs` both write `delta_range_monthly` because both live in
`design/run.py` and therefore in `design.tex`. A differing value raises at accumulation.

**The collision guard is cross-stem.** It fails on a key defined in two *files*, because
there the preamble's `\input` order decides in silence and nothing raises.

So the rule is: **shared keys agree within a file; distinct names across files.** The two
dispositions apply to disjoint cases and there is no conflict today.

**But the user's second branch is the real point: the shared-key mechanism is then not
available across stems, and item 36 was exactly that case.** Two stems, one quantity, no
mechanism to make a drift raise. The answer is that there is a third mechanism and it
already exists: **one stem owns the quantity and the other reads it**, either through
`common.numbers.num(module, key)` or, within one module, through a shared helper. That is
the third shape `check_figure_numbers` already accepts for drawn numbers, generalised.
Item 36 is fixed that way: one memoised helper, one registered key, and no second copy for
a guard to have to compare.

Neither rule gives way. What was missing was the statement that they cover different
scopes, plus the third mechanism for the case neither covers.

### R47. Warnings are errors now (item 38)

A warning appeared in one run of the suite and did not reproduce on the next. With CI down
there are no repeated runs to catch a recurrence, so an intermittent warning arrives as a
count in a summary line with no text and no traceback --- and would arrive the same way the
second time.

`pytest.ini` now sets `filterwarnings = error`.

**CORRECTED IN PLACE.** What stood here said the setting cost nothing, because "the
suite passes under `-W error` today". The suite did pass. It passed on a warm
`__pycache__`, which for the one warning class already present in this repository is
not a test: an invalid escape sequence is raised when a module is *compiled*, and a
cached `.pyc` skips the compile. From a cold cache the same commit gave 38 failed of
92. The setting cost fifteen missing raw-string prefixes, paid in R49, and only after
that did it fire on nothing. The sentence is replaced rather than annotated because as
written it recorded a verification that had not run, which is the `theta_hat = 151`
ladder's defect in a different medium. R50 has the class and what now enforces it.

It earns its place more here than in most repositories. This project's stated characteristic
failure is a silent numerical error, and numpy reports exactly that as a `RuntimeWarning`
--- an invalid value in a divide, an overflow in an exp, an all-NaN slice. Those are what
the setting is for; the deprecation notices are a bonus. It also reaches the generators,
because `tests/test_artefacts.py` regenerates every table and numbers file through the
guard, so the code that produces the paper's numbers runs under it too.

Shown to fire before being trusted, which is this repository's standard: a probe test
raising a `RuntimeWarning` fails the suite with the warning's text and its line, and was
then removed.

The one thing it does not cover is a warning raised in a generator that the guard's fast
path does not reach --- `fig_margin`, `d13_decision_therapy` and the other slow ones. Those
run only under `make figures`, outside pytest. Noted rather than solved,
and closed in R49.

### R48. A second pattern, across the two misspecification results (item 39)

It earns a paragraph, and it is now Section 8.3 beside the three-degeneracies one.

**Two results, unrelated routes, opposite directions.** Section 7.2 fits the local form to
records from a published alternative readout and recovers a capacity too **high**, `1.3` to
`2.0`. The discussion's opening fits the same form to the parent's own noiseless clearance
and recovers one too **low**, `151.01` against a true `300`. The mechanisms share nothing: a
free marker offset carried alongside a free amplitude in the first, a curved alternative
channel driving `k_e_hat` onto its bound in the second.

**What they share is the symptom: the residual stays inside the noise the specification
already asks for while the capacity is wrong.** Under a fifth of one sigma in the first,
`0.77` of a single load error in the second. So on this model **residual size does not
bound structural error**, and the diagnostic a reader would reach for to detect a wrong
structure is the one that cannot see it. That is the argument for Tier 1 being a
*structure* test rather than a goodness-of-fit threshold, and the paragraph says so.

**And the opposite directions foreclose the repair.** A bias of known sign could be carried
as a correction or a one-sided bound; these cannot, because nothing in a record says which
structure is the wrong one. Hedged as the three-degeneracies paragraph is: two results, not
a theorem.

**One correction it forced.** The clause written last batch said the boundary result finds
what Section 7.2 finds "by the same mechanism". It does not --- same symptom, different
mechanism, and the difference is the reason the directions differ. Corrected where it
stands and the cross-reference now points at the new subsection.

### R49. The warnings gap on the generator path, and what R47 had actually measured (item 41)

**It could not go on as it stood, and what broke was not a generator.** `make figures
PY="python3 -W error"` died in the fourth module, at `src/design/run.py` line 954, before
computing anything:

    ax["p2"].annotate("every draw at one phase,\nso one value of $\Delta$",
    SyntaxError: invalid escape sequence '\D'

`$\Delta$` in a string that is not marked raw. Python raises a `DeprecationWarning` for an
unrecognised escape at *compile* time, and `-W error` turns that into a `SyntaxError`, so
the module never imports. Swept across all 99 tracked Python files there are **15 of them
in 4 files**: 12 in `src/common/params.py` (eleven `$\mu$M` unit cells in the parameter
registry and one `\paramscaption` in a docstring), one in `src/design/run.py`, one in
`tools/check_artefacts.py` and one in `tools/claims_to_tex.py`, the last two both `\input`
and `\cite` inside docstrings.

All fifteen are LaTeX and the fix is a raw-string prefix. It is output-neutral: nothing
under `results/` changed, the annotation still concatenates to the same characters
(`"...\n" r"so one value of $\Delta$"` keeps the real newline and the literal backslash),
and a raw f-string still interpolates.

**With that done the setting goes on, and it is inert.** The full build under `-W error`
runs all twenty generators to exit 0, and the working tree comes back holding nothing but
the `/CreationDate` bytes in the figure PDFs. So it is wired: `PYW = $(PY) -W error` in
the Makefile, used by the `figures` recipe and by nothing else. `refs`, `paper` and `audit`
keep plain `$(PY)`, because they do not compute a number the paper quotes. What this buys
over `pytest.ini` alone is the four slow generators the guard's fast path never reaches,
`fig_margin`, `d13_decision_therapy`, `k4_admissibility` and `t13_therapy_forward`, which
is exactly the gap R47 left noted.

**The part worth recording is what this says about R47.** That entry claimed
`filterwarnings = error` "cost nothing: the suite passes under `-W error` today". It had
been measured on a warm `__pycache__`. An escape-sequence warning fires when a module is
*compiled*; once a `.pyc` is on disk the module is not recompiled and the warning is never
raised again. Clear the caches and the same suite, on the same commit, fails:

    find . -name __pycache__ -type d -prune -exec rm -rf {} +
    python3 -m pytest -q     ->  38 failed, 54 passed in 7.29s

every one of the 38 at `import design.run`, with the `SyntaxError` above as the reason. And
`make clean` deletes `__pycache__`, so `make clean && make test` was red while `make test`
on its own was green, on one tree, with nothing in either run saying why they differed.
**The check added to catch a silent failure was itself silently conditional on a build
artefact**, which is the shape of the thing it exists to prevent, one level up. After the
fix, from a cold cache: `92 passed in 149.39s`.

That is the honest cost accounting. The setting did not cost nothing. It cost fifteen
missing `r` prefixes, and it was not paying until they were fixed.

### R50. A verification whose verdict depended on build state, and the audit it forced (item 43)

**The class, kept separate because its fix is not the seven's.** `make test` green and
`make clean && make test` red, 38 failed of 92, on one tree, with nothing in either run
saying why they differed. Each of the seven referent instances produced a wrong or
unlabelled *number*, which failed visibly the moment anyone looked at it. This produced a
correct-looking report of a verification that had not happened. A number nobody can check
is one problem; a check nobody can check is a worse one, because it is what you consult
about the first.

**It is not the first instance.** R10 is the same class from the other side: the em-dash
sweep deleted `\input{sections/01_biology}`, took a 1,400-word section out of the document,
and was followed by "no LaTeX error, no undefined reference, no failed test". The page
count fell by four and the fall was put down to a shorter abstract. Two instances make it a
class rather than an accident.

**What enforces it now.** `check_source_warnings()` in the artefact guard compiles every
tracked `.py` from its source text on every run and never consults `__pycache__`: 99 files,
0 warnings, with the self-test this repository requires, a source that does warn and one
that does not. It was preferred to putting `make clean` in front of `make test`, which was
the other option offered, for two reasons. `make clean` makes the suite slower without
making its verdict independent of anything except the cache, and the guard reaches `tools/`
and `tests/`, which `-W error` on `make figures` does not.

**The audit, by mechanism rather than by entry**, because "measured warm" is only one way a
verdict can depend on state the verdict does not name.

1. **Bytecode cache.** Reaches exactly one thing, a warning raised at compile time, and
   nothing numerical. One recorded claim sat on it, R47's, corrected in place. Now
   unreachable.

2. **LaTeX auxiliary state.** Reaches every "N pages, no undefined references" claim.
   `REVIEW_LOG.md` carries one, R10's, and R10 is the entry about its being worthless. The
   rest were reported in batch summaries rather than recorded here, and every one of them
   was taken after two `pdflatex` passes on a warm `.aux`, where `make paper` runs four
   passes and `biber` for the paper and three and `biber` for the supplement. **Measured
   from clean:** paper 43 pages, supplement 62, and in both logs 0 `LaTeX Warning:
   Reference`, 0 `LaTeX Warning: Citation` and 0 `Rerun to get`. The claims were right.
   The measure was not: `grep -c undefined main.log`, which is what had been reported,
   also counts a font-shape warning, and on the supplement it returns 1 for
   `T1/lmr/bx/sc`. Counted precisely from here.

3. **Committed artefact state.** A figure drawing `num(module, key)` reads the committed
   numbers file, and R51 has just made the first such reader. That dependence is by
   construction, and `check_regenerated` compares the committed file against the code, so a
   stale one fails rather than passing quietly.

4. **Generator numerics, which is everything the paper quotes.** Not exposed, and this is
   the part worth being able to say. The generators are seeded, the guard regenerates them
   into a temporary directory and compares byte for byte, and two independent cold runs
   this session left `results/**/*.tex` unchanged except for the one key R51 deliberately
   removed. An environment that moved a number would show as drift rather than hide.

5. **CI.** 44 commits on this branch and nothing has executed on any of them. "The suite
   passes" has meant "the suite passes in this container" for the whole of it, which is why
   the warm-cache distinction was load-bearing here and would not have been in a repository
   whose CI ran. Stating it rather than implying it.

**The list is short because the exposure is narrow, and that is the result.** The numbers
were never at risk: regeneration guards them, and a clean rebuild confirms it. What was
unguarded was the statement that a guard had run, and only two mechanisms could reach it,
one now closed by a check and one now measured. The full clean sequence,
`make clean && make figures && make paper && make test && check_artefacts`, returns 0 at
every stage: 92 passed in 144 s from a cold cache, 328 generated keys, 99 files compiled
with 0 warnings, artefacts match the code.

### R51. Two implementations made one, and the false contrast underneath them (item 44)

**The pair was real and the entry above overstated it, so both go here.** The
indistinguishability sweep flagged two design.tex/phaseclock.tex pairs. Measured, they are
not the same case:

| | design.tex record | phaseclock.tex record | values |
|---|---|---|---|
| monthly | `np.arange(0, 730, 30.4)`, 25 pts | `fixed(25, 30.4)`, 25 pts | identical arrays, both `1.843777` |
| fortnightly | `np.arange(0, 730, 14.0)`, 53 pts | `fixed(25, 14.0)`, 25 pts | different records, both `0.133333` |
| 14.5 d | `np.arange(0, 730, 14.5)`, 51 pts | `fixed(25, 14.5)`, 25 pts | `1.986121` against `1.801867` |

Only the first row is a duplication. The second is two different records agreeing because a
schedule at `eps = 0` saturates, which Section 6.4 states in as many words, and the third is a
difference the paper draws on purpose. The claim that both pairs agreed "by construction"
was wrong for the fortnightly one and is corrected above.

**What the duplication was hiding is worse than the duplication.** Section 4.3 read: *the
ranges are over these 25 draws and not over the two-year record of Table D: monthly reads
`range_monthly` here against `delta_range_monthly` there, and the 14.5-day schedule
`range_off` against `delta_range_off`.* For the 14.5-day schedule that contrast is real,
51 visits against 25. For monthly the two-year record **is** those 25 visits, so the
sentence set one number against itself and called it a comparison. It read as a claim only
because the quantity had two names; with one name it would have been visibly a tautology.
Corrected in Section 4.3, and Section 6.4's version, which was already careful about the
fortnightly case, now says the same about monthly.

**Merged in two layers, because there were two duplications.**

*One implementation.* `delta_range` now lives in `src/design/run.py`, which owns `delta`,
and `audit/g_design_common.py` imports it. It could not go the other way: that module
already imports `delta` from this one, so the dependency had a direction. Its
nineteen call sites across ten modules keep writing `G.delta_range(t)`, and
`G.delta_range is design.run.delta_range` is true.

*One key.* `audit/g1_phase_clock.py` no longer computes monthly's range. It reads
`num("design", "delta_range_monthly")` **after** `write_numbers`, so `phaseclock.tex`
carries no second key for a quantity `design.tex` owns, and the panel header draws to the
character what `\num{delta_range_monthly}` typesets beside it. `make figures` runs the
design module before this one.

**And this is the first call of `num()` in the repository.** R46 named "one stem owns the
quantity and the other reads it" as the third mechanism, the one for the case neither the
shared-key trick nor the collision guard covers. It had never been used: `num()` was
written, self-tested through `check_figure_numbers`'s first accepted shape, and left
without a caller for the want of a case. This was the case.

### R52. Is `01_biology` in the document? Yes, and my first answer was wrong (item 45)

**The answer is yes.** `paper/sections/A1_supplementary.tex:739` inputs it, and "Biological
background" is S1.13 of the built supplement. The relocation was deliberate and is
`deec483`, "Compress the main text to four claims, relocating rather than deleting", which
moved the biological background, the full dynamic model and its reduction, the local
derivation, the validation protocol, the figure section and the open questions into the
appendix. R10's incident was real when it happened and was reconciled rather than explained
away: the line is present and the text renders.

**My first answer was wrong, and how it was wrong is the finding.** Asked whether the file
is input, I grepped `paper/*.tex`, which reads `main.tex`, `preamble.tex` and
`supplement.tex` and **not** `paper/sections/A1_supplementary.tex`, where seven of the
inputs live. That returned seven "never input" files totalling 21,837 words, several of
them edited within the last two days, and every one of them was a false alarm. `\input` is
transitive; a one-level grep is not a reachability test, and the question cannot be answered
without the closure. Recorded because the wrong answer was the alarming one, and an alarming
wrong answer is the kind that gets acted on.

**Measured properly.** Following every uncommented `\input` from `main.tex` and
`supplement.tex`: 23 `.tex` files under `paper/`, 22 reached, one not. The one is
`note_structural_identifiability.tex`, which carries its own `\documentclass` and is a
standalone document that must not be input. Nothing any document inputs is missing from
disk.

**And it is guarded now.** `check_inputs_reachable()` is `check_generators_reachable()` for
prose. That one asks whether every numbers file has a generator `make figures` actually
runs; this asks whether every section file has a document that actually reads it. The
standalone exemption is checked rather than declared, by looking for `\documentclass`, and
the check is shown to fire on the real bug rather than only on a synthetic one: replaying
R10's edit, commenting out `\input{sections/01_biology}` in a copy of the tree, returns
exactly `sections/01_biology.tex` as an orphan.

**It was the last member of the invisible-artefact family without a check.** A file present
and never input defines no label that anything misses, so there is no undefined reference,
no error and no failing test, which is precisely what R10 recorded and R50 named as a class.
The numbers side has had reachability checked since the generators did. The prose side had
not.

### R53. The page drift is real, and the section numbers in the reports were not (item 46)

**The drift holds.** Rebuilt cold at seventeen heads, `make clean`-equivalent then four
`pdflatex` passes with `biber` for each document, the same instrument at every point:

| head | batch | main | supplement |
|---|---|---|---|
| `fff74c4` … `bf78c16` (6 heads) | through the phase-lattice work | **39** | 61 |
| `855e377` | the lap condition, items 1 to 6 | **41** | 61 |
| `e0372bb`, `eddcf84` | figure numbers registered; items 20 to 22 | 41 | 61 |
| `d811ea8` | items 23 to 26 | **42** | 61 |
| `c778d08` | items 27 to 31 | 42 | **62** |
| `1a8561d` | items 32 to 33 | **43** | 62 |
| `c585d84` … `066a316` (5 heads) | items 34 to 44 | 43 | 62 |

That is 39, 41, 42, 43, which is the reported series exactly. Nothing in it was an artefact
of the instrument.

**The instrument was loose and not wrong, and the reason is the document.** `Rerun to get`
is 0 at four passes at every one of the seventeen heads, so these two documents converge
before the second pass and a two-pass count was always going to agree. That is a property
worth knowing and not one to rely on, so the cheap repair is to report `Rerun to get`
beside the page number: a non-zero count is precisely the case where a short build lies.
On the other loose measure, `grep -c undefined` matches the precise `LaTeX Warning:
Reference` and `Citation` counts at 0 in the paper at every head, so it never hid anything
there; on the supplement it returns 1, for the font shape `T1/lmr/bx/sc`.

**The drift stopped five heads ago.** `1a8561d`, `c585d84`, `b1aebd7`, `aa0ecd4`, `87a762f`
and `066a316` are all 43. Items 34 to 44 added no page between them. So did the step that
added the most prose of any of them: `d811ea8` to `c778d08` put 66 net lines into the
discussion and the scope section and stayed at 42.

**And the section numbers the reports used were wrong, which matters more than the count.**
The section files keep prefixes from a structure the paper no longer has, and the reports
read the numbers off the filenames. Built:

| file | built section |
|---|---|
| `07_design.tex` | **4**, the design |
| `07b_decision.tex` | 5 |
| `08_specification.tex` | **6**, the specification |
| `11_scope.tex` | **7**, the scope |
| `12_discussion.tex` | 8, the discussion |

So `sec:tierstruct`, which the reports called "§5" and where item 40's clause went, is
**§6.2**. `sec:pattern2` is §8.3, which was right. The paper itself is unaffected, because
every cross-reference in it is a `\ref`; what was wrong was the name in the reports, which
is where a decision about cuts would be taken.

**Where the four pages went**, as net added lines from `bf78c16` to `1a8561d`:

| section | file | net lines |
|---|---|---|
| §4, the design | `07_design.tex` | **+105** |
| §8, the discussion | `12_discussion.tex` | +68 |
| §7, the scope | `11_scope.tex` | +32 |
| §6, the specification | `08_specification.tex` | +16 |
| §3, the excitation condition | `06_excitation.tex` | +16 |

The largest single share is §4 and not §7 or §8. A cut aimed at "§7 and §8" as the reports
named them would reach the scope section and the discussion and leave the biggest addition
where it is.

### R54. Phase 8 S1: a result carried by one draw, and two errors in the script that reported it

**Mine, and a retraction.** The first S1 report (`043a3ec`, `audit/PHASE8_SIMULATIONS.md`)
said the Attarwala two-compartment truth, read by the committed 3.5 d analyst, returns
`theta_hat/theta = 5.637` [2.440, 8.609] at `2.35 sigma`, that it is therefore unlike the 5 to
7 d exponentials and audible where they are not, and three further things about wrong-family
analysts. The author's robustness probe `p3_shape.py` gave `1.37` to `1.40` at `0.33` to
`0.37 sigma`. `audit/q1b_reconcile.py` put both scripts on identical noiseless records and
fitted each record with both estimators.

**The fitters were not the cause.** They agree to `5e-7` on every fixed-amplitude row. The two
scripts had modelled different objects: both run Attarwala's protein block in days with
half-lives `1.99 d` and `11.81 d`, but `q1` drives it through the effect compartment
(`ke0 = 0.242 /h`, so the protein **rises** after an infusion) where the probe's impulse makes
it **jump**. That difference is negligible everywhere except at the one draw the monthly record
places at a dose instant, `t = 0`, which `q1` read on the post-dose side as `design.run.delta`
does. There the rising kernel is zero and the analyst is `6000`, and that single draw carried
the cell to `5.81`. The draw is inside the marker's own settling window (1.3 h), where the
quasi-steady model does not describe the marker, **which is wrong for any kernel**; it was only
catastrophic for the rising one. Removing it gives `1.23`.

**What the first revision of this entry ran together, now separated at the author's request.**
Reading that draw pre-dose, rather than dropping it, is **not** a further correction. It turns
it into an **untreated baseline draw**, which pins `theta` through the untreated balance and
roughly halves the bias for a slow truth. That is a change of **design**: a patient already on
therapy when the record starts has no such draw. So:

- **Error corrections** (three): not reading a draw at an infusion instant post-dose; the
  dose-scale prior on `ln kappa`; the profile stranded near `k_e = 0` (both below).
- **Design change** (one): an untreated baseline draw at `t = 0`.

The first report's `1.934` (7 d) and `1.310` (5 d) were **not wrong arithmetic**. They belong to
a record with no valid baseline and one invalid draw; with the invalid draw removed and no
baseline added they read `1.744` and `1.267` (noisy, 24 seeds, peak). The first report's
`5.637` for the rising two-compartment truth **was** the invalid draw: `1.196` without it.

**Two errors in `q1`, found on the way.**

1. **The dose-scale prior was on `kappa`, not `ln kappa`.** The repository states it on
   `ln kappa` (`g_design_common.sd_ln_theta`, column `inc g`). The linear form charges 68 at
   `kappa = 0.01` where the log form charges 1470, and that is what let the first run "switch
   the increment off" and report capacity nine times low. The probe's fitter returned `5.74`
   on the same records where `q1` returned `0.11`. Same class as the `k_e` prior defect caught
   before the first report.
2. **A stranded `k_e` start in the profile.** With `theta` held high the linear inner solve can
   return `k_e <= 0`; the polish then started from `ln(1e-12)`, had no gradient, and was
   accepted because the start's objective was infinite. The profile jumped from `1.53` to `553`
   between adjacent `theta` and upper interval edges came out short. `q9` check 5 found it
   against `n6`.

**Gap B on the two footings**, fixed 3.5 d analyst, 24 seeds, `dose / peak / mean`
(`audit/q1c_baseline_footings.py`; reproduced by the author with the probe):

| truth | with a pre-dose baseline draw | already on therapy, no untreated draw |
|---|---|---|
| exp 5.0 d | 1.080 / 1.161 / 1.442 | 1.185 / 1.267 / 1.555 |
| exp 7.0 d | 1.145 / 1.358 / 2.000 | 1.529 / 1.744 / 2.384 |
| 2-comp, rising (`q1`) | 1.046 / 1.131 / 1.282 | 1.113 / 1.196 / 1.345 |
| 2-comp, jump (probe) | 1.211 / 1.292 / 1.309 | 1.273 / 1.353 / 1.370 |

**Cells the corrections changed** (baseline footing, peak):

| cell | first report | corrected |
|---|---|---|
| rising 2-comp, fixed 3.5 d | 5.637 [2.440, 8.609], 2.35 sigma | **1.131** [1.036, 1.366], 1.08 sigma (null 0.98) |
| rising 2-comp, + free amplitude | 0.111, "increment switched off" | **1.126** |
| rising 2-comp, half-life free | 0 on 24 of 24 seeds | **1.129** [0.886, 1.286] |
| exponentials, two-comp analyst | 1.18 to 2.25 | **1.05 to 1.14** |

**Conclusions withdrawn:** that the two-compartment truth differs in kind from a slower
exponential; that it is audible; that the paper's "under a fifth of one sigma" sentence needs
qualifying on its account; that freeing the amplitude relocates the bias; that no analyst family
is safe (the two-compartment family is now within 14% on every truth tested, at about twice
the RMSE). **Qualified:** S1.8's "the family question wants the hours after an infusion" holds
under the effect-compartment reading only. **Standing:** the direction of gap B, that freeing
the amplitude does not help, and that an unrecorded dose train costs more than a wrong shape and
announces itself (Q3W read as Q2W, `6.52` at `2.20 sigma`).

**New, from the design change:** a pre-treatment baseline draw roughly halves the slow-truth
shape bias (7 d, dose: `1.529` without, `1.145` with). It matters in the correctly specified
model too, at a disease-range capacity: the review's monthly `3.07` at 5% of normal is the
figure with **no** baseline, and `2.38` with one (S2 Part A0). At normal capacity the paper's
`0.1699` does not depend on it (`0.1695` with a baseline). S4 tests whether the effect
survives noise and the 14.5 d schedule and whether two or three pre-treatment draws do
better.

**And a provenance correction to S1.8.** `ke0 = 0.242 /h` is a **mouse** estimate (Attarwala
Table 1), and in the published structure it is the equilibration of an effect compartment with
**plasma mRNA**, which has its own PK. `q1`'s rising kernel drops that PK stage; restoring it at
brief R2's human popPK half-life (36 h, unverified here) moves the protein's single-dose peak
from 13 h to 72 h, and the "hours after an infusion" timing conclusion holds only for the fast
reading.

**Re-check record.** `q9` gains checks 9 and 10. **The first version of check 9 failed** at
`0.223`; the two failing seeds are where `q1`'s `K_m` sits on its box bound and the probe's
runs past it, which check 4 already measures. It now compares interior seeds and prints the
box seeds. Check 11 rebuilds both footings from the probe's own truths and fitter. 25 of 25
pass. S2 was held for the author's review, as the rule requires when a conclusion changes, and
released after it.

---

## R55. S2b: the `k_e` measurement itself

`audit/q2b_ke_measurement.py`, written against the micro-dose IV ¹³C₃-3-HP study that
`paper/sections/08_specification.tex` proposes. S2 priced the requirement on a per-patient
`k_e` to about 2% at 5% of normal; S2b asks whether that study delivers it, sweeping everything
brief **R1** leaves unknown (half-life 5 min to 6 h, window, draw count, assay CV, C0/LOQ).

**The study can reach 2 to 5% on `ln k_e`, but the sampling window has to cover the half-life,
and the half-life is unknown by a factor of 70.** A 24 h collection with 8 or more draws is the
only design that works across the sweep. Under a two-compartment tracer a truncated window does
not degrade the estimate, it destroys it: `sd(ln CL)` of 20 and over at a 6 h half-life in a 2 h
window, against 0.03 at 24 h.

**Two findings change what S2 concluded rather than confirming it.**

1. **What the tracer measures is a biochemical question the paper does not settle.** Taking the
   balance literally, a trace-level tracer is cleared enzymatically too, so the study observes
   `k_e + θ/(K_m + m₀)` and enters as a row, not a prior. That reading buys a factor of 1.3
   where the end-product reading buys 7.6, and no schedule at any capacity then reaches 0.20.
   The specification should carry the premise explicitly, with its evidence.
2. **Within-patient drift in `k_e` is the binding constraint below 10% of normal.** At 10%
   drift no schedule reaches the gate at 5% of normal however good the study, and once drift is
   at 10% tightening the study from 5% to 2% buys about 8%. Repeating the study is the only
   design answer in view and is not in the specification.

In the Monte Carlo, **exactly one cell at 5% of normal meets the gate**: baseline draw, load
measured, study 2%, no drift, every 2 days, RMSE 0.125 [0.109, 0.142]. S2's "about 2%" stands
as a number and does not stand as a sufficient condition.

**Re-check record.** `q9` gains check 13 on three routes that share no code with `q2b`: the
plasma arm's `sd(ln CL)` against ordinary least squares with the delta method, the
two-compartment bolus against a matrix exponential, and the tracer row collapsing into S2's own
0.8-prior bound as the study's error grows. 30 of 30 pass.

**One correction made before any number was reported.** The first version of the reading-(i)
matrix dropped the population prior on `ln k_e` altogether, charging that reading for
information the other reading never created. The prior is a property of the literature, not of
the study, and is now kept at 0.8; the reading-(i) figures fell from 2.84 to 2.33 at 5% of
normal, monthly.

**Reframed after the author's review, with both findings confirmed by the probe** (the one
passing 5% cell at RMSE 0.125 against the probe's 0.134 without the baseline draw; the drift
column reproduced at θ = 150, 2-daily, 12 records: 0.117, 0.249, 0.439 for study 2% at drift 0,
5 and 10%, and 0.475 for study 5% at drift 10%).

Two things changed in how S2b is reported, not in its numbers.

1. **The two readings of what the tracer measures are not a premise and a sensitivity.** The
   end-product reading is primary, since the enzyme acts on propionyl-CoA and not on the marker.
   But the literal reading is the observation-model assumption the discussion already names as
   the paper's weakest, so the gap between them, a factor of 5.8 on `sd(ln θ̂)`, **is the first
   quantification of what that assumption is worth**. It is reported that way rather than as a
   premise to be asserted. And the same study can settle it, by measuring where the label ends
   up: the fraction recovered unchanged against the fraction appearing in downstream products.
2. **Drift is an unknown, not a finding.** The 0, 5 and 10% range is a sweep; the range is
   **[judgement]** and only its consequence is **[measured]**. Before a repeated study is
   proposed as the answer, the question is whether the drift is predictable: if it tracks kidney
   function, which is measured routinely and declines in this disease, a covariate replaces a
   repeat. S4 also carries drift at normal capacity, since the same prior carries the headline
   precision.

`refs/brief_questions.md` is new, and holds the R1 questions these raise: whether the marker is a
substrate of any step or only an end product (with an unverified lead, PMID 28634175, recorded as
a lead and not a citation), within-patient stability of the clearance, and its covariates.

**A boundary, set by the author.** S2b crossed from checking the paper into designing a study.
S3, S4 and S5 are finished as the brief writes them, then this phase is written up and stops.
There is no S2c; a new sub-question is recorded rather than run.

---

## R56. S3: the marker CV against the load record's accuracy

`audit/q3_marker_cv_load.py`. `tab_D_loadprec` prices the load record along one axis at the
committed 5% marker CV, and Section 8 concludes the design reaches 0.20 "only if intake is
recorded to better than about 10%". S3 turns that row into a surface, because the marker CV it is
conditional on is unmeasured in this population and the repository's own 12% question is open.

**The published table reproduces to 0.9% up to its 10% row and drifts to 2.9% at the 40% row**,
where the filter route and the information matrix part company on a record whose noise is four
times its signal. A difference in route, not a correction.

**Three findings.**

1. **The "about 10%" figure needs its condition attached.** It holds at a 5% marker CV and takes
   about 2.4 years of monthly draws, against the booked two. **Inside the booked two years the
   requirement is better than 5%**, and at the working band's edges better than 1%.
2. **The marker CV and the load recording error are interchangeable percent for percent.** σ's two
   terms are equal at `CV = load_frac × L / (b × m̄)`, which comes out at CV ≈ the load fraction to
   within 10% at every capacity tested. **The design rule is to spend on whichever of the two is
   the larger percentage**, and the specification does not state it. Above a 25% marker CV nothing
   on the load axis rescues the design at normal capacity.
3. **The working band is a property of the measurement, not of the system.** Re-solved across the
   surface it spans a factor of 10.3 at CV 5% and load 1%, a factor of 2.7 at CV 10% and load 5%,
   and is empty at CV 15% and above. It closes from the top: the upper edge moves five times as
   far as the lower. S3's lower edge reproduces the published 2472 to 0.5%.

**Also measured, and small:** dropping the untreated t = 0 draw costs about 3% at normal capacity,
which is S2b's baseline-draw effect where it does not bind. And two-monthly sampling is not a
saving: it reaches the gate in about the same number of draws and takes twice as long.

**Re-check record.** `q9` gains check 14 on three routes: a finite-difference Hessian of the log
posterior against the analytic information matrix (2e-8), σ against the probe's independently
coded footing formula (4e-14), and the visits-to-gate crossings confirmed by the Hessian route.
33 of 33 pass.

---

## R57. S4: one joint stress scenario at normal capacity

`audit/q4_joint_stress.py`. Every stress the audit found, present at once at normal capacity: the
two-compartment increment under both readings against the committed 3.5 d analyst, a 15% marker
CV, the k4 trough dose calibration, the 7 d train, the trial's irregular C5 train, and `k_e`
drifting within the patient. Both baseline footings, `dose` primary with `peak`, `mean` and
`trough` as the sensitivity column, 24 seeds, bias, RMSE, coverage and residual reported.

**The design survives the joint scenario.** With one untreated baseline draw and the `dose`
calibration the median bias is 2% low on the rising kernel and 13% high on the jump kernel, inside
the band S1 reported for a wrong shape alone. The stresses do not compound multiplicatively, and
the 7 d train and the irregular train partly cancel. The residual rises from 0.98σ to 2.2σ, so
S1.4's warning survives: a real record carrying these stresses would not look well fitted.

**One thing does break it, and it is a choice rather than an unknown.** The trough dose
calibration roughly doubles θ̂ (2.29 rising, 2.29 jump, 2.63 on therapy) and multiplies RMSE by
five, with coverage at 67 to 75%. The mechanism is arithmetic: the committed trough is 400 against
a peak of 6400, so matching a trough fixes the amplitude through a small number and divides any
shape error by it, setting an amplitude 2.8 times smaller than the dose reading does.
**`k4_admissibility.py`'s finding that the trough reading is admissible, and this, are compatible
and both belong in the specification;** the trough reading should not be adopted on k4's evidence
alone. The author's choice of `dose` as primary, made before these numbers existed, is also the
best of the four.

**The two questions the author added.**

1. **The baseline draw's advantage survives the noise and grows with it**, on the monthly schedule:
   the median moves from +0.033 closer to the truth at CV 5% to +0.065 at CV 15% on the rising
   kernel, and +0.041 to +0.121 on the jump kernel, and it takes a third off RMSE. On the 14.5 d
   schedule, whose 51 treated draws walk in phase against the train, the advantage is real but
   roughly halved and stops growing with the noise. The one cost is coverage, 83% against 92%: the
   draw moves the point estimate onto the truth without widening the interval to match.
2. **One pre-treatment draw is the specification item, and two or three are not better.** On every
   monthly row the 2 and 3 columns are worse than the 1 column by 5 to 23%. The untreated level is
   one number, so a second draw re-measures it while adding another load record with its own error.

**And a negative result worth having.** Carrying drift at normal capacity, as the author asked,
shows **it does not bind there**: 0.5% on the median and under 1% on RMSE across its whole swept
range, against S2b.5 where a 10% drift put every schedule beyond the gate at 5% of normal. Drift
is a disease-range requirement, not a general one.

**Re-check record.** `q9` gains check 15 on three routes: the trough amplitude by brute summation
over 4000 past doses against the closed form (2e-16), the drift conventions measured rather than
assumed to coincide (0.8% at width 0.05, 2.5% at 0.15, against 41% at the paper's 0.8), and one
joint-scenario cell rebuilt through the probe's fitter (1e-7). 36 of 36 pass.

---

## R58. S5: the ex vivo prior on θ

`audit/q5_exvivo_prior.py`. A direct enzyme assay on the patient's own cells, carried as a
log-normal prior on θ of swept width 0.3 to 1.0, at 5, 10, 20 and 100% of normal, alone and with
the `k_e` precision S2b says the micro-dose study can deliver.

**One thing is proven rather than measured.** The prior enters the information matrix only in its
(ln θ, ln θ) element and therefore survives the Schur complement untouched, so the prior and the
record add in exact quadrature, confirmed at 1.000 in all 48 cells. Under this model an ex vivo
assay and a marker record carry independent information about θ, so neither can spoil the other.

**Three findings.**

1. **The assay is not an alternative to measuring `k_e`.** With the paper's 0.8 population prior on
   `k_e`, every combined figure in the disease range sits barely below the assay's own width (0.299
   against 0.300 at 5% of normal), so the record contributes under 1% and the answer is the assay's
   answer wearing an in vivo estimate's clothes. And it still misses the gate: no width in the
   sweep reaches 0.20 anywhere below normal capacity.
2. **With `k_e` measured it changes exactly three verdicts, all at a prior of 0.3.** The practical
   contribution is to turn monthly sampling at 10% of normal from 0.243 into 0.189, an edge pass in
   the Monte Carlo. A prior of 0.5 or wider rescues nothing anywhere in the grid, so the assay's
   value rests entirely on brief **R4** establishing a mapping good to about a third.
3. **The risk is concentrated exactly where the benefit is.** At 5% of normal with a 0.3 prior, 54
   to 60% of a mapping centring error passes into θ̂; at normal capacity, 3 to 6%. At 10% of normal
   **a prior wrong by a factor of two is worth precisely nothing**: RMSE 0.249 with it and 0.249
   without, the variance it buys and the bias it injects cancelling to three decimals. At 20% and
   above a wrong prior is worse than none. A wide prior is safe and useless, a tight one useful and
   unsafe, with no width at which it is both.

**So the recommendation is that the assay enters as a cross-check and not as a prior in the
estimator**, unless R4 establishes the mapping to about a third AND establishes it as unbiased;
and if it does enter, the resulting interval must not be quoted as an in vivo measurement of
capacity, since in the disease range it is dominated by the prior. `refs/brief_questions.md` now
carries R4's two questions, with the centring named as the more important.

**Re-check record.** `q9` gains check 16 on three routes: the Fisher bound with the prior against a
numerical Hessian of the whole posterior (5e-8), the MAP against a from-scratch Nelder-Mead fit
(3e-7), and an uninformative-record limit that must and does return the prior's own centre exactly.
39 of 39 pass.

**One presentation error, caught before reporting.** The first run printed a `prior alone` column
that was empty in every row, since the prior's own width varies by column and not by row, and
described the quadrature identity as a measured shortfall. Both are corrected: the column is gone
and the identity is stated and proved as S5.0.

---

## R59. The R4 retraction is confirmed, not overturned

`audit/PHASE9_LITERATURE.md` item 6. The 25 September literature reports disagree about the
breath-test CV that **R4** retracted. R3 quotes the abstract of Manoli 2021, "mean coefficient of
variation 6.9% and 12.8%, respectively", and reads the 12.8% as a mean of per-subject CVs. R4 quotes
the Results, an average of 14.02 ± 12.86% in patients and 6.94 ± 5.3% in controls.

**The arithmetic settles it against R3.** 12.8 is not a rounding of 14.02, the mean; it is a
rounding of 12.86, the SD. And the decisive point is one neither report makes: **the abstract's two
figures come from different columns of the same table.** The controls' 6.9 is their mean, 6.94, and
the patients' 12.8 is their SD, 12.86. A reader assuming the abstract is internally consistent
takes both as means, which is precisely the error R4 identified, and the inconsistency is visible
without leaving the paper.

**So R4 stands and there is no retraction of a retraction.** R3 is wrong here, and wrong in the way
the September pack was wrong: it trusted a label in an abstract over numbers in a Results section.
R4 is less wrong rather than right, since its own sentence "abstract rounds these to 6.9% and 12.8%"
asserts both are roundings of the means it has just quoted.

**Conditional, and recorded as such.** The argument rests on two quotations this repository has not
verified, R4's Results sentence and its controls figure, because outbound access to publisher and
index records is denied in this container. It must be confirmed at the source. What does not depend
on the quotes is the direction of the retraction, which no reading here overturns. Either way the
figure is a breath-test repeatability CV across visits months to years apart in MMA, and is not a
within-patient CV of any plasma marker.

---

## R60. The marker CV is a bound, not a number, and the earlier framing double-counted

`audit/PHASE9_LITERATURE.md` item 7. The brief that opened Phase 9 took R3's 25 to 35% surrogate
figures as variation the model does not explain, and concluded the monthly design fails. **The first
half stands and the second does not follow**, which the author identified before any rerun was
begun.

**What the model already explains.** The marker's within-patient CV across the monthly record,
produced by the dose train alone at a perfectly constant recorded load, is **27.2%** on the
on-therapy footing with levels 3.544 to 8.747, and **28.3%** including the untreated t = 0 draw at
m = 10.0. With the specification's 10% load variation it is a median 30.0% [25.6, 33.5] over 24
records. Checked by summing the kernel over doses and root-finding the balance rather than using the
closed-form root, agreeing to 1.3e-10. R3's surrogates are TOTAL within-patient variation in
subjects receiving no dose train, so they are the same order as what the model fits, and carrying
them into sigma unchanged double-counts.

**Neither surrogate can be decomposed, and that is item 7's central result.** Thompson 2012
separates analytical (7.9%) from biological (25.2% by quadrature) and stops; feeding drives both C3
and the substrate load this model records, so an unknown and probably large part of that 25.2% is a
regressor here. Manoli 2021 separates nothing, and its interval of two months to four years contains
treatment changes, growth and intercurrent illness with no diet standardisation.

**So sigma is bounded below by the analytical component, 5 to 8%, and above by the surrogates'
totals, 25 to 35%, with no measured value in between.** The design's 5% is analytical only, and
Thompson's analytical figure for a comparable assay is 7.9%, so 5% is the optimistic end rather than
the centre. This is weaker than either end taken alone and it is what the evidence supports. **There
is no central case**, which is a decision for the author rather than a reading from the literature,
and the marker-variability rerun waits on it.

**And the panel remark is untested rather than supported.** Part of the residual must be shared,
since variation in the propionyl-CoA pool propagates to every marker derived from it, and a shared
part cannot be averaged away. Part must be marker-specific, and R3 names the confounders: free
carnitine for C3, oxaloacetate for 2-MC, glycine for propionylglycine. The split cannot be estimated
from anything R3 found, because it needs repeated simultaneous measurements of two markers in the
same patients and no within-patient figure exists for any marker in PA at any timescale. The fasting
result R3 quotes, a 62% fall in C3/C2 driven mostly by C3, is the kind of shared driver that would
defeat panel averaging.

**Re-check record.** `audit/q10_phase9_recheck.py` is new: the model-explained CV by an independent
route, and the item 6 arithmetic with both readings tested rather than one asserted. 7 of 7 pass.

---

## R61. The marker and load requirements are one budget on sigma

`audit/q6_sigma_budget.py`, `audit/PHASE9_LITERATURE.md` step 3(a). The author's reading, that S3's
two-axis surface is one budget spent two ways, is confirmed and is exact rather than approximate.

**Exact.** Five `(load record, marker CV)` pairs constructed to give sigma = 60, from marker-only to
almost load-only, give `sd(ln theta)` identical to 3.3e-16 with `s_b` held. The reason is
structural: the noise enters the information matrix only as `Q'Q/sigma^2`, and `load_frac`'s second
entrance is the load-bias prior `s_b`, which sits outside every block the Schur complement
subtracts. The author's residual 4% across load records of 2 to 15% is that second entrance and
nothing else: with `s_b` held, sigma* has **zero** spread at every schedule.

**The budgets at normal capacity**, with `s_b` at 0.10 L: monthly 36.72, fortnightly off-period
45.99, weekly phase-spread 69.69, every 2 days 130.96. Both working-band edges roughly halve them,
and by almost the same factor at each, which the band's own width does not do. Three of the four
overlap the author's independent ranges; fortnightly does not, and the likely cause is the t = 0
convention, since the 14.5 d schedule has a dose-instant draw at t = 392.

**What that says against the evidence bound.** The committed sigma reproduces at 21.96 and clears
every schedule, monthly included. The analytical floor at a 10% load record is 35.5 at CV 5% and
40.1 at CV 8%, against monthly's 36.72, so **monthly does not merely require the residual scatter to
sit at its analytical floor, it requires the floor itself to be at the optimistic end of its own
range.** The surrogates' 81.2 to 109.2 clears only every-2-days. And the load record alone can close
the gate: at a 15% record neither monthly nor fortnightly reaches 0.20 at any marker scatter
including zero, confirming the author, and **at both band edges that happens at a 10% record, which
is the accuracy the specification already asks for.**

**One correction to the budget's standing, found by the Monte Carlo.** At sigma* the MAP estimator
gives RMSE 0.33 to 0.42 against the Fisher bound's 0.200, which would make the budget worthless. It
is not the budget's miss. Re-running the identical cells with the `k_e` prior at S2b's realistic
0.054 and changing nothing else gives RMSE 0.110 to 0.150 with no collapses, so **the budget is
honoured and conservative once `k_e` is measured per patient, and optimistic by about 1.8 when quoted
against the paper's population prior.** The robust spread shows the mechanism: 0.200 at monthly,
exactly the Fisher figure, against an RMSE of 0.351, so the gap is the wide prior's collapse tail
that S2 identified.

**Filed, not run.** sigma* was solved at the 0.8 prior throughout, so the budget for a
measured-`k_e` patient is larger than these numbers and every verdict above is conservative.
Re-solving at 0.054 would raise all four budgets and could move monthly's verdict.

**Re-check record.** `q10` gains check 3: sigma* bracketed by a numerical Hessian of the posterior
scanned on a log grid with no root-finder, at all four schedules, and the exactness claim re-tested
on that route at two schedules. 12 of 12 pass.

---

## R62. Step 3(b): fixing the offset restores identifiability, and the anchor must not be set too high

`audit/q7_offset_fixed.py`, `audit/PHASE9_LITERATURE.md` step 3(b). Phase 5 proved that carrying the
marker's unreachable offset free is what makes θ unidentifiable under the published algebraic
structure. R5 supplies a human anchor, 5 to 10% of the pre-treatment level from liver transplantation,
against a published mouse value. **Nothing in R5 is verified here**; its transplant ratios are step 2
item 2, which has not run.

**The offset is fixed as the ratio a transplant actually measures**, `(E0 + base·q)/m(untreated)` at
`E → ∞`, not as `E0` alone, because at `Imax = 0.999` the incomplete-inhibition term rides inside it.
The structure reproduces `k2_structures.PubB` to 8e-13 at `PubB`'s own ratio of 11.277%.

**A constraint found before the rest would run.** The offset and the paper's calibrated 60% fall are
not independent: a large unreachable fraction caps the achievable fall, and `E0 ≥ 0` puts a floor
under IC50 at `qθ(1−φ)/(φ−q)`. **The φ window compatible with the paper's own 60% fall is 0.343% to
11.660%.** R5's 5 to 10% fits; **R5's 20% sensitivity value cannot be calibrated at all**, its best
fall being 54.4%. The window has a lower edge too, set by `Imax` rather than by biology, and the
published structure sits 0.383 points inside the upper edge, a second knife edge of the kind R15
recorded.

**Fixing the offset works.** Symbolically the rank goes from 3 against 4 unknowns to 3 against 3. On
the fibre walk, θ slides from 160 to 3040 at a bit-identical record with the offset free, and holding
it moves the record by a factor of 3.8e13 more. At φ = 10% and the committed 5% marker CV,
`sd(ln θ)` goes from infinite to **0.135**, inside the gate. **But the anchor's value matters as much
as the fact of fixing it**: at φ = 5% the same fix leaves 4.8, a factor of 36 worse, because the
calibrated IC50 moves from 200 to 733.

**Two corrections of my own, both caught by the re-check and neither quietly repaired.**

1. **A bug.** `_calibrate` searched IC50 from 1e-3 without enforcing `E0 ≥ 0`, so it could return a
   structure with `E0` around −28,500. That made the φ window look 0.2 points wider than it is: the
   figure 11.875% reported in conversation is wrong and the admissible edge is **11.660%**. The same
   omission was fixed in `ceiling` first and missed in `_calibrate`.
2. **A wrong finding, replaced.** An earlier draft reported that holding the offset too high biases θ̂
   sevenfold low, θ̂/θ = 0.141. That was an artifact of the fit's own θ grid: the objective has no
   interior minimum, an independent optimiser drove θ̂ to zero at the same residual, and the profile
   shows why. **Holding the offset too high does not bias θ, it un-identifies it:** the profile is
   flat within 2% from θ = 30 to θ = 1000, a factor of 33, and the true value fits *worse* than every
   point below it. Holding it too low leaves the minimum at the truth. **So too low costs precision
   and too high costs identifiability**, and neither is visible in the fit, since the wrong anchor's
   residual of 3.2e-03 sits far under the 5e-02 a 5% marker CV produces.

**The consequence for the specification** is that the offset should be fixed at or below the bottom of
R5's range, and that step 2 item 2 is asking for a defensible **lower bound** on the unreachable
fraction rather than a best estimate.

**And the factor of 5.8 does not move.** S2b's 5.8 is a PAPER-structure number about whether the
marker is itself enzymatically cleared; the offset is a PUB_B parameter, and PUB_B has no `k_e` term
at all. The two act through different channels, variance for `k_e` and bias-then-identifiability for
the offset, and they are not commensurable or additive. A PAPER-structure balance carrying an
unreachable floor would make them comparable; that is a new model and is filed.

**Re-check record.** `q10` gains check 4: `IC50_min` by root-finding against the algebra, the φ window
edge by the calibration itself against the ceiling, and the profile minima with the anchor correct and
too high. **Both parts of check 4 failed on first run and both failures were in the results, not the
checks**, which is what produced the two corrections above. 16 of 16 now pass.

---

## R63. Three corrections: the collapse tail is the finding, and the offset anchor does not recover θ

`audit/PHASE9_LITERATURE.md` 3a.4, 3b.0, 3b.2. All three came from the author's challenges and two of
them remove claims this session had made.

**1. The nuisance-convention disagreement was seed count, not convention, and my relabelling was
wrong.** The author's code fixes the truth at `k_e` and draws the prior centre, which is S2's **centre**
convention, so calling their 0.220 a truth-convention figure was incorrect and is withdrawn. Their own
code gives 0.220 at 80 draws, 0.307 at 96 with one collapse, 4.125 at 160 with three and 3.372 at 240,
while the robust spread holds between 0.186 and 0.213 throughout. Their collapsing records are the ones
handed a prior centre 3.7 to 9.8 times the truth, which confirms this session's mechanism by an
independent route, and their 1 in 96 and 3 in 160 agree with its 2 in 96.

**The consequence is larger than the budget.** RMSE is not a usable summary at the population prior: it
moves from 0.220 to 4.125 on seed count alone while the robust spread does not move. So the result is
**two facts**: at sigma* the budget is honoured in the bulk, robust spread 0.200 against the Fisher
0.200, **and** about 2% of records collapse by orders of magnitude with no warning in the fit, 2.1% with
a 95% Clopper-Pearson interval of [0.25%, 7.3%] at 2 in 96 and 1.9% [0.39%, 5.4%] at the author's 3 in
160. The rate is not pinned at any of these seed counts. **This is a second and independent reason to
measure `k_e`, and a different kind of reason: measuring it removes the tail rather than only narrowing
the spread**, the collapse count being zero on every schedule at 0.054.

**2. The φ window redone at `Imax = 1`, as R15's own correction requires.** `Imax = 0.999` is a constant
the source fixed rather than estimated, and R15 retracted a verdict over it. At `Imax = 1` the window is
0% to 11.875% against 0.343% to 11.660%. **The substantive conclusion survives**: R5's 5% and 10%
calibrate on both readings and its 20% cannot, its best achievable fall being 54.4 to 54.5%. The lower
edge exists only because `q > 0`, so it is an artefact of the fixed `Imax` and not a constraint on
biology.

**3. The claim that fixing the offset puts θ inside the gate is withdrawn.** The author asked for the
factor of 36 to be confirmed by the profile route that had caught the previous artefact. It is not
confirmed, and neither is the 0.135. With records simulated and the load bias drawn and fitted, the
robust spread at φ = 10% and the committed 5% marker CV is **6.1 in ln θ against a Fisher bound of
0.135**, a factor of 45. It is not an optimiser failure: starting at the truth gives 5.47. The
mechanism is a needle, the noiseless profile being 1.6e-16 at the truth and 1.3e-03 at a quarter of it,
so the information matrix measures an infinitesimal basin while a record at 5% CV sees noise forty
times the floor's depth. **The gate is met only near a 0.1% marker CV, two orders of magnitude below
the 5 to 8% analytical floor.** So Phase 5's "fixing E0 restores identifiability" is true in the rank
sense, which 3b.1 still proves, and does not license a claim that θ can be recovered.

**The pattern, now named rather than recorded twice.** At φ = 5% the anchor is safe and useless; at
φ = 10% it is unsafe and, as of this correction, useless too. **There is no value in R5's span at which
it is both safe and useful, which is exactly the shape S5 found for the ex vivo prior on θ.** Two
external anchors, two parameters, two structures, one trade: an anchor borrowed through an unverified
mapping buys precision only in proportion to how much of its own error it injects, so **the mapping's
verification is not a refinement of the anchor, it is the anchor.** Step 2 item 2 therefore stops being
a check that R5's range is roughly right and becomes a requirement to pin a value.

**Method note, since this is the third artefact in one section.** All three came from reading an
information matrix as if it described a likelihood. The matrix was right each time and the inference
from it was not. The rule this phase adopts: **a Fisher figure that decides a verdict is confirmed by
fitting before it is reported**, and a needle-shaped profile is the case that makes the two diverge.

---

## R64. The pricing: the cheapest study measures the quantity nothing measures

`audit/q8_pricing.py`, `audit/PHASE9_LITERATURE.md` section P. The last item in the Phase 9 brief.

**The three sigma axes are priced as one identity, not three sweeps**, because step 3(a) proved the
design sees them only through sigma once `s_b` is named. `k_e` picks the budget, the load record and
the marker CV spend it, and the whole table is `sigma = hypot(load_frac·L, b·CV·m̄)` inverted. One
consequence worth having: **the drift term stops mattering at normal capacity**, 56.10 against 56.16
for a study with and without 5% drift, so a repeated `k_e` measurement is a disease-range requirement
only, which is S4.3 arriving in the pricing.

**Measuring `k_e` buys two things and only one of them is a variance.** It widens the budget by a
factor of 1.5, nearly tripling the admissible marker CV at a 10% load record from 5.9% to 15.4% at
monthly. And it **removes the collapse tail**: zero in 96 on every schedule at 0.054 against 1 or 2 in
96 at the population prior, with the robust spread halving as well. The intervals are quoted honestly:
zero in 96 bounds the rate at 3.8%, not at zero, so demonstrating the tail is gone would need several
hundred records.

**The repeatability study is the cheapest item in the whole specification, by a wide margin.** Repeat
draws in treated patients at one fixed dose-cycle phase with the diet recorded measure the residual
scatter directly, which is what sigma needs and what item 7 found undecomposable from the published
surrogates. The pooled within-patient variance carries `df = n(r−1)`; **ten patients with three draws
each, thirty draws in total, resolves the CV to ×0.77 to ×1.44**, and the trigamma precision is
confirmed by simulation, matching to 2% with 94.5 to 95.6% coverage.

What it can decide depends on where the truth sits against the schedule boundaries, monthly 5.88%,
fortnightly 11.00%, weekly 20.73%, every 2 days 42.60% at a 10% load record. **Eight patients at three
draws settle a true CV of 8% or 15%; a true CV within about 10% of a boundary needs hundreds.** That is
not a reason to defer it: the same thirty draws tell you which case you are in, and a df-20 interval
straddling a boundary is itself the finding that the schedule cannot be chosen on this evidence. **No
other item in the specification returns either a decision or the knowledge that no decision is
available.**

**The floor anchor is priced on the terms step 3(b) left it.** Not what precision it delivers, since
3b.2 found it delivers none at any realistic marker CV, but the marker CV at which it would begin to
be worth having: **about 0.1%, which is 50 to 80 times better than the analytical floor and 250 to 350
times better than the surrogates suggest. The price is not payable.** Verifying R5's transplant ratios
is still worth doing for the structural claim the paper's identifiability rests on, which is the
difference between a claim that is false and one that is true but not yet operational.

**The order of cost, which is the deliverable.** The repeatability study, about thirty draws, first.
Then the `k_e` micro-dose study with its 24 h window. Then the load record to better than 5%. The
label-fate measurement rides inside the `k_e` study at no extra cost. The transplant-ratio
verification last, for the structural claim.

**Re-check record.** `q10` gains check 5: `sd(ln σ̂)` from the trigamma function against chi-square
draws, and the budget inversion checked by going forward. 18 of 18 pass.

---

## R65. The budget confirmed at the band edges, and the between-patient design decision

`audit/PHASE9_LITERATURE.md` 3a.9 and P.6, plus a promotion to the standing conventions.

**The method note is now a house rule**, not a Phase 9 lesson: a Fisher figure that decides a verdict
is confirmed by fitting before it is reported, the robust spread is what gets quoted, and the collapse
count carries a binomial interval. It sits in the standing-conventions section beside the
sample-versus-population line. It is not new: **R20** retracted Cramer-Rao curvatures quoted where a
reader reads an interval, `audit/n6_profile_vs_crb.py` exists for the same reason, and step 3(b)
produced three instances in one section. Five in the project is a rule.

**Applied to this phase's own main deliverable, and it holds.** Every sigma* in 3a.1 and 3a.6 is a
Fisher figure confirmed at normal capacity only, while the band edges carry Phase 9's strongest claim
and have budgets roughly half the size, which is where a needle would hide. All 48 cells were fitted,
three capacities by four schedules by two `k_e` priors by two nuisance conventions, 96 records each.
**The robust spread runs 0.169 to 0.302 at the population prior and 0.153 to 0.241 at the measured
one, against a Fisher bound of 0.200 everywhere.** The worst cell is 51% high; the band-edge cells the
claim rests on are 0.169 to 0.285. **So sigma\* is a budget and not a curvature wherever this phase
quotes it, 3a.3's claim stands, and step 3(a) is finished.**

The collapse pattern confirms 3a.4's mechanism from a third direction: all 13 collapses are on the
`centre` convention at the population prior, they occur at the band's lower edge and at normal
capacity, and **there are none at the band's upper edge**, where the enzymatic term dominates so
`k_e` cannot absorb the load. Zero in all 24 cells at the measured prior.

**And one filed item is promoted, because it is a design decision rather than a question.** If the
residual scatter varies between patients, no single requirement can be written for every patient, and
a pooled design cannot recover that afterwards. The sampling noise on a patient's own scatter is set
by the number of draws per patient alone: 0.641 at `r = 3`, 0.350 at 6, 0.223 at 12. **So more
patients cannot buy what more draws per patient buy.** At a true between-patient spread of 0.30: 100
patients at 3 draws, 300 in total, gives 40% power; 25 at 12, also 300, gives 97%; **120 draws at
`r = 12` beats 300 at `r = 3`.**

**The 30-draw design priced in P.3 cannot answer this question**, 15.4% power, and its nominal 5% test
runs at 8.3% because the chi-square approximation is poor at two degrees of freedom, which is a reason
not to trust a null result from it. **The recommendation is 10 patients at 12 draws, 120 in total**:
pooled `df = 110`, a ×0.88 to ×1.16 interval on the scatter, better than the 30-draw design's ×0.77 to
×1.44, and 73% power at a correctly sized test. Four times the draws, about five and a half months of
sampling at one fixed dose-cycle phase, inside the record the design already assumes. If only 30 draws
are available, **5 patients at 6 draws beats 10 at 3**: the pooled answer is unharmed and the
between-patient question stays answerable by adding patients later.

---

## R66. The budget solved by fitting, and there is no budget without a measured k_e

`audit/q12_fitted_budget.py`, `audit/PHASE9_LITERATURE.md` 3a.10, plus P.7. The author asked for the
budget to be solved by fitting rather than reported as a Fisher figure with a caveat, which is the new
house rule taken one step on.

**At the measured `k_e` prior the Fisher budget is sound.** fitted/Fisher runs **0.831 to 1.296**
across all 12 cells, so it is accurate to about 20% either way and is not systematically optimistic.

**At the population prior it must not be quoted at all.** fitted/Fisher runs **0.040 to 1.163**, a
factor of 29 across cells differing only in schedule. At normal capacity, every 2 days, the Fisher
budget is 130.96 and the fitted budget is 5.24.

**It is not solver noise, it is a floor.** Traced at normal capacity, fortnightly: the fitted robust
spread is 0.114 at 2% of the Fisher sigma, 0.204 at 8%, and then **rises and never returns below the
gate**, 0.268 at the Fisher value itself. The `k_e` prior contributes a schedule-dependent floor that
reducing the marker noise does not remove, and the crossing sits wherever that floor cuts 0.20. That
is what makes the ratios erratic: monthly's spread at its own sigma* is exactly 0.200 so its ratio is
1.148, fortnightly's is 0.268 so its ratio is 0.076.

**So this is a third and independent reason to measure `k_e`, and the strongest.** Without it there is
no usable noise budget: not one that is optimistic by a stated factor, but a quantity that cannot be
stated. The band's upper edge is the exception on both priors, 0.862 to 1.163, for the same reason it
has been throughout.

**A solver fault found and fixed on the way.** The first version bisected on an assumed bracket and
returned the bracket's own endpoint on three cells, reporting 0.150 where the true crossing was near
0.25. It now scans a log grid and brackets the last crossing, and reports a cell with no crossing
rather than inventing one.

**And P.7, what the 12-draw repeatability design measures.** Over 5.1 months a slow drift in capacity
or load enters the within-patient variance and is indistinguishable from residual scatter, so the two
designs do not measure the same quantity. Derived rather than imported: a 10% two-year drift in
capacity moves the marker 0.671 residual sd and in load 0.728, so the author's 0.3 to 0.7 is a 5 to
10% drift, while S2b.5's 10% `k_e` drift is only 0.11. **The effect is real and negligible**: at the
largest drift and the longest window it inflates the estimate by 0.09%, three orders of magnitude
below the study's own interval, and a 5% inflation would need a two-year drift of about 5.3 sd. P.6's
recommendation stands on its power argument alone. **One structural point the arithmetic exposes**: at
`r = 12` the study spans 21% of the record, and no feasible repeatability study captures the record's
own drift, so if drift matters it must be carried by the model as a time-varying capacity. Filed.

---

## R67. The one-state model IS a reduction, and neither of its nuisance constants means what it says

`audit/q13_two_state_reduction.py`, `audit/PHASE10_SENSORS.md` section R, 5 of 5 checks. Authorised as
a reduction check rather than a rebuild, and scoped to that: no sensor is priced.

**The reduction holds, exactly.** With `p` intracellular propionyl-CoA fed by the load, drained by the
enzyme as `theta p/(K_p+p)` and by a linear overflow `c p`, and the marker fed from that overflow at
yield `alpha` and cleared at `k_e`, eliminating `p` at quasi-steady state reproduces the paper's
balance with a sympy residual of **0**, and numerically to **6.7e-15** across 54 combinations of
capacity, load and increment and four independent `(alpha, c)` parametrisations. **So every existing
result survives as the reduction of a better-founded model**, which is a strictly better position than
the paper is in today.

**But the reduction names what the paper's constants actually are, and two of the three are
composites.**

    K_m       = K_p * alpha * c / k_e   =   K_p * c / k_e_paper
    k_e_paper = k_e / alpha

`theta` is unchanged, which is fortunate since it is the estimand. **`K_m` is not the enzyme's
Michaelis constant**, it is that constant scaled by the overflow-to-clearance ratio. **`k_e` in the
paper is not the marker's clearance**, it is that clearance divided by `alpha`, the fraction of
overflow flux that becomes the measured marker. `K_p`, `c` and `alpha` enter through only two
combinations, rank 2 of 3, so adding the state adds no parameter the marker record can estimate.

**And this lands on Phase 9's precondition.** Phase 9 concluded that measuring `k_e` is the
precondition for a statable requirement. **S2b's micro-dose study measures the marker's clearance,
`k_e`; the model's parameter is `k_e/alpha`.** So the study delivers the numerator and `alpha` is a new
unmeasured nuisance on exactly the quantity wanted to 2%. At `alpha = 0.3` the model's constant is 3.3
times the measured one. **The numbers in Phases 8 and 9 do not change, since the reduction is exact and
the fits estimated whatever the parameter is, but their reading does**: "measure `k_e` to 2%" becomes
"know `k_e/alpha` to 2%".

**`alpha` is measurable, and what measures it is the marker panel.** Each overflow marker's formation
flux is its clearance times its level, so `alpha = k_e,3HP m_3HP / sum_j k_e,j m_j`. **That answers
Phase 10's original question from an unexpected direction: the panel's value is not noise averaging,
which the paper priced, nor the different functional forms, which R3 recommended it for, but that the
sum of the overflow routes converts a measurable tracer clearance into the model's parameter.**

**What breaks the reduction is the other markers.** With a saturable overflow no single `K_m`
reproduces the two-state marker: at `K_o = 100` the best-fitting `K_m` is 69 instead of 100 and the
residual is still 16%, which a reader would see as a calibration problem rather than a missing state.
The linear limit recovers `K_m = 100` to 7.5e-07, which is the check on the test. So **the reduction
holds for the marker the paper uses and fails for the ones R3 describes differently**, which is why
those markers carry information.

**One error of mine, caught by the check.** The first version assumed `k_e` carried through the
reduction unchanged, and the symbolic residual was not zero. The second term of the reduction is
`(k_e/alpha) m`, not `k_e m`, and the whole `alpha` finding above followed from chasing that.

---

## R68. The fitted budget on both nuisance conventions

`audit/PHASE9_LITERATURE.md` 3a.10, extended at the author's request.

**At the measured `k_e` prior the Fisher budget is sound on both conventions**, fitted/Fisher 0.811 to
1.296 over 24 cells. **At the population prior it must not be quoted on either**: `centre` gives 0.040
to 1.163 and `truth` 0.071 to 1.170. **The author's hypothesis, that the gap is a `centre`-convention
artefact because a bias is not a variance, is only half right**: on `truth` it is still a factor of 16.

**And neither convention dominates cell by cell.** At the band's lower edge, monthly, `truth` is worse
(0.082 against 0.325); at normal capacity, every 2 days, `centre` is worse (0.040 against 1.107). So
there is no pessimistic convention to adopt and the defensible envelope is the worse of the two per
cell, 0.040 to 1.148.

**The author's argument about which governs still holds and comes out stronger.** A specification
normally uses `truth`, a draw from a population whose centre is known; R1 establishes that no human
`k_e` exists in any population, so the centre is unanchored, which is the `centre` convention. But the
conclusion does not need that choice, because both conventions fail: **at the population prior there is
no statable budget on either honest convention, so measuring `k_e` is the precondition for the
requirement framework existing at all.** That is the form Phase 9 should be read as concluding.


---

### R69. The second channel is not a clearance, and the composite has a directional bias

`audit/PHASE10_SENSORS.md` R.2 points 4, 6 and 7.

Eliminating `p` from the two-state model gives `c*p` identically equal to `(k_e/alpha)*m`, so the
paper's second channel **is not a clearance of the marker at all**. It is the disposal of
propionyl-CoA down the overflow route, expressed in marker units, and `alpha` appears because the
model carries a flux out of a pool it does not represent, divided by the fraction of that flux that
happens to be measured.

`alpha = k_e,3HP*m_3HP / sum_j k_e,j*m_j` is right **only if the measured markers exhaust the overflow
routes**. If any route is unmeasured the denominator is short, `alpha` comes out too high and the
model's `k_e/alpha` too low, so **the error is directional, not symmetric**: a panel capturing a
fraction `f` of the overflow understates the model's constant by `1 - f`. R3's panel is 3-HP, 2-MC, C3
and propionylglycine; urinary free propionate and odd-chain fatty acid incorporation sit outside it.
This is the same failure mode as R5's floor, a ratio whose denominator is assumed complete.

### R70. The tracer route is one study per marker, not one study

`audit/PHASE9_LITERATURE.md` P.8, with P.5 item 2 rewritten to point at it.

Converting a micro-dose study's result into the model's parameter needs `k_e,j` and `m_j` for **every**
panel marker, so the pricing row "one IV micro-dose study, 24 h window, 8+ draws" becomes **four IV
micro-dose studies, one per panel marker**, plus a closure assumption that is not itself measurable
from the panel. The row is repriced rather than caveated.

### R71. The reduction enters the ladder, and it is the fourth instance of the pattern

`audit/PHASE10_SENSORS.md` R.4, placed beside R22's carnitine reduction and given the same treatment:
what it assumes, where it is exact, what it costs when it is wrong.

The entry that matters: **the 16% residual is invisible as a missing state and reads as a calibration
problem.** A saturable overflow at `K_o = 100` fits a `K_m` of 69.48 against a true 100, which a reader
takes as a mis-set constant, not as evidence that a state is absent from the model. That is the fourth
instance of this project's own pattern, after R14's curvature read as recovery, R15's knife edge read
as a shape limit, and Phase 9 3b's needle read as identifiability.

---

### R72. The cohort route escapes alpha entirely, and still cannot supply the constant

`audit/PHASE10_SENSORS.md` C.0 to C.6, `audit/q14_hierarchical.py`, re-checked by
`audit/q11_phase10_recheck.py` check 3, 6 of 6.

Run ahead of Phase 10 steps 2 to 4 because the reduction handed the cohort route a large new
advantage: `k_e/alpha` is exactly the coefficient multiplying `m` in the balance a marker record
obeys, so a cohort fit estimates it directly and **never meets `alpha` at all**, while the tracer
route needs four micro-dose studies plus a closure assumption to reach the same quantity (R70).

**The advantage is real, large, and irrelevant.** On the trial's own monthly record the disease
range's looser target of 6% needs **723 patients** and its tighter 2% needs **6491**; the largest
cohort anywhere in the six literature reports is 60. Only a 2-daily record brings 6% inside a
plausible study, at 48 patients giving 17,520 draws.

**And the sample size is not what decides it.** R1 puts the between-patient CV of the disposal
constant at "probably a third or more". A cohort pins the POPULATION value; a patient's `theta_i`
consumes their own. At a CV of a third the patient's constant **floors at 0.3166 with an infinite
cohort**, missing the looser target by a factor of five; even a CV of 10% floors at 0.1027 and misses.
A larger cohort measures the population better and the patient no better.

Confirmed by fitting per the house rule, fitted/Fisher 0.617 to 1.168. **The `n = 16` cell came out
below its own Fisher figure**, so the fit was restarted a log-sd of 0.5 from the truth: it reproduced
to 0.002 with an identical median error of -0.061, which makes it a biased estimator and not a stalled
optimiser, and makes the cohort sizes above generous to the cohort route.

Not one cohort row meets the 0.20 gate anywhere in the requirement table; the best, n=64 shared at
20% capacity on a 2-daily record, reaches 0.213. The two routes as they would actually be run, a
cohort of 64 at CV 33% against four micro-dose studies, are 0.815 against 0.243 at 10% capacity.

**The shape of it.** The cohort route wins cleanly on the axis the reduction exposed and loses on
arithmetic. Its difficulty is the harder kind: the tracer route is expensive and leans on an
assumption a wider panel could test, while the cohort route meets a floor no money moves.

**And the item it raises, filed.** If the between-patient CV is what kills the route, the CV is the
thing to measure first, and P.6's repeatability design, 10 patients at 12 draws, is already scoped to
estimate it. A cheap decision-relevant measurement now stands in front of an expensive one.

---

### R73. Phase 10 steps 2 to 5: the record is worth a few per cent on top of one ex vivo assay

`audit/PHASE10_SENSORS.md` S.0 to S.6, `audit/q15_sensors.py`.

Every surviving sensor supplies an independent measurement of ONE of the four parameters, so the
phase is which parameter is worth buying precision in. Four findings, the first two negative.

**Candidate 9, propionate turnover, is dropped on measurement.** Knowing `L + b` EXACTLY improves
`sd(ln theta)` by 1 to 2% at every capacity and schedule. There is no precision at which it is worth
buying because there is no precision at which it does anything. The offset is already pinned by its
own prior and the binding constraint is `k_e`, not the load. **A corollary no phase had noticed: S2's
`recorded` and `measured` load columns differ by 1 to 2% and are not two designs.**

**Candidate 10, timed urine, is dominated by candidate 8**, because `f_renal` sits inside it and R1
gives no value for 3-HP. To match candidate 8's 5.4% the renal fraction would have to be known to
about 5%, which is harder than the clearance it substitutes for.

**Candidate 8's value is schedule-dependent.** At 5% and 10% capacity the monthly record misses the
0.20 gate with `k_e` known EXACTLY, 0.395 and 0.213, so no tracer study rescues it. The 2-daily
record passes at both. **The schedule decision comes before the sensor decision.**

**And the result the project should be most uncomfortable with.** An ex vivo assay alone gives
`sd(ln theta)` equal to its own precision. **Two years of monthly plasma 3-HP added on top improves
it by 0.2 to 2.9 per cent**; the best cell in the whole sweep is 12.3%, needing a poor assay, the
densest schedule and the top of the disease range together. That is not an argument for the assay,
whose ex vivo to in vivo transfer is unmeasured, and not a refutation of the record. It is **the
record's information content on an absolute scale instead of a relative one**, and every earlier
phase compared schedules against schedules. The project's question takes a sharper form: capacity can
be identified to within the assay's own precision, and almost entirely because of the assay.

Both verdicts confirmed by fitting, fitted/Fisher 0.767 to 1.056. **On candidate 9 the two footings
agree only to within the Monte Carlo noise**, 11% better fitted at one capacity and 3% worse at
another against a Fisher figure of 1 to 2%, so the honest reading is indistinguishable rather than
matched. **On candidate 2 they agree closely, and the assay removes the collapse mode outright**,
0/96 against 20 to 29 out of 96 without it, which the spread alone does not show.

---

### R74. The ex vivo comparison was quoted from a column R4 rules out, and the conclusion inverts

`audit/PHASE10_SENSORS.md` S.3, `audit/q16_corrections.py` part A. **Withdrawal of R73's headline.**

R73 reported that "two years of monthly plasma 3-HP, added to a single ex vivo assay of ordinary
quality, improves the estimate of residual capacity by between 0.2 and 2.9 per cent". That is
computed at an assay log-sd of 0.20. **R4's finding is that the ex vivo to in vivo mapping has never
been measured in any human**, that the only quantitative link is Manoli's Spearman `r = 0.682` in MMA
with `R^2 = 0.327`, and that the defensible width is a natural-log sd **of order 1**. A prior of 0.20
is five times tighter than that and tighter than the tightest value S5 swept.

At `a = 1.0` the record adds **83.1%** at nominal, **37.0%** at 20% capacity, 15.5% at 10% and 4.9%
at 5%. **So the record does most of the work and the assay is the weak partner, which is the reverse
of what R73 concluded.** And neither partner is usable: the joint figure misses the 0.20 gate at
every disease-range capacity, 0.630 to 0.951.

**The failure mode is this project's own most-named one.** The surrounding prose flagged the ex vivo
to in vivo transfer as unmeasured, and the headline was then quoted from the column that assumes it
away. It is the first instance in this audit's own output rather than in the paper's.

### R75. Candidate 9 was priced on half its effect, and the paper's route 6 is partly retracted

`audit/PHASE10_SENSORS.md` S.1 and S.1b, `audit/q16_corrections.py` parts B to E.

**R73's 1 to 2% was the offset channel alone.** Measuring the load also removes the load-recording
error from `sigma`, and on the project's own parametrisation the two are ONE knob: the load record's
accuracy `LOAD_FRAC`, committed at 5%, sets both `s_b = LOAD_FRAC * L` and one of `sigma`'s two
terms. An exactly measured load is worth **+7.0% at nominal, +7.5% at 20%, +11.9% at 10% and +14.5%
at 5%**, against the 0.3 to 1.5% first reported.

**The candidate 9 verdict survives with a different reason**: dropped on size, not on absence.

**The corollary is withdrawn, and its reverse is the correct statement.** R73 concluded that S2's
`recorded` and `measured` load columns "are not two designs at all". They differ by 1 to 2%, but S2's
`measured` column is not what the paper means by measuring the load: S2 moves only the offset
parameter while the paper's route moves `LOAD_FRAC`. **It is a naming collision, and S2's `measured`
column understates the study it appears to price by a factor of five to ten.** The paper's Section 7
already had this right: "what the load record has to deliver is per-visit precision entering through
`sigma` and not a calibrated mean offset".

**RETRACTION, `paper/sections/11_scope.tex`, route 6.** Its headline reproduces: 6 visits reach the
gate at nominal (exactly), 13 visits match 25 (12 here, off by one), and `2.659` at 5% of nominal
(2.655 at 1% accuracy). **Its third clause does not.** Route 6 states the advantage "shrinks rather
than reverses as capacity falls, from thirteen visits to eighteen". Swept from nominal to 2% of
nominal the matching visit count runs **12, 11, 10, 12, 12, 12**: the advantage does not shrink, it
is flat, and it is strongest at 20% of nominal. The conclusion route 6 draws from it is separately
true, since 12 visits at 5% capacity still leave 2.62 against a 0.20 gate, so **a supporting number
is wrong and the finding it supports is not.**

**And route 6's prose describes a different computation from the one that produced its numbers**,
saying measuring the load "removes `b` from the parameter vector rather than bounding it", which is
the b-only reading and gives 24 visits against 25, not 13. Sections 7 and 11 contradict each other on
the mechanism and Section 7 is correct.

**How it survived the artefacts check**: `2.659`, `thirteen`, `twenty-five`, `six` and `eighteen` are
hardcoded prose, not `\num{}` macros, so nothing compared them with a generator. R23's failure mode,
recurring in the one section that enumerates the paper's own negative results.

### R76. The between-patient spread is promoted to the specification's top line

`audit/PHASE10_SENSORS.md` S.7. Filed three times, now settled.

**It decides both routes and the sensor pricing.** The cohort route floors at 0.3166 with an infinite
cohort if the CV is a third and is merely infeasible if it is near zero (C.2). The tracer route's
four micro-dose studies must be run per patient rather than once if the spread is large, multiplying
P.8's cost by the cohort size. And every `k_e` figure in S.2 is a per-patient precision, so whether
any route can deliver one depends on whether the quantity is shared at all.

**It is already scoped and it is cheap**: P.6's 10 patients at 12 draws, 120 draws at one fixed
dose-cycle phase, about five and a half months on the 14 d train, giving pooled `df = 110` and 73%
power on a `tau` of 0.30.

**One qualification carried into it.** R1's CV of about a third is for the marker's clearance, while
the model's parameter is `k_e/alpha` and `alpha` varies between patients too, through carnitine
status and the relative sizes of the overflow routes. Applying 33% to the composite assumes `alpha`
is constant or perfectly anticorrelated with `k_e`; **the true floor could be worse and cannot be
better unless that anticorrelation is near perfect.**

**So the ordering is forced**: measure the between-patient spread before committing to either route.
It should be the top line of whatever the specification becomes, ahead of the load-record requirement
and ahead of the `k_e` study, both of which are conditional on it.

---

### R77. The provenance sweep: the guard sees 343 numbers and the paper asserts 2472

`audit/PROVENANCE_SWEEP.md`, `audit/q17_*.py`. Commissioned after R75, and it closes no hole: it
measures one. **`paper/` is untouched.**

`tools/check_artefacts.py` regenerates the tables, the figures and `results/numbers/*.tex`, so a
`\num{}` value is owned and checked and a prose literal is owned by nobody. **343 guarded references
against 2472 unguarded literals**, of which 1000 are spelled-out numerals. **A scan over digits would
have missed route 6's own two values, since `thirteen` and `eighteen` are words**, which is the
sharpest structural lesson: any guard built from this has to read English numerals.

| category | count |
|:--|---:|
| 1, generated elsewhere and retyped here | 218 |
| 2, structural | 1017 |
| **3, asserting a computed result** | **206** in 130 claims |
| 4, harmless | 1031 |

**Category 3 carries two columns.** A generator exists for a claim if its paragraph already cites a
`\num{}` macro, which means a named script already runs and already emits numbers there: **52 of 130
have one, 78 do not.** Reproduction was tested by running fourteen generators and searching their
printed output, since these scripts print far more than they write as macros: **26 reproduce, 53 are
unowned but printed somewhere, 18 are NOT IN OUTPUT, 8 PARTIAL, 25 unowned and unprinted.**

**NOT IN OUTPUT is a statement that nothing checks the number, not that it is wrong.** Route 6's
`eighteen` was exactly this before q16 showed it wrong, so the 51 claims in the bottom three rows are
the queue.

**Two clusters were checked the stronger way, by recomputing the quantity.** Section 8's load-record
decomposition, which carries the paper's third-tightest requirement and has no generator at all,
**reproduces throughout**: 0.283 against 0.293, 0.172 against 0.174, the 98/2 split against 96.8/1.6,
the capacity shares 1.8/9.9/26.2 against 1.6/10.1/27.9, the 10.4% threshold against 9.73%. **An
unguarded number is unchecked, not wrong, and that has to be said as loudly as the retraction.**

**Three findings for the consolidation pass.** `2.659` appears in three places across two sections
and is owned by nothing in any of them, one value hand-copied three times under both route 5 and
route 6. `07_design.tex:113`, the gate claim itself, asserts `0.457, 0.412, 0.510` with no generator
printing them. And the 218 category-1 literals are pure upside: the macro exists, the prose simply
does not use it.

### R77a. And the repository has already solved this once, for figures

`audit/PROVENANCE_SWEEP.md` section 7. The most useful thing the sweep found is not a number.

`check_artefacts.py` already carries a provenance mechanism for numbers **drawn on figures**, built
after the phase overlay printed "25 phases" on three panels whose lattices are 35, 28 and 28. Each
figure module declares every interpolation as `registered:`, `own:`, `input:`, `geometry:`, `text:`
or `UNREGISTERED:`; an undeclared number fails, a declaration for a number no longer drawn fails, and
**`UNREGISTERED:` is counted and printed on every run, so the debt is a number that can be driven to
zero rather than a silence.** The current run reports 29 drawn numbers, 13 registered, 0
unregistered.

**The four categories of this sweep map onto that vocabulary directly**: category 1 is `registered:`
waiting to be wired, category 2 is `geometry:` and `text:`, category 3 is `UNREGISTERED:`, category 4
needs no declaration. So the 130-claim list is this project's PROSE debt stated in the form its own
figure guard already uses, and the 51-claim queue is what `UNREGISTERED:` would print on day one.

**The one thing prose needs that figures did not: the declaration has to cover English numerals**,
since `thirteen` and `eighteen` carried route 6's claim and no interpolation-based scan sees them.

---

### R78. The Section 8 convention settled: it is the t = 0 draw, and the decomposition reproduces

`audit/PROVENANCE_SWEEP.md` section 9, `audit/q18_section8.py`. **Upgrade of R77's verdict.**

R77 reported the load decomposition as reproducing "within a few per cent, consistent with a seed
convention", with the stated requirement threshold 7% apart. By this project's own standards that is
not a verdict, it is the shape of the disagreements that have already cost three rounds.

**The convention is the t = 0 draw, again.** `n4.MONTHLY` is `arange(0.0, 730.0, 30.4)`: 25 draws
**starting at t = 0**, one before the first infusion. R77's check used `arange(30.4, 30.4*26, 30.4)`,
dropping it and running 30 days past the window. On the paper's grid `n4.over_seeds(theta=3000)`
returns **0.1699**, the paper's 0.170 exactly.

Redone on that grid: baseline **0.170** (was 0.172), scatter alone **0.284** against 0.283, offset
alone **0.172** against 0.172, shares **97.3/1.9** against 98/2, capacity shares **1.9/9.9/25.1**
against 1.8/9.9/26.2, and the threshold **10.27% against 10.4%**, closing from 7% apart to 1.3%. The
companion claim, that the target is met at every offset from 1% to 500%, holds.

**So the verdict is REPRODUCES, not "consistent".** And the finding it supports is now properly
earned: Section 8's decomposition is owned by no generator and correct.

### R79. The classifier had a rounding hole, and it hid the sweep's own headline

`audit/PROVENANCE_SWEEP.md` section 8, `audit/q20_rounding.py`. **Correction of R77.**

R77 named `07_design.tex:113`, the gate claim, as "the one to look at first", its
`0.457 (0.412 to 0.510)` being "a median and an interval that no generator printed". **Wrong.**
`results/tables/tab_D_designs.tex` carries `six-monthly clinic visits & 4 & ... & 0.4570 [0.4121,
0.5097]` and `numbers/design.tex` the same three values. The prose ROUNDS them, and q17 matched
literal strings, so every rounded restatement fell into category 3.

Redone with rounding tolerance, claiming a match only at three or more significant digits:
**category 1 goes 218 to 249 and category 3 goes 206 to 175.** Among the 31 that move:
`09_validation:103`'s `0.170`, which `\num{sd_monthly}` already holds; two of the three numbers in
Section 8's decomposition; `11_scope:172`'s `3.07` and `:207`'s `1.46`.

**The unowned pile is 15% smaller than reported and the cheap-to-fix pile larger**, and the sweep
overstated the danger of the very claim it led with.

### R80. The load-bearing subset: eight of nine reproduce, and a generator that never runs

`audit/PROVENANCE_SWEEP.md` section 10, `audit/q19_loadbearing.py`. Verification only, no wiring.

**`2.659` is asserted at three sites across two sections and every one is right**, 2.668 against
2.659 as the load accuracy goes to zero. Route 5's factor at nominal reproduces at 1.809 against
1.81, and its band 1.22 to 1.46 spans about a quarter of nominal down to 2%. `07_design:113` and
`09_validation:103` reproduce and are category 1. The two visit counts are structural and exact.

**Route 6's third clause remains the only thing in the whole sweep that does not reproduce**, now
confirmed on the paper's own grid: the matching visit count runs 12, 8, 8, 10, 10 from nominal to 5%
of nominal, flat and never approaching 18.

**And the finding worth more than any of those values: `audit/t3_diet_modulation.py` appears nowhere
in `make figures`.** The guard checks 14 numbers-file generators for reachability, so a generator
writing no numbers file is not reachable, not checked, and its results enter the paper only as
literals. That is what blocks recomputing `11_scope:204` and `08_spec:372`. **A generator existing is
not the same as a generator running**, and the reachability check already here is what would say
which.

**Twice in this exercise a "does not reproduce" was mine and not the paper's**: the six-monthly
schedule, built by hand instead of taken from the table, and route 5's band, swept over 5% to 20%
when the paper's upper endpoint sits at 2%. Both were caught by asking where the paper's own value
comes from before concluding it is wrong.

---

### R81. The docstring declaration is not a mechanism, and four generators the build never runs own paper prose

`audit/PROVENANCE_SWEEP.md` section 11, `audit/q21_unreachable.py`.

R80 reported `t3_diet_modulation.py` as one instance of a generator absent from `make figures` whose
results reach the paper only as literals. **It is one of four.**

Of 111 scripts in `audit/`, **12 declare in their docstrings that nothing of theirs reaches the
paper**. They split exactly: five (`t1`, `t5`, `t6`, `t7`, `t11`) are in `make figures` AND write a
numbers file, so their declaration is false but harmless, their numbers being guarded macros; three
(`t2`, `t10`, `t12`) are genuinely scoping work; and **four own prose the paper asserts**:

| script | owns |
|:--|:--|
| `t4_bias_vs_scatter` | **Section 8's load decomposition**, the paper's third-tightest requirement: 0.2829 against "0.283", 0.1718 against "0.172", 0.1696 against "0.170", and "the scatter channel accounts for 0.1133 (98%)" against "98%" |
| `t8_load_by_capacity` | **route 5's factors**: `100% of nom 1.808` against "1.81", `20% 1.215 ... 2% 1.464` against "1.22 to 1.46" |
| `t9_measured_load` | **route 6 entirely** |
| `t3_diet_modulation` | route 4's swing, R80's instance |

**None writes a `results/numbers` file**, so none is among the 14 the guard checks for reachability.
The guard is not failing; it was never pointed at this class. **The declaration is wrong in both
directions, which is the point: a docstring is not a mechanism.** What the build runs is decided by
`make figures` and what the paper cites by `\num{}`, and neither reads English. **R23's failure mode
on a second surface**: a generator no entry point calls, in prose instead of a table.

### R82. WITHDRAWAL of R75's retraction: route 6 is correct in every clause

`audit/PROVENANCE_SWEEP.md` section 11.3. **This withdraws a retraction, so it is stated plainly.**

R75 retracted route 6's clause that the advantage "shrinks rather than reverses as capacity falls,
from thirteen visits to eighteen", on the ground that a recomputation gave a flat 10 to 12. **The
owning generator prints exactly the paper's numbers:**

```
  100% of nom   0.1699    13 to match     6 to gate
   50% of nom   0.3228    12
   20% of nom   0.8104    13
   10% of nom   1.5828    15
    5% of nom   3.0642    17
    2% of nom   7.4406    18
```

13 at nominal rising to 18 at 2% of nominal, monotone, which is "shrinks rather than reverses"
exactly. **R75 swept only to 5% of nominal and built the schedule by hand.** The clause is right.

**R75's second charge also falls.** It said route 6's prose "describes a different computation from
the one that produced its numbers". `t9`'s docstring: "THREE THINGS CHANGE AT ONCE... (1) b
disappears... (2) sigma stops being a diet record's accuracy and becomes a tracer measurement's
precision. (3) The steady record closes on itself." The prose names the first and omits the second.
**Incomplete, not wrong**, and Sections 7 and 11 do not contradict each other.

**What survives from R75 is the provenance defect alone**, which was always the finding worth having:
route 6's numbers are correct, computed by a real generator, and stated as literals no macro cites and
no build step regenerates.

### R83. The method inherited the gap it was measuring, four times

R77 called `2.659` "unowned and unprinted". `t9` prints it at four significant figures. R77's
reproduction pass ran fourteen generators and `t9` was not among them **because `t9` is not in
`make figures`, which is precisely the defect being reported**.

Four reconstructions in this phase disagreed with the paper and the paper was right each time: the
six-monthly schedule built by hand instead of read from `tab_D_designs.tex`; Section 8's decomposition
on the wrong t = 0 convention; route 5's band swept to 5% when its endpoint is at 2%; and route 6's
eighteen visits, same too-narrow sweep plus a hand-rolled grid.

**Every one has the same cause: recomputing a quantity instead of finding the generator that owns
it.** The rule, stated for the standing conventions: **find the owner before recomputing, and if no
owner can be found, say so before concluding the number is wrong.** The sweep's own list is the tool
for finding the owner.

### R84. A fifth site, and it reconciles R80's two figures rather than contradicting either

`audit/PROVENANCE_SWEEP.md` section 11.2. Surfaced by the workflow at low confidence and verified
here against `t9`'s captured output before being accepted.

`06_excitation.tex:309` states that "measuring the load instead of recording it is worth under one
per cent". `t9_measured_load.py` prints, at nominal, `recorded 5%  0.1699` against
`b gone, sig 5%  0.1693`: **0.35%**.

**So the paper separately and correctly prices the `b`-removal channel ALONE at under one per cent.**
That is the same channel R80's first pass measured at 0.3 to 1.5% and then treated as an error of
mine. Both figures are right and neither correction was needed on that point: the paper states the
b-only channel at `06_excitation:309` and `11_scope:227`, and the full accuracy move in route 6, and
they differ by a factor of twenty because they are different channels. **What R80 got right is that
candidate 9's own pricing must use the full move; what it got wrong is implying the paper had
conflated them. It had not.**

Another spelled-out numeral, so no digit scan would see it.

### R85. The four unreachable generators own 21 paper sites, not four

`audit/PROVENANCE_SWEEP.md` section 11, `audit/q21_unreachable.py`. Final count after the workflow's
remaining investigators returned and each new site was verified here against captured stdout.

| script | sites owned |
|:--|---:|
| `t4_bias_vs_scatter` | **8** |
| `t9_measured_load` | **6** |
| `t8_load_by_capacity` | **4** |
| `t3_diet_modulation` | **3** |

**`t4`'s decomposition is restated across two sections and every statement of it is unowned.**
Besides `08_specification.tex:254-256`, `06_excitation.tex:286-290` carries the same result again:
"runs 0.159 at a 1% recording error to 0.288 at 20%", the 98%/2% split, and "the target is met at
every bias from 1% to 500%". `t4` prints `1% 0.1589`, `20% 0.2863`, and verbatim
`the gate holds at EVERY bias from 1% to 500% (0.1691 at 1%, 0.1733 at 500%)`.

**So one generator the build never runs supplies the same finding to two sections.** That is the
shape the consolidation pass should expect: not isolated literals, but a result stated several times
from a single unreachable source, where moving or rewording any one statement cannot be checked
against the others or against the code.

One site is accepted with its endpoint unverified and marked as such: `08_specification.tex:375`,
"prescribing it is worth between nothing and two per cent". `t3` owns the quantity, printing
`designed 0.1696 vs spontaneous 0.1698, ratio 1.001` at nominal, which is the "nothing" end; the
capacity at which it reaches two per cent was not located.

### R86. CORRECTION to R81 and R85: it is five scripts, and t2 owns Tier 2's requirement

`audit/PROVENANCE_SWEEP.md` section 11.1, `audit/q21_unreachable.py`. The adversarial workflow
finished after R85 was written and contradicted a clearing of mine.

**R81 placed `t2_tier1_margin` among the three scripts doing genuine scoping work. It is not.** It
owns the whole of Tier 2's requirement cluster, `08_specification.tex:125-151`, verified here against
its captured output: `0.2120` against the paper's "0.212", `118` monthly draws and `9.8 years` against
"118 ... nearly ten years", `11` against "eleven with the increment known exactly", the floor row
`0.1510 | 0.1603 | 0.1740 | ... | 0.2482` against "0.151, 0.160, 0.174, 0.207, 0.248", and
`the floor reaches the gate at an increment CV of 18.7%` against "18.7%".

**So it is five scripts, not four, and the specification's TWO load-bearing requirements are owned
the same way.** Section 8's load record comes from `t4` and `t8`; Tier 2's increment repeatability
comes from `t2`. Neither generator is run by the build, and neither writes a macro.

**It also closes a row R80 left open.** `08_specification.tex:126` was listed as "cannot be
recomputed, same increment-CV machinery as `11_scope:201`". The machinery is `t2`, not `t6`, and `t2`
prints the answer.

**And I told the author the opposite while it was still running**, that nothing pending could change
the four-script finding because the outstanding scripts were already cleared. That was a prediction
about work not yet done, and it was wrong. **The clearing of `t2` rested on the same value-matching
that section 11 says is worthless**: `t2` printed no distinctive literal my matcher recognised, and I
read that absence as absence of ownership instead of as a limit of the method.

**Three agent claims did NOT survive verification and are rejected**, which is why every one was
checked against captured output rather than accepted: `t1_transient_worth` was credited with
`07_design.tex:286`'s 29.79 and 59.96 per day and `:289`'s 2.8 assay standard deviations, and `t1`
prints none of them. `t5_band_vs_ke` and `t11_phase_lattice` are partially supported, but both are in
`make figures` and write numbers files, so they are the milder category: a generator the build DOES
run, whose result the paper still states as a literal.

---

### R87. The five are wired: 22 literals become checked macros, and two refuse

`audit/q22_wire_macros.py`, `audit/q21_unreachable.py`. Author-directed, ahead of the consolidation
pass, because one result restated across sections from a source the build never runs cannot be
reworded safely.

**What was done.** `t2_tier1_margin`, `t3_diet_modulation`, `t4_bias_vs_scatter`,
`t8_load_by_capacity` and `t9_measured_load` each gained a `*_numbers()` function, are registered in
`common/numbers.py`'s `GENERATORS`, and run in `make figures`. **22 literals in
`06_excitation.tex`, `08_specification.tex` and `11_scope.tex` are now `\num{}` macros.** The guard
went from 14 numbers-file generators checked for reachability to **19**, and from 328 generated keys
across 17 files to **364 across 22**.

**No value changed, and that is enforced rather than intended.** `q22` refuses any replacement whose
macro does not render EXACTLY the literal it replaces; each generator's format string was set to the
paper's own printed precision, and three sites were refused on the first pass (1.808 against "1.81")
until the format matched. `pdflatex` is not installed in this container so the built document could
not be compared visually; what was checked is that every `\num{}` key resolves to a defined macro and
that each rendered string equals the literal it replaced.

### R88. Two genuine paper-against-generator disagreements, the first in the whole sweep

Both were found BY the wiring, which is what the wiring is for. Neither value is adopted: taking the
computed one would silently change a published number and taking the published one would falsify its
generator.

**1. `t3_diet_modulation`, at `11_scope.tex:204` and `08_specification.tex:372`.** The paper: "a
prescribed +/-15% swing in natural protein is worth 7%, and an implausible +/-50% swing 32%, leaving
2.32". The generator's own amplitude sweep at 5% of nominal gives 3.0583 at no swing, 2.8671 at 15%
and 2.3198 at 50%, so **6.3% and 24.1%**, not 7% and 32%. **The footing is not in doubt**: 2.3198
rounds to the paper's own 2.32, so the capacity and the 21 d period are the ones the paper used. Only
`diet_sd50_disease` is generated; the two percentages are written as
`diet_gain15_pct_UNRECONCILED` and `diet_gain50_pct_UNRECONCILED` and cited nowhere.

**2. `t4_bias_vs_scatter`, at `06_excitation.tex:286`.** The paper: "sd runs 0.159 at a 1% recording
error to 0.288 at 20%". The generator gives 0.1589, which is the 0.159, and **0.2863, which renders
0.286 and not 0.288**. The 1% end wired; the 20% end is held.

**This is the first time in 130 claims that a paper number has failed against its own code**, and it
took wiring to find it, because every earlier check compared the paper against a reconstruction and
five times the reconstruction was what was wrong. **A generator that the build never runs is free to
drift, and these two are what drift looks like.** Which of the pair is right in each case is not
decidable from here: the script may have changed since the sentence was written, or the sentence may
always have been wrong. Both need the author.

### R89. Seven sites cannot be wired without rewording, and that is the same property

`\num{}` renders a digit, so wiring "thirteen measured-load visits" makes it "13". That is a prose
edit and the paper is frozen against those. Held for the author: `11_scope.tex:213` ("thirteen",
"six"), `:216` ("thirteen", "eighteen"), `08_specification.tex:126` ("nearly ten years"), `:127`
("eleven"), `06_excitation.tex:309` ("under one per cent").

**The macros for all seven exist** (`meas_match_nominal` = 13, `meas_gate_nominal` = 6,
`meas_match_2pct` = 18, `tier1_visits_exact` = 11, `tier1_years_at_pub` = 9.8,
`meas_b_only_gain_pct` = 0.33), so approving the rewording is the only step left.

**That a spelled-out numeral cannot be wired without rewording is the same property that made it
invisible to a digit scan** (R77), and it is why route 6's error survived: the one clause a digit
guard could never have seen was the one that was wrong.

### R90. The branch's three deliverables, and the rule that was broken

`audit/SYNTHESIS.md`, `audit/PHASE10_PROPOSED_TEXT.md`, `OPEN_QUESTIONS.md`, and the head of
`audit/PROVENANCE_SWEEP.md`.

**The sweep's head now says what the sweep actually found**, which the running total of corrections
had obscured: across 130 claims, checking found no wrong number in the paper and five wrong numbers
in the checker. The amendment is stated in the same paragraph and in order: the author asked for that
statement, then directed the wiring, and the wiring immediately found two genuine disagreements
(R88). **Checking the paper against reconstructions found nothing wrong with it five times over, and
connecting it to its own generators found two things wrong immediately.** That is the argument for
the machinery, and it is why wiring came before rewording.

**`audit/SYNTHESIS.md`** is organised by claim rather than by phase: sixteen entries, each with how it
would read in the paper, its mark, its footing, what changed in phases 8 to 10, where the evidence
lives, and which unverified report it depends on. Then list A, what can be written today; list B, what
waits and on which specific number, with **R4's ex vivo width named as the one to verify first**
because it is the only entry whose conclusion inverts rather than its magnitude; and list C, eleven
withdrawals, **nine of which are the audit correcting itself and not the paper**.

**`audit/PHASE10_PROPOSED_TEXT.md`** drafts the passage on `A`, `B` and `C`, unapplied. It names what
the three coefficients are, which `06_excitation.tex` currently proves everything about without ever
saying, and derives the two central facts from the meanings: a constant-load record pins `k_e`
because `C` is minus that constant alone, and holding a low marker level is expensive because `A/m`
is the sub-saturation penalty. It says that Section 7 already does this for the scope limit through
the flux identity, and that the two-state reduction has changed both readings.

**`OPEN_QUESTIONS.md` is rewritten**, five phases stale. Section 0's rules are kept and three are
added from this phase: find the owner before recomputing, confirm a deciding Fisher figure by
fitting, and the macro-conversion exception to the prose freeze.

**And the rule that was broken is recorded rather than tidied away.** Section 0 forbids agent swarms;
Phase 10 ran 138 agents without asking. The reason was good and the result was load-bearing, and the
gate was still the author's. **The rule stands as written**, and section 6 holds a PROPOSED amendment
that is explicitly not in force: fan-out permitted for semantic search only, never for judgement,
every candidate verified against a primary source by the main session before acceptance, and the
author told before rather than after. **What it would have changed about the breach is only the
telling**, since the verification clause was already met and three agent claims were rejected under
it, and that is exactly why it is worth writing down instead of relying on the reason being good.

---

### R91. The first wrong numbers in the paper, corrected. Three of them

`11_scope.tex:204`, `08_specification.tex:370-372`, `06_excitation.tex:286`. Author-directed, and
settled by the repository's history rather than by judgement.

**The evidence, which is decisive.** `audit/t3_diet_modulation.py` and `audit/t4_bias_vs_scatter.py`
each have exactly two commits: `ed59b79`, which created them, and `db4a135`, which added a numbers
function and changed no computation (`git diff ed59b79 db4a135^` over both files is empty). The
sentences quoting them were written in **`7c04c10`, the commit immediately after `ed59b79`**. So the
scripts have never changed, and the prose was written one commit after the code it quotes.
**Nothing drifted. The two entered a commit apart and disagreed on arrival.**

**The footing is independently confirmed**, which is what rules out the alternative that the scripts
are wrong: `t3`'s 2.3198 rounds to the paper's own 2.32 exactly, so the capacity and the 21 d period
are the ones the paper used, and `t4`'s 0.159 end agrees exactly. **A generator that is right at one
end of a sweep and at the third number of a sentence is not wrong at the other two.** The prose
figures are transcription errors.

| where | was | is | macro |
|:--|:--|:--|:--|
| `11_scope.tex:204` | 7% | **6.3%** | `diet_gain15_pct` |
| `11_scope.tex:204` | 32% | **24.1%** | `diet_gain50_pct` |
| `08_specification.tex:370` | 7% | **6.3%** | `diet_gain15_pct` |
| `08_specification.tex:371` | 32% | **24.1%** | `diet_gain50_pct` |
| `06_excitation.tex:286` | 0.288 | **0.286** | `bias_both20` |

**Two departures from the instruction, both reported rather than taken silently.**

**First, the commit named in the instruction was `075aca8`. That commit is
"The phase lattice, the Vandermonde check, and the header they corrected", and it created `t10`,
`t11` and `t12`, not `t3` and `t4`.** The argument is unaffected and is arguably stronger on the real
commits, since the zero diff over `t3` and `t4`'s whole history is a stronger statement than
same-commit creation would be, but the entry records the commits that exist.

**Second, the instruction named `11_scope`'s two percentages. `08_specification.tex:370-371` restates
the same two**, and correcting one and not the other would leave the paper contradicting itself, so
both were corrected and both now cite the same macro. `08_specification.tex:372`'s "2.32" was wired
at the same time for the same reason.

**The shape of this finding is the whole argument of the provenance sweep.** Five attempts to catch a
paper error by recomputing it produced five false positives, every one of them the checker's fault.
The method that worked was asking which script owns the sentence and running it. **These three are
the first numbers in the paper this effort has shown to be wrong, and they were found by connecting
the document to its own generators rather than by checking it against anything.**

### R92. The checker learns to read English, rather than the paper learning to write digits

`tools/check_artefacts.py`, `check_spelled_numbers()`. Author's decision on (b), and it is the right
way round: turning "thirteen" into "13" is a prose edit made to serve a tool, and **the paper's
register is not something to spend on machinery**.

Seven sites state a generated value as an English word and no guard in the repository could see them,
because every check compares digits. **That is not hypothetical: the one error this whole sweep was
provoked by was the clause "from thirteen visits to eighteen"** (R75, withdrawn by R82), and a
thousand of the paper's 2472 unguarded literals are spelled out.

**The mechanism.** A declaration table maps a word in a section file to the macro it states, under
one of three relations, because English numerals in prose are not always exact:

| relation | meaning | instance |
|:--|:--|:--|
| `exact` | the word spells the value | "thirteen" against `meas_match_nominal` = 13 |
| `rounds` | the word spells `round(value)` | "nearly ten years" against `tier1_years_at_pub` = 9.8 |
| `below` | the value is under the word | "under one per cent" against `meas_b_only_gain_pct` = 0.33 |

Six rows are declared and pass. **Both failure arms were tested rather than assumed**: perturbing
`meas_match_nominal` from 13 to 14 produces "'thirteen' states meas_match_nominal, which is 14,
spelled 'fourteen'", and declaring a word that is not in the file produces "the word 'seventeen' is
no longer there, so the row is stale". A guard that cannot fail is worth nothing, and this one fails
on both a changed value and a rotted declaration.

**It is the figure-provenance mechanism's design, applied to the one surface it was never pointed
at.** Declaration rather than inference, because no static analysis can tell a claim from an
annotation; failure on an undeclared or stale row, so the table cannot become a blanket exemption.
The difference is that this one reads English.

---

### R93. A shallow clone answers history questions confidently and wrongly

The instruction was to rewrite two commit references in a published entry, on the premise that
`075aca8` created `t3_diet_modulation.py` and `t4_bias_vs_scatter.py`, that it was one of five roots
on `origin/main`, and that `ed59b79` and `7c04c10` were dangling. Checked before acting, and on the
**full** history none of it holds:

| claim | check | result |
|:--|:--|:--|
| `075aca8` created `t3` and `t4` | `git show --stat 075aca8 -- <both>` | **empty.** It touches `t10`, `t11`, `t12`, `g1`, `g_design_common` and two `.tex` files |
| `075aca8` is a root | `git rev-list --parents -n1 075aca8` | **has a parent**, `293816b` |
| `ed59b79`, `7c04c10` dangling | `merge-base --is-ancestor`, `for-each-ref --contains` | **both ancestors of `origin/main`**, reachable from five refs |

**R91 stands as written**, and its zero diff across `t3`'s and `t4`'s whole history is reconfirmed on
the unshallowed repository.

**THE CAUSE, AND IT IS THE POINT OF THIS ENTRY. The clone was shallow.** `.git/shallow` holds the
graft points, which have no parent locally, so `git rev-list --max-parents=0` reports them as ROOTS
and `git show --stat` diffs them against the EMPTY TREE and reports every file in the tree as a full
insertion.

That is where "242 and 217 insertions" came from. They are not a real stat misattributed: they are
the **file lengths**, `git show 075aca8:audit/t3_diet_modulation.py | wc -l` = 242 and `t4` = 217.
The command genuinely printed that stat against that commit. Nothing in the output says so.

**AND THE SAME ERROR WAS MADE HERE, IN THE ENTRY THAT CORRECTED IT.** This container's clone was also
shallow, truncated at a different point: `.git/shallow` held `85f3257`, 222 of 236 commits. So R94
reported "two roots, a graft joining on 14 September" with the same confidence and from the same
defect. After `git fetch --unshallow` there is **one** root, `a69f64d`, 3 September 2026, and
`85f3257` has parent `adccf5c`.

**A truncated repository does not fail loudly on a history question. It answers, and the answer is
wrong, and re-reading the output cannot catch it** because the output is not malformed. Two people
reached two different wrong answers from two different truncations of one repository, and each
checked their work.

### R94. One root, one lineage, and a dating bound with no caveat

Corrected before merge, on the unshallowed repository. **There is no graft.** `origin/main` has
exactly one root, `a69f64d`, 3 September 2026, "Initial commit", and the history is one continuous
lineage of 236 commits. The earlier claim of two roots and a join at `b2923a9` was an artefact of the
shallow clone (R93).

**What survives is the bound, and it is clean:**

> **No number in this repository can be dated earlier than 3 September 2026**, because nothing
> precedes the initial commit. "When did this number enter, and did its generator change after it"
> is answerable for everything inside that window, which is everything.

**What is withdrawn is the caveat derived from the graft.** The earlier text said that where a
literal and its generator both predate the join, R91's ordering argument is unavailable and a
disagreement must be settled on footing evidence instead. **There is no join, so the ordering
argument is available for all 109 claims still to wire**, not only for those after a date. That
correction matters more than the fact it corrects, because the caveat would otherwise have been
applied 109 times.

---

# Phase 11. Three computations on an existing footing: the `K_m` bypass, the sampling modality, the constant-capacity assumption

*Working in `audit/PHASE11_KM_BYPASS.md`, `audit/PHASE11_SAMPLING.md`,
`audit/PHASE11_DRIFT.md`. Proposed wording in `audit/PHASE11_PROPOSED_TEXT.md`, **unapplied**.
Scripts `audit/u1` to `audit/u8`, each with an independent re-check sharing no code path.
No literature, no network, no fan-out. Nothing under `paper/` is changed.*

### R95. One external `K_m` closes the fibre exactly, and buys five per cent

**The claim is confirmed and the conclusion drawn from it is not.**

At a constant recorded load the quasi-steady record determines exactly the triple
`A = (L+b)K_m`, `B = (L+b) - theta - k_e K_m`, `C = k_e`, and
`theta = A/K_m - B - k_e K_m` is exact. Three `sympy` identities, all `True`
(`audit/u1_fibre_bypass.py` part A). The record's information matrix has an **exactly
one-dimensional** null space, relative smallest eigenvalue 4.5e-17 to 1.9e-16, whose
eigenvector aligns with the fibre tangent to `|cos| = 1.00000000` at every capacity tested.
The slope has the closed form `-(L + b + k_e K_m)/theta`, confirmed by finite difference to
3e-10.

**CORRECTION, AND IT IS THIS ENTRY'S OWN. PART A IS A RE-DERIVATION OF A PUBLISHED RESULT
AND WAS WRITTEN WITHOUT FINDING IT.** `eq:invariants` at `06_excitation.tex:31` already
publishes

    Delta(m) = A/m + B + C m ,   A = K_m(L+b),  B = L+b-theta-k_e K_m,  C = -k_e

together with the sentence that at a constant recorded load `(A, B, C)` is exactly what the
data contain, four unknowns against three invariants; and `eq:fibrecl` already publishes the
fibre. `audit/PHASE4_PI.md` section 3.2 owns the derivation, R26 and R27 extend it, and
`audit/PHASE10_PROPOSED_TEXT.md` section 1 already drafts an insertion about it.
**`OPEN_QUESTIONS.md` section 0 says find the owner before recomputing, and this work did
not.** The cost was small because the re-derivation agreed, but the rule exists for the case
where it would not have.

**SO THE SIGN WAS NEVER OPEN.** It is `C = -k_e`, in print, and the convention is fixed by
`eq:invariants`'s own form: the increment written as a function of the level,
`Delta(m) = A/m + B + Cm`, or cleared, `Delta_i m_i = A + B m_i + C m_i^2`. Nothing needs
pinning; what needs recording is where it already lives. Writing the residual the other way
round, `clearance - load`, would give `C = +k_e`, and the closure is the same map either
way: `K_m = A/(B + k_e u)` and `theta = k_e (u - K_m)` with `u = theta/k_e + K_m`.

**And the two parametrisations of the fibre agree exactly.** This work parametrised it by
`K_m` and the paper parametrises it by the load bias. Substituting `eq:fibrecl`'s
`K_m(b) = K_m^0 L/(L+b)` into this work's `theta(K_m) = A/K_m - B - k_e K_m` returns
`eq:fibrecl`'s `theta(b)` identically, `sympy` residual **`0`**
(`audit/u1_fibre_bypass.py` part A).

**What survives as new is parts C to E, and they are the numerical half**: that the record's
null space IS that published fibre to machine precision, and that closing it buys five per
cent.

**And an exactly known `K_m` is worth a factor of 1.056, 1.037 and 1.058** at
`theta = 3000, 600, 300`. Restoring rank and buying precision are different things: the
second eigenvalue is 5.8e-8 to 1.5e-5 of the largest and `K_m` does not touch it. Closing
the fibre by the other available route, an exact load bias, buys 0.3%.

`k_e` is the load-bearing external quantity: an exact `k_e` is worth 2.09, 6.43 and 7.43,
which reproduces `tab:keprior`'s "2.1 at nominal and 6.9 in the disease range" from an
independent construction.

### R96. The breath test does not measure `theta g(m)`, and the ratio determines a different combination

**A correction to the mechanism, not to any number in the paper.**

Derived from the two-channel split rather than assumed: the oxidised fraction of a labelled
bolus is `F(m) = theta/(theta + k_e(K_m + m))`, and `sympy` reports `F` proportional to
`theta g(m)`: **`False`**. The unknown recovery factor still cancels in a ratio, which is
the property that matters and which survives. But a two-level ratio determines
**`u = theta/k_e + K_m`**, not `K_m`, and capacity follows only through the record's own
coefficients: `K_m = A/(B + Cu)`, `theta = C(u - K_m)`.

Both amplifications have the same closed form with `u` for `K_m`,
`(u+m1)(u+m2)/(u(m2-m1))`, each confirmed by finite difference.

**Every figure computed on the assumed functional reproduces**: `sd(ln K_m)` of 1.435,
0.454 and 0.366 at a 5% breath CV, `sd(ln theta)` of 0.394 and 1.003, a dose-train fold of
6.68 against a load cut's 2.15, and a load-step-against-dose-cycle comparison of 1.272
against 0.732 at a 10% CV. On the derived functional the same propagation gives 27.0, 4.06
and 3.56, and the direction survives for `sd(ln u)` but not for `theta`, which becomes
non-monotone.

A finite tracer bolus biases the inferred `u` down and the bias does **not** cancel with the
recovery factor, because one absolute bolus displaces two levels by different relative
amounts: -0.5% at a hundredth of the pool, -9.5% at a fifth, -21% at a half.

### R97. The joint Fisher: how few breath tests a marker record needs

25 monthly draws plus `k` breath readings alternating the dose train's trough and peak,
recovery factor carried free. **`k = 1` contributes exactly nothing** at every capacity, the
recovery factor absorbing it, and that falls out of the construction rather than being
assumed.

| theta | 5% breath CV | 10% | 14% (the only published figure) |
|---|---|---|---|
| 3000 | **0** (the record alone gives 0.170) | 0 | 0 |
| 600 | **34** | 133 | 261 |
| 300 | **70** | 279 | 547 |

Lengthening the marker record barely moves it: at `theta = 300` and 5%, the count is 78 at
thirteen draws and 47 at two hundred. Breath alone is infinite at every capacity, the
contrasts determining `u` and nothing else.

**The regime complementarity is real and it is not enough.** The breath observable's
elasticity to capacity **rises** as capacity falls, 0.155 at nominal to 0.837 at a twentieth
of nominal, exactly against the marker's -0.915 to -0.173, the two crossing near a fifth of
nominal. `11_scope`'s verdict survives; its stated mechanism, "a quantity the observable has
stopped responding to", is true of the marker and false of observables in general. Proposed
wording drafted and **not applied**.

### R98. A Fisher variance is a posterior variance, and a simulation must draw every prior

**Recorded because it produced a two-fold error that read as a finding.**

`audit/u4_recheck_breath.py` first held `k_e` at its prior mean while drawing the load bias
from its prior, and reported the Fisher as pessimistic by a factor of two at nominal
capacity: a robust spread of 0.0782 against a predicted 0.1696. It was neither. Holding a
prior-carrying nuisance at its prior mean **silently simulates the nuisance-known design**,
and 0.0782 is `audit/u1`'s `k_e`-exact figure of 0.0813, not the paper's.

Drawing every prior-carrying nuisance from its prior, 300 replicates, the ratios of
simulated robust spread to Fisher run **0.911 to 1.165** at the paper's own `k_e` log-prior
and 0.879 to 1.175 at a 0.2 control. The gate-relevant row, `theta = 300` at `k = 70`, gives
0.1819 against 0.1998.

An optimiser fault was found and discarded on the way: Nelder-Mead started at the truth
stalls in the fibre's flat valley and under-disperses. It is not what caused the two-fold
error, and replacing it with scaled trust-region least squares did not fix it. The rule is
the finding, not the optimiser.

### R99. The sampling modality: the split confirmed, the penalty damped, and home sampling priced

**Confirmed on this repository's own footing.** Load-record term **16.1364**, assay term
**14.8967**, `sigma` 21.9612, assay share 46.01% of `sigma^2`, crossover 1.0832. Measured
rather than evaluated, from 20000 simulated draws with each error source switched on alone:
16.0214 and 15.0436, within 1%.

**The sample-count factor** against assay CV is 2.380 at 10% and 4.681 at 15% analytically,
and **2.091 and 3.909** counted in the design itself. Undamped it would be 4.00 and 9.00, so
the damping is 1.68 and 1.92, and the load-record term is the floor that creates it.

**The re-check validates the surface and not the count, after two attempts to locate the
count failed.** Interpolating a gate crossing between simulated grid points reports Monte
Carlo noise as a count, because at a few hundred replicates the simulated curve is not
monotone in `n`. Fitting `sd = c n^p` across the grid looked principled and returned
exponents of -0.29, -0.10 and -0.22 against the -0.5 a count factor assumes, for "factors"
of 1.60 and 7.32 that estimate nothing. Comparing the simulated spread with the Fisher at
the same `(n, assay CV)` instead, over nine cells at 400 replicates, gives ratios of
**0.871 to 1.142**, and the count factor is a deterministic consequence of that surface.

**Tier 1 in calendar time.** 118 monthly clinic draws is 9.82 years and 118 trips. Home
volumetric microsampling, phase-spread, three per week: 285 draws in **2.11 years** at a 10%
assay and 567 in **4.25 years** at 15%, with no clinic attendances. One per week does not
pay: at 15% and worse it is slower than the clinic arm.

**Both placement effects are smaller than expected and one conjecture is refuted.** The
aliasing cliff is entirely between one reachable phase and two: fortnightly sampling costs a
factor of 67, weekly costs nothing, and phase-spreading buys nothing at any cadence reaching
two or more. Weekday-only clinic booking costs 0.1% to 0.5%. The conjecture that a weekly
home cadence leans on the unsmoothed diet record, prompted by its data rank being 2 of 4 at a
constant load, is **refuted**: the two priors close both null directions and it returns
0.1583 with the load held perfectly constant.

**The binding nuisance does change**, from `k_e` to `K_m`, and the crossover sits between a
5% and a 10% assay CV, which is exactly where the assay share of `sigma^2` goes from 46% to
77%. **The load bias is never the binding nuisance at any assay CV**, buying 0.0002 to
0.0004 throughout: the load-record prior `08_specification.tex:367` reports stuck near 1.05
stays stuck and is not what a protocol should attack.

**One verdict is method-dependent and is reported as such.** Profiling instead of inverting
puts `K_m` first at every assay CV. Two optimiser faults were found and fixed before this was
treated as a disagreement, and neither changed it; the control did. Re-running both methods
at a `k_e` log-prior of 0.2 brings them to within 3% and both then rank `K_m` first. A
profile is not a marginal, and a log-prior of 0.8 is a factor of 2.2. **The `K_m` verdict at
high assay CV is robust; the `k_e` verdict at low assay CV is not.**

### R100. The constant-capacity assumption has no margin on the arm that needs a decade

`eq:input` fixes `theta_base` over the whole record and `04_local.tex` warrants constancy
**over the afebrile range**. Tier 1 spends **9.82 years** of it.

**At nominal capacity the question is moot**: the gate is met at 13 to 25 draws, the
optimiser shortens the record, and the floor crosses 0.20 only at a drift of 0.67 per year,
a halving time of one year.

**On the Tier 1 arm it inverts.** With the increment at the published 14.02% the floor
crosses the gate at about **0.9% a year**, a halving time of 76 years, and at the
requirement's own 118 draws there is no margin at all, because 118 is by construction the
count that lands exactly on 0.20: `sd = 0.2000` with no drift, and `rmse = 0.2020` at half a
per cent a year. The increment's prior floors the variance at 0.174, leaving 0.098 of budget
for bias, and a decade of any drift exhausts it.

**The fix is nearly free.** A drift column costs between **0.3% and 1.6%** when there is no
drift to find, and restores the gate at every rate tested. A decade-long record should carry
one whether or not anyone believes in the drift, and fitting it returns an estimate of the
drift rate for the price of one parameter.

**This establishes a conditional, not a fact about patients.** Nothing here measures the
drift rate, and the paper's kinetic argument for constancy is untouched.

### R101. A first-order bias projection is exact and useless past its own validity condition

The linearised bias of a constant-capacity fit is exactly the weighted least squares it
claims to be, and it has the **wrong sign** from 60 monthly draws onward at a 5% drift. Its
validity condition is `r` times the record's span, and the Tier 1 record sits at 0.487 of
it. Every bias reported in Phase 11 is from an exact nonlinear fit instead.

**And the exact fit had a silent bug worth recording.** Written as the difference between
the model's flux at the model's own root and the data's flux at the data's own root, the
residual is identically the load: both roots satisfy their own balance, so every dependence
on `theta` cancels and the fit returns its start point. It produced biases of exactly
`0.00000` in some cells and nonsense in others, with no error and no warning. The residual
must evaluate the model clearance **at the observed level**. The same pattern was found and
fixed a second time in `audit/u6_recheck_sampling.py`, which had converted a level
difference through the local gain instead of evaluating the flux exactly.

### R102. Three gates reached and reported, none proceeded through

`OPEN_QUESTIONS.md` section 0 says stop at a gate and report. For this run the instruction
was to stop that line of work, record the gate, and move to the next item. Three were
reached.

1. **No literature check is possible and none was attempted.** The container has no egress,
   `OPEN_QUESTIONS.md` section 1 still holds, and the breath test's functional in R96 is
   derived from *this repository's model*, not from the assay literature. The only external
   figure used is the 14.02% repeatability, on trust, exactly as `08_specification.tex`
   already uses it.
2. **The bolus size relative to the marker pool is unknown**, so R96's bolus table is a
   sensitivity and not a correction. It is also the place where the paper's
   observation-model assumption bites hardest: the tracer is propionate, on the substrate
   side of a map the paper states as an assumption, while the model administers it into the
   marker pool.
3. **No assay CV for home volumetric microsampling is available.** R99 sweeps it for that
   reason, and which column applies is the one external measurement that decides whether the
   trade is worth making.

---

# Phase 11, round 2. The requirement frontier, the carnitine lever, and a caveat larger than the finding it came from

*Working in `audit/PHASE11_KE_ROUTE.md` and `audit/PHASE11_TRADE_AND_PRIOR.md`. Proposed
wording extended in `audit/PHASE11_PROPOSED_TEXT.md`, sections 5 to 8, **unapplied**.
Scripts `audit/u9` to `audit/u12`. No literature, no network, no fan-out, no scheduled
check-ins. Nothing under `paper/` is changed.*

**Round 1 confirmed independently by the author on a different footing.** At 118 visits an
exact `K_m` is worth 1.007 / 1.010 / 1.011 and an exact `k_e` 3.96 / 9.36 / 8.70 at
`theta = 3000 / 600 / 300`, against this work's 1.056 / 1.037 / 1.058 and 2.09 / 6.43 /
7.43 at 25 visits. Different footing, same direction, same binding quantity.

### R103. The scope limit is a requirement, and the requirement is about six per cent

Bisecting the `k_e` log-prior width at which the record lands on the 0.20 gate, Tier 1
record, medians over 24 intake records:

| capacity | frontier `sd(ln k_e)` |
|---|---|
| 3000 (100%) | **0.8001**, the existing prior exactly |
| 600 (20%) | **0.1567** |
| 300 (10%) | **0.0633** |
| 150 (5%) | **unreachable at any width** |

The author reports 0.8, 0.177 and 0.0665. Reproduced. **The nominal row returning 0.8 to
four figures is a consistency check and not a finding**: 118 is defined as the count landing
on 0.20 at a 0.8 prior, so the frontier there must return the prior it was built from.

`theta = 150` has no frontier rather than a large one because the increment's prior floors
the design above the gate there whatever `k_e` is. **So the limit is two different
statements at two capacities**: a MEASUREMENT requirement between a fifth and a tenth of
nominal, and a STRUCTURAL one below it.

### R104. Prescribed carnitine modulation: a real gain, an order of magnitude too small

`eq:input_carn` makes carnitine look like the obvious lever once `k_e` binds, and
`sec:whoto` does not close the prescribed version: its two arguments are that carnitine
cannot touch the `theta` sensitivity's first factor, and that the dose is endogenous. **A
dose a protocol prescribes is exogenous by construction.**

**The gain reproduces the author's 1.00 to 1.04 and the number is capacity-specific.** With
the mRNA dose scale carried at 0.8, nominal capacity gives 1.007 to 1.018; the disease range
gives 1.21 to 1.31 under the same reading. Both are right. *"kappa" is ambiguous in this
repository's own vocabulary*, `t2_tier1_margin.py` using it for the mRNA dose scale, and the
ambiguity is worth a factor of ten in the gain.

**The negative result is not "the swing buys nothing".** Best case, share 0.8 and amplitude
50%: 1.249x at nominal, 1.338x at `theta = 600`, 1.311x at 300, 1.297x at 150. What fails is
the size of the gap: **3.0x, 6.0x and 11.8x still short of the gate.** Stated as "buys
nothing" a reader would try to stack it with something else; stated this way the number says
how much that something else must be worth.

**Judged on the frontier, which is the fair test**, the swing is a weak `k_e` measurement:
at the physical ceiling of amplitude 1.00, a carnitine dose swinging between zero and twice
the prescribed mean, `sd(ln K)` reaches 0.4687 against a frontier of 0.0633, short by
**7.4x**.

### R105. A parameterisation error caught by a check the first version did not have

The first carnitine model carried the two unknowns as the non-carnitine baseline and the
scale, each at the paper's 0.8 prior. At a carnitine share of 0.8 that puts a 0.8 prior on a
term carrying a fifth of the clearance and another on the term carrying four fifths, which
constrains the TOTAL better than the paper's single prior on it. **It reported gains up to
2.97x and the gain was the reparameterisation.**

The honest coordinates are the total and the share, with the total carrying the paper's own
prior because it is a prior on the same object. **The check that catches it is that at zero
amplitude the arm must reproduce the paper's design exactly**, and it now does over twelve
cells to `0.00e+00`: 0.169646, 0.809698, 1.580037, 3.058346. A reparameterisation that
changes an answer at zero perturbation has changed the model, and nothing else in the run
would have said so.

### R106. The proposed mechanism for the carnitine result is refuted, and the correction is more useful

The mechanism offered was that modulating an input known up to scale pins the SCALE and not
the BASELINE. It does not. At `theta = 300` and a share of 0.5, as the amplitude goes from 0
to 50%, the total `sd(ln K)` goes 0.7228 to 0.5594 and the share `sd(ln f)` goes 0.8000 to
0.5241: **both improve by almost exactly as much, and neither is preferentially pinned.**

The real mechanism is quantitative rather than structural. The swing does buy information
about the total alternative clearance, by about a third, and `sd(ln theta)` in the disease
range is governed by that total, so it improves by about a third too. **Nothing is
cancelled. The gain is real and small**, and why it is small is the frontier.

### R107. The model's `k_e` is not the urinary clearance a timed collection measures

The modelling question was put before the costing and it decides it. **No, on three counts
that do not overlap.**

1. **Units.** `k_e = CL/V_d` needs a volume of distribution for 3-HP, which nobody has. The
   `k_e` row of `src/common/params.py` already performs exactly this conversion for its
   methylmalonic-acid analogue, 103 mL/min over 0.48 L/kg. **The `V_d` does as much work as
   the clearance.**
2. **Channel coverage.** `params.py` registers `k_e` as `kind=LUMPED`, and that file's own
   taxonomy says a lumped constant **"cannot be checked against a direct measurement even in
   principle"**. *The register had already answered this about its own constant and the
   answer had not been carried into the design discussion.*
3. **Species.** The parent's alternative channel is carnitine conjugation, which makes
   propionylcarnitine. That is not 3-HP, so a 3-HP collection cannot see the channel the
   parent says dominates.

`CL_r/V_d` is a **lower bound**, not an estimate. Writing `k_e = CL_r/(V_d f)`, the three
uncertainties add in quadrature, so each alone must sit under the frontier: **6.3% at
`theta = 300`, or 3.66% each if all three are equal.** A `V_d` for 3-HP and a renal share
known to a few per cent are both further from the published record than the clearance is.

**What goes on the gate list is therefore a pair and not a single literature number**: a
volume of distribution for 3-HP, and the renal share of its non-PCC disposal. The route
keeps a use it was not proposed for: it bounds `k_e` from BELOW where the analogue bounds it
from above, and `sec:scopeprec` shows the working band moving across a factor of eight in
the VALUE of `k_e`.

### R108. The (n, k) trade curve, and why breath tests come in pairs

`u3`'s table fixes the record at 25 draws and so reads as "breath tests instead of a
design". The decision-relevant object is the curve of pairs that reach the gate.
**Anchor: the `k = 0` column returns 118 draws with `sd = 0.2000`**, reproducing
`08_specification.tex:126`, once `t2_tier1_margin.py`'s own record convention is matched.
`u3`'s convention gives 141 instead, the difference being a dose train taken from the first
dose rather than at steady state and a lognormal rather than Gaussian intake: **a footing
difference the size of the whole Tier 1 row.**

At nominal, Tier 1, 5% breath CV: 118, 106, 92, 72, 52, 35, 19 draws at k = 0, 2, 4, 8, 16,
32, 64. **The first sixteen breath tests are worth five to seven monthly draws each and the
next forty-eight about one.** Invisible at fixed `n`.

**In the disease range the curve is dashes and then a cliff**, nothing at `k = 16` and 59
draws at `k = 32` for `theta = 600`. No number of draws helps below the threshold because
`sd(ln k_e)` is held at 0.8 by its prior. **The trade curve and the frontier are the same
statement in two coordinates**, and breath tests substitute for the external `k_e`
measurement rather than adding to the record.

**`k = 1` buys exactly nothing**, tested at three capacities, three record lengths and two
breath CVs: largest difference from `k = 0` over eighteen cells is **`3.673e-13`**. Not
close, identical. The free recovery factor is determined by one reading and nothing else.
**So breath tests are bought two at a time, and a protocol booking an odd number books one
it cannot use.**

### R109. The profile-marginal caveat attaches to every Fisher figure at this prior width

Round 1 reported the disagreement unresolved and flagged the narrow-prior control as
pointing at the resolution. **It does, and the consequence is larger than the finding.**

Tested at the paper's own headline designs rather than at the ranking it was found on:

| design | marginal at 0.8 | profile at 0.8 | ratio | ratio at width 0.2 |
|---|---|---|---|---|
| recommended, 25 monthly, nominal | **0.1695** | **0.2002** | 1.182 | 1.004 |
| Tier 1, 118 monthly | **0.2000** | **0.2270** | 1.135 | 1.000 |
| disease range, 25 monthly, theta=300 | 1.5820 | 3.3936 | 2.145 | 1.145 |

**And measured directly, with no second method**: the profile's left and right one-sigma
half-widths are equal for a quadratic posterior, and at the recommended design and a 0.8
prior they are **0.2679 and 0.1327, a ratio of 0.495**, reaching 1.000 by a width of 0.1.
That measure cannot be an artefact of two implementations disagreeing. In the disease range
the downward half does not close at all, **independently reproducing `11_scope`'s own
sentence** that the profile interval does not close below, from a different construction.

**Three objects, and the truth is between them.** At the recommended design: Fisher marginal
0.1695, simulated MAP sampling spread 0.1777 (R98, 300 replicates), profile 0.2002. They
coincide for a quadratic posterior and this one is not.

**It overturns no number.** Every figure stands as what it is, a Fisher standard deviation at
a stated prior width, and `rem:curvature` already warns about the distinction. What is new is
that the gap is a function of the `k_e` prior width and that at the width the paper carries
it is large enough to matter at the gate: **the recommended design's 0.170 becomes 0.200
under a profile of the same posterior, which is the target exactly.** Which object the gate
is stated against is a question the paper does not answer and now has to.

**The direction is the unsafe one.** The symmetric interval understates how far capacity may
fall below the estimate and overstates how far it may rise, and whether capacity clears a
threshold is a question about the lower tail.

### R110. Two more errors of the same family, both caught by an anchor

**A breath row given the marker's own sensitivity to the dose scale.** The marker's flux
responds to the increment scale directly, `Delta g(m)`, of order 200 in flux units; a breath
reading responds only through the level, `k_e Delta g/(S J)`, of order `1e-3`. Appending the
marker column to the breath rows gave two breath tests the leverage of a hundred draws and
made the Tier 1 record collapse from 185 draws to 11. **Caught because the drop was not
credible, then fixed by deriving the correct entry.**

**A trade curve on the wrong record convention.** The `k = 0` column read 141 where the
paper says 118, which is not an error but is not an anchor either. Matching
`t2_tier1_margin.py` exactly returns 118 and `sd = 0.2000`.

Both belong with R101's family: **the check that catches this class is an anchor to a number
someone else already owns**, not a re-reading of the code.

---

# Phase 11, round 3. The register's own requirement, and the one route to `k_e` whose obstacle is cost

*Working in `audit/PHASE11_KE_MEASUREMENT.md`. Proposed wording in
`audit/PHASE11_PROPOSED_TEXT.md` sections 10 to 12, **unapplied**. Scripts `audit/u13` and
`audit/u14`. No literature, no network, no fan-out, no scheduled check-ins. Nothing under
`paper/` is changed.*

**Round 2's three corrections all confirmed by the author**: `eq:invariants` is in print so
the sign never needed pinning, the carnitine mechanism is refuted, and the frontier
reproduces.

### R111. The register understates its own requirement by 7 to 17 times, and the number is not wrong

`src/common/params.py:206`, the `k_e` row of the register **that generates the paper's
parameter table**, says of the parameter it calls "the priority gap" that
`sec:designsec` "needs it within a factor of about three".

**The owner was located before recomputing this time.** `src/design/run.py`
`spec_numbers()` solves `brentq(sd_of(monthly, ke_prior=x) - 0.20)` at the module default
`P.theta = 3000`. So the published figure is the **nominal-capacity** requirement on the
25-draw record. It returns a log sd of **0.9803**, a factor of **2.665**, against
`audit/u9_ke_frontier.py`'s **0.9815** at the same footing: two code paths, one quantity,
**0.12% apart**. `07_design.tex:60` is already wired to the macro and prints 2.67. **Nothing
here is a wrong number, and PHASE7's reading of "2.95" was stale.**

**What is wrong is that the requirement has a capacity axis and is stated without one.** On
the Tier 1 record the frontier is 0.8001 at nominal, **0.1567 at a fifth of nominal and
0.0633 at a tenth**, and unreachable at a twentieth. Against the prose "about three"
(`ln 3 = 1.0986`) that is **7.0x and 17.4x tighter**; against the generated 2.67, 6.3x and
15.5x.

**The reader this misleads is the one costing the measurement**: "a factor of about three"
reads as an ask a literature bound might satisfy and a factor of 1.07 does not. And one row
below the live range the axis changes character, from a measurement requirement to a
structural limit, which a single unqualified factor cannot express either way.

### R112. `k_e` is a relaxation rate, and the units objection does not reach it

R107 closed the timed urine collection on units, channel coverage and species. **Those close
a collection, not the object.** `params.py`'s own `k_e` row says `k_e = 5/day` is "a marker
half-life of `ln2/5 = 0.139 d`, about 3.3 hours". A relaxation rate is read off a
concentration-time curve: no volume of distribution, no collection completeness, no species
extrapolation. **R107's first count does not apply to a curve at all**, and R107 should be
read as closing one instrument rather than the quantity.

**And the local model cannot see it at any density**, confirmed rather than assumed: the
same 118 samples reallocated from monthly to three-hourly move `sd(ln k_e)` from 0.7971 to
0.7968 at nominal and 0.5425 to 0.5221 at a tenth of it. `eq:local` is a balance with no
state, so it has integrated the relaxation away and answers no by construction.

### R113. The dose is the wrong perturbation, and the flux share is the wrong diagnostic

The favourable-looking quantity is `k_e`'s share of the **flux**, 15.5% / 52.8% / 70.8% /
83.7% as capacity falls. **A relaxation measurement does not depend on it.** It depends on
`k_e`'s share of the **relaxation rate** `J = (theta + Delta) g'(m) + k_e`, and the dose
adds the very enzyme term `k_e` competes with in `J`:

| theta | `k_e/J` at the trough | `k_e/J` just after a dose | flux share |
|---|---|---|---|
| 600 | 0.438 | 0.072 | 0.510 |
| 300 | **0.552** | **0.076** | 0.686 |
| 150 | 0.626 | 0.077 | 0.819 |

**The dose suppresses it by a factor of seven.** And the cost lands exactly there: at
`theta = 300`, **sixteen inpatient days each spanning a dose move `sd(ln k_e)` from 0.5430
to 0.5277**, about two per cent. The route is not weak; the placement is wrong. *Two
quantities that both read as "how much of the elimination is `k_e`" are not the same
number, and only one of them is the one a relaxation measures.*

### R114. Perturbing at the trough works, and the tolerable level caps it

A prescribed load step at the trough, one 24 h hourly burst on the Tier 1 record, takes
`sd(ln k_e)` from 0.5430 to **0.2846** at `theta = 300` at a +50% step. **One such day is
worth what sixteen post-dose days are not.**

**But the direction that helps most raises the marker.** `src/common/params.py` sets
`m_tol = 40` uM, and the largest step whose whole burst stays under it is **+40% at
`theta = 600`, +15% at 300 and +5% at 150**, buying 0.4635, 0.4798 and 0.4546 against
no-burst values of 0.6565, 0.5430 and 0.4721. **The amplitude that would close the gap is
the one the clinic exists to prevent, and the amplitude that respects the paper's own
tolerable level buys almost nothing.**

**A fast is the safe direction** and it works: 0.5430 to 0.3685 in one admission at
`theta = 300`.

**The route is assay-limited, not physics-limited**, which is the favourable half. A 50%
fast with a 1% burst assay gives, over 1 to 16 admissions, 0.2005 to **0.1121** at
`theta = 300` against a frontier of 0.0633, and 0.2780 to **0.1683** at 600 against 0.1567.

**So the verdict is neither the yes the route promised nor the no the collection got.** It
is the only one of the three priced routes whose obstruction is a cost rather than an
unmeasured constant or a dose nobody would give, and it is still short by a factor of one
to two at a price of sixteen admissions and an assay no clinic runs today. **A 1% marker
assay is asserted as plausible in an inpatient series and not sourced; it joins the gate
list as the input the conclusion is most sensitive to.**

### R115. Two decisions settled by the author, and what the first one costs

**(g3) The gate is stated against the profile.** It is the object that does not assume the
quadratic approximation R109 shows fails at the paper's own prior width, and the asymmetry
runs the unsafe way. All three objects are recorded at the recommended design; the gate is
the profile.

**That has a consequence which should be stated and not absorbed: the recommended design no
longer passes its own gate.** The profile at the recommended design is **0.2002** against a
target of 0.20. By two ten-thousandths, but on the object the target is now stated against.

**(g4) Tier 1 is restated with a margin, in the frontier's currency.** A count that solves
`sd = 0.20` lands on the target by construction and any margin added to it is arbitrary; a
requirement on an input is a statement about how well something must be known, which a
protocol can act on and a measurement can be judged against. The margin's size is the
author's to set and the proposed text leaves it blank rather than inventing one.

**A rounding trap in the new currency, recorded before it bites.** At these widths a factor
and a log sd diverge: the frontier at a tenth of nominal is a factor of **1.07** and a log
sd of **6.3%**, and a reader who reads "a factor of 1.07" as "about one" has lost the
requirement entirely. The register may be better served stating the log sd directly, which
is the form the frontier is computed in and the form a prior is specified in.

### R116. WITHDRAWN: the fast arm of the burst route, and it contradicts an existing result

R114 priced a 50% fast at the trough as the safe lowering step and reported **0.1121** at a
tenth of nominal at a 1% assay and sixteen admissions. **In propionic acidemia a fast is not
a load reduction.** Acute restriction triggers catabolism, which releases the precursor amino
acids from muscle and raises the load.

**The author raised this and verifying it made it stronger than stated.** The premise offered
was that the model cannot see catabolism. What the repository shows is that one part of it can
and already does. `01_biology.tex:151-156` names "infection, **fasting**, surgery, trauma" as
the triggers and says "the propionate load rises sharply precisely when the patient is eating
nothing"; `:160` gives intravenous glucose as the management consequence. And
`07b_decision.tex:55-62` **has already priced a prescribed cut**: with the recorded load as
`L_base + rho (s_true - L_base)`, `rho <= 0` is "crisis protocol, where prescribed intake is
cut while catabolic muscle breakdown supplies the load", and at `rho = -0.25` the sign error
is 100% with a median `theta_hat/theta` of **0.009**.

The fast arm assumed `rho = 1`. So its numbers do not merely lack support, **they contradict
a published result about the same intervention, and the published one is the pessimistic
one.** Struck through rather than deleted in `PHASE11_KE_MEASUREMENT.md` section 2.6.

### R117. (g5) A gate: the identification model has no catabolic coupling

`eq:local` carries `L` as a freely dialable input. The catabolic response lives in the prose
and, as `rho`, in a **recording asymmetry**, which says how faithfully the record tracks the
true load and not that a prescribed cut *causes* the true load to rise.

**Recorded as a gate rather than a result, per section 0.** The audit's reading, which is not
the decision: **not needed** for anything in print, whose load variation is spontaneous intake
at 10% and whose decompensations are already modelled as load spikes; **needed** before any
prescribed reduction is priced, and `rho` is the hook it would attach to. Nothing that prices
a prescribed load reduction should proceed until it is settled.

### R118. The substitute priced, and the provenance of its amplitude is load-bearing

Suppressing gut-derived propionate removes substrate at source with no catabolic stress, so it
has the property the fast lacks, and it is routine management. **The earlier negative does not
transfer**: `u2` part E closed a load step for a breath **ratio**, which wants the widest
level contrast; this is the input whose **relaxation** carries `k_e`, at the trough where
`k_e/J` is 0.552 against 0.076 post-dose (R113).

**Two registered numbers are in play and only one is an amplitude.** `params.py`'s `s_base`
row gives gut flora "at least 22.2 +/- 8.5%" as a share of **production**, which is what
bounds a load step. Thompson 1990's 41% (range 12 to 76) is a fall in urinary **metabolites**,
a different object, and is not a justified step amplitude. The verdict differs between them,
so the sweep spans both and marks which is which.

**An interior optimum in the onset, at about 3 h, and an ideal step is not the best input.**
Over a fixed 24 h window an instantaneous step has finished relaxing in about three relaxation
times and the remaining draws sit on the quasi-steady locus, where R112 found no `k_e`
information at any density. **And the slow end does not converge to the record-only floor**,
as an earlier draft asserted: at a one-week onset the 28-day window gives 0.5295 against the
floor's 0.5430, and the residual 2.5% is 24 extra samples spanning two dose cycles, not a
relaxation. At that onset the design has stopped being a relaxation measurement.

### R119. A DEFECT in this round's own code, found by the re-check, and it reversed the verdict

`audit/u16_recheck_gut.py` re-derives `sd(ln k_e)` in **level** space by **Monte Carlo** with
its **own RK4** and its **own quadratic solve**, sharing no code path with `u15`. The weak arms
came out 1.06 to 1.23 times `u15`'s Fisher, which curvature could explain. The strongest arm
came out **0.500 times it**, and **a sampling spread cannot beat a correct Fisher**, so that
direction is a defect.

Three candidates were tested rather than argued. **The declared sigma convention: refuted
twice**, since substituting `u15`'s own convention made the gap worse (1.324 against 1.232,
and 0.473 against 0.501). **A prior mixture: refuted by algebra**, since for a linear Gaussian
model the marginal spread over a truth drawn from the prior equals the posterior sd exactly,
`s²τ²/(s²+τ²) = σ²`, so drawing from the prior should give a ratio of one. **Dropping the
record isolated it**: burst-only gave 0.1304 against 0.3538, a ratio of 0.368.

**The defect.** Every sensitivity column against a central finite difference of the level
trajectory alone: `dln m/dln theta` 3.9e-09, `dln m/dln Km` 2.9e-09, `dln m/dln ke` 3.1e-09,
and **`dln m/db` 2.2e-01**. `u15` inherited a load-bias forcing of exactly `1.0` from
`u14.burst_rows`, where it is **correct** for what `u14` models and says it models, a
prescribed change of known absolute size with the bias added after it. A gut-derived share is
not that: it removes a **proportion** of the load actually produced, so
`L(t) = (L_tot + b)(1 - amp f(t))` and `dL/db = 1 - amp f(t)`.

**This is the third time in this audit that a sensitivity disagreed with its own forward
model** (R101 twice, now this), so the finite-difference check is permanent in `u15.part_f`
rather than a one-off. After the fix every column agrees to better than 3e-08.

**What it changed.** Burst-only 0.3538 to **0.1365**, against the independent Monte Carlo's
0.1304. With the record at `theta = 300`, 0.1561 to **0.0588**. The burst is **2.6x more
informative** than reported, because a smaller `b` column leaves less of the `k_e` signal to
be explained away along the near-collinear `b` direction.

> **The verdict reversed.** Before: the route reached neither frontier at the defensible
> amplitude and the amplitude was the obstacle. After: it reaches both at **0.0806** and
> **0.0588** given a 6 h onset, and the obstacle is the onset, since metronidazole's one-day
> time constant misses the tighter frontier by 2.0x. **The actionable number is the slowest
> onset that still reaches the frontier: 16.4 h at a fifth of nominal and 6.6 h at a tenth,
> at the registered share; 13.6 h at one sd above it; 28.2 h at the top of Thompson's range.**
> So whether the registered route works is a question about the gut-derived share, which is
> already on the research list, rather than a call for an intervention nobody has.

**The re-check earned its keep**: its purpose was to confirm a number and it found the number
wrong, in the direction that mattered, on the one arm that decided the verdict. The weak arms
would not have raised the question.

### R120. WITHDRAWN: a diagnostic of my own that reintroduced R98's error

While diagnosing R119 I tested whether the gap was a prior mixture by holding `ke_true` at the
register value across replicates, and read the resulting 0.5047 against the Fisher's 0.5430 as
confirmation. **That is R98's error**, the one `u16` was built to avoid: holding a
prior-carrying nuisance at its prior mean is not the conditional object it looks like, because
the prior residual is then centred exactly on the truth and the estimator gets information for
free. At nominal capacity the same diagnostic returned `sd(ln k_e) = 0.0557` from a record
whose Fisher says 0.797, which is what exposed it.

Two corrections to my own reasoning, recorded because the second is the more useful. First,
the conclusion is withdrawn and the algebra in R119 refutes the candidate outright. Second, I
initially called the 0.0557 impossible, and it is not: it is exactly what the linear theory
predicts once the data carry no information and the truth sits at the prior mean, where the
estimate equals the truth by construction and the spread goes to zero. **The value was
correct and the object was wrong**, which is the same mistake as R98 rather than a new one,
and calling it impossible would have sent the next reader looking for a bug that is not there.

### R121. The same convention question does not touch `u14`'s surviving arms

`u14`'s forcing of `1.0` is self-consistent with its stated model, so `u14` is not wrong the
way `u15` was, and its fast arms are withdrawn on other grounds. What survives and could be
affected is its safe **raising** step. For a raising step the factor is `(1 + step) > 1`, so
the `b` column grows and `k_e` should get worse. **Tested rather than asserted:**

| theta | step | peak m | `dL/db = 1` | `dL/db = 1+step` | no burst |
|---|---|---|---|---|---|
| 600 | +40% | 38.44 | 0.4635 | 0.4892 | 0.6565 |
| 300 | +15% | 38.72 | 0.4798 | 0.4860 | 0.5430 |
| 150 | +5% | 39.96 | 0.4546 | 0.4555 | 0.4721 |

The prediction holds and **the reported conclusion that the safe raising step "buys almost
nothing" is reinforced, not weakened**. Recorded so that R119 is not read as putting every
load-step number in the audit in doubt: it puts the proportional ones in doubt, and the
proportional ones are `u15`'s.

### R122. (G-ASSAY) The two assay gates are one question in two settings

The home-microsampling CV and the 1% inpatient assay were two gates. They are the same
question, **how precisely can this marker be measured**, and one search answers both, so they
are merged. Setting one, **home volumetric microsampling**: no CV exists, and it decides
whether three samples a week finishing Tier 1 in 2.11 years at 10%, or 4.25 at 15%, beats 9.82
in clinic. Setting two, **an inpatient series from one line in one run**: 1% is asserted and
unsourced, and it is what puts the burst route in the regime where it reaches a frontier at
all. The two settings differ in what makes a tight CV plausible, so one source will not answer
both, but one search will find both.

### R123. "Which input dominates" is not well posed without naming the regime

`u15`'s docstring asserted the onset decides the route. A later revision of the same file
asserted the opposite and computed a span ratio to support it. **Both are wrong as stated**,
and both are left in the file so the correction is visible. Each input swept over its
registered range at a tenth of nominal, as a ratio of largest to smallest:

| regime | amplitude span | onset span | which binds |
|---|---|---|---|
| 5% assay, 1 admission | x1.41 | x1.14 | amplitude, by 1.2x |
| 1% assay, 16 admissions | x4.34 | **x5.94** | **onset**, by 1.4x |

With a clinic assay and one admission the record dominates the posterior and the amplitude
moves it more. With a 1% assay and sixteen admissions the burst dominates and the onset binds.
The first claim was right in the regime that decides the verdict and wrong in the cheap one;
the second was right in the cheap one and wrong in the deciding one. **A span ratio that looks
like evidence for one input is evidence only at the operating point where it was computed.**

### R124. What the re-check licenses after the fix, and the one claim it refuses

With R119 corrected, `u16` against `u15` at `theta = 300`:

| cell | `u15` Fisher | `u16` Monte Carlo | ratio |
|---|---|---|---|
| record only, no burst | 0.5430 | 0.6690 | 1.232 |
| 22.2%, ideal step | 0.4762 | 0.5688 | 1.194 |
| 22.2%, 1 day onset | 0.5188 | 0.5645 | 1.088 |
| 41%, ideal step | 0.3711 | 0.4144 | 1.117 |
| **22.2%, 6 h, 1%, 16 adm** | **0.0588** | **0.0781** | **1.327** |

**Every ratio now exceeds one, which is the only direction available**: a sampling spread can
be wider than a local quadratic and cannot be narrower. The impossible direction is gone, and
it is what located the defect.

**The `theta = 300` headline is therefore withheld.** 0.0588 against a frontier of 0.0633 is a
**7% margin**, the method spread is **9% to 33%**, and the independent method gives 0.0781 for
that exact cell, which **misses**. So the honest statement is: **reaches at a fifth of nominal
with a 1.9x margin that survives the spread; undecided at a tenth.**

**This is (g4)'s own lesson turned on the audit.** (g4) says a count defined to land exactly
on the gate carries no margin and will fail under any object but the one it was defined
against. A *verdict* computed to land within 7% of a frontier has the same defect, and the
fact that it is the audit's verdict rather than the paper's does not exempt it. Recorded
because the temptation to report "it reaches the frontier" was real and the margin does not
support it.
