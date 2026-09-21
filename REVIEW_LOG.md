# Review log

Where the mathematics of Section 3 stands, so that the next session does not
have to re-derive it. Every model number below is produced by `src/` and asserted
in `tests/test_models.py`; none is typed into the paper by hand. Numbers taken
from or derived off the literature are not model numbers and do not live here:
they live in `paper/claims.csv`, against the source and with the verdict, and
the exact Poisson intervals in `894ae04` are recorded there as computed rather
than quoted.

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
qualification that the reduced separatrix is a lower bound on the basin rather than the
basin. `alt_flux_C`'s docstring called this an adiabatic elimination of a fast state and
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

*Section 5.3 said the figure draws three numbers it does not draw.* The sentence gave the
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
only here. Section 3.2 said `prop:rank`'s null vector is the fibre's tangent at `b = 0`
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

Section 4.3 now gives the first its `n`, its span and its lap count, and says in one
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
dosing guidance: Section 5.6 already says a protocol written from this work is not a rule
for escalating a patient, and a forward simulation at an **assumed** capacity is weaker
evidence than that, not stronger. What it can inform is a protocol's missed-dose
contingency, which is a document-level decision.

Needs no identifiability at all, which is the point. `theta` is swept as a percentage of
the model's normal (100, 20, 10, 5, 2, 1) and never estimated or cited; the objective is
marker exposure and not the sign of `theta - theta_crit`, whose threshold Section 4.6
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
