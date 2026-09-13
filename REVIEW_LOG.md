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
