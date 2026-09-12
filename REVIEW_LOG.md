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
