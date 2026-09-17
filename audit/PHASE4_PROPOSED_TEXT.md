# Proposed replacement text. NOT APPLIED.

**Nothing in this file has been written into `paper/`.** It is a proposal, held at the gate.
The paper and the supplement are frozen and stay frozen until the author signs off on the
wording below. Written after the adversarial stress test of `audit/PHASE4_STRESS.md`, so the
scope here is narrower than a first draft had it in two places and wider in two others, and
both changes are flagged where they occur.

Fifteen passages. Passages 9 and 10 were added after Phase 5, from the primary sources
supplied in September 2026; they are flagged NEW where they occur, as passages 5 to 8 were.
Passages 11 to 15 were added after the September 2026 verification pass
(`audit/PHASE6_REVIEWS_SEP2026.md`), once network egress made the primary sources readable
rather than merely citable. Passage 11 carries retraction `R17`; passage 15 is the second half
of the same retraction and is kept separate because it moves in the opposite direction.

---

## 1. `paper/sections/06_excitation.tex:71-75`, two clauses only

**Not** the whole span. `06_excitation.tex:66`, "flat only to first order", is a true and
reproducible statement about `prop:rank`'s null *vector* (local power 2.0000 to 2.0002 over a
decade) and stays.

### Current

> The combination is therefore recovered through curvature rather than through
> slope, and the consequence is a price rather than a prohibition: the recoverable
> displacement scales as the square root of the noise level, where a first-order
> direction would scale linearly. Halving the assay CV buys a factor of $\sqrt2$ along this
> direction and a factor of $2$ along every other.

### Proposed

> That is a statement about the null \emph{vector}, and the degeneracy is a \emph{curve}.
> Solving the balance for the increment that produces a level $m$ gives, with no remainder,
> \begin{equation}
> \Delta(m) \;=\; \frac{A}{m} \;+\; B \;+\; C\,m,
> \qquad A = K_m(L+b), \quad B = L + b - \theta - k_e K_m, \quad C = -k_e,
> \label{eq:invariants}
> \end{equation}
> so the increment that produces a given level depends on the parameters only through
> $(A,B,C)$. Since $A>0$ and $C<0$ the map is strictly decreasing from $(0,\infty)$ onto
> $\mathbb{R}$, hence a bijection, so the level a known $\Delta$ produces is a function of
> $(A,B,C)$ alone; and three distinct levels determine $(A,B,C)$, the relevant determinant
> being $-\prod_{i<j}(m_i-m_j)/(m_1m_2m_3)$. So on this locus, at a constant recorded load
> and a record visiting at least three levels, $(A,B,C)$ is exactly what the data contain.
> Four unknowns against three invariants leave a one-parameter family of exact
> observational twins,
> \begin{equation}
> K_m(b) = \frac{K_m^{0}\,L}{L+b},
> \qquad \theta(b) = \theta^{0} + b + k_e\bigl(K_m^{0} - K_m(b)\bigr),
> \qquad k_e \text{ unchanged},
> \label{eq:fibrecl}
> \end{equation}
> which is the integral curve of the null direction and along which the marker itself is
> unchanged. Proposition~\ref{prop:rank}'s null vector is its tangent at $b=0$, so the
> $\varepsilon^4$ above measures the curvature of \eqref{eq:fibrecl} against that tangent.
> The combination is therefore not recovered at second order. On this locus it is not
> recovered at all, and halving the assay CV buys nothing along it.
>
> Two things close it from outside the model, and the paper uses both. A bound on the load
> bias closes it directly: \eqref{eq:fibrecl} is admissible for every $b > -272.7$, and at
> the $5\%$ load record of Section~\ref{sec:spec} it leaves $\theta$ open by a factor of
> $1.028$. And a recorded load that \emph{varies} closes it from the data, because $L_k$ and
> the constant column separate; that is the last clause of Proposition~\ref{prop:rank} and
> the reason Section~\ref{sec:spec} asks for ordinary variability to be recorded rather than
> smoothed. Where the data do the work the price is ordinary: the half-width scales as
> $\sigma^{0.98}$, not as $\sqrt\sigma$.
>
> A singular Fisher matrix on this locus does not contradict the structural verdict of
> Section~\ref{sec:model}, because the two are statements about different records: the
> structural result reads a transient, from which every time derivative of the output is
> available, and the clinic samples a locus on which $\dot m$ has already vanished. That
> sentence is the whole reconciliation, and is not stated again.

### Why

Both retracted clauses are false, and the second is now refuted quantitatively rather than
structurally. "Global injectivity of the map" is dropped from the reconciliation sentence
because on this locus at a constant load the map is exactly non-injective; the Julia verdict
is about the transient, which the section's own opening already says.

---

## 2. `paper/sections/06_excitation.tex:497-516`, `rem:second_order`

### Proposed

> Proposition~\ref{prop:rank} is a statement about the design matrix. The information filter
> of \eqref{eq:infofilter} is a first-order object, so when it meets a rank-deficient design
> it fills the deficient direction from the prior and reports the result as though the data
> had supplied it. At a constant recorded load that is exactly what happens: the load-record
> prior gives a $95\%$ half-width of $\pm31.6$ in $b$, and the likelihood supplies nothing,
> because \eqref{eq:fibrecl} is an exact fibre.
>
> A diagnostic that reads the power of the noiseless $\chi^2$ along the smallest
> eigendirection does not detect this, and an earlier version of this remark drew the
> opposite conclusion from one. Moving along the straight eigendirection departs from the
> curved fibre at second order, so $\chi^2 \propto \varepsilon^4$ whether or not the data
> constrain anything. Two further properties of that measurement are worth recording,
> because both are easy to mistake for evidence. It is \emph{parametrisation dependent},
> being second order: the same direction walked in $(\theta, K_m)$ rather than in
> $(\ln\theta, \ln K_m)$ changes the constant by a factor of three. And it is strongly
> \emph{design dependent}: the half-width runs $36$, $49$, $69$ across two-daily, weekly and
> monthly sampling, so no single value of it can speak for every covariance this note
> reports. The diagnostic that does work costs more and is the one to run: re-optimise the
> remaining parameters at each fixed value of the suspect coordinate. That \emph{profile} is
> flat, and its optimum traces \eqref{eq:fibrecl}.
>
> What does the work instead is the last clause of Proposition~\ref{prop:rank}. At the
> day-to-day variability the specification asks for, the fibre closes and the profile becomes
> finite, though widely spread over realisations of the intake. Across eight simulated
> records the $95\%$ half-width in $b$ runs $2100$ to $4200$ at a load coefficient of
> variation of $0.5\%$, $420$ to $690$ at $2\%$ and $80$ to $130$ at $10\%$. Those remain
> wider than the prior's $\pm31.6$,
> so the prior stays the tighter of the two; what changes is that the direction is no longer
> exactly unidentified. Either way $\operatorname{sd}(\ln\hat\theta)$ barely moves: on the
> accessible monthly design the load-bias prior is worth $1.4\%$ of it.
>
> This is the second time in this paper that a linearisation has misrepresented an estimand;
> the first is Remark~\ref{rem:no_delta}. They do not point the same way, and the difference
> is the lesson. There the linearisation \emph{understated} what the record contained. Here a
> curvature measured along a ray was read as the record containing something, when the exact
> statement is that on this locus it contains nothing in that direction. A diagnostic is a
> model too, and a power law measured along the wrong curve answers a question the design
> cannot pose.

### Why, and what changed since the first draft

The `±31.6` against `±36.9` comparison goes, but **not** for the reason first given. The
first draft charged that the two were computed on different designs, with `±54.2` on the
monthly design. That charge was wrong: `χ²/ε⁴` is parametrisation dependent, `prop:rank` is
stated in log parameters, and in those coordinates the two-daily design reproduces the
printed pair to 0.5%. The surviving charges are that `±36.9` is a ray rather than a profile,
and that one design-dependent number is being asked to defend covariances on three designs.

---

## 3. `paper/sections/07_design.tex:141`, the mechanism clause only

### Current

> ... because along the near-null direction of Proposition~\ref{prop:rank} the gradient is
> third order and a started-away optimiser does not climb back.

### Proposed

> ... because the direction of Proposition~\ref{prop:rank} is an exactly flat fibre at a
> constant load and a long shallow valley at the variability a real record carries, so a
> started-away optimiser reports the spread it started with rather than the spread the data
> imply.

### Why

"The gradient is third order" is true of the straight line and says nothing about the fibre,
where every derivative vanishes. The replacement names the mechanism the audit observed,
which is start dependence, and which Phase 1's variable-projection fit is designed to remove.
Narrower than the first draft, which called the existing sentence false.

---

## 4. `paper/sections/07_design.tex:417-419`, `rem:no_delta`

### Current

> This is the first of the paper's two linearisation failures; the second is
> Remark~\ref{rem:second_order}, where a first-order filter reports a curvature-held
> direction as prior-held. They point the same way.

### Proposed

> This is the first of the paper's two linearisation failures; the second is
> Remark~\ref{rem:second_order}, where a curvature measured along a ray was read as the data
> pinning a direction that an exact fibre leaves free. They do not point the same way, and
> Remark~\ref{rem:second_order} says why.

### Why

Exactly inverted as it stands, and the false parallel goes with it.

---

## 5. `paper/sections/06_excitation.tex:371`, a step in the proof. NEW

This is independent of everything above and should be fixed whatever is decided about the
rest.

### Current

> The increment is eliminated by the balance itself: $(\theta+\Delta)\,u = L - k_e m$ ...

### Proposed

> The increment is eliminated by the balance itself: $(\theta+\Delta)\,u = L + b - k_e m$ ...

with the consequent changes carried all the way through, which is what makes this a fix
rather than a patch: \eqref{eq:nullvec} becomes
$n = (L + b + k_e K_m,\; -K_m,\; 0,\; L + b)$, the normalisation line reads "dividing
through by $\theta(L+b)$" rather than by $\theta L$, and \eqref{eq:rank} becomes
\[
c_\theta = \frac{1 + K_m k_e/(L+b)}{\theta}, \qquad c_{K_m} = -\frac{1}{L+b}.
\]
Only then is the fibre the integral curve of the null direction, which is how the proposed
Section 6 text describes it. Doing half of this would leave the new text asserting something
the amended proposition still would not support.

### Why

The balance carries the bias and the elimination drops it, so `eq:nullvec` as printed is the
`b = 0` value of a field that depends on `b`. It annihilates the Jacobian at the truth and
fails to at `b = 31.6`. Carrying `b` through, the field's integral curve is the exact fibre,
with symbolic difference zero; integrating the printed vector instead gives
`K_m = K_0 exp(-b/L)`, which leaves the fibre at second order. This is the reason the exact
fibre is not derivable from the paper as written.

---

## 6. `paper/sections/06_excitation.tex`, inside `rem:rankvsprec`. NEW

The first draft of this document called `rem:rankvsprec` correct as it stands. It is not.

### Current

> At any real day-to-day variation the data fill it instead and the same $0.1696$ is earned
> rather than assumed.

### Proposed

> At the variability the specification asks for the data fill it instead, though not at once:
> the same table's $0.5\%$ and $2\%$ rows still return $0.86$ and $0.26$ prior-free against
> the $0.1696$ the filter reports, so what the prior is holding is given back gradually and
> is not fully given back until the record carries several per cent.

### Why

The claim is false at two of the four rows the remark itself prints. Its reading of the
zero-CV row, that the prior is filling the deficient direction, is correct and stays.

---

## 7. `paper/sections/07_design.tex:25` and `A1_supplementary.tex:516`. NEW

Added by the stress test, and treated the way passage 3 was: replace the mechanism clause,
keep the empirical conclusion.

### Current

> $K_m$ lies inside the exact degeneracy, so a prior on it cannot help ... This is a reason,
> not a convenience.

### Proposed

> $K_m$ moves \emph{along} the degeneracy rather than across it, so a prior on it does bound
> the load bias, through $b = L(K_m^{0}/K_m - 1)$: a log-prior of standard deviation $s$
> confines $b$ to $[L(e^{-1.96s}-1),\, L(e^{1.96s}-1)]$. What that shows is how tight such a
> prior would have to be. At $s = 0.8$ it leaves $b$ in $[-255, +1225]$, forty times the load
> record's own $\pm31.6$, and matching the record would need $s = 0.048$, a factor of $1.10$
> on $K_m$. No published figure is within reach of that, and
> Remark~\ref{rem:no_km_orig} forbids importing one in any case. This is a reason, not a
> convenience.

### Why

The mechanism as stated is inverted, and the fibre is what shows it: a $K_m$ prior is
precisely a prior on position along the fibre. With the load prior removed and $k_e$ free, a
$K_m$ log-prior alone takes the filter from singular to $\operatorname{sd}(\ln\hat\theta) =
0.2820$ at $s = 0.8$ and $0.1854$ at $s = 0.2$ on the monthly design. So it helps; the
conclusion survives because the help needed is unobtainable, which is a stronger argument
than the one in place.

---

## 8. `src/common/params.py`, the `D` row's justification. NEW

From Phase 2 item 4, `audit/h1_admissibility.py` check 8 and `audit/h2_recheck.py` check 3.
This is the **only** failing check of the twenty, and it is a wording, not a constraint.

### Current

> `ASSUMED`, "chosen so the marker falls by about half at peak effect, the order of the
> reductions reported in the phase 1/2 interim analysis~\cite{trial_phase12}. Only the
> excursion it produces is used. ..."

### Proposed

> `ASSUMED`, "chosen so the marker falls by **about two thirds ($62$--$66\%$, depending on
> whether the fall is read from the untreated level or the pre-dose trough, and on the first
> dose or the steady-state peak)** at peak effect, the order of the reductions reported in
> the phase 1/2 interim analysis~\cite{trial_phase12}. Only the excursion it produces is
> used. ..."

The rest of the row is unchanged, including the recoverability sentence.

### Why

Computed in exact rationals, the marker falls from $10.0000$ to $3.5099$ at the first dose
($64.90\%$), from $10.0000$ to $3.3638$ at the steady-state peak ($66.36\%$), and from the
pre-dose trough $8.9106$ to $3.3638$ ($62.25\%$). **Every reading is larger than half**, the
smallest by a factor of $1.24$. It is not a rounding artefact: `h2` check 3 recomputes all
three symbolically.

Two reasons this is a wording fix and not a substantive one. The row's own next sentence,
"Only the excursion it produces is used", is the honest reading: nothing in the design turns
on the fall being any particular fraction, and the excursion is what is actually measured, as
$\Delta/\theta$ running $2.1333$ at peak to $0.1333$ at trough. And the "from untreated"
readings are anchored on a definition rather than a measurement, since $m(\Delta=0)=10$
exactly by the construction of $L$.

So the order of magnitude is right and only the word is wrong. **If the author prefers to
leave the row as it stands, nothing downstream breaks**; the proposal is offered because a
reader who checks the number will find it does not match, and a justification that does not
survive checking costs more than the word saves.

---

## 9. `paper/sections/A1_supplementary.tex`, after `\S`"The load-bearing part of the form". NEW

From Phase 5 steps 2, 3 and 4: `audit/k2_structures.py`, `audit/k5_identifiability.py`,
re-checked in `audit/k3_recheck.py` and `audit/k6_recheck.py`.

The supplement already ends that paragraph with

> The obstruction survives a change of structure; this parametrisation of it does not.

That sentence is now measured rather than asserted, and what the measurement says is sharper
and less comfortable than the sentence. **Proposed: keep the sentence and add the remark
below immediately after it.**

### Proposed

> \begin{remark}[$\theta$ is identified relative to a structure]
> \label{rem:structure_relative}
> Both halves of the preceding sentence can be checked, by generating records from each
> published structure and fitting \eqref{eq:A} to them.
>
> \paragraph{What survives.} The obstruction does. A record whose samples all sit at one
> dose-cycle phase leaves $K_m$ unbounded under \emph{every} structure tested, and breaking
> the period lock bounds it under every structure whose marker is algebraic in the enzyme.
> So does the design advice: the ranking of schedules by $\mathrm{sd}(\ln\theta)$ is
> unchanged, and the two recommended designs agree across structures to about $3\%$. What
> does \emph{not} survive is any hope that the recommended schedules would also warn of the
> substitution: against a turnover marker they do, leaving two to three $\sigma$ of residual,
> but against the algebraic one no schedule tested raises a signal at all.
>
> \paragraph{What does not.} The value. Fitting \eqref{eq:A} to records generated by the
> translational model's algebraic $I_{\max}$ structure~\cite{pkpd_translational} returns
> $\hat\theta$ high by a factor of $1.3$ to $2.0$ at this paper's own dose, and by a factor
> of $4$ to $20$ at a dose calibrated to the falls the phase 1/2 interim analysis actually
> reports~\cite{trial_phase12}. The residual of that fit can be under a fifth of one
> $\sigma$, in which case \emph{nothing in the fit signals the error}: the misspecification
> is then invisible in exactly the diagnostic a reader would use to detect it. The direction
> of the bias is robust to how that structure is coupled to the substrate load and to which
> estimator is used; its \emph{magnitude} is not, and is several times larger once this
> paper's own priors are removed, so the figures above are the favourable case.
>
> \paragraph{Why, exactly.} Not estimator failure. Under that structure the marker is
> $E_0 + \mathrm{base}\,(1 - I_{\max}E/(\mathrm{IC}_{50}+E))$ with $E = \theta + \Delta$,
> which is a M\"obius function of $\Delta$,
> \begin{equation}
> m(\Delta) \;=\; \frac{P\Delta + U}{\Delta + V},
> \qquad P = E_0 + \mathrm{base}\,q,\quad U = P\theta + (E_0+\mathrm{base})\mathrm{IC}_{50},
> \quad V = \mathrm{IC}_{50} + \theta,
> \label{eq:moebius}
> \end{equation}
> with $q = 1 - I_{\max}$. A record therefore carries three numbers about that map, and the
> structure has four free parameters at a fixed $I_{\max}$. Eliminating
> $(E_0, \mathrm{base}, \mathrm{IC}_{50})$ leaves \emph{no} relation on $\theta$: for any
> $\theta \in (0,V)$, setting
> $\mathrm{IC}_{50} = V-\theta$, $\mathrm{base} = (U-PV)/((1-q)(V-\theta))$ and
> $E_0 = P - q\,\mathrm{base}$ reproduces the record exactly. Under that structure $\theta$
> is not imprecise; it is \emph{not identifiable}, at any design and at any number of visits.
> So the estimator was reporting a point on a fibre along which the data are silent, and
> \eqref{eq:A}'s parametrisation chose which point.
>
> \paragraph{And whether it is identifiable at all is itself an assumption.} That structure
> was fitted to a fixed-diet animal colony and carries \emph{no substrate load}, so where a
> load enters it is a modelling choice and not a transcription. Of four plausible choices,
> three leave $\theta$ unidentifiable exactly as above and one restores it: a load that lifts
> only the suppressible pool $\mathrm{base}$, separating it from the unaffected floor $E_0$.
> Fixing $E_0$ to zero, or measuring $\mathrm{IC}_{50}$ externally, restores it too. What
> creates the deficit in the first place is carrying a \emph{free marker offset alongside a
> free marker amplitude}, which is what a phenomenological readout has and a flux balance such
> as \eqref{eq:A} does not. So the statement a reader should take is not that $\theta$ is
> unidentifiable under this structure, but that \emph{whether it is identifiable turns on an
> assumption nobody has made}.
>
> \paragraph{And it is not generic.} The other published structure, a turnover marker whose
> production the therapy inhibits~\cite{dose_selection_poster,mrna_moa}, \emph{is}
> identifiable in $\theta$, at a constant load and with no external information. The
> difference is two-fold and it is worth naming, because it says what to look for in any
> third structure: its marker is a \emph{state}, so the transient pins its own rate constant
> separately from the steady map, which buys one number; and it carries no free marker
> scale, where the algebraic structure carries two ($E_0$ and $\mathrm{base}$) and spends
> one. Section~\ref{sec:excitation}'s own model sits between the two: $\theta$ is
> identifiable, but only once the recorded load varies.
>
> \paragraph{What follows for the reader.} No claim about the \emph{value} of $\theta$ can
> rest on goodness of fit, because a well-fitting wrong structure is exactly the case that
> goodness of fit cannot see. The check has to come from outside the longitudinal record,
> which is what the tiered specification of Section~\ref{sec:spec} is for: Tier~1 is a direct
> measurement rather than an inference from the record, and it is the only anchor in this
> paper that a structural bias in $\hat\theta$ could not pass through unnoticed. Everything
> this paper claims about \emph{designs} is structure-robust on the evidence above.
> Everything it might be read to claim about a \emph{number} is not.
> \end{remark}

### And one pointer, `paper/sections/12_discussion.tex`, "What would change the answer"

Add a fourth item:

> \item \textbf{An external anchor on $\theta$.} $\theta$ is identified relative to the model
> structure, and Remark~\ref{rem:structure_relative} shows a published alternative structure
> under which it is not identifiable at all, fitted at a residual too small to notice. Tier~1
> of Section~\ref{sec:spec} is the only anchor here that sits outside the longitudinal
> record. A second one, independent of it, would turn every statement about the value of
> $\theta$ from conditional into checkable.

### Why

Three reasons this belongs in the paper rather than in the audit.

1. **The supplement already makes the claim**, in one sentence, with nothing behind it. The
   remark supplies what is behind it, and the answer is not the comfortable one: the half of
   the sentence that says the parametrisation does not survive is doing much more work than
   its phrasing suggests.
2. **The failure mode is invisible in the standard diagnostic.** A reader who fits this model
   and sees a residual of $0.8\sigma$ on twelve points and four parameters will conclude the
   model is adequate. It is exactly the case where that conclusion is wrong, and the paper
   should say so in the place where a reader looks for it.
3. **It strengthens the paper's own recommendation.** Tier~1's justification in
   Section~\ref{sec:tier1} is currently about the precision of $\Delta$. This adds a second,
   independent reason for the same requirement, and the second reason does not depend on the
   $12\%$ figure being right.

**Scope, stated so it can be narrowed.** The remark asserts (a) the M\"obius reduction and
its consequence, which is proven symbolically and by exact rational interpolation, confirmed
by the numerical rank of the sensitivity map and by walking the fibre explicitly over a
$990$-fold range of $\theta$ at a bit-identical record; (b) the bias figures, which are
measured over $8$ diet records per cell, are calibration-dependent with both branches carried,
and are the \emph{favourable} end of a range that widens to $5$--$9$ without this paper's own
priors; (c) the schedule-ranking robustness, measured on five schedules and three load
variabilities; and (d) the insertion dependence, measured over four couplings of the load plus
the uncoupled case. If the author wants a shorter version, (a) and (d) carry the argument,
(c) is the reassuring half, and (b) can be cut to a footnote, though (b) is the part a referee
will remember.

**One thing the remark deliberately does not claim.** The three-number count behind (a) is
elementary within a stated class (structures whose fixed-load $(\Delta,m)$ graph is cut out by
a polynomial with at most four monomials) and is verified for the three structures at issue.
It is \emph{not} offered as a general theorem, and the remark should not be edited into one.

---

## 10. Five annotation rows that read mouse estimates as human. NEW

From Phase 5 step 0, `audit/k1_sources.py` check 2, against the primary sources supplied in
September 2026.

**No value moves.** `60.1`, `0.0075`, `0.0987` and `51%` are all annotation prose;
`M_BAR` stays `10`, `LAMBDA` stays `ln 2/3.5`, `ASSAY_CV` stays `5`. `tab_params.tex` is
generated from `params.py` by `make figures`, so each item below is one edit, not two.

**The fact all five turn on.** Attarwala's PD parameters are fitted in `Pcca-/- (A138T)`
mice, `N = 111`, and the paper's own scaling sentence covers PK only:

> "The preclinical PK/PD model described earlier was extrapolated to humans based on
> allometric scaling of model parameters. Based on the PCCA/B mRNA PK model described
> earlier, estimated allometric exponent was 0.63 for CL12 and CL32 and 1.1 for CL23 and
> CL20. Allometric exponent for the volume of distribution parameters (ie, V and V2) was
> fixed to 1."

Every parameter named there is a clearance or a volume. The companion paper says the same in
one clause: *"Each preclinical PK/PD model was extrapolated to humans based on allometric
scaling of PK parameters"*~\cite{mrna_moa}. So `E0`, `base`, `IC50`, `Imax`, `Slope`, `kdeg`,
`kq`, `ke0` and both inter-individual variabilities are murine.

### 10a. `src/common/params.py:60-63`, the `$\bar m$` row

**Current**

> "On a real scale: the translational model puts the untreated plasma 3-HP baseline at
> $60.1$~$\mu$M and the fully corrected floor near $0.1$~$\mu$M~\cite{pkpd_translational}, so
> a managed patient sits between the two and $\bar m$ is read in $\mu$mol/L"

**Proposed**

> "On a real scale: the translational model of PCC-deficient \textbf{mice} puts the untreated
> plasma 3-HP baseline at $60.1$~$\mu$M and the fully corrected floor at
> $E_0 = 0.0987$~$\mu$M~\cite{pkpd_translational}. \textbf{Neither figure was scaled to
> human}: that paper scales its PK parameters allometrically and not its PD ones. The nearest
> \emph{human} anchor is the newborn dried-blood-spot $99.9$th centile of
> $38.13$~$\mu$mol/L~\cite{dbs_3hp_cutoff}, within a factor of $1.6$ of the murine baseline,
> so the \emph{scale} carries across even though the number is not a patient's. $\bar m$ is
> read in $\mu$mol/L."

### 10b. `src/common/params.py:107-113`, the `$\lambda$` row

**Current**

> "\textbf{Corroborated}: the translational PK/PD model estimates a PCC-protein degradation
> rate of $0.0075$~h$^{-1}$~\cite{pkpd_translational}, a half-life of $3.85$~d, so this is
> within $10\%$ of a published value rather than free."

**Proposed**

> "\textbf{Partly corroborated, and only early}: the translational PK/PD model of
> PCC-deficient \textbf{mice} estimates $k_{\mathrm{deg}} = 0.0075$~h$^{-1}$, but that is one
> rate constant inside a \emph{two-compartment} protein model whose own stated protein
> half-life is $\sim 7$~d~\cite{pkpd_translational}. The protein declines biexponentially, so
> its half-life is a function of time: the model's local half-life passes $3.5$~d at about
> five days after a dose, reaches $9.0$~d by the end of a fortnightly interval, and tends to
> $11.8$~d. So $\lambda$ tracks the published decline over the first week and decays
> \emph{faster} than it thereafter, which is the direction that makes the inter-dose
> excursion look larger than it is. The $3.85$~d previously quoted here is $\ln 2/k_{\deg}$,
> which is not a half-life of the protein."

### 10c. `src/common/params.py:135-136`, the `CV` row

**Current**

> "Biological variability is far larger, the translational model carrying $51\%$ between
> subjects on the 3-HP baseline~\cite{pkpd_translational}"

**Proposed**

> "Biological variability is far larger, the translational model carrying $51.1\%$ between
> \textbf{mice of one inbred strain} on the 3-HP baseline~\cite{pkpd_translational}, which
> bounds between-\emph{patient} variability from below rather than above. The phase 1/2
> interim analysis reports standard deviations of $24$ to $57$ percentage points on the
> per-cohort mean 3-HP change, on cohorts of three or four~\cite{trial_phase12}; that figure
> confounds assay, within-patient and between-patient variation and bounds each from above."

### 10d. `src/common/params.py:82-83`, the `$k_e$` row

Not a species error. A structure error, and the only one of the five that names a route that
does not exist.

**Current**

> "The nearest route is the first-order 3-HP elimination implicit in the published
> indirect-response model~\cite{pkpd_translational}"

**Proposed**

> "The nearest route is the first-order marker loss of the \emph{turnover} model of the
> markers~\cite{dose_selection_poster,mrna_moa}, which publishes no value for it. The
> translational model~\cite{pkpd_translational} does \textbf{not} supply one: its marker is
> an algebraic function of the enzyme, $\mathrm{3\text{-}HP} = E_0 + \mathrm{base}(1 -
> I_{\max}\mathrm{PCC}/(\mathrm{IC}_{50}+\mathrm{PCC}))$, with no marker state and no
> elimination constant anywhere in it."

`paper/claims.csv` C46 already records this ("no biomarker turnover constant at all"), so
this brings the row into line with a note the repository already holds. Supplementary
Table S1(C) of `pkpd_translational` is the equation.

### 10e. `paper/refs.sources:102-116`, the `pkpd_translational` annotation

**Current**

> "Semi-mechanistic translational PK/PD model for mRNA-3927, fitted largely to preclinical
> data and allometrically scaled to human."
> ...
> "PROVENANCE: these figures are transcribed from the author's parameter-grounding review of
> September 2026, not read off the paper in this repository's own pass. The citekey is
> confirmed to be the right article; the numbers await a check against the primary source."

**Proposed**

> "Semi-mechanistic translational PK/PD model for mRNA-3927. The PK half is fitted across
> mouse, rat and monkey and scaled allometrically to human (exponents 0.63 on CL12/CL32,
> 1.1 on CL23/CL20, 1 on the volumes). The PD half is fitted in Pcca-/- (A138T) mice, N=111,
> and is NOT scaled: every PD figure below is a murine estimate.
>   PCC-protein degradation rate kdeg = 0.0075 /h. This is ONE rate constant in a
>     two-compartment protein model, not the protein's half-life; the paper states the
>     latter as ~7 d. kq = 0.00474 /h, ke0 = 0.242 /h.
>   plasma 3-HP disease baseline 60.1 umol/L (14.4% RSE), IIV 51.1%; floor E0 = 0.0987.
>   IC50 for 3-HP 37.5 ug/g PCC (12.2% RSE), IIV 73.8%.
>   Slope 55.6 (1.93% RSE), printed as (mg/g/h)/(mg/mL). It is the LINEAR coefficient from
>     effect-compartment mRNA to PCC synthesis rate, Ksyn = Ce*Slope, NOT a Hill exponent.
>   Imax FIXED at 0.999, not estimated (Table 1, footnote m).
>   The marker equation is ALGEBRAIC in the enzyme (Suppl Table S1(C)), so there is no
>     first-order 3-HP elimination in this model and k_e cannot be extracted from it.
>   TWO UNIT MISPRINTS, on the paper's own internal evidence. 'mg/g' in the IC50 rows is
>     ug/g: Fig. 2B's axis is labelled 'PCC protein, ug/g' and runs to ~700, so mg/g would
>     put the IC50 a thousandfold above every fitted point. 'mM/L' in the marker rows is
>     umol/L: the same text gives the 3-HP LOQ as 25 'mM/L', and Baek et al. Fig. 2e plots
>     3-HP in uM, 0 to 200, in the same mouse strain. ('mM/L' is not a unit in any case.)
> PROVENANCE: read off the primary source (full text, Supplementary Table S1 and
> Supplementary Figures S1-S3) in this repository's own pass, September 2026."

### Why, and what this does not claim

The provenance warning was the repository's own, and it can be lifted: all six figures are
confirmed against Table 1, and the two unit readings the repository had already guessed
correctly are now confirmed from the paper's own figure axes rather than inferred.

**These are annotation corrections and nothing downstream moves.** No conclusion in the paper
depends on `60.1`, on `51%`, or on `3.85` d; they appear only as provenance for assumed
constants. What changes is whether a reader who checks them finds what the row says. Three of
the five (10a, 10c, 10e) say "mouse" where the row now says nothing; one (10b) replaces a
claim of corroboration with a bounded one; and one (10d) removes a named extraction route
that the cited paper does not contain.

**10b is the one with a downstream consequence, and it is not applied here either.** If
$\lambda$ were moved to $\ln 2/7$, a monthly record needs $30$ visits rather than $11$ to
reach $\mathrm{sd}(\ln\theta) = 0.20$. The row above keeps $\lambda = \ln 2/3.5$ and states
the qualification; moving the constant is a separate decision and is not proposed.

---

## Also needing a note, outside the paper

- `audit/PHASE2C.md:141-150`, where `±31.6` against `±36.9` first appears with the inverted
  reading.
- `audit/PHASE3_CHECKS.md:301-302`, the `sqrt(noise)` sentence.
- `audit/e5_structural.py` `q3`. Its computation is correct and its interpretation is what is
  retracted, so the fix is to append what the `eps^2` measurement does and does not
  establish, not to change any number.
- `audit/e9_rank_vs_precision.py`, whose printed commentary carries the `rem:rankvsprec`
  sentence corrected in passage 6.
- `audit/PHASE5_ROBUSTNESS.md` carries the measurements behind passages 9 and 10, and
  `audit/PHASE4_SUMMARY.md` now separates what is structure-robust from what is not.
- `audit/g_design_common.py` and every Phase 3 figure are computed at `lambda = ln 2/3.5`
  and `D = 6000`. If either constant moves, those figures move with it; the sensitivity
  table in `PHASE5_ROBUSTNESS.md` says by how much, and nothing has been regenerated.

---

## 11. `paper/sections/08_specification.tex:48-50` and `:118-120`. The load recording error. NEW

Carries retraction `R17`. **This is the largest single change proposed in this file**, because
it converts a sentence that reads as robustness into a requirement the design can fail.

### Current, at `:48-50`

> Two rows already report their own insensitivity: the load record gives $0.170$ whether it is
> kept to $1\%$ or to $20\%$, and $K_m$ needs nothing at all.

### Current, at `:118-120`

> \paragraph{The load record: ordinary precision, and keep the variation.} Recording intake
> twenty times more precisely buys about six per cent, so precision is not the binding
> variable.

### Proposed, at `:48-50`

> One row reports its own insensitivity: $K_m$ needs nothing at all. The load row reports two
> different quantities and only one of them is insensitive. How much the record \emph{varies}
> changes nothing ($0.170$ whether the day-to-day coefficient of variation is $1\%$ or $20\%$);
> how well it is \emph{recorded} is the second-tightest requirement in the table.

### Proposed, at `:118-120`, replacing the first sentence and keeping the rest of the paragraph

> \paragraph{The load record: the recording error binds, and keep the variation.} Two different
> things are called the load record and they behave oppositely. Making the record more
> \emph{variable} buys nothing, and making it twenty times less variable costs nothing. Making
> it more \emph{accurate} is the second-tightest requirement here, after the increment:
> \begin{center}
> \small
> \begin{tabular}{@{}lcc@{}}
> \toprule
> \textbf{substrate load recorded to} & $\sigma$ & $\operatorname{sd}(\ln\hat\theta)$ \\
> \midrule
> $1\%$  & $15.2$  & $\mathbf{0.157}$ \\
> $5\%$  & $22.0$  & $\mathbf{0.166}$ \\
> $10\%$ & $35.5$  & $\mathbf{0.192}$ \\
> $20\%$ & $66.2$  & $0.264$ \\
> $40\%$ & $130.0$ & $0.402$ \\
> \bottomrule
> \end{tabular}
> \end{center}
> so the design reaches $0.20$ only if intake is recorded to better than about $12\%$. That
> is not a comfortable margin. European practice records what patients are \emph{prescribed}
> rather than what they eat~\cite{pa_diet_survey}, no study in propionic acidemia reports
> repeated as-eaten intake at all, and the one study in any protein-restricted inborn error
> that has compared the two found actual natural protein running at $167\%$ of
> prescribed, a median of $40\%$ above it~\cite{pku_prescribed_actual}. At $40\%$ this design
> returns $0.402$, twice its own target. \textbf{The feasibility of everything downstream rests
> on a measurement of intake that is currently not made.}

### Why

`spec_numbers()` computes $\sigma$ once at a $5\%$ load error and reuses it on every row, so the
published "$1\%$ or $20\%$" rows sweep `load_cv`, the record's variation, and never touch the
recording error, which enters only through $\sigma$. Full working and the seed spread are in
`audit/m1_loadbias_prior.py` and `REVIEW_LOG.md` `R17`. **No constant is moved**: `LOAD_FRAC`
stays at `0.05` and the headline stays at $0.166$. What changes is that the table now says what
that $0.05$ is buying.

Two citations are needed and both are registered: `pa_diet_survey` (Daly et al. 2017,
PMID 29021961) and `pku_prescribed_actual` (Pinto et al., *Nutrients* 2023;15(23):4903).

---

## 12. `paper/sections/08_specification.tex`, a new paragraph after the load record. NEW

### Current

Nothing. The paper discusses cyclical intestinal antibiotics in
`A1_supplementary.tex` and in the "three readings" of `07_design.tex`, in both cases as a
\emph{pre-therapy instrument}, and never as an input that coexists with the dose train.

### Proposed

> \paragraph{Record the antibiotic weeks, and one reason that is not the reason you expect.}
> Most patients on this therapy are also on cyclical intestinal antibiotics. The consensus
> regimen is one to two weeks of metronidazole alternating with two to three weeks
> off~\cite{pa_guidelines}, and the measured effect on propionate metabolite excretion in four
> PA and five MMA children is a mean fall of $41\%$, range $12$ to $76\%$, with plasma
> propionate falling from $45.0$ to $25.1\,\mu$mol/L~\cite{metronidazole_pa}. So the
> non-enzymatic load is not a constant in these patients. It is close to a square wave of
> order $40\%$ amplitude on a three- to five-week period, running alongside a two- or
> three-week dose train and a monthly draw.
>
> Recorded, it is free and almost inert: it changes
> $\operatorname{sd}(\ln\hat\theta)$ from $0.166$ to between $0.167$ and $0.173$ at every period
> tried, which is what Section~\ref{sec:instruments} predicts for a perturbation an order of
> magnitude below $K_m$. \emph{Un}recorded, the same cycle lands in the residual and takes the
> design to between $0.297$ and $0.313$, missing the target by half. The instruction is
> therefore one line and it costs nothing: \textbf{record which weeks the patient was on
> antibiotics}, in the same record as intake and for the same reason.
>
> One row of that calculation is worth stating explicitly because it inverts the warning of
> Remark~\ref{rem:aliasing}. If the antibiotic period happens to equal the sampling interval,
> the unrecorded cycle becomes \emph{harmless} ($0.167$), because every draw then lands at the
> same phase, the swing is a constant offset, and the load-bias column absorbs it. Aliasing
> destroys the dose train, which is the \emph{signal}; aliasing rescues an unrecorded
> nuisance input. The two are not the same phenomenon and a reader who has just met the first
> will guess the wrong sign for the second.

### Why

**Both sources were already registered here, with correct annotations**: `pa_guidelines` and
`metronidazole_pa` (Thompson et al., *Eur J Pediatr* 1990;149(11):792-6, PMID 2226555, which is
a different paper from the Thompson 1990 *Metabolism* partitioning study also held as
`propionate_sources`). What is new is not the sources but the **conjunction**: nobody had put
the guideline's intermittent schedule next to the trial's effect size and asked what a periodic
`b` does to a design built on a periodic `Delta`. The `pa_guidelines` annotation has been
extended with the regimen, which was an open verification item.
`audit/m2_antibiotic_cycle.py` carries the sweep, 24 seeds, and cross-checks against the
independent $\sigma$ sweep of passage 11 (its unrecorded case runs at $\sigma \approx 68$, and
the $\sigma$ sweep gives $0.264$ at $\sigma = 66.2$). **What this passage does not claim.** The square wave is the shape the regimen implies. Nobody
has measured the shape, and no PA cohort reports antibiotic status alongside markers.

---

## 13. `paper/sections/08_specification.tex:73-75`. The PA breath-test cohort. NEW

### Current

> And the cohort is methylmalonic acidemia, not propionic acidemia, where the block sits at PCC
> itself.

### Proposed

> The cohort is also methylmalonic acidemia rather than propionic acidemia, where the block
> sits at PCC itself. That last point cuts less deeply than it appears, and in the
> favourable direction. The assay was first applied to this disease in five patients with
> propionic acidemia among twenty subjects, where recovery of $^{13}$CO$_2$ was significantly
> below normal and, notably, \emph{``well correlated with clinical severity in PA, but did not
> correlate in MMA''}~\cite{breath_test_original}. The measurement Tier~1 asks for is therefore
> not an extrapolation from the wrong disease. What is missing in propionic acidemia is
> narrower and more specific: a \emph{repeatability} figure, which that study does not report
> and which nobody has published since.

### Why

Barshop et al., *Pediatr Res* 1991;30(1):15-22 (PMID 1909779), verified. The current sentence
is true of `breath_test_mma` and reads as a statement about the literature, which it is not:
there are PA data, and the severity correlation in them is better than in MMA. Correcting it
strengthens Tier 1's rationale rather than weakening it, and it keeps the gap that Tier 1
actually fills stated precisely. Registered as `breath_test_original`.

---

## 14. `paper/sections/07_design.tex:253-259`, `rem:transient`. A second untested assumption. NEW

### Current

> Against that: the fast dynamics of a plasma overflow metabolite is not the fast dynamics of
> the mitochondrial substrate, and Remark~\ref{sec:marker} already names the monotone-tracking
> assumption as the weakest link in the chain. A transient interrogates exactly the part of
> the observation model that assumption is least able to support, at steady state the
> marker only has to track the substrate monotonically, whereas during a transient it has to
> track it with the right \emph{time constant}, and nothing in the literature establishes
> that.

### Proposed

> Against that, two things, and they are different objections. The first is the one
> Remark~\ref{sec:marker} already names: the fast dynamics of a plasma overflow metabolite is
> not the fast dynamics of the mitochondrial substrate, and at steady state the marker only has
> to track the substrate monotonically, whereas during a transient it has to track it with the
> right \emph{time constant}. Nothing in the literature establishes that.
>
> The second is about the plasma compartment itself, and it is the one piece of human evidence
> that bears on the $0.81$~h constant directly. After intermittent haemodialysis is stopped in
> an acute decompensation, the organic-acid pool does not stay drawn down: it refills, from
> continued extrahepatic production and from redistribution out of tissue, which is the
> standard argument for preferring continuous therapy~\cite{pa_guidelines}. A slow compartment
> behind the plasma pool is exactly the structure that would break the separation of $104$ this
> section is built on, and it would do so without violating monotone tracking at all. The
> observation is qualitative and it is made about crisis management rather than about steady
> state, so it does not move a number here. It does say that the two assumptions are
> independent, and that a transient design would be testing both at once.

### Why

The September 2026 clinical review's conclusion that capacity changes act over days does **not**
follow from its own evidence (post-transplant decay is consistent with a fast marker tracking a
graft whose function ramps over a week, and the review's Di Campli entry shows exactly that
ramp). The dialysis rebound is the part of that review that does bear on the question, and it
bears on a different assumption from the one the remark currently names. Barshop's suggestion
that PA and MMA differ in the \emph{``distribution of propionate pools''}~\cite{breath_test_original}
points the same way from the other direction.

**No number changes.** `rem:transient` remains an open direction and not a recommendation.

---

## 15. `paper/sections/08_specification.tex:120-126`. The other half of `R17`. NEW

The second correction inside the load-record paragraph, moving the opposite way to passage 11.

### Current

> By Proposition~\ref{prop:rank} a
> perfectly steady load makes the design exactly rank deficient, and that null direction must
> be filled by something; as the specification stands it is filled by the $5\%$ prior on the
> load bias, and only real day-to-day variation lets the \emph{data} fill it instead. So
> ``record intake as eaten, do not smooth it'' is the difference between a result the data
> support and one a prior supports.

### Proposed

> By Proposition~\ref{prop:rank} a perfectly steady load makes the \emph{untreated} design
> exactly rank deficient. Under the dose train it does not, because the increment varies even
> when the load does not: at a constant recorded load with the prior on the load bias removed
> entirely, $\operatorname{sd}(\ln\hat\theta)$ is $0.266$ rather than unbounded. The prior is
> worth a factor of $1.6$ there and is not holding the result up. So ``record intake as eaten,
> do not smooth it'' buys a real but bounded improvement, and it is the \emph{accuracy} of the
> record, not its variability, that decides whether the design works at all.

### Why

Measured in `audit/m1_loadbias_prior.py`: at a constant recorded load, `0.1659` at the published
$5\%$ prior, `0.1911` at $40\%$, `0.2659` with the prior flat. The published claim was checked
against the varying-load case and never against the constant-load case it is a statement about.
This is a **softening of the paper's own argument** and is proposed for that reason: the
sentence as written claims more for the instruction than the code supports, and passage 11
gives the instruction a better justification anyway.
