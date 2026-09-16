# Proposed replacement text. NOT APPLIED.

**Nothing in this file has been written into `paper/`.** It is a proposal, held at the gate.
The paper and the supplement are frozen and stay frozen until the author signs off on the
wording below. Written after the adversarial stress test of `audit/PHASE4_STRESS.md`, so the
scope here is narrower than a first draft had it in two places and wider in two others, and
both changes are flagged where they occur.

Six passages.

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
> and $\{1/m,\,1,\,m\}$ are independent, so $(A,B,C)$ is everything a constant-load record of
> pairs $(m_k,\Delta_k)$ on this locus contains. Four unknowns against three invariants leave
> a one-parameter family of exact observational twins,
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
> finite: a $95\%$ half-width in $b$ of about $2100$ at a load coefficient of variation of
> $0.5\%$, $509$ at $2\%$ and $96$ at $10\%$. Those remain wider than the prior's $\pm31.6$,
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

with the consequent change to \eqref{eq:nullvec}, which becomes
$n = (L + b + k_e K_m,\; -K_m,\; 0,\; L + b)$.

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

## Also needing a note, outside the paper

- `audit/PHASE2C.md:141-150`, where `±31.6` against `±36.9` first appears with the inverted
  reading.
- `audit/PHASE3_CHECKS.md:301-302`, the `sqrt(noise)` sentence.
- `audit/e5_structural.py` `q3`. Its computation is correct and its interpretation is what is
  retracted, so the fix is to append what the `eps^2` measurement does and does not
  establish, not to change any number.
