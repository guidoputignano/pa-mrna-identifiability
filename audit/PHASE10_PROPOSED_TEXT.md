# Proposed replacement text, Phase 10. NOT APPLIED.

**Nothing in this file has been written into `paper/`.** It is a proposal held at the gate, in the
style of `audit/PHASE4_PROPOSED_TEXT.md`. The paper and supplement are frozen and stay frozen until
the author signs off on the wording.

One passage. It says what `A`, `B` and `C` ARE, which the paper currently proves everything about
without ever stating.

---

## 1. `paper/sections/06_excitation.tex:33-40`, an insertion after `\eqref{eq:invariants}`

### The problem

Equation `\eqref{eq:invariants}` introduces

```
Delta(m) = A/m + B + C m,    A = K_m(L+b),  B = L + b - theta - k_e K_m,  C = -k_e
```

and then proves, correctly and in one paragraph, that the map is a bijection, that three distinct
levels determine the triple, that the system is Vandermonde in the marker levels, and that four
unknowns against three invariants leave the fibre. **Every one of those is a statement about the
triple and none is a statement about what the triple MEANS.** A reader finishes the paragraph knowing
that `(A,B,C)` is what the data contain and not knowing what any of the three is.

That matters because the paper's two central facts both follow from the meanings in a line each, and
are currently derived by algebra instead: that a constant-load record pins `k_e` outright, and that
holding a low marker level is expensive in capacity.

### Proposed insertion, after "...so the increment depends on the parameters only through $(A,B,C)$."

> **The three coefficients are not opaque, and naming them gives the rest of this section for free.**
> $C$ is minus the non-enzymatic rate constant, nothing more: a bare parameter read straight off the
> linear term. $A/m$ is the \emph{sub-saturation penalty}, the extra increment the therapy must supply
> because the enzyme is working below saturation; it vanishes as $m$ grows past $K_m$ and grows
> without bound as the marker falls. $B$ is the \emph{shortfall}: the total load the patient carries,
> less the capacity they have, less the non-enzymatic flux evaluated at the half-saturation
> concentration.
>
> Read that way three things follow immediately. **A constant-load record pins $k_e$ outright**,
> because $C$ is $-k_e$ alone and the linear term is the one a record at several levels determines
> best; the degeneracy therefore cannot involve $k_e$, and \eqref{eq:fibrecl} duly leaves it
> unchanged. **Holding a low marker level is expensive in capacity**, because $A/m$ is the penalty
> term and management drives $m$ down: the better the patient is controlled, the larger the fraction
> of the increment that is spent on sub-saturation rather than on excitation. And **the fibre is the
> shortfall's ambiguity**: $A$ and $B$ both move with $b$, so a load bias can be traded against
> capacity and half-saturation at fixed $(A,B,C)$, which is \eqref{eq:fibrecl} in words.

### Two things the passage must also say, and why

**First, Section 7 already does exactly this, for a different quantity.** The scope limit is derived
there through the flux identity, $\theta = S(L+b)/g(m)$ with $S = 1 - k_e m/(L+b)$, which is the same
move: give a composite a name and read the conclusion off it. The proposed passage should say so,
because two independent instances of one technique read as a method and one reads as a trick.

**Second, and this is the part that cannot wait for the consolidation pass: the two-state reduction
has changed what two of these three symbols mean.** `audit/PHASE10_SENSORS.md` R.1 to R.3 proves that
the paper's one-state model is the exact reduction of a two-state model carrying intracellular
propionyl-CoA, and that under that reduction

- the paper's `k_e` is `k_e/alpha`, the disposal of propionyl-CoA down the overflow route expressed
  in marker units, **not the marker's clearance**, so `C = -k_e` is minus an overflow coefficient
  divided by the measured fraction; and
- the paper's `K_m` is `K_p c / k_e_paper`, **not the enzyme's Michaelis constant**, so the
  "half-saturation concentration" in `B` and the saturation scale in `A/m` are composites.

**So the passage must name the coefficients in the reduced model's own terms**, or it will teach a
reader a reading the reduction has already withdrawn. The wording above is deliberately written to
survive that: "non-enzymatic rate constant" and "half-saturation concentration" are true of the
one-state model as stated, and a single following sentence can carry the reduction:

> Under the two-state reduction of Appendix~\ref{sec:reduction} these last two are composites: the
> rate constant is the overflow coefficient divided by the fraction of that flux the assay measures,
> and the half-saturation concentration is the propionyl-CoA constant rescaled by the same ratio.
> Neither substitution changes any statement above, because all of them are statements about
> $(A,B,C)$ and the reduction leaves that triple invariant. It changes what a reader may import from
> in vitro enzymology, which is nothing.

### Status

**[judgement]** on the wording. **[proven]** on the three readings, which are the definitions of `A`,
`B` and `C` rearranged, and on the reduction, `audit/q13_two_state_reduction.py`, 5 of 5.

**Blocked on one author decision**: whether the two-state reduction goes into the paper at all. If it
does not, the last block must be dropped and the passage carries a reading the audit knows to be
incomplete. That is a reason to decide the reduction's status before this passage is applied, not a
reason to delay the rest of it.
