# Proposed: the constancy assumption, stated with its failure threshold

*Drafted, not applied. `paper/` is frozen for prose. This is PHASE4 passage 14's surviving content,
which R162 kept when the other three orphans were withdrawn: **the paper assumes capacity is constant
across a 730-day record, says so nowhere, and already contains the result showing the assumption is
load-bearing.***

---

## Why this one survived the reading and the other three did not

R162 read all four orphaned passages against the built PDF. Three were made or made moot by the
rewrite that removed their anchors. This one was not, and the searches are the evidence:

| searched, unbounded, in the 45-page build | hits |
|---|---|
| `"over days"` | **0** |
| `"act over"` | **0** |
| `dialysis` | **1**, inside a reference title only |
| `drift` | **1**, and it is between-visit **assay** drift, a different quantity |

**Every precision in the paper is computed on a monthly record over two years. Nothing in the
document says that capacity is assumed not to move over those two years.**

## What makes it worth stating rather than merely true

An unstated assumption is worth one sentence. **An unstated assumption with a known failure
threshold is worth a remark**, and this one has a threshold already computed in
`audit/PHASE11_DRIFT.md`.

**And the threshold has two values, for two different questions.** Both are in that document and
conflating them would misstate the result:

| arm | what it asks | crossing |
|---|---|---|
| **Tier 1 as written**, `n = 118` | does the requirement in the specification survive drift? | **fails at `r = 0.005`**, 0.5% a year: rmse **0.2020** against a gate of 0.20 |
| **schedule re-optimised**, `n*` free to 240 | can any monthly schedule recover it? | crosses at `r ≈ 0.009`, **and the re-check says read it as "of order one per cent a year, not as a sharp threshold at 0.0091"** |

**The 118-draw row is the one that belongs in the paper**, because 118 is the specification's own
number and `PHASE11_DRIFT.md` records why it is brittle: *"118 was **defined** as the count that
lands on the gate"*, so **the Tier 1 requirement has no margin against drift at all.** At
`r = 0.005` the standard deviation is still 0.2006 and it is the **bias**, −0.0241, that pushes the
rmse over.

**What it does not establish, and the proposed text must not overclaim.** `PHASE11_DRIFT.md` §6:
*"It does not establish that capacity drifts. Nothing here measures `r`, and the paper's own kinetic
argument for constancy is untouched."* The statement is a **conditional**.

---

## The sharpening that makes it a ceiling rather than a price

**The excess over the gate at `r = 0.005` is mostly bias, and the arithmetic says by how much.**
In variance terms, `sd^2 - gate^2 = 0.000240` and `bias^2 = 0.000581`, so **the bias is 70.7% of
the excess.** That much says drift dominates. It does not yet say the requirement cannot buy its
way out.

**What makes it a ceiling is where the bias comes from.** Round 1 established that the bias tracks

```
    bias  ~  - r * span / 2
```

and at `r = 0.005` over 118 monthly draws the span is 3587 days, or 9.82 years, giving
**-0.02455 against the reported -0.0241**. The bias is set by the record's **SPAN**, not by its
**COUNT**.

**Two consequences, and together they close the exit.**
- **More draws inside a fixed window cannot help.** They shrink the variance and leave the bias
  exactly where it is.
- **Lengthening the window to fit more draws makes it worse**, because the span is what the bias
  scales with.

**So the sentence is not that no number of draws reduces the second term, it is that no SCHEDULE
reduces it.** That is the strongest statement available from this result, and it is what makes 118
a boundary rather than a price: the requirement cannot be bought out of a drift by any
rearrangement of when the blood is taken.

---

## The proposed passage

*Placement settled: **beside the Tier 1 requirement**, not with the scope conditions. The finding
is specific to that number. 118 is defined as the count landing exactly on the gate, so it carries
no margin and drift is the first thing to spend it; as a scope condition it would read as a
general caveat about a general assumption. And the delivery pass found this paper's characteristic
weakness to be qualifications that are present, correct and late, so a qualification belongs where
its number is.*

> \begin{remark}[Capacity is assumed constant across the record, and the requirement has no margin
> against it]
> \label{rem:constancy}
> Every precision in this paper is computed on a record spanning two years, and each assumes
> $\theta$ does not move over that span. The assumption is not incidental to the design. A
> drifting capacity \emph{biases} the fit rather than merely widening it, so the quantity that
> meets the gate is $\mathrm{rmse} = \sqrt{\operatorname{sd}^2 + \mathrm{bias}^2}$.
>
> \textbf{The requirement has no margin against it.} The $\num{tier1_n}$-draw count of
> Section~\ref{sec:spec} is \emph{defined} as the number landing exactly on
> $\operatorname{sd}(\ln\hat\theta) = \num{gate}$, so any bias at all carries it over. At an
> exponential drift of $\num{drift_r_pct}\%$ per year the standard deviation is still
> $\num{drift_sd}$ and the bias is $\num{drift_bias}$, giving
> $\mathrm{rmse} = \num{drift_rmse}$: \textbf{the requirement fails on a drift whose halving time
> is $\num{drift_halflife}$ years.}
>
> \textbf{And no schedule recovers it}, which is what makes the count a boundary rather than a
> price. The bias scales as $-r\,T/2$ in the record's span $T$ and not in its number of draws, so
> sampling more densely inside a fixed window shrinks the variance and leaves the bias where it
> was, while lengthening the window to fit more draws increases it. Allowing the schedule to
> re-optimise moves the crossing only to of order one per cent a year.
>
> \textbf{This is a conditional and not a measurement.} Nothing here establishes that capacity
> drifts, no rate is measured, and the kinetic argument for constancy in
> Section~\ref{sec:intrinsic} is untouched. What it establishes is that the specification as
> written cannot absorb a drift it does not currently exclude, and that a decade-long study should
> record whatever would let $r$ be estimated afterwards rather than assumed.
> \end{remark}

## The macros, which are to be wired

**Decision taken: wire them.** They are computed, they can move, and R143's degenerate-macro
exception does not apply, since their content is not that they stay fixed.

| macro | value | owner |
|---|---|---|
| `tier1_n` | 118 | `src/design/run.py` *(already generated)* |
| `gate` | 0.20 | `src/design/run.py` *(check whether a key already exists before adding one)* |
| `drift_r_pct` | 0.5 | `audit/u7_drift.py`, the `n = 118` table |
| `drift_sd` | 0.2006 | same |
| `drift_bias` | −0.0241 | same |
| `drift_rmse` | 0.2020 | same |
| `drift_halflife` | 138.6 | same, `ln 2 / r` |

**A generator must be written or extended before any of these is cited**, per R141: a macro whose
value is typed by hand is worse than a literal, because the guard then vouches for it. `u7_drift.py`
computes all five drift figures today and does not emit a numbers file.

---

## What is left for the author

**Nothing on placement or wiring; both are settled above.** What remains is the prose freeze itself:
whether this passage, as worded, goes in.
