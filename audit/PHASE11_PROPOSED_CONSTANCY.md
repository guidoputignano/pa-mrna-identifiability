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

## The proposed passage

*For `07_design.tex`, near the Tier 1 requirement, or `11_scope.tex` beside the other scope
conditions. Placement is the author's; the content is the same either way.*

> \begin{remark}[Capacity is assumed constant across the record]
> \label{rem:constancy}
> Every precision in this paper is computed on a record spanning two years, and every one of them
> assumes that $\theta$ does not move over that span. The assumption is not incidental to the
> design: a drifting capacity biases the fit rather than merely widening it, so the relevant
> quantity is $\mathrm{rmse} = \sqrt{\mathrm{sd}^2 + \mathrm{bias}^2}$ and no number of draws
> reduces the second term.
>
> \textbf{The requirement has no margin against it.} The $\num{tier1_n}$-draw count of
> Section~\ref{sec:spec} is defined as the number that lands exactly on
> $\operatorname{sd}(\ln\hat\theta) = 0.20$, so any bias at all carries it over. At an exponential
> drift of $0.5\%$ per year the standard deviation is still $0.201$ but the bias is $-0.024$, and
> the root mean square is $0.202$: the requirement fails on a drift whose halving time is more
> than a century. Allowing the schedule to lengthen in response moves the crossing only to of
> order one per cent a year.
>
> \textbf{This is a conditional and not a measurement.} Nothing here establishes that capacity
> drifts, no rate is measured, and the kinetic argument for constancy in
> Section~\ref{sec:intrinsic} is untouched. What it establishes is that the specification as
> written cannot absorb a drift it does not currently exclude, and that a study running for a
> decade should record whatever would let $r$ be estimated afterwards rather than assumed.
> \end{remark}

**Macros this passage would need**, none of which exist yet and all of which have generators:
`tier1_n` (118, `src/design/run.py`), and, if the numbers are wired rather than left as prose,
the three from `audit/u7_drift.py`'s `n = 118` table: sd 0.2006, bias −0.0241, rmse 0.2020 at
`r = 0.005`.

---

## What the author decides

1. **Whether the passage lands at all**, which is the prose freeze.
2. **Where**: beside the Tier 1 requirement, where the failure bites, or in `11_scope.tex` with the
   other scope conditions, where it reads as one more limit on reach.
3. **Whether the three numbers are wired or left as prose.** They are single-arm figures from one
   table rather than headline results, and R143's rule applies: a macro is for a number that could
   move.
