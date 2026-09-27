# Phase 11, round 3. The register understates its own requirement, and the one route to `k_e` whose obstacle is cost

*Proposals only. Nothing under `paper/` is changed by this file.*

Scripts: `audit/u13_ke_requirement_audit.py`, `audit/u14_dynamic_burst.py`.

---

## 1. `params.py`'s `k_e` row understates its own requirement by 7 to 17 times

### 1.1 The claim, and where it lives

`src/common/params.py:206`, the `k_e` row of the register **that generates the paper's
parameter table**, says of the parameter it itself calls "the priority gap":

> Section~\ref{sec:designsec} needs it within a factor of about three, and no independent
> measurement in humans supplies it.

`SUMMARY.md` carries the same figure in its specification row.

### 1.2 The number is right, and it is right at one capacity

**The owner was located before recomputing this time.** `src/design/run.py` `spec_numbers()`
solves `brentq(sd_of(monthly, ke_prior=x, seed=s) - 0.20)` over seeds. `monthly` is the
25-draw record and `sd_of` runs at the module default `P.theta = 3000`. So the published
figure is **the nominal-capacity requirement on the recommended record**.

| | log sd | as a factor |
|---|---|---|
| `spec_numbers()`, the owner | **0.9803** (0.9761 to 0.9865) | **2.665** |
| `audit/u9_ke_frontier.py`, same capacity and record | **0.9815** | 2.669 |
| difference | 0.0012 | **0.12%** |

Two code paths, one quantity. `07_design.tex:60` is already wired to the macro and prints
2.67 with its seed range; the register's "about three" is that value rounded, well inside
the repository's own convention. **Nothing here is a wrong number.**

### 1.3 But the requirement has a capacity axis and is stated without one

Tier 1 record, the arm that needs a decade:

| theta | % of nominal | frontier log sd | as a factor | vs "about three" | vs the generated 2.67 |
|---|---|---|---|---|---|
| 3000 | 100% | 0.8001 | 2.226 | 1.4x | 1.2x |
| **600** | **20%** | **0.1567** | **1.170** | **7.0x** | **6.3x** |
| **300** | **10%** | **0.0633** | **1.065** | **17.4x** | **15.5x** |
| 150 | 5% | *unreachable at any width* | | | |

> **The correction.** Not that 2.67 is wrong, but that it is stated as *the* requirement for
> the parameter the paper calls its priority gap, with no capacity qualifier, in the row a
> reader consults to cost the measurement. At the capacities where the clinical question
> lives the same quantity is **7 to 17 times tighter**.

**The reader this misleads is the one costing the measurement.** "Within a factor of about
three" reads as a modest ask that a literature bound or an analogue might satisfy. A factor
of 1.07 does not read that way at all.

**And one row below the live range the axis changes character.** At a twentieth of nominal
there is no frontier, because the increment's own prior floors the design above the gate
whatever `k_e` is. So the requirement is a *measurement* one between a fifth and a tenth of
nominal and a *structural* limit below that, and a single unqualified factor cannot express
either. Proposed wording is in `audit/PHASE11_PROPOSED_TEXT.md` §5.

---

## 2. The `k_e` measurement question, reopened

### 2.1 What `u9` closed, and what it did not

`audit/u9` closed the **timed urine collection** on three counts: units, channel coverage,
species. Those close a *collection*. They do not close the *object*.

**`params.py`'s own `k_e` row says `k_e = 5/day` is "a marker half-life of `ln2/5 = 0.139 d`,
about 3.3 hours".** So `k_e` is a **relaxation rate**, and a relaxation rate is read off a
concentration-time curve. No volume of distribution, no collection completeness, no species
extrapolation. **The units objection does not apply to it at all**, because nothing is
divided by a volume.

### 2.2 The quasi-steady record cannot see it, at any density

Confirmed rather than asserted. The same 118 samples reallocated from monthly to three-hourly:

| theta | monthly | weekly | daily | 3-hourly |
|---|---|---|---|---|
| 3000, `sd(ln theta)` | 0.1543 | 0.1522 | 0.1540 | 0.1540 |
| 3000, `sd(ln k_e)` | 0.7971 | 0.7988 | 0.7971 | **0.7968** |
| 300, `sd(ln k_e)` | 0.5425 | 0.7723 | 0.5532 | **0.5221** |

`eq:local` is a balance with no state, so it has integrated the relaxation away and answers
no by construction. **The question has to go to the dynamic model of Section 3 or nowhere.**

### 2.3 The obvious burst fails, and the diagnostic that looked favourable is the wrong one

`k_e`'s share of the **flux** is 15.5% / 52.8% / 70.8% / 83.7% as capacity falls, so the
measurement should get easier exactly where the frontier tightens. **A relaxation
measurement does not depend on the flux share.** It depends on `k_e`'s share of the
**relaxation rate** `J = (theta + Delta) g'(m) + k_e`, and the dose adds the very enzyme
term `k_e` competes with in `J`:

| theta | at the trough: `k_e/J` | just after a dose: `k_e/J` | flux share |
|---|---|---|---|
| 3000 | 0.149 | 0.054 | 0.154 |
| 600 | 0.438 | 0.072 | 0.510 |
| **300** | **0.552** | **0.076** | 0.686 |
| 150 | 0.626 | 0.077 | 0.819 |

**The dose suppresses `k_e`'s share of the relaxation rate by a factor of seven** while the
flux share barely moves. And the cost lands exactly there: at `theta = 300`, **sixteen
inpatient days each spanning a dose move `sd(ln k_e)` from 0.5430 to 0.5277**, about two per
cent. *The route is not weak; the placement is wrong.*

### 2.4 Perturb at the trough instead

One 24 h hourly burst, Tier 1 record, a prescribed load step at the trough:

| theta | no burst | fast −50% | fast −30% | load +30% | load +50% | post-dose | frontier |
|---|---|---|---|---|---|---|---|
| 600 | 0.6565 | 0.5046 | 0.5701 | 0.5251 | 0.4072 | 0.6546 | 0.1567 |
| **300** | **0.5430** | 0.3685 | 0.4380 | 0.3869 | **0.2846** | 0.5399 | **0.0633** |
| 150 | 0.4721 | — | — | 0.3204 | 0.2339 | 0.4678 | n/a |

**One load step at the trough is worth what sixteen post-dose days are not.**

### 2.5 The tolerable level caps the useful amplitude

The direction that helps most raises the marker, and `src/common/params.py` sets
`m_tol = 40` uM. The largest prescribed step whose whole burst stays under it:

| theta | trough m | max safe step | peak m | `sd(ln k_e)` | no burst | frontier |
|---|---|---|---|---|---|---|
| 600 | 24.80 | **+40%** | 38.44 | 0.4635 | 0.6565 | 0.1567 |
| 300 | 31.23 | **+15%** | 38.72 | 0.4798 | 0.5430 | 0.0633 |
| 150 | 35.64 | **+5%** | 39.96 | 0.4546 | 0.4721 | n/a |

**The amplitude that would close the gap is the one the clinic exists to prevent, and the
amplitude that respects the paper's own tolerable level buys almost nothing.** At a tenth of
nominal the safe raising step takes 0.5430 to 0.4798.

What is left is the other direction. **A fast lowers the marker and is safe.**

### 2.6 The safe route, and what it actually reaches

A 50% fast at the trough, 24 hourly draws per admission. The burst's assay CV is separated
from the record's because an inpatient series drawn from one line and analysed in one run is
where a tighter assay is plausible if it is plausible anywhere.

**theta = 600, frontier 0.1567**

| burst assay CV | 1 admission | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| 5% | 0.5046 | 0.4456 | 0.3874 | 0.3377 | 0.2993 |
| 2% | 0.3539 | 0.3121 | 0.2780 | 0.2481 | 0.2192 |
| **1%** | 0.2780 | 0.2481 | 0.2192 | 0.1922 | **0.1683** |

**theta = 300, frontier 0.0633**

| burst assay CV | 1 admission | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| 5% | 0.3685 | 0.3181 | 0.2751 | 0.2416 | 0.2158 |
| 2% | 0.2523 | 0.2246 | 0.2005 | 0.1769 | 0.1534 |
| **1%** | 0.2005 | 0.1769 | 0.1534 | 0.1313 | **0.1121** |

**The route is assay-limited, not physics-limited.** That is the favourable half of the
answer: the single input that moves it most is the one an inpatient setting is best placed
to improve. At a 1% burst assay one admission already beats the frontier at `theta = 600`
when the step is a +50% load rather than a fast, and the safe version needs sixteen.

---

## 3. The three routes, priced against one frontier

| route | verdict | obstacle |
|---|---|---|
| timed urine collection | **blocked** | not precision: a `V_d` for 3-HP and the renal share of non-PCC disposal, each to a few per cent (R107) |
| prescribed carnitine modulation | **open but weak** | at the physical ceiling, `sd(ln k_e) = 0.469` against 0.063 (R104) |
| **sub-daily burst on a load step** | **open, short by 1.1 to 1.8** | a 1% assay and sixteen inpatient days |

**The burst is the only one of the three whose obstruction is a cost rather than an
unmeasured quantity or a dose nobody would give**, and the only one for which the units
objection never arises, because a relaxation rate is read off a curve and never divided by a
volume. That is a materially different answer from the one the collection got, and it is the
reason the question was worth reopening.

**It still does not reach the frontier.** Best safe combination — a 50% fast, a 1% burst
assay, sixteen inpatient days — gives **0.1683 against 0.1567 at a fifth of nominal** and
**0.1121 against 0.0633 at a tenth**.

---

## 4. What is not settled

- **The prescribed load step is treated as exactly known.** A delivered step with its own
  error degrades everything in §2.4 onward, and nothing here bounds that error.
- **Safety is assumed, not established**, as `t3_diet_modulation.py` assumes it for diet.
  §2.5 uses the paper's own `m_tol` as the constraint, which is a convention and not a
  clinical limit; a real protocol would need a clinical one.
- **A 1% marker assay is asserted as plausible in an inpatient series and not sourced.** It
  is the input the conclusion is most sensitive to, so it joins the gate list.
- **The burst inherits the observation-model assumption** the paper names as its first
  limit, and inherits it through the relaxation rather than the level.
- Nothing here is checked against a primary source; the container has no egress.
