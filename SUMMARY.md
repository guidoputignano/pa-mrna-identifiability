# One page, for a reader who has not followed the work

*Assembled entirely from the paper (`paper/main.pdf`; derivations and background are in
`paper/supplement.pdf`). Every number here appears there with its derivation.
Frozen with the paper at `b289ef8`. For the technical overview of the repository see
`README.md` — but read the note at the foot of this page first.*

---

## The question

Propionic acidemia is treated, in trial, with an mRNA enzyme replacement. The clinical
decision it creates is a **dosing interval**, and that interval is set by a quantity nobody
measures: the patient's residual enzymatic capacity, written **θ**. This note asks a narrow
question — given a clinical record of the kind the disease actually generates, can θ be
identified, and does identifying it decide anything?

## The four claims

**1. θ is identifiable; the obstruction is a sampling regime, not a structure; and the
therapy is what removes it.** The model is globally structurally identifiable, so the
information is in the system. What a monthly record loses is the *transient* — the marker
settles in about 0.8 h against a therapeutic half-life of 3.5 days. On the quasi-steady
locus that remains, the obstruction is a rank-one identity whose null direction is capacity,
affinity and the unrecorded load, with the alternative clearance provably absent from it.
The excitation that closes it comes from the dose train the trial already administers, at
**42×** the marker band's contribution.

**2. Identifying θ does not decide the dosing interval.** The decision is usually posed as
the sign of a margin against a critical capacity. That threshold is not measured: it is
defined through two chosen constants, and across their defensible range it moves by a factor
of **24** and the sign flips for the same patient with perfect data.

**3. The threshold framing is itself doubtful.** The only endpoint evidence available
reports a smooth dose–response and no threshold at all. Recorded as a finding, not adopted
as a replacement estimand.

**4. Capacity stops being recoverable where it matters most.** Above roughly three to four
times the managed marker level, two patients differing tenfold in residual capacity are
indistinguishable — which is the severely affected population.

## The specification

| | |
|---|---|
| **Who** | Patients receiving a therapy that **varies** over the sampling window. The dose train is the design's *condition of existence*, not an improvement on it — see below. |
| **How many** | **14** to confirm the design achieves its stated precision (±20% relative); **96** to bound the fraction of patients it works for (±10 points). Fourteen is inside a phase 1/2 trial of sixteen. |
| **How many draws** | **Monthly draws over two years (n ≈ 25)** → sd(ln θ̂) = 0.170. Thirteen two-monthly draws also pass, at 0.187. Ordinary outpatient blood work. |
| **Which days** | **Not on dosing days, nor any multiple of the dose interval.** Every draw then lands at the same phase and the precision goes from 0.17 to **3.51**. Moving the schedule **half a day** recovers it. This is the one schedule a protocol would naturally write. |
| **The load record** | Ordinary precision — twenty-fold better buys six per cent. But **record intake as eaten and do not smooth it**: a perfectly steady load is the exactly singular case. This is an instruction to omit a data-cleaning step, not to add a measurement. |
| **Known beforehand** | The non-enzymatic clearance within a factor of about **three**, which the literature supplies. The delivered increment characterised to about **12%**. **Nothing is assumed about the Michaelis constant** — a prior twenty times looser leaves every number unchanged. |

**Tier 1 comes first.** One measurement nobody has published — the same-day repeatability of
a 1-¹³C-propionate breath test across a single dose. Until it reports, no protocol can be
written.

## The qualification that matters

**The 0.170 is one patient.** Across a plausibly heterogeneous cohort the median achieved
precision is **0.192** and **just over half** meet a 0.20 target. The design confirms
identifiability in a cohort; it does not deliver the stated precision to every patient in
one.

**And the shortfall is systematic, not random.** The strongest correlate of the achieved
precision is the increment over capacity. The dose is delivered in absolute terms; the
capacity it is measured against is not. So a patient with **more** residual capacity
receives less *relative* excitation and is identified less well — by arithmetic, not by
chance. θ̂ should always be reported with its own interval, never as a guarantee attached to
an individual record.

**Why the therapy cannot be substituted.** On the same footing, spontaneous load variability
alone gives **1.77**, and adding recorded carnitine perturbations makes it *worse* — carnitine
acts on the wrong channel and moves the patient deeper into the linear regime to buy a small
gain. Reaching the target without a dose train would need day-to-day intake varying by
**229%**, or daily venepuncture for two years, which still lands at 0.67.

## What was withdrawn

Nine claims were retracted during the work and are kept in `REVIEW_LOG.md` with their causes
rather than deleted. A reader of the commit history needs these:

- **The seven-sample single dose cycle does not pass.** It passed only with the alternative clearance held known.
- **A 52.5% decision figure is one seed**, not a finding (43.9% ± 3.5% over 24 seeds). *That retraction was itself retracted once*, on a check that sampled the wrong seeds.
- **A sign-error plateau, a correlation reversal and a breath-test repeatability figure** are all withdrawn — the last inherited from a source whose abstract contradicts its own Results.
- **"A constant increment leaves only a lumped sum identifiable" is wrong**; the increment is a *known* input. Caught by an external tool, not by re-running the project's own code.
- **A correction to a filter-exactness claim is itself withdrawn** — the original claim was right and the note already argued it.
- **Two published precision figures (0.054, 0.023) are withdrawn**: they reproduce only on a dense design with the alternative clearance held exactly known, a footing the note disowns, and they were computed with the dose train present in a passage about what load variation buys.
- **"Information in the degenerate direction saturates" is withdrawn.** It grows with sample size; it never *arrives*, because the prior holds that direction at every sample size a study could reach.

Three standing rules came out of these. No script states an expected result in a comment.
A stated mechanism is a claim carrying the same burden as a number — retractions included.
And every quoted precision carries its footing in the same breath: which design, which
priors, and whether the therapy channel was present.

---

### A caution about `README.md`

`README.md` was written before the corrections above and **has not been updated**. It still
presents as findings several numbers this page lists as withdrawn — among them the
seven-sample design, a 53% sign error, a 0.004 precision figure and its derived "factor of
122", a cohort of ~3,100, a sign-error plateau near 6%, and a 12.8% breath-test
repeatability. **Read the paper or this page in preference to it.** The discrepancy is
recorded here rather than repaired, because the document set was deliberately frozen.
