# Phase 5. Robustness to model structure

**No paper edits, no merge.** Step 0 is a gate; this file currently contains Step 0 only.

---

## Step 0 (GATE). What the repository holds for the two published models

Everything below was found in the repository. Nothing was fetched and nothing is inferred
from memory.

### `pkpd_translational` — Attarwala, Lumley, Liang, Ivaturi, Senn (2023)

*Translational Pharmacokinetic/Pharmacodynamic Model for mRNA-3927, an Investigational
Therapeutic for the Treatment of Propionic Acidemia*, Nucleic Acid Therapeutics 33(2):141-147,
`doi:10.1089/nat.2022.0036`. Full CrossRef metadata in `paper/refs.bib:71`.

**Text held: two verbatim sentences**, in `paper/claims.csv` C46:

> "PCC protein expression rate was modeled as a linear function of plasma PCCA/B mRNA
> concentration using a 2-compartment indirect response PD model with an empirical effect
> compartment"
>
> "PCC protein-mediated reduction of 2-MC, 3-HP, and C3/C2 ratios was modeled using a direct
> sigmoidal maximum inhibition (Imax) model"

**Numbers held**, all in `paper/refs.sources:104-112`:

| quantity | value | unit as recorded |
|---|---|---|
| PCC-protein degradation rate `kdeg` | 0.0075 | /h (half-life 3.85 d) |
| plasma 3-HP disease baseline | 60.1 (14.4% RSE) | umol/L |
| between-subject variability on that baseline | 51.1 | % |
| fully corrected floor `E0` | ~0.0987 | umol/L |
| `IC50` for 3-HP | 37.5 | ug/g PCC |
| "slope" | 55.6 | **no unit recorded** |

**One structural statement about the parameter table**, C46's note: Table 1 "carries `E0`,
`base`, `IC50` and a fixed `Imax` and no biomarker turnover constant at all."

**PROVENANCE WARNING, the repository's own** (`paper/refs.sources:113-116`):

> "these figures are transcribed from the author's parameter-grounding review of September
> 2026, **not read off the paper in this repository's own pass**. The citekey is confirmed to
> be the right article; **the numbers await a check against the primary source**."

**What is missing.**

1. **The equations.** Not one differential or algebraic equation is held. The two sentences
   name a structure; they do not write it down.
2. **The PK half.** No mRNA/LNP disposition model, no two-compartment parameters, and no
   definition of the "empirical effect compartment" or its transfer to PCC protein.
3. **The `Imax` value.** C46 says it is fixed; the value is not recorded anywhere.
4. **What "slope 55.6" is.** No unit, and three readings are open: a Hill exponent on the
   sigmoidal `Imax`; the linear slope from plasma mRNA to PCC expression rate; or a
   steepness parameter on the effect compartment. A Hill exponent of 55.6 would be a step
   function, so the reading changes the model qualitatively. **This is the single most
   blocking gap.**
5. **Units on `IC50`.** "37.5 ug/g PCC" does not say ug of what per g of what, nor in which
   compartment, so it cannot be put into the paper's units.
6. **Whether the three markers share parameters** or each carries its own set, and whether
   `base` is per-marker.
7. **The residual-error model** and whether the 51.1% is on the baseline only.

### `dose_selection_poster` — Gao, He, Li, Liang, Luo, Natenshon, Patel, Real, Van, Ward, Zhang

*Model-Informed Dose Selection for mRNA-3927 in Participants with Propionic Acidemia (PA)*,
International Society of Pharmacometrics, `doi:10.70534/JVQU8233`. Full CrossRef metadata in
`paper/refs.bib:73`. **It is a conference abstract.**

**Text held: one verbatim sentence**, in `paper/claims.csv` C46:

> "An indirect inhibition Kin model was used to describe the PD response of 3-HP and 2-MC."

**Facts held** (`paper/refs.sources:119-121`, and claims C17 / C18): 0.6 mg/kg biweekly,
selected by simulation in 1000 virtual participants aged 0-36 y over body-weight
distributions, with 3-HP and 2-MC as the response; body weight is a covariate; **genotype is
not**.

**What is missing: everything structural beyond the model-class name.** The repository says so
itself, in C46's note:

> "it carries no PD equations and no biomarker parameters"

So there is no `kin`, no `kout`, no inhibition function, no `Imax` or `IC50`, no baseline, no
units, and no random-effect structure.

---

## What this means for Steps 1 to 4

**Steps 3 and 4 do not need the missing material.** Identifiability of a structure, and the
theorem about when a quasi-steady rational balance becomes linear in a reparametrisation, are
questions about the **form** of the equations. The forms are named unambiguously enough to
write down: (A) turnover with inhibited production and first-order loss; (B) indirect response
on an enzyme state with an algebraic saturating map to the marker.

**Step 2 does need it**, or an explicit substitute. Generating data from A and B and fitting
our model requires numbers. With what is held I can only supply them myself, and that would
make the misspecification test a test of **my** choice of parameters rather than of the
published models.

**One further gap, not about the references.** `StructuralIdentifiability.jl` is **not
installed in this container** (`julia` is absent; `symbolic/julia/` holds the two scripts from
earlier phases but they cannot be run here). Step 3's transient case can be done in `sympy` by
resultants and Groebner elimination, which is what Phase 0 used and what `f5`/`f6` already do,
but it will not be the same tool the paper cites. Say if you want Julia installed.

---

## The ask

**Stopping here, as the gate requires.** To proceed properly I need:

1. **`pkpd_translational` full text**, and above all **the model-equations section and Table 1**
   with units. The blocking item is what "slope 55.6" denotes.
2. **`dose_selection_poster`**: the poster itself if anything exists beyond the abstract. If
   the abstract is all there is, say so and I will treat (A) as *a named structure with no
   published parameters*.

**An alternative, if you would rather not wait.** I can run Steps 1, 3 and 4 now on the two
structures as named, and run Step 2 with A and B **calibrated to this paper's own operating
point** rather than to published values, flagged in every table as my parameterisation and not
theirs. That still answers the question Step 2 asks, which is whether the design rules and the
`theta` estimate survive a change of structure. It does not answer whether they survive *the
published models*, and the report would say so in those words.

Tell me which, and I will proceed.
