# Identifiability of residual enzymatic capacity in propionic acidemia

A control-theoretic study of individualised dosing for an mRNA enzyme-replacement
therapy in a monogenic metabolic disease.

**The clinical decision** is a dosing interval. **The quantity that determines it** is the
patient's residual enzymatic capacity `theta`, which is never measured. **The question** is
how much of it can be known, from data it is ethical to collect, before the decision has to
be made.

---

## Start here

| if you want | read |
|---|---|
| the result, on one page | **[`SUMMARY.md`](SUMMARY.md)** — four claims, the specification, the coverage qualification, what was withdrawn |
| the argument and the proofs | `paper/main.pdf`, 27 pages (`make paper`) |
| the background, derivations, validation protocol and open questions | `paper/supplement.pdf`, 47 pages, numbered S1, S2, ... and cited from the paper by section number |
| what this work got wrong, and why | **[`REVIEW_LOG.md`](REVIEW_LOG.md)** — nine retractions, each kept with its cause |
| how the repository is laid out | the rest of this file |

**This README was itself rewritten against the current state.** An earlier version
presented as findings several results the work has since withdrawn. If you have a cached
copy, discard it and read `SUMMARY.md`.

---

## The model, at three resolutions

One model. The letters A, B and C survive only as directory names: two of the three were
the same model read at different resolutions, and the third was a hypothesis about one term.

```
dm/dt = L(t) + L0 - (theta + Delta(t)) * g(m) - alpha(m),    g(m) = m/(Km + m)
```

`g` is the **Michaelis-Menten saturation** and `theta` **absorbs** the maximal velocity: it
*is* the patient's residual `V_max`, not a fraction of a population one. `Delta(t)` is the
therapeutic increment, known from the dosing schedule. `L0` is the **unrecorded** part of
the propiogenic load, carried as an unknown.

| resolution | obtained by | what it is for |
|---|---|---|
| **the dynamic model** | as written | the trajectory, the excitation condition |
| **the equilibrium reduction** | `dm/dt = 0` | the fold, `theta_crit`, the safety margin |
| **the local identification model** | `alpha(m) -> ke*m`, over the visited band | the fibre, the margin, the design |

**The third is not a reduction of the first**, and the check is one line. The tangent
intercept of `theta*g + alpha` at the operating point is `a = theta*g^2 + [alpha - m*alpha']`,
so the local model's identity `a = theta*g(m_bar)^2` holds only when the alternative channel
is a straight line through the origin. At the illustrative constants the second bracket is
374 against the enzymatic term's 61. The local model is therefore the *parent with one term
replaced*, and it is retained because that replacement is the **most favourable** case.

A fourth structure, `alpha(m) = ke*h(m)` with a ceiling instead of a fold, survives as one
remark. It is the null against which the fold is testable, and it prices biological fidelity
twice over: a factor of about 450 in condition number, and — because both channels then share
a functional form — an *exact* two-point degeneracy in which capacity and alternative
clearance exchange. That is a scope limit on this work rather than support for it.

---

## What the note claims

The four claims in full, with their numbers, are in [`SUMMARY.md`](SUMMARY.md). In brief:

1. **`theta` is identifiable; the obstruction is a sampling regime, not a structure; and the
   therapy is what removes it.** The model is globally structurally identifiable, so the
   information is in the system — what a monthly record loses is the transient. On the
   quasi-steady locus that remains, the obstruction is a rank-one identity whose null
   direction is capacity, affinity and the unrecorded load, with the alternative clearance
   provably absent from it. The excitation that closes it is the dose train the trial already
   administers, at **42x** the marker band's contribution.
2. **Identifying `theta` does not decide the dosing interval.** The threshold it would be
   compared against is stipulated, not measured: across its defensible range it moves by a
   factor of **24** and the sign flips for the same patient with perfect data.
3. **The threshold framing is itself doubtful.** The only endpoint evidence available reports
   a smooth dose-response and no threshold. Recorded as a finding, not adopted as a
   replacement estimand.
4. **Capacity stops being recoverable where it matters most** — above roughly three to four
   times the managed marker level, which is the severely affected population.

**Two results a reader should not miss**, both in `SUMMARY.md`:

- **The therapy is the design's condition of existence, not an improvement on it.**
  Spontaneous load variability alone gives 1.77 against a 0.20 target, and adding carnitine
  perturbations makes it *worse*. There is no second specification for untreated patients.
- **The stated precision is not delivered to every patient.** Median 0.192 across a
  heterogeneous cohort, with just over half meeting the target, and the shortfall is
  systematic rather than random.

Supporting results that survive unchanged: with `Km` unknown — which is this work's own
situation — the two input channels stop being informationally substitutable by six orders of
magnitude in `det G` and by a factor of 46 to 807 in `sd(ln theta)`; and the assay term
dominates the load term whenever the load is recorded better than about its own CV, the two
balancing at `CV = 1.08 x (load error)`.

---

## Layout

```
paper/          LaTeX source, split by section
  main.tex        build target
  sections/       one file per section
  note_structural_identifiability.tex   standalone companion note
src/
  common/models.py      the kinetics and the three clearance structures
  common/params.py      every constant, with a source column
  common/dynamics.py    the sensitivity system and the observability Gramian
  common/excitation.py  phi, E, E_max, and the design quantities
  model_C/run.py  the dynamic model and its equilibrium reduction
  model_A/run.py  the local identification model: the fibre, the margin, the floor
  model_B/run.py  the ceiling hypothesis (one remark in the paper)
  design/run.py   the excitation figure, the covariance collapse, the decision
                  characteristic, the excursion table, the data specification
symbolic/       structural identifiability, three independent tools
  julia/          StructuralIdentifiability.jl
  maple/          SIAN
  daisy/          DAISY
  sympy/          observational equivalence, worked by hand alongside the tools
audit/          the check scripts behind every gated decision, one file per check
tools/
  fetch_refs.py     builds paper/refs.bib from CrossRef and PubMed
  claims_to_tex.py  renders paper/claims.csv as the claim register appendix
results/
  figures/        every figure in the paper, regenerated by `make figures`
  tables/         LaTeX tables generated from code
tests/          consistency checks asserting the claims in the paper
```

No figure in the paper is drawn by hand. Everything under `results/` is regenerated from
`src/`.

---

## Running

```bash
pip install -r requirements.txt

make figures     # regenerate all figures and tables (seconds)
make refs        # fetch paper/refs.bib, render the claim register
make paper       # build paper/main.pdf (runs refs first)
make test        # assert the structural claims still hold
make audit       # claim register and citation integrity
make all
```

`make test` passes 76 tests; `make paper` builds 73 pages; `make audit` reports 46 claims,
none contradicted, none unsourced, four matches-with-drift each carrying a note saying what
could not be confirmed. See `REFERENCES.md`: no citation is typed by hand, and no factual
claim enters the paper without a row in `paper/claims.csv` saying whether the source actually
says it.

Structural identifiability needs Julia and is slower:

```bash
julia -e 'using Pkg; Pkg.add("StructuralIdentifiability")'
julia symbolic/julia/run_identifiability.jl
julia symbolic/julia/run_channel_exchange.jl
```

**These scripts record what the tool returned, not what was expected of it.** An earlier
version of this file, and of the scripts themselves, wrote the expected answer beside each
model; two of those expectations were wrong, and each would have read as a confirmation of
itself. The rule the repository now follows is that no script states an expected result in a
comment: an expected value goes in an assertion that can fail, or nowhere.

---

## Known issues

**`make audit` is not deterministic, and this is unfixed by choice.** It resolves every DOI
and PMID over the network on each run, so an external service's availability sits inside a
required gate. A transient failure on a single lookup fails the whole check even though the
reference is already cached — observed once in CI, where the same commit passed in one run
and failed in another three seconds apart. The fix is to fall back to the cached entry on a
transient fetch failure and fail only on genuinely unresolvable references. It is deferred
to after review rather than done now, because the repository is frozen for reading.

---

## Open items

The paper's own list, in its discussion, of what would change the answer.

1. **A same-day paired repeatability figure for the breath test.** It does not exist, and it
   is what Tier 1 of the specification asks for. Everything the therapy channel buys is
   conditional on the delivered increment being characterisable, and only a between-visit
   upper bound is published. **Until it reports, no protocol can be written.**
2. **An endpoint model that does not require a threshold.** The natural history study's own
   hazard formulation is the shape to build on. Doing so properly needs the reference
   concentration and the model specification that publication does not report. This is the
   route by which the worst-illness amplitude leaves the decision boundary the design cannot
   resolve.
3. **A marker that is the substrate, or an assay of the substrate.** That would remove the
   observation-model assumption on which the whole construction rests, and it is the single
   change that would most enlarge what a record can support.

---

## Caveats a reader should hold

Parameter values are illustrative, not fitted. Their role is to fix the qualitative structure
and the orders of magnitude, and the design results are reported relative to the patient's own
operating level.

`m` in the model is the measured plasma marker; the enzyme's actual substrate is
intracellular propionyl-CoA, which is never measured. The map between them is an assumption,
and it is the first thing the discussion names as a limit on what any of this licenses.

`k_e` is treated as a population parameter on the grounds that it does not involve PCC, but
carnitine is also administered therapeutically at patient-specific doses. That dose is
adjusted in response to the patient's state, so it is endogenous and using it as a
perturbation would need an instrument — where the dose train is fixed by protocol and
exogenous by construction.

The depletion mechanism in the dynamic model is a hypothesis. A plain two-state depletion
model gives saturation, not decline; decline requires that replenishment itself be impaired
as the metabolite accumulates. The check is whether free carnitine collapses during
decompensation, which is measured routinely.

## Licence

The code in this repository is released under the MIT licence; see [LICENSE](LICENSE).
