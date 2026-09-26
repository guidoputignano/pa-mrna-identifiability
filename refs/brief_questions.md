# Questions for the next round of literature briefs

**Nothing in this file is a finding, a citation, or a verified fact.** It is a list of questions
that simulation work has raised and that no document in this repository answers. Each is filed
against the brief that would cover it, using the labels the Phase 8 task uses: **R1** for the
marker's non-enzymatic clearance `k_e` and its variability, **R2** for the increment time
course, **R3** for the marker's measurement CV, **R4** for the ex vivo mapping to θ.

**The rule for a lead.** A PMID, DOI or author named here has **not** been retrieved, opened or
checked in this repository. It is a pointer for whoever runs the brief, not a reference. It does
not go into `paper/refs.sources` or any prose until it has been resolved against a publisher
record, which is the same standard `refs/literature_pending.md` states.

**Do not edit `refs/reviews_2026-09/`.** Those documents are stored verbatim as received. New
questions live here.

---

## R1. The marker's non-enzymatic clearance

Raised by `audit/PHASE8_SIMULATIONS.md` S2b, which prices the whole disease-range requirement on
a per-patient measurement of `k_e`.

1. **Is 3-HP itself a substrate of any step in the pathway, or only an end product?** The model
   writes the enzymatic term as a clearance of the marker. Taken literally, a trace-level tracer
   would be cleared enzymatically too, and S2b.4 prices that reading at a factor of 5.8 on
   `sd(ln θ̂)`: it is the difference between a study worth running and one that is not. The
   biochemistry, not the model, has to answer it.
   - **A lead, not registered and not verified here:** a rat liver study on the interrelations of
     3-hydroxypropionate and propionate metabolism, said to be relevant to disorders of
     propionyl-CoA metabolism, **PMID 28634175**. Retrieve and read it; do not cite it from this
     entry.
   - **What would settle it experimentally**, inside the micro-dose study the specification
     already proposes: the fraction of label recovered unchanged against the fraction appearing
     in downstream products.
2. **Within-patient stability of the marker's clearance over months to a year.** S2b.5 makes
   this the binding constraint below 10% of normal, and it is swept at 0, 5 and 10% because
   nothing measures it. A repeatability study, or serial clearance measurements in the same
   patients, would replace the sweep.
3. **Does that drift have a covariate?** If it tracks kidney function, which is measured
   routinely and declines in this disease, a covariate replaces a repeated study. Anything on
   marker clearance against eGFR, measured GFR, or renal tubular secretion would decide it.
4. **Still open from S2b.0 to S2b.3, and unchanged:** a human plasma 3-HP half-life, clearance,
   volume of distribution, and the renal fraction of its elimination. The half-life is unknown by
   a factor of about 70 across the sweep, and it is what the tracer study's sampling window must
   be chosen against.

---

## R4. The ex vivo mapping to θ

Raised by `audit/PHASE8_SIMULATIONS.md` S5, which carries a direct enzyme assay on the patient's
own cells as a log-normal prior on θ and finds its value and its danger in the same cells.

1. **How wide is the mapping** from an ex vivo activity, per milligram of cell protein in
   fibroblasts or lymphocytes under assay conditions, to the whole-body maximal velocity the model
   calls θ? S5 sweeps sd(ln θ) from 0.3 to 1.0 and finds that **only 0.3 changes any verdict**, so
   the question has a sharp answer-or-nothing quality: a mapping good to about a third is worth
   one cell of the requirement map, and anything wider is worth none.
2. **Is the mapping CENTRED?** S5.3 makes this the more important of the two, and it is the one
   nobody asks. At 5% of normal with a 0.3 prior, 54 to 60% of any centring error passes straight
   into θ̂, and at 10% a prior wrong by a factor of two is worth exactly nothing (RMSE 0.249 with
   it and 0.249 without). A width without a centring claim cannot be used in an estimator.
3. **What would establish the centring**, rather than the width: paired measurements in the same
   patients, ex vivo activity against any in vivo flux measure (a labelled propionate oxidation
   study is the obvious candidate), across more than one cell type.
