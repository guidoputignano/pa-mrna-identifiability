# The genotype-prior document: first verification findings, and what 118 actually rests on

*Two questions from the author, answered before the R4 queue runs. Nothing is imported from the new
document into the repository by this file.*

---

## 1. The first verification finding against the new document, and it is a repeat

**The document states:** *"test–retest CV 6.9% (controls)/12.8% (patients)"* for Manoli 2021.

**That is the abstract's figure and it is wrong.** `REVIEW_LOG.md` R4 retracted it: Manoli's Results
give **14.02 ± 12.86%**, range 0.39 to 50.84, and the abstract reports **12.86, the standard
deviation, as though it were the mean**. R59 added the argument that settles it without leaving the
paper: **the abstract's two figures come from different columns of the same table.** The controls'
6.9 is their mean (6.94); the patients' 12.8 is their SD (12.86). A reader assuming the abstract is
internally consistent takes both as means, which is exactly the error.

**This is the third time this repository has met this figure and the second time it has been imported
from outside.** `refs/reviews_2026-09/README.md` records the September 2026 pack repeating it, and
counts it among that pack's seven citation errors. `08_specification.tex` already prints the
corrected figure with the retraction stated, and `PHASE7_TABLE2_PROVENANCE.md` D10 records that
`params.py` and `tab_params.tex` still carry the retracted one.

**What this is evidence of, which is the author's reason for asking.** The document's Manoli
paragraph is otherwise accurate and careful: it correctly identifies the study as **MMA and not PA**,
correctly notes the enzyme is **MUT and not PCC**, gives the right n (57 + 16, 83 tests), and draws
the right conclusion that no PA study pairs in-vitro activity with in-vivo flux. **The error is
narrow and it is a copy from the abstract rather than the Results.** That is a specific failure mode,
not general unreliability: where the document quotes an abstract it inherits the abstract's errors.
**Every figure it gives should be checked against the source's Results section**, and its qualitative
findings of absence, which are its main contribution, are not touched by this.

---

## 2. What 118 rests on, which is the load-bearing question

**Why it matters.** `u30_wire_drift.py`'s footing assertion established that Tier 1 needs 118 draws
**because `kappa_cv = 0.1402`**. Remove that and the count lands at a different number. So the
headline design figure is sized by Manoli's increment CV, and its verification status is
load-bearing.

**The repository reports it two ways and the author asked which is right.**

### The answer is neither branch the author offered, and it is better than one and worse than the other

**SYNTHESIS is NOT stale, and it is right.** `SYNTHESIS.md:122` and `:255` list the 14.02% as
**unverified**, and `REVIEW_LOG.md` **R59 says so explicitly and in terms**:

> *"Conditional, and recorded as such. The argument rests on two quotations this repository has not
> verified, R4's Results sentence and its controls figure, because outbound access to publisher and
> index records is denied in this container. It must be confirmed at the source."*

**So this repository has never read Manoli's Results.** It cannot have: there is no egress, and an
unbounded search of the whole repository finds no supplied copy of the paper. The primary sources the
author did supply (`k1_sources.py`, Phase 5 step 0) cover `LAMBDA` and `D_AMP`, **not the breath
test**.

**But it is not "the same weak footing" either, and that is the part worth having.** The figure has
**two independent arrivals**:

| | source | what it carries |
|---|---|---|
| **R4 / `PHASE2D_CHECKS.md`** | earlier phase | three verbatim Results passages, the subject-level causes (eGFR 57→43, VCO2 146.0→159.2), and a check of the published Correction's DOI |
| **`refs/reviews_2026-09/05`** | supplied September 2026, **later** | the same Results sentence with section pointers, the controls figure 6.94 ± 5.3%, the Correction's DOI and PMID, and the discrepancy **independently flagged** |

The repository's own header on that review file records the judgement: *"This is the one document in
the set whose retrieval evidently worked: it carries verbatim quotations and section pointers, and it
**independently found** the Manoli abstract-versus-Results discrepancy that R4 had already
retracted."*

**The internal argument is SUGGESTIVE, NOT ARITHMETIC, and an earlier version of this file
overstated it.** R59 and I both wrote that 12.8 is "a rounding of 12.86". **It is not: 12.86 to one
decimal is 12.9.** So the abstract's 12.8 is a truncation, or a differently computed figure, and not
a correct rounding of the Results SD.

**The direction still looks right**, because 12.8 is 0.06 from 12.86 and 1.22 from 14.02, and the
SD-as-mean story explains the confusion where nothing explains 12.8 as a mean. **But it is a
plausibility argument rather than an identity**, so it does not settle the question on its own.
**That makes the retrieval more necessary, not less.**

### So, stated precisely

**The 14.02% is CORROBORATED BUT NOT PRIMARY-VERIFIED.** Two independent secondaries agree on the
verbatim Results text, one of them with section pointers and a Correction check, and the retraction's
direction is internally checkable. That is stronger than "one unverified quotation" and weaker than
"read at primary". **SYNTHESIS's "unverified" is the correct label** and should not be relaxed.

**And where R4 got its quotations is not recorded anywhere**, which is the genuine gap. R4 predates
the September pack, the container has no egress, and no supplied copy exists. **A quotation presented
as *"I take the primary source over the repo's register"* with no record of how it was obtained is
exactly the provenance failure this audit exists to catch**, and it is in the audit's own files rather
than the paper's.

---

## 3. What follows

1. **`kappa_cv = 0.1402`, and therefore 118, should carry its status.** Not "unverified" in a way that
   implies doubt about the number, and not silence. The accurate sentence is that two independent
   secondary readings agree on the Results figure and no primary read has been made in this
   repository.
2. **This is the first item for the R4 queue when egress exists**, ahead of the rest, because it sizes
   the specification's headline count. One retrieval of Manoli's Results settles it.
3. **The new document's Manoli figure must not be imported**, and the pattern it reveals — abstract
   quoted where Results were needed — should be applied to every other number in it.
4. **`params.py:133-134` and `tab_params.tex` still print the retracted 12.8%** (`PHASE7` D10). That
   is a live inconsistency inside the repository, independent of the new document, and it is now the
   third record of it.
