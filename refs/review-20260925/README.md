# Literature reports supplied 25 September 2026, stored verbatim

**Nothing in this directory is a finding of this repository.** These are six reports commissioned
by the author and committed here **exactly as received**, under the opaque filenames they arrived
with, including their own errors, their own confidence labels and their own caveats. They are raw
input.

**Do not cite anything from these files.** They are not sources. They are summaries of sources.
A number moves from here into `paper/`, into `src/common/params.py` or into a conclusion only
after its primary source has been opened and the value read in it, which is the same standard
`refs/literature_pending.md` states and `refs/reviews_2026-09/` was held to. **No constant has
moved on any of this, and none may.**

**They were not available when Phase 8 ran.** `audit/PHASE8_SIMULATIONS.md` carries four of the
questions these reports address as swept unknowns, named there as briefs R1 to R4, because no
document answered them. Where a report replaces a sweep with a number, that is a change to Phase
8's footing and not a confirmation of it. What Phase 8 asked for is in `refs/brief_questions.md`.

---

## The mapping

Nothing in the filenames says which report is which, and nothing in the repository recorded it
before this file. The mapping below was checked by reading each file's own title and scope, not
taken on trust.

| file | label | subject | its own stated target path |
|---|---|---|---|
| `compass_artifact_wf-b9d7ba6c-...md` | **R1** | the non-enzymatic clearance `k_e`: can it be pinned down well enough to identify θ | `01_ke_clearance.md` |
| `compass_artifact_wf-8d9f61a7-...md` | **R2** | the increment time course Δ(t) | `refs/reviews_2026-10/02_increment_time_course.md` |
| `compass_artifact_wf-f8c2cfcf-...md` | **R3** | markers, within-patient variability, and the case for a second marker | `refs/reviews_2026-10/03_markers.md` |
| `compass_artifact_wf-71022447-...md` | **R4** | the ex vivo assay as an anchor on θ | `refs/reviews_2026-10/04_ex_vivo_anchor.md` |
| `compass_artifact_wf-74a958c0-...md` | **R5** | the extrahepatic floor of the marker, anchoring `E0` | `refs/reviews_2026-10/05_extrahepatic_floor.md` |
| `compass_artifact_wf-6f9f8eb2-...md` | **R6** | the endpoint model: can a 3-HP to decompensation rule be written | none stated |
| `ideas.md` | n/a | **not a report.** A one-line placeholder, contents `files of research` | n/a |

The target paths are the reports' own, under `refs/reviews_2026-10/`, a directory that does not
exist. The files live here instead, under the names they arrived with, because renaming a document
stored verbatim is an edit.

---

## Three process warnings, each verified here against the files themselves

**1. Three reports claim a branch and a pull request that do not exist.** R4's header names branch
`reviews/2026-10-ex-vivo-anchor` and its closing line says "Open the PR and stop". R5's header says
"Work performed on branch `review/extrahepatic-floor`; PR opened". R2's header says "Branch:
`review/increment-time-course`. Status: PR opened, awaiting merge." `git ls-remote --heads origin`
on 26 September 2026 returns exactly two branches, `main` and `claude/modest-shannon-7qmfq1`.
**No branch of any of those three names has ever existed on the remote, so nothing those reports
describe as committed was committed, and no pull request was opened.** Every claim they make about
the state of the repository is therefore a claim about work that was never done. Read their
findings on their own merits and treat their process statements as false.

**2. R3 states outright that it could not read this repository.** Its own process caveat, section
1: the repository "could not be cloned or read, so the specific files named in the task ... could
not be verified verbatim, and a git branch/PR could not be opened", and its reconstruction of the
pending-reference table "is a reasonable reconstruction and should be reconciled with the real
`refs/literature_pending.md` before commit". **So every statement R3 makes about what this
repository contains, assumes or already records is unverified, including anything it says about
`params.py`, `REVIEW_LOG.md` or the other reviews.** Its literature findings do not depend on that
access and are not tarred by it; its claims about the repository are.

**3. The September 2026 pack contained seven citation errors, and the same standard applies here.**
`refs/reviews_2026-09/README.md` records what checking the previous five reviews found: seven
citation errors, including a DOI resolving to an unrelated neuroscience paper, a design paper
attributed to the wrong author, and a retracted breath-test CV repeated a second time; several
claims were wrong against the primary source. Three of those five reviews had stated that their own
retrieval tools failed. **Nothing about this pack earns a lower bar.** These reports do label their
claims VERIFIED, SECONDARY or NOT FOUND, and R1 in particular states that species and population
are given "because the repository has previously mistaken a mouse estimate for a human one and an
SD for a mean", which is a good sign and not a substitute for checking.

A fourth thing, not a warning but worth recording: **R3 and R4 give conflicting readings of the
same figure**, the breath-test repeatability in Manoli 2021. R3 says the 12.8% is the mean of
per-subject CVs and that the ±12.86% elsewhere in the paper is the SD of that distribution. R4
quotes the 12.8% as an intrasubject CV without that distinction. `REVIEW_LOG.md` R4 retracted the
figure on the grounds that the abstract reported an SD as a mean. **That conflict is unresolved
here and is not adjudicated in this file.** It has to be settled at the source, because if the
repository's retraction is wrong then a retraction is itself retracted.

---

## What has and has not been done

**Done:** the six files are committed verbatim, the mapping is recorded, and the three process
claims above are checked against the files and against the remote.

**Not done:** no citation in any of the six has been resolved against its primary source, no
number has been carried into a conclusion, no constant has moved, and none of the decisions these
reports bear on has been taken. Which numbers are worth checking, and their verdicts, will go in
`audit/PHASE9_LITERATURE.md` when that work is authorised. Eighty-odd citations are cited across
the six reports and the plan is not to check them all: only the ones that would move a conclusion.
