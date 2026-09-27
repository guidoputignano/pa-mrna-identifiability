# The supplement pass

*`paper/supplement.pdf`, 64 pages, 36,979 words, built at `7c0f521` with 0 errors and 0 unresolved
references. **Part 1 was written with the supplement open and nothing else**, per the author's
instruction, including `main.pdf`. Parts 2 and 3 follow, with the main document open.*

**Judged as a reference, not as a read-through**, per the author: can a reader who arrives at one
place, from one pointer, get what they came for without reading around it.

---

# Part 1. The cold read

## 1.1 What a reader consulting it would take away

**That it is two documents bound as one.** `S1` is titled *"material moved out of the main text"*,
and it contains **two parallel treatments of the same material**: a short set at `S1.1`–`S1.12`
(pages 3 to 12) and a long set at `S1.13`–`S1.20` (pages 13 to 52). The contents page shows the
duplication in its headings:

| the short set | the long set |
|---|---|
| `S1.5` What the dynamic model would have to establish (p6) | `S1.15.7` What would have to be established (p25) |
| `S1.6` The testable discriminant, and the observable it has to be stated in (p7) | `S1.15.8` Testable discriminant, and the observable it has to be stated in (p25) |
| `S1.7` What the accessible input costs (p8) | `S1.16.5` **What the accessible input costs** (p32) |
| `S1.8` The local reduction's known defect (p9) | `S1.16.8` Known defect (p37) |

Two of those headings are **verbatim identical**. I checked two pairs against their bodies and they
are the same material differently worded, not stubs: `S1.7` and `S1.16.5` both price the carnitine
route, both carry `0.170` and `1.91 to 2.02`, and they are numbered **Remark 4** and **Remark 17**.
`S1.6` and `S1.15.8` both give the fold-versus-ceiling discriminant, with the longer one adding the
free-carnitine restatement.

**Beyond that duplication the document is better made than the main text in two respects and worse
in one.** `S2` (provenance) and `S3` (reproducibility) are the strongest writing in the project.
`S1.20` (questions since answered) does, unprompted, the thing this entire audit has been enforcing:
it preserves superseded text verbatim and states item by item what replaced it. What is worse is
navigation, below.

## 1.2 What is wrong, garbled or stale

### (a) Ten pointers send the reader to the document they are already in. This is the largest defect

Counted with fixed strings over the whole text, no result limit:

| form | count |
|---|---|
| `"Appendix S1"` | **6** |
| `"supplement §S1"` | **3** |
| `"the supplement's local-model section"` | **1** |
| **total** | **10** |

**Every one occurs inside Appendix S1.** Examples, verbatim:

> *"`S1.15.7` What would have to be established. Four things, set out in **Appendix S1**"*
> *"Two qualifications, and **Appendix S1** carries them in full"*
> *"Three assumptions, stated here and warranted in **Appendix S1**"*
> *"Cost bounded in **supplement §S1**: over `Kc` from 25 to 800..."*

**For a reference reader these are dead ends dressed as live pointers.** Three of them sit in Table
S5's last column, which is the column that discharges every assumed row: a reader consulting one row
for its cost bound is told the bound is in "supplement §S1", which is the document in their hands,
forty pages earlier, with no subsection named.

**The likely cause is visible in S1's own title.** These sentences read correctly when they sat in
the main text and "Appendix S1" was a forward reference. The block was moved into the appendix and
the cross-references were not rewritten. It is the compression, not carelessness.

### (b) 58 cross-document references do not say they are cross-document

`"Section 3"` ×15, `"Section 6"` ×17, `"Section 4"` ×11, `"Section 7"` ×8, `"Table 9"` ×8,
`"Table 2"` ×4, `"Table 7"` ×3, `"Section 5"` ×3, `"Section 8"` ×1. Against that, the main text is
named explicitly **5 times** (`"the main text"` ×3, `"of the main"` ×2).

**None is broken** — nothing resolves to `??` or `0`, so the `xr` mechanism works. The defect is
that the supplement's own sections are `S1.x` and `S2`–`S4`, so a bare "Section 6" must mean the
main document, and the reader has to infer that 58 times. In a document that is consulted rather
than read, that inference is made cold each time.

### (c) The contents page collides from S1.10 onward

Every entry with a two-digit subsection number runs its number into its title: **"S1.10The Saudi
context"**, "S1.11Where the linear alternative channel...", "S1.12Why the priors are...",
"S1.13Biological background", and every `S1.1x.y` entry likewise. A `\numberline` box too narrow for
two digits. **It is the first page a reference reader sees.**

### (d) S1.19 poses as open the question S1.20 answers, with no pointer between them

`S1.19` (p47) opens: *"Three quantities in the formulation are not measured... `mtol` (clinical
convention), `speak` (varies with the precipitating event), and the map from intracellular substrate
to plasma biomarker. θcrit is therefore an operational definition... **Whether the decision built on
it is invariant to that arbitrariness is the central question.**"*

`S1.20` (p50) opens: *"Well-posedness, and whether the decision is invariant to it. Asked above as
the central question. **Answered, and in the negative.** ...across their defensible range it moves by
a factor of 24 while the sign of the margin flips for the same patient with perfect data."*

**Read in order this is coherent.** Consulted at `S1.19`, it tells the reader the central question is
open when the next section says it is closed. `S1.19` carries no forward pointer.

### (e) S1.1 repeats its own heading as the first words of its body

*"`S1.1` The enzymatic term is Michaelis-Menten kinetics" / "**The enzymatic term is
Michaelis-Menten kinetics.** The function g is the Michaelis-Menten law..."*

### (f) One specific staleness suspicion, flagged for Part 3 rather than asserted

`S3` states: *"the decision operating characteristic uses **11** for the figure and **12, 13 and 21**
for the numbers quoted in Section 5 and Table 9."* That describes single-seed sourcing for the
decision numbers. **Whether the main document still quotes numbers that way is a cross-document
check**, so it is not decided here.

`S3` also states *"Everything under `results/` is regenerated from `src/`"*. **`results/` is now also
written by generators under `audit/`**, so that sentence is narrower than the truth it asserts.

## 1.3 What I went looking for and could not find

**A pointer from the supplement back to where each of its sections is cited from.** `S1` is material
moved out of the main text, and a reference reader arrives at, say, `S1.16.3` from somewhere. Nothing
in the section says where. The traffic is one-way: the main text points in, and no section says what
points at it. In a 64-page reference that is the single thing that would most help the reader the
framing describes.

**A statement of which of the two parallel treatments is authoritative.** With `S1.7` and `S1.16.5`
both pricing the carnitine route, nothing tells a reader which to cite. `S1.7` points forward to
`S1.16`, which suggests the long one, but it is inference.

**An index or a symbol table reachable from the front.** `Table S1: Notation` exists, on page 16,
which is inside `S1.14`. A reader consulting one row for one symbol has no way to know that from the
contents page, which does not list tables.

**And one absence I checked rather than assumed.** I expected `S2`'s countable claims about its own
table to be wrong, because a count above a table is the classic stale statement. They are right:
*"four rows lean on estimates fitted in mice and two on a rat in vitro preparation"*, and the table
carries exactly 4 `mouse` tags and 2 `rat, in vitro` tags. **I could not verify *"Nineteen of the
twenty-four rows are assumed"* from the extracted text**, because `pdftotext` does not preserve the
column, and I am recording that as unverified rather than guessing it.

## 1.4 The thing worth saying most plainly

**The supplement already contains, as its central open question, the finding this effort reached at
R163 and treated as new.** `S1.19` names `mtol`, `speak` and the substrate-to-biomarker map as the
three unmeasured quantities that make θcrit an operational definition, and `S1.20` quantifies the
consequence: **a factor of 24 in θcrit, and a sign flip in the margin for the same patient with
perfect data.**

R163 wrote θcrit's provenance as a precondition on the decidability result. **Applied to
`|mu|/sd(mu)` that is new. The underlying arbitrariness and its size were already in print, and
better quantified than I put them.** They were in the document nobody had read.

**That is the pass's own justification.** The defect the author predicted lives here does live here,
ten times over in the pointers. But the more expensive thing was a result the project already owned
and had lost track of.

---

# Part 2. Every pointer from the main document into the supplement

**Enumerated unbounded across all main-document sections, comments excluded: 29 pointers.**

| target as written | count | what it resolves to |
|---|---|---|
| `supplement, \S S1` | **19** | **pages 3 to 52. A 49-page section.** |
| `supplement \S S2` | **6** | pages 52 to 61 |
| named by content (`supplement's local-model section`, `supplement's fibre calculation`) | 2 | locatable |
| loose (`the supplement`, `supplement, which`) | 2 | the document |

## 2.1 The structural verdict, which is uniform and therefore not a table

**Not one pointer names a subsection.** There is no `\S S1.16.3` anywhere. All 19 `\S S1` pointers
resolve to a section that exists, contains the material, and is **49 pages long**.

**For the standard the author set, that is a failure of all 19.** A reader who arrives at
*"the argument is in the supplement, \S S1"* has been told the argument is somewhere in 49 pages,
in a document with a contents page whose entries collide from `S1.10` onward and which does not list
tables. Per-pointer verdicts would be 19 identical rows.

**And it explains the self-pointers of Part 1.** The whole project uses `S1` as an undifferentiated
destination, so a sentence inside `S1` saying *"set out in Appendix S1"* did not read as circular to
anyone writing it.

## 2.2 Promises tested individually. Two are broken, one was mine

Where a pointer promises something specific, I checked whether the supplement contains it, searching
for the **substance** rather than the main text's phrasing.

| pointer | promise | verdict |
|---|---|---|
| `07_design:98`, `:264` | *"repeats the table with one further column: the increment dropped from the sensitivity while still being used to generate the levels"* | **KEPT.** The material is in `S1` |
| `08_spec:381` | the load-prior figure sits at 1.05 and does not move | **KEPT** |
| `08_spec:472` | the excursion requirement and the null-genotype anchor | **KEPT** |
| `07b:93` | the margin theorem of the fibre calculation | **KEPT** |
| `07_design:620` | the validation protocol | **KEPT** (`S1.17`) |
| `02_model:93` | carnitine's second exit | **KEPT** |
| `12_disc:164` | the channel-exchange result | **KEPT** |
| **`07_design:451`** | *"The supplement, \S S1, **releases each assumption in turn**"* | **BROKEN** |
| **`06_excitation:146`** | *"The **one-line proof** is in the supplement, \S S1"* | **BROKEN** |
| **`07_design:504`** | *"The evidence for each row is in the supplement, \S S1"* | **WRONG TARGET, AND IT IS MINE** |

**`07_design:451` is broken and the chain is worth seeing.** The main text sends the reader to `S1`
for the release of the three assumptions. `S1.16.2` **states** them with motivations and opens
*"Three assumptions, stated here and warranted in **Appendix S1**"* — a self-pointer. And nothing in
the supplement releases them: searched unbounded for `releases each`, `release each`, `relaxed`,
`relaxing`, `released`, `one at a time` and `each assumption`, **all zero**; the single `in turn` is
about carnitine against enzymatic increments. **So the path is main text → \S S1 → "Appendix S1" →
nothing.**

**`06_excitation:146` is broken.** It promises the one-line proof of `prop:straddle`, *"the
excitation is maximised straddling `Km`"*. The supplement's only three `straddl` hits are a
dependency **diagram** cell reading *"Prop. 2: E is maximal straddling Km"*, a design row named
`straddling Km`, and an unrelated *"it straddles the target"*. **The proposition is referred to and
not proved there.**

**`07_design:504` is a defect I created in this effort, one pass ago.** R158 found that this caption
promised evidence in `\S S1` that was absent, and fixed it by adding the full admissibility table to
the supplement. **I put it in `\S S4` and left the caption pointing at `\S S1`.** So the promise
became satisfiable and the pointer stayed wrong. **Repointed to `\S S4` in this pass.** Repairing my
own breakage, on the same footing the author granted for the doubled clause: a section number that is
factually wrong is neither a claim nor a value.

---

# Part 3. Statements made in both documents

**The load-bearing shared figures, counted with fixed strings over both extracted texts.**

| figure | main | supplement | verdict |
|---|---|---|---|
| `14.02` / `12.86` | 1 / 1 | 1 / 1 | **AGREE** (the doubled clause fixed at R174 was the last disagreement here) |
| `82%` / `627%`, the band | 5 / 3 | 2 / 2 | **AGREE** |
| `factor of 24`, θcrit's movement | 1 | 1 | **AGREE** |
| `118` | 5 | 6 | **AGREE** |
| `1.79`, the no-dose-train footing | 2 | 2 | **AGREE** |
| **`0.1699`, the headline result** | **16** | **5** | — |
| **`0.170`, the same quantity** | **0** | **3** | **DISAGREE** |

## The one real disagreement, and it is the paper's central number

**The main document writes the headline design result as `0.1699` sixteen times and never as
`0.170`. The supplement writes it as `0.1699` five times and as `0.170` three times.**

The three are: *"On the footing that gives the dose train `sd(ln θ̂) = 0.170`"* (twice, in the two
duplicated carnitine remarks `S1.7` and `S1.16.5`), and *"two years reach `sd(ln θ̂) = 0.170`
(Section 4)"*.

**No number is wrong**: `0.170` is a correct rounding of `0.1699`. **But the supplement is
internally inconsistent and inconsistent with the main text about the precision of the one figure
the whole design turns on**, and it is inconsistent in the duplicated sections, so the duplication is
what carries it.

**Why it is not cosmetic.** R142 was a footing conflict decided at the **fourth** digit, `0.1699`
against `0.1695`, and this session settled the profile-marginal question on that digit. A document
that quotes the same quantity to three digits in one place and four in another is one transcription
away from re-opening exactly that. **Reported, not fixed: a rounding of a headline value is a value,
and the freeze protects values.**

---

# What this pass found, in order of what it cost to find

1. **Ten self-pointers** sending the reader to the document they are holding, three of them in the
   column that discharges every assumed constant. Found by reading.
2. **Two broken promises** from the main text, `07_design:451` and `06_excitation:146`, in a class
   with one prior instance. Found by enumeration, as the author predicted.
3. **One mis-pointed caption of my own**, created one pass ago while fixing the prior instance of
   the same class. Fixed here.
4. **The headline number in two precisions**, the supplement disagreeing with the main text and with
   itself. Found by enumeration.
5. **`S1.19` posing as open the question `S1.20` answers**, and with it the discovery that the
   supplement already owned, quantified better than R163 put it, the θcrit arbitrariness this effort
   treated as new. Found by reading, and no enumeration would have found it.

**The read and the enumeration found different classes, and neither substitutes for the other.** The
author's prediction that the enumeration would be higher yield is right by count, two broken promises
and a precision conflict against one structural observation. **It is wrong by consequence**: the
thing worth most was the rediscovery in 1.4, and only reading found it.

---

# Part 4. The contradiction with S1.20, resolved by reading its derivation

**The author stopped the defect list for this, and it was the right call: the resolution is a third
answer, and it supersedes both documents' summaries as well as mine.**

## 4.1 What the derivation actually establishes

The derivation is `S1.16.6`, page 33, equations (S35)–(S39). Read rather than recomputed.

**Step 1. The margin is better conditioned than the parameter it is built from, by one order.**
Parametrising the unidentified direction by `Km` and substituting the fitted pair `(a,b)`:

> *"The `Km²` terms cancel identically. The margin is **linear** along the unidentified direction
> where θ is **quadratic** in it, so the decision quantity is better conditioned than the parameter
> it is built from, **by one order**."*

**Step 2. Therefore the sign is constant along the fibre unless its endpoints disagree.**

> *"What survives is the sign, and **because µ is linear the sign is constant along the whole fibre
> unless the two endpoints disagree.** They disagree only for `speak` inside [(S38)]."*

**Step 3. Outside that window the decision is determined by the data even though θ is not.**

> *"Outside that window the decision is determined by the data even though θ is not, and at the
> illustrative values **θ is unknown there to a factor of 169**."*

**Step 4. And the window is narrow at the registered convention.**

> *"**The window is narrow where it matters.** At `mtol = 40` it is `(3.25, 3.77)` times the
> patient's baseline load, **14.8% wide** about its centre, and an illness of the amplitude the
> equilibrium reduction draws..., half again the baseline load, **sits well below it with µ > 0 at
> both ends**. Its width... closes as that convention approaches the band the patient is actually
> observed in: **1.2% at `mtol = 15`, 14.8% at 40, 67% at 160**."*

## 4.2 So which claim survives

**The paper's own derivation says the sign is ROBUST outside a narrow window, and that the illness
amplitude the model draws sits outside it with µ > 0 at both ends.** That is closer to the margin
document's conclusion than to `S1.20`'s.

| statement | verdict against the derivation |
|---|---|
| `S1.20` and `01_question`: *"the sign of the margin flips for the same patient with perfect data"* | **TRUE ONLY INSIDE THE WINDOW, and both state it without the window.** As written it reads unconditional. It is the strongest fair sentence about the 14.8% interior and a misleading one about everything else |
| my margin document: *"the sign of µ is robust in the disease range"* | **CONCLUSION RIGHT, ARGUMENT WRONG.** I supported it by `|mu|` being large at low capacity. The paper's reason is better and was already in print: `Km²` cancels, µ is linear where θ is quadratic, so the sign is constant along the fibre except in a computable window |
| the author's proposed retreat to *"relative decidability only"* | **NOT REQUIRED by this derivation.** The ordering survives, and so does the classification outside the window |

**The salvageable claim is therefore stronger than the author's fallback and weaker than S1.20's
headline:** the sign of µ is determined by the data outside a window in `speak` that is **14.8% wide
at `mtol = 40`**, and undetermined inside it, where no width of band settles the decision because
every point of the fibre fits equally well and they disagree on the sign. **The window's width is
itself the result, and it is computable per patient from `(a,b)`** — which is what makes this a
specification rather than a caveat.

## 4.3 The factor of 24 has no derivation in either document

**Searched exhaustively, fixed strings, both extracted texts.** `"factor of 24"` occurs **twice
total**: `01_question` in the main text and `S1.20` in the supplement. **Neither derives it**, and
`01_question` sources it to *"(supplement, §S1)"*, the 49-page pointer.

Every other `24` in either document is unrelated: 24 intake records, `Km = 124`, 24.5%, 24.1%,
`24%–647%`. The θcrit sensitivities that **are** derived are much smaller and are for the
**intrinsic** θcrit: **1.56** over a 32-fold `Kc` range, **1.37** over an eight-fold `Mi` range,
**1.6%** for `ρ0`. The figure the window derivation reports is **169**, and it is the uncertainty in
**θ**, not in θcrit.

**This is the third broken promise and the most serious**, because `S1.20` calls the claim it
supports *"the second of the paper's four claims and the reason it ends in a data specification
rather than a dosing rule."* **A load-bearing claim, asserted twice, sourced to a section whose
actual derivation says something different and narrower.**

## 4.4 What this does to R163

**R163's provenance paragraph is superseded by the document it was rediscovering**, and its
reasoning was weaker. It wrote θcrit's dependence on two conventions as a condition on the
decidability result. `S1.16.6` had already established the sharper thing: **the dependence is
confined to a computable window, and is a property of the fibre's endpoints rather than of the
conventions' ranges.** Recorded as superseded rather than corrected, since applying it to
`|mu|/sd(mu)` remains new.

---

# Part 5. The duplication, anatomised for the author's decision

**Three questions were asked. The third changes the decision.**

## 5.1 Which set is later: the SHORT one

**The long set is not supplement material at all. It is the original main text, displaced.**
`paper/sections/A1_supplementary.tex` lines **739–745** input seven files:

```
\input{sections/01_biology}    \input{sections/02_system}
\input{sections/03_dynamic}    \input{sections/04_local}
\input{sections/09_validation} \input{sections/10_results}
\input{sections/A1_open}
```

**None of those is inputted by `main.tex`.** They carry the original section numbering, and the
commit that did it names the operation: `f328ed4`, *"Rewrite Section 3 at half length, derivations to
the supplement"*.

So `A1_supplementary.tex` lines 1–738 are the **later, purpose-written supplement** (`S1.1`–`S1.12`),
and lines 739–745 pull in the **earlier originals** (`S1.13`–`S1.20`).

## 5.2 Which the pointers land in: neither, and that is the same finding

All 19 pointers say `\S S1`, which spans both. **The pointers cannot distinguish the sets.** What a
reader following one actually meets first is the short set, since it occupies pages 3 to 12.

## 5.3 Does either contain what the other lacks: YES, BOTH WAYS. It is a merge, not a deletion

**The overlap is four subsections, not two sets.** `S1.5`–`S1.8` against `S1.15.7`, `S1.15.8`,
`S1.16.5`, `S1.16.8`. Everything else in each set is unique:

| unique to the short set | unique to the long set |
|---|---|
| `S1.2` relation to the two published PK/PD models | `S1.13` biological background (3 subsections) |
| `S1.3` derivations for the excitation condition | `S1.14` what is known, measured, unknown, and Table S1 (notation) |
| `S1.4` two conditionings, kept separate | `S1.15` the dynamic model, 8 subsections, including the intrinsic θcrit |
| `S1.9` invariance under monotone rescaling | `S1.16` the local model, 8 subsections, **including the window derivation of Part 4** |
| `S1.10` the Saudi context | `S1.17` the validation protocol |
| `S1.11.1` what separates two channels in general | `S1.18` the figures |
| `S1.12` why the priors are what they are | `S1.19`, `S1.20` open and answered questions |

**And even within the four overlapping pairs neither side is a subset.** `S1.7` carries the therapy's
forty-two-fold excitation and the reach argument that `S1.16.5` lacks; `S1.16.5` carries Proposition
4's informational-versus-statistical distinction and the confounded-by-indication argument that
`S1.7` lacks. **Deleting either loses material.**

**So the decision is narrower and harder than "which set goes":** four pairs to merge, each losing
nothing, and the remaining ~19 subsections to leave alone. The one clean deletion available is the
**duplicate remark numbering** — one argument should not be both Remark 4 and Remark 17.
