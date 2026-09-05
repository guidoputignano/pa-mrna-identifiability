# References and the claim register

Two artefacts, one rule: **no citation is typed by hand**, and **no factual
claim goes into the paper without a row in the register**.

```
paper/refs.sources        citekey -> DOI or PMID. The only file you edit.
paper/refs.manual.bib     entries with no DOI. Hand-typed, so hand-checked.
paper/refs.bib            GENERATED. Do not edit.
paper/refs.todo           GENERATED. Identifiers that failed to resolve.
paper/claims.csv          one row per factual assertion in the paper.
results/tables/tab_claims.tex   GENERATED from claims.csv.
tools/fetch_refs.py       pulls BibTeX from CrossRef and PubMed.
tools/claims_to_tex.py    renders the register as a LaTeX table.
```

## Why a fetcher instead of a `.bib`

Bibliographies written from memory acquire wrong volumes, wrong years and
authors who were never on the paper. Here every field in `refs.bib` comes
from the registrar. `refs.sources` holds identifiers only, which are short
enough to check by eye against the record. Anything that fails to resolve
lands in `refs.todo` and is never written into the bibliography, so a
broken identifier shows up as a missing citation at build time rather than
as a plausible-looking but wrong entry.

## Build

```bash
pip install requests
python tools/fetch_refs.py          # writes paper/refs.bib
python tools/claims_to_tex.py       # writes results/tables/tab_claims.tex
```

Set `MAILTO` at the top of `fetch_refs.py` to a real address. CrossRef
serves anonymous requests from a slower pool.

Both scripts take `--check`, which exits non-zero when something is
unresolved. Wire that into CI so an unsourced claim fails the build.

`claims_to_tex.py --check` also closes the loop in both directions between
the register and the paper:

- **forwards**, every citekey a row names has an entry in `refs.bib`, so a
  row cannot point at a source the bibliography does not carry;
- **backwards**, every `\cite{}` key in `paper/**/*.tex` has at least one
  row, so a citation cannot enter the paper without a claim recorded
  against it.

An unauthorised citation fails the build exactly as an unresolved claim
does. A mistyped citekey trips both directions at once, which is the
signature to look for.

## LaTeX wiring

In `paper/preamble.tex`, after `hyperref`:

```latex
\usepackage[backend=biber, style=numeric-comp, sorting=none,
            maxbibnames=6, minbibnames=3, giveninits=true]{biblatex}
\addbibresource{refs.bib}
```

At the end of `paper/main.tex`, before `\end{document}`:

```latex
\printbibliography[title={References}]

\appendix
\section{Claim register}
Every factual assertion in this note, its source, and whether the source
says what the assertion says. Rows are ordered worst first.
\begin{adjustwidth}{-1.4cm}{-1.4cm}
\footnotesize
\input{../results/tables/tab_claims}
\end{adjustwidth}
```

Then in `Makefile`, replace the `paper` target:

```makefile
refs:
	$(PY) tools/fetch_refs.py
	$(PY) tools/claims_to_tex.py

paper: refs
	cd paper && $(LATEX) -interaction=nonstopmode main.tex >/dev/null
	cd paper && biber main >/dev/null
	cd paper && $(LATEX) -interaction=nonstopmode main.tex >/dev/null
	cd paper && $(LATEX) -interaction=nonstopmode main.tex >/dev/null
	@echo ">> built $(PAPER).pdf"

audit:
	$(PY) tools/fetch_refs.py --check
	$(PY) tools/claims_to_tex.py --check
```

`biber` needs the `.bcf`, so the first `pdflatex` pass has to run before
it. If you would rather stay on `bibtex`, drop `backend=biber` and use
`natbib` instead; the generated `.bib` works with either.

The register is emitted as a `longtable`, not a `tabular`. Thirty rows of
wrapped prose do not fit on one page and a `tabular` cannot break, so it
runs off the bottom without complaining. `longtable` builds its alignment
at the outer level and positions itself with `\LTleft`/`\LTright`, which
means it does not inherit the indentation of the enclosing `adjustwidth`:
the column widths are fractions of `\textwidth` and must sum to about
0.88, since four inter-column gaps add roughly another 0.10.

## Adding a reference

1. Find the DOI or PMID on the record itself, not from a search snippet.
2. Add a line to `paper/refs.sources` with a comment saying what the source
   actually reports, in its own units.
3. Run `python tools/fetch_refs.py`.
4. Add the row to `paper/claims.csv` before writing the sentence that
   cites it.

If there is no DOI, put the entry in `refs.manual.bib` with a `note`
recording which fields still need checking. Check first that there really
is none: `dose_selection_poster` was assumed to have no usable identifier
and turned out to resolve through CrossRef, which supplied an eleven-author
list, the publisher and the year that a hand-written entry would have had
to leave marked "unverified". Anything in `refs.manual.bib` that later
acquires an identifier moves to `refs.sources`.

## What the fetcher normalises

Three transformations, all of them re-spellings of what the registrar
returned rather than additions to it.

- **Author names from PubMed.** esummary returns `Al-Hamed MH`. BibTeX
  reads an unpunctuated name as *First Last*, making the family name `MH`
  and rendering it `A.-H. MH`. Names are emitted as `Al-Hamed, M. H.`
  instead. A name whose final token is not an initials group is braced, so
  BibTeX treats it as one literal unit rather than guessing.
- **Month macros.** CrossRef emits `month=July`. BibTeX's month macros are
  three letters, so `July` is undefined, biber warns and drops the field.
  It is mapped onto `jul`.
- **Unicode.** Greek letters go to math mode, dashes and quotes to ASCII.
  Anything non-ASCII that survives and is not a Latin-1 letter is reported
  at the end of the run, so a new source cannot quietly introduce a
  codepoint the build may not accept.

## Verdicts in `claims.csv`

| verdict | meaning |
|---|---|
| `matches` | the source says what the paper says, in the same units and under the same definition |
| `matches-with-drift` | the number is right but the definition, population or qualifier differs. **This is the dangerous one.** |
| `contradicted` | the source says something else |
| `unsourced` | no source found yet |

`matches-with-drift` exists because this literature reuses words for
different things. A decompensation event is hospitalisation-coded in the
claims studies and investigator-judged in the natural history study.
Early versus late onset is cut at three months in one paper and at
28 days in others. A rate moved across those boundaries is not the same
rate. Any number entering the paper carries the definition that produced
it, or it does not enter.

## Current state

30 claims: 9 match, 7 drift, 7 contradicted, 7 unsourced.

The seven contradicted rows are not typos. Four of them change the
argument:

- **C08** the "70 percent experience no event" figure is the wrong way
  round. Open item 2 uses it as motivation.
- **C19** 2-MC, 3-HP and C3 are treated as one quantity in Table 1, but
  only two of the three carry a risk signal.
- **C23** the `m << Km` inequality compares plasma C3 with an
  intracellular propionyl-CoA constant. The degeneracy result rests on it.
- **C24** the measured biomarker does not rise through a decompensation,
  which is what both Model B and Model C predict.

Resolve these before extending the models. C23 and C24 are the two that
decide whether the current state variable is the right one.
