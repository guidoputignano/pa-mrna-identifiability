#!/usr/bin/env python3
"""
Render paper/claims.csv as results/tables/tab_claims.tex.

The point of putting the audit in the paper is that a claim whose
verdict is 'contradicted' or 'unsourced' becomes visible rather than
quietly surviving a revision. Rows are sorted worst-first.

Usage:
    python tools/claims_to_tex.py
    python tools/claims_to_tex.py --check    # exit 1 if anything is unresolved
"""

import argparse
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "paper" / "claims.csv"
OUT = ROOT / "results" / "tables" / "tab_claims.tex"
BIB = ROOT / "paper" / "refs.bib"
TEXDIR = ROOT / "paper"

ORDER = {"contradicted": 0, "unsourced": 1, "matches-with-drift": 2, "matches": 3}
MARK = {
    "matches": r"\textcolor{green!45!black}{matches}",
    "matches-with-drift": r"\textcolor{orange!80!black}{drift}",
    "contradicted": r"\textbf{\textcolor{red!70!black}{contradicted}}",
    "unsourced": r"\textcolor{red!70!black}{unsourced}",
}


def esc(s):
    """
    Escape a register cell for LaTeX.

    claims.csv is plain prose, so every LaTeX special has to be neutralised,
    not just the four that happened to appear first. A source value quoting
    "R^2" put the table into math mode and broke the build, and because the
    appendix only renders under the draft option, `make paper` stayed green
    while `make paper-draft` failed. Backslash goes first, or the escapes
    introduced below would themselves be escaped.
    """
    for a, b in [("\\", r"\textbackslash{}"),
                 ("{", r"\{"), ("}", r"\}"),
                 ("&", r"\&"), ("%", r"\%"), ("_", r"\_"), ("#", r"\#"),
                 ("$", r"\$"),
                 ("^", r"\textasciicircum{}"), ("~", r"\textasciitilde{}")]:
        s = s.replace(a, b)
    return s


def citekeys(row):
    """A row may name more than one source; C17 and C29 do."""
    return [k.strip() for k in row["citekey"].split(",") if k.strip()]


# A claim-id comment must be the last thing on its line. LaTeX discards
# everything after `%` including the newline, so `... % C02 Some text` silently
# deletes "Some text" and joins what follows straight onto the citation. The
# source looks correct and the loss is visible only in the PDF, which is how one
# such line survived until a full read-through. `make audit` now fails on it.
CLAIM_COMMENT = re.compile(r"(?<!\\)%\s*C\d+\b(?P<rest>.*)$")
TRAILING_IDS = re.compile(r"^(?:[,;]?\s*C\d+\b)*\s*$")


def check_claim_comments():
    """Every `% Cxx` marker sits at end of line, with nothing after it."""
    problems = []
    for tex in sorted(TEXDIR.rglob("*.tex")):
        for n, line in enumerate(tex.read_text().splitlines(), 1):
            m = CLAIM_COMMENT.search(line)
            if m and not TRAILING_IDS.match(m.group("rest")):
                problems.append(
                    f"{tex.relative_to(ROOT)}:{n} has text after a claim id, "
                    f"which LaTeX discards: {m.group('rest').strip()!r}")
    return problems


def cross_check(rows):
    r"""
    Both directions between the register and the paper.

    Forwards: every citekey the register names has an entry in refs.bib, so a
    row cannot point at a source the bibliography does not carry.

    Backwards: every \cite key in paper/**/*.tex has at least one row, so a
    citation cannot enter the paper without a claim recorded against it. An
    unauthorised citation fails the build exactly as an unresolved claim does.

    Returns a list of problems, empty when both hold.
    """
    problems = []

    declared = {k for r in rows for k in citekeys(r)}
    if BIB.exists():
        in_bib = set(re.findall(r"@\w+\s*\{\s*([^,\s]+)\s*,", BIB.read_text()))
        for key in sorted(declared - in_bib):
            owners = ", ".join(r["id"] for r in rows if key in citekeys(r))
            problems.append(f"citekey {key!r} ({owners}) is not in refs.bib")
    else:
        problems.append(f"{BIB.relative_to(ROOT)} missing: run tools/fetch_refs.py")

    for tex in sorted(TEXDIR.rglob("*.tex")):
        for group in re.findall(r"\\cite\{([^}]*)\}", tex.read_text()):
            for key in (k.strip() for k in group.split(",")):
                if key and key not in declared:
                    problems.append(
                        f"{tex.relative_to(ROOT)} cites {key!r}, "
                        f"which no row in claims.csv authorises")
    return problems


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    rows = list(csv.DictReader(SRC.open()))
    rows.sort(key=lambda r: (ORDER.get(r["verdict"], 9), r["id"]))

    # longtable, not tabular: thirty rows of wrapped prose do not fit on one
    # page, and a tabular cannot break, so it would silently run off the
    # bottom.
    #
    # Widths are fractions of \textwidth. longtable builds its alignment at
    # the outer level and positions itself with \LTleft/\LTright, so it does
    # not inherit the indentation of an enclosing adjustwidth: sizing against
    # the local line width there gives a table wider than the block it is
    # drawn in, whose last column runs off the page.
    head = (r"\textbf{} & \textbf{Claim} & \textbf{Source} & "
            r"\textbf{Verdict} & \textbf{What the source says} \\")
    lines = [
        r"% Generated by tools/claims_to_tex.py. DO NOT EDIT.",
        r"% Edit paper/claims.csv and re-run.",
        # Fractions sum to 0.88, not 1: four inter-column gaps of
        # 2\tabcolsep add about 0.10\textwidth on top. The widths given in
        # the task sheet summed to 0.945 and overflowed, which showed up as
        # the verdict column running into the one beside it.
        r"\begin{longtable}{@{}p{0.03\textwidth} p{0.36\textwidth} "
        r"p{0.05\textwidth} p{0.115\textwidth} p{0.325\textwidth}@{}}",
        r"\toprule", head, r"\midrule", r"\endfirsthead",
        r"\toprule", head, r"\midrule", r"\endhead",
        r"\midrule",
        r"\multicolumn{5}{@{}r@{}}{\footnotesize\itshape continued overleaf} \\",
        r"\endfoot",
        r"\bottomrule", r"\endlastfoot",
    ]
    for r in rows:
        cite = rf"\cite{{{r['citekey']}}}" if r["citekey"] else r"\emph{none}"
        lines.append(
            f"{r['id']} & {esc(r['claim'])} & {cite} & "
            f"{MARK.get(r['verdict'], r['verdict'])} & "
            f"{esc(r['source_value'] or r['note'])} \\\\"
        )
    lines += [r"\end{longtable}"]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n")

    counts = {v: sum(1 for r in rows if r["verdict"] == v) for v in ORDER}
    print(f"wrote {OUT.relative_to(ROOT)}  ({len(rows)} claims)")
    for v, n in counts.items():
        print(f"  {v:20s} {n}")

    problems = cross_check(rows) + check_claim_comments()
    for msg in problems:
        print(f"  ! {msg}")

    unresolved = counts["contradicted"] + counts["unsourced"]
    if args.check and (unresolved or problems):
        why = []
        if unresolved:
            why.append(f"{unresolved} claims unresolved")
        if problems:
            why.append(f"{len(problems)} register/paper mismatch"
                       f"{'es' if len(problems) > 1 else ''}")
        sys.exit("\n" + "; ".join(why))


if __name__ == "__main__":
    main()
