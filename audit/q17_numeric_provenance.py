"""
Every numeric literal in paper/ that the artefact guard cannot see.

WHY. `tools/check_artefacts.py` regenerates the tables, the figures AND
`results/numbers/*.tex`, and compares each with what is committed. So a number
written `\\num{name}` is owned by a generator and guarded. A number typed
directly into the prose is owned by nobody, and nothing has ever compared it
with anything. R75 found what that costs: route 6's `2.659`, `thirteen`, `six`
and `eighteen` are literals, and one of them had been wrong since it was
written.

This is the sweep. It does not fix anything and does not touch paper/.

THE FOUR CATEGORIES, as the author set them:

  1  GENERATED-ELSEWHERE: the same quantity exists in results/numbers/*.tex or
     a committed table, and the prose retypes it instead of citing it.
  2  STRUCTURAL: an equation constant, a section or item count, a year, a
     version, a float specifier. Nothing generates it and nothing should.
  3  ASSERTING-COMPUTED: prose-only, and it states the outcome of a
     computation. THIS IS THE LIST THE AUTHOR ASKED FOR, because each member
     is a claim no generator owns and no check compares.
  4  HARMLESS: prose-only and asserting nothing computed. Page furniture,
     rhetorical quantities, ranges quoted from a citation.

Category 3 is reported with, for each member, whether a generator exists that
could own it and whether the value still reproduces. The reproduction checks
themselves live in `q17_reproduce.py`, because they are slow and this file is a
scanner.

Run:  python audit/q17_numeric_provenance.py
"""
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper"
#: every macro whose ARGUMENT is a generated value, so the literal inside is
#: owned by a generator and not by the prose
OWNED = ("num",)
#: LaTeX control sequences whose numeric arguments are typesetting, not claims
TYPESET = re.compile(
    r"\\(?:"
    r"(?:v|h)space\*?|(?:set|add)length|(?:top|bottom|text|line|column|paper)"
    r"(?:width|height|sep|skip)|includegraphics|scalebox|resizebox|"
    r"(?:small|med|big)skip|arraystretch|tabcolsep|extrarowheight|"
    r"baselinestretch|rule|raisebox|hspace|kern|mskip|multicolumn|cline|"
    r"cmidrule|midrule|toprule|bottomrule|label|ref|eqref|cite|citep|citet|"
    r"input|usepackage|documentclass|geometry|definecolor|color|textcolor|"
    r"pagenumbering|setcounter|refstepcounter|thepage|footnotemark"
    r")\b")
NUMBER = re.compile(r"(?<![A-Za-z0-9_@.\\])(\d+(?:[.,]\d+)*)(?![A-Za-z0-9_@])")


def strip_comments(line):
    """LaTeX comments, honouring the escaped percent."""
    out, esc = [], False
    for ch in line:
        if esc:
            out.append(ch)
            esc = False
            continue
        if ch == "\\":
            out.append(ch)
            esc = True
            continue
        if ch == "%":
            break
        out.append(ch)
    return "".join(out)


def mask_owned(line):
    """Blank out the ARGUMENTS of guarded macros, so their literals do not
    count as prose. Brace-matched, because arguments nest."""
    for name in OWNED:
        i = 0
        pat = "\\" + name + "{"
        while True:
            j = line.find(pat, i)
            if j < 0:
                break
            k, depth = j + len(pat), 1
            while k < len(line) and depth:
                depth += (line[k] == "{") - (line[k] == "}")
                k += 1
            line = line[:j] + " " * (k - j) + line[k:]
            i = j + 1
    return line


def mask_typeset(line):
    """Blank a typesetting control sequence and the rest of its argument run,
    so lengths and float specifiers do not enter the census."""
    out = line
    for m in list(TYPESET.finditer(line)):
        k = m.end()
        #  swallow following bracket/brace groups and bare length arguments
        while k < len(out) and out[k] in " ":
            k += 1
        while k < len(out) and out[k] in "[{":
            close = "]" if out[k] == "[" else "}"
            depth, k2 = 1, k + 1
            while k2 < len(out) and depth:
                if out[k2] == out[k]:
                    depth += 1
                elif out[k2] == close:
                    depth -= 1
                k2 += 1
            k = k2
        m2 = re.match(r"\s*-?[\d.]+\s*(pt|em|ex|cm|mm|in|bp|sp|\\[a-z]+)?", out[k:])
        if m2:
            k += m2.end()
        out = out[:m.start()] + " " * (k - m.start()) + out[k:]
    return out


def sweep():
    rows = []
    files = sorted(PAPER.glob("*.tex")) + sorted((PAPER / "sections").glob("*.tex"))
    for f in files:
        for n, raw in enumerate(f.read_text().splitlines(), 1):
            line = mask_typeset(mask_owned(strip_comments(raw)))
            for m in NUMBER.finditer(line):
                rows.append({"file": str(f.relative_to(ROOT)), "line": n,
                             "value": m.group(1), "ctx": raw.strip()[:160]})
    return rows


def main():
    rows = sweep()
    print("=" * 100)
    print("  NUMERIC LITERALS IN paper/ THAT THE ARTEFACT GUARD CANNOT SEE")
    print("=" * 100)
    print(f"\n  {len(rows)} unguarded literals across "
          f"{len(set(r['file'] for r in rows))} files.")
    print("\n  By file:")
    c = Counter(r["file"] for r in rows)
    for f, k in c.most_common():
        print(f"    {k:5d}  {f}")
    #  and how many ARE guarded, for the ratio
    tot = sum(len(re.findall(r"\\num\{", p.read_text()))
              for p in list(PAPER.glob("*.tex"))
              + list((PAPER / "sections").glob("*.tex")))
    print(f"\n  Guarded (\\num) references: {tot}. Unguarded literals: "
          f"{len(rows)}. Ratio {len(rows) / max(tot, 1):.2f} unguarded per "
          f"guarded.")
    out = ROOT / "audit" / "q17_literals.txt"
    with out.open("w") as fh:
        for r in sorted(rows, key=lambda r: (r["file"], r["line"])):
            fh.write(f"{r['file']}:{r['line']}\t{r['value']}\t{r['ctx']}\n")
    print(f"\n  Full census written to {out.relative_to(ROOT)} for "
          f"classification.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
