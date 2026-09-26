"""
For each category-3 claim: does a generator exist that could own it?

THE CRITERION, and it is deliberately checkable rather than interpretive. A
generator "exists" for a literal if the SAME PARAGRAPH already cites at least
one `\\num{}` macro. That macro is written by a named audit script into a named
`results/numbers/*.tex` file, so that script is already loaded, already runs in
`make audit`, and already emits numbers into this paragraph. Adding one more is
an edit to the generator, not a new generator.

Where no `\\num{}` appears in the paragraph the answer is `none in paragraph`,
which is the harder case: the claim sits in prose no generator reaches, and
somebody would have to decide which script should own it.

Run:  python audit/q17_owners.py
"""
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NUMDIR = ROOT / "results" / "numbers"
CAT3 = ROOT / "audit" / "q17_category3.txt"


def macro_owner():
    """macro name -> the numbers file that defines it."""
    own = {}
    for p in sorted(NUMDIR.glob("*.tex")):
        for m in re.finditer(r"csname paper@num@([A-Za-z0-9_]+)\\endcsname",
                             p.read_text()):
            own[m.group(1)] = p.stem
    return own


def paragraphs(text):
    """(start_line, end_line) for each blank-line-delimited paragraph."""
    out, start = [], 1
    lines = text.splitlines()
    for n, line in enumerate(lines, 1):
        if not line.strip():
            if n > start:
                out.append((start, n - 1))
            start = n + 1
    if start <= len(lines):
        out.append((start, len(lines)))
    return out


def main():
    own = macro_owner()
    claims = defaultdict(list)
    for line in CAT3.read_text().splitlines():
        loc, val, why, ctx = line.split("\t", 3)
        f, n = loc.rsplit(":", 1)
        claims[(f, int(n))].append(val)
    print("=" * 100)
    print("  DOES A GENERATOR EXIST FOR EACH CATEGORY-3 CLAIM")
    print("=" * 100)
    cache = {}
    rows = []
    for (f, n), vals in sorted(claims.items()):
        text = cache.setdefault(f, (ROOT / f).read_text())
        para = next((p for p in paragraphs(text) if p[0] <= n <= p[1]),
                    (n, n))
        block = "\n".join(text.splitlines()[para[0] - 1:para[1]])
        macros = re.findall(r"\\num\{([A-Za-z0-9_]+)\}", block)
        owners = sorted({own.get(m, "?") for m in macros})
        rows.append({"file": f, "line": n, "vals": vals,
                     "owners": owners,
                     "ctx": " ".join(block.split())[:200]})
    have = [r for r in rows if r["owners"]]
    print(f"\n  {len(rows)} claims, {sum(len(r['vals']) for r in rows)} "
          f"literals.")
    print(f"  {len(have)} sit in a paragraph that already cites a generated "
          f"macro, so a generator exists and could own them.")
    print(f"  {len(rows) - len(have)} do not, and would need an owner chosen.")
    print("\n  Which generator is already present, by claim count:")
    c = Counter(o for r in have for o in r["owners"])
    for k, v in c.most_common():
        print(f"    {v:5d}  results/numbers/{k}.tex")
    out = ROOT / "audit" / "q17_claims.tsv"
    with out.open("w") as fh:
        fh.write("file\tline\tvalues\tgenerator_in_paragraph\tcontext\n")
        for r in rows:
            fh.write(f"{r['file']}\t{r['line']}\t{','.join(r['vals'])}\t"
                     f"{'+'.join(r['owners']) or 'NONE'}\t{r['ctx']}\n")
    print(f"\n  Written to {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
