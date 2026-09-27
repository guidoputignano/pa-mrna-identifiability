"""
The PHASE4 three-way reconciliation the author asked for: applied, stale,
orphaned, per passage, as a reading rather than a decision.

WHY A SCRIPT RATHER THAN A READING. The author's framing was that this is
*"yours to do, and it is a reading, not a decision"*, and that he decides only
the orphans. A reading of fifteen passages against a 45-page document is exactly
where an assertion drifts from what the files say, so the classification is
computed and the evidence for each verdict is printed with it.

HOW A PASSAGE IS CLASSIFIED. Each passage in `PHASE4_PROPOSED_TEXT.md` has a
`### Current` block, the text it would replace, and a `### Proposed` block, the
replacement. Both are matched against the live section files:

  APPLIED   the PROPOSED text is in the paper. The edit already happened, by
            whatever route, and the passage is spent.
  STALE     the CURRENT text is in the paper and the proposed text is not. The
            passage still has its target and could still be applied. This is the
            only category that is actionable as written.
  ORPHANED  NEITHER is in the paper. The target has been rewritten or removed
            since the passage was drafted, so there is nothing to replace and
            nothing to compare against. **These are the author's to decide**,
            because withdrawing a passage is a judgement about whether its point
            still needs making.
  BOTH      both match, which would mean the passage is partly applied. Reported
            separately rather than folded into one of the three, because it is
            the one state that needs a human to look at the diff.

MATCHING, AND WHY IT IS A SHINGLE RATHER THAN AN EXACT COMPARE. The `Current`
blocks quote the paper's prose with line breaks and LaTeX normalised away by
whoever drafted them, so an exact string compare fails on whitespace and reports
every passage orphaned. Text is therefore reduced to lowercase alphanumeric
words and matched on the longest run of consecutive words found in the paper.
The run length is printed, so a weak match is visible as a weak match instead of
being rounded to a verdict. `MIN_RUN` is the threshold and it is stated rather
than tuned: twelve consecutive words is far beyond coincidence in this corpus
and far below the length of any quoted block.

Nothing here records an expected result in advance, and nothing is decided: the
orphan list is output, not action.

Run:  python audit/u27_phase4_reconcile.py
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "audit" / "PHASE4_PROPOSED_TEXT.md"
SECTIONS = sorted((ROOT / "paper" / "sections").glob("*.tex"))

#: consecutive normalised words that count as a match. Stated, not tuned.
MIN_RUN = 12


def words(text):
    """Lowercase alphanumeric words, LaTeX and punctuation discarded."""
    text = re.sub(r"\\[A-Za-z]+\*?", " ", text)
    text = re.sub(r"[^A-Za-z0-9 ]+", " ", text)
    return [w for w in text.lower().split() if w]


def longest_run(needle, haystack):
    """Longest run of consecutive `needle` words appearing in `haystack`."""
    if not needle:
        return 0
    hay = " " + " ".join(haystack) + " "
    best = 0
    for i in range(len(needle)):
        #  extend greedily from i while the phrase is still present
        j = i + best
        while j < len(needle):
            phrase = " " + " ".join(needle[i:j + 1]) + " "
            if phrase in hay:
                best = j - i + 1
                j += 1
            else:
                break
    return best


def target_is_register(title, body):
    """
    Does this passage target `src/common/params.py` rather than the paper?

    THIS DISTINCTION WAS MISSING AND IT COST TWO FALSE ORPHANS. The first
    version searched every passage against `paper/sections` and reported six
    orphans. Two of them target the register: passage 8 rewrites the `D` row's
    justification, and passage 10 five annotation rows whose own text says
    "`tab_params.tex` is generated from `params.py`". Their text is absent from
    the paper because it was never meant to be there. **A wrong haystack gives a
    confident wrong verdict**, the same defect as reading a threshold off the
    wrong design (REVIEW_LOG.md R153), and here it would have handed the author
    two decisions that were not his to make.

    Title and body are both consulted, because passage 10's title names no file
    at all.
    """
    blob = (title + " " + body).lower()
    return "params.py" in blob or "tab_params" in blob


def passages():
    """Each numbered passage with its title, Current block and Proposed block."""
    raw = DOC.read_text(errors="replace")
    parts = re.split(r"\n## (\d+)\. ", raw)
    out = []
    for i in range(1, len(parts), 2):
        num, body = parts[i], parts[i + 1]
        title = body.splitlines()[0].strip()

        def block(name):
            m = re.search(rf"### {name}[^\n]*\n(.*?)(?=\n### |\n## |\Z)",
                          body, re.S)
            if not m:
                return ""
            #  only the quoted lines are the text itself
            return "\n".join(l.lstrip("> ").rstrip()
                             for l in m.group(1).splitlines()
                             if l.lstrip().startswith(">"))
        out.append((int(num), title, block("Current"), block("Proposed"),
                    target_is_register(title, body)))
    return out


if __name__ == "__main__":
    hay = words("\n".join(p.read_text(errors="replace") for p in SECTIONS))
    reg = words((ROOT / "src" / "common" / "params.py").read_text(errors="replace"))
    print("=" * 100)
    print("PHASE4: APPLIED, STALE, ORPHANED, per passage")
    print("=" * 100)
    print(f"\n  {len(SECTIONS)} section files, {len(hay)} words searched, "
          f"{len(reg)} words in the register, "
          f"match threshold {MIN_RUN} consecutive words\n")

    rows, tally = [], {}
    for num, title, cur, prop, is_reg in passages():
        wc, wp = words(cur), words(prop)
        if not wc and not wp:
            #  A PASSAGE WITH NO QUOTED BLOCKS IS NOT A TEXT REPLACEMENT.
            #  Passage 10 is a list of five annotation edits described in prose,
            #  with no Current/Proposed pair to match. Calling it orphaned would
            #  report an absence the script never looked for, which is the
            #  bounded-search error in miniature.
            v = "NO-BLOCKS"
        else:
            h = reg if is_reg else hay
            rc = longest_run(wc, h) if wc else 0
            rp = longest_run(wp, h) if wp else 0
            has_c, has_p = rc >= MIN_RUN, rp >= MIN_RUN
            if is_reg:
                v = ("REGISTER-APPLIED" if has_p else
                     "REGISTER-PENDING" if has_c else "REGISTER-UNMATCHED")
            elif has_c and has_p:
                v = "BOTH"
            elif has_p:
                v = "APPLIED"
            elif has_c:
                v = "STALE"
            else:
                v = "ORPHANED"
        rc = rc if (wc or wp) else 0
        rp = rp if (wc or wp) else 0
        tally[v] = tally.get(v, 0) + 1
        rows.append((v, num, title, len(wc), rc, len(wp), rp))

    ORDER = ("STALE", "ORPHANED", "BOTH", "APPLIED", "REGISTER-PENDING",
             "REGISTER-APPLIED", "REGISTER-UNMATCHED", "NO-BLOCKS")
    print(f"  {'verdict':>18} {'#':>3}  {'cur':>5} {'match':>6}   "
          f"{'prop':>5} {'match':>6}  target")
    print("  " + "-" * 96)
    for v in ORDER:
        for _, num, title, nc, rc, np_, rp in [x for x in rows if x[0] == v]:
            print(f"  {v:>18} {num:>3}  {nc:>5} {rc:>6}   {np_:>5} {rp:>6}  "
                  f"{title[:40]}")

    print()
    print("=" * 100)
    print("  " + "   ".join(f"{k} {tally.get(k, 0)}" for k in ORDER
                            if tally.get(k, 0)))
    print("=" * 100)
    orphans = [r[1] for r in rows if r[0] == "ORPHANED"]
    stale = [r[1] for r in rows if r[0] == "STALE"]
    print(f"""
  STALE ({len(stale)}) is the only category actionable as written: the target text is
  still in the paper and the passage still replaces it. {stale if stale else ''}

  ORPHANED ({len(orphans)}) is the author's: {orphans}. Neither the text the passage
  would replace nor the text it proposes is in the paper, so it cannot be applied
  and cannot be checked against anything. Withdrawing it is a judgement about
  whether its point still needs making, which is not a reading.

  REGISTER rows target src/common/params.py and are searched there instead. They
  are not orphans and they are not the author's to decide on this table.

  A LOW `match` IS NOT NOTHING. It means the block and the haystack share a few
  words, which is what two texts on one subject do. The threshold separates that
  from a real match and is printed so the margin stays visible.""")
