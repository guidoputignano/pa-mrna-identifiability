"""
Step 1 of the landing pass: does each proposed passage still apply to the paper
as it is today?

WHY THIS RUNS AS A CHECK AND NOT AS A READING. `audit/PHASE4_PROPOSED_TEXT.md`
was drafted several phases ago. Since then Phase 10 corrected three numbers in
`paper/` and converted 22 bare literals to `\\num{}` macros, and later rounds
moved more. **A proposal whose "Current" text has since changed will apply
cleanly and wrongly**, and nothing in the output will look malformed. That is the
same failure mode as the shallow clone (R93) and the bounded search (R125), and
`audit/LANDING_PLAN.md` names it as the concentrated risk of landing. So it gets
an explicit PASS or FAIL per passage rather than a careful read.

WHAT IT CHECKS. Every `## N.` section of a proposed-text file names its target in
backticks in the heading and quotes the text to be replaced under `### Current`.
The check is whether that quoted text is still present, VERBATIM up to
whitespace, in the named file. Whitespace is normalised because the proposals
re-wrap lines and LaTeX does not care; nothing else is normalised, so a changed
digit or a changed macro is a FAIL.

VERDICTS, and the second one is why this script exists in its present form. A
first version reported only present-or-absent and returned ten absences, which
reads as ten stale proposals. Checking whether the PROPOSED text is present
instead showed the opposite: the passages had ALREADY BEEN APPLIED, so their
Current text is gone because it was replaced. An absence alone cannot tell those
two apart, and they call for opposite actions.

  PENDING   Current present, Proposed absent. Applies as drafted.
  APPLIED   Proposed present, Current absent. Already in the paper; the
            proposed-text file is stale as a to-do list, not the paper.
  STALE     neither present. The paper moved underneath it; must be rebased.
  PARTIAL   both present. Ambiguous, and a human must look.
  NO-ANCHOR the section quotes no Current text, usually because it is an
            insertion rather than a replacement. This check cannot speak to it.
  NO-TARGET the heading names no file this script can resolve.

Nothing here records an expected result in advance.

Run:  python audit/u20_staleness.py
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

FILES = ["audit/PHASE4_PROPOSED_TEXT.md",
         "audit/PHASE10_PROPOSED_TEXT.md",
         "audit/PHASE11_PROPOSED_TEXT.md"]

#: a heading's backticked target, e.g. `paper/sections/06_excitation.tex:71-75`
TARGET = re.compile(r"`([A-Za-z0-9_./]+\.(?:tex|py))(?::[\d\-]+)?`")


def norm(s):
    """Collapse whitespace. Nothing else is normalised."""
    return re.sub(r"\s+", " ", s).strip()


def sections(text):
    """(number, heading, body) for each '## N.' section."""
    out, cur = [], None
    for line in text.splitlines():
        m = re.match(r"^## (\d+)\.\s*(.*)$", line)
        if m:
            if cur:
                out.append(cur)
            cur = [m.group(1), m.group(2), []]
        elif cur is not None:
            cur[2].append(line)
    if cur:
        out.append(cur)
    return [(n, h, "\n".join(b)) for n, h, b in out]


def block(body, name):
    """The text under '### <name>', with blockquote markers stripped."""
    m = re.search(r"^### " + name + r"\s*$(.*?)(?:^### |\Z)", body, re.S | re.M)
    if not m:
        return None
    lines = [re.sub(r"^> ?", "", ln) for ln in m.group(1).splitlines()
             if ln.strip().startswith(">")]
    out = "\n".join(lines).strip()
    return out or None


def current_block(body):
    return block(body, "Current")


def probes(txt, width=55, want=8):
    """
    Several short fragments that survive later WIRING, for asking whether a
    passage reached the paper.

    A single long verbatim probe does not work, and finding that out is the
    point of this function. PHASE4 section 1's proposed text IS in the paper,
    and a 140-character probe through its middle misses, because Phase 10
    afterwards replaced the literal $1.028$ with $\num{bias_theta_factor}$ and
    added a trailing clause. Applied text keeps being edited, so the probe has
    to be short and has to avoid the parts most likely to have been rewritten:
    anything with a digit, a macro or a math delimiter.
    """
    t = norm(txt)
    #  candidate windows on sentence-ish boundaries, then filtered
    cands = [c.strip() for c in re.split(r"(?<=[.;:]) ", t) if len(c.strip()) >= width]
    clean = [c for c in cands
             if "$" not in c and "\\" not in c and not re.search(r"\d", c)]
    pool = clean or cands
    step = max(1, len(pool) // want)
    return [p[:width * 2] for p in pool[::step]][:want]


def check(path, quoted):
    """Is the quoted text present verbatim up to whitespace?"""
    f = ROOT / path
    if not f.exists():
        return None, f"target file missing: {path}"
    hay = norm(f.read_text(errors="replace"))
    return norm(quoted) in hay, None


if __name__ == "__main__":
    print("=" * 96)
    print("STALENESS RE-CHECK of every proposed passage against paper/ as it is today")
    print("=" * 96)
    tally = {"PENDING": 0, "APPLIED": 0, "STALE": 0, "PARTIAL": 0,
             "NO-ANCHOR": 0, "NO-TARGET": 0}
    failures, anchorless = [], []
    for rel in FILES:
        src = (ROOT / rel).read_text(errors="replace")
        print(f"\n  {rel}")
        print("  " + "-" * 92)
        for num, head, body in sections(src):
            targets = TARGET.findall(head)
            cur = current_block(body)
            label = f"§{num}"
            short = re.sub(r"[`*]", "", head)[:54]
            if not targets:
                tally["NO-TARGET"] += 1
                print(f"    {label:>5}  {'NO-TARGET':>9}  {short}")
                continue
            if cur is None:
                tally["NO-ANCHOR"] += 1
                anchorless.append(f"{rel} {label}")
                print(f"    {label:>5}  {'NO-ANCHOR':>9}  {short}"
                      f"   targets {', '.join(targets)}")
                continue
            prop = block(body, "Proposed")
            cur_in = any(check(t, cur)[0] for t in targets)
            prop_in, hits, total = False, 0, 0
            if prop:
                ps = probes(prop)
                total = len(ps)
                for t in targets:
                    f = ROOT / t
                    if not f.exists():
                        continue
                    hay = norm(f.read_text(errors="replace"))
                    h = sum(1 for pr in ps if pr in hay)
                    hits = max(hits, h)
                prop_in = hits > 0
            if cur_in and not prop_in:
                tally["PENDING"] += 1
                v = "PENDING"
            elif prop_in and not cur_in:
                tally["APPLIED"] += 1
                v = "APPLIED"
            elif cur_in and prop_in:
                tally["PARTIAL"] += 1
                v = "PARTIAL"
                failures.append(f"{rel} {label}: {short.strip()} (both texts present)")
            else:
                tally["STALE"] += 1
                v = "STALE"
                failures.append(f"{rel} {label}: {short.strip()}")
            ev = f"probes {hits}/{total}" if total else "no Proposed block"
            print(f"    {label:>5}  {v:>9}  {short}   {ev}")
    print()
    print("=" * 96)
    print("  " + "   ".join(f"{k} {v}" for k, v in tally.items()))
    print("=" * 96)
    if failures:
        print("\n  NEED ATTENTION (stale, or both texts present):")
        for f in failures:
            print(f"    {f}")
    if anchorless:
        print("\n  NO ANCHOR, so this check cannot clear them and a human must look:")
        for a in anchorless:
            print(f"    {a}")
    print("""
  PENDING means only that the text to be replaced is still there. It does not
  mean the passage is correct, nor that its numbers are wired. STALE, PARTIAL
  and NO-ANCHOR all block application.""")
    sys.exit(0)
