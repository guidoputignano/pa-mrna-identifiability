#!/usr/bin/env python3
"""
Every generated artefact, against what is committed.

The paper's tables and figures are produced by code and committed as files.
Nothing checked that the two agreed, and they stopped agreeing: the design
table was written by a function no entry point called, so it froze before the
seed work and printed single-seed numbers while the text beside it quoted
medians over 24 records (REVIEW_LOG.md R23). This is the check that would
have caught it the day it drifted.

Two modes, because regenerating a table costs seconds and regenerating the
decision characteristic costs minutes.

    python tools/check_artefacts.py
        Regenerate every table into a scratch directory and compare. About
        25 seconds, so `make test` runs it on every commit. Writes nothing
        into results/, because a failing check must not destroy the evidence
        it just found.

    python tools/check_artefacts.py --working-tree
        Compare everything in results/ against git HEAD, figures included.
        Run it AFTER `make figures`, which is what CI does. PDF creation and
        modification dates are stripped before comparing, because matplotlib
        stamps the current time into every PDF it writes and the whole tree
        would otherwise report as changed on every run -- which is its own
        hazard, since a dirty tree then says nothing about drift.

Exit status is 1 if anything differs, so either mode can gate a build.

One table is exempt from byte equality and compared number by number
instead, because its cells are simulated and a simulated number does not
reproduce bit for bit on another machine. What the exemption is, how narrow
it is, and what measurement sets its tolerance are all at SIMULATED below.
"""
import argparse
import io
import re
import subprocess
import sys
import tempfile
from contextlib import redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

#: matplotlib writes the wall clock into every PDF. Nothing about the drawing
#: changes when it does, so the comparison ignores these two fields.
PDF_DATES = re.compile(rb"/(CreationDate|ModDate)\s*\([^)]*\)")

#: Tables whose cells are Monte Carlo estimates, with the relative tolerance
#: their comparison allows.
#:
#: A simulated number is not byte-reproducible across machines. numpy's exp
#: takes a different SIMD path on a different CPU, which moves a least-squares
#: trajectory by an ulp, which flips the sign of a replicate that sat on the
#: decision boundary. Measured here by perturbing the synthetic record by a
#: relative 1e-13 -- a few ulp, the size of the effect -- a cell of the
#: decision table moves by up to 0.04 percentage points, which is enough to
#: print 2.9 where another machine prints 3.0. CI observed exactly that: two
#: runners on the SAME commit, minutes apart, one passing and one failing
#: (runs 35472135298 and 35472137934).
#:
#: So for these, and only these, the rule is weaker but still sharp: the first
#: column and every non-numeric character must match exactly, the numbers must
#: line up one for one, and each must agree to within the tolerance below.
#: Five per cent is well inside the seed-to-seed spread of every cell in this
#: table (6 to 11 per cent of its own value), so nothing meaningful about the
#: estimate escapes the check, and it is orders of magnitude inside the drift
#: the check exists to catch -- the design table froze at a factor of 20 off.
SIMULATED = {"tab_D_decision.tex": 0.05}

#: a signed decimal or scientific literal, wherever it sits in the markup
NUMBER = re.compile(r"-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


def cells_agree(old, new, rel):
    """
    True if two cells say the same thing and their numbers agree to `rel`.

    The text around the numbers has to be identical -- a changed word, a
    changed column, a number that appeared or vanished, all fail.
    """
    if NUMBER.sub("#", old) != NUMBER.sub("#", new):
        return False
    a = [float(x) for x in NUMBER.findall(old)]
    b = [float(x) for x in NUMBER.findall(new)]
    return all(abs(x - y) <= rel * max(abs(x), abs(y), 1e-30)
               for x, y in zip(a, b))


def simulated_match(old, new, rel):
    """
    Compare a table whose cells are simulated.

    Structure is compared exactly: the same lines, the same number of columns
    on each, the same row labels. Only the value columns get the tolerance.
    """
    a, b = old.splitlines(), new.splitlines()
    if len(a) != len(b):
        return False
    for x, y in zip(a, b):
        if x == y:
            continue
        fx, fy = x.split("&"), y.split("&")
        #  a line with no columns is rule, header or note: it must match
        if len(fx) != len(fy) or len(fx) == 1:
            return False
        #  the row label is not a simulated quantity
        if fx[0] != fy[0]:
            return False
        if not all(cells_agree(u, v, rel) for u, v in zip(fx[1:], fy[1:])):
            return False
    return True


def same_table(name, old, new):
    """Byte equality, unless `name` is a table this file declares simulated."""
    if old == new:
        return True
    rel = SIMULATED.get(name)
    return rel is not None and simulated_match(old, new, rel)


def tolerated(path, head):
    """
    True if a declared-simulated table differs only within its tolerance.

    `head` is the committed copy as raw bytes, which is how git hands it over.
    Anything that is not one of these tables, or does not decode as text, is
    not tolerated and falls through to the ordinary comparison.
    """
    rel = SIMULATED.get(path.name)
    if rel is None:
        return False
    try:
        return simulated_match(head.decode(), path.read_text(), rel)
    except UnicodeDecodeError:
        return False


def canonical(path_or_bytes):
    b = (path_or_bytes.read_bytes() if isinstance(path_or_bytes, Path)
         else path_or_bytes)
    return PDF_DATES.sub(b"", b)


def regenerate_tables(out):
    """Write every committed table into `out`, and return the names written."""
    from common.style import output_to
    import common.params as params
    import model_A.run as A
    import model_C.run as C
    import design.run as D

    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "n4", ROOT / "audit" / "n4_theta_scope.py")
    n4 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(n4)

    with output_to(out), redirect_stdout(io.StringIO()):
        params.write_params_table()
        params.write_params_table(compact=True)
        C.table_shape()
        A.table_design()
        D.table_designs()
        D.table_excursion()
        #  the quasi-steady check contributes the settling time a protocol
        #  needs, so it runs here too or the numbers file comes out short
        D.quasi_steady_check()
        #  the dose-train figure contributes the increment's peak and trough
        #  and the phase coverage of each schedule, which the captions quote
        D.fig_dose_train()
        D.table_factors()
        D.table_load_precision()
        D.table_instruments()
        D.table_spec()
        #  the scope convention lives in the audit script that checks it,
        #  so the tables that use it are written there rather than copied
        n4.table_elasticity()
        n4.table_scope()
        n4.table_ke_prior()
        #  and the numbers the text quotes in prose, which drift the same way
        from common.numbers import write_numbers, NUMDIR
        import common.numbers as _num
        _num.NUMDIR = Path(out) / "numbers"
        write_numbers("design", D.PAPER_NUMBERS, D.PAPER_FORMATS)
        _num.NUMDIR = NUMDIR
    return sorted(p.name for p in (Path(out) / "tables").glob("*.tex"))


def check_regenerated(live=None):
    """
    Compare the regenerated tables with the committed ones.

    `live` overrides the directory holding the committed copies, which is
    what the guard's own self-test uses: it must be able to perturb a table
    and see the check reject it WITHOUT editing anything in results/, or a
    process killed mid-test leaves a corrupted number in the repository.
    """
    live = Path(live) if live else ROOT / "results" / "tables"
    with tempfile.TemporaryDirectory() as tmp:
        written = regenerate_tables(tmp)
        bad, checked = [], 0
        for name in written:
            new = Path(tmp) / "tables" / name
            old = live / name
            if not old.exists():
                bad.append((name, "regenerated but not committed"))
                continue
            checked += 1
            if not same_table(name, old.read_text(), new.read_text()):
                a, b = old.read_text().splitlines(), new.read_text().splitlines()
                n = sum(1 for x, y in zip(a, b) if x != y) + abs(len(a) - len(b))
                bad.append((name, f"differs on {n} line(s)"))
        #  the generated numbers file, compared the same way
        nums = Path(tmp) / "numbers" / "design.tex"
        live = ROOT / "results" / "numbers" / "design.tex"
        if not nums.exists():
            bad.append(("numbers/design.tex", "not regenerated"))
        elif not live.exists():
            bad.append(("numbers/design.tex", "regenerated but not committed"))
        else:
            checked += 1
            if live.read_text() != nums.read_text():
                bad.append(("numbers/design.tex", "differs"))
        print(f"  {checked} generated tables and number files checked "
              f"against the code")
        # tab_claims is rendered from paper/claims.csv by claims_to_tex.py,
        # which `make refs` runs before every paper build, so it cannot go
        # stale the way a figure module's table can.
        for name, why in bad:
            print(f"    {name}: {why}")
        return not bad


def check_working_tree():
    """
    Compare results/ against HEAD after a regeneration.

    Tables and the numbers file are text the code writes, so any difference
    is drift and fails.

    Figures are not. A .pdf or .png is a rendering, and matplotlib does not
    render byte-identically across versions: a runner with a newer
    matplotlib, or different fonts, produces different bytes from identical
    data. Failing on that would report library drift as artefact drift, and
    the first thing it did was exactly that. So a figure difference is
    reported and not fatal; what IS fatal for a figure is being missing or
    empty, which means its generator did not run.

    What this does not catch is a figure whose DATA changed while its table
    did not. Nothing cheap does: the fix for that is to make the number the
    figure turns on a generated number, where the guard above sees it.
    """
    out = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain",
                          "results"], capture_output=True, text=True).stdout
    paths = [line[3:].strip() for line in out.splitlines() if line.strip()]
    real, dates_only, rendered, simulated = [], [], [], []
    for rel in paths:
        p = ROOT / rel
        is_figure = "/figures/" in rel
        if not p.exists():
            real.append((rel, "deleted"))
            continue
        head = subprocess.run(["git", "-C", str(ROOT), "show", f"HEAD:{rel}"],
                              capture_output=True).stdout
        if not head:
            (rendered if is_figure else real).append(
                (rel, "not committed"))
        elif canonical(p) == canonical(head):
            dates_only.append(rel)
        elif tolerated(p, head):
            simulated.append(rel)
        elif is_figure:
            rendered.append((rel, "renders differently"))
        else:
            real.append((rel, "content differs"))
    missing = []
    for stem in sorted({p.stem for p in (ROOT / "results/figures").glob("*.pdf")}):
        f = ROOT / "results/figures" / f"{stem}.pdf"
        if not f.exists() or f.stat().st_size == 0:
            missing.append((str(f.relative_to(ROOT)), "missing or empty"))
    real += missing
    print(f"  {len(paths)} artefact(s) differ from HEAD by bytes;"
          f" {len(dates_only)} by their embedded dates alone")
    for rel in simulated:
        print(f"  {rel}: simulated cells within tolerance, structure "
              f"unchanged (see SIMULATED)")
    if rendered:
        print(f"  {len(rendered)} figure(s) render differently, which is not"
              f" checked here (see this function's docstring):")
        for rel, why in rendered[:8]:
            print(f"    {rel}: {why}")
        if len(rendered) > 8:
            print(f"    ... and {len(rendered) - 8} more")
    for rel, why in real:
        print(f"    {rel}: {why}")
    return not real


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--working-tree", action="store_true",
                    help="compare results/ against HEAD instead of "
                         "regenerating the tables")
    args = ap.parse_args()
    ok = check_working_tree() if args.working_tree else check_regenerated()
    print("  ARTEFACTS MATCH THE CODE" if ok else "  ARTEFACTS HAVE DRIFTED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
