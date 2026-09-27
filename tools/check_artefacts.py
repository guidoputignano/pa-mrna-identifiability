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
import ast
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

#: one generated macro line, as common/numbers.py writes it
_MACRO_LINE = re.compile(r"\\csname paper@num@([A-Za-z0-9_]+)\\endcsname\{(.*)\}$")

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
        D.quasi_steady_check(tol=0.01)
        #  the dose-train figure contributes the increment's peak and trough
        #  and the phase coverage of each schedule, which the captions quote
        D.fig_dose_train()
        #  the two figures that contribute generated numbers rather than a
        #  table: the excitation legend's E for each band, and the protocol
        #  figure's range bars, which write the SAME delta_range_* keys
        #  tab_D_designs does. Without them the regenerated design.tex comes
        #  out short and the comparison fails on the guard's own fast path.
        D.fig_excitation()
        D.fig_protocol()
        D.table_factors()
        D.table_load_precision()
        D.table_instruments()
        D.table_spec()
        #  the scope convention lives in the audit script that checks it,
        #  so the tables that use it are written there rather than copied
        n4.table_scope()
        n4.table_ke_prior()
        #  and the numbers the text quotes in prose, which drift the same way
        from common.numbers import write_numbers, NUMDIR
        import common.numbers as _num
        _num.NUMDIR = Path(out) / "numbers"
        write_numbers("design", D.PAPER_NUMBERS, D.PAPER_FORMATS)
        #  the scope result: the band and the enzymatic share, which the
        #  abstract, Section 2.5 and Section 6.2 all quote
        n4.scope_numbers()
        #  the K_m sweep, quoted in two sections; 0.6 s
        spec_n1 = importlib.util.spec_from_file_location(
            "n1", ROOT / "audit" / "n1_km_cost.py")
        n1 = importlib.util.module_from_spec(spec_n1)
        spec_n1.loader.exec_module(n1)
        n1.km_numbers()
        #  the reduction's off-equilibrium cost, which the k row of
        #  the constants registry quotes as that constant's bound
        C.reduction_cost_numbers()
        #  the reduction's off-equilibrium cost, which the k row
        #  of the constants registry quotes as that constant's bound
        C.reduction_cost_numbers()
        #  the m_tol window, which the registry's m_tol row quotes as that
        #  convention's cost; it is produced by the fibre-margin figure
        A.fig_margin_fibre()
        #  the fibre figures' two asserted numbers, which are cheap; the
        #  margin figure's four are not (48 s), so results/numbers/margin.tex
        #  is covered by --working-tree after `make figures` instead, as
        #  decision.tex and profile.tex are.
        A.fig_load_error()
        A.fig_fibre()
        A.fibre_numbers()
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
        #  decision.tex and profile.tex are written by the slow generators and
        #  are covered by --working-tree after `make figures`, not here
        for stem in ("design", "scope", "kmcost", "local", "reduction",
                     "fibre"):
            nums = Path(tmp) / "numbers" / f"{stem}.tex"
            here = ROOT / "results" / "numbers" / f"{stem}.tex"
            if not nums.exists():
                bad.append((f"numbers/{stem}.tex", "not regenerated"))
            elif not here.exists():
                bad.append((f"numbers/{stem}.tex",
                            "regenerated but not committed"))
            else:
                checked += 1
                if here.read_text() != nums.read_text():
                    bad.append((f"numbers/{stem}.tex", "differs"))
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


def _bodies(text):
    """{function name: its source text}, for module-level defs only."""
    out, name, buf = {}, None, []
    for line in text.splitlines(keepends=True):
        m = re.match(r"def (\w+)\(", line)
        if m:
            if name:
                out[name] = "".join(buf)
            name, buf = m.group(1), [line]
        elif name is not None:
            if line and not line[0].isspace() and line.strip():
                out[name] = "".join(buf)
                name, buf = None, []
            else:
                buf.append(line)
    if name:
        out[name] = "".join(buf)
    return out


def _reachable_from_main(text, main_body):
    """
    Every module-level function the __main__ block reaches, transitively.

    A literal scan of the __main__ block is not enough: n1's km_numbers() and
    n4's scope_numbers() are called from main(), and design/run.py's
    spec_numbers() from table_spec(). Those are reachable and the first
    version of this check called all three broken. What it must catch is a
    generator NOTHING in the build reaches, which is a graph question.
    """
    bodies = _bodies(text)
    seen, stack = set(), [fn for fn in bodies if f"{fn}(" in main_body]
    while stack:
        fn = stack.pop()
        if fn in seen:
            continue
        seen.add(fn)
        body = bodies.get(fn, "")
        stack.extend(g for g in bodies if g not in seen and f"{g}(" in body)
    return seen


def check_generators_reachable():
    """
    Every *_numbers() function must be reachable from `make figures`.

    THREE INSTANCES OF ONE BUG, found in a single session and all invisible:

      * src/design/run.py defined quasi_steady_check() and its __main__ never
        called it, so `make figures` wrote a design.tex five keys short;
      * audit/t1_transient_worth.py defined transient_numbers() and nothing
        called it, so the build ran the file and left transient.tex untouched;
      * that file was not in the Makefile at all.

    All three survived because THIS GUARD'S FAST PATH called the generators
    directly while the build did not, so both paths were green and they were
    green about different things. A generator the build cannot reach is a
    generator that cannot catch drift, which is the one job it has.

    The check is static: find every `def *_numbers()` and every module that
    writes a numbers file, and require each to be named by a Makefile recipe
    line AND to be called from that module's own `__main__` block.
    """
    mk = (ROOT / "Makefile").read_text()
    figures = mk.split("figures:", 1)[1].split("\n\n", 1)[0] if "figures:" in mk else ""
    problems = []
    seen = 0
    for py in sorted(list((ROOT / "src").rglob("*.py"))
                     + list((ROOT / "audit").glob("*.py"))):
        text = py.read_text()
        if "write_numbers(" not in text:
            continue
        rel = str(py.relative_to(ROOT))
        if py.name == "numbers.py":
            continue
        seen += 1
        if rel not in figures:
            problems.append(f"{rel} writes a numbers file and no `make "
                            f"figures` recipe line runs it")
            continue
        main = text.split('if __name__ == "__main__":', 1)
        if len(main) < 2:
            problems.append(f"{rel} writes a numbers file and has no "
                            f"__main__ block, so the build cannot reach it")
            continue
        reached = _reachable_from_main(text, main[1])
        for fn in re.findall(r"^def (\w*numbers)\(", text, re.M):
            if fn not in reached:
                problems.append(f"{rel} defines {fn}() and nothing "
                                f"__main__ runs ever calls it")
        #  a module whose numbers come from helpers rather than a
        #  *_numbers() function: require the write itself to be reachable
        if not re.search(r"^def \w*numbers\(", text, re.M):
            owner = [fn for fn in _bodies(text)
                     if "write_numbers(" in _bodies(text)[fn]]
            if owner and not any(fn in reached for fn in owner):
                problems.append(f"{rel} calls write_numbers() only from "
                                f"{', '.join(owner)}, which __main__ never "
                                f"reaches")
    #  the same standard C2 of audit/n2_table2_recheck.py holds itself to: a
    #  check that cannot be shown to fail is not evidence. These are the three
    #  shapes actually found, written out.
    probe_ok = []
    probe = ('def foo_numbers():\n    write_numbers("x", {})\n\n'
             'def main():\n    print(1)\n\n'
             'if __name__ == "__main__":\n    main()\n')
    probe_ok.append("foo_numbers" not in
                    _reachable_from_main(probe, probe.split(
                        'if __name__ == "__main__":', 1)[1]))
    probe2 = probe.replace("    print(1)", "    foo_numbers()")
    probe_ok.append("foo_numbers" in
                    _reachable_from_main(probe2, probe2.split(
                        'if __name__ == "__main__":', 1)[1]))
    probe3 = ('def foo_numbers():\n    write_numbers("x", {})\n\n'
              'if __name__ == "__main__":\n    foo_numbers()\n')
    probe_ok.append("foo_numbers" in
                    _reachable_from_main(probe3, probe3.split(
                        'if __name__ == "__main__":', 1)[1]))
    if not all(probe_ok):
        problems.append("the reachability check's own self-test failed: "
                        f"{probe_ok}")
    print(f"  {seen} numbers-file generator(s) checked for reachability from "
          f"`make figures` (self-test: "
          f"{'PASS' if all(probe_ok) else 'FAIL'})")
    for why in problems:
        print(f"    {why}")
    return not problems



#: matplotlib calls that put a string on a figure. `set_xticklabels` and
#: `set_yticklabels` are deliberately NOT here: a tick is a coordinate, not a
#: claim, and the rule below is about claims. They are counted and reported so
#: the exclusion is visible rather than assumed.
DRAWS_TEXT = {"set_title", "suptitle", "text", "figtext", "annotate",
              "set_xlabel", "set_ylabel", "supxlabel", "supylabel",
              "set_label", "bar_label"}
TICK_TEXT = {"set_xticklabels", "set_yticklabels"}

#: the classifications a module may give a number it draws. Everything drawn
#: must carry one of these, and `registered:` is the only one the machine can
#: check; the rest are declarations a person made and signed.
DRAWN_KINDS = ("registered:", "own:", "input:", "geometry:", "text:",
               "UNREGISTERED:")


def _text_helpers():
    """
    This repository's own text-drawing wrappers, found rather than listed.

    corner(), over_x(), beside_y() and at_point() all end in ax.text or
    ax.annotate, so a number passed to one of them is drawn exactly as if it
    had been passed to matplotlib. A hardcoded list of matplotlib names misses
    every one of them, and did: the first version of this check saw 23 drawn
    numbers and the four helpers were hiding 8 more, among them theta_crit and
    the safety margin, which are two of the figures' load-bearing claims.
    """
    sty = (ROOT / "src" / "common" / "style.py").read_text()
    return set(re.findall(r"^def (\w+)\(ax[^)]*\btext\b", sty, re.M))


def _interpolations(node):
    """
    Every expression whose VALUE is formatted into a string inside `node`.

    Covers the three ways this repository builds text: f-strings, `%` on a
    literal (which design/run.py uses for the dose-train trough) and
    "...".format(). A number reaching a drawn string by any other route is
    outside this check; see check_figure_numbers.__doc__.
    """
    out = []
    for sub in ast.walk(node):
        if isinstance(sub, ast.JoinedStr):
            out += [ast.unparse(v.value) for v in sub.values
                    if isinstance(v, ast.FormattedValue)]
        elif (isinstance(sub, ast.BinOp) and isinstance(sub.op, ast.Mod)
              and isinstance(sub.left, ast.Constant)
              and isinstance(sub.left.value, str)):
            r = sub.right
            out += [ast.unparse(e)
                    for e in (r.elts if isinstance(r, ast.Tuple) else [r])]
        elif (isinstance(sub, ast.Call) and isinstance(sub.func, ast.Attribute)
              and sub.func.attr == "format"
              and isinstance(sub.func.value, ast.Constant)):
            out += [ast.unparse(a) for a in sub.args]
    return out


def _drawn(text, helpers):
    """
    {expression source: first line it is drawn on}, plus the tick-label ones.

    A `label=` keyword on any call counts: a legend is drawn text, and
    design/run.py draws the excitation number E in one.
    """
    drawn, ticks = {}, {}
    for node in ast.walk(ast.parse(text)):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        nm = (fn.attr if isinstance(fn, ast.Attribute) else
              fn.id if isinstance(fn, ast.Name) else None)
        into = (drawn if nm in DRAWS_TEXT or nm in helpers else
                ticks if nm in TICK_TEXT else None)
        if into is not None:
            for e in _interpolations(node):
                into.setdefault(e, node.lineno)
        for kw in node.keywords:
            if kw.arg == "label":
                for e in _interpolations(kw.value):
                    drawn.setdefault(e, node.lineno)
    return drawn, ticks


def _registered_exprs(text):
    """
    The source text of every expression this module registers as a number.

    Collected from the two ways a value reaches results/numbers/: a dict
    literal handed straight to write_numbers(), and a `.update({...})` on the
    module's accumulating numbers dict, which is how design/run.py and
    model_A/run.py build theirs across a dozen functions.
    """
    tree = ast.parse(text)
    books, out = set(), set()
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                and node.func.id == "write_numbers"):
            if len(node.args) > 1 and isinstance(node.args[1], ast.Name):
                books.add(node.args[1].id)
            for a in node.args[1:2]:
                if isinstance(a, ast.Dict):
                    out |= {ast.unparse(v) for v in a.values}
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "update"
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id in books
                and node.args and isinstance(node.args[0], ast.Dict)):
            out |= {ast.unparse(v) for v in node.args[0].values}
    return books, out


def _registers(text, expr):
    """
    True if `expr` draws a registered number rather than recomputing one.

    Three shapes count, and all three are checked rather than declared.

      num("design", "delta_range_off")   reads the committed file. This is
                                         the cross-module case: the figure
                                         draws what \num{} typesets.
      N["key"] where write_numbers(_, N) the drawn value and the registered
                                         value are one object.
      the same expression registered      `trough` is drawn by the dose-train
      elsewhere in this module            figure and registered four lines
                                         later as train_trough_frac. One
                                         computation, two uses.

    The third is text identity and is therefore narrow on purpose. A figure
    that draws np.interp(M_TOL, mts, wid) while the register holds
    margin_window(a0, b0, 40.0) is drawing a DIFFERENT computation of the
    same quantity, which is exactly the drift this check is for, so it does
    not pass -- it has to be made one computation first.
    """
    if re.search(r"\bnum\(", expr):
        return True
    books, regd = _registered_exprs(text)
    if expr in regd:
        return True
    m = re.match(r"^(\w+)\[", expr)
    return bool(m) and m.group(1) in books


def _built_elsewhere(text, helpers):
    """
    How many interpolated strings this module builds away from a drawing call.

    THE LIMIT OF THE CHECK, counted rather than assumed. A number reaches
    drawn text two ways: formatted inside the call, which _drawn() sees, or
    formatted into a name and drawn later, which it cannot see without
    following the value. g5_conditions.py builds its legend strings into a
    list first; so, harmlessly, does every module that formats a table cell.

    Console output and table cells are subtracted, because those have their
    own guards -- stdout is not an artefact and a table cell is compared
    line by line by check_regenerated(). What is left is the residue this
    check does not cover, and it is printed so the gap is a number.
    """
    tree = ast.parse(text)
    drawn_at, _ = _drawn(text, helpers)
    covered, total = set(), []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        nm = (fn.attr if isinstance(fn, ast.Attribute) else
              fn.id if isinstance(fn, ast.Name) else None)
        if nm in ("print", "write_numbers", "write_table", "append", "update",
                  "format", "join", "warn") or nm in DRAWS_TEXT \
                or nm in TICK_TEXT or nm in helpers:
            for sub in ast.walk(node):
                if isinstance(sub, ast.JoinedStr):
                    covered.add(id(sub))
        for kw in node.keywords:
            if kw.arg == "label":
                for sub in ast.walk(kw.value):
                    if isinstance(sub, ast.JoinedStr):
                        covered.add(id(sub))
    for node in ast.walk(tree):
        if (isinstance(node, ast.JoinedStr) and id(node) not in covered
                and any(isinstance(v, ast.FormattedValue)
                        for v in node.values)):
            total.append(node)
    return len(total)


def _collisions(sources):
    """{key: (stem, value)} and the collisions, from {stem: file text}."""
    seen, problems = {}, []
    for stem in sorted(sources):
        for line in sources[stem].splitlines():
            m = _MACRO_LINE.search(line)
            if not m:
                continue
            key, val = m.group(1), m.group(2)
            if key in seen:
                other, oval = seen[key]
                problems.append(
                    f"{key!r} is defined in both {other}.tex and {stem}.tex "
                    f"({oval!r} against {val!r}); the preamble's input order "
                    f"would decide which the paper uses")
            else:
                seen[key] = (stem, val)
    return seen, problems


def check_key_collisions():
    r"""
    No generated key may be defined in two numbers files.

    THE ONE MECHANICAL CASE of this repository's characteristic failure, a
    number that does not say what it is a number of (REVIEW_LOG.md). Most of
    that class is semantic and no check can see it: whether "the managed
    level" means 10 or 81.93 is a question about a phrase. But its crudest
    form IS visible, and is the one that would be worst if it happened.

    The preamble \input{}s every file in results/numbers/. A key defined in
    two of them does not collide, error or warn: the LAST \input silently
    wins, and every \num{} of that key in the paper takes whichever file the
    preamble happened to list second. Two different quantities would then be
    one macro, and the build would be green. There is no collision today; this
    is here so that there is still none the day a second module wants a name
    like theta_crit, which margin.tex already owns and which boundary.tex
    came within one key of taking.
    """
    numdir = ROOT / "results" / "numbers"
    sources = {f.stem: f.read_text() for f in sorted(numdir.glob("*.tex"))}
    seen, problems = _collisions(sources)
    #  the same standard the other checks hold themselves to: a check that
    #  cannot be shown to fail is not evidence. The detector is run against a
    #  synthetic pair that DOES collide and a pair that does not, so the
    #  green result above is a result and not an absence.
    line = r"\expandafter\def\csname paper@num@%s\endcsname{%s}"
    clash = _collisions({"a": line % ("k", "1"), "b": line % ("k", "2")})[1]
    clean = _collisions({"a": line % ("k", "1"), "b": line % ("j", "2")})[1]
    probe_ok = len(clash) == 1 and not clean
    print(f"  {len(seen)} generated key(s) checked for cross-file collisions "
          f"across {len(list(numdir.glob('*.tex')))} numbers files "
          f"(self-test: {'PASS' if probe_ok else 'FAIL'})")
    for why in problems:
        print(f"    {why}")
    return not problems


#: value coincidences already looked at and confirmed as distinct quantities.
#: Each needs a reason, because the point of the check is that somebody read it.
DECLARED_VALUE_TWINS = {
    #  --- confirmed DISTINCT quantities that happen to coincide ---
    ("sd_tier1_with_kappa", "sim_spread_recommended"):
        "Tier 1's marginal and the recommended design's simulated spread both "
        "land on 0.2002. Different designs, different objects, coincident to "
        "four places (REVIEW_LOG.md R145).",
    ("E_info", "bias_alone20"):
        "An excitation measure and a load-bias offset, both 0.172. Unrelated "
        "quantities from different generators; coincidence at three places.",
    ("bias_sweep", "trans_noprior_nominal"):
        "A load-bias sweep endpoint and the transient's prior-free nominal, "
        "both 0.1733. Unrelated quantities from different generators.",

    #  --- confirmed REDUNDANT. Declared so the build is green while the
    #  finding stays visible in the code, NOT because they are acceptable.
    #  Collapsing them changes generated output and possibly the paper, which
    #  is the author's call; recorded at REVIEW_LOG.md R147 and in
    #  OPEN_QUESTIONS.md as an open item. ---
    ("bias_both20", "bias_loose"):
        "REDUNDANT, PENDING COLLAPSE. audit/t4_bias_vs_scatter.py assigns BOTH "
        "keys the same variable `both` (lines 114 and 119), so these are one "
        "number under two names and can be edited apart while this guard "
        "reports a match. Exactly the class this check exists for, found on "
        "its first run, in code predating it.",
    ("sd_load_1pct", "sd_load_steady"):
        "DIFFERENT INPUTS, SAME ANSWER, AND THAT IS THE POINT. A load record "
        "good to 1% and one with no variability at all both give 0.1696, "
        "which is the insensitivity src/design/run.py:1686 states by "
        "contrasting 1% against 20%. Coincidence here is the result.",
    ("sd_load_1pct", "trans_base_nominal"):
        "Same design at nominal capacity reached by two generators, 0.1696. "
        "Distinct scripts, one underlying quantity; see the entry above.",
    ("sd_load_steady", "trans_base_nominal"):
        "Same design at nominal capacity reached by two generators, 0.1696.",
    ("loadcap_perfect_5pct", "meas_exact_5pct"):
        "REDUNDANT, PENDING COLLAPSE. audit/t8_load_by_capacity.py computes "
        "over(150.0, 1e-6, 1e-6) and audit/t9_measured_load.py computes "
        "over(150.0, measured=True, prec=1e-6), both 2.659: the same quantity "
        "from two generators. PROVENANCE_SWEEP.md separately records 2.659 as "
        "appearing in three places in the paper and 'owned by nothing in any "
        "of them', so this is that value seen from the generator side.",
}


def _value_twins(sources, min_decimals=3):
    """
    Distinct keys whose emitted values are identical to full precision.

    WHY THIS IS NOT THE COLLISION CHECK ABOVE. That one catches one key in two
    files, where the preamble silently picks a winner. This catches the
    opposite and quieter case: TWO keys carrying the SAME value. The guard
    compares keys, so a redundant key passes it, and two names for one number
    can be edited apart later while the guard still reports a match. That
    happened: `sd_recommended_postdose` was emitted beside the existing
    `sd_monthly`, both 0.1699 by different code, and it was caught by reading
    rather than by any check (REVIEW_LOG.md R146).

    A shared value is not itself wrong. `0.20` and `25` coincide all over a
    paper and mean nothing by it, which is why only values carrying at least
    `min_decimals` decimals are considered: coincidence at four places is
    worth a human's glance, coincidence at two is arithmetic.

    AND THE CHECK MUST NOT TRY TO DECIDE, because the two cases it separates
    look identical from here. `sd_load_1pct` and `sd_load_steady` agree to four
    places, and they are a load record good to 1% against one with no
    variability at all: TWO DIFFERENT INPUTS GIVING THE SAME ANSWER, which is
    an invariance `src/design/run.py:1686` already states by contrasting 1%
    against 20%. A redundant key and a real invariance are the same bytes in
    this file. So the check reports both and adjudicates neither, and the
    declaration carries a human's reason: for a redundancy, that it is pending
    collapse; for an invariance, that the coincidence IS the result. A check
    that guessed between them would suppress findings half the time and
    manufacture them the other half.
    """
    def _stem(k):
        """A median and its seed range are one quantity under three keys."""
        for suf in ("_hi", "_lo"):
            if k.endswith(suf):
                return k[: -len(suf)]
        return k

    vals = {}
    for stem in sorted(sources):
        for line in sources[stem].splitlines():
            m = _MACRO_LINE.search(line)
            if not m:
                continue
            key, val = m.group(1), m.group(2)
            if "." not in val:
                continue
            if len(val.split(".")[-1]) < min_decimals:
                continue
            vals.setdefault(val, []).append((key, stem))
    twins = []
    for val, ks in sorted(vals.items()):
        if len(ks) < 2:
            continue
        for i in range(len(ks)):
            for j in range(i + 1, len(ks)):
                a, b = sorted([ks[i][0], ks[j][0]])
                #  A MEDIAN AND ITS RANGE ARE NOT A DUPLICATE. `x`, `x_hi` and
                #  `x_lo` carry one quantity and its seed spread, and a tight
                #  spread makes all three print the same digits. The first
                #  version of this check reported those as twins and buried
                #  its real finds under them: 13 reports, of which 9 were one
                #  median-and-range cluster. Characterising the false-positive
                #  mode rather than handing over the noisy list is the same
                #  lesson as u21's SUSPECT column (REVIEW_LOG.md R141).
                if _stem(a) == _stem(b):
                    continue
                #  declarations are by STEM, so one entry covers a key and
                #  its _hi/_lo range variants rather than needing three
                if tuple(sorted((_stem(a), _stem(b)))) in DECLARED_VALUE_TWINS:
                    continue
                twins.append((a, ks[i][1], b, ks[j][1], val))
    return twins


def check_value_twins():
    r"""
    Two keys carrying the same value must be declared, or collapsed.

    A duplicate that AGREES today is the only kind that can silently diverge
    tomorrow, which is why this fails rather than warns: confirming a pair
    means adding it to DECLARED_VALUE_TWINS with a reason, and that is a
    thirty-second act that leaves a record.
    """
    numdir = ROOT / "results" / "numbers"
    sources = {f.stem: f.read_text() for f in sorted(numdir.glob("*.tex"))}
    twins = _value_twins(sources)
    line = r"\expandafter\def\csname paper@num@%s\endcsname{%s}"
    probe_hit = _value_twins({"a": line % ("k", "0.1234"),
                              "b": line % ("j", "0.1234")})
    probe_miss = _value_twins({"a": line % ("k", "0.20"),
                               "b": line % ("j", "0.20")})
    probe_ok = len(probe_hit) == 1 and not probe_miss
    print(f"  {len(sources)} numbers file(s) checked for undeclared value twins, "
          f"{len(DECLARED_VALUE_TWINS)} declared "
          f"(self-test: {'PASS' if probe_ok else 'FAIL'})")
    for a, sa, b, sb, val in twins:
        print(f"    {a!r} ({sa}.tex) and {b!r} ({sb}.tex) are both {val}: "
              f"collapse them, or declare the pair with a reason")
    return not twins and probe_ok


def check_numbers_loaded():
    r"""
    Every generated numbers file must be \input by the document, and every
    \num{} key the paper uses must live in a file that IS.

    THE GAP THIS CLOSES, AND IT WAS OPEN THE WHOLE TIME (REVIEW_LOG.md R152).
    The other checks establish that a key is produced by a generator the build
    reaches, and that no two files define one key. **None of them asks whether
    the preamble loads the file the key lives in.** It did not, for six files:
    `bias`, `tier1`, `loadcap`, `diet`, `measured` and `frontier`. Twenty-four
    keys the paper prints were therefore undefined at compile time and raised
    \PackageError, twenty-nine times, while this script reported
    ARTEFACTS MATCH THE CODE.

    **It is the exact failure mode the guard exists to prevent, one level up.**
    A number with no generator is what `check_regenerated` catches. A number
    with a generator nobody loads is invisible to every check here, because
    each one stops at the generated file and never asks whether LaTeX reads it.
    `frontier` is Phase 11's own wiring, so the pass that added nineteen
    macros shipped them into a file the document does not open, and the guard
    certified that too.

    Two directions, because they fail differently. An unloaded file whose keys
    the paper never cites is silent today and a trap tomorrow; an unloaded file
    whose keys it does cite is already broken. Both are reported.
    """
    numdir = ROOT / "results" / "numbers"
    generated = {f.stem for f in sorted(numdir.glob("*.tex"))}
    #  which files the document actually opens, read from every .tex in paper/
    loaded = set()
    for tex in sorted((ROOT / "paper").rglob("*.tex")):
        for m in re.finditer(r"\\input\{[^}]*results/numbers/([A-Za-z0-9_]+)\}",
                             tex.read_text(errors="replace")):
            loaded.add(m.group(1))
    unloaded = sorted(generated - loaded)
    missing_files = sorted(loaded - generated)

    #  which keys the paper cites, and where each one lives
    cited = set()
    for tex in sorted((ROOT / "paper").rglob("*.tex")):
        cited |= set(re.findall(r"\\num\{([A-Za-z0-9_]+)\}",
                                tex.read_text(errors="replace")))
    home = {}
    for stem in generated:
        for line in (numdir / f"{stem}.tex").read_text(errors="replace").splitlines():
            m = _MACRO_LINE.search(line)
            if m:
                home[m.group(1)] = stem
    undefined = sorted(k for k in cited
                       if home.get(k) is not None and home[k] not in loaded)

    print(f"  {len(generated)} numbers file(s) generated, {len(loaded)} loaded "
          f"by the document, {len(cited)} \\num{{}} key(s) cited")
    for stem in unloaded:
        n = sum(1 for k, v in home.items() if v == stem)
        cites = sum(1 for k in cited if home.get(k) == stem)
        print(f"    results/numbers/{stem}.tex is GENERATED BUT NEVER \\input "
              f"({n} keys, {cites} of them cited by the paper): add it to "
              f"paper/preamble.tex")
    for stem in missing_files:
        print(f"    paper/ inputs results/numbers/{stem}, which no generator "
              f"writes")
    for k in undefined:
        print(f"    \\num{{{k}}} is cited but lives in the unloaded "
              f"{home[k]}.tex, so it raises \\PackageError at compile time")
    return not unloaded and not missing_files and not undefined


#: generated artefacts the document deliberately does not use, each with the
#: reason a human wrote. Same discipline as DECLARED_VALUE_TWINS: an unused
#: artefact is not always a defect, but the difference is a judgement and the
#: judgement belongs in the file rather than in whoever last looked.
DECLARED_UNUSED_ARTEFACTS = {
    "tab_S_admissible":
        "The FULL admissibility table, three columns including an EVIDENCE "
        "column. The supplement inputs tab_S_admissible_compact instead, which "
        "fits two constraint/verdict pairs across and drops that column. Both "
        "come from audit/h1_admissibility.py:440,450. Kept generated because "
        "the evidence column is the part a reader would want and the compact "
        "table is the space-constrained choice, not the better one. OPEN FOR "
        "THE AUTHOR: whether the supplement should carry the full version.",
}


def check_generated_artefacts_used():
    r"""
    The same question as `check_numbers_loaded`, asked of the other two things
    `results/` generates: is every generated TABLE \input somewhere, and is
    every generated FIGURE actually included.

    WHY THIS EXISTS, and it is the author's instruction after R152: the defect
    found there is one file type over. `check_inputs_reachable` asks "does the
    document read this?" of section files, `check_numbers_loaded` asks it of
    numbers files, and nothing asked it of tables or figures. **A generated
    artefact nobody reads is the same failure whatever its extension**, and the
    reason R152 went unnoticed for weeks is that the question was being asked
    of one artefact class at a time, as each one happened to cause trouble,
    rather than of all of them at once.

    Both directions, because they fail differently and only one is loud.
      GENERATED BUT UNUSED  costs a build step and misleads a reader of
        `make figures` into thinking the paper shows it. Silent.
      USED BUT NOT GENERATED  is a broken `\input` or a missing image, which
        LaTeX reports, but only once someone compiles. Since nothing in this
        repository compiled the paper for the whole of this audit, "LaTeX would
        catch it" was not a defence (R152).

    Figures are matched on basename without extension, because the document
    cites `fig_x.pdf` while `make figures` writes both `.pdf` and `.png`, and
    a check keyed on the exact filename would report every `.png` as unused.
    """
    texts = {}
    for tex in sorted((ROOT / "paper").rglob("*.tex")):
        texts[tex] = tex.read_text(errors="replace")
    blob = "\n".join(texts.values())

    #  ---- tables -------------------------------------------------------
    tdir = ROOT / "results" / "tables"
    gen_t = {f.stem for f in sorted(tdir.glob("*.tex"))} if tdir.is_dir() else set()
    used_t = set(re.findall(r"\\input\{[^}]*tables/([A-Za-z0-9_]+)\}", blob))
    t_unused = sorted(gen_t - used_t - set(DECLARED_UNUSED_ARTEFACTS))
    t_missing = sorted(used_t - gen_t)

    #  ---- figures ------------------------------------------------------
    fdir = ROOT / "results" / "figures"
    gen_f = ({f.stem for f in fdir.glob("*.pdf")} |
             {f.stem for f in fdir.glob("*.png")}) if fdir.is_dir() else set()
    used_f = {Path(m).stem for m in
              re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]*)\}", blob)}
    #  only judge names that look like this repository's generated artefacts,
    #  so a logo or an externally supplied image is not reported as missing
    f_unused = sorted(n for n in gen_f - used_f
                      if n.startswith(("fig_", "tab_"))
                      and n not in DECLARED_UNUSED_ARTEFACTS)
    f_missing = sorted(n for n in used_f - gen_f
                       if n.startswith(("fig_", "tab_")))

    #  a self-test, because a check that silently stopped checking is the
    #  failure mode R152 was: feed it a name nothing uses and it must report it
    probe = sorted({"tab_probe_unused"} - used_t - set(DECLARED_UNUSED_ARTEFACTS))
    print(f"  {len(DECLARED_UNUSED_ARTEFACTS)} declared unused "
          f"(self-test: {'PASS' if probe else 'FAIL'})")
    print(f"  {len(gen_t)} generated table(s), {len(t_unused)} unused; "
          f"{len(gen_f)} generated figure file(s) covering "
          f"{len({n for n in gen_f})} name(s), {len(f_unused)} unused")
    for n in t_unused:
        print(f"    results/tables/{n}.tex is GENERATED BUT NEVER \\input: "
              f"use it, or stop generating it")
    for n in t_missing:
        print(f"    paper/ inputs results/tables/{n}, which no generator writes")
    for n in f_unused:
        print(f"    results/figures/{n} is GENERATED BUT NEVER INCLUDED: "
              f"use it, or stop generating it")
    for n in f_missing:
        print(f"    paper/ includes {n}, which no generator writes")
    return not (t_unused or t_missing or f_unused or f_missing)


def check_figure_numbers():
    """
    Every number a figure DRAWS must be classified, and a claim must be
    registered.

    THE BUG THIS EXISTS FOR. The phase overlay drew "25 phases" as a panel
    header on three of five schedules whose lattices are 35, 28 and 28. The
    header was min(q, N) and N was the visit count, so three panels printed
    their own sample size back at the reader as though it were a property of
    the schedule. Nothing caught it: the number was computed inside the
    figure module at draw time, so it was in no table, in no numbers file, and
    in nothing the artefact guard compares. A table cell that wrong fails
    check_regenerated() the same day.

    THE RULE. A number a figure draws as a claim about the design must come
    from a registered generator, exactly as a table number does. A number that
    is an axis tick, a geometry fact about the drawing, or a restatement of
    the setting the panel was drawn at, need not.

    WHAT IS ENFORCEABLE AND WHAT IS NOT. The machine can find every number
    drawn at a call site and it can verify that an expression reads the
    register. It CANNOT tell a claim from an axis annotation: `{v:.2f}` in a
    title is the same syntax either way, and no amount of static analysis
    decides which one it is. So the rule is enforced as a DECLARATION. Each
    figure module carries a DRAWN dict mapping the source text of every
    interpolation it draws to one of DRAWN_KINDS, and:

      * a number drawn and not declared FAILS -- so nothing new slips in;
      * a declaration for a number no longer drawn FAILS -- so the list
        cannot rot into a blanket exemption;
      * a `registered:` declaration is checked, by _registers(), against what
        the expression actually does;
      * `UNREGISTERED:` is a claim with no generator yet. It is counted and
        printed on every run, so the debt is a number that can be driven to
        zero rather than a silence.

    What is NOT enforceable: whether a declaration is honest. Writing
    `geometry:` over a claim passes. That is a review question, and the
    declaration is what gives a reviewer something to review.

    Also outside the check: text assembled away from the drawing call, as
    g5_conditions.py builds its legend strings into a list first. Those are
    reported by count at the end and are not classified.
    """
    helpers = _text_helpers()
    problems, debt, n_drawn, n_ticks, n_reg = [], [], 0, 0, 0
    elsewhere = []
    for py in sorted(list((ROOT / "audit").glob("*.py"))
                     + list((ROOT / "src").rglob("*.py"))):
        text = py.read_text()
        if "save(" not in text and "savefig(" not in text:
            continue
        rel = str(py.relative_to(ROOT))
        drawn, ticks = _drawn(text, helpers)
        n_drawn += len(drawn)
        n_ticks += len(ticks)
        n_away = _built_elsewhere(text, helpers)
        if n_away:
            elsewhere.append(f"{rel}: {n_away}")
        if not drawn:
            continue
        ns = {}
        for node in ast.walk(ast.parse(text)):
            if (isinstance(node, ast.Assign) and len(node.targets) == 1
                    and isinstance(node.targets[0], ast.Name)
                    and node.targets[0].id == "DRAWN"):
                ns = ast.literal_eval(node.value)
        if not ns:
            problems.append(f"{rel} draws {len(drawn)} interpolated "
                            f"number(s) and declares no DRAWN dict")
            continue
        for expr, line in sorted(drawn.items(), key=lambda kv: kv[1]):
            kind = ns.get(expr)
            if kind is None:
                problems.append(f"{rel}:{line} draws {expr!r} and DRAWN does "
                                f"not classify it")
            elif not kind.startswith(DRAWN_KINDS):
                problems.append(f"{rel}:{line} classifies {expr!r} as "
                                f"{kind!r}, which is not one of "
                                f"{DRAWN_KINDS}")
            elif kind.startswith("registered:"):
                n_reg += 1
                if not _registers(text, expr):
                    problems.append(
                        f"{rel}:{line} declares {expr!r} registered, but it "
                        f"neither calls num() nor indexes the dict this "
                        f"module writes")
            elif kind.startswith("UNREGISTERED:"):
                debt.append(f"{rel}:{line} {expr} -- {kind.split(':', 1)[1].strip()}")
        for expr in sorted(set(ns) - set(drawn)):
            problems.append(f"{rel} DRAWN classifies {expr!r} and nothing "
                            f"draws it")
    #  the same standard the reachability check holds itself to: a check that
    #  cannot be shown to fail is not evidence. These are the three shapes
    #  that matter -- an undeclared number, a false `registered:`, and a
    #  helper-drawn number the matplotlib name list would miss.
    probe = ('import x\ndef f(ax, v, N):\n'
             '    ax.set_title(f"factor {v:.2f}")\n'
             '    at_point(ax, (0, 0), f"q = {N[\'k\']}")\n'
             '    save(fig, "p")\n')
    pd, pt = _drawn(probe, {"at_point"})
    probe_ok = [set(pd) == {"v", "N['k']"},
                not _registers(probe, "v"),
                not _registers(probe, "N['k']"),
                _registers(probe + 'write_numbers("m", N)\n', "N['k']"),
                _registers(probe + 'write_numbers("m", {"a": v})\n', "v"),
                _registers(probe + 'B = {}\nwrite_numbers("m", B)\n'
                           'B.update({"a": v})\n', "v"),
                _registers(probe, 'num("design", "sigma")')]
    if not all(probe_ok):
        problems.append(f"the figure-number check's own self-test failed: "
                        f"{probe_ok}")
    if elsewhere:
        print(f"  {sum(int(e.split(': ')[1]) for e in elsewhere)} "
              f"interpolated string(s) in figure modules are built away from "
              f"any drawing call and are NOT classified: "
              f"{'; '.join(elsewhere)}")
    print(f"  {n_drawn} number(s) drawn on figures checked for provenance "
          f"({n_reg} registered, {len(debt)} unregistered claim(s), "
          f"{n_ticks} tick-label interpolation(s) out of scope) "
          f"(self-test: {'PASS' if all(probe_ok) else 'FAIL'})")
    for d in debt:
        print(f"    debt: {d}")
    for why in problems:
        print(f"    {why}")
    return not problems

_INPUT_RE = re.compile(r"\\(?:input|include)\{([^}]*)\}")


def _uncommented(text):
    """`text` with LaTeX comments removed, respecting \\%."""
    out = []
    for line in text.splitlines():
        i, cut = 0, len(line)
        while i < len(line):
            if line[i] == "\\":
                i += 2
                continue
            if line[i] == "%":
                cut = i
                break
            i += 1
        out.append(line[:cut])
    return "\n".join(out)


def _tex_name(target):
    """An \\input target as a key: path from the compile directory, with .tex."""
    return target if target.endswith(".tex") else target + ".tex"


def _input_closure(sources, roots):
    """
    Every file the documents actually read, and every target that resolves to
    nothing.

    `sources` maps a name to its text, keyed the way LaTeX resolves an
    \\input: against the directory pdflatex runs in, which is paper/, and NOT
    against the including file. sections/A1_supplementary.tex says
    `\\input{sections/01_biology}` and means paper/sections/01_biology.tex.
    """
    seen, stack, missing = set(), list(roots), []
    while stack:
        name = stack.pop()
        if name in seen:
            continue
        if name not in sources:
            missing.append(name)
            continue
        seen.add(name)
        for m in _INPUT_RE.finditer(_uncommented(sources[name])):
            stack.append(_tex_name(m.group(1)))
    return seen, missing


def check_inputs_reachable():
    r"""
    Every .tex in paper/ is read by a document, or is a document itself.

    THE BUG THIS EXISTS FOR, and it is the one check in this file whose
    absence is recorded in the log before the check. The em-dash sweep
    swallowed the line `\input{sections/01_biology}` and took a 1,400-word
    section out of the build. Nothing failed: no LaTeX error, no undefined
    reference, no failed test, because a section that is not input defines no
    label that anything misses. The page count fell by four and the fall was
    put down to a shorter abstract (REVIEW_LOG.md R10, R50).

    It is check_generators_reachable() for prose. That one asks whether every
    numbers file has a generator `make figures` actually runs; this asks
    whether every section file has a document that actually reads it. A file
    present and never input is as invisible as a number never generated.

    The exemption is standalone documents, and it is CHECKED rather than
    declared: a file nothing inputs passes only if it carries its own
    \documentclass, which is what note_structural_identifiability.tex is.

    THE LIMIT, stated rather than assumed. This follows \input and \include
    and nothing else, so it sees a file that left the build and not a passage
    that did. Prose commented out inside a file the documents DO read, or put
    behind a conditional that is never true, is invisible to it exactly as
    01_biology was before. R10's other two defects were of that kind: two
    sentences pushed behind a comment marker, which no check here catches and
    which is why R10's rule is that a sweep is read or it is not done.
    """
    paper = ROOT / "paper"
    sources = {str(p.relative_to(paper)): p.read_text()
               for p in sorted(paper.rglob("*.tex"))}
    roots = ["main.tex", "supplement.tex"]
    seen, missing = _input_closure(sources, roots)
    #  targets outside paper/ are the generated numbers and tables, which
    #  check_regenerated already owns; they are required to exist, not to be
    #  reachable from anywhere else.
    outside = [n for n in missing if n.startswith("..")]
    missing = [n for n in missing if not n.startswith("..")]
    for n in outside:
        if not (paper / n).exists():
            missing.append(n)
    orphans = []
    for name in sorted(set(sources) - seen):
        if r"\documentclass" not in sources[name]:
            orphans.append(name)
    #  the same standard as the other checks: the detector is run against a
    #  tree with an orphan and a tree without one, and against a dangling
    #  \input, so a clean sweep is a result and not an absence.
    probe = {"main.tex": r"\input{a}" + "\n% \\input{c}\n", "a.tex": "", "b.tex": ""}
    p_seen, p_missing = _input_closure(probe, ["main.tex"])
    clean_seen, _ = _input_closure({"main.tex": r"\input{a}", "a.tex": ""},
                                   ["main.tex"])
    _, dangling = _input_closure({"main.tex": r"\input{gone}"}, ["main.tex"])
    probe_ok = (set(probe) - p_seen == {"b.tex"} and not p_missing
                and not set({"main.tex", "a.tex"}) - clean_seen
                and dangling == ["gone.tex"])
    standalone = len(set(sources) - seen) - len(orphans)
    print(f"  {len(sources)} .tex file(s) in paper/ checked for reachability "
          f"from main.tex and supplement.tex: {len(seen)} read, "
          f"{standalone} standalone document(s), {len(orphans)} orphan(s), "
          f"{len(missing)} dangling input(s) (self-test: "
          f"{'PASS' if probe_ok else 'FAIL'})")
    for name in orphans:
        print(f"    paper/{name} is never input and is not a document")
    for name in missing:
        print(f"    \\input{{{name[:-4]}}} resolves to no file")
    return not orphans and not missing and probe_ok


def _compile_warnings(sources):
    """{path: [(line, category, message)]} from {path: source text}."""
    import warnings
    out = {}
    for path in sorted(sources):
        with warnings.catch_warnings(record=True) as got:
            warnings.simplefilter("always")
            try:
                compile(sources[path], path, "exec")
            except SyntaxError as e:
                out.setdefault(path, []).append(
                    (e.lineno, "SyntaxError", str(e.msg)))
                continue
        for w in got:
            out.setdefault(path, []).append(
                (w.lineno, w.category.__name__, str(w.message)))
    return out


def _tracked_python():
    """Every .py file under version control, as {path: source text}."""
    try:
        out = subprocess.run(["git", "ls-files", "*.py"], cwd=ROOT,
                             capture_output=True, text=True, check=True).stdout
        names = out.split()
    except (OSError, subprocess.CalledProcessError):
        names = [str(p.relative_to(ROOT)) for p in ROOT.rglob("*.py")
                 if ".git" not in p.parts]
    return {n: (ROOT / n).read_text() for n in sorted(names)}


def check_source_warnings():
    r"""
    Every tracked module must compile from source without a warning.

    THE BUG THIS EXISTS FOR, and its shape is not the shape of the others.
    `pytest.ini` sets `filterwarnings = error`, and that was recorded as
    having cost nothing because the suite passed under it. The suite passed
    on a warm `__pycache__`. A warning about a string literal is raised when
    a module is COMPILED; once a `.pyc` is on disk the module is not
    recompiled and the warning is never raised again. Fifteen invalid escape
    sequences sat in four files while the same tree returned `92 passed` or
    `38 failed, 54 passed` depending on nothing but whether `make clean` had
    run first, and neither run said why (REVIEW_LOG.md R49, R50).

    THE RULE IT ENFORCES is not about escapes. It is that a check's verdict
    must not depend on build state. This one compiles the source text itself
    on every run and never consults the cache, so it cannot be warm. It
    covers `tools/` and `tests/` as well as the generators, which `-W error`
    on `make figures` does not reach.
    """
    sources = _tracked_python()
    found = _compile_warnings(sources)
    #  the same standard the other checks hold themselves to. The detector is
    #  run against a source that DOES warn and one that does not, so a clean
    #  sweep is a result rather than an absence.
    bad = 'x = "\\d"\n'
    good = 'x = r"\\d"\n'
    probe_ok = (len(_compile_warnings({"p": bad}).get("p", [])) == 1
                and not _compile_warnings({"g": good}))
    print(f"  {len(sources)} tracked python file(s) compiled from source, "
          f"{sum(len(v) for v in found.values())} warning(s) "
          f"(self-test: {'PASS' if probe_ok else 'FAIL'})")
    for path in sorted(found):
        for line, cat, msg in found[path]:
            print(f"    {path}:{line}: {cat}: {msg}")
    return not found and probe_ok


# ==========================================================================
#  SPELLED-OUT NUMERALS
# ==========================================================================
#: A paper number written as an English word is invisible to every guard in
#: this file, because all of them compare digits. That is not hypothetical: the
#: one error this repository's provenance sweep was provoked by was the clause
#: "from thirteen visits to eighteen" (REVIEW_LOG.md R75, withdrawn by R82),
#: and a thousand of the paper's 2472 unguarded literals are spelled out.
#:
#: The alternative was to reword those sentences into digits so the existing
#: machinery could see them. The author declined, and rightly: the paper's
#: register is not something to spend on a tool's convenience. So the tool
#: learns to read English instead.
#:
#: Each row declares that a WORD in a section file states a GENERATED value,
#: under one of three relations, because English numerals in prose are not
#: always exact:
#:
#:    exact      the word spells the value        "thirteen" <- 13
#:    rounds     the word spells round(value)     "ten" <- 9.8, "nearly ten years"
#:    below      the value is under the word      "one" <- 0.33, "under one per cent"
#:
#: A row whose word is no longer in the file FAILS, so the table cannot rot
#: into a blanket exemption, and a row whose relation stops holding FAILS,
#: which is the whole point.
SPELLED = (
    ("11_scope.tex", "thirteen", "meas_match_nominal", "exact"),
    ("11_scope.tex", "eighteen", "meas_match_2pct", "exact"),
    ("11_scope.tex", "six", "meas_gate_nominal", "exact"),
    ("08_specification.tex", "eleven", "tier1_visits_exact", "exact"),
    ("08_specification.tex", "ten", "tier1_years_at_pub", "rounds"),
    ("06_excitation.tex", "one", "meas_b_only_gain_pct", "below"),
)

_ONES = ("zero one two three four five six seven eight nine ten eleven twelve "
         "thirteen fourteen fifteen sixteen seventeen eighteen nineteen").split()
_TENS = {20: "twenty", 30: "thirty", 40: "forty", 50: "fifty", 60: "sixty",
         70: "seventy", 80: "eighty", 90: "ninety"}


def _spell(n):
    """An integer in English, over the range a design quantity occupies."""
    n = int(n)
    if n < 0 or n > 999:
        return None
    if n < 20:
        return _ONES[n]
    if n < 100:
        t, r = divmod(n, 10)
        return _TENS[t * 10] + ("-" + _ONES[r] if r else "")
    h, r = divmod(n, 100)
    return _ONES[h] + " hundred" + (" " + _spell(r) if r else "")


def check_spelled_numbers(numdir=None):
    """
    Every declared English numeral against the macro it states.

    This is the guard's one check that does not compare digits, and it exists
    because the paper is allowed to write "thirteen".
    """
    numdir = (ROOT / "results" / "numbers") if numdir is None else Path(numdir)
    vals = {}
    for f in sorted(Path(numdir).glob("*.tex")):
        for m in re.finditer(
                r"csname paper@num@([A-Za-z0-9_]+)\\endcsname\{([^}]*)\}",
                f.read_text()):
            vals[m.group(1)] = m.group(2)
    bad = []
    for fn, word, key, rel in SPELLED:
        src = ROOT / "paper" / "sections" / fn
        text = _uncommented(src.read_text()) if src.exists() else ""
        if not re.search(r"(?<![-\w])" + re.escape(word) + r"(?![-\w])",
                         text, re.I):
            bad.append(f"{fn}: the word '{word}' is no longer there, so the "
                       f"row for {key} is stale")
            continue
        if key not in vals:
            bad.append(f"{fn}: '{word}' declares {key}, which no generator "
                       f"writes")
            continue
        try:
            v = float(vals[key])
        except ValueError:
            bad.append(f"{key} is not numeric: {vals[key]!r}")
            continue
        if rel == "exact":
            want = _spell(round(v))
            if want != word.lower():
                bad.append(f"{fn}: '{word}' states {key}, which is {vals[key]}, "
                           f"spelled '{want}'")
        elif rel == "rounds":
            want = _spell(round(v))
            if want != word.lower():
                bad.append(f"{fn}: '{word}' rounds {key} = {vals[key]}, which "
                           f"rounds to '{want}'")
        elif rel == "below":
            bound = _ONES.index(word.lower()) if word.lower() in _ONES else None
            if bound is None or not v < bound:
                bad.append(f"{fn}: '{word}' claims {key} = {vals[key]} is "
                           f"below {bound}, and it is not")
    for b in bad:
        print("  SPELLED: " + b)
    print(f"  {len(SPELLED)} spelled-out numeral(s) checked against the macros "
          f"they state (self-test: {'PASS' if _spell(13) == 'thirteen' and _spell(18) == 'eighteen' and _spell(25) == 'twenty-five' else 'FAIL'})")
    return not bad


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--working-tree", action="store_true",
                    help="compare results/ against HEAD instead of "
                         "regenerating the tables")
    args = ap.parse_args()
    ok = check_working_tree() if args.working_tree else check_regenerated()
    #  run in BOTH modes: the reachability question is about the build, and
    #  the fast path is precisely the place where a build-unreachable
    #  generator looks fine.
    ok = check_generators_reachable() and ok
    #  and the same question for a number a figure DRAWS, which until this
    #  existed was the one number in the repository with no provenance at all
    ok = check_figure_numbers() and ok
    #  and the one mechanical case of the referent failure: two files
    #  defining one key, which \input order would settle in silence
    ok = check_key_collisions() and ok
    #  and its quieter opposite: two keys carrying ONE value. The guard
    #  compares keys, so a redundant key passes it, and a duplicate that
    #  agrees today is the only kind that can silently diverge tomorrow
    #  (REVIEW_LOG.md R146, R147)
    ok = check_value_twins() and ok
    #  and the level none of the above reached: a key can be generated,
    #  reachable, collision-free and STILL undefined in the document, because
    #  nothing here asked whether the preamble inputs the file it lives in.
    #  Six files were in that state and the paper raised 29 LaTeX errors while
    #  this script printed ARTEFACTS MATCH THE CODE (REVIEW_LOG.md R152)
    ok = check_numbers_loaded() and ok
    #  and the same question one file type over, which is where R152's defect
    #  would have been found earlier: a generated table nobody inputs and a
    #  generated figure nobody includes are the same failure as a numbers file
    #  nobody loads (the author's instruction after R152)
    ok = check_generated_artefacts_used() and ok
    #  and the class the cache hid: a check whose verdict depended on
    #  whether `make clean` had run. This one compiles the source itself.
    ok = check_source_warnings() and ok
    #  and the prose half of the reachability question: a section file that
    #  no document inputs is as invisible as a number no generator writes
    ok = check_inputs_reachable() and ok
    #  and the half no digit-based check can see: a number the paper spells
    #  out in English. The one error this repository's sweep was provoked by
    #  was exactly that (REVIEW_LOG.md R92)
    ok = check_spelled_numbers() and ok
    print("  ARTEFACTS MATCH THE CODE" if ok else "  ARTEFACTS HAVE DRIFTED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
