"""
Numbers the paper quotes, written by the code that computes them.

A number typed into the .tex is a number nothing checks. The design table
drifted for a month because no entry point regenerated it (REVIEW_LOG.md
R23), and a literal in a sentence is the same hazard with no file to diff.

So every quantity a module computes and the paper quotes is written here as
a macro, and the text says \num{sd_monthly} rather than 0.170. The file is
generated, so tools/check_artefacts.py compares it like any other artefact
and a number that moves fails a check instead of sitting in a sentence.

Keys keep the names the code uses, underscores and all, because a key that
has to be transliterated into a LaTeX control sequence is a key that will be
transliterated wrongly. The preamble defines

    \\newcommand{\\num}[1]{\\csname paper@num@#1\\endcsname}

and \\csname accepts any character tokens, so nothing has to be renamed.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NUMDIR = ROOT / "results" / "numbers"


#: which file owns each generated numbers file.
#:
#: The header used to be built from the output stem, so every file in
#: results/numbers/ named a generator that does not exist: design.tex said
#: `src/design.py` for `src/design/run.py`, scope.tex said `src/scope.py` for
#: `audit/n4_theta_scope.py`, and so on for all seven. In a repository whose
#: rule is that no number reaches the paper without a row saying where it came
#: from, a provenance line pointing at a missing file is worse than none.
#:
#: It is a table and not the caller's own `__file__` because the artefact
#: guard regenerates every one of these into a temporary directory and
#: compares the result byte for byte. Read off the caller's frame, the header
#: would say `tools/check_artefacts.py` on the guard's copy and the generator
#: on the committed one, and the guard would fail on its own provenance line.
#: A stem with no entry here is a stem nobody has claimed, and it says so.
GENERATORS = {
    "design": "src/design/run.py",
    "local": "src/model_A/run.py",
    "scope": "audit/n4_theta_scope.py",
    "kmcost": "audit/n1_km_cost.py",
    "decision": "audit/d13_decision_therapy.py",
    "profile": "audit/n6_profile_vs_crb.py",
    "transient": "audit/t1_transient_worth.py",
    "reduction": "src/model_C/run.py",
    "kesweep": "audit/t5_band_vs_ke.py",
    "cohort": "audit/t6_cohort_with_kappa.py",
    "structure": "audit/t7_structure_transient.py",
    "lattice": "audit/t11_phase_lattice.py",
    "phaseclock": "audit/g1_phase_clock.py",
    "fibre": "src/model_A/run.py",
    "margin": "src/model_C/run.py",
    "therapy": "audit/t13_therapy_forward.py",
}


class DuplicateKey(Exception):
    """A key given two different values."""


class NumberBook(dict):
    """
    A dict of generated numbers that refuses to change a value silently.

    THE BUG THIS EXISTS FOR. `track_pct` was the marker's departure AFTER a
    tolerance-dependent settling window, keyed as though it were
    tolerance-independent. Three quasi_steady_check() calls in one build each
    ran .update() on it and whichever went last decided the value; the
    artefact guard and `make figures` called that function a different number
    of times in a different order and wrote 0.9861 and 4.9745 for the same
    name. Nothing flagged it, because last-write-wins is silent.

    A key with two values is worse than a stale one. A stale number is wrong
    everywhere and can be caught by regenerating; a duplicated one is right on
    the path that happened to run last and wrong on every other, so the two
    paths have to disagree before anyone notices.

    The check is at ACCUMULATION and not at the write, because that is where
    the collision happens: by the time write_numbers() sees the dict it holds
    one value and the other is gone. Rewriting a key with the SAME value is
    allowed, since a generator legitimately called twice at the same settings
    must not fail.
    """

    def __setitem__(self, key, value):
        old = self.get(key, _MISSING)
        if old is not _MISSING and not _same(old, value):
            raise DuplicateKey(_explain({key: (old, value)}))
        dict.__setitem__(self, key, value)

    def update(self, other=(), **kw):
        clashes = {}
        items = dict(other, **kw) if other else dict(kw)
        for k, v in items.items():
            old = self.get(k, _MISSING)
            if old is not _MISSING and not _same(old, v):
                clashes[k] = (old, v)
        if clashes:
            raise DuplicateKey(_explain(clashes))
        dict.update(self, items)

    def setdefault(self, key, default=None):
        if key not in self:
            self[key] = default
        return self[key]


def _explain(clashes):
    lines = "\n".join(f"    {k}: {a!r} then {b!r}"
                       for k, (a, b) in sorted(clashes.items()))
    return (f"{len(clashes)} generated key(s) given two different values:\n"
            f"{lines}\n"
            "  A generated key must name one quantity. If two callers compute\n"
            "  different things, they need different keys; see\n"
            "  quasi_steady_check's tolerance suffix for the pattern.")


#: every (module, key) this process has written, and what it wrote.
#:
#: `track_pct` was the departure AFTER a tolerance-dependent settling window
#: and was keyed as though it were tolerance-independent, so the three
#: quasi_steady_check() calls in one build each overwrote it and whichever ran
#: last decided the value. The guard and `make figures` called it a different
#: number of times in a different order and wrote 0.9861 and 4.9745 for the
#: same name. Nothing flagged it, because last-write-wins is silent.
#:
#: A key with two values is worse than a stale one: a stale number is wrong
#: everywhere and a duplicated one is wrong only on some paths. So writing an
#: existing key with a DIFFERENT value is an error, and writing it with the
#: same value is fine -- a generator legitimately called twice at the same
#: settings must not fail.
_WRITTEN = {}


def _check_duplicates(module, values):
    """
    The same rule across write_numbers() calls, for a module written twice.

    NumberBook catches the collision inside one accumulating dict. This
    catches a module whose generator is run twice in one process and produces
    a different number the second time, which a fresh dict each time would
    otherwise hide.
    """
    clashes = {}
    for k, v in values.items():
        seen = _WRITTEN.get((module, k), _MISSING)
        if seen is not _MISSING and not _same(seen, v):
            clashes[k] = (seen, v)
        _WRITTEN[(module, k)] = v
    if clashes:
        raise DuplicateKey(f"into results/numbers/{module}.tex: "
                           + _explain(clashes))


class _Missing:
    pass


_MISSING = _Missing()


def _same(a, b):
    if isinstance(a, float) and isinstance(b, float):
        if a != a and b != b:
            return True
        return a == b or abs(a - b) <= 1e-12 * max(abs(a), abs(b), 1.0)
    return a == b


def reset_duplicate_check():
    """Forget what this process has written. For the guard's own self-test."""
    _WRITTEN.clear()


def write_numbers(module, values, formats=None, digits=4):
    """
    Write `values` as \\num-able macros into results/numbers/<module>.tex.

    A value may be a float, an int or a preformatted string. Non-finite
    values are written as \\infty, because a design that does not identify
    anything is a result and not a missing number.

    `formats` gives a format string per key, because the paper quotes the
    same quantity at the precision its sentence needs: the aliasing cost is
    "a factor of 27" and the monthly precision is "0.170", and a file that
    emitted 26.9776 and 0.1699 for both would force the text to round by
    hand, which is the habit this file exists to remove.
    """
    formats = formats or {}
    _check_duplicates(module, values)
    NUMDIR.mkdir(parents=True, exist_ok=True)
    src = GENERATORS.get(module, f"an unregistered generator for {module!r}")
    lines = [f"%  generated by {src} -- do not edit",
             "%  every entry is quoted in the paper through \\num{<key>}"]
    for k in sorted(values):
        v = values[k]
        if isinstance(v, str):
            out = v
        elif v is None:
            continue
        elif v != v or v in (float("inf"), float("-inf")):
            out = r"\infty"
        elif isinstance(v, int):
            out = str(v)
        elif k in formats:
            out = formats[k].format(v)
        else:
            out = f"{v:.{digits}f}"
        lines.append(r"\expandafter\def\csname paper@num@%s\endcsname{%s}"
                     % (k, out))
    path = NUMDIR / f"{module}.tex"
    path.write_text("\n".join(lines) + "\n")
    rel = path.relative_to(ROOT) if path.is_relative_to(ROOT) else path
    print(f"  wrote {rel}  ({len(lines) - 2} numbers)")
    return path


#: a \expandafter\def line as write_numbers() emits it
_MACRO = re.compile(r"\\csname paper@num@([A-Za-z0-9_]+)\\endcsname\{(.*)\}$")


def num(module, key):
    """
    The registered value of a generated number, as the paper prints it.

    For FIGURE TEXT. A table cell and a \\num{} in a sentence both come from
    results/numbers/, and until this existed a number drawn onto a figure did
    not: it was recomputed at draw time from whatever the figure module had
    to hand, so the figure and the sentence beside it could disagree and
    nothing would say so. The phase overlay did exactly that -- it printed
    min(q, N) as "phases" on three of five panels while the lattice count in
    the text was 35, 28 and 28 (REVIEW_LOG.md, the phase-lattice entry).

    The STRING is returned, not the float, and that is the point: the figure
    then draws character for character what \\num{key} typesets, so the two
    cannot differ by a rounding either. A caller wanting arithmetic should
    compute it and register the result, not parse this back.

    Raises rather than falling back, because a silent fallback is how a
    figure goes back to drawing its own number.
    """
    path = NUMDIR / f"{module}.tex"
    if not path.exists():
        raise KeyError(
            f"results/numbers/{module}.tex does not exist; "
            f"{GENERATORS.get(module, 'its generator')} has not run. "
            f"`make figures` runs the generators in dependency order.")
    for line in path.read_text().splitlines():
        m = _MACRO.search(line)
        if m and m.group(1) == key:
            return m.group(2)
    raise KeyError(f"{module}.tex has no key {key!r}; it has "
                   f"{', '.join(sorted(registered(module)))}")


def registered(module):
    """Every key results/numbers/<module>.tex defines."""
    path = NUMDIR / f"{module}.tex"
    if not path.exists():
        return set()
    return {m.group(1) for m in
            (_MACRO.search(ln) for ln in path.read_text().splitlines()) if m}
