"""
The generated tables in results/ against the code that writes them.

The paper prints tables that are generated and committed. Nothing checked
that the committed file still matched its generator, and one of them stopped
matching: the design table was written by a function no entry point called,
so it froze before the seed work and printed single-seed numbers while the
text beside it quoted medians over 24 records (REVIEW_LOG.md R23). Table 3
then disagreed with Table 5 about the same design, in the same document.

This is that check. It regenerates every table into a scratch directory and
compares, so it writes nothing into results/ and a failure leaves the
evidence intact. Figures are not regenerated here because the decision
characteristic is simulated and costs minutes; CI runs `make figures` and
then `tools/check_artefacts.py --working-tree`, which compares the whole of
results/ against HEAD with the PDF date fields stripped.

Run:  pytest tests/test_artefacts.py -q
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tools"))

import numpy as np                                            # noqa: E402

from common.models import ParamsA, sensitivity_A              # noqa: E402


def test_generated_tables_match_the_code():
    """Every committed table is what its generator produces today."""
    import check_artefacts
    assert check_artefacts.check_regenerated(), (
        "a committed table no longer matches its generator; run "
        "`make figures` and read the diff before committing anything else"
    )


def test_check_artefacts_fails_when_a_table_drifts(tmp_path):
    """The guard itself: a perturbed table must fail the check.

    A check that cannot fail is worse than no check, because it reads as
    evidence. This perturbs one committed number, asserts the check rejects
    it, and puts the file back.
    """
    import shutil
    import check_artefacts
    live = ROOT / "results" / "tables"
    copy = tmp_path / "tables"
    shutil.copytree(live, copy)
    target = copy / "tab_D_designs.tex"
    original = target.read_text()
    target.write_text(original.replace("0.1699", "0.1698", 1))
    assert target.read_text() != original, "the perturbation did not apply"
    assert not check_artefacts.check_regenerated(live=copy)
    #  and nothing in the repository was touched to find that out
    assert (live / "tab_D_designs.tex").read_text() == original


def test_three_distinct_levels_identify_the_local_model():
    """What the orphaned rank table used to assert, as an assertion.

    `table_rank()` wrote results/tables/tab_A_rank.tex on every build and no
    document ever showed it. The claim it carried is Proposition 1's: the
    sensitivity matrix reaches full rank at three distinct marker levels and
    not before, so a fourth adds nothing to the rank.
    """
    p = ParamsA()
    for levels, want in (([10.0], 1),
                         ([10.0, 12.0], 2),
                         ([10.0, 12.0, 14.0], 3),
                         ([10.0, 50.0, 200.0], 3),
                         ([10.0, 12.0, 14.0, 16.0], 3)):
        S = np.array([sensitivity_A(m, p) for m in levels])
        assert np.linalg.matrix_rank(S, tol=1e-10) == want, levels


#  ------------------------------------------------------------------ the
#  exemption for simulated cells. The decision table's cells come out of a
#  nonlinear least-squares fit to 400 synthetic records, twelve times over.
#  That does not reproduce bit for bit on another machine, and CI proved it:
#  two runners on the same commit, one printing 2.9 where the other printed
#  3.0. The guard tolerates that and nothing else, and these are the tests
#  that say where the line is.

def _decision_table():
    return (ROOT / "results" / "tables" / "tab_D_decision.tex").read_text()


def test_the_measured_machine_noise_is_tolerated():
    """The exact difference CI saw between two runners must not fail."""
    import check_artefacts
    old = _decision_table()
    new = old.replace(r"$42.3\% \pm 2.9$", r"$42.3\% \pm 3.0$", 1)
    assert new != old, "the table no longer contains the cell this tests"
    assert check_artefacts.same_table("tab_D_decision.tex", old, new)


def test_the_tolerance_does_not_apply_to_other_tables():
    """Only the declared table is exempt; every other one stays byte-exact."""
    import check_artefacts
    old = _decision_table()
    new = old.replace(r"$42.3\% \pm 2.9$", r"$42.3\% \pm 3.0$", 1)
    assert not check_artefacts.same_table("tab_D_designs.tex", old, new)


def test_drift_past_the_tolerance_still_fails():
    """A number that moved by more than the simulation's own noise."""
    import check_artefacts
    old = _decision_table()
    for was, now in ((r"$19.1\% \pm 2.1$", r"$22.0\% \pm 2.1$"),
                     (r"$4900$", r"$7000$")):
        new = old.replace(was, now, 1)
        assert new != old, f"the table no longer contains {was}"
        assert not check_artefacts.same_table("tab_D_decision.tex", old, new)


def test_structure_is_still_compared_exactly():
    """Everything that is not a value in a value column fails on any change."""
    import check_artefacts
    old = _decision_table()
    cases = {
        "row label": (r"monthly across dose cycles", r"weekly across dose cycles"),
        "words in a cell": (r"n/a", r"n/b"),
        "the note": (r"twelve seeds", r"eleven seeds"),
        "the header": (r"\textbf{cohort $N$}", r"\textbf{cohort $M$}"),
    }
    for what, (was, now) in cases.items():
        new = old.replace(was, now, 1)
        assert new != old, f"the table no longer contains {was} ({what})"
        assert not check_artefacts.same_table("tab_D_decision.tex", old, new), what
    #  a dropped row is a changed line count
    lines = old.splitlines()
    dropped = "\n".join(lines[:4] + lines[5:]) + "\n"
    assert not check_artefacts.same_table("tab_D_decision.tex", old, dropped)


def test_the_working_tree_mode_uses_the_same_rule(tmp_path):
    """CI runs --working-tree, which compares a file on disk with git's bytes.

    The unit tests above exercise the comparison; this one exercises the path
    that reaches it from a working tree, because that is the one CI takes.
    """
    import check_artefacts
    head = _decision_table()
    on_disk = tmp_path / "tab_D_decision.tex"
    on_disk.write_text(head.replace(r"$42.3\% \pm 2.9$", r"$42.3\% \pm 3.0$", 1))
    assert check_artefacts.tolerated(on_disk, head.encode())
    on_disk.write_text(head.replace(r"$42.3\% \pm 2.9$", r"$42.3\% \pm 9.9$", 1))
    assert not check_artefacts.tolerated(on_disk, head.encode())
    #  a table that is not declared simulated gets no tolerance at all
    other = tmp_path / "tab_D_designs.tex"
    other.write_text(head.replace(r"$42.3\% \pm 2.9$", r"$42.3\% \pm 3.0$", 1))
    assert not check_artefacts.tolerated(other, head.encode())
