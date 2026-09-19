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
