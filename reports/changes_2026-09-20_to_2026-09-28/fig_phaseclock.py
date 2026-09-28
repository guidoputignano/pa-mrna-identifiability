"""
Figure 1 of the change report: the paper's Figure 6 as it stood in each draft.

Both drafts include results/figures/fig_E1_phaseclock.pdf as Figure 6 (the
20 Sept source at paper/sections/07_design.tex:300, the 27 Sept one at :386).
The file was overwritten after 20 Sept, so the 20 Sept panel cannot be read off
the working tree. This script regenerates each version with its own generator,
audit/g1_phase_clock.py as committed at that revision (the last commits on or
before 20 and 27 Sept, UTC), in a throwaway git
worktree, and copies the result here. Nothing in the main working tree is
touched.

Run:  python reports/changes_2026-09-20_to_2026-09-28/fig_phaseclock.py
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "fig"


def git(*args, cwd=ROOT):
    return subprocess.run(["git", *args], cwd=cwd, check=True,
                          capture_output=True, text=True).stdout.strip()


OLD = git("rev-list", "-1", "--before=2026-09-20T23:59:59Z", "HEAD")
NEW = git("rev-list", "-1", "--before=2026-09-27T23:59:59Z", "HEAD")


def build(rev, name):
    with tempfile.TemporaryDirectory() as tmp:
        wt = Path(tmp) / "wt"
        git("worktree", "add", "--detach", "-q", str(wt), rev)
        try:
            subprocess.run([sys.executable, "audit/g1_phase_clock.py"], cwd=wt,
                           check=True, capture_output=True, text=True)
            shutil.copy(wt / "results/figures/fig_E1_phaseclock.pdf", OUT / f"{name}.pdf")
        finally:
            git("worktree", "remove", "--force", str(wt))
    print(f"{name}.pdf  from audit/g1_phase_clock.py @ {rev[:7]}")


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    build(OLD, "phaseclock_0920")
    build(NEW, "phaseclock_0927")
