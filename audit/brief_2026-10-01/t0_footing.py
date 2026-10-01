"""
The repository facts REPORT.md quotes, printed from their sources so that every
number in the report can be traced to an out/*.txt file.

Nothing here is computed beyond reading constants, one information-matrix
agreement check, and the marker range a record visits.
"""
import re

import numpy as np

import core as C
from design.run import D_AMP, DOSE_INTERVAL, LAMBDA, LOAD_CV, LOAD_FRAC, ASSAY_CV


def main():
    print("=" * 92)
    print("  REPOSITORY FOOTING QUOTED BY REPORT.md")
    print("=" * 92)
    print(f"  design.run: D_AMP = {D_AMP:g}, DOSE_INTERVAL = {DOSE_INTERVAL:g} d, increment half-life "
          f"{np.log(2)/LAMBDA:.2f} d")
    print(f"  design.run: LOAD_CV = {LOAD_CV:.2f} (day-to-day intake variation, recorded), "
          f"LOAD_FRAC = {LOAD_FRAC:.2f}, ASSAY_CV = {ASSAY_CV:.2f}")
    print(f"  core: monthly record of {C.MONTHLY.size} draws; SEEDS = {list(C.SEEDS)[0]}..{list(C.SEEDS)[-1]}")

    #  the brief's ln-m information matrix against the flux-space one
    worst = 0.0
    for th in (150.0, 3000.0, 20000.0):
        for s in (0, 3):
            r = C.record(th, seed=s)
            F1, F2 = C.fim(r, noise="perdraw"), C.fim_lnm(r)
            worst = max(worst, float(np.abs(F1 - F2).max() / np.abs(F1).max()))
    print(f"  ln-m FIM by finite differences of the steady-state solve, against the flux-space FIM "
          f"(per-draw noise): worst relative {worst:.1e}")

    #  the enzyme's saturation g = m/(K_m + m), untreated and under therapy
    for th, lab in ((C.TH0, "normal"), (C.TH_5, "5%")):
        L, m0 = C.N4.anchor(th, C.KM0, C.KE, "A")
        gs = []
        for s in C.SEEDS:
            r = C.record(th, seed=s)
            gs.append(r["m"] / (C.KM0 + r["m"]))
        gs = np.concatenate(gs)
        print(f"  saturation g at {lab:>6}: untreated {m0/(C.KM0+m0):.3f} (m = {m0:.2f}); on therapy, "
              f"over the 24 monthly records: {gs.min():.3f} to {gs.max():.3f} (median {np.median(gs):.3f})")

    #  the parent model's fold, as the repository's own generated numbers state it
    txt = (C.ROOT / "results" / "numbers" / "boundary.tex").read_text()
    m = re.search(r"fold_load_ratio\\endcsname\{([0-9.]+)\}", txt)
    print(f"  results/numbers/boundary.tex: fold_load_ratio = {m.group(1) if m else 'NOT FOUND'}")

    #  the breath-test repeatability line the report cites
    lines = (C.ROOT / "audit" / "PHASE4_SUMMARY.md").read_text().splitlines()
    print(f"  audit/PHASE4_SUMMARY.md:208: {lines[207].strip()}")
    print("  (n6_profile_vs_crb.py's own output, run unchanged in this session, is out/n6_profile_vs_crb.txt)")


if __name__ == "__main__":
    main()
