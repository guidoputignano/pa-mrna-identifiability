"""
The generator for the constancy passage's numbers.

WHY THIS EXISTS BEFORE THE PASSAGE IS CITED, per R141 and the author's
instruction. `audit/u7_drift.py` computes every drift figure the proposed
constancy remark quotes and emits no numbers file, so the passage's macros would
have no generator. **A macro whose value is typed by hand is worse than a
literal, because the guard then vouches for it.** So the generator comes first
and the passage cites it, not the other way round.

EVERY VALUE IS COMPUTED BY CALLING THE FUNCTION THAT OWNS IT. Nothing here is
tabulated from `PHASE11_DRIFT.md`, which is a findings document rather than a
source; `u7_drift.sd_bias` is the owner and is called directly at the same
`n` and `r` the passage names.

THE FOOTING IS THE TIER 1 REQUIREMENT AS WRITTEN, `n = 118`, and that is the
whole point of the passage. `PHASE11_DRIFT.md` reports two crossings for two
different questions: the fixed 118-draw requirement, which fails at
`r = 0.005`, and a schedule free to re-optimise, which crosses near
`r = 0.009`. **Quoting one for the other would misstate the result**, so this
file emits only the first and names it in the key.

THE HALF-LIFE IS DERIVED, not asserted: `ln 2 / r`. It is the figure that makes
the sentence land, because a requirement failing on a drift whose halving time
is longer than a century is a requirement with no margin.

Run:  python audit/u30_wire_drift.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

from common.numbers import write_numbers

import u7_drift as U7

#: the Tier 1 count, and the drift rate the passage quotes
TIER1_N = 118
R = 0.005
#: THE FOOTING, AND OMITTING IT WAS CAUGHT BY THIS FILE'S OWN CHECK. Tier 1
#: needs 118 draws precisely BECAUSE the delivered increment is characterised
#: only to the published 14.02%, and PHASE11_DRIFT.md's Tier 1 table carries
#: that throughout. The first version of this generator called `med_rmse`
#: without it and the no-drift anchor came out 0.1543 instead of the 0.20 that
#: 118 is defined to land on. The number looked plausible and was on the wrong
#: design. With it, r = 0 reproduces 0.2000 exactly.
KAPPA_CV = 0.1402
#: what the no-drift anchor MUST return, since 118 is defined by it
GATE = 0.20
#: days per month and per year, for the span the bias scales with
MONTH_D = 30.4
YEAR_D = 365.25


def drift_group():
    """
    sd, bias and rmse at the Tier 1 count under a 0.5%/yr drift, plus the span
    the bias is set by.

    `med_rmse` returns the median over seeds with its range, so the figures
    carry a spread rather than being one draw, which is this repository's
    standing rule.
    """
    sd, bias, rmse, lo, hi = U7.med_rmse(TIER1_N, R, kappa_cv=KAPPA_CV)
    span_d = TIER1_N * MONTH_D
    span_y = span_d / YEAR_D
    out = {
        "drift_r_pct": 100.0 * R,
        "drift_sd": sd,
        "drift_bias": bias,
        "drift_rmse": rmse,
        "drift_rmse_lo": lo,
        "drift_rmse_hi": hi,
        "drift_halflife": float(np.log(2.0) / R),
        "drift_span_yr": span_y,
        #  the ceiling argument's own check: the bias tracks -r T / 2, so the
        #  predicted value is emitted beside the computed one and a reader can
        #  see they agree without running anything
        "drift_bias_span_pred": -R * span_y / 2.0,
    }
    #  and the no-drift anchor, which is what makes the count a boundary: 118
    #  is DEFINED as the draw count landing on the gate, so its zero-drift rmse
    #  must be the gate itself
    sd0, bias0, rmse0, _, _ = U7.med_rmse(TIER1_N, 0.0, kappa_cv=KAPPA_CV)
    if abs(rmse0 - GATE) > 5e-4:
        raise AssertionError(
            f"footing not matched: the no-drift rmse at n = {TIER1_N} is "
            f"{rmse0:.4f} against the {GATE} that count is DEFINED by. "
            f"Nothing is emitted, because a drift figure computed on a design "
            f"that is not Tier 1 would be quoted as though it were.")
    out["drift_rmse_nodrift"] = rmse0
    return out


if __name__ == "__main__":
    print("  computing the Tier 1 drift figures from u7_drift.sd_bias ...")
    N = drift_group()
    fmt = {k: "{:.4f}" for k in N}
    fmt.update({"drift_r_pct": "{:.1f}",
                "drift_halflife": "{:.1f}",
                "drift_span_yr": "{:.2f}"})
    write_numbers("drift", N, fmt)
    print()
    for k in sorted(N):
        print(f"    {k:<26} {N[k]}")
    print(f"""
  THE CHECKS A READER SHOULD MAKE, stated rather than left implicit.
    the no-drift rmse is {N['drift_rmse_nodrift']:.4f}, and 118 is defined as
    the count landing on 0.20, so this confirms the footing;
    the bias is {N['drift_bias']:.4f} against -r T / 2 =
    {N['drift_bias_span_pred']:.4f}, which is the ceiling argument: the bias is
    set by the record's SPAN of {N['drift_span_yr']:.2f} years and not by its count.""")
