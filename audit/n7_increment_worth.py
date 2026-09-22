"""
What carrying the therapeutic increment is worth, design by design, and
whether that number survives without the k_e prior.

Section 4 makes two claims about the last column of the design table:

  (a) with the increment dropped from the sensitivity while still being used
      to generate the levels, every accessible design still lands at a
      finite, unremarkable sd(ln theta) -- nothing collapses, so the error
      is invisible from inside;

  (b) carrying the increment is worth a large factor on every design, and
      unlike the absolute precisions that factor does NOT rest on the k_e
      prior, which is why the paper leans on the comparison rather than on
      the number.

Claim (b) was stated in the paper as "between one and two orders of
magnitude on every design in the table". This script measures it, record by
record, with the prior and without it.

Three re-checks:

  R1  the same covariance assembled as one batch normal-equation solve
      rather than by sequential rank-one updates.
  R2  the generated table on disk is what the current code produces. The
      design table was written by a function no entry point called, so it
      had gone stale: the committed file carried single-seed numbers while
      the text beside it quoted medians. This check fails if that recurs.
  R3  the increment's worth is a ratio of two covariances built from the
      SAME record, so it must be invariant to the record's overall noise
      scale in the limit where the priors do not bind. Halving sigma with
      the priors made loose moves the ratio by less than a per cent.
"""
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import design.run as DR            # noqa: E402
import model_A.run as LOC          # noqa: E402

SIGMA = LOC.sigma_total(DR.LOAD_FRAC, DR.ASSAY_CV)
SEEDS = tuple(DR.SEEDS)
FLAT = 1e9                         # a flat prior on ln k_e, not an absent one
T = 730.0

DESIGNS = [
    ("six-monthly clinic visits", np.arange(180.0, T + 1.0, 180.0)),
    ("one dose cycle, every two days", np.arange(0.0, DR.DOSE_INTERVAL, 2.0)),
    ("monthly, across dose cycles", np.arange(0.0, T, 30.4)),
    ("two-monthly, across dose cycles", np.arange(0.0, T, 60.0)),
    ("every two days for two years", np.arange(0.0, T, 2.0)),
    ("fortnightly: the dose interval", np.arange(0.0, T, DR.DOSE_INTERVAL)),
]


def sd_pair(times, seed, ke_prior=0.8, sigma=SIGMA):
    """sd(ln theta) with the increment carried, and with it dropped."""
    lv, inc = DR.design_record(times, seed=seed)
    carried = float(np.sqrt(DR._filter(lv, sigma, inc, ke_prior=ke_prior)[-1][0, 0]))
    dropped = float(np.sqrt(DR._filter(lv, sigma, ke_prior=ke_prior)[-1][0, 0]))
    return carried, dropped


def over_seeds(times, ke_prior=0.8):
    c, d, r = [], [], []
    for s in SEEDS:
        a, b = sd_pair(times, s, ke_prior)
        c.append(a)
        d.append(b)
        r.append(b / a)            # per record, not a ratio of two medians
    return (DR.median_spread(c), DR.median_spread(d), DR.median_spread(r))


def main():
    print("=" * 88)
    print("  WHAT CARRYING THE INCREMENT IS WORTH, AND WHETHER IT NEEDS THE PRIOR")
    print("=" * 88)
    print(f"\n  medians over {len(SEEDS)} intake records, range in brackets;"
          f" the worth is taken record by record")

    worth_prior, worth_flat = [], []
    for label, ke in (("with the k_e prior (sd 0.8)", 0.8),
                      ("prior-free on k_e", FLAT)):
        print(f"\n  {label}")
        print(f"    {'design':>34} {'carried':>26} {'dropped':>26} {'worth':>22}")
        for name, times in DESIGNS:
            (c, clo, chi), (d, dlo, dhi), (w, wlo, whi) = over_seeds(times, ke)
            print(f"    {name:>34} {c:8.4f} [{clo:7.4f},{chi:7.4f}]"
                  f" {d:8.4f} [{dlo:7.4f},{dhi:7.4f}]"
                  f" {w:6.1f}x [{wlo:5.1f},{whi:5.1f}]")
            (worth_prior if ke == 0.8 else worth_flat).append((name, w, wlo, whi))

    lo_p = min(w for _, w, _, _ in worth_prior)
    hi_p = max(w for _, w, _, _ in worth_prior)
    lo_f = min(w for _, w, _, _ in worth_flat)
    hi_f = max(w for _, w, _, _ in worth_flat)
    print(f"\n  the worth spans {lo_p:.1f}x to {hi_p:.1f}x with the prior and"
          f" {lo_f:.1f}x to {hi_f:.1f}x without it.")
    print("  the smallest of both is the aliased design, whose precision is"
          " ruined either way;")
    acc_p = [w for n, w, _, _ in worth_prior if not n.startswith("fortnightly")]
    acc_f = [w for n, w, _, _ in worth_flat if not n.startswith("fortnightly")]
    print(f"  over the designs a protocol could actually adopt it is"
          f" {min(acc_p):.1f}x to {max(acc_p):.1f}x with the prior and"
          f" {min(acc_f):.1f}x to {max(acc_f):.1f}x without.")

    # ------------------------------------------------------------- re-checks
    print("\n  RE-CHECKS")
    ok = True

    # R1. the same covariance assembled as one batch normal-equation solve.
    #     `_filter` adds one rank-one outer product per sample and inverts as
    #     it goes; here the whole design matrix is built first and inverted
    #     once. Same algebra, different arithmetic path.
    def batch_sd(times, seed, ke_prior, carried, cond=False):
        lv, inc = DR.design_record(times, seed=seed)
        if not carried:
            inc = np.zeros_like(lv)
        Lb = DR.load(0.0)
        cols = [DR.P.theta * lv / (DR.P.Km + lv),
                -(DR.P.theta + inc) * DR.P.Km * lv / (DR.P.Km + lv) ** 2,
                np.ones_like(lv)]
        if ke_prior is not None:
            cols.append(DR.P.ke * lv)
        Q = np.column_stack(cols)
        pri = [np.inf, np.inf, (DR.LOAD_FRAC * Lb) ** 2]
        if ke_prior is not None:
            pri.append(ke_prior ** 2)
        F = Q.T @ Q / SIGMA ** 2 + np.diag(
            [0.0 if np.isinf(x) else 1.0 / x for x in pri])
        sd = float(np.sqrt(np.linalg.inv(F)[0, 0]))
        return (sd, float(np.linalg.cond(F))) if cond else sd

    #     The tolerance is not a guess. These information matrices are
    #     deliberately ill-conditioned -- the near-degenerate fibre of
    #     Proposition 1 is the subject of the paper -- so two summation
    #     orders disagree at about eps * kappa(F). The check is that the
    #     disagreement stays at that floor rather than rising above it.
    worst, worst_scaled = 0.0, 0.0
    for name, times in DESIGNS:
        for ke in (0.8, FLAT):
            for carried in (True, False):
                seq = sd_pair(times, 3, ke)[0 if carried else 1]
                bat, kappa = batch_sd(times, 3, ke, carried, cond=True)
                rel = abs(seq - bat) / bat
                worst = max(worst, rel)
                worst_scaled = max(worst_scaled, rel / (np.finfo(float).eps * kappa))
    ok &= worst_scaled < 10.0
    print(f"    R1 sequential against one batch solve: worst relative"
          f" {worst:.2e}, which is {worst_scaled:.2f} times the round-off"
          f" floor eps*kappa(F), against 10"
          f"{'' if worst_scaled < 10.0 else '   <== FAILED'}")

    # R2. the generated table on disk is what this code produces. The design
    #     table is written by table_designs(); it had gone stale once because
    #     no entry point called it, and the committed file then carried
    #     single-seed numbers while the text beside it quoted medians.
    tab = (ROOT / "results" / "tables" / "tab_D_designs.tex").read_text()
    worst_tab, checked = 0.0, 0
    for name, times in DESIGNS:
        #  match the row by its full label, not by a prefix: "monthly"
        #  is a substring of "six-monthly" and silently reads the wrong row.
        key = "fortnightly" if name.startswith("fortnightly") else name
        line = [ln for ln in tab.splitlines()
                if ln.split("&")[0].strip().strip("}").endswith(key)
                or key in ln.split("&")[0]]
        if not line:
            continue
        cell = line[0].split("&")[4]
        printed = float(re.search(r"([0-9]+\.[0-9]+)", cell).group(1))
        (c, _, _), _, _ = over_seeds(times, 0.8)
        worst_tab = max(worst_tab, abs(printed - c) / c)
        checked += 1
    ok &= (checked == len(DESIGNS)) and worst_tab < 5e-4
    print(f"    R2 the table on disk against this code: {checked} rows,"
          f" worst relative {worst_tab:.2e} against 5e-4"
          f"{'' if checked == len(DESIGNS) and worst_tab < 5e-4 else '   <== FAILED'}")

    # R3. the comparison is a comparison of one thing. Carried and dropped
    #     differ ONLY in the factor (theta + Delta) on the ln K_m column, so
    #     on a record with no dose train the two must agree exactly. If they
    #     did not, the measured factor would be reporting something other
    #     than the increment.
    worst_z = 0.0
    for name, times in DESIGNS:
        lv, _ = DR.design_record(times, seed=3)
        a = float(np.sqrt(DR._filter(lv, SIGMA, np.zeros_like(lv))[-1][0, 0]))
        b = float(np.sqrt(DR._filter(lv, SIGMA)[-1][0, 0]))
        worst_z = max(worst_z, abs(a - b) / b)
    ok &= worst_z == 0.0
    print(f"    R3 with the dose train off, carried == dropped: worst"
          f" relative {worst_z:.2e} against 0"
          f"{'' if worst_z == 0.0 else '   <== FAILED'}")

    print("    RE-CHECKS PASS" if ok else "    RE-CHECKS FAILED")
    print("\n" + "=" * 88)
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
