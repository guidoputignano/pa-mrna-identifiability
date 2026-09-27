"""
The k_e requirement frontier, and whether a timed urine collection can meet it.

WHY THIS REPLACES THE K_m BYPASS. audit/u1_fibre_bypass.py found that closing
the fibre buys about five per cent and that k_e is the load-bearing external
quantity, worth a factor of 2 to 9 depending on footing. That turns the scope
limit into a REQUIREMENT rather than a verdict: capacity is recoverable in the
disease range if and only if the non-PCC clearance is known well enough, and
the number that matters is how well.

PART A states the frontier: the k_e log-prior width at which the record
reaches sd(ln theta) = 0.20, per capacity, on two records. It is solved by
bisection on the width rather than read off a sweep, because the quantity
wanted is a threshold.

PART B gives the curve behind it, so the frontier is visible as a crossing
rather than asserted as a number.

PART C asks the modelling question the frontier raises, and it is NOT
rhetorical: is the model's k_e the same object as the urinary clearance a
timed collection measures? The answer decides whether the frontier is
reachable by a routine test or not reachable at all, and this script does not
assume they are equal. Three obstructions are priced separately.

FOOTING. As audit/PHASE11_KM_BYPASS.md section 0. theta = 3000 nominal with
the RECORDED LOAD HELD FIXED at L = 322.727 umol/L/d and the operating level,
sigma and the load-bias prior re-solved at each capacity, which is
Table tab:scope's convention. K_m = 100, k_e = 5/d, D = 6000 on a 14 d train
with a 3.5 d half-life, load recorded to 5% with a day-to-day intake CV of
10%, load-bias prior 0.05 L, gate 0.20, medians over 24 intake records with
the seed range printed. The Tier 1 record is 118 monthly draws with the
increment's scale carried at the published 14.02%, which is
audit/t2_tier1_margin.py's construction and the arm that needs a decade.

Nothing here records an expected result in advance.

Run:  python audit/u9_ke_frontier.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from design.run import (ASSAY_CV, DOSE_INTERVAL, D_AMP, LAMBDA, LOAD_CV,
                        LOAD_FRAC, load)

KM, KE = 100.0, 5.0
LB = float(load(0.0))
GATE = 0.20
MONTH = 30.4
YEAR = 365.25
KE_PRIOR = 0.8
PUBLISHED_KAPPA_CV = 0.1402
T_DOSE = 50.0 * DOSE_INTERVAL
PEAK = D_AMP / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))
SEEDS = range(24)
TIER1_N = 118


def op_level(theta):
    B = KE * KM + theta - LB
    return (-B + np.sqrt(B ** 2 + 4.0 * KE * LB * KM)) / (2.0 * KE)


def sigma_at(theta, assay_cv=ASSAY_CV):
    mbar = op_level(theta)
    slope = theta * KM / (KM + mbar) ** 2 + KE
    return float(np.hypot(LOAD_FRAC * LB, slope * assay_cv * mbar))


def levels(theta, n, seed, step=MONTH):
    t = step * np.arange(n, dtype=float)
    tau = np.mod(t - T_DOSE, DOSE_INTERVAL)
    inc = PEAK * np.exp(-LAMBDA * tau)
    Ls = LB * (1.0 + LOAD_CV *
               np.random.default_rng(seed).standard_normal(t.size))
    B = KE * KM + theta + inc - Ls
    return (-B + np.sqrt(B ** 2 + 4.0 * KE * Ls * KM)) / (2.0 * KE), inc


def sd_ln_theta(theta, m, inc, sigma, ke_prior=KE_PRIOR, kappa_cv=None):
    gm = m / (KM + m)
    cols = [theta * gm,
            -(theta + inc) * KM * m / (KM + m) ** 2,
            np.ones_like(m),
            KE * m]
    pri = [np.inf, np.inf, (LOAD_FRAC * LB) ** 2,
           ke_prior ** 2 if ke_prior is not None else np.inf]
    if kappa_cv is not None:
        cols.append(inc * gm)
        pri.append(kappa_cv ** 2)
    Q = np.stack(cols, axis=1)
    F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri]) \
        + Q.T @ Q / sigma ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def med(theta, n, ke_prior, kappa_cv=None):
    sig = sigma_at(theta)
    vals = np.array([sd_ln_theta(theta, *levels(theta, n, s), sig,
                                 ke_prior=ke_prior, kappa_cv=kappa_cv)
                     for s in SEEDS], float)
    f = vals[np.isfinite(vals)]
    if not f.size:
        return np.inf, np.inf, np.inf
    return float(np.median(f)), float(f.min()), float(f.max())


def frontier(theta, n, kappa_cv=None, lo=1e-4, hi=5.0, tol=1e-6):
    """
    The k_e log-prior width at which the record lands exactly on the gate.

    Bisected on the width. sd(ln theta) is increasing in the width, so the
    bracket is [tight, loose]; if the loose end already meets the gate the
    existing prior suffices and there is no requirement to state.
    """
    if med(theta, n, hi, kappa_cv)[0] <= GATE:
        return None, "the existing prior already suffices"
    if med(theta, n, lo, kappa_cv)[0] > GATE:
        return np.nan, "unreachable: even an exact k_e misses the gate"
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if med(theta, n, mid, kappa_cv)[0] <= GATE:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), ""


# ======================================================================
def part_a():
    print("=" * 78)
    print("PART A.  the k_e requirement frontier")
    print("=" * 78)
    print("\n  The k_e log-prior width at which the record reaches the 0.20")
    print("  gate. Below the stated width the gate is met, above it is not.")
    print("  'suffices' means the paper's own 0.8 prior already reaches it.\n")
    arms = [("Tier 1: 118 monthly, increment at 14.02%", TIER1_N,
             PUBLISHED_KAPPA_CV),
            ("the recommended record: 25 monthly, increment exact", 25, None)]
    out = {}
    for name, n, kap in arms:
        print(f"  {name}")
        print(f"    {'theta':>8} {'% of nominal':>13} {'required sd(ln k_e)':>20} "
              f"{'as a %':>9} {'factor vs 0.8':>14}")
        print("    " + "-" * 70)
        for theta in (3000.0, 600.0, 300.0, 150.0):
            s, note = frontier(theta, n, kap)
            out[(name, theta)] = s
            if s is None:
                print(f"    {theta:8.0f} {100 * theta / 3000:12.1f}% "
                      f"{'0.8 suffices':>20} {'':>9} {'':>14}")
            elif np.isnan(s):
                print(f"    {theta:8.0f} {100 * theta / 3000:12.1f}% "
                      f"{'unreachable':>20} {'':>9} {'':>14}")
            else:
                print(f"    {theta:8.0f} {100 * theta / 3000:12.1f}% "
                      f"{s:20.4f} {100 * s:8.2f}% {KE_PRIOR / s:14.1f}")
        print()
    print("  READING. This is the scope limit restated as a requirement. It")
    print("  is not 'capacity is unrecoverable in the severe patient'. It is")
    print("  'capacity is recoverable there if and only if the non-PCC")
    print("  clearance is known to the width in this table', which is a")
    print("  statement a reader can act on and the previous one is not.")
    return out


def part_b():
    print()
    print("=" * 78)
    print("PART B.  the curve behind the frontier")
    print("=" * 78)
    print("\n  sd(ln theta) against the k_e log-prior width, Tier 1 record,")
    print("  medians over 24 intake records. The frontier is where each")
    print("  column crosses 0.20.\n")
    widths = (0.8, 0.4, 0.2, 0.1, 0.05, 0.02, 0.01, 1e-4)
    print(f"    {'sd(ln k_e)':>11} " +
          "".join(f"{f'theta={t:.0f}':>22}" for t in (3000, 600, 300)))
    print("    " + "-" * 78)
    for s in widths:
        cells = []
        for theta in (3000.0, 600.0, 300.0):
            m, lo, hi = med(theta, TIER1_N, s, PUBLISHED_KAPPA_CV)
            cells.append(f"{m:.4f} [{lo:.4f},{hi:.4f}] ")
        print(f"    {s:11.4f} " + "".join(f"{c:>22}" for c in cells))
    print("\n  The theta = 3000 column never crosses because it starts below")
    print("  the gate; the other two cross, and where they cross is PART A.")


# ======================================================================
def part_c(front):
    print()
    print("=" * 78)
    print("PART C.  is the model's k_e the urinary clearance a collection measures")
    print("=" * 78)
    print("""
  THE QUESTION IS NOT RHETORICAL AND THE ANSWER IS NO, ON THREE COUNTS.
  A timed urine collection with paired plasma is the standard way to measure
  a clearance, it is routine, and it aims at the binding quantity. What it
  returns is CL_r, a renal clearance of the MARKER, in volume per time. The
  model's k_e is a first-order rate constant on the whole non-PCC route, in
  reciprocal time. Three separate things stand between them.

  (i) UNITS. k_e = CL / V_d, so the conversion needs a volume of
      distribution for 3-HP, which is not measured either. This repository
      already performs exactly this conversion for its own analogue: the
      k_e row of src/common/params.py takes methylmalonic acid at
      103 mL/min and divides by a volume of distribution of 0.48 L/kg to
      reach a rate. The V_d is doing as much work as the clearance.

  (ii) CHANNEL COVERAGE. src/common/params.py registers k_e with
      kind=LUMPED, and that file's own taxonomy says a lumped constant
      "absorbs several distinct mechanisms into one number" and therefore
      "cannot be checked against a direct measurement even in principle".
      A urinary collection sees the renal channel. Any non-renal disposal of
      3-HP sits inside k_e and outside the measurement.

  (iii) SPECIES. The parent model's alternative channel is carnitine
      conjugation, which esterifies propionyl-CoA to propionylcarnitine.
      That is a different molecule from 3-HP. So the channel the parent says
      dominates the alternative route does not remove the marker at all, and
      a 3-HP collection cannot see it even in principle.

  WHAT THAT LEAVES. CL_r / V_d is a LOWER BOUND on k_e, not an estimate of
  it, and the frontier is a requirement on PRECISION rather than on value.
  Writing k_e = CL_r / (V_d * f) with f the renal share of the non-PCC
  route, the three uncertainties add in quadrature:

      sd(ln k_e)^2 = sd(ln CL_r)^2 + sd(ln V_d)^2 + sd(ln f)^2

  so each term alone must be below the frontier. The table below reads that
  requirement off part A.
""")
    for name, theta in [(k[0], k[1]) for k in front
                        if k[1] in (600.0, 300.0)
                        and k[0].startswith("Tier 1")]:
        s = front[(name, theta)]
        if s is None or (isinstance(s, float) and np.isnan(s)):
            continue
        print(f"  theta = {theta:.0f}  ({100 * theta / 3000:.0f}% of nominal), "
              f"frontier sd(ln k_e) = {s:.4f}")
        print(f"    each of the three terms alone must be under "
              f"{100 * s:.2f}%")
        print(f"    and if all three were EQUAL, each must be under "
              f"{100 * s / np.sqrt(3):.2f}%")
        print()
    print("  So the modelling question decides the costing, and it decides it")
    print("  against the route as a way of REACHING the disease-range")
    print("  frontier: a volume of distribution for 3-HP known to a few per")
    print("  cent, and a renal share of the non-PCC route known to a few per")
    print("  cent, are both further from the published record than the")
    print("  clearance itself.")
    print()
    print("  WHAT THE ROUTE IS STILL WORTH, and it is not nothing:")
    print("    - it bounds k_e from BELOW, where the methylmalonic-acid")
    print("      analogue in params.py bounds it from above, and")
    print("      Section sec:scopeprec shows the working band moving across a")
    print("      factor of eight in the VALUE of k_e, which is a different")
    print("      question from its precision and one a bound does address;")
    print("    - the renal share f is itself measurable, as the fraction of")
    print("      an administered tracer recovered in urine, which is the same")
    print("      collection;")
    print("    - nothing here needs new technology, which is the property")
    print("      that made it worth asking about.")
    print()
    print("  IT IS NOT ON THE GATE LIST AS A LITERATURE NUMBER. The gate")
    print("  list wants 'the achievable precision on a timed collection',")
    print("  and that is only the first of three terms. What goes on the")
    print("  list is the pair: a volume of distribution for 3-HP, and the")
    print("  renal share of its non-PCC disposal. Neither was located here,")
    print("  and this container cannot check.")


if __name__ == "__main__":
    front = part_a()
    part_b()
    part_c(front)
