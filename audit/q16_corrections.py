"""
The two Phase 10 corrections, computed rather than transcribed.

CORRECTION 1. S.3's ex vivo comparison was run at an assay width of 0.20. R4's
finding is that the ex vivo to in vivo mapping has never been measured in any
human, that the only quantitative link is Manoli's Spearman r of 0.682 in MMA
(R-squared 0.327), and that the defensible width is a natural-log sd of order 1,
about threefold per standard deviation. 0.20 is five times tighter than that and
tighter than the tightest value S5 swept, so the headline was quoted from a
column nobody can currently deliver. Part A redoes it with 1.0 primary.

CORRECTION 2. S.1 priced candidate 9 as removing `b` from the parameter vector.
Measuring the load does TWO things: it removes the offset AND it removes the
load-recording error from `sigma`. Part B separates them, and Part C reconciles
the result with the paper's own route 6 in Section 11, which claims thirteen
measured-load visits match twenty-five recorded ones, about 28% in sd. Three
numbers now exist for one route and they cannot all be right.

Note on route 6's provenance: `2.659`, `thirteen`, `twenty-five` and `eighteen`
are hardcoded in the prose of `paper/sections/11_scope.tex`, not `\\num{}` macros,
so the artefacts check never compared them with a generator. R23's failure mode.

Nothing here records an expected result in advance.

Run:  python audit/q16_corrections.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import n4_theta_scope as N4
import q1_increment_shape as Q
import q2_requirement_map as S2

KM, KE, TH_NOM = S2.KM, S2.KE, S2.TH_NOM
GATE = S2.GATE
#: nominal plus the disease range, since route 6 is a claim about nominal
CAPS = (1.00, 0.20, 0.10, 0.05)
#: R4's supported width, and the optimistic one S.3 was quoted from
ASSAYS = (1.0, 0.50, 0.30, 0.20)
_K = Q.kernel_exp(Q.TH_HALF0)
_TN = np.arange(0.0, 900.0, 14.0)
_WN = np.ones(int(np.ceil(900.0 / 14.0)))


def banner(t):
    print("\n" + "=" * 100)
    print("  " + t)
    print("=" * 100)


def inc_at(times):
    return Q.D_AMP * Q.assemble(_K, _TN, _WN, Q.at_draw(times))


def sd_theta(theta, times, seed, load_frac=None, free_b=True,
             ke_sd=S2.KE_SDS[0]):
    """
    sd(ln theta) on the PROJECT's own load parametrisation.

    `LOAD_FRAC` is the load record's ACCURACY, distinct from the 10% load
    VARIABILITY that S2 draws the loads with. The paper's Section 8 states the
    two behave oppositely, and the accuracy is ONE knob driving BOTH channels:
    `sigma = hypot(load_frac * L, b * assay_cv * m_bar)` and the offset prior
    `s_b = load_frac * L`. Setting it to 0 is measuring the load exactly, which
    leaves `sigma` at its assay term rather than at zero.

    `free_b` False instead removes the offset from the parameter vector while
    leaving `sigma` untouched. That is what S.1 priced, and it is only half of
    what measuring the load does.
    """
    lf = N4.LOAD_FRAC if load_frac is None else load_frac
    L, m_bar = N4.anchor(theta, KM, KE, "A")
    pA = N4.ParamsA(theta=theta, Km=KM, ke=KE)
    sigma = float(N4.LOC.sigma_total(lf, N4.ASSAY_CV, p=pA, m_bar=m_bar))
    s_b = lf * L
    Ls = S2.loads_for(times.size, seed, L)
    inc = inc_at(times)
    m = S2.level(theta, inc, Ls)
    g = m / (KM + m)
    cols = [theta * g, -(theta + inc) * KM * m / (KM + m) ** 2,
            np.ones_like(m), KE * m]
    keep, pri = [0, 1], [0.0, 0.0]
    if free_b and s_b > 0:
        keep.append(2)
        pri.append(1.0 / s_b ** 2)
    keep.append(3)
    pri.append(1.0 / ke_sd ** 2)
    Qm = np.stack([cols[i] for i in keep], axis=1)
    F = np.diag(pri) + Qm.T @ Qm / sigma ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 1e-14 * w[-1]:
        return np.inf
    return float(np.sqrt(np.linalg.inv(F)[0, 0]))


def cell(theta, times, **kw):
    v = np.array([sd_theta(theta, times, s, **kw) for s in S2.SEEDS])
    v = v[np.isfinite(v)]
    return float(np.median(v)) if v.size else np.inf


def monthly(n):
    """n monthly draws, the paper's own step, no baseline."""
    return np.arange(30.4, 30.4 * (n + 1) - 1e-9, 30.4)


# ==========================================================================
def part_a():
    banner("A. THE EX VIVO COMPARISON AT R4's SUPPORTED WIDTH")
    print("""  An ex vivo assay of log-sd `a` measures ln theta directly, so the joint precision
  is the quadrature sum: 1/sd^2 = 1/record^2 + 1/a^2. `record adds` is how much the
  two-year monthly record improves on the assay used alone. R4 supports a = 1.0;
  a = 0.20 is what S.3 was quoted from.""")
    print(f"\n    {'capacity':>9s} {'record alone':>13s}" +
          "".join(f"{f'adds, a={a:.2f}':>16s}" for a in ASSAYS))
    for cap in CAPS:
        r = cell(TH_NOM * cap, monthly(25))
        cells = []
        for a in ASSAYS:
            j = 1.0 / np.sqrt(1.0 / r ** 2 + 1.0 / a ** 2)
            cells.append(f"{100 * (1 - j / a):15.1f}%")
        print(f"    {100 * cap:8.0f}% {r:13.3f}" + "".join(cells))
    print("""
  And the joint figure itself, against the 0.20 gate, which is what decides whether
  either partner is usable:""")
    print(f"\n    {'capacity':>9s} {'assay alone':>12s} {'joint, a=1.0':>14s} "
          f"{'joint, a=0.20':>14s}")
    for cap in CAPS:
        r = cell(TH_NOM * cap, monthly(25))
        j1 = 1.0 / np.sqrt(1.0 / r ** 2 + 1.0)
        j2 = 1.0 / np.sqrt(1.0 / r ** 2 + 25.0)
        print(f"    {100 * cap:8.0f}% {1.0:12.3f} {j1:14.3f} {j2:14.3f}")


def part_b():
    banner("B. CANDIDATE 9's TWO CHANNELS, SEPARATED")
    print("""  Measuring the load removes the offset from the parameter vector AND removes the
  load-recording error from `sigma`. S.1 priced only the first. On the project's
  own knob the two are ONE parameter, the load record's accuracy `LOAD_FRAC`,
  whose committed value is 5%: it sets the offset prior `s_b = LOAD_FRAC * L` and
  one of `sigma`'s two terms. Measuring the load exactly sends it to 0, which
  leaves `sigma` at its assay term, not at zero.

  Improvement in sd(ln theta) against the committed design, 25 monthly draws.""")
    print(f"\n    {'capacity':>9s} {'committed':>10s} {'b dropped only':>15s}"
          + "".join(f"{f'accuracy {100*r:.0f}%':>17s}" for r in (0.02, 0.01, 0.0)))
    for cap in CAPS:
        th = TH_NOM * cap
        t = monthly(25)
        base = cell(th, t)
        cells = [f"{f'{100 * (1 - cell(th, t, free_b=False) / base):+.1f}%':>15s}"]
        for r in (0.02, 0.01, 0.0):
            v = cell(th, t, load_frac=r)
            cells.append(f"{f'{100 * (1 - v / base):+.1f}%':>17s}")
        print(f"    {100 * cap:8.0f}% {base:10.3f}" + "".join(cells))
    print("""
  `b dropped only` is what S.1 reported and it is the smaller half everywhere.
  The accuracy columns are what candidate 9 actually buys.""")


def part_c():
    banner("C. THE PAPER'S ROUTE 6 RECONCILED")
    print("""  Section 11 route 6: "thirteen measured-load visits match twenty-five recorded
  ones and six reach the target" at nominal, and the advantage shrinks "from
  thirteen visits to eighteen" as capacity falls. Checked on both readings: the
  `b`-only one its prose states, and the accuracy one its 2.659 implies.""")
    for cap, lab in ((1.00, "nominal"), (0.05, "5% of nominal")):
        th = TH_NOM * cap
        target = cell(th, monthly(25))
        print(f"\n  {lab.upper()}: 25 committed monthly draws give "
              f"sd(ln theta) = {target:.4f}")
        for kw, name in ((dict(free_b=False), "b dropped only (route 6's prose)"),
                         (dict(load_frac=0.0), "load measured exactly")):
            hit = next((n for n in range(2, 60)
                        if cell(th, monthly(n), **kw) <= target), None)
            print(f"    {name:34s} n = {hit if hit else '>59':>4} measured visits "
                  f"match 25   (at n=25: {cell(th, monthly(25), **kw):.4f})")
        for kw, name in ((dict(free_b=False), "b dropped only"),
                         (dict(load_frac=0.0), "load measured exactly")):
            hit = next((n for n in range(2, 400)
                        if cell(th, monthly(n), **kw) <= GATE), None)
            print(f"    {name:34s} n = {hit if hit else '>399':>4} measured visits "
                  f"reach the 0.20 gate")
    print("""
  Route 6 claims 13 against 25 at nominal, six reaching the gate, and 13 to 18 as
  capacity falls.""")


def part_d():
    banner("D. ROUTE 6's CHECKABLE FIGURE: 2.659 AT 5% OF NOMINAL")
    print("""  Routes 5 and 6 both quote 2.659 as what a perfect load record leaves at 5% of
  nominal. It is reproducible here, so it says which channels the paper's own
  computation moved.""")
    th = TH_NOM * 0.05
    t = monthly(25)
    for kw, name in (({}, "committed, LOAD_FRAC = 5%"),
                     (dict(free_b=False), "b dropped, sigma untouched"),
                     (dict(load_frac=0.01), "accuracy 1%"),
                     (dict(load_frac=0.0), "load measured exactly")):
        print(f"    {name:32s} {cell(th, t, **kw):8.4f}")
    print("""
  The row reading 2.659 is the computation the paper ran.""")


def part_e():
    banner("E. ROUTE 6's THIRD CLAUSE: DOES THE ADVANTAGE SHRINK BELOW NOMINAL")
    print("""  Route 6 says the advantage "shrinks rather than reverses as capacity falls, from
  thirteen visits to eighteen", and uses that to argue the route opens nothing
  below nominal. Swept across the whole range.""")
    print(f"\n    {'capacity':>9s} {'25 committed give':>18s} "
          f"{'measured visits to match':>26s}")
    for cap in (1.00, 0.50, 0.20, 0.10, 0.05, 0.02):
        th = TH_NOM * cap
        tgt = cell(th, monthly(25))
        n = next((n for n in range(2, 60)
                  if cell(th, monthly(n), load_frac=0.0) <= tgt), None)
        print(f"    {100 * cap:8.1f}% {tgt:18.4f} "
              f"{n if n else '>59':>26}")
    print("""
  If the column is flat, the advantage does not shrink and route 6's third clause
  is wrong, whatever its conclusion. The conclusion is separately checkable: 12
  visits at 5% of nominal still leave sd(ln theta) far outside the 0.20 gate, so
  "opens nothing below nominal" can be right while its supporting number is not.""")


def main():
    print("=" * 100)
    print("  THE TWO PHASE 10 CORRECTIONS")
    print("=" * 100)
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
