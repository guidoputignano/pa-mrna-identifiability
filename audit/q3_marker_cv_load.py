"""
S3. The marker CV against the load record's accuracy, as a two-way table.

`results/tables/tab_D_loadprec` prices the load record's accuracy along ONE
axis, at the committed 5% marker CV, and moves both channels the accuracy feeds
(the ordinate noise sigma and the load-bias prior) together. The paper then
says the design reaches 0.20 "only if intake is recorded to better than about
10%". That sentence is conditional on a marker CV nobody has measured in this
population: brief R3 puts it somewhere between an assay's own repeatability and
the biological variation of an organic-acid marker in an outpatient, and the
repository's own 12% question (REVIEW_LOG R4) is unresolved. So the requirement
is a SURFACE, not a row, and this script measures the surface.

    marker CV        5, 10, 15, 20, 25, 30%          (brief R3, swept)
    load recorded to 1, 5, 10, 20, 40%               (tab_D_loadprec's own grid)

Both channels move with the load accuracy, as in the published table; the
marker CV moves only sigma, which is what a marker CV is.

THE CELL IS THE NUMBER OF VISITS TO 0.20, which is the quantity a clinic can
act on, with sd at the booked two-year record beside it. Schedules: monthly and
two-monthly.

CAPACITY. Normal, theta = 3000, and the working band's two edges, theta = 2472
and 18814, which Section 6.2 reports and `n5_scope_contribution.py` records as
82% and 627% of normal. The band is itself a function of the marker CV and the
load accuracy, so Part E re-solves both edges across the surface rather than
holding them. The DISEASE range is S2's, not S3's, and is not repeated here.

THE t = 0 CONVENTION, stated once. Dose-day draws are read PRE-DOSE
(`q1.DOSE_DAY_DRAW`), and the schedules here are ON THERAPY: the first draw is
one dose interval in, so no draw is an untreated baseline. That is S2's primary
footing. Part A's reproduction row is the exception, and is flagged there: the
published table's monthly record starts at t = 0, which on the pre-dose reading
is an untreated baseline draw.

Nothing here records an expected result in advance.

Run:  python audit/q3_marker_cv_load.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from common.models import ParamsA
import model_A.run as LOC

import n4_theta_scope as N4
import q1_increment_shape as Q
import q2_requirement_map as S2

MARKER_CVS = (0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
LOAD_FRACS = (0.01, 0.05, 0.10, 0.20, 0.40)
CAPS = {"band low, 2472": 2472.0, "normal, 3000": 3000.0,
        "band high, 18814": 18814.0}
SCHED = {"monthly": 30.4, "two-monthly": 60.8}
SEEDS = S2.SEEDS
GATE = S2.GATE
KE_PRIOR = 0.8          #: the paper's own population prior, held across S3
MAX_VISITS = 240        #: 20 years monthly; beyond this the cell reads "never"


def banner(t):
    print("\n" + "=" * 104)
    print("  " + t)
    print("=" * 104)


# ==========================================================================
#  The footing at a given marker CV and load accuracy.
# ==========================================================================
def footing(theta, cv, load_frac):
    """(L, m_bar, sigma, s_b). Convention A, as S2 and n4 use: the recorded
    load stays where theta = 3000 put it and the marker rises."""
    L, m_bar = N4.anchor(theta, S2.KM, S2.KE, "A")
    p = ParamsA(theta=theta, Km=S2.KM, ke=S2.KE)
    sigma = float(LOC.sigma_total(load_frac, cv, p=p, m_bar=m_bar))
    return float(L), float(m_bar), sigma, load_frac * float(L)


def times_on_therapy(step, n):
    """`n` draws, the first one dose interval in, none untreated."""
    return step * np.arange(1.0, n + 1.0)


def fisher_sd(theta, times, seed, cv, load_frac):
    """sd(ln theta-hat) in (ln theta, ln K_m, b, ln k_e), load RECORDED."""
    L, _, sigma, s_b = footing(theta, cv, load_frac)
    Ls = S2.loads_for(times.size, seed, L)
    inc = S2.increment(times)
    m = S2.level(theta, inc, Ls)
    g = m / (S2.KM + m)
    Qm = np.stack([theta * g, -(theta + inc) * S2.KM * m / (S2.KM + m) ** 2,
                   np.ones_like(m), S2.KE * m], axis=1)
    F = (np.diag([0.0, 0.0, 1.0 / s_b ** 2, 1.0 / KE_PRIOR ** 2])
         + Qm.T @ Qm / sigma ** 2)
    w = np.linalg.eigvalsh(F)
    if w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def cell(theta, times, cv, load_frac, seeds=SEEDS):
    v = np.array([fisher_sd(theta, times, s, cv, load_frac) for s in seeds])
    v = v[np.isfinite(v)]
    if not v.size:
        return np.inf, np.inf, np.inf
    return float(np.median(v)), float(v.min()), float(v.max())


def visits_to_gate(theta, step, cv, load_frac, cap=MAX_VISITS):
    """
    The smallest number of draws whose MEDIAN over the seed set meets the gate.

    sd falls with n but not always smoothly, because each added draw also adds a
    load record, so the search is a scan rather than a bisection. Returns None
    if the gate is not met by `cap` draws.
    """
    lo, hi = 2, cap
    #  a geometric scan first, so a cell that never meets the gate costs little
    n = 4
    while n <= cap:
        if cell(theta, times_on_therapy(step, n), cv, load_frac)[0] <= GATE:
            hi = n
            break
        lo = n
        n *= 2
    else:
        return None
    for k in range(lo + 1, hi + 1):
        if cell(theta, times_on_therapy(step, k), cv, load_frac)[0] <= GATE:
            return k
    return hi


# ==========================================================================
def part_a():
    banner("A. THE PUBLISHED ROW, REPRODUCED BEFORE THE SURFACE IS BUILT")
    print("""  tab_D_loadprec, at the committed 5% marker CV, theta = 3000, its own monthly
  record `arange(0, 730, 30.4)`, 24 intake records. NOTE THE CONVENTION: that
  record starts at t = 0, which on the pre-dose reading is an UNTREATED BASELINE
  draw. The rest of S3 is on therapy and has no such draw, so this row is the
  only one comparable with the published table.""")
    pub = {0.01: 0.1590, 0.05: 0.1699, 0.10: 0.1991, 0.20: 0.2876,
           0.40: 0.5049}
    monthly0 = np.arange(0.0, 730.0, 30.4)
    print(f"\n    {'load to':>8s} {'sigma':>8s} {'q3, published record':>26s} "
          f"{'tab_D_loadprec':>15s} {'rel':>8s}"
          f"   {'q3, on therapy':>22s}")
    for lf in LOAD_FRACS:
        sig = footing(3000.0, 0.05, lf)[2]
        c = cell(3000.0, monthly0, 0.05, lf)
        ot = cell(3000.0, times_on_therapy(30.4, monthly0.size), 0.05, lf)
        print(f"    {100 * lf:7.0f}% {sig:8.1f}   {c[0]:.4f} "
              f"[{c[1]:.4f},{c[2]:.4f}] {pub[lf]:15.4f} "
              f"{100 * (c[0] / pub[lf] - 1):+7.2f}%   {ot[0]:.4f} "
              f"[{ot[1]:.4f},{ot[2]:.4f}]")
    print("\n  And what the marker CV alone does to sigma at theta = 3000, "
          "load recorded to 10%:")
    print(f"    {'marker CV':>10s} {'sigma':>8s} {'assay term':>12s} "
          f"{'load term':>11s}")
    L, m_bar = N4.anchor(3000.0, S2.KM, S2.KE, "A")
    p = ParamsA(theta=3000.0, Km=S2.KM, ke=S2.KE)
    b = LOC.fit_ab(p, m_bar)[1]
    for cv in MARKER_CVS:
        sig = footing(3000.0, cv, 0.10)[2]
        print(f"    {100 * cv:9.0f}% {sig:8.1f} {b * cv * m_bar:12.1f} "
              f"{0.10 * L:11.1f}")


def part_b():
    banner("B. THE TWO-WAY SURFACE: sd(ln theta-hat) AT THE BOOKED TWO-YEAR "
           "RECORD")
    print("  Fisher, 24 intake records, median [min, max]. `*` at or below "
          f"{GATE:.2f}. On therapy, so")
    print("  monthly is 24 draws and two-monthly 12, both over two years.")
    out = {}
    for cname, th in CAPS.items():
        for sn, step in SCHED.items():
            n = int(730.0 // step)
            t = times_on_therapy(step, n)
            print(f"\n  {cname}, {sn}, n = {n}")
            print(f"    {'load to':>8s}" +
                  "".join(f"{f'CV {int(100 * c)}%':>22s}" for c in MARKER_CVS))
            for lf in LOAD_FRACS:
                cells = []
                for cv in MARKER_CVS:
                    c = cell(th, t, cv, lf)
                    out[(cname, sn, lf, cv)] = c
                    cells.append(S2.fmt(c, 22))
                print(f"    {100 * lf:7.0f}%" + "".join(cells))
    return out


def part_c():
    banner("C. THE HEADLINE: VISITS TO 0.20")
    print("  The smallest number of draws whose median over 24 intake records "
          "meets the gate.")
    print(f"  `never` means not by {MAX_VISITS} draws. In brackets, the "
          "calendar time that many draws")
    print("  takes on that schedule, in years.")
    out = {}
    for cname, th in CAPS.items():
        for sn, step in SCHED.items():
            print(f"\n  {cname}, {sn}")
            print(f"    {'load to':>8s}" +
                  "".join(f"{f'CV {int(100 * c)}%':>18s}" for c in MARKER_CVS))
            for lf in LOAD_FRACS:
                cells = []
                for cv in MARKER_CVS:
                    n = visits_to_gate(th, step, cv, lf)
                    out[(cname, sn, lf, cv)] = n
                    cells.append(f"{'never':>18s}" if n is None else
                                 f"{n:12d} ({n * step / 365.25:3.1f}y)")
                print(f"    {100 * lf:7.0f}%" + "".join(cells))
    return out


def part_d(out_b, out_c):
    banner("D. WHICH AXIS BINDS")
    print("""  Part C's ratio of visits is unusable where the slower cell never meets the gate,
  which is most of the surface, so the comparison is made three ways instead.
  (a) THE TWO-YEAR FRONTIER: the loosest setting on each axis that still meets
      the gate inside the booked 24 monthly draws, and the loosest that meets it
      at all by 240 draws.
  (b) sd AT A FIXED RECORD, 24 monthly draws, so both axes are compared at the
      same cost: what each axis multiplies sd by across its own sweep.
  (c) THE CROSSOVER, the marker CV at which sigma's assay term equals its load
      term, which is where the two axes trade places.""")
    for cname in CAPS:
        print(f"\n  (a) {cname}, monthly. `-` means no setting on that row "
              "qualifies.")
        print(f"      {'marker CV':>10s} {'loosest load, <=24 draws':>26s} "
              f"{'loosest load, <=240':>22s}")
        for cv in MARKER_CVS:
            in24 = [lf for lf in LOAD_FRACS
                    if out_b[(cname, "monthly", lf, cv)][0] <= GATE]
            ever = [lf for lf in LOAD_FRACS
                    if out_c[(cname, "monthly", lf, cv)] is not None]
            f24 = f"{100 * max(in24):.0f}%" if in24 else "-"
            fev = f"{100 * max(ever):.0f}%" if ever else "-"
            print(f"      {100 * cv:9.0f}% {f24:>26s} {fev:>22s}")
    print("\n  (b) sd at 24 monthly draws, the ratio across each axis:")
    print(f"      {'capacity':20s} {'load 40% / load 1%, CV held':>30s} "
          f"{'CV 30% / CV 5%, load held':>28s}")
    for cname in CAPS:
        g = out_b
        rl = [g[(cname, "monthly", 0.40, cv)][0] / g[(cname, "monthly", 0.01, cv)][0]
              for cv in MARKER_CVS]
        rc = [g[(cname, "monthly", lf, 0.30)][0] / g[(cname, "monthly", lf, 0.05)][0]
              for lf in LOAD_FRACS]
        print(f"      {cname:20s} {f'{min(rl):.2f} to {max(rl):.2f}':>30s} "
              f"{f'{min(rc):.2f} to {max(rc):.2f}':>28s}")
    print("""
  (c) the crossover. sigma = hypot(load_frac * L, b * CV * m_bar), so the two
  terms are equal at CV = load_frac * L / (b * m_bar):""")
    for cname, th in CAPS.items():
        L, m_bar = N4.anchor(th, S2.KM, S2.KE, "A")
        p_ = ParamsA(theta=th, Km=S2.KM, ke=S2.KE)
        b = LOC.fit_ab(p_, m_bar)[1]
        print(f"    {cname:18s} " + "  ".join(
            f"load {int(100 * lf):2d}%: CV {100 * lf * L / (b * m_bar):5.1f}%"
            for lf in LOAD_FRACS))


def part_e():
    banner("E. THE WORKING BAND ACROSS THE SURFACE")
    print("""  Section 6.2's band, theta in [2472, 18814], is quoted at one marker CV and one
  load accuracy. Both edges are re-solved here on the booked monthly two-year
  record, on therapy, by scanning theta and bracketing each crossing of the
  gate. `-` means the gate is not met at any theta on the scan.""")
    grid = np.geomspace(200.0, 3e5, 120)
    print(f"\n    {'load to':>8s}" +
          "".join(f"{f'CV {int(100 * c)}%':>24s}" for c in MARKER_CVS))
    t = times_on_therapy(30.4, 24)
    for lf in LOAD_FRACS:
        cells = []
        for cv in MARKER_CVS:
            vals = np.array([cell(th, t, cv, lf)[0] for th in grid])
            ok = vals <= GATE
            if not ok.any():
                cells.append(f"{'-':>24s}")
                continue
            i = np.where(ok)[0]
            cells.append(f"[{grid[i[0]]:7.0f},{grid[i[-1]]:8.0f}]".rjust(24))
        print(f"    {100 * lf:7.0f}%" + "".join(cells))
    print("""
  The scan is on a 120-point log grid from 200 to 3e5, so each edge is located
  to about 6%; it is the band's MOVEMENT across the surface that is the result,
  not a sharper edge than Section 6.2's.""")


def main():
    print("=" * 104)
    print("  S3. THE MARKER CV AGAINST THE LOAD RECORD'S ACCURACY")
    print("=" * 104)
    part_a()
    out_b = part_b()
    out_c = part_c()
    part_d(out_b, out_c)
    part_e()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
