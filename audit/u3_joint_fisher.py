"""
The joint Fisher precision on ln theta from a marker record PLUS k breath
measurements, and the smallest k that reaches the 0.20 gate.

THE QUESTION, as posed. Not whether a breath test alone identifies theta, and
not what a two-point ratio is worth on its own. How few breath measurements a
MARKER RECORD needs, which is a joint-information question and is not answered
by propagating a ratio through the fibre. audit/u2_breath_functional.py does
the propagation and states its two biases; this script replaces it.

THE PARAMETER VECTOR is (ln theta, ln K_m, b, ln k_e, ln c). The fifth entry
is the breath test's unknown subject-specific recovery factor, carried FREE
with no prior. That is what makes c cancel: a single breath reading then
contributes exactly nothing, because ln c absorbs it, and k readings
contribute through their k-1 contrasts. Nothing in this script assumes the
cancellation; it falls out of carrying ln c as a free parameter.

THE TWO BLOCKS.
  marker   residual (theta + Delta_i) g(m_i) + k_e m_i - (L_i + b), ordinate
           noise sigma from model_A.run.sigma_total re-solved at each
           capacity. Columns are g_design_common.columns's, and this script
           reproduces that module's sd(ln theta) as its anchor check.
  breath   ln y_j = ln c + ln F(m_j), F = theta/(theta + k_e(K_m + m_j)),
           derived in audit/u2_breath_functional.py part A, with a relative
           CV per reading.

TWO CONVENTIONS ON THE BREATH VISIT'S LEVEL, both reported, because they
bracket the answer and the difference is not small.
  measured  m_j is treated as a known covariate, so F depends on the
            parameters only explicitly. This is the author's own ratio
            convention, "two KNOWN marker levels". It ignores the assay error
            on m_j.
  modelled  m_j is the quasi-steady level the parameters imply, so the
            implicit derivative dm/dparam enters. This is the full likelihood
            and it is the primary number.

PLACEMENT. Breath measurements are attached to visits the marker record
already contains, so they cost an assay and not a blood draw. They are placed
at the dose cycle's trough and peak in alternation, which audit/u2 part E
establishes is the widest lever available.

FOOTING. theta = 3000 nominal with the RECORDED LOAD HELD FIXED at
L = 322.727 umol/L/d and the operating level, sigma and the load-bias prior
re-solved at each capacity, which is Table tab:scope's convention. K_m = 100,
k_e = 5/d, k_e log-prior 0.8, load-bias prior 0.05 L, D = 6000 on a 14 d train
with a 3.5 d half-life. Marker record: 25 monthly draws unless stated. Load
CV 0 for the headline table and 10% with 24 seeds for the spread.

Nothing here records an expected result in advance.

Run:  python audit/u3_joint_fisher.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from design.run import P, delta, load, DOSE_INTERVAL, LAMBDA, D_AMP, LOAD_FRAC, ASSAY_CV

TH, KM, KE = P.theta, P.Km, P.ke
LB = float(load(0.0))
T = DOSE_INTERVAL
GATE = 0.20
KE_PRIOR = 0.8

DEC = float(np.exp(-LAMBDA * T))
DELTA_PEAK = D_AMP / (1.0 - DEC)
DELTA_TROUGH = DELTA_PEAK * DEC


# ---------------------------------------------------------------- model
def level(theta, inc, Ltot, Km=KM, ke=KE):
    inc = np.asarray(inc, dtype=float)
    Bq = ke * Km + theta + inc - Ltot
    return (-Bq + np.sqrt(Bq ** 2 + 4.0 * ke * Ltot * Km)) / (2.0 * ke)


def sigma_at(theta, mbar, load_frac=LOAD_FRAC, assay_cv=ASSAY_CV):
    slope = theta * KM / (KM + mbar) ** 2 + KE
    return float(np.hypot(load_frac * LB, slope * assay_cv * mbar))


def marker_rows(theta, m, inc):
    """d(residual)/d(ln theta, ln K_m, b, ln k_e, ln c), in flux units."""
    g = m / (KM + m)
    return np.stack([theta * g,
                     -(theta + inc) * KM * m / (KM + m) ** 2,
                     np.ones_like(m),
                     KE * m,
                     np.zeros_like(m)], axis=1)


def breath_rows(theta, m, inc, modelled=True):
    """
    d(ln y)/d(ln theta, ln K_m, b, ln k_e, ln c) for the derived functional
    F = theta / (theta + k_e (K_m + m)).

    `modelled` adds the implicit derivative of the quasi-steady level with
    respect to each parameter; `modelled=False` treats m as a known covariate.
    """
    m = np.asarray(m, dtype=float)
    inc = np.asarray(inc, dtype=float)
    S = theta + KE * (KM + m)
    d_lt = 1.0 - theta / S
    d_lK = -KE * KM / S
    d_b = np.zeros_like(m)
    d_lk = -KE * (KM + m) / S
    d_lc = np.ones_like(m)
    if modelled:
        dF_dm = -KE / S
        J = (theta + inc) * KM / (KM + m) ** 2 + KE       # dPhi/dm
        dm_dlt = -(theta * m / (KM + m)) / J
        dm_dlK = +((theta + inc) * KM * m / (KM + m) ** 2) / J
        dm_db = +1.0 / J
        dm_dlk = -(KE * m) / J
        d_lt = d_lt + dF_dm * dm_dlt
        d_lK = d_lK + dF_dm * dm_dlK
        d_b = d_b + dF_dm * dm_db
        d_lk = d_lk + dF_dm * dm_dlk
    return np.stack([d_lt, d_lK, d_b, d_lk, d_lc], axis=1)


def sd_joint(theta, m, inc, sigma, mb, incb, breath_cv, modelled=True,
             ke_prior=KE_PRIOR, load_prior=True, km_prior=None):
    """sd(ln theta) from the marker block, the breath block and the priors."""
    pri = [np.inf,
           km_prior ** 2 if km_prior is not None else np.inf,
           (LOAD_FRAC * LB) ** 2 if load_prior else np.inf,
           ke_prior ** 2 if ke_prior is not None else np.inf,
           np.inf]
    F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri])
    if len(m):
        Q = marker_rows(theta, m, inc)
        F = F + Q.T @ Q / sigma ** 2
    if len(mb):
        R = breath_rows(theta, mb, incb, modelled=modelled)
        F = F + R.T @ R / breath_cv ** 2
    else:
        #  no breath data means no ln c to estimate. Carrying the column with
        #  a flat prior and no row would make the matrix singular and report
        #  the marker record as uninformative, which is the opposite of true.
        F = F[:4, :4]
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def breath_design(k):
    """k readings alternating trough, peak, trough, ... at the train's ends."""
    if k <= 0:
        return np.array([]), np.array([])
    incs = np.array([DELTA_TROUGH if j % 2 == 0 else DELTA_PEAK
                     for j in range(k)], dtype=float)
    return incs, incs


def smallest_k(theta, m, inc, sigma, breath_cv, modelled=True, kmax=4000,
               gate=GATE, **kw):
    """
    The fewest breath readings that bring sd(ln theta) to the gate.

    Searched from k = 0 upward rather than bisected, because sd(ln theta) is
    not monotone in k at k = 1: one reading carries no information at all,
    the free recovery factor absorbing it exactly, so a bisection on a
    monotone assumption would step over the k = 0 answer.
    """
    for k in range(0, kmax + 1):
        incs, _ = breath_design(k)
        mb = level(theta, incs, LB) if k else np.array([])
        if sd_joint(theta, m, inc, sigma, mb, incs, breath_cv,
                    modelled=modelled, **kw) <= gate:
            return k
    return None


def record(theta, times, load_cv=0.0, seed=0):
    times = np.asarray(times, dtype=float)
    inc = np.asarray(delta(times), dtype=float)
    if load_cv:
        rng = np.random.default_rng(seed)
        s2 = np.log1p(load_cv ** 2)
        Ls = LB * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
    else:
        Ls = np.full(times.size, LB)
    return level(theta, inc, Ls), inc


# ======================================================================
def anchor():
    """Reproduce g_design_common's own number, as the one shared footing."""
    print("=" * 74)
    print("ANCHOR.  this script's marker block against g_design_common")
    print("=" * 74)
    import g_design_common as G
    t = 30.4 * np.arange(25, dtype=float)
    mine = []
    theirs = []
    for cv, seeds in ((0.0, [0]), (0.10, range(8))):
        for s in seeds:
            m, inc = record(TH, t, load_cv=cv, seed=s)
            mine.append(sd_joint(TH, m, inc, G.SIG, [], [], 0.05))
            mg, ig, _ = G.record(t, load_cv=cv, seed=s)
            theirs.append(G.sd_ln_theta(mg, ig))
    mine, theirs = np.array(mine), np.array(theirs)
    gap = np.abs(mine - theirs).max()
    print(f"\n  9 comparisons, constant load and 10% load CV over 8 seeds")
    print(f"  largest absolute difference: {gap:.3e}")
    print(f"  this script at the paper's own design: {mine[0]:.4f}")
    print(f"  g_design_common at the same:           {theirs[0]:.4f}")
    return gap


def part_a():
    print()
    print("=" * 74)
    print("PART A.  where the breath row points, against where the record is")
    print("=" * 74)
    print("\n  the record's information has an exact null direction (the")
    print("  fibre, audit/u1) and a second, nearly-null one. What decides")
    print("  whether breath helps is whether its row has weight on the")
    print("  SECOND. Cosines in the scaled coordinates (ln theta, ln K_m,")
    print("  b/L, ln k_e), with ln c projected out of the breath contrast.\n")
    t = 30.4 * np.arange(25, dtype=float)
    Sc = np.diag([1.0, 1.0, LB, 1.0])
    for theta in (3000.0, 600.0, 300.0):
        mbar = float(level(theta, 0.0, LB))
        sig = sigma_at(theta, mbar)
        m, inc = record(theta, t, 0.0)
        Q = marker_rows(theta, m, inc)[:, :4] @ Sc
        Fm = Q.T @ Q / sig ** 2
        w, V = np.linalg.eigh(Fm)
        incs, _ = breath_design(2)
        mb = level(theta, incs, LB)
        R = breath_rows(theta, mb, incs, modelled=True)
        contrast = (R[1] - R[0])[:4] @ Sc          # ln c drops out of it
        contrast = contrast / np.linalg.norm(contrast)
        print(f"  theta = {theta:6.0f}   lambda: "
              f"{np.array2string(w, precision=2)}")
        for j, nm in ((0, "null (fibre)"), (1, "second weakest")):
            print(f"      |cos(breath contrast, {nm:>14})| = "
                  f"{abs(float(V[:, j] @ contrast)):.4f}")
        print(f"      breath contrast = "
              f"{np.array2string(contrast, precision=4)}")
    print("\n  READING. Weight on the null direction is what closes the fibre.")
    print("  Weight on the second is what buys precision, and audit/u1 part D")
    print("  showed an exact K_m (pure null direction) buys about 5%.")


def part_b():
    print()
    print("=" * 74)
    print("PART B.  the joint precision against k, and the smallest k")
    print("=" * 74)
    t = 30.4 * np.arange(25, dtype=float)
    for cv in (0.05, 0.10, 0.14):
        print(f"\n  breath CV {cv:.0%} per reading, 25 monthly marker draws,"
              f" constant recorded load")
        print(f"    {'theta':>7} {'conv':>9} " +
              "".join(f"{f'k={k}':>9}" for k in (0, 1, 2, 4, 8, 16, 32, 64)) +
              f"{'min k':>8}")
        print("    " + "-" * 96)
        for theta in (3000.0, 600.0, 300.0):
            mbar = float(level(theta, 0.0, LB))
            sig = sigma_at(theta, mbar)
            m, inc = record(theta, t, 0.0)
            for modelled in (True, False):
                vals = []
                for k in (0, 1, 2, 4, 8, 16, 32, 64):
                    incs, _ = breath_design(k)
                    mb = level(theta, incs, LB) if k else np.array([])
                    vals.append(sd_joint(theta, m, inc, sig, mb, incs, cv,
                                         modelled=modelled))
                kmin = smallest_k(theta, m, inc, sig, cv, modelled)
                tag = "modelled" if modelled else "measured"
                print(f"    {theta:7.0f} {tag:>9} " +
                      "".join(f"{v:9.3f}" for v in vals) +
                      (f"{kmin:8d}" if kmin is not None else f"{'none':>8}"))


def part_c():
    print()
    print("=" * 74)
    print("PART C.  the same at the specification's 10% load CV, 24 seeds")
    print("=" * 74)
    print("\n  medians with the seed range in brackets, breath CV 10%,")
    print("  full likelihood ('modelled') only.\n")
    t = 30.4 * np.arange(25, dtype=float)
    cv = 0.10
    for theta in (3000.0, 600.0, 300.0):
        mbar = float(level(theta, 0.0, LB))
        sig = sigma_at(theta, mbar)
        print(f"  theta = {theta:.0f}")
        for k in (0, 2, 4, 8, 16, 32):
            incs, _ = breath_design(k)
            mb = level(theta, incs, LB) if k else np.array([])
            vv = np.array([
                sd_joint(theta, *record(theta, t, load_cv=0.10, seed=s),
                         sig, mb, incs, cv)
                for s in range(24)], float)
            print(f"    k = {k:3d}   {np.median(vv):8.4f} "
                  f"[{vv.min():.4f}, {vv.max():.4f}]"
                  + ("   <- gate" if np.median(vv) <= GATE else ""))
        print()


def part_bb():
    print()
    print("=" * 74)
    print("PART B2.  the smallest k against the marker record's own length")
    print("=" * 74)
    print("\n  monthly draws, constant recorded load, full likelihood.")
    print("  k = 0 means the record already meets the gate. A dash means the")
    print("  gate is not reached at 4000 breath readings.\n")
    for cv in (0.05, 0.10, 0.14):
        print(f"  breath CV {cv:.0%}")
        print(f"    {'n visits':>9} {'years':>6} " +
              "".join(f"{f'theta={t:.0f}':>12}" for t in (3000, 600, 300)))
        print("    " + "-" * 56)
        for n in (13, 25, 50, 100, 200):
            cells = []
            for theta in (3000.0, 600.0, 300.0):
                mbar = float(level(theta, 0.0, LB))
                sig = sigma_at(theta, mbar)
                m, inc = record(theta, 30.4 * np.arange(n, dtype=float), 0.0)
                kk = smallest_k(theta, m, inc, sig, cv)
                cells.append("-" if kk is None else str(kk))
            print(f"    {n:9d} {n * 30.4 / 365.25:6.1f} " +
                  "".join(f"{c:>12}" for c in cells))
        print()


def part_d():
    print()
    print("=" * 74)
    print("PART D.  breath alone, and the marker record alone, side by side")
    print("=" * 74)
    print("\n  'breath alone' drops the marker block entirely, so the record")
    print("  supplies nothing: that is the regime where the author's ratio")
    print("  propagation lives, and it needs A, B and C from somewhere.")
    print("  Reported as the joint number with n = 0 marker draws and the")
    print("  same priors, which is why it is finite at all.\n")
    t = 30.4 * np.arange(25, dtype=float)
    cv = 0.10
    print(f"    {'theta':>7} {'record only':>13} {'breath k=8 only':>17} "
          f"{'joint k=8':>11} {'joint k=32':>11}")
    print("    " + "-" * 64)
    for theta in (3000.0, 600.0, 300.0):
        mbar = float(level(theta, 0.0, LB))
        sig = sigma_at(theta, mbar)
        m, inc = record(theta, t, 0.0)
        rec = sd_joint(theta, m, inc, sig, [], [], cv)
        i8, _ = breath_design(8)
        mb8 = level(theta, i8, LB)
        i32, _ = breath_design(32)
        mb32 = level(theta, i32, LB)
        bre = sd_joint(theta, [], [], sig, mb8, i8, cv)
        j8 = sd_joint(theta, m, inc, sig, mb8, i8, cv)
        j32 = sd_joint(theta, m, inc, sig, mb32, i32, cv)
        print(f"    {theta:7.0f} {rec:13.3f} {bre:17.3f} {j8:11.3f} "
              f"{j32:11.3f}")


if __name__ == "__main__":
    g = anchor()
    part_a()
    part_b()
    part_bb()
    part_c()
    part_d()
    print()
    print("=" * 74)
    print(f"ANCHOR verdict: largest gap to g_design_common = {g:.3e}")
    print("=" * 74)
