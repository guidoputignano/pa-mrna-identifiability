"""
S5. What an ex vivo enzyme assay on the patient's own cells would be worth.

Residual PCC activity can be measured directly, in cultured fibroblasts or
lymphocytes, and that measurement is what the clinical literature on propionic
acidemia actually reports. It is not the same quantity as the in vivo theta the
paper estimates: it is an activity per milligram of cell protein in one cell
type under assay conditions, and mapping it to a whole-body maximal velocity
carries the cell-type mix, the assay's own conditions and the protein content.
Brief R4 is the mapping, and it is UNRESOLVED. So the assay is carried here as
a LOG-NORMAL PRIOR on theta of a stated width, swept over the whole plausible
range of how bad that mapping is:

    sd(ln theta) = 0.3   the mapping is good to about a third
                   0.5
                   0.7
                   1.0   the mapping is worth little more than an order of
                         magnitude either way
    none            the record alone, which is S2's column

at 5, 10, 20 and 100% of normal capacity, ALONE and combined with the k_e
precision S2b says the micro-dose study can actually deliver:

    0.8      no study: the paper's own population prior on ln k_e
    0.054    study at 2% with a 5% within-patient drift, in quadrature, which
             S2b.5 calls the realistic case
    0.02     study at 2% with no drift, which S2b calls the optimistic one

AND THE RISK THE PRIOR CARRIES, which is the whole reason a sweep of its width
is not enough. A prior borrowed through an unresolved mapping may be CENTRED
WRONG, and a wrongly centred prior pulls the estimate with it in exact
proportion to how tight it is. Part D offsets the prior's centre by 0, 0.3 and
0.7 in ln theta and measures the bias, which is the quantity that decides
whether the assay should enter the estimator at all or stay a cross-check.

THE t = 0 CONVENTION AND THE FOOTING. As S2: dose-day draws pre-dose,
schedules on therapy, with the `baseline` footing adding one untreated draw at
t = 0. Convention A for the footing, `dose` for the amplitude.

Nothing here records an expected result in advance.

Run:  python audit/q5_exvivo_prior.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import least_squares

import q1_increment_shape as Q
import q2_requirement_map as S2

CAPS = (0.05, 0.10, 0.20, 1.00)
TH_SDS = (0.3, 0.5, 0.7, 1.0, None)          # None = the record alone
KE_CASES = (("no study, 0.8", 0.8),
            ("realistic, 0.054", float(np.hypot(0.02, 0.05))),
            ("optimistic, 0.02", 0.02))
SCHEDS = ("monthly", "weekly 7.3 d", "every 2 d")
OFFSETS = (0.0, 0.3, 0.7)
SEEDS = S2.SEEDS
GATE = S2.GATE


def banner(t):
    print("\n" + "=" * 104)
    print("  " + t)
    print("=" * 104)


# ==========================================================================
#  Fisher, with a prior on ln theta.
# ==========================================================================
def fisher_sd(theta, times, seed, ke_sd, load, th_sd):
    L, _, sigma, s_b = S2.footing(theta)
    Ls = S2.loads_for(times.size, seed, L)
    inc = S2.increment(times)
    m = S2.level(theta, inc, Ls)
    g = m / (S2.KM + m)
    cols = [theta * g, -(theta + inc) * S2.KM * m / (S2.KM + m) ** 2,
            np.ones_like(m), S2.KE * m]
    pri = [0.0 if th_sd is None else 1.0 / th_sd ** 2, 0.0,
           1.0 / s_b ** 2, None if ke_sd == 0 else 1.0 / ke_sd ** 2]
    keep = [0, 1]
    if load == "recorded":
        keep.append(2)
    if ke_sd != 0:
        keep.append(3)
    Qm = np.stack([cols[i] for i in keep], axis=1)
    F = np.diag([pri[i] for i in keep]) + Qm.T @ Qm / sigma ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def cell(theta, times, ke_sd, load, th_sd, seeds=SEEDS):
    v = np.array([fisher_sd(theta, times, s, ke_sd, load, th_sd)
                  for s in seeds])
    v = v[np.isfinite(v)]
    if not v.size:
        return np.inf, np.inf, np.inf
    return float(np.median(v)), float(v.min()), float(v.max())


# ==========================================================================
#  The MAP, with a prior on ln theta that may be offset.
# ==========================================================================
def mc_fit(m, inc, s, sigma, s_b, ke_sd, ke_centre, load, th_sd, th_centre):
    """S2's MAP with a log-normal prior on theta added, centred at
    `th_centre`. The grid stage carries the prior through a delta-method
    linearisation and the polish carries it exactly."""
    free_b = load == "recorded"
    free_ke = ke_sd != 0
    g = m[None, :] / (S2.KM_GRID[:, None] + m[None, :])
    cols, sd, mu = [g], [None if th_sd is None else th_sd * th_centre], \
        [th_centre]
    if free_ke:
        cols.append(np.broadcast_to(m, g.shape))
        sd.append(ke_sd * ke_centre)
        mu.append(ke_centre)
    if free_b:
        cols.append(np.ones_like(g))
        sd.append(s_b)
        mu.append(0.0)
    y = s[None, :] - inc[None, :] * g
    if not free_ke:
        y = y - ke_centre * m[None, :]
    X = np.stack(cols, axis=-1)
    p, J = Q._linear_solve(X, y, sigma, sd, mu)
    j = int(np.argmin(J))
    th0 = max(float(p[j, 0]), 1e-6)
    ke0 = float(p[j, 1]) if free_ke else ke_centre
    ke0 = ke0 if ke0 > 1e-3 * ke_centre else ke_centre
    b0 = float(p[j, -1]) if free_b else 0.0
    q0 = [np.log(th0), np.log(S2.KM_GRID[j])]
    lb, ub = [-np.inf, 0.0], [np.inf, np.log(1000.0)]
    if free_ke:
        q0.append(np.log(ke0)); lb.append(-np.inf); ub.append(np.inf)
    if free_b:
        q0.append(b0); lb.append(-np.inf); ub.append(np.inf)

    def res(q):
        th, K = np.exp(q[0]), np.exp(q[1])
        i = 2
        ke = np.exp(q[i]) if free_ke else ke_centre
        i += free_ke
        b = q[i] if free_b else 0.0
        r = [((th + inc) * m / (K + m) + ke * m + b - s) / sigma]
        if th_sd is not None:
            r.append([(q[0] - np.log(th_centre)) / th_sd])
        if free_ke:
            r.append([(q[2] - np.log(ke_centre)) / ke_sd])
        if free_b:
            r.append([b / s_b])
        return np.concatenate(r)

    q0 = np.clip(q0, np.array(lb) + 1e-9, np.array(ub) - 1e-9)
    out = least_squares(res, q0, bounds=(lb, ub), method="trf",
                        max_nfev=6000, xtol=1e-13, ftol=1e-13, gtol=1e-13)
    return float(np.exp(out.x[0]))


def mc_cell(theta, times, ke_sd, load, th_sd, offset=0.0, reps=96, seed0=0):
    """
    Error of ln theta_hat, with the ex vivo prior carried and possibly offset.

    The k_e prior is on S2b's `centre` convention, a per-patient measurement.
    The theta prior's centre is the truth times exp(offset), which is what a
    biased ex vivo mapping would supply, plus its own draw from the prior's own
    width, which is what its stated sd means.
    """
    L, _, sigma, s_b = S2.footing(theta)
    inc = S2.increment(times)
    err = []
    for r in range(reps):
        sd_ = seed0 + r
        rng = np.random.default_rng(10_007 * sd_ + 3)
        Ls = S2.loads_for(times.size, sd_, L)
        z = rng.normal(0.0, ke_sd) if ke_sd else 0.0
        ke_t, ke_c = S2.KE, S2.KE * np.exp(z)
        b_t = rng.normal(0.0, s_b) if load == "recorded" else 0.0
        m = S2.level(theta, inc, Ls + b_t, ke=ke_t)
        s = Ls + rng.normal(0.0, sigma, times.size)
        if th_sd is None:
            th_c = theta
        else:
            th_c = theta * np.exp(offset + rng.normal(0.0, th_sd))
        th = mc_fit(m, inc, s, sigma, s_b, ke_sd, ke_c, load, th_sd, th_c)
        err.append(np.log(th) - np.log(theta))
    e = np.array(err)
    boot = np.random.default_rng(99).choice(e, (400, e.size))
    rb = np.sqrt(np.mean(boot ** 2, axis=1))
    return dict(rmse=float(np.sqrt(np.mean(e ** 2))),
                bias=float(np.median(e)),
                rmse_lo=float(np.percentile(rb, 5)),
                rmse_hi=float(np.percentile(rb, 95)),
                robust=float((np.percentile(e, 84)
                              - np.percentile(e, 16)) / 2),
                collapse=int(np.sum(np.abs(e) > S2.COLLAPSE)), n=e.size)


# ==========================================================================
def part_a():
    banner("A. WHAT THE PRIOR IS WORTH, ALONE AND WITH THE RECORD")
    print("""  sd(ln theta-hat), Fisher, median over 24 load records, `*` at or below 0.20.
  Each `prior v` column is the prior AND the record together; what the assay
  delivers on its own is v itself, so a column whose entries exceed its own
  heading means the record is adding nothing. `record alone` is S2's column.
  On therapy, load recorded.""")
    out = {}
    for kname, ks in KE_CASES:
        print(f"\n  k_e: {kname}")
        for sn in SCHEDS:
            t = S2.times_for(sn, "no baseline")
            print(f"    {sn:13s}" +
                  "".join(f"{f'prior {v}' if v else 'record alone':>22s}"
                          for v in TH_SDS))
            for cap in CAPS:
                th = S2.TH_NOM * cap
                cells = []
                for v in TH_SDS:
                    c = cell(th, t, ks, "recorded", v)
                    out[(kname, sn, cap, v)] = c
                    cells.append(S2.fmt(c, 22))
                print(f"    {100 * cap:12.0f}%" + "".join(cells))
    return out


def part_b(out):
    banner("B. THE PRIOR AND THE RECORD, SEPARATED")
    print("""  THE TWO ADD IN QUADRATURE EXACTLY, and that is an identity rather than a
  measurement. The prior enters the information matrix only in its (ln theta,
  ln theta) element, and the marginal precision of ln theta is the SCHUR
  COMPLEMENT F00 - F0r Frr^-1 Fr0, in which the added term survives untouched
  because it sits outside every block the complement subtracts. So

      1 / sd(both)^2  =  1 / sd(record)^2  +  1 / prior^2

  holds to machine precision for every cell, and the `ratio` column below is a
  check on the arithmetic, not a finding. What it means is the substantive part:
  under this model an ex vivo assay and a marker record carry INDEPENDENT
  information about theta, so neither can spoil the other and the weaker one is
  never wasted. It also means the combination can be read off without running
  anything, which the columns confirm.""")
    for kname, ks in KE_CASES:
        print(f"\n  k_e: {kname}, monthly, load recorded")
        print(f"    {'cap':>5s} {'prior':>7s} {'record':>8s} {'both':>8s} "
              f"{'quadrature':>11s} {'ratio':>7s} {'meets':>6s}")
        for cap in CAPS:
            rec = out[(kname, "monthly", cap, None)][0]
            for v in (0.3, 0.5, 0.7, 1.0):
                both = out[(kname, "monthly", cap, v)][0]
                quad = 1.0 / np.sqrt(1.0 / v ** 2 + 1.0 / rec ** 2)
                print(f"    {100 * cap:4.0f}% {v:7.2f} {rec:8.3f} "
                      f"{both:8.3f} {quad:11.3f} {both / quad:7.3f} "
                      f"{'yes' if both <= GATE else 'no':>6s}")


def part_c(out):
    banner("C. WHERE THE PRIOR CHANGES A VERDICT")
    print("""  The cells S2 could not reach, re-asked with the prior. A cell counts as
  RESCUED if the record alone misses the gate and the record with the prior meets
  it, on the Fisher bound, and it is then re-run through the MAP estimator, 96
  records each, because S2 found the bound optimistic near the gate.""")
    rows = []
    for kname, ks in KE_CASES:
        for sn in SCHEDS:
            for cap in CAPS:
                rec = out[(kname, sn, cap, None)][0]
                if rec <= GATE:
                    continue
                for v in (0.3, 0.5, 0.7, 1.0):
                    if out[(kname, sn, cap, v)][0] <= GATE:
                        rows.append((kname, ks, sn, cap, v,
                                     out[(kname, sn, cap, v)][0], rec))
                        break
    print(f"  {len(rows)} cells are rescued on the bound.\n")
    print(f"    {'k_e':18s} {'sched':13s} {'cap':>5s} {'prior':>6s} "
          f"{'record':>7s} {'Fisher':>7s} {'RMSE [90%]':>22s} {'bias':>7s} "
          f"{'coll':>5s} {'verdict':>8s}")
    for kname, ks, sn, cap, v, fis, rec in rows:
        t = S2.times_for(sn, "no baseline")
        r = mc_cell(S2.TH_NOM * cap, t, ks, "recorded", v)
        print(f"    {kname:18s} {sn:13s} {100 * cap:4.0f}% {v:6.2f} "
              f"{rec:7.3f} {fis:7.3f} {r['rmse']:6.3f} "
              f"[{r['rmse_lo']:5.3f},{r['rmse_hi']:6.3f}] {r['bias']:+7.3f} "
              f"{r['collapse']:5d} {S2.verdict(r):>8s}")
    return rows


def part_d():
    banner("D. THE RISK: A PRIOR CENTRED WRONG")
    print("""  A prior borrowed through an unresolved mapping may be centred wrong, and a
  tight prior transmits that offset into the estimate almost undiminished. The
  prior's centre is offset by 0, 0.3 and 0.7 in ln theta, which is a mapping
  wrong by a factor of 1, 1.35 and 2. The MAP, 96 records, monthly, load
  recorded, k_e at the realistic 0.054. `pass-through` is the median bias divided
  by the offset: 1 means the estimate moved the whole way with the prior.""")
    for cap in CAPS:
        th = S2.TH_NOM * cap
        t = S2.times_for("monthly", "no baseline")
        ks = KE_CASES[1][1]
        print(f"\n    {100 * cap:.0f}% of normal")
        print(f"      {'prior':>6s}" +
              "".join(f"{f'offset {o}':>34s}" for o in OFFSETS))
        for v in (0.3, 0.5, 0.7, 1.0, None):
            cells = []
            for o in OFFSETS:
                r = mc_cell(th, t, ks, "recorded", v, offset=o)
                pt = "-" if (v is None or o == 0) else f"{r['bias'] / o:.2f}"
                cells.append(f"bias {r['bias']:+6.3f} RMSE {r['rmse']:5.3f} "
                             f"pass {pt:>5s}".rjust(34))
            lab = f"{v:.1f}" if v is not None else "none"
            print(f"      {lab:>6s}" + "".join(cells))


def main():
    print("=" * 104)
    print("  S5. THE EX VIVO PRIOR ON THETA")
    print("=" * 104)
    out = part_a()
    part_b(out)
    part_c(out)
    part_d()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
