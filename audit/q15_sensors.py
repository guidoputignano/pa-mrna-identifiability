"""
Phase 10 steps 2 to 4: the five candidates that survive step 1, priced.

Step 1 dropped five of ten candidates because the model has no propionyl-CoA
state. The five that remain are candidate 2 (ex vivo PCC activity), 3 (the plasma
3-HP record itself), 8 (IV 13C3-3-HP), 9 (propionate turnover) and 10 (timed
urine). This prices them.

THE FRAMING THAT MAKES THEM COMPARABLE. The record's own information about
`(ln theta, ln K_m, b, ln k_e)` is a fixed 4x4 block. Every one of the four
remaining sensors does the same kind of thing to it: it supplies an independent
measurement of ONE of those four parameters, which adds a precision to one
diagonal entry. So the phase is not "which sensor is best" but **which parameter
is worth buying precision in, and at what price**, and the four sensors happen to
buy three different parameters:

    candidate 2,  ex vivo PCC activity  ->  ln theta   directly
    candidate 8,  IV 13C3-3-HP          ->  ln k_e     (which is k_e/alpha, R.2)
    candidate 9,  propionate turnover   ->  b          (it measures L + b)
    candidate 10, timed urine           ->  ln k_e     partially, f_renal * k_e

CANDIDATE 9 IS THE ROW NO PHASE HAS PRICED, and the framing shows why it was
missed. S2's requirement map has a `recorded` load column and a `measured` load
column, and treats them as two different studies. **They are the two endpoints of
candidate 9**: `recorded` is a turnover measurement of infinite width, `measured`
is one of zero width. Everything in between was never computed, and the sensor
lives there.

Nothing here records an expected result in advance.

Run:  python audit/q15_sensors.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import least_squares

import q2_requirement_map as S2

KM, KE, TH_NOM = S2.KM, S2.KE, S2.TH_NOM
GATE = S2.GATE
CAPS = (0.05, 0.10, 0.20)
SCHEDS = ("monthly", "weekly 7.3 d", "every 2 d")
SEEDS = S2.SEEDS
#: the population prior on ln k_e, which is what the record has without any sensor
KE_PRIOR = 0.8
#: fractional precisions a sensor might plausibly deliver
PRECS = (np.inf, 0.50, 0.30, 0.20, 0.10, 0.05, 0.02, 0.0)
_INC = {}


def banner(t):
    print("\n" + "=" * 104)
    print("  " + t)
    print("=" * 104)


def inc_of(times):
    k = (times.size, float(times[0]), float(times[-1]))
    if k not in _INC:
        _INC[k] = S2.increment(times)
    return _INC[k]


def sd_theta(theta, times, seed, sd_th=np.inf, sd_ke=KE_PRIOR, sd_b=None):
    """
    sd(ln theta) from the record plus whatever the sensors supply.

    `sd_th` is candidate 2's precision on ln theta, `sd_ke` the prior on ln k_e
    that candidates 8 and 10 tighten, `sd_b` the prior on the offset that
    candidate 9 tightens. `sd_b = None` means the record's own `s_b`, which is
    S2's `recorded` column; `sd_b = 0` is its `measured` column.
    """
    if sd_th == 0.0:
        return 0.0              # candidate 2 exact: theta is measured, not fitted
    L, _, sigma, s_b = S2.footing(theta)
    Ls = S2.loads_for(times.size, seed, L)
    inc = inc_of(times)
    m = S2.level(theta, inc, Ls)
    g = m / (KM + m)
    cols = [theta * g, -(theta + inc) * KM * m / (KM + m) ** 2,
            np.ones_like(m), KE * m]
    #  candidate 9 supplies an INDEPENDENT measurement of L + b, so its precision
    #  adds to the record's own offset prior rather than replacing it
    if sd_b is None:
        pb = 1.0 / s_b ** 2
    elif sd_b <= 0.0:
        pb = np.inf
    else:
        pb = 1.0 / s_b ** 2 + 1.0 / (sd_b * L) ** 2
    keep = [0, 1]
    pri = [0.0 if not np.isfinite(sd_th) else 1.0 / sd_th ** 2, 0.0]
    if np.isfinite(pb):
        keep.append(2)
        pri.append(pb)
    if np.isfinite(sd_ke) and sd_ke > 0:
        keep.append(3)
        pri.append(1.0 / sd_ke ** 2)
    elif sd_ke == 0:
        pass                        # k_e exactly known: the column drops out
    else:
        keep.append(3)
        pri.append(0.0)
    Qm = np.stack([cols[i] for i in keep], axis=1)
    F = np.diag(pri) + Qm.T @ Qm / sigma ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 1e-14 * w[-1]:
        return np.inf
    return float(np.sqrt(np.linalg.inv(F)[0, 0]))


def cell(theta, times, **kw):
    v = np.array([sd_theta(theta, times, s, **kw) for s in SEEDS])
    v = v[np.isfinite(v)]
    return float(np.median(v)) if v.size else np.inf


def mark(v):
    return f"{v:8.3f}{'*' if v <= GATE else ' '}"


# ==========================================================================
def step_2():
    banner("STEP 2. EACH SENSOR ALONE, AGAINST ITS OWN PRECISION")
    print("""  sd(ln theta), median over 24 load seeds, `*` at or inside the 0.20 gate. The
  first row of each block is the record with no sensor at all. Each sensor is then
  swept over the precision it might deliver, everything else left as the record
  has it.""")
    for cap in CAPS:
        th = TH_NOM * cap
        print(f"\n  CAPACITY {100 * cap:.0f}% OF NOMINAL")
        print(f"    {'sensor':38s} {'precision':>10s}" +
              "".join(f"{s:>12s}" for s in SCHEDS))
        base = {s: cell(th, S2.times_for(s, "no baseline")) for s in SCHEDS}
        print(f"    {'no sensor, record alone':38s} {'':>10s}" +
              "".join(mark(base[s]) + "   " for s in SCHEDS))
        specs = [("candidate 2, ex vivo PCC activity", "sd_th", True),
                 ("candidate 8/10, k_e sensor", "sd_ke", True),
                 ("candidate 9, propionate turnover", "sd_b", False)]
        for lab, key, is_log in specs:
            for p in PRECS:
                if key == "sd_ke" and not np.isfinite(p):
                    p = KE_PRIOR              # the no-sensor value, not infinity
                    tag = f"{p:.2f} (prior)"
                elif not np.isfinite(p):
                    tag = "none"
                elif p == 0.0:
                    tag = "exact"
                else:
                    tag = f"{100 * p:.0f}%"
                kw = {key: p}
                cells = "".join(
                    mark(cell(th, S2.times_for(s, "no baseline"), **kw)) + "   "
                    for s in SCHEDS)
                print(f"    {lab:38s} {tag:>10s}" + cells)
            print()


def step_3():
    banner("STEP 3. THE COMBINATIONS, AND WHICH ARE REDUNDANT")
    print("""  Two sensors at a plausible precision each, against the same two alone. A pair
  whose joint result equals the better of the two alone is REDUNDANT: the second
  sensor bought nothing. Monthly and every 2 d, at 5 and 10% capacity.

  Precisions used: candidate 2 at 20% (an ex vivo assay's repeatability),
  candidate 8 at 5.4% (P.8's four micro-dose studies, alpha's closure assumed
  free), candidate 9 at 10% (a turnover study at the load record's own CV).""")
    P2, P8, P9 = 0.20, 0.054, 0.10
    combos = [("record alone", {}),
              ("c2 only, 20%", dict(sd_th=P2)),
              ("c8 only, 5.4%", dict(sd_ke=P8)),
              ("c9 only, 10%", dict(sd_b=P9)),
              ("c2 + c8", dict(sd_th=P2, sd_ke=P8)),
              ("c2 + c9", dict(sd_th=P2, sd_b=P9)),
              ("c8 + c9", dict(sd_ke=P8, sd_b=P9)),
              ("c2 + c8 + c9", dict(sd_th=P2, sd_ke=P8, sd_b=P9))]
    for cap in (0.05, 0.10):
        th = TH_NOM * cap
        print(f"\n  CAPACITY {100 * cap:.0f}%")
        print(f"    {'combination':22s}" +
              "".join(f"{s:>14s}" for s in ("monthly", "every 2 d")))
        for lab, kw in combos:
            cells = "".join(
                mark(cell(th, S2.times_for(s, "no baseline"), **kw)) + "     "
                for s in ("monthly", "every 2 d"))
            print(f"    {lab:22s}" + cells)


def step_4():
    banner("STEP 4. WHAT EACH SENSOR BUYS PER UNIT OF PRECISION")
    print("""  The elasticity: how far sd(ln theta) falls when a sensor goes from its
  no-sensor value to a good one. This is the ranking that matters, because the
  sensors cost very different amounts and a large fall from a cheap sensor beats a
  small fall from four tracer studies. Monthly, the trial's own record.""")
    print(f"\n    {'cap':>5s} {'sensor, and the move priced':44s} "
          f"{'from':>8s} {'to':>8s} {'factor':>8s}")
    for cap in CAPS:
        th = TH_NOM * cap
        t = S2.times_for("monthly", "no baseline")
        base = cell(th, t)
        rows = [("candidate 2, ex vivo: none to 20%", dict(sd_th=0.20)),
                ("candidate 2, ex vivo: none to 10%", dict(sd_th=0.10)),
                ("candidate 8, k_e: 0.80 prior to 5.4%", dict(sd_ke=0.054)),
                ("candidate 9, turnover: none to 10%", dict(sd_b=0.10)),
                ("candidate 9, turnover: none to exact", dict(sd_b=0.0))]
        for lab, kw in rows:
            v = cell(th, t, **kw)
            print(f"    {100 * cap:4.0f}% {lab:44s} {base:8.3f} {v:8.3f} "
                  f"{base / v:8.2f}")
        print()


def step_4b():
    banner("4b. WHAT THE TWO-YEAR RECORD IS WORTH ON TOP OF THE EX VIVO ASSAY")
    print("""  Step 2 shows sd(ln theta) tracking candidate 2's own precision almost exactly,
  which says the record is contributing little. This measures how little. An ex
  vivo assay ALONE gives sd(ln theta) equal to its own precision, with no record
  at all. The record's marginal value is the gap between that and the joint
  figure. It is the question the whole project is about, asked the other way
  round.""")
    print(f"\n    {'cap':>5s} {'assay':>7s} {'assay alone':>12s}" +
          "".join(f"{s:>26s}" for s in SCHEDS))
    for cap in CAPS:
        th = TH_NOM * cap
        for a in (0.30, 0.20, 0.10):
            cells = []
            for sc in SCHEDS:
                v = cell(th, S2.times_for(sc, "no baseline"), sd_th=a)
                cells.append(f"{v:14.4f} ({100 * (1 - v / a):5.1f}%)")
            print(f"    {100 * cap:4.0f}% {100 * a:6.0f}% {a:12.4f}"
                  + "".join(cells))
    print("""
  The bracket is how much the two-year record improves on the assay by itself. A
  figure near zero means the record adds nothing to a patient who has had the
  assay, which is a statement about the record, not about the assay.""")


def step_4c():
    banner("4c. CANDIDATE 10's OWN FLOOR, WHICH IS NOT THE ASSAY")
    print("""  Timed urine measures `f_renal * k_e`, not `k_e`. Its precision on ln k_e is the
  assay's precision and the uncertainty in `f_renal` in quadrature, so `f_renal`
  sets a floor no urine collection improves on. R1 gives no `f_renal` for 3-HP;
  the rows below are what different states of knowledge about it would imply.""")
    print(f"\n    {'f_renal known to':28s} {'floor on ln k_e':>16s}" +
          "".join(f"{f'{100 * c:.0f}% cap, monthly':>22s}" for c in CAPS))
    for lab, sf in (("a factor of 3 (95% interval)", np.log(3.0) / 4.0),
                    ("+/- 20%", 0.20),
                    ("a factor of 2 (95% interval)", np.log(2.0) / 4.0),
                    ("+/- 5%", 0.05),
                    ("exactly, which R1 does not", 0.0)):
        floor = float(np.hypot(sf, 0.02))     # a 2% assay, the best plausible
        cells = "".join(
            f"{cell(TH_NOM * c, S2.times_for('monthly', 'no baseline'), sd_ke=floor):22.3f}"
            for c in CAPS)
        print(f"    {lab:28s} {floor:16.3f}" + cells)
    print("""
  Gate is 0.20. The comparison is against candidate 8, which measures `k_e`
  directly and has no `f_renal` in it at all.""")


COLLAPSE = S2.COLLAPSE


def mc_cell(theta, times, sd_th=np.inf, sd_ke=KE_PRIOR, sd_b=None, reps=96,
            seed0=0):
    """
    The fitted spread of ln theta_hat with the sensors present.

    The nuisance convention is `centre`, S2's per-patient measurement one: the
    patient's true `k_e` is `KE` and the prior CENTRE is drawn, because a sensor
    of stated precision is a measurement that lands somewhere. Candidate 2 enters
    the same way, as a drawn measurement of `ln theta` with the stated precision.
    Spread is reported robustly (1.4826 x MAD) with the collapse count, per the
    standing conventions.
    """
    L, _, sigma, s_b = S2.footing(theta)
    inc = inc_of(times)
    rng = np.random.default_rng(90210 + seed0)
    eff_b = s_b if sd_b is None else (
        0.0 if sd_b <= 0 else 1.0 / np.sqrt(1.0 / s_b ** 2
                                            + 1.0 / (sd_b * L) ** 2))
    errs, coll = [], 0
    for r in range(reps):
        Ls = S2.loads_for(times.size, seed0 * 1000 + r, L)
        b_t = float(rng.normal(0.0, s_b))
        m = S2.level(theta, inc, Ls + b_t)
        obs = ((theta + inc) * m / (KM + m) + KE * m + b_t
               + rng.normal(0.0, sigma, times.size))
        #  each sensor is a measurement that LANDS somewhere, so its centre is
        #  drawn about the truth at its own stated precision
        ke_c = KE * np.exp(rng.normal(0.0, sd_ke)) if sd_ke > 0 else KE
        th_c = (theta * np.exp(rng.normal(0.0, sd_th))
                if np.isfinite(sd_th) else None)
        b_c = (b_t + rng.normal(0.0, sd_b * L)
               if sd_b is not None and sd_b > 0 else None)
        best, bJ = None, np.inf
        for K0 in (30.0, 100.0, 300.0):
            def res(q):
                th, K, ke, b = np.exp(q[0]), np.exp(q[1]), np.exp(q[2]), q[3]
                out = [((th + inc) * m / (K + m) + ke * m + b - obs) / sigma,
                       np.array([(q[2] - np.log(ke_c)) / max(sd_ke, 1e-9)])]
                if eff_b > 0:
                    ref = 0.0 if b_c is None else b_c
                    out.append(np.array([(b - ref) / eff_b]))
                else:
                    out.append(np.array([(b - b_t) / (1e-6 * L)]))
                if th_c is not None:
                    out.append(np.array([(q[0] - np.log(th_c)) / sd_th]))
                return np.concatenate(out)
            q0 = np.array([np.log(theta if th_c is None else th_c),
                           np.log(K0), np.log(ke_c), 0.0])
            o = least_squares(res, q0, method="lm", xtol=1e-13, ftol=1e-13,
                              max_nfev=40000)
            J = float(np.sum(o.fun ** 2))
            if J < bJ:
                bJ, best = J, o
        e = float(best.x[0] - np.log(theta))
        if abs(e) > COLLAPSE:
            coll += 1
        errs.append(e)
    a = np.array(errs)
    rob = float(1.4826 * np.median(np.abs(a - np.median(a))))
    return rob, coll, reps


def step_5():
    banner("STEP 5. THE TWO VERDICTS CONFIRMED BY FITTING")
    print("""  The standing house rule. Two verdicts in steps 2 to 4 rest on Fisher figures,
  so both are refitted before being reported: that candidate 9 buys essentially
  nothing, and that the record adds almost nothing on top of candidate 2. Spread
  is robust (1.4826 x MAD), with the collapse count beside it, 96 records per
  cell, `centre` nuisance convention.""")
    print(f"\n    {'cap':>5s} {'configuration':34s} {'Fisher':>8s} "
          f"{'fitted':>8s} {'fitted/Fisher':>14s} {'collapses':>10s}")
    for cap in (0.10, 0.20):
        th = TH_NOM * cap
        t = S2.times_for("monthly", "no baseline")
        rows = [("record alone", {}),
                ("record + c9 exact turnover", dict(sd_b=0.0)),
                ("record + c2 at 20%", dict(sd_th=0.20)),
                ("record + c2 at 20% + c9 exact",
                 dict(sd_th=0.20, sd_b=0.0))]
        for lab, kw in rows:
            f = cell(th, t, **kw)
            rob, c, n = mc_cell(th, t, **kw)
            print(f"    {100 * cap:4.0f}% {lab:34s} {f:8.3f} {rob:8.3f} "
                  f"{rob / f:14.3f} {f'{c}/{n}':>10s}", flush=True)
        print()
    print("""  If the c9 rows match the rows above them under fitting as they do under Fisher,
  the verdict stands on both footings and candidate 9 is dropped on measurement,
  not on a bound.""")


def main():
    print("=" * 104)
    print("  PHASE 10 STEPS 2 TO 4: THE FIVE SURVIVING CANDIDATES, PRICED")
    print("=" * 104)
    step_2()
    step_3()
    step_4()
    step_4b()
    step_4c()
    step_5()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
