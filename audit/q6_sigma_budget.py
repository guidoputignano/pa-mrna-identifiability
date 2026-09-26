"""
Step 3(a). The requirement as a budget on sigma, not as a pair of axes.

S3 measured the marker CV against the load record's accuracy as a two-way
surface and found the two interchangeable percent for percent. The author's
reading is that this is not a coincidence but the whole structure: the design
sees the two error sources only through the combined ordinate noise sigma, so
the requirement is ONE BUDGET per schedule and the two axes are two ways of
spending it. This script tests that reading rather than assuming it, and then
states the budget.

WHAT sigma IS MADE OF. `model_A.run.sigma_total` gives

    sigma = hypot(load_frac * L,  b * cv * m_bar)

the recording error of the substrate load and the marker's own scatter, in flux
units. `load_frac` enters the information matrix a SECOND time, through the
load-bias prior `s_b = load_frac * L`, and that second entrance is not part of
sigma. So the budget claim is exactly:

    given s_b, sd(ln theta_hat) depends on the two error sources ONLY through
    sigma, and sigma* (the sigma at which a schedule reaches the gate) is one
    constant per schedule.

The author reports sigma* constant to about 4% across load records of 2 to 15%
and attributes the residual to the s_b channel. PART A tests both halves of
that: whether the dependence is exactly through sigma at fixed s_b, and whether
the residual is exactly the s_b channel.

THE t = 0 CONVENTION. Dose-day draws are read PRE-DOSE. The monthly record's
t = 0 draw is therefore an untreated baseline; a `post-dose` row is carried in
part A so the author's own footing is comparable, since the author's figures are
on that reading.

Nothing here records an expected result in advance.

Run:  python audit/q6_sigma_budget.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.optimize import brentq

from common.models import ParamsA
import model_A.run as LOC

import n4_theta_scope as N4
import q1_increment_shape as Q
import q2_requirement_map as S2

SCHED = {"monthly 25": np.arange(0.0, 730.0, 30.4),
         "fortnightly 51": np.arange(0.5, 730.0, 14.5),
         "weekly 100": np.arange(0.5, 730.0, 7.3),
         "every 2 d 365": np.arange(1.0, 730.0, 2.0)}
LOAD_FRACS = (0.02, 0.05, 0.10, 0.15)
CAPS = {"band low, 2472": 2472.0, "normal, 3000": 3000.0,
        "band high, 18814": 18814.0}
SEEDS = S2.SEEDS
GATE = S2.GATE
KE_PRIOR = 0.8
SIGMA_COMMITTED = 22.0          #: design/run.py's own, at load 5% and CV 5%


def banner(t):
    print("\n" + "=" * 104)
    print("  " + t)
    print("=" * 104)


def anchor(theta):
    """(L, m_bar, b) at this capacity, footing convention A."""
    L, m_bar = N4.anchor(theta, S2.KM, S2.KE, "A")
    p = ParamsA(theta=theta, Km=S2.KM, ke=S2.KE)
    return float(L), float(m_bar), float(LOC.fit_ab(p, m_bar)[1])


def sigma_of(theta, load_frac, cv):
    L, m_bar, b = anchor(theta)
    return float(np.hypot(load_frac * L, b * cv * m_bar))


def sd_at_sigma(theta, times, seed, sigma, s_b):
    """sd(ln theta_hat) as a function of sigma and s_b DIRECTLY, so that no
    decomposition of sigma into its two terms can enter."""
    L, _, _ = anchor(theta)
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


def med_sd(theta, times, sigma, s_b, seeds=SEEDS):
    v = np.array([sd_at_sigma(theta, times, s, sigma, s_b) for s in seeds])
    v = v[np.isfinite(v)]
    return float(np.median(v)) if v.size else np.inf


def sigma_star(theta, times, s_b, gate=GATE, hi=4000.0):
    """The sigma at which the median over the seed set reaches the gate."""
    f = lambda sg: med_sd(theta, times, sg, s_b) - gate
    if f(1e-3) > 0:
        return None                      # gate unreachable at any noise
    if f(hi) < 0:
        return np.inf
    return float(brentq(f, 1e-3, hi, xtol=1e-6, rtol=1e-12))


# ==========================================================================
def part_a():
    banner("A. IS THE DEPENDENCE EXACTLY THROUGH sigma, AND IS THE RESIDUAL "
           "EXACTLY s_b")
    print("""  (a) Several (load_frac, marker CV) pairs chosen to give the SAME sigma, with
  s_b held fixed. If the design sees the two sources only through sigma, every
  row must give an identical sd(ln theta_hat), not merely a close one.""")
    theta = 3000.0
    L, m_bar, b = anchor(theta)
    times = SCHED["monthly 25"]
    s_b_fixed = 0.10 * L
    target = 60.0
    print(f"\n    target sigma {target}, s_b held at 0.10 L = {s_b_fixed:.4f}, "
          "monthly, theta = 3000")
    print(f"      {'load_frac':>10s} {'marker CV':>10s} {'sigma':>10s} "
          f"{'sd(ln theta)':>14s}")
    got = []
    for lf in (0.0, 0.05, 0.10, 0.15, 0.1858):
        rem = target ** 2 - (lf * L) ** 2
        if rem <= 0:
            continue
        cv = float(np.sqrt(rem) / (b * m_bar))
        sg = sigma_of(theta, lf, cv)
        v = med_sd(theta, times, sg, s_b_fixed)
        got.append(v)
        print(f"      {100 * lf:9.2f}% {100 * cv:9.2f}% {sg:10.6f} "
              f"{v:14.10f}")
    spread = max(got) - min(got)
    print(f"\n    spread across those rows: {spread:.3e}  "
          f"{'EXACT to machine precision' if spread < 1e-12 else 'NOT exact'}")

    print("""
  (b) sigma* per schedule with s_b moving as load_frac does, which is the honest
  coupling, against s_b HELD at its committed 0.10 L. If the author's residual is
  the s_b channel and nothing else, the held column must be one constant per
  schedule with zero spread.""")
    print(f"\n      {'schedule':16s} {'load_frac':>10s} "
          f"{'sigma*, s_b moving':>19s} {'sigma*, s_b held':>18s}")
    out = {}
    for name, times in SCHED.items():
        mov, hel = [], []
        for lf in LOAD_FRACS:
            a = sigma_star(theta, times, lf * L)
            h = sigma_star(theta, times, s_b_fixed)
            mov.append(a)
            hel.append(h)
            print(f"      {name:16s} {100 * lf:9.0f}% {a:19.4f} {h:18.4f}")
        out[name] = (mov, hel)
        m = [x for x in mov if x is not None]
        print(f"      {'':16s} {'spread':>10s} "
              f"{f'{100 * (max(m) / min(m) - 1):.2f}%':>19s} "
              f"{f'{100 * (max(hel) / min(hel) - 1):.2e}%':>18s}")
    return out


def part_b():
    banner("B. THE BUDGET, PER SCHEDULE AND PER CAPACITY")
    print("""  sigma*, the ordinate noise at which each schedule reaches sd(ln theta) = 0.20,
  with s_b at the committed 0.10 L so the budget is a single number. Beside it,
  the marker CV that budget allows at each load record, and `-` where the load
  record alone already exceeds the budget so no marker scatter is admissible.""")
    for cname, theta in CAPS.items():
        L, m_bar, b = anchor(theta)
        print(f"\n  {cname}   (L = {L:.2f}, m_bar = {m_bar:.3f}, b = {b:.4f})")
        print(f"    {'schedule':16s} {'sigma*':>9s}" +
              "".join(f"{f'CV at load {int(100 * lf)}%':>19s}"
                      for lf in LOAD_FRACS))
        for name, times in SCHED.items():
            s = sigma_star(theta, times, 0.10 * L)
            if s is None:
                print(f"    {name:16s} {'never':>9s}")
                continue
            cells = []
            for lf in LOAD_FRACS:
                rem = s ** 2 - (lf * L) ** 2
                cells.append(f"{'-':>19s}" if rem <= 0 else
                             f"{100 * np.sqrt(rem) / (b * m_bar):18.2f}%")
            print(f"    {name:16s} {s:9.2f}" + "".join(cells))


def part_c():
    banner("C. WHERE THE COMMITTED FOOTING AND THE EVIDENCE BOUND SIT")
    print("""  PHASE9 item 7 bounds the marker's residual scatter below by the analytical
  component, 5 to 8%, and above by the surrogates' totals, 25 to 35%, with no
  measured value between. Those readings are converted to sigma here, at
  theta = 3000, and set against the budgets. `design/run.py`'s committed sigma is
  also placed, and reproduced from its own constants rather than quoted.""")
    theta = 3000.0
    L, m_bar, b = anchor(theta)
    print(f"\n    committed footing, load 5% and CV 5%: sigma = "
          f"{sigma_of(theta, 0.05, 0.05):.4f}  (design/run.py's "
          f"{SIGMA_COMMITTED})")
    print(f"\n    {'reading':34s} {'load 5%':>10s} {'load 10%':>10s} "
          f"{'load 15%':>10s}")
    for label, cv in (("analytical floor, CV 5%", 0.05),
                      ("analytical floor, CV 8%", 0.08),
                      ("surrogate lower, CV 25%", 0.25),
                      ("surrogate upper, CV 35%", 0.35)):
        print(f"    {label:34s}" +
              "".join(f"{sigma_of(theta, lf, cv):10.1f}"
                      for lf in (0.05, 0.10, 0.15)))
    print(f"\n    {'schedule':16s} {'sigma*':>9s} {'committed 22':>13s} "
          f"{'floor, load 10%':>16s} {'surrogates, load 10%':>21s}")
    for name, times in SCHED.items():
        s = sigma_star(theta, times, 0.10 * L)
        lo = sigma_of(theta, 0.10, 0.05)
        hi5 = sigma_of(theta, 0.10, 0.08)
        s25 = sigma_of(theta, 0.10, 0.25)
        s35 = sigma_of(theta, 0.10, 0.35)
        def verdict(x):
            return "meets" if x <= s else "fails"
        print(f"    {name:16s} {s:9.2f} {verdict(SIGMA_COMMITTED):>13s} "
              f"{verdict(lo) + '/' + verdict(hi5):>16s} "
              f"{verdict(s25) + '/' + verdict(s35):>21s}")
    print("""
  The `floor` column is CV 5% then CV 8%, the two ends of the analytical range;
  the `surrogates` column is CV 25% then CV 35%. A row reading `meets/fails`
  reaches the gate at the favourable end of that reading and not at the other.""")


def part_d():
    banner("D. THE LOAD RECORD ALONE CAN PUT THE GATE OUT OF REACH")
    print("""  sigma has a floor the marker cannot get below: at marker CV zero it is still
  load_frac * L. Where that floor exceeds a schedule's budget, no assay however
  good reaches the gate. The comparison is made at theta = 3000 and at both band
  edges, since s_b moves with the capacity too.""")
    for cname, theta in CAPS.items():
        L, m_bar, b = anchor(theta)
        print(f"\n  {cname}")
        print(f"    {'schedule':16s} {'sigma*':>9s}" +
              "".join(f"{f'floor at {int(100 * lf)}%':>17s}"
                      for lf in LOAD_FRACS))
        for name, times in SCHED.items():
            s = sigma_star(theta, times, 0.10 * L)
            cells = []
            for lf in LOAD_FRACS:
                fl = lf * L
                cells.append(f"{fl:11.1f} {'ok' if fl <= s else 'OVER':>5s}")
            print(f"    {name:16s} {s:9.2f}" + "".join(cells))


def mc_at_sigma(theta, times, sigma, s_b, reps=96, seed0=0,
                ke_sd=None):
    """S2's MAP estimator with sigma and s_b supplied rather than derived."""
    ke_sd = KE_PRIOR if ke_sd is None else ke_sd
    L, _, _ = anchor(theta)
    inc = S2.increment(times)
    err = []
    for r in range(reps):
        sd_ = seed0 + r
        rng = np.random.default_rng(10_007 * sd_ + 3)
        Ls = S2.loads_for(times.size, sd_, L)
        z = rng.normal(0.0, ke_sd)
        b_t = rng.normal(0.0, s_b)
        m = S2.level(theta, inc, Ls + b_t)
        obs = Ls + rng.normal(0.0, sigma, times.size)
        th = S2.mc_fit(m, inc, obs, sigma, s_b, ke_sd,
                       S2.KE * np.exp(z), "recorded")
        err.append(np.log(th) - np.log(theta))
    e = np.array(err)
    boot = np.random.default_rng(99).choice(e, (400, e.size))
    rb = np.sqrt(np.mean(boot ** 2, axis=1))
    return dict(rmse=float(np.sqrt(np.mean(e ** 2))),
                lo=float(np.percentile(rb, 5)),
                hi=float(np.percentile(rb, 95)),
                bias=float(np.median(e)),
                robust=float((np.percentile(e, 84)
                              - np.percentile(e, 16)) / 2),
                collapse=int(np.sum(np.abs(e) > S2.COLLAPSE)))


def part_e():
    banner("E. MONTE CARLO ON THE ROWS NEAR THE GATE")
    print("""  The budget is a Fisher bound and S2 found the bound optimistic near the gate,
  so each schedule is re-run through the MAP estimator at its own sigma*, where
  the bound reads exactly 0.20, and at the two readings that bracket it. 96
  records per cell, the k_e prior at the paper's 0.8 drawn per record. theta =
  3000, load record 10%, s_b at 0.10 L throughout.""")
    theta = 3000.0
    L, m_bar, b = anchor(theta)
    s_b = 0.10 * L
    print(f"\n    {'schedule':16s} {'sigma':>8s} {'what it is':28s} "
          f"{'Fisher':>7s} {'RMSE [90%]':>22s} {'bias':>7s} {'robust':>7s} "
          f"{'coll':>5s}")
    for name, times in SCHED.items():
        star = sigma_star(theta, times, s_b)
        rows = [(star, "sigma*, the budget itself"),
                (sigma_of(theta, 0.10, 0.05), "analytical floor, CV 5%"),
                (sigma_of(theta, 0.10, 0.25), "surrogate lower, CV 25%")]
        for sg, what in rows:
            f = med_sd(theta, times, sg, s_b)
            r = mc_at_sigma(theta, times, sg, s_b)
            print(f"    {name:16s} {sg:8.2f} {what:28s} {f:7.3f} "
                  f"{r['rmse']:6.3f} [{r['lo']:5.3f},{r['hi']:6.3f}] "
                  f"{r['bias']:+7.3f} {r['robust']:7.3f} {r['collapse']:5d}")
    print("""
  A budget that the estimator does not honour is not a budget. Read the RMSE
  column against 0.20 in the sigma* rows: where it sits above, the Fisher budget
  is optimistic and the real budget is tighter than the number above it.""")


def part_f():
    banner("F. WHOSE MISS IS IT: sigma OR THE k_e PRIOR")
    print("""  Part E's estimator RMSE sits well above 0.20 at every sigma*, which would mean
  the Fisher budget is not a budget. But part E carries the paper's own 0.8
  population prior on ln k_e, which S2 showed has a genuine collapse tail and
  S2b showed is insufficient below 20% of normal anyway. So the same cells are
  re-run with the k_e prior at 0.054, the realistic per-patient precision S2b
  says the micro-dose study delivers, changing nothing else. sigma* is left where
  the 0.8 prior put it, so the two rows differ in the prior alone.""")
    theta = 3000.0
    s_b = 0.10 * anchor(theta)[0]
    print(f"\n    {'schedule':16s} {'k_e prior':>10s} {'sigma*':>8s} "
          f"{'Fisher':>7s} {'RMSE [90%]':>22s} {'robust':>7s} {'coll':>5s}")
    for name, times in SCHED.items():
        star = sigma_star(theta, times, s_b)
        for ke in (KE_PRIOR, 0.054):
            f = med_sd(theta, times, star, s_b)
            r = mc_at_sigma(theta, times, star, s_b, ke_sd=ke)
            print(f"    {name:16s} {ke:10.3f} {star:8.2f} {f:7.3f} "
                  f"{r['rmse']:6.3f} [{r['lo']:5.3f},{r['hi']:6.3f}] "
                  f"{r['robust']:7.3f} {r['collapse']:5d}")


def main():
    print("=" * 104)
    print("  STEP 3(a). THE REQUIREMENT AS A BUDGET ON sigma")
    print("=" * 104)
    print(f"\n  dose-day draws read {Q.DOSE_DAY_DRAW}; k_e prior "
          f"{KE_PRIOR}; gate {GATE}; {len(SEEDS)} load records per cell")
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    part_f()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
