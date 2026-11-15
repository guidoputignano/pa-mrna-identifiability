"""
Round 20 (close): the within-patient contrast, simulated rather than linearised.

Result 3 (t22, t24, t30, t42) gives the contrast between two windows of 104 weekly draws
at 5 % of normal as a Fisher value, sd(ln th2 - ln th1) = 0.133. That number has never been
checked by fitting simulated record pairs. Here it is, on t24's exact footing:

  * two windows, 104 draws each over 730 d, on therapy (D = 6000), 10 % load variation,
    per-draw noise (5 % assay CV, 5 % load-record scatter: core.sigmas 'perdraw');
  * K_m, the load bias b and k_e shared by the two windows (hard-tied, as in t24), theta
    free in each; t24 adds the k_e prior (sd 0.8 on ln k_e) and the b prior (5 % of L) to
    BOTH tied copies, so the shared parameter carries each prior twice: sd 0.8/sqrt 2 and
    0.05 L/sqrt 2 here;
  * joint least squares with theta >= 1e-3, multi-start;
  * window 2 of a changed patient: the same load series generator, theta2 = ratio x theta1,
    and the marker recomputed at theta2 (core.level), so the patient's load scale is unchanged.

Questions: (1) does the fitted log ratio have the Fisher sd, and (2) does the difference
theta2_hat - theta1_hat? (3) What does each detect, and how often does each raise a false
alarm?

Run: python -W error t56_contrast_mc.py   (writes out/t56_contrast_mc.txt via stdout)
"""
import numpy as np
from scipy.optimize import least_squares

import core as C
import t24_contrast_stress as T24

TH = C.TH_5
N, WIN = 104, 730.0
LB = np.log(1e-3)
KEP = C.KE_PRIOR / np.sqrt(2.0)
BFRAC = 0.05 / np.sqrt(2.0)
REPS = 400


def fit_pair(rs, ys, sgs, L, starts):
    sb = BFRAC * L

    def res(x):
        lt1, lt2, lk, b, lke = x
        out = []
        for r, y, sg, lt in zip(rs, ys, sgs, (lt1, lt2)):
            m, inc = r["m"], r["inc"]
            pred = (np.exp(lt) + inc) * m / (np.exp(lk) + m) + np.exp(lke) * m + b
            out.append((y - pred) / sg)
        out.append([b / sb, (lke - np.log(C.KE)) / KEP])
        return np.concatenate(out)

    best = None
    for x0 in starts:
        s = least_squares(res, x0, method="trf", bounds=([LB, LB, -np.inf, -np.inf, -np.inf], [np.inf] * 5),
                          xtol=1e-12, ftol=1e-12, max_nfev=4000)
        if best is None or s.cost < best.cost:
            best = s
    return best.x


def pair(k, ratio, rng):
    times = np.linspace(0.0, WIN, N)
    r1 = C.record(TH, seed=2 * k, times=times)
    r2 = dict(C.record(TH, seed=2 * k + 1, times=times))
    if ratio != 1.0:
        r2["theta"] = TH * ratio
        r2["m"] = C.level(r2["theta"], C.KM0, C.KE, r2["loads"], r2["inc"])
    rs = [r1, r2]
    sgs = [C.sigmas(r, "perdraw") for r in rs]
    ys = [(r["theta"] + r["inc"]) * r["m"] / (r["Km"] + r["m"]) + r["ke"] * r["m"] + rng.normal(0.0, sg)
          for r, sg in zip(rs, sgs)]
    tr = np.array([np.log(TH), np.log(TH * ratio), np.log(C.KM0), 0.0, np.log(C.KE)])
    starts = [tr] + [tr + np.array([d, d, 0.0, 0.0, -0.15 * d]) for d in (-3.0, -1.5, -0.5, 0.5, 1.0)]
    x = fit_pair(rs, ys, sgs, r1["L"], starts)
    return x


def main():
    print("=" * 92)
    print("  ROUND 20 (close). THE WITHIN-PATIENT CONTRAST, SIMULATED (5 % of normal, 2 x 104 weekly draws)")
    print("=" * 92)
    fims = [T24.contrast(N, WIN, 2 * k, 2 * k + 1)[0] for k in range(REPS)]
    sd_ln = float(np.median(fims))
    sd_d = sd_ln * TH
    print(f"  Fisher (t24) on the same {REPS} no-change pairs: sd(ln th2 - ln th1) median {sd_ln:.3f};"
          f" linearised at theta = {TH:.0f}, that is sd(th2 - th1) = {sd_d:.1f} uM/d")
    print(f"  tests: log ratio rejects when |ln th2 - ln th1| > 1.96 x {sd_ln:.3f};"
          f" difference rejects when |th2 - th1| > 1.96 x {sd_d:.1f} uM/d")
    print(f"\n  {'true ratio':>10} {'at floor':>9} {'ln-ratio: median':>17} {'robust sd':>10} {'2.5-97.5 %':>17}"
          f" {'outside +-1.96 sd':>18} | {'diff: mean':>10} {'sd':>6} | {'reject: ln':>10} {'diff':>6}")
    for ratio in (1.0, 1.2, 0.8, 1.45, 0.70):
        rng = np.random.default_rng(int(1000 * ratio) + 7)
        xs = np.array([pair(k, ratio, rng) for k in range(REPS)])
        lr = xs[:, 1] - xs[:, 0]
        df = np.exp(xs[:, 1]) - np.exp(xs[:, 0])
        floor = np.mean(np.min(xs[:, :2], axis=1) <= LB + 1e-6)
        q = np.percentile(lr, [2.5, 25, 50, 75, 97.5])
        out = np.mean(np.abs(lr - np.log(ratio)) > 1.96 * sd_ln)
        rej_ln = np.mean(np.abs(lr) > 1.96 * sd_ln)
        rej_d = np.mean(np.abs(df) > 1.96 * sd_d)
        print(f"  {ratio:>10.2f} {floor:>9.2f} {q[2]:>+8.3f} (true {np.log(ratio):+.3f}) {(q[3]-q[1])/1.349:>8.3f}"
              f" [{q[0]:+6.2f}, {q[4]:+6.2f}] {out:>18.2f} | {df.mean():>+10.1f} {df.std():>6.1f} |"
              f" {rej_ln:>10.2f} {rej_d:>6.2f}")
    print(f"""
  Reading. At no change the middle of the fitted log ratio matches the Fisher value, but its
  tails do not: whenever one window's theta goes to the floor (the record cannot rule out zero
  capacity, t1 1d), the other follows it towards zero and their ratio is arbitrary. The
  difference theta2_hat - theta1_hat has the Fisher sd and is nearly unbiased with or without a
  change, because k_e and the load bias shift both windows' estimates by the same AMOUNT
  (d theta = -d k_e (K_m + m) at a fixed level), not by the same factor. So what the record
  measures is the change in uM/d; the 0.133 is that change's sd divided by the true 150. A test
  on the difference is calibrated; one on the log ratio raises false alarms. With sd
  {sd_d:.0f} uM/d the change detected with 80 % power (two-sided 5 %) is 2.80 x {sd_d:.0f} =
  {2.80 * sd_d:.0f} uM/d either way, {100 * 2.80 * sd_d / TH:.0f} % of a 5 % patient's capacity.""")


if __name__ == "__main__":
    main()
