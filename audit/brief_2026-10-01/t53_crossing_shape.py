"""
Round 19: does the 21-day knife-edge survive Model C's unmeasured shape constants?

The 3.5452 d crossing (R-N) is where the 21-day trough increment equals Model C's
intrinsic critical capacity, theta_crit_int = 150.775 at the repo's constants. That
threshold depends on three constants of the carnitine channel that nobody has
measured: K_c (models.py: "ARBITRARY. There is no measurement behind it"), M_i and
rho_0. The repository already sweeps them (model_C/run.py: shape_sensitivity,
collapse_check), re-deriving R0 each time so the alternative route still peaks at
500. It never carried that sweep into the dosing crossing. This script does.

Same-rate trough convention as t44/t48/t50: Delta_trough(tau) = D14 (tau/14) x/(1-x),
x = exp(-lambda tau), lambda = ln2 / half-life.

Two readings:
  * the repo's one-at-a-time sweep, each constant over the range it lists;
  * a joint sweep over the same ranges collapse_check uses (log-uniform K_c and M_i,
    uniform rho_0). The share of sampled sets in which a half-life clears the line
    is a SENSITIVITY summary over the admissible box, not a probability: nothing
    measured puts a distribution on these constants.
"""
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
import model_C.run as MC  # noqa: E402

D14, TAU0 = 6000.0, 14.0
MEASURED, CENTRAL = 3.6, 4.3


def trough(tau, thalf):
    lam = np.log(2.0) / thalf
    x = np.exp(-lam * tau)
    return D14 * (tau / TAU0) * x / (1.0 - x)


def crossing(thc, tau=21):
    """Half-life at which Delta_trough(tau) = theta_crit."""
    return brentq(lambda t: trough(tau, t) - thc, 0.3, 30.0)


def main():
    base = MC.theta_crit_intrinsic(MC.PC, MC.S_BASE)
    print("=" * 92)
    print("  ROUND 19. THE 21-DAY CROSSING UNDER MODEL C'S UNMEASURED SHAPE CONSTANTS")
    print("=" * 92)
    print(f"""
  repo constants (K_c=100, M_i=300, rho_0=0.5): theta_crit_int = {base:.3f}
    crossing half-life: 14-day {crossing(base, 14):.3f} d   15-day {crossing(base, 15):.3f} d   21-day {crossing(base):.4f} d
""")
    print("  ONE AT A TIME (model_C/run.py shape_sensitivity, R0 re-derived each time):")
    print(f"    {'constant':>8} {'value':>7} {'theta_crit':>10} {'21-d crossing':>13} {'3.6 d clears?':>14} {'4.3 d clears?':>14}")
    for name, rec in MC.shape_sensitivity().items():
        for v, tc_cal, _tc_raw, _mg in rec:
            if not np.isfinite(tc_cal):
                continue
            c = crossing(tc_cal)
            print(f"    {name:>8} {v:>7} {tc_cal:>10.1f} {c:>12.3f}d {('yes' if MEASURED > c else 'NO'):>14}"
                  f" {('yes' if CENTRAL > c else 'NO'):>14}")

    rng = np.random.default_rng(0)
    tcs = []
    for _ in range(300):
        q = replace(MC.PC, Kc=float(10 ** rng.uniform(1.4, 2.9)),
                    Mi=float(10 ** rng.uniform(2.0, 3.1)),
                    rho0=float(rng.uniform(0.35, 1.0)))
        q = replace(q, R0=MC.R0_for_peak(q))
        t = MC.theta_crit_intrinsic(q, MC.S_BASE)
        if np.isfinite(t):
            tcs.append(t)
    tcs = np.array(tcs)
    c14 = np.array([crossing(t, 14) for t in tcs])
    c15 = np.array([crossing(t, 15) for t in tcs])
    c21 = np.array([crossing(t) for t in tcs])
    print(f"""
  JOINT SWEEP over collapse_check's ranges (K_c 25-794, M_i 100-1259, rho_0 0.35-1.0;
  {len(tcs)} sampled sets with a finite threshold, seed 0):
    theta_crit_int           {tcs.min():6.1f} .. {tcs.max():6.1f}   (median {np.median(tcs):.1f})
    14-day crossing          {c14.min():6.3f} .. {c14.max():6.3f} d
    15-day crossing          {c15.min():6.3f} .. {c15.max():6.3f} d
    21-day crossing          {c21.min():6.3f} .. {c21.max():6.3f} d (median {np.median(c21):.3f})
    sampled sets in which the half-life clears the 21-day line:
        measured 3.6 d       {np.mean(MEASURED > c21):.0%}
        liver estimate 4.3 d {np.mean(CENTRAL > c21):.0%}
    sampled sets in which 3.6 d clears the 14- and 15-day lines: {np.mean(MEASURED > c14):.0%} / {np.mean(MEASURED > c15):.0%}
""")
    # the eight corners of the same box: the extremes, which 300 random draws only approach
    import itertools
    print("  THE EIGHT CORNERS OF THE BOX (R0 re-derived at each):")
    print(f"    {'K_c':>6} {'M_i':>6} {'rho_0':>5} {'theta_crit':>10} {'14-d':>6} {'15-d':>6} {'21-d':>6}  3.6 d clears 21 d?")
    corners = []
    for kc, mi, r0 in itertools.product((10 ** 1.4, 10 ** 2.9), (10 ** 2.0, 10 ** 3.1), (0.35, 1.0)):
        q = replace(MC.PC, Kc=float(kc), Mi=float(mi), rho0=float(r0))
        q = replace(q, R0=MC.R0_for_peak(q))
        t = MC.theta_crit_intrinsic(q, MC.S_BASE)
        row = (kc, mi, r0, t, crossing(t, 14), crossing(t, 15), crossing(t))
        corners.append(row)
        print(f"    {kc:>6.1f} {mi:>6.0f} {r0:>5.2f} {t:>10.1f} {row[4]:>6.3f} {row[5]:>6.3f} {row[6]:>6.3f}  "
              f"{'yes' if MEASURED > row[6] else 'NO'}")
    cr = np.array(corners)
    print(f"""    box: theta_crit {cr[:,3].min():.1f}..{cr[:,3].max():.1f}; worst 14-day {cr[:,4].max():.3f} d, 15-day {cr[:,5].max():.3f} d;
         21-day {cr[:,6].min():.3f}..{cr[:,6].max():.3f} d; 3.6 d clears at {int((MEASURED > cr[:,6]).sum())} of 8 corners, 4.3 d at {int((CENTRAL > cr[:,6]).sum())} of 8.
""")
    lo, hi = MC.collapse_check()
    print(f"  For comparison, the repo's own collapse_check (12 draws) spans theta_crit {lo:.1f}..{hi:.1f},"
          f" a factor {hi/lo:.2f};\n  the 300 draws above span a factor {tcs.max()/tcs.min():.2f}, the box corners"
          f" {cr[:,3].max()/cr[:,3].min():.2f}. Twelve draws understate the joint spread.")
    from scipy.stats import norm
    print("\n  R-P's 'about 65-70 % likely to hold' recomputed at the crossing (centre 4.3 d; R-P's wide band")
    print("  sd(ln) 0.53 and tight band sd(ln) 0.35):")
    for lab, c in (("repo constants", crossing(base)), ("largest sampled crossing", c21.max()),
                   ("worst box corner", cr[:,6].max())):
        pw = 1 - norm.cdf(np.log(c / CENTRAL) / 0.53)
        pt = 1 - norm.cdf(np.log(c / CENTRAL) / 0.35)
        print(f"    {lab:<26} crossing {c:.3f} d -> P(liver half-life above it) = {pw:.0%} wide, {pt:.0%} tight")
    # The paper's own ranges (results/tables/tab_C_shape.tex: K_c 25-800, M_i 150-1200,
    # rho_0 0.35-1), which its remark on the shape constants sweeps jointly and reports as
    # "inside a factor of 1.4". Also the paper's largest admissible increment at 5 % of
    # normal (11_scope.tex: D = 2672) beside the booked D = 6000.
    D_ADM = 2672.0
    def cross_D(t, tau, D):
        return brentq(lambda h: trough(tau, h) * D / D14 - t, 0.3, 40.0)
    print("\n  THE PAPER'S OWN RANGES (tab_C_shape: K_c 25-800, M_i 150-1200, rho_0 0.35-1), corners,")
    print("  R0 re-derived; crossings at the booked D = 6000 and at the largest admissible D = 2672:")
    print(f"    {'K_c':>5} {'M_i':>5} {'rho_0':>5} {'theta_crit':>10}   D=6000: 14-d  15-d  21-d    D=2672: 14-d  15-d  21-d")
    pc = []
    for kc, mi, r0 in itertools.product((25.0, 800.0), (150.0, 1200.0), (0.35, 1.0)):
        q = replace(MC.PC, Kc=kc, Mi=mi, rho0=r0)
        q = replace(q, R0=MC.R0_for_peak(q))
        t = MC.theta_crit_intrinsic(q, MC.S_BASE)
        row = [t] + [cross_D(t, tau, D) for D in (D14, D_ADM) for tau in (14, 15, 21)]
        pc.append(row)
        print(f"    {kc:>5.0f} {mi:>5.0f} {r0:>5.2f} {t:>10.1f}           {row[1]:.2f}  {row[2]:.2f}  {row[3]:.2f}"
              f"            {row[4]:.2f}  {row[5]:.2f}  {row[6]:.2f}")
    pc = np.array(pc)
    nom = [cross_D(base, tau, D) for D in (D14, D_ADM) for tau in (14, 15, 21)]
    print(f"""    nominal constants                    {nom[0]:.2f}  {nom[1]:.2f}  {nom[2]:.2f}            {nom[3]:.2f}  {nom[4]:.2f}  {nom[5]:.2f}
    corners: theta_crit {pc[:,0].min():.1f}..{pc[:,0].max():.1f}, a factor {pc[:,0].max()/pc[:,0].min():.2f} (the paper's remark: 'inside a factor of 1.4')
    ranges over the corners, D=6000: 14-d {pc[:,1].min():.2f}-{pc[:,1].max():.2f}, 15-d {pc[:,2].min():.2f}-{pc[:,2].max():.2f}, 21-d {pc[:,3].min():.2f}-{pc[:,3].max():.2f}
                             D=2672: 14-d {pc[:,4].min():.2f}-{pc[:,4].max():.2f}, 15-d {pc[:,5].min():.2f}-{pc[:,5].max():.2f}, 21-d {pc[:,6].min():.2f}-{pc[:,6].max():.2f}
    3.6 d clears the 21-day line (D=6000) at {int((MEASURED > pc[:,3]).sum())} of 8 corners; 4.3 d at {int((CENTRAL > pc[:,3]).sum())} of 8""")
    for lab, c in (("21-d, nominal, D=6000", nom[2]), ("21-d, worst corner, D=6000", pc[:,3].max()),
                   ("21-d, nominal, D=2672", nom[5]), ("14-d, worst corner, D=2672", pc[:,4].max()),
                   ("15-d, worst corner, D=2672", pc[:,5].max())):
        pw = 1 - norm.cdf(np.log(c / CENTRAL) / 0.53); pt = 1 - norm.cdf(np.log(c / CENTRAL) / 0.35)
        print(f"    P(liver half-life above the {lab} crossing {c:.2f} d) = {pw:.0%} wide, {pt:.0%} tight")
    print(f"""
  READING. The 14- and 15-day conclusions do not depend on Model C's unmeasured
  shape: even at the worst corner of the box their crossings are about {cr[:,4].max():.2f} and {cr[:,5].max():.2f} d,
  under the measured 3.6 d. (They still need the liver half-life above those values.)
  The 21-day conclusion at 3.6 d does depend on the shape. The crossing moves from about
  {cr[:,6].min():.2f} to {cr[:,6].max():.2f} d across the box ({c21.min():.2f}-{c21.max():.2f} d over the 300 sampled sets), which is
  wider than the 0.055 d (1.5 %) by which 3.6 d clears the line at K_c = 100. One constant
  at a time, 3.6 d fails at K_c <= 50 or at M_i = 150; jointly it clears in 80 % of the
  sampled sets and at 6 of the 8 corners. The liver estimate of 4.3 d clears the line
  everywhere in the box. So the knife-edge has two blades: the true liver half-life, which
  the data requests address, and the carnitine channel's shape, which no data in the
  repository addresses.
""")


if __name__ == "__main__":
    main()
