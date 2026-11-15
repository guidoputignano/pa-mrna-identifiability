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
  {len(tcs)} sets with a finite threshold, seed 0):
    theta_crit_int           {tcs.min():6.1f} .. {tcs.max():6.1f}   (median {np.median(tcs):.1f})
    14-day crossing          {c14.min():6.3f} .. {c14.max():6.3f} d
    15-day crossing          {c15.min():6.3f} .. {c15.max():6.3f} d
    21-day crossing          {c21.min():6.3f} .. {c21.max():6.3f} d (median {np.median(c21):.3f})
    sets in which the half-life clears the 21-day line:
        measured 3.6 d       {np.mean(MEASURED > c21):.0%}
        liver estimate 4.3 d {np.mean(CENTRAL > c21):.0%}
    sets in which 3.6 d clears the 14- and 15-day lines: {np.mean(MEASURED > c14):.0%} / {np.mean(MEASURED > c15):.0%}

  READING. The 14- and 15-day conclusions do not depend on Model C's unmeasured
  shape: their crossings stay below 3.0 and 3.2 d, far under the measured 3.6 d.
  The 21-day conclusion at 3.6 d does. The crossing itself moves from about 3.3 to
  4.0 d across the admissible constants, which is wider than the 0.055 d (1.5 %)
  by which 3.6 d clears the line at K_c = 100. At K_c <= 50, or M_i = 150, the
  measured 3.6 d falls on the needs-theta side. The liver estimate of 4.3 d clears
  the line for every sampled set. So the knife-edge has two blades: the true liver
  half-life, and the carnitine channel's shape constants. The first is a data
  request; the second is a model assumption no data in the repository can fix.
""")


if __name__ == "__main__":
    main()
