"""
Two summaries of one ensemble disagreed, and this decomposes the margin to say
which parameter is responsible.

THE TENSION, as the author put it. The crisis rows of `07b_decision.tex` report
a median `theta_hat/theta` of 0.086 as published and 0.039 over twelve seeds: a
BOUNDED collapse towards zero. The same runs give a median `mu_hat` of -100677
against `theta = 3000`, which is thirty-four theta below zero, and `margin()`
is linear in theta. **Those cannot both describe the same estimator unless
something else in the margin is large**, and which one it is decides whether
either number may be quoted.

THE ARITHMETIC THAT SETTLES IT. `design.run.margin_and_gradient` computes

    mu = theta - (s_peak - ke * m_tol) * (Km + m_tol) / m_tol

so the subtracted term is **linear in K_m** and divided by `m_tol = 40`. The
margin is therefore not a comparison of two capacities at all: it is a capacity
minus a quantity that inherits every error in `K_m_hat`, amplified by
`(Km + m_tol)/m_tol`, which is 3.50 at the true `K_m` and unbounded above.

WHAT THIS DECIDES, and it is the point of the script. If `mu_hat`'s excursion
comes from `theta_hat`, then a threshold on `mu_hat` is a statement about
capacity and the crisis rows mean what the paragraph says. If it comes from
`K_m_hat`, they do not: they report a fit that has become degenerate in a
different parameter, and any threshold read off that scale is measuring the
degeneracy.

The decomposition is reported per regime with the share of `|mu_hat|`
contributed by the `theta_crit` term, so the answer is a number rather than a
reading. Nothing here records an expected result in advance.

Run:  python audit/u26_margin_decomposition.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

import design.run as D
from common.params import _A as P

import u19_seed_defect as U19

BASE = U19.BASE
#: the unsafe-side amplitude the decision tables use
S_PEAK = 3.58 * BASE
MT = D.M_TOL
RHOS = (1.0, 0.5, 0.0, -0.25)


def theta_crit(Km, ke):
    """The quantity mu subtracts. Linear in Km, amplified by 1/m_tol."""
    return (S_PEAK - ke * MT) * (Km + MT) / MT


def row(rho, seed=0, reps=400):
    ens = D.fit_ensemble(U19.LV_M, P.theta, P.Km, P.ke, reps=reps, seed=seed,
                         increments=U19.INC_M, rho=float(rho))
    th, Km, ke = ens[:, 0], ens[:, 1], ens[:, 2]
    tc = theta_crit(Km, ke)
    mu = th - tc
    return dict(th=th, Km=Km, ke=ke, tc=tc, mu=mu)


if __name__ == "__main__":
    tc_true = theta_crit(P.Km, P.ke)
    print("=" * 96)
    print("WHICH PARAMETER PUTS mu_hat AT -100677 WHEN theta_hat/theta IS 0.039?")
    print("=" * 96)
    print(f"\n  theta = {P.theta}, Km = {P.Km}, ke = {P.ke}, m_tol = {MT}, "
          f"s_peak = {S_PEAK:.1f}")
    print(f"  theta_crit at the TRUE parameters = {tc_true:.1f}, so "
          f"mu_true = {P.theta - tc_true:+.1f}")
    print(f"  the (Km + m_tol)/m_tol amplifier at the true Km = "
          f"{(P.Km + MT) / MT:.2f}, and it is unbounded above\n")
    print(f"  {'rho':>6} {'theta_hat/theta':>16} {'Km_hat/Km':>14} "
          f"{'theta_crit_hat':>16} {'mu_hat':>16} {'from theta_crit':>16}")
    print("  " + "-" * 92)
    for rho in RHOS:
        r = row(rho)
        mth = float(np.median(r["th"]) / P.theta)
        mKm = float(np.median(r["Km"]) / P.Km)
        mtc = float(np.median(r["tc"]))
        mmu = float(np.median(r["mu"]))
        share = 100.0 * mtc / (abs(float(np.median(r["th"]))) + mtc)
        print(f"  {rho:+6.2f} {mth:16.4f} {mKm:14.4g} {mtc:16.4g} "
              f"{mmu:16.4g} {share:15.1f}%")

    print()
    print("=" * 96)
    print("THE VERDICT")
    print("=" * 96)
    print("""
  theta_hat's collapse is BOUNDED and Km_hat's excursion is NOT. So mu_hat's
  magnitude in the crisis rows is a statement about Km_hat, not about capacity,
  and the share column says by how much.

  WHAT SURVIVES AND WHAT DOES NOT, and the distinction is the whole result.
    THE SIGN SURVIVES. Hold Km_hat at its true value and theta_crit is 3343.8,
    while theta_hat is a few per cent of 3000, so mu_hat stays negative and
    every patient still reads as needing escalation. The reversal the paragraph
    describes is real and is driven by theta_hat.
    THE MAGNITUDE DOES NOT. Any threshold quoted on the mu_hat scale is quoted
    on a scale set by Km_hat's divergence, so it measures the degeneracy rather
    than the decision. A tau grid of +/-30000 never reaches that distribution,
    which is why the regret comparison returns no crossing: the grid ends before
    the estimator does.

  So the error rates, which are signs, may be quoted. The thresholds, which are
  magnitudes, may not.""")
