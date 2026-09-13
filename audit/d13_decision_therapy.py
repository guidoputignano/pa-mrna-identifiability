"""
The decision characteristic on the therapy-channel design, both regimes.

Everything earlier in 1.9 was computed on the fixed-point design (nine levels
from m_bar to 5 m_bar, no increment). This asks the same question of the design
the specification now recommends: monthly sampling across dose cycles, with the
increment carried and intake varying day to day.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
import design.run as D

P = D.P
base = D.clearance_A(D.M_BAR, P)


def main():
    print(__doc__)
    old = np.geomspace(D.M_BAR, 5 * D.M_BAR, 9)
    t_m = np.arange(0.0, 730.0, 30.4)
    lv_m, inc_m = D.design_record(t_m)
    cases = [("fixed point, n=9 (the old design)", old, None),
             ("monthly across dose cycles, n=25", lv_m, inc_m)]

    print("=" * 92)
    print("SIGN ERROR either side of the boundary, 12 seeds x 400 replicates")
    print("=" * 92)
    print(f"{'design':>36} {'kappa':>6} | {'treat unneeded':>22} {'fail to treat':>22}")
    for nm, lv, inc in cases:
        for kcv in ((0.0,) if inc is None else (0.0, 0.117)):
            a, b = [], []
            for sd_ in range(12):
                a.append(D.sign_error(lv, P.theta, P.Km, P.ke, reps=400, seed=sd_,
                                      s_peak=2.98 * base, increments=inc,
                                      kappa_cv=kcv)["p"])
                b.append(D.sign_error(lv, P.theta, P.Km, P.ke, reps=400, seed=sd_,
                                      s_peak=3.58 * base, increments=inc,
                                      kappa_cv=kcv)["p"])
            a, b = np.array(a), np.array(b)
            lab = nm if kcv == 0.0 else ""
            print(f"{lab:>36} {kcv:6.3f} | {a.mean():8.1%} +/- {a.std():4.1%} "
                  f"[{a.min():5.1%},{a.max():5.1%}] {b.mean():8.1%} +/- {b.std():4.1%} "
                  f"[{b.min():5.1%},{b.max():5.1%}]")

    print()
    print("=" * 92)
    print("THE COHORT, on each design")
    print("=" * 92)
    for nm, lv, inc in cases:
        ps = [D.sign_error(lv, P.theta, P.Km, P.ke, reps=400, seed=sd_,
                           s_peak=3.6 * base, increments=inc)["p"] for sd_ in range(8)]
        p = float(np.mean(ps))
        N = D.cohort_size(p)
        print(f"  {nm:>36}: p(err) near the boundary = {p:5.1%} +/- {np.std(ps):.1%}"
              f"   N = {('%d' % N) if np.isfinite(N) else 'infinite'}")

    print()
    print("=" * 92)
    print("CRISIS REGIME on the therapy-channel design")
    print("=" * 92)
    print(f"{'recording fraction rho':>24} | {'treat unneeded':>16} {'fail to treat':>15}"
          f" {'median th/th':>14}")
    for rho in (1.0, 0.5, 0.0, -0.25):
        ens = D.fit_ensemble(lv_m, P.theta, P.Km, P.ke, reps=400, seed=11,
                             increments=inc_m, rho=rho)
        e_s = D.sign_error(lv_m, P.theta, P.Km, P.ke, reps=400, seed=11,
                           s_peak=2.98 * base, increments=inc_m, rho=rho)["p"]
        e_u = D.sign_error(lv_m, P.theta, P.Km, P.ke, reps=400, seed=11,
                           s_peak=3.58 * base, increments=inc_m, rho=rho)["p"]
        print(f"{rho:24.2f} | {e_s:15.1%} {e_u:14.1%} "
              f"{np.median(ens[:,0])/P.theta:14.3f}")


if __name__ == "__main__":
    main()
