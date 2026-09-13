"""
C4 -- regenerate the A1 and A3 numbers PHASE2A.md states, from committed code.

A1: the ablation table (the 1.8 filter's assumptions peeled off one at a time),
    the placement ratio on the footing A1 recommends, and the excursion curve.
A3: the excitation available at the fold ranges the literature supports against
    the ones the repo hardcodes, and what the decision curve returns there.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
import design.run as D, model_A.run as LOC
from common.excitation import excitation, excitation_max

P = D.P
Lb = D.load(0.0)
SIG = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)


def sd_theta(levels, inc=None, ke_free=False, ke_prior=0.8, flat_load=False,
             sigma=SIG):
    """sd(ln theta) with each of the 1.8/1.9 assumptions switchable."""
    pri = [0.8 ** 2, 0.8 ** 2,
           np.inf if flat_load else (D.LOAD_FRAC * Lb) ** 2]
    if ke_free:
        pri.append(ke_prior ** 2)
    Finv = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri])
    lv = np.asarray(levels, dtype=float)
    dd = (np.zeros_like(lv) if inc is None
          else np.broadcast_to(np.asarray(inc, dtype=float), lv.shape))
    for m, d in zip(lv, dd):
        q = [P.theta * m / (P.Km + m), -(P.theta + d) * P.Km * m / (P.Km + m) ** 2, 1.0]
        if ke_free:
            q.append(P.ke * m)
        q = np.array(q)
        Finv = Finv + np.outer(q, q) / sigma ** 2
    w = np.linalg.eigvalsh(Finv)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(Finv)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def main():
    print("=" * 90)
    print("A1  the ablation table  (Delta DROPPED throughout, as PHASE2A.md states it)")
    print("=" * 90)
    designs = {
        "dense 2-daily x 2y (n=365)": D.quasi_steady_level(np.arange(0.0, 730.0, 2.0)),
        "straddling K_m (n=365)": np.geomspace(P.Km / 2, P.Km * 2, 365),
        "straddling K_m (n=9)": np.geomspace(P.Km / 2, P.Km * 2, 9),
        "band [m_bar, 5 m_bar] (n=9)": np.geomspace(D.M_BAR, 5 * D.M_BAR, 9),
    }
    print(f"{'design':>30} {'as published':>13} {'+fit k_e':>10} "
          f"{'+flat load':>11} {'+both':>9}")
    for name, lv in designs.items():
        a = sd_theta(lv)
        b = sd_theta(lv, ke_free=True)
        c = sd_theta(lv, flat_load=True)
        d = sd_theta(lv, ke_free=True, flat_load=True)
        print(f"{name:>30} {a:13.4f} {b:10.4f} {c:11.4f} {d:9.4f}")

    dn, pl = designs["dense 2-daily x 2y (n=365)"], designs["straddling K_m (n=365)"]
    print(f"\n  placement ratio, published footing (k_e known): "
          f"{sd_theta(dn)/sd_theta(pl):.1f}x")
    print(f"  placement ratio, A1's recommended footing (k_e fitted, 0.8 log-prior): "
          f"{sd_theta(dn, ke_free=True)/sd_theta(pl, ke_free=True):.2f}x")
    inc = np.asarray(D.delta(np.arange(0.0, 730.0, 2.0)), dtype=float)
    print(f"  placement ratio, 2C footing (Delta carried, k_e fitted): "
          f"{sd_theta(dn, inc, ke_free=True)/sd_theta(pl, inc, ke_free=True):.2f}x")

    print()
    print("=" * 90)
    print("A1  the excursion curve  (PHASE2A.md: 'the curve does not plateau, it turns')")
    print("=" * 90)
    tops = (2.0, 5.0, 10.0, 20.0, 40.0, 80.0, 160.0, 320.0)
    top, err, f_at, curve = D.excursion_for_decision(
        target=0.05, mu_frac=0.30, reps=300, seed=31, n=9, tops=tops)
    print(f"  asked at an illness of {f_at:.4f}x baseline load (mu = -0.30 theta)")
    print("  " + "  ".join(f"{t:>6.0f}" for t, _ in curve))
    print("  " + "  ".join(f"{100*e:>5.1f}%" for _, e in curve))
    lo = int(np.argmin([e for _, e in curve]))
    print(f"  minimum {100*err:.1f}% at {curve[lo][0]:.0f} m_bar, then "
          f"{'RISING' if curve[-1][1] > err else 'flat'} to "
          f"{100*curve[-1][1]:.1f}% at {curve[-1][0]:.0f} m_bar")
    print(f"  first top reaching the 5% target: {top}")
    print(f"  PAPER PRINTS (08_specification.tex): 36 / 31 / 18 / 10 / 6 / 7 %")

    print()
    print("=" * 90)
    print("A3  excitation available, against the 0.50 the decision needs")
    print("=" * 90)
    print(f"{'source':>42} {'fold':>6} {'E':>9} {'E_max':>9} {'frac of max':>12} "
          f"{'frac of 0.50':>13}")
    rows = [(lab, lo, hi) for lab, lo, hi, _ in D.excursion_sources()]
    rows += [("*Haijes 2020 urine 3-HP, AMD/no-AMD (P=NS)", D.M_BAR, 3.18 * D.M_BAR),
             ("*the '20x' the specification asks for", D.M_BAR, 20.0 * D.M_BAR)]
    for lab, lo, hi in rows:
        r = hi / lo
        E = excitation(lo, hi, P.Km)
        Em = excitation_max(r)
        print(f"{lab:>42} {r:6.2f} {E:9.4f} {Em:9.4f} {E/Em:12.2%} {E/0.50:13.2%}")
    print("  (rows marked * are not in excursion_sources(); the 8x row is the one A3 refutes)")

    print()
    print("  steady load needed at Delta = 0 to hold each excursion top:")
    for k in (2.0, 5.0, 8.0, 10.0, 20.0, 40.0, 80.0):
        m = k * D.M_BAR
        print(f"    {k:5.1f} m_bar -> {D.clearance_A(m, P)/Lb:6.2f}x baseline load")
    base = D.clearance_A(D.M_BAR, P)
    grid = np.linspace(1.0, 20.0, 8001)
    mus = np.array([D.margin_and_gradient(P.theta, P.Km, P.ke, s_peak=f * base)[0]
                    for f in grid])
    f_at = float(np.interp(0.30 * P.theta, -mus, grid))
    B = P.ke * P.Km + P.theta - f_at * base
    m_at = (-B + np.sqrt(B ** 2 + 4.0 * P.ke * f_at * base * P.Km)) / (2.0 * P.ke)
    print(f"\n  the illness amplitude excursion_for_decision itself fixes: "
          f"{f_at:.4f}x load  ->  marker at {m_at/D.M_BAR:.2f} m_bar")


if __name__ == "__main__":
    main()
