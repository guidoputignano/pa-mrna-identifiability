"""
2E.1 -- what actually binds the decision?

The sign error on mu = theta - theta_crit is decomposed into its sources by
recomputing mu_hat from the SAME fitted ensemble four ways:

    full    mu_hat = theta_hat - theta_crit(K_m_hat, k_e_hat)   as reported
    (a)     mu_hat = theta_hat - theta_crit(K_m,     k_e    )   theta only
    (b)     mu_hat = theta     - theta_crit(K_m_hat, k_e_hat)   threshold only
    none    mu     = theta     - theta_crit(K_m,     k_e    )   no error

One correction to the question as posed. The 511-3285 band (a factor of 6.4) is
the band on the MARGIN -- separatrix minus attractor, in metabolite units --
and not on theta_crit. Over the same M_i sweep theta_crit_intrinsic moves by a
factor of 1.37, and over K_c by 1.56. So the premise "theta_crit is known only
to an order of magnitude" is not what the 2B.3 band says. The hypothesis is
still worth testing, at the right magnitude, and with the convention constants
s_peak and m_tol included, which are arguably the looser ones.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
import design.run as D

P = D.P
base = D.clearance_A(D.M_BAR, P)


def mu_of(theta, Km, ke, s_peak, m_tol=None):
    m_tol = D.M_TOL if m_tol is None else m_tol
    return D.margin_and_gradient(theta, Km, ke, m_tol, s_peak)[0]


def decompose(lv, inc, s_peak, reps=400, seeds=range(8), m_tol=None):
    """Sign error under each combination of error sources."""
    out = {k: [] for k in ("full", "a_theta_only", "b_threshold_only", "none")}
    for sd in seeds:
        ens = D.fit_ensemble(lv, P.theta, P.Km, P.ke, reps=reps, seed=sd,
                             increments=inc)
        mu_true = mu_of(P.theta, P.Km, P.ke, s_peak, m_tol)
        sgn = np.sign(mu_true)
        full = np.array([mu_of(t, K, k, s_peak, m_tol) for t, K, k in ens])
        a = np.array([mu_of(t, P.Km, P.ke, s_peak, m_tol) for t, _, _ in ens])
        b = np.array([mu_of(P.theta, K, k, s_peak, m_tol) for _, K, k in ens])
        out["full"].append(np.mean(np.sign(full) != sgn))
        out["a_theta_only"].append(np.mean(np.sign(a) != sgn))
        out["b_threshold_only"].append(np.mean(np.sign(b) != sgn))
        out["none"].append(0.0)
    return {k: (float(np.mean(v)), float(np.std(v))) for k, v in out.items()}


def main():
    print(__doc__)
    t_m = np.arange(0.0, 730.0, 30.4)
    lv, inc = D.design_record(t_m)

    print("=" * 84)
    print("(a) AGAINST (b): which error source drives the sign error?")
    print("=" * 84)
    print("Monthly design, 8 seeds x 400 replicates, at the two amplitudes the")
    print("decision table quotes.")
    print()
    print(f"{'amplitude':>12} {'mu_true':>9} | {'full':>14} {'theta only (a)':>16} "
          f"{'threshold only (b)':>20}")
    for f, lab in ((2.98, "safe side"), (3.58, "unsafe side")):
        r = decompose(lv, inc, f * base)
        mt = mu_of(P.theta, P.Km, P.ke, f * base)
        print(f"{f:11.2f}x {mt:9.0f} | "
              f"{100*r['full'][0]:6.1f}% +/-{100*r['full'][1]:4.1f} "
              f"{100*r['a_theta_only'][0]:8.1f}% +/-{100*r['a_theta_only'][1]:4.1f} "
              f"{100*r['b_threshold_only'][0]:12.1f}% +/-{100*r['b_threshold_only'][1]:4.1f}")

    print()
    print("  If (a) and (b) are each much smaller than `full`, the two error")
    print("  sources are ANTI-correlated and the margin theorem is working.")
    print("  If (b) ~ full, the boundary is what is not located.")

    print()
    print("=" * 84)
    print("(c) THE DISTANCE FROM THE BOUNDARY")
    print("=" * 84)
    print(f"{'illness':>9} {'mu_true':>10} {'|mu|/theta':>11} | {'full':>14} "
          f"{'(a)':>12} {'(b)':>12}")
    for f in (2.0, 2.6, 2.98, 3.28, 3.58, 4.07, 5.0, 6.5):
        r = decompose(lv, inc, f * base, seeds=range(4))
        mt = mu_of(P.theta, P.Km, P.ke, f * base)
        print(f"{f:8.2f}x {mt:10.0f} {abs(mt)/P.theta:11.3f} | "
              f"{100*r['full'][0]:10.1f}% {100*r['a_theta_only'][0]:11.1f}% "
              f"{100*r['b_threshold_only'][0]:11.1f}%")

    print()
    print("=" * 84)
    print("(b') THE CONVENTION CONSTANTS, which are not in the simulation at all")
    print("=" * 84)
    print("theta_crit_local = (s_peak - k_e m_tol)(K_m + m_tol)/m_tol depends on")
    print("two quantities nobody measures: the illness amplitude s_peak and the")
    print("tolerable level m_tol. Their effect on the THRESHOLD, at truth:")
    print()
    print(f"{'m_tol':>8} | " + "  ".join(f"s={f:.2f}x" for f in (2.98, 3.28, 3.58)))
    for mt_ in (20.0, 40.0, 80.0, 160.0):
        cells = [D.theta_crit_local(P.ke, P.Km, m_tol=mt_, s_peak=f * base)
                 for f in (2.98, 3.28, 3.58)]
        print(f"{mt_:8.0f} | " + "  ".join(f"{c:9.0f}" for c in cells))
    lo = min(D.theta_crit_local(P.ke, P.Km, m_tol=mt_, s_peak=f * base)
             for mt_ in (20.0, 40.0, 80.0, 160.0) for f in (2.98, 3.28, 3.58))
    hi = max(D.theta_crit_local(P.ke, P.Km, m_tol=mt_, s_peak=f * base)
             for mt_ in (20.0, 40.0, 80.0, 160.0) for f in (2.98, 3.28, 3.58))
    print()
    print(f"  threshold range over those choices: {lo:.0f} to {hi:.0f}, "
          f"a factor of {hi/lo:.1f}")
    print(f"  against theta = {P.theta:.0f}, so the SIGN of mu is decided by the")
    print(f"  convention over that whole range.")
    print()
    print("  For comparison, the model-C shape constants move theta_crit_intrinsic")
    print("  by a factor of 1.37 (M_i) and 1.56 (K_c). The 6.4-fold band of 2B.3 is")
    print("  on the MARGIN in metabolite units, not on theta_crit.")


if __name__ == "__main__":
    main()
