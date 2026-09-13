"""
3.0.2 and 3.0.4 -- the rate estimand.

THE MAP. The natural history study reports that a 50% lower plasma 3-HP is
associated with a decompensation rate ratio of 0.7. Read as a log-linear
dose-response over the observed range,

    rate(m) = rate_ref * (m / m_ref)^beta ,   beta = -log2(0.70) = 0.5146 ,

so the rate is a POWER of the marker level with an exponent below one. The
model supplies the level: at a recorded load L the patient sits where
(theta + Delta) g(m) + k_e m = L. So the chain is

    (theta_hat, K_m_hat, k_e_hat)  ->  m_hat(L)  ->  rate_hat .

3.0.2 asks whether the cancellation that made mu stable survives this chain.
3.0.4 asks what precision on the delivered increment a stated precision on the
rate requires.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
import design.run as D

P = D.P
BETA = -np.log2(0.70)


def level_at(theta, Km, ke, L, dlt=0.0):
    """Quasi-steady marker level at a recorded load, for any parameter set."""
    B = ke * Km + theta + dlt - L
    return (-B + np.sqrt(B ** 2 + 4.0 * ke * L * Km)) / (2.0 * ke)


def rate_of(m, m_ref, rate_ref=1.0):
    return rate_ref * (m / m_ref) ** BETA


def main():
    print(__doc__)
    print(f"  beta = -log2(0.70) = {BETA:.4f}")
    L0 = D.load(0.0)
    t = np.arange(0.0, 730.0, 30.4)
    lv, inc = D.design_record(t)
    m_ref = float(level_at(P.theta, P.Km, P.ke, L0))
    print(f"  reference level at baseline load: {m_ref:.4f}")

    ens = D.fit_ensemble(lv, P.theta, P.Km, P.ke, reps=400, seed=11, increments=inc)
    th, Km, ke = ens[:, 0], ens[:, 1], ens[:, 2]
    print(f"  ensemble: {len(ens)} fits, theta_hat spans "
          f"{th.min():.3g} to {th.max():.3g}, a factor of {th.max()/th.min():.3g}")

    print()
    print("=" * 86)
    print("3.0.2  DOES THE CANCELLATION SURVIVE?")
    print("=" * 86)
    print("  sd of the log of each quantity across the ensemble, at three loads.")
    print(f"{'load':>22} | {'sd(ln theta_hat)':>17} {'sd(ln m_hat)':>13} "
          f"{'sd(ln rate_hat)':>16} | {'ratio th/rate':>14}")
    for lab, Lf in (("the patient's own load", 1.0),
                    ("1.5x baseline", 1.5),
                    ("3.0x baseline", 3.0)):
        L = Lf * L0
        mh = np.array([level_at(a, b, c, L) for a, b, c in zip(th, Km, ke)])
        rh = rate_of(mh, m_ref)
        s_th = float(np.std(np.log(th), ddof=1))
        s_m = float(np.std(np.log(mh), ddof=1))
        s_r = float(np.std(np.log(rh), ddof=1))
        print(f"{lab:>22} | {s_th:17.4f} {s_m:13.4f} {s_r:16.4f} | "
              f"{s_th/max(s_r,1e-12):14.1f}")

    print()
    print("  The fitted parameters are the ones that REPRODUCE the record, so at the")
    print("  load the record was taken at they all predict nearly the same level: the")
    print("  fibre is exactly the set that agrees there. The cancellation is therefore")
    print("  stronger than the one the margin enjoyed, and it is exact at the operating")
    print("  point rather than approximate. It degrades as the load is extrapolated.")

    print()
    print("=" * 86)
    print("3.0.2b  AGAINST THE MARGIN, on the same ensemble")
    print("=" * 86)
    base = D.clearance_A(D.M_BAR, P)
    mu_h = np.array([D.margin_and_gradient(a, b, c, D.M_TOL, 3.58 * base)[0]
                     for a, b, c in zip(th, Km, ke)])
    mu_t = D.margin_and_gradient(P.theta, P.Km, P.ke, D.M_TOL, 3.58 * base)[0]
    mh0 = np.array([level_at(a, b, c, L0) for a, b, c in zip(th, Km, ke)])
    rh0 = rate_of(mh0, m_ref)
    print(f"    sign error on mu at 3.58x:                 "
          f"{100*np.mean(np.sign(mu_h) != np.sign(mu_t)):.1f}%")
    print(f"    relative sd of the rate at the same design: "
          f"{100*np.std(rh0, ddof=1)/np.mean(rh0):.2f}%")
    q = np.percentile(rh0, [2.5, 97.5])
    print(f"    95% of the ensemble's rate estimates lie in "
          f"[{q[0]:.4f}, {q[1]:.4f}] x the reference rate")
    print("    (the sign is a coin toss on the same fits that pin the rate to a "
          "few per cent)")

    print()
    print("=" * 86)
    print("3.0.4  WHAT PRECISION ON THE INCREMENT DOES A RATE TARGET NEED?")
    print("=" * 86)
    print("  Because rate ~ m^beta with beta = 0.515, a relative error in the level")
    print("  is COMPRESSED by beta into the rate: 1% on m is 0.51% on the rate.")
    print()
    print(f"{'target on the rate':>22} {'implied target on m':>21} "
          f"{'implied on ln theta':>21}")
    for tgt in (0.05, 0.10, 0.20, 0.30):
        on_m = tgt / BETA
        # d ln m / d ln theta at the operating point, from the balance
        h = 1e-6
        dlnm = (np.log(level_at(P.theta * np.exp(h), P.Km, P.ke, L0))
                - np.log(level_at(P.theta * np.exp(-h), P.Km, P.ke, L0))) / (2 * h)
        print(f"{tgt:21.0%} {on_m:20.1%} {on_m/abs(dlnm):20.1%}")
    h = 1e-6
    dlnm = (np.log(level_at(P.theta * np.exp(h), P.Km, P.ke, L0))
            - np.log(level_at(P.theta * np.exp(-h), P.Km, P.ke, L0))) / (2 * h)
    print()
    print(f"  d ln m / d ln theta at the operating point = {dlnm:.4f}")
    print(f"  so the level is {1/abs(dlnm):.1f}x LESS sensitive than theta, and the")
    print(f"  rate is a further {1/BETA:.2f}x less sensitive than the level.")
    print(f"  Net: a 20% target on the rate needs {0.20/BETA/abs(dlnm):.0%} on theta.")


if __name__ == "__main__":
    main()


def does_the_rate_need_the_channel():
    """
    The question 3.0.2 raises by answering well: if the rate is pinned even
    where theta is unidentified, does the therapy channel still matter?
    """
    print()
    print("=" * 86)
    print("3.0.2c  IF THE RATE IS PINNED ANYWAY, WHAT IS THE THERAPY CHANNEL FOR?")
    print("=" * 86)
    L0 = D.load(0.0)
    m_ref = float(level_at(P.theta, P.Km, P.ke, L0))
    designs = {
        "six-monthly (n=4)": np.arange(180.0, 731.0, 180.0),
        "monthly (n=25)": np.arange(0.0, 730.0, 30.4),
        "fortnightly = dose interval (n=53)": np.arange(0.0, 730.0, D.DOSE_INTERVAL),
    }
    print(f"{'design':>36} {'increment':>10} | {'sd(ln theta)':>13} "
          f"{'sd(ln rate) at L0':>18} {'at 3 L0':>10}")
    for nm, t in designs.items():
        lv, inc = D.design_record(t)
        for carry, lab in ((True, "carried"), (False, "dropped")):
            ens = D.fit_ensemble(lv, P.theta, P.Km, P.ke, reps=250, seed=11,
                                 increments=inc if carry else None)
            th, Km, ke = ens[:, 0], ens[:, 1], ens[:, 2]
            m0 = np.array([level_at(a, b, c, L0) for a, b, c in zip(th, Km, ke)])
            m3 = np.array([level_at(a, b, c, 3 * L0) for a, b, c in zip(th, Km, ke)])
            s_th = float(np.std(np.log(th), ddof=1))
            s_r0 = BETA * float(np.std(np.log(m0), ddof=1))
            s_r3 = BETA * float(np.std(np.log(m3), ddof=1))
            print(f"{(nm if carry else ''):>36} {lab:>10} | {s_th:13.3f} "
                  f"{s_r0:18.4f} {s_r3:10.4f}")
    print()
    print("  Read the last two columns. If the rate is pinned with the increment")
    print("  DROPPED as well, the therapy channel is not what makes the rate")
    print("  estimable, and the paper must say so.")


def kappa_for_a_rate_target():
    """3.0.4, by propagating the Fisher covariance rather than simulating."""
    print()
    print("=" * 86)
    print("3.0.4b  THE INCREMENT PRECISION A RATE TARGET NEEDS, by propagation")
    print("=" * 86)
    import model_A.run as LOC
    from scipy.optimize import brentq
    L0 = D.load(0.0)
    SIG = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)
    t = np.arange(0.0, 730.0, 30.4)
    lv, inc = D.design_record(t)

    # gradient of ln m_hat at the operating load, in the same coordinates
    h = 1e-6
    g = []
    for i, (a, b, c) in enumerate(((P.theta, P.Km, P.ke),)):
        pass
    def lnm(th, km, k, bias):
        return np.log(level_at(th, km, k, L0 + bias))
    base_args = (P.theta, P.Km, P.ke, 0.0)
    g = np.array([
        (lnm(P.theta*np.exp(h), P.Km, P.ke, 0.0) - lnm(P.theta*np.exp(-h), P.Km, P.ke, 0.0))/(2*h),
        (lnm(P.theta, P.Km*np.exp(h), P.ke, 0.0) - lnm(P.theta, P.Km*np.exp(-h), P.ke, 0.0))/(2*h),
        (lnm(P.theta, P.Km, P.ke, h*L0) - lnm(P.theta, P.Km, P.ke, -h*L0))/(2*h)/L0,
        (lnm(P.theta, P.Km, P.ke*np.exp(h), 0.0) - lnm(P.theta, P.Km, P.ke*np.exp(-h), 0.0))/(2*h),
    ])
    print(f"  gradient of ln m at the operating load, (ln th, ln Km, bias, ln ke):")
    print(f"    {np.array2string(g, precision=5)}")

    def sd_rate(kap=None, ke_prior=0.8):
        pri = [np.inf, np.inf, (D.LOAD_FRAC*L0)**2, ke_prior**2]
        gg = list(g)
        if kap is not None:
            pri.append(kap**2); gg.append(0.0)
        F = np.diag([0.0 if np.isinf(x) else 1.0/x for x in pri])
        for m, d in zip(lv, inc):
            gm = m/(P.Km+m)
            q = [P.theta*gm, -(P.theta+d)*P.Km*m/(P.Km+m)**2, 1.0, P.ke*m]
            if kap is not None:
                q.append(d*gm)
            q = np.array(q); F = F + np.outer(q, q)/SIG**2
        w = np.linalg.eigvalsh(F)
        if w[0] <= 0 or w[0] <= 1e-14*w[-1]:
            return np.inf
        C = np.linalg.inv(F)
        gg = np.array(gg)
        return float(BETA*np.sqrt(gg @ C @ gg))

    print()
    print(f"{'increment known to':>22} {'sd(ln rate)':>13} {'against 20% target':>20}")
    for kap in (None, 0.05, 0.117, 0.20, 0.40, 0.80, 2.0):
        v = sd_rate(kap)
        lab = "exactly" if kap is None else f"{kap:.1%}"
        print(f"{lab:>22} {v:13.4f} {'PASS' if v < 0.20 else 'fail':>20}")
    for tgt in (0.05, 0.10, 0.20):
        try:
            c = brentq(lambda x: sd_rate(x) - tgt, 1e-3, 20.0, xtol=1e-5)
            print(f"  a {tgt:.0%} target on the rate needs the increment to {c:.1%}")
        except ValueError:
            v0 = sd_rate(20.0)
            print(f"  a {tgt:.0%} target on the rate: "
                  f"{'met at ANY increment precision' if v0 < tgt else 'unreachable'}"
                  f"  (flat-kappa value {v0:.4f})")


does_the_rate_need_the_channel()
kappa_for_a_rate_target()
