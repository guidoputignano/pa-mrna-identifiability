"""
2D.6 -- the decision asymmetry in both regimes, and the regret-equalised
threshold in each.

THE TWO REGIMES.  fit_ensemble regresses the RECORDED load on the marker. In
the managed regime the record tracks the true load to its recording error, and
that is what the paper has always assumed. Under crisis protocol it does not:
natural protein is stopped while catabolic muscle breakdown supplies the load,
so recorded intake falls towards zero exactly as the marker climbs.

One parameter covers both. With a recording fraction rho,

    s_rec = [ L_base + rho (s_true - L_base) ] (1 + eta) + assay noise ,

rho = 1 is the managed regime and reproduces fit_ensemble exactly; rho = 0 is
the recorded load pinned at baseline while the true load rises; rho < 0 is
intake actively cut below baseline, which is what the protocol prescribes.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
from scipy.optimize import least_squares
import design.run as D

P = D.P
_CACHE = {}


def fit_regime(levels, rho, reps=400, seed=0, sigma_load=None, cv=None):
    sigma_load = D.LOAD_FRAC if sigma_load is None else sigma_load
    cv = D.ASSAY_CV if cv is None else cv
    lv = np.asarray(levels, dtype=float)
    key = (tuple(np.round(lv, 9)), float(rho), int(reps), int(seed))
    if key in _CACHE:
        return _CACHE[key]
    rng = np.random.default_rng(seed)
    b_loc = P.theta * P.Km / (P.Km + D.M_BAR) ** 2 + P.ke
    s_true = P.theta * lv / (P.Km + lv) + P.ke * lv
    L_base = D.clearance_A(D.M_BAR, P)
    s_rec0 = L_base + rho * (s_true - L_base)
    out = []
    for _ in range(int(reps)):
        s_obs = (s_rec0 * (1.0 + rng.normal(0.0, sigma_load, lv.size))
                 + rng.normal(0.0, b_loc * cv * D.M_BAR, lv.size))

        def res(q):
            t, K, k = np.exp(np.clip(q[:3], -25.0, 25.0))
            return t * lv / (K + lv) + k * lv + q[3] - s_obs
        start = np.concatenate([np.log([P.theta, P.Km, P.ke])
                                + rng.normal(0.0, 0.5, 3), [0.0]])
        try:
            r = least_squares(res, start, method="lm", max_nfev=3000)
        except Exception:
            continue
        out.append(np.exp(np.clip(r.x[:3], -25.0, 25.0)))
    arr = np.array(out)
    _CACHE[key] = arr
    return arr


def mu_hats(ens, s_peak):
    return np.array([D.margin_and_gradient(t, K, k, D.M_TOL, s_peak)[0]
                     for t, K, k in ens])


def err_at(ens, s_peak, tau):
    """Error rate of the rule 'escalate if mu_hat < tau' at this patient."""
    mu_t = D.margin_and_gradient(P.theta, P.Km, P.ke, D.M_TOL, s_peak)[0]
    mh = mu_hats(ens, s_peak)
    if mu_t < 0:                      # truly needs escalation
        return float(np.mean(mh >= tau)), "fail to treat", mu_t
    return float(np.mean(mh < tau)), "treat unneeded", mu_t


def main():
    print(__doc__)
    base = D.clearance_A(D.M_BAR, P)
    lv9 = np.geomspace(D.M_BAR, 5 * D.M_BAR, 9)

    # the illness amplitudes that put the patient at +/- the same |mu|
    grid = np.linspace(1.0, 20.0, 8001)
    mus = np.array([D.margin_and_gradient(P.theta, P.Km, P.ke, s_peak=f * base)[0]
                    for f in grid])
    f_lo = float(np.interp(-339.0, -mus, grid))     # mu = +339, safe side
    f_hi = float(np.interp(+339.0, -mus, grid))     # mu = -339, unsafe side
    print(f"  illness amplitudes at mu = +/-339: {f_lo:.4f}x and {f_hi:.4f}x baseline")
    print(f"  (the published table quotes 2.98x and 3.58x; mu there is "
          f"{D.margin_and_gradient(P.theta,P.Km,P.ke,s_peak=2.98*base)[0]:+.0f} and "
          f"{D.margin_and_gradient(P.theta,P.Km,P.ke,s_peak=3.58*base)[0]:+.0f})")
    f_lo, f_hi = 2.98, 3.58          # use the published amplitudes exactly

    print("=" * 82)
    print("THE ASYMMETRY, at |mu| = 339 either side of the boundary, tau = 0")
    print("=" * 82)
    print(f"{'regime':>34} {'treat unneeded':>16} {'fail to treat':>15} {'ratio':>9}")
    rows = {}
    for lab, rho in (("managed (rho = 1, as published)", 1.0),
                     ("partial recording (rho = 0.5)", 0.5),
                     ("crisis, intake pinned (rho = 0)", 0.0),
                     ("crisis, protein stopped (rho = -0.25)", -0.25)):
        ens = fit_regime(lv9, rho, reps=400, seed=21)
        e_safe, _, _ = err_at(ens, f_lo * base, 0.0)
        e_uns, _, _ = err_at(ens, f_hi * base, 0.0)
        rows[lab] = (rho, ens, e_safe, e_uns)
        r = (e_uns / e_safe) if e_safe > 0 else np.inf
        print(f"{lab:>34} {e_safe:15.1%} {e_uns:14.1%} {r:9.1f}")
    print()
    print("  rho = 1 'fails towards declaring patients safe': the costly error dominates.")
    print("  rho <= 0 reverses it: the fit is driven to low theta, so every patient")
    print("  reads as needing escalation and the cheap error dominates instead.")

    print()
    print("=" * 82)
    print("REGRET-EQUALISED THRESHOLD, per regime")
    print("=" * 82)
    print("  Rule: escalate if mu_hat < tau. Regret is equalised where")
    print("  C_fail * P(fail to treat) = C_unneeded * P(treat unneeded) at +/-|mu| = 339.")
    print(f"{'regime':>34} {'C ratio 1:1':>13} {'10:1':>10} {'30:1':>10} {'100:1':>10}")
    for lab, (rho, ens, _, _) in rows.items():
        cells = []
        for R in (1.0, 10.0, 30.0, 100.0):
            taus = np.linspace(-30000.0, 30000.0, 3001)
            gap = np.array([R * err_at(ens, f_hi * base, tau)[0]
                            - err_at(ens, f_lo * base, tau)[0] for tau in taus])
            sgn = np.where(np.diff(np.sign(gap)) != 0)[0]
            if sgn.size:
                cells.append(float(taus[sgn[0]]) / P.theta)
            elif gap[0] > 0:
                cells.append(np.inf)     # fail-to-treat regret always dominates
            else:
                cells.append(-np.inf)
        def fmt(c):
            if np.isposinf(c):
                return "escalate all"
            if np.isneginf(c):
                return "escalate none"
            return f"{c:.3f}"
        print(f"{lab:>34} " + " ".join(f"{fmt(c):>12}" for c in cells))
    print("  (tau* quoted as a fraction of theta; tau* > 0 means escalate EARLIER")
    print("   than the nominal boundary, tau* < 0 means later)")

    print()
    print("=" * 82)
    print("WHAT THE FIT DOES TO theta_hat IN EACH REGIME (the mechanism)")
    print("=" * 82)
    print(f"{'regime':>34} {'median theta_hat':>17} {'/ theta':>9} {'IQR in decades':>16}")
    for lab, (rho, ens, _, _) in rows.items():
        th = ens[:, 0]
        q1, q3 = np.percentile(th, [25, 75])
        print(f"{lab:>34} {np.median(th):17.1f} {np.median(th)/P.theta:9.3f} "
              f"{np.log10(q3/q1):16.2f}")


if __name__ == "__main__":
    main()
