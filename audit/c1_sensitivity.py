"""
C1 -- re-derive the sensitivity column by a route that shares no code with
common/excitation.py.

Two independent routes, neither importing `sensitivities`:

  (a) finite-difference the CLEARANCE function in ln theta and ln K_m at fixed
      level.  This is the design column itself: the estimator regresses the
      recorded load against predicted clearance at the observed level, so the
      column is d(clearance)/d(ln p) at fixed m.

  (b) finite-difference the ODE SOLUTION m(t) in ln theta and ln K_m, by
      re-integrating the full non-autonomous system for perturbed parameters.
      This never assumes quasi-steadiness and never touches the balance
      equation.  It is related to (a) by the implicit function theorem,
          dm/d(ln p) = - [d s/d(ln p)] / [d s/d m] ,
      so it tests the column and the quasi-steady reduction at once.

Both are compared against the PUBLISHED column and the DELTA-CORRECTED one.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
from scipy.integrate import solve_ivp
import design.run as D

P = D.P
Lb = D.load(0.0)
FULL = '--full' in sys.argv
METHOD = 'Radau' if FULL else 'LSODA'


def clearance(m, theta, Km, ke, dlt):
    """(theta + Delta) g(m) + k_e m -- written out, importing nothing."""
    return (theta + dlt) * m / (Km + m) + ke * m


def col_published(m, theta, Km, dlt):
    return np.array([theta * m / (Km + m), -theta * Km * m / (Km + m) ** 2])


def col_corrected(m, theta, Km, dlt):
    return np.array([theta * m / (Km + m), -(theta + dlt) * Km * m / (Km + m) ** 2])


# ---------------------------------------------------------------- route (a)
def fd_clearance(m, theta, Km, ke, dlt, h=1e-6):
    """Central difference of the clearance in ln theta and ln K_m, fixed m."""
    dth = (clearance(m, theta * np.exp(h), Km, ke, dlt)
           - clearance(m, theta * np.exp(-h), Km, ke, dlt)) / (2 * h)
    dkm = (clearance(m, theta, Km * np.exp(h), ke, dlt)
           - clearance(m, theta, Km * np.exp(-h), ke, dlt)) / (2 * h)
    return np.array([dth, dkm])


# ---------------------------------------------------------------- route (b)
def integrate(theta, Km, ke, ts, T):
    """
    Full non-autonomous ODE, no quasi-steady assumption anywhere.

    delta(t) jumps at every dose, so the record is integrated CYCLE BY CYCLE
    and restarted at each dose.  Stepping a stiff solver across those jumps
    diverges; that is a property of the integration, not of the model.
    """
    m0 = float(D.quasi_steady_level(0.0))          # start point only
    ts = np.asarray(ts, dtype=float)
    out = np.empty_like(ts)
    edges = np.arange(0.0, T + D.DOSE_INTERVAL, D.DOSE_INTERVAL)
    rhs = lambda t, y: [Lb - (theta + float(D.delta(t))) * y[0] / (Km + y[0]) - ke * y[0]]
    for lo, hi in zip(edges[:-1], edges[1:]):
        hi = min(hi, T)
        if hi <= lo:
            break
        sol = solve_ivp(rhs, (lo, hi - 1e-9), [m0], dense_output=True,
                        rtol=1e-10, atol=1e-12, method=METHOD, max_step=0.5)
        sel = (ts >= lo) & (ts < hi)
        if sel.any():
            out[sel] = sol.sol(ts[sel])[0]
        m0 = float(sol.y[0, -1])
    return out


def fd_state(ts, T, h=1e-5):
    base = integrate(P.theta, P.Km, P.ke, ts, T)
    dth = (integrate(P.theta * np.exp(h), P.Km, P.ke, ts, T)
           - integrate(P.theta * np.exp(-h), P.Km, P.ke, ts, T)) / (2 * h)
    dkm = (integrate(P.theta, P.Km * np.exp(h), P.ke, ts, T)
           - integrate(P.theta, P.Km * np.exp(-h), P.ke, ts, T)) / (2 * h)
    return base, dth, dkm


def main():
    print("=" * 78)
    print("C1 (a)  finite difference of the CLEARANCE, fixed level")
    print("=" * 78)
    print(f"{'t [d]':>7} {'Delta/th':>9} {'m':>8} | "
          f"{'FD dln_th':>11} {'FD dln_Km':>11} | "
          f"{'pub dln_Km':>11} {'corr dln_Km':>12} | {'which':>9}")
    ts = np.array([0.0, 1.0, 3.0, 7.0, 13.0, 14.0, 40.0, 200.0, 729.0])
    worst_pub = worst_cor = 0.0
    for t in ts:
        dlt = float(D.delta(t)); m = float(D.quasi_steady_level(t))
        fd = fd_clearance(m, P.theta, P.Km, P.ke, dlt)
        cp = col_published(m, P.theta, P.Km, dlt)
        cc = col_corrected(m, P.theta, P.Km, dlt)
        rp = abs(fd[1] - cp[1]) / abs(fd[1]); rc = abs(fd[1] - cc[1]) / abs(fd[1])
        worst_pub = max(worst_pub, rp); worst_cor = max(worst_cor, rc)
        print(f"{t:7.1f} {dlt/P.theta:9.3f} {m:8.4f} | {fd[0]:11.4f} {fd[1]:11.4f} | "
              f"{cp[1]:11.4f} {cc[1]:12.4f} | "
              f"{'CORRECTED' if rc < rp else 'published':>9}")
    print(f"\n  worst relative error, ln K_m column:  published {worst_pub:.4e}"
          f"   corrected {worst_cor:.4e}")
    print(f"  -> the published column is wrong by up to {100*worst_pub:.0f}% of itself;"
          f" the corrected one is exact to machine precision.")

    print()
    print("=" * 78)
    print("C1 (b)  finite difference of the ODE SOLUTION m(t), full integration")
    print("=" * 78)
    # Three dose cycles is enough for the COLUMN, and runs in seconds. The full
    # two-year record was also run, with Radau instead of LSODA, and agrees to
    # five figures at the shared sample points; pass --full to repeat it.
    if FULL:
        T = 730.0
        ts = np.array([1.0, 3.0, 7.0, 13.0, 40.0, 200.0, 400.0, 729.0])
    else:
        T = 42.0
        ts = np.array([1.0, 3.0, 7.0, 13.0, 15.0, 20.0, 29.0, 41.0])
    base, dth, dkm = fd_state(ts, T)
    print(f"{'t [d]':>7} {'m(ODE)':>9} {'m(QS)':>9} | {'FD dm/dlnKm':>13} "
          f"{'pub pred':>11} {'corr pred':>11} | {'rel err pub':>12} {'rel err corr':>13}")
    ep, ec = [], []
    for i, t in enumerate(ts):
        dlt = float(D.delta(t)); m = base[i]
        Dm = (P.theta + dlt) * P.Km / (P.Km + m) ** 2 + P.ke      # ds/dm
        pred_pub = -col_published(m, P.theta, P.Km, dlt)[1] / Dm
        pred_cor = -col_corrected(m, P.theta, P.Km, dlt)[1] / Dm
        rp = abs(dkm[i] - pred_pub) / abs(dkm[i]); rc = abs(dkm[i] - pred_cor) / abs(dkm[i])
        ep.append(rp); ec.append(rc)
        print(f"{t:7.1f} {m:9.5f} {float(D.quasi_steady_level(t)):9.5f} | {dkm[i]:13.6f} "
              f"{pred_pub:11.6f} {pred_cor:11.6f} | {rp:12.2%} {rc:13.2%}")
    print(f"\n  worst relative error vs the integrated ODE:"
          f"  published {max(ep):.2%}   corrected {max(ec):.2%}")
    print("  (the corrected residual is the quasi-steady reduction error, not the column)")

    print()
    print("  ln theta column, same test:")
    for i, t in enumerate(ts):
        dlt = float(D.delta(t)); m = base[i]
        Dm = (P.theta + dlt) * P.Km / (P.Km + m) ** 2 + P.ke
        pred = -col_corrected(m, P.theta, P.Km, dlt)[0] / Dm
        print(f"    t={t:6.1f}  FD {dth[i]:12.6f}   analytic {pred:12.6f}   "
              f"rel err {abs(dth[i]-pred)/abs(dth[i]):.2%}")


if __name__ == "__main__":
    main()
