"""
Item 10. Is the one-compartment Michaelis-Menten specification right?

The honest constraint first: the repository has NO patient data. Every record is
SIMULATED from the one-compartment MM model, so residuals from the MM fit are
white by construction -- that tests the estimator, not the model. A real
specification test needs real data.

What is testable here is the DIAGNOSTIC's POWER: if real data came from a
different mechanism (a second compartment, a saturable non-enzymatic clearance,
or a load that does not match its record), would the residuals from an MM fit
show it, at the repository's noise and sampling? This simulates each
misspecified data-generating process, fits the MM model, and reports the
residual signature and the power to detect it on the monthly record against a
dense excursion.
"""
import numpy as np
from scipy.optimize import least_squares

import core as C
import dyn as D

RNG = np.random.default_rng(31102026)
L0, M0 = C.record(C.TH_5, seed=0)["L"], None


# ------------------------------------------------------------- the MM fit
def mm_fit(loads, inc, m_obs):
    """MAP fit of (ln th, ln Km, b, ln ke) to (recorded load, marker) pairs."""
    sb = C.LOAD_FRAC * L0
    sg = np.full(m_obs.size, 0.05 * L0)

    def res(x):
        lt, lk, b, lke = x
        th, K, ke = np.exp(lt), np.exp(lk), np.exp(lke)
        pred = (th + inc) * m_obs / (K + m_obs) + ke * m_obs + b
        return np.r_[(loads - pred) / sg, b / sb, (lke - np.log(C.KE)) / C.KE_PRIOR]
    x0 = np.array([np.log(C.TH_5), np.log(C.KM0), 0.0, np.log(C.KE)])
    s = least_squares(res, x0, method="lm", xtol=1e-12, ftol=1e-12, max_nfev=4000)
    lt, lk, b, lke = s.x
    th, K, ke = np.exp(lt), np.exp(lk), np.exp(lke)
    pred = (th + inc) * m_obs / (K + m_obs) + ke * m_obs + b
    return loads - pred            # residuals in flux (= recorded-load) units


def slope_t(x, r):
    """t-statistic of the OLS slope of residual r on regressor x (centred)."""
    x = x - x.mean()
    if np.allclose(x, 0):
        return 0.0
    b = (x @ r) / (x @ x)
    dof = max(x.size - 2, 1)
    s2 = ((r - b * x) ** 2).sum() / dof
    se = np.sqrt(s2 / (x @ x))
    return 0.0 if se == 0 else b / se


# ---------------------------------------------- misspecified steady states
def m_true(dgp, load, inc, **kw):
    """solve the TRUE balance of the named DGP at one (load, inc)."""
    th, K, ke = C.TH_5, C.KM0, C.KE
    if dgp == "mm":
        return C.m_ss(th, K, ke, load, inc)
    if dgp == "sat_ne":          # k_e m -> Vne m/(Kne+m), matched at m_ref
        Kne = kw["Kne"]
        Vne = ke * (Kne + kw["m_ref"])          # so Vne m/(Kne+m) ~ ke m near m_ref
        # solve (th+inc) m/(K+m) + Vne m/(Kne+m) = load
        m = C.m_ss(th, K, ke, load, inc)
        for _ in range(60):
            f = (th + inc) * m / (K + m) + Vne * m / (Kne + m) - load
            fp = ((th + inc) * K / (K + m) ** 2 + Vne * Kne / (Kne + m) ** 2)
            m = max(m - f / fp, 1e-6)
        return m
    raise ValueError(dgp)


def monthly_power(dgp, nmc=300, **kw):
    """power of each residual test on the monthly record under a DGP."""
    rej = np.zeros(3)
    sig = np.zeros(3)
    for j in range(nmc):
        r = C.record(C.TH_5, seed=j % 24)
        loads, inc, t = r["loads"], r["inc"], r["times"]
        mt = np.array([m_true(dgp, loads[i], inc[i], **kw) for i in range(loads.size)])
        m_obs = mt * np.exp(RNG.normal(0, 0.05, mt.size))      # 5% assay
        res = mm_fit(loads, inc, m_obs)
        ts = np.array([slope_t(m_obs, res), slope_t(t, res), slope_t(loads, res)])
        rej += np.abs(ts) > 1.96
        sig += ts
    return rej / nmc, sig / nmc


def excursion_power(dgp, nmc=200, dt_min=5.0, assay=0.05, **kw):
    """power on a within-record load ramp over 48 h, sampling step dt_min minutes."""
    tt = np.arange(0.0, 2.0 + 1e-9, dt_min / 60 / 24)
    loads = L0 * (0.6 + 0.9 * tt / tt[-1])      # 0.6L -> 1.5L, within fold at the top
    inc = np.zeros_like(tt)
    rej = np.zeros(3)
    for j in range(nmc):
        mt = np.array([m_true(dgp, loads[i], 0.0, **kw) for i in range(loads.size)])
        m_obs = mt * np.exp(RNG.normal(0, assay, mt.size))
        res = mm_fit(loads, inc, m_obs)
        ts = np.array([slope_t(m_obs, res), slope_t(tt, res), slope_t(loads, res)])
        rej += np.abs(ts) > 1.96
    return rej / nmc, tt.size


def main():
    print("=" * 96)
    print("  ITEM 10.  RESIDUAL DIAGNOSTICS: IS THE MM SPECIFICATION RIGHT?")
    print("=" * 96)
    print("\n  CAVEAT FIRST: all data is simulated from the one-compartment MM model,")
    print("  so MM residuals are white by construction. This is a POWER analysis of")
    print("  the diagnostic, not a test of the model. A real test needs real data.")

    print("\n  (0) SANITY: residuals under the TRUE MM model are structureless.")
    p, s = monthly_power("mm", nmc=300)
    print(f"      monthly record, test rejection rate at 5% (should be ~5%):")
    print(f"        vs marker {100*p[0]:.1f}%, vs time {100*p[1]:.1f}%, vs load {100*p[2]:.1f}%")

    print("\n  (1) SATURABLE NON-ENZYMATIC CLEARANCE (k_e m -> Vne m/(Kne+m)).")
    print("      monthly record (n=25): nearly blind.")
    for Kne in (150.0, 60.0):
        pm, sm = monthly_power("sat_ne", Kne=Kne, m_ref=M0 or 20.0)
        print(f"        Kne = {Kne:.0f}: power vs marker {100*pm[0]:.0f}%, vs time {100*pm[1]:.0f}%,"
              f" vs load {100*pm[2]:.0f}%")
    print("      Dense excursion: the residual-slope test is SIZE-INVALID at large N --")
    print("      a noisy nonlinear regressor (errors-in-variables) inflates the false")
    print("      positive rate. Compare the TRUE-MM rejection (should be ~5%) at each N:")
    for dt, lab in ((5.0, "5-min"), (60.0, "hourly")):
        mm_rej, N = excursion_power("mm", dt_min=dt)
        print(f"        {lab} sampling (N={N}):  TRUE-MM false positive"
              f" vs time {100*mm_rej[1]:.0f}%, vs load {100*mm_rej[2]:.0f}%  <- the null")
        for Kne in (150.0, 60.0):
            pe, _ = excursion_power("sat_ne", Kne=Kne, m_ref=M0 or 20.0, dt_min=dt)
            print(f"           Kne={Kne:.0f}: vs time {100*pe[1]:.0f}%, vs load {100*pe[2]:.0f}%")
    print("      So at the dense N the '100%' is mostly the EIV artifact; at the coarse,")
    print("      size-valid N only STRONG saturation (Kne=60) genuinely rejects above the")
    print("      null, while MILD saturation (Kne=150) is still undetected.")

    print("\n  (2) LOAD MISMATCH: a 10%/20% linear drift of true load off its record,")
    print("      beyond the constant bias b the model allows.")
    for dr in (0.10, 0.20):
        rej = np.zeros(3)
        for j in range(300):
            r = C.record(C.TH_5, seed=j % 24)
            loads, inc, t = r["loads"], r["inc"], r["times"]
            Ltrue = loads * (1 + dr * (t / t[-1] - 0.5))     # drift, zero-mean so b can't absorb
            mt = np.array([C.m_ss(C.TH_5, C.KM0, C.KE, Ltrue[i], inc[i]) for i in range(loads.size)])
            m_obs = mt * np.exp(RNG.normal(0, 0.05, mt.size))
            res = mm_fit(loads, inc, m_obs)
            rej += np.abs(np.array([slope_t(m_obs, res), slope_t(t, res), slope_t(loads, res)])) > 1.96
        rej /= 300
        print(f"      drift {int(dr*100)}%: monthly power vs marker {100*rej[0]:.0f}%,"
              f" vs time {100*rej[1]:.0f}%, vs load {100*rej[2]:.0f}%")

    print("\n  (3) SECOND COMPARTMENT: a slow exchange pool. At STEADY STATE the")
    print("      input-output balance is unchanged, so the monthly record is BLIND")
    print("      to it; it only shows in a transient as a second relaxation phase.")
    tt = np.arange(0.0, 2.0 + 1e-9, 5 / 60 / 24)
    exW = D.Experiment([(-1, 0.5 * L0), (0, 0.0)], L0=0.5 * L0)
    m1, _ = D.simulate(D.params(ex=exW), exW, tt, sens=False)
    # a two-pool washout: m exchanges with a slow pool at rate kx
    def twopool(kx, frac):
        th, ke, K = C.TH_5, C.KE, C.KM0
        m, q = m1[0], frac * m1[0]
        out = [m]
        for i in range(1, tt.size):
            dt = tt[i] - tt[i - 1]
            f = exW.u(tt[i]) + 0.5 * L0 - th * m / (K + m) - ke * m - kx * (m - q)
            m = max(m + dt * f, 1e-6); q = max(q + dt * kx * (m - q), 1e-6)
            out.append(m)
        return np.array(out)
    for kx in (0.5, 2.0):
        m2 = twopool(kx, 1.0)
        # fit a single-compartment washout by scaling the MM prediction's shape is
        # not needed: report the extra curvature as the gap in log-marker decay
        d1 = np.diff(np.log(m1)); d2 = np.diff(np.log(m2))
        print(f"      exchange rate kx = {kx}/d: peak |d ln m| gap single vs two-pool"
              f" = {np.max(np.abs(d2 - d1)) / np.max(np.abs(d1)) * 100:.0f}% of the decay rate"
              f" (0% at steady state; grows in the transient)")

    print("\n  VERDICT (power of residual diagnostics, by mechanism and design):")
    print("   - The residual-vs-MARKER test is degenerate here: the MM fit is itself")
    print("     a marker model, so it removes marker-linear structure by construction")
    print("     (0% even under a real misspecification). The vs-TIME and vs-LOAD")
    print("     tests carry the power.")
    print("   - SECOND COMPARTMENT: the monthly record is BLIND (identical steady-state")
    print("     balance); it shows only in a transient's second decay phase.")
    print("   - SATURABLE non-enzymatic clearance: the monthly record is nearly BLIND")
    print("     (narrow on-therapy m-range looks linear; power ~5% at n=25). A load")
    print("     EXCURSION carries the signal, BUT the residual-slope test is size-invalid")
    print("     at large N (EIV inflates the true-MM false positive to ~100%). At a")
    print("     size-valid N only STRONG saturation is genuinely detected; MILD")
    print("     saturation (Kne=150) is not. The excursion is necessary, not sufficient.")
    print("   - LOAD MISMATCH (a time-varying drift): the monthly record CAN catch it,")
    print("     via residual-vs-time (88-100% power at a 10-20% drift). The constant")
    print("     bias b absorbs only a constant offset, not a trend.")
    print("   - So the load-record assumption is partly checkable now; a second")
    print("     compartment and saturable clearance are not, short of a transient or")
    print("     excursion. And with NO real data, this remains a power analysis, not")
    print("     a test of the specification.")


if __name__ == "__main__":
    # fix a baseline marker for the saturable-NE matching point
    globals()["M0"] = float(C.record(C.TH_5, seed=0)["m"].mean())
    main()
