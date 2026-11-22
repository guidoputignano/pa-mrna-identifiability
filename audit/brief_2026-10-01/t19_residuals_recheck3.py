"""
Independent adversarial recheck of ITEM 10 (out/t19_residuals.txt).

Imports ONLY core/dyn primitives, NOT t19_residuals.py. Every piece of the
pipeline -- the MM MAP fit, the residual-slope t-tests, the misspecified
data-generating processes, and the two-pool washout -- is re-derived here from
the stated model, so this is a from-scratch reconstruction, not a re-run.

Model fitted (flux / recorded-load ordinate):
    pred_i = (theta + inc_i) m_i/(K + m_i) + k_e m_i + b
MAP objective: sum_i ((load_i - pred_i)/sg)^2 + (b/sb)^2 + ((ln ke - ln KE)/0.8)^2
    sg = 0.05 * L0     (assay CV mapped to flux at anchor load)
    sb = LOAD_FRAC*L0  (load-bias prior)
Residual (flux units) = load - pred.
Test: OLS slope t of residual on a centred regressor; reject at |t|>1.96.
"""
import numpy as np
from scipy.optimize import least_squares, brentq

import core as C
import dyn as D

RNG = np.random.default_rng(7)             # my own seed, deliberately different

TH, K0, KE = C.TH_5, C.KM0, C.KE
L0 = float(C.record(TH, seed=0)["L"])
M0 = float(C.record(TH, seed=0)["m"].mean())
SG = 0.05 * L0
SB = C.LOAD_FRAC * L0


# ---------------------------------------------------------- my own MM MAP fit
def mm_fit(loads, inc, m_obs):
    """Return flux-unit residuals load - pred at the MAP estimate."""
    def res(x):
        lt, lk, b, lke = x
        th, K, ke = np.exp(lt), np.exp(lk), np.exp(lke)
        pred = (th + inc) * m_obs / (K + m_obs) + ke * m_obs + b
        return np.concatenate([(loads - pred) / SG,
                               [b / SB, (lke - np.log(KE)) / C.KE_PRIOR]])
    x0 = np.array([np.log(TH), np.log(K0), 0.0, np.log(KE)])
    s = least_squares(res, x0, method="lm", xtol=1e-12, ftol=1e-12, max_nfev=5000)
    lt, lk, b, lke = s.x
    th, K, ke = np.exp(lt), np.exp(lk), np.exp(lke)
    pred = (th + inc) * m_obs / (K + m_obs) + ke * m_obs + b
    return loads - pred


def slope_t(x, r):
    """t-stat of the OLS slope of r on centred x (my own derivation)."""
    x = x - x.mean()
    sxx = x @ x
    if sxx <= 0:
        return 0.0
    b = (x @ r) / sxx
    dof = max(x.size - 2, 1)
    s2 = ((r - b * x) ** 2).sum() / dof
    se = np.sqrt(s2 / sxx)
    return 0.0 if se == 0 else b / se


# -------------------------------------------------- misspecified steady states
def m_mm(load, inc):
    return float(C.m_ss(TH, K0, KE, load, inc))


def m_sat(load, inc, Kne, m_ref):
    """TRUE balance with saturable non-enzymatic clearance, matched at m_ref.
       (theta+inc) m/(K+m) + Vne m/(Kne+m) = load, Vne = KE*(Kne+m_ref)."""
    Vne = KE * (Kne + m_ref)
    f = lambda m: (TH + inc) * m / (K0 + m) + Vne * m / (Kne + m) - load
    # bracket: f(0)=-load<0, f grows unboundedly
    hi = max(load, 1.0)
    while f(hi) < 0:
        hi *= 2.0
    return brentq(f, 1e-9, hi, xtol=1e-12, rtol=1e-12)


# ----------------------------------------------------------------- power loops
def monthly_power(dgp_m, nmc=300):
    """dgp_m(load, inc) -> true marker. Return (rej_rates, mean_t) for the 3 tests."""
    rej = np.zeros(3)
    sig = np.zeros(3)
    for j in range(nmc):
        r = C.record(TH, seed=j % 24)
        loads, inc, t = r["loads"], r["inc"], r["times"]
        mt = np.array([dgp_m(loads[i], inc[i]) for i in range(loads.size)])
        m_obs = mt * np.exp(RNG.normal(0, 0.05, mt.size))
        res = mm_fit(loads, inc, m_obs)
        ts = np.array([slope_t(m_obs, res), slope_t(t, res), slope_t(loads, res)])
        rej += np.abs(ts) > 1.96
        sig += ts
    return rej / nmc, sig / nmc


def excursion_power(dgp_m, nmc=120):
    """Dense within-record load sweep 0.6L -> 1.5L, 5-min over 2 d, inc=0."""
    tt = np.arange(0.0, 2.0 + 1e-9, 5 / 60 / 24)
    loads = L0 * (0.6 + 0.9 * tt / tt[-1])
    inc = np.zeros_like(tt)
    mt = np.array([dgp_m(loads[i], 0.0) for i in range(loads.size)])
    rej = np.zeros(3)
    for j in range(nmc):
        m_obs = mt * np.exp(RNG.normal(0, 0.05, mt.size))
        res = mm_fit(loads, inc, m_obs)
        ts = np.array([slope_t(m_obs, res), slope_t(tt, res), slope_t(loads, res)])
        rej += np.abs(ts) > 1.96
    return rej / nmc


def drift_power(dr, nmc=300):
    """True load = recorded load * (1 + dr*(t/T - 0.5)): zero-mean linear drift."""
    rej = np.zeros(3)
    for j in range(nmc):
        r = C.record(TH, seed=j % 24)
        loads, inc, t = r["loads"], r["inc"], r["times"]
        Ltrue = loads * (1 + dr * (t / t[-1] - 0.5))
        mt = np.array([m_mm(Ltrue[i], inc[i]) for i in range(loads.size)])
        m_obs = mt * np.exp(RNG.normal(0, 0.05, mt.size))
        res = mm_fit(loads, inc, m_obs)
        ts = np.array([slope_t(m_obs, res), slope_t(t, res), slope_t(loads, res)])
        rej += np.abs(ts) > 1.96
    return rej / nmc


# ------------------------------------------------------- second compartment
def twopool_gap(kx, frac=1.0):
    """Peak |d ln m| gap (single vs two-pool washout) as % of max single-pool rate.
       Independent integrator (RK4, finer grid) for the exchange pool."""
    tt = np.arange(0.0, 2.0 + 1e-9, 5 / 60 / 24)
    exW = D.Experiment([(-1, 0.5 * L0), (0, 0.0)], L0=0.5 * L0)
    m1, _ = D.simulate(D.params(ex=exW), exW, tt, sens=False)

    def rhs(m, q, t):
        u = exW.u(t)
        dm = u + 0.5 * L0 - TH * m / (K0 + m) - KE * m - kx * (m - q)
        dq = kx * (m - q)
        return dm, dq

    # RK4 on a 10x finer grid, sampled back to tt
    nsub = 10
    m = m1[0]; q = frac * m1[0]
    out = [m]
    for i in range(1, tt.size):
        h = (tt[i] - tt[i - 1]) / nsub
        tcur = tt[i - 1]
        for _ in range(nsub):
            k1m, k1q = rhs(m, q, tcur)
            k2m, k2q = rhs(m + 0.5 * h * k1m, q + 0.5 * h * k1q, tcur + 0.5 * h)
            k3m, k3q = rhs(m + 0.5 * h * k2m, q + 0.5 * h * k2q, tcur + 0.5 * h)
            k4m, k4q = rhs(m + h * k3m, q + h * k3q, tcur + h)
            m += h * (k1m + 2 * k2m + 2 * k3m + k4m) / 6
            q += h * (k1q + 2 * k2q + 2 * k3q + k4q) / 6
            m = max(m, 1e-9); q = max(q, 1e-9)
            tcur += h
        out.append(m)
    m2 = np.array(out)
    d1 = np.diff(np.log(m1)); d2 = np.diff(np.log(m2))
    gap = np.max(np.abs(d2 - d1)) / np.max(np.abs(d1)) * 100
    ss_gap = abs(d2[-1] - d1[-1]) / np.max(np.abs(d1)) * 100  # near steady state
    return gap, ss_gap


def main():
    print("=" * 84)
    print("  RECHECK3 of ITEM 10 -- residual-diagnostic power (independent rebuild)")
    print("=" * 84)
    print(f"  L0={L0:.3f}  M0(m_ref)={M0:.3f}  SG={SG:.3f}  SB={SB:.3f}")

    print("\n(a) SANITY: true MM data -> residuals white, reject ~ nominal 5%")
    p, s = monthly_power(lambda L, inc: m_mm(L, inc), nmc=400)
    print(f"    monthly reject @5%: vs marker {100*p[0]:.1f}%, vs time {100*p[1]:.1f}%, "
          f"vs load {100*p[2]:.1f}%  (mean t: mk {s[0]:+.2f} tm {s[1]:+.2f} ld {s[2]:+.2f})")

    print("\n(b) SATURABLE NON-ENZYMATIC CLEARANCE: monthly (blind) vs excursion (detects)")
    for Kne in (150.0, 60.0):
        pm, sm = monthly_power(lambda L, inc, Kne=Kne: m_sat(L, inc, Kne, M0), nmc=300)
        pe = excursion_power(lambda L, inc, Kne=Kne: m_sat(L, inc, Kne, M0), nmc=120)
        print(f"    Kne={Kne:.0f}:")
        print(f"      monthly:   vs marker {100*pm[0]:.0f}%, vs time {100*pm[1]:.0f}%, "
              f"vs load {100*pm[2]:.0f}%  (mean residual-vs-marker t={sm[0]:+.2f})")
        print(f"      excursion: vs marker {100*pe[0]:.0f}%, vs time {100*pe[1]:.0f}%, "
              f"vs load {100*pe[2]:.0f}%")

    print("\n(c) LOAD DRIFT (zero-mean, b cannot absorb): monthly residual-vs-time power")
    for dr in (0.10, 0.20):
        rej = drift_power(dr, nmc=300)
        print(f"    drift {int(dr*100)}%: vs marker {100*rej[0]:.0f}%, "
              f"vs time {100*rej[1]:.0f}%, vs load {100*rej[2]:.0f}%")

    print("\n(d) SECOND COMPARTMENT: invisible at steady state, grows in transient")
    for kx in (0.5, 2.0):
        gap, ss = twopool_gap(kx)
        print(f"    kx={kx}/d: peak |d ln m| gap = {gap:.0f}% of decay rate; "
              f"near-steady-state gap = {ss:.1f}%")


if __name__ == "__main__":
    main()
