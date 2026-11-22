"""
Adversarial independent recheck of ITEM 9 (out/t13_conditioning.txt).

Imports ONLY core primitives (and numpy/scipy). It does NOT import
t13_conditioning. Every number is rebuilt from first principles:

  * the per-draw loads are reconstructed from the stated RNG recipe,
  * each marker m(t_k) is re-solved with an INDEPENDENT bracket root finder
    (scipy brentq) on the mass balance, not with core.m_ss/level,
  * the sensitivity rows, the per-draw flux sigma, the priors and the
    information matrix are all re-written here and the sd(ln theta) is
    recomputed by hand.

core.record is exercised as the object under test for clause (a): we check
that what it returns is a per-draw load array and a marker that solves the
balance at each draw's own (load_k, inc_k).
"""
import numpy as np
from scipy.optimize import brentq

import core as C


def balance_resid(m, theta, Km, ke, inc, L):
    """f(m) = (theta + inc) m/(Km+m) + ke m - L ; the mass balance."""
    return (theta + inc) * m / (Km + m) + ke * m - L


def solve_balance_independent(theta, Km, ke, loads, inc):
    """Solve the balance per draw with brentq; shares no algebra with core."""
    loads = np.atleast_1d(np.asarray(loads, float))
    inc = np.broadcast_to(np.asarray(inc, float), loads.shape)
    out = np.empty_like(loads)
    for k in range(loads.size):
        # m is positive; upper bracket from ke*m <= L -> m <= L/ke, pad it
        hi = loads[k] / ke + loads[k] / max(theta + inc[k], 1e-12) + 1.0
        out[k] = brentq(balance_resid, 1e-12, hi,
                        args=(theta, Km, ke, inc[k], loads[k]),
                        xtol=1e-14, rtol=1e-15, maxiter=500)
    return out


def rebuild_loads(L, seed, cv, n):
    """Reconstruct record()'s lognormal per-draw load array from scratch."""
    if not cv:
        return np.full(n, L)
    rng = np.random.default_rng(seed)
    s2 = np.log1p(cv ** 2)
    return L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), n))


def flux_rows_independent(theta, Km, ke, m, inc):
    """d(clearance)/d(ln theta, ln K_m, b, ln k_e), rebuilt from the formulas."""
    g = m / (Km + m)
    col_lntheta = theta * g                       # theta * dC/dtheta
    col_lnKm = -(theta + inc) * Km * m / (Km + m) ** 2
    col_b = np.ones_like(m)
    col_lnke = ke * m                             # ke * dC/dke
    return np.stack([col_lntheta, col_lnKm, col_b, col_lnke], axis=1)


def sigma_perdraw_independent(theta, Km, ke, m, inc, L,
                              load_frac=0.05, assay_cv=0.05):
    """Per-draw flux sigma: quadrature of 5% load scatter and 5% assay CV."""
    slope = (theta + inc) * Km / (Km + m) ** 2 + ke     # C'(m)
    return np.hypot(load_frac * L, slope * assay_cv * m)


def fim_independent(theta, Km, ke, m, inc, L,
                    b_prior=None, ke_prior=0.8,
                    load_frac=0.05, assay_cv=0.05):
    """Information matrix in (ln theta, ln K_m, b, ln k_e), rebuilt by hand."""
    Q = flux_rows_independent(theta, Km, ke, m, inc)
    sg = sigma_perdraw_independent(theta, Km, ke, m, inc, L, load_frac, assay_cv)
    F = (Q / sg[:, None] ** 2).T @ Q
    sb = load_frac * L if b_prior is None else b_prior
    pri = np.zeros(4)
    pri[2] = sb ** -2            # bias b prior
    pri[3] = ke_prior ** -2      # ln k_e prior
    return F + np.diag(pri)


def sd_lntheta(F):
    return float(np.sqrt(np.linalg.inv(F)[0, 0]))


def main():
    out = []
    p = out.append
    p("=" * 92)
    p("  ITEM 9 INDEPENDENT RECHECK  (rebuilt from first principles)")
    p("=" * 92)

    theta, Km, ke = C.TH_5, C.KM0, C.KE
    r = C.record(theta, seed=0)                 # the object under test
    L = r["L"]
    n = r["loads"].size
    p(f"  record: theta={theta:g}  L={L:.4f}  draws={n}")
    p(f"  analytic anchor-A L = 3000*10/110 + 5*10 = "
      f"{3000*10/110 + 5*10:.4f}")
    p("")

    # ---- (a) per-draw load array + marker solves balance at each draw -------
    p("  (a) PER-DRAW LOAD AND MARKER SOLVE THE BALANCE")
    loads = r["loads"]
    inc = r["inc"]
    p(f"      loads is an array shape {loads.shape}, varies: "
      f"min {loads.min():.2f}  median {np.median(loads):.2f}  "
      f"max {loads.max():.2f}  sd {loads.std():.3f} ({100*loads.std()/L:.2f}% of L)")
    # reconstruct loads from the RNG recipe, fully independent
    loads_rebuilt = rebuild_loads(L, 0, 0.10, n)
    p(f"      loads reconstructed from RNG recipe matches record: "
      f"max|diff| = {np.abs(loads_rebuilt - loads).max():.2e}")
    p(f"      marker m varies draw to draw: min {r['m'].min():.4f}  "
      f"median {np.median(r['m']):.4f}  max {r['m'].max():.4f}  "
      f"(unique={np.unique(np.round(r['m'],6)).size}/{n})")
    # independent brentq solve of the balance at each draw's own load & inc
    m_indep = solve_balance_independent(theta, Km, ke, loads, inc)
    p(f"      brentq re-solve of (theta+inc)g(m)+ke*m=load_k vs record m: "
      f"max|diff| = {np.abs(m_indep - r['m']).max():.2e}")
    # does it instead solve at ONE constant L? test the alternative hypothesis
    m_if_constL = solve_balance_independent(theta, Km, ke, np.full(n, L), inc)
    p(f"      [skeptic] if it used ONE constant L, m would differ from record "
      f"by up to {np.abs(m_if_constL - r['m']).max():.4f} (so it does NOT)")
    # confirm the increment is per-draw too
    p(f"      increment Delta(t_k) is per-draw: min {inc.min():.2f}  "
      f"median {np.median(inc):.2f}  max {inc.max():.2f}")
    p("")

    # ---- (b) bias b is a single constant column, prior 5% of L -------------
    p("  (b) THE BIAS b IS A SINGLE CONSTANT COLUMN, PRIOR 5% OF L")
    Q_core = C.flux_rows(r)
    Q_mine = flux_rows_independent(theta, Km, ke, r["m"], inc)
    p(f"      core.flux_rows col 2 unique values: "
      f"{np.unique(np.round(Q_core[:,2], 12))}  (all ones = "
      f"{np.allclose(Q_core[:,2], 1.0)})")
    p(f"      my rebuilt flux_rows matches core to max|diff| = "
      f"{np.abs(Q_core - Q_mine).max():.2e}")
    p(f"      cols 0,1,3 DO vary (they are per-draw): "
      f"n_unique = {[np.unique(np.round(Q_core[:,j],6)).size for j in (0,1,3)]}")
    p(f"      bias prior sd = LOAD_FRAC * L = 0.05 * {L:.4f} = "
      f"{0.05*L:.4f} flux units  (core.LOAD_FRAC={C.LOAD_FRAC})")
    p("")

    # ---- (c) the three sd(ln theta) numbers --------------------------------
    p("  (c) sd(ln theta): per-draw record, load-flat, bias-pinned")
    # base: per-draw record
    F_base = fim_independent(theta, Km, ke, r["m"], inc, L)
    base = sd_lntheta(F_base)
    # flat: load_cv = 0 -> loads all L, re-solve markers independently
    loads_flat = rebuild_loads(L, 0, 0.0, n)
    m_flat = solve_balance_independent(theta, Km, ke, loads_flat, inc)
    F_flat = fim_independent(theta, Km, ke, m_flat, inc, L)
    flat = sd_lntheta(F_flat)
    # tight: bias pinned, b_prior = 1e-6 * L
    F_tight = fim_independent(theta, Km, ke, r["m"], inc, L, b_prior=1e-6 * L)
    tight = sd_lntheta(F_tight)

    p(f"      per-draw record         : {base:.4f}   (claim ~3.06)")
    p(f"      load held flat (cv=0)   : {flat:.4f}   (claim ~3.14)")
    p(f"      bias b pinned (L known) : {tight:.4f}   (claim ~3.03)")
    p("")
    p(f"      ordering check: flat ({flat:.3f}) > base ({base:.3f}) > "
      f"tight ({tight:.3f})  -> {flat > base > tight}")
    p(f"      claim 'flat slightly WORSE than base': "
      f"{flat > base} (delta = +{flat-base:.4f})")

    # cross-check against core.fim / core.sd_th (should agree exactly)
    p("")
    p("  CROSS-CHECK vs core.fim/core.sd_th (must match my hand-built FIM):")
    r_flat = C.record(theta, seed=0, load_cv=0.0)
    core_base = C.sd_th(C.fim(r, noise="perdraw"))
    core_flat = C.sd_th(C.fim(r_flat, noise="perdraw"))
    core_tight = C.sd_th(C.fim(r, noise="perdraw", b_prior=1e-6 * L))
    p(f"      base : mine {base:.6f}  core {core_base:.6f}  "
      f"|diff| {abs(base-core_base):.2e}")
    p(f"      flat : mine {flat:.6f}  core {core_flat:.6f}  "
      f"|diff| {abs(flat-core_flat):.2e}")
    p(f"      tight: mine {tight:.6f}  core {core_tight:.6f}  "
      f"|diff| {abs(tight-core_tight):.2e}")

    # round to 2 decimals as the task asks
    p("")
    p("  ROUNDED TO 2 DECIMALS:")
    p(f"      base={base:.2f}  flat={flat:.2f}  tight={tight:.2f}")
    p(f"      claim base=3.06 flat=3.14 tight=3.03")

    text = "\n".join(out)
    print(text)
    return dict(base=base, flat=flat, tight=tight)


if __name__ == "__main__":
    main()
