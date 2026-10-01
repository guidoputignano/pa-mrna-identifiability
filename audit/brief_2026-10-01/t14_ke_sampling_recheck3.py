"""
Independent re-derivation of ITEM 12 (out/t14_ke_sampling.txt).

Does NOT import t14_ke_sampling. Uses only core primitives to build the
record (C.record), then rebuilds the (ln theta, ln K_m, b, ln k_e) FIM from
scratch with my own flux-row algebra, my own per-draw sigma, my own prior
block, and my own inversion. Cross-checks against C.fim / C.sd_th as a
separate control.

record: n draws on np.linspace(0,730,n), theta = 5% of normal, dose train on,
intake load_cv 0.10, median of sd(ln theta) over the 24 intake seeds, with a
ln-k_e prior of sd = ke_sd (and the load-bias prior at 5% of L, per core).
"""
import numpy as np

import core as C

HORIZON = 730.0
SEEDS = list(C.SEEDS)                       # 0..23  -> 24 seeds
LF, AC = C.LOAD_FRAC, C.ASSAY_CV            # 0.05, 0.05


# --------------------------------------------------- independent FIM builder
def my_fim(r, ke_sd):
    """
    Rebuild the information matrix in (ln theta, ln K_m, b, ln k_e) from the
    record fields, independently of core.fim. 'perdraw' noise.
    """
    th, K, ke = r["theta"], r["Km"], r["ke"]
    m, inc, L = r["m"], r["inc"], r["L"]
    # flux sensitivities per draw (my own transcription of the algebra)
    q_th = th * m / (K + m)
    q_km = -(th + inc) * K * m / (K + m) ** 2
    q_b = np.ones_like(m)
    q_ke = ke * m
    Q = np.stack([q_th, q_km, q_b, q_ke], axis=1)
    # per-draw clearance slope and perdraw flux sigma
    slope = (th + inc) * K / (K + m) ** 2 + ke
    sg = np.hypot(LF * L, slope * AC * m)
    F = (Q / sg[:, None] ** 2).T @ Q
    # priors: load bias at 5% of L (index 2), ln k_e at ke_sd (index 3)
    F[2, 2] += (LF * L) ** -2
    F[3, 3] += ke_sd ** -2
    return F


def my_sd_th(F):
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    return float(np.sqrt(np.linalg.inv(F)[0, 0]))


def record_n(n, seed):
    times = np.linspace(0.0, HORIZON, int(n))
    return C.record(C.TH_5, seed=seed, times=times)


def sd_at(ke_sd, n, use_core=False):
    """median over 24 seeds of sd(ln theta)."""
    vals = []
    for s in SEEDS:
        r = record_n(n, s)
        if use_core:
            v = C.sd_th(C.fim(r, noise="perdraw", ke_prior=ke_sd))
        else:
            v = my_sd_th(my_fim(r, ke_sd))
        vals.append(v)
    vals = np.array(vals, float)
    vals = vals[np.isfinite(vals)]
    return float(np.median(vals))


def crossing(ke_sd, gate, ns, vals):
    """smallest n at which median sd <= gate, log-log interpolated; else None."""
    below = np.where(vals <= gate)[0]
    if not below.size:
        return None
    i = below[0]
    if i == 0:
        return float(ns[0])
    x0, x1 = np.log(ns[i - 1]), np.log(ns[i])
    y0, y1 = np.log(vals[i - 1]), np.log(vals[i])
    xc = x0 + (np.log(gate) - y0) * (x1 - x0) / (y1 - y0)
    return float(np.exp(xc))


def main():
    ke_sds = [1e-4, 0.05, 0.10, 0.20, 0.40, 0.80]
    ns = np.unique(np.round(np.geomspace(24, 1460, 18)).astype(int))
    print("grid n (2yr):", list(ns))
    print("draws/yr     :", [round(n / 2) for n in ns])
    print()

    # full table, independent builder
    table = {}
    print("sd(ln theta) median/24 seeds  [INDEPENDENT FIM]")
    print("ke_sd \\ n   " + "".join(f"{n:>7d}" for n in ns))
    for ke_sd in ke_sds:
        row = np.array([sd_at(ke_sd, n) for n in ns])
        table[ke_sd] = row
        lab = "exact" if ke_sd < 1e-3 else f"{ke_sd:.2f}"
        print(f"{lab:>9}  " + "".join(f"{v:7.3f}" for v in row))
    print()

    # control: core.fim path, a few cells, to confirm my builder matches core
    print("CONTROL  (my_fim vs core.fim), max |diff| over spot cells:")
    diffs = []
    for ke_sd in ke_sds:
        for n in (ns[0], ns[len(ns) // 2], ns[-1]):
            a = sd_at(ke_sd, n)
            b = sd_at(ke_sd, n, use_core=True)
            diffs.append(abs(a - b))
    print(f"   max|my - core| = {max(diffs):.2e}")
    print()

    # (a),(b) crossings
    print("CROSSINGS (n over 2yr; /yr in brackets):")
    for gate in (0.30, 0.20):
        print(f"  gate {gate:.2f}:")
        for ke_sd in ke_sds:
            nc = crossing(ke_sd, gate, ns, table[ke_sd])
            lab = "exact" if ke_sd < 1e-3 else f"{ke_sd:.2f}"
            if nc is None:
                print(f"    {lab:>6}: not reached by {ns[-1]} (min sd over grid "
                      f"= {table[ke_sd].min():.3f} at n={ns[np.argmin(table[ke_sd])]})")
            else:
                print(f"    {lab:>6}: n = {nc:7.1f} ({nc/2:6.1f}/yr)")
    print()

    # (c) achievable precision row, min/max across 24..daily(=1460)
    r020 = table[0.20]
    print(f"(c) ke_sd=0.20 row: min={r020.min():.3f} (n={ns[np.argmin(r020)]}), "
          f"max={r020.max():.3f} (n={ns[np.argmax(r020)]})")
    # value at 'daily' ~ 365/yr = 730 over 2yr, and feasible ceiling 182/yr=364
    print(f"    ke_sd=0.20 at n=730 (daily): {sd_at(0.20, 730):.3f}")
    print(f"    ke_sd=0.20 at n=364 (182/yr ceiling): {sd_at(0.20, 364):.3f}")
    print()

    # (d) feasibility box: exact n = 2*rate
    print("FEASIBILITY BOX (n = 2*rate, exact):")
    for rate, lab in ((26, "biweekly 26/yr"), (52, "weekly 52/yr"),
                      (182, "every-other-day 182/yr")):
        print(f"   ke_sd=0.20 {lab:>26}: sd = {sd_at(0.20, 2 * rate):.3f}")
    print(f"   ke_sd=exact every-other-day 182/yr: sd = {sd_at(1e-4, 2 * 182):.3f}")
    print()

    # (d) best feasible point: min over feasible box (ke>=0.20, draws<=ceiling)
    #     and the absolute floor of the ke=0.20 row (daily, 2x ceiling)
    print("(d) FEASIBLE-REGION FLOOR (achievable ke_sd>=0.20):")
    print(f"    best feasible (ke=0.20, 182/yr ceiling) sd = {sd_at(0.20, 364):.3f}")
    print(f"    even ke=0.20 at 1460 draws (730/yr) sd = {r020[-1]:.3f}  (> 0.20? "
          f"{r020[-1] > 0.20})")
    print(f"    -> 0.20 contour reached by any feasible point? "
          f"{'NO' if sd_at(0.20, 364) > 0.20 else 'YES'}")
    print()

    # spot checks to 2 decimals against the published table
    print("SPOT CHECKS vs published table (2 dp):")
    published = {
        (1e-4, 24): 0.408, (1e-4, 1460): 0.112, (1e-4, 309): None,
        (0.05, 24): 0.463, (0.05, 1460): 0.252,
        (0.10, 24): 0.599, (0.10, 707): 0.455, (0.10, 1460): 0.416,
        (0.20, 24): 0.963, (0.20, 1460): 0.598,
        (0.40, 436): 1.096, (0.80, 24): 3.096, (0.80, 1460): 0.743,
    }
    for (ke_sd, n), pub in published.items():
        if pub is None:
            continue
        mine = sd_at(ke_sd, n)
        ok = abs(round(mine, 2) - round(pub, 2)) <= 0.01
        lab = "exact" if ke_sd < 1e-3 else f"{ke_sd:.2f}"
        print(f"   ke={lab:>5} n={n:>5}: mine={mine:.3f} pub={pub:.3f}  "
              f"{'OK' if ok else 'MISMATCH'}")


if __name__ == "__main__":
    main()
