"""
The minimal defensible footing, and the dose-precision requirement on it.

audit/d11_priorfree.py shows that of the four priors the filter carries, only
the one on ln k_e does any work: dropping the ln K_m prior changes
sd(ln theta) by nothing at all (0.1474 either way, and still 0.1474 at a K_m
prior of sd 3.0, a factor of twenty), dropping the ln theta prior changes it to
0.1499, dropping the load-bias prior to 0.1483, and dropping the ln k_e prior
to 0.6451.

So the paper should claim exactly one prior -- k_e within a factor of about two,
which the literature supports -- and let theta and K_m float. That also settles
the objection rem:no_km would otherwise raise. Every number below is on that
footing.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
from scipy.optimize import brentq
import design.run as D, model_A.run as LOC

P = D.P
L0 = D.load(0.0)
SIG = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)
GATE = 0.20


def record(t, load_cv, seed=3):
    rng = np.random.default_rng(seed)
    dd = np.asarray(D.delta(t), dtype=float)
    Ls = (L0 * (1.0 + load_cv * rng.standard_normal(t.size)) if load_cv
          else np.full(t.size, L0))
    B = P.ke * P.Km + P.theta + dd - Ls
    lv = (-B + np.sqrt(B ** 2 + 4.0 * P.ke * Ls * P.Km)) / (2.0 * P.ke)
    return lv, dd


def sd(lv, dd, ke_prior=0.8, kap=None, load_frac=D.LOAD_FRAC, carry=True):
    """theta and K_m FREE. Priors only on ln k_e and the load bias."""
    pri = [np.inf, np.inf, ke_prior ** 2, (load_frac * L0) ** 2]
    if kap is not None:
        pri.append(kap ** 2)
    F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri])
    for m, d in zip(lv, dd):
        g = m / (P.Km + m)
        eff = P.theta + (d if carry else 0.0)
        q = [P.theta * g, -eff * P.Km * m / (P.Km + m) ** 2, P.ke * m, 1.0]
        if kap is not None:
            q.append(d * g)
        q = np.array(q)
        F = F + np.outer(q, q) / SIG ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def main():
    print(__doc__)
    designs = {
        "six-monthly (n=4)": np.arange(180.0, 731.0, 180.0),
        "one dose cycle (n=7)": np.arange(0.0, D.DOSE_INTERVAL, 2.0),
        "monthly x 2y (n=24)": np.arange(0.0, 730.0, 30.4),
        "fortnightly x 2y (n=52)": np.arange(0.0, 730.0, 14.0),
        "weekly x 2y (n=105)": np.arange(0.0, 730.0, 7.0),
        "two-daily x 2y (n=365)": np.arange(0.0, 730.0, 2.0),
    }
    print("=" * 88)
    print("THE DESIGN TABLE on the minimal footing (only ln k_e and the load are prior)")
    print("=" * 88)
    print(f"{'design':>24} | {'Delta dropped':>14} {'Delta carried':>14} | "
          f"{'ratio':>7} {'gate':>6}")
    for nm, t in designs.items():
        lv, dd = record(t, 0.10)
        a, b = sd(lv, dd, carry=False), sd(lv, dd, carry=True)
        f = lambda v: "singular" if not np.isfinite(v) else f"{v:.4f}"
        r = a / b if np.isfinite(a) else np.inf
        print(f"{nm:>24} | {f(a):>14} {f(b):>14} | {r:7.1f} "
              f"{'PASS' if b < GATE else 'fail':>6}")

    print()
    print("=" * 88)
    print("HOW WELL THE DELIVERED INCREMENT HAS TO BE KNOWN, on this footing")
    print("=" * 88)
    print(f"{'design':>24} | " + " ".join(f"{c:>9}" for c in
          ("exact", "5%", "10%", "16%", "30%")) + f" | {'gate needs':>11}")
    for nm, t in designs.items():
        lv, dd = record(t, 0.10)
        cells = [sd(lv, dd, kap=k) for k in (None, 0.05, 0.10, 0.159, 0.30)]
        try:
            crit = brentq(lambda x: sd(lv, dd, kap=x) - GATE, 1e-3, 3.0, xtol=1e-5)
            lab = f"{crit:.1%}"
        except ValueError:
            lab = "any" if cells[-1] < GATE else "unreachable"
        print(f"{nm:>24} | " + " ".join(f"{c:9.4f}" for c in cells) + f" | {lab:>11}")

    print()
    print("=" * 88)
    print("AND AGAINST HOW WELL k_e IS KNOWN  (two-daily, load CV 10%, Delta exact)")
    print("=" * 88)
    t = designs["two-daily x 2y (n=365)"]
    lv, dd = record(t, 0.10)
    for kep in (0.1, 0.3, 0.5, 0.8, 1.2, 2.0, 1e3):
        v = sd(lv, dd, ke_prior=kep)
        lab = "flat" if kep >= 1e3 else f"{kep:.1f} (factor {np.exp(kep):.1f})"
        print(f"    ln k_e prior sd {lab:>22}:  sd(ln theta) = {v:8.4f}"
              f"  {'PASS' if v < GATE else 'fail'}")
    try:
        c = brentq(lambda x: sd(lv, dd, ke_prior=x) - GATE, 1e-3, 3.0, xtol=1e-5)
        print(f"\n    the gate needs k_e to within a log-sd of {c:.3f}, "
              f"a factor of {np.exp(c):.2f}")
    except ValueError:
        print("\n    no k_e precision reaches the gate")


if __name__ == "__main__":
    main()
