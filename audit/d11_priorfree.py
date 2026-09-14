"""
How much of the gate is the data and how much is the prior?

The information filter of Section 7 carries log-priors of 0.8 on ln theta,
ln K_m and ln k_e, and 5% on the load bias. A log-prior of 0.8 is a factor of
2.2, and Remark rem:no_km forbids importing K_m from in vitro work -- so a
number that depends on assuming K_m within a factor of two is a number this
note is not entitled to. This strips the priors off and asks what the record
alone delivers, on each design and against day-to-day variation in intake.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
import design.run as D, model_A.run as LOC

P = D.P
L0 = D.load(0.0)
SIG = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)


def record(t, load_cv, seed=3):
    rng = np.random.default_rng(seed)
    dd = np.asarray(D.delta(t), dtype=float)
    Ls = (L0 * (1.0 + load_cv * rng.standard_normal(t.size)) if load_cv
          else np.full(t.size, L0))
    B = P.ke * P.Km + P.theta + dd - Ls
    lv = (-B + np.sqrt(B ** 2 + 4.0 * P.ke * Ls * P.Km)) / (2.0 * P.ke)
    return lv, dd


def sd(lv, dd, priors, carry=True):
    """4 parameters throughout: ln theta, ln K_m, ln k_e, load bias."""
    F = (np.linalg.inv(np.diag([0.8 ** 2, 0.8 ** 2, 0.8 ** 2, (D.LOAD_FRAC * L0) ** 2]))
         if priors else np.zeros((4, 4)))
    for m, d in zip(lv, dd):
        eff = P.theta + (d if carry else 0.0)
        q = np.array([P.theta * m / (P.Km + m),
                      -eff * P.Km * m / (P.Km + m) ** 2, P.ke * m, 1.0])
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
        "weekly x 2y (n=105)": np.arange(0.0, 730.0, 7.0),
        "two-daily x 2y (n=365)": np.arange(0.0, 730.0, 2.0),
    }
    print("=" * 92)
    print("PRIOR-FREE sd(ln theta), increment carried, four parameters")
    print("=" * 92)
    print(f"{'design':>24} | " + " ".join(f"{c:>10}" for c in
          ("CV 0%", "CV 2%", "CV 5%", "CV 10%", "CV 20%")) + f" | {'with priors':>12}")
    for name, t in designs.items():
        cells = []
        for cv in (0.0, 0.02, 0.05, 0.10, 0.20):
            lv, dd = record(t, cv)
            v = sd(lv, dd, priors=False)
            cells.append("singular" if not np.isfinite(v) else f"{v:.4f}")
        lv, dd = record(t, 0.10)
        wp = sd(lv, dd, priors=True)
        print(f"{name:>24} | " + " ".join(f"{c:>10}" for c in cells) + f" | {wp:12.4f}")

    print()
    print("  A design is worth reporting only where the prior-free column is finite and")
    print("  small: that is where the RECORD identifies theta rather than the prior.")

    print()
    print("=" * 92)
    print("WHAT THE INCREMENT IS WORTH, PRIOR-FREE  (the same designs, load CV 10%)")
    print("=" * 92)
    print(f"{'design':>24} {'increment dropped':>20} {'increment carried':>20} {'ratio':>9}")
    for name, t in designs.items():
        lv, dd = record(t, 0.10)
        a = sd(lv, dd, priors=False, carry=False)
        b = sd(lv, dd, priors=False, carry=True)
        r = a / b if np.isfinite(a) and np.isfinite(b) and b > 0 else np.nan
        f = lambda v: "singular" if not np.isfinite(v) else f"{v:.4f}"
        print(f"{name:>24} {f(a):>20} {f(b):>20} {r:9.1f}")

    print()
    print("=" * 92)
    print("THE GATE, PRIOR-FREE, against how well the increment is known")
    print("=" * 92)
    t = designs["two-daily x 2y (n=365)"]
    for cv in (0.05, 0.10, 0.20):
        lv, dd = record(t, cv)
        row = []
        for kap in (None, 0.05, 0.159, 0.30):
            pri = [np.inf] * 4 + ([kap ** 2] if kap else [])
            F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri])
            for m, d in zip(lv, dd):
                g = m / (P.Km + m)
                q = [P.theta * g, -(P.theta + d) * P.Km * m / (P.Km + m) ** 2,
                     P.ke * m, 1.0] + ([d * g] if kap else [])
                q = np.array(q)
                F = F + np.outer(q, q) / SIG ** 2
            w = np.linalg.eigvalsh(F)
            row.append(np.inf if w[0] <= 1e-14 * w[-1]
                       else float(np.sqrt(np.linalg.inv(F)[0, 0])))
        print(f"  load CV {cv:4.0%}:  Delta exact {row[0]:.4f}   "
              f"kappa 5% {row[1]:.4f}   kappa 15.9% {row[2]:.4f}   kappa 30% {row[3]:.4f}")


if __name__ == "__main__":
    main()


def which_prior():
    """Which of the four priors is doing the work?"""
    print()
    print("=" * 92)
    print("WHICH PRIOR IS LOAD-BEARING?  Two-daily x 2y, load CV 10%, increment carried.")
    print("=" * 92)
    t = np.arange(0.0, 730.0, 2.0)
    lv, dd = record(t, 0.10)
    names = ("ln theta", "ln K_m", "ln k_e", "load bias")

    def run(active):
        pri = [0.8 ** 2, 0.8 ** 2, 0.8 ** 2, (D.LOAD_FRAC * L0) ** 2]
        F = np.diag([1.0 / p if a else 0.0 for p, a in zip(pri, active)])
        for m, d in zip(lv, dd):
            q = np.array([P.theta * m / (P.Km + m),
                          -(P.theta + d) * P.Km * m / (P.Km + m) ** 2,
                          P.ke * m, 1.0])
            F = F + np.outer(q, q) / SIG ** 2
        w = np.linalg.eigvalsh(F)
        return (np.inf if w[0] <= 1e-14 * w[-1]
                else float(np.sqrt(np.linalg.inv(F)[0, 0])))

    print(f"  all four priors on:                      {run([1,1,1,1]):8.4f}")
    print(f"  none:                                    {run([0,0,0,0]):8.4f}")
    print()
    print("  dropping ONE prior at a time:")
    for i, nm in enumerate(names):
        a = [1, 1, 1, 1]; a[i] = 0
        print(f"    without the {nm:>10} prior:            {run(a):8.4f}")
    print()
    print("  keeping ONE prior at a time:")
    for i, nm in enumerate(names):
        a = [0, 0, 0, 0]; a[i] = 1
        v = run(a)
        print(f"    only the {nm:>10} prior:               "
              f"{'singular' if not np.isfinite(v) else f'{v:8.4f}'}")
    print()
    print("  the combination rem:no_km permits (no prior on ln K_m):")
    print(f"    theta + k_e + load, K_m free:            {run([1,0,1,1]):8.4f}")
    print(f"    k_e + load only, theta and K_m free:     {run([0,0,1,1]):8.4f}")
    print()
    print("  and how a WEAKER K_m prior behaves (theta, k_e, load priors kept):")
    for kmp in (0.4, 0.8, 1.2, 2.0, 3.0):
        pri = [0.8 ** 2, kmp ** 2, 0.8 ** 2, (D.LOAD_FRAC * L0) ** 2]
        F = np.linalg.inv(np.diag(pri))
        for m, d in zip(lv, dd):
            q = np.array([P.theta * m / (P.Km + m),
                          -(P.theta + d) * P.Km * m / (P.Km + m) ** 2,
                          P.ke * m, 1.0])
            F = F + np.outer(q, q) / SIG ** 2
        print(f"    ln K_m prior sd {kmp:4.1f} (a factor of {np.exp(kmp):5.1f}):  "
              f"sd(ln theta) = {np.sqrt(np.linalg.inv(F)[0,0]):.4f}")


which_prior()
