"""
How well the delivered enzymatic increment has to be known, on each footing.

The dose amplitude enters as a free scale kappa: Delta_i = kappa * delta_i with
delta_i read off the dosing record and kappa the per-patient conversion from
dose to enzymatic capacity.  Nobody measures kappa directly; a prior on it is
what a PK/PD characterisation of the product buys.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
from scipy.optimize import brentq
import design.run as D, model_A.run as LOC

P = D.P
Lb = D.load(0.0)
SIG = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)
GATE = 0.20


def sd_theta(lv, dd, kap_prior=None, ke_free=True, load_frac=D.LOAD_FRAC):
    """kap_prior=None means Delta known exactly (no kappa column)."""
    pri = [0.8 ** 2, 0.8 ** 2, (load_frac * Lb) ** 2]
    if ke_free:
        pri.append(0.8 ** 2)
    if kap_prior is not None:
        pri.append(kap_prior ** 2)
    Finv = np.linalg.inv(np.diag(pri))
    for m, d in zip(lv, dd):
        g = m / (P.Km + m)
        q = [P.theta * g, -(P.theta + d) * P.Km * m / (P.Km + m) ** 2, 1.0]
        if ke_free:
            q.append(P.ke * m)
        if kap_prior is not None:
            q.append(d * g)
        q = np.array(q)
        Finv = Finv + np.outer(q, q) / SIG ** 2
    w = np.linalg.eigvalsh(Finv)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(Finv)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def main():
    print(__doc__)
    for label, t in (("two-daily x 2 years (n=365)", np.arange(0.0, 730.0, 2.0)),
                     ("one dose cycle (n=7)", np.arange(0.0, D.DOSE_INTERVAL, 2.0)),
                     ("monthly, 2 years (n=24)", np.arange(0.0, 730.0, 30.4))):
        lv = D.quasi_steady_level(t)
        dd = np.asarray(D.delta(t), dtype=float)
        print("=" * 78)
        print(f"{label}")
        print("=" * 78)
        for ke_free in (False, True):
            head = "k_e known (published footing)" if not ke_free else \
                   "k_e fitted, 0.8 log-prior (the footing A1 recommends)"
            print(f"  {head}")
            print(f"    Delta known exactly:  sd(ln theta) = "
                  f"{sd_theta(lv, dd, None, ke_free):.4f}")
            for kp in (0.05, 0.10, 0.15, 0.20, 0.30, 0.50):
                v = sd_theta(lv, dd, kp, ke_free)
                flag = "PASS" if v < GATE else "fail"
                print(f"    kappa known to {kp:4.0%}:      sd(ln theta) = "
                      f"{v:8.4f}   {flag}")
            try:
                crit = brentq(lambda x: sd_theta(lv, dd, x, ke_free) - GATE,
                              1e-3, 2.0, xtol=1e-5)
                print(f"    -> the gate (0.20) needs kappa to {crit:.1%}")
            except ValueError:
                v0 = sd_theta(lv, dd, 1e-3, ke_free)
                print(f"    -> the gate is {'met at any kappa precision' if v0 < GATE else 'UNREACHABLE at any kappa precision'}"
                      f" (sd -> {v0:.4f} as the prior tightens)")
        print()

    print("=" * 78)
    print("What a paired breath test buys, if Delta is measured rather than inferred")
    print("=" * 78)
    print("  A 1-13C-propionate breath test measures theta + Delta directly, so a pair")
    print("  bracketing a dose gives Delta as a difference. Reported repeatability is a")
    print("  CV of 6.9% in controls and 12.8% in patients (breath_test_mma).")
    lv = D.quasi_steady_level(np.arange(0.0, 730.0, 2.0))
    dd = np.asarray(D.delta(np.arange(0.0, 730.0, 2.0)), dtype=float)
    for lab, cv in (("controls, single pair", 0.069),
                    ("patients, single pair", 0.128),
                    ("patients, difference of a pair", 0.128 * np.sqrt(2)),
                    ("patients, 2 pairs averaged", 0.128 * np.sqrt(2) / np.sqrt(2)),
                    ("patients, 4 pairs averaged", 0.128 * np.sqrt(2) / 2.0)):
        for ke_free in (True,):
            v = sd_theta(lv, dd, cv, ke_free)
            print(f"    {lab:>32}  CV {cv:5.1%}  ->  sd(ln theta) = {v:.4f}  "
                  f"{'PASS' if v < GATE else 'fail'}")


if __name__ == "__main__":
    main()
