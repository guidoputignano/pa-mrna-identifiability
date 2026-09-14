"""
2D.10 -- the second channel. How well must k_e be known, and are the two
channels substitutes or complements?

A fibre closes when an independent channel is known. Two candidate channels:

  Delta -- the delivered enzymatic increment, known to a relative precision
           `kap` by characterising the product (a PK/PD sub-study, or paired
           breath tests bracketing a dose).
  k_e   -- the non-enzymatic clearance of the marker, known to a relative
           precision `kep` in logs by some independent measurement.

Both are swept on the SAME footing: Delta carried in the sensitivity, the
module's 5% load-record prior, 0.8 log-prior on ln theta and ln K_m, and the
gate at sd(ln theta) < 0.20.

SUBSTITUTES would mean either channel alone reaches the gate.
COMPLEMENTS would mean neither alone does.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
from scipy.optimize import brentq
import design.run as D, model_A.run as LOC

P = D.P
L = D.load(0.0)
SIG = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)
GATE = 0.20
FLAT = 1e3          # a log-prior of sd 1000 is flat for any practical purpose


def sd_theta(lv, dd, kap=None, kep=0.8, load_frac=D.LOAD_FRAC, carry_delta=True):
    """
    kap  = relative prior sd on the dose-to-capacity scale; None = Delta known.
    kep  = prior sd on ln k_e; FLAT = k_e unknown, small = k_e known.
    """
    pri = [0.8 ** 2, 0.8 ** 2, (load_frac * L) ** 2, kep ** 2]
    if kap is not None:
        pri.append(kap ** 2)
    Finv = np.linalg.inv(np.diag(pri))
    for m, d in zip(lv, dd):
        g = m / (P.Km + m)
        eff = P.theta + (d if carry_delta else 0.0)
        q = [P.theta * g, -eff * P.Km * m / (P.Km + m) ** 2, 1.0, P.ke * m]
        if kap is not None:
            q.append(d * g)
        q = np.array(q)
        Finv = Finv + np.outer(q, q) / SIG ** 2
    w = np.linalg.eigvalsh(Finv)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(Finv)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def crit(f, lo, hi):
    try:
        return brentq(f, lo, hi, xtol=1e-5)
    except ValueError:
        return None


def main():
    print(__doc__)
    t = np.arange(0.0, 730.0, 2.0)
    lv, dd = D.quasi_steady_level(t), np.asarray(D.delta(t), dtype=float)

    print("=" * 76)
    print("THE 2x2.  Gate is sd(ln theta) < 0.20, two-daily x 2 years.")
    print("=" * 76)
    print(f"{'':>34} {'k_e UNKNOWN':>14} {'k_e known 10%':>15} {'k_e known exactly':>18}")
    for lab, kap in (("Delta UNKNOWN (scale flat)", FLAT),
                     ("Delta known to 20%", 0.20),
                     ("Delta known to 16%", 0.16),
                     ("Delta known to 10%", 0.10),
                     ("Delta known exactly", None)):
        row = [sd_theta(lv, dd, kap, kep) for kep in (FLAT, 0.10, 1e-4)]
        cells = "  ".join(f"{v:12.4f}{'*' if v < GATE else ' '}" for v in row)
        print(f"{lab:>34}  {cells}")
    print("   (* passes the gate)")

    print()
    print("=" * 76)
    print("EACH CHANNEL ALONE -- how precise does it have to be?")
    print("=" * 76)
    c1 = crit(lambda x: sd_theta(lv, dd, x, 0.8) - GATE, 1e-3, 5.0)
    print(f"  Delta alone (k_e at the module's own 0.8 log-prior):"
          f"  needs the dose scale to {c1:.1%}" if c1 else
          "  Delta alone: no precision reaches the gate")
    c1b = crit(lambda x: sd_theta(lv, dd, x, FLAT) - GATE, 1e-3, 5.0)
    print(f"  Delta alone (k_e genuinely unknown):"
          f"  {'needs the dose scale to %.1f%%' % (100*c1b) if c1b else 'NEVER reaches the gate'}")
    c2 = crit(lambda x: sd_theta(lv, dd, FLAT, x) - GATE, 1e-4, 5.0)
    print(f"  k_e alone (dose scale genuinely unknown):"
          f"  {'needs ln k_e to %.1f%%' % (100*c2) if c2 else 'NEVER reaches the gate'}")
    c2b = crit(lambda x: sd_theta(lv, dd, None, x) - GATE, 1e-4, 5.0)
    print(f"  k_e with Delta known exactly:"
          f"  {'needs ln k_e to %.1f%%' % (100*c2b) if c2b else 'any precision passes'}")
    print(f"     (for reference, k_e unknown + Delta exact = "
          f"{sd_theta(lv, dd, None, FLAT):.4f}; k_e at 0.8 + Delta exact = "
          f"{sd_theta(lv, dd, None, 0.8):.4f})")

    print()
    print("=" * 76)
    print("THE TRADE-OFF CURVE -- how much k_e precision buys how much Delta slack")
    print("=" * 76)
    print(f"{'prior on ln k_e':>18} | {'dose scale the gate then needs':>32}")
    for kep in (FLAT, 0.80, 0.50, 0.30, 0.20, 0.10, 0.05, 1e-4):
        c = crit(lambda x: sd_theta(lv, dd, x, kep) - GATE, 1e-3, 5.0)
        lab = "flat" if kep >= FLAT else (f"{kep:.0%}" if kep > 1e-3 else "exact")
        base = sd_theta(lv, dd, None, kep)
        print(f"{lab:>18} | {('%.1f%%' % (100*c)) if c else 'unreachable':>20}"
              f"   (Delta exact -> {base:.4f})")

    print()
    print("=" * 76)
    print("IS k_e A REAL CHANNEL AT ALL?  What it is worth with Delta DROPPED")
    print("=" * 76)
    print("  (i.e. on the pre-2C footing, where the therapy channel is not used)")
    for kep in (FLAT, 0.8, 0.1, 1e-4):
        v = sd_theta(lv, dd, None, kep, carry_delta=False)
        lab = "flat" if kep >= FLAT else (f"{kep:.0%}" if kep > 1e-3 else "exact")
        print(f"    ln k_e prior {lab:>6}:  sd(ln theta) = {v:.4f}"
              f"  {'PASS' if v < GATE else 'fail'}")


if __name__ == "__main__":
    main()
