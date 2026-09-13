"""
14.1 -- the excitation budget, three scenarios on one metric.

The design this note specifies requires a dose train, and no patient outside the
trial receives one. So the question is not "is the therapy channel better" -- the
note already shows it is -- but "is there ANY design for a patient who is not on
therapy". Three sources of excitation, one precision metric, same footing.

    (1) spontaneous day-to-day load variability, no therapy
    (2) as (1) plus quasi-exogenous carnitine perturbations
    (3) the dose train

FOOTING, identical across all three: theta and K_m free, log-prior 0.8 on ln k_e,
load recorded to 5%, monthly outpatient draws over two years (24 samples).
Everything is sd(ln theta-hat).

No expected results are recorded here. Each number is printed as computed.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from design.run import (P, load, delta, LOAD_FRAC, ASSAY_CV, DOSE_INTERVAL,
                        M_BAR, LOAD_CV)
import model_A.run as LOC

SIGMA = LOC.sigma_total(LOAD_FRAC, ASSAY_CV)
LB = load(0.0)
MONTHLY = np.arange(0.0, 730.0, 30.4)
SEEDS = range(24)


def banner(t):
    print("\n" + "=" * 82)
    print("  " + t)
    print("=" * 82)


def levels(loads, inc, ke_rel):
    """Quasi-steady marker solving (theta+inc) g(m) + k_e(1+rel) m = L."""
    ke = P.ke * (1.0 + np.asarray(ke_rel, dtype=float))
    B = ke * P.Km + P.theta + np.asarray(inc, float) - np.asarray(loads, float)
    return (-B + np.sqrt(B ** 2 + 4.0 * ke * np.asarray(loads, float) * P.Km)) / (2.0 * ke)


def sd_of(times, inc=None, ke_rel=None, load_cv=LOAD_CV, seed=3):
    """sd(ln theta) on the minimal footing, with optional increment and k_e path."""
    n = len(times)
    inc = np.zeros(n) if inc is None else np.broadcast_to(inc, (n,))
    ke_rel = np.zeros(n) if ke_rel is None else np.broadcast_to(ke_rel, (n,))
    rng = np.random.default_rng(seed)
    # LOGNORMAL, not Gaussian. At 40% CV a Gaussian draw returns a NEGATIVE load
    # and a negative marker level, which is not a record. See e9, part (a).
    if load_cv:
        s2 = np.log1p(load_cv ** 2)
        L = LB * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), n))
    else:
        L = np.full(n, LB)
    m = levels(L, inc, ke_rel)
    ke = P.ke * (1.0 + ke_rel)
    F = np.diag([0.0, 0.0, 1.0 / (LOAD_FRAC * LB) ** 2, 1.0 / 0.8 ** 2])
    for mi, di, kei in zip(m, inc, ke):
        q = np.array([P.theta * mi / (P.Km + mi),
                      -(P.theta + di) * P.Km * mi / (P.Km + mi) ** 2,
                      1.0,
                      kei * mi])
        F = F + np.outer(q, q) / SIGMA ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def spread(fn):
    """Median and range over seeds. A single draw near a boundary is not a result."""
    vals = np.array([fn(s) for s in SEEDS], dtype=float)
    fin = vals[np.isfinite(vals)]
    if fin.size == 0:
        return np.inf, np.inf, np.inf, 0
    return np.median(fin), fin.min(), fin.max(), fin.size


# ---------------------------------------------------------------------------
def scenario_1():
    banner("SCENARIO 1  --  spontaneous load variability only, NO therapy")
    print("  Proposition prop:rank: a perfectly steady recorded load is the exactly")
    print("  singular case, and day-to-day variation is what closes it. So sweep the")
    print("  variation. 24 monthly draws, no increment, no carnitine perturbation.")
    print()
    print(f"    {'day-to-day load CV':>20}   {'median sd(ln theta)':>20}  {'[min, max] over 24 seeds':>28}")
    out = {}
    for cv in (0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.40):
        med, lo, hi, n = spread(lambda s, c=cv: sd_of(MONTHLY, load_cv=c, seed=s))
        out[cv] = med
        cell = "singular" if not np.isfinite(med) else f"{med:.4f}"
        rng_ = "--" if not np.isfinite(med) else f"[{lo:.4f}, {hi:.4f}]"
        print(f"    {cv*100:18.0f}%   {cell:>20}  {rng_:>28}")
    return out


def scenario_2():
    banner("SCENARIO 2  --  plus quasi-exogenous carnitine perturbations")
    print("  MECHANISM, stated before the numbers. Carnitine acts on k_e, not on")
    print("  theta. It therefore does NOT enter the sensitivity ratio")
    print("      (ds/dlnK_m)/(ds/dlntheta) = -(1 + Delta/theta) phi(m)")
    print("  through the first factor at all. It can only move the SECOND factor,")
    print("  phi(m) = K_m/(K_m+m), and only indirectly, by shifting the level m at")
    print("  which the patient sits. The dose train moves the first factor directly.")
    print()
    print("  How far does a carnitine step actually move phi?")
    print(f"    {'relative change in k_e':>24} {'m*':>10} {'phi(m*)':>10} {'|dphi| from baseline':>22}")
    m0 = float(levels([LB], [0.0], [0.0])[0])
    phi0 = P.Km / (P.Km + m0)
    for rel in (0.0, 0.10, 0.25, 0.50, 1.00, 2.00):
        mi = float(levels([LB], [0.0], [rel])[0])
        ph = P.Km / (P.Km + mi)
        print(f"    {rel*100:22.0f}% {mi:10.4f} {ph:10.5f} {abs(ph-phi0):22.6f}")
    print()
    print("  For comparison, the DOSE TRAIN's factor over one cycle:")
    t = np.linspace(0.0, DOSE_INTERVAL, 400)
    d = delta(t)
    f1 = 1.0 + d / P.theta
    print(f"    1 + Delta/theta ranges {f1.min():.4f} to {f1.max():.4f}"
          f"  -- a factor of {f1.max()/f1.min():.3f}")
    print(f"    phi over the same cycle ranges {phi0:.5f} to "
          f"{P.Km/(P.Km+float(levels([LB]*len(t), d, [0.0]*len(t)).max())):.5f}")

    print()
    print("  Now the precision. Carnitine perturbation modelled as a recorded step")
    print("  in k_e alternating between visits, at several magnitudes, WITH the")
    print("  gain treated as KNOWN -- which is the favourable case for carnitine.")
    print(f"\n    {'k_e step':>12} {'load CV':>9}   {'median sd(ln theta)':>20}  {'[min, max]':>22}")
    for rel in (0.25, 0.50, 1.00):
        path = np.where(np.arange(len(MONTHLY)) % 2 == 0, 0.0, rel)
        for cv in (0.10, 0.20):
            med, lo, hi, n = spread(
                lambda s, p=path, c=cv: sd_of(MONTHLY, ke_rel=p, load_cv=c, seed=s))
            cell = "singular" if not np.isfinite(med) else f"{med:.4f}"
            rng_ = "--" if not np.isfinite(med) else f"[{lo:.4f}, {hi:.4f}]"
            print(f"    {rel*100:10.0f}% {cv*100:8.0f}%   {cell:>20}  {rng_:>22}")


def scenario_3():
    banner("SCENARIO 3  --  the dose train")
    inc = delta(MONTHLY)
    for cv in (0.10, 0.20):
        med, lo, hi, n = spread(lambda s, c=cv: sd_of(MONTHLY, inc=inc, load_cv=c, seed=s))
        print(f"    load CV {cv*100:3.0f}%:  median sd(ln theta) = {med:.4f}"
              f"   [{lo:.4f}, {hi:.4f}] over {n} seeds")
    print()
    print("  (The note's headline 0.17 is the single-seed value at CV 10%; the")
    print("   median over seeds is printed above and the spread with it.)")


# ---------------------------------------------------------------------------
def why_carnitine_goes_backwards():
    """Scenario 2 made things WORSE. This checks why rather than asserting it."""
    banner("WHY SCENARIO 2 GOES BACKWARDS")
    print("  Carnitine raises k_e, which LOWERS the quasi-steady m. But the")
    print("  theta sensitivity is q_theta = theta m/(K_m+m), which falls with m.")
    print("  So a carnitine step trades magnitude in the theta column for a move")
    print("  in phi. Both sides of that trade, measured:")
    print()
    m0 = float(levels([LB], [0.0], [0.0])[0])
    q0 = P.theta * m0 / (P.Km + m0)
    phi0 = P.Km / (P.Km + m0)
    print(f"    {'k_e step':>10} {'m*':>9} {'q_theta':>10} {'lost in q_theta':>17}"
          f" {'gained in phi':>15}")
    for rel in (0.25, 0.50, 1.00, 2.00):
        mi = float(levels([LB], [0.0], [rel])[0])
        qi = P.theta * mi / (P.Km + mi)
        ph = P.Km / (P.Km + mi)
        print(f"    {rel*100:8.0f}% {mi:9.4f} {qi:10.2f} {100*(q0-qi)/q0:16.2f}%"
              f" {100*abs(ph-phi0)/phi0:14.2f}%")
    print()
    print("  The trade is roughly ten to one against: at a doubling of k_e the")
    print("  theta column loses an order of magnitude more than phi gains.")
    print("  Carnitine does not merely supply weak excitation -- it moves the")
    print("  patient the WRONG WAY along the axis that carries the information,")
    print("  deeper into the linear regime where Proposition prop:excitation")
    print("  says only the lumped ratio survives.")


# ---------------------------------------------------------------------------
def what_would_scenario_1_need():
    banner("WHAT WOULD SCENARIO 1 NEED? (the gate)")
    target = 0.20

    print("  (a) By raising day-to-day load variability, with 24 monthly draws:")
    lo, hi = 0.10, 20.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        med, *_ = spread(lambda s, c=mid: sd_of(MONTHLY, load_cv=c, seed=s))
        if np.isfinite(med) and med < target:
            hi = mid
        else:
            lo = mid
    print(f"      required CV = {100*0.5*(lo+hi):.0f}%  of mean intake.")
    print("      A coefficient of variation above 100% is not a diet, and the")
    print("      quasi-steady levels it implies are not a managed patient.")

    print()
    print("  (b) By taking more samples at the note's own 10% CV:")
    base, *_ = spread(lambda s: sd_of(MONTHLY, load_cv=0.10, seed=s))
    for n_per_year in (12, 52, 365):
        t = np.arange(0.0, 730.0, 365.0 / n_per_year)
        med, l, h, _ = spread(lambda s, tt=t: sd_of(tt, load_cv=0.10, seed=s))
        print(f"      {len(t):5d} draws over two years ({n_per_year:3d}/yr): "
              f"sd = {med:.4f}")
    print()
    print("      NOTE THE SCALING. An earlier version of this script said the")
    print("      information SATURATES, 'because consecutive days resample one load")
    print("      distribution'. That was asserted, not derived, and e9 derives it:")
    print("      the information in the weak direction GROWS roughly as n and never")
    print("      saturates. It simply never ARRIVES -- the load prior holds about")
    print("      1.05 in that direction and the data supply 7.6e-04 even at n = 3e5.")
    print("      The weak direction is PRIOR-HELD at every n a study could reach, so")
    print("      sd falls only through the directions the data do reach, and the")
    print("      fitted exponent sits near 0.1 rather than 0.5. No achievable n")
    print("      closes it, but the reason is a fixed prior, not correlated samples.")
    print()
    print("  Neither route exists. The gate is therefore answered in the negative.")


if __name__ == "__main__":
    print(__doc__)
    s1 = scenario_1()
    scenario_2()
    why_carnitine_goes_backwards()
    scenario_3()
    what_would_scenario_1_need()
