"""
What actually buys the information when the increment is carried?

A first guess -- that a varying Delta breaks the sub-saturation degeneracy of
Proposition 2(1) -- does not survive the arithmetic. On an idealised band with
Delta sweeping monotonically with m, carrying it makes sd(ln K_m) much better
(0.5870 -> 0.2089) and sd(ln theta) slightly WORSE (0.5760 -> 0.6647).

Yet on the real dose train at the same operating point, carrying it takes
sd(ln theta) from 0.5020 to 0.1475. So the gain is not "Delta varies". The
test below isolates what it is.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
import design.run as D, model_A.run as LOC

P = D.P
Lb = D.load(0.0)
SIG = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)


def sd(lv, dd, carry=True):
    pri = [0.8 ** 2, 0.8 ** 2, (D.LOAD_FRAC * Lb) ** 2, 0.8 ** 2]
    F = np.linalg.inv(np.diag(pri))
    for m, d in zip(lv, dd):
        eff = P.theta + (d if carry else 0.0)
        q = np.array([P.theta * m / (P.Km + m),
                      -eff * P.Km * m / (P.Km + m) ** 2, 1.0, P.ke * m])
        F = F + np.outer(q, q) / SIG ** 2
    C = np.linalg.inv(F)
    return float(np.sqrt(C[0, 0])), float(np.sqrt(C[1, 1]))


def main():
    print(__doc__)
    t = np.arange(0.0, 730.0, 2.0)
    lv = D.quasi_steady_level(t)
    dd = np.asarray(D.delta(t), dtype=float)

    # How many DISTINCT increments does each level see?
    key = np.round(lv, 6)
    uniq = {}
    for k, d in zip(key, dd):
        uniq.setdefault(k, set()).add(round(d, 6))
    revisits = np.array([len(v) for v in uniq.values()])

    print("=" * 84)
    print("THE REAL DESIGN")
    print("=" * 84)
    a, b = sd(lv, dd, False); c, e = sd(lv, dd, True)
    print(f"  n = {len(lv)},  {len(uniq)} distinct levels,  "
          f"levels revisited at {revisits.mean():.1f} distinct increments on average")
    print(f"  corr(m, Delta) over the record = {np.corrcoef(lv, dd)[0,1]:+.4f}")
    print(f"  increment dropped: sd(ln theta) = {a:.4f}   carried: {c:.4f}")

    print()
    print("=" * 84)
    print("THE SAME LEVELS, WITH THE INCREMENT REASSIGNED")
    print("=" * 84)
    print("  Each row keeps the identical set of marker levels and the identical set")
    print("  of increments, and only changes which increment goes with which level.")
    order = np.argsort(lv)
    variants = {
        "as the dose train delivers them": dd,
        "sorted: Delta monotone increasing in m": np.sort(dd)[np.argsort(order)],
        "sorted: Delta monotone decreasing in m": np.sort(dd)[::-1][np.argsort(order)],
        "shuffled at random": np.random.default_rng(0).permutation(dd),
        "held at its mean (no variation at all)": np.full_like(dd, dd.mean()),
    }
    print(f"  {'assignment':>42} {'corr(m,D)':>10} {'sd(ln th)':>10} {'sd(ln Km)':>10}")
    for lab, d in variants.items():
        x, y = sd(lv, d, True)
        print(f"  {lab:>42} {np.corrcoef(lv, d)[0,1]:+10.4f} {x:10.4f} {y:10.4f}")

    print()
    print("  The dose train and the random shuffle are the same design in every respect")
    print("  except the pairing, and the monotone assignments are the same again. What")
    print("  separates them is whether ONE LEVEL IS VISITED AT SEVERAL INCREMENTS.")

    print()
    print("=" * 84)
    print("THE CLEANEST FORM OF THE STATEMENT")
    print("=" * 84)
    print("  Take a SINGLE marker level and vary only the increment. This isolates the")
    print("  mechanism, and it is an illustration, NOT a design -- see the caveat below.")
    m0 = float(np.median(lv))
    for k in (1, 2, 3, 5, 9):
        d = np.linspace(0.0, 2.0 * P.theta, k)
        d = np.tile(d, int(np.ceil(365 / k)))[:365]
        x, y = sd(np.full(365, m0), d, True)
        print(f"    {k:2d} distinct increments at the single level m = {m0:.2f}:"
              f"  sd(ln theta) = {x:8.4f}   sd(ln K_m) = {y:8.4f}")
    print()
    print("  One level and one increment is rank one and nothing is learnt. Two")
    print("  increments make it rank two and theta moves to 0.166 with no excursion at")
    print("  all: the therapy channel, stated without reference to the band.")
    print()
    print("  THE CAVEAT, and it is why this is not a design. At a single level k_e m and")
    print("  the constant load-bias column are proportional, so the matrix is rank two of")
    print("  four however many increments are used, and two directions are carried by")
    print("  priors rather than by data. Prior-free the design is singular. Put the")
    print("  dose scale in as a fifth parameter at the 15.9% the gate requires and it")
    print("  gives 0.2209, which FAILS. The real dose-train design gives 0.1475 and")
    print("  passes, because it varies the level as well and so separates k_e from the")
    print("  load bias. The increment is what identifies theta; the level variation is")
    print("  what keeps the nuisance parameters apart.")


if __name__ == "__main__":
    main()
