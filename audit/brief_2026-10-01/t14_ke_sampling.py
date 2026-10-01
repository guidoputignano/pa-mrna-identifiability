"""
Item 12. Replace the 'k_e exact' premise with a trade-off surface.

sd(ln theta_hat) at 5% of normal as a function of
    - the external precision on k_e, sd(ln k_e), and
    - the number of blood draws per year,
over a fixed two-year programme (so n_draws = 2 * draws_per_year), each draw on
the booked dose train, median over the 24 intake seeds.

Contours at 0.20 and 0.30. The achievable region uses this project's own
findings: per-patient k_e is reachable to about 0.2 log units (item 15), and
every-other-day venipuncture (~182/yr) is the clinical ceiling, weekly (52/yr)
already demanding, monthly (12/yr) routine.
"""
import numpy as np

import core as C

HORIZON = 730.0  # days, two years, as in the monthly record


def record_n(n, seed):
    """n draws spread uniformly over the two-year window, dose train on."""
    times = np.linspace(0.0, HORIZON, int(n))
    return C.record(C.TH_5, seed=seed, times=times)


def sd_at(ke_sd, n):
    def one(s):
        return C.sd_th(C.fim(record_n(n, s), noise="perdraw", ke_prior=ke_sd))
    return C.over_seeds(one)[0]


def crossing_n(ke_sd, gate, ns):
    """smallest n in the grid at which the median sd <= gate (interp), or None."""
    vals = np.array([sd_at(ke_sd, n) for n in ns])
    below = np.where(vals <= gate)[0]
    if not below.size:
        return None, vals
    i = below[0]
    if i == 0:
        return float(ns[0]), vals
    #  log-log interpolate the crossing between ns[i-1] and ns[i]
    x0, x1 = np.log(ns[i - 1]), np.log(ns[i])
    y0, y1 = np.log(vals[i - 1]), np.log(vals[i])
    xc = x0 + (np.log(gate) - y0) * (x1 - x0) / (y1 - y0)
    return float(np.exp(xc)), vals


def main():
    print("=" * 96)
    print("  ITEM 12.  sd(ln theta) AT 5% OVER (sd(ln k_e), draws per year), TWO-YEAR PROGRAMME")
    print("=" * 96)

    ke_sds = [1e-4, 0.05, 0.10, 0.20, 0.40, 0.80]   # 1e-4 ~ 'exact'
    ns = np.unique(np.round(np.geomspace(24, 1460, 18)).astype(int))  # 12..730 /yr
    print(f"  draws over two years: {list(ns)}")
    print(f"  (draws per year = n/2: {[round(n/2) for n in ns]})\n")

    print("  sd(ln theta), median over 24 seeds:")
    head = "  sd(ln k_e) \\ n  " + "".join(f"{n:>7d}" for n in ns)
    print(head)
    table = {}
    for ke_sd in ke_sds:
        row = [sd_at(ke_sd, n) for n in ns]
        table[ke_sd] = np.array(row)
        lab = "exact" if ke_sd < 1e-3 else f"{ke_sd:.2f}"
        print(f"  {lab:>12}   " + "".join(f"{v:7.3f}" for v in row))

    print("\n  CONTOUR CROSSINGS (draws over two years; draws/yr in brackets):")
    for gate in (0.30, 0.20):
        print(f"    sd = {gate:.2f}:")
        for ke_sd in ke_sds:
            nc, _ = crossing_n(ke_sd, gate, ns)
            lab = "exact" if ke_sd < 1e-3 else f"{ke_sd:.2f}"
            if nc is None:
                print(f"      sd(ln k_e) {lab:>5}: not reached by {ns[-1]} draws "
                      f"({ns[-1]//2}/yr)")
            else:
                print(f"      sd(ln k_e) {lab:>5}: n = {nc:7.0f}  ({nc/2:6.0f}/yr)")

    print("\n  FEASIBILITY BOX (this project's own numbers):")
    print("    achievable k_e precision:  sd(ln k_e) >= 0.20  (item 15: ~0.2 log units)")
    print("    clinically feasible draws: <= 52/yr weekly (demanding), <= 26/yr routine;")
    print("                               182/yr (every other day) is the stated ceiling.")
    #  the decisive query: best feasible sd, i.e. k_e at 0.20 and draws at the ceiling
    for rate, lab in ((26, "biweekly 26/yr"), (52, "weekly 52/yr"),
                      (182, "every-other-day 182/yr (ceiling)")):
        v = sd_at(0.20, 2 * rate)
        print(f"    sd(ln k_e)=0.20, {lab:>32}: sd(ln theta) = {v:.3f}")
    #  and with k_e exact at the ceiling, for reference
    print(f"    sd(ln k_e)=exact, every-other-day 182/yr        : "
          f"sd(ln theta) = {sd_at(1e-4, 2*182):.3f}")
    print("\n  FEASIBLE POINT ON THE 0.20 CONTOUR?  Read it off the box above.")


if __name__ == "__main__":
    main()
