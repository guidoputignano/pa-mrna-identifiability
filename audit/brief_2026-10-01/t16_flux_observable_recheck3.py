"""
ADVERSARIAL RECHECK of ITEM 11 (out/t16_flux_observable.txt).

Independent of t16_flux_observable.py. Imports ONLY the shared primitives
(core, t5_value, t12_flux_mc). Every FIM here is rebuilt from scratch:
  - the flux reading row (d (c theta g)/d params) is reconstructed from m_bar,
    not read from t16;
  - the 5x5 (ln theta, ln Km, b, ln ke, ln c) FIM is assembled by hand and
    inverted directly;
  - the single-calibrated reading at Y_eff is assembled by hand too;
  - the (s,L) share sensitivity is finite-differenced through the balance solve
    (core.m_ss) with my own step, independent of t16.share_row.

Claims under test:
  (a) sd(ln theta) over the (Y,Z) grid EQUALS a single calibrated reading at
      Y_eff = hypot(Y,Z) to machine precision; 0.20 contour hypot(Y,Z) <= Y*
      with Y* ~ 19.3% (FIM); cells (10%,10%)->~0.150, (20%,exact)->~0.206.
  (b) (s,L) reparam reaches 0.20 on essentially the same locus (crossing in Y
      ~19-20%, weakly dependent on X in 1-5%).
"""
import numpy as np

import core as C
import t5_value as T5
import t12_flux_mc as T12

CV = C.ASSAY_CV          # 0.05 assay CV
TH = C.TH_5              # 150
KM = C.KM0               # 100
GATE = C.GATE            # 0.20
SEEDS = C.SEEDS


# ----------------------------------------------------------------------
#  (a)  INDEPENDENT reconstruction of the flux reading and the 5x5 FIM
# ----------------------------------------------------------------------
def my_flux_row(r):
    """d (theta g(m0)) / d(ln theta, ln Km, b, ln ke) rebuilt from m_bar."""
    m0 = r["m_bar"]
    K = r["Km"]
    th = r["theta"]
    g = m0 / (K + m0)
    F = th * g
    d_lnth = F                                  # d(th g)/d ln th = th g
    d_lnKm = -th * K * m0 / (K + m0) ** 2        # d(th g)/d ln Km
    eps = K / (K + m0)                           # d ln g / d ln m = Km/(Km+m)
    return np.array([d_lnth, d_lnKm, 0.0, 0.0]), F, eps


def my_fim5(r, Y, Z):
    """5x5 FIM in (ln th, ln Km, b, ln ke, ln c); marginal sd(ln theta)."""
    F4 = C.fim(r)                                # record (incl. baseline draw)
    row4, F, eps = my_flux_row(r)
    sdf = F * np.hypot(Y, eps * CV)              # absolute reading sd
    F5 = np.zeros((5, 5))
    F5[:4, :4] = F4
    row5 = np.r_[row4, F]                         # d(c th g)/d ln c = F at c=1
    F5 += np.outer(row5, row5) / sdf ** 2
    F5[4, 4] += Z ** -2                           # Gaussian prior on ln c, sd Z
    return float(np.sqrt(np.linalg.inv(F5)[0, 0]))


def my_single(r, Yeff):
    """single calibrated flux reading at Yeff, 4x4; marginal sd(ln theta)."""
    F4 = C.fim(r)
    row4, F, eps = my_flux_row(r)
    sdf = F * np.hypot(Yeff, eps * CV)
    return C.sd_th(F4 + np.outer(row4, row4) / sdf ** 2)


def med5(Y, Z):
    return float(np.median([my_fim5(T5.base_record(s, baseline=True), Y, Z)
                            for s in SEEDS]))


def med_single(Yeff):
    return float(np.median([my_single(T5.base_record(s, baseline=True), Yeff)
                            for s in SEEDS]))


# ----------------------------------------------------------------------
#  (b)  INDEPENDENT share sensitivity by finite-differencing the balance
# ----------------------------------------------------------------------
def my_share(th, K, ke, b, Lrec):
    """fraction oxidised s = theta g(m0)/(L-b), m0 the untreated balance root."""
    m0 = C.m_ss(th, K, ke, Lrec - b)
    return th * m0 / (K + m0) / (Lrec - b)


def my_share_row(r, h=1e-6):
    """s and d s/d(ln th, ln Km, b, ln ke) by central differences of the solve."""
    th, K, ke, Lrec = r["theta"], r["Km"], r["ke"], r["L"]
    base = np.array([np.log(th), np.log(K), 0.0, np.log(ke)])
    s0 = my_share(th, K, ke, 0.0, Lrec)
    g = np.zeros(4)
    for j in range(4):
        hj = h * (Lrec if j == C.I_B else 1.0)
        xp, xm = base.copy(), base.copy()
        xp[j] += hj
        xm[j] -= hj

        def sh(x):
            return my_share(np.exp(x[0]), np.exp(x[1]), np.exp(x[3]), x[2], Lrec)
        g[j] = (sh(xp) - sh(xm)) / (2 * hj)
    return s0, g


def my_sd_sL(Y, X, h=1e-6):
    """record + one s-reading to Y% + load to X% (bias prior X*L); median."""
    def one(seed):
        r = T5.base_record(seed, baseline=True)
        s0, g = my_share_row(r, h=h)
        return C.sd_th(C.fim(r, b_prior=max(X, 1e-7) * r["L"],
                             extra=[(g, Y * s0)]))
    return C.over_seeds(one)[0]


def crossing_in_Y(f, lo=0.01, hi=0.60):
    if f(lo) > GATE:
        return None
    if f(hi) <= GATE:
        return hi
    for _ in range(50):
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if f(mid) <= GATE else (lo, mid)
    return lo


# ----------------------------------------------------------------------
def main():
    print("=" * 96)
    print("  RECHECK3 of ITEM 11 -- independent rebuild of the FIMs")
    print("=" * 96)

    Ys = [0.05, 0.10, 0.15, 0.175, 0.20, 0.25]
    Zs = [1e-7, 0.05, 0.10, 0.15, 0.20]

    # -- (a.1) the (Y,Z) grid and the hypot collapse --------------------
    print("\n  (a) 5x5 marginal sd(ln theta), median over 24 seeds (rows Z, cols Y):")
    print("        Z \\ Y   " + "".join(f"{y*100:7.1f}%" for y in Ys))
    worst = 0.0
    worst_at = None
    grid5 = {}
    for Z in Zs:
        rowvals = []
        for Y in Ys:
            Zeff = 0.0 if Z < 1e-6 else Z
            v5 = med5(Y, Z) if Z >= 1e-6 else med_single(Y)  # Z=exact: no ln c
            yeff = float(np.hypot(Y, Zeff))
            ref = med_single(yeff)
            gap = abs(v5 - ref)
            if gap > worst:
                worst, worst_at = gap, (Y, Z)
            rowvals.append(v5)
            grid5[(round(Y, 4), round(Z, 7))] = v5
        lab = "exact" if Z < 1e-6 else f"{Z:.2f}"
        print(f"      {lab:>7}   " + "".join(f"{v:7.3f} " for v in rowvals))
    print(f"\n      hypot collapse: worst |5x5 - single(Y_eff)| over grid = {worst:.2e}"
          f"  at (Y,Z)={worst_at}")

    # cross-check the hand-built 5x5 against the allowed primitive t12.fim_sd
    worstp = 0.0
    for Z in Zs:
        for Y in Ys:
            mine = med5(Y, Z) if Z >= 1e-6 else med_single(Y)
            prim = T12.fim_seeds(Y=Y, cal_sd=(None if Z < 1e-6 else Z))
            worstp = max(worstp, abs(mine - prim))
    print(f"      cross-check vs t12.fim_sd primitive: worst |mine - t12| = {worstp:.2e}")

    # -- (a.2) the single-reading crossing Y* ---------------------------
    ystar = crossing_in_Y(lambda y: med_single(y), lo=0.05, hi=0.40)
    print(f"\n      single calibrated-reading crossing Y* = {100*ystar:.3f}%"
          f"   (sd at Y* = {med_single(ystar):.4f})")
    print(f"        bracket: sd(Y=19.0%)={med_single(0.19):.4f},"
          f" sd(Y=19.33%)={med_single(0.1933):.4f}, sd(Y=19.7%)={med_single(0.197):.4f}")

    # verify the contour hypot(Y,Z) <= Y*: for several Z find the Y crossing
    print("\n      contour test -- for each Z, find Y where sd=0.20, report hypot(Y,Z):")
    for Z in (0.0, 0.05, 0.10, 0.15):
        if Z == 0.0:
            yc = crossing_in_Y(lambda y: med_single(y), lo=0.01, hi=0.40)
        else:
            yc = crossing_in_Y(lambda y: med5(y, Z), lo=0.01, hi=0.40)
        hyp = np.hypot(yc, Z) if yc is not None else float("nan")
        print(f"        Z={Z:4.2f}:  Y*={100*yc:6.3f}%   hypot(Y*,Z)={100*hyp:6.3f}%")

    # -- (a.3) the two named cells --------------------------------------
    c_10_10 = med5(0.10, 0.10)
    c_20_ex = med_single(0.20)
    print(f"\n      named cells:  (Y=10%,Z=10%) -> {c_10_10:.3f}  (claim ~0.150)")
    print(f"                    (Y=20%,Z=exact) -> {c_20_ex:.3f}  (claim ~0.206)")
    # also reproduce the whole brief table for Z in {exact,.05,.10,.15,.20}
    print("      brief-table reproduction (cols Y=5,10,15,17.5,20,25):")
    for Z in Zs:
        lab = "exact" if Z < 1e-6 else f"{Z:.2f}"
        vals = [grid5[(round(Y, 4), round(Z, 7))] for Y in Ys]
        print(f"        {lab:>7} " + "".join(f"{v:6.3f}" for v in vals))

    # -- (b) the (s,L) reparameterisation -------------------------------
    print("\n  (b) (s,L) reparam: sd(ln theta), median over seeds (rows X=load prec):")
    Xs = [1e-7, 0.01, 0.02, 0.05]
    print("        X \\ Y   " + "".join(f"{y*100:7.1f}%" for y in Ys))
    for X in Xs:
        rowvals = [my_sd_sL(Y, X) for Y in Ys]
        lab = "exact" if X < 1e-6 else f"{X*100:.0f}%"
        print(f"      L {lab:>5}  " + "".join(f"{v:7.3f} " for v in rowvals))

    print("\n      0.20 crossing in Y per load precision X (my own FD, h=1e-6):")
    cross = {}
    for X, lab in ((1e-7, "L exact"), (0.01, "L to 1%"), (0.02, "L to 2%"),
                   (0.05, "L to 5%")):
        yc = crossing_in_Y(lambda Y: my_sd_sL(Y, X))
        cross[lab] = yc
        if yc is None:
            print(f"        {lab:>8}: never reaches 0.20 in [1%,60%]")
        else:
            print(f"        {lab:>8}: s to Y <= {100*yc:.2f}%   "
                  f"(hypot(Y,X)={100*np.hypot(yc, X):.1f}%)")

    # FD robustness: redo the crossings with a coarser step
    print("\n      FD robustness -- crossings recomputed with h=1e-5:")
    for X, lab in ((1e-7, "L exact"), (0.05, "L to 5%")):
        yc = crossing_in_Y(lambda Y: my_sd_sL(Y, X, h=1e-5))
        print(f"        {lab:>8}: s to Y <= {100*yc:.2f}%")

    # is the X-dependence weak?  spread of the crossing over X in [0,5%]
    vals = [cross[k] for k in ("L exact", "L to 1%", "L to 2%", "L to 5%")]
    print(f"\n      crossing spread over X in [0,5%]: "
          f"{100*max(vals):.2f}% - {100*min(vals):.2f}% = {100*(max(vals)-min(vals)):.2f} pts"
          f"  -> {'WEAK' if (max(vals)-min(vals)) < 0.01 else 'NOT weak'}")

    # does X enter exactly like a quadrature calibration?  check hypot(Y,X)=const
    hyps = [np.hypot(cross[k], x) for k, x in
            (("L exact", 0.0), ("L to 1%", 0.01), ("L to 2%", 0.02), ("L to 5%", 0.05))]
    print(f"      hypot(Y*,X) across X: "
          + ", ".join(f"{100*h:.1f}%" for h in hyps)
          + f"  (rises by {100*(hyps[-1]-hyps[0]):.2f} pts -> NOT exact quadrature)")


if __name__ == "__main__":
    main()
