"""
Task 5. Value of additional measurements, at theta = 5% of normal.

Every measurement enters the information matrix of the repository's monthly
record (n4's footing) as an extra row or a changed prior:

  (a) k_e known to a factor F: the ln k_e prior sd becomes ln F (F is a 1-sd
      factor, the convention under which the repository's 0.8 is "a factor of
      2.2"; the 95% reading, ln F / 1.96, is printed beside it).
  (b) the load MEASURED to X% instead of recorded to 5%: the bias prior becomes
      X% of L and the per-visit load scatter in sigma becomes X%.
  (c) a direct observation of the enzymatic flux F_enz = theta g(m) to Y%.
      Physically this has to be taken BEFORE therapy: on therapy any PCC flux
      measurement sees (theta + Delta) g(m) with Delta 3 to 40 times theta, so
      the patient's own flux is not separable. It is therefore modelled as one
      untreated visit carrying a marker draw and the flux reading. Its row is
      d F_enz / d(ln th, ln Km, b, ln ke) = (theta g, -theta Km m/(Km+m)^2, 0, 0):
      zero in k_e and in b. Its sd is Y% of F_enz in quadrature with the 5%
      assay error on m carried through g, eps * 5%.
  (d) combinations.

The baseline marker draw that (c) needs is priced on its own, so that the
flux reading is not credited with what the extra draw buys.
"""
import numpy as np

import core as C
from plotstyle import plt, C1, C2, C3, C4, INK, MUTED, finish, save, plain

TH = C.TH_5


def base_record(seed, baseline=False):
    r = C.record(TH, seed=seed)
    if not baseline:
        return r
    #  one extra untreated draw at the recorded load, before the first dose
    L, m0 = r["L"], r["m_bar"]
    for k, v in (("loads", L), ("inc", 0.0), ("m", m0), ("times", -1.0)):
        r[k] = np.r_[v, r[k]]
    return r


def flux_row(r):
    m0 = r["m_bar"]
    g = m0 / (C.KM0 + m0)
    F = TH * g
    return np.array([F, -TH * C.KM0 * m0 / (C.KM0 + m0) ** 2, 0.0, 0.0]), F, C.eps_of(m0)


def sd5(ke_F=None, load_X=None, flux_Y=None, baseline=False, n_flux=1):
    def one(s):
        r = base_record(s, baseline=baseline or flux_Y is not None)
        kw = {}
        if ke_F is not None:
            kw["ke_prior"] = np.log(ke_F) if ke_F > 1 else 1e-4
        if load_X is not None:
            kw["load_frac"] = max(load_X, 1e-4)
        extra = []
        if flux_Y is not None:
            row, F, eps = flux_row(r)
            sd = F * np.hypot(flux_Y, eps * C.ASSAY_CV) / np.sqrt(n_flux)
            extra.append((row, sd))
        return C.sd_th(C.fim(r, extra=extra, **kw))
    return C.over_seeds(one)[0]


def crossing(f, lo, hi):
    """Largest measurement sd in [lo, hi] at which sd(ln theta) <= 0.20, or None."""
    if f(lo) > C.GATE:
        return None
    if f(hi) <= C.GATE:
        return hi
    for _ in range(50):
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if f(mid) <= C.GATE else (lo, mid)
    return lo


def main():
    print("=" * 92)
    print("  TASK 5. VALUE OF ADDITIONAL MEASUREMENTS AT 5% OF NORMAL")
    print("=" * 92)
    base = sd5()
    bl = sd5(baseline=True)
    print(f"\n  current design: {base:.4f};  plus one untreated baseline draw: {bl:.4f}")

    Fs = (1.05, 1.1, 1.2, 1.3, 1.5, 1.75, 2.0)
    Xs = (0.01, 0.02, 0.05, 0.10, 0.20)
    Ys = (0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
    print("\n  (a) k_e known to a factor F (1 sd);  [95% reading: F at 1.96 sd]")
    A = {}
    for F in Fs:
        A[F] = sd5(ke_F=F)
        a95 = sd5(ke_F=np.exp(np.log(F) / 1.96))
        print(f"      F = {F:4.2f}  sd(ln k_e) = {np.log(F):.4f}  ->  {A[F]:.4f}   [{a95:.4f}]")
    print(f"      exact k_e: {sd5(ke_F=1.0):.4f}")
    print("\n  (b) load measured to X% (bias prior and per-visit scatter both X%)")
    Bv = {}
    for X in Xs:
        Bv[X] = sd5(load_X=X)
        print(f"      X = {100*X:4.0f}%  ->  {Bv[X]:.4f}")
    print(f"      exact load: {sd5(load_X=0.0):.4f}")
    print("\n  (c) one direct enzymatic-flux reading at an untreated visit, Y%")
    Cv = {}
    for Y in Ys:
        Cv[Y] = sd5(flux_Y=Y)
        print(f"      Y = {100*Y:4.0f}%  ->  {Cv[Y]:.4f}   (three readings averaged: "
              f"{sd5(flux_Y=Y, n_flux=3):.4f})")

    print("\n  (d) combinations")
    combos = [
        ("k_e 1.2 + load 2%", dict(ke_F=1.2, load_X=0.02)),
        ("k_e 1.05 + load 1%", dict(ke_F=1.05, load_X=0.01)),
        ("k_e exact + load exact", dict(ke_F=1.0, load_X=0.0)),
        ("flux 20% + k_e 2.0", dict(flux_Y=0.20, ke_F=2.0)),
        ("flux 20% + k_e 1.5", dict(flux_Y=0.20, ke_F=1.5)),
        ("flux 30% + k_e 1.5", dict(flux_Y=0.30, ke_F=1.5)),
        ("flux 30% + load 2%", dict(flux_Y=0.30, load_X=0.02)),
        ("flux 20% + k_e 1.2 + load 2%", dict(flux_Y=0.20, ke_F=1.2, load_X=0.02)),
        ("flux 30% + k_e 1.2 + load 2%", dict(flux_Y=0.30, ke_F=1.2, load_X=0.02)),
    ]
    D = {}
    for lab, kw in combos:
        D[lab] = sd5(**kw)
        print(f"      {lab:32} -> {D[lab]:.4f}")

    print("\n  RANKING: the loosest measurement that alone reaches sd(ln theta) = 0.20")
    ya = crossing(lambda v: sd5(ke_F=np.exp(v)), 1e-6, 1.0)
    yb = crossing(lambda v: sd5(load_X=v), 1e-6, 0.5)
    yc = crossing(lambda v: sd5(flux_Y=v), 1e-3, 1.0)
    print(f"      (a) k_e:            {'never, not even exact' if ya is None else f'sd(ln k_e) <= {ya:.4f}'}")
    print(f"      (b) load:           {'never, not even exact' if yb is None else f'X <= {100*yb:.2f}%'}")
    print(f"      (c) enzymatic flux: {'never' if yc is None else f'Y <= {100*yc:.1f}% (one reading)'}")
    yac = crossing(lambda v: sd5(flux_Y=0.30, ke_F=np.exp(v)), 1e-6, 2.0)
    if yac is None:
        print(f"      (c)+(a): with the flux at 30%, no k_e precision reaches the gate; exact k_e "
              f"gives {sd5(flux_Y=0.30, ke_F=1.0):.4f}")
    else:
        print(f"      (c)+(a): with the flux at 30%, k_e needs sd(ln k_e) <= {yac:.4f} "
              f"(factor {np.exp(yac):.2f})")
    # elasticity of sd(ln th) to each measurement's own sd, at a practical value
    def el(f, v):
        return (np.log(f(v * 1.1)) - np.log(f(v / 1.1))) / (2 * np.log(1.1))
    print("      elasticity d ln sd(ln th) / d ln (measurement sd), at a practical value:")
    print(f"        (a) at F = 1.2:   {el(lambda v: sd5(ke_F=np.exp(v)), np.log(1.2)):.3f}")
    print(f"        (b) at X = 5%:    {el(lambda v: sd5(load_X=v), 0.05):.3f}")
    print(f"        (c) at Y = 15%:   {el(lambda v: sd5(flux_Y=v), 0.15):.3f}")

    # ---------------------------------------------------- complete deficiency
    print("\n  COMPLETE DEFICIENCY. The repository holds no patient records: every 'record' is")
    print("  simulated (README: parameter values are illustrative, not fitted). So no theta ~ 0")
    print("  record exists to read k_e from. What one WOULD give, if such a patient were at an")
    print("  untreated steady state (model A): k_e = (L - b)/m, so its error is the anchor's")
    print("  load bias (5% recorded) and assay error averaged over its own draws, and it")
    print("  transfers to another patient only through the between-patient sd omega:")
    for n_anchor in (1, 25):
        se_anchor = np.hypot(C.LOAD_FRAC, C.ASSAY_CV / np.sqrt(n_anchor))
        cells = []
        for w in (0.0, 0.1, 0.2, 0.33):
            sk = np.hypot(se_anchor, w)
            cells.append(f"w={w}: {sd5(ke_F=np.exp(sk)):.3f}")
        print(f"      anchor with {n_anchor:2d} draws (k_e to {se_anchor:.3f}): " + ", ".join(cells))
    print("  And in the parent model (models.ParamsC) a theta = 0 patient has NO steady state:")
    print("  R0 is set so the alternative route alone peaks at 500 against a baseline load of")
    print("  600. A complete-deficiency patient is not at a steady state to be read.")

    # ---------------------------------------------------------------- figure
    fig, ax = plt.subplots(1, 3, figsize=(13.0, 3.8), sharey=True)
    for a, (xs, ys, lab, col, fmt) in zip(ax, (
            ([100 * (F - 1) for F in Fs], [A[F] for F in Fs], "(a) k_e known to a factor F", C1, "F − 1, %"),
            ([100 * X for X in Xs], [Bv[X] for X in Xs], "(b) load measured to X%", C2, "X, %"),
            ([100 * Y for Y in Ys], [Cv[Y] for Y in Ys], "(c) untreated enzymatic flux to Y%", C3, "Y, %"))):
        a.plot(xs, ys, "-o", color=col, ms=4, lw=2)
        a.axhline(C.GATE, color=INK, ls="--", lw=1)
        a.axhline(base, color=MUTED, ls=":", lw=1)
        a.set_title(lab)
        a.set_xlabel(fmt)
        a.set_yscale("log")
        finish(a); plain(a)
    ax[0].text(ax[0].get_xlim()[0] + 1, base * 0.82, f"current {base:.2f}", fontsize=8, color=MUTED)
    ax[0].text(ax[0].get_xlim()[0] + 1, 0.21, "target 0.20", fontsize=8, color=INK)
    ax[0].axhline(sd5(ke_F=1.0), color=C1, lw=0.8, ls="-.")
    ax[0].text(60, sd5(ke_F=1.0) * 0.85, "k_e exact", fontsize=8, color=INK)
    ax[0].set_ylabel(r"sd(ln $\hat\theta$) at 5% of normal")
    save(fig, "t5_value")


if __name__ == "__main__":
    main()
