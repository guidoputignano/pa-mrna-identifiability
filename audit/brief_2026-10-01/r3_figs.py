"""
Round-3 figures:
  r3_flux_contour   item 11: sd(ln theta) over (reading precision Y, calibration Z), 0.20 contour
  r3_ke_sampling    item 12: sd(ln theta) over (sd(ln k_e), draws/year), 0.20 & 0.30 contours
  r3_asymptote      item 13: sd(ln theta) vs n with k_e exact -- no floor, ~1/sqrt(n)
  r3_excursion      item 15: sd(ln theta) vs load multiple, admissible region shaded
"""
import numpy as np
import matplotlib.pyplot as plt

import plotstyle as P
import core as C
import t5_value as T5
import t14_ke_sampling as T14
import dyn as D
import t7_design as T7

SEEDS8 = range(8)


def sd_ke(ke_sd, n):
    def one(s):
        return C.sd_th(C.fim(T14.record_n(n, s), noise="perdraw", ke_prior=ke_sd))
    return C.over_seeds(one, seeds=SEEDS8)[0]


# ---------------------------------------------------------------- figure 1
def fig_flux():
    Ys = np.linspace(0.03, 0.30, 46)
    Zs = np.linspace(0.0, 0.26, 44)
    Z = np.array([[T5.sd5(flux_Y=float(np.hypot(y, z))) for y in Ys] for z in Zs])
    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    im = ax.pcolormesh(Ys * 100, Zs * 100, np.clip(Z, 0, 0.6), cmap="Blues_r",
                       shading="gouraud", vmin=0.05, vmax=0.45)
    cs = ax.contour(Ys * 100, Zs * 100, Z, levels=[0.20], colors=[P.INK], linewidths=2)
    ax.clabel(cs, fmt={0.20: "sd = 0.20"}, fontsize=9)
    ax.contour(Ys * 100, Zs * 100, Z, levels=[0.15, 0.30], colors=[P.MUTED],
               linewidths=0.8, linestyles="--")
    fig.colorbar(im, ax=ax, label="sd(ln θ̂) at 5% of normal")
    ax.set_xlabel("reading precision Y  (% relative)")
    ax.set_ylabel("calibration precision on c,  Z = sd(ln c)  (%)")
    ax.set_title("Flux observable O = c·θ·g(m): target needs √(Y²+Z²) ≤ 19%")
    ax.text(21, 6, "√(Y²+Z²)=19%", color=P.INK, fontsize=8, rotation=-38)
    P.finish(ax)
    P.save(fig, "r3_flux_contour")


# ---------------------------------------------------------------- figure 2
def fig_ke_sampling():
    kes = np.geomspace(0.02, 0.8, 26)
    dpy = np.geomspace(12, 365, 24)          # draws per year
    ns = np.round(2 * dpy).astype(int)       # over two years
    Z = np.array([[sd_ke(k, n) for n in ns] for k in kes])
    fig, ax = plt.subplots(figsize=(5.8, 4.2))
    im = ax.pcolormesh(dpy, kes, np.clip(Z, 0, 1.5), cmap="Blues_r",
                       shading="gouraud", vmin=0.1, vmax=1.2)
    cs = ax.contour(dpy, kes, Z, levels=[0.20, 0.30], colors=[P.INK, P.MUTED],
                    linewidths=[2, 1])
    ax.clabel(cs, fmt={0.20: "0.20", 0.30: "0.30"}, fontsize=9)
    ax.set_xscale("log"); ax.set_yscale("log")
    # achievable / feasible box: k_e no better than ~0.2 log units, draws <= ~52/yr
    ax.axhspan(0.2, 0.8, color=P.C2, alpha=0.10)
    ax.axvspan(12, 52, color=P.C3, alpha=0.08)
    ax.axhline(0.2, color=P.C2, lw=1.2, ls=":")
    ax.text(14, 0.52, "achievable k_e precision\n(no better than ~0.2 log units)",
            color=P.C2, fontsize=7.5, va="center")
    ax.text(16, 0.028, "feasible\nsampling\n(≤ weekly)", color=P.C3, fontsize=7.5)
    fig.colorbar(im, ax=ax, label="sd(ln θ̂) at 5% of normal")
    ax.set_xlabel("blood draws per year  (two-year programme)")
    ax.set_ylabel("external precision on k_e,  sd(ln k_e)")
    ax.set_title("No feasible point on the 0.20 contour")
    P.plain(ax); P.finish(ax)
    P.save(fig, "r3_ke_sampling")


# ---------------------------------------------------------------- figure 3
def fig_asymptote():
    ns = np.unique(np.round(np.geomspace(25, 100000, 22)).astype(int))

    def sd(n, bX=None):
        kw = dict(noise="perdraw", ke_prior=1e-4)
        if bX is not None:
            kw["b_prior"] = bX * C.record(C.TH_5, seed=0)["L"]

        def one(s):
            r = C.record(C.TH_5, seed=s, times=np.linspace(0, 730, int(n)))
            return C.sd_th(C.fim(r, **kw))
        return C.over_seeds(one, seeds=SEEDS8)[0]
    a = np.array([sd(n) for n in ns])
    b = np.array([sd(n, bX=1e-7) for n in ns])
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    ax.loglog(ns, a, "-o", color=P.C1, ms=4, lw=1.8, label="load bias b at its 5% prior")
    ax.loglog(ns, b, "-s", color=P.C3, ms=4, lw=1.8, label="load measured exactly (b pinned)")
    guide = a[0] * np.sqrt(ns[0] / ns)
    ax.loglog(ns, guide, color=P.MUTED, lw=1.0, ls="--", label="1/√n guide")
    ax.axhline(0.20, color=P.C2, lw=1.2, ls=":")
    ax.text(ns[-1], 0.205, "target 0.20", color=P.C2, fontsize=8, ha="right", va="bottom")
    ax.set_xlabel("number of draws over two years")
    ax.set_ylabel("sd(ln θ̂) at 5%,  k_e exact")
    ax.set_title("With k_e exact the fall has no floor (≈ 1/√n)")
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    P.plain(ax); P.finish(ax)
    P.save(fig, "r3_asymptote")


# ---------------------------------------------------------------- figure 4
def fig_excursion():
    fs = np.array([1.177, 1.3, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0])
    tt = np.arange(0, 2.0 + 1e-9, 5 / 60 / 24)
    L = D.LTOT
    sds = []
    for f in fs:
        ex = D.Experiment([(-1, 0.5 * L), (0, 0.5 * L + (f - 1) * L)], L0=0.5 * L)
        Fd = D.fim(D.params(ex=ex), ex, tt)
        sds.append(T7.joint_record(Fd, {"ke": 0.8, "Km": T7.KM_REC, "L0": 0.01 * L}))
    sds = np.array(sds)
    fig, ax = plt.subplots(figsize=(5.8, 4.0))
    ax.axvspan(0.9, 1.177, color=P.C3, alpha=0.14)
    ax.text(1.03, sds.max() * 0.9, "admissible\n(load ≤ 1.177×)", color="#157f58",
            fontsize=8, ha="center")
    ax.plot(fs, sds, "-o", color=P.C1, ms=5, lw=1.9)
    ax.axhline(0.20, color=P.C2, lw=1.3, ls="--")
    ax.text(7.7, 0.215, "target 0.20", color=P.C2, fontsize=8, ha="right")
    # mark the crossing
    xc = np.interp(np.log(0.20), np.log(sds[::-1]), fs[::-1])
    ax.plot([xc], [0.20], "v", color=P.INK, ms=8)
    ax.text(xc, 0.255, f"crosses at\nf ≈ {xc:.1f}×", color=P.INK, fontsize=8, ha="center")
    ax.set_xlabel("sustained load excursion,  multiple of baseline load")
    ax.set_ylabel("sd(ln θ̂) at 5%  (dense transient + record, generous priors)")
    ax.set_title("Every excursion reaching 0.20 is past the safety threshold")
    ax.set_ylim(0.10, sds.max() * 1.08)
    P.finish(ax)
    P.save(fig, "r3_excursion_safety")


if __name__ == "__main__":
    fig_flux()
    fig_ke_sampling()
    fig_asymptote()
    fig_excursion()
