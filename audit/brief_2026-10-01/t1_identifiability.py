"""
Task 1. Identifiability analysis on the repository's footing.

  1a  sd(ln theta_hat) from 1% to 1000% of normal: the U-shape and its minimum
  1b  the mechanism: correlations of ln theta with each nuisance, and the
      variance each nuisance carries, against capacity
  1c  the analytic elasticity d ln m / d ln theta = -s/(s eps + 1 - s) against
      a Richardson central difference of the steady-state solve; eps at normal
  1d  the profile likelihood at s = 0.16 (5% of normal), its actual shape

Footing: core.py's docstring. Every number is a median over the 24 intake
records unless it says otherwise.
"""
import numpy as np
from scipy.optimize import minimize_scalar

import core as C
from plotstyle import plt, C1, C2, C3, C4, INK, MUTED, finish, save

PCTS = np.geomspace(1.0, 1000.0, 61)


def scan_curve():
    rows = []
    for pct in PCTS:
        th = C.TH0 * pct / 100
        rows.append((pct, *C.sd_theta(th)))
    return np.array(rows)


def u_minimum():
    f = lambda lp: C.sd_theta(C.TH0 * np.exp(lp) / 100)[0]
    res = minimize_scalar(f, bounds=(np.log(100), np.log(500)), method="bounded",
                          options=dict(xatol=1e-4))
    return float(np.exp(res.x)), float(res.fun)


def band_edges(times=C.MONTHLY):
    f = lambda th: C.sd_theta(th, rec_kw=dict(times=times))[0]
    return (C.N4.crossing(f, 3000.0, 300.0), C.N4.crossing(f, 6000.0, 30000.0))


def mechanism(pcts=(1, 2, 5, 10, 20, 50, 82, 100, 150, 200, 300, 627, 1000)):
    """Correlations, and sd(ln theta) with each nuisance removed, per capacity."""
    out = []
    for pct in pcts:
        th = C.TH0 * pct / 100
        cs, cond = [], {k: [] for k in ("ke", "Km", "b", "ke+b", "all")}
        for s in C.SEEDS:
            F = C.fim(C.record(th, seed=s))
            P = C.cov(F)
            cs.append([C.corr(P, 0, C.I_KE), C.corr(P, 0, C.I_KM),
                       C.corr(P, 0, C.I_B), np.sqrt(P[3, 3]), np.sqrt(P[1, 1]),
                       P[0, 3] / P[3, 3]])
            for k, fix in (("ke", [3]), ("Km", [1]), ("b", [2]), ("ke+b", [2, 3]),
                           ("all", [1, 2, 3])):
                keep = [i for i in range(4) if i not in fix]
                cond[k].append(np.sqrt(np.linalg.inv(F[np.ix_(keep, keep)])[0, 0]))
        cs = np.median(np.array(cs), axis=0)
        out.append(dict(pct=pct, s=C.share(th), sd=C.sd_theta(th)[0],
                        c_ke=cs[0], c_km=cs[1], c_b=cs[2], post_ke=cs[3],
                        post_km=cs[4], beta=cs[5],
                        **{f"sd_{k}": float(np.median(v)) for k, v in cond.items()}))
    return out


def elasticity_check():
    """Analytic -s/(s eps + 1 - s) against a Richardson difference of the solve."""
    rows, worst = [], 0.0
    for pct in (1, 5, 10, 20, 50, 100, 200, 500, 1000):
        th = C.TH0 * pct / 100
        L, m = C.N4.anchor(th, C.KM0, C.KE, "A")
        s, eps = C.share(th), C.eps_of(m)
        ana = -s / (s * eps + 1 - s)
        f = lambda lt: np.log(C.m_ss(np.exp(lt), C.KM0, C.KE, L))
        d = [(f(np.log(th) + h) - f(np.log(th) - h)) / (2 * h) for h in (1e-3, 5e-4)]
        num = (4 * d[1] - d[0]) / 3
        worst = max(worst, abs(ana - num))
        #  eps recovered from the numerical elasticity and the share alone
        eps_fit = (s / abs(num) - (1 - s)) / s
        rows.append((pct, s, m, eps, ana, num, eps_fit, (1 - s) / s))
    # the same identity per draw under the dose train: the share and eps are
    # then each draw's own, and the increment enters the denominator
    worst_d = 0.0
    for th in (150.0, 3000.0):
        r = C.record(th, seed=0)
        m, inc, L = r["m"], r["inc"], r["loads"]
        g = m / (C.KM0 + m)
        s_i, d_i, e_i = th * g / L, inc * g / L, C.eps_of(m)
        ana = -s_i / ((s_i + d_i) * e_i + 1 - s_i - d_i)
        f = lambda lt: np.log(C.m_ss(np.exp(lt), C.KM0, C.KE, L, inc))
        num = (f(np.log(th) + 1e-5) - f(np.log(th) - 1e-5)) / 2e-5
        worst_d = max(worst_d, float(np.abs(ana - num).max()))
    return rows, worst, worst_d


def profiles(theta=C.TH_5, seeds=C.SEEDS, span=(-8.0, 3.5), n=116):
    lt = np.log(theta) + np.linspace(*span, n)
    curves, his, los, mles, lo_open = [], [], [], [], 0
    for s in seeds:
        r = C.record(theta, seed=s)
        y, sg = C.noisy(r, s)
        c, _ = C.profile(r, y, sg, lt)
        c = c - c.min()
        curves.append(c)
        inside = lt[c <= 1.92]
        his.append(np.exp(inside.max()))
        los.append(np.exp(inside.min()))
        lo_open += inside.min() <= lt[0] + 1e-12
        mles.append(np.exp(lt[np.argmin(c)]))
    r = C.record(theta, seed=0)
    y0 = (r["theta"] + r["inc"]) * r["m"] / (r["Km"] + r["m"]) + r["ke"] * r["m"]
    c0, x0 = C.profile(r, y0, C.sigmas(r), lt)
    return dict(lt=lt, curves=np.array(curves), his=np.array(his),
                los=np.array(los), mles=np.array(mles), lo_open=lo_open,
                expected=c0 - c0.min(), ke_path=np.exp(x0[:, 2]),
                km_path=np.exp(x0[:, 0]), b_path=x0[:, 1], L=r["L"])


def main():
    print("=" * 92)
    print("  TASK 1. IDENTIFIABILITY ON THE REPOSITORY'S FOOTING")
    print("=" * 92)

    # ---- 1a
    curve = scan_curve()
    pmin, sdmin = u_minimum()
    lo, hi = band_edges()
    lo_d, _ = (C.N4.crossing(lambda th: C.sd_theta(th, rec_kw=dict(times=C.DENSE))[0],
                             3000.0, 100.0), None)
    print("\n  1a. sd(ln theta_hat), median [min, max] over 24 records")
    for pct, med, mn, mx in curve[::6]:
        print(f"      {pct:7.1f}% {med:9.4f} [{mn:8.4f}, {mx:8.4f}]")
    print(f"      minimum {sdmin:.4f} at {pmin:.0f}% of normal; band where sd <= 0.20: "
          f"{100*lo/C.TH0:.1f}% to {100*hi/C.TH0:.1f}%")
    print(f"      no therapy (D = 0): 5% -> {C.sd_theta(C.TH_5, rec_kw=dict(D=0))[0]:.3f}, "
          f"normal -> {C.sd_theta(C.TH0, rec_kw=dict(D=0))[0]:.3f};"
          f" constant load and no therapy: {C.sd_theta(C.TH_5, rec_kw=dict(D=0, load_cv=0))[0]}")
    print(f"      per-draw ln-m noise instead of n4's single sigma: 5% -> "
          f"{C.sd_theta(C.TH_5, noise='perdraw')[0]:.4f}, normal -> "
          f"{C.sd_theta(C.TH0, noise='perdraw')[0]:.4f}")

    # ---- 1b
    mech = mechanism()
    print("\n  1b. Mechanism. Correlations from inv(FIM), and sd(ln theta) with a nuisance known")
    print(f"      {'pct':>5} {'s':>6} {'(1-s)/s':>8} {'beta':>7} {'c(th,ke)':>9} {'c(th,Km)':>9} "
          f"{'c(th,b)':>8} {'sd':>8} {'|ke':>7} {'|Km':>7} {'|b':>7} {'|ke,b':>7} {'|all':>7}"
          f" {'ke share of var':>16}")
    for d in mech:
        print(f"      {d['pct']:5.0f} {d['s']:6.3f} {(1-d['s'])/d['s']:8.3f} {d['beta']:7.3f} "
              f"{d['c_ke']:9.4f} {d['c_km']:9.4f} {d['c_b']:8.4f} {d['sd']:8.4f} "
              f"{d['sd_ke']:7.4f} {d['sd_Km']:7.4f} {d['sd_b']:7.4f} {d['sd_ke+b']:7.4f} "
              f"{d['sd_all']:7.4f} {1 - (d['sd_ke']/d['sd'])**2:15.1%}")
    print("      beta = Cov(ln th, ln ke)/Var(ln ke) from the posterior: the effective")
    print("      amplification, to set against the brief's (1-s)/s.")

    # ---- 1c
    rows, worst, worst_d = elasticity_check()
    print("\n  1c. d ln m / d ln theta, analytic -s/(s eps + 1 - s) against Richardson")
    print(f"      {'pct':>5} {'s':>6} {'m':>7} {'m/Km':>6} {'eps':>6} {'analytic':>9} "
          f"{'numeric':>9} {'eps from numeric':>17}")
    for pct, s, m, eps, ana, num, ef, amp in rows:
        print(f"      {pct:5.0f} {s:6.3f} {m:7.2f} {m/C.KM0:6.3f} {eps:6.3f} {ana:9.5f} "
              f"{num:9.5f} {ef:17.5f}")
    print(f"      worst |analytic - numeric| untreated {worst:.2e}; per draw under the dose "
          f"train, with the increment in the denominator, {worst_d:.2e}")

    # ---- 1d
    pr = profiles()
    crb = C.sd_theta(C.TH_5)[0]
    print(f"\n  1d. Profile likelihood at theta = 150 (s = {C.share(C.TH_5):.3f}), 24 records,")
    print(f"      grid ln theta in ln150 + [-8, +3.5], i.e. {np.exp(pr['lt'][0]):.3f} to "
          f"{np.exp(pr['lt'][-1]):.0f}")
    print(f"      CRB sd {crb:.3f} implies a 95% factor of {np.exp(2*1.96*crb):,.0f}")
    print(f"      95% upper end: median {np.median(pr['his']):.0f} "
          f"({100*np.median(pr['his'])/C.TH0:.1f}% of normal), range "
          f"{pr['his'].min():.0f} to {pr['his'].max():.0f}")
    print(f"      lower end reaches the grid floor theta = {np.exp(pr['lt'][0]):.3f} in "
          f"{pr['lo_open']} of {len(pr['his'])} records")
    print(f"      MLE: median {np.median(pr['mles']):.1f}, range {pr['mles'].min():.2f} "
          f"to {pr['mles'].max():.0f}; {np.mean(pr['mles'] < 15):.0%} of records put it "
          f"below 15 (a tenth of the truth)")
    ex = pr["expected"]
    lt = pr["lt"]
    print("      noiseless record (the expected profile), drop from the minimum:")
    for f in (0.01, 0.1, 0.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0):
        j = np.argmin(np.abs(lt - np.log(C.TH_5 * f)))
        print(f"        theta = {f:5.2f} x truth: {ex[j]:7.3f}   k_e on the profile "
              f"{pr['ke_path'][j]:.3f}, K_m {pr['km_path'][j]:.1f}, b/L {pr['b_path'][j]/pr['L']:+.3f}")
    up = lt[(ex <= 1.92)].max()
    print(f"      expected-profile 95% upper end {np.exp(up):.0f} ({100*np.exp(up)/C.TH0:.1f}%"
          f" of normal); below the truth the expected profile never exceeds "
          f"{ex[lt < np.log(C.TH_5)].max():.3f}")

    # ---- figure
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.0))
    a = ax[0]
    a.fill_between(curve[:, 0], curve[:, 2], curve[:, 3], color=C1, alpha=0.18, lw=0)
    a.plot(curve[:, 0], curve[:, 1], color=C1, lw=2, label="median over 24 records")
    a.axhline(0.20, color=INK, lw=1, ls="--")
    a.text(1.1, 0.215, "target 0.20", color=INK, fontsize=8)
    a.axvspan(100 * lo / C.TH0, 100 * hi / C.TH0, color=C3, alpha=0.10, lw=0)
    a.text(np.sqrt(lo * hi) / C.TH0 * 100, 6, f"band\n{100*lo/C.TH0:.0f}–{100*hi/C.TH0:.0f}%",
           ha="center", fontsize=8, color=MUTED)
    a.plot([5], [crb], "o", color=C2, ms=6, zorder=5)
    a.annotate(f"5%: {crb:.2f}", (5, crb), (9, 5.5), fontsize=8, color=INK,
               arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
    a.plot([pmin], [sdmin], "o", color=C1, ms=6, zorder=5)
    a.annotate(f"min {sdmin:.3f}\nat {pmin:.0f}%", (pmin, sdmin), (12, 0.125),
               fontsize=8, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("capacity, % of normal"); a.set_ylabel(r"sd(ln $\hat\theta$)")
    a.set_title("(a) monthly record, 25 draws, therapy on")
    finish(a)

    b = ax[1]
    P = np.array([d["pct"] for d in mech])
    b.plot(P, [d["c_ke"] for d in mech], "-o", color=C1, ms=4, lw=2)
    b.plot(P, [d["c_km"] for d in mech], "-s", color=C2, ms=4, lw=2)
    b.plot(P, [d["c_b"] for d in mech], "-^", color=C3, ms=4, lw=2)
    b.text(1.2, -0.90, r"corr(ln$\theta$, ln$k_e$)", color=C1, fontsize=8.5)
    b.text(160, 0.78, r"corr(ln$\theta$, ln$K_m$)", color=C2, fontsize=8.5, ha="right")
    b.text(1.2, -0.04, r"corr(ln$\theta$, $b$)", color=C3, fontsize=8.5)
    b.axhline(-1, color=MUTED, lw=0.6, ls=":")
    b.set_xscale("log"); b.set_ylim(-1.05, 1.05)
    b.set_xlabel("capacity, % of normal"); b.set_ylabel("posterior correlation")
    b.set_title("(b) what ln θ is confounded with")
    finish(b)

    c = ax[2]
    x = np.exp(pr["lt"]) / C.TH_5
    for k, cv in enumerate(pr["curves"]):
        c.plot(x, cv, color=C1, lw=0.7, alpha=0.35,
               label="24 noisy records" if k == 0 else None)
    c.plot(x, ex, color=INK, lw=2, label="noiseless record")
    q = 0.5 * (np.log(x) / crb) ** 2
    c.plot(x, q, color=C2, lw=1.5, ls="--", label=f"CRB quadratic, sd {crb:.2f}")
    c.axhline(1.92, color=MUTED, lw=0.8, ls=":")
    c.text(0.012, 2.15, "95% cut (1.92)", fontsize=8, color=MUTED)
    c.axvline(np.median(pr["his"]) / C.TH_5, color=C3, lw=1.2)
    c.text(np.median(pr["his"]) / C.TH_5 * 0.7, 8,
           f"95% upper end\nmedian {np.median(pr['his']):.0f}\n= {100*np.median(pr['his'])/C.TH0:.0f}% of normal",
           fontsize=8, color=INK, ha="right")
    c.set_xscale("log"); c.set_ylim(-0.3, 14)
    c.set_xlabel(r"$\theta$ / true $\theta$  (truth 150 = 5% of normal)")
    c.set_ylabel(r"$-\Delta\log$ profile likelihood")
    c.set_title("(c) profile at s = 0.16")
    c.legend(fontsize=7.5, loc="upper left", frameon=False)
    finish(c)
    save(fig, "t1_identifiability")


if __name__ == "__main__":
    main()
