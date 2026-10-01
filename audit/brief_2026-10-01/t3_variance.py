"""
Task 3. Var(ln theta_hat) split into a 1/n part and an n-independent part.

THE SPLIT, exact for a Gaussian posterior (law of total variance):

    Var(ln th) = Var(ln th | ln ke)  +  beta^2 Var(ln ke),
                 beta = Cov(ln th, ln ke) / Var(ln ke)

The first term is what the record determines with k_e held: it falls with the
number of draws. The second is k_e's posterior variance carried through the
amplification beta; it is n-independent for as long as the record teaches
nothing about k_e, so Var(ln ke) stays at the prior's 0.64. Both are checked
rather than assumed: the first against 1/n, the second against a fit of
sd(n) = sqrt(floor^2 + c^2/n) over record length (the repository's g2 form).

Then the 15x sampling question: the band's lower edge with draws every 2 days
for the same two years (365 against 25), and the edge in the limit of
infinitely many draws with k_e still held only by its prior.
"""
import numpy as np
from scipy.optimize import curve_fit

import core as C
from plotstyle import plt, C1, C2, C3, C4, INK, MUTED, finish, save, plain

EDGE = 2472.0          # n4's lower band edge, 82% of normal


def split(theta, times=C.MONTHLY, ke_prior=C.KE_PRIOR):
    """Median over records of (total, given-ke, ke part, posterior Var ln ke, beta)."""
    out = []
    for s in C.SEEDS:
        P = C.cov(C.fim(C.record(theta, seed=s, times=times), ke_prior=ke_prior))
        beta = P[0, 3] / P[3, 3]
        cond = P[0, 0] - P[0, 3] ** 2 / P[3, 3]
        out.append((P[0, 0], cond, beta ** 2 * P[3, 3], P[3, 3], beta))
    return np.median(np.array(out), axis=0)


def sd_vs_n(theta, ns, ke_prior=C.KE_PRIOR):
    """Monthly records of n draws (prefixes of n4's GRID, so longer records)."""
    return np.array([C.sd_theta(theta, rec_kw=dict(times=C.GRID[:n]),
                                ke_prior=ke_prior)[0] for n in ns])


def main():
    print("=" * 92)
    print("  TASK 3. VARIANCE DECOMPOSITION")
    print("=" * 92)
    print("\n  A. The exact split, monthly 25 draws, median over 24 records")
    print(f"     {'pct':>5} {'s':>6} {'sd':>8} {'Var':>9} {'Var|ke':>9} {'beta^2 Var ke':>14} "
          f"{'n-indep share':>14} {'Var ke post':>12} {'beta':>8} {'-(1-s)/s':>9}")
    rows = {}
    for pct in (5, 10, 20, 50, 82.4, 100, 200, 627):
        th = C.TH0 * pct / 100
        tot, cond, kep, vke, beta = split(th)
        rows[pct] = (tot, cond, kep, vke, beta)
        s = C.share(th)
        print(f"     {pct:5.1f} {s:6.3f} {np.sqrt(tot):8.4f} {tot:9.5f} {cond:9.5f} {kep:14.5f} "
              f"{kep/tot:14.1%} {vke:12.4f} {beta:8.4f} {-(1-s)/s:9.4f}")
    tot, cond, kep, vke, beta = rows[82.4]
    s = C.share(EDGE)
    print(f"\n     At the band edge (theta = {EDGE:.0f}, s = {s:.3f}): {kep/tot:.1%} of the variance is "
          f"the k_e term.")
    print(f"     Implied Var(ln k_e) = (beta^2 Var ke)/beta^2 = {vke:.4f} (sd {np.sqrt(vke):.4f});"
          f" via the brief's (1-s)/s instead: {kep/((1-s)/s)**2:.4f} (sd "
          f"{np.sqrt(kep/((1-s)/s)**2):.4f}). The prior is 0.8^2 = 0.64.")

    # ---- is the conditional part really 1/n, and is the k_e part really flat?
    print("\n  B. Each term against the number of draws (monthly, longer records), theta = 2472")
    ns = np.array([6, 12, 25, 50, 100, 200, 400])
    print(f"     {'n':>5} {'sd':>8} {'Var|ke * n':>11} {'beta^2 Var ke':>14} {'Var ke post':>12}")
    terms = []
    for n in ns:
        tot, cond, kep, vke, beta = split(EDGE, times=C.GRID[:n])
        terms.append((n, tot, cond, kep, vke))
        print(f"     {n:5d} {np.sqrt(tot):8.4f} {cond*n:11.5f} {kep:14.5f} {vke:12.4f}")
    terms = np.array(terms)

    # ---- the repository's own g2 functional form
    ns_fit = np.array([8, 12, 16, 25, 36, 50, 75, 100, 150, 200, 300, 400])
    sd = sd_vs_n(EDGE, ns_fit)
    f = lambda n, fl, c: np.sqrt(fl ** 2 + c ** 2 / n)
    (fl, c), _ = curve_fit(f, ns_fit, sd, p0=(0.18, 0.4))
    print(f"\n     Fit sd(n) = sqrt(floor^2 + c^2/n) over n = 8..400: floor {fl:.4f}, c {c:.4f};"
          f" at n = 25 the floor is {fl**2/f(25, fl, c)**2:.1%} of the variance")
    print(f"     worst fit residual {np.abs(sd - f(ns_fit, fl, c)).max():.4f}")

    # ---- the 15x question
    lo25 = C.N4.crossing(lambda th: C.sd_theta(th)[0], 3000.0, 300.0)
    lo365 = C.N4.crossing(lambda th: C.sd_theta(th, rec_kw=dict(times=C.DENSE))[0], 3000.0, 300.0)

    def floor_only(th):
        # n -> infinity with k_e held at its prior: only the beta^2 * 0.64 term,
        # beta from a very long record (n4 GRID, 600 monthly draws)
        P = C.cov(C.fim(C.record(th, seed=0, times=C.GRID)))
        return np.sqrt((P[0, 3] / P[3, 3]) ** 2 * C.KE_PRIOR ** 2)
    lo_inf = C.N4.crossing(floor_only, 3000.0, 300.0)
    print("\n  C. The lower band edge against sampling effort (gate 0.20)")
    print(f"     monthly, 25 draws:                    {lo25:6.0f} = {100*lo25/C.TH0:5.1f}% of normal")
    print(f"     every 2 days, 365 draws (14.6x):      {lo365:6.0f} = {100*lo365/C.TH0:5.1f}% "
          f"(moves {100*(1-lo365/lo25):.1f}%)")
    print(f"     infinitely many, k_e at its prior:    {lo_inf:6.0f} = {100*lo_inf/C.TH0:5.1f}%")
    print(f"     the share s at those edges: {C.share(lo25):.3f}, {C.share(lo365):.3f}, "
          f"{C.share(lo_inf):.3f}")
    v25, v365 = split(lo365)[3], split(lo365, times=C.DENSE)[3]
    print(f"     dense beats the k_e-held limit because 365 draws DO teach k_e a little: at theta = "
          f"{lo365:.0f}")
    print(f"     the posterior Var(ln k_e) is {v25:.4f} after 25 monthly draws and {v365:.4f} after 365 "
          f"dense ones ({100*(1-v365/v25):.1f}% less).")
    print("     and the edge with a tighter k_e, monthly 25 draws:")
    for kp in (0.4, 0.2, 0.1):
        e = C.N4.crossing(lambda th: C.sd_theta(th, ke_prior=kp)[0], 3000.0, 30.0)
        print(f"       ln k_e prior {kp:4.2f} (factor {np.exp(kp):.2f}): edge {e:6.0f} = "
              f"{100*e/C.TH0:5.1f}% of normal")

    # ---- figure
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.9))
    a = ax[0]
    a.plot(terms[:, 0], terms[:, 2], "-o", color=C1, ms=4, lw=2)
    a.plot(terms[:, 0], terms[:, 3], "-s", color=C2, ms=4, lw=2)
    a.plot(terms[:, 0], terms[:, 1], "-", color=INK, lw=1.2)
    a.text(300, terms[-1, 2] * 1.6, "Var(ln θ | k_e)\n∝ 1/n", color=INK, fontsize=8, ha="center")
    a.text(7, terms[0, 3] * 1.25, "β² Var(ln k_e): flat", color=INK, fontsize=8)
    a.text(9, terms[1, 1] * 1.05, "total", color=INK, fontsize=8)
    a.axhline(0.04, color=MUTED, ls="--", lw=1)
    a.text(250, 0.043, "gate 0.20²", fontsize=8, color=MUTED, ha="center")
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("monthly draws n (record length grows with n)")
    a.set_ylabel(r"variance of ln $\hat\theta$")
    a.set_title(f"(a) at the band edge, θ = {EDGE:.0f} (82%)")
    finish(a); plain(a)
    b = ax[1]
    P = np.geomspace(1, 1000, 41)
    shares = []
    for pct in P:
        tot, cond, kep, vke, beta = split(C.TH0 * pct / 100)
        shares.append(kep / tot)
    b.plot(P, shares, color=C1, lw=2)
    b.axvline(82.4, color=MUTED, lw=0.8, ls=":")
    b.text(82.4 * 0.9, 0.3, f"band edge\n{rows[82.4][2]/rows[82.4][0]:.0%} n-independent",
           fontsize=8, ha="right", color=INK)
    b.set_xscale("log"); b.set_ylim(0, 1.02)
    b.set_xlabel("capacity, % of normal")
    b.set_ylabel("share of Var(ln θ) held by k_e")
    b.set_title("(b) n-independent share, monthly 25 draws")
    finish(b); plain(b)
    save(fig, "t3_variance")


if __name__ == "__main__":
    main()
