"""
Task 4. Set identification: what the balance and non-negativity alone imply.

THE ALGEBRA. At an untreated steady state theta g(m) + k_e m = L, with every
flux non-negative and g <= 1:

    k_e m <= L                     =>  k_e in [0, L/m]
    theta = (L - k_e m)/g(m)       =>  theta in [ (L - k_e,hi m)/g ,  (L - k_e,lo m)/g ]
    g(m) = m/(K_m + m)             =>  upper end = (L - k_e,lo m)(K_m + m)/m

so the one-sided bound is

    theta_max(data) = (L - k_e,lo m) (K_m,hi + m) / m ,

increasing in K_m. NON-NEGATIVITY ALONE (k_e,lo = 0) THEREFORE BOUNDS THETA
ONLY IF K_m IS BOUNDED; with K_m free it does not bound theta at all. The
lower end is 0 whenever k_e,hi m >= L, which at 5% of normal it is for any
k_e above L/m = 5.97 /d, i.e. 1.2 times the nominal 5.

Under the dose train each draw gives its own version,

    theta <= (L_i - b - k_e,lo m_i)(K_m + m_i)/m_i - Delta_i ,

and the record identifies K_m to about 6% (Task 1), which is what can supply
K_m,hi. Part B evaluates this on noiseless records, part C on noisy ones as a
constrained profile: maximise theta subject to the data, k_e >= 0 and NO k_e
prior, and compare with n6's prior-driven 24%.
"""
import numpy as np

import core as C
from plotstyle import plt, C1, C2, C3, C4, INK, MUTED, finish, save, plain

KE_LO_95 = C.KE * np.exp(-1.96 * C.KE_PRIOR)      # 95% lower end of the 0.8 prior


def theta_max(L, m, Km_hi, ke_lo):
    return (L - ke_lo * m) * (Km_hi + m) / m


def main():
    print("=" * 92)
    print("  TASK 4. SET IDENTIFICATION")
    print("=" * 92)
    th = C.TH_5
    L, m = C.N4.anchor(th, C.KM0, C.KE, "A")
    print(f"\n  A. Untreated steady state of the 5% patient: L = {L:.2f}, m = {m:.2f}, "
          f"L/m = {L/m:.3f} /d")
    print(f"     theta_max = (L - k_e,lo m)(K_m,hi + m)/m, as % of normal (truth 5%)")
    kes = (("k_e >= 0 (non-negativity only)", 0.0),
           (f"k_e >= {KE_LO_95:.2f} (95% low end of the 0.8 prior)", KE_LO_95),
           (f"k_e >= {C.KE*np.exp(-C.KE_PRIOR):.2f} (1 sd low end)", C.KE * np.exp(-C.KE_PRIOR)))
    kms = (100.0, 107.0, 150.0, 200.0, 500.0, 1e4)
    print(f"     {'':52}" + "".join(f"{'K_m<=' + format(k, 'g'):>11}" for k in kms))
    for lab, kl in kes:
        print(f"     {lab:52}" + "".join(
            f"{100*theta_max(L, m, k, kl)/C.TH0:10.1f}%" for k in kms))
    print("     K_m unbounded: no finite bound in any row.")
    print(f"     Lower end: 0 for any k_e,hi >= L/m = {L/m:.2f}; the 0.8 prior's 95% upper end "
          f"is {C.KE*np.exp(1.96*C.KE_PRIOR):.1f}, so the lower bound is 0.")

    # ---- B: per-draw bounds on noiseless dose-train records, K_m known
    print("\n  B. Under the dose train, noiseless, K_m known: theta <= min_i [(L_i - k_e,lo m_i)/g(m_i) - Delta_i]")
    for lab, kl in kes[:2]:
        vals = []
        for s in C.SEEDS:
            r = C.record(th, seed=s)
            g = r["m"] / (C.KM0 + r["m"])
            b = (r["loads"] - kl * r["m"]) / g - r["inc"]
            vals.append((b.min(), r["m"][np.argmin(b)], r["inc"][np.argmin(b)], b.max()))
        v = np.array(vals)
        print(f"     {lab:52}: tightest draw {np.median(v[:,0]):6.0f} = "
              f"{100*np.median(v[:,0])/C.TH0:4.1f}% (at m = {np.median(v[:,1]):.1f}, Delta = "
              f"{np.median(v[:,2]):.0f}); loosest {100*np.median(v[:,3])/C.TH0:.1f}%")
    print("     The tightest draws are the high-Delta ones, where the bound is L/g - Delta: a")
    print("     difference of two numbers near 6500 whose own error is a few hundred (5% on m")
    print("     alone moves L/g by 5%), so on real data the per-draw bound there is noise.")

    # ---- C: the constrained profile on noisy records
    print("\n  C. Noisy records: the 95% upper end of the profile, by what closes it")
    lt = np.log(th) + np.linspace(-3.0, 4.0, 71)
    cases = (("n6 footing: k_e prior 0.8, b prior 5%", dict()),
             ("NO k_e prior, k_e > 0, b prior 5%", dict(ke_prior=None)),
             ("NO k_e prior, K_m prior 0.07 as well", dict(ke_prior=None, km_prior=0.07)),
             ("k_e prior 0.8, K_m prior 0.07", dict(km_prior=0.07)))
    res = {}
    for lab, kw in cases:
        his, kemin = [], []
        for s in C.SEEDS:
            r = C.record(th, seed=s)
            y, sg = C.noisy(r, s)
            c, x = C.profile(r, y, sg, lt, **kw)
            c -= c.min()
            ok = c <= 1.92
            j = np.where(ok)[0].max()
            his.append(np.exp(lt[j]) if j < lt.size - 1 else np.inf)
            kemin.append(np.exp(x[j, 2]))
        his = np.array(his)
        fin = his[np.isfinite(his)]
        res[lab] = his
        print(f"     {lab:42}: median {np.median(his):7.0f} = {100*np.median(his)/C.TH0:5.1f}%"
              f" of normal; {np.sum(~np.isfinite(his))} of 24 unbounded on the grid; "
              f"k_e at the bound {np.median(kemin):.3f}")

    # noiseless versions of the same, the expected profile
    print("     noiseless record (expected profile), same four footings:")
    r = C.record(th, seed=0)
    y0 = (r["theta"] + r["inc"]) * r["m"] / (r["Km"] + r["m"]) + r["ke"] * r["m"]
    for lab, kw in cases:
        c, x = C.profile(r, y0, C.sigmas(r), lt, **kw)
        c -= c.min()
        j = np.where(c <= 1.92)[0].max()
        tag = ("beyond grid" if j == lt.size - 1 else
               f"{np.exp(lt[j]):6.0f} = {100*np.exp(lt[j])/C.TH0:4.1f}%")
        print(f"       {lab:42}: {tag}; k_e there {np.exp(x[j,2]):.3f}, K_m {np.exp(x[j,0]):.1f}")

    # ---- figure: the analytic bound against K_m,hi, two k_e floors
    fig, ax = plt.subplots(figsize=(5.6, 3.9))
    K = np.geomspace(60, 2000, 200)
    for (lab, kl), col in zip(kes[:2], (C1, C2)):
        ax.plot(K, 100 * theta_max(L, m, K, kl) / C.TH0, color=col, lw=2)
    ax.text(1100, 100 * theta_max(L, m, 1100, 0) / C.TH0 * 1.08, "k_e ≥ 0", color=INK,
            fontsize=8, ha="center")
    ax.text(1100, 100 * theta_max(L, m, 1100, KE_LO_95) / C.TH0 * 0.80,
            f"k_e ≥ {KE_LO_95:.2f}\n(prior's 95% low end)", color=INK, fontsize=8, ha="center")
    ax.axhline(100 * np.median(res[cases[0][0]]) / C.TH0, color=MUTED, ls="--", lw=1)
    ax.text(62, 100 * np.median(res[cases[0][0]]) / C.TH0 * 0.86,
            f"n6 profile bound, {100*np.median(res[cases[0][0]])/C.TH0:.0f}%", fontsize=8, color=MUTED)
    ax.axvline(100, color=MUTED, lw=0.8, ls=":")
    ax.text(103, 400, "true K_m", fontsize=8, color=MUTED)
    ax.axhline(5, color=C3, lw=1)
    ax.text(62, 5.6, "truth, 5%", fontsize=8, color=INK)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"assumed upper bound on $K_m$")
    ax.set_ylabel(r"$\theta_{max}$, % of normal")
    ax.set_title("non-negativity bound, untreated 5% patient")
    finish(ax); plain(ax)
    save(fig, "t4_set")


if __name__ == "__main__":
    main()
