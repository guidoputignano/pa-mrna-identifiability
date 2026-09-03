"""
Model A: saturating enzyme, linear alternative clearance.

    dm/dt = s(t) - theta*g(m) - ke*m ,   g(m) = m/(Km+m)

Question answered: can theta be identified from clinical data?

Produces
    fig_A1_channels          the two clearance channels and the linear regime
    fig_A2_conditioning      condition number against excursion of m
    fig_A3_indistinguishable two patients with identical measurements
    fig_A4_population        deficiency of the population system against N
    fig_A5_load_error        amplification of misspecification in s
    tab_A_rank               rank against number of distinct levels

Run:  python src/model_A/run.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq, least_squares

from common.models import (ParamsA, g, clearance_A, steady_state,
                           sensitivity_A, condition_number)
from common.style import save, write_table, GREEN, RED, BLUE, AMBER, GREY, MUTED


P = ParamsA()


# =====================================================================
def fig_channels():
    """The two clearance channels. In the low-m band only their sum is seen."""
    m = np.linspace(0, 600, 600)
    enz = P.theta * g(m, P.Km)
    alt = P.ke * m
    tot = enz + alt

    fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.1))

    ax[0].axvspan(0, 40, color=BLUE, alpha=0.06)
    ax[0].plot(m, enz, color=GREEN, label=r"enzymatic  $\theta\,g(m)$")
    ax[0].plot(m, alt, color=AMBER, label=r"alternative  $k_e m$")
    ax[0].axhline(P.theta, color=GREEN, ls=":", lw=0.9)
    ax[0].text(560, P.theta + 60, r"$\theta$", color=GREEN, fontsize=9)
    ax[0].text(6, 2700, r"$m \ll K_m$", color=BLUE, fontsize=8)
    ax[0].set_xlabel(r"metabolite $m$")
    ax[0].set_ylabel("elimination flux")
    ax[0].set_title("the two channels")
    ax[0].legend(loc="lower right")

    # two parameter sets that agree in the operating regime and differ outside
    from scipy.optimize import brentq as _bq
    m_op = 10.0
    keB = _bq(lambda k: clearance_A(m_op, ParamsA(1000.0, P.Km, k))
              - clearance_A(m_op, P), 1e-6, 1e4)
    P2 = ParamsA(1000.0, P.Km, keB)

    ax[1].axvspan(0, 25, color=BLUE, alpha=0.06)
    ax[1].plot(m, clearance_A(m, P), color=GREEN,
               label=r"A: $\theta=3000,\ k_e=5.0$")
    ax[1].plot(m, clearance_A(m, P2), color=RED, ls="--",
               label=rf"B: $\theta=1000,\ k_e={keB:.1f}$")
    ax[1].plot([m_op], [clearance_A(m_op, P)], "o", ms=5, color=MUTED, zorder=4)
    ax[1].text(11, clearance_A(m_op, P) - 400,
               "identical here", color=MUTED, fontsize=7.5)
    ax[1].text(3, 3400, "operating\nregime", color=BLUE, fontsize=7.5)
    ax[1].set_xlim(0, 300)
    ax[1].set_ylim(0, 4200)
    ax[1].set_xlabel(r"metabolite $m$")
    ax[1].set_title("two patients, same total clearance where they live")
    ax[1].legend(loc="lower right")

    save(fig, "fig_A1_channels")


# =====================================================================
def fig_conditioning():
    """Condition number against the relative excursion of m."""
    spreads = np.logspace(-2, 1.3, 60)
    conds = []
    m0 = 10.0
    for sp in spreads:
        levels = m0 * np.array([1.0, 1 + sp, 1 + 2 * sp])
        conds.append(condition_number(sensitivity_A, P, levels))
    conds = np.array(conds)

    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    ax.loglog(spreads * 100, conds, color=BLUE)
    ref = conds[0] * (spreads / spreads[0]) ** -2
    ax.loglog(spreads * 100, ref, color=GREY, ls="--", lw=1.0,
              label=r"slope $-2$")
    ax.axvspan(1, 10, color=RED, alpha=0.07)
    ax.text(1.3, conds[-1] * 3, "clinically\nstable patient", color=RED, fontsize=7.5)
    ax.set_xlabel("relative excursion of $m$  [%]")
    ax.set_ylabel("condition number")
    ax.set_title("conditioning degrades quadratically as the patient is stabilised")
    ax.legend()
    ax.grid(True, which="both", alpha=0.4)
    save(fig, "fig_A2_conditioning")


# =====================================================================
def fig_indistinguishable():
    """Two patients with different biology and identical measurements."""
    thA, keA = 3000.0, 5.0
    pA = ParamsA(thA, P.Km, keA)
    s = clearance_A(10.0, pA)
    thB = 1000.0
    keB = brentq(lambda k: clearance_A(10.0, ParamsA(thB, P.Km, k)) - s, 1e-6, 1e4)
    pB = ParamsA(thB, P.Km, keB)

    deltas = np.linspace(0, 6000, 120)
    mA = [steady_state(clearance_A, ParamsA(thA + d, P.Km, keA), s) for d in deltas]
    mB = [steady_state(clearance_A, ParamsA(thB + d, P.Km, keB), s) for d in deltas]
    mA, mB = np.array(mA), np.array(mB)

    fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.1))

    ax[0].plot(deltas, mA, color=GREEN, label=r"A: $\theta=3000,\ k_e=5.0$")
    ax[0].plot(deltas, mB, color=RED, ls="--", label=r"B: $\theta=1000,\ k_e=23.2$")
    ax[0].set_xlabel(r"therapeutic increment $\Delta$")
    ax[0].set_ylabel(r"steady-state $m$")
    ax[0].set_title("response to a selective input on the enzymatic channel")
    ax[0].legend()

    rel = 100 * np.abs(mA - mB) / mA
    ax[1].plot(deltas, rel, color=BLUE)
    ax[1].axhline(5, color=GREY, ls=":", lw=1.0)
    ax[1].text(150, 5.3, "typical assay variability", color=GREY, fontsize=7.5)
    ax[1].set_xlabel(r"therapeutic increment $\Delta$")
    ax[1].set_ylabel("relative difference  [%]")
    ax[1].set_ylim(0, 8)
    ax[1].set_title("separation stays below assay noise")

    save(fig, "fig_A3_indistinguishable")
    return pA, pB, s


# =====================================================================
def fig_population(N_max=30, seed=0):
    """Deficiency of the steady-state population system, and anchor recovery."""
    Ns = np.arange(1, N_max + 1)
    unknowns = Ns + 2
    equations = Ns

    fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.1))

    ax[0].plot(Ns, unknowns, color=RED, label="unknowns  $N+2$")
    ax[0].plot(Ns, equations, color=GREEN, label="equations  $N$")
    ax[0].fill_between(Ns, equations, unknowns, color=RED, alpha=0.08)
    ax[0].annotate("deficiency 2\nfor every $N$", xy=(20, 21), xytext=(9, 27),
                   color=RED, fontsize=8,
                   arrowprops=dict(arrowstyle="->", color=RED, lw=0.9))
    ax[0].set_xlabel("cohort size $N$")
    ax[0].set_ylabel("count")
    ax[0].set_title("enlarging the cohort does not help")
    ax[0].legend(loc="upper left")

    # anchored recovery
    rng = np.random.default_rng(seed)
    N = 12
    thetas = np.concatenate([[0.0], rng.uniform(300, 4000, N - 1)])
    loads = rng.uniform(200, 500, N)
    ms = np.array([steady_state(clearance_A, ParamsA(t, P.Km, P.ke), s)
                   for t, s in zip(thetas, loads)])

    ke_hat = loads[0] / ms[0]                      # anchor (i): null genotype
    i = 1                                          # anchor (ii): two loads
    s2 = loads[i] * 1.6
    m2 = steady_state(clearance_A, ParamsA(thetas[i], P.Km, P.ke), s2)

    def resid(v):
        th, km = v
        return [loads[i] - ke_hat * ms[i] - th * ms[i] / (km + ms[i]),
                s2 - ke_hat * m2 - th * m2 / (km + m2)]

    th_hat, Km_hat = least_squares(resid, [1000.0, 50.0],
                                   bounds=([0, 1e-3], [1e6, 1e5])).x
    rec = [(loads[j] - ke_hat * ms[j]) * (Km_hat + ms[j]) / ms[j]
           for j in range(1, N)]

    ax[1].scatter(thetas[1:], rec, s=26, color=GREEN, zorder=3,
                  label="recovered")
    lim = [0, 4200]
    ax[1].plot(lim, lim, color=GREY, ls="--", lw=1.0, label="identity")
    ax[1].set_xlim(lim); ax[1].set_ylim(lim)
    ax[1].set_xlabel(r"true $\theta_i$")
    ax[1].set_ylabel(r"recovered $\theta_i$")
    ax[1].set_title("with a null genotype and one repeated observation")
    ax[1].legend(loc="upper left")

    save(fig, "fig_A4_population")
    print(f"    ke recovered {ke_hat:.4f} (true {P.ke}), "
          f"Km recovered {Km_hat:.3f} (true {P.Km})")


# =====================================================================
def fig_load_error():
    """How an error in the assumed substrate load propagates into theta."""
    theta, Km, ke = 1500.0, P.Km, P.ke
    p = ParamsA(theta, Km, ke)
    m = 20.0
    s = clearance_A(m, p)
    eps = np.linspace(-0.3, 0.3, 200)
    th_hat = np.array([(s * (1 + e) - ke * m) * (Km + m) / m for e in eps])
    rel = 100 * (th_hat - theta) / theta

    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    ax.plot(eps * 100, rel, color=BLUE)
    ax.plot(eps * 100, eps * 100, color=GREY, ls="--", lw=1.0,
            label="no amplification")
    ax.axhline(0, color=GREY, lw=0.6)
    ax.axvline(0, color=GREY, lw=0.6)
    slope = np.polyfit(eps * 100, rel, 1)[0]
    ax.set_xlabel(r"error in assumed load $s$  [%]")
    ax.set_ylabel(r"resulting error in $\theta$  [%]")
    ax.set_title(f"amplification factor  $\\approx$ {slope:.2f}")
    ax.legend()
    save(fig, "fig_A5_load_error")
    return slope


# =====================================================================
def table_rank():
    rows = []
    for levels in ([10.0], [10.0, 12.0], [10.0, 12.0, 14.0],
                   [10.0, 50.0, 200.0], [10.0, 12.0, 14.0, 16.0]):
        S = np.array([sensitivity_A(m, P) for m in levels])
        rows.append([len(levels), np.linalg.matrix_rank(S, tol=1e-10),
                     "yes" if np.linalg.matrix_rank(S, tol=1e-10) == 3 else "no"])
    write_table("tab_A_rank",
                ["distinct levels of $m$", "rank of $S$", "all three identifiable"],
                rows,
                note="three distinct levels are necessary and generically sufficient")


if __name__ == "__main__":
    print("MODEL A")
    fig_channels()
    fig_conditioning()
    fig_indistinguishable()
    fig_population()
    slope = fig_load_error()
    table_rank()
    print(f"  load-error amplification factor: {slope:.3f}")
