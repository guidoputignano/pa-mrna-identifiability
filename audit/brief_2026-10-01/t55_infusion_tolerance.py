"""
Round 20 (continued): how late can an infusion be? A specification on the delay.

mRNA-3927 is given as an intravenous infusion every 14 days (paper, 01_biology). The
clinical form of the dosing question is how late an infusion visit can be, and what a
missed infusion costs, alone and when the patient misses it because they are ill.

Model: the parent model's two-state system (src/common/models.py rhs_C, k = 1) with the
delivered increment Delta(t) added to the enzyme channel (as in t54, part B):
    dm/dt = s(t) - (theta + Delta(t)) g(m) - k c m/(K_c + m)
    dc/dt = R0 - k c m/(K_c + m) - rho0 k (1 + m/M_i) c
Fortnightly doses of size D settle the cycle (8 doses); the next is given d days late
(or skipped, a 28-day gap), then doses resume on time. Load s = 600 (baseline), or the
paper's registered illness, 1.5 x baseline for 14 days (common/params ILLNESS_FOLD and
model_C's illness window), starting on the day the dose was due.

Three criteria, all stated (the dynamic ones at k = 1, the model default, and k = 10,
near the reduced limit the paper uses elsewhere; k is ASSUMED in common/params.py):
  * trough condition (the model's stability test): theta + Delta > theta_crit_int at the
    end of the gap; slack d_q = ln(Delta_tr(14) / (theta_crit - theta)) / lambda;
  * relative peak: the metabolite's peak in the gap stays below twice the patient's usual
    on-therapy peak; slack d_2x by bisection on d. This is a convention, chosen because
    the model does not define a harmful metabolite level (the paper's m_tol = 40 uM is a
    convention on the local model's marker scale, which the parent's on-therapy cycle
    already crosses);
  * separatrix: the gap peak stays below the tipping point (unstable equilibrium) of the
    parent model at the patient's usual trough capacity, beyond which the metabolite keeps
    rising on its own until the next infusion, or below the model's validity limit
    m ~ 5.5 M_i (paper 03_dynamic: 1650 uM nominal, 825 at the worst corner), whichever is
    lower. Model-intrinsic, but it also depends on k.
The 5 % row sits 0.8 uM/d below theta_crit_int, so its trough-condition slack is sensitive
(31.5 d at theta = 150, 21 d at 145); the dynamic criteria are not.

Settings: nominal carnitine constants and the worst corner of the paper's ranges
(K_c 25, M_i 150, rho_0 1.0; R0 re-derived), booked D = 6000 and the paper's largest
admissible D = 2672 at 5 % of normal (11_scope).

Writes out/t55_infusion_tolerance.txt (via stdout) and pi_note/fig/pi_tolerance.pdf/png.
Run: python -W error t55_infusion_tolerance.py
"""
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
sys.path.insert(0, str(HERE))
from common.models import ParamsC, g  # noqa: E402
from common.params import ILLNESS_FOLD  # noqa: E402
import model_C.run as MC  # noqa: E402

D0, D_ADM, TAU = 6000.0, 2672.0, 14.0
S_BASE = 600.0
K = 1.0
TH_NORMAL = 3000.0
ILL_DAYS = MC.ILLNESS[1] - MC.ILLNESS[0]          # 14 d, the paper's illustrative window
CAPS = [(0.0, "0 %"), (30.0, "1 %"), (75.0, "2.5 %"), (150.0, "5 %")]


def shape(worst):
    if not worst:
        return MC.PC
    q = replace(MC.PC, Kc=25.0, Mi=150.0, rho0=1.0)
    return replace(q, R0=MC.R0_for_peak(q))


def thc_of(p):
    return MC.theta_crit_intrinsic(p, MC.S_BASE)


def trough(tau, thalf, D):
    x = np.exp(-np.log(2.0) / thalf * tau)
    return D * (tau / TAU) * x / (1.0 - x)


def quasi_slack(theta, thalf, D, thc):
    need = thc - theta
    if need <= 0:
        return np.inf
    return np.log(trough(TAU, thalf, D) / need) / (np.log(2.0) / thalf)


def simulate(theta, thalf, D, p, delay=0.0, skip=False, illness=False, n_pre=8, n_post=4, k=1.0):
    """Usual cycle peak, peak in the long gap (to 3 d after the late dose), last-cycle peak."""
    lam = np.log(2.0) / thalf
    pp = replace(p, theta=float(theta))
    pre = [TAU * k for k in range(n_pre)]
    due = pre[-1] + TAU
    late = due + (TAU if skip else delay)
    doses = np.array(pre + [late + TAU * k for k in range(n_post)])
    ill = (due, due + ILL_DAYS) if illness else None

    def cap(t):
        past = doses[doses <= t + 1e-12]
        return theta + D * np.exp(-lam * (t - past)).sum()

    def load(t):
        return S_BASE * ILLNESS_FOLD if (ill and ill[0] <= t <= ill[1]) else S_BASE

    def rhs(t, y):
        m, c = y
        conj = k * c * m / (pp.Kc + m)
        delta = pp.rho0 * k * (1.0 + m / pp.Mi)
        return [load(t) - cap(t) * g(m, pp.Km) - conj, pp.R0 - conj - delta * c]

    edges = sorted(set(list(doses) + [doses[-1] + TAU] + ([ill[0], ill[1]] if ill else [])))
    y = [20.0, 500.0]
    T, M = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        sol = solve_ivp(rhs, (a, b), y, method="LSODA", rtol=1e-8, atol=1e-8, max_step=0.05)
        T.append(sol.t)
        M.append(sol.y[0])
        y = sol.y[:, -1]
    T, M = np.concatenate(T), np.concatenate(M)
    usual = M[(T >= pre[-2]) & (T < pre[-1])].max()
    gap = M[(T >= pre[-1]) & (T <= late + 3.0)].max()
    last = M[T >= doses[-2]].max()
    return usual, gap, last


def separatrix(theta, thalf, D, p):
    """Unstable equilibrium (tipping point) at the patient's usual trough capacity, load 600."""
    from common.models import clearance_C_quasistatic, equilibria_nonmonotone
    capq = theta + trough(TAU, thalf, D)
    pp = replace(p, theta=float(capq))
    roots, stab = equilibria_nonmonotone(clearance_C_quasistatic, pp, S_BASE, np.linspace(1e-3, 20000, 40000))
    un = [r for r, st in zip(roots, stab) if not st]
    return un[0] if un else np.inf


def dyn_slack(theta, thalf, D, p, crit="2x", k=1.0, hi=28.0, tol=0.2):
    """Days of delay before the gap peak reaches 2x the usual peak, or the separatrix."""
    if crit == "2x":
        level = 2.0 * simulate(theta, thalf, D, p, delay=0.0, k=k)[0]
    else:
        level = min(separatrix(theta, thalf, D, p), 5.5 * p.Mi)   # paper 03_dynamic: valid below ~5.5 M_i
    if simulate(theta, thalf, D, p, delay=hi, k=k)[1] < level:
        return np.inf
    lo_, hi_ = 0.0, hi
    while hi_ - lo_ > tol:
        mid = 0.5 * (lo_ + hi_)
        if simulate(theta, thalf, D, p, delay=mid, k=k)[1] < level:
            lo_ = mid
        else:
            hi_ = mid
    return lo_


def fmt(d):
    if not np.isfinite(d):
        return ">28 d"
    return "none" if d < 0.05 else f"{d:.1f} d"


def main():
    print("=" * 92)
    print("  ROUND 20. HOW LATE CAN AN INFUSION BE? (two-state parent model, fortnightly infusions)")
    print("=" * 92)
    pn, pw = shape(False), shape(True)
    print(f"  theta_crit_int: nominal {thc_of(pn):.1f}, worst corner {thc_of(pw):.1f} uM/d; illness"
          f" {ILLNESS_FOLD:.2f} x baseline for {ILL_DAYS:.0f} d")

    print("\n  1. DAYS OF SLACK. Four numbers per cell: trough condition | peak below 2x usual at k = 1 /")
    print("     k = 10 | peak below the separatrix at the usual trough capacity at k = 1 / k = 10.")
    print("     k (carnitine conjugation rate) is ASSUMED in the repo (params.py); k = 1 is the model's")
    print("     default, k = 10 is near the reduced (quasi-steady carnitine) limit the paper uses elsewhere.")
    settings = [("booked dose, nominal constants", D0, pn), ("booked dose, worst constants", D0, pw),
                ("admissible dose, nominal constants", D_ADM, pn), ("admissible dose, worst constants", D_ADM, pw)]
    for lab, D, p in settings:
        thc = thc_of(p)
        print(f"\n    {lab} (D = {D:.0f}):")
        print(f"      {'half-life':>9} {'capacity':>9} {'trough':>8} | {'2x, k=1':>8} {'2x, k=10':>9} |"
              f" {'sep, k=1':>9} {'sep, k=10':>10}")
        for th in (3.5, 4.3):
            for theta, c in CAPS:
                q = quasi_slack(theta, th, D, thc)
                qd = "fails" if q < 0 else fmt(q)
                vals = [dyn_slack(theta, th, D, p, crit=cr, k=kk) for cr in ("2x", "sep") for kk in (1.0, 10.0)]
                print(f"      {th:>8.1f}d {c:>9} {qd:>8} | {fmt(vals[0]):>8} {fmt(vals[1]):>9} |"
                      f" {fmt(vals[2]):>9} {fmt(vals[3]):>10}")

    print("\n  2. A SKIPPED INFUSION (28-day gap), and an infusion DELAYED DURING ILLNESS (illness starts on the")
    print("     day the infusion is due; 14 d late = postponed until recovery, the same gap as a skip).")
    print("     Half-life 3.5 / 4.3 d.")
    print("     Peak metabolite (uM) and its fold over the usual on-therapy peak; 'back' = the last")
    print("     cycle's peak once infusions resume, against the usual peak.")
    for (lab, D, p), kk in [(st, 1.0) for st in settings[:3]] + [(settings[0], 10.0)]:
        print(f"\n    {lab} (D = {D:.0f}), k = {kk:g}:")
        print(f"      {'capacity':>10} {'half-life':>9} {'usual':>7} {'skipped':>15} {'ill, on time':>15}"
              f" {'ill, 3 d late':>15} {'ill, 7 d late':>15} {'ill, 14 d late':>15} {'back':>11}")
        for theta, c in CAPS:
            for th in (3.5, 4.3):
                u = simulate(theta, th, D, p, k=kk)[0]
                sk = simulate(theta, th, D, p, skip=True, k=kk)
                ills = [simulate(theta, th, D, p, delay=d, illness=True, k=kk) for d in (0.0, 3.0, 7.0, ILL_DAYS)]
                cell = lambda r: f"{r[1]:7.0f} ({r[1]/u:3.1f}x)"
                print(f"      {c:>10} {th:>8.1f}d {u:7.1f} {cell(sk):>15} " + " ".join(f"{cell(r):>15}" for r in ills)
                      + f" {ills[-1][2]:5.0f}/{u:<5.0f}")

    # figure: slack against half-life, booked dose, nominal constants
    import matplotlib
    matplotlib.use("Agg")
    import plotstyle as P
    from plotstyle import plt
    plt.rcParams.update({"font.size": 7.5, "axes.labelsize": 8, "axes.titlesize": 8,
                         "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7})
    hs = np.linspace(3.0, 6.0, 13)
    cols = [P.C2, P.C4, P.C3, P.C1]
    fig, ax = plt.subplots(1, 2, figsize=(6.9, 2.8), sharey=True)
    thc = thc_of(pn)
    print("\n  3. FIGURE DATA (booked dose, nominal constants): slack in days by half-life")
    print(f"      {'half-life':>9} " + " ".join(f"{'q '+c:>9} {'2x '+c:>9} {'k10 '+c:>9}" for _, c in CAPS))
    rows = []
    for th in hs:
        q = [quasi_slack(theta, th, D0, thc) for theta, _ in CAPS]
        d2 = [dyn_slack(theta, th, D0, pn) for theta, _ in CAPS]
        d10 = [dyn_slack(theta, th, D0, pn, k=10.0) for theta, _ in CAPS]
        rows.append((q, d2, d10))
        print(f"      {th:>8.2f}d " + " ".join(f"{fmt(a):>9} {fmt(b):>9} {fmt(c):>9}" for a, b, c in zip(q, d2, d10)))
    for k, ((theta, c), col) in enumerate(zip(CAPS, cols)):
        qv = np.array([min(r[0][k], 28.0) for r in rows])
        dv = np.array([min(r[1][k], 28.0) for r in rows])
        dv10 = np.array([min(r[2][k], 28.0) for r in rows])
        ax[0].plot(hs, np.maximum(qv, 0.0), color=col, lw=1.6, label=f"{c} of normal")
        ax[1].plot(hs, dv, color=col, lw=1.6)
        ax[1].plot(hs, dv10, color=col, lw=1.1, ls="--")
    ax[1].text(5.95, 1.0, "solid: k = 1;  dashed: k = 10", fontsize=6.5, color=P.MUTED, ha="right")
    for a, ttl in zip(ax, ("(a) trough stays above the critical capacity",
                           "(b) metabolite peak stays below 2x usual")):
        a.axvline(3.6, color=P.MUTED, lw=0.8, ls=":")
        a.axvline(4.3, color=P.MUTED, lw=0.8, ls="--")
        a.set_xlim(3.0, 6.0); a.set_ylim(0, 28.5)
        a.set_yticks([0, 7, 14, 21, 28])
        a.set_xlabel("PCC protein half-life (d)")
        a.set_title(ttl, loc="left")
        P.finish(a)
    ax[0].set_ylabel("days an infusion can be late")
    ax[0].text(3.62, 26.5, "3.6 d\nfibroblast", fontsize=6.5, color=P.MUTED, va="top")
    ax[0].text(4.32, 26.5, "4.3 d\nlikely liver", fontsize=6.5, color=P.MUTED, va="top")
    hd, lb = ax[0].get_legend_handles_labels()
    fig.tight_layout(w_pad=1.2, rect=(0, 0.08, 1, 1))
    fig.legend(hd, [x.replace(" of normal", "") for x in lb], loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.02), title="residual capacity, % of normal", title_fontsize=7)
    out = HERE / "pi_note" / "fig"
    fig.savefig(out / "pi_tolerance.pdf", bbox_inches="tight")
    fig.savefig(out / "pi_tolerance.png", dpi=170, bbox_inches="tight")
    print("\n  wrote pi_note/fig/pi_tolerance.pdf (capped at 28 d = a whole missed infusion)")


if __name__ == "__main__":
    main()
