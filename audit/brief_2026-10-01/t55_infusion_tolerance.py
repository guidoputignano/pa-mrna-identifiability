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
(or skipped, a 28-day gap), then doses resume on time. The gap peak is the highest
metabolite level from the last on-time dose to 3 d after the late one, extended during an
illness to the next dose after the late one. Peaks above the model's validity limit (5.5 M_i: 1650 uM
nominal, 825 at the worst corner) are flagged with '!': they mean only that the model's
range is left. Load s = 600 (baseline), or the
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


def corner(Kc, Mi, rho0):
    """A corner of the paper's carnitine ranges (tab_C_shape), R0 re-derived as the repo does."""
    q = replace(MC.PC, Kc=Kc, Mi=Mi, rho0=rho0)
    return replace(q, R0=MC.R0_for_peak(q))


CORNERS = [(kc, mi, r0) for kc in (25.0, 800.0) for mi in (150.0, 1200.0) for r0 in (0.35, 1.0)]


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


def simulate(theta, thalf, D, p, delay=0.0, skip=False, illness=False, n_pre=8, n_post=4, k=1.0,
             traj=False):
    """Usual cycle peak, peak in the long gap (to 3 d after the late dose; during an illness, to
    the next dose after it), last-cycle peak."""
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
    # The gap peak runs from the last on-time dose to 3 d after the late dose and, during an
    # illness, to the next dose after the late one: with an on-time or slightly late dose the
    # highest level comes at the end of the illness or, at the smaller dose, after it, just before
    # that next dose (earlier versions stopped at late + 3 d, then at the illness's end).
    end = late + TAU if ill else late + 3.0
    gap = M[(T >= pre[-1]) & (T <= end)].max()
    last = M[T >= doses[-2]].max()
    if traj:
        return T - pre[-2], M, doses - pre[-2], (None if not ill else (ill[0] - pre[-2], ill[1] - pre[-2]))
    return usual, gap, last


def separatrix(theta, thalf, D, p):
    """Unstable equilibrium (tipping point) at the patient's usual trough capacity, load 600."""
    from common.models import clearance_C_quasistatic, equilibria_nonmonotone
    capq = theta + trough(TAU, thalf, D)
    pp = replace(p, theta=float(capq))
    roots, stab = equilibria_nonmonotone(clearance_C_quasistatic, pp, S_BASE, np.linspace(1e-3, 20000, 40000))
    un = [r for r, st in zip(roots, stab) if not st]
    return un[0] if un else np.inf


def dyn_slack(theta, thalf, D, p, crit="2x", k=1.0, hi=28.0, tol=0.05):
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
    return 0.5 * (lo_ + hi_)          # within tol/2 of the crossing


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
    print("     cycle's peak once infusions resume, against the usual peak. '!' = above the model's")
    print("     validity limit 5.5 M_i (1650 uM nominal, 825 uM worst corner): outside the model.")
    for (lab, D, p), kk in [(st, 1.0) for st in settings[:3]] + [(settings[0], 10.0)]:
        print(f"\n    {lab} (D = {D:.0f}), k = {kk:g}:")
        print(f"      {'capacity':>10} {'half-life':>9} {'usual':>7} {'skipped':>15} {'ill, on time':>15}"
              f" {'ill, 3 d late':>15} {'ill, 7 d late':>15} {'ill, 14 d late':>15} {'back':>11}")
        for theta, c in CAPS:
            for th in (3.5, 4.3):
                u = simulate(theta, th, D, p, k=kk)[0]
                sk = simulate(theta, th, D, p, skip=True, k=kk)
                ills = [simulate(theta, th, D, p, delay=d, illness=True, k=kk) for d in (0.0, 3.0, 7.0, ILL_DAYS)]
                vl = 5.5 * p.Mi
                cell = lambda r: f"{r[1]:7.0f} ({r[1]/u:3.1f}x){'!' if r[1] > vl else ' '}"
                print(f"      {c:>10} {th:>8.1f}d {u:7.1f} {cell(sk):>15} " + " ".join(f"{cell(r):>15}" for r in ills)
                      + f" {ills[-1][2]:5.0f}/{u:<5.0f}")

    print("\n  5. LEAST FAVOURABLE CONSTANTS, per column: the 8 corners of the paper's carnitine ranges")
    print("     (K_c 25/800, M_i 150/1200, rho_0 0.35/1; R0 re-derived), no own enzyme, k = 1. The corner")
    print("     with the highest theta_crit (25, 150, 1) is the worst for the trough condition but not")
    print("     always for the doubling, so each column takes its own minimum.")
    for D in (D0, D_ADM):
        for th in (4.3, 3.5):
            rows5 = []
            for kc, mi, r0 in CORNERS:
                p = corner(kc, mi, r0)
                rows5.append((kc, mi, r0, dyn_slack(0.0, th, D, p), quasi_slack(0.0, th, D, thc_of(p)),
                              simulate(0.0, th, D, p)[0]))
            print(f"\n    D = {D:.0f}, half-life {th} d:   K_c    M_i  rho0   2x (d)   trough (d)   usual peak (uM)")
            for kc, mi, r0, d2, q, u in rows5:
                print(f"    {'':28s}{kc:5.0f} {mi:6.0f} {r0:5.2f} {fmt(d2):>8} {('fails' if q < 0 else fmt(q)):>11}"
                      f" {u:12.1f}")
            d2m = min(r[3] for r in rows5)
            qm = min(r[4] for r in rows5)
            print(f"    {'':28s}minimum: 2x {fmt(d2m)}, trough {'fails on time' if qm < 0 else fmt(qm)}")

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
    fig.savefig(out / "pi_tolerance_slack.pdf", bbox_inches="tight")
    fig.savefig(out / "pi_tolerance_slack.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("\n  wrote pi_note/fig/pi_tolerance_slack.pdf (capped at 28 d late = two missed infusions; a skipped one is 14 d late)")
    fig_clear(pn)


def fig_clear(pn, thalf=4.3):
    """The figure for the note: (a) what a late infusion looks like; (b) peak against days late."""
    import plotstyle as P
    from plotstyle import plt
    fig, ax = plt.subplots(1, 2, figsize=(6.9, 2.9), gridspec_kw=dict(width_ratios=[1.15, 1]))
    a = ax[0]
    print(f"\n  4. FIGURE 5 DATA (booked dose, nominal constants, k = 1, half-life {thalf} d)")
    blues = plt.cm.Blues
    runs = [(None, "skipped", blues(0.95), "-"), (7.0, "7 days late", blues(0.72), "-"),
            (3.0, "3 days late", blues(0.5), "-"), (0.0, "on time", P.INK, "-")]
    u0 = simulate(0.0, thalf, D0, pn)[0]
    a.axhline(2.0 * u0, color=P.INK, lw=0.8, ls=":", zorder=0)
    a.text(69, 2.0 * u0 + 3, "twice the usual peak", ha="right", va="bottom", fontsize=6.5, color=P.INK)
    for d, lab, col, ls in runs:
        T, M, doses, _ = simulate(0.0, thalf, D0, pn, delay=d or 0.0, skip=d is None, traj=True)
        keep = (T >= 0) & (T <= 70)
        a.plot(T[keep], M[keep], color=col, lw=1.1 if d == 0.0 else 1.6, ls=ls, label=lab)
        cyc = [M[(T > doses[i] + 7.0) & (T < doses[i + 1])].max() / u0 for i in (8, 9)]   # pre-dose peaks
        print(f"    no own enzyme, {lab:>12}: peak in the gap {M[(T >= 14) & (T <= 45)].max():7.1f} uM;"
              f" after resuming, cycle peaks {cyc[0]:.3f}x then {cyc[1]:.3f}x usual")
    for t0 in (0, 14):
        a.axvline(t0, color=P.GRID, lw=0.8, zorder=0)
    a.axvline(28, color=P.MUTED, lw=0.8, ls=":", zorder=0)
    a.text(28.7, 222, "infusion due", fontsize=6.5, color=P.MUTED, va="top", ha="left")
    a.set_xlim(0, 70)
    a.set_ylim(0, 230)
    a.set_xlabel("days")
    a.set_ylabel("metabolite (\u00b5M)")
    a.set_title("(a) a patient with no enzyme of their own, not ill", loc="left")
    a.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="none", framealpha=1.0,
             title="day-28 infusion", title_fontsize=6.5, bbox_to_anchor=(0.0, 0.93))
    P.finish(a)

    b = ax[1]
    days = np.arange(0, 15, 1.0)
    vl = 5.5 * pn.Mi
    for theta, cap_lab, ls in ((0.0, "no own enzyme", "-"), (75.0, "2.5 % of normal", "--")):
        u = simulate(theta, thalf, D0, pn)[0]
        well = [simulate(theta, thalf, D0, pn, delay=d)[1] / u for d in days]
        ill = [simulate(theta, thalf, D0, pn, delay=d, illness=True)[1] / u for d in days]
        b.plot(days, well, color=P.C1, ls=ls, lw=1.5)
        ok = np.array(ill) * u <= vl          # stop the curve at the model's validity limit
        last = int(np.argmax(~ok)) if (~ok).any() else len(days)
        b.plot(days[:last], np.array(ill)[:last], color=P.C2, ls=ls, lw=1.5)
        if last < len(days):
            b.plot(days[last - 1], ill[last - 1], marker="o", ms=3, color=P.C2)
            print(f"    {cap_lab:>16}: ill curve stops after {days[last-1]:.0f} d late (next point"
                  f" {ill[last]*u:.0f} uM > validity limit {vl:.0f} uM)")
        print(f"    {cap_lab:>16}: fold over usual by days late 0..14, not ill: "
              + " ".join(f"{v:.1f}" for v in well))
        print(f"    {cap_lab:>16}: fold over usual by days late 0..14, ill:     "
              + " ".join(f"{v:.1f}" for v in ill))
    b.axhline(2.0, color=P.INK, lw=0.8, ls=":")
    b.text(14, 1.93, "peak doubled", ha="right", va="top", fontsize=6.5, color=P.INK)
    b.text(4.4, 12, "ill during\nthe delay", ha="right", va="center", fontsize=7, color=P.C2)
    b.text(9.6, 3.6, "not ill", ha="center", va="bottom", fontsize=7, color=P.C1)
    b.text(13.9, 5.9, "14 days late\n= skipped", ha="right", va="bottom", fontsize=6.3, color=P.MUTED)
    b.text(0.25, 52, "solid: no enzyme of their own\ndashed: 2.5 % of normal", ha="left", va="top",
           fontsize=6.3, color=P.MUTED)
    b.set_yscale("log")
    b.set_ylim(0.9, 60)
    b.set_yticks([1, 2, 5, 10, 20, 50])
    b.set_yticklabels(["1", "2", "5", "10", "20", "50"])
    b.text(10.5, 44, "curve stops at\nthe model's\nvalidity limit", fontsize=6.3, color=P.MUTED, va="top")
    b.set_xlim(0, 14)
    b.set_xticks([0, 2, 4, 6, 8, 10, 12, 14])
    b.set_xlabel("days the infusion is late")
    b.set_ylabel("peak / usual peak")
    b.set_title("(b) how much higher the peak goes", loc="left")

    P.finish(b)
    fig.tight_layout(w_pad=1.5)
    out = HERE / "pi_note" / "fig"
    fig.savefig(out / "pi_tolerance.pdf", bbox_inches="tight")
    fig.savefig(out / "pi_tolerance.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  wrote pi_note/fig/pi_tolerance.pdf")


if __name__ == "__main__":
    main()
