"""
Task 7. Experiment design for a transient at theta = 5% of normal.

  A  the washout: perturbation, sampling interval, duration, against what is
     known beforehand (a knowledge ladder), alone and joined to the monthly
     record. Is sd(ln theta) <= 0.20 reachable, and with what?
  B  the tension in tau, for the efficient estimator (levels, FIM) and for the
     finite-difference estimator (derivatives), as a fraction of 1/lambda
  C  the excursion: how sd(ln theta) and the conditioning of [1, m, m^2] move
     with how far m is swept, and whether crossing K_m is what matters
  D  the drug as the forcing input: an indirect rise to a peak at 13 h (the
     Attarwala two-compartment kernel's peak, audit/q1_increment_shape.py),
     3.5 d decay (design.run.LAMBDA), no substrate load, k_e untouched, with
     the delivered scale kappa known, known to 12% (the specification's Tier 1)
     or free; booked and admissible amplitudes (n4: the ceiling at 5% of
     normal allows a peak increment of 2850, the booked train peaks at 6400)

Assumptions this needs and the repository does not supply are dyn.py's (the
split of the load into known u and unrecorded L0 = 0.5 L, constant through
the experiment); L0's share is swept in part A.
"""
import numpy as np
from scipy.optimize import brentq

import core as C
import dyn as D
from plotstyle import plt, C1, C2, C3, C4, INK, MUTED, finish, save, plain

L = D.LTOT
H = 1.0 / 24
KD = np.log(2) / 3.5
KR = brentq(lambda k: np.log(k / KD) / (k - KD) - 13 / 24, KD * 1.01, 200)
TP = np.log(KR / KD) / (KR - KD)
NORM = np.exp(-KD * TP) - np.exp(-KR * TP)
FR = C.fim(C.record(C.TH_5, seed=0), ke_prior=None)     # the monthly record, data + b prior
KM_REC = 0.062                                         # sd(ln K_m) the record delivers (Task 1)

#: the knowledge ladder. Each rung adds one thing to the one above.
LADDER = (("nothing", {}),
          ("k_e prior 0.8", {"ke": 0.8}),
          ("+ K_m from the record (0.062)", {"ke": 0.8, "Km": KM_REC}),
          ("+ L0 to 5% of L", {"ke": 0.8, "Km": KM_REC, "L0": 0.05 * L}),
          ("+ L0 to 1% of L", {"ke": 0.8, "Km": KM_REC, "L0": 0.01 * L}),
          ("k_e exact, K_m record, L0 1%", {"ke": 1e-4, "Km": KM_REC, "L0": 0.01 * L}))
IDX = {"ke": 1, "Km": 2, "L0": 3, "kap": 4}


def kernel(Dpk):
    return lambda t: 0.0 if t < 0 else Dpk * (np.exp(-KD * t) - np.exp(-KR * t)) / NORM


def with_priors(F, pri):
    F = F.copy()
    for k, sd in pri.items():
        F[IDX[k], IDX[k]] += sd ** -2
    return F


def joint_record(Fd, pri):
    """The transient joined to the monthly record. z = (th, ke, Km, L0, b[, kap])."""
    n = Fd.shape[0] + 1
    Z = np.zeros((n, n))
    idd = [0, 1, 2, 3] + ([5] if Fd.shape[0] > 4 else [])
    Z[np.ix_(idd, idd)] += Fd
    Z[np.ix_([0, 2, 4, 1], [0, 2, 4, 1])] += FR
    pri = {k: v for k, v in pri.items() if k != "Km"}     # the record supplies K_m itself
    for k, sd in pri.items():
        j = 5 if k == "kap" else IDX[k]
        Z[j, j] += sd ** -2
    return D.sd_lnth(Z)


def washout(r_out=0.5, up=None, T=48, tau=1.0, L0frac=0.5):
    if up is None:
        steps = [(-1, (1 - L0frac) * L), (0, (1 - L0frac - r_out) * L)]
    else:
        steps = [(-1, (1 - L0frac) * L), (0, (1 - L0frac) * L + (up - 1) * L), (T / 48, 0.0)]
    ex = D.Experiment(steps, L0=L0frac * L)
    t = np.arange(0, T / 24 + 1e-9, tau / 24)
    p = D.params(ex=ex)
    return D.fim(p, ex, t), ex, p, t


def drug(Dpk, T=72, tau=0.25, kfree=False):
    ex = D.Experiment([(-1, 0.5 * L)], L0=0.5 * L, delta=kernel(Dpk), kappa_free=kfree)
    t = np.arange(0, T / 24 + 1e-9, tau / 24)
    p = D.params(ex=ex)
    return D.fim(p, ex, t), ex, p, t


def main():
    print("=" * 92)
    print("  TASK 7. EXPERIMENT DESIGN AT 5% OF NORMAL")
    print("=" * 92)

    # ------------------------------------------------------------------ A
    print("\n  A. The washout, hourly, 48 h, 5% CV. sd(ln theta) down the knowledge ladder")
    designs = (("withdraw 25% of the load", dict(r_out=0.25)),
               ("withdraw 50% (all the diet)", dict(r_out=0.5)),
               ("load x2 for 24 h, then diet out", dict(up=2)),
               ("load x3 for 24 h, then diet out", dict(up=3)),
               ("load x6 for 24 h, then diet out", dict(up=6)))
    print(f"     {'':34}" + "".join(f"{lab[:13]:>14}" for lab, _ in LADDER) + f"{'+ record':>10}")
    Aw = {}
    for lab, kw in designs:
        F, ex, p, t = washout(**kw)
        m, _ = D.simulate(p, ex, t, sens=False)
        row = [D.sd_lnth(with_priors(F, pri)) for _, pri in LADDER]
        jr = joint_record(F, LADDER[4][1])
        Aw[lab] = row + [jr]
        print(f"     {lab:34}" + "".join(f"{v:14.3f}" for v in row) + f"{jr:10.3f}"
              f"   m {m.min():.0f}-{m.max():.0f}")
    print("     '+ record' = the washout joined to the 25-draw monthly record, k_e prior 0.8,")
    print("     L0 to 1% of L; the record supplies K_m and its own load-bias prior.")

    print("\n     The half-load washout over (tau, T), rung '+ L0 to 1% of L':")
    taus = (0.083, 0.25, 0.5, 1, 2, 4)
    print(f"     {'T':>6}" + "".join(f"{f'tau {x:g} h':>12}" for x in taus))
    for T in (12, 24, 48, 96):
        row = [D.sd_lnth(with_priors(washout(T=T, tau=x)[0], LADDER[4][1])) for x in taus]
        print(f"     {T:4d} h" + "".join(f"{v:12.3f}" for v in row))
    print("     The x6 excursion over (tau, T), same rung:")
    for T in (24, 48, 96):
        row = [D.sd_lnth(with_priors(washout(up=6, T=T, tau=x)[0], LADDER[4][1])) for x in taus]
        print(f"     {T:4d} h" + "".join(f"{v:12.3f}" for v in row))
    print("     Sensitivity to the assumed unrecorded share L0/L, x6 design, rung '+ L0 1%':")
    for fr in (0.3, 0.5, 0.7):
        F, *_ = washout(up=6, L0frac=fr)
        print(f"       L0 = {fr:.1f} L: {D.sd_lnth(with_priors(F, LADDER[4][1])):.3f}")

    # ------------------------------------------------------------------ B
    print("\n  B. The tension in tau")
    F, ex, p, t = washout()
    lam0 = D.lam(p, ex, D.steady(p, ex))
    m_inf = float(C.m_ss(C.TH_5, C.KM0, C.KE, 0.5 * L))
    lam1 = D.lam(p, ex, m_inf)
    print(f"     washout: lambda = {lam0:.3f} /d at m0 = 54.0 and {lam1:.3f} at the new level "
          f"{m_inf:.1f}; 1/lambda = {24/lam0:.2f} to {24/lam1:.2f} h")
    print("     (i) efficient estimator, a FIXED number of samples N = 24, rung '+ L0 1%':")
    rows = []
    for tl in (0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0):
        tau_h = tl * 24 / lam0
        T = 23 * tau_h
        ex2 = D.Experiment([(-1, 0.5 * L), (0, 0.0)], L0=0.5 * L)
        tt = np.arange(0, 24) * tau_h / 24
        v = D.sd_lnth(with_priors(D.fim(D.params(ex=ex2), ex2, tt), LADDER[4][1]))
        rows.append((tl, v))
        print(f"        tau = {tl:4.2f}/lambda ({tau_h:5.2f} h, T = {T:5.1f} h): {v:.3f}")
    tl_opt = min(rows, key=lambda r: r[1])[0]
    print(f"        optimum near tau = {tl_opt}/lambda: N samples should span ~{23*tl_opt:.0f}/lambda")
    print("     (ii) fixed duration 48 h: denser is never worse for the efficient estimator")
    print("          (part A's rows fall monotonically as tau shrinks).")
    print("     (iii) the finite-difference derivative psi: MSE over the transient of the")
    print("           forward and central differences, noise 5% CV, exponential relaxation")
    tt_f = np.linspace(0, 4 / lam0, 400)
    A0, m1 = 54.02 - m_inf, m_inf
    res = []
    for tl in np.geomspace(0.01, 3, 80):
        tau = tl / lam0
        m_t = m1 + A0 * np.exp(-lam0 * tt_f)
        sig = 0.05 * m_t
        d_true = -lam0 * A0 * np.exp(-lam0 * tt_f)
        fwd_bias = (A0 * np.exp(-lam0 * (tt_f + tau)) - A0 * np.exp(-lam0 * tt_f)) / tau - d_true
        cen_bias = (A0 * np.exp(-lam0 * (tt_f + tau / 2)) - A0 * np.exp(-lam0 * (tt_f - tau / 2))) / tau - d_true
        noise = 2 * sig ** 2 / tau ** 2
        res.append((tl, np.sqrt(np.mean(fwd_bias ** 2 + noise)) / np.sqrt(np.mean(d_true ** 2)),
                    np.sqrt(np.mean(cen_bias ** 2 + noise)) / np.sqrt(np.mean(d_true ** 2))))
    res = np.array(res)
    jf, jc = np.argmin(res[:, 1]), np.argmin(res[:, 2])
    print(f"           forward: optimum tau = {res[jf,0]:.3f}/lambda, relative RMS error {res[jf,1]:.2f}")
    print(f"           central: optimum tau = {res[jc,0]:.3f}/lambda, relative RMS error {res[jc,2]:.2f}")
    print("           a relative error near 1 means the derivative is not measured at all,")
    print("           which is why every FD row of Task 6's bench A fails.")

    # ------------------------------------------------------------------ C
    print("\n  C. The excursion: sd(ln theta), and cond of the [1, m, m^2] columns over the path")
    print(f"     {'load factor':>12} {'m range':>14} {'m_max/K_m':>10} {'cond[1,m,m^2]':>14} "
          f"{'cond FIM':>10} {'rung k_e+K_m':>13} {'rung +L0 1%':>12}")
    Cx = []
    for up in (1.25, 1.5, 2, 3, 4, 6, 8, 12):
        F, ex, p, t = washout(up=up)
        m, _ = D.simulate(p, ex, t, sens=False)
        X = np.stack([np.ones_like(m), m, m ** 2], axis=1)
        cX = np.linalg.cond(X / np.linalg.norm(X, axis=0))
        a = D.sd_lnth(with_priors(F, LADDER[2][1]))
        b = D.sd_lnth(with_priors(F, LADDER[4][1]))
        Cx.append((up, m.max() / m.min(), m.max() / C.KM0, cX, a, b))
        print(f"     {up:12.2f} {m.min():6.1f}-{m.max():6.1f} {m.max()/C.KM0:10.2f} {cX:14.1f} "
              f"{np.linalg.cond(F):10.2e} {a:13.3f} {b:12.3f}")
    Cx = np.array(Cx)

    # ------------------------------------------------------------------ D
    print("\n  D. The drug as the forcing input: first dose, every 15 min for 72 h")
    print(f"     kernel: rise rate {KR:.2f} /d, peak at {24*TP:.1f} h, decay half-life 3.5 d")
    print(f"     {'peak Delta':>10} {'Delta/th':>9} {'m min':>6}" + "".join(
        f"{lab[:13]:>14}" for lab, _ in LADDER) + f"{'+record':>9} {'kap 12%+rec':>12} {'kap free+rec':>13}")
    Dd = {}
    for Dpk in (150, 450, 1500, 2850, 6400):
        F, ex, p, t = drug(Dpk)
        m, _ = D.simulate(p, ex, t, sens=False)
        row = [D.sd_lnth(with_priors(F, pri)) for _, pri in LADDER]
        jr = joint_record(F, LADDER[4][1])
        Fk, *_ = drug(Dpk, kfree=True)
        jk = joint_record(Fk, {**LADDER[4][1], "kap": 0.12})
        jfree = joint_record(Fk, LADDER[4][1])
        Dd[Dpk] = row + [jr, jk, jfree]
        tag = "" if Dpk <= 2850 else "  (booked; NOT admissible at 5%)"
        print(f"     {Dpk:10d} {Dpk/C.TH_5:9.1f} {m.min():6.1f}" + "".join(f"{v:14.3f}" for v in row)
              + f"{jr:9.3f} {jk:12.3f} {jfree:13.3f}{tag}")
    print("     The burst with k_e EXACT and the record, peak 2850:")
    F, *_ = drug(2850)
    print(f"       L0 unknown {joint_record(F, {'ke': 1e-4}):.3f};  L0 to 1% of L "
          f"{joint_record(F, {'ke': 1e-4, 'L0': 0.01 * L}):.3f}")
    print("     How well k_e must be known for burst + record to reach 0.20, L0 unknown:")
    Fk, *_ = drug(2850, kfree=True)
    for lab, FF, extra in (("kappa known", F, {}), ("kappa to 12%", Fk, {"kap": 0.12})):
        f = lambda v: joint_record(FF, {"ke": v, **extra})
        lo, hi = 1e-4, 0.8
        if f(lo) > C.GATE:
            print(f"       {lab}: not reachable")
            continue
        for _ in range(50):
            mid = np.sqrt(lo * hi)
            lo, hi = (mid, hi) if f(mid) <= C.GATE else (lo, mid)
        print(f"       {lab}: sd(ln k_e) <= {lo:.4f} (k_e to {100*lo:.1f}%)")
    #  and the record alone at the same k_e, for comparison
    Fr = C.fim(C.record(C.TH_5, seed=0), ke_prior=0.05)
    print(f"     the monthly record alone with k_e to 5%: {C.sd_th(Fr):.3f}; burst + record "
          f"with k_e to 5%: {joint_record(F, {'ke': 0.05}):.3f}")

    # ---------------------------------------------------------------- figure
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 3.9))
    a = ax[0]
    labs = [lab for lab, _ in LADDER]
    for (lab, _), col in zip(designs, (C1, C2, C3, C4, INK)):
        a.plot(range(len(LADDER)), Aw[lab][:len(LADDER)], "-o", color=col, ms=4, lw=1.8, label=lab)
    a.axhline(C.GATE, color=INK, ls="--", lw=1)
    a.set_yscale("log")
    a.set_xticks(range(len(LADDER)), ["none", "k_e 0.8", "+K_m", "+L0 5%", "+L0 1%", "k_e exact"],
                 rotation=30, fontsize=8)
    a.set_ylabel(r"sd(ln $\hat\theta$), 5% of normal")
    a.text(0, 0.215, "target 0.20", fontsize=8, color=INK)
    a.set_title("(a) washout, hourly 48 h: what has to be known")
    a.legend(fontsize=7, frameon=False, loc="upper right")
    finish(a); plain(a)
    b = ax[1]
    b.plot(res[:, 0], res[:, 1], color=C1, lw=2)
    b.plot(res[:, 0], res[:, 2], color=C2, lw=2)
    b.plot([res[jf, 0]], [res[jf, 1]], "o", color=C1, ms=5)
    b.plot([res[jc, 0]], [res[jc, 2]], "o", color=C2, ms=5)
    b.annotate(f"forward difference\nmin {res[jf,1]:.2f} at τ = {res[jf,0]:.2f}/λ",
               (res[jf, 0], res[jf, 1]), (0.02, 0.35), fontsize=8,
               arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
    b.annotate(f"central difference\nmin {res[jc,2]:.2f} at τ = {res[jc,0]:.2f}/λ",
               (res[jc, 0], res[jc, 2]), (0.02, 0.12), fontsize=8,
               arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
    b.set_ylim(0.08, 30)
    b.set_xscale("log"); b.set_yscale("log")
    b.set_xlabel("sampling interval τ × λ")
    b.set_ylabel("relative RMS error of dm/dt")
    b.set_title("(b) the tension exists for derivatives only")
    finish(b); plain(b)
    c = ax[2]
    c.plot(Cx[:, 2], Cx[:, 4], "-o", color=C1, ms=4, lw=2, label="k_e 0.8 + K_m from record")
    c.plot(Cx[:, 2], Cx[:, 5], "-s", color=C2, ms=4, lw=2, label="… + L0 to 1% of L")
    c.axhline(C.GATE, color=INK, ls="--", lw=1)
    c.axvline(1.0, color=MUTED, ls=":", lw=1)
    c.text(1.03, 3, "m = K_m", fontsize=8, color=MUTED)
    c.text(5.5, 0.21, "target 0.20", fontsize=8, color=INK)
    c.set_xscale("log"); c.set_yscale("log")
    c.set_xlabel(r"peak marker / $K_m$ reached by the load step")
    c.set_ylabel(r"sd(ln $\hat\theta$)")
    c.set_title("(c) excursion: crossing K_m")
    c.legend(fontsize=7.5, frameon=False, loc="lower left")
    finish(c); plain(c)
    save(fig, "t7_design")


if __name__ == "__main__":
    main()
