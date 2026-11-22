"""
Task 2. The joint fit across records, and what borrowing k_e buys.

WHAT THE 24 RECORDS ARE. design.run.SEEDS = range(24) are 24 draws of the
day-to-day intake sequence for ONE simulated patient (theta, K_m = 100,
k_e = 5 fixed). n4 scores each with its own information matrix at the truth
and reports the median and range. Nothing is fitted and nothing is shared, so
yes: "one at a time". But they are not patients, so "refit jointly with k_e
shared across patients" cannot be done on them. Pooling them is reported once,
in part A, for what it is: one patient eating 24 different two-year diets.

THE COHORT QUESTION, posed properly (parts B-D). N patients, each with own
ln theta_i and own load-record bias b_i (a diary bias is a property of the
patient's diary). K_m (PCC's Michaelis constant) and k_e are population
constants, with between-patient variability in k_e allowed:

    ln k_e,i = ln k_e,pop + eta_i ,   eta_i ~ N(0, omega^2)
    ln k_e,pop ~ N(ln 5, sigma_pop^2),  sigma_pop^2 + omega^2 = 0.8^2

so that any one patient's marginal prior on ln k_e is the repository's 0.8,
and N = 1 reproduces n4 exactly. omega = 0 is "k_e genuinely shared"; omega =
0.8 is "the 0.8 is all between-patient variability", where the population
mean is known and pooling cannot help. Which of the two the repository's 0.8
is, the repository does not say (it calls it a population prior).

Part E checks the information matrix by actually fitting: a joint MAP over
N patients at 5% of normal with k_e and K_m shared, replicated.
"""
import numpy as np
from scipy.optimize import least_squares

import core as C
from plotstyle import plt, C1, C2, C3, C4, INK, MUTED, finish, save, plain

POP = C.KE_PRIOR
NS = (1, 2, 4, 8, 16, 32, 64, 128, 256, 1024)
OMEGAS = (0.0, 0.1, 0.2, 0.33, 0.5, 0.8)


def data_fim(theta, seed, times=C.MONTHLY):
    """One patient's data information in (ln th, ln Km, b, ln ke), b prior in."""
    return C.fim(C.record(theta, seed=seed, times=times), ke_prior=None)


def cohort(thetas, seeds, omega, times=C.MONTHLY, target=0):
    """
    Joint posterior for a cohort. Shared g = (ln Km, ln ke_pop); local
    l_i = (ln th_i, b_i, eta_i) (eta_i dropped when omega = 0).
    Returns (sd ln theta of patient `target`, sd ln ke_pop, sd ln Km, sd of
    patient target's own ln ke).
    """
    sp2 = POP ** 2 - omega ** 2
    G = np.zeros((2, 2))
    G[1, 1] = 1.0 / sp2 if sp2 > 0 else 1e12   # omega = 0.8: pop mean known
    locs = []
    for th, s in zip(thetas, seeds):
        F = data_fim(th, s, times)
        # map (ln th, ln Km, b, ln ke) -> (ln th, b, eta | ln Km, ln ke_pop)
        if omega > 0:
            T = np.zeros((4, 5))
            T[0, 0] = 1; T[2, 1] = 1; T[3, 2] = 1          # local: th, b, eta
            T[1, 3] = 1; T[3, 4] = 1                       # shared: Km, ke_pop
            Fx = T.T @ F @ T
            Fx[2, 2] += 1.0 / omega ** 2
            nl = 3
        else:
            T = np.zeros((4, 4))
            T[0, 0] = 1; T[2, 1] = 1; T[1, 2] = 1; T[3, 3] = 1
            Fx = T.T @ F @ T
            nl = 2
        A, B, D = Fx[:nl, :nl], Fx[:nl, nl:], Fx[nl:, nl:]
        Ai = np.linalg.inv(A)
        G += D - B.T @ Ai @ B
        locs.append((Ai, B, nl))
    Sg = np.linalg.inv(G)
    Ai, B, nl = locs[target]
    # Var(local) = Ai + Ai B Sg B' Ai
    W = Ai @ B
    Vl = Ai + W @ Sg @ W.T
    sd_th = np.sqrt(Vl[0, 0])
    if nl == 3:
        # ln ke_i = ln ke_pop + eta_i
        cov_eta_pop = -(W @ Sg)[2, 1]
        v_ke = Sg[1, 1] + Vl[2, 2] + 2 * cov_eta_pop
    else:
        v_ke = Sg[1, 1]
    return float(sd_th), float(np.sqrt(Sg[1, 1])), float(np.sqrt(Sg[0, 0])), float(np.sqrt(v_ke))


def disease_cohort(n, frac=0.05):
    return [C.TH0 * frac] * n, list(range(n))


def mixed_cohort(n):
    """q14's composition: 5, 8, 10, 15, 20% of normal, cycled; target = a 5%."""
    caps = (0.05, 0.08, 0.10, 0.15, 0.20)
    return [C.TH0 * caps[i % 5] for i in range(n)], list(range(n))


# ------------------------------------------------------------ Monte Carlo
def joint_map(thetas, recs, ys, sgs, rng, share_km=True):
    """
    Joint MAP: shared (ln Km, ln ke), per-patient (ln th_i, b_i), with the
    priors read as noisy pseudo-observations (centres drawn about the truth),
    which is the frequentist reading under which the Bayesian FIM is the
    sampling covariance. Returns (ln th_hat per patient, ln ke_hat).
    """
    n = len(recs)
    ke_c = np.log(C.KE) + rng.normal(0, POP)
    b_c = np.array([rng.normal(0, C.LOAD_FRAC * r["L"]) for r in recs])
    sb = np.array([C.LOAD_FRAC * r["L"] for r in recs])
    nd = recs[0]["m"].size

    def unpack(x):
        return x[0], x[1], x[2:2 + n], x[2 + n:]

    def res(x):
        lk, lke, lth, b = unpack(x)
        K, ke = np.exp(lk), np.exp(lke)
        out = []
        for i, r in enumerate(recs):
            m, inc = r["m"], r["inc"]
            pred = (np.exp(lth[i]) + inc) * m / (K + m) + ke * m + b[i]
            out.append((ys[i] - pred) / sgs[i])
        out.append([(lke - ke_c) / POP])
        out.append((b - b_c) / sb)
        return np.concatenate(out)

    def jac(x):
        lk, lke, lth, b = unpack(x)
        K, ke = np.exp(lk), np.exp(lke)
        J = np.zeros((n * nd + 1 + n, 2 + 2 * n))
        for i, r in enumerate(recs):
            m, inc = r["m"], r["inc"]
            th = np.exp(lth[i])
            rows = slice(i * nd, (i + 1) * nd)
            J[rows, 0] = (th + inc) * K * m / (K + m) ** 2 / sgs[i]
            J[rows, 1] = -ke * m / sgs[i]
            J[rows, 2 + i] = -th * m / (K + m) / sgs[i]
            J[rows, 2 + n + i] = -1.0 / sgs[i]
        J[n * nd, 1] = 1.0 / POP
        J[n * nd + 1:, 2 + n:] = np.diag(1.0 / sb)
        return J

    lo = np.r_[np.log(5.0), np.log(1e-3), np.full(n, np.log(1e-3)), np.full(n, -np.inf)]
    hi = np.r_[np.log(5000.0), np.log(500.0), np.full(n, np.log(1e5)), np.full(n, np.inf)]
    best = None
    for th0, ke0 in ((np.log(thetas), ke_c), (np.full(n, np.log(600.0)), ke_c - 0.5),
                     (np.full(n, np.log(30.0)), ke_c + 0.5)):
        x0 = np.r_[np.log(C.KM0), np.clip(ke0, lo[1] + 1e-6, hi[1] - 1e-6),
                   np.clip(th0, lo[2] + 1e-6, hi[2] - 1e-6), np.zeros(n)]
        s = least_squares(res, x0, jac=jac, bounds=(lo, hi), method="trf",
                          x_scale="jac", xtol=1e-10, ftol=1e-10, max_nfev=3000)
        if best is None or s.cost < best.cost:
            best = s
    _, lke, lth, _ = unpack(best.x)
    return lth, lke


def mc(n, reps=120, frac=0.05, seed0=0):
    rng = np.random.default_rng(seed0)
    th = C.TH0 * frac
    e_th, e_ke = [], []
    for k in range(reps):
        recs, ys, sgs = [], [], []
        for i in range(n):
            r = C.record(th, seed=1000 * k + i)
            y, sg = C.noisy(r, 1000 * k + i)
            recs.append(r); ys.append(y); sgs.append(sg)
        lth, lke = joint_map([th] * n, recs, ys, sgs, rng)
        e_th.append(lth[0] - np.log(th))
        e_ke.append(lke - np.log(C.KE))
    e_th, e_ke = np.array(e_th), np.array(e_ke)
    robust = lambda e: 0.5 * (np.percentile(e, 84) - np.percentile(e, 16))
    return dict(rmse=float(np.sqrt(np.mean(e_th ** 2))), robust=float(robust(e_th)),
                bias=float(np.median(e_th)), ke_rmse=float(np.sqrt(np.mean(e_ke ** 2))),
                off4=float(np.mean(np.abs(e_th) > np.log(4.5))))


def main():
    print("=" * 92)
    print("  TASK 2. JOINT FIT: WHAT BORROWING k_e ACROSS A COHORT BUYS")
    print("=" * 92)

    # ---------------------------------------------------------------- A
    print("\n  A. The 24 records, as the repository treats them, and pooled")
    for th in (C.TH0, C.TH_5):
        one = C.sd_theta(th)
        Fs = [C.fim(C.record(th, seed=s), ke_prior=None) for s in C.SEEDS]
        Fp = sum(Fs)
        Fp[C.I_KE, C.I_KE] += POP ** -2
        print(f"     theta = {th:5.0f}: one at a time, median {one[0]:.4f} [{one[1]:.4f}, "
              f"{one[2]:.4f}];  all 24 pooled as ONE patient (600 draws, own b per record "
              f"but one b would be similar): {C.sd_th(Fp):.4f}")

    # ---------------------------------------------------------------- B
    print("\n  B. A disease cohort: N patients all at 5% of normal, monthly records.")
    print("     Target patient's sd(ln theta), and sd of the shared ln k_e,pop")
    print(f"     {'N':>5} " + " ".join(f"{'w=' + str(w):>14}" for w in OMEGAS))
    B = {}
    for n in NS:
        th, sd = disease_cohort(n)
        cells = []
        for w in OMEGAS:
            a = cohort(th, sd, w)
            B[(n, w)] = a
            cells.append(f"{a[0]:6.3f} ({a[1]:5.3f})")
        print(f"     {n:5d} " + " ".join(f"{c:>14}" for c in cells))
    print("     cell: sd(ln theta_1) (sd ln k_e,pop). w = omega, between-patient sd of ln k_e.")
    print(f"     N = 1 reproduces n4's 3.0673 for record 0: {B[(1, 0.0)][0]:.4f} vs "
          f"{C.sd_th(C.fim(C.record(C.TH_5, seed=0))):.4f}")
    print("     The patient's OWN ln k_e sd, which is the quantity that limits theta_i:")
    print(f"     {'N':>5} " + " ".join(f"{'w=' + str(w):>8}" for w in OMEGAS))
    for n in (1, 16, 64, 1024):
        print(f"     {n:5d} " + " ".join(f"{B[(n, w)][3]:8.3f}" for w in OMEGAS))

    # ---------------------------------------------------------------- C
    print("\n  C. q14's composition (5/8/10/15/20% cycled), target the first 5% patient")
    print(f"     {'N':>5} {'sd ln th, w=0':>14} {'w=0.33':>8} {'sd ke_pop w=0':>14}")
    Cm = {}
    for n in (1, 4, 16, 64, 128, 1024):
        th, sd = mixed_cohort(n)
        a, b = cohort(th, sd, 0.0), cohort(th, sd, 0.33)
        Cm[n] = (a, b)
        print(f"     {n:5d} {a[0]:14.3f} {b[0]:8.3f} {a[1]:14.4f}")
    print("     q14 printed (PHASE10_SENSORS.md:535-541) 1.578 at n=16 and 0.969 at n=64,")
    print("     passing only k_e's sd into a record whose K_m stays free; this passes both.")

    # ---------------------------------------------------------------- D
    print("\n  D. The ceiling on borrowing: patient sd(ln theta) with k_e EXACT, single record")
    for frac in (0.05, 0.10, 0.20):
        th = C.TH0 * frac
        F = lambda s: C.fim(C.record(th, seed=s))
        keep = [0, 1, 2]
        v = np.median([np.sqrt(np.linalg.inv(F(s)[np.ix_(keep, keep)])[0, 0]) for s in C.SEEDS])
        print(f"     {100*frac:3.0f}% of normal: {v:.4f}  (and with K_m exact too: "
              f"{np.median([1/np.sqrt(F(s)[0,0] - F(s)[0,2]**2/F(s)[2,2]) for s in C.SEEDS]):.4f})")

    # ---------------------------------------------------------------- E
    print("\n  E. Monte Carlo: joint MAP over N patients at 5%, k_e and K_m shared (omega = 0)")
    print(f"     {'N':>4} {'FIM sd':>8} {'MC RMSE':>8} {'MC robust':>10} {'median err':>11} "
          f"{'>4.5x off':>10} {'ke RMSE':>8} {'FIM ke':>8}")
    E = {}
    for n, reps in ((1, 200), (4, 150), (16, 100), (64, 60)):
        r = mc(n, reps=reps)
        E[n] = r
        f = B[(n, 0.0)]
        print(f"     {n:4d} {f[0]:8.3f} {r['rmse']:8.3f} {r['robust']:10.3f} {r['bias']:+11.3f} "
              f"{r['off4']:10.1%} {r['ke_rmse']:8.3f} {f[1]:8.3f}")

    # ---------------------------------------------------------------- figure
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.0))
    a = ax[0]
    cols = {0.0: C1, 0.2: C3, 0.33: C4, 0.5: C2}
    for w, col in cols.items():
        y = [B[(n, w)][0] for n in NS]
        a.plot(NS, y, "-o", color=col, ms=3.5, lw=2)
        a.text(NS[-1] * 1.15, y[-1], f"ω = {w}", color=INK, fontsize=8, va="center")
    a.plot(list(E), [E[n]["robust"] for n in E], "D",
           color=INK, ms=5, mfc="white", label="joint MAP, robust spread (ω = 0)")
    a.axhline(0.20, color=INK, ls="--", lw=1)
    a.axhline(0.388, color=MUTED, ls=":", lw=1)
    a.text(1.05, 0.40, "k_e exact, one record: 0.39", fontsize=8, color=MUTED, va="bottom")
    a.text(1.05, 0.205, "target 0.20", fontsize=8, color=INK, va="bottom")
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(0.8, 4000)
    a.set_xlabel("patients in the cohort, all at 5% of normal")
    a.set_ylabel(r"one patient's sd(ln $\hat\theta$)")
    a.set_title("(a) borrowing k_e: per-patient precision")
    a.legend(fontsize=7.5, frameon=False, loc="upper right")
    finish(a)
    b = ax[1]
    for w, col in cols.items():
        y = [B[(n, w)][3] for n in NS]
        b.plot(NS, y, "-o", color=col, ms=3.5, lw=2)
        b.text(NS[-1] * 1.15, y[-1], f"ω = {w}", color=INK, fontsize=8, va="center")
    b.set_xscale("log"); b.set_yscale("log"); b.set_xlim(0.8, 4000)
    b.set_xlabel("patients in the cohort, all at 5% of normal")
    b.set_ylabel("sd of the patient's own ln k_e")
    b.set_title("(b) what the cohort teaches about k_e")
    finish(b)
    for x in ax:
        plain(x)
    save(fig, "t2_cohort")


if __name__ == "__main__":
    main()
