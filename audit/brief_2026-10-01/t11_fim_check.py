"""
Follow-up item 6. Does the information matrix mislead, and in which direction?

THE CLAIM UNDER TEST. t2_cohort's Monte Carlo (part E) put the FIM below the
MC spread at N = 1 (FIM 3.08 against an MC robust spread 6.6) and above it at
N = 64 (FIM 0.81 against MC RMSE 0.49, robust 0.46). REPORT.md says the cause
is not isolated. This script tests candidate causes; it does not assume one.

THE ESTIMATOR. t2_cohort.joint_map is a dense trust-region fit over all
2 + 2N parameters (its time per fit is printed in Part 1.0). Here the SAME
objective (same residuals, same priors read as pseudo-observations, same
bounds) is minimised by variable projection: for fixed (ln K_m, ln k_e) the
model is linear in every patient's (theta_i, b_i), and the bound on theta_i is
a single box constraint on a 2x2 convex quadratic, so each patient's
(theta_i, b_i) is solved exactly. What is left is a 2-D search over
(ln K_m, ln k_e), done on a grid over the parameter box followed by local
refinement from the best grid cells. The 2-D grid is also the global-optimum
check. Part 1.0 shows it reproduces t2's printed numbers digit for digit and
agrees with joint_map itself.

PART 1 (item 6). (a) dependence on the arbitrary lower bound theta >= 1e-3;
(b) MC at N in {1, 2, 4, 8, 16, 32, 64}: the constrained MAP, the same draws
refitted with theta_i unconstrained (real-valued, may go negative), with the
prior centres held at the truth, and a linear-Gaussian surrogate; global
optimum checks; profile-likelihood intervals for theta_1 at N = 64; (c) the
surrogate with 2000 replicates, where MC and FIM must agree.

PART 2. Task 5's top scenarios at 5% (t5_value.sd5 footing: one simulated
patient, monthly record, therapy on), checked by profile likelihood over the
24 intake seeds and by Monte Carlo MAP. Every measurement that the FIM treats
as information (prior on k_e, load bias, flux reading) is simulated as a noisy
measurement whose centre is drawn about the truth, which is the reading under
which data-plus-prior information is the Fisher information.

Run from this directory:  python3 -W error t11_fim_check.py > out/t11_fim_check.txt
"""
import os
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")     # one BLAS thread: the machine is shared, and the
                                       # small dense solves here run ~100x faster this way
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, brentq

import core as C
import t2_cohort as T2
import t5_value as T5

T_START = time.time()
HERE = Path(__file__).resolve().parent
TH = C.TH_5
LN_TH = np.log(TH)
POP = C.KE_PRIOR
CHI95 = 3.841458820694124
LK_LO, LK_HI = np.log(5.0), np.log(5000.0)          # t2.joint_map's bounds on K_m
LKE_LO, LKE_HI = np.log(1e-3), np.log(500.0)        # and on k_e
TH_HI = 1e5
NS = (1, 2, 4, 8, 16, 32, 64)
REPS = 200
BOUNDS = (1e-2, 1e-3, 1e-5)


def banner(s):
    print("\n" + "=" * 100)
    print("  " + s)
    print("=" * 100)


def robust(e):
    e = np.asarray(e).ravel()
    return float(0.5 * (np.percentile(e, 84) - np.percentile(e, 16)))


def rmse(e):
    e = np.asarray(e).ravel()
    return float(np.sqrt(np.mean(e ** 2)))


def boot_ci(e, fn, B=2000, seed=0):
    """Percentile bootstrap 95% interval of fn over replicates (rows)."""
    e = np.asarray(e)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, e.shape[0], (B, e.shape[0]))
    v = np.array([fn(e[i]) for i in idx])
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def elapsed():
    return f"[{time.time() - T_START:.0f} s]"


# =====================================================================
#  The MAP by variable projection
# =====================================================================
class Problem:
    """
    N patients sharing (ln K_m, ln k_e); each has its own (theta_i, b_i).
    Residuals, exactly as t2_cohort.joint_map writes them:
        data   (y_ij - (theta_i + inc_ij) g(m_ij) - k_e m_ij - b_i) / sigma_ij
        prior  (ln k_e - kec) / ke_sd ;  (b_i - bc_i) / sb_i
    plus, optionally, one direct flux reading per patient,
        flux   (z_i - theta_i g(m0_i)) / sdf_i.
    theta_i is box-constrained to [th_lo, TH_HI] (th_lo = -inf: unconstrained);
    fix0 holds patient 0's theta fixed (for profiling).
    """

    def __init__(self, M, INC, Y, SG, sb, bc, kec, ke_sd, flux=None,
                 th_lo=1e-3, fix0=None):
        self.M, self.INC, self.Y = (np.atleast_2d(np.asarray(a, float)) for a in (M, INC, Y))
        SG = np.atleast_2d(np.asarray(SG, float))
        self.SG = np.broadcast_to(SG, self.M.shape)
        self.W = 1.0 / self.SG ** 2
        self.SW = self.W.sum(1)
        self.sb, self.bc = np.asarray(sb, float).ravel(), np.asarray(bc, float).ravel()
        self.kec, self.ke_sd = float(kec), float(ke_sd)
        self.flux = None if flux is None else tuple(np.asarray(a, float).ravel() for a in flux)
        self.lo, self.fix0 = th_lo, fix0
        self.n = self.M.shape[0]

    def _theta_b(self, a11, a12, a22, r1, r2):
        th = np.clip((a22 * r1 - a12 * r2) / (a11 * a22 - a12 ** 2), self.lo, TH_HI)
        if self.fix0 is not None:
            th = np.array(th, copy=True)
            th[..., 0] = self.fix0
        return th, (r2 - a12 * th) / a22

    def _stats(self, K):
        """Per-patient sufficient statistics at K (any leading shape of K)."""
        K = np.asarray(K, float)[..., None, None]
        G = self.M / (K + self.M)
        Yp = self.Y - self.INC * G
        W, M = self.W, self.M
        S = lambda X: (W * X).sum(-1)
        st = dict(YY=S(Yp * Yp), YM=S(Yp * M), MM=S(M * M), GY=S(G * Yp), GM=S(G * M),
                  Y=S(Yp), Mm=S(M), GG=S(G * G), G=S(G))
        if self.flux is not None:
            z, sdf, m0 = self.flux
            g0 = m0 / (K[..., 0] + m0)
            st.update(fg=g0, fz=z, fw=1.0 / sdf ** 2)
        return st

    def _cost_from_stats(self, st, ke, lke):
        """Profiled cost; st arrays (..., N), ke and lke broadcastable to (..., 1)."""
        a11 = st["GG"] + (st["fw"] * st["fg"] ** 2 if self.flux is not None else 0.0)
        a12 = st["G"]
        a22 = self.SW + 1.0 / self.sb ** 2
        r1 = st["GY"] - ke * st["GM"] + (st["fw"] * st["fg"] * st["fz"] if self.flux is not None else 0.0)
        r2 = st["Y"] - ke * st["Mm"] + self.bc / self.sb ** 2
        a11, r1, r2 = np.broadcast_arrays(a11, r1, r2)
        th, b = self._theta_b(a11, a12, a22, r1, r2)
        ZZ = st["YY"] - 2 * ke * st["YM"] + ke ** 2 * st["MM"]
        dat = (ZZ - 2 * th * (st["GY"] - ke * st["GM"]) - 2 * b * (st["Y"] - ke * st["Mm"])
               + th ** 2 * st["GG"] + 2 * th * b * st["G"] + b ** 2 * self.SW)
        if self.flux is not None:
            dat = dat + st["fw"] * (st["fz"] - th * st["fg"]) ** 2
        pri = ((b - self.bc) / self.sb) ** 2
        return 0.5 * (dat + pri).sum(-1) + 0.5 * ((lke[..., 0] - self.kec) / self.ke_sd) ** 2

    def inner(self, lk, lke):
        K, ke = np.exp(lk), np.exp(lke)
        G = self.M / (K + self.M)
        Z = self.Y - self.INC * G - ke * self.M
        W = self.W
        a11 = (W * G * G).sum(1)
        a12 = (W * G).sum(1)
        a22 = self.SW + 1 / self.sb ** 2
        r1 = (W * G * Z).sum(1)
        r2 = (W * Z).sum(1) + self.bc / self.sb ** 2
        if self.flux is not None:
            z, sdf, m0 = self.flux
            g0 = m0 / (K + m0)
            a11 = a11 + g0 ** 2 / sdf ** 2
            r1 = r1 + g0 * z / sdf ** 2
        th, b = self._theta_b(a11, a12, a22, r1, r2)
        return th, b, G, Z

    def resid(self, x):
        lk, lke = x
        th, b, G, Z = self.inner(lk, lke)
        R = [((Z - th[:, None] * G - b[:, None]) / self.SG).ravel()]
        if self.flux is not None:
            z, sdf, m0 = self.flux
            R.append((z - th * m0 / (np.exp(lk) + m0)) / sdf)
        R.append([(lke - self.kec) / self.ke_sd])
        R.append((b - self.bc) / self.sb)
        return np.concatenate(R)

    def cost(self, x):
        r = self.resid(x)
        return 0.5 * float(r @ r)

    def grid_cost(self, lks, lkes):
        st = self._stats(np.exp(lks))                        # (A, N)
        st = {k: (v[:, None, :] if np.ndim(v) == 2 else v) for k, v in st.items()}   # (A, 1, N)
        ke = np.exp(lkes)[None, :, None]                      # (1, B, 1)
        return self._cost_from_stats(st, ke, np.broadcast_to(lkes[None, :, None], (1, lkes.size, 1)))

    def fit(self, n_lk=41, n_lke=61, top=6, extra=(), full_box=False, width=6.0):
        lks = np.linspace(LK_LO, LK_HI, n_lk)
        if full_box:
            a, b = LKE_LO, LKE_HI
        else:
            a = max(LKE_LO, self.kec - width * self.ke_sd)
            b = min(LKE_HI, self.kec + width * self.ke_sd)
        lkes = np.linspace(a, b, n_lke)
        Cg = self.grid_cost(lks, lkes)
        idx = np.argsort(Cg.ravel())[:top]
        starts = ([np.array([lks[i // n_lke], lkes[i % n_lke]]) for i in idx]
                  + [np.asarray(e, float) for e in extra])
        lo, hi = np.array([LK_LO, LKE_LO]), np.array([LK_HI, LKE_HI])
        best = None
        for x0 in starts:
            x0 = np.clip(x0, lo + 1e-9, hi - 1e-9)
            r = least_squares(self.resid, x0, bounds=(lo, hi), method="trf",
                              xtol=1e-12, ftol=1e-12, gtol=1e-12, diff_step=1e-7)
            if best is None or r.cost < best.cost:
                best = r
        th, b, _, _ = self.inner(*best.x)
        return dict(lk=float(best.x[0]), lke=float(best.x[1]), th=th, b=b,
                    cost=float(best.cost), gridmin=float(Cg.min()))


# ------------------------------------------------------------- cohorts
def cohort_draw(n, k, rng, th=TH):
    """t2_cohort.mc's replicate k, drawn in exactly the same order."""
    recs, ys, sgs = [], [], []
    for i in range(n):
        r = C.record(th, seed=1000 * k + i)
        y, sg = C.noisy(r, 1000 * k + i)
        recs.append(r); ys.append(y); sgs.append(sg)
    ke_c = np.log(C.KE) + rng.normal(0, POP)
    b_c = np.array([rng.normal(0, C.LOAD_FRAC * r["L"]) for r in recs])
    return recs, ys, sgs, ke_c, b_c


def cohort_problem(recs, ys, sgs, kec, bc, **kw):
    return Problem([r["m"] for r in recs], [r["inc"] for r in recs], ys, sgs,
                   [C.LOAD_FRAC * r["L"] for r in recs], bc, kec, POP, **kw)


def jac_full(recs, sgs, th, lk, lke, drop=()):
    """
    Whitened Jacobian of every residual in (ln K_m, ln k_e, ln theta_i, b_i),
    priors included. drop: patients whose ln theta column is removed (theta at
    an active bound is not a free parameter).
    """
    n, nd = len(recs), recs[0]["m"].size
    K, ke = np.exp(lk), np.exp(lke)
    J = np.zeros((n * nd + 1 + n, 2 + 2 * n))
    for i, r in enumerate(recs):
        m, inc, s = r["m"], r["inc"], sgs[i]
        rows = slice(i * nd, (i + 1) * nd)
        J[rows, 0] = -(th[i] + inc) * K * m / (K + m) ** 2 / s
        J[rows, 1] = ke * m / s
        J[rows, 2 + i] = th[i] * m / (K + m) / s
        J[rows, 2 + n + i] = 1.0 / s
    J[n * nd, 1] = 1.0 / POP
    sb = np.array([C.LOAD_FRAC * r["L"] for r in recs])
    J[n * nd + 1:, 2 + n:] = np.diag(1.0 / sb)
    keep = [j for j in range(2 + 2 * n) if not (2 <= j < 2 + n and (j - 2) in drop)]
    return J[:, keep], keep


def prior_precision(recs):
    n = len(recs)
    P = np.zeros(2 + 2 * n)
    P[1] = POP ** -2
    P[2 + n:] = np.array([C.LOAD_FRAC * r["L"] for r in recs]) ** -2
    return np.diag(P)


# =====================================================================
#  PART 1.0: the estimator is t2's
# =====================================================================
def parse_t2_table_E():
    txt = (HERE / "out" / "t2_cohort.txt").read_text()
    blk = txt.split("E. Monte Carlo")[1]
    rows = {}
    for line in blk.splitlines():
        f = line.split()
        if len(f) == 8 and f[0].isdigit():
            rows[int(f[0])] = f
    return rows


def part1_0():
    banner("PART 1.0  THE ESTIMATOR: variable projection reproduces t2_cohort's joint MAP")
    printed = parse_t2_table_E()
    print("  (i) t2's table E rerun with this estimator (same seeds, same draws, same bounds):")
    print(f"     {'N':>4} {'reps':>5} {'MC RMSE':>8} {'MC robust':>10} {'median err':>11} "
          f"{'>4.5x off':>10} {'ke RMSE':>8}   t2 printed (RMSE, robust, median, >4.5x, ke)   match")
    for n, reps in ((1, 200), (4, 150), (16, 100), (64, 60)):
        rng = np.random.default_rng(0)
        e, ek = [], []
        for k in range(reps):
            recs, ys, sgs, kec, bc = cohort_draw(n, k, rng)
            f = cohort_problem(recs, ys, sgs, kec, bc).fit()
            e.append(np.log(f["th"][0]) - LN_TH)
            ek.append(f["lke"] - np.log(C.KE))
        e, ek = np.array(e), np.array(ek)
        mine = [f"{rmse(e):.3f}", f"{robust(e):.3f}", f"{np.median(e):+.3f}",
                f"{np.mean(np.abs(e) > np.log(4.5)):.1%}", f"{rmse(ek):.3f}"]
        t2 = printed[n]
        theirs = [t2[2], t2[3], t2[4], t2[5], t2[6]]
        ok = mine == theirs
        assert ok, (n, mine, theirs)
        print(f"     {n:4d} {reps:5d} {mine[0]:>8} {mine[1]:>10} {mine[2]:>11} {mine[3]:>10} "
              f"{mine[4]:>8}   {', '.join(theirs):44}   {'yes' if ok else 'NO'}")

    print("\n  (ii) against t2_cohort.joint_map itself on the same replicate (it draws the same")
    print("       prior centres from an identically seeded generator):")
    print(f"     {'N':>4} {'reps':>5} {'max |d ln theta_1|':>19} {'max |d ln k_e|':>15} {'t2 s/fit':>9}")
    for n, reps in ((1, 20), (2, 20), (4, 20), (8, 20), (16, 20), (32, 20), (64, 20)):
        ra, rb = np.random.default_rng(0), np.random.default_rng(0)
        d1, d2, tt = 0.0, 0.0, 0.0
        for k in range(reps):
            recs, ys, sgs, kec, bc = cohort_draw(n, k, ra)
            t0 = time.time()
            lth, lke = T2.joint_map([TH] * n, recs, ys, sgs, rb)
            tt += time.time() - t0
            f = cohort_problem(recs, ys, sgs, kec, bc).fit()
            d1 = max(d1, abs(lth[0] - np.log(f["th"][0])))
            d2 = max(d2, abs(lke - f["lke"]))
        print(f"     {n:4d} {reps:5d} {d1:19.2e} {d2:15.2e} {tt / reps:9.2f}")
    print(f"  {elapsed()}")


# =====================================================================
#  PART 1 (a, b): the Monte Carlo
# =====================================================================
def mc_cohort(n, reps=REPS):
    rng = np.random.default_rng(0)
    out = {key: [] for key in ("e3", "e2", "e5", "eu", "ec", "ke3", "keu", "fim", "fimke",
                               "fim_fixc", "gn", "lin", "linke", "nb3", "act", "same",
                               "rel3")}
    for k in range(reps):
        recs, ys, sgs, kec, bc = cohort_draw(n, k, rng)
        mk = lambda **kw: cohort_problem(recs, ys, sgs, kec, bc, **kw)
        f3 = mk().fit()
        f2 = mk(th_lo=1e-2).fit()
        f5 = mk(th_lo=1e-5).fit()
        fu = mk(th_lo=-np.inf).fit()
        fc = cohort_problem(recs, ys, sgs, np.log(C.KE), np.zeros(n)).fit()
        out["e3"].append(np.log(f3["th"]) - LN_TH)
        out["rel3"].append(f3["th"] / TH - 1)
        out["e2"].append(np.log(f2["th"]) - LN_TH)
        out["e5"].append(np.log(f5["th"]) - LN_TH)
        out["eu"].append(fu["th"] / TH - 1)                  # linear: may be < -1
        out["ec"].append(np.log(fc["th"][0]) - LN_TH)
        out["ke3"].append(f3["lke"] - np.log(C.KE))
        out["keu"].append(fu["lke"] - np.log(C.KE))
        out["nb3"].append(int(np.sum(f3["th"] <= 1e-3 * (1 + 1e-9))))
        act = bool(np.any(fu["th"] < 1e-3))                    # constraint would bind
        out["act"].append(act)
        if not act:
            out["same"].append(abs(np.log(f3["th"][0]) - np.log(fu["th"][0])))
        # information at the truth, on this replicate's own records
        J0, _ = jac_full(recs, sgs, np.full(n, TH), np.log(C.KM0), np.log(C.KE))
        F0 = J0.T @ J0
        P0 = np.linalg.inv(F0)
        out["fim"].append(np.sqrt(np.diag(P0)[2:2 + n]))
        out["fimke"].append(np.sqrt(P0[1, 1]))
        S = P0 @ (F0 - prior_precision(recs)) @ P0          # centres held at the truth
        out["fim_fixc"].append(np.sqrt(S[2, 2]))
        # linear-Gaussian surrogate: one Gauss-Newton step from the truth on the same draws
        sb = np.array([C.LOAD_FRAC * r["L"] for r in recs])
        y0 = np.array([(TH + r["inc"]) * r["m"] / (C.KM0 + r["m"]) + C.KE * r["m"] for r in recs])
        w = np.concatenate([((np.array(ys) - y0) / np.array(sgs)).ravel(),
                            [(kec - np.log(C.KE)) / POP], bc / sb])
        d = P0 @ (J0.T @ w)
        out["lin"].append(d[2:2 + n])
        out["linke"].append(d[1])
        # Gauss-Newton information at the MAP (active bounds removed)
        drop = [i for i in range(n) if f3["th"][i] <= 1e-3 * (1 + 1e-9)]
        if 0 not in drop:
            J1, keep = jac_full(recs, sgs, f3["th"], f3["lk"], f3["lke"], drop=drop)
            out["gn"].append(np.sqrt(np.linalg.inv(J1.T @ J1)[keep.index(2), keep.index(2)]))
    return {k: np.array(v) for k, v in out.items()}


def part1_ab():
    banner("PART 1 (a, b)  THE COHORT MONTE CARLO, N = 1 ... 64, 200 replicates each")
    print("  Footing: t2_cohort.mc exactly (omega = 0, all patients at 5%, priors as pseudo-")
    print("  observations with centres drawn about the truth, seeds 1000k + i). Patient 0 is")
    print("  t2's target; 'pooled' uses all N exchangeable patients of every replicate.")
    R, FIMT2 = {}, {}
    for n in NS:
        R[n] = mc_cohort(n)
        FIMT2[n] = T2.cohort(*T2.disease_cohort(n), 0.0)
        print(f"     N = {n:2d} done {elapsed()}")

    print("\n  (a) THE ARBITRARY LOWER BOUND ON theta. Patient 0, same replicates, three bounds")
    print(f"     {'N':>3} {'bound':>7} {'RMSE':>7} {'robust':>7} {'median':>7} {'at bound':>9}"
          f" {'pooled RMSE':>12} {'pooled robust':>14} {'ln(bound/theta)':>16}")
    for n in NS:
        for bd, key in zip(BOUNDS, ("e2", "e3", "e5")):
            e = R[n][key]
            atb = np.mean(e[:, 0] <= np.log(bd / TH) + 1e-6)
            print(f"     {n:3d} {bd:7.0e} {rmse(e[:, 0]):7.3f} {robust(e[:, 0]):7.3f} "
                  f"{np.median(e[:, 0]):+7.3f} {atb:9.1%} {rmse(e):12.3f} {robust(e):14.3f} "
                  f"{np.log(bd / TH):16.2f}")

    print("\n  (b-i) FIM against the constrained MAP (bound 1e-3, t2's), patient 0")
    print(f"     {'N':>3} {'FIM t2':>7} {'FIM own':>8} {'RMSE':>7} {'95% CI':>15} {'robust':>7}"
          f" {'95% CI':>15} {'median':>7} {'mean':>7} {'at bnd':>7} {'>4.5x':>6} {'robust/FIM':>10}")
    for n in NS:
        r = R[n]
        e0 = r["e3"][:, 0]
        fo = float(np.sqrt(np.mean(r["fim"][:, 0] ** 2)))
        ci_r = boot_ci(e0, rmse)
        ci_b = boot_ci(e0, robust)
        print(f"     {n:3d} {FIMT2[n][0]:7.3f} {fo:8.3f} {rmse(e0):7.3f} [{ci_r[0]:5.2f}, {ci_r[1]:5.2f}]"
              f" {robust(e0):7.3f} [{ci_b[0]:5.2f}, {ci_b[1]:5.2f}] {np.median(e0):+7.3f} {np.mean(e0):+7.3f}"
              f" {np.mean(e0 <= np.log(1e-3 / TH) + 1e-6):7.1%}"
              f" {np.mean(np.abs(e0) > np.log(4.5)):6.1%} {robust(e0) / fo:10.2f}")
    print("     FIM t2: t2_cohort.cohort on seeds 0..N-1 (t2's column); FIM own: rms over the")
    print("     replicates of the FIM on each replicate's own records. 95% CI: bootstrap over replicates.")

    print("\n  (b-ii) pooled over all N patients, constrained MAP, and the same quantity in LINEAR theta")
    print("     (theta_hat/theta - 1, which is bounded below by -1 and so does not depend on the bound)")
    print(f"     {'N':>3} {'FIM own':>8} {'RMSE ln':>8} {'robust ln':>10} {'median ln':>10} {'cohorts w/ a':>13}"
          f" {'RMSE lin':>9} {'mean lin':>9} {'median lin':>11}")
    print(f"     {'':>3} {'':>8} {'':>8} {'':>10} {'':>10} {'pt at bound':>13}")
    for n in NS:
        r = R[n]
        fo = float(np.sqrt(np.mean(r["fim"] ** 2)))
        print(f"     {n:3d} {fo:8.3f} {rmse(r['e3']):8.3f} {robust(r['e3']):10.3f} {np.median(r['e3']):+10.3f}"
              f" {np.mean(r['nb3'] > 0):13.1%} {rmse(r['rel3']):9.3f} {np.mean(r['rel3']):+9.3f}"
              f" {np.median(r['rel3']):+11.3f}")

    print("\n  (b-iii) THE SAME DRAWS WITHOUT THE POSITIVITY CONSTRAINT, AND A LINEAR-GAUSSIAN SURROGATE")
    print("     unconstrained: theta_i real-valued (no bound); its error is theta_hat/theta - 1 and the")
    print("     FIM's linear-theta sd equals its ln-theta sd at the truth. surrogate: one Gauss-Newton")
    print("     step from the truth on the same noise and prior draws, exactly linear-Gaussian.")
    print(f"     {'N':>3} {'FIM own':>8} {'uncon sd':>9} {'uncon rob':>10} {'uncon mean':>11} {'theta<0':>8}"
          f" {'surr sd':>8} {'surr rob':>9} {'constraint':>11} {'max|d| where':>13}")
    print(f"     {'':>3} {'':>8} {'pooled':>9} {'pooled':>10} {'pooled':>11} {'pooled':>8}"
          f" {'pooled':>8} {'pooled':>9} {'binds':>11} {'it does not':>13}")
    for n in NS:
        r = R[n]
        fo = float(np.sqrt(np.mean(r["fim"] ** 2)))
        same = r["same"]
        print(f"     {n:3d} {fo:8.3f} {np.std(r['eu']):9.3f} {robust(r['eu']):10.3f} {np.mean(r['eu']):+11.3f}"
              f" {np.mean(r['eu'] < -1):8.1%} {np.std(r['lin']):8.3f} {robust(r['lin']):9.3f}"
              f" {np.mean(r['act']):11.1%} {(same.max() if same.size else float('nan')):13.1e}")
    print("     constraint binds: replicates where the unconstrained optimum has some theta_i < 1e-3;")
    print("     in every other replicate the constrained MAP must equal the unconstrained one.")

    print("\n  (b-iv) k_e, the shared parameter: FIM against the MC")
    print(f"     {'N':>3} {'FIM sd':>7} {'+1.96 sd':>9} {'MC RMSE':>8} {'P2.5':>7} {'P50':>7} {'P97.5':>7}"
          f" {'uncon RMSE':>11} {'uncon P2.5':>11} {'uncon P97.5':>12} {'surr sd':>8}")
    for n in NS:
        r = R[n]
        fk = float(np.sqrt(np.mean(r["fimke"] ** 2)))
        p = np.percentile(r["ke3"], [2.5, 50, 97.5])
        pu = np.percentile(r["keu"], [2.5, 97.5])
        print(f"     {n:3d} {fk:7.3f} {1.96 * fk:+9.3f} {rmse(r['ke3']):8.3f} {p[0]:+7.3f} {p[1]:+7.3f} {p[2]:+7.3f}"
              f" {rmse(r['keu']):11.3f} {pu[0]:+11.3f} {pu[1]:+12.3f} {np.std(r['linke']):8.3f}")

    print("\n  (b-v) LIKE FOR LIKE? Prior centres drawn about the truth (as above) against centres held")
    print("     AT the truth (core.profile's and t1's convention). With centres held, the right linear")
    print("     comparison is the sandwich inv(F) (F - P_prior) inv(F), not inv(F).")
    print(f"     {'N':>3} {'FIM':>7} {'sandwich':>9} {'MC robust, drawn':>17} {'MC robust, held':>16} {'median held':>12}")
    for n in NS:
        r = R[n]
        fo = float(np.sqrt(np.mean(r["fim"][:, 0] ** 2)))
        fs = float(np.sqrt(np.mean(r["fim_fixc"] ** 2)))
        print(f"     {n:3d} {fo:7.3f} {fs:9.3f} {robust(r['e3'][:, 0]):17.3f} {robust(r['ec']):16.3f}"
              f" {np.median(r['ec']):+12.3f}")

    print("\n  (b-vi) the Gauss-Newton information at the MAP instead of at the truth, patient 0")
    print(f"     {'N':>3} {'FIM at truth':>13} {'GN at MAP: median':>18} {'P10':>7} {'P90':>7} {'n':>5}")
    for n in NS:
        r = R[n]
        g = r["gn"]
        print(f"     {n:3d} {np.sqrt(np.mean(r['fim'][:, 0] ** 2)):13.3f} {np.median(g):18.3f}"
              f" {np.percentile(g, 10):7.3f} {np.percentile(g, 90):7.3f} {g.size:5d}")
    print(f"  {elapsed()}")
    return R, FIMT2


def part1_global():
    banner("PART 1 (b)  IS THE MAP THE GLOBAL OPTIMUM? Dense full-box grid and 20 starts")
    print("  Standard fit: 41 x 61 grid in (ln K_m, ln k_e), k_e within 6 prior sd of its centre,")
    print("  6 local refinements. Check: 161 x 241 grid over the WHOLE box of t2's bounds")
    print("  (K_m 5..5000, k_e 1e-3..500), 20 refinements. (theta_i, b_i) are exact for any (K_m, k_e).")
    print(f"     {'N':>3} {'reps':>5} {'max cost gain':>14} {'gains > 1e-6':>13} {'max |d ln theta_1|':>19}")
    for n in NS:
        rng = np.random.default_rng(0)
        gain, cnt, dth = 0.0, 0, 0.0
        reps = 30
        for k in range(reps):
            recs, ys, sgs, kec, bc = cohort_draw(n, k, rng)
            p = cohort_problem(recs, ys, sgs, kec, bc)
            f = p.fit()
            g = p.fit(n_lk=161, n_lke=241, top=20, full_box=True)
            gain = max(gain, f["cost"] - g["cost"])
            cnt += (f["cost"] - g["cost"]) > 1e-6
            dth = max(dth, abs(np.log(f["th"][0]) - np.log(g["th"][0])))
        print(f"     {n:3d} {reps:5d} {gain:14.2e} {cnt:13d} {dth:19.2e}")
    print(f"  {elapsed()}")


# ------------------------------------------------------------- profiles
def pl_interval(make, lt_grid, cmap, lt_map):
    """
    95% profile-likelihood interval in ln theta_0 from a grid then brentq.
    make(fix) returns a Problem with patient 0's theta fixed. Returns
    (lo, hi, lo_open, hi_open, profile array).
    """
    thr = CHI95 / 2.0
    prev = None
    prof = []
    for lt in lt_grid:
        p = make(np.exp(lt))
        f = p.fit(top=3, extra=() if prev is None else (prev,))
        prev = np.array([f["lk"], f["lke"]])
        prof.append(f["cost"])
    prof = np.array(prof)
    cmin = min(cmap, prof.min())
    dev = prof - cmin - thr

    def pf(lt):
        return make(np.exp(lt)).fit(top=3)["cost"] - cmin - thr

    inside = np.where(dev <= 0)[0]
    if inside.size == 0:
        j = int(np.argmin(np.abs(lt_grid - lt_map)))
        inside = np.array([j])
    i0, i1 = inside.min(), inside.max()
    lo_open, hi_open = i0 == 0, i1 == lt_grid.size - 1
    lo = lt_grid[0] if lo_open else brentq(pf, lt_grid[i0 - 1], lt_grid[i0], xtol=1e-4)
    hi = lt_grid[-1] if hi_open else brentq(pf, lt_grid[i1], lt_grid[i1 + 1], xtol=1e-4)
    gaps = np.any(dev[i0:i1 + 1] > 0)
    return lo, hi, lo_open, hi_open, gaps, prof


def part1_profile(ncoh=12):
    banner("PART 1 (b)  PROFILE LIKELIHOOD OF theta_1 IN THE JOINT FIT, N = 64")
    n = 64
    grid = np.r_[np.linspace(-11.5, -3.5, 9), np.arange(-3.4, 2.55, 0.1)] + LN_TH
    print(f"  {ncoh} cohorts (t2's replicates k = 0..{ncoh - 1}); every other parameter, including the other")
    print("  63 thetas, re-optimised at each fixed theta_1. 95%: drop of chi2_1(0.95)/2 in -log L.")
    print("  ln units relative to the truth. FIM interval: MAP +- 1.96 FIM sd (own records).")
    print(f"     {'k':>3} {'MAP':>7} {'PL lo':>7} {'PL hi':>7} {'PL width':>9} {'FIM width':>10}"
          f" {'GN@MAP width':>13} {'PL/FIM':>7} {'up/down':>8} {'covers':>7} {'lin half':>9} {'FIM sd':>7}")
    rng = np.random.default_rng(0)
    rows = []
    for k in range(ncoh):
        recs, ys, sgs, kec, bc = cohort_draw(n, k, rng)
        f = cohort_problem(recs, ys, sgs, kec, bc).fit()
        lt_map = np.log(f["th"][0])
        J0, _ = jac_full(recs, sgs, np.full(n, TH), np.log(C.KM0), np.log(C.KE))
        fim_sd = np.sqrt(np.linalg.inv(J0.T @ J0)[2, 2])
        drop = [i for i in range(n) if f["th"][i] <= 1e-3 * (1 + 1e-9)]
        J1, keep = jac_full(recs, sgs, f["th"], f["lk"], f["lke"], drop=drop)
        gn_sd = np.sqrt(np.linalg.inv(J1.T @ J1)[keep.index(2), keep.index(2)])
        make = lambda fix: cohort_problem(recs, ys, sgs, kec, bc, fix0=fix)
        lo, hi, lo_o, hi_o, gaps, _ = pl_interval(make, grid, f["cost"], lt_map)
        w = hi - lo
        cov = lo <= LN_TH <= hi
        lin = (np.exp(hi) - np.exp(lo)) / (3.92 * TH)
        rows.append((w, 3.92 * fim_sd, 3.92 * gn_sd, (hi - lt_map) / (lt_map - lo), cov, lo_o, hi_o, gaps,
                     lin, fim_sd))
        print(f"     {k:3d} {lt_map - LN_TH:+7.3f} {lo - LN_TH:+7.2f}{'<' if lo_o else ' '}{hi - LN_TH:+7.2f}"
              f"{'>' if hi_o else ' '}{w:8.3f} {3.92 * fim_sd:10.3f} {3.92 * gn_sd:13.3f} {w / (3.92 * fim_sd):7.2f}"
              f" {(hi - lt_map) / (lt_map - lo):8.2f} {'yes' if cov else 'no':>7} {lin:9.3f} {fim_sd:7.3f}"
              f"{'  (gap)' if gaps else ''}")
    rows = np.array(rows, dtype=float)
    print(f"     median PL width {np.median(rows[:, 0]):.3f}, FIM width {np.median(rows[:, 1]):.3f}, "
          f"ratio {np.median(rows[:, 0] / rows[:, 1]):.2f} (range {np.min(rows[:, 0] / rows[:, 1]):.2f}"
          f"-{np.max(rows[:, 0] / rows[:, 1]):.2f}); GN-at-MAP width {np.median(rows[:, 2]):.3f}; "
          f"upper/lower half-width median {np.median(rows[:, 3]):.2f}")
    print(f"     linear theta: PL (theta_hi - theta_lo)/(3.92 theta) median {np.median(rows[:, 8]):.3f} against the FIM sd "
          f"{np.median(rows[:, 9]):.3f} (ratio {np.median(rows[:, 8] / rows[:, 9]):.2f})")
    print("     '<' marks an interval still open at the grid floor ln theta - ln 150 = -11.5 (theta = 0.0015).")
    print(f"     covers the truth in {int(rows[:, 4].sum())} of {ncoh}; open at the grid floor in "
          f"{int(rows[:, 5].sum())}, at the top in {int(rows[:, 6].sum())}; disjoint pieces in {int(rows[:, 7].sum())}")
    print(f"  {elapsed()}")
    return rows


def part1_c(reps=2000):
    banner("PART 1 (c)  VALIDATING THE MC MACHINERY: the linear-Gaussian surrogate, 2000 replicates")
    print("  Same draws (t2's records, noise, prior centres) fed through the model LINEARISED at the")
    print("  truth, with no bounds. Its estimate is exactly Gaussian with covariance inv(FIM), so the")
    print("  whitened error (error / FIM sd on that replicate's records) must have sd 1.")
    print(f"     {'N':>3} {'reps':>5} {'FIM rms':>8} {'MC sd':>7} {'whitened sd th_1':>17} {'95% range for 1':>16}"
          f" {'whitened sd k_e':>16} {'mean, whitened':>15}")
    out = {}
    for n, k0, k1, seed in ((1, 0, reps, 0), (16, 0, reps, 0), (64, 0, reps, 0), (64, reps, 4 * reps, 1)):
        rng = np.random.default_rng(seed)
        z, zk, e, fs = [], [], [], []
        for k in range(k0, k1):
            recs, ys, sgs, kec, bc = cohort_draw(n, k, rng)
            J0, _ = jac_full(recs, sgs, np.full(n, TH), np.log(C.KM0), np.log(C.KE))
            P0 = np.linalg.inv(J0.T @ J0)
            sb = np.array([C.LOAD_FRAC * r["L"] for r in recs])
            y0 = np.array([(TH + r["inc"]) * r["m"] / (C.KM0 + r["m"]) + C.KE * r["m"] for r in recs])
            w = np.concatenate([((np.array(ys) - y0) / np.array(sgs)).ravel(),
                                [(kec - np.log(C.KE)) / POP], bc / sb])
            d = P0 @ (J0.T @ w)
            z.append(d[2] / np.sqrt(P0[2, 2])); zk.append(d[1] / np.sqrt(P0[1, 1]))
            e.append(d[2]); fs.append(np.sqrt(P0[2, 2]))
        z, zk = np.array(z), np.array(zk)
        nr = k1 - k0
        half = 1.96 / np.sqrt(2 * nr)
        if (n, k0) == (64, reps):
            zz = np.r_[out[(64, 0)][2], z]
            hz = 1.96 / np.sqrt(2 * zz.size)
            print(f"     {n:3d} {nr:5d} {np.sqrt(np.mean(np.square(fs))):8.3f} {np.std(e):7.3f} {np.std(z):17.3f}"
                  f"  [{1 - half:.3f}, {1 + half:.3f}] {np.std(zk):16.3f} {np.mean(z):+15.3f}"
                  f"   (replicates {k0}..{k1 - 1}, new prior-centre stream)")
            print(f"     {n:3d} {zz.size:5d} {'':8} {'':7} {np.std(zz):17.3f}  [{1 - hz:.3f}, {1 + hz:.3f}]"
                  f"   (both N = 64 blocks together)")
            out[(64, "all")] = (float(np.std(zz)), hz, zz)
        else:
            print(f"     {n:3d} {nr:5d} {np.sqrt(np.mean(np.square(fs))):8.3f} {np.std(e):7.3f} {np.std(z):17.3f}"
                  f"  [{1 - half:.3f}, {1 + half:.3f}] {np.std(zk):16.3f} {np.mean(z):+15.3f}")
            out[(n, k0)] = (float(np.std(z)), half, z)
    print("  The 95% range is the Gaussian sampling range of an sd estimated from that many replicates.")
    print(f"  {elapsed()}")
    return out


# =====================================================================
#  PART 2: Task 5's scenarios
# =====================================================================
SCEN = [
    ("current design (reference)", dict()),
    ("(c) flux reading to 10%", dict(flux_Y=0.10)),
    ("(c) flux reading to 20%", dict(flux_Y=0.20)),
    ("(a) k_e exact", dict(ke_F=1.0)),
    ("(a) k_e to 5%", dict(ke_F=1.05)),
    ("(d) k_e exact + load exact", dict(ke_F=1.0, load_X=0.0)),
    ("(d) flux 20% + k_e 1.2 + load 2%", dict(flux_Y=0.20, ke_F=1.2, load_X=0.02)),
]


def scen_settings(ke_F=None, load_X=None, flux_Y=None):
    """t5_value.sd5's reading of each measurement."""
    ke_sd = POP if ke_F is None else (np.log(ke_F) if ke_F > 1 else 1e-4)
    lf = C.LOAD_FRAC if load_X is None else max(load_X, 1e-4)
    return ke_sd, lf, flux_Y


def scen_record(seed, kw):
    return T5.base_record(seed, baseline=kw.get("flux_Y") is not None)


def scen_fim(seed, kw):
    """Exactly t5_value.sd5's one(s)."""
    ke_sd, lf, Y = scen_settings(**kw)
    r = scen_record(seed, kw)
    extra = []
    if Y is not None:
        row, F, eps = T5.flux_row(r)
        extra.append((row, F * np.hypot(Y, eps * C.ASSAY_CV)))
    return C.sd_th(C.fim(r, extra=extra, ke_prior=ke_sd, load_frac=lf)), r, extra


def scen_jac_fim(r, kw):
    """The MAP's own model, differentiated at the truth: must equal C.fim."""
    ke_sd, lf, Y = scen_settings(**kw)
    sg = C.sigmas(r, "repo", load_frac=lf)
    th, K, ke, m, inc = r["theta"], r["Km"], r["ke"], r["m"], r["inc"]
    rows = [np.stack([th * m / (K + m), -(th + inc) * K * m / (K + m) ** 2, np.ones_like(m), ke * m], 1) / sg[:, None]]
    if Y is not None:
        m0 = r["m_bar"]
        F = th * m0 / (K + m0)
        sdf = F * np.hypot(Y, C.eps_of(m0) * C.ASSAY_CV)
        rows.append(np.array([[F, -th * K * m0 / (K + m0) ** 2, 0, 0]]) / sdf)
    rows.append(np.array([[0, 0, 1.0 / (lf * r["L"]), 0]]))
    rows.append(np.array([[0, 0, 0, 1.0 / ke_sd]]))
    J = np.concatenate(rows)
    return J.T @ J


def scen_problem(r, kw, rng, fix0=None):
    """One simulated measurement set for this scenario; centres drawn about the truth."""
    ke_sd, lf, Y = scen_settings(**kw)
    e = rng.standard_normal(26)
    ef, eke, eb = rng.standard_normal(3)
    e = e if r["m"].size == 26 else e[1:]
    sg = C.sigmas(r, "repo", load_frac=lf)
    y = (r["theta"] + r["inc"]) * r["m"] / (r["Km"] + r["m"]) + r["ke"] * r["m"] + sg * e
    flux = None
    if Y is not None:
        m0 = r["m_bar"]
        F = TH * m0 / (C.KM0 + m0)
        sdf = F * np.hypot(Y, C.eps_of(m0) * C.ASSAY_CV)
        flux = (np.array([F + sdf * ef]), np.array([sdf]), np.array([m0]))
    sb = lf * r["L"]
    return Problem(r["m"][None, :], r["inc"][None, :], y[None, :], sg[None, :], [sb], [sb * eb],
                   np.log(C.KE) + ke_sd * eke, ke_sd, flux=flux, fix0=fix0)


def part2(reps=240):
    banner("PART 2  TASK 5's TOP SCENARIOS AT 5% (t5_value.sd5 footing), FIM vs PROFILE vs MC MAP")
    print("  One simulated patient (theta = 150, K_m = 100, k_e = 5), monthly record, dose train on.")
    print("  Flux scenarios add the untreated baseline draw and one direct flux reading (t5's row).")
    print("  'exact' is t5's 1e-4. Every measurement is simulated with its centre drawn about the")
    print("  truth (pseudo-observation). MC: 240 replicates = 10 per intake seed 0..23, common random")
    print("  numbers across scenarios. Profile: one measurement set per intake seed 0..23.")

    print("\n  (i) the FIM: t5_value.sd5 recomputed, and the MAP's own model differentiated at the truth")
    print(f"     {'scenario':34} {'sd5':>8} {'here':>8} {'own-Jacobian FIM, max rel diff':>31}")
    FIM = {}
    for lab, kw in SCEN:
        vals, dmax = [], 0.0
        for s in C.SEEDS:
            v, r, extra = scen_fim(s, kw)
            vals.append(v)
            ke_sd, lf, _ = scen_settings(**kw)
            F1 = C.fim(r, extra=extra, ke_prior=ke_sd, load_frac=lf)
            F2 = scen_jac_fim(r, kw)
            dmax = max(dmax, float(np.max(np.abs(F1 - F2) / (np.abs(F1) + 1e-300))))
        med = float(np.median(vals))
        ref = T5.sd5(**kw)
        assert abs(med - ref) <= 1e-12 * ref, (lab, med, ref)
        FIM[lab] = (med, np.array(vals))
        print(f"     {lab:34} {ref:8.4f} {med:8.4f} {dmax:31.1e}")

    print("\n  (ii) Monte Carlo MAP, ln theta_hat - ln theta")
    print(f"     {'scenario':34} {'FIM':>7} {'FIM rms':>8} {'RMSE':>7} {'95% CI':>15} {'robust':>7} {'95% CI':>15}"
          f" {'median':>7} {'at bnd':>7} {'RMSE/FIM':>9} {'rob/FIM':>8} {'lin RMSE':>8} {'lin mean':>8} {'min err':>8}")
    MC = {}
    for lab, kw in SCEN:
        e, fims = [], []
        for j in range(reps):
            s = j % 24
            r = scen_record(s, kw)
            p = scen_problem(r, kw, np.random.default_rng(70_000 + j))
            f = p.fit()
            e.append(np.log(f["th"][0]) - LN_TH)
            fims.append(FIM[lab][1][s])
        e = np.array(e)
        lin = np.exp(e) - 1.0
        frms = float(np.sqrt(np.mean(np.square(fims))))
        ci_r, ci_b = boot_ci(e, rmse), boot_ci(e, robust)
        MC[lab] = dict(rmse=rmse(e), robust=robust(e), med=float(np.median(e)), ci_r=ci_r, ci_b=ci_b,
                       frms=frms, atb=float(np.mean(e <= np.log(1e-3 / TH) + 1e-6)),
                       lin=rmse(lin), lin_mean=float(np.mean(lin)), emin=float(e.min()),
                       ci_l=boot_ci(lin, rmse))
        m = MC[lab]
        print(f"     {lab:34} {FIM[lab][0]:7.4f} {frms:8.4f} {m['rmse']:7.4f} [{ci_r[0]:6.4f},{ci_r[1]:6.4f}]"
              f" {m['robust']:7.4f} [{ci_b[0]:6.4f},{ci_b[1]:6.4f}] {m['med']:+7.4f} {m['atb']:7.1%}"
              f" {m['rmse'] / frms:9.3f} {m['robust'] / frms:8.3f} {m['lin']:8.4f} {m['lin_mean']:+8.4f}"
              f" {m['emin']:+8.2f}")
    print("     lin: theta_hat/theta - 1, which the bound cannot reach below -1. min err: the worst ln error.")
    print(f"  {elapsed()}")

    print("\n  (iii) profile likelihood, 95% interval for ln theta (k_e, K_m, b re-optimised), 24 seeds")
    print("       ln units relative to the truth; FIM width = 2 x 1.96 x the seed's own FIM sd")
    print(f"     {'scenario':34} {'med lo':>7} {'med hi':>7} {'med width':>10} {'FIM width':>10} {'PL/FIM':>7}"
          f" {'range':>12} {'down/1.96':>10} {'up/1.96':>8} {'lin half':>9} {'covers':>7} {'open':>5} {'med MAP':>8}")
    PL = {}
    grid = np.r_[np.linspace(-11.5, -3.5, 9), np.arange(-3.4, 2.05, 0.1)] + LN_TH
    for lab, kw in SCEN[1:]:
        rows = []
        for s in C.SEEDS:
            r = scen_record(s, kw)
            f = scen_problem(r, kw, np.random.default_rng(90_000 + s)).fit()
            lt_map = np.log(f["th"][0])
            make = lambda fix: scen_problem(r, kw, np.random.default_rng(90_000 + s), fix0=fix)
            lo, hi, lo_o, hi_o, gaps, _ = pl_interval(make, grid, f["cost"], lt_map)
            fw = 3.92 * FIM[lab][1][s]
            rows.append((lo - LN_TH, hi - LN_TH, hi - lo, fw, (lt_map - lo) / 1.96, (hi - lt_map) / 1.96,
                         (np.exp(hi) - np.exp(lo)) / (3.92 * TH), lo <= LN_TH <= hi, lo_o or hi_o or gaps,
                         lt_map - LN_TH))
        rows = np.array(rows, dtype=float)
        ratio = rows[:, 2] / rows[:, 3]
        PL[lab] = dict(width=float(np.median(rows[:, 2])), ratio=float(np.median(ratio)),
                       half=float(np.median(rows[:, 2])) / 3.92, cov=int(rows[:, 7].sum()),
                       down=float(np.median(rows[:, 4])), up=float(np.median(rows[:, 5])),
                       lin=float(np.median(rows[:, 6])), nopen=int(rows[:, 8].sum()))
        print(f"     {lab:34} {np.median(rows[:, 0]):+7.3f} {np.median(rows[:, 1]):+7.3f} {np.median(rows[:, 2]):10.4f}"
              f" {np.median(rows[:, 3]):10.4f} {np.median(ratio):7.3f} {ratio.min():5.3f}-{ratio.max():5.3f}"
              f" {np.median(rows[:, 4]):10.4f} {np.median(rows[:, 5]):8.4f} {np.median(rows[:, 6]):9.4f}"
              f" {int(rows[:, 7].sum()):4d}/24 {int(rows[:, 8].sum()):5d} {np.median(rows[:, 9]):+8.3f}")
    print("     down/1.96, up/1.96: (MAP - lower end)/1.96 and (upper end - MAP)/1.96, each to set against")
    print("     the FIM sd. lin half: (theta_hi - theta_lo)/(3.92 theta), the interval's half-width in")
    print("     linear theta. open: interval still open at the grid floor (ln theta - ln 150 = -11.5,")
    print("     theta = 0.0015) or top (+2.0), or split in pieces.")
    print("     The current design's profile is in t1 (1d): no lower end at all, so it is not repeated.")
    print(f"  {elapsed()}")

    print("\n  (iv) SUMMARY: ranking, ratio to the FIM, and whether the FIM is a lower bound")
    labs = [l for l, _ in SCEN[1:]]
    fimv = np.array([FIM[l][0] for l in labs])
    mcr = np.array([MC[l]["rmse"] for l in labs])
    mcb = np.array([MC[l]["robust"] for l in labs])
    linr = np.array([MC[l]["lin"] for l in labs])
    plh = np.array([PL[l]["half"] for l in labs])
    pll = np.array([PL[l]["lin"] for l in labs])
    rk = lambda v: np.argsort(np.argsort(v)) + 1
    print(f"     {'scenario':34} {'FIM':>7} {'PL w/3.92':>10} {'PL lin':>7} {'MC RMSE':>8} {'MC robust':>10}"
          f" {'lin RMSE':>9}   ranks: FIM  PL  PLlin  RMSE  rob  lin")
    for i, l in enumerate(labs):
        print(f"     {l:34} {fimv[i]:7.4f} {plh[i]:10.4f} {pll[i]:7.4f} {mcr[i]:8.4f} {mcb[i]:10.4f} {linr[i]:9.4f}"
              f"   {rk(fimv)[i]:9d} {rk(plh)[i]:3d} {rk(pll)[i]:6d} {rk(mcr)[i]:5d} {rk(mcb)[i]:4d} {rk(linr)[i]:4d}")
    print("     PL w/3.92: median ln-theta profile width / 3.92; PL lin: median linear-theta half-width.")
    for nm, v in (("profile width", plh), ("profile, linear", pll), ("MC RMSE", mcr), ("MC robust", mcb),
                  ("linear RMSE", linr)):
        print(f"     ranking under {nm:16} identical to the FIM's: {'yes' if np.all(rk(fimv) == rk(v)) else 'no'}")
    for i, l in enumerate(labs):
        for j, l2 in enumerate(labs):
            if i < j and (fimv[i] - fimv[j]) * (mcr[i] - mcr[j]) < 0:
                print(f"     order flips between FIM and MC RMSE: '{l}' ({fimv[i]:.4f} / {mcr[i]:.4f}) vs "
                      f"'{l2}' ({fimv[j]:.4f} / {mcr[j]:.4f})")
    print("\n     Is the FIM sd a lower bound? Each MC measure / the FIM rms, with its bootstrap 95% CI;")
    print("     'below' means the whole CI is under 1 (the FIM is NOT a lower bound for that measure),")
    print("     'above' the whole CI over 1, 'at' a CI containing 1.")
    word = lambda lo, hi: "below" if hi < 1 else ("above" if lo > 1 else "at")
    for l in labs:
        m = MC[l]
        f = m["frms"]
        parts = []
        for nm, v, ci in (("ln RMSE", m["rmse"], m["ci_r"]), ("robust", m["robust"], m["ci_b"]),
                          ("lin RMSE", m["lin"], m["ci_l"])):
            parts.append(f"{nm} {v / f:.3f} [{ci[0] / f:.3f}, {ci[1] / f:.3f}] {word(ci[0] / f, ci[1] / f)}")
        print(f"     {l:34} " + "; ".join(parts))
    gate = [l for i, l in enumerate(labs) if fimv[i] <= C.GATE]
    gate_mc = [l for i, l in enumerate(labs) if mcr[i] <= C.GATE]
    gate_pl = [l for i, l in enumerate(labs) if plh[i] <= C.GATE]
    print(f"     at or under the 0.20 gate by FIM: {gate}")
    print(f"     by MC RMSE: {gate_mc}")
    print(f"     by profile half-width: {gate_pl}")
    ref = MC[SCEN[0][0]]
    print(f"     reference, the current design: FIM {FIM[SCEN[0][0]][0]:.4f}, MC RMSE {ref['rmse']:.3f}, "
          f"robust {ref['robust']:.3f}, {ref['atb']:.1%} of fits at the theta >= 1e-3 bound")
    return FIM, MC, PL


# =====================================================================
def verdict1(R, FIMT2, prof, surr):
    banner("PART 1  WHAT THE NUMBERS ABOVE SAY (every qualitative word below is computed)")
    cmp = lambda a, b, tol: "within" if abs(a / b - 1) <= tol else ("above" if a > b else "below")
    r1, r64 = R[1], R[64]
    f1 = float(np.sqrt(np.mean(r1["fim"] ** 2)))
    f64 = float(np.sqrt(np.mean(r64["fim"][:, 0] ** 2)))
    f64p = float(np.sqrt(np.mean(r64["fim"] ** 2)))
    tol = 1.96 / np.sqrt(2 * REPS)              # 95% sampling band of an sd estimated from REPS draws
    rb = [robust(r1[k][:, 0]) for k in ("e2", "e3", "e5")]
    rm = [rmse(r1[k][:, 0]) for k in ("e2", "e3", "e5")]
    print(f"  N = 1. Bounds 1e-2 / 1e-3 / 1e-5: robust {rb[0]:.2f} / {rb[1]:.2f} / {rb[2]:.2f}, "
          f"RMSE {rm[0]:.2f} / {rm[1]:.2f} / {rm[2]:.2f};")
    print(f"     at the 1e-3 bound: {np.mean(r1['e3'][:, 0] <= np.log(1e-3 / TH) + 1e-6):.1%} of fits. "
          f"Robust spread {'rises' if rb[2] > rb[1] > rb[0] else 'does not rise monotonically'} as the bound "
          f"is lowered (x {rb[2] / rb[0]:.2f} from 1e-2 to 1e-5).")
    su = float(np.std(r1["eu"]))
    print(f"     Unconstrained theta (same draws): sd(theta_hat/theta - 1) = {su:.3f} against the FIM's {f1:.3f} "
          f"({cmp(su, f1, tol)} the +-{tol:.0%} band), theta_hat < 0 in {np.mean(r1['eu'] < -1):.1%}.")
    e64, e64p = r64["e3"][:, 0], r64["e3"]
    t2_64 = parse_t2_table_E()[64]
    print(f"  N = 64. Patient-0 RMSE over {e64.size} replicates {rmse(e64):.3f} (t2 printed {t2_64[2]} from 60 replicates); pooled over "
          f"the 64 patients {rmse(e64p):.3f};")
    print(f"     pooled RMSE with bounds 1e-2 / 1e-3 / 1e-5: {rmse(r64['e2']):.3f} / {rmse(e64p):.3f} / "
          f"{rmse(r64['e5']):.3f}; cohorts with a patient at the 1e-3 bound: {np.mean(r64['nb3'] > 0):.1%}.")
    print(f"     So the ln-theta RMSE is {cmp(rmse(e64p), f64p, tol)} the FIM ({f64p:.3f}) and "
          f"{'depends' if rmse(r64['e5']) > 1.05 * rmse(r64['e2']) else 'does not depend'} on the bound.")
    rc, rcp = robust(e64), robust(e64p)
    print(f"     Robust spread: patient 0 {rc:.3f}, pooled {rcp:.3f}: {cmp(rcp, f64p, tol)} the FIM; "
          f"median error {np.median(e64p):+.3f} (pooled), linear-theta RMSE {rmse(r64['rel3']):.3f} "
          f"({cmp(rmse(r64['rel3']), f64p, tol)} the FIM).")
    su64, ru64, sl64 = float(np.std(r64["eu"])), robust(r64["eu"]), float(np.std(r64["lin"]))
    print(f"     Without the positivity constraint: sd {su64:.3f}, robust {ru64:.3f} ({cmp(su64, f64p, tol)} / "
          f"{cmp(ru64, f64p, tol)} the FIM); linear surrogate {sl64:.3f} ({cmp(sl64, f64p, tol)}).")
    k3, ku = np.percentile(r64["ke3"], 97.5), np.percentile(r64["keu"], 97.5)
    fk = 1.96 * float(np.sqrt(np.mean(r64["fimke"] ** 2)))
    print(f"     The constraint changes the estimate in {np.mean(r64['act']):.1%} of cohorts. k_e's 97.5% point: "
          f"{k3:+.3f} constrained, {ku:+.3f} unconstrained, FIM {fk:+.3f}; 2.5% point "
          f"{np.percentile(r64['ke3'], 2.5):+.3f} / {np.percentile(r64['keu'], 2.5):+.3f} / {-fk:+.3f}.")
    w = prof[:, 0] / prof[:, 1]
    print(f"     Profile likelihood, N = 64, {prof.shape[0]} cohorts: width / FIM width median {np.median(w):.2f} "
          f"(range {w.min():.2f}-{w.max():.2f}); upper/lower half-width median {np.median(prof[:, 3]):.2f}; "
          f"open at the floor in {int(prof[:, 5].sum())}.")
    keys = [(1, 0), (16, 0), (64, 0), (64, "all")]
    lab = {(1, 0): "N=1", (16, 0): "N=16", (64, 0): "N=64", (64, "all"): "N=64 all blocks"}
    print("  MACHINERY. Linear-Gaussian surrogate, whitened sd of theta_1: "
          + ", ".join(f"{lab[k]} {surr[k][0]:.3f} ({'inside' if abs(surr[k][0] - 1) <= surr[k][1] else 'OUTSIDE'}"
                      f" +-{surr[k][1]:.3f})" for k in keys) + ".")
    print("  WHICH BOUND. With the prior centres drawn about a fixed truth, inv(FIM) is the Cramer-Rao bound")
    print("     of the augmented experiment for UNBIASED estimators (the surrogate above is that case). Mean")
    print(f"     error in linear theta at N = 64: constrained MAP {np.mean(r64['rel3']):+.3f}, unconstrained "
          f"{np.mean(r64['eu']):+.3f}; at N = 1: {np.mean(r1['rel3']):+.3f} / {np.mean(r1['eu']):+.3f}.")
    print("     A van Trees (Bayesian) bound would need the truth drawn from a proper prior on every")
    print("     parameter; here the truth is fixed and ln theta_i, ln K_m carry no prior, so none is computed.")


def figure(R, FIM, MC):
    from plotstyle import plt, C1, C2, C3, C4, INK, MUTED, finish, save, plain
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw=dict(width_ratios=[1.0, 1.15]))
    a = ax[0]
    ns = np.array(NS, float)
    fimv = [np.sqrt(np.mean(R[n]["fim"] ** 2)) for n in NS]
    rob = [robust(R[n]["e3"]) for n in NS]
    rm = [rmse(R[n]["e3"]) for n in NS]
    rlo = [rmse(R[n]["e2"]) for n in NS]
    rhi = [rmse(R[n]["e5"]) for n in NS]
    unc = [np.std(R[n]["eu"]) for n in NS]
    ring = dict(mec="white", mew=1.5, ms=7, ls="none", zorder=3)
    a.plot(ns, fimv, "-", color=C1, lw=2, label="FIM sd (each replicate's own records)", zorder=2)
    a.vlines(ns * 1.06, rlo, rhi, color=C3, lw=1.2, zorder=2)
    a.plot(ns * 1.06, rm, "o", color=C3, label="MAP RMSE in ln θ (bar: bound 1e-2 to 1e-5)", **ring)
    a.plot(ns, rob, "D", color=C2, label="MAP robust spread (16-84%) in ln θ", **ring)
    a.plot(ns / 1.06, unc, "s", color=C4, label="no positivity constraint: sd of θ̂/θ − 1", **ring)
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xticks(ns); a.set_xlim(0.8, 80)
    a.set_xlabel("patients in the cohort, all at 5% of normal (k_e shared)")
    a.set_ylabel("spread of patient's estimate")
    a.set_title("(a) cohort fit: the FIM against the joint MAP")
    a.set_ylim(0.3, 60)
    a.legend(fontsize=7.5, frameon=False, loc="upper right")
    finish(a); plain(a)
    from matplotlib.ticker import NullFormatter
    a.xaxis.set_minor_formatter(NullFormatter())

    b = ax[1]
    labs = [l for l, _ in SCEN[1:]]
    y = np.arange(len(labs))[::-1]
    fv = [FIM[l][0] for l in labs]
    mr = [MC[l]["rmse"] for l in labs]
    mb = [MC[l]["robust"] for l in labs]
    lo = [MC[l]["ci_r"][0] for l in labs]
    hi = [MC[l]["ci_r"][1] for l in labs]
    b.hlines(y - 0.12, lo, hi, color=C3, lw=1.2, zorder=2)
    b.plot(mr, y - 0.12, "o", color=C3, label="MAP RMSE in ln θ (bar: bootstrap 95%)", **ring)
    b.plot(mb, y + 0.12, "D", color=C2, label="MAP robust spread", **ring)
    b.plot(fv, y, "|", color=C1, ms=16, mew=2.5, ls="none", label="FIM sd (median, 24 seeds)", zorder=4)
    b.axvline(C.GATE, color=INK, ls="--", lw=1)
    b.text(C.GATE * 1.02, y.max() + 0.55, "target 0.20", fontsize=8, color=INK, va="bottom")
    b.set_yticks(y)
    b.set_yticklabels([l.replace("(d) ", "").replace("(c) ", "").replace("(a) ", "") for l in labs], fontsize=8)
    b.set_xscale("log")
    b.set_ylim(-0.7, len(labs) - 0.2)
    b.set_xlabel(r"sd of ln $\hat\theta$, one patient at 5% of normal")
    b.set_title("(b) Task 5's scenarios: FIM against the MAP, 240 replicates")
    b.legend(fontsize=7.5, frameon=False, loc="lower right")
    finish(b); plain(b)
    fig.tight_layout()
    save(fig, "t11_fim_check")


def main():
    print("=" * 100)
    print("  FOLLOW-UP ITEM 6. THE FIM AGAINST THE LIKELIHOOD AND THE MONTE CARLO")
    print("=" * 100)
    print("  Repository footing (n4 via core): theta = 150 (5%), K_m = 100, k_e = 5, L = 322.73,")
    print("  25 monthly draws, dose train known, sigma = n4's, k_e prior 0.8, b prior 0.05 L.")
    print("  There is no patient data: every 'patient' below is simulated.")
    part1_0()
    R, FIMT2 = part1_ab()
    part1_global()
    prof = part1_profile()
    surr = part1_c()
    verdict1(R, FIMT2, prof, surr)
    FIM, MC, PL = part2()
    figure(R, FIM, MC)
    print(f"\n  total runtime {time.time() - T_START:.0f} s")


if __name__ == "__main__":
    main()
