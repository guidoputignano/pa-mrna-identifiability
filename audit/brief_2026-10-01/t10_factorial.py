"""
Follow-up item 3. TB-A against TB-B, and whether 5% of normal can reach 0.20.

  (a) The two benches of t6_estimators.bench differ in four ways at once:
        capacity    5% (theta = 150)              20% (theta = 600)
        sampling    every 1 h                     every 15 min
        assay CV    5%                            1%
        excursion   diet (half the load) out      load x6 for 24 h, then out
      Everything else as bench(): T = 48 h, L0 = 0.5 L, the same recorded load
      L = 322.73 at both capacities (convention A). The full 2^4 factorial of
      the information-matrix sd(ln theta_hat) (dyn.fim), with no priors as in
      t6, and again on the rung "k_e prior 0.8 + K_m from the record (0.062)".
      Main effects, interactions and Shapley shares in log10 sd. Then a Monte
      Carlo of t6_estimators.est_nls (t6's three starts) on the one-factor
      path cells: A, A with each single factor at its B value, B, and B with
      each single factor back at its A value.
  (b) At 5% (theta = 150, fixed): the cheapest design in a grid that reaches
      sd(ln theta_hat) <= 0.20, or none. Axes: sampling interval {5, 15, 30,
      60 min}, duration {24, 48, 72 h}, assay CV {5, 2, 1, 0.5%}, excursion
      (11 kinds), the knowledge rungs of t7_design.LADDER[:5], and the transient
      alone or joined to the monthly record (t7_design.joint_record).
  (c) The cheapest passing designs checked by Monte Carlo (est_nls with no
      priors; a MAP with the priors as pseudo-observations otherwise, and the
      monthly record's own likelihood when joined).

THE COST ORDERING, lexicographic, most important first:
  1  excursion safety class. Class 0 adds no propiogenic substrate: drug bursts
     at admissible peaks (n4's capacity ceiling at 5% allows a peak increment of
     2850) and diet withdrawals. Classes 1..6 are load increases x1.5 ... x8 in
     that order. The paper says decompensations "are triggered by anything that
     raises substrate load" (paper/sections/01_biology.tex:150-151), and
     audit/u15_gut_suppression.py withdraws the fast because catabolism raises
     the true load. The diet withdrawal here is u17's glucose-protected kind;
     that L0 stays constant under it is unverified (dyn.py).
  2  knowledge that must be bought. Free: no prior, the k_e prior 0.8 (the
     repository's), K_m from the record, and joining the record (the patient
     already has it). Bought: L0 to 5% of L (cost 1), L0 to 1% (cost 2, a
     research-grade measurement of the unrecorded load).
  3  assay CV better than the clinical 5%: 5% < 2% < 1% < 0.5%.
  4  the number of samples N = T/tau + 1.
  then, as tie-breaks only: shorter duration; the simpler analysis (fewer
  priors, alone before joined); within class 0, drug 450 < 1500 < 2850 <
  withdraw 25% < withdraw 50%.
The ordering is a judgement. The full frontier is printed so that any other
ordering can be read off it, and the cheapest design under a second ordering
(the clinical 5% assay first) is printed too.

ASSUMPTIONS beyond dyn.py's (L0 = 0.5 L, constant through the experiment):
  * the drug's delivered scale kappa is known (t7_design.drug's default); the
    cost of kappa at 12% and of kappa free is printed for the drug designs
  * joined designs keep the monthly record at its own clinical 5% assay CV
    (the record already exists) whatever the transient's assay
  * a better assay CV is treated as a property of the assay. Under the model
    (independent assay noise only) n replicates at 5% are worth one reading at
    5%/sqrt(n), so N (0.05/cv)^2 "5%-assay equivalents" is also printed. That
    equivalence fails if part of the 5% is pre-analytical or biological.

Run from this directory:  python3 -W error t10_factorial.py > out/t10_factorial.txt
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import itertools
import multiprocessing as mp
import time
from math import factorial

import numpy as np
from scipy.optimize import least_squares

import core as C
import dyn as D
import t6_estimators as T6
import t7_design as T7

T_START = time.time()
L = T7.L
GATE = C.GATE
KM_REC = T7.KM_REC
IDX = T7.IDX
RUNGS = T7.LADDER[:5]
RUNG_SHORT = ("none", "k_e 0.8", "+K_m rec", "+L0 5%", "+L0 1%")
KCOST = (0, 0, 0, 1, 2)
MODES = ("alone", "joined")
REPS_A, REPS_A_REV, REPS_C = 30, 20, 40


def banner(s):
    print("\n" + "=" * 100)
    print("  " + s)
    print("=" * 100)


def workers():
    try:
        la = os.getloadavg()[0]
    except OSError:
        la = 2.0
    return 3 if la < 1.6 else 2


# =========================================================== (a) factorial
FACTORS = ("capacity", "sampling", "assay CV", "excursion")
LEVELS = {"capacity": (C.TH_5, 0.20 * C.TH0),
          "sampling": (1.0, 0.25),                 # hours
          "assay CV": (0.05, 0.01),
          "excursion": ("diet out", "x6 24 h")}
LEV_TXT = {"capacity": ("5%", "20%"), "sampling": ("1 h", "15 min"),
           "assay CV": ("5%", "1%"), "excursion": ("diet out", "x6 24h")}
CELLS = list(itertools.product((0, 1), repeat=4))


def fact_cell(x, T=48):
    th = LEVELS["capacity"][x[0]]
    tau = LEVELS["sampling"][x[1]]
    cv = LEVELS["assay CV"][x[2]]
    if x[3] == 0:
        steps = [(-1, 0.5 * L), (0, 0.0)]
    else:
        steps = [(-1, 0.5 * L), (0, 5.5 * L), (T / 48, 0.0)]
    ex = D.Experiment(steps, L0=0.5 * L)
    t = np.arange(0, T / 24 + 1e-9, tau / 24)
    p = D.params(theta=th, ex=ex)
    return th, ex, t, cv, p


def effects(y):
    """Main effects, all interactions (B minus A in coded units) of a 2^4 table."""
    X = 2 * np.array(CELLS) - 1
    out = {}
    for k in range(1, 5):
        for S in itertools.combinations(range(4), k):
            c = np.prod(X[:, list(S)], axis=1)
            out[S] = float(np.mean(y[c > 0]) - np.mean(y[c < 0]))
    return out


def shapley(y):
    """Shapley shares of y(B) - y(A) over the four factors, v(S) = y(S at B) - y(A)."""
    pos = {c: i for i, c in enumerate(CELLS)}
    v = lambda S: y[pos[tuple(int(j in S) for j in range(4))]] - y[pos[(0, 0, 0, 0)]]
    n = 4
    phi = np.zeros(n)
    for i in range(n):
        rest = [j for j in range(n) if j != i]
        for k in range(n):
            for S in itertools.combinations(rest, k):
                w = factorial(k) * factorial(n - k - 1) / factorial(n)
                phi[i] += w * (v(set(S) | {i}) - v(set(S)))
    return phi


def nls_task(args):
    """One est_nls replicate on a factorial cell (t6's three starts)."""
    x, rep, seed = args
    th, ex, t, cv, p = fact_cell(x)
    m, _ = D.simulate(p, ex, t, sens=False)
    rng = np.random.default_rng(seed)
    y = m * np.exp(rng.normal(0, cv, m.size))
    L00 = p[3]
    ps = [p, p + np.array([0.7, -0.3, 0.5, -0.2 * L00]),
          p - np.array([0.7, -0.3, 0.5, -0.2 * L00])]
    t0 = time.time()
    try:
        xh = T6.est_nls(y, t, ex, ps)
    except Exception:          # every start failed
        xh = np.full(4, np.nan)
    return x, rep, np.asarray(xh, float), time.time() - t0


def mc_stats(lth_hat, lth):
    e = np.asarray(lth_hat, float) - lth
    ok = np.isfinite(e)
    e = e[ok]
    return dict(n=int(e.size), fail=int((~ok).sum()), bias=float(np.median(e)),
                rmse=float(np.sqrt(np.mean(e ** 2))),
                robust=float(0.5 * (np.percentile(e, 84) - np.percentile(e, 16))),
                off2=float(np.mean(np.abs(e) > np.log(2))),
                in20=float(np.mean(np.abs(e) <= GATE)), e=e)


def boot_rmse(e, seed, nb=2000):
    rng = np.random.default_rng(seed)
    b = np.sqrt(np.mean(rng.choice(e, (nb, e.size)) ** 2, axis=1))
    return float(np.percentile(b, 5)), float(np.percentile(b, 95))


def part_a():
    banner("(a) TB-A AGAINST TB-B: THE 2^4 FACTORIAL OF THE FIM sd(ln theta_hat)")
    # the cells reproduce t6's benches exactly
    for nm, x in (("A", (0, 0, 0, 0)), ("B", (1, 1, 1, 1))):
        b = T6.bench(nm)
        th, ex, t, cv, p = fact_cell(x)
        Fb = D.fim(b["p"], b["ex"], b["t"], cv=b["cv"])
        Fc = D.fim(p, ex, t, cv=cv)
        rel = float(np.abs(Fb - Fc).max() / np.abs(Fb).max())
        assert rel < 1e-12 and b["t"].size == t.size and b["cv"] == cv and b["th"] == th
        print(f"  cell {nm} reproduces t6_estimators.bench('{nm}'): max |dF|/max|F| = {rel:.1e}, "
              f"sd(ln theta) no priors {D.sd_lnth(Fc):.3f}")
    sdk = {}
    for th in (C.TH_5, 0.20 * C.TH0):
        sdk[th] = C.over_seeds(lambda s: float(np.sqrt(np.linalg.inv(
            C.fim(C.record(th, seed=s)))[1, 1])))
    print(f"  sd(ln K_m) the monthly record gives (k_e prior 0.8; median [min, max] over 24 seeds):"
          f"\n     5%: {sdk[C.TH_5][0]:.4f} [{sdk[C.TH_5][1]:.4f}, {sdk[C.TH_5][2]:.4f}];"
          f"  20%: {sdk[0.2*C.TH0][0]:.4f} [{sdk[0.2*C.TH0][1]:.4f}, {sdk[0.2*C.TH0][2]:.4f}]."
          f"  The rung uses t7's {KM_REC} at both capacities.")
    pri2 = {"ke": 0.8, "Km": KM_REC}
    Y = {"no priors": [], "k_e 0.8 + K_m 0.062": []}
    Y20 = []
    rows = []
    for x in CELLS:
        th, ex, t, cv, p = fact_cell(x)
        F = D.fim(p, ex, t, cv=cv)
        s0 = D.sd_lnth(F)
        s2 = D.sd_lnth(T7.with_priors(F, pri2))
        Y["no priors"].append(np.log10(s0))
        Y["k_e 0.8 + K_m 0.062"].append(np.log10(s2))
        if x[0] == 1:
            Y20.append(np.log10(D.sd_lnth(T7.with_priors(F, {"ke": 0.8, "Km": sdk[0.2 * C.TH0][0]}))) - np.log10(s2))
        rows.append((x, t.size, s0, s2))
    for k in Y:
        Y[k] = np.array(Y[k])
    print("\n  The 16 cells (0 = TB-A level, 1 = TB-B level)")
    print(f"     {'capacity':>9} {'sampling':>9} {'assay CV':>9} {'excursion':>10} {'N':>5} "
          f"{'sd, no priors':>14} {'sd, k_e+K_m':>12}")
    for x, n, s0, s2 in rows:
        tag = "   <- TB-A" if x == (0, 0, 0, 0) else ("   <- TB-B" if x == (1, 1, 1, 1) else "")
        print(f"     {LEV_TXT['capacity'][x[0]]:>9} {LEV_TXT['sampling'][x[1]]:>9} "
              f"{LEV_TXT['assay CV'][x[2]]:>9} {LEV_TXT['excursion'][x[3]]:>10} {n:5d} "
              f"{s0:14.4f} {s2:12.4f}{tag}")
    print(f"  capacity-matched K_m sd at 20% ({sdk[0.2*C.TH0][0]:.4f} instead of {KM_REC}) moves "
          f"log10 sd by at most {np.abs(Y20).max():.4f} over the 8 cells at 20%")

    print("\n  Effects in log10 sd(ln theta_hat): B level minus A level (coded -1/+1),")
    print("  so -1.0 means the B setting divides sd by 10. 10^effect is the factor on sd.")
    names = [f[:4] for f in ("capacity", "sampling", "assayCV", "excursion")]
    E = {k: effects(v) for k, v in Y.items()}
    hdr = f"     {'term':44}" + "".join(f"{k:>24}" for k in Y)
    print(hdr)
    for S in E["no priors"]:
        lab = " x ".join(FACTORS[i] for i in S)
        print(f"     {lab:44}" + "".join(f"{E[k][S]:+14.4f} (x{10**E[k][S]:7.3g})" for k in Y))
    # CV only scales the data's information when there are no priors
    cv_main = E["no priors"][(2,)]
    cv_int = max(abs(v) for S, v in E["no priors"].items() if 2 in S and len(S) > 1)
    assert abs(cv_main + np.log10(5)) < 1e-7 and cv_int < 1e-7
    print(f"     check: with no priors the assay CV only rescales the information, so its effect is "
          f"-log10(5) = {-np.log10(5):.4f}\n     (deviation {abs(cv_main + np.log10(5)):.1e}) and its "
          f"interactions vanish (max |.| {cv_int:.1e}), up to the inversion of\n     near-singular "
          f"matrices; asserted < 1e-7")
    for k, y in Y.items():
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        ss = {S: 4 * v ** 2 for S, v in E[k].items()}
        main = sum(v for S, v in ss.items() if len(S) == 1) / ss_tot
        two = sum(v for S, v in ss.items() if len(S) == 2) / ss_tot
        high = sum(v for S, v in ss.items() if len(S) > 2) / ss_tot
        assert abs(main + two + high - 1) < 1e-9
        dAB = y[CELLS.index((1, 1, 1, 1))] - y[CELLS.index((0, 0, 0, 0))]
        odd = sum(v for S, v in E[k].items() if len(S) in (1, 3))
        assert abs(dAB - odd) < 1e-9
        phi = shapley(y)
        assert abs(phi.sum() - dAB) < 1e-9
        print(f"\n  [{k}] variance of log10 sd over the 16 cells: main effects {main:.1%}, "
              f"two-way {two:.1%}, three- and four-way {high:.1%}")
        print(f"     log10 sd(B) - log10 sd(A) = {dAB:+.4f} (sd {10**y[0]:.3f} -> {10**y[-1]:.4f})"
              f" = sum of main effects {sum(E[k][(i,)] for i in range(4)):+.4f}"
              f" + three-way {sum(v for S, v in E[k].items() if len(S) == 3):+.4f}")
        print("     Shapley shares of that change (mean over the 24 orders of switching A -> B):")
        print("       " + ", ".join(f"{FACTORS[i]} {phi[i]:+.4f} ({phi[i]/dAB:.1%})" for i in range(4)))
        i0 = CELLS.index((0, 0, 0, 0))
        i1 = CELLS.index((1, 1, 1, 1))
        fromA = [y[CELLS.index(tuple(int(j == i) for j in range(4)))] - y[i0] for i in range(4)]
        fromB = [y[i1] - y[CELLS.index(tuple(int(j != i) for j in range(4)))] for i in range(4)]
        print("     one factor switched alone, from A (A+f minus A) and as the last switch into B"
              " (B minus B-f):")
        print("       " + ", ".join(f"{FACTORS[i]} {fromA[i]:+.3f} / {fromB[i]:+.3f}" for i in range(4)))
    return Y, E


def part_a_mc(Y):
    banner("(a) MONTE CARLO: est_nls (t6's three starts, no priors) ON THE ONE-FACTOR PATH CELLS")
    path = [("A", (0, 0, 0, 0))] + [(f"A, {FACTORS[i]} at B", tuple(int(j == i) for j in range(4)))
                                    for i in range(4)] + [("B", (1, 1, 1, 1))] + \
           [(f"B, {FACTORS[i]} at A", tuple(int(j != i) for j in range(4))) for i in range(4)]
    tasks = []
    for ci, (nm, x) in enumerate(path):
        nrep = REPS_A if ci < 6 else REPS_A_REV
        for r in range(nrep):
            tasks.append((x, r, 31_000 + 1000 * ci + r))
    #  slow cells first, for load balance
    order = sorted(range(len(tasks)), key=lambda k: -(4 * (tasks[k][0][1] == 1) + (tasks[k][0][0] == 0)))
    nw = workers()
    t0 = time.time()
    with mp.get_context("fork").Pool(nw) as pool:
        res = pool.map(nls_task, [tasks[k] for k in order], chunksize=1)
    print(f"  {len(tasks)} fits on {nw} worker processes in {time.time()-t0:.0f} s. Replicates: "
          f"{REPS_A} per requested cell (A, A+one factor, B) and {REPS_A_REV} per reverse cell"
          f"\n  (B with one factor back at A), against t6's 80, to hold the runtime;"
          f" RMSE from {REPS_A} replicates is uncertain by about +-{1/np.sqrt(2*REPS_A):.0%}.")
    byx = {}
    for x, r, xh, dt in res:
        byx.setdefault(x, []).append((r, xh, dt))
    print(f"\n     {'cell':30} {'reps':>4} {'FIM sd':>9} {'MC RMSE':>8} {'90% CI (bootstrap)':>19} "
          f"{'robust':>7} {'bias':>7} {'off >2x':>7} {'|err|<=.2':>9} {'s/fit':>6}")
    out = {}
    for ci, (nm, x) in enumerate(path):
        th, ex, t, cv, p = fact_cell(x)
        fsd = D.sd_lnth(D.fim(p, ex, t, cv=cv))
        rr = sorted(byx[x], key=lambda z: z[0])
        st = mc_stats([z[1][0] for z in rr], np.log(th))
        lo, hi = boot_rmse(st["e"], 77 + ci)
        out[x] = (fsd, st)
        print(f"     {nm:30} {st['n']:4d} {fsd:9.3f} {st['rmse']:8.3f}   [{lo:6.3f}, {hi:6.3f}]  "
              f"{st['robust']:7.3f} {st['bias']:+7.3f} {st['off2']:7.0%} {st['in20']:9.0%} "
              f"{np.mean([z[2] for z in rr]):6.1f}" + (f"  ({st['fail']} failed)" if st['fail'] else ""))
    xa, xb = (0, 0, 0, 0), (1, 1, 1, 1)
    print("\n  Attribution by single switches, log10 of the MC robust half-width (RMSE in brackets)"
          " against the FIM's log10 sd:")
    print(f"     {'factor':12} {'A+f vs A: MC':>22} {'FIM':>8} {'B vs B-f: MC':>22} {'FIM':>8}")
    for i in range(4):
        xf = tuple(int(j == i) for j in range(4))
        xr = tuple(int(j != i) for j in range(4))
        lr = lambda u, v: float(np.log10(u / v)) if u > 0 and v > 0 else np.nan
        a = lr(out[xf][1]["robust"], out[xa][1]["robust"])
        ar = lr(out[xf][1]["rmse"], out[xa][1]["rmse"])
        b = lr(out[xb][1]["robust"], out[xr][1]["robust"])
        br = lr(out[xb][1]["rmse"], out[xr][1]["rmse"])
        fa = np.log10(out[xf][0] / out[xa][0])
        fb = np.log10(out[xb][0] / out[xr][0])
        print(f"     {FACTORS[i]:12} {a:+10.3f} [{ar:+7.3f}] {fa:+8.3f} {b:+10.3f} [{br:+7.3f}] {fb:+8.3f}")
    return out


# =========================================================== (b) the 5% grid
EXC = (  # name, safety class, within-class rank, kind, argument
    ("drug burst, peak 450", 0, 0, "drug", 450),
    ("drug burst, peak 1500", 0, 1, "drug", 1500),
    ("drug burst, peak 2850", 0, 2, "drug", 2850),
    ("withdraw 25% of the load", 0, 3, "out", 0.25),
    ("withdraw 50% (the diet)", 0, 4, "out", 0.5),
    ("load x1.5, then out", 1, 0, "up", 1.5),
    ("load x2, then out", 2, 0, "up", 2),
    ("load x3, then out", 3, 0, "up", 3),
    ("load x4, then out", 4, 0, "up", 4),
    ("load x6, then out", 5, 0, "up", 6),
    ("load x8, then out", 6, 0, "up", 8))
TAUS = (5, 15, 30, 60)          # minutes
TS = (24, 48, 72)               # hours
CVS = (0.05, 0.02, 0.01, 0.005)


def grid_fim(kind, arg, T, tau_min, kfree=False):
    """Unit-CV information of one transient at 5% (t7_design.washout / drug)."""
    if kind == "drug":
        F, ex, p, t = T7.drug(arg, T=T, tau=tau_min / 60, kfree=kfree)
    elif kind == "out":
        F, ex, p, t = T7.washout(r_out=arg, T=T, tau=tau_min / 60)
    else:
        F, ex, p, t = T7.washout(up=arg, T=T, tau=tau_min / 60)
    assert t.size == round(T * 60 / tau_min) + 1
    return F * D.CV ** 2, ex, p, t


def joint_matrix(Fd, pri):
    """t7_design.joint_record, returning the matrix instead of its sd."""
    n = Fd.shape[0] + 1
    Z = np.zeros((n, n))
    idd = [0, 1, 2, 3] + ([5] if Fd.shape[0] > 4 else [])
    Z[np.ix_(idd, idd)] += Fd
    Z[np.ix_([0, 2, 4, 1], [0, 2, 4, 1])] += T7.FR
    pri = {k: v for k, v in pri.items() if k != "Km"}
    for k, sd in pri.items():
        j = 5 if k == "kap" else IDX[k]
        Z[j, j] += sd ** -2
    return Z


def nsamp(T, tau):
    return round(T * 60 / tau) + 1


def arank(r, j):
    return 5 * j + r


def part_b():
    banner("(b) AT 5% OF NORMAL: THE CHEAPEST DESIGN THAT REACHES sd(ln theta_hat) <= 0.20")
    shape = (len(EXC), len(TS), len(TAUS), len(CVS), len(RUNGS), 2)
    SD = np.full(shape, np.nan)
    for e, (nm, cls, rk, kind, arg) in enumerate(EXC):
        for a, T in enumerate(TS):
            for b, tau in enumerate(TAUS):
                F1, *_ = grid_fim(kind, arg, T, tau)
                for c, cv in enumerate(CVS):
                    F = F1 / cv ** 2
                    for r, (_, pri) in enumerate(RUNGS):
                        SD[e, a, b, c, r, 0] = D.sd_lnth(T7.with_priors(F, pri))
                        SD[e, a, b, c, r, 1] = T7.joint_record(F, pri)
    n_all = SD.size
    print(f"  {n_all} grid points: {len(EXC)} excursions x {len(TS)} durations x {len(TAUS)} intervals"
          f" x {len(CVS)} CVs x {len(RUNGS)} rungs x alone/joined")
    # ----- reproductions of t7's printed values
    eu4 = [nm for nm, *_ in EXC].index("load x4, then out")
    eu6 = [nm for nm, *_ in EXC].index("load x6, then out")
    ed = [nm for nm, *_ in EXC].index("drug burst, peak 2850")
    i48, i60, i5 = TS.index(48), TAUS.index(60), CVS.index(0.05)
    v4, v6 = SD[eu4, i48, i60, i5, 4, 0], SD[eu6, i48, i60, i5, 4, 0]
    vd, vdj = SD[ed, TS.index(72), TAUS.index(15), i5, 4, 0], SD[ed, TS.index(72), TAUS.index(15), i5, 4, 1]
    assert round(v4, 3) == 0.205 and round(v6, 3) == 0.176
    assert round(vd, 3) == 1.481 and round(vdj, 3) == 1.260
    print(f"  reproduces t7_design (hourly, 48 h, 5% CV, k_e 0.8 + K_m 0.062 + L0 to 1%): load x4 "
          f"{v4:.4f} (misses), x6 {v6:.4f} (passes);\n  drug 2850 every 15 min for 72 h, same rung: "
          f"alone {vd:.3f}, joined {vdj:.3f} (t7 printed 0.205, 0.176, 1.481, 1.260; asserted)")
    # the same x4 elsewhere in the grid
    print("  the x4 'miss' is specific to hourly over 48 h; on the same rung (L0 to 1%, 5% CV, alone):")
    for a, T in enumerate(TS):
        print(f"     T = {T:2d} h: " + ", ".join(f"{tau:2d} min {SD[eu4, a, b, i5, 4, 0]:.3f}"
                                         for b, tau in enumerate(TAUS)))

    # ----- the frontier: per (excursion, rung, mode, CV) the minimal N
    def best_cell(e, c, r, j):
        cand = []
        for a, T in enumerate(TS):
            for b, tau in enumerate(TAUS):
                v = SD[e, a, b, c, r, j]
                if v <= GATE:
                    cand.append((nsamp(T, tau), T, v, tau))
        if not cand:
            return None, float(np.nanmin(SD[e, :, :, c, r, j]))
        return min(cand), None

    print("\n  THE FRONTIER: for each (excursion, rung, alone/joined) and assay CV, the fewest samples")
    print("  N that reach 0.20 in the grid, as 'N tau/T sd', or 'none' with the best sd the grid gives.")
    print("  Joined: the record supplies K_m, so '+K_m rec' joined equals 'k_e 0.8' joined.")
    print(f"     {'excursion':26} {'rung':9} {'':6}" + "".join(f"{'CV ' + format(100*cv, 'g') + '%':>21}"
                                                       for cv in CVS))
    FR_ = {}
    for e, (nm, *_ ) in enumerate(EXC):
        for j in range(2):
            for r in range(len(RUNGS)):
                cells = []
                for c in range(len(CVS)):
                    bc, worst = best_cell(e, c, r, j)
                    FR_[(e, r, j, c)] = bc
                    if bc is None:
                        cells.append(f"none ({worst:.3g})")
                    else:
                        n, T, v, tau = bc
                        cells.append(f"{n:4d} {tau:2d}m/{T:2d}h {v:.3f}")
                print(f"     {nm[:26]:26} {RUNG_SHORT[r]:9} {MODES[j]:6}" + "".join(f"{s:>21}" for s in cells))
        print()

    # ----- compact: per (excursion, knowledge cost, CV), best analysis at that cost or below
    print("  COMPACT FRONTIER: fewest samples per (excursion, knowledge bought, CV); any analysis whose")
    print("  knowledge cost is at most the column's (free = none / k_e 0.8 / K_m from the record /")
    print("  joined to the record). Entry 'N tau/T sd [analysis]'.")
    print(f"     {'excursion':26} {'knowledge':10}" + "".join(f"{'CV ' + format(100*cv, 'g') + '%':>30}" for cv in CVS))
    KLAB = ("free", "L0 to 5%", "L0 to 1%")
    for e, (nm, *_ ) in enumerate(EXC):
        for kc in range(3):
            cells = []
            for c in range(len(CVS)):
                cand = []
                for r in range(len(RUNGS)):
                    if KCOST[r] > kc:
                        continue
                    for j in range(2):
                        bc = FR_[(e, r, j, c)]
                        if bc is not None:
                            n, T, v, tau = bc
                            cand.append((n, T, arank(r, j), v, tau, r, j))
                if not cand:
                    best = np.nanmin(SD[e, :, :, c, [r for r in range(5) if KCOST[r] <= kc], :])
                    cells.append(f"none ({best:.3g})")
                else:
                    n, T, _, v, tau, r, j = min(cand)
                    cells.append(f"{n:4d} {tau:2d}m/{T:2d}h {v:.3f} [{RUNG_SHORT[r][:7]},{MODES[j][0]}]")
            print(f"     {nm[:26]:26} {KLAB[kc]:10}" + "".join(f"{s:>30}" for s in cells))

    # ----- every passing design, ranked
    passing = []
    for e, (nm, cls, rk, kind, arg) in enumerate(EXC):
        for a, T in enumerate(TS):
            for b, tau in enumerate(TAUS):
                for c in range(len(CVS)):
                    for r in range(len(RUNGS)):
                        for j in range(2):
                            v = SD[e, a, b, c, r, j]
                            if v <= GATE:
                                passing.append(dict(e=e, cls=cls, rk=rk, T=T, tau=tau, c=c, r=r, j=j,
                                                    kc=KCOST[r], N=nsamp(T, tau), sd=float(v)))
    print(f"\n  {len(passing)} of {n_all} grid points pass. Passing grid points by safety class and assay CV:")
    print(f"     {'class':>5} {'excursions':36}" + "".join(f"{'CV ' + format(100*cv, 'g') + '%':>10}" for cv in CVS))
    for cls in range(7):
        names = ("drug 450/1500/2850, withdraw 25/50%" if cls == 0 else
                 [nm.split(",")[0] for nm, c_, *_ in EXC if c_ == cls][0])
        cnt = [sum(1 for d in passing if d["cls"] == cls and d["c"] == c) for c in range(len(CVS))]
        print(f"     {cls:5d} {names:36}" + "".join(f"{n:10d}" for n in cnt))

    # ----- the drug designs with the delivered scale kappa uncertain
    print("\n  The drug bursts with the delivered scale kappa known to 12% (the specification's Tier 1)"
          " or free:\n  fewest samples reaching 0.20 at the best analysis of each knowledge cost, or the best"
          " sd in the grid")
    print(f"     {'burst':12} {'kappa':6} {'knowledge':10}" + "".join(f"{'CV ' + format(100*cv, 'g') + '%':>22}" for cv in CVS))
    for e, (nm, cls, rk, kind, arg) in enumerate(EXC):
        if kind != "drug":
            continue
        Fk = {(T, tau): grid_fim(kind, arg, T, tau, kfree=True)[0] for T in TS for tau in TAUS}
        for klab, kx in (("12%", {"kap": 0.12}), ("free", {})):
            for kc in range(3):
                cells = []
                for c, cv in enumerate(CVS):
                    cand, bestv = [], np.inf
                    for (T, tau), F1 in Fk.items():
                        F = F1 / cv ** 2
                        for r, (_, pri) in enumerate(RUNGS):
                            if KCOST[r] > kc:
                                continue
                            for j in range(2):
                                v = (D.sd_lnth(T7.with_priors(F, {**pri, **kx})) if j == 0
                                     else T7.joint_record(F, {**pri, **kx}))
                                bestv = min(bestv, v)
                                if v <= GATE:
                                    cand.append((nsamp(T, tau), T, v, tau))
                    if cand:
                        n, T, v, tau = min(cand)
                        cells.append(f"{n:4d} {tau:2d}m/{T:2d}h {v:.3f}")
                    else:
                        cells.append(f"none ({bestv:.3g})")
                print(f"     {str(arg):12} {klab:6} {KLAB[kc]:10}" + "".join(f"{x:>22}" for x in cells))
    key1 = lambda d: (d["cls"], d["kc"], d["c"], d["N"], d["T"], arank(d["r"], d["j"]), d["rk"])
    key2 = lambda d: (d["c"], d["cls"], d["kc"], d["N"], d["T"], arank(d["r"], d["j"]), d["rk"])

    def desc(d):
        return (f"{EXC[d['e']][0]:26} every {d['tau']:2d} min for {d['T']:2d} h (N = {d['N']:3d}), "
                f"CV {100*CVS[d['c']]:g}%, {RUNG_SHORT[d['r']]} {MODES[d['j']]}: sd {d['sd']:.4f}"
                f"  [{d['N'] * (0.05 / CVS[d['c']]) ** 2:,.0f} 5%-assay equivalents]")

    for lab, key in (("ORDERING (i), the brief's: safety > knowledge bought > assay CV > N", key1),
                     ("ORDERING (ii), the clinical 5% assay first: assay CV > safety > knowledge > N", key2)):
        srt = sorted(passing, key=key)
        print(f"\n  {lab}")
        print("     the cheapest design for each excursion, in that order:")
        seen = []
        for d in srt:
            if d["e"] in seen:
                continue
            seen.append(d["e"])
            print(f"       {len(seen):2d}. " + desc(d))
        for e, (nm, *_ ) in enumerate(EXC):
            if e not in seen:
                print(f"           {nm:26} none in the grid")

    # ----- the multi-objective Pareto set over (class, knowledge cost, CV, N)
    vec = lambda d: (d["cls"], d["kc"], d["c"], d["N"])
    pts = {}
    for d in sorted(passing, key=key1):
        pts.setdefault(vec(d), d)
    V = list(pts)
    par = [v for v in V if not any(all(w[i] <= v[i] for i in range(4)) and w != v for w in V)]
    par.sort()
    print(f"\n  THE PARETO SET over (safety class, knowledge bought, assay CV, N): {len(par)} non-dominated"
          f" cost vectors out of {len(V)}")
    print(f"     {'class':>5} {'knowledge':>10} {'CV':>5} {'N':>4}   representative design")
    for v in par:
        d = pts[v]
        print(f"     {v[0]:5d} {KLAB[v[1]]:>10} {100*CVS[v[2]]:4g}% {v[3]:4d}   " + desc(d))
    return SD, passing, key1, key2, desc


# ========================================================== (c) Monte Carlo
LO = np.array([np.log(1e-2), np.log(1e-3), np.log(1.0), -1e4, -1e4])
HI = np.array([np.log(1e5), np.log(1e3), np.log(1e5), 1e5, 1e4])
REC = C.record(C.TH_5, seed=0)                    # the record T7.FR is built on
REC_SG = C.sigmas(REC)
REC_SB = C.LOAD_FRAC * REC["L"]


def rec_clear(p):
    th, ke, K = np.exp(p[0]), np.exp(p[1]), np.exp(p[2])
    mr, inc = REC["m"], REC["inc"]
    return (th + inc) * mr / (K + mr) + ke * mr


def rec_jac(p):
    th, ke, K = np.exp(p[0]), np.exp(p[1]), np.exp(p[2])
    mr, inc = REC["m"], REC["inc"]
    J = np.zeros((mr.size, 5))
    J[:, 0] = th * mr / (K + mr)
    J[:, 1] = ke * mr
    J[:, 2] = -(th + inc) * K * mr / (K + mr) ** 2
    J[:, 4] = 1.0
    return J / REC_SG[:, None]


def est_map(y, t, ex, cv, pri, cen, p0s, rec=None):
    """
    MAP through the ODE: ln y on ln m with sd cv, diagonal priors read as
    pseudo-observations (centres cen), and, if joined, the monthly record's
    flux ordinate with its load-bias prior. p = (ln th, ln ke, ln Km, L0[, b]).
    """
    ly = np.log(y)
    npar = 4 + (rec is not None)
    idx = sorted(pri)
    memo = {}

    def sim(p):
        k = p[:4].tobytes()
        if k not in memo:
            memo.clear()
            memo[k] = D.simulate(p[:4], ex, t, rtol=1e-7)
        return memo[k]

    def res(p):
        m, _ = sim(p)
        r = [(np.log(np.maximum(m, 1e-12)) - ly) / cv,
             np.array([(p[j] - cen[j]) / pri[j] for j in idx])]
        if rec is not None:
            r += [(rec_clear(p) + p[4] - rec["y"]) / REC_SG, np.array([(p[4] - rec["cb"]) / REC_SB])]
        return np.concatenate(r)

    def jac(p):
        m, S = sim(p)
        J = np.zeros((m.size, npar))
        J[:, :4] = S / m[:, None] / cv
        P = np.zeros((len(idx), npar))
        for k, j in enumerate(idx):
            P[k, j] = 1.0 / pri[j]
        blocks = [J, P]
        if rec is not None:
            B = np.zeros((1, npar))
            B[0, 4] = 1.0 / REC_SB
            blocks += [rec_jac(p), B]
        return np.vstack(blocks)

    lo, hi = LO[:npar], HI[:npar]
    best = None
    for p0 in p0s:
        try:
            s = least_squares(res, np.clip(p0, lo + 1e-6, hi - 1e-6), jac=jac, bounds=(lo, hi),
                              method="trf", x_scale="jac", xtol=1e-10, ftol=1e-10, max_nfev=200)
        except Exception:
            continue
        if best is None or s.cost < best.cost:
            best = s
    return (None if best is None else best.x), res, jac


def design_setup(d):
    nm, cls, rk, kind, arg = EXC[d["e"]]
    F1, ex, p, t = grid_fim(kind, arg, d["T"], d["tau"])
    cv = CVS[d["c"]]
    pri = {IDX[k]: v for k, v in RUNGS[d["r"]][1].items()}
    if d["j"] == 1:
        pri.pop(IDX["Km"], None)               # the record supplies K_m itself
    return F1, ex, p, t, cv, pri


def map_task(args):
    d, rep, seed = args
    F1, ex, p, t, cv, pri = design_setup(d)
    m, _ = D.simulate(p, ex, t, sens=False)
    rng = np.random.default_rng(seed)
    y = m * np.exp(rng.normal(0, cv, m.size))
    cen = {j: p[j] + sd * rng.normal() for j, sd in pri.items()}
    rec = None
    if d["j"] == 1:
        rec = dict(y=rec_clear(p) + rng.normal(0, REC_SG), cb=REC_SB * rng.normal())
    L00 = p[3]
    dp = np.array([0.7, -0.3, 0.5, -0.2 * L00])
    t0 = time.time()
    if not pri and rec is None:
        try:
            xh = T6.est_nls(y, t, ex, [p, p + dp, p - dp])
        except Exception:
            xh = np.full(4, np.nan)
        return d["id"], rep, np.asarray(xh, float), time.time() - t0
    base = [p, p + dp, p - dp]
    pc = p.copy()
    for j, c in cen.items():
        pc[j] = c
    base.append(pc)
    if rec is not None:
        base = [np.r_[q, 0.0] for q in base]
    xh, _, _ = est_map(y, t, ex, cv, pri, cen, base, rec)
    if xh is None:
        xh = np.full(4 + (rec is not None), np.nan)
    return d["id"], rep, np.asarray(xh, float), time.time() - t0


def part_c(passing, key1, key2, desc):
    banner("(c) MONTE CARLO OF THE CHEAPEST PASSING DESIGNS")
    sel = []
    seen = set()
    for d in sorted(passing, key=key1):
        if d["e"] not in seen:
            seen.add(d["e"])
            sel.append(dict(d, why="ordering (i), cheapest for its excursion"))
        if len(sel) == 3:
            break
    same = lambda u, v: all(u[k] == v[k] for k in ("e", "T", "tau", "c", "r", "j"))
    d2 = min(passing, key=key2)
    if not any(same(d2, s) for s in sel):
        sel.append(dict(d2, why="ordering (ii): the cheapest at the clinical 5% assay"))
    d3 = min([d for d in passing if d["kc"] == 0], key=key2)
    if not any(same(d3, s) for s in sel):
        sel.append(dict(d3, why="ordering (ii) with free knowledge only: the cheapest at the clinical"
                                " 5% assay that buys nothing"))
    for i, d in enumerate(sel):
        d["id"] = i
    print("  The designs checked, chosen by the orderings of part (b):")
    for d in sel:
        print(f"     [{d['id']}] " + desc(d) + f"\n         ({d['why']})")
    # the estimator's information at the truth equals the FIM the grid used
    print("  check: J'J of each estimator's residuals at the truth (noise-free, centres at truth)"
          " against the grid's FIM")
    for d in sel:
        F1, ex, p, t, cv, pri = design_setup(d)
        Fd = F1 / cv ** 2
        if d["j"] == 0:
            Ft = T7.with_priors(Fd, RUNGS[d["r"]][1])
            pp = p
            rec = None
        else:
            Ft = joint_matrix(Fd, RUNGS[d["r"]][1])
            assert abs(D.sd_lnth(Ft) - T7.joint_record(Fd, RUNGS[d["r"]][1])) < 1e-12
            pp = np.r_[p, 0.0]
            rec = dict(y=rec_clear(p), cb=0.0)
        m, _ = D.simulate(p, ex, t, sens=False)
        _, res, jac = est_map(m, t, ex, cv, pri, {j: p[j] for j in pri}, [], rec)
        J = jac(pp)
        rel = float(np.abs(J.T @ J - Ft).max() / np.abs(Ft).max())
        r0 = float(np.abs(res(pp)).max())
        assert rel < 1e-6 and r0 < 1e-3
        print(f"     [{d['id']}] max |J'J - F| / max|F| = {rel:.1e}; largest residual at truth {r0:.1e} noise sd"
              f" (ODE rtol 1e-7 against the FIM's 1e-9);"
              f"\n         FIM sd {D.sd_lnth(Ft):.4f} (grid {d['sd']:.4f}); asserted < 1e-6 and < 1e-3")
    for d in sel:
        nm, cls, rk, kind, arg = EXC[d["e"]]
        if kind != "drug":
            continue
        Fk, *_ = grid_fim(kind, arg, d["T"], d["tau"], kfree=True)
        Fk = Fk / CVS[d["c"]] ** 2
        pri = RUNGS[d["r"]][1]
        f = (lambda F, q: D.sd_lnth(T7.with_priors(F, q))) if d["j"] == 0 else \
            (lambda F, q: T7.joint_record(F, q))
        print(f"     [{d['id']}] the drug's delivered scale kappa: known {d['sd']:.4f}; to 12% "
              f"{f(Fk, {**pri, 'kap': 0.12}):.4f}; free {f(Fk, pri):.4f}")
    tasks = [(d, r, 52_000 + 1000 * d["id"] + r) for d in sel for r in range(REPS_C)]
    tasks.sort(key=lambda z: -z[0]["N"])
    nw = workers()
    t0 = time.time()
    with mp.get_context("fork").Pool(nw) as pool:
        res = pool.map(map_task, tasks, chunksize=1)
    print(f"\n  {len(tasks)} fits ({REPS_C} replicates per design) on {nw} worker processes in "
          f"{time.time()-t0:.0f} s. Estimator: est_nls (t6, three starts) where the analysis has no"
          f"\n  priors; otherwise est_map, with t6's three starts plus one at the prior centres. Each"
          f" replicate draws the\n  transient's noise, the prior centres (truth + sd N(0,1)) and,"
          f" when joined, the record's noise and its load-bias centre.")
    print(f"\n     {'design':6} {'FIM sd':>8} {'MC RMSE':>8} {'90% CI (bootstrap)':>19} {'robust':>7} "
          f"{'bias':>7} {'off >2x':>7} {'|err|<=.2':>9} {'MC/FIM':>7} {'s/fit':>6}")
    out = []
    for d in sel:
        rr = sorted([z for z in res if z[0] == d["id"]], key=lambda z: z[1])
        st = mc_stats([z[2][0] for z in rr], np.log(C.TH_5))
        lo, hi = boot_rmse(st["e"], 991 + d["id"])
        ke = np.array([z[2][1] for z in rr]) - np.log(C.KE)
        out.append((d, st, lo, hi))
        print(f"     [{d['id']}]    {d['sd']:8.4f} {st['rmse']:8.4f}   [{lo:6.3f}, {hi:6.3f}]  "
              f"{st['robust']:7.4f} {st['bias']:+7.3f} {st['off2']:7.0%} {st['in20']:9.0%} "
              f"{st['rmse']/d['sd']:7.2f} {np.mean([z[3] for z in rr]):6.1f}"
              + (f"  ({st['fail']} failed)" if st['fail'] else "")
              + f"   ln k_e RMSE {np.sqrt(np.nanmean(ke ** 2)):.3f}")
        print(f"            target 0.20: MC RMSE {'<=' if st['rmse'] <= GATE else '>'} 0.20, bootstrap"
              f" 90% interval {'entirely below' if hi <= GATE else ('entirely above' if lo > GATE else 'straddles')} 0.20")
    return out


def main():
    print("=" * 100)
    print("  FOLLOW-UP ITEM 3. TB-A VERSUS TB-B, AND WHETHER 5% OF NORMAL CAN REACH 0.20")
    print("=" * 100)
    print(f"  footing: theta = {C.TH_5:g} (5%) and {0.2*C.TH0:g} (20%), K_m = {C.KM0:g}, k_e = {C.KE:g} /d,"
          f" L = {L:.2f} at both (convention A),\n  L0 = 0.5 L unrecorded and constant, target sd(ln theta_hat)"
          f" <= {GATE}. Time in days inside the model; hours and minutes in the labels.")
    Y, E = part_a()
    part_a_mc(Y)
    SD, passing, key1, key2, desc = part_b()
    part_c(passing, key1, key2, desc)
    print(f"\n  total runtime {time.time() - T_START:.0f} s")


if __name__ == "__main__":
    main()
