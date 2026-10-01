"""
Follow-up item 1. The 0.39 "floor": sd(ln theta_hat) with k_e known exactly.

  (a) Var(ln theta_hat | k_e exact) split over four factors: the load bias b
      (prior 0.05 L), K_m (no prior), and the two parts of the per-visit sigma,
      the 5% load-record scatter and the 5% assay CV carried through C'(m_bar).
      Shapley values over the four factors, the Mobius (Harsanyi) dividends that
      show where the interactions sit, and the order-dependent figures (each
      factor removed first, removed last, and along two explicit orders).
      At 5%, 10%, 20%, 50% and 100% of normal.
  (b) n-dependence at 5% with k_e exact: longer monthly records and denser
      sampling within two years, the fit sd(n) = sqrt(floor^2 + c^2/n), and the
      same with b exact, with the load measured exactly, and with K_m exact.
  (c) the minimal combination reaching 0.20 at 5%: k_e precision x load-bias
      precision x load scatter x draws.
  (d) the k_e variance share 1 - Var(ln th | k_e)/Var(ln th) against capacity.

Footing: core.py's docstring (n4's, convention A). "Removing" a factor means:
  b        b known exactly: its coordinate deleted from the information matrix
  K_m      K_m known exactly: its coordinate deleted
  scatter  the load-record term of sigma set to zero; the b prior is passed
           separately and stays at 0.05 L (core.fim's load_frac otherwise
           moves both)
  assay    the assay CV set to zero in sigma
k_e is removed throughout (its coordinate deleted), which is t1's "|ke".
With both sigma terms at zero the record is exact and Var = 0: that is not
computed through a singular matrix but asserted from the rank of the
remaining sensitivity columns, and checked against a vanishing sigma.

Every number is a median over the 24 intake records (core.SEEDS) unless it
says otherwise.
"""
import itertools
from math import factorial

import numpy as np
from scipy.optimize import curve_fit

import core as C
from design.run import delta, D_AMP, LAMBDA, DOSE_INTERVAL
from plotstyle import plt, C1, C2, C3, C4, INK, MUTED, finish, save, plain

TH5 = C.TH_5
EDGE = 2472.0                      # n4's lower band edge (t3's EDGE)
FACTORS = ("b", "K_m", "scatter", "assay")
NF = len(FACTORS)
SUBSETS = [frozenset(s) for k in range(NF + 1)
           for s in itertools.combinations(FACTORS, k)]
GATE = C.GATE
LAST_DOSE = 199 * DOSE_INTERVAL    # design.run.delta sums doses k = 0..199


# ------------------------------------------------------------- the record
def delta_cont(t, D=D_AMP):
    """
    The booked dose train continued indefinitely, in closed form. A variant of
    design.run.delta, which stops after 200 doses (day 2786): every dose k with
    k*14 <= t, summed as a geometric series.
    """
    t = np.asarray(t, float)
    n = np.floor(t / DOSE_INTERVAL)
    n = np.where((n + 1) * DOSE_INTERVAL <= t, n + 1, n)
    n = np.where(n * DOSE_INTERVAL > t, n - 1, n)
    q = np.exp(-LAMBDA * DOSE_INTERVAL)
    return D * np.exp(-LAMBDA * (t - n * DOSE_INTERVAL)) * (1 - q ** (n + 1)) / (1 - q)


def rec(theta, times, seed, cont=False):
    """
    core.record, with the increment supplied here so that it is computed once
    per schedule rather than once per record. cont=True continues the dose
    train past day 2786; otherwise the repository's train is used (identical
    to delta_cont before day 2786, n4's truncated GRID_INC slice after it).
    """
    times = np.asarray(times, float)
    r = C.record(theta, D=0, seed=seed, times=times)
    if cont or times.max() < LAST_DOSE:
        r["inc"] = delta_cont(times)
    else:
        r["inc"] = C.N4.delta_at(times, C.N4.D0)
    r["m"] = C.N4.level(theta, r["Km"], r["ke"], r["loads"], r["inc"])
    return r


# ------------------------------------------------------------ the variance
def var_lnth(r, ke_sd=0.0, km_exact=False, b_frac=C.LOAD_FRAC,
             scatter=C.LOAD_FRAC, assay=C.ASSAY_CV):
    """
    Var(ln theta_hat) on core.fim's information matrix. ke_sd = 0 and
    b_frac = 0 mean exact (coordinate deleted); ke_sd = None means no prior.
    b_frac is the bias prior in units of L; scatter and assay are sigma's
    two terms, set independently.
    """
    fixed = []
    if ke_sd == 0:
        fixed.append(C.I_KE)
    if km_exact:
        fixed.append(C.I_KM)
    if b_frac == 0:
        fixed.append(C.I_B)
    keep = [i for i in range(4) if i not in fixed]
    if scatter == 0 and assay == 0:
        Q = C.flux_rows(r)[:, keep]
        assert np.linalg.matrix_rank(Q) == len(keep), "noise-free record not full rank"
        return 0.0
    F = C.fim(r, ke_prior=(ke_sd if ke_sd else None),
              b_prior=(b_frac * r["L"] if b_frac else np.inf),
              load_frac=scatter, assay_cv=assay)
    P = C.cov(F[np.ix_(keep, keep)])
    return np.inf if P is None else float(P[0, 0])


def removed_kw(S):
    """var_lnth keywords for the factor set S removed (k_e exact throughout)."""
    return dict(ke_sd=0.0, km_exact="K_m" in S,
                b_frac=0.0 if "b" in S else C.LOAD_FRAC,
                scatter=0.0 if "scatter" in S else C.LOAD_FRAC,
                assay=0.0 if "assay" in S else C.ASSAY_CV)


def med(v):
    return float(np.median(np.asarray(v, float)))


# ------------------------------------------------------------- the games
def shapley(V):
    """
    Shapley values of the elimination game: phi_i is factor i's reduction of
    Var, V(S) - V(S + i), averaged over all orders of removal. They sum to
    V(none) - V(all).
    """
    phi = {}
    for f in FACTORS:
        tot = 0.0
        for S in SUBSETS:
            if f in S:
                continue
            w = factorial(len(S)) * factorial(NF - len(S) - 1) / factorial(NF)
            tot += w * (V[S] - V[S | {f}])
        phi[f] = tot
    return phi


def dividends(V):
    """
    Mobius dividends of the switch-on game w(T) = V(all factors outside T
    removed), i.e. the variance with only the factors in T present. Main
    effects are the singletons; the rest are interactions. They sum to V(none).
    """
    full = frozenset(FACTORS)
    w = {T: V[full - T] for T in SUBSETS}
    d = {}
    for T in SUBSETS:
        d[T] = sum((-1) ** (len(T) - len(S)) * w[S] for S in SUBSETS if S <= T)
    return d


def order_figs(V):
    """Each factor's reduction along every one of the 4! removal orders."""
    out = {f: [] for f in FACTORS}
    for order in itertools.permutations(FACTORS):
        S = frozenset()
        for f in order:
            out[f].append(V[S] - V[S | {f}])
            S = S | {f}
    return out


def lab(S):
    return "{" + ", ".join(f for f in FACTORS if f in S) + "}" if S else "{}"


# ------------------------------------------------------------------ (a)
def part_a():
    print("\n" + "-" * 92)
    print("  (a) DECOMPOSITION OF Var(ln theta_hat | k_e exact), monthly 25 draws")
    print("-" * 92)
    # checks first: the conditional variance is the inverse of the sub-matrix,
    # and var_lnth on the defaults is core.fim's and n4's number
    worst_cond, worst_n4 = 0.0, 0.0
    for s in C.SEEDS:
        r = C.record(TH5, seed=s)
        F = C.fim(r)
        P = np.linalg.inv(F)
        a = P[0, 0] - P[0, 3] ** 2 / P[3, 3]
        worst_cond = max(worst_cond, abs(var_lnth(r) / a - 1))
        worst_n4 = max(worst_n4, abs(np.sqrt(var_lnth(r, ke_sd=C.KE_PRIOR))
                                     / C.N4.sd_ln_theta(TH5, seed=s) - 1))
    assert worst_cond < 1e-9 and worst_n4 < 1e-9
    print(f"  checks: inv(F[keep,keep]) = P00 - P03^2/P33 to {worst_cond:.1e} (relative); "
          f"with the k_e prior, n4.sd_ln_theta to {worst_n4:.1e}")
    r0 = C.record(TH5, seed=0)
    lim = [np.sqrt(var_lnth(r0, scatter=e * C.LOAD_FRAC, assay=e * C.ASSAY_CV))
           for e in (1e-2, 1e-4)]
    print(f"  vanishing sigma (seed 0, both terms scaled by 1e-2, 1e-4): sd = {lim[0]:.3e}, "
          f"{lim[1]:.3e}; taken as 0 at zero noise")
    r = C.record(TH5, seed=0)
    sl = TH5 * C.KM0 / (C.KM0 + r["m_bar"]) ** 2 + C.KE
    print(f"  sigma's two terms at 5%: load scatter 0.05 L = {C.LOAD_FRAC * r['L']:.3f}, "
          f"assay C'(m_bar) 0.05 m_bar = {sl * C.ASSAY_CV * r['m_bar']:.3f}, "
          f"total {C.sigmas(r)[0]:.3f} (flux units)")

    results = {}
    for pct in (5, 10, 20, 50, 100):
        th = C.TH0 * pct / 100
        per = {S: [] for S in SUBSETS}
        phis = []
        for s in C.SEEDS:
            r = C.record(th, seed=s)
            Vs = {S: var_lnth(r, **removed_kw(S)) for S in SUBSETS}
            for S in SUBSETS:
                per[S].append(Vs[S])
            ph = shapley(Vs)
            assert abs(sum(ph.values()) - (Vs[frozenset()] - Vs[frozenset(FACTORS)])) < 1e-12
            phis.append([ph[f] / Vs[frozenset()] for f in FACTORS])
        V = {S: med(per[S]) for S in SUBSETS}
        results[pct] = dict(V=V, phis=np.array(phis), sd_med=med(np.sqrt(per[frozenset()])))
    pcts = list(results)

    print("\n  The 16 subsets: sd(ln theta_hat) with k_e exact and the listed factors removed")
    print("  (sqrt of the median variance over 24 records)")
    print(f"    {'removed':32}" + "".join(f"{p:>9d}%" for p in pcts))
    for S in SUBSETS:
        print(f"    {lab(S):32}" + "".join(f"{np.sqrt(results[p]['V'][S]):10.4f}" for p in pcts))
    print(f"    {'median sd, nothing removed':32}" + "".join(f"{results[p]['sd_med']:10.4f}" for p in pcts))

    print("\n  Shapley shares of Var(ln theta | k_e) (game on the median variances; sums to 100%)")
    print("  and, in parentheses, the median [min-max] of each record's own Shapley share")
    print(f"    {'pct':>4} {'Var':>9}  " + "".join(f"{f:>28}" for f in FACTORS))
    shap = {}
    for p in pcts:
        V = results[p]["V"]
        ph = shapley(V)
        tot = V[frozenset()]
        assert abs(sum(ph.values()) - tot) < 1e-12
        shap[p] = {f: ph[f] / tot for f in FACTORS}
        ps = results[p]["phis"]
        cells = [f"{ph[f]/tot:6.1%} ({np.median(ps[:, j]):5.1%} [{ps[:, j].min():5.1%}-{ps[:, j].max():5.1%}])"
                 for j, f in enumerate(FACTORS)]
        print(f"    {p:4d} {tot:9.5f}  " + "".join(f"{c:>28}" for c in cells))
    print("    the two sigma terms together: " + ", ".join(
        f"{p}%: {shap[p]['scatter'] + shap[p]['assay']:.1%}" for p in pcts))

    print("\n  Order-dependent figures: each factor's reduction of Var as a share of Var(none)")
    print("    removed FIRST = V({}) - V({i});  removed LAST = V(all but i) - V(all);")
    print("    [min, max] over all 24 removal orders; Shapley = the mean over the 24 orders")
    for p in pcts:
        V = results[p]["V"]
        tot = V[frozenset()]
        of = order_figs(V)
        print(f"    {p}% of normal (Var {tot:.5f}, sd {np.sqrt(tot):.4f})")
        print(f"      {'factor':9} {'first':>8} {'last':>8} {'min':>8} {'max':>8} {'Shapley':>8}")
        for f in FACTORS:
            first = V[frozenset()] - V[frozenset({f})]
            rest = frozenset(FACTORS) - {f}
            last = V[rest] - V[frozenset(FACTORS)]
            print(f"      {f:9} {first/tot:8.1%} {last/tot:8.1%} {min(of[f])/tot:8.1%} "
                  f"{max(of[f])/tot:8.1%} {np.mean(of[f])/tot:8.1%}")
        for order in (("b", "K_m", "scatter", "assay"), ("assay", "scatter", "K_m", "b"),
                      ("K_m", "b", "assay", "scatter")):
            S, steps = frozenset(), []
            for f in order:
                steps.append(f"{f} {(V[S] - V[S | {f}])/tot:.1%}")
                S = S | {f}
            print("      order " + " -> ".join(steps))

    print("\n  Mobius dividends of the switch-on game (variance with ONLY the factors in T present),")
    print("  as shares of Var(none). Singletons are main effects, the rest interactions; they sum")
    print("  to 100%. b and K_m have zero main effect: without per-visit noise there is no variance,")
    print("  whatever b's prior and however unknown K_m is.")
    print(f"    {'T (present)':32}" + "".join(f"{p:>9d}%" for p in pcts))
    D = {p: dividends(results[p]["V"]) for p in pcts}
    for p in pcts:
        assert abs(sum(D[p].values()) - results[p]["V"][frozenset()]) < 1e-12
    for T in SUBSETS[1:]:
        print(f"    {lab(T):32}" + "".join(
            f"{D[p][T]/results[p]['V'][frozenset()]:10.1%}" for p in pcts))
    return results, shap


# ------------------------------------------------------------------ (b)
SCEN = (("k_e exact", dict()),
        ("k_e + b exact (scatter kept)", dict(b_frac=0.0)),
        ("k_e + load measured exactly", dict(b_frac=0.0, scatter=0.0)),
        ("k_e + K_m exact", dict(km_exact=True)),
        ("k_e + K_m + b exact", dict(km_exact=True, b_frac=0.0)))
NS = np.array([6, 8, 12, 16, 25, 36, 50, 75, 100, 150, 200, 300, 400, 600])
NS_LONG = np.array([1200, 2400, 6000])
DTS = (30.4, 14.5, 7.3, 2.0, 1.0)
NSCAN = 800


def floor_fit(n, sd):
    f = lambda n, fl, c: np.sqrt(fl ** 2 + c ** 2 / n)
    (fl, c), _ = curve_fit(f, n, sd, p0=(0.1, 1.0), bounds=([0, 0], [np.inf, np.inf]))
    return fl, c, float(np.abs(sd - f(n, fl, c)).max()), f


def sd_design(times, kw, cont):
    return med([np.sqrt(var_lnth(rec(TH5, times, s, cont), **kw)) for s in C.SEEDS])


def b_post_ratio(times, cont):
    """Posterior sd(b) / prior sd(b), k_e exact, everything else at the footing."""
    out = []
    for s in C.SEEDS:
        r = rec(TH5, times, s, cont)
        F = C.fim(r, ke_prior=None)
        P = np.linalg.inv(F[np.ix_([0, 1, 2], [0, 1, 2])])
        out.append(np.sqrt(P[2, 2]) / (C.LOAD_FRAC * r["L"]))
    return med(out)


def part_b():
    print("\n" + "-" * 92)
    print("  (b) n-DEPENDENCE AT 5% OF NORMAL WITH k_e EXACT")
    print("-" * 92)
    t = np.arange(0.0, LAST_DOSE, 0.37)
    err = float(np.max(np.abs(delta_cont(t) / delta(t) - 1)))
    assert err < 1e-10
    print(f"  delta_cont (the dose train continued) equals design.run.delta to {err:.1e} (relative) "
          f"before day {LAST_DOSE:.0f}.")
    n_tr = int(np.sum(C.GRID <= LAST_DOSE))
    print(f"  NOTE. design.run.delta stops after 200 doses, the last on day {LAST_DOSE:.0f}. On n4's GRID")
    print(f"  that is monthly draw number {n_tr} (t = {C.GRID[n_tr-1]:.1f} d); draw {n_tr+1} sits "
          f"{C.GRID[n_tr]-LAST_DOSE:.1f} d after it (increment {C.N4.GRID_INC[n_tr]:.1f}), and from")
    print(f"  draw {n_tr+2} on the increment is below {C.N4.GRID_INC[n_tr+1]:.2f}, i.e. off therapy. So "
          f"prefixes of core.GRID longer than ~{n_tr+1} draws")
    print("  add UNTREATED draws. Both versions are shown: 'n4 GRID' as-is, and 'therapy continued'.")
    worst = 0.0
    for times in (C.MONTHLY, C.DENSE, C.GRID, np.arange(0.0, 730.0, 7.3)):
        for s in (0, 7, 23):
            a, b = rec(TH5, times, s), C.record(TH5, seed=s, times=times)
            worst = max(worst, float(np.max(np.abs(a["m"] / b["m"] - 1))),
                        float(np.max(np.abs(a["inc"] - b["inc"]) / (1 + b["inc"]))))
    assert worst < 1e-9
    print(f"  rec() reproduces core.record (monthly, every 2 d, the 600-draw GRID, every 7.3 d) "
          f"to {worst:.1e}")
    rA, rB = C.record(TH5, seed=3, times=C.GRID[:50]), C.record(TH5, seed=3, times=C.GRID[:25])
    print(f"  prefix property (the 25-draw record is the first 25 draws of the 50-draw one): "
          f"{bool(np.allclose(rA['loads'][:25], rB['loads']))}")

    tabs = {}
    for cont, tag in ((False, "n4 GRID as-is (untreated after draw ~93)"),
                      (True, "therapy continued")):
        ns = NS if not cont else np.r_[NS, NS_LONG]
        print(f"\n  (i) longer monthly records, {tag}: sd(ln theta_hat)")
        print(f"    {'n':>5}" + "".join(f"{nm:>30}" for nm, _ in SCEN))
        rows = np.array([[sd_design(np.arange(n) * 30.4, kw, cont) for _, kw in SCEN] for n in ns])
        for n, row in zip(ns, rows):
            print(f"    {n:5d}" + "".join(f"{v:30.4f}" for v in row))
        print(f"    {'n*Var':>5} at n = 25, 100, 600" + (", 6000:" if cont else ":"))
        for j, (nm, _) in enumerate(SCEN):
            sel = [25, 100, 600] + ([6000] if cont else [])
            vals = [rows[list(ns).index(k), j] ** 2 * k for k in sel]
            print(f"      {nm:30} " + "  ".join(f"{v:7.3f}" for v in vals))
        tabs[cont] = (ns, rows)

    # local slopes on the continued record
    ns, rows = tabs[True]
    print("\n  local slope d ln sd / d ln n, therapy continued (-0.5 is 1/sqrt(n); 0 is a floor)")
    pairs = [(25, 50), (50, 100), (100, 200), (200, 400), (400, 600), (600, 1200), (1200, 2400), (2400, 6000)]
    print(f"    {'n1-n2':>10}" + "".join(f"{nm:>30}" for nm, _ in SCEN))
    il = list(ns)
    for a, b in pairs:
        sl = [np.log(rows[il.index(b), j] / rows[il.index(a), j]) / np.log(b / a) for j in range(len(SCEN))]
        print(f"    {f'{a}-{b}':>10}" + "".join(f"{v:30.3f}" for v in sl))

    print("\n  what the data teach about b (k_e exact): posterior sd(b) / prior sd(b) = 0.05 L")
    for n in (6, 25, 100, 365, 600, 2400, 6000):
        times = np.arange(0.0, 730.0, 2.0) if n == 365 else np.arange(n) * 30.4
        lab_ = "365 every 2 d" if n == 365 else f"{n} monthly, therapy continued"
        print(f"    {lab_:32} {b_post_ratio(times, True):.4f}")

    # dense within two years
    print("\n  (ii) denser sampling within two years (therapy on throughout): sd(ln theta_hat)")
    print(f"    {'dt (d)':>7} {'n':>5}" + "".join(f"{nm:>30}" for nm, _ in SCEN))
    dn, drows = [], []
    for dt in DTS:
        tt = np.arange(0.0, 730.0, dt)
        row = [sd_design(tt, kw, False) for _, kw in SCEN]
        dn.append(tt.size)
        drows.append(row)
        print(f"    {dt:7.1f} {tt.size:5d}" + "".join(f"{v:30.4f}" for v in row))
    dn, drows = np.array(dn), np.array(drows)

    print("\n  fits of sd(n) = sqrt(floor^2 + c^2/n) (least squares on sd, floor >= 0)")
    print(f"    {'scenario':30} {'design':34} {'floor':>8} {'c':>8} {'worst res':>10}")
    fits = {}
    for j, (nm, _) in enumerate(SCEN):
        for key, (n, sd), dlab in (
                ("grid", (tabs[False][0], tabs[False][1][:, j]), "monthly n4 GRID, n 6-600"),
                ("cont", (NS, tabs[True][1][:len(NS), j]), "monthly therapy cont., n 6-600"),
                ("dense", (dn, drows[:, j]), "2 years, dt 30.4-1 d, n 25-730")):
            fl, c, wr, f = floor_fit(n.astype(float), sd)
            fits[(j, key)] = (fl, c, f)
            print(f"    {nm:30} {dlab:34} {fl:8.4f} {c:8.4f} {wr:10.4f}")
    print("\n  out of sample: the therapy-continued fit over n 6-600, extrapolated, against the FIM")
    print(f"    {'scenario':30} " + "".join(f"{'n='+str(n)+' fit/FIM':>22}" for n in NS_LONG))
    for j, (nm, _) in enumerate(SCEN):
        fl, c, f = fits[(j, "cont")]
        cells = [f"{f(n, fl, c):.4f}/{tabs[True][1][len(NS) + k, j]:.4f}" for k, n in enumerate(NS_LONG)]
        print(f"    {nm:30} " + "".join(f"{x:>22}" for x in cells))

    print("\n  asymptotic coefficient c_inf = sqrt(n Var) at n = 6000, therapy continued:")
    for j, (nm, _) in enumerate(SCEN):
        print(f"    {nm:30} c_inf = {np.sqrt(6000) * tabs[True][1][-1, j]:.4f}")

    # draws needed to reach the gate
    print("\n  draws needed to reach sd(ln theta_hat) <= 0.20 (median over 24 records)")
    print("  monthly, therapy continued: records are nested prefixes, so each record's variance")
    print("  falls monotonically with n and a bisection is exact. Within two years, at spacing")
    print(f"  730/n d: records are not nested and the dose phases alias, so every n from 2 to {NSCAN} is")
    print("  scanned; 'first' is the smallest passing n, 'stays' the n from which every larger n passes.")
    need = {}
    ns2 = np.arange(2, NSCAN + 1)
    scan = np.zeros((ns2.size, len(SCEN)))
    for i, n in enumerate(ns2):
        tt = np.arange(n) * (730.0 / n)
        rs = [rec(TH5, tt, s) for s in C.SEEDS]
        for j, (_, kw) in enumerate(SCEN):
            scan[i, j] = med([np.sqrt(var_lnth(r, **kw)) for r in rs])
    for j, (nm, kw) in enumerate(SCEN):
        f = lambda n: sd_design(np.arange(n) * 30.4, kw, True)
        lo, hi = 2, 6000
        if f(lo) <= GATE:
            nm_month = lo
        else:
            while hi - lo > 1:
                mid = (lo + hi) // 2
                lo, hi = (lo, mid) if f(mid) <= GATE else (mid, hi)
            nm_month = hi
        ok = scan[:, j] <= GATE
        first = int(ns2[np.argmax(ok)]) if ok.any() else None
        fail = np.where(~ok)[0]
        stays = (int(ns2[fail[-1] + 1]) if fail.size and fail[-1] + 1 < ns2.size
                 else (int(ns2[0]) if not fail.size else None))
        need[nm] = (nm_month, first, stays)
        fn = lambda k: "none" if k is None else f"{k} (every {730/k:.2f} d, sd {scan[k - 2, j]:.4f})"
        print(f"    {nm:30} monthly: n = {nm_month:4d} ({nm_month*30.4/365.25:5.1f} years);  two years: "
              f"first n = {fn(first)}, stays from n = {fn(stays)}")
    return tabs, (dn, drows), need


# ------------------------------------------------------------------ (c)
KE_LEV = (("exact", 0.0), ("1%", np.log(1.01)), ("2%", np.log(1.02)), ("5%", np.log(1.05)))
B_LEV = (("exact", 0.0), ("1%", 0.01), ("2%", 0.02), ("5%", 0.05))
SC_LEV = (("0", 0.0), ("1%", 0.01), ("2%", 0.02), ("5%", 0.05))
DRAWS = (("730 daily", np.arange(0.0, 730.0, 1.0)),
         ("365 every 2 d", np.arange(0.0, 730.0, 2.0)),
         ("100 every 7.3 d", np.arange(0.0, 730.0, 7.3)),
         ("25 monthly", np.arange(0.0, 730.0, 30.4)))
_REC = {}


def recs(times_key, times):
    if times_key not in _REC:
        _REC[times_key] = [rec(TH5, times, s) for s in C.SEEDS]
    return _REC[times_key]


def sd_req(dkey, ke_sd, b_frac, scatter, km_exact=False):
    rs = recs(dkey, dict(DRAWS)[dkey])
    return med([np.sqrt(var_lnth(r, ke_sd=ke_sd, b_frac=b_frac, scatter=scatter,
                                 km_exact=km_exact)) for r in rs])


def loosest(f, lo, hi):
    """Largest x in [lo, hi] with f(x) <= GATE (f increasing), or None if f(lo) > GATE."""
    if f(lo) > GATE:
        return None
    if f(hi) <= GATE:
        return hi
    for _ in range(40):
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if f(mid) <= GATE else (lo, mid)
    return lo


def part_c():
    print("\n" + "-" * 92)
    print("  (c) REQUIREMENTS FOR sd(ln theta_hat) <= 0.20 AT 5% OF NORMAL (therapy on, 2-year record)")
    print("-" * 92)
    print("  k_e precision = sd(ln k_e) = ln(1 + p); load-bias precision = the b prior as a fraction")
    print("  of L; load scatter = the per-visit load-record term of sigma; the 5% assay CV is kept.")
    print("  '*' marks a cell at or below 0.20.")
    grid = {}
    for dkey, _ in DRAWS:
        print(f"\n    {dkey}:  rows k_e precision / load scatter, columns load-bias precision")
        print(f"      {'k_e':>6} {'scatter':>8}" + "".join(f"{'b '+bl:>11}" for bl, _ in B_LEV))
        for kl, kv in KE_LEV:
            for sl, sv in SC_LEV:
                cells = []
                for bl, bv in B_LEV:
                    v = sd_req(dkey, kv, bv, sv)
                    grid[(dkey, kl, bl, sl)] = v
                    cells.append(f"{v:9.4f}{'*' if v <= GATE else ' '}")
                print(f"      {kl:>6} {sl:>8} " + " ".join(cells))

    # Pareto-loosest passing combinations on the grid
    dk = [d for d, _ in DRAWS]           # index grows with looseness (fewer draws)
    kk = [k for k, _ in KE_LEV]
    bk = [b for b, _ in B_LEV]
    sk = [s for s, _ in SC_LEV]
    passing = [(dk.index(d), kk.index(k), bk.index(b), sk.index(s))
               for (d, k, b, s), v in grid.items() if v <= GATE]
    pareto = [p for p in passing
              if not any(all(q[i] >= p[i] for i in range(4)) and q != p for q in passing)]
    print("\n  The loosest passing combinations on that grid (Pareto set: nothing passes that is")
    print("  looser in one requirement and no tighter in any other):")
    print(f"    {'draws':>16} {'k_e':>6} {'b prior':>8} {'scatter':>8} {'sd':>8}")
    for p in sorted(pareto, key=lambda p: (-p[0], -p[1], -p[2], -p[3])):
        key = (dk[p[0]], kk[p[1]], bk[p[2]], sk[p[3]])
        print(f"    {key[0]:>16} {key[1]:>6} {key[2]:>8} {key[3]:>8} {grid[key]:8.4f}")

    print("\n  Outside the grid: 25 monthly draws with K_m ALSO known exactly (a fifth requirement)")
    print(f"      {'k_e':>6} {'scatter':>8}" + "".join(f"{'b '+bl:>11}" for bl, _ in B_LEV))
    for kl, kv in KE_LEV:
        for sl, sv in SC_LEV:
            cells = []
            for bl, bv in B_LEV:
                v = sd_req("25 monthly", kv, bv, sv, km_exact=True)
                cells.append(f"{v:9.4f}{'*' if v <= GATE else ' '}")
            print(f"      {kl:>6} {sl:>8} " + " ".join(cells))

    print("\n  Continuous requirements, by bisection (None = not reached even at the tightest end):")
    print(f"    {'draws':>16} {'k_e':>6} {'sd at b 5%, sc 5%':>18} {'b needed (sc 5%)':>19} "
          f"{'sc needed (b 5%)':>19} {'both at X':>19}")
    req = {}
    for dkey, _ in DRAWS:
        for kl, kv in KE_LEV:
            base = sd_req(dkey, kv, 0.05, 0.05)
            xb = loosest(lambda x: sd_req(dkey, kv, x, 0.05), 1e-4, 0.2)
            if xb is None and sd_req(dkey, kv, 0.0, 0.05) <= GATE:
                xb = "exact"
            xs = loosest(lambda x: sd_req(dkey, kv, 0.05, x), 1e-4, 0.2)
            if xs is None and sd_req(dkey, kv, 0.05, 0.0) <= GATE:
                xs = "zero"
            xx = loosest(lambda x: sd_req(dkey, kv, x, x), 1e-4, 0.2)
            if xx is None and sd_req(dkey, kv, 0.0, 0.0) <= GATE:
                xx = "exact"
            req[(dkey, kl)] = (base, xb, xs, xx)
            fmt = lambda x: ("None" if x is None else x if isinstance(x, str)
                             else ">= 20% (range end)" if x >= 0.2 else f"{100*x:.2f}%")
            print(f"    {dkey:>16} {kl:>6} {base:18.4f} {fmt(xb):>19} {fmt(xs):>19} {fmt(xx):>19}")

    print("\n  The loosest k_e precision (as ln(1+p), reported as p) for a given load specification:")
    print(f"    {'draws':>16} " + "".join(f"{lbl:>22}" for lbl in
                                        ("b 5%, sc 5% (now)", "b 1%, sc 5%", "b 1%, sc 1%", "load exact")))
    kreq = {}
    for dkey, _ in DRAWS:
        cells = []
        for bv, sv in ((0.05, 0.05), (0.01, 0.05), (0.01, 0.01), (0.0, 0.0)):
            x = loosest(lambda v: sd_req(dkey, v, bv, sv), 1e-4, 0.8)
            if x is None:
                cells.append("exact misses" if sd_req(dkey, 0.0, bv, sv) > GATE else "exact only")
            else:
                cells.append(f"{100*(np.exp(x)-1):.2f}%")
            kreq[(dkey, bv, sv)] = x
        print(f"    {dkey:>16} " + "".join(f"{c:>22}" for c in cells))
    return grid, req, kreq


# ------------------------------------------------------------------ (d)
def part_d():
    print("\n" + "-" * 92)
    print("  (d) THE k_e VARIANCE SHARE 1 - Var(ln th | k_e)/Var(ln th), monthly 25 draws")
    print("-" * 92)
    pcts = np.unique(np.r_[np.geomspace(1.0, 1000.0, 31), 5.0, 100.0 * EDGE / C.TH0, 100.0])
    print(f"    {'pct':>8} {'theta':>8} {'s':>6} {'sd':>8} {'sd|ke':>8} {'share (median)':>15} "
          f"{'[min, max]':>16} {'1-(sd|ke/sd)^2':>15}")
    rows = []
    for pct in pcts:
        th = C.TH0 * pct / 100
        tot, cond = [], []
        for s in C.SEEDS:
            r = C.record(th, seed=s)
            tot.append(var_lnth(r, ke_sd=C.KE_PRIOR))
            cond.append(var_lnth(r))
        tot, cond = np.array(tot), np.array(cond)
        sh = 1 - cond / tot
        sd_t, sd_c = med(np.sqrt(tot)), med(np.sqrt(cond))
        alt = 1 - (sd_c / sd_t) ** 2
        rows.append((pct, th, med(sh), sh.min(), sh.max(), sd_t, sd_c, alt))
        print(f"    {pct:8.2f} {th:8.1f} {C.share(th):6.3f} {sd_t:8.4f} {sd_c:8.4f} {med(sh):15.2%} "
              f"[{sh.min():6.2%},{sh.max():6.2%}] {alt:15.2%}")
    rows = np.array(rows)
    i5 = int(np.argmin(np.abs(rows[:, 0] - 5.0)))
    ie = int(np.argmin(np.abs(rows[:, 1] - EDGE)))
    user = 1 - (0.39 / 3.07) ** 2
    print(f"\n    the user's arithmetic 1 - (0.39/3.07)^2 = {user:.4f}")
    print(f"    at 5%:  median per-record share {rows[i5, 2]:.4f}; from the median sds "
          f"1 - ({rows[i5, 6]:.4f}/{rows[i5, 5]:.4f})^2 = {rows[i5, 7]:.4f}")
    print(f"    at the band edge theta = {rows[ie, 1]:.0f} ({rows[ie, 0]:.1f}% of normal): median share "
          f"{rows[ie, 2]:.4f}; from the median sds {rows[ie, 7]:.4f}; sd there {rows[ie, 5]:.4f}")
    assert abs(rows[i5, 2] - user) < 0.005
    print(f"    |user's figure - median share at 5%| = {abs(rows[i5, 2] - user):.4f} (asserted < 0.005)")
    up = np.diff(rows[:, 2])
    print(f"    the share rises between neighbouring capacities in {int(np.sum(up > 0))} of {up.size} "
          f"steps, by at most {up.max():.2e} (between {rows[np.argmax(up), 0]:.2f}% and "
          f"{rows[np.argmax(up) + 1, 0]:.2f}%); its largest fall is {-up.min():.4f}")
    for target in (0.9, 0.5, 0.1):
        k = int(np.argmax(rows[:, 2] < target))
        print(f"    the share first falls below {target:.0%} between {rows[k-1, 0]:.1f}% and {rows[k, 0]:.1f}% of normal")
    return rows, i5, ie


# ------------------------------------------------------------------ figure
def figure(shap, tabs, dense, rows, i5, ie):
    fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.2))
    a = ax[0]
    pcts = list(shap)
    cols = (C1, C2, C3, C4)
    y = np.arange(len(pcts))[::-1]
    left = np.zeros(len(pcts))
    for f, col in zip(FACTORS, cols):
        w = np.array([shap[p][f] for p in pcts])
        a.barh(y, w, left=left, color=col, height=0.62, edgecolor="white", linewidth=2, label=f)
        for yi, l0, wi in zip(y, left, w):
            if wi > 0.07:
                a.text(l0 + wi / 2, yi, f"{wi:.0%}", ha="center", va="center", fontsize=7.5,
                       color="white" if col in (C1, C2) else INK)
        left += w
    a.set_yticks(y, [f"{p}%" for p in pcts])
    a.set_xlim(0, 1)
    a.set_xlabel("Shapley share of Var(ln θ | k_e exact)")
    a.set_ylabel("capacity, % of normal")
    a.set_title("(a) what holds the k_e-exact variance")
    a.legend(fontsize=7.5, frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.17))
    finish(a)

    b = ax[1]
    ns, rws = tabs[True]
    dn, drows = dense
    for j, col in zip(range(4), cols):
        b.plot(ns, rws[:, j], "-", color=col, lw=2, label=SCEN[j][0])
        b.plot(dn, drows[:, j], "o", ms=8, mfc="white", mec=col, mew=1.6)
    b.axhline(GATE, color=INK, ls="--", lw=1)
    b.text(2000, GATE * 1.06, "target 0.20", fontsize=8, color=INK, ha="center")
    b.axvline(25, color=MUTED, ls=":", lw=0.8)
    b.set_xscale("log"); b.set_yscale("log")
    b.set_xlabel("number of draws n")
    b.set_ylabel(r"sd(ln $\hat\theta$) at 5% of normal")
    b.set_title("(b) lines: monthly, therapy continued; circles: 2 years, denser")
    b.legend(fontsize=7, frameon=False, loc="lower left")
    finish(b); plain(b)

    c = ax[2]
    c.plot(rows[:, 0], rows[:, 2], "-", color=C1, lw=2)
    c.fill_between(rows[:, 0], rows[:, 3], rows[:, 4], color=C1, alpha=0.18, lw=0)
    for i, txt in ((i5, "5%"), (ie, "band edge")):
        c.plot(rows[i, 0], rows[i, 2], "o", ms=8, color=C1, mec="white", mew=2)
        c.annotate(f"{txt}: {rows[i, 2]:.1%}", (rows[i, 0], rows[i, 2]), xytext=(8, -16),
                   textcoords="offset points", fontsize=8, color=INK)
    c.set_xscale("log"); c.set_ylim(0, 1.02)
    c.set_xlabel("capacity, % of normal")
    c.set_ylabel("1 − Var(ln θ | k_e) / Var(ln θ)")
    c.set_title("(c) k_e's share of Var(ln θ), monthly 25 draws")
    finish(c); plain(c)
    save(fig, "t8_floor")


def main():
    print("=" * 92)
    print("  FOLLOW-UP ITEM 1. THE 0.39 FLOOR: sd(ln theta_hat) AT 5% WITH k_e KNOWN EXACTLY")
    print("=" * 92)
    results, shap = part_a()
    tabs, dense, need = part_b()
    part_c()
    rows, i5, ie = part_d()
    figure(shap, tabs, dense, rows, i5, ie)


if __name__ == "__main__":
    main()
