"""
Follow-up: where does ONE untreated, calibrated flux reading reach sd(ln theta_hat) = 0.20,
by Monte Carlo rather than by the information matrix?

THE QUESTION. Task 5 (t5_value.sd5) found by FIM that one direct, calibrated, off-therapy
reading of the enzymatic flux F_enz = theta g(m0) at an untreated visit (plus that visit's
marker draw), added to the monthly record at 5% of normal, reaches the 0.20 gate at a flux
precision Y <= 19.3%. t11_fim_check's Monte Carlo covered only Y = 10% and 20% (240
replicates on the repository's 24 intake seeds), so the crossing by Monte Carlo was not
located. This script locates it.

FOOTING (all imported, nothing restated): core / n4 convention A, theta = 150 (5% of 3000),
K_m = 100, k_e = 5, L = 322.73; 25 monthly draws with the booked dose train (known
increment); n4's sigma; ln k_e ~ N(ln 5, 0.8^2), b ~ N(0, (0.05 L)^2), nothing on ln K_m or
ln theta. The flux visit adds t5_value.base_record's untreated baseline draw and one reading
with t5's sd, F_enz hypot(Y, eps 5%). There is no patient data: every record is simulated,
for ONE simulated patient.

THE ESTIMATOR is t11_fim_check.Problem (imported, not copied): the MAP by variable
projection, theta >= 1e-3, K_m in [5, 5000], k_e in [1e-3, 500], a 41 x 61 grid in
(ln K_m, ln k_e) then 6 local refinements. Part 0 reruns it on t11's own draws and must
reproduce t11's printed RMSEs digit for digit; Part 4 checks it against a dense full-box grid.

THE MONTE CARLO. Every measurement the FIM counts as information is simulated as a noisy
pseudo-observation about the truth: the 26 marker-ordinate draws, the flux reading, the k_e
prior centre and the load-bias prior centre (and the recovery-factor prior centre in
variant ii). FRESH seeds, none used by t11: replicate j draws its intake sequence from
SeedSequence([71212026, 1, j]) and its measurement errors from SeedSequence([71212026, 2, j]);
t11 used integer seeds (intake 0..23 and 1000k+i, default_rng(70000+j), default_rng(90000+s),
default_rng(0)) and its recheck SeedSequence([20261001, ...]). Each replicate has its own
intake record. Common random numbers: the same replicate's standard-normal errors are used at
every Y and in both variants, so differences between scenarios are paired.

VARIANTS AT Y = 15%.
  (i)  three readings averaged: t5's definition, sd = F_enz hypot(Y, eps 5%) / sqrt(3)
       (t5 divides the marker-assay term by sqrt(3) too; the alternative with that term
       shared is printed as a FIM sensitivity only).
  (ii) unknown recovery: the reading is z = c F_enz + noise with a multiplicative
       calibration factor c, ln c ~ N(ln c_true, 0.10^2) read as a pseudo-observation whose
       centre is drawn about c_true = 1. The MAP is then over (ln K_m, ln k_e, ln c) with
       (theta, b) still exact by projection, since for fixed c the reading is
       (z / c) ~ theta g(m0) with sd sdf / c.

Run from this directory:  python3 -W error t12_flux_mc.py > out/t12_flux_mc.txt
"""
import os
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import re
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import kstest, norm

import core as C
import t5_value as T5
import t11_fim_check as T11

T_START = time.time()
HERE = Path(__file__).resolve().parent
TH = C.TH_5
LN_TH = np.log(TH)
POP = C.KE_PRIOR
CV = C.ASSAY_CV
LF = C.LOAD_FRAC
TAG = 71212026                       # this script's seed label
YS = (0.10, 0.125, 0.15, 0.175, 0.20, 0.225)   # 0.225 added beyond the brief's grid, for bracketing
Y_BRIEF = (0.10, 0.125, 0.15, 0.175, 0.20)
Y_VAR = 0.15
CAL_SD = 0.10
YR = (0.10, 0.125, 0.15, 0.175)      # recovery variant: 15% is the brief's; the others are an extension
REPS = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
B = 2000                             # bootstrap resamples
N_GLOBAL = int(sys.argv[2]) if len(sys.argv) > 2 else 100
LN_BOUND = np.log(1e-3 / TH)


def elapsed():
    return f"[{time.time() - T_START:.0f} s]"


def banner(s):
    print("\n" + "=" * 100)
    print("  " + s)
    print("=" * 100)


def rmse(e):
    e = np.asarray(e)
    return float(np.sqrt(np.mean(e ** 2)))


def robust(e):
    e = np.asarray(e)
    return float(0.5 * (np.percentile(e, 84) - np.percentile(e, 16)))


# =====================================================================
#  The FIM for each scenario
# =====================================================================
def flux_sd(r, Y, nread=1, shared_assay=False):
    row, F, eps = T5.flux_row(r)
    if shared_assay:
        return row, F, F * np.hypot(Y / np.sqrt(nread), eps * CV)
    return row, F, F * np.hypot(Y, eps * CV) / np.sqrt(nread)


def fim_sd(r, Y, nread=1, cal_sd=None, shared_assay=False):
    """sd(ln theta) from the record (with the baseline draw) plus the flux reading(s)."""
    F4 = C.fim(r)
    row, F, sdf = flux_sd(r, Y, nread, shared_assay)
    if cal_sd is None:
        return C.sd_th(F4 + np.outer(row, row) / sdf ** 2)
    F5 = np.zeros((5, 5))
    F5[:4, :4] = F4
    row5 = np.r_[row, F]                       # d(c theta g)/d ln c = F at c = 1
    F5 += np.outer(row5, row5) / sdf ** 2
    F5[4, 4] += cal_sd ** -2
    return float(np.sqrt(np.linalg.inv(F5)[0, 0]))


def fim_seeds(**kw):
    """Median over the repository's 24 intake seeds, t5_value.sd5's convention."""
    return float(np.median([fim_sd(T5.base_record(s, baseline=True), **kw) for s in C.SEEDS]))


# =====================================================================
#  The calibrated-reading MAP: (ln K_m, ln k_e, ln c) outer, (theta, b) exact
# =====================================================================
class CalProblem:
    def __init__(self, P, z, sdf, m0, lc0, sc=CAL_SD):
        self.P, self.z, self.sdf, self.m0, self.lc0, self.sc = P, z, sdf, m0, lc0, sc

    def set_c(self, lc):
        c = np.exp(lc)
        self.P.flux = (np.array([self.z / c]), np.array([self.sdf / c]), np.array([self.m0]))

    def resid(self, x):
        self.set_c(x[2])
        return np.r_[self.P.resid(x[:2]), (x[2] - self.lc0) / self.sc]

    def fit(self, offs=np.linspace(-3.0, 3.0, 7), n_lk=41, n_lke=61, top=6, full_box=False):
        P = self.P
        lks = np.linspace(T11.LK_LO, T11.LK_HI, n_lk)
        if full_box:
            a, b = T11.LKE_LO, T11.LKE_HI
        else:
            a = max(T11.LKE_LO, P.kec - 6.0 * P.ke_sd)
            b = min(T11.LKE_HI, P.kec + 6.0 * P.ke_sd)
        lkes = np.linspace(a, b, n_lke)
        cands = []
        for o in offs:
            lc = self.lc0 + self.sc * o
            self.set_c(lc)
            Cg = (P.grid_cost(lks, lkes) + 0.5 * o ** 2).ravel()
            for i in np.argsort(Cg)[:top]:
                cands.append((float(Cg[i]), lks[i // n_lke], lkes[i % n_lke], lc))
        cands.sort(key=lambda t: t[0])
        lo = np.array([T11.LK_LO, T11.LKE_LO, self.lc0 - 2.0])
        hi = np.array([T11.LK_HI, T11.LKE_HI, self.lc0 + 2.0])
        best = None
        for _, a0, b0, c0 in cands[:top]:
            x0 = np.clip(np.array([a0, b0, c0]), lo + 1e-9, hi - 1e-9)
            r = least_squares(self.resid, x0, bounds=(lo, hi), method="trf",
                              xtol=1e-12, ftol=1e-12, gtol=1e-12, diff_step=1e-7)
            if best is None or r.cost < best.cost:
                best = r
        self.set_c(best.x[2])
        th, bb, _, _ = P.inner(*best.x[:2])
        return dict(lk=float(best.x[0]), lke=float(best.x[1]), lc=float(best.x[2]), th=th, b=bb,
                    cost=float(best.cost), gridmin=float(cands[0][0]))


# =====================================================================
#  One replicate
# =====================================================================
def draws(j):
    r = T5.base_record(np.random.SeedSequence([TAG, 1, j]), baseline=True)   # fresh intake
    rng = np.random.default_rng(np.random.SeedSequence([TAG, 2, j]))
    e = rng.standard_normal(r["m"].size)       # 26 marker-ordinate errors
    ef = rng.standard_normal(3)                # up to three flux-reading errors
    eke, eb, ec = rng.standard_normal(3)       # k_e, load-bias and recovery prior centres
    return r, e, ef, eke, eb, ec


def problem(r, e, ef, eke, eb, Y, nread=1):
    sg = C.sigmas(r, "repo")
    y = (r["theta"] + r["inc"]) * r["m"] / (r["Km"] + r["m"]) + r["ke"] * r["m"] + sg * e
    _, F, sdf1 = flux_sd(r, Y)
    z = F + sdf1 * np.mean(ef[:nread])         # mean of nread readings, each sd sdf1
    sdf = sdf1 / np.sqrt(nread)
    sb = LF * r["L"]
    P = T11.Problem(r["m"][None, :], r["inc"][None, :], y[None, :], sg[None, :], [sb], [sb * eb],
                    np.log(C.KE) + POP * eke, POP,
                    flux=(np.array([z]), np.array([sdf]), np.array([r["m_bar"]])))
    return P, z, sdf


KC = "recovery sd .10"


def scen_keys():
    return ([("1 reading", Y, 1, None) for Y in YS] + [("3 readings", Y_VAR, 3, None)]
            + [(KC, Y, 1, CAL_SD) for Y in YR])


def label(key):
    kind, Y, _, _ = key
    return f"{kind}, Y = {100 * Y:.1f}%"


def run_replicate(j, keys):
    r, e, ef, eke, eb, ec = draws(j)
    out = {}
    for key in keys:
        _, Y, nread, cal = key
        P, z, sdf = problem(r, e, ef, eke, eb, Y, nread)
        if cal is None:
            f = P.fit()
            lc = np.nan
        else:
            f = CalProblem(P, z, sdf, r["m_bar"], cal * ec, cal).fit()
            lc = f["lc"]
        at_box = (min(abs(f["lk"] - T11.LK_LO), abs(f["lk"] - T11.LK_HI),
                      abs(f["lke"] - T11.LKE_LO), abs(f["lke"] - T11.LKE_HI)) < 1e-6)
        out[key] = (np.log(f["th"][0]) - LN_TH, fim_sd(r, Y, nread, cal), lc, f["lke"] - np.log(C.KE),
                    f["lk"] - np.log(C.KM0), float(at_box))
    return out, np.r_[e, ef, eke, eb, ec]


# =====================================================================
#  PART 0: the estimator is t11's
# =====================================================================
def part0():
    banner("PART 0  THE ESTIMATOR: t11_fim_check.Problem rerun on t11's OWN Task 5 draws")
    print("  t11 Part 2: 240 replicates, intake seed j % 24, errors from default_rng(70000 + j).")
    txt = (HERE / "out" / "t11_fim_check.txt").read_text()
    for Y in (0.10, 0.20):
        kw = dict(flux_Y=Y)
        e = []
        for j in range(240):
            r = T11.scen_record(j % 24, kw)
            f = T11.scen_problem(r, kw, np.random.default_rng(70_000 + j)).fit()
            e.append(np.log(f["th"][0]) - LN_TH)
        e = np.array(e)
        m = re.search(rf"\(c\) flux reading to {int(100 * Y)}%\s+(\S+)\s+(\S+)\s+(\S+)\s+\[\S+\]\s+(\S+)\s+\[\S+\]\s+(\S+)",
                      txt.split("(ii) Monte Carlo MAP")[1])
        printed = (m.group(3), m.group(4), m.group(5))
        mine = (f"{rmse(e):.4f}", f"{robust(e):.4f}", f"{np.median(e):+.4f}")
        assert mine == printed, (Y, mine, printed)
        print(f"     Y = {100 * Y:4.1f}%: RMSE {mine[0]}, robust {mine[1]}, median {mine[2]};"
              f"  t11 printed {', '.join(printed)}: {'identical' if mine == printed else 'DIFFERENT'}")
    print("  The same estimator is used below, on fresh draws.")
    print(f"  {elapsed()}")


# =====================================================================
#  PART 1: FIM values
# =====================================================================
def part1():
    banner("PART 1  THE FIM, on t5_value.sd5's footing (median over the repository's 24 intake seeds)")
    print(f"     {'scenario':34} {'FIM here':>9} {'t5.sd5':>8} {'match':>6}")
    F = {}
    for Y in YS:
        a = fim_seeds(Y=Y)
        ref = T5.sd5(flux_Y=Y)
        assert abs(a - ref) <= 1e-12 * ref
        F[("1 reading", Y, 1, None)] = a
        print(f"     {'1 reading, Y = %.1f%%' % (100 * Y):34} {a:9.4f} {ref:8.4f} {'yes':>6}")
    a3 = fim_seeds(Y=Y_VAR, nread=3)
    ref3 = T5.sd5(flux_Y=Y_VAR, n_flux=3)
    assert abs(a3 - ref3) <= 1e-12 * ref3
    F[("3 readings", Y_VAR, 3, None)] = a3
    print(f"     {'3 readings averaged, Y = 15%':34} {a3:9.4f} {ref3:8.4f} {'yes':>6}")
    a3s = fim_seeds(Y=Y_VAR, nread=3, shared_assay=True)
    print(f"     {'  (marker-assay term shared)':34} {a3s:9.4f} {'-':>8} {'':>6}  sensitivity only, not simulated")
    for Y in YR:
        ac = fim_seeds(Y=Y, cal_sd=CAL_SD)
        yeff = float(np.hypot(Y, CAL_SD))
        refc = T5.sd5(flux_Y=yeff)
        F[(KC, Y, 1, CAL_SD)] = ac
        print(f"     {'recovery sd(ln c) .10, Y = %.1f%%' % (100 * Y):34} {ac:9.4f} {refc:8.4f} "
              f"{'yes' if abs(ac - refc) <= 1e-9 * refc else 'NO':>6}  (t5.sd5 at Y_eff = hypot(Y, 10%) = {100 * yeff:.2f}%;"
              f" 5x5 FIM with ln c)")
    yc = T5.crossing(lambda v: T5.sd5(flux_Y=v), 1e-3, 1.0)
    ycc = T5.crossing(lambda v: fim_seeds(Y=v, cal_sd=CAL_SD), 1e-3, 1.0)
    yc3 = T5.crossing(lambda v: T5.sd5(flux_Y=v, n_flux=3), 1e-3, 1.0)
    scc = T5.crossing(lambda v: fim_seeds(Y=Y_VAR, cal_sd=v), 1e-4, 1.0)
    print(f"  FIM crossing of 0.20: one reading Y <= {100 * yc:.2f}%; three readings Y <= {100 * yc3:.2f}%;")
    print(f"     one reading with recovery known to 10%: Y <= {100 * ycc:.2f}%; at Y = 15% the recovery must be")
    print(f"     known to sd(ln c) <= {scc:.4f}")
    print(f"  {elapsed()}")
    return F, yc


# =====================================================================
#  PART 2-3: the Monte Carlo
# =====================================================================
def boot_stats(e, idx):
    E = e[idx]
    P16, P84 = np.percentile(E, [16, 84], axis=1)
    return dict(rmse=np.sqrt(np.mean(E ** 2, 1)), robust=0.5 * (P84 - P16), median=np.median(E, 1),
                lin=np.sqrt(np.mean((np.exp(E) - 1.0) ** 2, 1)), atb=np.mean(E <= LN_BOUND + 1e-6, 1))


def point_stats(e):
    return dict(rmse=rmse(e), robust=robust(e), median=float(np.median(e)),
                lin=rmse(np.exp(e) - 1.0), atb=float(np.mean(e <= LN_BOUND + 1e-6)))


def ci90(v):
    return float(np.percentile(v, 5)), float(np.percentile(v, 95))


def crossing(ys, v, gate=C.GATE):
    """First upward crossing of gate by piecewise-linear interpolation. (Y, flag)."""
    ys, v = np.asarray(ys), np.asarray(v)
    if v[0] > gate:
        return np.nan, "below grid"
    for i in range(len(ys) - 1):
        if v[i] <= gate < v[i + 1]:
            return float(ys[i] + (gate - v[i]) * (ys[i + 1] - ys[i]) / (v[i + 1] - v[i])), "ok"
    return np.nan, "above grid"


def crossing_quad(ys, v, gate=C.GATE):
    c = np.polyfit(ys, np.asarray(v) - gate, 2)
    rts = [float(x.real) for x in np.roots(c) if abs(x.imag) < 1e-12 and ys[0] - 0.05 <= x.real <= ys[-1] + 0.05]
    return min(rts) if rts else np.nan


def n_crossings(v, gate=C.GATE):
    s = np.sign(np.asarray(v) - gate)
    return int(np.sum(s[1:] != s[:-1]))


def part2(F_seed, yc_fim):
    keys = scen_keys()
    banner(f"PART 2  MONTE CARLO MAP, {REPS} replicates per scenario, fresh seeds, common random numbers")
    res = {k: [] for k in keys}
    Z = []
    for j in range(REPS):
        out, z = run_replicate(j, keys)
        for k in keys:
            res[k].append(out[k])
        Z.append(z)
        if (j + 1) % 500 == 0:
            print(f"     {j + 1} replicates done {elapsed()}")
    res = {k: np.array(v) for k, v in res.items()}
    Z = np.array(Z)

    print("\n  (a) are the draws representative? standard-normal errors pooled over replicates")
    names = ["26 marker errors"] + ["flux reading 1", "flux reading 2", "flux reading 3", "k_e centre",
                                    "load-bias centre", "recovery centre"]
    cols = [Z[:, :26].ravel()] + [Z[:, 26 + i] for i in range(6)]
    print(f"     {'stream':20} {'n':>7} {'mean':>8} {'sd':>7} {'median':>8} {'frac > 0':>9} {'KS p':>7}")
    for nm, v in zip(names, cols):
        print(f"     {nm:20} {v.size:7d} {np.mean(v):+8.4f} {np.std(v):7.4f} {np.median(v):+8.4f}"
              f" {np.mean(v > 0):9.3f} {kstest(v, 'norm').pvalue:7.3f}")

    rng = np.random.default_rng(np.random.SeedSequence([TAG, 9]))
    idx = rng.integers(0, REPS, (B, REPS))
    BS, PT = {}, {}
    for k in keys:
        e = res[k][:, 0]
        BS[k] = boot_stats(e, idx)
        PT[k] = point_stats(e)

    print("\n  (b) ln theta_hat - ln theta. [..] = percentile-bootstrap 90% interval over replicates "
          f"(B = {B}, paired across scenarios).")
    print("      FIM sd5: t5_value.sd5 (median, 24 repository seeds); FIM rms: over this MC's own records.")
    print(f"     {'scenario':28} {'FIM sd5':>7} {'FIM rms':>7} {'MC RMSE':>7} {'90% CI':>15} {'robust':>7}"
          f" {'90% CI':>15} {'median':>7} {'90% CI':>17} {'lin RMSE':>8} {'90% CI':>15} {'at bnd':>6}"
          f" {'RMSE/FIM':>8} {'rob/FIM':>7} {'P(RMSE<=.2)':>11}")
    FR = {}
    for k in keys:
        p, b = PT[k], BS[k]
        fr = float(np.sqrt(np.mean(res[k][:, 1] ** 2)))
        FR[k] = fr
        c1, c2, c3, c4 = ci90(b["rmse"]), ci90(b["robust"]), ci90(b["median"]), ci90(b["lin"])
        print(f"     {label(k):28} {F_seed[k]:7.4f} {fr:7.4f} {p['rmse']:7.4f} [{c1[0]:.4f},{c1[1]:.4f}]"
              f" {p['robust']:7.4f} [{c2[0]:.4f},{c2[1]:.4f}] {p['median']:+7.4f} [{c3[0]:+.4f},{c3[1]:+.4f}]"
              f" {p['lin']:8.4f} [{c4[0]:.4f},{c4[1]:.4f}] {p['atb']:6.2%} {p['rmse'] / fr:8.3f}"
              f" {p['robust'] / fr:7.3f} {np.mean(b['rmse'] <= C.GATE):11.3f}")
    print(f"     at bnd: fraction of fits at theta = 1e-3 (ln err {LN_BOUND:.2f}). P(RMSE<=.2): bootstrap share.")
    print("\n  (c) further diagnostics per scenario")
    print(f"     {'scenario':28} {'at bnd n':>8} {'at box':>7} {'mean err':>9} {'min err':>8} {'max err':>8}"
          f" {'|err|>ln2':>9} {'ln k_e RMSE':>11} {'ln K_m RMSE':>11} {'ln c RMSE':>9}")
    for k in keys:
        e = res[k][:, 0]
        lc = res[k][:, 2]
        print(f"     {label(k):28} {int(np.sum(e <= LN_BOUND + 1e-6)):8d} {np.mean(res[k][:, 5]):7.2%}"
              f" {np.mean(e):+9.4f} {e.min():+8.3f} {e.max():+8.3f} {np.mean(np.abs(e) > np.log(2)):9.2%}"
              f" {rmse(res[k][:, 3]):11.4f} {rmse(res[k][:, 4]):11.4f}"
              f" {(f'{rmse(lc):9.4f}' if np.all(np.isfinite(lc)) else '        -')}")
    print("     at box: K_m or k_e at its box bound. ln c RMSE: estimated recovery factor against c_true = 1.")
    nz = int(np.sum([np.sum(res[k][:, 0] <= LN_BOUND + 1e-6) for k in keys]))
    print(f"     fits at the theta bound, all scenarios: {nz}; with zero of n fits, the 95% upper bound on the"
          f" rate is 3/n = {3 / REPS:.2%}")
    print(f"  {elapsed()}")
    return keys, res, BS, PT, FR


def lognormal_bench(s):
    """
    ln-scale RMSE and robust half-width of ln(1 + s eps), eps ~ N(0, 1), floored at the MAP's
    theta bound: what an estimator would show if theta_hat/theta - 1 were exactly N(0, s^2)
    (linear-Gaussian in theta at the FIM sd). Deterministic quadrature, no sampling.
    """
    x = np.linspace(-40.0, 40.0, 800_001)
    w = np.exp(-0.5 * x ** 2)
    w /= w.sum()
    le = np.log(np.maximum(1.0 + s * x, 1e-3 / TH))
    zq = float(norm.ppf(0.84))
    return float(np.sqrt(np.sum(w * le ** 2))), float(0.5 * (np.log(1 + zq * s) - np.log(1 - zq * s)))


def flo(t, y0):
    return f"{100 * t:.2f}" if np.isfinite(t) else f"below {100 * y0:g}"


def fhi(t, y1):
    return f"{100 * t:.2f}" if np.isfinite(t) else f"above {100 * y1:g}"


def cross_block(title, kset, BS, PT):
    ys = np.array([k[1] for k in kset])
    print(f"  {title}")
    rows = {}
    for nm, key in (("MC RMSE", "rmse"), ("MC robust half-width", "robust"), ("linear-theta RMSE", "lin")):
        v = np.array([PT[k][key] for k in kset])
        yhat, flag = crossing(ys, v)
        yq = crossing_quad(ys, v) if ys.size >= 3 else np.nan
        vb = np.stack([BS[k][key] for k in kset], 1)
        cb = [crossing(ys, vb[i]) for i in range(B)]
        yb = np.array([c[0] for c in cb])
        fl = [c[1] for c in cb]
        ok = np.isfinite(yb)
        nx = np.array([n_crossings(vb[i]) for i in range(B)])
        print(f"     {nm:22} values at Y = {', '.join(f'{100 * y:g}%' for y in ys)}: {' '.join(f'{x:.4f}' for x in v)}")
        if ok.sum() < 20:
            print(f"     {'':22} crossing: {flag}; resamples inside the grid {np.mean(ok):.1%}")
            rows[nm] = (yhat, np.nan, np.nan, np.nan, flag, float(np.mean(ok)), ys[0], ys[-1])
            continue
        lo, hi = ci90(yb[ok])
        lo95, hi95 = np.percentile(yb[ok], [2.5, 97.5])
        # resamples that do not cross inside the grid lie beyond it: fold them into the interval as such
        #  (sentinels -1 and 9 stand for "below" and "above" the grid; 'nearest' never interpolates them)
        ya = np.where(ok, yb, np.where(np.array(fl) == "below grid", -1.0, 9.0))
        lo_all = float(np.percentile(ya, 5, method="nearest"))
        hi_all = float(np.percentile(ya, 95, method="nearest"))
        lo_all = lo_all if ys[0] <= lo_all <= ys[-1] else np.nan
        hi_all = hi_all if ys[0] <= hi_all <= ys[-1] else np.nan
        rows[nm] = (yhat, lo_all, hi_all, float(np.std(yb[ok])), flag, float(np.mean(ok)), ys[0], ys[-1])
        print(f"     {'':22} crossing Y = {100 * yhat:.2f}% ({flag}); bootstrap 90% [{flo(lo_all, ys[0])}, {fhi(hi_all, ys[-1])}]%"
              f" (crossing resamples only: 90% [{100 * lo:.2f}, {100 * hi:.2f}]%, 95% [{100 * lo95:.2f}, {100 * hi95:.2f}]%,"
              f" sd {100 * np.std(yb[ok]):.2f} pp); quadratic fit over the grid {100 * yq:.2f}%")
        print(f"     {'':22} resamples: crossing inside the grid {np.mean(ok):.1%}, below {fl.count('below grid')},"
              f" above {fl.count('above grid')}; more than one crossing {int(np.sum(nx > 1))}")
    return rows


def part3(keys, res, BS, PT, FR, F_seed, yc_fim):
    banner("PART 3  THE CROSSING OF 0.20 BY MONTE CARLO")
    print("  Piecewise-linear interpolation between grid points; uncertainty from the paired bootstrap")
    print("  (every resample re-interpolates all Y on the same replicate indices). The 90% interval counts")
    print("  resamples that do not cross inside the grid as lying beyond it.")
    k1 = [k for k in keys if k[0] == "1 reading"]
    ys = np.array([k[1] for k in k1])
    rows = cross_block("(a) ONE CALIBRATED READING, Y = 10 ... 22.5%", k1, BS, PT)
    fr = np.array([FR[k] for k in k1])
    fs = np.array([F_seed[k] for k in k1])
    print(f"     {'FIM rms (own records)':22} values {' '.join(f'{x:.4f}' for x in fr)}; crossing "
          f"{100 * crossing(ys, fr)[0]:.2f}% (interpolated)")
    print(f"     {'FIM sd5 (24 seeds)':22} values {' '.join(f'{x:.4f}' for x in fs)}; crossing "
          f"{100 * crossing(ys, fs)[0]:.2f}% interpolated, {100 * yc_fim:.2f}% exact (t5_value.crossing)")
    kb = [k for k in k1 if k[1] in Y_BRIEF]
    yb_ = np.array([k[1] for k in kb])
    for nm, key in (("MC RMSE", "rmse"), ("MC robust half-width", "robust")):
        v = np.array([PT[k][key] for k in kb])
        print(f"     brief's grid only (10-20%), {nm}: crossing {100 * crossing(yb_, v)[0]:.2f}% ({crossing(yb_, v)[1]})")

    print("\n  (b) BENCHMARK: what the ln-scale measures would be if theta_hat/theta - 1 were exactly")
    print("      N(0, FIM rms^2) (linear-Gaussian in theta, no estimator inefficiency); quadrature, not sampled")
    print(f"     {'Y':>6} {'FIM rms':>8} {'bench ln RMSE':>14} {'MC ln RMSE':>11} {'bench robust':>13} {'MC robust':>10}"
          f" {'MC lin RMSE':>12} {'lin RMSE/FIM':>13} {'90% CI':>15}")
    bl, br = [], []
    for k in k1:
        a, b = lognormal_bench(FR[k])
        bl.append(a)
        br.append(b)
        c = ci90(BS[k]["lin"] / FR[k])
        print(f"     {100 * k[1]:5.1f}% {FR[k]:8.4f} {a:14.4f} {PT[k]['rmse']:11.4f} {b:13.4f} {PT[k]['robust']:10.4f}"
              f" {PT[k]['lin']:12.4f} {PT[k]['lin'] / FR[k]:13.3f} [{c[0]:.3f},{c[1]:.3f}]")
    print(f"     benchmark crossings: ln RMSE {100 * crossing(ys, bl)[0]:.2f}%, robust {100 * crossing(ys, br)[0]:.2f}%")

    kr = [k for k in keys if k[0] == KC]
    print()
    rows_c = cross_block("(c) ONE READING WITH RECOVERY KNOWN TO 10% (sd ln c = 0.10), Y = 10 ... 17.5% "
                         "(15% is the brief's point; the others extend it)", kr, BS, PT)
    frc = np.array([FR[k] for k in kr])
    print(f"     {'FIM rms (own records)':22} values {' '.join(f'{x:.4f}' for x in frc)}; crossing "
          f"{100 * crossing([k[1] for k in kr], frc)[0]:.2f}% (interpolated)")
    rows.update({"recovery " + nm: v for nm, v in rows_c.items()})

    print("\n  (d) PAIRED DIFFERENCES at Y = 15%, against one calibrated reading (same replicates); bootstrap 90%")
    base = ("1 reading", Y_VAR, 1, None)
    for k in keys:
        if k[0] == "1 reading" or k[1] != Y_VAR:
            continue
        for nm, key in (("RMSE", "rmse"), ("robust", "robust")):
            d = BS[k][key] - BS[base][key]
            c = ci90(d)
            print(f"     {label(k):30} {nm:7} {PT[k][key]:.4f} - {PT[base][key]:.4f} = {PT[k][key] - PT[base][key]:+.4f}"
                  f"  [{c[0]:+.4f}, {c[1]:+.4f}]")
    print("\n  (e) DOES AN UNKNOWN RECOVERY ACT LIKE A LOOSER CALIBRATED READING? Recovery variant at Y against")
    print("      one calibrated reading at Y_eff = hypot(Y, 10%), interpolated on (a)'s grid; paired bootstrap 90%")
    for k in kr:
        yeff = float(np.hypot(k[1], CAL_SD))
        parts = []
        for nm, key in (("RMSE", "rmse"), ("robust", "robust")):
            v = np.array([PT[q][key] for q in k1])
            at = float(np.interp(yeff, ys, v))
            vb = np.stack([BS[q][key] for q in k1], 1)
            atb = np.array([np.interp(yeff, ys, vb[i]) for i in range(B)])
            c = ci90(BS[k][key] - atb)
            parts.append(f"{nm} {PT[k][key]:.4f} vs {at:.4f}, diff {PT[k][key] - at:+.4f} [{c[0]:+.4f}, {c[1]:+.4f}]")
        print(f"     Y = {100 * k[1]:4.1f}% (Y_eff {100 * yeff:5.2f}%): " + "; ".join(parts))
    for nm in ("MC RMSE", "MC robust half-width"):
        y1 = rows[nm][0]
        yd = np.sqrt(max(y1 ** 2 - CAL_SD ** 2, 0.0))
        print(f"     {nm}: crossing with recovery known to 10%, derived as sqrt(Y*^2 - 0.10^2) from (a): {100 * yd:.2f}%;"
              f" simulated directly in (c): {100 * rows['recovery ' + nm][0]:.2f}%")
    print(f"  {elapsed()}")
    return rows


# =====================================================================
#  PART 4: global optimum checks
# =====================================================================
def part4():
    banner(f"PART 4  IS THE MAP THE GLOBAL OPTIMUM? first {N_GLOBAL} replicates, dense full-box search")
    print("  Standard fit as above. Check: 161 x 241 grid over the whole (K_m, k_e) box, 20 refinements;")
    print("  for the recovery variant also 21 values of ln c over +-5 prior sd. (theta, b) exact throughout.")
    print(f"     {'scenario':28} {'reps':>5} {'max cost gain':>14} {'gains > 1e-6':>13} {'max |d ln theta|':>17}")
    for key in [("1 reading", 0.10, 1, None), ("1 reading", 0.15, 1, None), ("1 reading", 0.20, 1, None),
                ("1 reading", 0.225, 1, None), ("3 readings", Y_VAR, 3, None), (KC, Y_VAR, 1, CAL_SD),
                (KC, 0.175, 1, CAL_SD)]:
        _, Y, nread, cal = key
        gain, cnt, dth = -np.inf, 0, 0.0
        for j in range(N_GLOBAL):
            r, e, ef, eke, eb, ec = draws(j)
            P, z, sdf = problem(r, e, ef, eke, eb, Y, nread)
            if cal is None:
                f = P.fit()
                g = P.fit(n_lk=161, n_lke=241, top=20, full_box=True)
            else:
                f = CalProblem(P, z, sdf, r["m_bar"], cal * ec, cal).fit()
                g = CalProblem(P, z, sdf, r["m_bar"], cal * ec, cal).fit(
                    offs=np.linspace(-5, 5, 21), n_lk=161, n_lke=241, top=20, full_box=True)
            gain = max(gain, f["cost"] - g["cost"])
            cnt += (f["cost"] - g["cost"]) > 1e-6
            dth = max(dth, abs(np.log(f["th"][0]) - np.log(g["th"][0])))
        print(f"     {label(key):28} {N_GLOBAL:5d} {gain:14.2e} {cnt:13d} {dth:17.2e}")
    # the recovery MAP collapses to the calibrated one as sd(ln c) -> 0
    d = 0.0
    for j in range(N_GLOBAL):
        r, e, ef, eke, eb, ec = draws(j)
        P, z, sdf = problem(r, e, ef, eke, eb, Y_VAR)
        f0 = P.fit()
        P2, z2, sdf2 = problem(r, e, ef, eke, eb, Y_VAR)
        f1 = CalProblem(P2, z2, sdf2, r["m_bar"], 0.0, 1e-6).fit()
        d = max(d, abs(np.log(f0["th"][0]) - np.log(f1["th"][0])))
    print(f"     recovery MAP with sd(ln c) = 1e-6 and centre 0 against the calibrated MAP, Y = 15%: "
          f"max |d ln theta| = {d:.2e}")
    print(f"  {elapsed()}")


# =====================================================================
def part5(keys, PT, rows):
    banner("PART 5  AGAINST t11 AND ITS RECHECK, at Y = 10% and 20% (as printed in out/)")
    t11 = (HERE / "out" / "t11_fim_check.txt").read_text().split("(ii) Monte Carlo MAP")[1]
    rc = (HERE / "out" / "t11_fim_check_recheck.txt").read_text()
    for Y in (0.10, 0.20):
        k = ("1 reading", Y, 1, None)
        m = re.search(rf"\(c\) flux reading to {int(100 * Y)}%\s+\S+\s+\S+\s+(\S+)\s+\[(\S+),(\S+)\]\s+(\S+)", t11)
        m2 = re.search(rf"flux {int(100 * Y)}%\s+FIM.*?ln RMSE (\S+) \[(\S+), (\S+)\].*?robust (\S+)", rc)
        print(f"     Y = {100 * Y:4.1f}%: here RMSE {PT[k]['rmse']:.4f}, robust {PT[k]['robust']:.4f} ({REPS} reps);"
              f" t11 RMSE {m.group(1)} [95% {m.group(2)}, {m.group(3)}], robust {m.group(4)} (240 reps);"
              f" recheck RMSE {m2.group(1)} [95% {m2.group(2)}, {m2.group(3)}], robust {m2.group(4)} (400 reps)")


def summary(keys, PT, BS, FR, F_seed, rows, yc):
    banner("SUMMARY (every number repeated from above)")
    for nm in ("MC RMSE", "MC robust half-width", "linear-theta RMSE"):
        y, lo, hi, sd, flag, frac, y0, y1 = rows[nm]
        print(f"  one calibrated reading, {nm:21} crosses 0.20 at Y = {100 * y:.2f}%  (bootstrap 90% {flo(lo, y0)}% to"
              f" {fhi(hi, y1)}%, sd {100 * sd:.2f} pp)")
    print(f"  one calibrated reading, FIM (t5_value.sd5)        crosses 0.20 at Y = {100 * yc:.2f}%")
    for nm in ("MC RMSE", "MC robust half-width"):
        y, lo, hi, sd, flag, frac, y0, y1 = rows["recovery " + nm]
        if np.isfinite(y):
            print(f"  recovery known to 10%,  {nm:21} crosses 0.20 at Y = {100 * y:.2f}%  (bootstrap 90% {flo(lo, y0)}% to"
                  f" {fhi(hi, y1)}%; resamples crossing inside 10-17.5%: {frac:.1%})")
        else:
            print(f"  recovery known to 10%,  {nm:21} crossing {flag} (10-17.5%)")
    for k in keys:
        if k[1] == Y_VAR:
            c = ci90(BS[k]["rmse"])
            print(f"  {label(k):30} FIM {F_seed[k]:.4f}; MC RMSE {PT[k]['rmse']:.4f} [{c[0]:.4f}, {c[1]:.4f}], robust"
                  f" {PT[k]['robust']:.4f}, median {PT[k]['median']:+.4f}, P(RMSE <= 0.20) {np.mean(BS[k]['rmse'] <= C.GATE):.3f}")


def main():
    print("=" * 100)
    print("  FOLLOW-UP: THE FLUX-READING CROSSING OF sd(ln theta_hat) = 0.20 BY MONTE CARLO")
    print("=" * 100)
    print("  Footing (core / n4, convention A): theta = 150 (5%), K_m = 100, k_e = 5, L = 322.73, 25 monthly")
    print("  draws with the booked dose train, n4's sigma, k_e prior 0.8, load-bias prior 0.05 L; plus one")
    print("  untreated visit (marker draw + flux reading). ONE simulated patient; there is no patient data.")
    print(f"  Y grid {[100 * y for y in YS]}% (22.5% added beyond the brief's 10-20% to bracket the crossing);")
    print(f"  {REPS} replicates per scenario; bootstrap B = {B}; seed label {TAG}.")
    part0()
    F_seed, yc = part1()
    keys, res, BS, PT, FR = part2(F_seed, yc)
    rows = part3(keys, res, BS, PT, FR, F_seed, yc)
    part4()
    part5(keys, PT, rows)
    summary(keys, PT, BS, FR, F_seed, rows, yc)
    print(f"\n  total runtime {time.time() - T_START:.0f} s")


if __name__ == "__main__":
    main()
