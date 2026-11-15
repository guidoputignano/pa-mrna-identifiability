"""
Round 20: two questions the owner asked about the dosing result.

  1. How late can a fortnightly dose be, and what does a fully skipped dose cost?
     A. Quasi-static: the trough condition theta + Delta(t) > theta_crit_int at the
        end of a delayed gap. Same-rate dose train (t44/t48/t50/t53):
        the increment just before an on-time dose is Delta_tr(14); a dose that is
        d days late lets it decay to Delta_tr(14) exp(-lambda d).
     B. Dynamic: the parent model's two-state system (src/common/models.py rhs_C,
        k = 1) with the dose train added to the enzyme channel,
            dm/dt = s - (theta + Delta(t)) g(m) - k c m/(K_c + m),
            dc/dt = R0 - k c m/(K_c + m) - rho0 k (1 + m/M_i) c,
        at the baseline load s = 600, nominal carnitine constants. Twelve
        fortnightly doses settle the cycle; the thirteenth is given d days late
        (or skipped, a 28-day gap); five more follow on time. Reported: the
        metabolite's usual peak, its peak in the long gap, and whether it returns
        to the usual cycle. Whether a given peak is clinically dangerous is not
        something the model defines (the paper treats the tolerable level as a
        convention); as a scale, the paper's illustrative 10 %-capacity patient
        (theta = 300) sits stably at 82 uM untreated at this load.

  2. Is sd(ln theta) = 0.44 (the k_e tracer plus the load turnover) enough to decide
     the interval per patient? A patient can be booked on a longer interval tau if
     theta > theta_need = theta_crit_int - Delta_tr(tau). With an estimate
     theta_hat = theta + s(theta) Z, where s(theta) = theta sd(ln theta) is the
     Fisher-information sd at that capacity (t26's machinery, recomputed at each
     capacity: precision is not one number), certify when theta_hat - 1.645 s > theta_need
     (one-sided 95 %). Population: theta log-uniform over 0.17-5 % of normal, as in R-J.
     FIM-based, not simulated fits (R-M's rule asks for simulated fits for decisions).

Run: python -W error t54_dose_timing.py
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from common.models import ParamsC, g  # noqa: E402
import model_C.run as MC  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parent))
import core as C  # noqa: E402

D0, TAU0 = 6000.0, 14.0
D_ADM = 2672.0                   # paper 11_scope: largest admissible increment at 5 %
TH_NORMAL = 3000.0
S = 600.0
K = 1.0
THC_NOM = MC.theta_crit_intrinsic(MC.PC, MC.S_BASE)
THC_WORST = 232.1                # t53: worst corner of the paper's ranges


def trough(tau, thalf, D=D0):
    x = np.exp(-np.log(2.0) / thalf * tau)
    return D * (tau / TAU0) * x / (1.0 - x)


def max_delay(theta, thalf, thc, D=D0):
    """Days a fortnightly dose can be late before theta + Delta < theta_crit."""
    need = thc - theta
    if need <= 0:
        return np.inf
    lam = np.log(2.0) / thalf
    return np.log(trough(14, thalf, D) / need) / lam


def sd_ln_at(theta, sd_ke, turnover, seeds=None):
    """t26.sd_theta at an arbitrary capacity (t26 fixes 5 % of normal)."""
    import core as CC
    seeds = CC.SEEDS if seeds is None else seeds

    def one(s):
        r = CC.record(theta, seed=s)
        L = r["L"]
        b_pr = CC.LOAD_FRAC * L
        if turnover is not None:
            b_pr = 1.0 / np.sqrt(1.0 / (CC.LOAD_FRAC * L) ** 2 + 1.0 / (turnover * L) ** 2)
        return CC.sd_th(CC.fim(r, noise="perdraw", ke_prior=sd_ke, b_prior=b_pr))
    return CC.over_seeds(one, seeds)[0]


def simulate(theta, thalf, D=D0, delay=0.0, skip=False, n_pre=12, n_post=5):
    lam = np.log(2.0) / thalf
    p = ParamsC(**{**MC.PC.__dict__, "theta": float(theta)})
    doses = [TAU0 * k for k in range(n_pre)]
    nxt = doses[-1] + (2 * TAU0 if skip else TAU0 + delay)
    doses = np.array(doses + [nxt + TAU0 * k for k in range(n_post)])

    def cap(t):
        past = doses[doses <= t + 1e-12]
        return theta + D * np.exp(-lam * (t - past)).sum()

    def rhs(t, y):
        m, c = y
        conj = K * c * m / (p.Kc + m)
        delta = p.rho0 * K * (1.0 + m / p.Mi)
        return [S - cap(t) * g(m, p.Km) - conj, p.R0 - conj - delta * c]

    y = [20.0, 500.0]
    T, M = [], []
    edges = list(doses) + [doses[-1] + TAU0]
    for a, b in zip(edges[:-1], edges[1:]):
        sol = solve_ivp(rhs, (a, b), y, method="LSODA", rtol=1e-8, atol=1e-8, max_step=0.05)
        T.append(sol.t)
        M.append(sol.y[0])
        y = sol.y[:, -1]
    T, M = np.concatenate(T), np.concatenate(M)
    usual = M[(T >= doses[n_pre - 2]) & (T < doses[n_pre - 1])].max()
    gap = M[(T >= doses[n_pre - 1]) & (T <= doses[n_pre] + 3.0)].max()
    final = M[T >= doses[-1]].max()
    return usual, gap, final


def main():
    print("=" * 92)
    print("  ROUND 20. LATE AND MISSED DOSES, AND WHETHER 0.44 DECIDES THE INTERVAL PER PATIENT")
    print("=" * 92)
    print(f"  theta_crit_int nominal {THC_NOM:.1f}, worst corner of the paper's ranges {THC_WORST:.1f}"
          f" uM/d; normal capacity {TH_NORMAL:.0f}")

    print("\n  A. QUASI-STATIC: days a FORTNIGHTLY dose can be late before the trough condition fails")
    print("     (inf = never fails: the patient's own capacity alone clears theta_crit)")
    caps = [(0.0, "0 %"), (30.0, "1 %"), (75.0, "2.5 %"), (150.0, "5 %")]
    for D in (D0, D_ADM):
        for thc, lab in ((THC_NOM, "nominal"), (THC_WORST, "worst corner")):
            print(f"\n    D = {D:.0f}, constants {lab}:")
            print(f"      {'half-life':>9} " + " ".join(f"{'theta '+c:>11}" for _, c in caps))
            for th in (3.5, 3.6, 4.3):
                row = []
                for theta, _ in caps:
                    d = max_delay(theta, th, thc, D)
                    row.append("inf" if not np.isfinite(d) else ("fails on time" if d < 0 else f"{d:.1f} d"))
                print(f"      {th:>8.1f}d " + " ".join(f"{r:>11}" for r in row))
    print("\n    A fully SKIPPED dose (28-day gap) keeps the trough condition only for theta above:")
    for D in (D0, D_ADM):
        for th in (3.5, 3.6, 4.3):
            lam = np.log(2.0) / th
            left = trough(14, th, D) * np.exp(-lam * 14)
            print(f"      D = {D:.0f}, half-life {th} d: increment left {left:6.1f} -> theta > "
                  f"{THC_NOM - left:6.1f} uM/d ({100*(THC_NOM-left)/TH_NORMAL:.1f} % of normal), nominal constants")

    print("\n  B. DYNAMIC (two-state parent model, load 600, nominal constants): metabolite m in uM")
    print("     usual = peak over a normal 14-day cycle; then the peak in the long gap; every run")
    print("     returns to its usual cycle once doses resume (final cycle peak shown).")
    for D in (D0, D_ADM):
        for th in (3.5, 4.3):
            print(f"\n    D = {D:.0f}, half-life {th} d:")
            print(f"      {'capacity':>12} {'usual':>7} {'1 d late':>9} {'3 d late':>9} {'7 d late':>9} {'skipped':>9} {'final':>7}")
            for theta, lab in caps:
                runs = [simulate(theta, th, D, delay=d) for d in (1.0, 3.0, 7.0)]
                sk = simulate(theta, th, D, skip=True)
                print(f"      {lab+' ('+str(int(theta))+')':>12} {runs[0][0]:7.1f} {runs[0][1]:9.1f} {runs[1][1]:9.1f}"
                      f" {runs[2][1]:9.1f} {sk[1]:9.1f} {sk[2]:7.1f}")
    p300 = ParamsC(**{**MC.PC.__dict__, "theta": 300.0})
    from common.models import clearance_C_quasistatic, equilibria_nonmonotone
    eq = equilibria_nonmonotone(clearance_C_quasistatic, p300, S, np.linspace(1e-3, 5000, 20000))[0]
    print(f"\n    Scale: the paper's illustrative patient (theta = 300, 10 % of normal) sits stably at"
          f" {eq[0]:.0f} uM untreated at this load (separatrix {eq[1]:.0f} uM).")

    print("\n  C. CAN AN ESTIMATE OF theta DECIDE A LONGER INTERVAL PER PATIENT?")
    print("     Precision is NOT one number: sd(ln theta) grows as capacity falls. Recomputed with")
    print("     t26's machinery (25 monthly draws, per-draw noise) at each capacity:")
    pcts = (0.17, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0)
    confs = (("record alone", C.KE_PRIOR, None), ("k_e tracer + load turnover", 0.054, 0.10),
             ("k_e and unrecorded load exact", 1e-6, 1e-6))
    tab = {}
    print(f"      {'% of normal':>11} " + " ".join(f"{c[0]:>30}" for c in confs))
    for pct in pcts:
        th = TH_NORMAL * pct / 100
        row = [sd_ln_at(th, ke, tu) for _, ke, tu in confs]
        tab[pct] = row
        print(f"      {pct:>11.2f} " + " ".join(f"{v:>21.3f} (+-{v*th:5.0f} uM/d)" for v in row))
    print("     In absolute terms the uncertainty is nearly the same at every capacity (sd(ln theta)")
    print("     x theta, in brackets): the likelihood is close to Gaussian in theta itself. So the")
    print("     model below is theta_hat = theta + s(theta) Z, with s = theta sd(ln theta) interpolated")
    print("     in log capacity and held constant beyond the computed range.")
    print("     Population: theta log-uniform over 0.17-5 % of normal (R-J). Certify if the one-sided")
    print("     95 % lower bound theta_hat - 1.645 s exceeds theta_need = theta_crit - Delta_tr(tau).")
    lp = np.log(np.array(pcts))
    rng = np.random.default_rng(0)
    n = 200_000
    theta = TH_NORMAL * 10 ** rng.uniform(np.log10(0.0017), np.log10(0.05), n)
    z = rng.standard_normal(n)
    lt = np.log(100 * theta / TH_NORMAL)
    cases = [("21 d, half-life 3.5 d, D 6000, nominal", THC_NOM - trough(21, 3.5)),
             ("21 d, half-life 3.5 d, D 6000, worst corner", THC_WORST - trough(21, 3.5)),
             ("21 d, half-life 4.3 d, D 2672, nominal", THC_NOM - trough(21, 4.3, D_ADM)),
             ("14 d, half-life 3.5 d, D 2672, worst corner", THC_WORST - trough(14, 3.5, D_ADM))]
    for lab, need in cases:
        safe = theta > need
        print(f"\n    {lab}: theta_need = {need:.1f} uM/d ({100*need/TH_NORMAL:.2f} % of normal);"
              f" truly safe {safe.mean():.0%} of the population")
        print(f"      {'configuration':>30} {'certified':>10} {'of the safe':>12} {'false-safe among unsafe':>24}"
              f" {'capacity certified 50 % / 80 % of the time':>44}")
        for k, (name, _, _) in enumerate(confs):
            s_abs = [tab[q][k] * TH_NORMAL * q / 100 for q in pcts]
            s_lin = np.interp(lt, lp, s_abs)
            est = theta + s_lin * z
            cert = est - 1.645 * s_lin > need
            fs = (cert & ~safe).sum() / max((~safe).sum(), 1)
            def s_of(th, s_abs=s_abs):
                return np.interp(np.log(100 * th / TH_NORMAL), lp, s_abs)
            t50 = brentq(lambda th: th - 1.645 * s_of(th) - need, 1e-3, 1e5)
            t80 = brentq(lambda th: th - (1.645 + 0.8416) * s_of(th) - need, 1e-3, 1e5)
            print(f"      {name:>30} {cert.mean():>10.0%} {cert[safe].mean():>12.0%} {fs:>24.1%}"
                  f" {100*t50/TH_NORMAL:>27.1f} % / {100*t80/TH_NORMAL:.1f} % of normal")
    print("""
  READING. (1) In the parent model a late dose is not a cliff edge. The trough condition is
  a quasi-static test; dynamically the metabolite rises gradually while the increment is
  low and falls back once the next dose arrives. How far it rises depends on the residual
  capacity, the half-life and the dose: see B. (2) The 0.44 of the two tracers is the
  precision at 5 % of normal only; in absolute terms the uncertainty is about the same at
  every capacity, so an estimate certifies a longer interval only for patients well above
  theta_need (see the lowest certified capacity). These are model statements on simulated
  patients; the decision numbers in C are FIM-based (R-M's rule asks for simulated fits).
""")


if __name__ == "__main__":
    main()
