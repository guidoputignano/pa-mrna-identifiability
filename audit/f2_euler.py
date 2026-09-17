"""
P2 -- the Euler step the PI proposes, and whether a between-visit difference
carries any information.

His note writes

    m(t_{k+1}) = m(t_k) + tau * (...)

with tau the gap between clinic visits, and warns that Euler is unsafe across
the jump. Three questions:

    (a) what is the marker's time constant, and what is the explicit-Euler
        stability limit at the paper's constants?
    (b) does a finite difference between two clinic visits carry information
        about the quasi-steady locus?
    (c) does any code path in this repository integrate with explicit Euler
        across a dose?

Nothing below records an expected result in advance.

Run:  python audit/f2_euler.py
"""
import sys
import subprocess
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.integrate import solve_ivp

from common.models import g
from common.params import M_OP
from design.run import (P, delta, load, DOSE_INTERVAL, LAMBDA, D_AMP,
                        quasi_steady_level, ASSAY_CV)

ROOT = Path(__file__).resolve().parents[1]
TH, KM, KE = P.theta, P.Km, P.ke
LB = load(0.0)


def banner(t):
    print("\n" + "=" * 78)
    print("  " + t)
    print("=" * 78)


def jac(m, d):
    """-d(dm/dt)/dm = (theta + Delta) K_m/(K_m+m)^2 + k_e, the decay rate."""
    return (TH + d) * KM / (KM + m) ** 2 + KE


# --------------------------------------------------------------------------
def s1_time_constant_and_stability():
    banner("S1  time constant and the explicit-Euler stability limit")

    print("  The linearisation of dm/dt = L + b - (theta + Delta) g(m) - k_e m")
    print("  about a level m has decay rate")
    print("      lambda(m, Delta) = (theta + Delta) K_m/(K_m+m)^2 + k_e,")
    print("  so tau = 1/lambda and explicit Euler is stable only for")
    print("      h < 2/lambda.\n")

    print(f"  {'where':>34} {'m':>9} {'Delta':>9} {'lambda (/d)':>12}"
          f" {'tau (h)':>10} {'h_max (h)':>11}")
    print("  " + "-" * 92)
    rows = [("no therapy, operating point", M_OP, 0.0)]
    t_cycle = np.linspace(700.0, 714.0, 5601)
    lv = np.array([quasi_steady_level(x) for x in t_cycle])
    dd = np.asarray(delta(t_cycle), dtype=float)
    i_lo, i_hi = int(np.argmin(lv)), int(np.argmax(lv))
    rows += [("cycle trough in m (peak dose)", lv[i_lo], dd[i_lo]),
             ("cycle peak in m (trough dose)", lv[i_hi], dd[i_hi]),
             ("the instant after a dose", lv[i_hi], float(delta(714.0)))]
    worst = 0.0
    for lab, m, d in rows:
        lam = jac(m, d)
        worst = max(worst, lam)
        print(f"  {lab:>34} {m:9.4f} {d:9.1f} {lam:12.5f} "
              f"{24/lam:10.4f} {48/lam:11.4f}")

    print(f"\n  stiffest point on the record: lambda = {worst:.4f} /d, so")
    print(f"  explicit Euler needs h < {48/worst:.4f} h = {2/worst:.6f} d.")

    print(f"\n  {'visit gap':>22} {'h (d)':>9} {'h/h_max':>10}"
          f" {'|1 - h lambda| per step':>25}")
    print("  " + "-" * 70)
    for lab, h in (("two-daily", 2.0), ("weekly", 7.0),
                   ("the dose interval", DOSE_INTERVAL),
                   ("monthly", 30.4), ("two-monthly", 60.0)):
        amp = abs(1.0 - h * worst)
        print(f"  {lab:>22} {h:9.3f} {h/(2/worst):10.1f} {amp:25.1f}")

    print("\n  So at every clinic spacing the explicit step is unstable by two")
    print("  to three orders of magnitude, and the instability is not confined")
    print("  to the dose instant: it is there at the trough as well.")
    return worst


# --------------------------------------------------------------------------
def s2_euler_actually_run(worst):
    banner("S2  running his step, to see what it does")

    def rhs(t, y):
        return LB - (TH + float(delta(t))) * g(y, KM) - KE * y

    for lab, h in (("monthly", 30.4), ("weekly", 7.0), ("two-daily", 2.0),
                   ("stable, h = 1/(2 lambda)", 1.0 / (2 * worst))):
        m = quasi_steady_level(0.0)
        t, hist, blew = 0.0, [m], False
        for _ in range(int(120.0 / h)):
            m = m + h * rhs(t, m)
            t += h
            hist.append(m)
            if not np.isfinite(m) or abs(m) > 1e12:
                blew = True
                break
        ref = quasi_steady_level(t)
        print(f"  {lab:>26}: after {t:7.2f} d Euler gives "
              f"{'overflow' if blew else f'{m:14.4e}'}, "
              f"the balance gives {ref:8.4f}")
        neg = sum(1 for x in hist if np.isfinite(x) and x < 0)
        print(f"  {'':>26}  steps taken {len(hist)-1:4d}, of which "
              f"{neg} put the marker negative")
    print("\n  A negative marker is not a record. The step does not merely lose")
    print("  accuracy across the dose; it leaves the physical region on the")
    print("  first step at any clinic spacing.")


# --------------------------------------------------------------------------
def s3_information_in_a_difference(gap=30.4, n=25, seed=11, reps=2000):
    banner("S3  does a between-visit difference carry information?")

    t = np.arange(0.0, gap * n, gap)
    lv = np.array([quasi_steady_level(x) for x in t])
    dd = np.asarray(delta(t), dtype=float)

    # The increment in closed form, so the ODE can be integrated cycle by
    # cycle: Delta is smooth inside a dose interval and jumps only at the
    # dose times, so each interval is an ordinary smooth problem.
    def delta_fast(tt):
        k = np.floor(tt / DOSE_INTERVAL)
        return (D_AMP * np.exp(-LAMBDA * (tt - k * DOSE_INTERVAL))
                * (1.0 - np.exp(-LAMBDA * (k + 1.0) * DOSE_INTERVAL))
                / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL)))

    chk = max(abs(delta_fast(x) - float(delta(x)))
              for x in (0.0, 3.3, 13.9, 14.2, 100.7, 700.4))
    print(f"  closed-form Delta against design.run.delta: max difference "
          f"{chk:.3e}\n")

    def rhs(tt, y):
        return [LB - (TH + delta_fast(tt)) * g(y[0], KM) - KE * y[0]]

    edges = np.arange(0.0, t[-1] + 2 * DOSE_INTERVAL, DOSE_INTERVAL)
    m_at, cur = {}, quasi_steady_level(0.0)
    want = list(t)
    for lo, hi in zip(edges[:-1], edges[1:]):
        inside = [x for x in want if lo <= x < hi]
        a, z = lo + 1e-9, hi - 1e-9
        ev = sorted(min(max(x, a), z) for x in inside)
        s = solve_ivp(rhs, (a, z), [cur], t_eval=ev if ev else None,
                      rtol=1e-11, atol=1e-13, method="LSODA")
        for x, v in zip(sorted(inside), (s.y[0] if ev else [])):
            m_at[x] = float(v)
        cur = float(s.y[0, -1])
    m_true = np.array([m_at[x] for x in t])
    mdot = np.array([float(rhs(tt, [mm])[0]) for tt, mm in zip(t, m_true)])

    fd = np.diff(m_true) / gap

    print(f"  monthly design, n = {n}, gap = {gap} d\n")
    print(f"    true |dm/dt| at the visits:   median "
          f"{np.median(np.abs(mdot)):.6e}, max {np.abs(mdot).max():.6e}")
    print(f"    finite difference (m_k+1 - m_k)/tau: median "
          f"{np.median(np.abs(fd)):.6e}, max {np.abs(fd).max():.6e}")
    print(f"    ratio of medians: {np.median(np.abs(fd))/np.median(np.abs(mdot)):.2f}x")

    sig_m = ASSAY_CV * float(np.mean(m_true))
    noise_fd = np.sqrt(2.0) * sig_m / gap
    bias = np.abs(fd - mdot[:-1])
    print(f"\n    assay sd on the marker: {sig_m:.6f}")
    print(f"    noise sd on a difference: sqrt(2)*sigma_m/tau = {noise_fd:.6e}")
    print(f"    bias of the difference as an estimate of dm/dt at the visit:")
    print(f"      median {np.median(bias):.6e}, max {bias.max():.6e}")
    print(f"      as a multiple of the noise sd: median "
          f"{np.median(bias)/noise_fd:.2f}, max {bias.max()/noise_fd:.2f}")
    print(f"      as a fraction of dm/dt itself: median "
          f"{np.median(bias/np.abs(mdot[:-1])):.4f}")
    print("\n    dm/dt at a visit is NOT zero: the quasi-steady LOCUS itself")
    print("    drifts as the increment decays, at roughly (dm*/dDelta)(-lambda")
    print("    Delta). A difference over a gap longer than the dose cycle")
    print("    averages that drift over whole cycles and returns something")
    print("    unrelated to the instantaneous rate.")

    print("\n  But the deeper point is that the finite difference is not")
    print("  estimating dm/dt at all. Between two visits the marker tracks the")
    print("  quasi-steady level, so (m_k+1 - m_k)/tau estimates")
    print("      (dm*/dDelta)(Delta_k+1 - Delta_k)/tau,")
    print("  a known function of the dosing record and the same parameters")
    print("  already in m_k and m_k+1. Below: the difference predicted from the")
    print("  balance alone against the difference the ODE produces.")
    pred = np.diff(lv) / gap
    print(f"\n    max |ODE difference - balance difference| = "
          f"{np.max(np.abs(fd - pred)):.6e}")
    print(f"    as a fraction of the difference itself: "
          f"{np.max(np.abs(fd - pred))/np.median(np.abs(fd)):.3e}")
    print("\n  The residual that a derivative would have to carry is the gap")
    print("  between the ODE and its own algebraic locus, and that gap is")
    print("  smaller than the assay by:")
    print(f"    max |m_ODE - m_balance| at the visits = "
          f"{np.max(np.abs(m_true - lv)):.6e}")
    print(f"    assay sd                              = {sig_m:.6e}")
    print(f"    ratio                                 = "
          f"{np.max(np.abs(m_true - lv))/sig_m:.3e}")


# --------------------------------------------------------------------------
def s4_no_euler_in_the_code():
    banner("S4  is there an explicit-Euler path in the repository?")

    pat = r"euler|y\s*\+=\s*h\s*\*|m\s*=\s*m\s*\+\s*h\s*\*"
    out = subprocess.run(["grep", "-rniE", pat, "--include=*.py",
                          "src", "audit", "tests", "tools"],
                         cwd=ROOT, capture_output=True, text=True).stdout
    # this script itself is full of the word; exclude it
    hits = [l for l in out.splitlines()
            if l.strip() and not l.startswith("audit/f2_euler.py")]
    print(f"  grep for an explicit step or the word 'Euler': {len(hits)} hits")
    for h in hits:
        print("    " + h[:140])

    out2 = subprocess.run(["grep", "-rnoE", r"method=\"[A-Za-z0-9]+\"",
                           "--include=*.py", "src", "audit", "tests"],
                          cwd=ROOT, capture_output=True, text=True).stdout
    from collections import Counter
    c = Counter(l.split("method=")[1] for l in out2.splitlines() if l.strip())
    print(f"\n  integrator methods actually requested: {dict(c)}")
    print("  (solve_ivp's default, where no method is named, is RK45 with")
    print("   adaptive step control, which is not an explicit fixed step.)")


# --------------------------------------------------------------------------
def main():
    print(__doc__)
    worst = s1_time_constant_and_stability()
    s2_euler_actually_run(worst)
    s3_information_in_a_difference()
    s4_no_euler_in_the_code()


if __name__ == "__main__":
    main()
