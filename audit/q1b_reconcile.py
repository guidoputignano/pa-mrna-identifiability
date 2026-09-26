"""
S1 reconciled against the 25 September robustness probes.

`audit/q1_increment_shape.py` reported theta_hat/theta = 5.64 [2.44, 8.61] at a
2.35 sigma residual for the Attarwala two-compartment truth. The probe
`audit/probes_2026-09-25/p3_shape.py` gives 1.37 to 1.40 at 0.33 to 0.37 sigma
on noiseless records. On the exponential truths the two agree, so the question
is whether they model the same two-compartment object, and if they do not,
which difference moves the estimate. This script answers that on IDENTICAL
noiseless records, fitted by BOTH estimators, so that an estimator error and a
modelling difference cannot be confused.

Three things are crossed.

    KERNEL   probe  impulse straight into PCC, the 2x2 protein block only.
             q1     the same block driven through Attarwala's effect
                    compartment, ke0 = 0.242 /h, so the protein RISES from zero
                    over a few hours after each infusion.
    DRAW     as run     a draw at the instant of an infusion is read on the
                        POST-dose side, which is what `design.run.delta` does
                        (`t >= t_k`) and so what the committed analyst does.
             pre-dose   a draw at a dose instant is taken BEFORE the infusion,
                        for the truth and the analyst alike. Clinical practice,
                        and the only reading under which a quasi-steady model
                        can be asked about that draw at all.
             drop       every draw at a dose instant removed.
    FITTER   probe  `probes_2026-09-25/fit.py`, and its p3b and p3c variants.
             q1     `q1_increment_shape.Analyst`, the global K_m scan with the
                    exact polish.

The monthly schedule `np.arange(0, 730, 30.4)` puts exactly one draw at a dose
instant in two years, t = 0, because 30.4 k is a multiple of 14 only at
k = 35. The probe's 14.5 d schedule, `np.arange(0.5, 730, 14.5)`, puts one at
t = 392, where 0.5 + 14.5 k is a multiple of 14 at k = 27; an earlier draft of
this docstring said it had none. The 2-daily schedule sits on odd days and
never meets a dose. `dose_instants` finds them rather than assuming.

Footing: the probe's. theta = 3000, noiseless, load CV 10%, the probe's own 8
load records (seeds 0 to 7), the paper's two priors. The amplitude of every
truth is matched on the steady-state peak of the nominal train, which is both
scripts' headline convention; section R4 varies it.

Nothing here records an expected result in advance.

Run:  python audit/q1b_reconcile.py
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "src"))
sys.path.insert(0, str(HERE / "probes_2026-09-25"))

import numpy as np
from scipy.optimize import least_squares

import q1_increment_shape as Q
import core as C
import fit as PF

SEEDS = range(8)
EPS = 1e-9
SCHED = {"monthly": np.arange(0.0, 730.0, 30.4),
         "14.5 d": np.arange(0.5, 730.0, 14.5),
         "2-daily": np.arange(1.0, 730.0, 2.0)}
FINE = np.arange(0.0, 730.0, 0.01)
SS = (FINE > 600) & (FINE < 614)
A_PK = C.train(FINE)[SS].max()


def banner(t):
    print("\n" + "=" * 92)
    print("  " + t)
    print("=" * 92)


# --------------------------------------------------------------------------
#  The truths, as functions of the draw times. `pre` evaluates just before
#  each time, so a draw at a dose instant sees the previous doses only.
# --------------------------------------------------------------------------
REF = Q.committed_ref()
TN, WN = Q.nominal_train()


def q1_truth(name, conv="peak"):
    k = Q.KERNELS[name]
    A = Q.amplitude(k, conv, REF, key=("recon", name))
    return lambda t: A * Q.assemble(k, TN, WN, t)


def probe_truth(name):
    if name == "biexp":
        b = C.train_biexp(FINE)[SS].max()
        return lambda t: C.train_biexp(t) * (A_PK / b)
    half = {"exp 5.0 d": 5.0, "exp 7.0 d": 7.0, "exp 3.5 d": 3.5}[name]
    pk = C.train(FINE, half=half)[SS].max()
    return lambda t: C.train(t, half=half) * (A_PK / pk)


def dose_instants(times, period=14.0, tol=1e-6):
    """Which draws fall at an instant of the nominal dose train."""
    ph = np.mod(times, period)
    return (ph < tol) | (period - ph < tol)


def draws(times, conv):
    """(times the truth and analyst are evaluated at, times kept)."""
    if conv == "as run":
        return times, times
    if conv == "pre-dose":
        return times - EPS, times
    keep = times[~dose_instants(times)]
    return keep, keep


def records(truth, times, conv):
    """The probe's noiseless records, with the truth read at the stated side."""
    t_eval, t_keep = draws(times, conv)
    out = []
    for s in SEEDS:
        m, y, sig, _ = PF.record(C.TH_NOM, t_keep, s, inc_true=truth(t_eval),
                                 noisy=False)
        out.append((m, y, sig, t_eval))
    return out


def summ(v):
    v = np.asarray(v, float)
    return f"{np.median(v):6.2f} [{v.min():5.2f},{v.max():5.2f}]"


# --------------------------------------------------------------------------
def part_r1():
    banner("R1. WHAT EACH SCRIPT'S TWO-COMPARTMENT KERNEL IS")
    k = Q.kernel_attarwala()
    A = np.array([[-(0.00474 + 0.0075), 0.00474], [0.00474, -0.00474]]) * 24.0
    w = np.linalg.eigvals(A)
    print(f"  q1     rates per day {np.round(k.w, 6)}")
    print(f"         half-lives    {np.round(np.log(2) / -k.w, 4)} d")
    print(f"         coefficients  {np.round(k.c, 4)}")
    print(f"  probe  half-lives    {np.round(np.sort(np.log(2) / -w), 4)} d")
    print("""
  Both take kdeg = 0.0075 /h and kq = 0.00474 /h from Attarwala Table 1 and
  multiply the rates by 24, so both run in DAYS. Both share the two protein
  half-lives, 1.99 d and 11.81 d, which PHASE5_ROBUSTNESS.md finding 3 gives.

  They differ in where the dose enters. The probe puts a unit impulse straight
  into PCC, so the protein JUMPS at the infusion. q1 puts it into the effect
  compartment, Ce(0) = 1 with ke0 = 0.242 /h and Slope = 55.6, so there is a
  third mode, 0.119 d, and PCC RISES from zero over a few hours. Neither has
  mRNA PK beyond that compartment.

  Where each amplitude convention is matched:
    q1     dose  the SINGLE-dose peak, i.e. the first cycle, of one impulse.
           peak  the steady-state peak of the infinite nominal train, in
                 closed form (ExpSum.steady_peak).
           mean  the steady-state cycle mean, which is the impulse area over
                 the interval, also closed form.
    probe  peak and mean on days 600 to 614 of a 730 d train from t = 0, which
           is steady state to within exp(-600 ln2 / 11.81), so the same
           object as q1's.
  So `peak` and `mean` are matched at steady state in both, and q1's `dose` is
  the only first-cycle convention in either script.""")
    print("\n  The truths at the first eight monthly draws, peak-matched:")
    t = SCHED["monthly"][:8]
    an = C.train(t)
    qa = q1_truth("Attarwala 2-comp")(t)
    pb = probe_truth("biexp")(t)
    print(f"    {'t (d)':>7s} {'phase':>6s} {'analyst':>8s} {'q1 truth':>9s} "
          f"{'probe truth':>12s}")
    for i in range(t.size):
        print(f"    {t[i]:7.1f} {t[i] % 14:6.2f} {an[i]:8.0f} {qa[i]:9.0f} "
              f"{pb[i]:12.0f}")


# --------------------------------------------------------------------------
def fit_fixed(recs, kappa, est):
    """The committed 3.5 d shape. kappa None: amplitude fixed; else its prior."""
    A = Q.make_analysts((TN, WN), REF, "peak",
                        which=("fixed 3.5 d", "fixed + free amp"))
    an = A["fixed + free amp" if kappa else "fixed 3.5 d"]
    th, rms = [], []
    for m, y, sig, t_eval in recs:
        if est == "probe":
            r = PF.fit(m, C.train(t_eval), y, sig, C.L_A, kappa_cv=kappa)
            th.append(r["theta"] / C.TH_NOM)
            rms.append(r["rms"])
        else:
            f = an.fit(m, y, times=t_eval)
            th.append(f["theta"] / C.TH_NOM)
            rms.append(Q.rms_residual(f["theta"], f["Km"], f["ke"], f["b"],
                                      f["kappa"], m, y, f["phi"]))
    return th, rms


def part_r2():
    banner("R2. THE FIXED 3.5 d ANALYST, BOTH KERNELS, BOTH FITTERS")
    print("  theta_hat/theta median [min, max] over the probe's 8 noiseless "
          "records, and the")
    print("  RMS residual in sigma. Identical records go to both fitters.\n")
    print(f"  {'kernel':8s} {'draw':10s} {'schedule':9s} {'amplitude':10s} "
          f"{'probe fitter':>28s} {'q1 fitter':>28s}")
    out = {}
    for kname, tr in (("probe", probe_truth("biexp")),
                      ("q1", q1_truth("Attarwala 2-comp"))):
        for conv in ("as run", "pre-dose", "drop"):
            for sname, times in SCHED.items():
                if conv != "as run" and not dose_instants(times).any():
                    continue          # no draw at a dose instant to move
                recs = records(tr, times, conv)
                for kap in (None, 0.12):
                    cells = []
                    for est in ("probe", "q1"):
                        th, rms = fit_fixed(recs, kap, est)
                        out[(kname, conv, sname, kap, est)] = (th, rms)
                        cells.append(f"{summ(th)} {np.median(rms):5.2f}s")
                    print(f"  {kname:8s} {conv:10s} {sname:9s} "
                          f"{'fixed' if kap is None else '12% prior':10s} "
                          f"{cells[0]:>28s} {cells[1]:>28s}")
    worst = max(abs(np.median(out[k[:4] + ("probe",)][0])
                    - np.median(out[k[:4] + ("q1",)][0]))
                for k in out if k[4] == "q1")
    print(f"\n  worst gap between the two fitters' medians, every row: "
          f"{worst:.2e}")
    return out


# --------------------------------------------------------------------------
def fit_hl_probe(m, t, y, sig, kappa_free):
    """p3b (kappa fixed at D = 6000) and p3c (kappa free, flat), verbatim forms."""
    bsd = C.LOAD_FRAC * C.L_A

    def res(p):
        lt, lk, b, lke, lhl = p[:5]
        ka = np.exp(p[5]) if kappa_free else 1.0
        th, K, ke, hl = np.exp([lt, lk, lke, lhl])
        pred = (th + ka * C.train(t, half=hl)) * m / (K + m) + ke * m + b
        return np.concatenate([(y - pred) / sig, [b / bsd],
                               [(lke - np.log(C.KE)) / 0.8]])
    best = None
    for dt in (0, -1, 1):
        for dh in (0, 0.5, -0.5, 1.0):
            p0 = [np.log(3000) + dt, np.log(C.KM), 0.0, np.log(C.KE),
                  np.log(3.5) + dh] + ([0.0] if kappa_free else [])
            r = least_squares(res, p0, x_scale="jac", max_nfev=6000,
                              xtol=1e-12, ftol=1e-12, gtol=1e-12)
            if best is None or r.cost < best.cost:
                best = r
    rms = float(np.sqrt(np.mean(best.fun[:m.size] ** 2)))
    return float(np.exp(best.x[0])) / C.TH_NOM, float(np.exp(best.x[4])), rms


def part_r3():
    banner("R3. THE WRONG-FAMILY CELLS")
    print("""  A single-exponential analyst with its half-life free, on the
  two-compartment truth. Three analysts, because the scripts did not free the
  same thing:
    p3b / dose   half-life free, dose amplitude FIXED at D = 6000 per dose.
                 This is probe p3b, and it is q1's `dose` convention.
    q1 / peak    half-life free, amplitude re-normalised to the committed
                 steady PEAK for every candidate half-life. This is what q1's
                 Part 1 ran, and it is NOT the probe's analyst.
    p3c          half-life and amplitude both free, flat. Probe p3c.
  Monthly, 8 noiseless records, both kernels, draw as run and pre-dose.""")
    out = {}
    print(f"\n  {'kernel':8s} {'draw':10s} {'analyst':12s} {'fitter':7s} "
          f"{'theta_hat/theta':>22s} {'half-life':>8s} {'rms':>6s}")
    for kname, tr in (("probe", probe_truth("biexp")),
                      ("q1", q1_truth("Attarwala 2-comp"))):
        for conv in ("as run", "pre-dose"):
            recs = records(tr, SCHED["monthly"], conv)
            for lab, conv_a, kfree in (("p3b / dose", "dose", False),
                                       ("q1 / peak", "peak", False),
                                       ("p3c", "dose", True)):
                a = Q.Analyst(lab, (TN, WN), "exp", free_kappa=kfree,
                              kappa_prior=None)
                a.set_reference(REF, conv_a)
                rows = {"q1": [], "probe": []}
                for m, y, sig, t_eval in recs:
                    f = a.fit(m, y, times=t_eval)
                    rows["q1"].append((f["theta"] / C.TH_NOM,
                                       float(np.exp(f["q"][0])),
                                       Q.rms_residual(f["theta"], f["Km"],
                                                      f["ke"], f["b"],
                                                      f["kappa"], m, y,
                                                      f["phi"])))
                    if conv_a == "dose":
                        rows["probe"].append(fit_hl_probe(m, t_eval, y, sig,
                                                          kfree))
                for est, v in rows.items():
                    if not v:
                        continue
                    v = np.array(v)
                    out[(kname, conv, lab, est)] = v
                    print(f"  {kname:8s} {conv:10s} {lab:12s} {est:7s} "
                          f"{summ(v[:, 0]):>22s} {np.median(v[:, 1]):8.2f} "
                          f"{np.median(v[:, 2]):6.2f}")
    return out


# --------------------------------------------------------------------------
def part_r4():
    banner("R4. DOES THE DOSE-INSTANT DRAW MOVE THE EXPONENTIAL TRUTHS TOO?")
    print("  The exponential truths JUMP at the infusion, as the analyst does, so "
          "the t = 0")
    print("  draw should matter much less there. Checked rather than assumed, "
          "q1 fitter,")
    print("  fixed 3.5 d analyst, monthly, 8 noiseless records, both "
          "amplitude conventions.\n")
    fa = Q.make_analysts((TN, WN), REF, "peak",
                         which=("fixed 3.5 d",))["fixed 3.5 d"]
    print(f"  {'truth':18s} {'convention':10s} "
          + "".join(f"{c:>24s}" for c in ("as run", "pre-dose", "drop")))
    for name in ("exp 2.5 d", "exp 5.0 d", "exp 7.0 d", "Attarwala 2-comp"):
        for tc in ("peak", "dose", "mean"):
            cells = []
            for conv in ("as run", "pre-dose", "drop"):
                recs = records(q1_truth(name, tc), SCHED["monthly"], conv)
                th = [fa.fit(m, y, times=t)["theta"] / C.TH_NOM
                      for m, y, s, t in recs]
                cells.append(summ(th))
            print(f"  {name:18s} {tc:10s}" + "".join(f"{c:>24s}" for c in cells))


def part_r6():
    banner("R6. DOES S1.8's FAMILY DISCRIMINATION SURVIVE THE PROBE'S KERNEL?")
    print("""  q1 Part 5 (c) found that telling the families apart wants a draw in
  the hours after an infusion, because there one shape jumps and the other
  rises. That rise IS the effect compartment, so the claim is re-run here on
  the probe's jump kernel, same method: both shapes given their best free
  amplitude, the best opposed pair of phases, a 3-sigma call, peak-matched.""")
    ref = Q.committed_ref()
    A = np.array([[-(0.00474 + 0.0075), 0.00474], [0.00474, -0.00474]]) * 24.0
    w, V = np.linalg.eig(A)
    c = np.linalg.solve(V, np.array([1.0, 0.0]))
    jump = Q.ExpSum(np.real(V[0] * c), np.real(w))
    kern = {"effect compartment (q1)": Q.kernel_attarwala(),
            "jump (probe)": jump}
    cvs = (0.05, 0.10, 0.20, 0.30, 0.50)
    for readout in ("delta", "capacity"):
        print(f"\n  readout = {readout}")
        print(f"    {'kernel':25s} {'earliest':>9s} {'best pair (d)':>15s} "
              f"{'sep/CV':>8s}" + "".join(f"{f'N@{int(100*v)}%':>8s}"
                                            for v in cvs))
        for kn, k in kern.items():
            for lo in (0.04, 0.25, 1.0, 2.0):
                g = np.round(np.arange(lo, Q.T_DOSE - lo + 1e-9, 0.02), 3)

                def ly(kk):
                    d = Q.steady_profile(kk, ref, "peak", g)
                    return np.log(d if readout == "delta" else Q.TH + d)
                d = ly(k) - ly(Q.kernel_exp(Q.TH_HALF0))
                j, i = int(np.argmax(d)), int(np.argmin(d))
                sep = float(abs(d[j] - d[i]) / np.sqrt(2.0))
                n = "".join(f"{int(2 * max(1, np.ceil((3 * v / sep) ** 2))):8d}"
                            for v in cvs)
                print(f"    {kn:25s} {lo:9.2f} "
                      f"{f'{g[j]:.2f}/{g[i]:.2f}':>15s} {sep:8.4f}" + n)


def part_r5():
    banner("R5. A FLAG, NOT A FINDING: THE SAME DRAW IN THE PAPER'S OWN HEADLINE")
    print("""  The committed design record is np.arange(0, 730, 30.4), and every
  design number in the paper, 0.1697 monthly among them, is computed on it with
  design.run.delta, which reads the t = 0 draw POST-dose. That draw is inside
  the marker's settling window, so the quasi-steady model the Fisher bound is
  taken from does not describe it. Here is what it is worth to the bound,
  through g_design_common's own sd_over_seeds, 24 load records at CV 10%:""")
    import g_design_common as G
    t = SCHED["monthly"]
    for lab, tt in (("as committed (t = 0 read post-dose)", t),
                    ("t = 0 read pre-dose (no increment yet)", None),
                    ("t = 0 dropped, n = 24", t[1:])):
        if tt is None:
            vals = []
            for s in range(24):
                m, inc, _ = G.record(t, 0.10, seed=s)
                inc = inc.copy()
                inc[0] = 0.0
                Ls = G.record(t, 0.10, seed=s)[2]
                B = G.KE * G.KM + G.TH + inc - Ls
                m = (-B + np.sqrt(B ** 2 + 4.0 * G.KE * Ls * G.KM)) / (2 * G.KE)
                vals.append(G.sd_ln_theta(m, inc))
            v = np.array(vals)
            lo, md, hi = v.min(), np.median(v), v.max()
        else:
            lo, md, hi = G.sd_over_seeds(tt, 0.10, seeds=range(24))
        print(f"    {lab:42s} sd(ln theta) {md:.4f} [{lo:.4f}, {hi:.4f}]")
    print("""  Reported so the author can decide whether it matters. Not acted on:
  the paper is frozen, and S1's question is misspecification, not the
  headline.""")


def main():
    print("=" * 92)
    print("  S1 RECONCILED AGAINST THE 25 SEPTEMBER PROBES")
    print("=" * 92)
    #  this script moves the dose-day draw itself, through `draws`, so q1's
    #  own convention is pinned to the post-dose reading for the duration and
    #  restored afterwards. Left at q1's default the analyst would be read
    #  pre-dose on the "as run" rows while the truth was not.
    with Q.rebind(Q, DOSE_DAY_DRAW="post-dose"):
        part_r1()
        part_r2()
        part_r3()
        part_r4()
        part_r5()
        part_r6()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
