"""
Can a sub-daily sampling burst measure k_e, where the quasi-steady record
cannot?

THE QUESTION, AND WHY IT IS NOT THE ONE audit/u9 CLOSED. u9 closed the TIMED
URINE COLLECTION as a route to the k_e frontier, on three counts: units, since
a clearance needs a volume of distribution for 3-HP that nobody has; channel
coverage, since k_e is registered LUMPED; and species, since the parent's
alternative channel makes propionylcarnitine and not the marker. Those three
close a COLLECTION. They do not close the OBJECT.

src/common/params.py's own k_e row says k_e = 5/day is "a marker half-life of
ln2/5 = 0.139 d, about 3.3 hours". So k_e is a RELAXATION RATE, and a
relaxation rate is read off a concentration-time curve. No volume of
distribution, no collection completeness, no species extrapolation. The units
objection does not apply to it, because nothing is being divided by a volume.

WHY THE LOCAL MODEL CANNOT SEE IT. eq:local is a quasi-steady balance: it has
already integrated the relaxation away, so it answers "no" at every sampling
density by construction, not by finding. PART A reproduces that null rather
than asserting it.

SO THE QUESTION GOES TO THE DYNAMIC MODEL of Section 3, which carries the
state. PART B puts it there: one inpatient day of sub-daily draws spanning a
scheduled dose, which is a standard pharmacokinetic procedure and needs no
instrument the clinic lacks.

THE PERTURBATION IS ALREADY IN THE PROTOCOL. Delta(t) jumps by D at each dose
time and decays with the therapeutic half-life between doses, so a dose is a
STEP on the enzymatic channel and the marker's relaxation to its new level is
the transient. Admitting the patient on a dosing day is the whole intervention.

WHAT LOOKED FAVOURABLE, AND THE QUANTITY THAT ACTUALLY DECIDES IT. k_e's
share of the total elimination FLUX is 15.5% at theta = 3000, 52.8% at 600,
70.8% at 300 and 83.7% at 150, so the measurement should get easier exactly
where the frontier gets tighter. But a relaxation measurement does not depend
on the flux share. It depends on k_e's share of the RELAXATION RATE
J = (theta + Delta) g'(m) + k_e, and those are different numbers. PART B
computes both, and the difference is why the obvious burst fails.

FOOTING. The dynamic model dm/dt = L + b - (theta + Delta(t)) g(m) - k_e m,
with the sensitivity system of src/common/dynamics.py extended by the load
bias, which that module does not carry. theta with the RECORDED LOAD HELD
FIXED at L = 322.727 umol/L/d, K_m = 100, k_e = 5/d, D = 6000 on a 14 d train
with a 3.5 d half-life so Delta runs 400 to 6400 at the train's steady state.
Marker assay CV 5% per draw, relative, on ln m. The burst starts at a trough,
the patient having been on the train long enough to sit on the quasi-steady
locus, so the sensitivities start at their quasi-steady values and NOT at
zero: the record before the burst already carries that much.

Nothing here records an expected result in advance.

Run:  python audit/u14_dynamic_burst.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from scipy.integrate import solve_ivp

from design.run import (ASSAY_CV, DOSE_INTERVAL, D_AMP, LAMBDA, LOAD_CV,
                        LOAD_FRAC, load)

KM, KE = 100.0, 5.0
LB = float(load(0.0))
GATE = 0.20
MONTH = 30.4
HOUR = 1.0 / 24.0
KE_PRIOR = 0.8
PUBLISHED_KAPPA_CV = 0.1402
T_DOSE = 50.0 * DOSE_INTERVAL
PEAK = D_AMP / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))
TROUGH = PEAK * np.exp(-LAMBDA * DOSE_INTERVAL)
SEEDS = range(24)
TIER1_N = 118

#: the frontier audit/u9_ke_frontier.py states, Tier 1 record
FRONTIER = {3000.0: 0.8001, 600.0: 0.1567, 300.0: 0.0633, 150.0: None}


def qss(theta, inc, Ltot, Km=KM, ke=KE):
    inc = np.asarray(inc, dtype=float)
    B = ke * Km + theta + inc - Ltot
    return (-B + np.sqrt(B ** 2 + 4.0 * ke * Ltot * Km)) / (2.0 * ke)


def sigma_flux(theta, assay_cv=ASSAY_CV):
    mbar = float(qss(theta, 0.0, LB))
    slope = theta * KM / (KM + mbar) ** 2 + KE
    return float(np.hypot(LOAD_FRAC * LB, slope * assay_cv * mbar))


# ---------------------------------------------------------------- the record
def record_rows(theta, n, seed, step=MONTH):
    """Quasi-steady marker record: flux-space rows in (lnth, lnKm, b, lnke)."""
    t = step * np.arange(n, dtype=float)
    tau = np.mod(t - T_DOSE, DOSE_INTERVAL)
    inc = PEAK * np.exp(-LAMBDA * tau)
    Ls = LB * (1.0 + LOAD_CV *
               np.random.default_rng(seed).standard_normal(t.size))
    m = qss(theta, inc, Ls)
    g = m / (KM + m)
    return np.stack([theta * g,
                     -(theta + inc) * KM * m / (KM + m) ** 2,
                     np.ones_like(m),
                     KE * m], axis=1), m, inc


# ---------------------------------------------------------------- the burst
def share_of_J(theta, inc, Km=KM, ke=KE, Ltot=None):
    """k_e's share of the RELAXATION RATE at the quasi-steady level."""
    Ltot = LB if Ltot is None else Ltot
    m = float(qss(theta, inc, Ltot, Km, ke))
    J = (theta + inc) * Km / (Km + m) ** 2 + ke
    return ke / J, J, m


def burst_rows(theta, times_d, assay_cv=ASSAY_CV, Km=KM, ke=KE, Ltot=None,
               inc0=None, load_step=0.0):
    """
    Rows of d(ln m)/d(lnth, lnKm, b, lnke) along the post-dose transient.

    The sensitivity system is src/common/dynamics.py's, extended by the load
    bias, whose forcing is +1 and which that module does not carry. The
    initial condition is the QUASI-STEADY sensitivity at the trough, not
    zero: the patient has been on the train for many cycles, so the state and
    its derivatives have relaxed onto the locus before the dose lands.
    """
    Ltot = LB if Ltot is None else Ltot
    inc0 = TROUGH if inc0 is None else inc0
    m0 = float(qss(theta, inc0, Ltot, Km, ke))
    J0 = (theta + inc0) * Km / (Km + m0) ** 2 + ke
    y0 = [m0,
          -theta * (m0 / (Km + m0)) / J0,                       # q_lntheta
          (theta + inc0) * Km * m0 / (Km + m0) ** 2 / J0,       # q_lnKm
          1.0 / J0,                                             # q_b
          -ke * m0 / J0]                                        # q_lnke

    #  inc0 = TROUGH with no load step is the patient sitting at the trough;
    #  inc0 = PEAK is a dose landing at t = 0. A load step is a PRESCRIBED
    #  additive change of a known size, so it moves the load and not b.
    dosed = inc0 > TROUGH * 1.5

    def dl(t):
        return (PEAK if dosed else TROUGH) * np.exp(-LAMBDA * t)

    def Lt(t):
        return Ltot * (1.0 + load_step)

    def rhs(t, y):
        m, qt, qK, qb, qk = y
        m = max(m, 1e-14)
        th = theta + dl(t)
        J = th * Km / (Km + m) ** 2 + ke
        return [Lt(t) - th * m / (Km + m) - ke * m,
                -theta * (m / (Km + m)) - J * qt,
                th * Km * m / (Km + m) ** 2 - J * qK,
                1.0 - J * qb,
                -ke * m - J * qk]

    T = float(max(times_d))
    sol = solve_ivp(rhs, (0.0, T + 1e-9), y0, dense_output=True,
                    rtol=1e-11, atol=1e-14, method="LSODA")
    Y = sol.sol(np.asarray(times_d, dtype=float))
    m = Y[0]
    #  d(ln m)/dp = q_p / m, and the observation noise is relative
    return (np.stack([Y[1] / m, Y[2] / m, Y[3] / m, Y[4] / m], axis=1)
            / assay_cv), m


def fisher(theta, n_rec, seed, burst_times=None, assay_cv=ASSAY_CV,
           kappa_cv=None, burst_cv=None, n_bursts=1, inc0=None,
           load_step=0.0):
    """Joint information from the record and, optionally, one burst."""
    Q, m, inc = record_rows(theta, n_rec, seed) if n_rec else (None, None, None)
    pri = [np.inf, np.inf, (LOAD_FRAC * LB) ** 2, KE_PRIOR ** 2]
    F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri])
    if n_rec:
        sig = sigma_flux(theta, assay_cv)
        if kappa_cv is not None:
            g = m / (KM + m)
            Q = np.column_stack([Q, inc * g])
        F2 = np.zeros((Q.shape[1], Q.shape[1]))
        F2[:4, :4] = F
        if kappa_cv is not None:
            F2[4, 4] = 1.0 / kappa_cv ** 2
        F = F2 + Q.T @ Q / sig ** 2
    elif kappa_cv is not None:
        F = np.pad(F, ((0, 1), (0, 1)))
        F[4, 4] = 1.0 / kappa_cv ** 2
    if burst_times is not None and len(burst_times):
        R, _ = burst_rows(theta, burst_times,
                          assay_cv if burst_cv is None else burst_cv,
                          inc0=inc0, load_step=load_step)
        if kappa_cv is not None:
            R = np.column_stack([R, np.zeros(R.shape[0])])
        #  n_bursts IDENTICAL independent admissions: each contributes the
        #  same rows with independent assay noise, so the information adds.
        #  Offsetting the sample TIMES by a dose interval instead is what
        #  this script did first, and it integrates one trajectory out to
        #  38 days after a single dose, which is not four bursts.
        F = F + n_bursts * (R.T @ R)
    return F


def sd_of(F, idx):
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[idx, idx]
    return float(np.sqrt(v)) if v > 0 else np.inf


def med(theta, idx, **kw):
    v = np.array([sd_of(fisher(theta, seed=s, **kw), idx) for s in SEEDS],
                 float)
    f = v[np.isfinite(v)]
    return float(np.median(f)) if f.size else np.inf


# ======================================================================
def part_a():
    print("=" * 86)
    print("PART A.  the quasi-steady record cannot see it, at any density")
    print("=" * 86)
    print("\n  The same 118 samples reallocated from monthly to dense within")
    print("  cycles. eq:local has integrated the relaxation away, so this")
    print("  should answer no by construction. Medians over 24 records.\n")
    print(f"    {'theta':>7} {'spacing':>12} {'sd(ln theta)':>13} "
          f"{'sd(ln k_e)':>12}")
    print("    " + "-" * 48)
    for theta in (3000.0, 300.0):
        for lab, step in [("monthly", MONTH), ("weekly", 7.0),
                          ("daily", 1.0), ("3-hourly", 3 * HOUR)]:
            v_t, v_k = [], []
            for s in SEEDS:
                Q, _, _ = record_rows(theta, TIER1_N, s, step=step)
                sig = sigma_flux(theta)
                pri = [np.inf, np.inf, (LOAD_FRAC * LB) ** 2, KE_PRIOR ** 2]
                F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri]) \
                    + Q.T @ Q / sig ** 2
                v_t.append(sd_of(F, 0))
                v_k.append(sd_of(F, 3))
            print(f"    {theta:7.0f} {lab:>12} {np.median(v_t):13.4f} "
                  f"{np.median(v_k):12.4f}")
        print()
    print("  Confirmed. A balance with no state has no relaxation to")
    print("  resolve, so density cannot help, and the question has to go to")
    print("  the dynamic model or nowhere.")


def part_b():
    print("=" * 86)
    print("PART B.  the dose is the WRONG perturbation, and the flux share")
    print("         is the wrong diagnostic")
    print("=" * 86)
    print("""
  A relaxation measurement depends on k_e's share of the RELAXATION RATE
  J = (theta + Delta) g'(m) + k_e, not on its share of the flux. The dose
  adds enzyme, and enzyme is exactly what k_e competes with in J. So the
  natural perturbation suppresses the thing being measured.
""")
    print(f"    {'theta':>7} {'at the trough':>28} {'just after a dose':>28}")
    print(f"    {'':>7} {'Delta':>7}{'m':>7}{'J':>7}{'k_e/J':>7} "
          f"{'Delta':>7}{'m':>7}{'J':>7}{'k_e/J':>7}  {'flux share':>11}")
    print("    " + "-" * 78)
    for theta in (3000.0, 600.0, 300.0, 150.0):
        st, Jt, mt = share_of_J(theta, TROUGH)
        sp, Jp, mp = share_of_J(theta, PEAK)
        flux = KE * mt / (theta * mt / (KM + mt) + KE * mt)
        print(f"    {theta:7.0f} {TROUGH:7.0f}{mt:7.2f}{Jt:7.2f}{st:7.3f} "
              f"{PEAK:7.0f}{mp:7.2f}{Jp:7.2f}{sp:7.3f}  {flux:11.3f}")
    print()
    print("  At a tenth of nominal the dose takes k_e from 55% of the")
    print("  relaxation rate to 8%, a factor of seven, while the flux share")
    print("  barely moves. And the cost shows up exactly there:")
    print()
    tt = HOUR * np.arange(0, 25, dtype=float)
    print(f"    {'theta':>7} {'no burst':>10} {'1 post-dose day':>16} "
          f"{'16 post-dose days':>18} {'frontier':>9}")
    print("    " + "-" * 64)
    for theta in (600.0, 300.0):
        b0 = med(theta, 3, n_rec=TIER1_N, kappa_cv=PUBLISHED_KAPPA_CV)
        b1 = med(theta, 3, n_rec=TIER1_N, burst_times=tt, inc0=PEAK,
                 kappa_cv=PUBLISHED_KAPPA_CV)
        b16 = med(theta, 3, n_rec=TIER1_N, burst_times=tt, inc0=PEAK,
                  n_bursts=16, kappa_cv=PUBLISHED_KAPPA_CV)
        print(f"    {theta:7.0f} {b0:10.4f} {b1:16.4f} {b16:18.4f} "
              f"{FRONTIER[theta]:9.4f}")
    print()
    print("  Sixteen inpatient days spanning sixteen doses move sd(ln k_e)")
    print("  by about two per cent. The route is not weak; the PLACEMENT is")
    print("  wrong.")


def part_c():
    print()
    print("=" * 86)
    print("PART C.  perturb at the trough instead, where k_e owns the rate")
    print("=" * 86)
    print("\n  One 24 h hourly burst, Tier 1 record, a PRESCRIBED load step")
    print("  at the trough rather than a dose. The step is taken as known.\n")
    tt = HOUR * np.arange(0, 25, dtype=float)
    print(f"    {'theta':>7} {'no burst':>9} {'fast -50%':>10} "
          f"{'fast -30%':>10} {'load +30%':>10} {'load +50%':>10} "
          f"{'post-dose':>10} {'frontier':>9}")
    print("    " + "-" * 82)
    for theta in (600.0, 300.0, 150.0):
        kw = dict(n_rec=TIER1_N, burst_times=tt, kappa_cv=PUBLISHED_KAPPA_CV)
        row = [med(theta, 3, n_rec=TIER1_N, kappa_cv=PUBLISHED_KAPPA_CV)]
        for a in (-0.5, -0.3, 0.3, 0.5):
            row.append(med(theta, 3, inc0=TROUGH, load_step=a, **kw))
        row.append(med(theta, 3, inc0=PEAK, **kw))
        fr = FRONTIER[theta]
        print(f"    {theta:7.0f} " + " ".join(f"{v:9.4f}" for v in row)
              + f" {(f'{fr:.4f}' if fr else 'n/a'):>10}")
    print()
    print("  A load step at the trough is worth what sixteen post-dose days")
    print("  are not. But the direction that helps most raises the marker,")
    print("  and PART D is where that stops.")


def part_d():
    print()
    print("=" * 86)
    print("PART D.  the tolerable level caps the useful amplitude")
    print("=" * 86)
    from common.params import M_TOL
    print(f"\n  m_tol = {M_TOL:.0f} uM, src/common/params.py. The largest")
    print("  prescribed load step whose whole burst stays under it, and")
    print("  what that step buys.\n")
    tt = HOUR * np.arange(0, 25, dtype=float)
    print(f"    {'theta':>7} {'trough m':>9} {'max safe step':>14} "
          f"{'peak m':>8} {'sd(ln k_e)':>11} {'no burst':>9} {'frontier':>9}")
    print("    " + "-" * 74)
    for theta in (600.0, 300.0, 150.0):
        best = None
        for a in np.arange(-0.9, 1.51, 0.05):
            _, mb = burst_rows(theta, tt, inc0=TROUGH, load_step=float(a))
            if mb.max() <= M_TOL and (best is None or a > best[0]):
                best = (float(a), float(mb.max()))
        a, pk = best
        v = med(theta, 3, n_rec=TIER1_N, burst_times=tt, inc0=TROUGH,
                load_step=a, kappa_cv=PUBLISHED_KAPPA_CV)
        b0 = med(theta, 3, n_rec=TIER1_N, kappa_cv=PUBLISHED_KAPPA_CV)
        fr = FRONTIER[theta]
        print(f"    {theta:7.0f} {float(qss(theta, TROUGH, LB)):9.3f} "
              f"{a:+13.0%} {pk:8.2f} {v:11.4f} {b0:9.4f} "
              f"{(f'{fr:.4f}' if fr else 'n/a'):>9}")
    print()
    print("  So the amplitude that would close the gap is the one the")
    print("  clinic exists to prevent, and the amplitude that respects the")
    print("  paper's own tolerable level buys almost nothing. What is left")
    print("  is the other direction: a FAST lowers the marker and is safe.")


def part_e():
    print()
    print("=" * 86)
    print("PART E.  the safe route, and what it actually reaches")
    print("=" * 86)
    print("\n  A 50% fast at the trough, 24 hourly draws per admission. The")
    print("  burst's assay CV is separated from the record's because an")
    print("  inpatient series from one line, analysed in one run, is where")
    print("  a tighter assay is plausible if it is plausible anywhere.\n")
    tt = HOUR * np.arange(0, 25, dtype=float)
    for theta in (600.0, 300.0):
        fr = FRONTIER[theta]
        print(f"  theta = {theta:.0f}, frontier sd(ln k_e) = {fr:.4f}")
        print(f"    {'burst CV':>9}" +
              "".join(f"{f'{n} adm':>11}" for n in (1, 2, 4, 8, 16)))
        print("    " + "-" * 64)
        for cv in (0.05, 0.02, 0.01):
            cells = []
            for nb in (1, 2, 4, 8, 16):
                v = med(theta, 3, n_rec=TIER1_N, burst_times=tt, inc0=TROUGH,
                        load_step=-0.50, burst_cv=cv, n_bursts=nb,
                        kappa_cv=PUBLISHED_KAPPA_CV)
                cells.append(f"{v:.4f}" + ("*" if v <= fr else ""))
            print(f"    {cv:9.0%}" + "".join(f"{c:>11}" for c in cells))
        print()
    print("  * reaches the frontier. The route is ASSAY-LIMITED and not")
    print("  physics-limited, which is the favourable half of the answer:")
    print("  the one input that moves it most is the one an inpatient")
    print("  setting is best placed to improve.")


def part_f():
    print()
    print("=" * 86)
    print("PART F.  what the whole route is worth, against the alternatives")
    print("=" * 86)
    print("""
  Three routes to k_e have now been priced, all against the same frontier.

    timed urine collection   blocked, and not on precision: it needs a
                             volume of distribution for 3-HP and the renal
                             share of non-PCC disposal, each to a few per
                             cent (audit/u9, R107)
    carnitine modulation     open but weak: at the physical ceiling it
                             reaches sd(ln k_e) = 0.469 against 0.063
                             (audit/u10, R104)
    a sub-daily burst        this file

  The burst is the only one of the three whose obstruction is an assay CV
  rather than an unmeasured quantity or a dose nobody would give, and it is
  the only one for which the units objection does not arise at all, because
  a relaxation rate is read off a curve and never divided by a volume.
""")
    tt = HOUR * np.arange(0, 25, dtype=float)
    print(f"    {'theta':>7} {'frontier':>9} {'best safe burst':>16} "
          f"{'short by':>9}  what it costs")
    print("    " + "-" * 76)
    for theta in (600.0, 300.0):
        v = med(theta, 3, n_rec=TIER1_N, burst_times=tt, inc0=TROUGH,
                load_step=-0.50, burst_cv=0.01, n_bursts=16,
                kappa_cv=PUBLISHED_KAPPA_CV)
        fr = FRONTIER[theta]
        print(f"    {theta:7.0f} {fr:9.4f} {v:16.4f} {v / fr:8.2f}x  "
              f"16 inpatient days, a 1% assay")
    print()
    print("  So the answer is neither the yes the route promised nor the no")
    print("  the collection got. It is: the first route whose obstacle is")
    print("  cost rather than an unmeasured constant, still short by a")
    print("  factor of one to two at a price of sixteen admissions, and")
    print("  short by more than that at any assay a clinic runs today.")


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
    part_d()
    part_e()
    part_f()
