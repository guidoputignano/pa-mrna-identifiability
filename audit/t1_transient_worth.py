"""
WHAT A POST-DOSE TRANSIENT IS WORTH, AND AT WHICH CAPACITIES.

Scoping only. Nothing here reaches the paper.

Section 4.1 says that within a few multiples of the marker's time constant
after an infusion the balance is not algebraic and the observation separates
k_e from K_m, and that a transient is therefore what would measure the one
parameter whose prior the design rests on. It does not say what that is worth,
and nobody has asked the question at a disease-range capacity, which is where
Section 6.2 says the estimate has already gone.

This asks it. Three capacities on the paper's own convention A, the recorded
load held at its nominal value while the marker rises:

    theta = 3000   100% of nominal, the illustrative patient
    theta =  300    10% of nominal
    theta =  150     5% of nominal, which is what the disease is

and at each, the recommended monthly record with n post-dose samples added,
log-spaced over the first few hours after one infusion.

TWO CALCULATIONS, AND THEY ANSWER DIFFERENT QUESTIONS.

(1) THE INFORMATION, in the paper's own local model. The transient rows are
    d m(t_j) / d(ln theta, ln K_m, b, ln k_e), obtained by integrating
    eq:local through a dose and differencing each parameter, with the initial
    condition re-solved inside each perturbed integration because a patient
    with a different theta sits at a different trough. This is the OPTIMISTIC
    bound: it assumes the local model describes the transient, which is
    exactly what is in question.

(2) WHETHER THAT TRANSIENT EXISTS AS MODELLED. Every steady-state number in
    the paper rests on the carnitine state being eliminable, which is an
    identity at a fixed point and nothing at all off one. During a transient
    dc/dt is not zero. So the two-state system is integrated against its own
    quasi-steady reduction through the same dose, and the gap between them is
    reported against the assay noise. If the gap is larger than the noise, the
    information in (1) is information about a shape the patient does not have.

WHAT THE TRANSIENT OBSERVATION MODEL ASSUMES BEYOND THE STEADY-STATE ONE.
Remark 2 of Section 2 requires only that the marker track the substrate
MONOTONICALLY: at a fixed point any monotone map carries the same information
about where the fixed point is. A transient requires more. The marker must
track the substrate with the RIGHT TIME CONSTANT, because the quantity being
read off is a rate of approach and not a level. 3-HP is an overflow metabolite
downstream of accumulated propionyl-CoA, with its own pool, its own
distribution volume and its own clearance; nothing published establishes that
its relaxation time equals the substrate's. If it does not, the transient
measures the marker's kinetics and not the enzyme's, and the estimator cannot
tell the two apart from a marker record alone.

Run:  python audit/t1_transient_worth.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import model_A.run as LOC
from common.models import ParamsA, ParamsC, g, alt_flux_C
from design.run import (ASSAY_CV, DOSE_INTERVAL, D_AMP, LAMBDA, LOAD_CV,
                        LOAD_FRAC, SEEDS, load, median_spread)

KM, KE = 100.0, 5.0
LB = load(0.0)
KE_PRIOR_LOG = 0.8
GATE = 0.20
HOUR = 1.0 / 24.0

#: a late dose instant, so the train is at its steady peak and trough rather
#: than still filling; the first dose of the record is not representative
T_DOSE = 50.0 * DOSE_INTERVAL
PEAK = D_AMP / (1.0 - np.exp(-LAMBDA * DOSE_INTERVAL))
TROUGH = PEAK - D_AMP


def op_level(theta, inc=0.0, Lb=LB):
    """The quasi-steady marker level at capacity `theta` and increment `inc`."""
    B = KE * KM + theta + inc - Lb
    return (-B + np.sqrt(B ** 2 + 4.0 * KE * Lb * KM)) / (2.0 * KE)


def sigma_of(theta):
    """The ordinate noise of the steady fit, re-solved at this capacity."""
    m_op = op_level(theta)
    return LOC.sigma_total(LOAD_FRAC, ASSAY_CV,
                           p=ParamsA(theta=theta, Km=KM, ke=KE), m_bar=m_op)


def sigma_marker(theta):
    """
    The assay error in MARKER units.

    `sigma_total` forms its assay term as (ds/dm) * cv * m_bar, so the error it
    is mapping is cv * m_bar on the marker. A transient sample is observed on
    the marker directly, so that is the sd it carries.
    """
    return ASSAY_CV * op_level(theta)


# ------------------------------------------------------------------ steady
def steady_block(theta, times, seed):
    """The paper's own information contribution from a monthly record."""
    t = np.asarray(times, float)
    tau = np.mod(t - T_DOSE, DOSE_INTERVAL)
    inc = PEAK * np.exp(-LAMBDA * tau)
    Ls = LB * (1.0 + LOAD_CV *
               np.random.default_rng(seed).standard_normal(t.size))
    B = KE * KM + theta + inc - Ls
    lv = (-B + np.sqrt(B ** 2 + 4.0 * KE * Ls * KM)) / (2.0 * KE)
    sig = sigma_of(theta)
    F = np.zeros((4, 4))
    for m, d in zip(lv, inc):
        gm = m / (KM + m)
        q = np.array([theta * gm, -(theta + d) * KM * m / (KM + m) ** 2,
                      1.0, KE * m])
        F += np.outer(q, q) / sig ** 2
    return F


# --------------------------------------------------------------- transient
def _traj(theta, Km, ke, b, t_eval):
    """
    m(t) through one dose, from the trough the patient was sitting at.

    The initial condition is the pre-dose quasi-steady level AT THESE
    parameters, so a perturbed patient starts where that patient would be.
    """
    Lb = LB + b
    B = ke * Km + theta + TROUGH - Lb
    m0 = (-B + np.sqrt(B ** 2 + 4.0 * ke * Lb * Km)) / (2.0 * ke)
    inc = lambda t: PEAK * np.exp(-LAMBDA * t)
    rhs = lambda t, y: [Lb - (theta + inc(t)) * (y[0] / (Km + y[0]))
                        - ke * y[0]]
    sol = solve_ivp(rhs, (0.0, float(t_eval[-1])), [m0], t_eval=t_eval,
                    method="LSODA", rtol=1e-10, atol=1e-12)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[0]


def transient_block(theta, t_post, rel=1e-5):
    """
    Information from post-dose marker samples, d m(t)/d(parameters).

    Central differences in (ln theta, ln K_m, b, ln k_e). The load bias b is a
    shift and is differenced absolutely, on the scale of its own prior.
    """
    te = np.asarray(t_post, float)
    base = dict(theta=theta, Km=KM, ke=KE, b=0.0)
    cols = []
    for name, step in (("theta", rel * theta), ("Km", rel * KM),
                       ("b", rel * LOAD_FRAC * LB), ("ke", rel * KE)):
        hi, lo = dict(base), dict(base)
        hi[name] = base[name] + step
        lo[name] = base[name] - step
        d = (_traj(t_eval=te, **hi) - _traj(t_eval=te, **lo)) / (2.0 * step)
        #  log-parameters for theta, K_m and k_e; b stays absolute
        cols.append(d * base[name] if name != "b" else d)
    Q = np.stack(cols, axis=-1)
    sm = sigma_marker(theta)
    return (Q.T @ Q) / sm ** 2


def prior_block():
    return np.diag([0.0, 0.0, 1.0 / (LOAD_FRAC * LB) ** 2,
                    1.0 / KE_PRIOR_LOG ** 2])


def sds(F):
    """(sd(ln theta), sd(ln k_e)) from an information matrix."""
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf, np.inf
    C = np.linalg.inv(F)
    return float(np.sqrt(C[0, 0])), float(np.sqrt(C[3, 3]))


def post_times(n, lo_h=0.25, hi_h=6.0):
    """n draws log-spaced over the first few hours after the infusion."""
    if n == 0:
        return np.zeros(0)
    if n == 1:
        return np.array([np.sqrt(lo_h * hi_h) * HOUR])
    return np.geomspace(lo_h, hi_h, n) * HOUR


# ------------------------------------------------------------------ parent
def parent_gap(theta):
    """
    How far the two-state system departs from its own reduction on a transient.

    Both are the PARENT family: the carnitine flux and its quasi-steady
    elimination. This isolates the carnitine lag, which is the reduction the
    transient breaks, and says nothing about the separate question of whether
    the local model's linear channel describes either.
    """
    p = ParamsC(theta=theta, Km=KM)
    te = np.linspace(0.0, 6.0 * HOUR, 400)
    inc = lambda t: PEAK * np.exp(-LAMBDA * t)

    #  start both at the pre-dose fixed point of the REDUCED system
    #  the reduced parent has TWO roots, an attractor and a separatrix
    #  (Section 2.3), and total clearance is negative at both ends of any
    #  wide bracket, so brentq over [0, 1e6] finds neither. Bracket the
    #  attractor between zero and the clearance maximum.
    f = lambda m: theta * g(m, KM) + alt_flux_C(m, p) - LB
    mm = np.geomspace(1e-4, 1e5, 40000)
    m_peak = float(mm[np.argmax(theta * g(mm, KM) + alt_flux_C(mm, p))])
    if f(m_peak) <= 0:
        return None          # no equilibrium at all: this is a decompensation
    m0 = brentq(f, 1e-9, m_peak, xtol=1e-12)
    u0 = m0 / (p.Kc + m0)
    c0 = p.R0 / (1.0 * (u0 + p.rho0 * (1.0 + m0 / p.Mi)))

    #  rhs_C with the therapeutic increment added to the enzymatic term, which
    #  rhs_C itself does not carry; k = 1, the registered value
    def two_state(t, y):
        m, c = y
        conj = c * m / (p.Kc + m)
        return [LB - (theta + inc(t)) * g(m, KM) - conj,
                p.R0 - conj - p.rho0 * (1.0 + m / p.Mi) * c]

    two = solve_ivp(two_state, (0.0, te[-1]), [m0, c0], t_eval=te,
                    method="LSODA", rtol=1e-10, atol=1e-12)
    red = solve_ivp(lambda t, y: [LB - (theta + inc(t)) * g(y[0], KM)
                                  - alt_flux_C(y[0], p)],
                    (0.0, te[-1]), [m0], t_eval=te, method="LSODA",
                    rtol=1e-10, atol=1e-12)
    if not (two.success and red.success):
        return None
    dm = np.abs(two.y[0] - red.y[0])
    return dict(m0=m0, peak_abs=float(dm.max()),
                peak_rel=float((dm / red.y[0]).max()),
                sigma_marker=sigma_marker(theta),
                in_sigmas=float(dm.max() / sigma_marker(theta)))


def transient_numbers():
    """
    The four figures Section 4.1 quotes, generated rather than typed.

    Only what the paper carries: what eight post-dose samples move at each
    capacity, what they move for k_e, the angle between the transient and the
    degeneracy, and the excursion in assay sigmas.
    """
    from common.numbers import write_numbers
    monthly = np.arange(0.0, 730.0, 30.4)
    N = {}
    for th, tag in ((3000.0, "nominal"), (150.0, "disease")):
        S = [steady_block(th, monthly, s) for s in SEEDS]
        tb = transient_block(th, post_times(8))
        base = median_spread([sds(prior_block() + F)[0] for F in S])[0]
        with8 = median_spread([sds(prior_block() + F + tb)[0] for F in S])[0]
        ke0 = median_spread([sds(prior_block() + F)[1] for F in S])[0]
        ke8 = median_spread([sds(prior_block() + F + tb)[1] for F in S])[0]
        N["trans_gain_" + tag] = 100.0 * (1.0 - with8 / base)
        N["trans_ke_gain_" + tag] = 100.0 * (1.0 - ke8 / ke0)
        #  the angle between the transient's best direction and the record's worst
        worst = np.linalg.eigh(prior_block() + S[0])[1][:, 0]
        worst = worst / np.linalg.norm(worst)
        te = post_times(8)
        base_p = dict(theta=th, Km=KM, ke=KE, b=0.0)
        cols = []
        for name, step in (("theta", 1e-5 * th), ("Km", 1e-5 * KM),
                           ("b", 1e-5 * LOAD_FRAC * LB), ("ke", 1e-5 * KE)):
            hi, lo = dict(base_p), dict(base_p)
            hi[name] = base_p[name] + step
            lo[name] = base_p[name] - step
            d = (_traj(t_eval=te, **hi) - _traj(t_eval=te, **lo)) / (2.0 * step)
            cols.append(d * base_p[name] if name != "b" else d)
        vt = np.linalg.svd(np.stack(cols, axis=-1), full_matrices=False)[2]
        best = vt[0] / np.linalg.norm(vt[0])
        N["trans_angle_" + tag] = float(
            np.degrees(np.arccos(min(abs(float(best @ worst)), 1.0))))
        m = _traj(th, KM, KE, 0.0, np.r_[1e-9, te])
        N["trans_sigmas_" + tag] = abs(m[0] - m[-1]) / sigma_marker(th)
    fmt = {k: ("{:.1f}" if "angle" in k or "sigmas" in k else "{:.1f}")
           for k in N}
    write_numbers("transient", N, fmt)
    return N


# -------------------------------------------------------------------- main
def main():
    print(__doc__)
    bar = "=" * 78
    monthly = np.arange(0.0, 730.0, 30.4)
    caps = [(3000.0, "100%"), (300.0, "10%"), (150.0, "5%")]

    print(bar)
    print("1. WHAT THE RECORD CARRIES, WITH AND WITHOUT A TRANSIENT")
    print(bar)
    print(f"  monthly record: {monthly.size} draws over two years, "
          f"medians over {len(SEEDS)} intake records")
    print(f"  transient: n draws log-spaced over 0.25 to 6 h after one "
          f"infusion, assay {100*ASSAY_CV:.0f}% of the operating level\n")
    hdr = (f"  {'theta':>8} {'m_op':>7} {'sigma':>7} {'sig_m':>7} "
           f"{'n':>4} {'sd(ln th)':>22} {'sd(ln ke)':>22}")
    for th, lab in caps:
        print(f"\n  --- theta = {th:.0f} ({lab} of nominal), "
              f"operating marker {op_level(th):.2f} "
              f"({op_level(th)/op_level(3000.0):.2f}x managed) ---")
        print(hdr)
        for n in (0, 3, 4, 5, 6, 8):
            tb = transient_block(th, post_times(n)) if n else np.zeros((4, 4))
            th_s, ke_s = [], []
            for s in SEEDS:
                a, b = sds(prior_block() + steady_block(th, monthly, s) + tb)
                th_s.append(a); ke_s.append(b)
            mt, lt, ht = median_spread(th_s)
            mk, lk, hk = median_spread(ke_s)
            flag = " *" if mt < GATE else "  "
            print(f"  {th:8.0f} {op_level(th):7.2f} {sigma_of(th):7.2f} "
                  f"{sigma_marker(th):7.3f} {n:4d} "
                  f"{mt:9.4f} [{lt:7.4f},{ht:7.4f}]{flag}"
                  f"{mk:9.4f} [{lk:7.4f},{hk:7.4f}]")

    print("\n" + bar)
    print("2. HOW MANY TRANSIENT DRAWS REACH THE GATE, IF ANY NUMBER DOES")
    print(bar)
    print("  n beyond about a dozen means the window repeated across dose")
    print("  cycles: the information scales the same way either.\n")
    for th, lab in caps:
        tb8 = transient_block(th, post_times(8))
        S = [steady_block(th, monthly, s) for s in SEEDS]
        def sd_at(lam):
            return median_spread([sds(prior_block() + F + lam * tb8)[0]
                                  for F in S])[0]
        if sd_at(1.0) < GATE:
            print(f"  theta = {th:6.0f} ({lab:>4}): already inside the gate "
                  f"on the steady record alone ({sd_at(1.0):.4f})")
            continue
        #  the block for n = 8 scaled by lam is the block for n = 8 lam draws
        #  at the same phases, which is the window repeated across cycles
        hi = 1e14
        if sd_at(hi) > GATE:
            print(f"  theta = {th:6.0f} ({lab:>4}): NOT reached at any n")
            continue
        lam = np.exp(brentq(lambda L: sd_at(np.exp(L)) - GATE,
                            0.0, np.log(hi), xtol=1e-4))
        ladder = "  ".join(f"n={n}: {sd_at(n / 8.0):.4f}"
                           for n in (8, 100, 1000, 10000))
        print(f"  theta = {th:6.0f} ({lab:>4}): n = {8 * lam:,.0f} draws")
        print(f"                      {ladder}")

    print("\n" + bar)
    print("3. WHAT THE TRANSIENT IS WORTH FOR k_e ALONE")
    print(bar)
    print("  k_e is the one load-bearing prior in the paper. Without the")
    print("  transient the record barely moves it off the prior's own 0.8.\n")
    for th, lab in caps:
        row = [f"  theta = {th:6.0f} ({lab:>4}):"]
        for n in (0, 3, 8, 32):
            tb = transient_block(th, post_times(n)) if n else np.zeros((4, 4))
            med = median_spread([sds(prior_block() + steady_block(th, monthly, s)
                                     + tb)[1] for s in SEEDS])[0]
            row.append(f"n={n:<4d} {med:7.4f}")
        print("  ".join(row))

    print("\n" + bar)
    print("4. DOES THAT TRANSIENT EXIST? THE CARNITINE STATE DOES NOT KEEP UP")
    print(bar)
    print("  Two-state parent against its own quasi-steady reduction, same")
    print("  constants, through the same dose. The reduction is an identity")
    print("  at a fixed point and carries no guarantee off one.\n")
    for th, lab in caps:
        r = parent_gap(th)
        if r is None:
            print(f"  theta = {th:6.0f} ({lab:>4}): no fixed point to start from")
            continue
        print(f"  theta = {th:6.0f} ({lab:>4}): m0 = {r['m0']:8.3f}, "
              f"worst gap {r['peak_abs']:8.4f} uM = {100*r['peak_rel']:6.2f}% "
              f"of the level = {r['in_sigmas']:6.2f} assay sigmas")

    print("\n" + bar)
    print("5. WHY IT CARRIES SO LITTLE: THE TRANSIENT IS ORTHOGONAL TO THE")
    print("   DIRECTION THE RECORD CANNOT SEE")
    print(bar)
    print("  The steady record's worst-determined eigendirection against the")
    print("  transient's best-determined one, in (ln theta, ln K_m, b, ln k_e).")
    print("  A transient that helped would lie ALONG the first, not across it.\n")
    for th, lab in caps:
        F = prior_block() + steady_block(th, monthly, 0)
        worst = np.linalg.eigh(F)[1][:, 0]
        worst = worst / np.linalg.norm(worst)
        te = post_times(8)
        base = dict(theta=th, Km=KM, ke=KE, b=0.0)
        cols = []
        for name, step in (("theta", 1e-5 * th), ("Km", 1e-5 * KM),
                           ("b", 1e-5 * LOAD_FRAC * LB), ("ke", 1e-5 * KE)):
            hi, lo = dict(base), dict(base)
            hi[name] = base[name] + step
            lo[name] = base[name] - step
            d = (_traj(t_eval=te, **hi) - _traj(t_eval=te, **lo)) / (2.0 * step)
            cols.append(d * base[name] if name != "b" else d)
        sv, vt = np.linalg.svd(np.stack(cols, axis=-1), full_matrices=False)[1:]
        best = vt[0] / np.linalg.norm(vt[0])
        cos = abs(float(best @ worst))
        print(f"  theta = {th:6.0f} ({lab:>4}): "
              f"|cos| = {cos:.4f}, {np.degrees(np.arccos(min(cos, 1.0))):5.2f} "
              f"degrees apart; the transient's own condition number is "
              f"{sv[0] / sv[-1]:.2g}")

    print("\n" + bar)
    print("6. THE MARKER EXCURSION THE TRANSIENT CARRIES, FOR SCALE")
    print(bar)
    print("  The transient is not weak. What is weak is its component along")
    print("  the direction the steady record is missing.\n")
    for th, lab in caps:
        te = post_times(8)
        m = _traj(th, KM, KE, 0.0, np.r_[1e-9, te])
        print(f"  theta = {th:6.0f} ({lab:>4}): the marker falls "
              f"{m[0]:7.3f} to {m[-1]:6.3f} uM over six hours, which is "
              f"{abs(m[0] - m[-1]) / sigma_marker(th):5.1f} assay sigmas")

    print("\n" + bar)
    print("R1. RE-CHECK: the transient block must vanish with no samples")
    print(bar)
    z = np.zeros((4, 4))
    a0 = sds(prior_block() + steady_block(3000.0, monthly, 0) + z)[0]
    a1 = sds(prior_block() + steady_block(3000.0, monthly, 0))[0]
    print(f"  n=0 against no block at all: {a0:.6f} vs {a1:.6f}, "
          f"identical: {abs(a0 - a1) < 1e-12}")

    print("\nR2. RE-CHECK: the steady block reproduces the paper at nominal")
    med, lo, hi = median_spread([sds(prior_block()
                                     + steady_block(3000.0, monthly, s))[0]
                                 for s in SEEDS])
    print(f"  sd(ln theta) monthly at theta = 3000: {med:.4f} "
          f"[{lo:.4f}, {hi:.4f}]  (the paper's Table 3 row: 0.1699 "
          f"[0.1684, 0.1709])")

    print("\nR3. RE-CHECK: finite-difference step independence")
    for rel in (1e-4, 1e-5, 1e-6):
        tb = transient_block(150.0, post_times(8), rel=rel)
        print(f"  rel = {rel:.0e}: sd(ln theta) = "
              f"{sds(prior_block() + steady_block(150.0, monthly, 0) + tb)[0]:.6f}")


if __name__ == "__main__":
    main()
    #  main() only PRINTS. transient_numbers() is what the paper reads, and
    #  nothing called it, so `make figures` ran this file and left
    #  results/numbers/transient.tex untouched -- the same gap the design
    #  module had with quasi_steady_check(). A generator the build does not
    #  reach is a generator that cannot catch drift.
    transient_numbers()
