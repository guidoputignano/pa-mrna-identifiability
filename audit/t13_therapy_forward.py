"""
Item 22 -- how the therapy behaves as a function of residual capacity.

FORWARD SIMULATION, NOT INFERENCE. Every other question in this repository is
"what can a record tell us about theta". This one is "given theta, what does
the marker do", which needs no identifiability at all and therefore sidesteps
the six sections that close it. Nothing below estimates anything.

THREE CONSTRAINTS ON THE SETUP, WHICH MATTER MORE THAN THE RUNS.

1. THETA IS SWEPT, NOT ESTIMATED AND NOT CITED. The literature gives a
   fraction -- residual PCC activity under 5% of control -- and never an
   absolute V_max, so it says where to sweep and not what the value is. The
   sweep is a percentage of the model's normal: 1, 2, 5, 10, 20%, with 100%
   as reference. EVERY NUMBER BELOW IS CONDITIONAL ON THE SWEPT VALUE and is
   printed with that percentage attached. Nothing here is a measurement of a
   patient's capacity and nothing should be quoted as one.

2. THE OBJECTIVE IS MARKER EXPOSURE, NOT THE SIGN OF theta - theta_crit.
   Section 4.6 shows that threshold moves by a factor of 24 across its own
   constants and flips sign for the same patient with perfect data, so any
   objective built on it inherits that. Exposure is the only quantity here
   with a published endpoint association:

       mean      time-averaged marker over a cycle, the exposure integral
       peak      the largest marker the cycle reaches
       above     the fraction of the cycle spent above the managed level

3. THE MODEL STOPS HOLDING BEFORE THE INTERESTING REGIME DOES. The parent
   \\eqref{eq:parent} has a clearance maximum, and above it there is no steady
   state to sit at; higher still it loses its stable equilibrium altogether.
   Both bounds are computed here from model_C rather than typed, and every
   result that crosses the first is MARKED rather than reported beside the
   rest. The missed-dose runs go there, and that is itself the finding.

WHAT THIS CANNOT BE READ AS. Section 5.6 already says a protocol written from
this work is not a rule for escalating an individual patient. A forward
simulation at an ASSUMED capacity is weaker evidence than that, not stronger.
What it can legitimately inform is a protocol's missed-dose contingency,
which is a document-level decision and not a patient-level one.

Nothing here records an expected result in advance.

Run:  python audit/t13_therapy_forward.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from common.models import g, ParamsC, clearance_C_quasistatic
from design.run import P, M_BAR, D_AMP, LAMBDA, DOSE_INTERVAL, load
import model_C.run as C

import importlib.util as _ilu
_spec = _ilu.spec_from_file_location(
    "n4", Path(__file__).resolve().parent / "n4_theta_scope.py")
N4 = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(N4)

L = float(load(0))
T = DOSE_INTERVAL
KM, KE = P.Km, P.ke
THETA_NORMAL = P.theta

#: the steady-state peak of the booked dose train
PEAK = D_AMP / (1.0 - np.exp(-LAMBDA * T))

#: the managed level, which is this repository's reference throughout:
#: `load()` is calibrated so the patient sits at M_BAR with the identification
#: increment off, and Section 7.1's ceiling is quoted as a multiple of it.
MANAGED = float(M_BAR)

#: the marker's limit as theta -> 0: the alternative channel does not
#: saturate, so the marker approaches L/k_e rather than rising without bound.
#: This is \eqref{eq:ceilinglimit}, and the paper quotes it as 6.45 x managed.
CEILING = L / KE

#: the capacity sweep, as a percentage of the model's normal
PCTS = (100, 20, 10, 5, 2, 1)


def validity_bounds():
    """
    Where the parent stops supporting a trajectory, computed not typed.

    Two multiples of the managed level, both properties of the parent's
    clearance curve rather than of the local model: the marker at which the
    clearance curve peaks, above which there is no steady state to sit at,
    and the upper equilibrium at baseline load, above which the stable one is
    gone. Section 8 quotes them as 3.4 and 12.6 times managed.
    """
    PC = ParamsC()
    _, m_peak = C.clearance_peak(PC)
    grid = np.geomspace(1e-3, 1e5, 400000)
    f = C.S_BASE - clearance_C_quasistatic(grid, PC)
    sign = np.sign(f)
    roots = [brentq(lambda m: C.S_BASE - clearance_C_quasistatic(m, PC),
                    grid[i], grid[i + 1])
             for i in np.where(np.diff(sign) != 0)[0]]
    attractor, upper = min(roots), max(roots)
    return m_peak / attractor, upper / attractor, attractor


def qs(theta, d):
    """The quasi-steady marker at a frozen increment `d`."""
    return brentq(lambda m: (theta + d) * g(m, KM) + KE * m - L, 1e-12, 1e9)


def _rhs(theta):
    return lambda t, u: [L - (theta + u[1]) * g(u[0], KM) - KE * u[0],
                         -LAMBDA * u[1]]


def steady_cycle(theta, amp=D_AMP, period=T, n=4001, tol=1e-11, itmax=300):
    """
    One steady-state cycle, as the fixed point of the cycle map.

    Delta is periodic and known in closed form on a steady train, so only the
    marker needs the fixed point. Iterating the cycle map is used rather than
    integrating forty cycles because it converges in two or three passes and
    because "forty cycles" is a guess about the slowest capacity in the sweep.
    """
    peak = amp / (1.0 - np.exp(-LAMBDA * period))
    m = qs(theta, peak)
    rhs = lambda t, u: [L - (theta + peak * np.exp(-LAMBDA * t))
                        * g(u[0], KM) - KE * u[0]]
    for it in range(itmax):
        s = solve_ivp(rhs, (0, period), [m], rtol=1e-11, atol=1e-13)
        nxt = float(s.y[0, -1])
        if abs(nxt - m) < tol * max(m, 1.0):
            break
        m = nxt
    s = solve_ivp(rhs, (0, period), [m], rtol=1e-11, atol=1e-13,
                  dense_output=True)
    tt = np.linspace(0.0, period, n)
    return tt, s.sol(tt)[0], peak * np.exp(-LAMBDA * tt)


def forward(theta, doses, amp=D_AMP, t_end=None, m0=None, d0=0.0, per_day=40,
            amps=None):
    """
    Integrate (marker, increment) with the increment JUMPING at each dose.

    The increment obeys its own equation, d(Delta)/dt = -lambda Delta, with a
    reset of +amp at each dose time, so a dose is an exact discontinuity and
    not a stiff forcing term the solver has to resolve. Sampling is a fixed
    number of points per DAY rather than per run, because the missed-dose
    scenarios are read off their peaks and a fixed total would resolve a long
    run more coarsely than a short one.
    """
    doses = np.asarray(doses, float)
    #  per-dose amplitude, so that TWO doses at one time are two doses. The
    #  first version keyed the jump off `np.isclose(doses, a)` after a
    #  np.unique() on the edges, which silently merged a doubled dose back
    #  into a single one: the "skipped, double at next" row came out
    #  bit-identical to "skipped" and looked like a real null result.
    amps = (np.full(doses.shape, float(amp)) if amps is None
            else np.asarray(amps, float))
    order = np.argsort(doses)
    doses, amps = doses[order], amps[order]
    t_end = float(doses[-1] + T if t_end is None else t_end)
    y = np.array([qs(theta, d0) if m0 is None else float(m0), float(d0)])
    rhs = _rhs(theta)
    edges = np.unique(np.concatenate([[0.0], doses[doses <= t_end], [t_end]]))
    ts, ms, ds = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        hit = np.isclose(doses, a)
        if hit.any():
            y = y + np.array([0.0, float(amps[hit].sum())])
        if b <= a:
            continue
        s = solve_ivp(rhs, (a, b), y, rtol=1e-10, atol=1e-12,
                      dense_output=True)
        tt = np.linspace(a, b, max(8, int(round(per_day * (b - a)))))
        uu = s.sol(tt)
        ts.append(tt), ms.append(uu[0]), ds.append(uu[1])
        y = s.y[:, -1]
    return np.concatenate(ts), np.concatenate(ms), np.concatenate(ds)


def exposure(t, m, ref=None):
    """mean, peak and the fraction of the window spent above `ref`."""
    ref = MANAGED if ref is None else ref
    span = t[-1] - t[0]
    return (float(np.trapezoid(m, t) / span), float(m.max()),
            float(np.trapezoid((m > ref).astype(float), t) / span))


def admissible_D(theta):
    """
    The largest booked dose with theta + Delta_peak <= normal.

    IMPORTED rather than rewritten, and the first draft here got it wrong.
    The check infers normal from the capacity as `theta / frac` with frac the
    reference residual activity of 5%, so the bound is Delta_peak <= 19 theta
    and it TIGHTENS as capacity falls. Rewriting it as a fixed
    normal = 3000 happens to agree at 5% of nominal, where theta/0.05 is
    exactly 3000, and disagrees everywhere else: it reported the booked dose
    inadmissible at NORMAL capacity, which is not what Section 7.1 says and
    is not what the check does. Using audit/n4_theta_scope.py's own function
    removes the second spelling along with the error.
    """
    return N4.d_ceiling(float(theta))


def admissibility(theta):
    """(theta + Delta_peak) / normal, the ratio the check reports."""
    return (theta + PEAK) / (theta / N4.REFERENCE_FRAC)


def bars(title):
    print("\n" + "=" * 84)
    print(title)
    print("=" * 84)


def flag(m_peak, r_peak, r_upper):
    """Mark a result whose trajectory leaves what the parent supports."""
    if m_peak > r_upper * MANAGED:
        return "  <<< PAST THE STABLE EQUILIBRIUM"
    if m_peak > r_peak * MANAGED:
        return "  <<< past the clearance maximum"
    return ""


def s0_setup(r_peak, r_upper, attractor):
    bars("0. THE SETUP, AND THE THREE CONSTRAINTS IT IS BUILT UNDER")
    print(f"  managed level (the reference throughout)   {MANAGED:8.3f}")
    print(f"  marker ceiling as theta -> 0, L/k_e        {CEILING:8.3f}"
          f"   = {CEILING/MANAGED:.2f} x managed")
    print(f"  booked dose D                              {D_AMP:8.1f}"
          f"   every {T:.0f} d, half-life {np.log(2)/LAMBDA:.2f} d")
    print(f"  its steady-state peak increment            {PEAK:8.1f}"
          f"   = {PEAK/THETA_NORMAL:.3f} x normal capacity")
    print()
    print("  WHERE THE MODEL STOPS HOLDING, from model_C rather than typed:")
    print(f"    parent attractor at baseline load        {attractor:8.3f}")
    print(f"    clearance maximum                        {r_peak:8.3f}"
          f" x managed  ->  marker {r_peak*MANAGED:.2f}")
    print(f"    upper equilibrium lost                   {r_upper:8.3f}"
          f" x managed  ->  marker {r_upper*MANAGED:.2f}")
    print("    Above the first there is no steady state to sit at. Every")
    print("    result crossing it is marked below rather than reported beside")
    print("    the rest.")


def s1_booked(r_peak, r_upper):
    bars("1. RUN (a). WHAT THE BOOKED SCHEDULE DOES, BY CAPACITY")
    print("  The same booked dose at every capacity, because a dose in mg/kg")
    print("  is fixed in absolute terms while theta is not. Delta_peak/theta")
    print("  is therefore the whole story, and it runs over a factor of 200")
    print("  across this sweep: the SAME INTERVENTION at the two ends is not")
    print("  the same intervention at all.\n")
    print(f"  {'theta':>6} {'theta':>8} {'D_pk/th':>8} {'untreated':>10} "
          f"{'mean':>8} {'peak':>8} {'trough':>8} {'fold':>6} {'above':>7} "
          f"{'(th+D)/norm':>12}")
    print(f"  {'%':>6} {'':>8} {'':>8} {'marker':>10} {'':>8} {'':>8} "
          f"{'':>8} {'':>6} {'managed':>7} {'':>12}")
    print("  " + "-" * 96)
    out = {}
    for pct in PCTS:
        th = THETA_NORMAL * pct / 100.0
        tt, m, d = steady_cycle(th)
        mean, pk, above = exposure(tt, m)
        adm = admissibility(th)
        out[pct] = (mean, pk, float(m.min()), above, adm)
        print(f"  {pct:6d} {th:8.1f} {PEAK/th:8.1f} {qs(th,0):10.3f} "
              f"{mean:8.3f} {pk:8.3f} {m.min():8.3f} {pk/m.min():6.2f} "
              f"{100*above:6.1f}% {adm:12.2f}"
              + flag(pk, r_peak, r_upper))
    print()
    print("  THE SWING GROWS AS CAPACITY FALLS, which is the opposite of what")
    print("  a fixed dose is usually assumed to do: at 100% the marker moves")
    print(f"  by a factor of {out[100][1]/out[100][2]:.2f} over a cycle and at 1% by "
          f"{out[1][1]/out[1][2]:.2f}, because the")
    print("  increment is a larger and larger multiple of what is left.")
    print()
    print("  ADMISSIBILITY, WHICH IS SECTION 7.1's CHECK READ AS A TRAJECTORY.")
    print("  (theta + Delta_peak)/normal above 1 means the booked dose drives")
    print("  total capacity past the enzyme's own normal V_max, which is not a")
    print("  therapy at that capacity but a different enzyme:")
    for pct in PCTS:
        th = THETA_NORMAL * pct / 100.0
        adm = admissibility(th)
        if adm > 1.0:
            print(f"    {pct:3d}% of nominal: {adm:.2f} x normal, largest "
                  f"admissible D = {admissible_D(th):.0f} against the booked "
                  f"{D_AMP:.0f}")
    return out


def s2_interval_dose(r_peak, r_upper):
    bars("2. RUN (b). INTERVAL AND DOSE AS A THERAPEUTIC TRADE")
    print("  Section 4.2 has the increment's trough against the interval and")
    print("  reads it as what IDENTIFICATION costs. The same sweep read as")
    print("  exposure asks a different question, and the two do not have to")
    print("  agree. Delivered rate is held fixed, D proportional to h, so a")
    print("  longer interval is a larger dose less often and not less drug.\n")
    HS = (7.0, 10.0, 14.0, 17.0, 21.0, 24.0, 28.0)
    best = {}
    for pct in PCTS:
        th = THETA_NORMAL * pct / 100.0
        print(f"  theta = {pct}% of normal ({th:.0f}):")
        print(f"    {'h  d':>6} {'D':>8} {'D_pk/th':>8} {'mean':>8} "
              f"{'peak':>8} {'trough':>8} {'D_trough':>9} {'above':>7}")
        rows = []
        for h in HS:
            amp = D_AMP * h / T
            tt, m, d = steady_cycle(th, amp=amp, period=h)
            mean, pk, above = exposure(tt, m)
            rows.append((h, amp, mean, pk, float(m.min()), float(d.min()),
                         above))
            print(f"    {h:6.0f} {amp:8.0f} {amp/(1-np.exp(-LAMBDA*h))/th:8.1f} "
                  f"{mean:8.3f} {pk:8.3f} {m.min():8.3f} {d.min():9.1f} "
                  f"{100*above:6.1f}%" + flag(pk, r_peak, r_upper))
        b = min(rows, key=lambda r: r[2])
        best[pct] = b[0]
        edge = " -- AT THE EDGE OF THE SWEPT RANGE, so not an optimum" \
            if b[0] in (HS[0], HS[-1]) else ""
        print(f"    least exposure at h = {b[0]:.0f} d "
              f"(mean {b[2]:.3f}){edge}\n")
    print(f"  DOES THE BEST INTERVAL DEPEND ON THETA? "
          f"{'yes' if len(set(best.values())) > 1 else 'NO'}: "
          + ", ".join(f"{p}% -> {h:.0f} d" for p, h in best.items()))
    print()
    print("  BUT THE HONEST READING IS THAT THERE IS NO INTERIOR OPTIMUM AT")
    print("  ALL, AND SAYING '7 DAYS' WOULD BE REPORTING THE EDGE OF A SWEEP")
    print("  AS A RESULT. Exposure falls monotonically as the interval")
    print("  shortens, towards the limit of continuous delivery at the same")
    print("  rate, where the increment settles at R/lambda:")
    R = D_AMP / T
    d_inf = R / LAMBDA
    print(f"    delivered rate R = {R:.1f}/d, so Delta -> R/lambda = "
          f"{d_inf:.1f} and the marker -> its value there\n")
    print(f"    {'theta':>6} {'h=28':>8} {'h=14':>8} {'h=7':>8} {'h=3.5':>8} "
          f"{'h=1':>8} {'continuous':>11}")
    for pct in PCTS:
        th = THETA_NORMAL * pct / 100.0
        cells = []
        for h in (28.0, 14.0, 7.0, 3.5, 1.0):
            tt, m, _ = steady_cycle(th, amp=D_AMP * h / T, period=h)
            cells.append(exposure(tt, m)[0])
        print(f"    {pct:5d}% " + " ".join(f"{c:8.3f}" for c in cells)
              + f" {qs(th, d_inf):11.3f}")
    print()
    print("  So within this model 'more often' is monotonically better on")
    print("  exposure and the interval question has no pharmacological")
    print("  answer, only a practical one. What DOES move with theta is the")
    print("  peak, and that is what the admissibility ceiling binds on, so")
    print("  the trade a protocol actually faces is peak against burden and")
    print("  not exposure against interval.")
    print()
    print("  And the dose at a fixed 14-day interval, which is the other axis:")
    print(f"    {'theta':>6} {'D':>8} {'mean':>8} {'peak':>8} {'trough':>8} "
          f"{'above':>7} {'admissible':>11}")
    for pct in (100, 5):
        th = THETA_NORMAL * pct / 100.0
        for amp in (1000.0, 2000.0, 2672.0, 4000.0, 6000.0, 12000.0):
            tt, m, d = steady_cycle(th, amp=amp)
            mean, pk, above = exposure(tt, m)
            ok = "yes" if amp <= admissible_D(th) else "NO"
            print(f"    {pct:5d}% {amp:8.0f} {mean:8.3f} {pk:8.3f} "
                  f"{m.min():8.3f} {100*above:6.1f}% {ok:>11}"
                  + flag(pk, r_peak, r_upper))
    return best


#: how long after the disturbance to watch, and how much train to carry
#: past it so the window is never reading the end of the integration
WINDOW = 10 * T
N_DOSES = 14


def steady_start(theta, amp=D_AMP, period=T):
    """
    The state that makes the train periodic from t = 0, dose included.

    Starting a missed-dose run from the untreated level and integrating
    thirty doses of burn-in costs four times the runtime and leaves a
    question about whether thirty was enough. Starting ON the cycle removes
    both: `steady_cycle` already solves for the periodic marker, and the
    increment just after a dose is the closed-form peak, so seeding with
    d0 = peak - amp makes the t = 0 dose land exactly on it.
    """
    tt, m, d = steady_cycle(theta, amp=amp, period=period)
    return float(m[0]), float(d[0] - amp)


def _on_grid(run, t1, n=8001):
    """A run resampled onto a common grid, so two runs can be subtracted."""
    ts, ms, ds = run
    grid = np.linspace(0.0, t1, n)
    return grid, np.interp(grid, ts, ms), np.interp(grid, ts, ds)


def s3_missed(r_peak, r_upper):
    bars("3. RUN (c). THE MISSED DOSE")
    print("  The most concrete of the three, and the one most likely to be")
    print("  asked of a protocol. A dose is disturbed deep inside the steady")
    print("  train and the trajectory is read against the undisturbed one.")
    print(f"  At a half-life of {np.log(2)/LAMBDA:.2f} d, a {2*T:.0f}-day gap decays the")
    print(f"  increment by a factor of {np.exp(LAMBDA*2*T):.0f}, so a skipped dose is")
    print("  most of an untreated cycle rather than a slightly smaller one.\n")
    base = list(np.arange(N_DOSES) * T)
    t_p = 0.0                     # the run starts ON the cycle, so dose 0 is it
    t1 = WINDOW + 2 * T
    rows = {}
    for pct in PCTS:
        th = THETA_NORMAL * pct / 100.0
        m0, d0 = steady_start(th)
        grid, m_u, _ = _on_grid(
            forward(th, base, t_end=t1, m0=m0, d0=d0, per_day=60), t1)
        w = grid <= WINDOW
        base_mean = float(np.trapezoid(m_u[grid <= T], grid[grid <= T]) / T)
        print(f"  theta = {pct}% of normal ({th:.0f}), undisturbed cycle mean "
              f"{base_mean:.3f}:")
        print(f"    {'scenario':>26} {'peak':>8} {'xmanaged':>9} "
              f"{'excess exposure':>16} {'rejoin d':>9}")
        rest = [d for d in base if d != t_p]
        for lab, doses, amps in (
                ("dose 1 d late", rest + [t_p + 1.0], None),
                ("dose 3 d late", rest + [t_p + 3.0], None),
                ("dose 7 d late", rest + [t_p + 7.0], None),
                ("skipped to the next dose", rest, None),
                ("skipped, double at next", rest,
                 [2.0 * D_AMP if d == t_p + T else D_AMP for d in rest])):
            run = forward(th, doses, t_end=t1, m0=m0, d0=d0, per_day=60,
                          amps=amps)
            _, m_p, d_p = _on_grid(run, t1)
            pk = float(m_p[w].max())
            #  excess exposure, as a multiple of one undisturbed cycle's worth
            excess = float(np.trapezoid(m_p[w] - m_u[w], grid[w])
                           / (base_mean * T))
            #  rejoin: the last time the two are more than 1% apart
            gap = np.abs(m_p - m_u) / np.maximum(m_u, 1e-12)
            after = np.where(w & (gap > 0.01))[0]
            rejoin = (grid[after[-1]] - t_p) if after.size else 0.0
            note = flag(pk, r_peak, r_upper)
            if amps is not None:
                d_pk = float(d_p.max())
                note += (f"   [double peaks at theta+Delta = "
                         f"{(th + d_pk)/(th/N4.REFERENCE_FRAC):.2f} x normal"
                         f"{', INADMISSIBLE' if th + d_pk > th/N4.REFERENCE_FRAC else ''}]")
            print(f"    {lab:>26} {pk:8.3f} {pk/MANAGED:9.2f} "
                  f"{excess:15.3f}x {rejoin:9.1f}" + note)
            rows[(pct, lab)] = (pk, excess, rejoin)
        print()
    print("  THE UNTREATED CEILING IS THE SCALE TO READ THESE AGAINST: "
          f"{CEILING:.2f},")
    print(f"  which is {CEILING/MANAGED:.2f} x managed and sits ABOVE the clearance")
    print(f"  maximum at {r_peak:.2f} x managed. So a long enough gap takes the")
    print("  marker into the regime where the parent says there is no steady")
    print("  state at all, and the honest reading of these rows is that the")
    print("  model can say a gap is bad and cannot say how bad.")
    return rows


def rechecks(r_peak, r_upper):
    """
    The adversarial pass, run BEFORE any of this is reported.

    The standing policy in REVIEW_LOG.md: a first report of a new measurement
    is provisional until this has run. It has already earned its place twice
    in this file. The admissibility check was first written with a fixed
    normal = 3000 and reported the booked dose inadmissible at NORMAL
    capacity, which is not what Section 7.1 says; and the doubled dose was
    silently merged into a single one by a np.unique() on the dose times, so
    "skipped, double at next" came out bit-identical to "skipped" and read as
    a real null result. Neither would have survived being asked to fail.
    """
    ok = []
    bars("R. THE ADVERSARIAL PASS")

    print("  R1. steady_cycle against a long forward run of the same train,")
    print("      which shares no code with the cycle-map fixed point.")
    for pct in (100, 5):
        th = THETA_NORMAL * pct / 100.0
        tt, m, _ = steady_cycle(th)
        ts, ms, _ = forward(th, np.arange(60) * T, per_day=200)
        late = ts >= 58 * T
        d_hi = abs(ms[late].max() - m.max()) / m.max()
        d_lo = abs(ms[late].min() - m.min()) / m.min()
        ok.append(max(d_hi, d_lo) < 2e-3)
        print(f"      theta {pct:3d}%: peak {d_hi:.2e}, trough {d_lo:.2e} "
              f"apart  {'PASS' if ok[-1] else 'FAIL'}")

    print("\n  R2. the constants, against what the paper already commits to.")
    for lab, got, want, tol in (
            ("clearance max / managed", r_peak, 3.4, 0.05),
            ("upper equilibrium / managed", r_upper, 12.6, 0.05),
            ("ceiling / managed", CEILING / MANAGED, 6.45, 0.01),
            ("largest admissible D at 5%",
             admissible_D(0.05 * THETA_NORMAL), 2672.0, 1.0),
            ("(th+D_pk)/normal at 5%",
             admissibility(0.05 * THETA_NORMAL), 2.18, 0.01),
            ("28-day decay of the increment",
             float(np.exp(LAMBDA * 2 * T)), 256.0, 1.0)):
        ok.append(abs(got - want) <= tol)
        print(f"      {lab:>30}: {got:12.4f} against {want:9.2f}  "
              f"{'PASS' if ok[-1] else 'FAIL'}")

    print("\n  R3. the increment channel against its closed form. The solver")
    print("      integrates dDelta/dt = -lambda Delta with resets; the train")
    print("      is also a geometric sum. They must agree off the dose")
    print("      times, where the output carries both sides of the jump.")
    ts, _, ds = forward(THETA_NORMAL, np.arange(12) * T, per_day=200)
    closed = np.zeros_like(ts)
    for k in range(12):
        on = ts >= k * T
        closed[on] += D_AMP * np.exp(-LAMBDA * (ts[on] - k * T))
    #  A dose time appears TWICE in the output, once ending the segment
    #  before the jump and once beginning the segment after it, which is the
    #  honest way to carry a discontinuity. The closed form has only the
    #  post-jump value there, so the comparison is made OFF the dose times.
    #  The first version of this check did not, and failed by exactly one
    #  dose amplitude, which reads as a broken integrator and was not.
    inner = ~np.any(np.isclose(ts[:, None], (np.arange(12) * T)[None, :]),
                    axis=1)
    worst = float(np.abs(ds[inner] - closed[inner]).max() / closed.max())
    ok.append(worst < 1e-6)
    print(f"      worst relative departure {worst:.3e} over the "
          f"{inner.sum()} samples off a dose time  "
          f"{'PASS' if ok[-1] else 'FAIL'}")

    print("\n  R4. REGRESSION for the bug above: two doses at ONE time must")
    print("      deliver twice the increment, not once.")
    _, _, d1 = forward(THETA_NORMAL, [0.0], t_end=1e-6, per_day=200)
    _, _, d2 = forward(THETA_NORMAL, [0.0, 0.0], t_end=1e-6, per_day=200)
    _, _, d3 = forward(THETA_NORMAL, [0.0], t_end=1e-6, per_day=200,
                       amps=[2 * D_AMP])
    ok.append(abs(d2[0] - 2 * d1[0]) < 1e-9 and abs(d3[0] - 2 * d1[0]) < 1e-9)
    print(f"      one {d1[0]:.1f}, two at one time {d2[0]:.1f}, "
          f"one at double {d3[0]:.1f}  {'PASS' if ok[-1] else 'FAIL'}")

    print("\n  R5. tolerance invariance: the headline peaks must not move when")
    print("      the solver is asked for more accuracy than the report needs.")
    th = 0.05 * THETA_NORMAL
    m0, d0 = steady_start(th)
    a = forward(th, [d for d in np.arange(N_DOSES) * T if d != 0.0],
                t_end=WINDOW, m0=m0, d0=d0, per_day=60)[1].max()
    b = forward(th, [d for d in np.arange(N_DOSES) * T if d != 0.0],
                t_end=WINDOW, m0=m0, d0=d0, per_day=400)[1].max()
    ok.append(abs(a - b) / b < 1e-3)
    print(f"      skipped-dose peak at 5%: {a:.4f} at 60/day against "
          f"{b:.4f} at 400/day, {abs(a-b)/b:.2e} apart  "
          f"{'PASS' if ok[-1] else 'FAIL'}")

    print("\n  R6. the rejoin time is not the window's edge in disguise:")
    print("      doubling the watch must not change it.")
    def rejoin(th, window):
        m0, d0 = steady_start(th)
        base = list(np.arange(int(window / T) + 4) * T)
        rest = [d for d in base if d != 0.0]
        t1 = window + 2 * T
        g1, mu, _ = _on_grid(forward(th, base, t_end=t1, m0=m0, d0=d0,
                                     per_day=60), t1)
        _, mp, _ = _on_grid(forward(th, rest, t_end=t1, m0=m0, d0=d0,
                                    per_day=60), t1)
        gap = np.abs(mp - mu) / np.maximum(mu, 1e-12)
        w = g1 <= window
        idx = np.where(w & (gap > 0.01))[0]
        return float(g1[idx[-1]]) if idx.size else 0.0
    r_a, r_b = rejoin(0.05 * THETA_NORMAL, WINDOW), \
        rejoin(0.05 * THETA_NORMAL, 2 * WINDOW)
    ok.append(abs(r_a - r_b) < 1.0)
    print(f"      5%, skipped: {r_a:.1f} d over {WINDOW:.0f} d of watch "
          f"against {r_b:.1f} d over {2*WINDOW:.0f}  "
          f"{'PASS' if ok[-1] else 'FAIL'}")

    print(f"\n  {sum(ok)} of {len(ok)} checks pass.")
    return all(ok)


def main():
    print(__doc__)
    r_peak, r_upper, attractor = validity_bounds()
    ok = rechecks(r_peak, r_upper)
    if not ok:
        print("\n  THE ADVERSARIAL PASS FAILED. Nothing below is reported.")
        return 1
    s0_setup(r_peak, r_upper, attractor)
    s1_booked(r_peak, r_upper)
    s2_interval_dose(r_peak, r_upper)
    s3_missed(r_peak, r_upper)
    print("\n" + "=" * 84)
    print("Every number above is conditional on a SWEPT capacity, is marked")
    print("where it leaves what the parent supports, and is not dosing")
    print("guidance. See this file's docstring.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
