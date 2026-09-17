"""
Phase 5, steps 1 and 2: simulators for the two published structures, and what
happens when the paper's model is fitted to records they generate.

A NAMING TRAP, stated first because it is easy to walk into. This repository
already has `src/model_A/`, which is the PAPER'S local model, and
`src/model_C/`, the parent. They have nothing to do with the two published
PK/PD structures. Throughout this file:

    PUB_A   the turnover structure: the marker is a STATE whose production the
            therapy inhibits and whose loss is first order. Named, with no
            equations and no parameter values, by Gao (ACoP2025) and by
            Baek et al. 2024 (Nat Commun 15:3804), Methods: "classic indirect
            response models with inhibition of kin ... incorporating
            PCC/hMUT enzyme-induced decreases in organic acid accumulation".
    PUB_B   Attarwala et al. 2023 (Nucleic Acid Ther 33:141-147), Suppl
            Table S1: the marker is an ALGEBRAIC direct-Imax function of the
            enzyme, with no state and no turnover constant of its own.
    PAPER   (theta + Delta) m/(K_m+m) + k_e m = L + b, the repository's own.

WHERE THE LOAD ENTERS IS OUR ASSUMPTION, NOT A TRANSCRIPTION.
Neither published structure has a substrate load in it at all. Both were
fitted to a fixed-diet animal colony, so the load never varied and nothing in
either paper says where it would go. We choose, and we say so in every table:

    PUB_A   the load is the marker's PRODUCTION rate, k_in = L + b. This is
            the least-committal choice available: k_in already IS a
            production rate and the load already IS a rate of substrate
            arrival, so the two are the same object.
    PUB_B   the marker is an instantaneous function of enzyme, so a load can
            only be imposed. We scale the whole marker by (L+b)/L0, i.e. at
            fixed enzyme the marker is proportional to substrate supply. This
            is the weakest assumption that gives the structure a load at all.

WHAT IS PUBLISHED AND WHAT IS OURS, per constant:
    Imax = 0.999                        PUBLISHED, Attarwala Table 1 note m,
                                        dimensionless, transfers unchanged.
    E0/base = 0.0987/60.1               PUBLISHED, Attarwala Table 1, a ratio,
                                        so it transfers unchanged.
    IC50 on OUR enzyme axis              OURS. Attarwala's 37.5 ug/g cannot be
                                        carried across: it is a liver protein
                                        concentration and theta is a flux
                                        capacity, and the untreated PA-mouse
                                        PCC level that would fix the ratio is
                                        not published. We calibrate it to the
                                        repository's own operating point.
    k_out (PUB_A only)                   OURS, and it follows from the
                                        calibration once IC50 is fixed.

THE CALIBRATION HAS TWO ROOTS. Matching one fall at one enzyme level is a
quadratic in IC50 and both roots reproduce the target exactly. We carry BOTH
as separate variants rather than picking one, so that no conclusion here rests
on a branch choice.

Nothing here records an expected result in advance.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy.integrate import solve_ivp

import g_design_common as G
import i1_estimator as I1
import design.run as DR
from design.run import delta

TH, KM, KE, LB = G.TH, G.KM, G.KE, G.LB
T, SIG = G.T, G.SIG
M_BAR = 10.0

IMAX = 0.999                 # Attarwala Table 1, footnote m
R_FLOOR = 0.0987 / 60.1      # Attarwala Table 1: E0 / 3-HPbase


# --------------------------------------------------------------------- shape
def delta_fast(t):
    """
    design.run.delta in closed form, for use inside an ODE right-hand side.

    `delta` sums 200 decaying doses in a Python loop at every call, which an
    integrator makes untenable. The same sum has a closed form: with
    r = exp(lambda T) and K = min(floor(t/T), 199),

        Delta(t) = D exp(-lambda t) (r^(K+1) - 1) / (r - 1),   t >= 0.

    _check_delta_fast() asserts the two agree to 1e-9 before anything uses it,
    so this is an optimisation and not a second model.
    """
    t = np.asarray(t, float)
    lam, D = DR.LAMBDA, DR.D_AMP
    r = np.exp(lam * T)
    K = np.clip(np.floor(t / T), 0.0, 199.0)
    out = D * np.exp(-lam * t) * (r ** (K + 1.0) - 1.0) / (r - 1.0)
    return np.where(t >= 0.0, out, 0.0)


def _check_delta_fast(n=4001):
    tt = np.linspace(0.0, 120.0 * T, n)
    a = np.asarray(delta(tt), float)
    b = delta_fast(tt)
    assert np.max(np.abs(a - b)) < 1e-9 * max(1.0, np.max(np.abs(a))), \
        "delta_fast disagrees with design.run.delta"
    return float(np.max(np.abs(a - b)))


def enzyme(t):
    """Our enzyme axis: the capacity the therapy acts through."""
    return TH + delta_fast(t)


def _inhib(E, ic50):
    """The published inhibition factor, 1 - Imax*E/(IC50+E)."""
    return 1.0 - IMAX * E / (ic50 + E)


def peak_enzyme():
    """The enzyme at the treated steady-state peak, on our axis."""
    return float(TH + delta(100 * T))


def _solve_shape(shape, target_fall, tol=1e-12):
    """
    Every IC50 at which `shape(c)` equals 1 - target_fall.

    `shape(c)` is the structure's OWN ratio marker(E_peak)/marker(E_untreated)
    at a constant load, so each structure is calibrated against itself and no
    constant is carried between them by accident.

    Found by bracketing sign changes on a log grid and bisecting, then
    VERIFIED by re-evaluating at the answer. Two roots are normal: the ratio
    is not monotone in IC50, tending to 1 at both ends. A target above the
    structure's ceiling returns an empty array, never a clipped one.
    """
    r = 1.0 - target_fall
    cs = np.geomspace(1e-6, 1e12, 900001)
    f = shape(cs) - r
    out = []
    for j in np.nonzero(np.sign(f[:-1]) * np.sign(f[1:]) < 0)[0]:
        lo, hi, flo = cs[j], cs[j + 1], f[j]
        for _ in range(300):
            mid = np.sqrt(lo * hi)
            fm = float(shape(mid)) - r
            if np.sign(fm) == np.sign(flo):
                lo, flo = mid, fm
            else:
                hi = mid
            if hi / lo - 1.0 < tol:
                break
        c = float(np.sqrt(lo * hi))
        assert abs(float(shape(c)) - r) < 1e-6, "calibration not verified"
        out.append(c)
    return np.sort(np.array(out))


def _ceiling(shape):
    """The largest fall `shape` can reach, maximised over IC50."""
    cs = np.geomspace(1e-6, 1e12, 900001)
    v = shape(cs)
    j = int(np.argmin(v))
    return float(1.0 - v[j]), float(cs[j])


# ------------------------------------------------------------------- PUB_B
class PubB:
    """Attarwala's structure: algebraic in the enzyme, load imposed by us."""

    name = "PUB_B (Attarwala, algebraic)"

    #: marker(E_hi)/marker(E_lo). The published floor ratio E0/base rides
    #: along and blunts the fall, so this is NOT the inhibition ratio.
    @staticmethod
    def shape(c, E_lo=None, E_hi=None):
        E_lo = TH if E_lo is None else E_lo
        E_hi = peak_enzyme() if E_hi is None else E_hi
        return ((R_FLOOR + _inhib(E_hi, c))
                / (R_FLOOR + _inhib(E_lo, c)))

    def __init__(self, ic50):
        self.ic50 = float(ic50)
        # base fixed so the untreated marker sits at the operating point
        self.base = M_BAR / (R_FLOOR + _inhib(TH, self.ic50))
        self.e0 = R_FLOOR * self.base

    def steady(self, E, load=None):
        load = LB if load is None else load
        return ((self.e0 + self.base * _inhib(np.asarray(E, float), self.ic50))
                * np.asarray(load, float) / LB)

    def marker(self, times, loads, b=0.0, **kw):
        return self.steady(enzyme(times), np.asarray(loads, float) + b)

    def half_life(self):
        return 0.0            # no state: the marker has no memory at all


# ------------------------------------------------------------------- PUB_A
class PubA:
    """Turnover with inhibition of k_in; the load IS k_in."""

    name = "PUB_A (turnover, inhibition of k_in)"

    #: marker(E_hi)/marker(E_lo) at steady state is the INHIBITION ratio, and
    #: k_out cancels out of it. The structure as named has no additive floor.
    @staticmethod
    def shape(c, E_lo=None, E_hi=None):
        E_lo = TH if E_lo is None else E_lo
        E_hi = peak_enzyme() if E_hi is None else E_hi
        return _inhib(E_hi, c) / _inhib(E_lo, c)

    def __init__(self, ic50):
        self.ic50 = float(ic50)
        # k_out fixed so the untreated steady state sits at the operating point
        self.kout = LB * _inhib(TH, self.ic50) / M_BAR

    def steady(self, E, load=None):
        load = LB if load is None else load
        return (np.asarray(load, float)
                * _inhib(np.asarray(E, float), self.ic50) / self.kout)

    def half_life(self):
        return float(np.log(2.0) / self.kout)

    def marker(self, times, loads=None, b=0.0, load_path=None, t_end=None):
        """
        Integrate  dm/dt = k_in(t)*inhib(E(t)) - k_out*m  and read m at
        `times`. k_in is piecewise constant on a DAILY grid, which is the
        resolution a diet record has; `load_path` supplies it. `loads` is
        accepted for signature parity and used only when no path is given.

        INTEGRATED INTERVAL BY INTERVAL between breakpoints. The right-hand
        side JUMPS at every dose time (the enzyme steps by D) and at every day
        boundary (the diet record steps), and an adaptive Runge-Kutta run
        straight through a jump loses its order there: a single solve_ivp call
        over the whole record disagreed with an independent quadrature by 2e-4
        relative on the slow branch and 5e-2 on the fast one. Restarting at
        each break removes that, and lets the steps inside a smooth stretch
        grow. k3_recheck.check_5 is what caught it.

        The run starts ten marker half-lives before the first visit, from the
        untreated steady state, so no visit is reading a start transient.
        """
        times = np.asarray(times, float)
        t0 = float(times.min()) - 10.0 * self.half_life()
        t_end = float(times.max()) + 1.0 if t_end is None else t_end
        days = np.arange(np.floor(t0) - 1.0, np.ceil(t_end) + 2.0)
        kin = (np.full(days.size, LB) if load_path is None
               else np.asarray(load_path(days), float)) + b

        def kin_at(t):
            j = int(np.clip(np.searchsorted(days, t, "right") - 1,
                            0, days.size - 1))
            return kin[j]

        doses = np.arange(0.0, t_end + T, T)
        brk = np.unique(np.concatenate([days, doses, times, [t0, t_end]]))
        brk = brk[(brk >= t0 - 1e-12) & (brk <= t_end + 1e-12)]

        out = {}
        m = kin_at(t0) * _inhib(float(enzyme(t0)), self.ic50) / self.kout
        out[round(float(brk[0]), 9)] = m
        for lo, hi in zip(brk[:-1], brk[1:]):
            if hi - lo <= 1e-12:
                continue
            k = kin_at(0.5 * (lo + hi))

            def rhs(t, y):
                return [k * _inhib(float(enzyme(t)), self.ic50)
                        - self.kout * y[0]]

            sol = solve_ivp(rhs, (lo, hi), [m], rtol=1e-11, atol=1e-14,
                            method="LSODA")
            m = float(sol.y[0, -1])
            out[round(float(hi), 9)] = m
        return np.asarray([out[round(float(t), 9)] for t in times], float)


# ------------------------------------------------------------------- PAPER
class Paper:
    """The repository's own model, as the generator, for the control row."""

    name = "PAPER (the repository's own)"

    def __init__(self, ic50=None):
        self.ic50 = np.nan

    def marker(self, times, loads, b=0.0, **kw):
        return I1.levels(np.asarray(times, float), np.asarray(loads, float),
                         b=b)

    def half_life(self):
        return 0.0             # quasi-steady by construction: no memory


# --------------------------------------------------------------------------
def daily_loads(load_cv, seed, lo, hi):
    """A lognormal daily diet record over [lo, hi], and a reader for it."""
    days = np.arange(np.floor(lo) - 1.0, np.ceil(hi) + 2.0)
    if not load_cv:
        vals = np.full(days.size, LB)
    else:
        rng = np.random.default_rng(seed)
        s2 = np.log1p(load_cv ** 2)
        vals = LB * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), days.size))

    def read(t):
        j = np.clip(np.searchsorted(days, np.asarray(t, float), "right") - 1,
                    0, days.size - 1)
        return vals[j]

    return read


# --------------------------------------------------------------------------
def daily_loads(load_cv, seed, lo, hi):
    """A lognormal DAILY diet record over [lo, hi], and a reader for it.

    The path is daily because that is the resolution a diet record has, and
    because it is what makes PUB_A structurally different: a turnover marker
    integrates the load history, so what it reads on a visit day is not the
    load on that day. PUB_B and the paper's model both read the day's load
    and nothing else.
    """
    days = np.arange(np.floor(lo) - 2.0, np.ceil(hi) + 3.0)
    if not load_cv:
        vals = np.full(days.size, LB)
    else:
        rng = np.random.default_rng(seed)
        s2 = np.log1p(load_cv ** 2)
        vals = LB * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), days.size))

    def read(t):
        j = np.clip(np.searchsorted(days, np.asarray(t, float), "right") - 1,
                    0, days.size - 1)
        return vals[j]

    return read


def generators(target_fall):
    """Every generator carried in the comparison, at one calibration target."""
    out = [("PAPER", Paper())]
    for cls, tag in ((PubA, "PUB_A"), (PubB, "PUB_B")):
        roots = _solve_shape(cls.shape, target_fall)
        if roots.size == 0:
            ceil, _ = _ceiling(cls.shape)
            out.append((f"{tag} (UNREACHABLE: ceiling {100*ceil:.2f}%)", None))
            continue
        for c, br in zip(roots, ("low", "high")):
            out.append((f"{tag} IC50 {c:.4g} ({br})", cls(c)))
    return out


def one_fit(gen, times, load_cv, seed, noiseless=False,
            load_frac=I1.LOAD_FRAC, cv=I1.CV):
    """
    Generate a record from `gen`, fit the PAPER'S model to it, and report.

    The ordinate and its noise are the repository's own, i1_estimator's
    `one_record` with eiv=False: the recorded load carries a relative error
    and the assay error is placed on the ordinate as b_loc*cv*m_bar.
    """
    times = np.asarray(times, float)
    path = daily_loads(load_cv, seed, times.min() - 400.0, times.max() + 2.0)
    L = path(times)
    m = gen.marker(times, L, load_path=path)
    inc = np.asarray(delta(times), float)
    if noiseless:
        s = L.copy()
    else:
        rng = np.random.default_rng(7_000_003 * seed + 11)
        s = (L * (1.0 + load_frac * rng.standard_normal(times.size))
             + I1.BLOC * cv * M_BAR * rng.standard_normal(times.size))
    fit = I1.varpro(m, inc, s)[0]
    fit = I1.polish(fit, m, inc, s)[0]
    r = ((fit["theta"] + inc) * m / (fit["Km"] + m)
         + fit["ke"] * m + fit["b"] - s)
    return dict(fit, rms=float(np.sqrt(np.mean(r ** 2)) / SIG),
                m=m, s=s, inc=inc)


def km_flatness(m, inc, s, rise=1.0):
    """
    The width of the K_m profile: the range of K_m over which the penalised
    objective sits within `rise` of its minimum.

    This is P5's fibre made measurable. At a constant recorded load the
    paper's own model is exactly singular in this direction and the width is
    the whole grid. Whether a record generated by a DIFFERENT structure is
    also flat is the question step 2 is really asking, because the paper's
    obstruction is claimed to survive a change of structure.
    """
    grid = np.geomspace(0.3, 3e4, 4001)
    J, _ = I1.varpro_profile(m, inc, s, grid)
    j = int(np.argmin(J))
    ok = J <= J[j] + rise
    lo, hi = grid[ok].min(), grid[ok].max()
    return float(lo), float(hi), float(grid[j]), bool(ok[0] or ok[-1])


def sd_at_fit(fit, m, inc):
    """
    sd(ln theta) the fitter would report: the information filter of
    g_design_common, but with the columns built at the FITTED parameters and
    on the GENERATED levels, not at the paper's own truth.

    Under misspecification there is no true theta to report a sd about. This
    is the number a person running the paper's method would write down, which
    is the quantity the design advice is stated in.
    """
    th, Km, ke = fit["theta"], fit["Km"], fit["ke"]
    g = m / (Km + m)
    Q = np.stack([th * g, -(th + inc) * Km * m / (Km + m) ** 2,
                  np.ones_like(m), ke * m], axis=1)
    F = Q.T @ Q / SIG ** 2 + np.diag([0.0, 0.0, 1.0 / I1.S_B ** 2, 1.0 / 0.8 ** 2])
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def _stat(v):
    v = np.asarray([x for x in v if np.isfinite(x)], float)
    if v.size == 0:
        return np.nan, np.nan, np.nan
    return float(np.median(v)), float(v.min()), float(v.max())


def step2(target_fall, n=12, seeds=range(8),
          cvs=(0.0, 0.02, 0.10), out=print):
    """The misspecification table at one calibration target."""
    gens = generators(target_fall)
    out("")
    out(f"  CALIBRATION TARGET: a {100*target_fall:.0f}% marker fall from the "
        f"untreated level")
    out(f"  to the treated steady-state peak (enzyme {TH:.0f} -> "
        f"{peak_enzyme():.0f}).")
    for cls, tag in ((PubA, "PUB_A"), (PubB, "PUB_B")):
        ceil, c_at = _ceiling(cls.shape)
        out(f"    {tag} ceiling on that fall, over ALL IC50: "
            f"{100*ceil:.2f}% (at IC50 {c_at:.4g})")
    for nm, g in gens:
        if g is None:
            out(f"    {nm}")
        elif hasattr(g, "kout"):
            out(f"    {nm}: k_out {g.kout:.4g} /d, marker half-life "
                f"{g.half_life():.3g} d")
        elif hasattr(g, "base"):
            out(f"    {nm}: base {g.base:.4g}, E0 {g.e0:.4g}, no marker state")
    out("")
    out("  Each cell: RMS residual of the paper's best fit, in units of the")
    out("  repository's own sigma, with NO observation noise at all, so it is")
    out("  structure and nothing else; then theta_hat/3000 from the same fit.")
    out("  Median over 8 daily diet records; [min, max] on the residual.")
    out("  A theta_hat of 3.3e-13 is the estimator's own lower clip: the fit")
    out("  wants theta <= 0 and cannot go there.")
    out("")
    head = f"      {'generator':30s}{'schedule':26s}"
    for c in cvs:
        head += f"{'CV ' + str(int(100 * c)) + '%':>26s}"
    out(head)
    for name, g in gens:
        if g is None:
            continue
        first = True
        for sname, fn in G.SCHEDULES.items():
            ts = fn(n)
            row = f"      {(name if first else ''):30s}{sname:26s}"
            first = False
            for cv in cvs:
                ss = [0] if cv == 0.0 else list(seeds)
                rows = [one_fit(g, ts, cv, s, noiseless=True) for s in ss]
                rm, rlo, rhi = _stat([r["rms"] for r in rows])
                tm, _, _ = _stat([r["theta"] / TH for r in rows])
                cell = (f"{rm:6.3f} {tm:9.3g}" if len(ss) == 1
                        else f"{rm:6.3f}[{rlo:5.2f},{rhi:6.2f}] {tm:9.3g}")
                row += f"{cell:>26s}"
            out(row)
        out("")


def step2_noisy(target_fall, n=12, seeds=range(8), cvs=(0.0, 0.02, 0.10),
                out=print):
    """The same matrix WITH the repository's observation noise."""
    gens = generators(target_fall)
    out("  With the repository's observation noise on the ordinate, which is")
    out(f"  what a real record has. sigma = {SIG:.4f}. Each cell: RMS/sigma,")
    out("  then sd(ln theta) as the FITTER would report it, both medians over")
    out("  8 seeds. The seeds carry the diet record AND the noise draw.")
    out("  'sing' means EVERY seed returned a singular information matrix at")
    out("  its own fit, so the fitter cannot report a precision for theta at")
    out("  all. That is a result, not a missing number.")
    out(f"  A well-specified fit of n={n} points for 4 parameters has an")
    out(f"  expected RMS/sigma of sqrt((n-4)/n) = {np.sqrt((n-4)/n):.3f}, which is the")
    out("  figure to read the residual column against.")
    out("")
    head = f"      {'generator':30s}{'schedule':26s}"
    for c in cvs:
        head += f"{'CV ' + str(int(100 * c)) + '%':>22s}"
    out(head)
    for name, g in gens:
        if g is None:
            continue
        first = True
        for sname, fn in G.SCHEDULES.items():
            ts = fn(n)
            row = f"      {(name if first else ''):30s}{sname:26s}"
            first = False
            for cv in cvs:
                rows = [one_fit(g, ts, cv, s) for s in seeds]
                rm, _, _ = _stat([r["rms"] for r in rows])
                sd, _, _ = _stat([sd_at_fit(r, r["m"], r["inc"]) for r in rows])
                sds = "    sing" if not np.isfinite(sd) else f"{sd:8.4f}"
                row += f"{f'{rm:7.3f} {sds}':>22s}"
            out(row)
        out("")


def step2_fibre(target_fall, n=12, out=print):
    """Does the fibre survive a change of structure, at a constant load?"""
    out("  THE FIBRE AT A CONSTANT RECORDED LOAD (load CV = 0, no noise).")
    out("  The width of the K_m profile within 1 unit of the minimum of the")
    out("  same penalised objective the estimator minimises. The paper's own")
    out("  model is EXACTLY singular here, so its width is the whole grid;")
    out("  'edge' means the flat region runs off the end of the scan.")
    out("")
    out(f"      {'generator':34s}{'schedule':26s}"
        f"{'K_m at min':>12s}{'flat from':>12s}{'to':>12s}  edge")
    for name, g in generators(target_fall):
        if g is None:
            continue
        for sname, fn in G.SCHEDULES.items():
            ts = fn(n)
            r = one_fit(g, ts, 0.0, 0, noiseless=True)
            lo, hi, kmin, edge = km_flatness(r["m"], r["inc"], r["s"])
            out(f"      {name:34s}{sname:26s}{kmin:12.4g}{lo:12.4g}"
                f"{hi:12.4g}  {'yes' if edge else 'no'}")
    out("")


def main():
    print("=" * 78)
    print("PHASE 5 STEP 1: THE TWO PUBLISHED STRUCTURES AS SIMULATORS")
    print("=" * 78)
    print(__doc__.strip())
    print()
    print(f"  delta_fast checked against design.run.delta: max abs difference "
          f"{_check_delta_fast():.3g}")
    print(f"  operating point m_bar = {M_BAR}, L0 = {LB:.4f}, "
          f"sigma = {SIG:.4f}")
    print()
    print("=" * 78)
    print("PHASE 5 STEP 2: FIT THE PAPER'S MODEL TO RECORDS THEY GENERATE")
    print("=" * 78)
    for tgt in (0.60, 0.30, 0.20):
        print("-" * 78)
        step2(tgt)
    for tgt in (0.60, 0.30):
        print("-" * 78)
        print(f"  NOISY PASS, {100*tgt:.0f}% calibration")
        print("-" * 78)
        step2_noisy(tgt)
    for tgt in (0.60, 0.30):
        print("-" * 78)
        print(f"  FIBRE, {100*tgt:.0f}% calibration")
        print("-" * 78)
        step2_fibre(tgt)


if __name__ == "__main__":
    main()
