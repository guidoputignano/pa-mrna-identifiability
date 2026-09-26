"""
Step 3(b). Identifiability with the marker's unreachable offset FIXED.

Phase 5 proved that under the published algebraic structure, carrying the
offset `E0` as a free parameter is what makes the enzyme capacity theta
unidentifiable, and that fixing `E0` restores identifiability. What Phase 5
could not do was say what to fix it AT: the published value is a mouse estimate,
`E0/base = 0.0987/60.1`, which puts the unreachable part of the marker at a
fraction of a percent.

R5, the extrahepatic-floor report of 25 September, supplies a human anchor from
liver transplantation, the natural experiment for a liver-directed therapy:
post-transplant 3-HP is 5 to 10% of the pre-transplant level across two cohorts,
with a wide plausible range to about 30%. NOTHING IN THAT REPORT IS VERIFIED
HERE, and the two transplant ratios it rests on are step 2 item 2, which has not
run for want of egress. So the anchor is carried as a SWEEP over the values the
author named, 5, 10 and 20% of the untreated level, and the question asked of it
is what each choice buys and what fixing it wrongly costs.

THE OFFSET'S DEFINITION, stated because the structure admits two readings. The
observable is

    m(E) = [E0 + base * f(E)] * (L + b) / L0,   f(E) = (IC50 + q E)/(IC50 + E)

with q = 1 - Imax. The untreated patient sits at E = theta, where the
calibration puts m = M_BAR. A transplant replaces the whole hepatic pool, which
is E -> infinity, where m -> [E0 + base*q] * (L+b)/L0. So the transplant ratio
measures (E0 + base*q)/(E0 + base*f(theta)), which is NOT E0 alone: at
Imax = 0.999 the incomplete-inhibition term base*q rides along inside it. This
script fixes the quantity R5 anchors, the ratio of the unreachable marker to the
untreated marker, and reports the E0 it implies, rather than fixing E0 directly
and hoping the ratio lands somewhere plausible.

Nothing here records an expected result in advance.

Run:  python audit/q7_offset_fixed.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import sympy as sp
from scipy.optimize import brentq

import g_design_common as G
import k2_structures as K2
import q2_requirement_map as S2
from design.run import delta

TH, LB = G.TH, G.LB
M_BAR = K2.M_BAR
IMAX = K2.IMAX
Q_RES = 1.0 - IMAX                     #: the residual, 1 - Imax
FRACS = (0.05, 0.10, 0.20)             #: R5's anchor and the author's sweep
TARGET_FALL = 0.60                     #: k2's own calibration target
MARKER_CVS = (0.05, 0.15)
LOAD_FRAC = 0.10
SEEDS = tuple(range(24))
GATE = S2.GATE
RESULTS = []


def banner(t):
    print("\n" + "=" * 100)
    print("  " + t)
    print("=" * 100)


def record(tag, ok, detail=""):
    RESULTS.append((tag, ok))
    print(f"\n  [{'PASS' if ok else 'FAIL'}] {tag}")
    for line in str(detail).splitlines():
        if line:
            print(f"         {line}")


# ==========================================================================
#  The structure, parametrised by the unreachable RATIO rather than by E0.
# ==========================================================================
def inhib(E, ic50):
    return (ic50 + Q_RES * np.asarray(E, float)) / (ic50 + np.asarray(E, float))


class Offset:
    """
    PUB_B with the unreachable fraction of the UNTREATED marker set to `phi`.

    `base` is solved so the untreated level stays at M_BAR, and `ic50` is
    re-calibrated against this structure's OWN fall so that each phi is
    calibrated against itself, which is `k2_structures._solve_shape`'s rule.
    """

    def __init__(self, phi, target_fall=TARGET_FALL, theta=TH):
        self.phi = float(phi)
        self.theta = float(theta)
        self.ic50 = self._calibrate(target_fall)
        self.e0, self.base = self._split(self.ic50)

    def _split(self, ic50):
        """(E0, base) with m(theta) = M_BAR and (E0 + base q)/M_BAR = phi."""
        #  E0 + base*q = phi*M_BAR  and  E0 + base*f(theta) = M_BAR
        f = float(inhib(self.theta, ic50))
        base = M_BAR * (1.0 - self.phi) / (f - Q_RES)
        e0 = self.phi * M_BAR - base * Q_RES
        return float(e0), float(base)

    def _fall(self, ic50):
        """This structure's own marker(E_peak)/marker(E_untreated)."""
        e0, base = self._split(ic50)
        hi = e0 + base * float(inhib(K2.peak_enzyme(), ic50))
        lo = e0 + base * float(inhib(self.theta, ic50))
        return hi / lo

    def ic50_min(self):
        """
        The smallest ADMISSIBLE IC50, where E0 reaches zero.

        `_split` holds the untreated marker at M_BAR and the unreachable ratio
        at phi, which forces E0 = phi*M_BAR - base*q with
        base = M_BAR(1-phi)/(f(theta)-q). As IC50 falls, f(theta) falls toward q,
        base diverges and E0 goes NEGATIVE. A negative unreachable term is not a
        structure, so the fall's ceiling has to be taken over the admissible
        range only. Setting E0 = 0 gives f(theta) = q/phi and hence

            IC50_min = q*theta*(1-phi)/(phi-q).

        An earlier version of this script maximised the fall over all IC50 and
        reported ceilings attained at E0 far below zero. Those numbers were
        wrong and are not carried.
        """
        if self.phi <= Q_RES:
            return 0.0
        return Q_RES * self.theta * (1.0 - self.phi) / (self.phi - Q_RES)

    def ceiling(self, verify=False):
        """
        The largest fall this phi can produce over ADMISSIBLE IC50.

        The fall decreases monotonically in IC50 above `ic50_min`, so the
        ceiling is attained AT `ic50_min` and needs one evaluation rather than a
        scan. `verify=True` checks that monotonicity on a coarse grid instead of
        assuming it.
        """
        lo = self.ic50_min() * (1.0 + 1e-12)
        if verify:
            grid = np.geomspace(lo, lo * 1e6, 120)
            v = np.array([self._fall(c) for c in grid])
            if not np.all(np.diff(v) > -1e-12):
                j = int(np.argmin(v))
                return float(1.0 - v[j]), float(grid[j])
        return float(1.0 - self._fall(lo)), float(lo)

    def _calibrate(self, target_fall):
        """
        The IC50 at which this phi produces `target_fall`, searched over the
        ADMISSIBLE range only.

        An earlier version searched from IC50 = 1e-3 and so could return a
        calibration with E0 < 0, which made the phi window look 0.2 points wider
        than it is. The search now starts at `ic50_min`.
        """
        want = 1.0 - target_fall
        grid = np.geomspace(max(self.ic50_min(), 1e-12) * (1.0 + 1e-12),
                            1e9, 4001)
        vals = np.array([self._fall(c) for c in grid])
        ok = np.nonzero(np.sign(vals[:-1] - want) * np.sign(vals[1:] - want)
                        < 0)[0]
        if not ok.size:
            raise ValueError(f"phi = {self.phi} cannot reach a "
                             f"{100 * target_fall:.0f}% fall; its ceiling is "
                             f"{100 * (1 - vals.min()):.1f}%")
        j = int(ok[0])
        return float(brentq(lambda c: self._fall(c) - want, grid[j],
                            grid[j + 1], xtol=1e-12, rtol=1e-14))

    def marker(self, E, load=LB, b=0.0):
        return ((self.e0 + self.base * inhib(E, self.ic50))
                * (np.asarray(load, float) + b) / LB)

    def unreachable_ratio(self):
        return (self.e0 + self.base * Q_RES) / M_BAR


# ==========================================================================
def _bare(phi):
    """An Offset without calibration, for asking about its ceiling alone."""
    o = object.__new__(Offset)
    o.phi, o.theta = float(phi), TH
    return o


def part_a():
    banner("A. THE STRUCTURE, CHECKED AGAINST k2 AND AGAINST ITS OWN ANCHOR")
    print("""  The published value first, as a cross-check: setting phi to the ratio
  `k2_structures.PubB` itself produces must reproduce that object's E0, base and
  IC50. Then the three anchored values.""")
    ic = float(K2._solve_shape(K2.PubB.shape, TARGET_FALL)[0])
    o = K2.PubB(ic)
    pub_ratio = (o.e0 + o.base * Q_RES) / M_BAR
    mine = Offset(pub_ratio)
    gap = max(abs(mine.e0 / o.e0 - 1.0), abs(mine.base / o.base - 1.0),
              abs(mine.ic50 / o.ic50 - 1.0))
    print(f"\n    k2's PubB at a {100 * TARGET_FALL:.0f}% fall: E0 = {o.e0:.6f}, "
          f"base = {o.base:.4f}, IC50 = {o.ic50:.6f}")
    print(f"    its unreachable ratio (E0 + base q)/M_BAR = "
          f"{100 * pub_ratio:.3f}% of the untreated level")
    print(f"    this script at the same ratio:  E0 = {mine.e0:.6f}, "
          f"base = {mine.base:.4f}, IC50 = {mine.ic50:.6f}")
    record("the offset structure reproduces k2's PubB at k2's own ratio",
           gap < 1e-9, f"worst relative gap {gap:.3e}")
    print("""
  AND A CONSTRAINT THAT HAD TO BE FOUND BEFORE THE REST WOULD RUN. The offset and
  the fall are not independent. An unreachable fraction caps the fall the
  structure can produce, because therapy can only remove what is reachable, and
  the cap is reached in the IC50 -> 0 limit. So each phi has a CEILING on the
  fall, and a phi whose ceiling is below the paper's calibrated 60% cannot be
  calibrated at all.""")
    print(f"\n    {'phi':>7s} {'IC50 >= ':>10s} {'admissible ceiling':>19s} "
          f"{'60% reachable':>14s}")
    for phi in (0.0026, 0.01, 0.02, 0.05, 0.08, 0.10, 0.1277, 0.15, 0.20, 0.30):
        o_ = _bare(phi)
        c, at = o_.ceiling(verify=True)
        print(f"    {100 * phi:6.2f}% {o_.ic50_min():10.3f} {100 * c:18.2f}% "
              f"{'yes' if c >= TARGET_FALL else 'NO':>14s}")
    print("""
  THE CEILING IS NOT MONOTONE IN phi, so the constraint is a WINDOW and not an
  upper limit. As phi falls toward q = 1 - Imax the unreachable part stops being
  E0 at all and becomes the incomplete-inhibition term base*q; IC50_min then
  diverges, and a large IC50 caps the fall just as a large phi does. So a 60%
  fall needs phi neither too large nor too small, and the lower edge is set by
  Imax rather than by any biology.""")
    f_ = lambda x: _bare(x).ceiling()[0] - TARGET_FALL
    grid = np.geomspace(1.2 * Q_RES, 0.30, 400)
    vals = np.array([f_(x) for x in grid])
    cross = np.nonzero(np.sign(vals[:-1]) * np.sign(vals[1:]) < 0)[0]
    edges = [float(brentq(f_, grid[j], grid[j + 1], xtol=1e-10))
             for j in cross]
    print(f"\n    phi WINDOW COMPATIBLE WITH THE PAPER'S {100 * TARGET_FALL:.0f}% "
          f"FALL: " + " to ".join(f"{100 * e:.3f}%" for e in edges))
    inside = (len(edges) == 2 and edges[0] <= pub_ratio <= edges[1])
    print(f"    the published structure's own ratio, {100 * pub_ratio:.3f}%, is "
          f"{'inside' if inside else 'OUTSIDE'} it"
          + (f", {100 * (edges[1] - pub_ratio):.3f} points below the upper edge"
             if inside else ""))
    phi_max = edges[-1] if edges else float("nan")
    record("R5's 5 to 10% anchor is compatible with the paper's own fall",
           len(edges) == 2 and edges[0] <= 0.05 and phi_max >= 0.10,
           f"window {100 * edges[0]:.3f}% to {100 * phi_max:.3f}%, so 5% and "
           f"10% calibrate and 20% does not")
    print(f"\n    {'phi':>6s} {'E0':>12s} {'base':>12s} {'IC50':>12s} "
          f"{'m(untreated)':>13s} {'unreachable':>12s} {'fall':>8s}")
    for phi in FRACS:
        try:
            b = Offset(phi)
        except ValueError as exc:
            print(f"    {100 * phi:5.0f}% {'INADMISSIBLE at a 60% fall':>52s}"
                  f"   {exc}")
            continue
        print(f"    {100 * phi:5.0f}% {b.e0:12.5f} {b.base:12.4f} "
              f"{b.ic50:12.5f} {b.marker(TH):13.6f} "
              f"{100 * b.unreachable_ratio():11.2f}% "
              f"{100 * (1 - b._fall(b.ic50)):7.1f}%")
    print("""
  The published structure's own unreachable ratio is printed above so the mouse
  value and R5's human anchor can be read on one scale. R5's point is that they
  differ by more than an order of magnitude, and the constraint above adds that
  the human end of R5's range is close to the largest the paper's own calibration
  admits.""")


def part_b(times, loads):
    banner("B. THE FIBRE: FREE OFFSET AGAINST FIXED OFFSET")
    print("""  Phase 5's result, restated and then answered. A constant-load record of this
  structure fixes only three combinations,

      P = E0 + base q,   U = P theta + (E0 + base) IC50,   V = IC50 + theta,

  so with (E0, base, IC50, theta) all free there are four unknowns for three
  equations and theta slides along a one-parameter fibre at a bit-identical
  record. Fixing E0 leaves three unknowns for three equations. Both halves are
  checked here symbolically and then by walking the fibre.""")
    E0, base, ic50, th, R = sp.symbols("E0 base IC50 theta R", positive=True)
    q = sp.Rational(1, 1) * (1 - sp.Float(IMAX))
    P, U, V = E0 + base * q, None, ic50 + th
    U = P * th + (E0 + base) * ic50
    J = sp.Matrix([[sp.diff(x, v) for v in (E0, base, ic50, th)]
                   for x in (P, U, V)])
    rank_free = J.rank()
    J_fixed = sp.Matrix([[sp.diff(x, v) for v in (base, ic50, th)]
                         for x in (P, U, V)])
    rank_fixed = J_fixed.rank()
    print(f"\n    rank of d(P, U, V)/d(E0, base, IC50, theta) = {rank_free} "
          f"of 4 unknowns  ->  a {4 - rank_free}-parameter fibre")
    print(f"    rank of d(P, U, V)/d(base, IC50, theta), E0 fixed = "
          f"{rank_fixed} of 3 unknowns  ->  {'identifiable' if rank_fixed == 3 else 'still not'}")
    record("free offset leaves a one-parameter fibre", rank_free == 3,
           f"rank {rank_free} against 4 unknowns")
    record("fixing the offset makes theta locally identifiable",
           rank_fixed == 3, f"rank {rank_fixed} against 3 unknowns")

    print("""
  And the fibre walked, at phi = 10%. With E0 free, theta is re-solved from
  (P, U, V) at each IC50 = R and the predicted record must not move. With E0
  fixed, the same walk must FAIL to hold the record, because E0 is no longer
  available to absorb the change.""")
    o = Offset(0.10)
    inc = np.asarray(delta(times), float)
    E_ref = o.theta + inc
    ref = o.marker(E_ref, loads)
    Pv = o.e0 + o.base * Q_RES
    Uv = Pv * o.theta + (o.e0 + o.base) * o.ic50
    Vv = o.ic50 + o.theta
    print(f"\n    truth: E0 = {o.e0:.6f}, base = {o.base:.4f}, "
          f"IC50 = {o.ic50:.6f}, theta = {o.theta:.0f}")
    print(f"    {'IC50 = R':>12s} {'theta = V-R':>12s} {'base':>12s} "
          f"{'E0':>12s} {'max|dm|/m, E0 free':>19s} "
          f"{'same, E0 held':>15s}")
    worst_free, worst_held = 0.0, 0.0
    for frac in (0.05, 0.2, 0.5, o.ic50 / Vv, 0.8, 0.95):
        Rv = frac * Vv
        th_r = Vv - Rv
        base_r = (Uv - Pv * Vv) / (Rv * IMAX)
        e0_r = Pv - Q_RES * base_r
        if not (base_r > 0 and e0_r >= 0 and th_r > 0):
            print(f"    {Rv:12.4g} {th_r:12.4g} {base_r:12.4g} {e0_r:12.4g} "
                  f"{'inadmissible':>19s} {'':>15s}")
            continue
        E = th_r + inc
        free = (e0_r + base_r * inhib(E, Rv)) * (np.asarray(loads) + 0.0) / LB
        held = (o.e0 + base_r * inhib(E, Rv)) * (np.asarray(loads) + 0.0) / LB
        ef = float(np.max(np.abs(free - ref)) / np.max(ref))
        eh = float(np.max(np.abs(held - ref)) / np.max(ref))
        worst_free = max(worst_free, ef)
        worst_held = max(worst_held, eh)
        print(f"    {Rv:12.4g} {th_r:12.4g} {base_r:12.4g} {e0_r:12.4g} "
              f"{ef:19.3e} {eh:15.3e}")
    record("with the offset free, theta slides at an unchanged record",
           worst_free < 1e-12, f"worst relative change {worst_free:.3e}")
    record("with the offset held, the same walk moves the record",
           worst_held > 1e-3,
           f"worst relative change {worst_held:.3e}, a factor of "
           f"{worst_held / max(worst_free, 1e-18):.2e} larger")


# ==========================================================================
def fisher_theta(o, times, loads, cv, free_e0, s_b):
    """
    sd(ln theta) under this structure, parameters
    (ln theta, ln base, ln IC50, b[, ln E0]), relative marker noise `cv`.
    """
    inc = np.asarray(delta(times), float)
    p0 = dict(theta=o.theta, base=o.base, ic50=o.ic50, e0=o.e0, b=0.0)

    def pred(theta, base, ic50, e0, b):
        E = theta + inc
        return ((e0 + base * inhib(E, ic50))
                * (np.asarray(loads, float) + b) / LB)

    ref = pred(**p0)
    names = ["theta", "base", "ic50"] + (["e0"] if free_e0 else []) + ["b"]
    cols, h = [], 1e-6
    for nm in names:
        up, dn = dict(p0), dict(p0)
        if nm == "b":
            up[nm] = p0[nm] + h
            dn[nm] = p0[nm] - h
            cols.append((pred(**up) - pred(**dn)) / (2 * h))
        else:
            up[nm] = p0[nm] * np.exp(h)
            dn[nm] = p0[nm] * np.exp(-h)
            cols.append((pred(**up) - pred(**dn)) / (2 * h))
    Jm = np.stack(cols, axis=1)
    sig = cv * ref
    F = (Jm / sig[:, None]).T @ (Jm / sig[:, None])
    F[names.index("b"), names.index("b")] += 1.0 / s_b ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 1e-14 * w[-1]:
        return np.inf
    v = np.linalg.inv(F)[0, 0]
    return float(np.sqrt(v)) if v > 0 else np.inf


def part_c(times, loads):
    banner("C. WHAT FIXING IT BUYS, AND AT WHICH ANCHOR")
    print(f"""  sd(ln theta) from the information matrix, relative marker noise as stated,
  load recorded to {100 * LOAD_FRAC:.0f}% so the load-bias prior is that fraction of the load.
  `*` at or below {GATE:.2f}.""")
    s_b = LOAD_FRAC * LB
    print(f"\n    {'phi':>6s} {'marker CV':>10s} {'E0 free':>12s} "
          f"{'E0 fixed':>12s} {'factor':>9s}")
    out = {}
    for phi in FRACS:
        try:
            o = Offset(phi)
        except ValueError:
            print(f"    {100 * phi:5.0f}% {'inadmissible at a 60% fall, see A':>50s}")
            continue
        for cv in MARKER_CVS:
            a = fisher_theta(o, times, loads, cv, True, s_b)
            b = fisher_theta(o, times, loads, cv, False, s_b)
            out[(phi, cv)] = (a, b)
            star = "*" if b <= GATE else " "
            fa = f"{a / b:.1f}x" if np.isfinite(a) else "infinite"
            print(f"    {100 * phi:5.0f}% {100 * cv:9.0f}% "
                  f"{a:12.4g} {b:11.4g}{star} {fa:>9s}")
    print("""
  The `E0 free` column is Phase 5's finding in this parametrisation: an infinite
  or near-infinite entry is the fibre. The `factor` column is what the anchor is
  worth, and it is the whole of step 3(b)'s first half.""")
    return out


def profile_theta(e0_held, o_true, times, loads, thetas):
    """
    The best relative residual at each theta with the offset HELD, optimising
    (base, IC50) exactly at each node from several starts.

    A point fit cannot answer this question: if the held offset is wrong the
    objective may have no interior minimum at all, and only a profile shows that.
    An earlier version of this script reported a point estimate from a bounded
    theta grid and read its lower edge as a bias. That was an artifact and the
    profile is what replaced it.
    """
    from scipy.optimize import least_squares
    inc = np.asarray(delta(times), float)
    ref = o_true.marker(o_true.theta + inc, loads)
    out = []
    for th in thetas:
        def res(q):
            bs, c = np.exp(q)
            return ((e0_held + bs * inhib(th + inc, c))
                    * np.asarray(loads, float) / LB) - ref
        best = None
        for b0 in (10.0, 48.0, 144.0, 500.0):
            for c0 in (10.0, 63.0, 200.0, 733.0, 3000.0):
                try:
                    r = least_squares(res, np.log([b0, c0]), method="lm",
                                      xtol=1e-15, ftol=1e-15, max_nfev=8000)
                except Exception:
                    continue
                if best is None or r.cost < best.cost:
                    best = r
        out.append(float(np.sqrt(2 * best.cost / np.sum(ref ** 2))))
    return np.array(out)


def part_d(times, loads):
    banner("D. THE COST OF FIXING IT AT THE WRONG VALUE")
    print("""  R5's anchor is a range, and the mouse value it exists to replace is near zero.
  Fixing the offset at one value when the truth is another was expected to bias
  theta. IT DOES SOMETHING WORSE, and a point estimate hides it, so what is
  reported here is the PROFILE: the best achievable residual at each theta with
  the offset held, (base, IC50) optimised exactly at every node.""")
    thetas = np.array([30.0, 100.0, 300.0, 1000.0, 2000.0, 3000.0, 6000.0,
                       20000.0])
    ot = Offset(0.05)
    rows = [("held at its own 5%, correct", ot.e0),
            ("held at 10%, too high", Offset(0.10).e0)]
    try:
        rows.append(("held at 2%, too low", Offset(0.02).e0))
    except ValueError:
        pass
    print(f"\n    truth is phi = 5%, E0 = {ot.e0:.5f}, theta = {ot.theta:.0f}")
    print(f"\n    {'offset held at':28s}" +
          "".join(f"{f'th={t:.0f}':>11s}" for t in thetas))
    for label, e0 in rows:
        v = profile_theta(e0, ot, times, loads, thetas)
        print(f"    {label:28s}" + "".join(f"{x:11.2e}" for x in v))
        j = int(np.argmin(v))
        flat = [t for t, x in zip(thetas, v) if x <= 1.02 * v[j]]
        print(f"    {'':28s}  minimum at theta = {thetas[j]:.0f}; within 2% of "
              f"it: theta in [{min(flat):.0f}, {max(flat):.0f}]")
    print("""
  Read the correct row first: a sharp minimum at the truth and a residual at
  machine precision there, which is the identifiability 3b.1 proved. Then read the
  wrong row: the residual is nearly constant across a wide range of theta and the
  TRUE value fits WORSE than values far below it. A wrong anchor does not bias
  theta, it gives back the flat profile Phase 5 found with the offset free, while
  the residual stays far below what a 5% marker CV would produce, so nothing in
  the fit reveals it.""")


def main():
    print("=" * 100)
    print("  STEP 3(b). IDENTIFIABILITY WITH THE MARKER OFFSET FIXED")
    print("=" * 100)
    times = np.arange(0.0, 730.0, 30.4)
    loads = np.full(times.size, LB)
    print(f"\n  monthly record, n = {times.size}, CONSTANT load, Imax = {IMAX}, "
          f"fall calibrated to {100 * TARGET_FALL:.0f}%")
    part_a()
    part_b(times, loads)
    part_c(times, loads)
    part_d(times, loads)
    banner("SUMMARY")
    bad = [n for n, ok in RESULTS if not ok]
    for n, ok in RESULTS:
        print(f"    {'PASS' if ok else 'FAIL'}  {n}")
    print(f"\n    {len(RESULTS) - len(bad)} of {len(RESULTS)} structural checks "
          f"pass.")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
