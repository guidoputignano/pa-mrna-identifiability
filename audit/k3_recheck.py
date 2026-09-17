"""
Phase 5, step 0-2: the independent re-check.

Every check here reaches the same number by a route that shares no code with
the script it is checking. Where k1 used numpy's eigensolver this uses exact
rational arithmetic; where k2 used a bisection on a grid this uses the closed
form of the same quadratic; where k2 integrated with solve_ivp this integrates
by variation of constants on a dose-aligned grid.

A check that fails is reported as a failure and not repaired in place.
"""
import sys
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import k1_sources as K1
import k2_structures as K2
import g_design_common as G
import i1_estimator as I1
import design.run as DR

getcontext().prec = 60
RESULTS = []


def check(name, ok, detail):
    RESULTS.append((name, bool(ok), detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    for line in detail.splitlines():
        print(f"         {line}")
    print()


# --------------------------------------------------------------------------
def check_1_eigenvalues():
    """k1's PCC eigenvalues, by exact rational arithmetic."""
    kq = Fraction(474, 100000)
    kdeg = Fraction(75, 10000)
    # A = [[-(kq+kdeg), kq], [kq, -kq]];  trace = -(2kq+kdeg), det = kq*kdeg
    tr, det = -(2 * kq + kdeg), kq * kdeg
    disc = tr * tr - 4 * det
    assert disc > 0
    sq = Decimal(disc.numerator).sqrt() / Decimal(disc.denominator).sqrt()
    trd = Decimal(tr.numerator) / Decimal(tr.denominator)
    lam_slow = (trd + sq) / 2          # the one nearer zero
    lam_fast = (trd - sq) / 2
    ln2 = Decimal(2).ln()
    t_slow = float(ln2 / (-lam_slow) / 24)
    t_fast = float(ln2 / (-lam_fast) / 24)

    A = np.array([[-(K1.KQ + K1.KDEG), K1.KQ], [K1.KQ, -K1.KQ]])
    w = np.sort(np.linalg.eigvals(A).real)[::-1]
    n_slow, n_fast = np.log(2) / -w[0] / 24, np.log(2) / -w[1] / 24

    ok = abs(t_slow - n_slow) < 1e-9 and abs(t_fast - n_fast) < 1e-9
    check("1. PCC eigenvalues by exact rational arithmetic", ok,
          f"exact    slow {t_slow:.6f} d, fast {t_fast:.6f} d\n"
          f"numpy    slow {n_slow:.6f} d, fast {n_fast:.6f} d\n"
          f"ln2/kdeg alone {float(ln2 / Decimal(kdeg.numerator) * Decimal(kdeg.denominator) / 24):.6f} d\n"
          f"So the paper's '~7 d' is neither: it sits between "
          f"{t_fast:.2f} d and {t_slow:.2f} d.")


def check_2_ic50_closed_form():
    """k2's bisected IC50 roots, by the closed form of the same quadratic."""
    q = 1.0 - K2.IMAX
    u = K2.R_FLOOR
    E_lo, E_hi = K2.TH, K2.peak_enzyme()
    lines, ok = [], True
    for tag, (A_of, B_of, C_of) in {
        "PUB_A": (lambda r: 1.0 - r,
                  lambda r: (E_lo + q * E_hi) - r * (q * E_lo + E_hi),
                  lambda r: q * E_hi * E_lo * (1.0 - r)),
        "PUB_B": (lambda r: (1.0 - r) * (u + 1.0),
                  lambda r: (u * (1.0 - r) * (E_hi + E_lo)
                             + (E_lo + q * E_hi) - r * (E_hi + q * E_lo)),
                  lambda r: (1.0 - r) * E_hi * E_lo * (u + q)),
    }.items():
        cls = K2.PubA if tag == "PUB_A" else K2.PubB
        for tgt in (0.20, 0.30, 0.60):
            r = 1.0 - tgt
            A, B, C = A_of(r), B_of(r), C_of(r)
            disc = B * B - 4 * A * C
            if disc < 0:
                closed = np.array([])
            else:
                closed = np.sort(np.array([(-B + np.sqrt(disc)) / (2 * A),
                                           (-B - np.sqrt(disc)) / (2 * A)]))
                closed = closed[closed > 0]
            grid = K2._solve_shape(cls.shape, tgt)
            good = (closed.size == grid.size
                    and np.allclose(closed, grid, rtol=1e-6))
            ok = ok and good
            lines.append(f"{tag} fall {100*tgt:2.0f}%  closed "
                         f"{np.array2string(closed, precision=4)}  "
                         f"grid {np.array2string(grid, precision=4)}  "
                         f"{'ok' if good else 'MISMATCH'}")
    check("2. IC50 calibration roots: closed form vs bisection", ok,
          "\n".join(lines))


def check_3_ceilings():
    """k2's ceilings, by the stationarity condition instead of a grid scan."""
    lines, ok = [], True
    for cls in (K2.PubA, K2.PubB):
        ceil, c_at = K2._ceiling(cls.shape)
        # the optimum by a derivative sign change on a bracketing bisection
        h = 1e-6
        d = lambda c: (cls.shape(c * (1 + h)) - cls.shape(c * (1 - h))) / (2 * c * h)
        lo, hi = 1e-3, 1e9
        for _ in range(300):
            mid = np.sqrt(lo * hi)
            if d(mid) < 0:
                lo = mid
            else:
                hi = mid
        c2 = float(np.sqrt(lo * hi))
        ceil2 = float(1.0 - cls.shape(c2))
        good = abs(ceil - ceil2) < 1e-7
        ok = ok and good
        lines.append(f"{cls.name:38s} grid {100*ceil:.4f}% at {c_at:.5g}; "
                     f"stationary {100*ceil2:.4f}% at {c2:.5g}  "
                     f"{'ok' if good else 'MISMATCH'}")
    lines.append("")
    lines.append("The repository's committed falls, for comparison:")
    for tag, inc in (("first-cycle peak", DR.D_AMP),
                     ("steady-state peak",
                      DR.D_AMP / (1 - np.exp(-DR.LAMBDA * G.T)))):
        lines.append(f"  {tag:18s} {100*K1.fall(inc):.2f}%")
    check("3. Structure ceilings on the achievable fall", ok, "\n".join(lines))


def check_4_D_calibration():
    """k1's calibrated D, by direct substitution into the balance."""
    lines, ok = [], True
    for reading in ("peak", "trough"):
        for tgt in (0.20, 0.30):
            D, mult = K1.calibrate_D(tgt, reading)
            inc = D * mult
            # solve the balance from scratch, by bisection on m, not by the
            # quadratic root formula k1 uses
            f = lambda m: ((K1.TH + inc) * m / (K1.KM + m)
                           + K1.KE * m - K1.LB)
            lo, hi = 1e-12, 1e6
            for _ in range(300):
                mid = 0.5 * (lo + hi)
                if f(mid) < 0:
                    lo = mid
                else:
                    hi = mid
            m_tr = 0.5 * (lo + hi)
            got = 1.0 - m_tr / K1.steady_level(0.0)
            good = abs(got - tgt) < 1e-8
            ok = ok and good
            lines.append(f"{reading:7s} target {100*tgt:2.0f}%  D {D:10.4f}  "
                         f"Delta {inc:9.4f}  m {m_tr:.8f}  fall "
                         f"{100*got:.6f}%  {'ok' if good else 'MISMATCH'}")
    lines.append("")
    lines.append("and the untreated level, which must be the operating point:")
    lo, hi = 1e-12, 1e6
    f0 = lambda m: K1.TH * m / (K1.KM + m) + K1.KE * m - K1.LB
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if f0(mid) < 0:
            lo = mid
        else:
            hi = mid
    m0 = 0.5 * (lo + hi)
    good0 = abs(m0 - K2.M_BAR) < 1e-9
    ok = ok and good0
    lines.append(f"  m(Delta=0) = {m0:.10f}  {'ok' if good0 else 'MISMATCH'}")
    check("4. D calibrated to a target fall, by bisection on the balance",
          ok, "\n".join(lines))


def check_5_puba_integration():
    """k2's solve_ivp path, by variation of constants on a dose-aligned grid."""
    lines, ok = [], True
    for tgt in (0.30, 0.60):
        for c in K2._solve_shape(K2.PubA.shape, tgt):
            g = K2.PubA(c)
            times = G.SCHEDULES["monthly (30.4 d)"](12)
            # (a) at a CONSTANT load and no therapy the answer is analytic
            old = DR.D_AMP
            try:
                DR.D_AMP = 0.0
                m_null = g.marker(times, None)
                exact = K2.LB * K2._inhib(K2.TH, c) / g.kout
            finally:
                DR.D_AMP = old
            e_null = float(np.max(np.abs(m_null - exact)) / exact)
            # (b) with therapy, against an exponential integrator
            t0 = float(times.min()) - 10.0 * g.half_life()
            t1 = float(times.max()) + 1.0
            days = np.arange(np.floor(t0) - 1.0, np.ceil(t1) + 2.0)
            kin = np.full(days.size, K2.LB)
            # a sub-grid that never straddles a dose time or a day boundary
            doses = np.arange(0.0, t1 + G.T, G.T)
            breaks = np.unique(np.concatenate(
                [days, doses[(doses >= t0) & (doses <= t1)], [t0, t1], times]))
            breaks = breaks[(breaks >= t0) & (breaks <= t1)]
            m, tprev = kin[0] * K2._inhib(float(K2.enzyme(t0)), c) / g.kout, t0
            rec = {}
            for tb in breaks[1:]:
                h = tb - tprev
                if h <= 0:
                    continue
                # sub-steps resolve BOTH the exponential kernel (1/k_out)
                # and the therapeutic decay (1/lambda), whichever is faster
                fine = min(0.02, 0.05 / g.kout, 0.05 / DR.LAMBDA)
                nsub = max(16, int(np.ceil(h / fine)) // 2 * 2)
                u = np.linspace(0.0, h, nsub + 1)
                j = int(np.clip(np.searchsorted(days, tprev + 1e-12, "right") - 1,
                                0, days.size - 1))
                # forcing on the interval OPEN AT THE TOP. enzyme(t) is
                # right-continuous and jumps at a dose time, so evaluating the
                # last quadrature node exactly AT tb would use the post-dose
                # enzyme for an interval that ends before the dose. Nudging it
                # inside takes the left limit, which is the one the integral
                # over (tprev, tb) wants. Without this the reference was wrong
                # by 2e-4 relative and the check failed against a correct
                # integrator.
                tq = np.minimum(tprev + u, tb - 1e-9)
                fu = (kin[j] * K2._inhib(K2.enzyme(tq), c)
                      * np.exp(-g.kout * (h - u)))
                w = np.ones(nsub + 1)
                w[1:-1:2], w[2:-1:2] = 4.0, 2.0
                integral = (h / nsub) / 3.0 * float(np.dot(w, fu))
                m = m * np.exp(-g.kout * h) + integral
                tprev = tb
                rec[round(tb, 9)] = m
            ref = np.array([rec[round(t, 9)] for t in np.sort(times)])
            got = g.marker(np.sort(times), None)
            e_ther = float(np.max(np.abs(got - ref) / ref))
            good = e_null < 1e-8 and e_ther < 1e-6
            ok = ok and good
            lines.append(
                f"fall {100*tgt:2.0f}%  IC50 {c:11.5g}  t1/2 {g.half_life():7.3g} d  "
                f"no-therapy rel err {e_null:.2e}  "
                f"Simpson/VOC rel err {e_ther:.2e}  "
                f"{'ok' if good else 'MISMATCH'}")
    check("5. PUB_A: solve_ivp vs variation of constants + Simpson", ok,
          "\n".join(lines))


def check_6_control_row():
    """The paper's own model must be recovered EXACTLY by the paper's fit."""
    lines, ok = [], True
    for sname, fn in G.SCHEDULES.items():
        ts = fn(12)
        for cv in (0.0, 0.10):
            r = K2.one_fit(K2.Paper(), ts, cv, 3, noiseless=True)
            dth = abs(r["theta"] / K2.TH - 1.0)
            dkm = abs(r["Km"] / K2.KM - 1.0)
            # at a constant recorded load the fibre is exactly flat, so K_m is
            # NOT recoverable and only theta's identity along the fibre is
            singular = (cv == 0.0)
            good = (r["rms"] < 1e-6) and (singular or (dth < 1e-4 and dkm < 1e-4))
            ok = ok and good
            lines.append(f"{sname:26s} CV {100*cv:4.0f}%  rms/sig {r['rms']:.2e}  "
                         f"theta_hat/theta {r['theta']/K2.TH:12.6g}  "
                         f"Km_hat {r['Km']:10.5g}  "
                         f"{'(fibre: K_m free)' if singular else ''}  "
                         f"{'ok' if good else 'MISMATCH'}")
    check("6. Control row: the paper's model fitted to its own data", ok,
          "\n".join(lines))


def main():
    print("=" * 78)
    print("PHASE 5: THE INDEPENDENT RE-CHECK")
    print("=" * 78)
    print()
    check_1_eigenvalues()
    check_2_ic50_closed_form()
    check_3_ceilings()
    check_4_D_calibration()
    check_5_puba_integration()
    check_6_control_row()
    n = sum(1 for _, ok, _ in RESULTS if ok)
    print("=" * 78)
    print(f"  {n} of {len(RESULTS)} checks pass.")
    for name, okk, _ in RESULTS:
        if not okk:
            print(f"    FAILED: {name}")
    print("=" * 78)


if __name__ == "__main__":
    main()
