"""
P3 -- the PI's multiplied-out linearisation, his constrained solve, and what
each of them is actually estimating.

His note multiplies the balance by (K_M + m) and writes

    K_M psi + psi m - L_0 + theta m + a K_M m + a m^2 = 0

with regressors [psi, -1, m, m, m^2], coefficients (K_M, L_0, theta, a K_M, a),
and the constraint  first * last = fourth.  psi is an UNKNOWN; Delta is the
therapy, known from the dosing record.

This script derives the paper's own multiplied-out balance, maps his symbols
onto it, tests two readings of psi, and then implements his solve exactly as
written so it can be compared with the true constrained least-squares
minimiser.

Nothing below records an expected result in advance.

Run:  python audit/f3_linearisation.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import sympy as sp
from scipy.optimize import least_squares

from common.models import g
from design.run import P, delta, load, LOAD_FRAC, ASSAY_CV
import model_A.run as LOC

TH0, KM0, KE0 = P.theta, P.Km, P.ke
LB = load(0.0)
SIG = LOC.sigma_total(LOAD_FRAC, ASSAY_CV)

MONTHLY = np.arange(0.0, 730.0, 30.4)


def banner(t):
    print("\n" + "=" * 78)
    print("  " + t)
    print("=" * 78)


def record(times, cv=0.0, seed=3, b=0.0):
    """Levels, increments and the per-sample RECORDED load."""
    inc = np.asarray(delta(times), dtype=float)
    if cv:
        rng = np.random.default_rng(seed)
        s2 = np.log1p(cv ** 2)
        L = LB * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), len(times)))
    else:
        L = np.full(len(times), LB)
    tot = L + b
    B = KE0 * KM0 + TH0 + inc - tot
    m = (-B + np.sqrt(B ** 2 + 4.0 * KE0 * tot * KM0)) / (2.0 * KE0)
    return m, inc, L


# --------------------------------------------------------------------------
def s1_the_papers_multiplied_out_balance():
    banner("S1  the paper's balance, multiplied out")

    th, km, ke, b, L, m, dl = sp.symbols(
        "theta K_m k_e b L m Delta", positive=True)
    lhs = sp.expand(sp.cancel(
        ((th + dl) * m / (km + m) + ke * m - (L + b)) * (km + m)))
    print("  Starting from  (theta + Delta) m/(K_m+m) + k_e m = L + b  and")
    print("  multiplying by (K_m + m):\n")
    print(f"    0 = {sp.collect(lhs, m)}")
    print()
    print("  Put the terms that are KNOWN sample by sample (Delta_k and m_k are")
    print("  both known: the increment from the dosing record, the marker from")
    print("  the assay) on one side:\n")
    resp = sp.expand(dl * m)
    rest = sp.expand(-(lhs - resp))
    print(f"    Delta m = {sp.collect(rest, m)}")
    print()
    print("  Grouping against the regressors (L, 1, m, m^2), with L the")
    print("  RECORDED load and everything else unknown:")
    print(f"    Delta m - L m = K_m L + K_m b + m (b - theta - k_e K_m) - k_e m^2")
    chk = sp.simplify(sp.expand(rest - L * m)
                      - (km * L + km * b + m * (b - th - ke * km) - ke * m ** 2))
    print(f"    check, difference from that form: {chk}")
    assert chk == 0
    print()
    print("  So the model is LINEAR in")
    print("      beta = (K_m, K_m b, b - theta - k_e K_m, -k_e)")
    print("  against the regressors (L_k, 1, m_k, m_k^2), with the response")
    print("  (Delta_k - L_k) m_k.  That is the second reader's beta exactly.")
    return dict(th=th, km=km, ke=ke, b=b, L=L, m=m, dl=dl)


# --------------------------------------------------------------------------
def s2_map_his_symbols(S):
    banner("S2  mapping his symbols onto the paper's")

    print("  His form:   K_M psi + psi m - L_0 + theta m + a K_M m + a m^2 = 0")
    print("  Paper's:    K_m (L+b) + (L+b) m - theta m - k_e K_m m - k_e m^2")
    print("                                                 - Delta m = 0")
    print()
    print("  Term by term:")
    rows = [
        ("K_M psi", "K_m (L+b)", "matches if psi is the TOTAL load L + b"),
        ("psi m", "(L+b) m", "matches, same reading of psi"),
        ("- L_0", "(nothing)", "the paper's balance has NO free constant "
                               "column;\n                     "
                               "the constant is K_m (L+b), already present"),
        ("+ theta m", "- theta m", "opposite sign, so his theta is -theta"),
        ("+ a K_M m", "- k_e K_m m", "matches with a = -k_e"),
        ("+ a m^2", "- k_e m^2", "matches with a = -k_e"),
        ("(absent)", "- Delta m", "HIS FORM HAS NO DOSE TERM"),
    ]
    print(f"    {'his term':>12}  {'paper term':>14}   comment")
    for h, p, c in rows:
        print(f"    {h:>12}  {p:>14}   {c}")

    print()
    print("  The single most consequential line is the last. Delta_k m_k is the")
    print("  entire therapy channel: it is the term that makes theta")
    print("  identifiable at all (PHASE2D_FINDINGS.md section 1, 42x the")
    print("  excitation the marker band supplies). Dropped, his design is the")
    print("  no-dose-train column of rem:rankvsprec.")
    print()
    print("  Delta_k m_k belongs on the KNOWN side: Delta_k is read off the")
    print("  dosing record and m_k off the assay, so it is a response, not a")
    print("  regressor with an unknown coefficient.")


# --------------------------------------------------------------------------
def s3_two_readings_of_psi(S):
    banner("S3  two readings of psi, and what each costs")

    th, km, ke, b, L, m, dl = (S[k] for k in ("th", "km", "ke", "b", "L", "m", "dl"))

    print("  READING 1: psi_k = L_k + b, the TRUE total load. L_k recorded and")
    print("  known, b an unknown constant bias.\n")
    print("    His products: K_M psi = K_M L_k + K_M b, and psi m = L_k m + b m.")
    print("    K_M L_k is (unknown x known): linear.")
    print("    K_M b     is (unknown x unknown): NOT linear in his parameters.")
    print("    b m       is (unknown x known): linear, but it lands on the SAME")
    print("              column as his theta m and his a K_M m.")
    print()
    print("    So his form is not linear in (K_M, L_0, theta, a). It becomes")
    print("    linear only after the reparametrisation beta above, in which")
    print("    K_M b is carried as one coefficient and never split.")
    print()
    print("    Where b appears: only inside psi. It never gets a column of its")
    print("    own. It reaches the fit as K_m b (the constant column) and as")
    print("    part of (b - theta - k_e K_m) (the m column).")
    print()
    print("    What his m coefficient really is: not theta, and not -theta")
    print("    either. It is (b - theta - k_e K_m): the load bias and the")
    print("    alternative channel are both inside it. Reading that coefficient")
    print("    as a capacity is a three-way confusion.")
    print()
    print("  READING 2: psi a single unknown CONSTANT load.\n")
    print("    Then the psi column is a constant column, and so is his -1")
    print("    column. K_M psi and L_0 multiply the same regressor, so only")
    print("    the combination (K_M psi - L_0) is identified; neither alone is.")
    print("    Adding L_0 to the model adds nothing the data can see.")
    print()
    print("  In BOTH readings, columns 3 and 4 of his design are literally the")
    print("  same regressor, m_k. theta and a K_M are not separately identified")
    print("  by the columns; only their sum is. His constraint")
    print("  first * last = fourth is what splits them, and it is correct as an")
    print("  algebraic identity: c1 c5 = K_M * a = a K_M = c4.")
    print()

    # rank of his design, both readings, constant and varying load
    print("  Rank of his 5-column design matrix [psi, -1, m, m, m^2], by")
    print("  reading and by whether the recorded load varies:\n")
    print(f"    {'reading':>34} {'load CV':>9} {'rank':>6} {'of':>4} "
          f"{'s_min/s_max':>13}")
    print("  " + "-" * 74)
    for cv in (0.0, 0.10):
        m_k, inc, L_k = record(MONTHLY, cv=cv)
        for lab, psi in (("psi = L_k (recorded, varying)", L_k),
                         ("psi = one unknown constant", np.ones_like(m_k))):
            X = np.column_stack([psi, -np.ones_like(m_k), m_k, m_k, m_k ** 2])
            Xn = X / np.linalg.norm(X, axis=0)
            sv = np.linalg.svd(Xn, compute_uv=False)
            print(f"    {lab:>34} {cv*100:8.1f}% "
                  f"{int(np.sum(sv > 1e-12 * sv[0])):6d} {X.shape[1]:4d} "
                  f"{sv[-1]/sv[0]:13.4e}")

    print()
    print("  At a constant recorded load, both readings give rank 3 of 5. The")
    print("  constraint removes one of the two deficiencies. The one left is")
    print("  the SAME object as prop:rank and as the constant-load fibre of")
    print("  f5_fibre.py: K_m and (L+b) enter only through their product.")


# --------------------------------------------------------------------------
_LAST_QUAD = (np.nan, np.nan, np.nan)


def which_partitions_can_be_formed(X):
    """
    His scheme profiles v out of the normal equations at a fixed scalar w.
    Columns 3 and 4 of his design are the same regressor, so the v-block is
    singular unless w is one of their two coefficients. This reports which of
    the five choices of w leaves an invertible v-block.
    """
    out = []
    for j in range(5):
        keep = [i for i in range(5) if i != j]
        Xv = X[:, keep]
        Xn = Xv / np.linalg.norm(Xv, axis=0)
        sv = np.linalg.svd(Xn, compute_uv=False)
        out.append((j, sv[-1] / sv[0]))
    return out


def his_solve(X, y, wcol=3):
    """
    His solve, implemented as written.

    v = the four coefficients other than w are profiled out by the conditional
    normal equations at fixed w; the constraint c1 c5 = c4 then becomes a
    scalar equation in w, which is a quadratic; the positive root is kept.

    `wcol` is which coefficient is held as w. Only wcol in {2, 3} (the two
    coefficients on the duplicated m column) leaves the v-block invertible,
    which is itself part of the finding.
    """
    keep = [i for i in range(5) if i != wcol]
    Xv, xw = X[:, keep], X[:, wcol]
    G = Xv.T @ Xv
    p = np.linalg.solve(G, Xv.T @ y)          # v(w) = p + q w
    q = -np.linalg.solve(G, Xv.T @ xw)

    def coef(i, w):
        """Coefficient i of the full five-vector, as a function of w."""
        if i == wcol:
            return w
        j = keep.index(i)
        return p[j] + q[j] * w

    # constraint c1 c5 - c4 = 0, a quadratic in w
    a1, b1 = (0.0, 1.0) if 0 == wcol else (p[keep.index(0)], q[keep.index(0)])
    a5, b5 = (0.0, 1.0) if 4 == wcol else (p[keep.index(4)], q[keep.index(4)])
    a4, b4 = (0.0, 1.0) if 3 == wcol else (p[keep.index(3)], q[keep.index(3)])
    A = b1 * b5
    B = a1 * b5 + a5 * b1 - b4
    C = a1 * a5 - a4
    global _LAST_QUAD
    _LAST_QUAD = (A, B, C)
    if abs(A) * max(abs(B), 1e-300) ** -1 < 1e-12:
        # A vanishes identically: the two m columns are the same regressor, so
        # v(w) moves only in the coefficient of the column w was taken from,
        # and q has a single nonzero entry. The quadratic is a linear equation
        # and the second "root" is an artefact of dividing by A.
        roots = np.array([-C / B]) if abs(B) > 0 else np.array([])
    else:
        disc = B * B - 4 * A * C
        roots = (np.array([]) if disc < 0 else
                 np.array([(-B + np.sqrt(disc)) / (2 * A),
                           (-B - np.sqrt(disc)) / (2 * A)]))
    pos = roots[roots > 0]
    cand = pos if pos.size else roots
    if not cand.size:
        return np.full(5, np.nan), roots
    # "keep the positive root" does not disambiguate when both are positive,
    # so the residual decides, and the ambiguity is itself reported.
    ss = [float(np.sum((np.array([coef(i, float(w)) for i in range(5)]) @ X.T
                        - y) ** 2)) for w in cand]
    w = float(cand[int(np.argmin(ss))])
    return np.array([coef(i, w) for i in range(5)]), roots


def true_constrained_ls(X, y, x0):
    """The actual minimiser of ||Xc - y||^2 subject to c1 c5 = c4."""
    def res(z):
        c1, c2, c3, c5 = z
        c = np.array([c1, c2, c3, c1 * c5, c5])
        return X @ c - y
    s = least_squares(res, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15)
    c1, c2, c3, c5 = s.x
    return np.array([c1, c2, c3, c1 * c5, c5]), float(np.sum(s.fun ** 2))


# --------------------------------------------------------------------------
def design_his(m_k, L_k):
    """
    His five columns, with the sign of the psi terms as the algebra requires.

    As literally written his equation is

        K_M psi + psi m - L_0 + theta m + a K_M m + a m^2 = 0,

    in which the psi terms and the clearance terms carry the SAME sign. They
    cannot: multiplying theta g(m) + a m = psi by (K_M + m) gives

        theta m + a K_M m + a m^2 - K_M psi - psi m = 0,

    so the psi terms enter negatively. Run as written, every coefficient comes
    back negative and his "positive root" rule selects nothing. The sign is
    corrected here, because the alternative is to report a transcription slip
    rather than test the method.
    """
    return np.column_stack([-L_k, -np.ones_like(m_k), m_k, m_k, m_k ** 2])


def response_his(m_k, inc, L_k, with_dose):
    """psi m, and -Delta m when the dose term is restored."""
    return L_k * m_k - (inc * m_k if with_dose else 0.0)


def s4_run_his_solve():
    banner("S4  his solve, run as written (sign of the psi terms corrected)")

    print("  The linear route only exists if psi is treated as KNOWN. With psi")
    print("  a free unknown, K_M psi is a product of unknowns and the problem")
    print("  is nonlinear; there is no set of normal equations to condition.")
    print("  So below psi := the recorded load L_k, which is known, and the")
    print("  unknown bias is carried where the algebra puts it.")
    print(design_his.__doc__)
    print(f"  At the truth the five coefficients are")
    print(f"    (K_M, L_0, theta, a K_M, a) = ({KM0:.1f}, 0.0, {TH0:.1f}, "
          f"{KE0*KM0:.1f}, {KE0:.1f}), and c1 c5 - c4 = "
          f"{KM0*KE0 - KE0*KM0:.1f}\n")

    m_k, inc, L_k = record(MONTHLY, cv=0.10)
    X0 = design_his(m_k, L_k)
    print("  First, which coefficient can be held as w at all. Holding w = c_j")
    print("  leaves the other four as the v-block; s_min/s_max of that block:")
    names = ["c1 = K_M", "c2 = L_0", "c3 = theta", "c4 = a K_M", "c5 = a"]
    for j, r in which_partitions_can_be_formed(X0):
        print(f"    w = {names[j]:>11}: v-block s_min/s_max = {r:10.3e}"
              f"   {'unformable' if r < 1e-12 else 'formable'}")
    print("  Only the two coefficients on the duplicated m column leave an")
    print("  invertible v-block, so w := c4 = a K_M below.\n")

    cases = [("load CV 10%, b_true = 0, noiseless", 0.10, 0.0, 0.0),
             ("load CV 10%, b_true = 32.3 (10% of L)", 0.10, 32.27, 0.0),
             ("load CV 10%, b_true = 0, marker CV 5%", 0.10, 0.0, ASSAY_CV),
             ("load CONSTANT, b_true = 0, noiseless", 0.0, 0.0, 0.0)]

    for lab, cv, b_true, mcv in cases:
        print(f"  ---- {lab} " + "-" * max(2, 56 - len(lab)))
        m_k, inc, L_k = record(MONTHLY, cv=cv, b=b_true)
        if mcv:
            m_k = m_k * (1.0 + mcv * np.random.default_rng(17)
                         .standard_normal(m_k.size))
        X = design_his(m_k, L_k)

        for with_dose, dlab in ((True, "with Delta_k m_k restored"),
                                (False, "his form, NO dose term")):
            y = response_his(m_k, inc, L_k, with_dose)
            try:
                c_his, roots = his_solve(X, y, wcol=3)
            except np.linalg.LinAlgError:
                print(f"    {dlab:>36}: normal equations SINGULAR, no solve")
                continue
            if not np.all(np.isfinite(c_his)):
                print(f"    {dlab:>36}: no finite root")
                continue
            c_tru, ss = true_constrained_ls(
                X, y, np.array([c_his[0], c_his[1], c_his[2], c_his[4]]))
            ss_his = float(np.sum((X @ c_his - y) ** 2))
            print(f"    {dlab}")
            print(f"      his scalar equation in w: A w^2 + B w + C with "
                  f"(A, B, C) = ({_LAST_QUAD[0]:.3e}, {_LAST_QUAD[1]:.6f}, "
                  f"{_LAST_QUAD[2]:.4f})")
            print(f"      roots: {np.array2string(roots, precision=4)}"
                  f"   positive roots: {int(np.sum(roots > 0))}")
            print(f"      {'':>20} {'truth':>14} {'his solve':>16} "
                  f"{'true constrained LS':>21}")
            truth = [KM0, 0.0, TH0, KE0 * KM0, KE0]
            for i, nm in enumerate(("c1 = K_M", "c2 = L_0", "c3 = theta",
                                    "c4 = a K_M", "c5 = a")):
                print(f"      {nm:>20} {truth[i]:14.4f} {c_his[i]:16.6f} "
                      f"{c_tru[i]:21.6f}")
            print(f"      {'residual SS':>20} {'':>14} {ss_his:16.6e} "
                  f"{ss:21.6e}")
            print(f"      {'constraint c1c5-c4':>20} {'':>14} "
                  f"{c_his[0]*c_his[4]-c_his[3]:16.3e} "
                  f"{c_tru[0]*c_tru[4]-c_tru[3]:21.3e}")
        print()


# --------------------------------------------------------------------------
def s5_is_it_singular_at_constant_load():
    banner("S5  is his solve singular at a constant load?")

    for cv in (0.0, 0.001, 0.01, 0.10):
        m_k, inc, L_k = record(MONTHLY, cv=cv)
        X = design_his(m_k, L_k)
        G = X[:, [0, 1, 2, 4]].T @ X[:, [0, 1, 2, 4]]
        w = np.linalg.eigvalsh(G)
        print(f"    load CV {cv*100:6.2f}%: cond = "
              f"{w[-1]/max(w[0], 1e-300):12.4e}   "
              f"{'SINGULAR' if w[0] <= 1e-12 * w[-1] else 'invertible'}")
    print()
    print("  (the block shown is the one his scheme CAN form, w = c4, so the")
    print("   duplicated m column is already out of it.)")
    print()
    print("  So his design is singular at a constant recorded load for two")
    print("  independent reasons at once: the psi column equals the constant")
    print("  column, and the two m columns are identical. His constraint")
    print("  repairs the second. The first is prop:rank, and no constraint on")
    print("  the coefficients repairs it -- only a varying load does.")


# --------------------------------------------------------------------------
def main():
    print(__doc__)
    S = s1_the_papers_multiplied_out_balance()
    s2_map_his_symbols(S)
    s3_two_readings_of_psi(S)
    s4_run_his_solve()
    s5_is_it_singular_at_constant_load()


if __name__ == "__main__":
    main()
