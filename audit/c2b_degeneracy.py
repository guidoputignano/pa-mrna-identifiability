"""
C2 (second half) -- what the corrected column does to the RANK.

Reinstating Delta does not double count.  It does something the audit did not
anticipate: it makes the three-column design matrix EXACTLY rank deficient.

THE IDENTITY.  Write u = m/(K_m+m), so K_m/(K_m+m) = 1-u.  On any record whose
levels satisfy the balance at a RECORDED LOAD L,

    (theta + Delta_i) u_i = L_i - k_e m_i ,

the corrected columns are

    q_theta = theta u ,    q_Km = -(theta+Delta)u(1-u) = -(L - k_e m)(1-u) ,
    q_b     = 1 .

Ask for constants c1, c2 with c1 q_theta + c2 q_Km = q_b.  Clearing (K_m+m),

    m (c1 theta + c2 K_m k_e - 1) = K_m (1 + c2 L) .

If L is CONSTANT both sides vanish identically at

    c2 = -1/L ,      c1 = (1 + K_m k_e / L) / theta ,

so q_b lies exactly in the span of the other two and the information matrix is
singular.  If L VARIES the first equation forces c2 = -1/L_i for every i, which
no constant satisfies, and the degeneracy is broken.

So the published column was not hiding a duplicate.  It was hiding a null
direction, and the null direction is the load bias.  Everything A5 reports rests
on the 5% load-record prior holding that direction down.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
import sympy as sp
import design.run as D, model_A.run as LOC

P = D.P
Lb = D.load(0.0)
SIG = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)


def symbolic():
    th, Km, ke, L, m, dl = sp.symbols('theta K_m k_e L m Delta', positive=True)
    # Delta eliminated through the balance, which is what generates the levels
    dl_sol = sp.solve(sp.Eq((th + dl) * m / (Km + m) + ke * m, L), dl)[0]
    q_th = th * m / (Km + m)
    q_Km = -(th + dl_sol) * Km * m / (Km + m) ** 2
    c1 = (1 + Km * ke / L) / th
    c2 = -1 / L
    resid = sp.simplify(c1 * q_th + c2 * q_Km - 1)
    print("  c1*q_theta + c2*q_Km - 1  simplifies to:", resid)
    assert resid == 0, "identity failed"
    print("  -> the load-bias column is EXACTLY in the span, for every m, at constant L")
    # and with L varying sample to sample it is not
    L1, L2 = sp.symbols('L_1 L_2', positive=True)
    print(f"  with L varying the same construction needs c2 = -1/L_i, "
          f"so it fails unless L_1 = L_2")


def fisher(lv, dd, Ls, use_delta, b_prior_frac=None, ke_free=False, ke_prior=0.8):
    pri = [0.8 ** 2, 0.8 ** 2]
    pri.append((b_prior_frac * Lb) ** 2 if b_prior_frac else np.inf)
    if ke_free:
        pri.append(ke_prior ** 2)
    Finv = np.diag([0.0 if np.isinf(p) else 1.0 / p for p in pri])
    for m, dl, L in zip(lv, dd, Ls):
        eff = P.theta + (dl if use_delta else 0.0)
        q = [P.theta * m / (P.Km + m), -eff * P.Km * m / (P.Km + m) ** 2, 1.0]
        if ke_free:
            q.append(P.ke * m)
        q = np.array(q)
        Finv = Finv + np.outer(q, q) / SIG ** 2
    w = np.linalg.eigvalsh(Finv)
    if w[0] <= 1e-9 * w[-1]:
        return np.inf, np.nan, w
    C = np.linalg.inv(Finv)
    return float(np.sqrt(C[0, 0])), float(C[0, 1] / np.sqrt(C[0, 0] * C[1, 1])), w


def levels_at(t, Ls):
    """Quasi-steady level solving the balance at a per-sample recorded load."""
    d = np.asarray(D.delta(t), dtype=float)
    B = P.ke * P.Km + P.theta + d - Ls
    return (-B + np.sqrt(B ** 2 + 4.0 * P.ke * Ls * P.Km)) / (2.0 * P.ke)


def main():
    print(__doc__)
    print("=" * 84)
    print("SYMBOLIC")
    print("=" * 84)
    symbolic()

    t = np.arange(0.0, 730.0, 2.0)
    dd = np.asarray(D.delta(t), dtype=float)
    const = np.full(t.size, Lb)
    lv = levels_at(t, const)

    print()
    print("=" * 84)
    print("NUMERIC -- constant recorded load, as the design chapter has it")
    print("=" * 84)
    for ud in (False, True):
        sd, r, w = fisher(lv, dd, const, ud, b_prior_frac=None)
        print(f"  Delta {'in ':>7} eigs {np.array2string(w, precision=3)}"
              if ud else f"  Delta {'dropped':>7} eigs {np.array2string(w, precision=3)}")
        print(f"           prior-free sd(ln theta) = "
              f"{'SINGULAR' if np.isinf(sd) else f'{sd:.4f}'}")

    print()
    print("  with the module's 5% load-record prior, which is what makes it finite:")
    for ud in (False, True):
        sd, r, _ = fisher(lv, dd, const, ud, b_prior_frac=D.LOAD_FRAC)
        print(f"    Delta {'in' if ud else 'dropped':>7}:  sd(ln theta) = {sd:.4f}   corr = {r:+.4f}")

    print()
    print("=" * 84)
    print("HOW HARD THE LOAD PRIOR HAS TO WORK  (this IS 2B.1, arrived at the other way)")
    print("=" * 84)
    print(f"{'load record':>12} | {'Delta dropped':>14} {'Delta in':>10} | "
          f"{'Delta dropped':>14} {'Delta in':>10}")
    print(f"{'':>12} | {'k_e known':>25} | {'k_e free (0.8 prior)':>25}")
    for f in (0.01, 0.02, 0.05, 0.10, 0.20, 0.50):
        row = []
        for kef in (False, True):
            for ud in (False, True):
                sd, _, _ = fisher(lv, dd, const, ud, b_prior_frac=f, ke_free=kef)
                row.append(sd)
        print(f"{f:11.0%} | {row[0]:14.4f} {row[1]:10.4f} | {row[2]:14.4f} {row[3]:10.4f}")

    print()
    print("=" * 84)
    print("DOES A VARYING LOAD BREAK THE DEGENERACY?  (day-to-day variation in intake)")
    print("=" * 84)
    print(f"{'load CV':>9} | {'prior-free, Delta in':>21} {'prior-free, Delta dropped':>26} | "
          f"{'5% prior, Delta in':>19}")
    rng = np.random.default_rng(0)
    for cv in (0.0, 0.02, 0.05, 0.10, 0.20, 0.40):
        Ls = Lb * (1.0 + cv * rng.standard_normal(t.size)) if cv else const
        lv2 = levels_at(t, Ls)
        a, _, _ = fisher(lv2, dd, Ls, True, b_prior_frac=None)
        b, _, _ = fisher(lv2, dd, Ls, False, b_prior_frac=None)
        c, _, _ = fisher(lv2, dd, Ls, True, b_prior_frac=D.LOAD_FRAC)
        f = lambda v: 'SINGULAR' if np.isinf(v) else f'{v:.4f}'
        print(f"{cv:8.0%} | {f(a):>21} {f(b):>26} | {f(c):>19}")


if __name__ == "__main__":
    main()
