"""
2D.0 -- the full singular spectrum of the design matrix under constant load.

The question the handoff asks: is the deficiency exactly rank one, or higher?

Q is built on the four columns the estimator actually uses at the footing the
paper claims -- (ln theta, ln K_m, ln k_e, load bias) -- so the rank question is
asked of the matrix that carries the result, not of the three-column reduction
in which the identity was first found.
"""
import sys, numpy as np
sys.path.insert(0, 'src')
import design.run as D, model_A.run as LOC

P = D.P
L = D.load(0.0)


def Q_of(lv, dd, cols=("theta", "Km", "ke", "b")):
    out = []
    for name in cols:
        if name == "theta":
            out.append(P.theta * lv / (P.Km + lv))
        elif name == "Km":
            out.append(-(P.theta + dd) * P.Km * lv / (P.Km + lv) ** 2)
        elif name == "ke":
            out.append(P.ke * lv)
        elif name == "b":
            out.append(np.ones_like(lv))
    return np.stack(out, axis=1)


def report(label, lv, dd, cols):
    Q = Q_of(lv, dd, cols)
    s = np.linalg.svd(Q, compute_uv=False)
    rank = int(np.linalg.matrix_rank(Q, tol=1e-10 * s[0]))
    print(f"  {label}")
    print(f"    columns {cols}")
    print(f"    singular values  {np.array2string(s, precision=6, suppress_small=False)}")
    print(f"    s_min/s_max = {s[-1]/s[0]:.3e}   numerical rank {rank} of {len(cols)}"
          f"   -> deficiency {len(cols)-rank}")
    return s, rank


def main():
    print(__doc__)
    t = np.arange(0.0, 730.0, 2.0)
    lv = D.quasi_steady_level(t)
    dd = np.asarray(D.delta(t), dtype=float)

    print("=" * 78)
    print("CONSTANT recorded load (the design chapter as coded)")
    print("=" * 78)
    report("three columns, as the C2 identity was stated", lv, dd, ("theta", "Km", "b"))
    print()
    s4, r4 = report("FOUR columns -- the footing the paper claims", lv, dd,
                    ("theta", "Km", "ke", "b"))
    print()
    report("three columns, k_e in place of the load bias", lv, dd, ("theta", "Km", "ke"))

    print()
    print("=" * 78)
    print("Is the second small direction real, or is it conditioning?")
    print("=" * 78)
    Q = Q_of(lv, dd)
    U, s, Vt = np.linalg.svd(Q, full_matrices=False)
    for i in range(len(s)):
        print(f"    s[{i}] = {s[i]:14.6e}   right vector (th,Km,ke,b) = "
              f"{np.array2string(Vt[i], precision=5, suppress_small=True)}")
    print()
    # the exact identity, expressed in the four-column space
    c = np.array([(1.0 + P.Km * P.ke / L) / P.theta, -1.0 / L, 0.0, -1.0])
    print(f"    the C2 identity vector, normalised: "
          f"{np.array2string(c/np.linalg.norm(c), precision=5, suppress_small=True)}")
    print(f"    ||Q c|| / ||Q|| = {np.linalg.norm(Q @ c)/np.linalg.norm(Q):.3e}")
    print(f"    |cos| against the smallest right singular vector: "
          f"{abs(np.dot(c/np.linalg.norm(c), Vt[-1])):.10f}")
    print()
    # a k_e-only identity?  k_e*m is not constant, so test explicitly
    print("    Is there a SECOND exact identity involving k_e? Least squares of each")
    print("    column on the others, residual as a fraction of the column norm:")
    for j, nm in enumerate(("ln theta", "ln K_m", "ln k_e", "load bias")):
        A = np.delete(Q, j, axis=1); y = Q[:, j]
        r = y - A @ np.linalg.lstsq(A, y, rcond=None)[0]
        print(f"      {nm:>10}: {np.linalg.norm(r)/np.linalg.norm(y):.3e}")

    print()
    print("=" * 78)
    print("VARYING recorded load -- does the deficiency close?")
    print("=" * 78)
    rng = np.random.default_rng(0)
    for cv in (0.0, 0.005, 0.01, 0.02, 0.05, 0.10):
        Ls = L * (1.0 + cv * rng.standard_normal(t.size)) if cv else np.full(t.size, L)
        B = P.ke * P.Km + P.theta + dd - Ls
        lv2 = (-B + np.sqrt(B ** 2 + 4.0 * P.ke * Ls * P.Km)) / (2.0 * P.ke)
        Q2 = Q_of(lv2, dd)
        s2 = np.linalg.svd(Q2, compute_uv=False)
        print(f"    load CV {cv:5.1%}:  s = {np.array2string(s2, precision=4)}"
              f"   s_min/s_max = {s2[-1]/s2[0]:.3e}"
              f"   rank {int(np.linalg.matrix_rank(Q2, tol=1e-10*s2[0]))}")


if __name__ == "__main__":
    main()
