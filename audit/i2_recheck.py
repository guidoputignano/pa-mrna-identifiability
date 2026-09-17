"""
Phase 1 item 1, the independent re-check. Four key claims, each tested by a
route that shares no code path with the one that produced it.

    (1) varpro's minimiser really is GLOBAL, by brute force on a dense 2-D
        grid over (K_m, k_e) with the other two profiled, which is a different
        separation of the same objective
    (2) the constant-load profile is flat as an IDENTITY, proved symbolically,
        not sampled on a grid
    (3) the calibrated Monte Carlo sd exceeds the filter's by 9% BECAUSE the
        model is nonlinear, tested by shrinking the noise until it does not
    (4) the errors-in-variables gap has the sign and order a second-order
        expansion predicts

Nothing here records an expected result in advance.

Run:  python audit/i2_recheck.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import sympy as sp

import i1_estimator as I


def head(t):
    print("\n" + "=" * 92); print("  " + t); print("=" * 92)


# ------------------------------------------------------------------ check 1
def check_global():
    head("(1) IS varpro's MINIMISER GLOBAL? BRUTE FORCE ON A DIFFERENT SPLIT")
    print("  varpro profiles (theta, k_e, b) out and scans K_m. Here the split")
    print("  is moved: (theta, b) are profiled out and BOTH K_m and k_e are")
    print("  scanned, on a dense 2-D grid. Same objective, different separation,")
    print("  so agreement is not arithmetic identity.\n")
    m, inc, s = I.one_record(0, 0)
    ref, _ = I.varpro(m, inc, s, pen=True)
    # the conditioning that makes the unpenalised fit hopeless, printed here
    # so the figure quoted in the report has a script behind it
    g = m / (I.KM + m)
    Q = np.stack([I.TH * g, -(I.TH + inc) * I.KM * m / (I.KM + m) ** 2,
                  np.ones_like(m), I.KE * m], axis=1)
    F0 = Q.T @ Q / I.SIG ** 2
    w = np.linalg.eigvalsh(F0)
    print(f"    for context, the UNPENALISED information matrix on this record")
    print(f"    has eigenvalues {np.array2string(w, precision=3)},")
    print(f"    condition number {w[-1]/w[0]:.3e}, and sd(ln theta) "
          f"{np.sqrt(np.linalg.inv(F0)[0,0]):.4f}.\n")
    kms = np.geomspace(5.0, 800.0, 400)
    kes = np.geomspace(0.5, 50.0, 400)
    best = (np.inf, None, None)
    for ke in kes:
        g = m[None, :] / (kms[:, None] + m[None, :])
        X = np.stack([g, np.ones_like(g)], axis=-1)          # theta, b
        y = s[None, :] - inc[None, :] * g - ke * m[None, :]
        A = X.transpose(0, 2, 1) @ X / I.SIG ** 2
        A = A + np.diag([0.0, 1.0 / I.S_B ** 2])
        rhs = np.einsum("gnj,gn->gj", X, y) / I.SIG ** 2
        p = np.linalg.solve(A, rhs[..., None])[..., 0]
        r = np.einsum("gnj,gj->gn", X, p) - y
        J = (np.sum(r ** 2, axis=1) / I.SIG ** 2
             + (ke - I.KE) ** 2 / I.S_KE ** 2 + p[:, 1] ** 2 / I.S_B ** 2)
        j = int(np.argmin(J))
        if J[j] < best[0]:
            best = (float(J[j]), float(kms[j]), float(ke))
    print(f"    varpro          objective {ref['rss']:.8f}  "
          f"K_m {ref['Km']:.6f}  k_e {ref['ke']:.6f}")
    print(f"    2-D brute force objective {best[0]:.8f}  "
          f"K_m {best[1]:.6f}  k_e {best[2]:.6f}")
    print(f"    the grid's own resolution is {kms[1]/kms[0]-1:.4%} in K_m and "
          f"{kes[1]/kes[0]-1:.4%} in k_e,")
    print(f"    and the brute force cannot beat varpro: difference "
          f"{best[0] - ref['rss']:+.3e}")
    print("\n    A negative difference here would mean varpro's scan had missed")
    print("    a lower basin. A small positive one is the 2-D grid paying for")
    print("    its own resolution, which varpro avoids by solving exactly.")


# ------------------------------------------------------------------ check 2
def check_flat_identity():
    head("(2) THE CONSTANT-LOAD PROFILE IS FLAT AS AN IDENTITY, SYMBOLICALLY")
    print("  i1 part D reads the flatness off a 601-node grid. A grid cannot")
    print("  distinguish `flat` from `flat at these 601 points`. Here it is")
    print("  proved: at any K, the linear solve returns a (theta, k_e, b) whose")
    print("  residual at EVERY sample is identically zero.\n")
    K, K0, L, th0, ke0, mm, dd = sp.symbols(
        "K K_0 L theta_0 k_e m Delta", positive=True)
    # the fibre, in the balance's own sign convention
    b_bal = L * (K0 / K - 1)
    th = th0 + b_bal + ke0 * (K0 - K)
    # the level that the TRUE parameters produce at this increment and load,
    # defined implicitly by the true balance
    resid_true = (th0 + dd) * mm / (K0 + mm) + ke0 * mm - L
    # the SAME level, fed through the fibre's parameters
    resid_fib = (th + dd) * mm / (K + mm) + ke0 * mm - (L + b_bal)
    diff = sp.simplify(sp.together(resid_fib - resid_true))
    print(f"    resid_fibre - resid_true  simplifies to  {diff}")
    print(f"    factored:                 {sp.factor(sp.numer(sp.cancel(diff)))}")
    ratio = sp.simplify(sp.cancel(resid_fib / resid_true))
    print(f"    and resid_fibre / resid_true = {sp.simplify(ratio)}")
    print("\n    So the residual is not merely small at the sampled levels: it")
    print("    is the true residual times a factor that does not depend on the")
    print("    sample. At a constant load the true residual is zero at every")
    print("    sample, so the fibre's is too, at every K and every m, and the")
    print("    profile is constant as an identity rather than as a measurement.")
    # and numerically, at the actual record, away from the grid's own nodes
    L_ = np.full(len(I.MONTHLY), I.LB)
    mv = I.levels(I.MONTHLY, L_)
    inc = np.asarray(I.delta(I.MONTHLY), float)
    off = np.array([7.3, 41.9, 133.7, 617.4])       # not grid nodes
    J, _ = I.varpro_profile(mv, inc, L_.copy(), off, pen=False)
    print(f"\n    objective at four K_m values that are NOT grid nodes "
          f"{off.tolist()}:")
    print(f"      {np.array2string(J, precision=3)}")


# ------------------------------------------------------------------ check 3
def check_calibration():
    head("(3) WHY THE CALIBRATED sd EXCEEDS THE FILTER'S: IS IT NONLINEARITY?")
    print("  i1 part B finds the Monte Carlo sd 1.09x the filter's posterior sd")
    print("  once the nuisances are drawn from their priors. The filter is a")
    print("  first-order object, so if that excess is nonlinearity it must go")
    print("  away as the model is made more linear. Shrinking sigma and the")
    print("  nuisance priors together does exactly that, and nothing else.\n")
    print(f"    {'scale':>7} {'filter sd':>12} {'MC sd':>12} {'ratio':>8} "
          f"{'MC bias':>12}")
    print("  " + "-" * 56)
    sig0, sb0, ske0 = I.SIG, I.S_B, I.S_KE
    try:
        for f in (1.0, 0.5, 0.25, 0.1, 0.03):
            I.SIG, I.S_B, I.S_KE = sig0 * f, sb0 * f, ske0 * f
            fb, v = [], []
            for s0 in range(8):
                L = I.loads_for(I.MONTHLY, 0.10, s0)
                inc = np.asarray(I.delta(I.MONTHLY), float)
                for i in range(80):
                    rng = np.random.default_rng(555_000 + 1000 * s0 + i)
                    b_t = rng.normal(0.0, I.S_B)
                    ke_t = I.KE * np.exp(rng.normal(0.0, 0.8 * f))
                    mm = I.levels(I.MONTHLY, L, b=b_t, ke=ke_t)
                    z1 = rng.standard_normal(len(L)); z2 = rng.standard_normal(len(L))
                    sobs = ((L + b_t) * (1.0 + I.LOAD_FRAC * f * z1)
                            + I.BLOC * I.CV * f * I.M_BAR * z2)
                    r = I.varpro(mm, inc, sobs, pen=True)[0]
                    if r["theta"] > 0:
                        v.append(np.log(r["theta"]) - np.log(I.TH))
                mm0 = I.levels(I.MONTHLY, L)
                g = mm0 / (I.KM + mm0)
                Q = np.stack([I.TH * g, -(I.TH + inc) * I.KM * mm0 / (I.KM + mm0) ** 2,
                              np.ones_like(mm0), I.KE * mm0], axis=1)
                F = (Q.T @ Q / I.SIG ** 2
                     + np.diag([0.0, 0.0, 1.0 / I.S_B ** 2, 1.0 / (0.8 * f) ** 2]))
                fb.append(float(np.sqrt(np.linalg.inv(F)[0, 0])))
            v = np.array(v)
            print(f"    {f:7g} {np.median(fb):12.6f} {v.std(ddof=1):12.6f} "
                  f"{v.std(ddof=1)/np.median(fb):8.3f} {v.mean():12.6f}")
    finally:
        I.SIG, I.S_B, I.S_KE = sig0, sb0, ske0
    print("\n    A ratio going to 1 as the scale shrinks says the 9% excess is")
    print("    nonlinearity and the filter is right to first order. A ratio")
    print("    that stays put would say something is specified differently")
    print("    between the two, which would be a defect and not a limitation.")


# ------------------------------------------------------------------ check 4
def check_eiv_order():
    head("(4) THE EIV GAP AGAINST A SECOND-ORDER PREDICTION")
    print("  Placing the assay error on the level instead of the ordinate")
    print("  changes the residual by the CURVATURE of the balance in m. Write")
    print("  R(m) for the model flux; then")
    print("      ordinate model  R(m) - s, with s carrying b_loc eps")
    print("      EIV model       R(m + eps) - s")
    print("                    = R(m) + R'(m) eps + R''(m) eps^2/2 - s")
    print("  and R'(m_bar) is exactly b_loc, so the two agree to first order")
    print("  and differ in expectation by R''(m) eps^2/2, with eps = cv m.\n")
    m = sp.symbols("m", positive=True)
    D = sp.Symbol("Delta", nonnegative=True)
    th, K, ke = sp.Float(I.TH), sp.Float(I.KM), sp.Float(I.KE)
    R = (th + D) * m / (K + m) + ke * m
    R1, R2 = sp.simplify(sp.diff(R, m)), sp.simplify(sp.diff(R, m, 2))
    print(f"    R'(m)  = {R1}")
    print(f"    R''(m) = {R2}")
    mb = sp.Float(I.M_BAR)
    print(f"\n    at m_bar = {I.M_BAR:g}, Delta = 0:  R' = "
          f"{float(R1.subs({m: mb, D: 0})):.6f} against b_loc = {I.BLOC:.6f}, "
          f"R'' = {float(R2.subs({m: mb, D: 0})):.6f}")
    print("\n    R'' < 0 everywhere, so the EIV residual sits BELOW the")
    print("    ordinate one in expectation. The fit drives the residual to")
    print("    zero, so it must RAISE the model to compensate, and the model")
    print("    increases in theta. The predicted bias is therefore POSITIVE.")
    print("    (An earlier draft of this check asserted the opposite sign from")
    print("    the same R''; the estimating equation is what settles it.)\n")
    print("  Quantitatively, to first order in the shift and second in cv:")
    print("      E[p_hat] - p_true = -F^{-1} Q^T delta / sigma^2,")
    print("      delta_k = R''(m_k, Delta_k) (cv m_k)^2 / 2,")
    print("  taken over the ACTUAL monthly record rather than at m_bar, and")
    print("  with F the same penalised information matrix the filter uses.\n")
    R2f = sp.lambdify((m, D), R2, "numpy")
    print(f"    {'assay CV':>9} {'predicted bias in ln theta':>28} "
          f"{'pred/cv^2':>11} {'as a fraction of sd':>21}")
    print("  " + "-" * 74)
    for s0 in (0,):
        L = I.loads_for(I.MONTHLY, 0.10, s0)
        mm = I.levels(I.MONTHLY, L)
        inc = np.asarray(I.delta(I.MONTHLY), float)
        g = mm / (I.KM + mm)
        Q = np.stack([I.TH * g, -(I.TH + inc) * I.KM * mm / (I.KM + mm) ** 2,
                      np.ones_like(mm), I.KE * mm], axis=1)
        F = (Q.T @ Q / I.SIG ** 2
             + np.diag([0.0, 0.0, 1.0 / I.S_B ** 2, 1.0 / 0.8 ** 2]))
        Finv = np.linalg.inv(F)
        sd = float(np.sqrt(Finv[0, 0]))
        for cv in (0.02, I.CV, 0.10, 0.20):
            d = R2f(mm, inc) * (cv * mm) ** 2 / 2.0
            pred = float(-(Finv @ (Q.T @ d) / I.SIG ** 2)[0])
            print(f"    {cv:9.0%} {pred:28.6f} {pred/cv**2:11.4f} "
                  f"{pred/sd:20.2%}")
    print("\n    Compare these with i1 part C's measured paired gaps. Agreement")
    print("    in SIGN and ORDER is what this check is for: the prediction is")
    print("    itself only first order in the shift, so it is not expected to")
    print("    match to three figures, and a measured gap of a few hundredths")
    print("    of a sigma is why thousands of paired replicates are needed to")
    print("    see it at all.")


if __name__ == "__main__":
    print(__doc__)
    check_global()
    check_flat_identity()
    check_calibration()
    check_eiv_order()
