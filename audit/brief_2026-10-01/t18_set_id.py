"""
Item 16. Set identification elevated: what IS recoverable when theta is not.

Three parts:
  (A) the one-sided bound on theta from the balance, non-negativity and K_m;
  (B) the reparameterisations the record DOES pin (K_m, the clearance slope
      b = theta g' + k_e, and the k_e-dominated combination), and the ones it
      does NOT (theta, a = theta g^2, the untreated flux, the fraction oxidised
      -- all proportional to theta);
  (C) the complete-deficiency check the brief asked for twice.
"""
import numpy as np

import core as C
import t4_set as T4

KE_LO_95 = C.KE * np.exp(-1.96 * C.KE_PRIOR)


def recoverable(r, h=1e-6):
    """sds of derived scalars from the monthly-record cov, median later."""
    P = C.cov(C.fim(r, noise="perdraw"))
    th, K, ke, Lr = r["theta"], r["Km"], r["ke"], r["L"]
    base = np.array([np.log(th), np.log(K), 0.0, np.log(ke)])

    def quants(x):
        t, k, b_, e = np.exp(x[0]), np.exp(x[1]), x[2], np.exp(x[3])
        m = C.m_ss(t, k, e, Lr - b_)
        g, gp = m / (k + m), k / (k + m) ** 2
        return dict(flux=t * g, share=t * g / (Lr - b_),
                    a=t * g * g, b=t * gp + e)   # Prop-1 affine coefficients
    keys = ["flux", "share", "a", "b"]
    q0 = quants(base)
    G = {q: np.zeros(4) for q in keys}
    for j in range(4):
        hj = h * (Lr if j == 2 else 1.0)
        xp, xm = base.copy(), base.copy()
        xp[j] += hj
        xm[j] -= hj
        qp, qm = quants(xp), quants(xm)
        for q in keys:
            G[q][j] = (np.log(qp[q]) - np.log(qm[q])) / (2 * hj)
    sd = lambda g: float(np.sqrt(g @ P @ g))
    sub = P[np.ix_([0, 3], [0, 3])]
    w, V = np.linalg.eigh(sub)
    return (np.sqrt(P[0, 0]), np.sqrt(P[1, 1]), sd(G["flux"]), sd(G["share"]),
            sd(G["a"]), sd(G["b"]), np.sqrt(w[0]), V[0, 0], V[1, 0])


def main():
    print("=" * 92)
    print("  ITEM 16.  SET IDENTIFICATION: WHAT IS RECOVERABLE WHEN THETA IS NOT")
    print("=" * 92)

    th = C.TH_5
    L, m = C.N4.anchor(th, C.KM0, C.KE, "A")
    print(f"\n  (A) THE ONE-SIDED BOUND ON THETA (set identification).")
    print(f"      Untreated balance theta g(m) + k_e m = L: L = {L:.2f}, m = {m:.2f},")
    print(f"      L/m = {L/m:.3f}/d. theta_max = (L - k_e,lo m)(K_m,hi + m)/m.")
    kms = (100.0, 150.0, 500.0, 1e4)
    print(f"      {'k_e floor':>34}" + "".join(f"{'Km<=' + format(k, 'g'):>12}" for k in kms))
    for lab, kl in (("k_e >= 0 (non-negativity)", 0.0),
                    (f"k_e >= {KE_LO_95:.2f} (prior 95% low)", KE_LO_95)):
        print(f"      {lab:>34}" + "".join(
            f"{100*T4.theta_max(L, m, k, kl)/C.TH0:11.1f}%" for k in kms))
    print(f"      With K_m from the record (pinned to ~6%), the usable bound is")
    print(f"      theta <= 30.7% (non-negativity) / 22.4% (k_e prior). The LOWER")
    print(f"      bound is 0: no record rules out a near-complete deficiency.")
    print(f"      => The defensible deliverable is a CEILING on capacity, not a point.")

    print(f"\n  (B) WHAT THE MONTHLY RECORD PINS, AND WHAT IT DOES NOT (median/24 seeds).")
    V = np.array([recoverable(C.record(th, seed=s)) for s in C.SEEDS])
    sth, sk, sf, ss, sa, sb, sc, v0, v1 = np.median(V, 0)
    print(f"      NOT recoverable (all proportional to theta):")
    print(f"        sd(ln theta)                       = {sth:.3f}")
    print(f"        sd(ln a),   a = theta g(m)^2       = {sa:.3f}")
    print(f"        sd(ln flux theta*g(m0), untreated) = {sf:.3f}")
    print(f"        sd(ln s),   s = fraction oxidised  = {ss:.3f}")
    print(f"      RECOVERABLE (independent of theta's ill-conditioning):")
    print(f"        sd(ln K_m)                         = {sk:.3f}   (~6%)")
    print(f"        sd(ln b),   b = theta g'(m)+k_e    = {sb:.3f}   (the relaxation")
    print(f"                    rate lambda at the operating point, ~{100*sb:.0f}%)")
    print(f"        best (ln theta, ln k_e) combination= {sc:.3f}   along")
    print(f"                    [{v0:+.2f} ln theta {v1:+.2f} ln k_e] (k_e-dominated)")
    print(f"      So the record identifies the LOCAL CLEARANCE CURVE (its slope b and")
    print(f"      K_m), hence k_e, but not the capacity theta: the enzyme's own")
    print(f"      curvature a = theta g^2 is swamped at this operating point. Neither")
    print(f"      the flux nor the fraction oxidised is a reparameterisation that")
    print(f"      escapes theta -- they ARE theta, up to the pinned g and L.")

    print(f"\n  (C) COMPLETE-DEFICIENCY CHECK (asked in Task 5, answered here).")
    ths = sorted({C.record(th, seed=s)["theta"] for s in C.SEEDS})
    print(f"      distinct theta across the 24 records: {ths} -- all at 5% of normal.")
    print(f"      NONE is a complete-deficiency (theta ~ 0) patient: the repository's")
    print(f"      records are 24 intake draws of ONE simulated 5% patient (README).")
    print(f"      If a theta ~ 0 patient existed at an untreated steady state, then")
    print(f"      theta g ~ 0 and the balance gives k_e = L/m directly, with no")
    print(f"      confounding -- a tight k_e prior at no procedural cost. But in the")
    print(f"      parent model (models.ParamsC) a theta = 0 patient has NO steady")
    print(f"      state (R0 is set so the alternative route alone peaks at 500 against")
    print(f"      a baseline load of 600), so there is no such record to read, and the")
    print(f"      counterfactual cannot be realised within this model.")


if __name__ == "__main__":
    main()
