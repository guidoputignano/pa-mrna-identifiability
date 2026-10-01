"""
Adversarial recheck (round 3) of ITEM 16 (out/t18_set_id.txt).

Independent of t18_set_id.py: imports ONLY core primitives (records, the
steady-state solve, SEEDS, the anchor and the constants). Everything that
t18 computes is rebuilt here from scratch:

  * theta_max from the balance algebra written out by hand;
  * the monthly-record FIM assembled from my own flux sensitivities and my
    own per-draw noise model (then cross-checked that it equals what core
    assembles, so the propagation below rests on a verified matrix);
  * the untreated steady state solved as the positive root of the quadratic
    (verified against core.m_ss), so the FD gradients of the derived scalars
    share no code path with core.lnm_rows_fd / t18.recoverable;
  * the distinct-theta / complete-deficiency check.

Run:  python3 -W error t18_set_id_recheck3.py   (output saved under out/)
"""
import numpy as np

import core as C

LF, AC = C.LOAD_FRAC, C.ASSAY_CV
KE_LO_95 = C.KE * np.exp(-1.96 * C.KE_PRIOR)   # 95% low end of the ln k_e prior
OUT = []


def log(s=""):
    OUT.append(s)
    print(s)


# ------------------------------------------------------------------ solvers
def my_mss(theta, Km, ke, L):
    """Positive root of theta*m/(Km+m) + ke*m = L (untreated balance),
    solved as a quadratic independently of core.level."""
    bq = ke * Km + theta - L
    return (-bq + np.sqrt(bq * bq + 4.0 * ke * L * Km)) / (2.0 * ke)


def theta_max(L, m, Km_hi, ke_lo):
    """Hand-written one-sided bound (L - ke_lo m)(Km_hi + m)/m."""
    return (L - ke_lo * m) * (Km_hi + m) / m


def my_fim(r):
    """Monthly-record information matrix in (ln theta, ln K_m, b, ln k_e),
    perdraw noise, b and k_e priors -- assembled from scratch."""
    th, K, ke, m, inc, L = (r["theta"], r["Km"], r["ke"],
                            r["m"], r["inc"], r["L"])
    q = np.stack([th * m / (K + m),
                  -(th + inc) * K * m / (K + m) ** 2,
                  np.ones_like(m),
                  ke * m], axis=1)
    slope = (th + inc) * K / (K + m) ** 2 + ke
    sig = np.hypot(LF * L, slope * AC * m)
    F = (q / sig[:, None] ** 2).T @ q
    pri = np.zeros(4)
    pri[2] = (LF * L) ** -2          # load-bias prior, 5% of L in flux units
    pri[3] = C.KE_PRIOR ** -2        # ln k_e prior, sd 0.8
    return F + np.diag(pri)


# --------------------------------------------------------- derived scalars
def quants(x, Lr):
    """flux, share, a=theta g^2, b=theta g'+k_e at untreated balance; the
    steady state is RE-SOLVED for each coordinate vector x."""
    t, k, b_, e = np.exp(x[0]), np.exp(x[1]), x[2], np.exp(x[3])
    m = my_mss(t, k, e, Lr - b_)
    g, gp = m / (k + m), k / (k + m) ** 2
    return dict(flux=t * g, share=t * g / (Lr - b_), a=t * g * g, b=t * gp + e)


def recover(r, h=1e-6):
    P = np.linalg.inv(my_fim(r))
    th, K, ke, Lr = r["theta"], r["Km"], r["ke"], r["L"]
    base = np.array([np.log(th), np.log(K), 0.0, np.log(ke)])
    keys = ["flux", "share", "a", "b"]
    G = {q: np.zeros(4) for q in keys}
    for j in range(4):
        hj = h * (Lr if j == 2 else 1.0)
        xp, xm = base.copy(), base.copy()
        xp[j] += hj
        xm[j] -= hj
        qp, qm = quants(xp, Lr), quants(xm, Lr)
        for q in keys:
            G[q][j] = (np.log(qp[q]) - np.log(qm[q])) / (2 * hj)
    sd = lambda g: float(np.sqrt(g @ P @ g))
    sub = P[np.ix_([0, 3], [0, 3])]
    w, V = np.linalg.eigh(sub)
    return (float(np.sqrt(P[0, 0])), float(np.sqrt(P[1, 1])),
            sd(G["flux"]), sd(G["share"]), sd(G["a"]), sd(G["b"]),
            float(np.sqrt(w[0])), float(V[0, 0]), float(V[1, 0]))


# --------------------------------------------------------------------- run
def main():
    log("=" * 88)
    log("  RECHECK3  ITEM 16  (independent of t18_set_id.py)")
    log("=" * 88)

    th = C.TH_5
    L, m = C.N4.anchor(th, C.KM0, C.KE, "A")

    # cross-checks that the rebuilt primitives match core
    log("\n  [0] primitive cross-checks (my code vs core):")
    log(f"      my_mss vs core.m_ss (untreated):  {my_mss(th, C.KM0, C.KE, L):.9f}"
        f"  {float(C.m_ss(th, C.KM0, C.KE, L)):.9f}")
    r0 = C.record(th, seed=0)
    log(f"      my_fim == core.fim(perdraw):      "
        f"{bool(np.allclose(my_fim(r0), C.fim(r0, noise='perdraw')))}")

    # ---- (A) the one-sided bound
    log("\n  (A) ONE-SIDED BOUND  theta_max = (L - k_e,lo m)(K_m,hi + m)/m")
    log(f"      L = {L:.2f}  m = {m:.2f}  L/m = {L/m:.3f}/d  (k_e,lo prior95 = {KE_LO_95:.3f})")
    kms = (100.0, 150.0, 500.0, 1e4)
    log(f"      {'k_e floor':>30}" + "".join(f"{'Km<=' + format(k, 'g'):>11}" for k in kms))
    for lab, kl in (("k_e >= 0", 0.0), (f"k_e >= {KE_LO_95:.2f}", KE_LO_95)):
        log(f"      {lab:>30}" + "".join(
            f"{100 * theta_max(L, m, k, kl) / C.TH0:10.1f}%" for k in kms))
    log(f"      claim: Km<=100 -> 30.7% (k_e>=0), 25.3% (k_e>=1.04).")
    log(f"      mine : Km<=100 -> {100*theta_max(L,m,100.0,0.0)/C.TH0:.1f}% (k_e>=0), "
        f"{100*theta_max(L,m,100.0,KE_LO_95)/C.TH0:.1f}% (k_e>={KE_LO_95:.2f}).")

    # ---- (B) recoverable vs not
    log("\n  (B) WHAT THE MONTHLY RECORD PINS (median over 24 seeds)")
    Vv = np.array([recover(C.record(th, seed=s)) for s in C.SEEDS])
    sth, sk, sf, ss, sa, sb, sc, v0, v1 = np.median(Vv, 0)
    log(f"      NOT recoverable:")
    log(f"        sd(ln theta)         = {sth:.3f}   (claim 3.045)")
    log(f"        sd(ln a=theta g^2)   = {sa:.3f}   (claim 3.131)")
    log(f"        sd(ln flux=theta g)  = {sf:.3f}   (claim 3.085)")
    log(f"        sd(ln share)         = {ss:.3f}   (claim 3.078)")
    log(f"      RECOVERABLE:")
    log(f"        sd(ln K_m)           = {sk:.3f}   (claim 0.062)")
    log(f"        sd(ln b=theta g'+ke) = {sb:.3f}   (claim 0.277)")
    log(f"        best (lnth,lnke) sd  = {sc:.3f}   (claim 0.084)")
    log(f"        combo direction      = [{v0:+.2f} ln theta {v1:+.2f} ln k_e]  (claim -0.22/-0.98)")
    log(f"      |k_e weight| {abs(v1):.2f} vs |theta weight| {abs(v0):.2f}  ->"
        f" {'k_e-dominated' if abs(v1) > abs(v0) else 'theta-dominated'}")
    # spread across seeds, to show it is not a one-seed fluke
    lo, hi = Vv.min(0), Vv.max(0)
    log(f"      ranges over seeds: sd(lnK_m) [{lo[1]:.3f},{hi[1]:.3f}]  "
        f"sd(ln b) [{lo[5]:.3f},{hi[5]:.3f}]  best [{lo[6]:.3f},{hi[6]:.3f}]  "
        f"sd(ln theta) [{lo[0]:.2f},{hi[0]:.2f}]")

    # ---- (C) complete-deficiency check
    log("\n  (C) COMPLETE-DEFICIENCY CHECK")
    ths = sorted({float(C.record(th, seed=s)["theta"]) for s in C.SEEDS})
    log(f"      distinct theta across 24 records: {ths}")
    log(f"      count of records: {len(C.SEEDS)}  all theta == 150 (5% of {C.TH0:.0f}): "
        f"{all(C.record(th, seed=s)['theta'] == 150.0 for s in C.SEEDS)}")
    log(f"      any theta ~ 0 (complete deficiency)? "
        f"{any(C.record(th, seed=s)['theta'] < 1e-6 for s in C.SEEDS)}")

    with open("out/t18_set_id_recheck3.txt", "w") as f:
        f.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    main()
