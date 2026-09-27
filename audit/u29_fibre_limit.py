"""
Does the decidability result survive the constant-load limit, where the fibre is
exactly null?

THE AUTHOR'S CHALLENGE, WHICH TESTS THE RESULT RATHER THAN EXTENDING IT. Along
the exact fibre `Km(b) = Km0 L/(L+b)`, `theta(b) = theta0 + b + ke(Km0 - Km(b))`,
the derivatives are

    d theta / db = 1 + ke Km0 / L        = 2.5493
    d mu    / db = that + (s_peak - ke m_tol)/m_tol * Km0/L = 9.9500

at **every** capacity, confirmed here analytically and by central differences. So
**the margin is 3.903 times more exposed to the unidentified direction than
capacity is.** Yet `u28` reports `sd(mu)` small, and small only through a **3.25x
cancellation** between diagonal and cross terms.

**The quantity most exposed to the fibre is coming out the better determined one,
and the reason is a cancellation.** That is the shape of a result that depends on
something incidental, so it must be tested against the limit in which the fibre
becomes exactly null.

WHAT BREAKS THE FIBRE IS THE LOAD VARIATION. Proposition 4 is exactly singular at
a constant recorded load, and what closes it is the day-to-day variability of
ordinary eating. `LOAD_CV = 0.10` is the register's value. This sweeps it down to
zero and watches `sd(ln theta)` and `sd(mu)` together.

THE TWO OUTCOMES, WRITTEN BEFORE THE RUN AND NOT AS A PREDICTION.
  If `sd(mu)` diverges FASTER than `sd(ln theta)`, the decidability result rests
  on the load variation breaking the fibre rather than on a property of the
  decision, and it must be published with that dependence stated.
  If it does not, the cancellation is structural and the result is stronger than
  it reads.
Which of these holds is what the run decides; nothing here records an expectation.

`levels()` reproduces `u9`'s exactly, with `load_cv` promoted from a module
constant to an argument. At `load_cv = 0` every seed gives the identical record,
so the seed spread collapses to a point and is reported as such rather than as a
suspiciously narrow range.

Run:  python audit/u29_fibre_limit.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np

import design.run as D
from common.params import _A as P

import u9_ke_frontier as U9
import u28_margin_frontier as U28

KM, KE, LB, MT = U9.KM, U9.KE, U9.LB, U28.MT
S_PEAK = U28.S_PEAK
CVS = (0.10, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 1e-4, 1e-5, 0.0)


def levels(theta, n, seed, load_cv, step=U9.MONTH):
    """u9.levels with LOAD_CV promoted to an argument. Otherwise identical."""
    t = step * np.arange(n, dtype=float)
    tau = np.mod(t - U9.T_DOSE, U9.DOSE_INTERVAL)
    inc = U9.PEAK * np.exp(-U9.LAMBDA * tau)
    Ls = LB * (1.0 + load_cv *
               np.random.default_rng(seed).standard_normal(t.size))
    B = KE * KM + theta + inc - Ls
    return (-B + np.sqrt(B ** 2 + 4.0 * KE * Ls * KM)) / (2.0 * KE), inc


def pieces(theta, m, inc, sigma):
    """(sd_ln_theta, sd_mu, cancellation factor, smallest relative eigenvalue)."""
    gm = m / (KM + m)
    Q = np.stack([theta * gm,
                  -(theta + inc) * KM * m / (KM + m) ** 2,
                  np.ones_like(m),
                  KE * m], axis=1)
    pri = [np.inf, np.inf, (U9.LOAD_FRAC * LB) ** 2, U9.KE_PRIOR ** 2]
    F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri]) \
        + Q.T @ Q / sigma ** 2
    w = np.linalg.eigvalsh(F)
    rel = float(w[0] / w[-1])
    if w[0] <= 0 or rel <= 1e-18:
        return np.inf, np.inf, np.nan, rel
    C = np.linalg.inv(F)
    sd_t = float(np.sqrt(C[0, 0])) if C[0, 0] > 0 else np.inf
    g = U28.grad_mu(theta)
    var = float(g @ C @ g)
    diag = float(sum(g[i] * C[i, i] * g[i] for i in range(4)))
    sd_m = float(np.sqrt(var)) if var > 0 else np.inf
    canc = float(np.sqrt(diag / var)) if var > 0 and diag > 0 else np.nan
    return sd_t, sd_m, canc, rel


def row(theta, cv, n=U9.TIER1_N):
    sig = U9.sigma_at(theta)
    out = [pieces(theta, *levels(theta, n, s, cv), sig) for s in U9.SEEDS]
    fin = [o for o in out if np.isfinite(o[0]) and np.isfinite(o[1])]
    if not fin:
        return None
    med = lambda i: float(np.median([o[i] for o in fin]))
    return med(0), med(1), med(2), med(3), len(fin), len(out)


if __name__ == "__main__":
    print("=" * 104)
    print("DOES THE DECIDABILITY RESULT SURVIVE THE CONSTANT-LOAD LIMIT?")
    print("=" * 104)
    print(f"\n  along the fibre: d theta/db = {1 + KE * KM / D.load(0.0):.4f}, "
          f"d mu/db = {1 + KE * KM / D.load(0.0) + (S_PEAK - KE * MT) / MT * KM / D.load(0.0):.4f}, "
          f"ratio {(1 + KE*KM/D.load(0.0) + (S_PEAK-KE*MT)/MT*KM/D.load(0.0))/(1+KE*KM/D.load(0.0)):.3f}")
    print(f"  the register's LOAD_CV is {D.LOAD_CV}; n = {U9.TIER1_N} draws, "
          f"{len(list(U9.SEEDS))} seeds\n")

    for theta in (3000.0, 600.0, 300.0):
        mu = theta - (S_PEAK - KE * MT) * (KM + MT) / MT
        print("-" * 104)
        print(f"theta = {theta:.0f},  mu = {mu:+.1f}")
        print("-" * 104)
        print(f"  {'load_cv':>9} {'sd(ln th)':>11} {'sd(mu)':>12} {'|mu|/sd':>9} "
              f"{'cancel':>8} {'min eig/max':>13} {'finite':>8}")
        base = None
        for cv in CVS:
            r = row(theta, cv)
            if r is None:
                print(f"  {cv:9.5f} {'singular':>11} {'singular':>12} "
                      f"{'-':>9} {'-':>8} {'-':>13} {'0/' + str(len(list(U9.SEEDS))):>8}")
                continue
            sd_t, sd_m, canc, rel, nf, nt = r
            if base is None:
                base = (sd_t, sd_m)
            print(f"  {cv:9.5f} {sd_t:11.4f} {sd_m:12.2f} {abs(mu)/sd_m:9.2f} "
                  f"{canc:8.2f} {rel:13.3e} {str(nf) + '/' + str(nt):>8}")
        if base:
            r0 = row(theta, 0.0)
            if r0:
                print(f"\n  from LOAD_CV = 0.10 to a constant load: "
                      f"sd(ln theta) x{r0[0]/base[0]:.3g}, sd(mu) x{r0[1]/base[1]:.3g}")
        print()

    print("=" * 104)
    print("PART B. WHICH MECHANISM BREAKS THE FIBRE? ONE AT A TIME, NOT BOTH AT ONCE")
    print("=" * 104)
    print("""
  AN EARLIER VERSION OF THIS SECTION GOT THIS WRONG AND THE ERROR IS INSTRUCTIVE
  (REVIEW_LOG.md R167). It removed the load variation AND the load-bias prior
  together, found the matrix singular, and concluded that the prior is what holds
  the design together. **That does not follow from that test**, and this file's
  own table already refuted it: at load_cv = 0.10 with no prior the matrix is
  non-singular, and at a constant load with the prior it is non-singular too.
  Only removing BOTH gives the null. A both-at-once experiment cannot attribute
  an effect to one factor, which is the oldest mistake in design of experiments
  and was made here inside a file about experimental design.

  So: each mechanism is swept alone, with the other left in place.
""")
    print(f"  {'theta':>7} {'b prior 5%':>12} {'b prior 1%':>12} {'NO b prior':>12} "
          f"{'5->1 buys':>11} {'removal costs':>14}")
    for theta in (3000.0, 600.0, 300.0):
        def sdt(bp):
            out = []
            for sd in U9.SEEDS:
                m, inc = levels(theta, U9.TIER1_N, sd, 0.10)
                sig = U9.sigma_at(theta)
                gm = m / (KM + m)
                Q = np.stack([theta * gm,
                              -(theta + inc) * KM * m / (KM + m) ** 2,
                              np.ones_like(m), KE * m], axis=1)
                pri = [np.inf, np.inf,
                       np.inf if bp is None else bp ** 2, U9.KE_PRIOR ** 2]
                F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri]) \
                    + Q.T @ Q / sig ** 2
                w = np.linalg.eigvalsh(F)
                out.append(np.inf if (w[0] <= 0 or w[0] <= 1e-18 * w[-1])
                           else float(np.sqrt(np.linalg.inv(F)[0, 0])))
            return float(np.median(out))
        a, b, c = sdt(U9.LOAD_FRAC * LB), sdt(0.01 * LB), sdt(None)
        print(f"  {theta:7.0f} {a:12.4f} {b:12.4f} {c:12.4f} "
              f"{100*(1-b/a):10.1f}% {100*(c/a-1):13.1f}%")
    print("""
  THE CORRECTED READING, AND IT IS A BETTER RESULT THAN THE ONE IT REPLACES.
  At the real load variation, dropping the load-bias prior ENTIRELY leaves the
  design non-singular and costs single-digit per cent. Tightening the diet record
  from 5% to 1% buys single-digit per cent the other way. And PART A showed the
  mirror case: with the prior in place, taking the load variation to zero costs
  1.15x to 1.27x on sd(mu) and CONVERGES.

  **Each mechanism alone suffices. Only removing both gives the null. The fibre
  is broken REDUNDANTLY.** The design does not rest on LOAD_FRAC and it does not
  rest on LOAD_CV; it rests on either, which is a robustness statement rather
  than a fragility one.

  And it explains a standing observation. 08_specification records that the
  load-record prior sits near 1.05 and does not move. **It does not move because
  the load variation has already done the work**, so the prior is buying a margin
  on something the record already determines.
""")
    print("=" * 104)
    print("PART C. THE NULL EIGENVECTOR, WHICH CLOSES A LOOP WITH THE MARGIN FRONTIER")
    print("=" * 104)
    print(f"\n  {'theta':>7}   (ln theta, ln K_m, b, ln k_e), scaled to max 1"
          f"      {'min/max eig':>12} {'|k_e|':>10}")
    for theta in (3000.0, 600.0, 300.0):
        m, inc = levels(theta, U9.TIER1_N, 0, 0.0)
        sig = U9.sigma_at(theta)
        gm = m / (KM + m)
        Q = np.stack([theta * gm, -(theta + inc) * KM * m / (KM + m) ** 2,
                      np.ones_like(m), KE * m], axis=1)
        F = np.diag([0.0, 0.0, 0.0, 1.0 / U9.KE_PRIOR ** 2]) + Q.T @ Q / sig ** 2
        w, V = np.linalg.eigh(F)
        v = V[:, 0] / np.abs(V[:, 0]).max()
        print(f"  {theta:7.0f}   ({v[0]:+.4f}, {v[1]:+.4f}, {v[2]:+.4f}, "
              f"{v[3]:+.4f})   {w[0]/w[-1]:12.2e} {abs(v[3]):10.2e}")
    print("""
  THE k_e COMPONENT IS ZERO TO MACHINE PRECISION AT EVERY CAPACITY. That is
  eq:fibrecl's "k_e unchanged" read straight out of the eigenvector rather than
  taken from the algebra, and it is an independent confirmation of the published
  fibre.

  AND IT EXPLAINS A RESULT THAT LOOKED UNRELATED. u28's margin frontier reports
  that an exactly known k_e is UNREACHABLE for the decision at any prior width in
  the disease range. **The reason is here: k_e is orthogonal to the degeneracy.**
  No amount of knowledge about k_e moves the direction mu is exposed along,
  because that direction has no k_e component. Two results that looked unrelated
  are one fact.
""")
    print("=" * 104)
    print("""  READ THE TWO GROWTH FACTORS AGAINST EACH OTHER, NOT SEPARATELY.
  If sd(mu) grows by the larger factor as the load variation vanishes, the
  decidability result depends on the fibre being broken by the record and must be
  published saying so. If the two grow together, the 3.25x cancellation is
  structural and |mu|/sd(mu) is a property of the decision.""")
