"""
Adversarial recheck #3 of ITEM 13 (out/t15_load_asymptote.txt).

Imports ONLY core primitives (and rebuilds the record from n4's footing via
core.N4 / core.level / core.delta_at). Does NOT import t15_load_asymptote.

Everything downstream of the record -- the sensitivity rows, the per-draw
flux-space sigma, the information matrix and its prior block -- is rebuilt here
from the model equations rather than taken from core.flux_rows / core.sigmas /
core.fim, so the check shares no FIM-assembly code with the script under test.

Record convention follows t15: times = linspace(0, 730, n), theta = 5% of
normal, k_e made "exact" by a tight ln-prior (default 1e-4), bias prior b and
per-draw load scatter separated so each can be pinned on its own.
"""
import numpy as np

import core as C

TH = C.TH_5
L0 = float(C.N4.anchor(TH, C.KM0, C.KE, "A")[0])   # recorded load 322.73
SEEDS = range(8)


# ----------------------------------------------- independent record (no core.record)
def rec(n, seed, D=C.D0, Km=C.KM0, ke=C.KE, load_cv=0.10):
    """A simulated intake record, built straight from n4's anchor/level/delta."""
    L, m_bar = C.N4.anchor(TH, Km, ke, "A")
    times = np.linspace(0.0, 730.0, int(n))
    rng = np.random.default_rng(seed)
    s2 = np.log1p(load_cv ** 2)
    loads = L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
    inc = C.delta_at(times, D) if D else np.zeros(times.size)
    m = C.level(TH, Km, ke, loads, inc)
    return dict(theta=TH, Km=Km, ke=ke, L=float(L), m_bar=float(m_bar),
                loads=loads, inc=inc, m=m, times=times)


# ----------------------------------------------- independent FIM assembly
def Qrows(r):
    th, K, ke, m, inc = r["theta"], r["Km"], r["ke"], r["m"], r["inc"]
    return np.stack([th * m / (K + m),
                     -(th + inc) * K * m / (K + m) ** 2,
                     np.ones_like(m),
                     ke * m], axis=1)


def slope(r):
    K = r["Km"]
    return (r["theta"] + r["inc"]) * K / (K + r["m"]) ** 2 + r["ke"]


def sigma_perdraw(r, lf, ac):
    return np.hypot(lf * r["L"], slope(r) * ac * r["m"])


def fim(r, ke_prior, b_prior, lf, ac):
    Q = Qrows(r)
    sg = sigma_perdraw(r, lf, ac)
    F = (Q / sg[:, None] ** 2).T @ Q
    pri = np.zeros(4)
    if b_prior:
        pri[2] = b_prior ** -2
    if ke_prior:
        pri[3] = ke_prior ** -2
    return F + np.diag(pri)


def sdth(F):
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf
    return float(np.sqrt(np.linalg.inv(F)[0, 0]))


def sd(n, ske=1e-4, b_prior=None, lf=C.LOAD_FRAC, ac=C.ASSAY_CV, seeds=SEEDS,
       bias_ties_lf=True):
    """
    Median sd(ln theta) over seeds.  b_prior semantics:
      - an explicit value is used as given (flux-unit sd of the bias),
      - None reproduces core.fim's default, b_prior = lf * L  (note: this ties
        the bias prior to the scatter knob -- the subtlety part (1b) probes).
    """
    if b_prior is None and bias_ties_lf:
        bp = lf * L0
    elif b_prior is None:
        bp = C.LOAD_FRAC * L0
    else:
        bp = max(b_prior, 1e-9)
    vals = []
    for s in seeds:
        vals.append(sdth(fim(rec(n, s), ske, bp, lf, ac)))
    vals = np.array(vals)
    vals = vals[np.isfinite(vals)]
    return float(np.median(vals)) if vals.size else np.inf


def main():
    out = []
    p = lambda s="": out.append(s)
    p("=" * 92)
    p("  ITEM 13 RECHECK #3 -- load-record term: BIAS vs SCATTER, and the n->inf limit")
    p("  (independent record from n4 footing; FIM rebuilt from model equations)")
    p("=" * 92)
    p(f"  L (recorded load) = {L0:.4f};  seeds = {list(SEEDS)}")

    # ------------------------------------------------------------------ (a) n->inf
    p("\n  (a) sd(ln theta) vs n, k_e EXACT (ske=1e-4), bias 5% prior.  Does it floor?")
    p(f"      {'n':>8}  {'sd (b 5%)':>10}  {'sd (b pinned)':>13}  {'sd*sqrt(n) b5%':>14}  {'slope dlnsd/dlnn':>16}")
    NS = (25, 100, 365, 730, 1460, 5000, 20000, 100000)
    prev = None
    tab = {}
    for n in NS:
        a = sd(n)                       # bias default = lf*L = 0.05*L
        b = sd(n, b_prior=1e-7 * L0)    # bias pinned
        tab[n] = a
        sl = "" if prev is None else f"{(np.log(a)-np.log(prev[1]))/(np.log(n)-np.log(prev[0])):>16.3f}"
        p(f"      {n:>8}  {a:>10.4f}  {b:>13.4f}  {a*np.sqrt(n):>14.3f}  {sl}")
        prev = (n, a)
    p(f"      claim check: n=25 ~0.38 -> {tab[25]:.4f};  n=730 ~0.15 -> {tab[730]:.4f};"
      f"  n=20000 ~0.03 -> {tab[20000]:.4f}")

    # push further out to see whether sd*sqrt(n) truly settles (no log blow-up)
    p("\n      further out (fewer seeds) to confirm sd*sqrt(n) settles, not climbs:")
    for n in (300000, 1000000):
        a = sd(n, seeds=range(3))
        p(f"      {n:>8}  sd={a:>10.5f}  sd*sqrt(n)={a*np.sqrt(n):>8.3f}")

    # ------------------------------------------------------------------ (b) decomposition
    p("\n  (b) DECOMPOSITION AT n=25 (k_e exact, K_m free): is BIAS or SCATTER binding?")
    full = sd(25, lf=C.LOAD_FRAC, b_prior=C.LOAD_FRAC * L0)      # scatter 5%, bias 5%
    scat_off_bias5 = sd(25, lf=1e-6, b_prior=C.LOAD_FRAC * L0)   # HONEST: scatter off, bias kept 5%
    scat5_bias_pin = sd(25, lf=C.LOAD_FRAC, b_prior=1e-7 * L0)   # scatter 5%, bias pinned
    both_off = sd(25, lf=1e-6, b_prior=1e-7 * L0)               # assay + K_m only
    # the script's as-implemented "no_scatter": lf=1e-6 with DEFAULT bias (=lf*L, pinned!)
    script_no_scatter = sd(25, lf=1e-6, b_prior=None)           # bias default tied to lf -> pinned
    p(f"      four honest corners (scatter x bias):")
    p(f"        scatter 5%, bias 5% prior   (full)          : {full:.4f}")
    p(f"        scatter OFF, bias 5% prior  (remove scatter): {scat_off_bias5:.4f}")
    p(f"        scatter 5%, bias pinned     (remove bias)   : {scat5_bias_pin:.4f}")
    p(f"        scatter OFF, bias pinned    (both off)      : {both_off:.4f}")
    p(f"      marginal value of removing SCATTER (keep bias 5%): {full - scat_off_bias5:+.4f}")
    p(f"      marginal value of removing BIAS    (keep scatter): {full - scat5_bias_pin:+.4f}")
    bind = "SCATTER" if (full - scat_off_bias5) > (full - scat5_bias_pin) else "BIAS"
    p(f"      -> more-binding single term at n=25: {bind}")
    p(f"      script's as-implemented 'remove scatter' (lf=1e-6, default bias): {script_no_scatter:.4f}")
    p(f"         (note: default bias prior = lf*L, so lf=1e-6 ALSO pins the bias;")
    p(f"          the script's 0.188 is the BOTH-OFF corner, not scatter-only,")
    p(f"          since scatter-only-with-bias-kept = {scat_off_bias5:.4f})")

    # ------------------------------------------------------------------ (c) measuring L
    p("\n  (c) COST OF MEASURING L (tighten bias prior to X%), k_e exact vs k_e nuisance 0.2")
    p(f"      {'X% on L':>10}  {'n=25 k_e exact':>14}  {'n=25 k_e=0.2':>13}  {'n=365 k_e exact':>15}")
    for X, lab in ((0.05, "5% (prior)"), (0.02, "2%"), (0.01, "1%"), (1e-7, "exact")):
        bp = max(X * L0, 1e-9)
        c_exact25 = sd(25, ske=1e-4, b_prior=bp)
        c_ke02_25 = sd(25, ske=0.2, b_prior=bp)
        c_exact365 = sd(365, ske=1e-4, b_prior=bp)
        p(f"      {lab:>10}  {c_exact25:>14.4f}  {c_ke02_25:>13.4f}  {c_exact365:>15.4f}")
    # explicit direction/magnitude summary for part (c)
    base25 = sd(25, ske=1e-4, b_prior=0.05 * L0)
    tight25 = sd(25, ske=1e-4, b_prior=0.01 * L0)
    base_ke02 = sd(25, ske=0.2, b_prior=0.05 * L0)
    tight_ke02 = sd(25, ske=0.2, b_prior=0.01 * L0)
    p(f"      k_e exact, n=25: 5% -> 1% moves {base25:.3f} -> {tight25:.3f} "
      f"(drop {base25-tight25:+.3f})")
    p(f"      k_e=0.2,  n=25: 5% -> 1% moves {base_ke02:.3f} -> {tight_ke02:.3f} "
      f"(drop {base_ke02-tight_ke02:+.3f}; ~nothing)")

    text = "\n".join(out)
    print(text)
    import pathlib
    pathlib.Path(C.HERE / "out" / "t15_load_asymptote_recheck3.txt").write_text(text + "\n")


if __name__ == "__main__":
    main()
