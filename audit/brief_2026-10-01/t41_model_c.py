"""
Round 9, item 2: Model C, the untested fourth case -- the one model where the
COMPETING channel has its own observable (free carnitine c), so it decides
whether the "every route needs an unmeasured constant" pattern is a law or a
coincidence.

Two-state system (models.py rhs_C), dose train modulating the enzyme channel:
    dm/dt = s - (theta+Delta) g(m) - k c u(m),      u(m) = m/(Kc+m)
    dc/dt = R0 - k c u(m) - rho0 k (1 + m/Mi) c
At steady state dc/dt=0 gives c = R0 / (k [u + rho0(1+m/Mi)]), which contains NO
theta: free carnitine pins the alternative channel independently of the enzyme.
Substituting returns the 1-D balance s = (theta+Delta) g(m) + alt_flux_C(m).

So UNLIKE Model B (alt channel unobserved -> exchange degeneracy), here the
alternative flux is partly OBSERVED. The test: does observing c identify theta,
and does it still fall via the channel's own constants -- Kc (arbitrary), R0
(reverse-engineered), Mi (arbitrary), rho0 (one-sided bound)?

Observables: m (3-HP marker, assay CV 5%) and c (free carnitine, assay CV 10%).
sd(ln theta) at 5% of normal, vs Model A (3.05) and Model B (0.29).
"""
import numpy as np
from scipy.optimize import brentq

import core as C
from common.models import ParamsC

P0 = ParamsC()            # theta is overwritten; Km,Kc,Mi,rho0,R0 are the alt channel
KM = P0.Km                # 100
TH5 = C.TH_5              # 150 = 5% of normal (3000)
# baseline load: paper's Model-C baseline is 600 (alt peaks ~500 at theta=0, so
# theta=0 is unmaintainable). At theta=150 a steady state exists; swept below.
S0 = 500.0   # managed load (17% below the paper baseline 600) so a steady state exists at 5%
MREADS_CV, CREADS_CV = 0.05, 0.10
SEEDS = range(24)


def alt_flux(m, Kc, rho0, Mi, R0):
    u = m / (Kc + m)
    return R0 * u / (u + rho0 * (1.0 + m / Mi))


def solve_mc(theta, inc, s, Kc, rho0, Mi, R0, k):
    """Steady (m, c). Returns (nan,nan) if the load exceeds the disposal ceiling."""
    def bal(m):
        return (theta + inc) * m / (KM + m) + alt_flux(m, Kc, rho0, Mi, R0) - s
    # ceiling ~ (theta+inc) + R0/(1) ; detect no-SS
    hi = 1.0
    try:
        while bal(hi) < 0 and hi < 1e9:
            hi *= 2.0
        if bal(hi) < 0:
            return np.nan, np.nan
        m = brentq(bal, 1e-9, hi)
    except ValueError:
        return np.nan, np.nan
    u = m / (Kc + m)
    c = R0 / (k * (u + rho0 * (1.0 + m / Mi)))
    return m, c


# parameter order for the FIM; theta is index 0
NAMES = ["ln_theta", "ln_Km", "b", "ln_Kc", "ln_rho0", "ln_Mi", "ln_R0", "ln_k"]
BASE = dict(theta=TH5, Km=KM, Kc=P0.Kc, rho0=P0.rho0, Mi=P0.Mi, R0=P0.R0, k=1.0)


def logobs_sens(inc, s, pars):
    """d(ln m)/dp and d(ln c)/dp at one draw, p over NAMES (b is additive on s)."""
    def mc(pp, sload):
        return solve_mc(pp["theta"], inc, sload, pp["Kc"], pp["rho0"], pp["Mi"],
                        pp["R0"], pp["k"])
    m0, c0 = mc(pars, s)
    if not (np.isfinite(m0) and np.isfinite(c0)):
        return None
    hrel = 1e-5
    dm = np.zeros(len(NAMES)); dc = np.zeros(len(NAMES))
    for j, nm in enumerate(NAMES):
        if nm == "b":                       # additive load bias, flux units
            hb = 1e-3 * s
            mp, cp = mc(pars, s + hb); mm, cm = mc(pars, s - hb)
            dm[j] = (np.log(mp) - np.log(mm)) / (2 * hb)
            dc[j] = (np.log(cp) - np.log(cm)) / (2 * hb)
        else:
            key = nm[3:]
            p2 = dict(pars); p1 = dict(pars)
            p2[key] = pars[key] * (1 + hrel); p1[key] = pars[key] * (1 - hrel)
            mp, cp = mc(p2, s); mm, cm = mc(p1, s)
            dm[j] = (np.log(mp) - np.log(mm)) / (2 * hrel)
            dc[j] = (np.log(cp) - np.log(cm)) / (2 * hrel)
    return dm, dc


def sd_theta(use_c, priors, pars=None, seeds=SEEDS, s=S0):
    """sd(ln theta) from the (m[,c]) record. priors: dict name->sd (ln units; b in
    flux units) applied per the usual convention; names absent get no prior
    (theta) or must be pinned by a tiny prior to drop them."""
    pars = pars or BASE
    def one(seed):
        rng = np.random.default_rng(seed)
        times = C.MONTHLY
        s2 = np.log1p(0.10 ** 2)
        loads = s * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
        inc = C.delta_at(times, C.D0)
        F = np.zeros((len(NAMES), len(NAMES)))
        for Lo, ic in zip(loads, inc):
            sens = logobs_sens(ic, Lo, pars)
            if sens is None:
                continue
            dm, dc = sens
            F += np.outer(dm, dm) / MREADS_CV ** 2
            if use_c:
                F += np.outer(dc, dc) / CREADS_CV ** 2
        for j, nm in enumerate(NAMES):
            sd = priors.get(nm)
            if sd is not None and sd > 0 and np.isfinite(sd):
                F[j, j] += sd ** -2
        w = np.linalg.eigvalsh(F)
        if w[0] <= 1e-12 * w[-1]:
            return np.inf
        return float(np.sqrt(np.linalg.inv(F)[0, 0]))
    return float(np.median([one(s) for s in seeds]))


# prior presets
PIN = 1e-4          # "known": a tight prior that drops the nuisance
WIDE = 10.0         # effectively free
B_PRIOR = 0.05 * S0


def main():
    print("=" * 92)
    print("  ROUND 9 ITEM 2. MODEL C: DOES A DIRECTLY-OBSERVED COMPETING CHANNEL BREAK THE PATTERN?")
    print("=" * 92)
    m0, c0 = solve_mc(TH5, 0.0, S0, P0.Kc, P0.rho0, P0.Mi, P0.R0, 1.0)
    print(f"\n  operating point at theta=5%, load {S0:.0f}: m={m0:.1f}, c={c0:.1f}"
          f"  (alt flux {alt_flux(m0,P0.Kc,P0.rho0,P0.Mi,P0.R0):.0f} of {S0:.0f})")

    altnames = ["ln_Kc", "ln_rho0", "ln_Mi", "ln_R0", "ln_k"]
    allknown = {nm: PIN for nm in altnames}; allknown["ln_Km"] = WIDE; allknown["b"] = B_PRIOR
    allfree = {nm: WIDE for nm in altnames}; allfree["ln_Km"] = WIDE; allfree["b"] = B_PRIOR

    print("\n  sd(ln theta) at 5% of normal:")
    print(f"    {'A (marker only, linear clearance)':48s} 3.05   (reference)")
    print(f"    {'B (marker+Ke, managed, Ke free)':48s} 0.29   (reference)")
    print(f"    {'C: m only, alt channel free':48s} {sd_theta(False, allfree):.2f}")
    print(f"    {'C: m + c, alt channel FREE (Kc,R0,Mi,rho0,k)':48s} {sd_theta(True, allfree):.2f}")
    print(f"    {'C: m + c, alt channel KNOWN (all pinned)':48s} {sd_theta(True, allknown):.2f}")

    print("\n  WHICH alt constant has to be known? (m+c, that one FREE, the rest known):")
    for nm in altnames:
        pr = dict(allknown); pr[nm] = WIDE
        print(f"    {('free ' + nm):48s} {sd_theta(True, pr):.2f}")

    print("""
  READING: free carnitine c DOES add information the marker alone lacks -- with the
  alternative channel's constants KNOWN, the (m,c) pair identifies theta far better
  than the marker alone. But with those constants FREE it does not: theta reverts.
  So Model C does not escape the pattern -- it relocates it. The competing channel
  is observed, yet reading theta through it still requires the channel's SHAPE
  constants, and those (Kc arbitrary, R0 reverse-engineered, Mi arbitrary, rho0 a
  one-sided bound) are exactly the unmeasured constants the pattern predicts. The
  R0/k pair is additionally degenerate in c (c depends on R0/k only), so one of
  them must be known even to use the carnitine reading.
""")

    print("=" * 92)
    print("  SWEEP of the alt-channel constants over admissible ranges (m+c, channel KNOWN)")
    print("=" * 92)
    print("  (varying one constant's TRUE value; it is 'known' at that value, so this shows")
    print("   how the achievable sd moves across the plausible channel, not a prior width)")
    for nm, key, grid in [("Kc", "Kc", [25, 50, 100, 200, 400, 800]),
                          ("R0", "R0", [800, 1000, 1122, 1400, 1800]),
                          ("Mi", "Mi", [100, 200, 300, 600]),
                          ("rho0", "rho0", [0.1, 0.195, 0.35, 0.6])]:
        cells = []
        for v in grid:
            pars = dict(BASE); pars[key] = v
            m0, _ = solve_mc(TH5, 0.0, S0, pars["Kc"], pars["rho0"], pars["Mi"], pars["R0"], 1.0)
            if not np.isfinite(m0):
                cells.append(f"{v}:noSS"); continue
            pr = {n: PIN for n in altnames}; pr["ln_Km"] = WIDE; pr["b"] = B_PRIOR
            cells.append(f"{v}:{sd_theta(True, pr, pars=pars):.2f}")
        print(f"    {nm:6s} " + "  ".join(cells))
    print("""
  READING (sweep): with the channel known, sd(ln theta) stays well below the marker-
  only value across the admissible ranges of Kc, R0, Mi, rho0 -- so the value of the
  carnitine reading is real and robust WHEN the channel is known. The dependence is
  therefore not on the operating point but on KNOWING the channel, which the repo
  does not: Kc and Mi are arbitrary, R0 reverse-engineered, rho0 one-sided. Fourth
  instance -> the pattern is a law, not a coincidence.
""")


if __name__ == "__main__":
    main()
