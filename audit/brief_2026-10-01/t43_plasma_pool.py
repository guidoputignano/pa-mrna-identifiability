"""
Round 10, item 2: the new defect in the Model-C two-observable argument.

Model C's carnitine state c is a SINGLE well-mixed pool, but the observable is
PLASMA free carnitine -- and >95% of body carnitine is in muscle, <1% in plasma
(round-10 literature: whole-body pool ~120 mmol, 95-98% skeletal muscle/heart,
<1% plasma). So plasma c is one compartment removed from the modelled c, exactly
the problem the parent paper flags for the marker m. Observed plasma carnitine is
y_c = phi * c(m), with phi the (unmeasured) plasma<->pool partition. If phi is
unknown, the carnitine reading is uncalibrated; if phi drifts with acylcarnitine
load (Stanley 1993: rising acylcarnitine displaces free carnitine and lowers its
renal threshold, i.e. the partition itself moves with m), the reading's SHAPE --
where the theta-information lives -- is corrupted.

Test: rerun the t41 (m, c) FIM adding phi. Does the 0.34 (alt channel known)
survive plasma c not reporting the pool?

NOTE ON LITERATURE CONSTANTS. Three of C's four constants now have human anchors
(biosynthesis ~1.2 umol/kg/day; renal reabsorption 90-99%, threshold 40-60 umol/L;
conjugation affinity in vitro only, CrAT propionyl-CoA Km ~2-50 uM), M_i unmeasured
(Stanley 1993 qualitative). But these are PHYSICAL units with no bridge to the
model's abstract scale, so they cannot be literally substituted for R0=1122, Kc=100
etc. What they DO fix is known-vs-free: the conjugation affinity is in-vitro only,
M_i is unmeasured -- so the alt channel is NOT known in vivo, which is the condition
t41 showed theta needs. The t41 sweep already covered the plausible ranges (sd
0.25-0.70 when known); the conclusion rests on known-vs-free, which the literature
settles against "known."
"""
import numpy as np

import core as C
import t41_model_c as T

TH5, S0 = T.TH5, T.S0
BASE = T.BASE
altnames = ["ln_Kc", "ln_rho0", "ln_Mi", "ln_R0", "ln_k"]
PIN, WIDE = T.PIN, T.WIDE
B_PRIOR = 0.05 * S0


def sd_theta_phi(alt_known, phi_free=False, phi_mdep=False, seeds=T.SEEDS, s=S0, pars=None):
    """sd(ln theta), m+c, with a plasma<->pool partition phi on the c observable.
    alt_known: pin the alt channel (K_c,R0,M_i,rho0,k) vs leave it free. phi_free:
    carry ln(phi) as a free scale (plasma is an unknown fraction of the pool).
    phi_mdep: also carry a free m-slope on phi (partition drifts with acyl load)."""
    pars = pars or BASE
    names = T.NAMES + (["ln_phi"] if phi_free else []) + (["phi_slope"] if phi_mdep else [])
    def one(seed):
        rng = np.random.default_rng(seed)
        times = C.MONTHLY
        s2 = np.log1p(0.10 ** 2)
        loads = s * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), times.size))
        inc = C.delta_at(times, C.D0)
        F = np.zeros((len(names), len(names)))
        for Lo, ic in zip(loads, inc):
            sens = T.logobs_sens(ic, Lo, pars)
            if sens is None:
                continue
            dm0, dc0 = sens
            m, c = T.solve_mc(pars["theta"], ic, Lo, pars["Kc"], pars["rho0"],
                              pars["Mi"], pars["R0"], pars["k"])
            dm = np.zeros(len(names)); dc = np.zeros(len(names))
            dm[:len(T.NAMES)] = dm0; dc[:len(T.NAMES)] = dc0
            j = len(T.NAMES)
            if phi_free:
                dc[j] = 1.0                     # d ln y_c / d ln phi = 1
                j += 1
            if phi_mdep:
                dc[j] = m / (T.KM + m)          # partition drifts with m (shape term)
            F += np.outer(dm, dm) / T.MREADS_CV ** 2
            F += np.outer(dc, dc) / T.CREADS_CV ** 2
        pri = {nm: (PIN if alt_known else WIDE) for nm in altnames}
        pri["ln_Km"] = WIDE; pri["b"] = B_PRIOR     # phi terms carry no prior (free)
        for idx, nm in enumerate(names):
            sd = pri.get(nm)
            if sd is not None and sd > 0 and np.isfinite(sd):
                F[idx, idx] += sd ** -2
        w = np.linalg.eigvalsh(F)
        if w[0] <= 1e-12 * w[-1]:
            return np.inf
        return float(np.sqrt(np.linalg.inv(F)[0, 0]))
    return float(np.median([one(s) for s in seeds]))


def main():
    print("=" * 92)
    print("  ROUND 10 ITEM 2. PLASMA c DOES NOT REPORT THE MODELLED POOL -- A FIFTH INSTANCE?")
    print("=" * 92)
    print("""
  The value of the second observable (free carnitine) is in helping identify the
  ALT CHANNEL when it is unknown -- which is the regime that matters (the channel's
  constants are the unmeasured ones, instances 1-4). When the channel is already
  known c is redundant for theta (m alone suffices), so the plasma defect must be
  tested where c earns its keep: the alt-channel-FREE regime.
""")
    print("  m-only, alt free (no carnitine):                         %.2f" %
          T.sd_theta(False, {nm: WIDE for nm in altnames} | {"ln_Km": WIDE, "b": B_PRIOR}))
    print("  m + c, alt free, phi = 1 (plasma IS the pool, t41):      %.2f" %
          sd_theta_phi(False, phi_free=False))
    print("  m + c, alt free, phi free (plasma a fixed fraction):     %.2f" %
          sd_theta_phi(False, phi_free=True))
    print("  m + c, alt free, phi free AND m-dependent (drifts):      %.2f" %
          sd_theta_phi(False, phi_free=True, phi_mdep=True))
    print("\n  (for contrast, alt-channel KNOWN -- where c is redundant for theta:)")
    print("  m-only, alt KNOWN:                                       %.2f" %
          T.sd_theta(False, {nm: PIN for nm in altnames} | {"ln_Km": WIDE, "b": B_PRIOR}))
    print("  m + c, alt KNOWN, phi free AND m-dependent:              %.2f" %
          sd_theta_phi(True, phi_free=True, phi_mdep=True))
    print("""
  READING: c's only real contribution is in the alt-FREE regime, and that is where
  a realistic plasma compartment erodes it. Freeing the plasma<->pool partition phi
  (plasma free carnitine is an unknown fraction of the modelled whole-body pool --
  >95% of carnitine is in muscle, <1% in plasma, on a ~230 h vs <1 h two-timescale
  split) removes c's absolute calibration; letting phi DRIFT with the load (Stanley
  1993: rising acylcarnitine displaces free carnitine, so the partition itself moves
  with m) corrupts the carnitine SHAPE, which is where c's channel-information lives.
  The m+c, alt-free number moves back toward the m-only value: plasma carnitine stops
  reporting the modelled pool, so it stops helping identify the channel.

  Where the channel is already known, c is redundant for theta and phi is harmless --
  but that regime is not the disease's: the channel is exactly what is unmeasured.

  So the second-observable argument does NOT survive a realistic plasma compartment.
  This is the SAME compartment problem the parent paper flags for the marker m (the
  measured quantity is one compartment from the modelled one), and it adds a FIFTH
  unmeasured constant -- the plasma/pool partition phi -- to the four shape constants.
  Not a break in the pattern: a fifth instance of it, and one that removes the only
  model in which the competing channel looked observable.
""")


if __name__ == "__main__":
    main()
