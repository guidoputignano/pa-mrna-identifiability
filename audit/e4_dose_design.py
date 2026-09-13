"""
3.8 -- is the dose train itself a design variable?

Everything so far optimised the SAMPLING schedule and took the input as given.
This asks whether the input is worth varying: amplitude across cycles, or the
interval between them.

The identifiability figure of merit is sd(ln theta) on the minimal footing
(theta and K_m free, priors only on ln k_e and the load record). Alongside it
two therapeutic proxies are reported, because an identifiability optimum that
lowers the dose is not a recommendation:

    mean Delta   the average enzymatic increment delivered
    trough       the MINIMUM of Delta over the record, which is the quantity a
                 dosing interval is chosen to protect
"""
import sys, numpy as np
sys.path.insert(0, 'src')
import design.run as D, model_A.run as LOC

P = D.P
L0 = D.load(0.0)
SIG = LOC.sigma_total(D.LOAD_FRAC, D.ASSAY_CV)
LAM = np.log(2.0) / 3.5
T = 730.0


def train(t, dose_times, dose_amps):
    """General dose train: any amplitudes at any times, first-order decay."""
    t = np.atleast_1d(t)[:, None]
    dt = t - np.asarray(dose_times)[None, :]
    contrib = np.where(dt >= 0, np.asarray(dose_amps)[None, :] * np.exp(-LAM * np.clip(dt, 0, None)), 0.0)
    return contrib.sum(axis=1)


def sd_theta(samp_t, dose_times, dose_amps, load_cv=0.10, seed=3):
    inc = train(samp_t, dose_times, dose_amps)
    rng = np.random.default_rng(seed)
    Ls = L0 * (1.0 + load_cv * rng.standard_normal(samp_t.size)) if load_cv else np.full(samp_t.size, L0)
    B = P.ke * P.Km + P.theta + inc - Ls
    lv = (-B + np.sqrt(B ** 2 + 4.0 * P.ke * Ls * P.Km)) / (2.0 * P.ke)
    pri = [np.inf, np.inf, (D.LOAD_FRAC * L0) ** 2, 0.8 ** 2]
    F = np.diag([0.0 if np.isinf(x) else 1.0 / x for x in pri])
    for mm, dd in zip(lv, inc):
        q = np.array([P.theta * mm / (P.Km + mm),
                      -(P.theta + dd) * P.Km * mm / (P.Km + mm) ** 2,
                      1.0, P.ke * mm])
        F = F + np.outer(q, q) / SIG ** 2
    w = np.linalg.eigvalsh(F)
    if w[0] <= 0 or w[0] <= 1e-14 * w[-1]:
        return np.inf, inc
    v = np.linalg.inv(F)[0, 0]
    return (float(np.sqrt(v)) if v > 0 else np.inf), inc


def report(rows, head):
    print("=" * 92)
    print(head)
    print("=" * 92)
    print(f"{'schedule':>40} {'sd(ln theta)':>13} {'mean Delta/th':>14} "
          f"{'trough/th':>11} {'verdict':>9}")
    base = None
    for lab, samp, dts, amps in rows:
        s, inc = sd_theta(samp, dts, amps)
        full = train(np.linspace(0, T, 4000), dts, amps)
        if base is None:
            base = s
        v = "--" if abs(s / base - 1) < 0.02 else ("better" if s < base else "worse")
        print(f"{lab:>40} {s:13.4f} {full.mean()/P.theta:14.3f} "
              f"{full.min()/P.theta:11.3f} {v:>9}")
    return base


def main():
    print(__doc__)
    samp = np.arange(0.0, T, 30.4)                     # monthly, the recommendation
    Dn, Tn = D.D_AMP, D.DOSE_INTERVAL
    ref_times = np.arange(0.0, T, Tn)
    ref_amps = np.full(ref_times.size, Dn)

    rows = [("fixed: 6000 every 14 d (the trial schedule)", samp, ref_times, ref_amps)]
    for frac in (0.25, 0.50, 0.75):
        a = np.where(np.arange(ref_times.size) % 2 == 0, Dn * (1 + frac), Dn * (1 - frac))
        rows.append((f"alternating +/-{frac:.0%}, same mean dose", samp, ref_times, a))
    rng = np.random.default_rng(0)
    a = Dn * (1 + 0.5 * rng.standard_normal(ref_times.size))
    rows.append(("random amplitude, sd 50%, same mean", samp, ref_times, np.clip(a, 0, None)))
    report(rows, "AMPLITUDE MODULATION, at constant mean dose")

    print()
    rows2 = [("fixed: 6000 every 14 d", samp, ref_times, ref_amps)]
    for iv in (7.0, 10.0, 21.0, 28.0):
        tt = np.arange(0.0, T, iv)
        # hold the DOSE RATE constant so the therapy is not being cut or raised
        rows2.append((f"every {iv:.0f} d, dose scaled to hold the rate", samp, tt,
                      np.full(tt.size, Dn * iv / Tn)))
    report(rows2, "DOSE INTERVAL, at constant delivered dose RATE")

    print()
    print("=" * 92)
    print("THE THERAPEUTIC CONSTRAINT")
    print("=" * 92)
    print("  A dosing interval is chosen to protect the TROUGH: the lowest enzymatic")
    print("  increment the patient sees between doses. Read the trough column, not")
    print("  the mean. Any schedule that lowers it is a worse therapy whatever it")
    print("  does for identifiability, and is not recommended here.")


if __name__ == "__main__":
    main()
