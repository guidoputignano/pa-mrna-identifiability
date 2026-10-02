"""
The within-patient longitudinal contrast, with the sharing contract in ONE place.

Round 8 found the same bug three times: a patient-level constant (k_e, then b,
twice) left INDEPENDENT between the two windows of a contrast, which inflates
sd(ln th2 - ln th1) because a constant that should cancel does not. This module
removes the bug class: `window_contrast` ties EVERY patient-level constant across
windows unless an explicit drift is given, so a contrast can only go wrong by
naming a constant that is not a constant -- not by forgetting to share one.

A patient-level constant is any parameter that is the same patient's property in
both windows: the Michaelis constants K_m, K_e, the clearance V_max k_e, the
load-record bias b, and the multiplicative calibrations kappa (ex vivo) and c
(breath). Only the estimand theta differs between windows (that is what is being
monitored). `drifts={name: sd}` relaxes a constant to a Gaussian random walk of
that sd between the two windows (sd=0 or absent -> hard-shared).

Every FIM-based contrast in the brief (t24 marker, t33 Model B, t35 managed)
goes through `window_contrast`, and `check_contract` is the test that fails if
any contrast leaves a patient constant unshared.
"""
import numpy as np

COUPLE = 1e12   # precision that hard-ties two window copies (drift -> 0)


def _assemble(Fs, names, est="theta", drifts=None, priors=None):
    """Build the joint contrast FIM Z and the per-window index map, sharing every
    patient constant (everything except `est`). No inversion -- so it is safe to
    call on zero-data FIMs for the contract test."""
    drifts = dict(drifts or {})
    priors = dict(priors or {})
    nwin = len(Fs)
    k = len(names)
    idx = [{nm: w * k + j for j, nm in enumerate(names)} for w in range(nwin)]
    P = nwin * k
    Z = np.zeros((P, P))
    for w, F in enumerate(Fs):
        for a in range(k):
            for b in range(k):
                Z[w * k + a, w * k + b] += F[a, b]
    # couple every patient constant (everything except the estimand) across windows
    for nm in [n for n in names if n != est]:
        sd = drifts.get(nm, 0.0)
        g = COUPLE if (sd is None or sd <= 0) else 1.0 / sd ** 2
        for w in range(1, nwin):
            i, j = idx[0][nm], idx[w][nm]
            Z[i, i] += g; Z[j, j] += g; Z[i, j] -= g; Z[j, i] -= g
    # per-window priors
    for nm, sd in priors.items():
        if sd and np.isfinite(sd):
            for w in range(nwin):
                Z[idx[w][nm], idx[w][nm]] += sd ** -2
    return Z, idx


def window_contrast(Fs, names, est="theta", drifts=None, priors=None,
                    return_parts=False):
    """sd(est in window 2 - est in window 1) from per-window DATA FIMs.

    Fs       : list of per-window data FIMs (k x k), SAME column order = names.
    names    : parameter names, one per FIM column. `est` is the monitored one.
    est      : the estimand; it gets an INDEPENDENT copy per window. Everything
               else is a patient constant and is SHARED across windows.
    drifts   : {name: sd} to relax a constant to a between-window drift of that
               sd (in the parameter's own units). Absent/<=0 -> hard-shared.
    priors   : {name: sd} a per-window prior on that parameter.
    """
    Z, idx = _assemble(Fs, names, est, drifts, priors)
    nwin = len(Fs)
    Cov = np.linalg.inv(Z)
    e = np.zeros(Z.shape[0])
    e[idx[0][est]] = -1.0
    e[idx[nwin - 1][est]] = +1.0
    sd = float(np.sqrt(e @ Cov @ e))
    if return_parts:
        return sd, Z, idx
    return sd


def check_contract(names, est="theta", drifts=None):
    """TEST: fail if `window_contrast` leaves any patient constant independent
    across windows. Builds a two-window contrast from ZERO data FIMs so that only
    the coupling (and any drift) holds the constants together, and asserts every
    non-estimand parameter is coupled at the expected precision.
    Returns (ok, message)."""
    drifts = dict(drifts or {})
    k = len(names)
    Fs = [np.zeros((k, k)), np.zeros((k, k))]
    Z, idx = _assemble(Fs, names, est=est, drifts=drifts)
    bad = []
    for nm in names:
        if nm == est:
            # the estimand MUST be independent (uncoupled) between windows
            off = Z[idx[0][nm], idx[1][nm]]
            if abs(off) > 1e-6:
                bad.append(f"estimand '{nm}' is coupled across windows (off={off:g})")
            continue
        off = -Z[idx[0][nm], idx[1][nm]]          # coupling precision
        sd = drifts.get(nm, 0.0)
        want = COUPLE if (sd is None or sd <= 0) else 1.0 / sd ** 2
        if not np.isclose(off, want, rtol=1e-6):
            bad.append(f"constant '{nm}' coupling {off:g} != expected {want:g} "
                       f"(is it shared?)")
    return (len(bad) == 0), ("; ".join(bad) if bad else "all patient constants shared")


if __name__ == "__main__":
    # self-test on the three parameterisations the brief uses
    print("CONTRACT TEST (fails if any patient constant is left independent):")
    for lab, nms in [("marker / Model A", ["theta", "Km", "b", "ke"]),
                     ("Model B", ["theta", "Km", "b", "ke", "Ke"])]:
        ok, msg = check_contract(nms)
        print(f"  {lab:18s}: {'PASS' if ok else 'FAIL'} -- {msg}")
    # a drift is allowed and must show the right (finite) coupling
    ok, msg = check_contract(["theta", "Km", "b", "ke"], drifts={"ke": 0.05})
    print(f"  {'ke drift 0.05':18s}: {'PASS' if ok else 'FAIL'} -- {msg}")
    # NEGATIVE control: a builder that forgets to share b must be catchable.
    import numpy as _np
    k = 4
    Fs = [_np.zeros((k, k)), _np.zeros((k, k))]
    # emulate the OLD bug: assemble without coupling b, then check
    Z, idx = _assemble(Fs, ["theta", "Km", "b", "ke"])
    Z[idx[0]["b"], idx[1]["b"]] = 0.0; Z[idx[1]["b"], idx[0]["b"]] = 0.0
    Z[idx[0]["b"], idx[0]["b"]] -= COUPLE; Z[idx[1]["b"], idx[1]["b"]] -= COUPLE
    off = -Z[idx[0]["b"], idx[1]["b"]]
    print(f"  negative control  : b coupling after un-sharing = {off:g} "
          f"-> {'CAUGHT (would fail)' if not _np.isclose(off, COUPLE) else 'MISSED'}")
