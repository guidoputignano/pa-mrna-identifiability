"""
Round-3 final. Three questions, then the study closes (the rest is empirical).

(1) DRUG RESPONSE: the drug is a KNOWN increment on the enzyme channel. What
    does a dose change do, can the contrast measure a drug response, and what is
    the protocol?
(2) ASSAY CV x DRAWS: can a less-invasive, noisier marker substitute more draws?
    Contour the contrast at 0.20.
(3) A FALSIFIABLE PREDICTION: the model's predicted OBSERVED within-patient
    marker variability, a number real serial-monitoring data can check.
"""
import numpy as np

import core as C
from n4_theta_scope import D0, delta_at, level


def data_fim(r, ac):
    return C.fim(r, noise="perdraw", ke_prior=0.0, b_prior=np.inf, km_prior=None, assay_cv=ac)


def contrast(n, s1, s2, ac=0.05, D1=D0, D2=D0, drift_ke=0.0, drift_b=0.0):
    r1 = C.record(C.TH_5, D=D1, seed=s1, times=np.linspace(0, 730, n))
    r2 = C.record(C.TH_5, D=D2, seed=s2, times=np.linspace(0, 730, n))
    L = r1["L"]
    F1, F2 = data_fim(r1, ac), data_fim(r2, ac)
    Z = np.zeros((7, 7)); i1, i2 = [0, 2, 3, 5], [1, 2, 4, 6]
    for a in range(4):
        for b in range(4):
            Z[i1[a], i1[b]] += F1[a, b]; Z[i2[a], i2[b]] += F2[a, b]

    def cp(i, j, sd):
        g = 1e12 if sd <= 0 else 1.0 / sd ** 2
        Z[i, i] += g; Z[j, j] += g; Z[i, j] -= g; Z[j, i] -= g
    cp(3, 4, drift_b * L); cp(5, 6, drift_ke)
    Z[5, 5] += 0.8 ** -2; Z[6, 6] += 0.8 ** -2
    Z[3, 3] += (0.05 * L) ** -2; Z[4, 4] += (0.05 * L) ** -2
    P = np.linalg.inv(Z); e = np.array([-1, 1, 0, 0, 0, 0, 0.])
    return float(np.sqrt(e @ P @ e))


def med(n, **kw):
    return float(np.median([contrast(n, 2 * s, 2 * s + 1, **kw) for s in range(6)]))


def sdlnm(n, load_cv, drug, assay, seeds=24):
    out = []
    for s in range(seeds):
        r = C.record(C.TH_5, D=(D0 if drug else 0.0), seed=s,
                     times=np.linspace(0, 730, n), load_cv=load_cv)
        m = r["m"]
        if assay:
            m = m * np.exp(np.random.default_rng(700 + s).normal(0, 0.05, m.size))
        out.append(np.std(np.log(m)))
    return float(np.median(out))


def main():
    print("=" * 92)
    print("  FINAL. DRUG RESPONSE, ASSAY-CV x DRAWS, AND A FALSIFIABLE PREDICTION")
    print("=" * 92)

    print("\n  (1) DRUG RESPONSE. The drug is a KNOWN per-draw increment Delta on the")
    print("      enzyme channel; at the booked dose theta+Delta ~ 14x theta (it dominates).")
    L, m0 = C.N4.anchor(C.TH_5, C.KM0, C.KE, "A")
    for sc, lab in ((0.0, "no drug"), (0.5, "half"), (1.0, "booked"), (2.0, "double")):
        inc = delta_at(C.MONTHLY, D0 * sc)
        m = level(C.TH_5, C.KM0, C.KE, L, inc)
        eff = (C.TH_5 + np.median(inc)) / C.TH_5
        print(f"        {lab:>6} dose: effective capacity theta+Delta = {eff:4.1f}x theta, "
              f"median marker {np.median(m):5.2f}")
    print("      A RECORDED dose change does NOT fake a theta change (Delta is modeled) --")
    print(f"      the contrast stays unbiased, at a precision cost: same dose "
          f"{med(104):.3f}, double dose on window 2 {med(104, D2=2*D0):.3f}. This is the")
    print("      GOOD DUAL of the antibiotic (which is unrecorded and biases, R-E).")
    print("      So the contrast is SPECIFIC to the patient's own capacity theta. Its")
    print("      drug-response use is detecting a change in theta ITSELF (a disease-")
    print("      modifying response), which it resolves at +45%/-30% (weekly) over a")
    print("      weeks window -- matching the contrast's validity. A dose change alone")
    print("      moves only the known Delta, so it needs no measurement; and if the")
    print("      delivered dose/kappa were UNCERTAIN it would leak into theta like the")
    print("      antibiotic, so dosing must be documented.")
    print("      PROTOCOL (drug-response monitoring): two windows of ~45-100 draws each")
    print("      (a ~6-week daily burst, or weekly for longer), SAME recorded dose and")
    print("      diet, antibiotic status fixed, gap = the expected response interval")
    print("      (weeks); read a >=45% change in theta at 80% power.")

    print("\n  (2) ASSAY CV x DRAWS PER WINDOW: contrast sd(d ln theta).")
    ns = [25, 50, 100, 200, 400, 800]
    print("      CV \\ draws " + "".join(f"{n:>7}" for n in ns))
    cvs = [0.05, 0.10, 0.15, 0.20, 0.30]
    tab = {}
    for cv in cvs:
        row = [med(n, ac=cv) for n in ns]
        tab[cv] = row
        print(f"      {cv:>7.0%}    " + "".join(f"{v:7.3f}" for v in row))
    print("      draws per window needed to reach 0.20, by assay CV (interp, ~ n prop CV^2):")
    for cv in cvs:
        xs, ys = np.log(ns), np.log(tab[cv])
        nc = np.exp(np.interp(np.log(0.20), ys[::-1], xs[::-1]))
        print(f"        CV {cv:>3.0%}: n ~ {nc:5.0f} draws/window")
    print("      VERDICT: variance ~ CV^2/n, so a less-invasive but noisier marker fails")
    print("      to substitute -- a 20% CV collection needs ~340 draws/window against ~45")
    print("      at 5%, and 30% CV needs ~760. Low-invasiveness at high CV is not")
    print("      deliverable; the protocol needs a ~5% assay.")

    print("\n  (3) A FALSIFIABLE PREDICTION: observed within-patient marker variability.")
    print("      The contrast's 0.13 rests on the shared nuisances cancelling, a property")
    print("      of the simulation. The CHECKABLE consequence is how much the marker")
    print("      ITSELF fluctuates within a patient on serial sampling. The model predicts")
    print("      serial within-patient sd(ln marker):")
    print(f"        dose-phase RANDOM sampling (full)        : {sdlnm(25,0.10,True,True):.3f}"
          f"  (CV ~{100*sdlnm(25,0.10,True,True):.0f}%) -- dominated by the 14-day dose cycle")
    print(f"        dose-SYNCHRONISED sampling (load+assay)  : {sdlnm(25,0.10,False,True):.3f}"
          f"  (CV ~{100*sdlnm(25,0.10,False,True):.0f}%)")
    print(f"          of which load (10%) alone {sdlnm(25,0.10,False,False):.3f}, "
          f"assay (5%) alone {sdlnm(25,0.0,False,True):.3f}")
    print("      THE TEST: on serial monitoring data sampled at a consistent time relative")
    print("      to dosing, the model predicts within-patient marker CV ~ 12%. If real")
    print("      patients fluctuate materially more than that (after removing dose-phase),")
    print("      the extra variance is either real theta movement (the 'stable theta'")
    print("      assumption fails) or uncancelled nuisance -- either way the cancellation")
    print("      is optimistic and the 0.13 contrast, the main result, is too. This is the")
    print("      one number that converts the study from simulated to testable.")


if __name__ == "__main__":
    main()
