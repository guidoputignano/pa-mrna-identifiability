"""
Round 19: can the on-therapy record settle WHICH SIDE of the 21-day crossing a
patient is on?

R-O and R-P said the record "reads lambda at 7-18 sigma" and therefore settles the
two sides of the 3.5452 d crossing "a fortiori". The 7-18 sigma is the separation
between a 3.6 d and a 7 d half-life, ln(7/3.6) = 0.665. The distance from the
measured 3.6 d to the crossing is ln(3.6/3.5452) = 0.015, about 43 times smaller.
So placing a value near 3.6 d on one side of the crossing is much HARDER than
telling 3.6 d from 7 d, not easier. This script prices how much harder.

Input: the per-patient sd(ln lambda) table that t46_timing.py prints (FIM, by
record length and cadence). It is parsed from out/t46_timing.txt, not retyped, so
the numbers here move if t46 is rerun.

Outputs, per patient:
  * z = ln(t_true / 3.5452) / sd(ln lambda) at the measured 3.6 d and at the
    round-17 central estimate 4.3 d
  * the 2-sigma "dead zone": true half-lives a single patient's record cannot
    place on either side of the crossing
  * the first record length at which a 4.3 d truth is placed at 2 sigma
  * how many patients a POOLED estimate would need to place a 3.6 d truth at
    2 sigma, if (and only if) lambda were the same in every patient

z is the EXPECTED separation (the log distance over the sd). One noisy record
lands at z + N(0, 1), so the chance that it actually clears 2 sigma on the right
side is Phi(z - 2): about 50 % when z = 2. The power columns report that chance,
and the pooled counts are given for 50 % and 80 % power.
"""
import re
from pathlib import Path

import numpy as np
from scipy.stats import norm

CROSS = 3.5452          # 21-day crossing (R-N, t48), Model C intrinsic theta_crit
MEASURED = 3.6          # Chandler & Ballard 1985, human fibroblast PCC
CENTRAL = 4.3           # round-17 compartment-anchored liver estimate (R-P, t50)
MOUSE = 7.0             # Attarwala 2023 effective, the comparator of R-N's 7-18 sigma
Z = 2.0

HERE = Path(__file__).resolve().parent


def t46_table():
    """
    {cadence: {months: sd(ln lambda)}} from t46. Full-precision values come from
    out/t46_timing_numbers.json (which t46 writes); the printed 3-decimal table in
    out/t46_timing.txt is the fallback.
    """
    js = HERE / "out" / "t46_timing_numbers.json"
    if js.exists():
        import json
        lam = json.loads(js.read_text())["lambda"]
        out = {"monthly": {}, "fortnightly": {}, "weekly": {}}
        for key, v in lam.items():
            cad, mo = key.split("/")
            out[cad][int(mo)] = float(v)
        return out
    txt = (HERE / "out" / "t46_timing.txt").read_text()
    block = txt.split("ITEM 2.")[1].split("spec reached by")[0]
    rows = re.findall(r"^\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$", block, re.M)
    if not rows:
        raise RuntimeError("could not parse the sd(ln lambda) table from out/t46_timing.txt")
    out = {"monthly": {}, "fortnightly": {}, "weekly": {}}
    for mo, a, b, c in rows:
        out["monthly"][int(mo)] = float(a)
        out["fortnightly"][int(mo)] = float(b)
        out["weekly"][int(mo)] = float(c)
    return out


CADENCE_DAYS = {"monthly": 30.4, "fortnightly": 14.0, "weekly": 7.0}


def sd_lnlam_at(thalf, months, cad_days, seeds=range(16)):
    """
    t46.sd_lnlam with the TRUE half-life as an argument. t46 evaluates the FIM at
    the repository's lambda = ln2/3.5 d; a patient whose true half-life is 3.6 d
    or 4.3 d has a slightly different record (the increment accumulates more) and
    so a slightly different sd. Same seeds, schedule jitter, missed doses and gap.
    """
    import core as C
    import t45_on_therapy as T45
    lam = np.log(2.0) / thalf
    T_days = months * 30.4
    n = max(4, int(round(T_days / cad_days)))
    tm = np.linspace(0.0, T_days, n)
    sk = dict(jitter=2.0, miss_p=0.08, years=T_days / 365.0)
    if T_days > 120:
        sk["gap"] = (0.4 * T_days, 28.0)
    out = []
    for sd_seed in seeds:
        dts, das = T45.schedule(seed=sd_seed, **sk)
        L, m_bar = C.N4.anchor(C.TH_5, T45.KM, T45.KE, "A")
        rng = np.random.default_rng(sd_seed)
        s2 = np.log1p(0.10 ** 2)
        loads = L * np.exp(rng.normal(-0.5 * s2, np.sqrt(s2), tm.size))
        inc, dinc = T45.inc_dinc(tm, dts, das, lam=lam)
        m = C.level(C.TH_5, T45.KM, T45.KE, loads, inc)
        r = dict(theta=C.TH_5, Km=T45.KM, ke=T45.KE, L=L, m_bar=m_bar, loads=loads,
                 inc=inc, dinc=dinc, m=m, times=tm)
        x = T45.sd_of(T45.fim5(r, with_lambda=True), 4)
        if np.isfinite(x):
            out.append(x)
    return float(np.median(out))


def truth_tables(tab):
    """sd(ln lambda) re-evaluated at each true half-life, plus a reproduction check."""
    check = max(abs(sd_lnlam_at(3.5, mo, CADENCE_DAYS[c]) - tab[c][mo]) for c in tab for mo in tab[c])
    at = {t: {c: {mo: sd_lnlam_at(t, mo, CADENCE_DAYS[c]) for mo in tab[c]} for c in tab}
          for t in (CROSS, MEASURED, CENTRAL)}
    return check, at


def main():
    tab = t46_table()
    check, at = truth_tables(tab)
    g_meas = np.log(MEASURED / CROSS)
    g_cent = np.log(CENTRAL / CROSS)
    g_mouse = np.log(MOUSE / MEASURED)
    print("=" * 92)
    print("  ROUND 19. CAN THE RECORD SETTLE WHICH SIDE OF THE 3.5452 d CROSSING? (per patient)")
    print("=" * 92)
    print(f"""
  log-distances:  3.6 d vs 7 d (R-N's 7-18 sigma)      ln = {g_mouse:.3f}
                  3.6 d vs the 3.5452 d crossing       ln = {g_meas:.4f}   ({g_mouse / g_meas:.0f}x smaller)
                  4.3 d vs the 3.5452 d crossing       ln = {g_cent:.3f}
  sd(ln lambda) per patient: t46's FIM (sd_lnlam), re-evaluated here at each TRUE half-life
  (t46 evaluates it at the repository's 3.5 d; reproduced exactly, max diff {check:.1e}).
  Dead zones use the sd at the crossing, z at 3.6 d the sd at 3.6 d, z at 4.3 d the sd at 4.3 d.
""")
    print(f"  {'cadence':>12} {'months':>6} {'sd @3.6d':>8} {'z @3.6d':>8} {'sd @4.3d':>8} {'z @4.3d':>8}"
          f" {'P(>=2s) @4.3d':>13} {'z 3.6 vs 7':>10}   2-sigma dead zone (true half-life, d)")
    for cad in ("monthly", "fortnightly", "weekly"):
        for mo in sorted(tab[cad]):
            s36, s43, sx = at[MEASURED][cad][mo], at[CENTRAL][cad][mo], at[CROSS][cad][mo]
            lo, hi = CROSS * np.exp(-Z * sx), CROSS * np.exp(Z * sx)
            print(f"  {cad:>12} {mo:>6} {s36:>8.3f} {g_meas / s36:>8.2f} {s43:>8.3f} {g_cent / s43:>8.2f}"
                  f" {norm.cdf(g_cent / s43 - Z):>13.0%} {g_mouse / s36:>10.1f}   {lo:.2f} - {hi:.2f}")
    best = min(min(v.values()) for v in at[MEASURED].values())
    print(f"""
  Best per-patient precision at a 3.6 d truth: sd(ln lambda) = {best:.4f}  ->  z = {g_meas / best:.2f}.
  A single patient's record never places a 3.6 d truth on a side of the crossing: even
  after 5 years of weekly draws it sits at {g_meas / best:.1f} sigma, inside the dead zone.
""")
    print("  A 4.3 d truth, one patient: first record length with expected z >= 2 (a coin flip that")
    print("  one record clears 2 sigma), and with an 80 % chance that it does (z >= 2.84):")
    for cad in ("monthly", "fortnightly", "weekly"):
        h50 = [mo for mo, s in sorted(at[CENTRAL][cad].items()) if g_cent / s >= Z]
        h80 = [mo for mo, s in sorted(at[CENTRAL][cad].items()) if norm.cdf(g_cent / s - Z) >= 0.80]
        f = lambda h: f"{h[0]} months" if h else "not within 60 months"
        print(f"    {cad:>12}: expected 2 sigma at {f(h50)};  80 % chance at {f(h80)}")
    print("\n  Pooling (ONLY if lambda is common to every patient, which nothing measured shows):")
    z80 = Z + norm.ppf(0.80)
    for cad, mo in (("monthly", 12), ("weekly", 12), ("weekly", 24)):
        s = at[MEASURED][cad][mo]
        n50 = int(np.ceil((Z * s / g_meas) ** 2))
        n80 = int(np.ceil((z80 * s / g_meas) ** 2))
        print(f"    {cad:>8} draws, {mo:>2} months: {n50:>4} patients (50 % power), {n80:>4} (80 % power)"
              f" to place a 3.6 d truth at 2 sigma")
    print(f"""
  READING. R-N's 7-18 sigma is right for its question (3.6 d or 7 d?). It does not carry
  to the crossing: the crossing is {g_mouse / g_meas:.0f}x closer (in log distance), so "a fortiori"
  runs the wrong way. Per patient the record settles the side only when the true liver
  half-life is clear of the crossing: near the 4.3 d central estimate, an 80 % chance
  within a year of monthly draws or six months of weekly ones; near the measured 3.6 d,
  never. A pooled trial estimate could do it with 23-77 patients at 50 % power or 46-155
  at 80 % (more at monthly draws), but only if lambda does not vary between patients,
  which is unmeasured. So the two emails still matter, and the on-therapy record is still
  the best route, but it settles the knife-edge only if the answer is not on the knife-edge.
""")


if __name__ == "__main__":
    main()
