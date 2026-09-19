"""
Phase 3, the design output: how many visits, and what precision they buy.

For each schedule, k_e prior and load variability, sd(ln theta_hat) against the
number of visits, with the 20% target marked, and a fit of

    sd(n) = sqrt(floor^2 + c^2/n)

whose two constants are the design answer: `c` is what visits buy and `floor`
is what they cannot. A design is worth extending only while c/sqrt(n) is still
above the floor.

Footing, stated once and true of every number below. Information filter in the
paper's coordinates (ln theta, ln K_m, b, ln k_e), load record 5%, sigma =
sigma_total(0.05, 0.05) = the repository's own. `k_e known` DROPS the column
rather than giving it a flat prior, which is _filter's convention and e9's. n
visits at the stated interval run from day 0, so the horizon grows with n.
Every load-varying figure is the median over 8 seeds, with the range reported.

Writes results/figures/fig_E2_visits.{pdf,png}.

Nothing here records an expected result in advance. The one number quoted from
outside is a second reader's estimate of the monthly floor and c at prior 0.8,
labelled where it is used.

Run:  python audit/g2_visits_vs_sd.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import warnings

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, OptimizeWarning

from common.style import save, figsize, INK, MUTED, RED, GREEN, BLUE, AMBER, GREY
import g_design_common as G

# a second reader's estimate for the monthly design at k_e prior 0.8, quoted
# so the fit below can be compared with it
REVIEWER_FLOOR, REVIEWER_C = 0.145, 0.45

PRIORS = [(0.8, "$k_e$ prior 0.8"), (0.3, "$k_e$ prior 0.3"),
          (0.05, "$k_e$ prior 0.05"), (None, "$k_e$ known")]
CVS = [(0.0, "load CV 0%"), (0.02, "load CV 2%"), (0.10, "load CV 10%")]
COLOUR = {"14 d (the dose interval)": RED,
          "14 d, already on therapy": MUTED, "14.5 d": AMBER,
          "monthly (30.4 d)": BLUE, "phase-spread": GREEN}
NS = np.array([4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20,
               25, 30, 36, 44, 52, 60])   # integer resolution to n=14,
                                          # where the target is crossed


def curve(schedule_fn, load_cv, ke_prior, ns=NS, seeds=range(8)):
    """Median sd(ln theta_hat) at each n, with the seed range."""
    lo, md, hi = [], [], []
    for n in ns:
        a, b, c = G.sd_over_seeds(schedule_fn(int(n)), load_cv,
                                  ke_prior=ke_prior, seeds=seeds)
        lo.append(a); md.append(b); hi.append(c)
    return np.array(lo), np.array(md), np.array(hi)


def fit_floor_c(ns, sd):
    """
    Fit sd(n) = sqrt(floor^2 + c^2/n).

    Fitted on the finite points only. Returns (floor, c, rms residual, n used)
    and nan if fewer than three points survive.
    """
    ok = np.isfinite(sd)
    if ok.sum() < 3:
        return np.nan, np.nan, np.nan, int(ok.sum())
    x, y = ns[ok].astype(float), sd[ok]
    f = lambda n, fl, c: np.sqrt(fl ** 2 + c ** 2 / n)
    # seed from the two-point solve: floor^2 = the large-n value
    p0 = [max(y.min() * 0.9, 1e-4), max(y[0] * np.sqrt(x[0]), 1e-3)]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            p, _ = curve_fit(f, x, y, p0=p0, maxfev=20000)
    except Exception:
        return np.nan, np.nan, np.nan, int(ok.sum())
    r = f(x, *p) - y
    return abs(float(p[0])), abs(float(p[1])), float(np.sqrt(np.mean(r ** 2))), int(ok.sum())


def figure():
    fig, axes = plt.subplots(len(CVS), len(PRIORS),
                             figsize=figsize(1.0, 0.80),
                             sharex=True, sharey=True)
    fig.set_layout_engine("none")
    fig.subplots_adjust(left=0.085, right=0.985, top=0.91, bottom=0.185,
                        wspace=0.10, hspace=0.18)

    for i, (cv, cvlab) in enumerate(CVS):
        for j, (kp, kplab) in enumerate(PRIORS):
            ax = axes[i, j]
            ax.axhline(G.GATE, color=GREY, ls=":", lw=1.0, zorder=1)
            for name, fn in G.SCHEDULES.items():
                lo, md, hi = curve(fn, cv, kp)
                good = np.isfinite(md)
                if good.any():
                    ax.plot(NS[good], md[good], color=COLOUR[name], lw=1.4,
                            zorder=3)
                    if cv:
                        ax.fill_between(NS[good], lo[good], hi[good],
                                        color=COLOUR[name], alpha=0.16, lw=0)
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_xlim(3.6, 66); ax.set_ylim(2e-2, 12)
            ax.set_xticks([4, 6, 10, 20, 40, 60])
            ax.set_xticklabels(["4", "6", "10", "20", "40", "60"])
            ax.set_yticks([0.03, 0.1, 0.3, 1.0, 3.0, 10.0])
            ax.set_yticklabels(["0.03", "0.1", "0.3", "1", "3", "10"])
            ax.minorticks_off()
            ax.grid(True, which="major", alpha=0.45)
            if i == 0:
                ax.set_title(kplab, fontsize=8)
            if j == 0:
                ax.set_ylabel(cvlab + "\n" + r"sd$(\ln\hat\theta)$", fontsize=7.5)
            if i == len(CVS) - 1:
                ax.set_xlabel("number of visits", fontsize=8)
            ax.tick_params(labelsize=6.5)

    axes[0, -1].text(0.98, G.GATE * 1.18, "0.20 target", ha="right",
                     transform=axes[0, -1].get_yaxis_transform(),
                     fontsize=6.5, color=GREY)
    handles = [plt.Line2D([], [], color=COLOUR[k], lw=1.6) for k in G.SCHEDULES]
    fig.legend(handles, list(G.SCHEDULES), loc="lower center", ncol=5,
               fontsize=7.0, bbox_to_anchor=(0.5, 0.012))
    #  the prior-conditional caveat belongs in the caption, where a
    #  reader meets it before the curves rather than under them.
    save(fig, "fig_E2_visits")


def saturation_n(fn, ns=NS, frac=0.95):
    """
    The first n at which the schedule's Delta/theta range reaches `frac` of
    its own maximum over the sweep.

    Below it the fit is describing two effects at once, more visits AND more
    increment range, and the one-over-root-n form cannot hold.

    Where n* lands is measured, not assumed, and it is larger than the phase
    step alone suggests: a schedule has not "covered the cycle" until it has
    sampled near BOTH ends of the increment, and the ends are where the dose
    train is changing fastest. A schedule whose n* sits at the end of the
    sweep has not saturated within it at all, and for that one the all-n fit
    is the only one there is; the table says so rather than printing nan.
    """
    r = np.array([G.delta_range(fn(int(n))) for n in ns])
    ok = np.where(r >= frac * r.max())[0]
    return int(ns[ok[0]]) if ok.size else int(ns[-1])


def table():
    print("\n  sd(n) = sqrt(floor^2 + c^2/n). floor is what visits cannot buy,")
    print("  c is what they can. Load-varying rows are fitted to the MEDIAN")
    print("  over 8 seeds. Two fits are given, because the form only describes")
    print("  a schedule once that schedule has covered the dose cycle: before")
    print("  then sd falls for two reasons at once and the fit absorbs both.")
    print("  n* is the first n at which the Delta/theta range reaches 95% of")
    print("  its own maximum, and the SATURATED fit uses n >= n*.\n")
    print("  `n at 0.20` is given twice: as the sweep OBSERVES it, the smallest")
    print("  n tried whose median sd is already at or below the target, and as")
    print("  the fit IMPLIES it. A starred implied value comes from the")
    print("  saturated fit, an unstarred one from the all-n fit because the")
    print("  saturated one had too few points.\n")
    print(f"  {'schedule':>26} {'k_e pr':>7} {'CV':>5} {'n*':>4} | "
          f"{'floor':>8} {'c':>7} {'rms':>9} | "
          f"{'floor*':>8} {'c*':>7} {'rms*':>9} | {'obs':>6} {'implied':>8}")
    print("  " + "-" * 116)
    rows = []
    for name, fn in G.SCHEDULES.items():
        nstar = saturation_n(fn)
        for kp, _ in PRIORS:
            for cv, _ in CVS:
                _, md, _ = curve(fn, cv, kp)
                fl, c, rms, k = fit_floor_c(NS, md)
                sel = NS >= nstar
                fl2, c2, rms2, _ = fit_floor_c(NS[sel], md[sel])
                if np.isnan(fl):
                    print(f"  {name:>26} {str(kp):>7} {cv:4.0%} {nstar:4d} | "
                          f"{'singular at every n' if k == 0 else 'unfittable':>58}")
                    continue
                # the saturated fit where it exists, the all-n fit where it
                # does not. Reading nan as "never" would report a schedule
                # that reaches the target in ten visits as one that never does.
                star = np.isfinite(fl2) and np.isfinite(c2)
                ufl, uc, tag = (fl2, c2, "*") if star else (fl, c, " ")
                need = (f"{uc**2/(G.GATE**2 - ufl**2):.0f}{tag}"
                        if ufl < G.GATE else "never")
                hit = NS[np.isfinite(md) & (md <= G.GATE)]
                obs = f"{int(hit[0])}" if hit.size else f">{int(NS[-1])}"
                f2 = lambda v: "     --" if not np.isfinite(v) else f"{v:7.4f}"
                print(f"  {name:>26} {str(kp):>7} {cv:4.0%} {nstar:4d} | "
                      f"{fl:8.4f} {c:7.4f} {rms:9.2e} | "
                      f"{f2(fl2):>8} {f2(c2):>7} "
                      f"{'       --' if not star else f'{rms2:9.2e}'} | "
                      f"{obs:>6} {need:>8}")
                rows.append((name, kp, cv, ufl, uc, need, obs))
    print("\n  Read the starred columns for floor and c, and the OBSERVED column")
    print("  for how many visits to book. Where rms and rms* differ by an order")
    print("  of magnitude the unstarred fit is describing cycle coverage rather")
    print("  than sample size.")
    print("\n  Where observed and implied disagree, implied is the optimistic")
    print("  one and the gap is the price of poor early coverage: the saturated")
    print("  fit answers `how many visits if the first n had the coverage the")
    print("  large-n design has`, and for a fixed interval they do not. The one")
    print("  schedule where the two agree is phase-spread, which is what it is")
    print("  for -- its coverage is built in at every n rather than accumulated.")
    return rows


def against_the_reviewer():
    print(f"\n  The second reader's estimate for the monthly design at k_e "
          f"prior 0.8 is floor {REVIEWER_FLOOR}, c {REVIEWER_C}.")
    fn = G.SCHEDULES["monthly (30.4 d)"]
    nstar = saturation_n(fn)
    for cv, lab in CVS:
        _, md, _ = curve(fn, cv, 0.8)
        fl, c, rms, _ = fit_floor_c(NS, md)
        sel = NS >= nstar
        fl2, c2, rms2, _ = fit_floor_c(NS[sel], md[sel])
        print(f"    {lab:>14}: all n     floor {fl:.4f} ({fl/REVIEWER_FLOOR:5.2f}x), "
              f"c {c:.4f} ({c/REVIEWER_C:5.2f}x), rms {rms:.2e}")
        print(f"    {'':>14}  n >= {nstar:<3d}  floor {fl2:.4f} ({fl2/REVIEWER_FLOOR:5.2f}x), "
              f"c {c2:.4f} ({c2/REVIEWER_C:5.2f}x), rms {rms2:.2e}")
    print("\n  The estimate is not attached to a load CV, so all three are shown.")


def caveat():
    print("\n  " + "-" * 76)
    print("  What the load-variation benefit assumes.")
    print("  Every gain from load CV above rests on the diet record capturing")
    print("  REAL day-to-day variation in intake. It is the recorded L_k that")
    print("  enters the design columns, so a record that is imputed from a")
    print("  prescription, or averaged to a weekly mean, carries no variation")
    print("  into the fit however much the patient's intake actually varies.")
    print("  On the numbers above that turns the load CV 10% rows into the")
    print("  load CV 0% rows. The paper already asks for this as a row of the")
    print("  specification (06_excitation.tex, after rem:rankvsprec); these")
    print("  figures price it.")


if __name__ == "__main__":
    print(__doc__)
    figure()
    table()
    against_the_reviewer()
    caveat()
