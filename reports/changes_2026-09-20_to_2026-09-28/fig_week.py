"""
The two figures for the weekly update.

week_cycle.pdf   where blood draws land in the dose cycle, booked two ways.
week_band.pdf    the spread of the estimate against how much enzyme the patient
                 has left, with the range it works over and the range this
                 disease sits in.

Both come from the repository's own code: the dose train and constants through
model_abc, and the capacity sweep through the module that generates the band
the paper quotes. No axis carries a symbol the update does not explain.

Run:  python reports/changes_2026-09-20_to_2026-09-28/fig_week.py
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "audit"))
import model_abc as M                                              # noqa: E402
from common.style import INK, MUTED, GREEN, RED, BLUE, GREY        # noqa: E402
import matplotlib.pyplot as plt                                    # noqa: E402

OUT = HERE / "fig"
PANEL = 3.3          # inches, so the two figures sit side by side on one page


def cycle_figure():
    """The same number of draws, booked two ways, over one dose cycle."""
    fig, axes = plt.subplots(2, 1, figsize=(PANEL, PANEL * 0.62), sharex=True)
    phase = np.linspace(0, M.T_DOSE, 400)
    curve = M.train_delta(phase) / M.DELTA_PEAK

    #  a dose-day draw is read just before the next infusion, as a clinic would
    panels = [(M.T_DOSE, M.T_DOSE - 0.02, "booked on dosing days", RED,
               "all %d draws land here" % M.VISITS_DRAWN),
              (M.T_DOSE + M.lap_step(), 0.0, "booked half a day later each time",
               GREEN, "the draws cover the cycle")]
    for ax, (step, start, title, colour, note) in zip(axes, panels):
        ax.plot(phase, curve, color=MUTED, lw=1.2, zorder=1)
        t = start + M.schedule(step, n=M.VISITS_DRAWN)
        ax.plot(np.mod(t, M.T_DOSE), M.train_delta(t) / M.DELTA_PEAK, "o",
                color=colour, ms=4, alpha=0.85, zorder=3)
        ax.set_title(title, fontsize=7.5, color=INK, pad=3)
        ax.text(0.97, 0.74, note, transform=ax.transAxes, ha="right",
                fontsize=6.8, color=colour)
        ax.set_xlim(-0.4, M.T_DOSE + 0.4)
        ax.set_ylim(-0.06, 1.14)
        ax.set_yticks([0, 1.0])
        ax.set_yticklabels(["none", "peak"], fontsize=6.8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.grid(axis="y", color=GREY, alpha=0.25, lw=0.5)
    axes[1].set_xticks([0, M.T_DOSE / 2, M.T_DOSE])
    axes[1].set_xticklabels(["dose", "a week later", "next dose"], fontsize=6.8)
    axes[0].set_ylabel("added enzyme", fontsize=7.2)
    axes[1].set_ylabel("added enzyme", fontsize=7.2)
    fig.tight_layout(h_pad=0.8)
    fig.savefig(OUT / "week_cycle.pdf")
    plt.close(fig)


def band_figure():
    """The estimate's spread against capacity, across the whole range."""
    import n4_theta_scope as N4
    pct = np.geomspace(1.0, 1000.0, 121)
    med, lo, hi = np.array([N4.over_seeds(theta=N4.TH0 * p / 100.0) for p in pct]).T
    #  read the band off the committed generated file rather than calling the
    #  generator, which would rewrite it as a side effect
    import re
    gen = (ROOT / "results/numbers/scope.tex").read_text()
    band = {k: float(v) for k, v in
            re.findall(r"paper@num@(band_lo_pct|band_hi_pct)\\endcsname\{([\d.]+)\}", gen)}
    lo_pct, hi_pct = band["band_lo_pct"], band["band_hi_pct"]
    dis_pct = 100.0 * M.THETA_DISEASE / M.THETA

    fig, ax = plt.subplots(figsize=(PANEL, PANEL * 0.62))
    ax.axvspan(lo_pct, hi_pct, color="#eef2f7", zorder=0)
    ax.axvspan(1, dis_pct, color="#fbe9e2", zorder=0)
    ax.fill_between(pct, lo, hi, color=BLUE, alpha=0.12, lw=0)
    ax.plot(pct, med, color=BLUE, lw=1.5)
    ax.axhline(0.20, color=INK, lw=0.7)
    ax.text(12.0, 0.215, "what we would need", color=INK, fontsize=6.4, va="bottom")
    ax.text(np.sqrt(lo_pct * hi_pct), 15.5, "the estimate works here",
            ha="center", va="top", fontsize=6.8, color=INK)
    ax.text(2.1, 0.30, "these\npatients", ha="center", va="bottom",
            fontsize=6.8, color="#8a3b22")
    for p_mark, va in ((dis_pct, "bottom"), (100.0, "top")):
        y = float(N4.over_seeds(theta=N4.TH0 * p_mark / 100.0)[0])
        ax.scatter([p_mark], [y], s=16, color=BLUE, edgecolor="white", lw=0.7, zorder=4)
        ax.text(p_mark * (1.25 if va == "bottom" else 1.35), y * (1.15 if va == "bottom" else 0.72), f"{y:.2f}",
                fontsize=6.6, color=INK, va=va)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1, 1000)
    ax.set_ylim(0.08, 20)
    ax.set_xticks([1, 5, 10, 100, 1000])
    ax.set_xticklabels(["1%", "5%", "10%", "100%", "1000%"], fontsize=6.8)
    ax.set_yticks([0.1, 0.2, 1, 10])
    ax.set_yticklabels(["0.1", "0.2", "1", "10"], fontsize=6.8)
    ax.grid(color=GREY, alpha=0.3, lw=0.5)
    ax.set_axisbelow(True)
    ax.set_xlabel("enzyme the patient still has, % of normal", fontsize=7.2)
    ax.set_ylabel("spread of the estimate,\nsmaller is better", fontsize=7.2)
    fig.tight_layout()
    fig.savefig(OUT / "week_band.pdf")
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    cycle_figure()
    band_figure()
    print("week_cycle.pdf, week_band.pdf")
