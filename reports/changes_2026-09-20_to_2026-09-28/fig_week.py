"""
One plain figure for the week summary: where blood draws land in the dose
cycle, and why that decides what the record can say.

Left panel: draws booked on dosing days. Every draw sees the same point of the
cycle, so the record contains one operating point repeated.
Right panel: draws booked a little off the dose interval. Over the study they
work their way round the cycle and see the whole swing.

Both panels use the repository's own dose train and constants through
model_abc. No axis carries a symbol the summary does not explain.

Run:  python reports/changes_2026-09-20_to_2026-09-28/fig_week.py
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import model_abc as M                                        # noqa: E402
from common.style import TEXTWIDTH_IN, INK, MUTED, GREEN, RED, GREY   # noqa: E402
import matplotlib.pyplot as plt                              # noqa: E402

OUT = HERE / "fig"


def week_figure():
    fig, axes = plt.subplots(1, 2, figsize=(TEXTWIDTH_IN * 0.86, TEXTWIDTH_IN * 0.30),
                             sharey=True)
    phase = np.linspace(0, M.T_DOSE, 400)
    curve = M.train_delta(phase) / M.DELTA_PEAK        # one cycle, as a fraction of the peak

    #  dose-day draws are read just before the next dose, as a clinic would
    panels = [(M.T_DOSE, M.T_DOSE - 0.02, "booked on dosing days", RED),
              (M.T_DOSE + M.lap_step(), 0.0, "booked half a day later each time", GREEN)]
    for ax, (step, start, title, colour) in zip(axes, panels):
        ax.plot(phase, curve, color=MUTED, lw=1.3, zorder=1)
        t = start + M.schedule(step, n=M.VISITS_DRAWN)
        ph = np.mod(t, M.T_DOSE)
        ax.plot(ph, M.train_delta(t) / M.DELTA_PEAK, "o", color=colour, ms=5,
                alpha=0.85, zorder=3)
        ax.set_title(title, fontsize=8.5, color=INK)
        ax.set_xlim(-0.4, M.T_DOSE + 0.4)
        ax.set_xticks([0, M.T_DOSE / 2, M.T_DOSE])
        ax.set_xticklabels(["dose", "a week later", "next dose"], fontsize=7.5)
        ax.set_ylim(-0.04, 1.08)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.grid(axis="y", color=GREY, alpha=0.25, lw=0.5)
    axes[0].set_ylabel("added enzyme,\nshare of its peak", fontsize=8)
    axes[0].set_yticks([0, 0.5, 1.0])
    axes[0].set_yticklabels(["none", "half", "peak"], fontsize=7.5)
    axes[0].annotate("all 25 draws\nland here", xy=(M.T_DOSE - 0.3, 0.06), xytext=(6.6, 0.74),
                     fontsize=7.5, color=RED, ha="center",
                     arrowprops=dict(arrowstyle="->", color=RED, lw=0.8))
    axes[1].annotate("the 25 draws cover\nthe whole cycle", xy=(7.5, 0.28), xytext=(6.0, 0.70),
                     fontsize=7.5, color=GREEN, ha="center",
                     arrowprops=dict(arrowstyle="->", color=GREEN, lw=0.8))
    fig.tight_layout(w_pad=1.4)
    fig.savefig(OUT / "week_cycle.pdf")
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    week_figure()
    print("week_cycle.pdf")
