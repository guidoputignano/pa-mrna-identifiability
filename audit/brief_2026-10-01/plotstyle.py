"""Figure style for this directory: the repository's rcParams, a validated palette."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import style as _repo_style  # noqa: F401  (applies the repo's rcParams)

#: categorical order, validated with the dataviz skill's validate_palette.js
#: (light surface): all checks pass, contrast WARN on slots 3 and 4, so every
#: series is also direct-labelled or named in a legend and tabulated.
C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, MUTED, GRID = "#1a1a1a", "#5a5a5a", "#d9d8d3"
FIG = Path(__file__).resolve().parent / "fig"
FIG.mkdir(exist_ok=True)


def finish(ax):
    ax.grid(True, color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote fig/{name}.png")


def plain(ax):
    """Plain-number tick labels on log axes, instead of 6x10^0."""
    from matplotlib.ticker import FuncFormatter, LogLocator
    f = FuncFormatter(lambda v, _: f"{v:g}")
    for axis in (ax.xaxis, ax.yaxis):
        if axis.get_scale() == "log":
            axis.set_major_formatter(f)
            axis.set_minor_formatter(FuncFormatter(
                lambda v, _: f"{v:g}" if f"{v:.0e}"[0] in "25" else ""))
