"""
One registry for every constant the note quotes, and where each one came from.

The rule this file enforces is that no number reaches the paper without a row
saying whether it was measured, derived from another number, or assumed. The
failure it is built against is specific and has happened three times in this
repository: ASSAY_CV was a percentage used as absolute flux, the load-error
amplification was a ratio quoted as a constant, and M_i was never written down
at all, so Figure 16 could not be reproduced from the paper. The first two were
caught by a reader; the third was caught only because someone tried to redraw
the figure.

Three source classes, and they mean different things:

    literature  a number read off a cited source, in that source's own units
    derived     fixed by a constraint elsewhere in the model, not free
    assumed     chosen. No measurement behind it. Its cost has to be bounded

An `assumed` row is not a defect on its own. An `assumed` row with no bound on
what it costs is.

`make figures` writes results/tables/tab_params.tex from this file, and
tests/test_models.py asserts that the values here are the values the model code
actually runs on, so the table cannot drift away from the simulation.
"""

from dataclasses import dataclass

from .models import ParamsA, ParamsC


LITERATURE, DERIVED, ASSUMED = "literature", "derived", "assumed"


@dataclass(frozen=True)
class Row:
    symbol: str          # as it appears in the paper, in LaTeX
    value: str           # rendered, with units where there are any
    meaning: str
    source: str          # one of the three classes above
    note: str            # what fixes it, or what its being assumed costs


_A = ParamsA()
_C = ParamsC()

# The operating point, the tolerable level and the illness fold live in
# model_A/run.py rather than in a dataclass, so they are named here against
# that module and checked against it in the tests.
M_OP, M_TOL, ILLNESS_FOLD = 10.0, 40.0, 1.50
ASSAY_CV = 5.0
S_BASE_C = 600.0
ALT_PEAK_C = 500.0


LOCAL = [
    Row(r"$\bar m$", f"{M_OP:.0f}", "operating point: the level management holds the patient at",
        ASSUMED, "illustrative. Every design quantity is reported relative to it, so its "
                 "absolute value cancels out of the conclusions"),
    Row(r"$\theta$", f"{_A.theta:.0f}", "residual capacity of the enzymatic channel",
        ASSUMED, "illustrative. It is the quantity the note says cannot be identified, so no "
                 "conclusion may depend on its value"),
    Row(r"$K_m$", f"{_A.Km:.0f}", "half-saturation of the enzymatic channel, at the plasma level",
        ASSUMED, r"illustrative, and \emph{not} the in vitro Michaelis constant of PCC "
                 r"($0.29$ mM for propionyl-CoA), which is a different molecule in a different "
                 r"compartment. The ratio $\bar m/K_m = 0.1$ is what the excitation condition "
                 r"is stated against"),
    Row(r"$k_e$", f"{_A.ke:.0f}", "alternative clearance rate constant",
        ASSUMED, "illustrative. Chosen so the enzymatic channel carries most of the flux, which "
                 "is the favourable case for identification"),
    Row(r"$\mtol$", f"{M_TOL:.0f}", "tolerable metabolite level",
        ASSUMED, r"a clinical convention, not a property of the patient. The width of the "
                 r"undetermined window runs $1.2\%$ to $67\%$ as it moves, so this is the most "
                 r"expensive assumption in the local model"),
    Row(r"$\spk/s_{\mathrm{base}}$", f"{ILLNESS_FOLD:.2f}", "illness amplitude, as a fold of the patient's own baseline load",
        ASSUMED, "stated as a fold rather than absolutely, because an absolute figure disagreed "
                 "between sections. Model C's own fold sits at 1.18"),
    Row(r"$\kappa$", "not assigned", "gain from recorded carnitine dose onto $k_e$",
        ASSUMED, r"never assigned a value. Only the product $\kappa u_c$ enters, as a lumped "
                 r"increment on $k_e$, and every figure varies that product directly. "
                 r"The ill-conditioning result of supplement \S S1 is conditional on it: an unknown $\kappa$ "
                 r"leaves the carnitine route two equations in four unknowns"),
    Row(r"$\lambda$", r"$\ln 2/3.5$ d$^{-1}$",
        "decay rate of the therapeutic increment",
        ASSUMED, r"a $3.5$-day half-life, so the increment has decayed sixteen-fold by the "
                 r"next dose and the inter-dose recovery is an excursion rather than a "
                 r"ripple. Nothing measures it. It was previously named and never assigned, "
                 r"which made Section~\ref{sec:collapse} unrunnable; its only role is to set "
                 r"how much excursion a fixed schedule creates, and "
                 r"Table~\ref{tab:excursion} reports the answer"),
    Row(r"$T_{\mathrm{dose}}$", "14 d", "interval between therapeutic doses",
        LITERATURE, r"the biweekly schedule selected for mRNA-3927 by model-informed dose "
                    r"selection~\cite{dose_selection_poster}. Adopted rather than assumed"),
    Row(r"$D$", "6000", "size of one therapeutic increment",
        ASSUMED, r"chosen so the marker falls by about half at peak effect, the order of the "
                 r"reductions reported in the phase 1/2 interim "
                 r"analysis~\cite{trial_phase12}. Only the excursion it produces is used"),
    Row(r"$\sigma$", "a family", "standard deviation of the ordinate of the fit",
        ASSUMED, r"deliberately not assigned. Every design conclusion is reported as a function "
                 r"of it, because nothing in this note bounds how well the substrate load is "
                 r"recorded. The assay contributes $b\,\sigma_m$ to the same total"),
    Row(r"$\mathrm{CV}$", f"{ASSAY_CV:.0f}\\%", "analytical variability of the marker assay",
        ASSUMED, r"illustrative, and a \emph{percentage}: it is never absolute flux. It enters "
                 r"the ordinate as $b\,\mathrm{CV}\,\bar m$"),
]

PARENT = [
    Row(r"$\theta$", f"{_C.theta:.0f}", "residual capacity, parent model",
        ASSUMED, r"illustrative, about $10\%$ of normal. Chosen to sit above "
                 r"$\thc^{\,\mathrm{int}}$ so the patient has a margin to measure"),
    Row(r"$K_m$", f"{_C.Km:.0f}", "half-saturation of the enzymatic channel",
        ASSUMED, "as above, and shared with the local model so the two can be compared"),
    Row(r"$K_c$", f"{_C.Kc:.0f}", "half-saturation of the conjugation step",
        ASSUMED, r"no measurement, and no reason a conjugation reaction should share an "
                 r"enzyme's half-saturation. Cost bounded in supplement \S S1"),
    Row(r"$M_i$", f"{_C.Mi:.0f}", "metabolite level at which acylcarnitine competition doubles the non-conjugative loss",
        ASSUMED, r"no measurement. Previously not written down anywhere, which is why "
                 r"the margin figure could not be reproduced from the paper. Cost bounded in "
                 r"supplement \S S1"),
    Row(r"$\rho_0$", f"{_C.rho0:.2f}", r"basal non-conjugative loss relative to conjugation, $\delta_0/k$",
        ASSUMED, r"no direct measurement. Bounded above by $1$, since the carnitine drain in "
                 r"this disease is dominated by propionylcarnitine, and below by the measured "
                 r"direction of free-carnitine excretion"),
    Row(r"$R_0$", f"{_C.R0:.1f}", "carnitine replenishment",
        DERIVED, r"\textbf{not free}. Fixed by requiring the alternative route alone to peak at "
                 rf"${ALT_PEAK_C:.0f}$ against a baseline load of ${S_BASE_C:.0f}$, so that a "
                 r"patient with $\theta = 0$ cannot be maintained. Changing $K_c$, $M_i$ or "
                 r"$\rho_0$ moves that peak and $R_0$ must be re-derived"),
    Row(r"$s_{\mathrm{base}}$", f"{S_BASE_C:.0f}", "baseline propiogenic load",
        ASSUMED, "illustrative. Sets the scale everything else in the parent model is quoted "
                 "against"),
    Row(r"$k$", "cancels", "conjugation rate constant",
        DERIVED, r"drops out of the quasi-steady reduction of supplement \S S1 exactly, so the "
                 r"reduced clearance does not depend on it and it is a free scale in the "
                 r"two-state system"),
]


def _tex(rows):
    out = []
    for r in rows:
        out.append(" & ".join([r.symbol, r.value, r.meaning,
                               r"\textsc{" + r.source + "}", r.note]) + r" \\")
    return out


def write_params_table():
    """Render both blocks into one table, source column included."""
    from .style import TABDIR, ROOT
    body = [r"\begin{tabular}{@{}p{0.07\textwidth} p{0.07\textwidth} "
            r"p{0.20\textwidth} p{0.11\textwidth} p{0.45\textwidth}@{}}",
            r"\toprule",
            r"\textbf{Symbol} & \textbf{Value} & \textbf{Meaning} & "
            r"\textbf{Source} & \textbf{What fixes it, or what assuming it costs} \\",
            r"\midrule",
            r"\multicolumn{5}{@{}l}{\emph{The local identification model}} \\",
            r"\addlinespace[2pt]"]
    body += _tex(LOCAL)
    body += [r"\addlinespace", r"\midrule",
             r"\multicolumn{5}{@{}l}{\emph{The parent model and its equilibrium reduction}} \\",
             r"\addlinespace[2pt]"]
    body += _tex(PARENT)
    body += [r"\bottomrule", r"\end{tabular}"]
    path = TABDIR / "tab_params.tex"
    path.write_text("\n".join(body))
    print(f"  wrote {path.relative_to(ROOT)}")
    return path


if __name__ == "__main__":
    write_params_table()
