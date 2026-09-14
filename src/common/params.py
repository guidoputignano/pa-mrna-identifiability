"""
One registry for every constant the paper quotes, and where each one came from.

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
        ASSUMED, r"illustrative, and every design quantity is reported relative to it, so "
                 r"its absolute value cancels out of the conclusions. On a real scale: the "
                 r"translational model puts the untreated plasma 3-HP baseline at "
                 r"$60.1$~$\mu$M and the fully corrected floor near "
                 r"$0.1$~$\mu$M~\cite{pkpd_translational}, so a managed patient sits between "
                 r"the two and $\bar m$ is read in $\mu$mol/L"),
    Row(r"$\theta$", f"{_A.theta:.0f}", "residual capacity of the enzymatic channel",
        ASSUMED, r"illustrative, and the quantity the paper says cannot be identified, so "
                 r"no conclusion may depend on its value. Residual PCC activity is published "
                 r"only as a fraction of control, under $5\%$ in PA, never as an absolute flux, "
                 r"so $\theta$ is an order-of-magnitude quantity and is confounded with "
                 r"$V_{\max}$ by construction. It is a \emph{capacity} and not a flux: at "
                 r"$\bar m/K_m = 0.1$ the channel runs at $9.1\%$ of $V_{\max}$, carrying $273$ "
                 r"against a load of $323$, so it is not comparable with a load figure directly"),
    Row(r"$K_m$", f"{_A.Km:.0f}", "half-saturation of the enzymatic channel, at the plasma level",
        ASSUMED, r"illustrative, and \emph{not} the in vitro Michaelis constant of PCC "
                 r"($0.29$ mM for propionyl-CoA~\cite{pcc_liver_kinetics}), a different molecule in a different "
                 r"compartment. The ratio $\bar m/K_m = 0.1$ is what the excitation condition "
                 r"is stated against"),
    Row(r"$k_e$", f"{_A.ke:.0f}", "alternative clearance rate constant",
        ASSUMED, r"illustrative, chosen so the enzymatic channel carries most of the flux, "
                 r"which is the favourable case for identification. \textbf{This is the "
                 r"priority gap}: Section~\ref{sec:designsec} needs it within a factor of about "
                 r"three, and no direct measurement exists. The nearest route is the "
                 r"first-order 3-HP elimination implicit in the published indirect-response "
                 r"model~\cite{pkpd_translational}"),
    Row(r"$\mtol$", f"{M_TOL:.0f}", "tolerable metabolite level",
        ASSUMED, r"a clinical convention, not a property of the patient. The width of the "
                 r"undetermined window runs $1.2\%$ to $67\%$ as it moves, so this is the most "
                 r"expensive assumption in the local model. No numeric plasma 3-HP target "
                 r"is published; the only quantitative anchor is that a $50\%$ lower 3-HP "
                 r"carries a decompensation rate ratio of $0.7$~\cite{natural_history}"),
    Row(r"$\spk/s_{\mathrm{base}}$", f"{ILLNESS_FOLD:.2f}", "illness amplitude, as a fold of the patient's own baseline load",
        ASSUMED, r"stated as a fold rather than absolutely, because an absolute figure "
                 r"disagreed between sections. Model C's own fold sits at 1.18. Muscle "
                 r"proteolysis during illness bounds this from above: 3-methylhistidine "
                 r"excretion rises by $227\%$ to $325\%$ in sepsis and trauma, so $1.50$ is a "
                 r"conservative value for a moderate intercurrent illness"),
    Row(r"$\kappa$", "not assigned", "gain from recorded carnitine dose onto $k_e$",
        ASSUMED, r"never assigned a value. Only the product $\kappa u_c$ enters, as a lumped "
                 r"increment on $k_e$, and every figure varies that product directly. "
                 r"The ill-conditioning result of supplement \S S1 is conditional on it: an unknown $\kappa$ "
                 r"leaves the carnitine route two equations in four unknowns. The "
                 r"\emph{direction} is established, carnitine supplementation raising "
                 r"propionylcarnitine excretion~\cite{carnitine_excretion_pa}, but no "
                 r"dose-response slope is published, so the magnitude cannot be assigned"),
    Row(r"$\lambda$", r"$\ln 2/3.5$ d$^{-1}$",
        "decay rate of the therapeutic increment",
        LITERATURE, r"a $3.5$-day half-life, so the increment has decayed sixteen-fold by the "
                 r"next dose and the inter-dose recovery is an excursion rather than a "
                 r"ripple. \textbf{Corroborated}: the translational PK/PD model estimates a "
                 r"PCC-protein degradation rate of $0.0075$~h$^{-1}$~\cite{pkpd_translational}, "
                 r"a half-life of $3.85$~d, so this is within $10\%$ of a published value "
                 r"rather than free. It was previously named and never assigned, "
                 r"which made Section~\ref{sec:collapse} unrunnable; its only role is to set "
                 r"how much excursion a fixed schedule creates, and "
                 r"Table~\ref{tab:excursion} reports the answer"),
    Row(r"$T_{\mathrm{dose}}$", "14 d", "interval between therapeutic doses",
        LITERATURE, r"the biweekly schedule selected for mRNA-3927 by model-informed dose "
                    r"selection~\cite{dose_selection_poster}; the phase 1/2 trial dosed every "
                    r"two or three weeks~\cite{trial_phase12}. Adopted rather than assumed"),
    Row(r"$D$", "6000", "size of one therapeutic increment",
        ASSUMED, r"chosen so the marker falls by about half at peak effect, the order of the "
                 r"reductions reported in the phase 1/2 interim "
                 r"analysis~\cite{trial_phase12}. Only the excursion it produces is used. "
                 r"It is recoverable in principle: the published dose-to-marker transfer "
                 r"function~\cite{pkpd_translational} turns a measured per-dose fall in 3-HP "
                 r"into a delivered increment"),
    Row(r"$\sigma$", "a family", "standard deviation of the ordinate of the fit",
        ASSUMED, r"deliberately not assigned. Every design conclusion is reported as a function "
                 r"of it, because nothing in this paper bounds how well the substrate load is "
                 r"recorded. The assay contributes $b\,\sigma_m$ to the same total"),
    Row(r"$\mathrm{CV}$", f"{ASSAY_CV:.0f}\\%", "analytical variability of the marker assay",
        ASSUMED, r"illustrative, and a \emph{percentage}: it is never absolute flux. It enters "
                 r"the ordinate as $b\,\mathrm{CV}\,\bar m$. This is the \emph{analytical} "
                 r"figure: the breath test reports repeatability of $6.9\%$ in controls and "
                 r"$12.8\%$ in patients~\cite{breath_test_mma}. Biological variability is far "
                 r"larger, the translational model carrying $51\%$ between subjects on the "
                 r"3-HP baseline~\cite{pkpd_translational}"),
]

PARENT = [
    Row(r"$\theta$", f"{_C.theta:.0f}", "residual capacity, parent model",
        ASSUMED, r"illustrative, about $10\%$ of normal, chosen to sit above "
                 r"$\thc^{\,\mathrm{int}}$ so the patient has a margin to measure. The "
                 r"published figure is a fraction rather than a flux: residual PCC activity "
                 r"is under $5\%$ of control in PA, and converting that to an absolute "
                 r"$V_{\max}$ has never been done, which would need a PA-specific "
                 r"$1$-$^{13}$C-propionate breath test~\cite{breath_test_mma}"),
    Row(r"$K_m$", f"{_C.Km:.0f}", "half-saturation of the enzymatic channel",
        ASSUMED, "as above, and shared with the local model so the two can be compared"),
    Row(r"$K_c$", f"{_C.Kc:.0f}", "half-saturation of the conjugation step",
        ASSUMED, r"no measurement, and no reason a conjugation reaction should share an "
                 r"enzyme's half-saturation. Cost bounded in supplement \S S1. The nearest "
                 r"measured analogue is the renal carnitine transporter, $K_m$ $0.08$ mM "
                 r"external and $1.2$ mM internal~\cite{octn2_renal_transport}, which is a "
                 r"different step; no half-saturation for the conjugation itself is published"),
    Row(r"$M_i$", f"{_C.Mi:.0f}", "metabolite level at which acylcarnitine competition doubles the non-conjugative loss",
        ASSUMED, r"no measurement. Previously not written down anywhere, which is why "
                 r"the margin figure could not be reproduced from the paper. Cost bounded in "
                 r"supplement \S S1. The mechanism is measured, acylcarnitines being the most "
                 r"effective inhibitors of that transporter~\cite{octn2_renal_transport}, but "
                 r"no propionylcarnitine inhibition constant is published, so the level at "
                 r"which competition doubles the loss is not available as a number"),
    Row(r"$\rho_0$", f"{_C.rho0:.2f}", r"basal non-conjugative loss relative to conjugation, $\delta_0/k$",
        ASSUMED, r"no direct measurement. Bounded above by $1$, since the carnitine drain in "
                 r"this disease is dominated by propionylcarnitine, and below by the measured "
                 r"direction of free-carnitine excretion: in PA and MMA the urinary excretion "
                 r"of free carnitine falls while short-chain acylcarnitine excretion "
                 r"rises~\cite{carnitine_excretion_pa}"),
    Row(r"$R_0$", f"{_C.R0:.1f}", "carnitine replenishment",
        DERIVED, r"\textbf{not free}. Fixed by requiring the alternative route alone to peak at "
                 rf"${ALT_PEAK_C:.0f}$ against a baseline load of ${S_BASE_C:.0f}$, so that a "
                 r"patient with $\theta = 0$ cannot be maintained. Changing $K_c$, $M_i$ or "
                 r"$\rho_0$ moves that peak and $R_0$ must be re-derived"),
    Row(r"$s_{\mathrm{base}}$", f"{S_BASE_C:.0f}", "baseline propiogenic load",
        ASSUMED, r"illustrative, and sets the scale everything else in the parent model is "
                 r"quoted against. Whole-body propionate turnover is published at roughly "
                 r"$250$ to $400$~$\mu$mol/kg/day, partitioned as about half from amino-acid "
                 r"catabolism, a fifth from gut fermentation and the rest from odd-chain fatty "
                 r"acids; a protein-restricted diet lowers the first two deliberately. The "
                 r"value here is not in those units and is not claimed to be"),
    Row(r"$k$", "not assigned", "conjugation rate constant, and the carnitine relaxation rate",
        ASSUMED, r"\textbf{This row previously read \textsc{derived}, ``a free scale in the "
                 r"two-state system''. The first half of that is kept and the second is "
                 r"withdrawn.} $k$ does drop out of the quasi-steady reduction of "
                 r"supplement~\S S1 exactly, so the reduced clearance does not depend on it. "
                 r"It does not follow that nothing does. $k$ sets the carnitine relaxation "
                 r"time $\tau_c = 1/[k(u+\rho)]$, and therefore whether the reduction is valid "
                 r"\emph{away} from equilibrium, which is the only place the two formulations "
                 r"differ at all. Nothing measures it, so it is assumed and not derived. Its "
                 r"cost is bounded rather than argued away, and the bound is not small: at "
                 r"$k=1$, the value the code runs, $\tau_c/\tau_m$ spans $0.52$ to $55$ across "
                 r"the sampled parameter region, so carnitine is comparable to or slower than "
                 r"the metabolite and the separation the reduction would need runs backwards. "
                 r"Published human kinetics say the same independently: whole-body carnitine "
                 r"turnover time is $38$ to $119$~h and even the rapidly equilibrating "
                 r"compartment distributes with a half-life of $2$ to $3$~h, against the "
                 r"marker's $0.81$~h. The resulting error scales as $1/k$: $3.4\%$ to $25.4\%$ "
                 r"on a transient, and up to $49.6\%$ on the basin boundary when carnitine "
                 r"starts away from the slow manifold. No equilibrium quantity is affected, "
                 r"because there the two formulations are one equation"),
]


def _tex(rows):
    out = []
    for r in rows:
        out.append(" & ".join([r.symbol, r.value, r.meaning,
                               r"\textsc{" + r.source + "}", r.note]) + r" \\")
    return out


def write_params_table():
    """Render both blocks as one longtable, source column included.

    A longtable rather than a tabular in a float: once every row carries its
    literature footing the table is taller than a page, and a tabular inside a
    table environment does not break -- it silently runs off the bottom and
    takes the rest of the table with it. The caption and the label are left to
    the including document, which defines \paramscaption and \paramslabel,
    because the paper and the supplement caption this table differently and
    the supplement imports the paper's labels with xr, so a hard-coded label
    would collide.
    """
    from .style import TABDIR, ROOT
    head = [r"\textbf{Symbol} & \textbf{Value} & \textbf{Meaning} & "
            r"\textbf{Source} & \textbf{What fixes it, or what assuming it costs} \\"]
    body = [r"\begin{longtable}{@{}p{0.07\textwidth} p{0.07\textwidth} "
            r"p{0.20\textwidth} p{0.11\textwidth} p{0.45\textwidth}@{}}",
            r"\caption{\paramscaption}\label{\paramslabel}\\",
            r"\toprule"] + head + [
            r"\midrule",
            r"\endfirsthead",
            r"\multicolumn{5}{@{}l}{\emph{Table \thetable, continued from the previous page}}\\",
            r"\toprule"] + head + [
            r"\midrule",
            r"\endhead",
            r"\midrule",
            r"\multicolumn{5}{r@{}}{\emph{continued on the next page}}\\",
            r"\endfoot",
            r"\bottomrule",
            r"\endlastfoot",
            r"\multicolumn{5}{@{}l}{\emph{The local identification model}} \\",
            r"\addlinespace[2pt]"]
    body += _tex(LOCAL)
    body += [r"\addlinespace", r"\midrule",
             r"\multicolumn{5}{@{}l}{\emph{The parent model and its equilibrium reduction}} \\",
             r"\addlinespace[2pt]"]
    body += _tex(PARENT)
    body += [r"\end{longtable}"]
    path = TABDIR / "tab_params.tex"
    path.write_text("\n".join(body))
    print(f"  wrote {path.relative_to(ROOT)}")
    return path


if __name__ == "__main__":
    write_params_table()
