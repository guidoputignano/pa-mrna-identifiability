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

The class says where a number came from. It does NOT say what the number is,
and Phase 7 found that readers were reading the two as the same thing. Every
row therefore carries a second tag, its KIND, which says what kind of object
the constant is:

    measurable  defined so that a single experiment could return it
    effective   an apparent constant at a level where the mechanism is not
    lumped      absorbs several distinct mechanisms into one number
    convention  a reporting or clinical choice, not a property of a patient
    calibration chosen to hit a stated target, and only that target is used
    derived     fixed by a constraint elsewhere
    eliminated  cancels out of the reduced curve. NOT the same as "does not
                matter": k cancels from alpha(m) and still sets the
                reduction's validity off equilibrium
    prior width the standard deviation of a prior, not a physical quantity

A row can be `literature` and `lumped` at once, and one is. The distinction
matters because an effective or lumped constant cannot be checked against a
direct measurement even in principle, so "no measurement exists" is a
different statement for those rows than for a measurable one.

A third tag names the SPECIES the supporting evidence comes from, where there
is any. Four rows in this table lean on estimates fitted in mice and two on a
rat in vitro preparation, and until Phase 7 the table did not say so.

UNITS. The paper's local model is a flux balance at the plasma level,

    (theta + Delta) m/(K_m + m) + k_e m  =  L + b,

so the unit system is fixed once the marker's unit is: m, K_m, M_i, K_c and
m_tol are concentrations in umol/L; theta, Delta, D, L, s_base, R_0 and sigma
are fluxes in umol/L/day; k_e, lambda and k are rates in /day; T_dose is in
days; and CV, the illness fold, rho_0 and kappa are dimensionless. Twelve rows
were dimensional and unitless before Phase 7. One consequence is worth stating
because it is nowhere else in the paper: k_e = 5/day is a marker half-life of
ln2/5 = 0.139 d, about 3.3 hours.

`make figures` writes results/tables/tab_params.tex from this file, and
tests/test_models.py asserts that the values here are the values the model code
actually runs on, so the table cannot drift away from the simulation.
"""

from dataclasses import dataclass

from .models import ParamsA, ParamsC


LITERATURE, DERIVED, ASSUMED = "literature", "derived", "assumed"


#: what kind of object the constant is, independent of where it came from
MEASURABLE, EFFECTIVE, LUMPED = "measurable", "effective", "lumped"
CONVENTION, CALIBRATION, ELIMINATED = "convention", "calibration", "eliminated"
PRIOR = "prior width"


@dataclass(frozen=True)
class Row:
    symbol: str          # as it appears in the paper, in LaTeX
    value: str           # rendered, with units where there are any
    meaning: str
    source: str          # one of the three classes above
    note: str            # what fixes it, or what its being assumed costs
    kind: str = ""       # measurable / effective / lumped / convention / ...
    species: str = ""    # where the supporting evidence comes from, if any


_A = ParamsA()
_C = ParamsC()

# The operating point, the tolerable level and the illness fold live in
# model_A/run.py rather than in a dataclass, so they are named here against
# that module and checked against it in the tests.
M_OP, M_TOL, ILLNESS_FOLD = 10.0, 40.0, 1.50
ASSAY_CV = 5.0
S_BASE_C = 600.0
ALT_PEAK_C = 500.0

# Two design constants that reached the paper's numbers without a row, which
# Phase 7 X2 recorded. They live as a default argument and a module constant in
# src/design/run.py; they are named here against that module and checked
# against it in the tests, exactly as M_OP and its neighbours are.
KE_PRIOR_LOG = 0.8      # sd of the log-normal prior on k_e
LOAD_FRAC = 0.05        # relative error in the RECORDED substrate load


LOCAL = [
    Row(r"$\bar m$", f"${M_OP:.0f}$~$\mu$M",
        "operating point: the level management holds the patient at",
        ASSUMED, r"illustrative, and every design quantity is reported relative to it, so "
                 r"its absolute value cancels out of the conclusions. \textbf{The scale "
                 r"quoted here is murine.} The translational PK/PD model puts the "
                 r"\emph{suppressible} plasma 3-HP pool at $60.1$~$\mu$M and the residual "
                 r"floor at $E_0 = 0.0987$~$\mu$M in \emph{PCC-deficient mice}; that paper "
                 r"scales its PK parameters allometrically and not its PD "
                 r"ones~\cite{pkpd_translational}, so neither is a human measurement. The "
                 r"nearest human anchor is the newborn dried-blood-spot $99.9$th centile, "
                 r"$38.13$~$\mu$mol/L~\cite{dbs_3hp_cutoff}, within a factor of $1.6$ of the "
                 r"murine figure, so the \emph{scale} carries across even though the number "
                 r"is not a patient's. \textbf{The illustrative point sits at the favourable "
                 r"end and nothing else in this table says so}: Section~\ref{sec:scope} puts "
                 r"the reach limit at three to four times the managed level, so a patient "
                 r"presenting untreated near $60$~$\mu$M is already outside the band in which "
                 r"capacity is recoverable",
        kind=CONVENTION, species="mouse"),
    Row(r"$\theta$", f"${_A.theta:.0f}$~$\mu$M/d",
        "capacity of the enzymatic channel, \textbf{at normal capacity}",
        ASSUMED, r"illustrative, and \textbf{a normal-capacity illustration rather than a "
                 r"patient}. That label is new and it resolves an ambiguity the row used to "
                 r"carry: \texttt{src/common/models.py} annotates this constant as "
                 r"$V_{\max}$, the \emph{maximal} enzymatic velocity, and annotates the "
                 r"parent model's $\theta = 300$ as ``$\sim 10\%$ of normal'', which puts "
                 r"normal at $3000$ and this block at $100\%$ of it. Residual PCC activity "
                 r"in PA is under $5\%$ of control~\cite{pcc_residual_activity}, so the "
                 r"disease sits at $\theta \le 150$, which is the row below. The two "
                 r"annotations are not in conflict once the label is stated; what was in "
                 r"conflict was reading an illustration of a healthy channel as a "
                 r"description of the patient the paper is about. "
                 r"\textbf{Some conclusions do depend on this value}: "
                 r"the identifiability results do not, but the "
                 r"supplement's local-model section (\S S1) flips a qualitative conclusion "
                 r"across it, the anchor attenuating at $3000$ and amplifying at the "
                 r"parent's $300$, and Section~\ref{sec:scope} now reports that the "
                 r"design's own precision does too. It "
                 r"remains an order-of-magnitude quantity and is confounded with $V_{\max}$ "
                 r"by construction: at $\bar m/K_m = 0.1$ the channel runs at $9.1\%$ of "
                 r"$V_{\max}$, carrying $273$ against a load of $323$",
        kind=CONVENTION),
    Row(r"$\theta_{\mathrm{PA}}$", r"$\le 150$~$\mu$M/d",
        "capacity in the disease, on the scale of the row above",
        LITERATURE, r"\textbf{This is the row that says who the paper is "
                    r"about.} Residual PCC activity in propionic acidemia is under "
                    r"$5\%$ of control~\cite{pcc_residual_activity}, which on the scale "
                    r"the row above fixes is $\theta \le 150$. \textbf{The design does not "
                    r"reach its target there, and the gap is large}: on the monthly record "
                    r"$\mathrm{sd}(\ln\hat\theta)$ is $0.1699$ at nominal, $1.5830$ at "
                    r"$10\%$ of nominal and $3.0673$ at $5\%$, against a target of $0.20$ "
                    r"(medians over $24$ intake records; \texttt{audit/n4\_theta\_scope.py}). "
                    r"The gate holds only for $\theta$ between $2472$ and $18814$, that is "
                    r"$82\%$ to $627\%$ of nominal, so the disease sits a factor of $16$ "
                    r"below the band. This is Section~\ref{sec:scope}'s own conclusion "
                    r"given a number: that section already says no schedule, excursion or "
                    r"precision recovers a quantity the observable has stopped responding "
                    r"to, and sampling every two days for two years instead of monthly "
                    r"moves the lower edge by $9\%$",
        kind=DERIVED, species="human"),
    Row(r"$K_m$", f"${_A.Km:.0f}$~$\mu$M",
        "half-saturation of the enzymatic channel, at the plasma level",
        ASSUMED, r"illustrative, and \emph{not} the in vitro Michaelis constant of PCC "
                 r"($0.29$ mM for propionyl-CoA~\cite{pcc_liver_kinetics}), a different "
                 r"molecule in a different compartment. The ratio $\bar m/K_m = 0.1$ is what "
                 r"the excitation condition is stated against. \textbf{Cost, which this row "
                 r"did not previously bound and which is the tightest in the local model "
                 r"after the increment}: on the paper's own monthly record, holding the "
                 r"patient managed at $\bar m$ and re-solving the load, the ordinate sd and "
                 r"the load-bias prior at each $K_m$, $\mathrm{sd}(\ln\hat\theta)$ runs "
                 r"$0.0507$ to $1.3608$ over $K_m$ from $5$ to $1000$~$\mu$M, a factor of "
                 r"$26.8$ (median of $24$ load records, load CV $10\%$; "
                 r"\texttt{audit/n1\_km\_cost.py}). It is $0.1699$ here and crosses the "
                 r"$0.20$ target at $K_m = 124$~$\mu$M, so the design has a factor of "
                 r"\textbf{$1.24$} of headroom on this constant; at the in vitro value the "
                 r"row declines to use, $290$~$\mu$M, it would be $0.4158$. "
                 r"\textbf{That headroom is smaller than it reads}: across the same $24$ "
                 r"records the spread is $0.1684$ to $0.1709$ at $K_m = 100$ and $0.1987$ "
                 r"to $0.2010$ at $K_m = 124$, so the crossing straddles the target and "
                 r"which side of it a given record falls on is not settled by the median",
        kind=EFFECTIVE),
    Row(r"$k_e$", f"${_A.ke:.0f}$~d$^{{-1}}$",
        "alternative clearance rate constant",
        ASSUMED, r"illustrative, chosen so the enzymatic channel carries most of the flux, "
                 r"which is the favourable case for identification. A rate of $5$~d$^{-1}$ is "
                 r"a marker half-life of $3.3$~h. \textbf{This is the priority gap}: "
                 r"Section~\ref{sec:designsec} needs it within a factor of about three, and "
                 r"\textbf{no source supplies it}. The nearest candidate is the "
                 r"first-order 3-HP elimination of the translational model; \emph{that model "
                 r"contains none}. Its 3-HP is a direct algebraic $I_{\max}$ function of PCC "
                 r"protein with no marker state and no elimination "
                 r"constant~\cite{pkpd_translational}, as \texttt{paper/claims.csv} C46 "
                 r"already recorded. The only published structure carrying a first-order "
                 r"marker loss is the turnover model of the dose-selection "
                 r"abstract~\cite{dose_selection_poster}, which publishes no value for it",
        kind=LUMPED),
    Row(r"$\mtol$", f"${M_TOL:.0f}$~$\mu$M", "tolerable metabolite level",
        ASSUMED, r"a clinical convention, not a property of the patient. The width of the "
                 r"undetermined window runs $1.2\%$ to $67\%$ as it moves, so this is the most "
                 r"expensive assumption in the local model. No numeric plasma 3-HP target "
                 r"is published; the only quantitative anchor is that a $50\%$ lower 3-HP "
                 r"carries a decompensation rate ratio of $0.7$~\cite{natural_history}",
        kind=CONVENTION),
    Row(r"$\spk/s_{\mathrm{base}}$", f"${ILLNESS_FOLD:.2f}$",
        "illness amplitude, as a fold of the patient's own baseline load",
        ASSUMED, r"stated as a fold rather than absolutely, because an absolute figure "
                 r"disagreed between sections. Model C's own fold sits at 1.18. Muscle "
                 r"proteolysis during illness bounds this from above: 3-methylhistidine "
                 r"excretion rises by $227\%$ to $325\%$ in sepsis and "
                 r"trauma~\cite{muscle_catabolism}, so $1.50$ is a conservative value for a "
                 r"moderate intercurrent illness. \textbf{The bound is looser than it looks}: "
                 r"3-methylhistidine indexes \emph{muscle} protein catabolism, amino-acid "
                 r"catabolism is about half of propionate production~\cite{propionate_sources} "
                 r"and only part of that is muscle, so the rise does not license an equal "
                 r"bound on total propiogenic load. The subjects are also adults, where the "
                 r"trial population's median age is $8.5$~years",
        kind=EFFECTIVE, species="human, adult"),
    Row(r"$\kappa$", "not assigned", "gain from recorded carnitine dose onto $k_e$",
        ASSUMED, r"never assigned a value. Only the product $\kappa u_c$ enters, as a lumped "
                 r"increment on $k_e$, and every figure varies that product directly. "
                 r"The ill-conditioning result of supplement \S S1 is conditional on it: an "
                 r"unknown $\kappa$ leaves the carnitine route two equations in four unknowns. "
                 r"The \emph{direction} is established, carnitine administration raising "
                 r"short-chain acylcarnitine excretion~\cite{coa_sequestration}, but no "
                 r"dose-response slope is published, so the magnitude cannot be assigned",
        kind=LUMPED, species="human"),
    Row(r"$\lambda$", r"$\ln 2/3.5$~d$^{-1}$",
        "decay rate of the therapeutic increment",
        LITERATURE, r"a $3.5$-day half-life, so the increment has decayed sixteen-fold by the "
                 r"next dose and the inter-dose recovery is an excursion rather than a "
                 r"ripple. \textbf{Partly corroborated, and only early.} The translational "
                 r"model estimates $k_{\deg} = 0.0075$~h$^{-1}$ for PCC "
                 r"protein~\cite{pkpd_translational}, but that is one rate constant inside a "
                 r"\emph{two-compartment} protein model whose own stated protein half-life is "
                 r"$\sim 7$~d, and it is fitted in mice, rats and monkeys and scaled "
                 r"allometrically to \emph{predict} human behaviour. The model's local "
                 r"half-life passes $3.5$~d at about five days after a dose and reaches "
                 r"$9.0$~d by the end of a fortnightly interval, so $\lambda$ tracks the "
                 r"published decline over the first week and decays \emph{faster} than it "
                 r"thereafter, which is the direction that makes the excursion look larger "
                 r"than it is. Its only role is to set how much excursion a fixed schedule "
                 r"creates, and Table~\ref{tab:excursion} reports the answer",
        kind=MEASURABLE, species="mouse, rat, monkey"),
    Row(r"$T_{\mathrm{dose}}$", "14~d", "interval between therapeutic doses",
        LITERATURE, r"the biweekly schedule selected for mRNA-3927 by model-informed dose "
                    r"selection~\cite{dose_selection_poster}; the phase 1/2 trial dosed every "
                    r"two or three weeks~\cite{trial_phase12}. Adopted rather than assumed",
        kind=MEASURABLE, species="human"),
    Row(r"$D$", r"$6000$~$\mu$M/d", "size of one therapeutic increment",
        ASSUMED, r"chosen so the marker falls by \textbf{about two thirds} at peak effect. "
                 r"The computed range is $62$--$66\%$ and it varies for two reasons, both "
                 r"conventions rather than measurements: \emph{which baseline} the fall is "
                 r"read from (the untreated level, or the patient's own pre-dose trough) and "
                 r"\emph{which peak} (the first dose, or the steady-state peak of the train). "
                 r"The three readings are $64.90\%$, $66.36\%$ and $62.25\%$. \textbf{The "
                 r"trial's own reported falls are smaller}: $-29$, $-3$, $-20$, $+5$ and "
                 r"$-20\%$ across the five cohorts after the third dose, with no dose "
                 r"trend~\cite{trial_phase12}, which suggests this value is too high. "
                 r"\textbf{The value matters, not merely its order}: over the range of $D$ "
                 r"those falls make defensible, a monthly record needs from $10$ to more than "
                 r"$60$ visits to reach $\mathrm{sd}(\ln\hat\theta) = 0.20$ "
                 r"(\texttt{audit/PHASE5\_ROBUSTNESS.md}), so an earlier version of this row "
                 r"saying only the excursion is used understated the dependence. It is "
                 r"recoverable in principle from a measured per-dose fall, but the trial "
                 r"reports no sampling time for its falls, so the reading is not determined",
        kind=CALIBRATION, species="mouse"),
    Row(r"$\sigma$", r"a family, $\mu$M/d",
        "standard deviation of the ordinate of the fit",
        ASSUMED, r"deliberately not assigned, because nothing in this paper bounds how well "
                 r"the substrate load is recorded. The assay contributes $b\,\sigma_m$ to the "
                 r"same total. \textbf{Not every design conclusion is reported as a "
                 r"conclusion is reported as a function of it. That is false of "
                 r"function of it}: \texttt{spec\_numbers()} computes one $\sigma$ at a recording error "
                 r"of $5\%$ and reuses it on every row, so that table is reported \emph{at} "
                 r"one $\sigma$, not as a function of it. Swept properly, "
                 r"$\mathrm{sd}(\ln\hat\theta)$ runs $0.1566$ at $1\%$ recording error to "
                 r"$0.2637$ at $20\%$, and the design meets its target only below roughly "
                 r"$12\%$ (median over $24$ simulated intake records at a load CV of "
                 r"$10\%$, spread $[0.1563, 0.1571]$ at $1\%$ and $[0.2563, 0.2690]$ at "
                 r"$20\%$; \texttt{audit/m1\_loadbias\_prior.py})",
        kind=LUMPED),
    Row(r"$\mathrm{CV}$", f"${ASSAY_CV:.0f}\\%$",
        "analytical variability of the marker assay",
        ASSUMED, r"illustrative, and a \emph{percentage}: it is never absolute flux. It enters "
                 r"the ordinate as $b\,\mathrm{CV}\,\bar m$. This is the \emph{analytical} "
                 r"figure. The nearest published repeatability is the "
                 r"$1$-$^{13}$C-propionate breath test, whose Results give "
                 r"$14.02 \pm 12.86\%$ in patients and $6.9\%$ in "
                 r"controls~\cite{breath_test_mma}; \textbf{the widely quoted $12.8\%$ is the "
                 r"standard deviation of the per-subject coefficients and not their mean, and "
                 r"is not their mean}. That figure is also "
                 r"\emph{between-visit}, over two months to four years, so it bounds the "
                 r"same-day repeatability Tier~1 asks for from above rather than measuring "
                 r"it, and the cohort is methylmalonic acidemia. Biological variability is far "
                 r"larger: the translational model carries $51\%$ between subjects on the "
                 r"3-HP baseline~\cite{pkpd_translational}, \emph{between mice of one inbred "
                 r"strain}, which bounds between-patient variability from below rather than "
                 r"above",
        kind=MEASURABLE, species="human (MMA); mouse"),
    Row(r"$s_{\ln k_e}$", f"${KE_PRIOR_LOG:.1f}$",
        "sd of the log-normal prior on $k_e$",
        ASSUMED, r"\textbf{The only prior that does any work.} "
                 r"A factor of $e^{0.8} = 2.23$ per standard deviation, which is wider "
                 r"than the factor of about three the gate needs, so the prior is "
                 r"deliberately weaker than the requirement. It lives as a default argument "
                 r"in \texttt{src/design/run.py} and is asserted against it in the tests. "
                 r"What it is grounded in is the open question, since by the $k_e$ row above "
                 r"no source supplies the constant it is a prior on. \textbf{What it is "
                 r"worth is an open question too}: Section~\ref{sec:designsec} states that "
                 r"removing it takes the same design from $0.1474$ to $0.6451$, and neither "
                 r"figure reproduces from this repository. On the densest design that "
                 r"section lists, sampling every two days for two years, "
                 r"$\mathrm{sd}(\ln\hat\theta)$ is $0.1500$ with the prior and $1.1096$ "
                 r"without, medians over $24$ simulated intake records with spreads "
                 r"$[0.1498, 0.1501]$ and $[1.0313, 1.1923]$, so the prior is worth a factor "
                 r"of $7.4$ rather than the $4.4$ the section prints "
                 r"(\texttt{audit/n3\_ke\_prior\_worth.py}, which also reports that every "
                 r"design figure in this paper is drawn at the one seed "
                 r"\texttt{design\_record} hard-codes)",
        kind=PRIOR),
    Row(r"$\sigma_L/L$", f"${100*LOAD_FRAC:.0f}\\%$",
        "relative error in the recorded substrate load",
        ASSUMED, r"\textbf{Sets $\sigma$} through "
                 r"$\sigma_{\mathrm{total}}$ and the load-bias prior through "
                 r"$s_b = \sigma_L/L \cdot L$. It is the recording error, distinct from the "
                 r"day-to-day variation of the load itself, and the two were conflated in "
                 r"Section~\ref{sec:spec} until R17. Measurable in principle by a weighed "
                 r"intake study; no such study exists in PA, where records are typically "
                 r"prescriptions rather than intake as eaten",
        kind=MEASURABLE),
]


PARENT = [
    Row(r"$\theta$", f"${_C.theta:.0f}$~$\mu$M/d", "residual capacity, parent model",
        ASSUMED, r"illustrative, chosen to sit above $\thc^{\,\mathrm{int}}$ so the patient "
                 r"has a margin to measure. On the scale \texttt{models.py} sets, where the "
                 r"local block's $3000$ is nominal, this is about $10\%$ of normal. "
                 r"\textbf{Residual PCC activity in PA is under $5\%$ of "
                 r"control}~\cite{pcc_residual_activity}, so even this block sits above the "
                 r"disease range, which is the open question flagged on the local $\theta$ "
                 r"row. That source is fibroblast enzymology rather than in vivo hepatic "
                 r"capacity; the $5\%$ reference adopted throughout the admissibility checks "
                 r"is a judgement informed by it and not a measurement of the quantity the "
                 r"model carries. Converting a fraction of control into an absolute "
                 r"$V_{\max}$ has never been done and would need a PA-specific "
                 r"$1$-$^{13}$C-propionate breath test~\cite{breath_test_mma}",
        kind=LUMPED),
    Row(r"$K_m$", f"${_C.Km:.0f}$~$\mu$M", "half-saturation of the enzymatic channel",
        ASSUMED, r"as above, and shared with the local model so the two can be compared. "
                 r"Sharing it across a tenfold difference in capacity assumes $K_m$ does not "
                 r"depend on residual activity, and \textbf{that assumption is supported}: "
                 r"four variant PCCs show no significant difference in $K_m$ for "
                 r"propionyl-CoA against wild type, although they carry about half the "
                 r"$k_{\mathrm{cat}}/K_m$~\cite{pcc_variants_kinetics}. The cost bound on the "
                 r"local $K_m$ row applies here too",
        kind=EFFECTIVE, species="human, in vitro"),
    Row(r"$K_c$", f"${_C.Kc:.0f}$~$\mu$M", "half-saturation of the conjugation step",
        ASSUMED, r"no measurement, and no reason a conjugation reaction should share an "
                 r"enzyme's half-saturation. Cost bounded in supplement \S S1: over $K_c$ "
                 r"from $25$ to $800$, a $32$-fold range, $\thc$ spans $186$ to $120$, a "
                 r"factor of $1.56$. The nearest measured analogue is the renal carnitine "
                 r"transporter, $K_m$ $0.08$ mM external and $1.2$ mM "
                 r"internal~\cite{octn2_renal_transport}, which is a different step and is "
                 r"\textbf{a rat preparation measured in vitro}; no half-saturation for the "
                 r"conjugation itself is published",
        kind=EFFECTIVE, species="rat, in vitro"),
    Row(r"$M_i$", f"${_C.Mi:.0f}$~$\mu$M",
        "metabolite level at which acylcarnitine competition doubles the non-conjugative loss",
        ASSUMED, r"no measurement. Previously not written down anywhere, which is why "
                 r"the margin figure could not be reproduced from the paper. \textbf{Its cost "
                 r"is bounded for the threshold and not for the margin}, and the two differ: "
                 r"over an eightfold range of $M_i$ supplement \S S1 moves $\thc$ by $1.37$ "
                 r"but the \emph{margin} by $6.44$, running $509$ to $3275$ where Figure~16 "
                 r"quotes $954$. The mechanism is measured, acylcarnitines being the most "
                 r"effective inhibitors of that transporter~\cite{octn2_renal_transport}, "
                 r"\textbf{in a rat preparation in vitro}, but no propionylcarnitine "
                 r"inhibition constant is published, so the level at which competition "
                 r"doubles the loss is not available as a number",
        kind=EFFECTIVE, species="rat, in vitro"),
    Row(r"$\rho_0$", f"${_C.rho0:.2f}$",
        r"basal non-conjugative loss relative to conjugation, $\delta_0/k$",
        ASSUMED, r"no direct measurement, and \textbf{the bound is one-sided}. Above it is a "
                 r"number: $\rho_0 < 1$, since the carnitine drain in this disease is "
                 r"dominated by propionylcarnitine. Below it is only a \emph{direction}, that "
                 r"in PA and MMA urinary free-carnitine excretion falls while short-chain "
                 r"acylcarnitine excretion rises~\cite{carnitine_excretion_pa}, which does "
                 r"not yield a number. A lower figure of ``about $0.35$'' quoted elsewhere in "
                 r"the code re-derives at $0.195$",
        kind=EFFECTIVE, species="human"),
    Row(r"$R_0$", f"${_C.R0:.1f}$~$\mu$M/d", "carnitine replenishment",
        DERIVED, r"\textbf{not free}. Fixed by requiring the alternative route alone to peak "
                 rf"at $\mathrm{{ALT}}_{{\mathrm{{peak}}}} = {ALT_PEAK_C:.0f}$ against a "
                 rf"baseline load of ${S_BASE_C:.0f}$, so that a patient with $\theta = 0$ "
                 r"cannot be maintained. Changing $K_c$, $M_i$ or $\rho_0$ moves that peak "
                 r"and $R_0$ must be re-derived, a failure the review log documents in three "
                 r"places. \textbf{The stored constant carries eight significant figures from "
                 r"two round inputs}; the digits past the fourth are an artefact of the "
                 r"derivation and carry no precision. It is left unrounded because every "
                 r"parent-model figure was produced with it",
        kind="derived"),
    Row(r"$\mathrm{ALT}_{\mathrm{peak}}$", f"${ALT_PEAK_C:.0f}$~$\mu$M/d",
        "peak of the alternative route alone, the constraint that fixes $R_0$",
        ASSUMED, r"\textbf{The constraint behind $R_0$}, which would otherwise be a "
                    r"\textsc{derived} row derived "
                 r"from a number with no row of its own, in a registry whose stated rule is "
                 rf"that no number reaches the paper without one. Set below the baseline load "
                 rf"of ${S_BASE_C:.0f}$ so that a patient with no residual capacity cannot be "
                 r"maintained by conjugation alone, which is the severity of null genotypes. "
                 r"The ratio to $s_{\mathrm{base}}$ is what matters and it is "
                 rf"${ALT_PEAK_C/S_BASE_C:.3f}$",
        kind=CALIBRATION),
    Row(r"$s_{\mathrm{base}}$", f"${S_BASE_C:.0f}$~$\mu$M/d", "baseline propiogenic load",
        ASSUMED, r"illustrative, and sets the scale everything else in the parent model is "
                 r"quoted against. The partition of whole-body propionate production is "
                 r"measured: amino-acid catabolism a mean of $51.7\%$ (range $24.5$ to "
                 r"$66.4\%$), gut flora at least $22.2 \pm 8.5\%$ by the fall under "
                 r"metronidazole, and roughly $30\%$ unaccounted for and attributed to "
                 r"odd-chain fatty acids~\cite{propionate_sources}; a protein-restricted diet "
                 r"lowers the first two deliberately. A whole-body turnover of roughly $250$ "
                 r"to $400$~$\mu$mol/kg/day has no source and is not used "
                 r"pending one. The value here is not in those units and is not claimed to be",
        kind=LUMPED, species="human"),
    Row(r"$k$", r"$1$~d$^{-1}$, cancels from the curve", "conjugation rate constant",
        ASSUMED, r"drops out of the quasi-steady reduction of supplement \S S1 exactly, so the "
                 r"reduced clearance does not depend on it. \textbf{It is not therefore a "
                 r"this row concluded from that it is a free scale and reclassified it "
                 r"\textsc{derived}}: $k$ alone sets the carnitine "
                 r"relaxation rate, $\tau_c = 1/(k[u(m)+\rho(m)])$, so it fixes the "
                 r"reduction's validity \emph{off} equilibrium even though it cancels from "
                 r"the curve, and nothing in this model fixes $k$ itself. \textbf{Cost}: "
                 r"starting the two-state system with carnitine at a pre-perturbation value "
                 r"and comparing against the reduced system, the peak relative error in $m$ "
                 r"is $21.6\%$ at the assumed $k = 1$, falling as $1/k$ for larger $k$ "
                 r"($8.5\%$ at $k=10$, $1.3\%$ at $k=100$) and saturating below it "
                 r"($26.9\%$ at $k=0.1$), because a relative error cannot exceed order one. "
                 r"Everything this paper computes from the parent is computed \emph{at} "
                 r"equilibrium, where the reduction is an identity and the cost is zero; the "
                 r"bound above is what a transient would cost and is the third reason "
                 r"Section~\ref{sec:discussion} excludes the decompensation route",
        kind=ELIMINATED, species="human"),
]


def _source_cell(r):
    """Class, kind and species in one column, one per line.

    Phase 7's first finding about the table's shape was that a single Source
    column reading `assumed` twenty times tells a reader nothing about what
    kind of constant each row holds, or that five of them rest on mice. The
    column carries all three now; a sixth column would not fit the page.
    """
    parts = [r"\textsc{" + r.source + "}"]
    if r.kind:
        parts.append(r"\emph{" + r.kind + "}")
    if r.species:
        parts.append(r"\emph{" + r.species + "}")
    return r"\newline ".join(parts)


#: the rows whose cost bound decides whether the paper's results hold, and so
#: the only notes the MAIN-TEXT table carries. Everything else's provenance
#: moves to the supplement, where a reader who wants it will look and a
#: reviewer counting pages will not.
MAIN_TEXT_NOTES = (r"$\theta$", r"$\theta_{\mathrm{PA}}$", r"$K_m$")


def _tex(rows, compact=False):
    """
    One LaTeX row per registry row.

    `compact` drops the provenance note from every row except the three whose
    cost bound is load-bearing, and points the rest at the supplement. It does
    NOT drop the kind or species tags: those are one word each, they are what
    stop a reader taking an assumed lumped murine constant for a measured
    human one, and they were the whole point of Phase 7.
    """
    out = []
    for r in rows:
        note = r.note
        if compact:
            note = note if r.symbol in MAIN_TEXT_NOTES else r"\S S2"
        out.append(" & ".join([r.symbol, r.value, r.meaning,
                               _source_cell(r), note]) + r" \\")
    return out


def write_params_table(compact=False, name=None):
    """Render both blocks as one longtable, source column included.

    A longtable rather than a tabular in a float: once every row carries its
    literature footing the table is taller than a page, and a tabular inside a
    table environment does not break -- it silently runs off the bottom and
    takes the rest of the table with it. The caption and the label are left to
    the including document, which defines \paramscaption and \paramslabel,
    because the paper and the supplement caption this table differently and
    the supplement imports the paper's labels with xr, so a hard-coded label
    would collide.

    `compact` writes the main-text version, which keeps every row and every
    tag and carries the note only for the three rows of MAIN_TEXT_NOTES. The
    full version goes in the supplement.
    """
    from .style import TABDIR, ROOT
    last = (r"\textbf{Cost, where it decides the result}" if compact
            else r"\textbf{What fixes it, or what assuming it costs}")
    widths = ((r"p{0.09\textwidth} p{0.10\textwidth} p{0.26\textwidth} "
               r"p{0.13\textwidth} p{0.34\textwidth}") if compact else
              (r"p{0.07\textwidth} p{0.07\textwidth} p{0.20\textwidth} "
               r"p{0.11\textwidth} p{0.45\textwidth}"))
    head = [r"\textbf{Symbol} & \textbf{Value} & \textbf{Meaning} & "
            r"\textbf{Source} & " + last + r" \\"]
    body = [r"\begin{longtable}{@{}" + widths + r"@{}}",
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
    body += _tex(LOCAL, compact)
    body += [r"\addlinespace", r"\midrule",
             r"\multicolumn{5}{@{}l}{\emph{The parent model and its equilibrium reduction}} \\",
             r"\addlinespace[2pt]"]
    body += _tex(PARENT, compact)
    body += [r"\end{longtable}"]
    path = TABDIR / (name or ("tab_params_compact.tex" if compact
                              else "tab_params.tex"))
    path.write_text("\n".join(body))
    print(f"  wrote {path.relative_to(ROOT)}")
    return path


if __name__ == "__main__":
    write_params_table()
    write_params_table(compact=True)
