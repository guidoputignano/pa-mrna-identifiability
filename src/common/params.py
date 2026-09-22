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
    #: one line saying what the value rests on: a citation, an explicit
    #: statement that nothing supplies it, or what it was chosen against.
    #: The main-text table carries this in place of the abridged cost column,
    #: whose cells were mostly a pointer at the supplement and told a reader
    #: nothing about the row they were reading.
    rests_on: str = ""


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
    Row(r"$\bar m$", rf"${M_OP:.0f}$~$\mu$M",
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
        kind=CONVENTION, species="mouse",
        rests_on=r"chosen for illustration, and every design quantity is relative to it. The scale is murine~\cite{pkpd_translational}, within $1.6$ of the human newborn cut-off~\cite{dbs_3hp_cutoff}"),
    Row(r"$\theta$", rf"${_A.theta:.0f}$~$\mu$M/d",
        r"capacity of the enzymatic channel, \textbf{at normal capacity}",
        ASSUMED, r"illustrative, and \textbf{a normal-capacity illustration rather than a "
                 r"patient}: \texttt{src/common/models.py} carries this constant as "
                 r"$V_{\max}$, the \emph{maximal} enzymatic velocity, and the parent "
                 r"model's $\theta = 300$ as $\sim 10\%$ of normal, which puts normal at "
                 r"$3000$ and this block at $100\%$ of it. Residual PCC activity in PA is "
                 r"under $5\%$ of control~\cite{pcc_residual_activity}, so the disease "
                 r"sits at $\theta \le 150$, which is the row below. "
                 r"\textbf{Some conclusions do depend on this value}: "
                 r"the identifiability results do not, but the "
                 r"supplement's local-model section (\S S1) flips a qualitative conclusion "
                 r"across it, the anchor attenuating at $3000$ and amplifying at the "
                 r"parent's $300$, and by Section~\ref{sec:scope} the design's own "
                 r"precision does too. It "
                 r"remains an order-of-magnitude quantity and is confounded with $V_{\max}$ "
                 r"by construction: at $\bar m/K_m = 0.1$ the channel runs at $9.1\%$ of "
                 r"$V_{\max}$, carrying $273$ against a load of $323$",
        kind=CONVENTION,
        rests_on=r"chosen for illustration as \emph{normal} capacity on the parent model's scale, against residual activity under $5\%$~\cite{pcc_residual_activity}"),
    Row(r"$\theta_{\mathrm{PA}}$", r"$\le 150$~$\mu$M/d",
        "capacity in the disease, on the scale of the row above",
        LITERATURE, r"\textbf{This is the row that says who the paper is "
                    r"about.} Residual PCC activity in propionic acidemia is under "
                    r"$5\%$ of control~\cite{pcc_residual_activity}, which on the scale "
                    r"the row above fixes is $\theta \le 150$. \textbf{The design does not "
                    r"reach its target there, and the gap is large}: on the monthly record "
                    r"$\mathrm{sd}(\ln\hat\theta)$ is $\num{sd_monthly}$ at nominal, "
                    r"$\num{sd_10pct}$ at $10\%$ of nominal and $\num{sd_5pct}$ at $5\%$, "
                    r"against a target of $0.20$ "
                    r"(medians over $24$ intake records; \texttt{audit/n4\_theta\_scope.py}). "
                    r"The gate holds only for $\theta$ between $\num{band_lo}$ and "
                    r"$\num{band_hi}$, that is $\num{band_lo_pct}\%$ to "
                    r"$\num{band_hi_pct}\%$ of nominal, so the disease sits a factor of "
                    r"$\num{disease_factor}$ "
                    r"below the band. This is Section~\ref{sec:scope}'s conclusion "
                    r"given a number: no schedule, excursion or precision recovers a "
                    r"quantity the observable has stopped responding to, and sampling "
                    r"every two days for two years instead of monthly moves the lower "
                    r"edge by $\num{band_dense_gain_pct}\%$",
        kind=DERIVED, species="human",
        rests_on=r"residual PCC activity in propionic acidemia is under $5\%$ of control~\cite{pcc_residual_activity}, which on the scale of the row above is $\theta \le 150$"),
    Row(r"$K_m$", rf"${_A.Km:.0f}$~$\mu$M",
        "half-saturation of the enzymatic channel, at the plasma level",
        ASSUMED, r"illustrative, and \emph{not} the in vitro Michaelis constant of PCC "
                 r"($0.29$ mM for propionyl-CoA~\cite{pcc_liver_kinetics}), a different "
                 r"molecule in a different compartment. The ratio $\bar m/K_m = 0.1$ is what "
                 r"the excitation condition is stated against. \textbf{Its cost is the "
                 r"tightest in the local model after the increment}: on the paper's own "
                 r"monthly record, holding the "
                 r"patient managed at $\bar m$ and re-solving the load, the ordinate sd and "
                 r"the load-bias prior at each $K_m$, $\mathrm{sd}(\ln\hat\theta)$ runs "
                 r"$\num{km_sweep_lo}$ to $\num{km_sweep_hi}$ over $K_m$ from $5$ to "
                 r"$1000$~$\mu$M, a factor of $\num{km_sweep_fold}$ (median of $24$ load "
                 r"records, load CV $10\%$; \texttt{audit/n1\_km\_cost.py}). It is "
                 r"$\num{sd_monthly}$ here and crosses the $0.20$ target at $K_m = "
                 r"\num{km_crossing}$~$\mu$M, so the design has a factor of "
                 r"\textbf{$\num{km_headroom}$} of headroom on this constant; at the in "
                 r"vitro value the row declines to use, $290$~$\mu$M, it would be "
                 r"$\num{km_invitro_sd}$. "
                 r"\textbf{That headroom is smaller than it reads}: across the same $24$ "
                 r"records the spread is $\num{sd_monthly_lo}$ to $\num{sd_monthly_hi}$ at "
                 r"$K_m = 100$ and $\num{km_crossing_sd_lo}$ to $\num{km_crossing_sd_hi}$ "
                 r"at the crossing, so it straddles the target and "
                 r"which side of it a given record falls on is not settled by the median",
        kind=EFFECTIVE,
        rests_on=r"chosen for illustration at $\bar m/K_m = 0.1$, and \emph{not} PCC's in vitro Michaelis constant, $0.29$~mM~\cite{pcc_liver_kinetics}, a different molecule in a different compartment"),
    Row(r"$k_e$", f"${_A.ke:.0f}$~d$^{{-1}}$",
        "alternative clearance rate constant",
        ASSUMED, r"chosen so the enzymatic channel carries most of the flux, which is the "
                 r"favourable case for identification. A rate of $5$~d$^{-1}$ is a marker "
                 r"half-life of $3.3$~h. \textbf{This is the priority gap}: "
                 r"Section~\ref{sec:designsec} needs it within a factor of about three, and "
                 r"\textbf{no independent measurement in humans supplies it}. The nearest "
                 r"candidate is the first-order 3-HP elimination of the translational model; "
                 r"\emph{that model contains none}. Its 3-HP is a direct sigmoidal "
                 r"$I_{\max}$ function of PCC protein with no marker state and no "
                 r"elimination constant, its only 3-HP parameters being $E_0 = 0.0987$ and a "
                 r"suppressible pool of $60.1$~$\mu$M/L fitted to the authors' own "
                 r"PCC-deficient mouse data, and the ${\sim}7$-day half-lives there "
                 r"belonging to PCC protein rather than to the "
                 r"marker~\cite{pkpd_translational}; \texttt{paper/claims.csv} C46 already "
                 r"recorded the structural fact. \textbf{That is not an oversight but the "
                 r"same gap, met by a published model of this disease and worked around the "
                 r"same way}: its authors record that their refinements were limited by the "
                 r"sparse literature on the interplay between 2-MC, 3-HP and the C3/C2 "
                 r"pathways. \textbf{And the two ways out are exclusive.} The one published "
                 r"structure carrying a first-order marker loss is the turnover model of the "
                 r"dose-selection abstract~\cite{dose_selection_poster}, which publishes no "
                 r"value for it; the other avoids needing the constant at all by making the "
                 r"marker algebraic in the enzyme, and Section~\ref{sec:structuredep} prices "
                 r"what that costs: under it $\theta$ is not identifiable at any design. "
                 r"Either a structure carries a marker state and needs a constant nobody has "
                 r"measured, or it drops the state and loses the parameter. "
                 r"\textbf{An order of magnitude is available, from the nearest measured "
                 r"analogue rather than from 3-HP.} Plasma clearance of methylmalonic acid "
                 r"during haemodialysis in a $\mathrm{mut}^0$ patient is $103$~mL/min at a "
                 r"volume of distribution of $0.48$~L/kg~\cite{ke_mma_analogue}, which at "
                 r"$\num{ke_mma_bw}$~kg is $\num{ke_mma}$~d$^{-1}$ and a half-life of "
                 r"$\num{ke_mma_half_h}$~h. The value carried here sits "
                 r"$\num{ke_mma_vs_carried_pct}\%$ above it, and the two coincide at "
                 r"$\num{ke_mma_cross_bw}$~kg, so for an adult heavier than that the "
                 r"illustrative rate is the faster of the two. \textbf{The caveats travel "
                 r"with the number}: that clearance is \emph{dialytic} and not endogenous, "
                 r"so it bounds the endogenous route from above rather than measuring it, "
                 r"and methylmalonic acid is dicarboxylic with strong tubular handling where "
                 r"3-HP is monocarboxylic (pK$_a$ $4.9$, MW $90$). It is a plausibility band "
                 r"and not a substitute. It is carried because "
                 r"Section~\ref{sec:scopeprec} shows the working band moving from "
                 r"$\num{band_lo_qtr_pct}\%$ to $\num{band_lo_dbl_pct}\%$ of nominal "
                 r"across a factor of eight in $k_e$: this row is load-bearing and should "
                 r"look like it. At the analogue's own $\num{ke_mma}$~d$^{-1}$ the band "
                 r"runs $\num{band_lo_mma_pct}\%$--$\num{band_hi_mma_pct}\%$ against "
                 r"$\num{band_lo_pct}\%$--$\num{band_hi_pct}\%$ at the carried value, so "
                 r"the correction moves the edge the favourable way and does not reach the "
                 r"disease, which sits at $5\%$",
        kind=LUMPED,
        rests_on=r"\textbf{no independent human measurement}: the translational model has no marker state and no elimination constant~\cite{pkpd_translational}, and the one structure with a first-order marker loss publishes no value~\cite{dose_selection_poster}. Bounded from above by the methylmalonic-acid analogue, $\num{ke_mma}$~d$^{-1}$~\cite{ke_mma_analogue}"),
    Row(r"$\mtol$", rf"${M_TOL:.0f}$~$\mu$M", "tolerable metabolite level",
        ASSUMED, r"a clinical convention, not a property of the patient. The width of the "
                 r"undetermined window runs $\num{mtol_window_lo_pct}\%$ to "
                 r"$\num{mtol_window_hi_pct}\%$ as it moves, so this is the most "
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
        kind=MEASURABLE, species="mouse, rat, monkey",
        rests_on=r"the translational model's PK, allometrically scaled from mouse, rat and monkey~\cite{pkpd_translational}"),
    Row(r"$T_{\mathrm{dose}}$", "14~d", "interval between therapeutic doses",
        LITERATURE, r"the biweekly schedule selected for mRNA-3927 by model-informed dose "
                    r"selection~\cite{dose_selection_poster}; the phase 1/2 trial dosed every "
                    r"two or three weeks~\cite{trial_phase12}. Adopted rather than assumed",
        kind=MEASURABLE, species="human",
        rests_on=r"the biweekly schedule of model-informed dose selection~\cite{dose_selection_poster}; the trial dosed every two or three weeks~\cite{trial_phase12}"),
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
        kind=CALIBRATION, species="mouse",
        rests_on=r"chosen for illustration so the marker falls by about two thirds at peak, against the trial's reported falls~\cite{trial_phase12}"),
    Row(r"$\sigma$", r"a family, $\mu$M/d",
        "standard deviation of the ordinate of the fit",
        ASSUMED, r"deliberately not assigned, because nothing in this paper bounds how well "
                 r"the substrate load is recorded. The assay contributes $b\,\sigma_m$ to the "
                 r"same total. \textbf{Not every design conclusion is reported as a "
                 r"function of it}: \texttt{spec\_numbers()} computes one $\sigma$ at a "
                 r"recording error of $5\%$ and reuses it on every row, so that table is "
                 r"reported \emph{at} one $\sigma$, not as a function of it. Swept "
                 r"properly, moving the accuracy in both places it enters, the ordinate "
                 r"noise and the load-bias prior, $\mathrm{sd}(\ln\hat\theta)$ runs "
                 r"$\num{sd_load_acc_1pct}$ at $1\%$ recording error to "
                 r"$\num{sd_load_acc_20pct}$ at $20\%$ and $\num{sd_load_acc_40pct}$ "
                 r"at $40\%$, and the design meets its target only below about "
                 r"$\num{load_accuracy_gate}\%$ "
                 r"(medians over $24$ simulated intake records at a load CV of $10\%$, "
                 r"spread $[\num{sd_load_acc_1pct_lo}, \num{sd_load_acc_1pct_hi}]$ at "
                 r"$1\%$ and $[\num{sd_load_acc_20pct_lo}, \num{sd_load_acc_20pct_hi}]$ "
                 r"at $20\%$; "
                 r"the generated table is \texttt{tab\_D\_loadprec}, which "
                 r"Section~\ref{sec:spec} prints and the artefact check regenerates. An "
                 r"earlier draft of this row quoted $0.1566$ to $0.2637$ and a gate of "
                 r"$12\%$, from a sweep that moved the accuracy in only one of the two "
                 r"places)",
        kind=LUMPED,
        rests_on=r"deliberately not assigned: nothing here bounds how well the load is recorded, so it is carried as a family indexed by the row below"),
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
        kind=MEASURABLE, species="human (MMA); mouse",
        rests_on=r"chosen for illustration. The nearest published repeatability, the breath test's, is between-visit and in methylmalonic acidemia~\cite{breath_test_mma}, so it bounds this from above rather than measuring it"),
    Row(r"$s_{\ln k_e}$", f"${KE_PRIOR_LOG:.1f}$",
        "sd of the log-normal prior on $k_e$",
        ASSUMED, r"\textbf{The only prior that does any work.} "
                 r"A factor of $e^{0.8} = 2.23$ per standard deviation, which is "
                 r"\emph{tighter} than the factor of $\num{ke_needed_factor}$ the gate "
                 r"needs, which is why "
                 r"Section~\ref{sec:designsec} can report that the gate passes: at the width "
                 r"the requirement names it would be met exactly and no more. It lives as a "
                 r"default argument in \texttt{src/design/run.py} and is asserted against it "
                 r"in the tests. What it is grounded in is the open question, since by the "
                 r"$k_e$ row above no source supplies the constant it is a prior on. What it "
                 r"is worth: on the densest design Section~\ref{sec:designsec} lists, "
                 r"sampling every two days for two years, "
                 r"$\mathrm{sd}(\ln\hat\theta)$ is $\num{sd_dense}$ with the prior and "
                 r"$\num{sd_dense_flat}$ "
                 r"without, medians over $24$ simulated intake records with spreads "
                 r"$[\num{sd_dense_lo}, \num{sd_dense_hi}]$ and $[\num{sd_dense_flat_lo}, "
                 r"\num{sd_dense_flat_hi}]$, a factor of $\num{ke_prior_worth}$ "
                 r"(\texttt{audit/n3\_ke\_prior\_worth.py})",
        kind=PRIOR,
        rests_on=r"chosen for illustration, tighter than the factor of $\num{ke_needed_factor}$ the gate needs, so the gate passes with it in place"),
    Row(r"$\sigma_L/L$", f"${100*LOAD_FRAC:.0f}\\%$",
        "relative error in the recorded substrate load",
        ASSUMED, r"\textbf{Sets $\sigma$} through "
                 r"$\sigma_{\mathrm{total}}$ and the load-bias prior through "
                 r"$s_b = \sigma_L/L \cdot L$. It is the recording error, distinct from the "
                 r"day-to-day variation of the load itself, and the two were conflated in "
                 r"Section~\ref{sec:spec} until R17. Measurable in principle by a weighed "
                 r"intake study; no such study exists in PA, where records are typically "
                 r"prescriptions rather than intake as eaten",
        kind=MEASURABLE,
        rests_on=r"chosen for illustration. Measurable by a weighed intake study; none exists in propionic acidemia, where records are prescriptions rather than intake as eaten"),
]


PARENT = [
    Row(r"$\theta$", rf"${_C.theta:.0f}$~$\mu$M/d", "residual capacity, parent model",
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
    Row(r"$K_m$", rf"${_C.Km:.0f}$~$\mu$M", "half-saturation of the enzymatic channel",
        ASSUMED, r"as above, and shared with the local model so the two can be compared. "
                 r"Sharing it across a tenfold difference in capacity assumes $K_m$ does not "
                 r"depend on residual activity, and \textbf{that assumption is supported}: "
                 r"four variant PCCs show no significant difference in $K_m$ for "
                 r"propionyl-CoA against wild type, although they carry about half the "
                 r"$k_{\mathrm{cat}}/K_m$~\cite{pcc_variants_kinetics}. The cost bound on the "
                 r"local $K_m$ row applies here too",
        kind=EFFECTIVE, species="human, in vitro"),
    Row(r"$K_c$", rf"${_C.Kc:.0f}$~$\mu$M", "half-saturation of the conjugation step",
        ASSUMED, r"no measurement, and no reason a conjugation reaction should share an "
                 r"enzyme's half-saturation. Cost bounded in supplement \S S1: over $K_c$ "
                 r"from $25$ to $800$, a $32$-fold range, $\thc$ spans $186$ to $120$, a "
                 r"factor of $1.56$. The nearest measured analogue is the renal carnitine "
                 r"transporter, $K_m$ $0.08$ mM external and $1.2$ mM "
                 r"internal~\cite{octn2_renal_transport}, which is a different step and is "
                 r"\textbf{a rat preparation measured in vitro}; no half-saturation for the "
                 r"conjugation itself is published",
        kind=EFFECTIVE, species="rat, in vitro"),
    Row(r"$M_i$", rf"${_C.Mi:.0f}$~$\mu$M",
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
    Row(r"$R_0$", rf"${_C.R0:.1f}$~$\mu$M/d", "carnitine replenishment",
        DERIVED, r"\textbf{not free}. Fixed by requiring the alternative route alone to peak "
                 rf"at $\mathrm{{ALT}}_{{\mathrm{{peak}}}} = {ALT_PEAK_C:.0f}$ against a "
                 rf"baseline load of ${S_BASE_C:.0f}$, so that a patient with $\theta = 0$ "
                 r"cannot be maintained. Changing $K_c$, $M_i$ or $\rho_0$ moves that peak "
                 r"and $R_0$ must be re-derived, a failure the review log documents in three "
                 r"places. \textbf{Changing the LOAD requires it too}, which is the case the "
                 r"row used to leave out: the constraint names a baseline load, so evaluating "
                 r"the parent at the local model's $L$ while holding $R_0$ makes a null "
                 r"genotype maintainable and $\theta_{\mathrm{crit}}$ undefined. That is not "
                 r"a technicality. \textbf{Any comparison of the two models at disease-range "
                 r"capacity turns on the convention chosen}: holding $R_0$ puts the parent's "
                 r"patient below the local model's level, re-deriving $R_0$ at the new load "
                 r"puts it above, and matching the alternative route to the local one at the "
                 r"operating point leaves several capacities with no equilibrium at all. The "
                 r"direction of the disagreement is a choice and not a result, so no such "
                 r"comparison is reported anywhere in this paper. What survives every "
                 r"convention is that the parent's clearance peak lies above "
                 r"$\sqrt{K_c M_i}$, which is $R_0$-free. \textbf{The stored constant carries eight significant figures from "
                 r"two round inputs}; the digits past the fourth are an artefact of the "
                 r"derivation and carry no precision. It is left unrounded because every "
                 r"parent-model figure was produced with it",
        kind="derived"),
    Row(r"$\mathrm{ALT}_{\mathrm{peak}}$", rf"${ALT_PEAK_C:.0f}$~$\mu$M/d",
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
    Row(r"$s_{\mathrm{base}}$", rf"${S_BASE_C:.0f}$~$\mu$M/d", "baseline propiogenic load",
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
    Row(r"$k$", r"$1$~d$^{-1}$",
        "conjugation rate constant, which cancels from the reduced curve",
        ASSUMED, r"drops out of the quasi-steady reduction of supplement \S S1 exactly, so the "
                 r"reduced clearance does not depend on it. \textbf{It is not therefore a "
                 r"free scale, which is what this row concluded from that, and which is "
                 r"why it was previously classed as derived rather than assumed}: $k$ "
                 r"alone sets the carnitine "
                 r"relaxation rate, $\tau_c = 1/(k[u(m)+\rho(m)])$, so it fixes the "
                 r"reduction's validity \emph{off} equilibrium even though it cancels from "
                 r"the curve, and nothing in this model fixes $k$ itself. \textbf{Cost}: "
                 r"stepping the substrate load to $\num{illness_fold}$ times baseline, "
                 r"which is the registered illness amplitude of the row above, with "
                 r"carnitine held at its pre-perturbation value, the peak relative error in "
                 r"$m$ against the reduced system is $\num{reduction_cost}\%$ at the "
                 r"assumed $k = 1$, falling as $1/k$ for larger $k$ "
                 r"($\num{reduction_cost_k10}\%$ at $k=10$, "
                 r"$\num{reduction_cost_k100}\%$ at $k=100$) and saturating below it "
                 r"($\num{reduction_cost_k0p1}\%$ at $k=0.1$), because a relative error "
                 r"cannot exceed order one. \textbf{The perturbation is named because "
                 r"without it the bound is not one}: an earlier version of this row gave "
                 r"four percentages and no perturbation, and a $10\%$ load step reproduces "
                 r"none of them, returning $\num{reduction_cost_step10}\%$ where this "
                 r"returns $\num{reduction_cost}\%$. \textbf{The cost is a function of "
                 r"the perturbation and not a constant}, so a perturbation this paper has "
                 r"not yet used has to be re-read rather than bounded by these: an "
                 r"isotope-dilution measurement of the load is itself a load perturbation "
                 r"and a small one by construction, and at $\num{tracer_fold}$ times "
                 r"baseline the reduction costs $\num{reduction_cost_tracer}\%$ at "
                 r"$k = 1$. Everything this paper computes from the "
                 r"parent is computed \emph{at} equilibrium, where the reduction is an "
                 r"identity and the cost is zero. The bound above is what leaving "
                 r"equilibrium costs, and \textbf{Section~\ref{sec:designsec}'s transient "
                 r"leaves it}, so it is load-bearing there as well as being the third "
                 r"reason Section~\ref{sec:discussion} excludes the decompensation route",
        kind=ELIMINATED, species="human",
        rests_on=r"cancels from the reduced clearance exactly, and nothing fixes it. Chosen for illustration; its cost off equilibrium is bounded in \S S2"),
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

#: the rows the MAIN TEXT reasons about, block by block. A constant earns its
#: place here by having its value quoted or depended on in the paper itself,
#: not by being registered: the operating point and the ratio the excitation
#: condition is stated against, the capacity rows, the two channels, the
#: dosing schedule the aliasing rule turns on, the noise the fit carries, the
#: one prior that does any work, and the load record. The parent block
#: contributes only `k`, which two passages of Section 2 cite by name for its
#: cost bound. Everything else is reference material and lives in the
#: supplement's full registry, which remains the authority: this is a
#: selection from it, not a second register.
MAIN_TEXT_ROWS = {
    "LOCAL": (r"$\bar m$", r"$\theta$", r"$\theta_{\mathrm{PA}}$", r"$K_m$",
              r"$k_e$", r"$\lambda$", r"$T_{\mathrm{dose}}$", r"$D$",
              r"$\sigma$", r"$\mathrm{CV}$", r"$s_{\ln k_e}$",
              r"$\sigma_L/L$"),
    "PARENT": (r"$k$",),
}


def _tex(rows, compact=False):
    """
    One LaTeX row per registry row.

    `compact` drops the provenance note from every row except those named in
    `keep_notes`, and points the rest at the supplement. It does NOT drop the
    kind or species tags: those are one word each, and they are what stop a
    reader taking an assumed lumped murine constant for a measured human one.

    """
    out = []
    for r in rows:
        #  The main-text table's last column used to be an abridged cost
        #  column whose cells were mostly a pointer at the supplement, which
        #  told a reader nothing about the row they were reading. It now says
        #  what the value rests on: a citation, or what the value was chosen
        #  against, or that the literature was searched and supplies nothing.
        note = r.rests_on if compact else r.note
        out.append(" & ".join([r.symbol, r.value, r.meaning,
                               _source_cell(r), note]) + r" \\")
    return out


def write_params_table(compact=False, name=None):
    r"""Render both blocks as one longtable, source column included.

    A longtable rather than a tabular in a float: once every row carries its
    literature footing the table is taller than a page, and a tabular inside a
    table environment does not break -- it silently runs off the bottom and
    takes the rest of the table with it. The caption and the label are left to
    the including document, which defines \paramscaption and \paramslabel,
    because the paper and the supplement caption this table differently and
    the supplement imports the paper's labels with xr, so a hard-coded label
    would collide.

    `compact` writes the main-text version: the rows of MAIN_TEXT_ROWS, with
    every tag, and a note only on the three rows of MAIN_TEXT_NOTES. The full
    registry goes in the supplement and stays the authority; the main-text
    table is a selection from it and says so in its own caption.
    """
    from .style import TABDIR, ROOT
    last = (r"\textbf{What the value rests on}" if compact
            else r"\textbf{What fixes it, or what assuming it costs}")
    #  \linewidth and not \textwidth: both tables are set inside an
    #  adjustwidth, so \textwidth is the wrong measure and the sum of the
    #  columns plus \tabcolsep ran past the rule by about 10pt
    #  the widths have to leave room for \tabcolsep on both sides of five
    #  columns, which is another 72pt on top of whatever they sum to; at
    #  0.89 of the line the table ran 46pt past the rule
    widths = ((r"p{0.075\linewidth} p{0.085\linewidth} p{0.22\linewidth} "
               r"p{0.11\linewidth} p{0.31\linewidth}") if compact else
              #  the symbol and value columns carry unbreakable math --
              #  $s_1 = s_2$, $3000$ $\mu$M/d -- which at 0.06 and 0.07 ran
              #  past their own column and took the whole row with them
              (r"p{0.08\linewidth} p{0.10\linewidth} p{0.15\linewidth} "
               r"p{0.085\linewidth} p{0.335\linewidth}"))
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
    local = ([r for r in LOCAL if r.symbol in MAIN_TEXT_ROWS["LOCAL"]]
             if compact else LOCAL)
    parent = ([r for r in PARENT if r.symbol in MAIN_TEXT_ROWS["PARENT"]]
              if compact else PARENT)
    body += _tex(local, compact)
    body += [r"\addlinespace", r"\midrule",
             r"\multicolumn{5}{@{}l}{\emph{The parent model and its equilibrium reduction}} \\",
             r"\addlinespace[2pt]"]
    body += _tex(parent, compact)
    if compact:
        body += [r"\addlinespace",
                 r"\multicolumn{5}{@{}p{0.92\textwidth}}{\emph{The remaining "
                 + str(len(LOCAL) + len(PARENT) - len(local) - len(parent))
                 + r" registered constants are in the supplement, \S S2, with "
                 r"the full provenance of every row in this table.}} \\"]
    body += [r"\end{longtable}"]
    path = TABDIR / (name or ("tab_params_compact.tex" if compact
                              else "tab_params.tex"))
    path.write_text("\n".join(body))
    #  a redirected write (tools/check_artefacts.py) lands outside the repo,
    #  where relative_to would raise rather than report
    shown = path.relative_to(ROOT) if path.is_relative_to(ROOT) else path
    print(f"  wrote {shown}")
    return path


if __name__ == "__main__":
    write_params_table()
    write_params_table(compact=True)
