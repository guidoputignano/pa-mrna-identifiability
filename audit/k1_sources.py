"""
Phase 5, step 0: the five source checks, and what two of them do to the
design numbers.

Nothing here records an expected result in advance. Every figure printed is
either read off a primary source in this pass (and named as such) or computed
from the repository's own code with one constant changed.

The two constants under test are `design.run.LAMBDA` and `design.run.D_AMP`.
Both are ASSUMED rows in src/common/params.py and both are corroborated in
their annotations against the same secondary reading. This script recomputes
the design numbers with each replaced, and with both replaced together.

Monkeypatching note: `design.run.delta` reads LAMBDA and D_AMP out of its own
module globals at call time, so rebinding the module attribute is enough and
no caller has to be touched. `g_design_common` imported the FUNCTION, not the
constants, so it follows along. The `_amp` and `_lam` context managers below
restore the originals so that a later block cannot inherit an earlier one's
setting -- an error worth guarding, because every number here is a ratio to
a number computed under a different setting.
"""
import sys
from contextlib import contextmanager
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

import design.run as DR
import g_design_common as G

TH, KM, KE, LB = G.TH, G.KM, G.KE, G.LB
T = G.T

# --------------------------------------------------------------------------
#  Primary-source values, read off the PDFs in this pass.
#  Attarwala H, Lumley M, Liang M, Ivaturi V, Senn J. Nucleic Acid Ther
#  2023;33(2):141-147, doi:10.1089/nat.2022.0036, Table 1 and Suppl Table S1.
# --------------------------------------------------------------------------
KDEG = 0.0075        # /h, PCC protein degradation, Table 1 (2.82% RSE)
KQ = 0.00474         # /h, PCC intercompartmental rate, Table 1 (6.7% RSE)
KE0 = 0.242          # /h, effect-compartment equilibration, Table 1 (3.35% RSE)
SLOPE = 55.6         # (mass/g/h)/(mass/mL), Table 1 (1.93% RSE)
IMAX = 0.999         # fixed, Table 1 footnote m
HP_BASE = 60.1       # "mM/L", Table 1 (14.4% RSE), IIV 51.1%
HP_E0 = 0.0987       # "mM/L", Table 1 (6.74% RSE)
HP_IC50 = 37.5       # "mg/g", Table 1 (12.2% RSE), IIV 73.8%
HP_LOQ = 25.0        # "mM/L", assay limit of quantification, text p.145
PAPER_T_HALF_D = 7.0  # "The model estimated a half-life of PCC protein of ~7 days."

# Extended Data Table 3 of Moderna's phase 1/2 interim analysis
# (Nature 2024;628:872-877, doi:10.1038/s41586-024-07266-7): 3-HP mean
# percentage change from baseline after the third dose, by cohort, with SD.
EDT3_3HP = {
    "1: 0.30 mg/kg Q3W (n=4)": (-29.0, 25.0),
    "2: 0.30 mg/kg Q2W (n=3)": (-3.0, 26.0),
    "3: 0.45 mg/kg Q2W (n=3)": (-20.0, 24.0),
    "4: 0.60 mg/kg Q2W (n=3)": (5.0, 57.0),
    "5: 0.90 mg/kg Q2W (n=3)": (-20.0, 41.0),
}

DBS_3HP_CUTOFF = 38.13   # umol/L, 99.9th centile, repository's dbs_3hp_cutoff


# --------------------------------------------------------------------------
def check_2_species():
    """Finding 2: which repository rows present mouse estimates as human."""
    print("=" * 74)
    print("2. SPECIES: the Table 1 PD estimates are MOUSE")
    print("=" * 74)
    print("  What the primary source says, read in this pass:\n")
    print("    - 'in PCC-deficient mice [(Strain Code: 4011020 (A138T)")
    print("       PCCA -/-); N = 111]'")
    print("    - the PD section is headed 'PK/PD model for plasma 2-MC, 3-HP,")
    print("       and C3/C2 ratio in PCC-deficient mice'")
    print("    - 'plasma levels of corresponding biomarkers that are")
    print("       susceptible to suppression by PCC protein expressed through")
    print("       mRNA-3927 in mice'")
    print("    - Suppl Table S1(A) scales on (BW/0.025)**cl, i.e. a 25 g")
    print("       reference animal")
    print("    - the ONLY scaling statement is: 'The preclinical PK/PD model")
    print("       described earlier was extrapolated to humans based on")
    print("       allometric scaling of model parameters. Based on the PCCA/B")
    print("       mRNA PK model described earlier, estimated allometric")
    print("       exponent was 0.63 for CL12 and CL32 and 1.1 for CL23 and")
    print("       CL20. Allometric exponent for the volume of distribution")
    print("       parameters (ie, V and V2) was fixed to 1.'\n")
    print("  Every parameter named in that sentence is a PK parameter. No PD")
    print("  parameter is scaled, and the companion paper says the same: 'Each")
    print("  preclinical PK/PD model was extrapolated to humans based on")
    print("  allometric scaling of PK parameters'.\n")
    print("  So E0, base, IC50, Imax, Slope, kdeg, kq and ke0 are MOUSE")
    print("  estimates, and the between-subject variabilities are between MICE.")
    print("  Supporting: the mouse model carries ~2% of normal liver PCC")
    print("  activity, which is a different residual capacity from a patient's.\n")
    print("  WHERE THE REPOSITORY PRESENTS THEM AS HUMAN")
    print("  (results/tables/tab_params.tex is generated from params.py by")
    print("  `make figures`, so each row below is one edit, not two):\n")
    for loc, quote, why in [
        ("src/common/params.py:60-63  -> tab_params.tex:19  (m-bar row)",
         "'the translational model puts the untreated plasma 3-HP baseline "
         "at 60.1 uM and the fully corrected floor near 0.1 uM, so a managed "
         "PATIENT sits between the two'",
         "reads a mouse baseline and a mouse floor as the scale a patient's "
         "marker is measured on. The floor 0.1 is E0 = 0.0987, also mouse."),
        ("src/common/params.py:107-113 -> tab_params.tex:26  (lambda row)",
         "'Corroborated: the translational PK/PD model estimates a "
         "PCC-protein degradation rate of 0.0075 /h, a half-life of 3.85 d, "
         "so this is within 10% of a published value'",
         "mouse, and see finding 3: 3.85 d is a rate constant inside the "
         "model, not the model's protein half-life, which the paper gives "
         "as ~7 d."),
        ("src/common/params.py:135-136 -> tab_params.tex:30  (CV row)",
         "'Biological variability is far larger, the translational model "
         "carrying 51% between subjects on the 3-HP baseline'",
         "51.1% IIV is between MICE (N = 111 of one inbred strain), so it "
         "is if anything a LOWER bound on between-patient variability, and "
         "the sentence uses it as the large figure."),
        ("src/common/params.py:82-83  -> tab_params.tex (k_e row)",
         "'The nearest route is the first-order 3-HP elimination implicit "
         "in the published indirect-response model'",
         "this one is not a species error but a structure error. Suppl "
         "Table S1(C) is ALGEBRAIC in the marker: 3-HP = E0 + base*(1 - "
         "Imax*PCC/(IC50+PCC)). There is no 3-HP turnover equation and no "
         "elimination constant anywhere in that model, so nothing can be "
         "extracted from it. claims.csv C46 already records this ('no "
         "biomarker turnover constant at all'). The model that DOES carry "
         "a first-order marker loss is the OTHER one, the inhibition-of-kin "
         "turnover model; it publishes no value for it."),
        ("paper/refs.sources:102-116  (the pkpd_translational annotation)",
         "'Semi-mechanistic translational PK/PD model for mRNA-3927, fitted "
         "largely to preclinical data and allometrically scaled to human' "
         "and the PROVENANCE note",
         "'allometrically scaled to human' is true of the PK and false of "
         "the PD, and the six figures beneath it are all PD. The provenance "
         "note ('the numbers await a check against the primary source') can "
         "be lifted: all six are confirmed against Table 1 in this pass."),
    ]:
        print(f"    {loc}")
        print(f"      says : {quote}")
        print(f"      issue: {why}\n")
    print("  Proposed corrections are listed in audit/PHASE5_ROBUSTNESS.md and")
    print("  are NOT applied here. No constant moves: 60.1, 0.0075 and 51% are")
    print("  all annotation text, not values the code uses. M_BAR stays 10.\n")


def check_5_model_identity():
    """Finding 5: what each simulator is, and what it must be called."""
    print("=" * 74)
    print("5. WHAT MODELS A AND B ARE")
    print("=" * 74)
    print("  MODEL A -- the turnover structure.")
    print("  Named, with no equations and no numbers, in two places:")
    print("    Gao (ACoP2025 abstract, dose_selection_poster), quoted in")
    print("      claims.csv C46: 'An indirect inhibition Kin model was used to")
    print("      describe the PD response of 3-HP and 2-MC.'")
    print("    Baek et al. 2024 (Nat Commun 15:3804), Methods, read in this")
    print("      pass: 'For PD biomarkers, classic indirect response models")
    print("      with inhibition of kin were developed for 3-HP and 2-MC/MMA,")
    print("      incorporating PCC/hMUT enzyme-induced decreases in organic")
    print("      acid accumulation.'")
    print("  Baek's Table 1 is a DOSE-PROJECTION table (LNP, starting dose,")
    print("  species, NOAEL, safety margin). It carries no PD parameter values")
    print("  at all, so it cannot supply numbers either.")
    print("  => Model A is a named structure with NO published parameters.")
    print("     Every number in our Model A is calibrated to OUR operating")
    print("     point and must be flagged as ours in every table.\n")
    print("  REGISTRATION: Baek et al. 2024 is ALREADY a source in this")
    print("  repository, as `mrna_moa` (doi:10.1038/s41467-024-47460-9,")
    print("  refs.sources:200-203). No new citekey is needed. What is missing")
    print("  from its annotation is the sentence above, which is the only")
    print("  place in the peer-reviewed literature where Model A's structure")
    print("  is stated; the repository currently attributes that structure to")
    print("  the conference abstract alone.\n")
    print("  MODEL B -- Attarwala's structure, Suppl Table S1(B) and (C):")
    print("      d(Ce)/dt  = ke0*(C13 - Ce)")
    print("      d(PCC)/dt = Ksyn + kq*(PCP - PCC) - kdeg*PCC,  Ksyn = Slope*Ce")
    print("      d(PCP)/dt = kq*(PCC - PCP)")
    print("      3-HP      = E0 + base*(1 - Imax*PCC/(IC50 + PCC))")
    print("  with the Table 1 values above. These are MOUSE values (finding 2)")
    print("  and the paper rescales NOTHING in this block to human, so our")
    print("  Model B carries them unscaled and says so.\n")
    print("  The structural point, which is what step 2 will test: Model B's")
    print("  marker is ALGEBRAIC in the enzyme. It has no marker state, so it")
    print("  cannot express a substrate load at all. Where L enters is OUR")
    print("  choice in both simulators, and it is an assumption, not a")
    print("  transcription. Step 1 states it explicitly.\n")


# --------------------------------------------------------------------------
@contextmanager
def _set(**kw):
    """Temporarily rebind design.run module constants."""
    old = {k: getattr(DR, k) for k in kw}
    for k, v in kw.items():
        setattr(DR, k, v)
    try:
        yield
    finally:
        for k, v in old.items():
            setattr(DR, k, v)


def steady_level(inc):
    """The quasi-steady marker at a given therapeutic increment."""
    B = KE * KM + TH + inc - LB
    return (-B + np.sqrt(B ** 2 + 4.0 * KE * LB * KM)) / (2.0 * KE)


def fall(inc):
    """Fractional fall of the marker below the untreated operating point."""
    return 1.0 - steady_level(inc) / steady_level(0.0)


# --------------------------------------------------------------------------
def check_1_units_and_form():
    """Finding 1: the transfer-function reading."""
    print("=" * 74)
    print("1. THE TRANSFER FUNCTION: units and functional form")
    print("=" * 74)
    print("  Read off Table 1 and Supplementary Table S1 in this pass.\n")

    print("  (a) The slope is LINEAR, not a Hill exponent.")
    print("      Suppl Table S1(B):  Ksyn = Ce * Slope")
    print("      Text p.144: 'the synthesis rate of PCC protein increased")
    print("      linearly with an increase in plasma PCC mRNA concentration")
    print("      through an effect compartment'.")
    print(f"      Slope = {SLOPE} with RSE 1.93%. There is no exponent in the")
    print("      PCC equation at all; the only saturating term in the whole")
    print("      model is the Imax expression on the MARKER, in part (C).\n")

    print("  (b) Imax is fixed, not estimated.")
    print(f"      Table 1 footnote m: 'Maximum inhibition (Imax) was fixed at {IMAX}.'")
    print("      So the marker's floor is E0, not a fitted maximum effect:")
    print(f"      3-HP -> E0 + base*(1 - {IMAX}) = {HP_E0:.4f} + "
          f"{HP_BASE * (1 - IMAX):.4f} = {HP_E0 + HP_BASE * (1 - IMAX):.4f}")
    print("      as PCC grows without bound.\n")

    print("  (c) The slope's unit as PRINTED is mg, not ug.")
    print("      Table 1 prints  Slope [(mg/g/h)/(mg/mL)].  The repository's")
    print("      refs.sources gives the IC50 in ug/g. The RATIO in the slope")
    print("      is scale-free if both masses carry the same prefix, so 55.6")
    print("      is unchanged either way -- but the printed unit is mg and the")
    print("      annotation should quote it as printed.\n")

    print("  (d) 'mg/g' in the IC50 row is a misprint for ug/g, on the")
    print("      paper's OWN evidence, not on ours:")
    print("      - Fig. 2B axis label, read in this pass: 'PCC protein, ug/g',")
    print("        with the plotted liver PCC running to about 700 on that axis.")
    print(f"      - The IC50s are printed 21.0 / {HP_IC50} / 32.1 'mg/g'. On the")
    print("        figure's own axis those sit inside the plotted range; read")
    print("        as mg/g they would be a THOUSANDFOLD above every point the")
    print("        model was fitted to, so no dose could reach half effect.")
    print("      Verdict: ug/g. The repository already has this right.\n")

    print("  (e) 'mM/L' in the marker rows is a misprint for umol/L.")
    print(f"      - 3-HP baseline printed {HP_BASE} 'mM/L'. As millimolar that")
    print("        is a plasma organic acid at 60 mM.")
    print(f"      - The same text gives the assay LOQ as {HP_LOQ} 'mM/L'. An LOQ at")
    print(f"        {HP_LOQ / HP_BASE:.2f} of the disease baseline is an ordinary")
    print("        assay; an LOQ at 25 mM is not an assay at all.")
    print(f"      - The repository's own population anchor, the newborn")
    print(f"        dried-blood-spot 99.9th centile, is {DBS_3HP_CUTOFF} umol/L,")
    print(f"        a factor {HP_BASE / DBS_3HP_CUTOFF:.2f} from {HP_BASE}.")
    print("      Verdict: umol/L. The repository already has this right.")
    print("      'mM/L' is in any case not a unit: mM is already per litre.\n")


def check_3_half_life():
    """Finding 3: 3.85 d or 7 d?"""
    print("=" * 74)
    print("3. THE PCC HALF-LIFE: which one corroborates lambda")
    print("=" * 74)
    A = np.array([[-(KQ + KDEG), KQ],
                  [KQ, -KQ]])
    w = np.linalg.eigvals(A)
    w = np.sort(w.real)[::-1]          # w[0] is the slow (terminal) root
    t_half_terminal = np.log(2.0) / (-w[0]) / 24.0
    t_half_fast = np.log(2.0) / (-w[1]) / 24.0
    t_half_kdeg = np.log(2.0) / KDEG / 24.0

    print("  Suppl Table S1(B) is TWO compartments for the protein:")
    print("      d(PCC)/dt = Ksyn + Kq*(PCP-PCC) - Kdeg*PCC")
    print("      d(PCP)/dt = Kq*(PCC-PCP)")
    print(f"  with Kdeg = {KDEG} /h and Kq = {KQ} /h, both from Table 1.\n")
    print("  (A note on S1 as printed: it gives 'd(PCC)/dt = Ksyn*Ce + ...'")
    print("  with 'Ksyn = Ce*Slope' on the next line, which makes the term")
    print("  quadratic in Ce and its units wrong. The text says the synthesis")
    print("  rate rises LINEARLY in mRNA, so the intended term is Ksyn alone.")
    print("  The eigenvalues below do not depend on which reading is taken:")
    print("  the input enters the forcing, not the matrix.)\n")
    # the model's own decline: a unit bolus into the effect compartment,
    # integrated, then the time from the PCC peak down to half of it.
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq

    def rhs(t, y):
        ce, pcc, pcp = y
        return [-KE0 * ce,
                SLOPE * ce + KQ * (pcp - pcc) - KDEG * pcc,
                KQ * (pcc - pcp)]

    sol = solve_ivp(rhs, (0.0, 4000.0), [1.0, 0.0, 0.0],
                    rtol=1e-10, atol=1e-14, dense_output=True)
    tt = np.linspace(0.0, 4000.0, 400001)
    pcc = sol.sol(tt)[1]
    ipk = int(np.argmax(pcc))
    tpk, ypk = tt[ipk], pcc[ipk]
    t_half_decay = brentq(lambda t: sol.sol(t)[1] - 0.5 * ypk,
                          tpk, 4000.0) - tpk

    # the LOCAL half-life along the decline, -ln2 / d(ln PCC)/dt. The profile
    # is biexponential, so this is not one number: it starts near the fast
    # phase and rises to the slow one. Where it passes the paper's ~7 d is
    # the checkable statement.
    post = tt >= tpk
    lg = np.log(np.clip(pcc[post], 1e-300, None))
    sl = np.gradient(lg, tt[post])
    local = -np.log(2.0) / sl / 24.0
    t_post = (tt[post] - tpk) / 24.0
    # near the peak the slope passes through zero and the local half-life
    # blows up, so the crossing is taken as the FIRST time after the profile
    # is genuinely declining (below 90% of peak) that it rises through 7 d.
    ok = np.isfinite(local) & (t_post > 0.0) & (pcc[post] < 0.9 * ypk)
    t_ok, l_ok = t_post[ok], local[ok]
    def _first_at(target):
        j = np.nonzero(l_ok >= target)[0]
        return float(t_ok[j[0]]) if j.size else float("nan")

    cross = _first_at(PAPER_T_HALF_D)
    cross_lam = _first_at(3.5)
    loc14 = float(np.interp(14.0, t_ok, l_ok))

    print(f"  ln2/Kdeg alone                  = {t_half_kdeg:.2f} d")
    print(f"  fast eigenvalue                 = {t_half_fast:.2f} d")
    print(f"  slow (terminal) eigenvalue      = {t_half_terminal:.2f} d")
    print(f"  peak-to-half of the simulated")
    print(f"    single-dose PCC profile       = {t_half_decay/24.0:.2f} d "
          f"(peak at {tpk/24.0:.2f} d)")
    print(f"  LOCAL half-life one dose")
    print(f"    interval (14 d) after the peak= {loc14:.2f} d")
    print(f"  the paper's own statement       = ~{PAPER_T_HALF_D:.0f} d\n")
    print("  NO single constant in Table 1 reproduces the paper's ~7 d. The")
    print("  simulated protein declines BIEXPONENTIALLY, which is what the")
    print("  text says the second compartment was added for, so its LOCAL")
    print("  half-life is a function of time, rising from the fast phase to")
    print("  the slow one:")
    for d in (1.0, 2.0, 5.0, 10.0, 14.0, 21.0, 40.0):
        print(f"      {d:5.0f} d after the peak   local half-life "
              f"{float(np.interp(d, t_ok, l_ok)):6.2f} d")
    print(f"  It passes the paper's {PAPER_T_HALF_D:.0f} d at {cross:.1f} d after the peak and")
    print(f"  tends to {t_half_terminal:.2f} d. So the ~7 d is a reading from the middle")
    print("  of that decline; it is consistent with the constants and is not")
    print("  an arithmetic disagreement with them. (The time from the peak to")
    print(f"  half the peak is a different quantity, {t_half_decay/24.0:.2f} d, because most")
    print("  of that drop happens in the fast phase.)\n")
    print(f"  What this does to lambda: params.py quotes {t_half_kdeg:.2f} d, calls it")
    print("  a 'half-life of PCC protein', and calls lambda's 3.5 d")
    print("  CORROBORATED as being within 10% of it.")
    print(f"  The comparison is defensible EARLY and not late. The local half-")
    print(f"  life passes 3.5 d at {cross_lam:.1f} d after the peak, so over roughly the")
    print("  first week lambda tracks the published model closely. Over the")
    print("  window lambda actually has to govern -- a whole 14 d dose")
    print(f"  interval -- the model's local half-life is {loc14:.2f} d, a ratio of")
    print(f"  {loc14/3.5:.1f}x, and the terminal figure is {t_half_terminal:.2f} d, a ratio of")
    print(f"  {t_half_terminal/3.5:.1f}x. A single-exponential Delta with a 3.5 d half-life")
    print("  therefore decays FASTER late in the interval than the published")
    print("  model's protein does, which is the direction that makes the")
    print("  excursion look larger than it is.")
    print(f"  What the annotation cannot support as written is that {t_half_kdeg:.2f} d is")
    print("  'a half-life of PCC protein'. It is one rate constant inside a")
    print("  two-compartment protein model whose own stated protein half-life")
    print(f"  is ~{PAPER_T_HALF_D:.0f} d.\n")
    return t_half_terminal, t_half_decay / 24.0, loc14


def check_4_calibration():
    """Finding 4: what fall does the trial actually report, and what D gives it."""
    print("=" * 74)
    print("4. THE DOSE SIZE: what fall the trial reports")
    print("=" * 74)
    print("  Nature 2024 Extended Data Table 3, read off the rendered page in")
    print("  this pass. 3-HP, mean % change from baseline after the third dose:\n")
    vals = []
    for k, (mu, sd) in EDT3_3HP.items():
        print(f"      cohort {k:26s}  {mu:+5.0f}%   (SD {sd:.0f})")
        vals.append(mu)
    v = np.array(vals)
    print(f"\n      n = {v.size} cohorts, mean of the cohort means {v.mean():+.1f}%,")
    print(f"      range {v.min():+.0f}% to {v.max():+.0f}%. Two of five cohorts fall by")
    print("      about 20-29%; one falls 3%; one RISES 5%.\n")
    print("  Against the repository's D_AMP = 6000, which params.py says was")
    print("  'chosen so the marker falls by about half at peak effect ... the")
    print("  order of the reductions reported in the phase 1/2 interim analysis':\n")
    with _set():
        peak_first = DR.D_AMP
        peak_steady = DR.D_AMP / (1.0 - np.exp(-DR.LAMBDA * T))
    trough = peak_steady * np.exp(-DR.LAMBDA * T)
    m_tr, m_pk = steady_level(trough), steady_level(peak_steady)
    print(f"      first-cycle peak   Delta = {peak_first:.0f}   fall {100*fall(peak_first):.2f}%"
          "   (from the untreated level)")
    print(f"      steady-state peak  Delta = {peak_steady:.0f}   fall {100*fall(peak_steady):.2f}%"
          "   (from the untreated level)")
    print(f"      steady trough->peak Delta {trough:.0f}->{peak_steady:.0f}   "
          f"fall {100*(1-m_pk/m_tr):.2f}%   (from the patient's own pre-dose)")
    print(f"      -> the three readings span "
          f"{100*min(fall(peak_first), fall(peak_steady), 1-m_pk/m_tr):.0f}-"
          f"{100*max(fall(peak_first), fall(peak_steady), 1-m_pk/m_tr):.0f}%, which is the "
          f"62-66% already\n         recorded unapplied in PHASE4_PROPOSED_TEXT.md.\n")
    print("  That is two to three times the reported falls, and the reported")
    print("  falls are not 'about half'.\n")
    print("  WHEN WAS 'AFTER THE THIRD DOSE' SAMPLED? The Nature paper does")
    print("  not say. Extended Data Table 3's caption is 'Biomarker percentage")
    print("  change from baseline after the third dose' with no time point;")
    print("  the Methods define the endpoint only as 'changes in blood 2-MC")
    print("  and 3-HP levels from baseline (predose levels) to post-dose")
    print("  levels'. The one sampling time the paper does name is a PREDOSE")
    print("  draw ('the predose sample for dose 3 (52 days after the first")
    print("  dose)'), which is a trough. Both readings are carried below.\n")
    return peak_first, peak_steady


def calibrate_D(target, reading, lam=None):
    """
    The D_AMP that makes the marker fall by `target` at the stated reading.

    reading = 'peak'   : sampled at the dose-3 peak,   Delta = D*sum_{j=0..2} q^j
    reading = 'trough' : sampled at the next predose,  Delta = D*sum_{j=1..3} q^j
    """
    lam = DR.LAMBDA if lam is None else lam
    q = np.exp(-lam * T)
    mult = {"peak": 1.0 + q + q ** 2, "trough": q + q ** 2 + q ** 3}[reading]
    # invert fall(inc) = target for inc, then divide by the multiplier
    lo, hi = 0.0, 1e7
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if fall(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi) / mult, mult


# --------------------------------------------------------------------------
def design_numbers(label, seeds=range(8)):
    """The three design quantities the reruns are reported on."""
    print(f"\n  --- {label} ---")
    print(f"      lambda = ln2/{np.log(2)/DR.LAMBDA:.2f} d,  D_AMP = {DR.D_AMP:.0f},  "
          f"steady peak Delta/theta = "
          f"{DR.D_AMP/(1-np.exp(-DR.LAMBDA*T))/TH:.2f}")
    print()
    print("      schedule                  Delta/theta range (n=12)   "
          "sd(ln theta), monthly n=12")
    print("                                                            "
          "CV 0     CV 2%              CV 10%")
    for name, fn in G.SCHEDULES.items():
        ts = fn(12)
        dr = G.delta_range(ts)
        row = f"      {name:24s}  {dr:22.4f}   "
        if name == "monthly (30.4 d)":
            outs = []
            for cv in (0.0, 0.02, 0.10):
                lo, md, hi = G.sd_over_seeds(ts, cv, seeds=seeds)
                outs.append(f"{md:.4f}" if cv == 0.0
                            else f"{md:.4f} [{lo:.4f},{hi:.4f}]")
            row += "  ".join(outs)
        print(row)

    print("\n      visits to sd(ln theta) <= 0.20, monthly record:")
    for cv in (0.0, 0.02, 0.10):
        hit = None
        for n in range(4, 61):
            _, md, _ = G.sd_over_seeds(G.SCHEDULES["monthly (30.4 d)"](n),
                                       cv, seeds=seeds)
            if md <= G.GATE:
                hit = n
                break
        print(f"          load CV {100*cv:4.0f}%   "
              f"{'n = ' + str(hit) if hit else 'not reached by n = 60'}")


def main():
    check_1_units_and_form()
    check_2_species()
    t_half_terminal, t_half_decay, loc14 = check_3_half_life()
    check_4_calibration()
    check_5_model_identity()

    print("=" * 74)
    print("THE RERUNS")
    print("=" * 74)
    print("  Every number below comes from g_design_common, unchanged, with")
    print("  design.run's LAMBDA and/or D_AMP rebound. Nothing in src/ is")
    print("  edited and no committed constant moves.\n")

    with _set():
        design_numbers("as committed (3.5 d half-life, D = 6000)")

    lam7 = np.log(2.0) / PAPER_T_HALF_D
    with _set(LAMBDA=lam7):
        design_numbers(f"finding 3: half-life {PAPER_T_HALF_D:.0f} d, D unchanged")

    print("\n  Finding 4 needs a D for each reading of 'after dose 3'. Solving")
    print("  fall(Delta at the sample) = target for D, at the committed lambda:\n")
    tab = {}
    for reading in ("peak", "trough"):
        for target in (0.20, 0.30):
            D, mult = calibrate_D(target, reading)
            tab[(reading, target)] = D
            print(f"      {reading:7s}  {100*target:.0f}% fall  ->  D = {D:8.1f}   "
                  f"(Delta at the sample = {D*mult:7.1f}, "
                  f"multiplier {mult:.4f})")

    for reading in ("peak", "trough"):
        for target in (0.20, 0.30):
            with _set(D_AMP=tab[(reading, target)]):
                design_numbers(f"finding 4: {100*target:.0f}% fall at the "
                               f"{reading}, half-life 3.5 d")

    print("\n  Both corrections together, at the reading that is hardest on the")
    print("  design (a 20% fall read at the peak) and at the softest of the")
    print("  four (a 30% fall read at the trough), with the 7 d half-life:\n")
    for reading, target in (("peak", 0.20), ("trough", 0.30)):
        D, mult = calibrate_D(target, reading, lam=lam7)
        with _set(LAMBDA=lam7, D_AMP=D):
            design_numbers(f"findings 3+4: {100*target:.0f}% fall at the "
                           f"{reading}, half-life {PAPER_T_HALF_D:.0f} d "
                           f"(D = {D:.1f})")

    print()
    print("=" * 74)
    print(f"(readings carried forward: terminal {t_half_terminal:.2f} d, "
          f"peak-to-half {t_half_decay:.2f} d, local at 14 d {loc14:.2f} d)")


if __name__ == "__main__":
    main()
