# Phase 7. Provenance and physical meaning of every row in Table 2

**Commit audited: `f71bf9af4442055b9d014b7dcff6aca5ff11c636`**, which is what `origin/main` and
local `main` both point at. `src/common/params.py` is byte-identical on that commit and on
this branch's head `5197106` (sha256 `f891996e985dd9339d09c7ad186340b7ed5fe9f2cfcf6c4b7b0cc825402d127c`
for both), because the branch's three commits touch no file under `src/`. There was no stale
copy. Nothing outside this file was changed and the paper was not touched.

The registry holds **20 rows: 16 `assumed`, 2 `literature` (`lambda`, `T_dose`), 2 `derived`
(`R_0`, `k`)**, counted from the dataclasses rather than from the caption.

## What this pass added to the record

Primary sources were read rather than cited. The decisive one is
`pkpd_translational` = Attarwala H, Lumley M, Liang M, Ivaturi V, Senn J,
*Nucleic Acid Ther* 2023;33(2):141-147, PMID 36577040, PMC10066765, read in full:

> "The translational PK/PD model was developed using **preclinical data available in mice with
> PA, Sprague Dawley rats, and cynomolgus monkeys** at dose levels ranging from 0.2 to 9 mg/kg
> ... The interspecies preclinical model was scaled allometrically to humans to **predict** the
> dose-response relationship."

> "PCC protein-mediated reduction of 2-MC, 3-HP, and C3/C2 ratios was modeled using a
> **direct sigmoidal maximum inhibition (Imax) model**"

> "Based on parameter estimates in Table 1, ~63.7%, 99.8%, and 86.9% of total plasma 2-MC,
> 3-HP, and C3/C2, respectively, was amenable to suppression ... **in PCC-deficient mice**"

Its Table 1, transcribed here from the full text: `k_e0` 0.242 /h (3.35% RSE), `slope` 55.6
(1.93%), **`k_deg` 0.0075 /h (2.82%)**, `k_q` 0.00474 /h (6.7%); and for 3-HP,
**`E0` 0.0987 uM/L (6.74%)**, **`base` 60.1 uM/L (14.4%) with IIV 51.1%**, `IC50` 37.5 ug/g
(12.2%, IIV 73.8%).

Three consequences, each used below:

1. **The source contains no human data.** Every number this repository takes from it is a
   preclinical estimate or an allometric prediction.
2. **It contains no marker turnover constant.** 3-HP is algebraic in PCC protein. The only
   `k_deg` in it is the **PCC protein** degradation rate.
3. **`refs.sources` has one structural mislabel.** It records "plasma 3-HP disease baseline
   60.1 umol/L". In the source, 60.1 is **`3-HP base`, the suppressible pool**; the untreated
   total is `E0 + base` = 60.2 uM. Numerically this is a 0.2% difference and harmless. The
   distinction is not harmless, because it is exactly the one `OPEN_QUESTIONS.md` decision (e)
   turns on, which already quotes Attarwala's definitions correctly. `refs.sources` is the file
   that is out of step, not `OPEN_QUESTIONS.md`.

**One good outcome.** `refs.sources:113-116` carries a PROVENANCE warning that its
`pkpd_translational` figures "are transcribed from the author's parameter-grounding review of
September 2026 ... the numbers await a check against the primary source." **That check is now
done and all seven numbers are correct as transcribed** (0.0075, 60.1, 0.0987, 51.1, 37.5,
55.6, 14.4 all appear in the full text). The warning can be discharged. What it did not catch,
and what matters more, is the species and the model structure.

---

## The table

`kind`: what sort of quantity it is. `absorbs`: what a lumped or effective quantity collapses
together. `depends`: whether a paper conclusion turns on the **value**, on the **order of
magnitude**, or on **neither**.

| # | Row | Value | Denotes | Kind | Absorbs | Source class declared / true | Species and context of source | Depends | Defect |
|---|---|---|---|---|---|---|---|---|---|
| L1 | `$\bar m$` | 10 | operating point the patient is held at | convention, and the unit anchor for the whole local model | the distribution volume, since `m` is a plasma concentration standing for a mitochondrial substrate | assumed / assumed | **mouse**: the 60.1 and 0.1 uM scale is Attarwala Table 1, PCC-deficient mice | order of magnitude (all design quantities are reported relative to it) | **D4, D5** |
| L2 | `$\theta$` | 3000 | residual capacity of the enzymatic channel | lumped | `V_max` of PCC, the functional hepatic PCC fraction, and the distribution volume (implied units umol/L/day, not enzyme units) | assumed / assumed | none. Residual activity is published only as a fraction of control | **value** | **D2** |
| L3 | `$K_m$` | 100 | half-saturation at the plasma level | effective | the whole mapping from plasma marker to mitochondrial substrate | assumed / assumed | human (`pcc_liver_kinetics`, *J Biol Chem* 1980;255:60-5, normal human liver), cited only to say the in vitro 0.29 mM is a **different** quantity | value (`m/K_m = 0.1` is what the excitation condition is stated against) | none. Well handled |
| L4 | `$k_e$` | 5 | alternative clearance rate constant | **lumped, and the catch-all** | every non-PCC route out of the marker pool. `tab:epistemic` names them, "carnitine conjugation and renal excretion"; this row does not | assumed / assumed | **the cited route does not exist in the cited source** | **value** (the design's binding prior) | **D1**, the most serious in the table |
| L5 | `$\mtol$` | 40 | tolerable metabolite level | convention, explicitly | a clinical judgement with no measurement behind it | assumed / assumed | none. `natural_history` supplies only a rate ratio, not a target | **value** (the undetermined window runs 1.2% to 67% as it moves) | none. The row says so itself |
| L6 | `$\spk/s_{\mathrm{base}}$` | 1.50 | illness amplitude, as a fold of baseline load | effective | every driver of load during illness: proteolysis, fasting, reduced intake, fever | assumed / assumed | **human, adults** (`muscle_catabolism`, *Metabolism* 1981;30:765-776, "adult normal subjects and during malnutrition, sepsis, and skeletal trauma") | order of magnitude, but a sign flips inside a narrow window of `s_peak` | **D8** |
| L7 | `$\kappa$` | not assigned | gain from recorded carnitine dose onto `k_e` | lumped | absorption, distribution, tissue carnitine loading, conjugation stoichiometry: a dose in mg/day mapped onto a rate constant in 1/day | assumed / not assigned | human but direction only (`carnitine_excretion_pa`, 2 PA + 2 MMA patients) | a **negative** conclusion depends on it being unknown, and the row says so | **none. The most honest row in the table** |
| L8 | `$\lambda$` | `ln2/3.5` d^-1 | decay rate of the therapeutic increment | measurable physical, in principle | nothing. It maps cleanly onto PCC protein turnover | **literature / literature-by-prediction** | **mouse, rat, monkey, allometrically scaled** (`k_deg` 0.0075 /h) | **value** (3.5 d gives n = 11 monthly visits, 7 d gives n = 30) | **D3, D6** |
| L9 | `$T_{\mathrm{dose}}$` | 14 d | interval between doses | measurable physical | nothing | literature / literature | human (`trial_phase12`), and the row **correctly** records "every two or three weeks" | value | none. **The best-sourced row in the table** |
| L10 | `$D$` | 6000 | size of one therapeutic increment | effective, a calibration choice | the entire dose-to-marker transfer: LNP delivery, translation, PCC assembly, and the marker's response to capacity | assumed / assumed | **mouse** transfer function, if recovered as the row proposes | **value**. `D` from 766 to 20977 across defensible readings gives n from >60 to 10 | **D7** |
| L11 | `$\sigma$` | a family | SD of the ordinate of the fit | lumped | the load recording error and the assay term, combined in quadrature by `sigma_total` | assumed / not assigned | none | value | **D9** |
| L12 | `$\mathrm{CV}$` | 5% | analytical variability of the assay | measurable physical | nothing | assumed / assumed | human (`breath_test_mma`, MMA not PA), plus a **mouse** IIV figure quoted beside it | order of magnitude | **D10, D3** |
| P1 | `$\theta$` | 300 | residual capacity, parent model | lumped, as L2 | as L2 | assumed / assumed | none | value | **D2** (the 5% / 10% half of it is `OPEN_QUESTIONS` (b)) |
| P2 | `$K_m$` | 100 | half-saturation, parent | effective, as L3 | as L3 | assumed / assumed | none cited, **although a supporting reference is already registered** | value | **D11**, constructive |
| P3 | `$K_c$` | 100 | half-saturation of the conjugation step | effective, model-internal | the whole saturation of carnitine conjugation | assumed / assumed | **rat, in vitro** (`octn2_renal_transport`), and the row does not say so | order of magnitude (`theta_crit` 1.56x over a 32-fold range at held calibration) | **D12** |
| P4 | `$M_i$` | 300 | level at which acylcarnitine competition doubles the non-conjugative loss | effective, model-internal | the whole competitive-inhibition mechanism at the transporter | assumed / assumed | **rat, in vitro** (`octn2_renal_transport`), mechanism only, and the row does not say so | order of magnitude for `theta_crit` (1.37x), **value for the margin (6.44x)** | **D12, D13** |
| P5 | `$\rho_0$` | 0.50 | basal non-conjugative loss relative to conjugation | effective, dimensionless ratio | the ratio `delta_0/k` | assumed / assumed | human, **direction only** (`carnitine_excretion_pa`) | order of magnitude (`theta_crit` 1.02x) | **D14** |
| P6 | `$R_0$` | 1122.0085 | carnitine replenishment | derived-internal | nothing, but it inherits everything its two inputs carry | derived / **derived from two assumptions, one unregistered** | none | value, and the row correctly warns it must be re-derived | **D15, D16** |
| P7 | `$s_{\mathrm{base}}$` | 600 | baseline propiogenic load | effective, sets the parent scale | the whole substrate supply | assumed / assumed | human turnover figures quoted **with no citation** | order of magnitude | **D17** |
| P8 | `$k$` | cancels | conjugation rate constant | **eliminated by the reduction** | nothing | derived / arguably a fourth class | none needed | neither | none. **Test-backed and sound** |

---

## Defects, with evidence

**D1. `k_e` cites a route that does not exist in the cited source. HIGH, and the most serious
finding here.** `src/common/params.py:76-79`:

> "The nearest route is the first-order 3-HP elimination implicit in the published
> indirect-response model~\cite{pkpd_translational}"

Attarwala models 3-HP as **"a direct sigmoidal maximum inhibition (Imax) model"** in PCC
protein. A direct algebraic function has no differential equation and therefore no elimination
rate constant. The indirect response in that paper is on **PCC protein synthesis**, not on the
marker. `paper/claims.csv` C46 already records this correctly ("no biomarker turnover constant
at all"), so `params.py` and `refs.sources:110-112` are on one side of an internal
contradiction and `claims.csv` is on the other, with the primary source agreeing with
`claims.csv`.

This is not a stray annotation. `k_e` is the parameter the design is most sensitive to:
`paper/sections/07_design.tex:40` records that its prior "is the only prior that does any work:
removing it takes the same design from `0.1474` to `0.6451`", and line 41 that "the gate needs
`k_e` within a factor of `2.95`, which is the one thing this paper asks the literature for".
And `paper/sections/08_specification.tex:129` asserts the design needs non-enzymatic clearance
"within a factor of about three, **which the literature supplies**; no measurement programme is
proposed for it."

**The literature does not supply it.** The registry row is honest that "no direct measurement
exists" and then names a route that is not there; the specification section claims the
literature supplies the number outright. The registry and the paper are wrong in two different
directions about the same parameter, and it is the one the whole design rests on.

**D2. The LOCAL `theta` is normal capacity, and a conclusion depends on its value. HIGH.**
`src/common/models.py:93` comments `ParamsC.theta = 300.0` as "residual V_max of PCC, **~10% of
normal**", which puts normal at 3000, and `ParamsA.theta = 3000.0` is commented simply
"V_max of PCC, maximal enzymatic velocity". `K_m = 100` is shared and the PARENT `K_m` row says
the sharing is "so the two can be compared", so the two blocks are on one scale. On that scale
the LOCAL illustrative patient sits at **100% of normal capacity** while its row calls the
value "residual capacity of the enzymatic channel" and justifies it by citing residual PCC
activity "under 5% in PA".

The row then says the value cannot matter:

> "the quantity the paper says cannot be identified, **so no conclusion may depend on its
> value**"

`paper/sections/04_local.tex:863-866` contradicts that directly:

> "Both scale as `1/theta` ... At `theta = 3000` they are `0.20` and `0.37` and the anchor
> **attenuates**; they pass one at `theta = 605` and `theta = 1100`; at the parent's own
> `theta = 300` they are `2.0` and `3.7` and the anchor **amplifies**. It is therefore worst
> for the severe patients the note is about, and **the favourable figure belongs to the mild
> end of the range**."

So a qualitative conclusion flips between the two blocks of the same table, the LOCAL value
sits on the favourable side of both crossings, and `paper/sections/08_specification.tex:226`
independently records that the coverage shortfall is `theta`-dependent in the same direction.
`OPEN_QUESTIONS.md` (b) covers the 5%-versus-10% half of this; **the LOCAL value being normal
rather than residual, and the "no conclusion may depend on its value" sentence, are not
covered anywhere.**

**D3. Five rows rest on an animal model presented as a human scale. HIGH, partly known.**
`m-bar`, `lambda`, `k_e`, `D` and `CV` all lean on `pkpd_translational`, which has no human
data. The `m-bar` row is the sharpest case because it sets the **units of the entire local
model**:

> "On a real scale: the translational model puts the untreated plasma 3-HP baseline at
> `60.1` uM and the fully corrected floor near `0.1` uM"

Both numbers are mouse PD estimates, and "on a real scale" reads as a human measurement.
`audit/PHASE4_PROPOSED_TEXT.md` passage 10 already proposes fixing five annotation rows that
"read mouse estimates as human" and is **unapplied**, so this is known. What is new is that the
species claim is now confirmed at the primary source rather than inferred from
`refs.sources`, and that 60.1 is the **suppressible pool** rather than the baseline.

**D4. `m-bar` has no unit in the Value column although it is the unit anchor. MEDIUM.** The
Value column reads "10". The note says the quantity "is read in umol/L". The only place the
paper states the marker's unit is `paper/sections/04_local.tex:312`, "The state `m` is a plasma
marker at micromolar", in a different section.

**D5. The illustrative patient is optimistically placed and no row says so. MEDIUM.**
`m-bar = 10` uM sits at about one sixth of the untreated total 3-HP of 60.2 uM in the cited
source, and the paper's own reach limit is that capacity stops being recoverable above roughly
three to four times the managed level, i.e. 30 to 40 uM. An untreated patient at 60 uM is
therefore already outside the recoverable band on the paper's own arithmetic. Nothing in the
table records that the illustrative operating point is at the favourable end.

**D6. `lambda` is declared `literature` here and `ASSUMED` in the code. MEDIUM.**
`src/common/params.py` gives the class as `LITERATURE` and calls the value "**Corroborated**"
and "within `10%` of a published value rather than free". `src/design/run.py:56` comments the
same quantity:

> `LAMBDA = np.log(2.0) / 3.5   # /day. ASSUMED. A 3.5-day half-life on the`

The registry and the code disagree about the epistemic class of one of only two `literature`
rows. On the evidence the code is closer to right: the corroborating `k_deg` is an allometric
estimate from mice, rats and monkeys. The **mechanism** mapping is sound, since the therapeutic
increment is to capacity and capacity is PCC protein, so `lambda` and `k_deg` are genuinely the
same quantity. Only the species and the prediction status are overstated.

**D7. `D`'s value matters more than the row implies. MEDIUM-HIGH.** The row says "Only the
excursion it produces is used", which is true and reads as reassurance.
`audit/PHASE5_ROBUSTNESS.md:727-731` shows what the excursion does:

| reading | `D` | fall at peak | monthly visits to 0.20 |
|---|---|---|---|
| committed | 6000 | 66.4% | **11** |
| 20% fall at peak | 766 | 20.0% | **>60** |
| 30% fall at peak | 1311 | 30.0% | **42** |
| 20% fall at trough | 12255 | 80.2% | **10** |
| 30% fall at trough | 20977 | 87.4% | **12** |

A 27-fold range of defensible `D` moves the headline design from 10 visits to more than 60.
That is dependence on the value, not on the order of magnitude, and the row does not say so.

**D8. `s_peak/s_base` bounds a whole-body load with a muscle-only index. MEDIUM.**
`muscle_catabolism` is human and adult, and 3-methylhistidine indexes **muscle protein**
catabolism. Per `propionate_sources` (Thompson 1990) amino-acid catabolism is about half of
propionate production, and only part of that comes from muscle, so a 227% to 325% rise in
3-methylhistidine does not license a 227% to 325% bound on total propiogenic load. The row uses
it only as an upper bound and then sets 1.50, far below, so **the conclusion is safe and the
bound is looser than it looks**. The subjects are also adults where the trial population's
median age is 8.5 years.

**D9. The `sigma` row's central claim is falsified by the code. HIGH, and a new extension of
R17.** The row says:

> "deliberately not assigned. **Every design conclusion is reported as a function of it**"

`REVIEW_LOG.md` R17 established that `src/design/run.py spec_numbers()` computes
`sigma = LOC.sigma_total(LOAD_FRAC, ASSAY_CV)` **once, at `LOAD_FRAC = 0.05`, and reuses that
one value on every row of the specification table**. So the design conclusions of Section 8 are
reported at one `sigma`, not as a function of it. R17 corrects the specification's prose; **it
does not touch the `sigma` row, which makes the same claim in the table that generates into the
paper.**

**D10. The `CV` row prints a retracted figure, and so does the generated table. HIGH, and
trivially fixable.** `src/common/params.py:133-134`:

> "the breath test reports repeatability of `6.9\%` in controls and `12.8\%` in patients"

`REVIEW_LOG.md` R4 retracted the 12.8%: Manoli's Results give `14.02 +/- 12.86%`, and the
abstract reports the standard deviation as though it were the mean. `results/tables/tab_params.tex`
**still prints both `6.9\%` and `12.8\%`**, while `paper/sections/08_specification.tex:72-73`
carries the correction in full. **Table 2 and Section 8 of the same paper currently disagree,
and Table 2 is on the retracted side.** `refs.sources`'s `breath_test_mma` annotation repeats
it too.

**D11. The PARENT `K_m` row is uncited where a supporting reference is already registered.
LOW, constructive.** The row's entire note is "as above, and shared with the local model so the
two can be compared". Sharing `K_m` across a tenfold difference in capacity assumes `K_m` is
independent of residual activity. `paper/refs.sources` already holds
`pcc_variants_kinetics` (PMID 15890657), annotated "Four variant PCCs: **NO significant
difference in Km for propionyl-CoA vs wild type**". That is direct support for the very
assumption the row leaves bare.

**D12. `K_c` and `M_i` cite a rat in vitro measurement without saying so. MEDIUM-HIGH, and
probably uncovered.** `refs.sources` annotates `octn2_renal_transport` as
"IN VITRO, RAT: the limitation is stated in the text". The `K_c` row says only "the nearest
measured analogue is the renal carnitine transporter ... which is a **different step**"; the
`M_i` row says only "the mechanism is measured, acylcarnitines being the most effective
inhibitors of that transporter". Neither names rat or in vitro. `PHASE4_PROPOSED_TEXT` passage
10 is about **mouse**; on its face it does not cover rat, so **this appears to be uncovered and
should be checked against passage 10 before being treated as new.**

**D13. `M_i`'s "cost bounded in supplement S1" is true of the threshold and false of the
margin. MEDIUM.** `REVIEW_LOG.md` section 0.1 reports, at a held calibration, `theta_crit`
moving 1.37x over an eightfold range of `M_i`, but the **margin** moving 6.44x and running 509
to 3275 over the same range while Figure 16 quotes 954 flat, and states in terms that "the
margin is **not** protected the way the threshold is". The row's flat "Cost bounded in
supplement S1" does not distinguish the two.

**D14. `rho_0`'s lower bound is a sign, not a number. MEDIUM.** The row reads as two-sided:
"Bounded above by `1` ... and below by the measured direction of free-carnitine excretion".
Above is a number; below is a direction. `REVIEW_LOG.md` section 0.1 separately records that
"the `rho_0` lower bound 'about 0.35' is an `R0`-held number; re-derived it is 0.195", so the
one lower figure that was quoted elsewhere did not survive.

**D15. `R_0` is derived from an unregistered constant. MEDIUM.** It is fixed by requiring the
alternative route to peak at `ALT_PEAK_C = 500` against `S_BASE_C = 600`.
`S_BASE_C` has a row. **`ALT_PEAK_C` does not.** A `derived` row in a registry whose stated
rule is that no number reaches the paper without a row is itself derived from a number with no
row. The row does warn correctly that changing `K_c`, `M_i` or `rho_0` forces a re-derivation,
which is the failure `REVIEW_LOG.md` section 0.1 documents in three places.

**D16. `R_0` is quoted to eight significant figures from two round inputs. LOW, and the most
visible.** `1122.0085` from 500 and 600.

**D17. `s_base` makes two literature claims with no citation. MEDIUM.** The row asserts
"Whole-body propionate turnover is published at roughly `250` to `400` umol/kg/day, partitioned
as about half from amino-acid catabolism, a fifth from gut fermentation and the rest from
odd-chain fatty acids" and carries **no `\cite`**. The partition is `propionate_sources`
(Thompson 1990), already registered and already annotated with exactly these figures, so the
citation is available for free. The turnover range is unattributed. There is also a tension
with `colonic_propionate`, annotated as "Endogenous colonic propionate rate of appearance
`0.27 +/- 0.09` umol/kg/min = ~389 umol/kg/day": if the gut is "a fifth" of 250 to 400, the gut
is 50 to 80, not 389. That source is Boets et al., *Nutrients* 2015;7:8916-8929,
"Quantification of in Vivo Colonic Short Chain Fatty Acid Production **from Inulin**", so its
figure is an inulin-stimulated colonic rate rather than a habitual whole-body one, which is
very likely the whole explanation. The row's closing disclaimer, "The value here is not in
those units and is not claimed to be", is honest and protects the constant itself.

---

## Cross-cutting findings

**X1. The two blocks are not at the same point on the saturation curve, although the table says
they are comparable.** Computed here: the LOCAL block sits at `m-bar = 10` with `K_m = 100`, so
`m/K_m = 0.10`. Solving the parent's clearance balance at `s_base = 600` puts its stable
attractor at **`m = 82`** and its separatrix at `m = 1032`, with the clearance peak at
`m = 276`, so the parent sits at **`m/K_m = 0.8`**. The PARENT `K_m` row says `K_m` is "shared
with the local model **so the two can be compared**", and the LOCAL `K_m` row says
"The ratio `m-bar/K_m = 0.1` **is what the excitation condition is stated against**".

An eightfold difference in saturation regime is not a detail here. `REVIEW_LOG.md` section 0.3
records that "moving the operating point from `m/K_m = 0.1` to `m/K_m = 2` drops the
collinearity to 0.982 and raises `det G` by fifteen orders of magnitude". The parent's own
operating point is most of the way along that interval. Sharing the number `K_m = 100` does not
put the two models in comparable regimes.

**X2. The registry does not meet its own stated rule.** The docstring says "no number reaches
the paper without a row saying whether it was measured, derived from another number, or
assumed". Constants a published result depends on that have **no row**:

| constant | where it lives | what depends on it |
|---|---|---|
| `ALT_PEAK_C = 500` | `src/common/params.py:53` | `R_0`, hence the whole parent calibration |
| `ParamsB.ke = 400`, `ParamsB.Ke = 300` | `src/common/models.py:81,83` | the condition number `4.8e7` and the factor of 450 quoted at `paper/sections/10_results.tex:217-218`, and `prop:exchange` |
| `ke_prior = 0.8` | `src/design/run.py:257`, a **default argument** | the only prior that does any work: `0.1474` to `0.6451` without it |
| `LOAD_FRAC = 0.05` | `src/design/run.py:44` | R17 shows this is the binding constraint on the whole design |
| `LOAD_CV = 0.10` | `src/design/run.py:89` | `prop:rank` and the excitation argument |
| `0.20` target | `src/design/run.py:366`, a literal | every "reaches the target" statement |
| the parent's operating point `m = 82` | computed, nowhere | every parent figure |

The most pointed of these is `ke_prior = 0.8`. The single most load-bearing number in the design
is a default argument in a function signature, in a repository whose registry exists so that no
such number can hide.

**X3. `k_e` carries two different dimensions under one row.** `ParamsA.ke = 5.0` is
"alternative clearance rate constant", units 1/day. `ParamsB.ke = 400.0` is commented
"**NOTE: now the `V_max` of the alternative route, a maximal flux and not a rate constant**",
units umol/L/day. The ceiling hypothesis is used in the paper, including the conditioning
result at `10_results.tex:218`. Table 2 has **one** `k_e` row, for the rate-constant meaning,
with no unit. `tab:epistemic` at `paper/sections/02_system.tex:155-157` carries both meanings in
one cell, "a rate constant in the local reduction, the `V_max` of the alternative route under
the ceiling hypothesis", which is accurate but means the symbol is dimensionally overloaded
across the paper.

**X4. No unit system is stated, and twelve rows carry no unit.** Only `lambda`, `T_dose` and
`CV` carry units in the Value column. Deriving the rest from the balance with `m` in micromolar:
the time unit is days (the marker time constant `1/((theta)g'(m-bar)+k_e) = 0.0336` d `= 0.81` h,
matching `paper/sections/07_design.tex:52`), `theta` and `s` and `D` are in umol/L/day, and
`k_e = 5` /day, which is **a non-enzymatic marker half-life of 3.3 hours**. That number appears
nowhere. The file's own docstring names a unit error as one of the three founding failures it
exists to prevent.

**X5. `ASSAY_CV` is defined three times in two unit conventions.**
`src/common/params.py:51` = `5.0`; `src/model_A/run.py:42` = `5.0`, commented "per cent";
`src/design/run.py:43` = `0.05`, commented "a fraction". Both comments insist the quantity is
"never flux units", which is the old failure being guarded, while the new hazard is that one
name now carries two scales an order of magnitude apart.

**X6. The `Source` column has three classes and needs four.** `theta` is **unidentifiable in
principle** and the paper is about that. `k_e` is merely **unmeasured**, and a measurement
would close it. `k` is **eliminated**: nothing derives it, it cancels out of the reduction
exactly, which `tests/test_models.py::test_the_two_state_system_reduces_to_the_quasi_steady_flux_symbolically`
verifies. All three are currently `assumed` or `derived`. The distinction between "no one has
measured this" and "this cannot be measured from these data" is the paper's whole subject and
the registry does not encode it.

**X7. The caption and the discussion disagree on the literature count, and the caption is
right.** `paper/sections/02_model.tex:262`: "Of the twenty registered constants, **two** are
taken from the literature and two are derived." `paper/sections/12_discussion.tex:7-8`:
"the constants of Table~\ref{tab:params}, of which **exactly one** is taken from the
literature." Counted from the dataclasses: 16 assumed, **2 literature**, 2 derived. The
discussion is wrong.

---

## Rows to change, in priority order

1. **`k_e`.** Strike the claim that the route exists in `pkpd_translational`. Say instead that
   no source contains a marker elimination constant: Attarwala's 3-HP is algebraic in PCC
   protein, and the only source with a first-order marker loss is `dose_selection_poster`,
   which `claims.csv` C46 already records as a conference abstract carrying "no PD equations
   and no biomarker parameters". Then reconcile `paper/sections/08_specification.tex:129`,
   which currently says the literature supplies the number. This is the parameter the design
   rests on and the claim is not supportable as written.
2. **The `CV` row's `12.8%`.** It is formally retracted in `REVIEW_LOG.md` R4 and still prints
   in `results/tables/tab_params.tex`. One edit, no argument required, and it removes a
   disagreement between Table 2 and Section 8 of the same paper.
3. **The LOCAL `theta` row.** Remove or qualify "no conclusion may depend on its value", which
   `04_local.tex:863-866` contradicts, and state what 3000 is on the shared scale. This is the
   row a reviewer is most likely to use to question the whole illustrative set.
4. **The `sigma` row.** "Every design conclusion is reported as a function of it" is false of
   Section 8 as coded. Align it with R17.
5. **Species flags on the five `pkpd_translational` rows**, which `PHASE4_PROPOSED_TEXT`
   passage 10 already drafts and which can now cite the abstract verbatim, and on the two
   `octn2_renal_transport` rows (rat, in vitro), which passage 10 may not cover.
6. **Units in the Value column**, and one sentence somewhere stating the unit system. Twelve
   rows are dimensional and unitless.
7. **`D`'s dependence.** Replace "Only the excursion it produces is used" with the range from
   `PHASE5_ROBUSTNESS`, which shows 10 to more than 60 visits.
8. **Register the missing constants** of X2, starting with `ke_prior = 0.8` and `LOAD_FRAC`.
9. **`M_i`'s cost claim** (threshold bounded, margin not), **`rho_0`'s one-sided bound**,
   **`R_0`'s precision and its unregistered input**, and **`s_base`'s two uncited claims**,
   which can cite `propionate_sources`, already registered.
10. **The PARENT `K_m` row** can cite `pcc_variants_kinetics` in support. This is the only item
    on the list that makes the paper stronger rather than more careful.

## What a pharmacometrics reviewer will challenge first, and why

1. **`k_e`, on two counts.** Borrowing a first-order elimination constant from a model that
   contains none is the kind of error this community is trained to catch, and the paper
   compounds it by telling the reader at `08_specification.tex:129` that the literature
   supplies the value. A reviewer who checks one citation will check this one, because the
   paper itself flags `k_e` as "the one thing this paper asks the literature for".
2. **`theta` and `K_m` estimability, and the units.** Reviewers in this field will immediately
   ask whether a `V_max`/`K_m` pair is separable from steady-state data at `m/K_m = 0.1`, which
   is far below half-saturation where the Michaelis form is nearly linear and the two
   parameters are nearly collinear. The paper's answer is good and is in fact its central
   result, but the table invites the question and then supplies `theta = 3000` with no unit,
   so the reviewer cannot check dimensional consistency without deriving the unit system.
3. **Mouse to human.** Five rows rest on a model whose abstract says it was built on mice, rats
   and monkeys and scaled allometrically to **predict** human behaviour. Allometric scaling of
   a clearance licenses a dose range; it does not license quoting a mouse `IIV` of 51% as human
   between-subject variability, or a mouse `E0` as "a real scale".
4. **Sixteen of twenty rows assumed, in one table with literature-cited rows under a single
   Source column.** A reviewer will read the table top to bottom and ask what is left that is
   measured. The answer is two rows, one of which is a dosing schedule.
5. **`R_0` at eight significant figures.** Small, certain to be noticed, and cheap to fix.
6. **Priors doing the work of data.** `07_design.tex:40` states plainly that removing the `k_e`
   prior takes the design from 0.1474 to 0.6451. A reviewer will pair that with finding 1 and
   ask what the prior is grounded in.

**What the same reviewer will find unusually good, and should be said plainly.** This registry
is better than the norm in this literature in four specific ways. It declares a source class on
every row rather than presenting all constants as if measured. It states what an assumption
**costs**, which almost no parameter table does. It is **generated from the code** and asserted
against the simulation in `tests/test_models.py`, so the table cannot drift from what was run.
And it contains rows that volunteer their own weakness: `kappa` is registered as "not assigned"
with an explicit statement of what that leaves undetermined, `sigma` is deliberately unassigned,
`mtol` says outright that it is a convention and not a property of the patient, and `k` records
that it cancels. The `T_dose` row is better sourced than the paper's own prose elsewhere,
because it records that the trial dosed "every two or three weeks" where
`src/design/run.py:52` and `tab:epistemic` both say "biweekly". Most of the defects above are
failures to apply this file's own standard consistently, not failures of standard.

## Method, and what was not done

Every claim above was checked against the file or the primary source and quoted. The Attarwala
full text was read from PMC10066765. The parent equilibrium in X1 was computed from
`clearance_C_quasistatic` rather than quoted. Two background reviewers, an outside
pharmacometrics lens and a completeness critic, were commissioned and had not returned when
this file was written; **anything they add is not yet folded in.**

**Nothing was changed but this file.** No `.tex`, no constant, no `refs.sources` entry, and no
row of `params.py`. The findings above are an audit, not a set of applied corrections, and the
priority list is a proposal for the author.

Two corrections to this session's own earlier notes, found while auditing:
`audit/PHASE6_REVIEWS_SEP2026.md` section 8 presents the trial's Q3W cohort 1 and the 52-day
versus nominal 42-day dose-3 interval as new; `audit/PHASE5_ROBUSTNESS.md:685-686` already
records both. And PHASE6 section 10 presents the content of the Manoli correction as new;
`REVIEW_LOG.md` R4 already records that it "restores a graphical abstract and changes no
numbers". What was new in PHASE6 was the registration of those corrections, not their content.

---
---

# Part 2. The outside reviewer's return, checked

The pharmacometrics lens commissioned in Part 1 returned after this file was committed. Every
claim below was **verified against the file before being accepted**; three of its claims did
not survive, and are recorded as such. The completeness critic has still not returned.

## What it found that Part 1 did not, and that checks out

**X8. D2 is understated. It is a scope objection, not a sensitivity objection. HIGH, and this
supersedes D2's framing.** `paper/sections/11_scope.tex:22-33` carries a table Part 1 did not
reach. Verified verbatim:

| `theta`, as % of nominal | marker, x managed | `abs(d ln m / d ln theta)` |
|---|---|---|
| 100% | 1.00 | **0.915** |
| 50% | 1.82 | 0.807 |
| 20% | 3.41 | 0.537 |
| 10% | 4.57 | 0.321 |
| **5%** | 5.40 | **0.173** |
| 2% | 6.00 | 0.072 |
| 1% | 6.22 | 0.036 |

The LOCAL block's `theta = 3000` **is the 100% row**, the top and most favourable line of the
paper's own scope table. Both `theta` rows state that residual PCC activity in PA is under 5%
of control, which is the `0.173` row, a fifth of the illustrative sensitivity and still
falling. So the headline `sd(ln theta-hat) = 0.170` is computed at a capacity no patient in
this indication has, and Section 11 already says in its own words that above about three to
four times the reference level "residual capacity is no longer recoverable".

The coverage study offered as the generalisation check does not close this.
`paper/sections/08_specification.tex:216-217` says the sweep is "log-uniform over `theta`
(nine-fold)" about the same 3000, so at its widest it reaches roughly a third of nominal and
**never enters the disease range at all**. Part 1's D2 said a conclusion flips between 3000 and
300. The correct and stronger statement is that the paper's own Section 11 shows the method
decaying to near zero across exactly the interval the patients occupy, and every headline is
computed outside it.

**X9. The main-text caption asserts the opposite of the table it captions. HIGH.**
`paper/sections/02_model.tex:265-266`, verified verbatim:

> "**No conclusion below depends on a value in this table**, and the design results are
> reported relative to the patient's own operating level."

`paper/sections/08_specification.tex:223-224`, verified verbatim:

> "What spreads it is `Delta/theta`, which correlates with the achieved `sd` at `+0.56` on
> logs, alongside **`K_m` at `+0.59`**."

So `K_m`'s value is the single strongest correlate of the achieved precision in the paper's own
coverage study, and `Delta/theta` is next. The caption is false for at least `theta`, `K_m` and
`D`. Note this is a **different and larger** claim than `rem:no_km`, which establishes only
that no *prior* on `K_m` is needed: `07_design.tex:36-39` shows the design returns `0.1474`
with the `K_m` prior, without it, and at a prior twenty times looser. "We need no prior on
`K_m`" and "nothing depends on `K_m`" are being conflated, and only the first is proved.

**X10. `ILLNESS_FOLD = 1.50` exceeds the load at which the parent model has no equilibrium.
MEDIUM-HIGH.** Confirmed by my own computation in X1: the parent's clearance peaks at
**706.4 at m = 276**, and `S_BASE_C = 600`, so the saddle-node sits at **1.177 x baseline
load**. The registered illness amplitude of 1.50 is 27% beyond it. The row
(`src/common/params.py:90-95`) discloses "Model C's own fold sits at 1.18" and in the same
sentence calls 1.50 "a conservative value for a moderate intercurrent illness". Conservative
for the local model's `theta_crit` arithmetic, and past the point of no return for the parent.
The two senses are not distinguished.

**X11. The same generated table is typeset twice, and only the supplement defines its Source
column.** `paper/main.tex` places `tab:params` via `02_model.tex`; `A1_supplementary.tex:547`
places the same `results/tables/tab_params.tex` as `tab:params_orig` via `02_system.tex:250-253`.
Only the supplement caption (`02_system.tex:251`) defines `literature`, `derived` and `assumed`.
A reader of the main paper meets a three-valued Source column that is never defined in the
document they are holding.

**X12. A second retracted number is still live, this time in the code.**
`src/common/models.py:109` comments `rho0` as "Bounded below at about **0.35** by the same
source". `REVIEW_LOG.md:59` retracted it: "lower bound 'about 0.35' is an `R0`-held number;
**re-derived it is 0.195**." So alongside D10's retracted `12.8%` in the `CV` row, the tree
carries a second figure the review log has corrected. Neither is in the `params.py` row text,
so neither reaches Table 2 directly, but both reach a reader of the source.

**X13. `paper/claims.csv` references four section files that no longer exist**:
`sections/03_model_A.tex`, `sections/05_model_C.tex`, `sections/06_results.tex`,
`sections/07_open.tex` (plus the sentinel `paper-wide`). The claim register has drifted from
the section layout, which matters because C46, the claim that settles D1, is one of the rows
pointing at a dead path.

## Three claims from the reviewer that did not survive checking

**R-a. It praises the `rho_0` lower bound of "about 0.35" as exemplary practice.** The *method*
is good and the praise is kept, but the *number* is retracted (X12), and it is not in the
`params.py` row at all, which states only a direction. The reviewer read `models.py:109` and
took a corrected figure for a current one.

**R-b. It states the coverage sweep is "centred on the same 3000".** The text says log-uniform
over a nine-fold range; "centred on" is an inference, not a quotation. The conclusion it draws
is nonetheless right, because even the bottom of a nine-fold range about 3000 is roughly a
third of nominal and the disease range is under 5%.

**R-c. It attributes "prior-free the design gives 1.03" to `07_design.tex:140-147`.** That span
says "does not identify it at all without", which supports the direction but is not the number
quoted. The figure `1.03` was not located and is not relied on here.

## What this does to the priority list

The Part 1 list stands with one change. **Item 3 moves to item 1** and is restated: re-anchor
the LOCAL `theta`, because Section 11 shows the method dies across the range the patients
occupy and every other number in the table is computed above it. Either move the illustrative
`theta` into the disease range and re-report, or relabel the row as the nominal non-disease
reference, add a disease-range row, re-centre the coverage sweep, and delete both "no
conclusion may depend on its value" (the row) and "No conclusion below depends on a value in
this table" (the main caption). The `k_e` route (D1) and the retracted `CV` figure (D10) follow
as before.

The reviewer's own ranking put this first for the same reason, independently.

## What this registry does better than the norm

The outside lens was explicit and unprompted on this, and it is worth recording against the
length of the defect list. The `kappa` row is the correct handling of an unidentifiable gain
and is rarely seen: it names what is established (a sign), what is not (a magnitude), and what
the non-assignment costs (a rank deficiency). The `K_m` row **refuses** to import a citable in
vitro constant on principled grounds and then demonstrates it did not need one, where the
common failure in this literature is the opposite. The `m_tol` row ranks itself as "the most
expensive assumption in the local model", and a table that ranks its own assumptions by cost is
unusual. The fifth column, "what fixes it, or what assuming it costs", does not exist in the
standard format. And several of the sharpest findings in this audit came out of the
repository's own records contradicting its own table, which is a system working mid-cycle
rather than a system failing.

---
---

# Part 3. The completeness critic's return, checked

The second commissioned lens returned 33 items. Verified as in Part 2. The ones below are new
and check out; three did not survive and are recorded.

## New and verified

**X14. `kappa` cites the wrong registered source, and the right one is already registered.
MEDIUM, and free to fix.** `src/common/params.py:101-102` says the direction "is established,
**carnitine supplementation raising** propionylcarnitine excretion~\cite{carnitine_excretion_pa}".
`refs.sources` annotates that key as describing the **untreated disease state**: "Two PA and two
MMA patients ... URINARY EXCRETION OF FREE CARNITINE DECREASED, while excretion of short-chain
acylcarnitines ... INCREASED." No supplementation. The supplementation experiment is registered
under a different key, `coa_sequestration`: "Carnitine **given** to an MMA patient ... raised
short-chain urinary acylcarnitines, and reduced methylmalonate and methylcitrate excretion."
So the row's claim is true and its citation points at a study that does not make it. Re-cite to
`coa_sequestration`, noting n = 1 and MMA.

This blunts Part 1's praise of `kappa` slightly. The row's **handling** of an unassignable gain
remains exemplary; its one citation is misattached.

**X15. Seven registered sources are cited nowhere in the paper, and four of them bear on rows
that say no source exists. MEDIUM-HIGH.** Counted here: `\cite` occurrences across
`paper/sections/` and `paper/*.tex` are **zero** for `dbs_3hp_cutoff`, `pa_guidelines`,
`propionate_sources`, `colonic_propionate`, `metronidazole_pa`, `muscle_catabolism` and
`mrna_moa` (`breath_test_mma` returns 7, so the count is working).

Four of those are the missing citations this audit has already found separately:
`propionate_sources` for `s_base`'s partition (D17), `muscle_catabolism` for the illness fold's
227% to 325% (D8, which quotes the figures with no `\cite` at all), `pa_guidelines` for
`m_tol`'s "no numeric plasma 3-HP target is published", and `dbs_3hp_cutoff`, whose 38.13 umol/L
99.9th-percentile cutoff sits beside `M_TOL = 40` with no row saying whether that is a
coincidence.

The most interesting is `metronidazole_pa`, annotated "Metronidazole 10-20 mg/kg/day in 5 MMA +
4 PA children. Mean decrease in urinary propionate metabolites 41% ... mean plasma propionate
45.0 -> 25.1 umol/l." That is a **human in vivo load perturbation with a measured plasma
response**, which is the closest published object to the excitation experiment this paper
designs, and it is cited nowhere. Whether it bounds `k_e` or the local slope is a question the
`k_e` row's "no direct measurement exists" should answer rather than pass over.

**X16. The claim register does not gate the table, and it extends D-nothing: this is new.**
`paper/sections/A2_reproducibility.tex:53-56` states that "Every factual assertion has a row in
`paper/claims.csv` ... `make audit` fails while any row is contradicted or unsourced".
`claims.csv` has **46 rows across 9 file targets, none of which is `params.py` or the parameter
table**. The table is injected from `results/tables/tab_params.tex`, so its cited assertions and
its three uncited ones (the `<5%` figure, the 250-400 umol/kg/day figure, the 227-325% figure)
never reach the gate at all. Four of the nine targets are also dead paths (X13).

**X17. The untreated patient sits at 93% of the null ceiling. HIGH, and it sharpens D5 and X8.**
Computed here from the registry: the load at the operating point is `L = 322.73`, so the
`theta -> 0` ceiling is `L/k_e = 64.55 = 6.45x m-bar`, which matches
`paper/sections/11_scope.tex:16` verbatim ("which at the illustrative values is `6.45` times the
managed level"). The untreated 3-HP the `m-bar` row itself quotes, 60.1 uM, is **6.0x m-bar,
which is 93% of that null ceiling**, and Section 11's table gives an elasticity of `0.072` at
6.00x. So on the row's own quoted scale an untreated patient is indistinguishable from a patient
with no capacity at all, and `m-bar = 10` is the favourable end of the interval. The row says
only "a managed patient sits between the two" and nothing about where.

**X18. The `<5% of control` figure is unsourced at the audited commit, and is the only empirical
statement about `theta` in the paper.** It appears at `params.py:67` and `:144` with no `\cite`,
and `refs.sources:177` shows the candidate `residual_carboxylase` commented out in a
not-resolved block. **Scope correction to the critic:** this is true at `f71bf9a` and is already
closed on this branch, where commit `5197106` registered `pcc_residual_activity` (Hsia et al.
1979, PMID 481943, verified: "Each of seven lines from patients with propionic acidemia had
clearly detectable PCC activity which was less than 5% of that in control lines"). The two
`theta` rows can cite it today, with the caveat that it is fibroblast enzymology rather than in
vivo hepatic capacity.

**X19. `K_c`'s stated cost basis spans a range the repository's own admissibility analysis
excludes.** `src/common/models.py:95-101` bounds the cost "over `Kc` from 25 to 800, a 32-fold
range". `src/model_C/run.py:615-621` records that the published figures for this "do not
reproduce under either" convention, and `REVIEW_LOG.md` section 0.1 gives the corrected reading:
the admissible window is two-sided, roughly `K_c` in `[100, 300]`, inside which `theta` moves by
1.41. So the cost is quoted over a range most of which is inadmissible, and `K_c = 100` sits at
the edge of what survives.

**X20. "Cost bounded in supplement S1" points nowhere locatable.** Three rows use it
(`params.py:99`, `:151`, `:159`, `:181`) for three different things. S1 is the whole supplement:
`A1_supplementary.tex:7` is its only `\section` and seven files are `\input` inside it. The cost
bounds actually live in `rem:shape` at `paper/sections/03_dynamic.tex:395-425`.

**X21. "Twenty registered constants" is seventeen.** The main caption
(`02_model.tex:263`) counts twenty, but `kappa` is "not assigned", `sigma` is "a family" and `k`
"cancels", so three rows carry no constant. Twenty **rows**, seventeen constants.

**X22. Six rows violate the file's own governing rule, and `K_m` is the one that matters.**
`params.py:19-20`: "An `assumed` row is not a defect on its own. An `assumed` row **with no bound
on what it costs** is." No cost is stated for `K_m` (both blocks), `CV`, `D`, the illness fold,
or the PARENT `theta`. `K_m` is the costly omission because the paper's central negative result
is stated at `m-bar/K_m = 0.1`: `07_design.tex:291-293` says the routine sources deliver
"between a fifth and a quarter of what their own fold range could deliver, **because they happen
an order of magnitude below `K_m`**". Move `K_m` and that verdict moves. The caption's protection,
that results are "reported relative to the patient's own operating level", does not reach it,
because the excitation measure depends on `m/K_m` and not on `m/m-bar`. **The cost of assuming
`K_m = 100`, expressed as the range of the excitation measure over a defensible `K_m` range,
does not exist anywhere in the repository.** That is the single largest computable gap this
audit found.

## Three critic claims that did not survive

**R-d. It says the capacity-ceiling test "passes only because it adopts the other reading".
False.** `tests/test_models.py:1533` carries
`@pytest.mark.parametrize("frac", [0.05, 0.10])` and the test is named
`test_the_capacity_ceiling_holds_at_both_readings_of_the_reference`, with a docstring stating
"**neither reading is chosen here: both must hold**". The test runs at both fractions and
asserts the ceiling at each. The critic read the body and not the decorator. The underlying
open question, 5% or 10%, is real and is `OPEN_QUESTIONS.md` (b); the test is not evidence for
it.

**R-e. Its `rho_0` item is a consistency defect only, and it says so itself.**
`src/model_C/run.py:663` still samples `rng.uniform(0.35, 1.0)` using the bound `REVIEW_LOG.md:59`
re-derived at 0.195. The critic recomputed `theta_crit` at 147.1 and 152.1 across the corrected
range against 149.7 to 152.1 across the stale one, so the conclusion does not move. Recorded as
the stale number it is (X12), not as a numerical error. Credit where due: the critic bounded its
own finding rather than inflating it.

**R-f. Its `sigma` units item is marked "COULD NOT SETTLE" and is left that way here.** It
observes that `sigma_total` is built in flux units while `A2_reproducibility.tex:47-49` describes
a log-scale observation model under which "`sigma` has unambiguous units", and could not
determine which construction the specification numbers use. Neither could this pass. It is
recorded as an open question, not a finding.

## What this does to the priority list, finally

No reordering. X22 joins the list as a new item: **compute what assuming `K_m = 100` costs**,
because it is the one gap here that is both load-bearing and cheap, and the repository's own
rule already requires it. X14, X15 and X18 are citation repairs that cost nothing and close four
rows that currently assert "no source exists" while the source sits registered and uncited.

The audit is complete. Both commissioned lenses have returned and been checked, three of their
claims were rejected, and nothing outside this file has been changed.
