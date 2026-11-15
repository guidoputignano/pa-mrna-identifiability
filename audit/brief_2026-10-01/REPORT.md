# Practical identifiability of enzyme capacity at 5% of normal

**Brief of 1 October 2026, worked on the repository's own footing, over nineteen rounds: the original eight tasks, a six-point follow-up, a nine-point round on conditioning, specification, and feasibility, and the rounds recorded in R-A to R-S. A plain-language account of the whole sequence, with its plots, is in `story_so_far.html`.** Every number comes from a script in this directory. `out/<script>.txt` holds each script's printed output and `fig/` its figures. The follow-up scripts (`t8`–`t11`, `t6c`, `t6d`) and the round-3 scripts (`t13`–`t19`) were re-derived by independent verifiers, by routes sharing no code with them; `t12` (the flux crossing) was not. Repository constants and the `n6` output this report quotes are printed in `out/t0_footing.txt` and `out/n6_profile_vs_crb.txt`. The footing is imported from `audit/n4_theta_scope.py`, not restated, and `core.fim` reproduces `n4.sd_ln_theta` to 1e-13. Claims that are assumed or not verified are collected in the [last section](#unverified-or-assumed).

> **Both headlines survive the parameters.** The two load-bearing conclusions — "no feasible design reaches 0.20 on the marker" and the ~19% flux contour — were stress-tested against K_m over a 4× range (50–200), assay CV 3–8%, load scatter 2–10%, the k_e-prior width, and the repo-vs-per-draw noise model: **both hold** (the feasibility sd stays ≥ 0.34; the flux crossing stays 18.8–20.2%) (R-B). The caveat that makes that test necessary: this is a **design study on a specified, simulated model, not an analysis of a cohort** — the repository holds no patient data (item 10; the brief's "24 patient records" was a mislabelling), so every number is conditional on the model and its constants. The saturation regime that governs which approximations hold is in R-A (at 5%, m/K_m ≈ 0.54, g ≈ 0.35 — sub-saturated).

> **How to read the numbers (round-14 audit).** A **log-scale sd** — every sd(ln θ), sd(ln λ), sd(ln k_e), the monitoring contrast, the interval *grid-resolution* tolerances — comes from the **FIM** and is reliable, tightening ~1/√n with the record. A **natural-scale or one-sided decision** built on the same posterior — a threshold crossing, a margin's lower bound — is **not** reliably FIM-propagated once sd(ln θ) is large (~3 at 5%): the log-normal posterior makes it tail-dominated, so it is read by **simulation** instead. Those decisions are the scenario-table Monte-Carlo column, the flux crossing (18.4%), round-13's eligibility, and round-14's spot-checks. The one place the two propagations were compared directly they disagreed by ~450 µM/d (R-M). The full FIM-vs-simulation classification, the two spot-checks (round-11's determined-subgroup and the cross-sectional headline — both survive), and a scan for hidden favourable cancellations are in **R-M**.

## The decision

> **Scope, stated first (round-7/8 finding).** Everything below is a statement about the **stable inter-crisis steady state**, computed on the parent paper's **Model A** (linear, unbounded non-PCC clearance; `sections/03_model_A.tex`). The round-6 literature shows non-PCC disposal in PA is **low-capacity** — patients accumulate and decompensate — and the repo's own Model C encodes this ("the alternative route alone peaks at 500 against a load of 600, so a patient with θ = 0 cannot be maintained"). Under a bounded, literature-consistent clearance a steady state exists only while the load fits under the disposal ceiling: θ > θ_crit = L − k_e. Because **L is actively managed** (protein restriction, gut suppression), a severe patient *does* have a steady state — **at their managed load**: at 5% capacity, untreated, they are below θ_crit and diverging, and a ~10% load reduction brings them above it (R-H, `t35`). So the correct scope statement is: **the analysis describes the managed, stable inter-crisis state; it does not describe the untreated or decompensating patient, and the margin μ = θ − θ_crit is thin — the record's own antibiotic/protein excursions (paper §3.5) can erase it, which is a decompensation.** Where a steady state exists, switching Model A → B moves the numbers (enzyme share 0.16 → 0.37; cross-sectional sd 3.05 → 0.29–0.92; contrast 0.13 → 0.08) but **not the sign of any verdict**, and B is *better* on both axes (the enzyme carries more of the load). The published figures inherit Model A's affine-clearance assumption — the paper states this explicitly (claim C23) — and B is not uniformly harder to identify: the marginal precision on θ *improves* under B even as the joint conditioning degrades (R-H, a correction to the paper's Remark 15).

**The general result (lead):** **only a direct, calibrated measurement of the enzyme reaches the target cross-sectionally. Inferring absolute capacity through the mass balance fails no matter how well the nuisances are measured** — the repo's full nuisance panel (an IV ¹³C₃ k_e tracer at 5.4% *and* an isotope-dilution load turnover at 10%) gives only **0.44** at 5%, and even k_e-exact + load-exact stays ~0.29 at monthly sampling. The two sensors that *do* reach target both read the enzyme itself, and both carry an **uncalibrated multiplicative factor** a single reading cannot separate from θ (R-G, item 1). The answer then splits by clinical use:

**1. Classifying absolute capacity — cross-sectional.** Reached only by a direct enzyme measurement, *and only if its calibration is supplied*:
- an **enzymatic-flux reading**, O = c·θ·g: reaches 0.20 at √(Y² + Z²) ≤ ~18–19% (reading precision Y, recovery calibration Z); with c uncalibrated it is worthless (R11).
- an **ex vivo PCC activity assay**, A = κ·θ: the repo's 20% is *repeatability only*, and κ (the cell-to-whole-body mapping) is unreported — structurally the breath's recovery factor. So the earlier "0.200" silently assumed κ known; it is knife-edge even then (25% repeatability → 0.25), and with κ **uncalibrated the route reverts to 3.0** (item 1). Same contour, √(Y² + Z²) ≤ 0.20.

So both direct routes need their multiplier (c or κ) calibrated for an *absolute* reading, and the repository supplies neither. No marker-reading or nuisance-sensor design reaches target at all.

> **This cross-sectional conclusion is model-dependent, but the escape it opens is a point, not a route (round-8, `t38`, `t40`).** Everything in the paragraph above is a **Model-A** result (linear, unbounded clearance), and the enzyme carries only a 16% share there, so θ is barely identified (3.05) and even perfect nuisance sensors floor at ~0.29. Under the **literature-consistent Model B** (saturating clearance), at a managed operating point the enzyme carries a larger share and the cross-sectional sd is **0.29 — missing 0.20 by ~1.5×, not ~15×.** Tighten the k_e prior to what the IV ¹³C₃ tracer delivers (5.4%) and one constructed point touches **~0.19**, at the gate. But on scrutiny that 0.19 is a **gate-hugging point, not a route**: swept across the plausible ceiling range under those same tightened conditions (`t40`), it only ever *hugs* the gate — at the tracer's k_e precision the whole ceiling sweep sits at **0.18–0.20, marginal everywhere and never comfortably below**, and at a looser k_e prior it slips off the gate (0.21–0.23) at higher ceilings; sd tracks the *ceiling* (the enzyme share), not the margin; and it holds at all only while **K_e is known to ~30%** (free K_e reopens the exchange degeneracy and sd blows up to ~9). So the honest reading is: Model A's low share makes the cross-sectional problem *look* worse than the mechanism requires, but under B the target is still **not reached by any route the data support** — it is touched only at a low, unmeasured ceiling with an unmeasured half-saturation assumed known.

> **The general negative result (round-8/9 — this is the report's conclusion).** The recurrence of "one unmeasured constant per route" is not a coincidence; it is forced by the structure. *An observable fixes absolute θ only in one of two ways, and each carries its own unmeasured constant:*
>
> - **Through the mass balance.** If the observable is (a function of) an accumulated metabolite state, its steady level is set by `load = θ·g(m) + [competing disposal]`. Attributing that balance to θ rather than to the competing channel requires the **competing channel's shape** — its V_max and half-saturation. These are k_e, K_e (Model B) and k, K_c, R₀, M_i, ρ₀ (Model C). None is measured in humans.
> - **By direct transduction.** If the observable is a scaled reading of the enzyme or its flux, `y = (const)·θ·(known function)`, recovering θ requires the **transduction constant**: κ (ex vivo) or c (breath/flux). Neither is measured in humans.
>
> **There is no third kind**, because θ enters the model only as the V_max of the flux θ·g(m): an observable depends on θ *either* through the metabolite state(s) it helps set (→ needs the competing channel to invert the balance) *or* directly through that flux (→ needs the proportionality to read it). An observable of an *upstream* pool (C3, 2-methylcitrate) needs **both** — a calibration to read the pool and the pool's own balance — so it is the product of the two kinds, not a third. **Checked against the repo's full ten-sensor panel** (`t26`/PHASE10): every one falls on a side (see R-G table). **And the decisive test — Model C, the only model where the competing channel has its own observable (free carnitine) — does not break it (`t41`, `t43`):** reading θ through it *still* needs the channel's shape constants (K_c/R₀/M_i/ρ₀ free → sd ~16; known → 0.34), and those are exactly the arbitrary (K_c, M_i), reverse-engineered (R₀) and one-sided-bounded (ρ₀) ones. Round-10 sharpened this: the carnitine observable **never earns its keep** — m-alone with the channel known already gives 0.34 (c is *redundant* for θ there), and in the channel-free regime c helps only marginally (20→16) and that help *erodes* once the realistic plasma→pool compartment gap is admitted (free carnitine is <1% of body carnitine, a different compartment from the modelled pool; its unknown, load-drifting partition is a **fifth** unmeasured constant, `t43`). So even a directly observed competing channel does not escape the pattern. **The conclusion: the accumulation marker cannot pin absolute capacity cross-sectionally without importing a constant unmeasured in humans — a calibration or a competing-channel shape — and across all three structural hypotheses (A/B/C) and the whole sensor panel, the repository and the human literature supply none.**

> **This is also why monitoring works (round-9, item 3; R-D).** κ, c, k_e, K_e and the load-bias b are *all patient-level constants*, and a within-patient contrast differences every one of them out. So the longitudinal result is **not a lucky exception** to the cross-sectional failures — it is immune to them *by the same mechanism that causes them*: the unmeasured constant each route must import is a patient constant, present identically in both windows, so it cancels in the change even though it blocks the absolute level.

> **The same theorem, derived independently from the biology (round-10 validation).** The round-10 literature reaches this conclusion from the other side, without the model: *at metabolic steady state the net carboxylation flux simply equals the rate of propiogenic-substrate disposal, so any concentration you measure reports a pool size set by inflow/outflow balance, not the throughput of the enzyme — and no concentration measurement can isolate the carboxylation step without a tracer.* That is the identity this brief proves from the sensor/FIM side (θ enters only as the flux's V_max, so a state read inverts the balance and a flux read needs a calibration). Two independent derivations, model and biology, land on the same negative result. And they close the one apparent escape — a direct flux read — the same way: **the 19% flux specification (R11) has no instrument.** The only observable truly proportional to carboxylation flux is a tracer-modelled V_PCC, which is rodent-validated only with no published human per-subject CV; and the ¹³C-propionate breath test is disqualified by its own criterion — an uncalibrated bicarbonate recovery (50–100% across states, not calibrated per subject) and a within-subject CV of 14.0 ± 12.9% (range 0.4–50.8%) in the one repeated-test series. So the flux route is not a pending measurement either: no human instrument for it exists.

*The **ratio route** removes the κ multiplier but not the deeper problem, and which it is turns on the estimand.* If θ is defined as A_patient / (healthy-control mean), the population mapping κ **cancels by construction**: a direct ratio read reaches **≈ 0.11–0.19** with the patient read to ~10% and 30–100 controls, no known-θ anchor (`t28`). **But that answers for *in vitro residual activity*, and the repository's θ is *whole-body capacity* — V_max in the mass balance (`models.py`).** The literature (user's round-6 research) is explicit that in vitro activity predicts in vivo capacity only weakly — the one same-subject study (Manoli 2021, in MMA) gives R² = 0.327, and the field calls the mapping "difficult or impossible to reduce to a single scaling factor." That introduces a **per-patient** scatter η that, unlike the population constant κ, **does not cancel in a ratio**. At the measured R² and a between-patient spread of ~0.8 log-units, η ≈ 0.66 — about **5× the ratio route's 0.135, far above the 0.20 gate** (R-G, item 1; `t29`). So the ex vivo ratio route reaches target **only if the clinical question is in vitro residual activity**; for the whole-body capacity the model carries, it does not, and the cross-sectional verdict reverts to **no inference-through-the-proxy route reaches 0.20** — the lone exception being a measurement that reads in-vivo capacity *directly* (the ¹³C-propionate breath oxidation test), bypassing the in-vitro mapping. The scatter figure is borrowed from MMA and is **unmeasured in PA**; the honest statement is that the ratio route silently assumes a tight mapping the only relevant data contradict.

**2. Monitoring within-patient change — longitudinal.** Here the marker reaches target *on its own*, and the calibration problem **disappears**: c, κ, k_e, b and K_m are patient constants, so they cancel in the contrast. sd(ln θ₂ − ln θ₁) is **≈ 0.13 for two weekly windows and ≈ 0.10 for two 6-month daily windows** — enough to detect a **+45% rise or −30% fall at 80% power** (weekly), needing no new assay and no calibration (R-D). It depends on draws-*per-window*, not calendar span, so short dense windows suffice. **Two large caveats from the round-6 literature**, which revise R-D: (i) the 0.13 is a model floor under *perfect* nuisance sharing and is **empirically untested** — the one cohort that collected these markers longitudinally (Moderna, n=50) found them too variable to use past baseline, and named real month-scale confounders (eGFR→k_e, carnitine dose→calibration) that the contrast's own stress tests show it tolerates only in small amounts; so the claim holds only for short intervals on a stable regimen with standardised sampling (`t30`). *Quantified against the only empirical serial-CV anchor (round-10): at a 5% per-draw CV the weekly contrast is 0.13, but at the one measured serial CV for C3 (26.4%, healthy adults) it is ~0.47 and would need ~590 draws per window for 0.20 — so monitoring is deliverable only if the true (unmeasured) PA serial CV is near the low end, <~10% (`t42`, R-D).* (ii) The **¹³C-propionate breath test is probably the better monitoring instrument**: it reads in-vivo capacity directly, its recovery c cancels in the contrast, and at a plausible ~8% repeatability a *single* reading per window matches 104 weekly blood draws (~0.13), a few readings beat it — and its reproducibility is the one such property the literature actually documents (`t31`, R-D).

**3. The one recommendation that costs nothing: record antibiotic status.** It is a precondition for *both* routes and is currently absent from the load record. Gut-flora propionate is antibiotic-suppressible (~40% load swing under metronidazole), so an unrecorded antibiotic course between two monitoring windows fakes a capacity change of **+1 to +3 ln units (a 3–20× apparent rise)**, destroying the contrast, and it collapses the cross-sectional estimate as well (R-E). **This is a safety point, not only a precision one: the fake is in the direction of apparent improvement, and antibiotics are started around clinical concern — so the failure mode is false reassurance exactly during deterioration.** Recording antibiotic status (and holding it fixed across a monitoring pair) is free and mandatory.

The cross-sectional routes, for completeness — everything but a *calibrated direct enzyme measurement* rests on an unobtainable premise, research-grade hardware, is clinically inadmissible, or is uncalibrated:

| route | verdict at 5% |
|---|---|
| **calibrated enzymatic-flux reading (O = c·θ·g) to Y ≤ ~18%** | **reaches 0.20** — only if absolutely calibrated. Contour **√(Y² + Z²) ≤ ~18–19%** (reading Y, recovery Z); flat at ~19% across capacities ~1–20% (`t16`/`t21`, R11/R-B). The repository's breath test reads a *fraction oxidised* with unknown recovery — off the contour until calibrated. |
| **ex vivo PCC activity (A = κ·θ)** | **reaches 0.20 only if κ is calibrated.** The repo's 20% is repeatability; with it, 0.200 (knife-edge — 25% → 0.25). κ (cell-to-whole-body mapping) is unreported, and **uncalibrated the route reverts to 3.0** — the same calibration problem as the breath test (`t27`, item 1). Contour √(Y² + Z²) ≤ 0.20. |
| nuisance sensors: IV ¹³C₃ k_e tracer (5.4%) + load turnover (10%) | **0.44 — does not reach 0.20**, even combined; even k_e-exact + load-exact ~0.29 (R-G). Inference through the balance fails however well the nuisances are known. |
| k_e known **exactly** + hundreds of draws | reaches 0.20 in principle (309 draws/2 yr) but **not feasibly**: k_e is achievable only to ~0.2 log units, and at that precision no sampling intensity reaches 0.20 (R12). |
| drug-forcing burst, peak 1500 | reaches 0.20 (0.131) only with a **research-grade 0.5% assay CV and κ known exactly** (`t10`). |
| load excursions large enough to reach 0.20 | **inadmissible** — every one sits past the parent model's crisis threshold of 1.177× baseline load (R15). |
| a cohort, a washout, a second-reading average, measuring the load | **none reach 0.20.** |

**k_e is the binding nuisance; the load matters only as a route to it.** *(This corrects the v2 "everything bottoms out at the load" framing, which the multi-observable panel showed is partly wrong — see the changed-claims table.)* k_e, not L, is what blocks cross-sectional θ: a **direct k_e sensor** (the IV ¹³C₃ tracer) at 5.4% takes the record to 0.45 **without touching L at all**, whereas measuring the load (turnover, 10%) does almost nothing (3.05 → 3.04). The load enters only indirectly, as the *route* to k_e through the safe transient (k_e's precision there is set by how well L₀ is known), and that route is worse than the direct sensor (turnover's 10% L supports only k_e ≈ 0.29). L is still the biggest term in the finite-n variance partition (62% at 25 draws) and the washout intercept, but neither of those is what decides the cross-sectional verdict. The reason no sensor suffices is deeper than the load: inference through the mass balance cannot separate θ from its clearance however well either nuisance is measured; only a direct enzyme read does, and only when calibrated.

**And measuring L to 1% is still not enough on its own.** If L₀ were known to 1%, the safe transient pins k_e to ~0.06; combined with daily sampling for two years that reaches only **0.33** — short of 0.20. Reaching the target by this route needs k_e within ~0.04 of exact (L₀ essentially exact) *and* ~300 draws, which is not feasible. So tightening the load helps, but does not open a marker-only path to 0.20.

**If no calibrated flux assay exists, the deliverable is set identification, not a point estimate** (R16): the record pins a one-sided ceiling **θ ≤ 30.7% of normal** (balance + non-negativity + K_m), tightening to **22–25%** with the k_e prior, and it pins K_m (~6%), the clearance slope b = θg′+k_e (~28%), and a k_e-dominated combination (~8%) — but not θ itself, nor the untreated flux or fraction oxidised, which all inherit θ's ~3.0. A defensible ceiling beats an unreachable target.

*(Monitoring, claim 2 above, is detailed with its stress tests — minimum programme, detectable change, drift and antibiotic breaking points — in R-D.)*

### The therapeutic decision the paper was written for — and what it actually needs (round 11, R-J)

All of the above prices how well θ can be *measured*. The clinical act the therapy actually forces is a different question: **which dosing interval to book.** That is a trough condition on capacity, not a measurement of θ. The therapeutic increment Δ(t) adds to capacity in (θ+Δ)·g(m); on a train of decaying doses it sags to a steady-state trough just before the next dose (paper §7: *"a dosing interval is chosen to protect the trough"*), and the safety constraint binds there:

> θ_base + Δ_trough(τ) > θ_crit.

This inequality is the whole decision, and it is worth seeing before the measurement questions that price it. The binding case is the longest interval anyone would book — the **21-day step** — against the parent paper's intrinsic critical capacity (θ_crit^int = 151 µM/d). The trough increment depends on one quantity, the **PCC protein half-life**, and the decision is a single crossing in it:

![the 21-day dosing step as a safety margin crossing zero at the PCC protein half-life](fig/r18_knife.png)

The y-axis is the **safety margin** μ = Δ_trough(21) − θ_crit, so the decision is a single **zero crossing**: μ > 0 (**right** of 3.5452 d) books the 21-day interval with **no knowledge of capacity** — the drug trough alone clears θ_crit for any residual θ ≥ 0 — and μ < 0 (**left**) needs θ, which is unmeasurable (the rest of this section). The repo's assumed 3.5 d sits just below zero (margin −8 µM/d); the one measured human value — **3.6 d** (Chandler & Ballard 1985) — sits just above (margin +10 µM/d). The two differ by **0.055 d — 1.5 % in half-life** *(at Model C's K_c = 100; across its unmeasured shape constants the crossing itself moves about 3.29–4.08 d over the corners of the admissible box, 3.30–3.97 d over 300 sampled sets, R-S)* (the only percentage on the figure is an x-axis, half-life difference, not the vertical margin gap). The compartment-anchored estimate (~4.3 d, R-P) with its 3–6 d bridging band straddles the crossing. That is the whole project in one image: **the decision does not need capacity — it needs this half-life, and the measured value sits ~1.5 % (in half-life) past the threshold that flips it.** Everything else in the report is how that figure was earned — the measurement questions above (why capacity itself is out of reach) and the interval arc that follows (R-J → R-P), which reduces the booking to λ and then measures it.

*Why this figure leads, and not the identifiability pencil of R-Q: the knife-edge's x-axis is the only quantity in this report that exists and has been measured, while the pencil's axes (θ, k_e) are exactly what the report concludes you cannot and need not measure.*

**Two questions, and the per-patient one is operative.** The intervals *are* far apart in capacity: at the repo's constants Δ_trough(14) = **400 µM/d** (the formula reproduces the repository's own dose train and its reported trough fraction 0.134·θ exactly), and adjacent grid choices differ by **287 µM/d (14 vs 10) and 600 (14 vs 7)** at the same delivered rate — so the *grid-resolution* tolerance is loose, sd(ln θ) ≈ 10–20 under the routine-load reading, which the record (3.05) clears. **But that is a fact about the grid, not about a patient.** For a per-patient booking the operative number is whether that individual's choice is pinned, and **a PA population straddles the single safety boundary θ_crit − Δ_trough(14).** Simulated over the whole posterior (the 3.05 already bundles the K_m fibre, the k_e prior, the load bias and the assay — so this *is* "doing it properly, not just on the fibre"): at the record precision per-patient classification is **61 % correct, 0 % unambiguous** (no patient is far enough from the boundary that no measurement could change the booking — barely above the ~52 % base rate), and only at the managed 0.29 does it reach **93 % correct with a 72 % determined subgroup.** So for an individual the interval choice is **not** settled at the record precision; the "no new measurement" reading is true of the grid spacing, not of the booking.

**Model C determines the decision by removing it — a weaker claim than it sounds.** The parent paper's intrinsic critical capacity (`prop:intrinsic`, no m_tol convention) is θ_crit^int(s_base = 600) = **150.8 µM/d**, *below* Δ_trough(14) = 400. The 100 %-determined result is therefore **not** "the boundary is well resolved" — it is "no patient is anywhere near the boundary," because the drug trough alone clears θ_crit for every residual capacity. That carries a **testable consequence**: Model C predicts the uniform fortnightly regimen is adequate for *every* PA patient, checkable against the trial's own marker reductions (the ~70 % C3 fall reported for mRNA-3927). If any patient is inadequately controlled at the booked interval, the intrinsic margin is too low or λ is shorter than assumed.

**λ is the sixth unmeasured constant, and its tolerance is the specification the report ends on.** *(⚠ Revised by R-N, round 15: λ is in fact **measured** in humans (3.6 d, Chandler 1985) and **transferable** as an intensive rate — so it is not a sixth *unmeasured* constant but the pattern's boundary. The specification framing stands, and sharpens: the measured value sits right on the thresholds.)* The decision does not turn on sd(ln θ) at all — it turns on λ, the human decay rate of the increment, which `design/run.py` marks **ASSUMED** ("nothing measures it") on an allometrically-scaled rodent PK model. Δ_trough(14) swings a factor of **42** (47 → 2000 µM/d) across a 2–7-day half-life, so "the drug term dominates θ_crit" is true at one end and false at the other. That makes the round-11 conclusion conditional on an unmeasured human constant — **the same pattern as κ, c, K_e, the Model-C shape constants and φ (R-I); λ is the sixth instance, and the first that gates the *decision* rather than the *estimate*.** As a specification, like the 19 % flux and the 30 % K_e:

| headline | holds while | half-life tolerance | λ tolerance |
|---|---|---|---|
| Model C: every patient clears τ = 14 | Δ_trough(14) > θ_crit^int = 151 | half-life > **2.62 d** | λ ≲ +34 % of assumed (~25 % margin) |
| routine-load coarse grid determined at the record 3.05 | req sd(14 vs 10) ≥ 3.05 | half-life > **3.32 d** | λ within **~5 %** of assumed (knife-edge) |

*(⚠ Round 19, R-S: 2.62 d and the ~25 % margin are at Model C's K_c = 100. At the worst corner of its admissible shape box the 14-day floor is 3.06 d, i.e. λ ≲ +14 % of the assumed 3.5 d, a ~13 % margin, and the 15-day floor is 3.21 d, ~8 %; `t53`.)*

So the fold-1.0 coarse-grid "no measurement needed" is *fragile* to λ (θ* ≈ 30 sits on the boundary); the robust version is Model C's ~25 % margin (~13 % at the worst corner of its shape box, R-S). **What would measure it:** the mRNA-3927 programme did dose selection by population PK/PD modelling, so the trial's own PK data should contain λ — the decay of PCC activity / the marker response between doses. The report's conclusion is thus not a tolerance on θ but a specification for a pharmacokinetic measurement. **And round 12 (R-K) shows that measurement is also *in the accumulated on-therapy record itself*:** two years of the dated dose schedule self-calibrate λ to sd(ln λ) ≈ 0.055 (monthly) — inside the Model-C spec and near the routine-load one — so the trial PK becomes a cross-check rather than a prerequisite. The decision's binding unknown is λ, and λ is recoverable *(⚠ round 19: to the 14- and 15-day specifications; for the 21-day step near the 3.5452 d crossing only conditionally, because one patient's record cannot place a value near 3.6 d on a side (R-R) and the crossing itself moves with Model C's unmeasured constants (R-S))*; what the passive record cannot do is individualise θ enough to re-choose the interval symmetrically (R-K).

### Supporting findings
- **The confounding.** ln θ and ln k_e are confounded (correlation −0.99 at 5%); k_e carries 98% of Var(ln θ̂). The likelihood already conditions on the per-draw recorded load L(t_k) and drug dose Δ(t_k), so the L–m covariation is inside the 3.07 (R9); the only load lever left is the systematic bias, worth little (R13).
- **No floor only when k_e is exact; otherwise there is a floor.** In one sentence: *with k_e exact the fall is 1/√n and crosses 0.20 near n ≈ 300 draws, but at the achievable sd(ln k_e) = 0.2 no sampling rate reaches even 0.30* (R12, R13). The 0.39 is not an n-floor in the first clause; the feasibility wall is the second. "k_e achievable to ~0.2" is itself not an external measurement — it is what a safe transient delivers when the unrecorded load L₀ is known to 5% (the repo's load prior); reaching the k_e ≈ 0.05 the target needs would require L₀ to ~1% (R-B).
- **Task 2 is closed as a negative result** (R14): at 24 patients with k_e shared the per-patient sd is 1.17 (still a factor 3.2 either way), and ω is unidentifiable from the available records under every likelihood tried.
- **The two-stage design is safe but short.** A mild, safe transient pins k_e to ~0.2 log units where it is nowhere near pinning θ, but k_e to 0.2 is not tight enough: the two-stage lands near 0.9, not 0.20 (R15).
- **The specification is untested, and largely untestable here.** All data is simulated from the one-compartment model; the monthly record is nearly blind to a second compartment and to saturable non-enzymatic clearance, though it can catch a time-varying load drift (R10). And the round-6 literature flags a concrete model-validity problem: the linear clearance k_e·m has the alternative route clearing ~84% of the load at 5% capacity, which conflicts with the PA literature (non-PCC disposal is low-capacity — patients accumulate); the repo's own Model B (saturating clearance) is the literature-consistent form (R-H).

## What changed across the rounds

Several round-1 claims were overturned by later rounds. The causes:

| claim (round) | corrected to | what caused the change |
|---|---|---|
| "the 0.39 is an n-independent floor" (R1 framing) | **not a floor**: ~1/√n with k_e exact, no positive asymptote (R13) | decomposing bias vs scatter and pushing n → 10⁵; both load terms are learned from the dose-modulated operating point |
| "the washout is the safe route to 0.20" (R1) | **the washouts that reach 0.20 add substrate past the 1.177× crisis threshold; inadmissible** (R15, F3) | making the safety threshold a hard constraint and sourcing it to `boundary.tex` / Model C |
| "the one-pass EKF is a bug" (R2 suspicion) | **linearisation error, not a bug**; dropped as unfair (F5) | the exact information identity at the truth + error growth with prior width |
| "the FIM misleads in both directions" (R2 suspicion) | **exact in linear θ; a lower bound only in ln θ where θ̂→0** (F6) | 150–300-replicate Monte Carlo against the FIM |
| best washout "0.237, fails" (R2, 40 replicates) | **0.199, straddles** (F3) | the verifier's 150-replicate re-run |
| "flux ≫ k_e ≫ load; flux reaches 0.20 alone" (R1/R2) | holds, **but only if absolutely calibrated**: √(Y²+Z²) ≤ ~19% (R11) | the generic-observable + (s, L) analysis, and the breath-test reparameterisation |
| "at n=25 the load scatter costs more than the bias" (R3 draft) | **reversed**: the bias is modestly more binding (R13) | the item-13 verifier caught a `core.fim` coupling that silently pinned the bias in the scatter-off run |
| "a dense excursion detects saturable clearance at ~100% power" (R3 draft) | **overstated**: that test is size-invalid at large N; only strong saturation is detected at a size-valid N (R10) | the item-10 verifier found the true MM model is rejected at ~100% too (errors-in-variables) |
| "no feasible design reaches 0.20 on the marker" (R3) | **split in two**: true for cross-sectional classification; **monitoring within-patient change reaches ~0.13** (R-D) | asking for the longitudinal contrast, where the shared nuisances cancel |
| "the calibrated flux reading is the only cross-sectional route" (R-B) | **also** a direct ex vivo PCC enzyme assay; the nuisance sensors (k_e tracer + load turnover) do **not** reach 0.20 even combined (R-G) | building the joint multi-observable FIM from the repo's own sensor panel — a gap in all prior rounds |
| "ex vivo PCC reaches 0.20 at its repo 20%" (R-G, v6) | **only if κ is calibrated**: the 20% is repeatability; κ (the mapping) is unreported, and uncalibrated the route reverts to 3.0 — the same problem as the breath test (item 1) | the κ question: κ is structurally the breath recovery factor, and the 20% covers reading repeatability only |
| "κ must be calibrated from known-θ patients" (v8, `t27`) | **κ need not be estimated at all** (**the user's ratio-route question**): if θ is defined as A_patient/mean(A_healthy), κ cancels by construction; the route reaches ~0.11–0.19 from a healthy control panel alone, no known-θ anchor (R-G step 1, `t28`). The anchor study I had costed is both unnecessary *and* practically unavailable (nothing measures θ to 20%) | the estimand is % of normal (scale-free), so the ratio cancels the population mapping — I had answered "how to calibrate κ" when the question was "need it be calibrated" |
| "the ex vivo ratio route reaches target, κ dissolved" (v9) | **only for an *in vitro activity* estimand; it fails for *whole-body capacity*, which is what the repo's θ is** (**the user's round-6 catch**): the ex vivo assay reads in-vitro activity, which predicts in-vivo θ only weakly (Manoli R²=0.327), leaving a per-patient scatter ~0.66 that does **not** cancel in a ratio — ~5× the 0.135 and above the gate (R-G step 2, `t29`). The cross-sectional verdict reverts: only a *direct in-vivo* read (the breath test) reaches target | v9 conflated in-vitro activity with whole-body capacity; the round-6 literature (Lehnert, Manoli 2021) shows the mapping is weak and unmeasured in PA, and per-patient scatter (unlike the population κ) does not cancel |
| "monitoring reaches ~0.13 on the marker (main result)" (R-D) | **a model floor under perfect nuisance sharing, empirically untested and probably optimistic** (**round-6 literature**): no within-patient CV is published; the Moderna n=50 cohort found the markers too variable to use past baseline; eGFR and carnitine confounders plausibly exceed the drift breakpoints (`t30`). And the **breath test is likely the better instrument** — it reads capacity directly, c cancels, one reading ≈ 104 draws (`t31`) | the 0.13 assumes the only non-shared variance is per-draw noise; the literature names real month-scale drift the contrast tolerates only in small amounts, and documents reproducibility for the breath test, not plasma C3 |
| "the analysis describes capacity across the disease range" (implicit, all rounds) | **it describes the *managed*, stable inter-crisis steady state only** (**the user's round-7/8 catch**): on a bounded clearance a 5% patient has a steady state only once the load is managed below θ+k_e (θ_crit = L − k_e moves with the managed load); untreated they diverge, a ~10% load cut restores stability, margin μ is thin (R-H, `t35`). Not "no steady state at 5%" — that overstated it | every brief number was on Model A; the ceiling is below the *untreated* load, but L is actively managed (protein restriction, gut suppression cut it 20–40%), so real patients are stable at their managed load |
| "Model B is uniformly harder to identify; A is the most favourable case" (parent paper, Remark 15 / Fig 13) | **holds for the joint conditioning, NOT for the marginal CRB on θ** (**a finding about the parent paper**): under B the enzyme share rises (0.16→0.37), so the *marginal, prior-regularised* precision on θ *improves* (3.05→0.29–0.92) even as the condition number degrades (the paper's direction) and the exact (θ,K_m)↔(k_e,K_e) exchange degeneracy appears (no-prior sd blows to 55). The two objects move oppositely — the relation Remark 3 left "not established" (R-H, `t36`) | the share gives θ more leverage while the ill-conditioning stays in the (prior-pinned) nuisance subspace; "fails a fortiori under B" is true for identifiability, not for the one parameter of interest |
| "the longitudinal contrast worsens under B, 0.133→0.167" (my v11) | **wrong — a bug; B's contrast is 0.088 (better)** (self-caught via the managed-load work): the Model-B contrast had failed to share the load-bias b between windows (b is a patient constant and cancels, as in R-D). Sharing it, B is better on *both* axes; there is no cross-sectional/longitudinal asymmetry | b not coupled between windows in the v11 Model-B contrast; the share mechanism moves both axes the same way |
| "no feasible design reaches target cross-sectionally" (all rounds, stated flatly) | **a Model-A statement; under Model B the gap is ~1.5× not ~15×, but the target is still not reached by any route the data support** (**the user's round-8 catch, then re-scrutinised**): A misses by ~15× (floor ~0.29 even with perfect nuisances); B at a managed point is 0.29. A tightened-k_e point touches ~0.19, but the sweep (`t40`) shows that holds only at a *low, unmeasured ceiling* (sd tracks the share, not the margin) and only with *K_e known to ~30%*, which no data fixes. So Model A's 16% share overstates the difficulty, yet the cross-sectional target remains unreached on the evidence | A's low enzyme share, not a deep obstruction, drives the *magnitude*; but reaching 0.20 under B needs a low ceiling and a known K_e — the third unmeasured constant (after κ, c) — so the sign of the conclusion survives the model change |
| "everything bottoms out at the unrecorded load L" (**the user's v2 framing**) | **k_e is the binding nuisance**; L matters only as the (inferior) route to k_e via the transient; a direct k_e sensor (5.4% → 0.45) beats it without touching L (item 3) | the IV ¹³C₃ tracer gives k_e directly; turnover (the load) adds ~nothing cross-sectionally |
| "Model C's free carnitine identifies θ (0.34) when the channel is known" (round-9, `t41`) | **slightly over-credited c; the 0.34 is from the channel-known *assumption*, not the carnitine reading** (self-caught via the plasma-compartment work): m-alone with the channel known already gives 0.34 (c redundant for θ), and in the channel-free regime c helps only marginally (20→16), eroded by the plasma-vs-pool partition to ~17.7 (R-I, `t43`). c never earns its keep | last round attributed to c what the channel-known prior supplied; the plasma defect (a fifth unmeasured constant) removes even the marginal channel-ID help |
| "identifying θ is the goal; the decision needs it" (implicit, all rounds) | **the therapeutic *interval* decision needs no new measurement of θ** (round-11, R-J, `t44`): it is a trough condition θ_base+Δ_trough>θ_crit in which the drug increment (Δ_trough(14)=400 µM/d) dwarfs a deep-PA residual θ; the coarse grid is determined at the record 3.05, and under Model C's convention-free intrinsic margin (θ_crit^int=150.8 < 400) it is determined for *every* patient at any precision. What it turns on instead is λ, which is unmeasured | separating the decision estimand from the parameter: θ is hard to estimate but the interval choice is a sign condition the drug term dominates |
| the dosing-trough formula uses a **fixed per-dose amplitude** D (the user's equation) | **the repository's §7 convention is same delivered rate**, D ∝ τ (so the mean increment D/(λτ) is held fixed) — a correction that reproduces §7's "87 % trough cut at 28 days" exactly and gives the narrower, conservative interval gaps (round-11, `t44`) | the paper compares schedules at matched delivered dose, not matched per-injection amplitude; both regimes keep the headline, fixed-D is looser |
| "the 3.05 record is the pre-therapy baseline" (round-12 premise, user) | **the 3.05 is *already* the on-therapy monthly record** (round-12, R-K, `t45`): the genuinely flat pre-therapy record is 8.56, and the 14-day dose train sampled monthly is 3.05 — the dose train is itself Prop-2 excitation, so most of the lift is bought before any events | the increment's phase-drift against monthly sampling already excites the operating point; `record()` carries `delta_at` by default |
| "λ must come from the trial's PK data (a prerequisite absent from the repo)" (round-11 close, R-J) | **the on-therapy record self-calibrates λ** to sd(ln λ) ≈ 0.055 monthly (→0.023 dense), inside the spec, so the trial PK is a *cross-check* not a prerequisite (round-12, R-K, `t45`) | the dated dose schedule's response to longer gaps identifies λ directly; it was never asked of the record before round 12 |
| "21-day lengthening is deliverable for a growing subgroup as the record lengthens (51–73 %)" (round-12, R-K `t45` item 3) | **not supportable from the passive record at any horizon** (round-13, R-L, `t46`): done per-patient with λ estimated from the same record, by NLS simulation, the margin estimate does not tighten with time (sd(μ̂) ≈ 140 µM/d flat 6 mo→5 yr); 21 d needs the active managed-load protocol | R-K used a representative, λ-fixed FIM sd for the whole cohort; FIM propagation is unreliable for this tail-dominated one-sided margin (delta vs MC disagree ~400 µM/d — the §7b skew pitfall), so simulation is required and it resolves to the pessimistic end |
| "λ is ASSUMED, from rodent PK scaled allometrically; the sixth *unmeasured* constant" (rounds 11–14) | **λ is measured in humans and transferable** (round-15, R-N, `t48`): Chandler & Ballard 1985 human fibroblast 0.19/d (~3.6 d); it transfers unscaled because protein k_deg is an intensive tissue-level rate, so allometry was never right — λ is the *boundary* of the unmeasured-constant pattern, not a sixth instance | a human measurement was found (two literature searches); and κ's transfer problem is extensive/tissue-weighted while k_deg's is intensive, so they have opposite verdicts from the same cell type |
| "λ half-life 3.5 d gives the decision comfortable margin" (implicit, rounds 11–14 at the assumed value) | **the measured human value (3.6 d) sits on the thresholds, and the correction runs toward them** (round-15, R-N, `t48`): 21-day safe-for-all flips at 3.55 d (round-13's 3.5-d "needs active protocol" → human 3.6-d "free for all", a ~1.5 % knife-edge); coarse grid clears 3.32 by 8 %; every interval headline flips inside the 3-fold measurement band (1.2–11.5 d), so λ is the dominant model uncertainty *(the "hepatic likely faster" clause reversed by R-P, round 17: mitochondrial-matrix comparator → hepatic ~4.3 d, at/above the thresholds)* | the old 3.5-d guess was luckily near the human value but sat just below the 21-day floor; the real hepatic value (unmeasured) is the single measurement the decision most needs, and the on-therapy record can read it (d′ ≈ 7–18 σ for 3.6 vs 7 d) |
| "human hepatic PCC likely turns over faster than fibroblast (≤ 3.6 d, toward the thresholds)" (round-15/16, R-N/R-O) | **reversed: hepatic is likely equal or slightly slower, central ~4.3 d, above the crossing** (round-17, R-P, `t50`): PCC is a mitochondrial-**matrix** enzyme, so the apt comparator is the liver-mitochondrial median (4.26 d, Kim 2012), not the whole-liver median (2.4–3.6 d); the fibroblast 3.6 d already sits in the slow tail of its own distribution. P(hepatic < crossing) falls from ~50 % to ~30 % | the earlier steer used the whole-liver median; the compartment-correct comparator runs the correction *away* from the thresholds (21-day step ~65–70 % free-for-all: 64 % on the wide band, 71 % on the tight, both at K_c = 100; 54–56 % / 56–59 % at the least favourable Model C shape constants, R-S), though the band still straddles |
| "the on-therapy record reads λ at 7–18 σ, so it settles the two sides of the 3.5452 d crossing a fortiori" (round-16/17, R-O/R-P) | **the a-fortiori runs the wrong way** (round-19, R-R, `t52`): the 7–18 σ is for 3.6 vs 7 d (ln gap 0.665); the crossing is 43× closer to 3.6 d (ln gap 0.015). Per patient the record places a value near 3.6 d on a side at only 0.2–0.6 σ, even after five years of weekly draws. It settles the side only if the true value is clear of the crossing (near the 4.3 d central estimate, an 80 % chance within a year of monthly draws or 6 months of weekly ones, against the line at K_c = 100; against the worst-corner line of R-S, 4.08 d, at most about 47 % even after five years of weekly draws). A pooled estimate needs 23–77 patients at 50 % power (46–155 at 80 %) and a λ shared by all of them, which is unmeasured. | overstated **positive**: the 7–18 σ answered a 3.6-vs-7 d question and was carried to a crossing 43× closer; caught by the round-19 completeness critic |
| "the measured 3.6 d clears the 21-day crossing by ~1.5 %" (round-15 onward, R-N, lead figure) | **true only at Model C's K_c = 100** (round-19, R-S, `t53`): K_c, M_i and ρ₀ are unmeasured. Over the corners of the repo's own admissible box the 21-day crossing moves about 3.29–4.08 d (3.30–3.97 d over 300 sampled sets). 3.6 d clears it in 80 % of sampled sets and at 6 of 8 corners; one constant at a time it fails at K_c ≤ 50 or M_i = 150. The 4.3 d liver estimate clears it everywhere in the box. The 14- and 15-day crossings stay ≤ 3.06 / 3.21 d at the worst corner, so those conclusions hold for any liver half-life above those values. | overstated **precision**: the crossing was computed at one setting of three unmeasured constants; found while answering the owner's question about which model the conclusions rest on |

**A note on the pipeline this table records.** The reversals run in *both* directions. Early rounds overstated **positives** — the ex vivo "0.200" (κ assumed known), the v9 ratio route (estimand conflated), "a k_e sensor reaches the gate" (a low-ceiling point), and "free carnitine identifies θ" (the channel-known prior, not c, round-9→10); one round overstated a **negative** — the v11 "contrast worsens under B" (a sharing bug); and round 12 overstated a **positive** — "21-day lengthening is deliverable for a growing subgroup" — which round 13 reversed on finding the one-sided margin does not resolve under simulation (R-L). That last reversal turned into a **method lesson** (R-M, round 14): a FIM is reliable for a log-scale sd but not for a natural-scale or one-sided decision on a heavy-tailed posterior; the two propagations of that margin disagreed by ~450 µM/d. The audit it prompted re-checked every conclusion that takes the same step — the two load-bearing ones (round-11's determined-subgroup, the cross-sectional headline) survive by simulation, and no hidden favourable correlation feeds an un-simulated optimistic claim. Adversarial verification, or a question that forced a re-check, has caught each error, in both directions, and the corrections (not just the conclusions) are logged here. That the errors go both ways is the reassuring part: the process is not biased toward a preferred answer, positive or negative — it is biased toward whatever survives an independent recomputation. The standing discipline (repo values only, verify every claim, keep this record, read log-scale sds from the FIM and natural-scale decisions from simulation) is what converts a sequence of confident-but-wrong intermediate claims into a conclusion that has been attacked from both sides and held, except where round 19 made the 21-day result at 3.6 d conditional. Round 19 added two reversals, one against R-O/R-P's text and one against the precision of the lead figure. Rounds 16–17 overstated a **positive**: that the on-therapy record, which reads λ at 7–18 σ for 3.6 vs 7 d, settles the side of the 3.5452 d crossing a fortiori; per patient it places a value near 3.6 d at only 0.2–0.6 σ, even after five years of weekly draws (R-R, `t52`), which a completeness critic run on `story_so_far.html` caught and an adversarial check confirmed. And rounds 15 onward stated the 21-day result with too much **precision**: the ~1.5 % by which 3.6 d clears the crossing holds only at Model C's K_c = 100, and over the corners of the repository's admissible shape constants the crossing moves 3.29–4.08 d (R-S, `t53`), found while answering the owner's question about which model the conclusions rest on.

## What the repository contains, and where the brief departs from it

| the brief says | the repository has | consequence |
|---|---|---|
| "~24 patient records" (and, in the follow-up, "the real cohort is 24 patients") | `design.run.SEEDS = range(24)`: **24 draws of one simulated patient's intake sequence**, all with k_e = 5. There is no patient data anywhere in the repository (README: "parameter values are illustrative, not fitted"). | Nothing can be "fit jointly across patients", and ω cannot be estimated from them: their between-record spread of k_e is 0 by construction (`t9`). There are no complete-deficiency records (`t5`). |
| "steady state only" | The 3.07 carries the **booked dose train**: Δ(t), D = 6000 every 14 d, 3.5 d half-life, known, acting on the enzyme channel. Draws are quasi-steady, and the recorded load varies by 10% day to day (`t0`). | Without the dose train the 5% figure is **8.59**. With constant load too, the FIM is singular. |
| FIM over (θ, k_e, K_m, L) of ln m | Flux-space FIM over (ln θ, ln K_m, b, ln k_e). b is a load-recording bias with a 5%·L prior; K_m has **no prior**; σ combines a 5% assay CV and a 5% load-record scatter. | Equivalent to first order. The finite-difference ln m route agrees to 3.9e-8 (`t0`); per-draw ln m noise gives 3.045 against 3.067. |
| "m/K_m ~ 10 at normal" | m = 10, K_m = 100: **m/K_m = 0.1**, ε = 0.909. | Untreated, the enzyme runs at 9% saturation at normal and 35% at 5%. Under therapy it runs at 3–10% at normal and 4–30% at 5% (`t0`). |
| "figure caption: median and range across records" | `n4_theta_scope.py` → `tab_S_scope`. Each seed is scored by its own FIM **at the truth**; nothing is fitted. | |
| "upper bound ruling out capacity above ~25%" | `n6_profile_vs_crb.py`: the median upper end of a **profile-likelihood** 95% interval, on 6 records with the k_e prior in force: 708, 24% (`out/n6_profile_vs_crb.txt`). | Prior-driven, not set identification (Task 4). |
| — | `design.run.delta` **stops after 200 doses** (day 2786), so draws on n4's GRID beyond about the 93rd are untreated. | Affects any "longer monthly record" result. `t8` uses a continued train; Task 3's n·Var figures are corrected below. |

The fixed parameters are normal θ = 3000, K_m = 100 and k_e = 5 /d, with L = 322.73 held as capacity falls (convention A). At 5% (θ = 150), m = 54.0 and s = θg/L = 0.163. Time is in days.

## Scenario table: sd(ln θ̂) at θ = 5% of normal

**For one patient, the FIM is right about the likelihood in linear θ, and is attained there in every single-patient scenario checked (`t11`, `t12`).** In ln θ it is optimistic wherever θ̂ can reach 0: on the current record, and with k_e exact or to 5%. There, ln θ Monte Carlo summaries depend on how θ̂ ≤ 0 is bounded. **In cohorts it errs both ways:** it understates the ln θ RMSE but overstates the robust spread, because the other patients' positivity constrains the shared k_e (F6). **Use the FIM for ranking. Treat its absolute values at 5% as optimistic for RMSE unless a Monte Carlo column says otherwise.**

Monte Carlo verdicts follow one rule: **pass** if the interval lies wholly below 0.20, **fail** if wholly above, **straddles** otherwise. Intervals are 90% bootstrap unless marked 95%.

| scenario | FIM | Monte Carlo (MAP) | source |
|---|---|---|---|
| **current design**: 25 monthly draws, therapy on, k_e prior 0.8 | **3.07** [2.66, 3.26] | RMSE 7.1 and robust 6.6, both set by the θ ≥ 1e-3 bound (35% of fits sit on it) | t1, t2, t11 |
| … with no therapy ("steady state only") | 8.59 | — | t1 |
| … plus one untreated baseline draw | 2.31 | — | t5 |
| … with k_e exact | 0.388 | robust 0.39–0.45; RMSE 0.64–1.49 across three independent Monte Carlos (heavy lower tail): fails | t8 recheck, t11, t11 recheck |
| … with k_e exact and the load measured exactly | 0.204 | 0.215–0.218 (t11: 95% [0.192, 0.240]): straddles | t8 recheck, t11 |
| k_e exact, **365 draws every 2 days** for 2 years | 0.195 | 0.219 (robust 0.199) | t8, t8 recheck |
| k_e exact, **730 daily draws** for 2 years | **0.151** | **0.151** | t8, t8 recheck |
| **joint cohort, 24 patients at 5%**, k_e shared exactly (ω = 0) | **1.17** | RMSE 2.34, robust 0.75 | t9 |
| joint cohort, 24 patients, ω = 0.33 | 1.82 | RMSE 3.90, robust 1.39 | t9 |
| joint cohort, 64 patients, ω = 0 | 0.81 (0.78 on the Monte Carlo's own records) | patient 0: RMSE 1.33, robust 0.52 (200 replicates) | t2, t11 |
| **washout**: diet withdrawn, hourly for 48 h, k_e prior 0.8 | **9.4** (383 with no priors) | — | t7 |
| washout + the monthly record + unrecorded load L₀ known to 1% | 0.62 | — | t7 |
| washout at 0.5% CV, every 5 min for 48 h, L₀ to 1% | 0.195 | 0.199 [0.179, 0.219] (150 replicates; t10's 40 gave 0.237): straddles | t10, t10 recheck |
| load ×3, every 5 min for 72 h, 5% CV, K_m from the record (inadmissible load) | 0.193 | 0.211 [0.191, 0.230] (150 replicates): straddles | t10, t10 recheck |
| load ×2, every 15 min for 72 h, 5% CV, L₀ to 1%, joined (inadmissible load) | 0.199 | 0.217 [0.201, 0.232] (300 replicates): **fails** | t10, t10 recheck |
| **drug burst** at the admissible peak (2850), every 15 min for 72 h, 0.5% CV, κ exact | 0.194 | **0.204 [0.185, 0.222]** (300 replicates): straddles | t10, t10 recheck |
| drug burst, peak 1500, every 5 min for 72 h (N = 865), 0.5% CV, κ exact | **0.124** | **0.131 [0.114, 0.149]** (150 replicates): passes | t10, t10 recheck |
| **direct enzymatic flux**: one untreated reading to 10% | **0.112** | **0.113 [0.110, 0.116]** (2000 replicates): passes | t5, t12 |
| … to 15% | 0.159 | 0.163 [0.158, 0.168]: passes | t5, t12 |
| … to 20% | 0.206 | 0.218 [0.211, 0.225]: **fails** | t5, t12 |
| … crossing of 0.20 | 19.3% | **18.4%** [17.9, 19.0] by RMSE; 19.35% by robust spread | t12 |
| … with the reading's recovery factor known only to 10% | crossing 16.5% | crossing **16.0%** [15.4, 16.6] | t12 |
| … three readings averaged, each to 15% | 0.097 | 0.097 | t12 |
| flux to 20% + k_e to 1.2 + load to 2% | 0.199 | 0.221, 95% [0.197, 0.244]: straddles | t5, t11 |
| k_e known to 5% (monthly record) | 0.445 | robust 0.49, RMSE 1.67: fails | t5, t11 |

## Which hypotheses held

| # | hypothesis | verdict |
|---|---|---|
| 1 | sd ≈ 0.17 at normal, ≈ 3.07 at 5%, band 82–627% | **Held exactly**: 0.1699, 3.0673 and 82.4–627.1%. It is on the dose-train footing, for one simulated patient. |
| 2 | U-shape with its minimum near 200–250% | **Held**: minimum 0.1226 at 212%. |
| 3 | corr(θ, k_e) → −1 as s falls | **Held**: −0.879 at normal, −0.991 to −0.992 at ≤10%. The high-capacity arm of the U is a different confounder: corr(θ, K_m) → +0.98. |
| 4 | amplification (1−s)/s = 5.1 at s = 0.16 | **Held** as the fixed-(m, L, K_m) derivative. The posterior's effective factor is β = −4.44 (−4.50 in the cohort limit with K_m known). |
| 5 | the binding constraint is k_e, not assay noise | **Partly wrong.** k_e is the dominant term, but exact k_e leaves 0.39. In that 0.39 the load record holds 62% and the assay 29% (`t8`). |
| 6 | elasticity −s/(sε + 1 − s) | **Held** to 3e-10. |
| 7 | m/K_m ≈ 10 at normal | **Wrong**: 0.1 (ε = 0.909). |
| 8 | the 24 records are fit one at a time | **They are not patients and are not fitted.** Each is scored by an FIM at the truth. |
| 9 | borrowing k_e across the cohort buys a lot for nothing | **Partly.** At N = 24 it buys 3.08 → 1.17 (FIM, for the same target patient) if k_e is genuinely shared, and 3.08 → 1.82 at ω = 0.33. It never reaches 0.20. |
| 10 | ~83% of the variance at the band edge is n-independent | **Held**: 84.3% exact, 83.9% from the sd(n) fit. The implied Var(ln k_e) is 0.638, which is the prior itself. |
| 11 | 15× sampling moves the edge only a few percent (82 → 75%) | **Held**: 82.4% → 75.1%. |
| 12 | set identification reproduces the ~25% upper bound | **No.** Non-negativity alone gives 30%, and only because the record pins K_m. The 22.4% (`t4`; n6's 24% on 6 records) is the k_e prior's lower tail. |
| 13 | direct flux > k_e > load, in value | **Held by FIM.** Monte Carlo confirms flux > k_e (`t11`, `t12`); no load-only scenario was run by Monte Carlo. |
| 14 | a complete-deficiency record gives k_e for free | **No such record exists.** In the parent model θ = 0 has no steady state. |
| 15 | λ → k_e because the enzyme is saturated | **Premise wrong**: the enzyme is not saturated. λ ≈ k_e at 5% only because θ is small, and θg′ is still 11–16% of λ. |
| 16 | FD regression: EIV attenuation and an MA(1) correlation of −0.5 | **MA(1) held**: −0.47 to −0.49 in the equation error. The bias is **not** classical EIV in m: noisy ψ sits on both sides and flips signs. **IV and bias-compensated LS do not fix it; generalised TLS with the MA(1) covariance does** (`t6c`). |
| 17 | the integral form is much better conditioned | **Only the 2-column saturated form, which cannot return θ.** My earlier "conditioning is not the issue" was wrong: the **noiseless** FD design is nearly collinear (cond 353 on TB-B, 4.8e4 on TB-A), and the noisy cond(X) of 26–104 hid that. |
| 18 | τ tension between derivative noise and aliasing | **At fixed duration, only for derivative estimators**: optimum 0.9–1.6/λ. For the efficient estimator at a fixed number of samples (N = 24) there is an interior optimum at about 0.3/λ. |
| 19 | conditioning depends strongly on the excursion | **Held**: cond[1, m, m²] goes from 76 to 18 as the peak m/K_m goes from 0.7 to 7.4. |
| F1 | 0.39 is an n-independent floor, so measuring L is mandatory | **No floor.** The fall is slow because of the load bias b (slope −0.22 to −0.33 against −0.5). An L measurement is necessary at monthly or weekly sampling, but **not sufficient** at 25 monthly draws (0.204; K_m must also be known). |
| F1 | the k_e share is 98% at 5% and ~83% at the band edge | **Held, as one curve**: 98.2–98.4% up to 10%, 93.5% at 50%, 84.3% at the band edge, 77% at 100%, 34% at 200%. |
| F2 | the plateau ≈ (1−s)/s · ω, so ω = 0.33 lands at ~1.55 | **Roughly.** Evaluated with (1−s)/s = 5.13, the formula gives 1.74 at ω = 0.33, against a computed 1.47: 14–26% high across ω = 0.2–0.5. Your ~1.55 lies closer. The posterior amplifies by β = 4.5, not 5.1, and the patient's own record shrinks k_e below ω. At ω = 0 the computed plateau is 0.38 (N → ∞, K_m known), or 0.43 at N = 1024, against your ~0.42. |
| F2 | estimate ω from the 24 real records | **Not possible: there are no real records.** 24 real monthly records could not estimate it either (confirmed under every likelihood). What would: each patient's own k_e measured to about sd 0.2, or every-2-day sampling in roughly 50–100 patients. |
| F3 | TB-A to TB-B: which factor? | **The excursion**: 63% of the change in log sd (no priors), then assay CV 19%, capacity 10% and sampling 8%. For NLS, the excursion is the one indispensable factor. |
| F3 | "×4 misses, ×6 passes" | **A sampling and duration edge, not an excursion threshold**: ×4 passes hourly over 72 h, or every 30 min over 48 h. |
| F4 | IV or TLS removes the FD bias | **IV and bias-compensated LS: no.** They need a CV of 0.01% or less. **Generalised TLS with the MA(1) covariance: yes**, matching NLS on TB-B. |
| F5 | the EKF degrading on TB-B is a bug | **One initialisation flaw, but not the cause.** P₀ ignored that m₀ is tied to the parameters, and tying it halves the TB-B RMSE (3.25 → 1.82). The dominant cause is one-pass linearisation error, which grows with the prior width: started at the truth it still gives 0.50 against a bound of 0.088. **Dropped** from the comparison. |
| F6 | the FIM misleads in both directions | **The FIM is right in linear θ.** "Optimistic at N = 1" compared it with a number set by the θ bound. "Pessimistic at N = 64" was a 60-replicate artefact. Only the core of the constrained MAP is narrower, and that is the positivity of the *other* patients' θ acting through the shared k_e (tested). |

---

## Follow-up results

### F1. The 0.39 is not a floor (`t8_floor.py`, rechecked)

![follow-up 1](fig/t8_floor.png)

**What holds the k_e-exact variance**, by Shapley share of Var(ln θ̂ | k_e) = 0.150 at 5%:

| factor | share | in your partition |
|---|---|---|
| load-record bias b | 29.2% | load record (b + scatter): **61.7%** |
| per-visit load-record scatter | 32.5% | |
| assay CV (per draw) | 29.1% | assay: **29.1%** |
| K_m | 9.3% | K_m: **9.3%** |

- **b and K_m have zero main effect.** With no per-visit noise the variance is 0. They act only through interaction with the noise, so order-dependent figures range widely: b takes 41% if removed first and 70% if removed after K_m.
- **Measuring the load exactly removes 72% of the variance** (0.388 → 0.204).
- Across 10–100% of normal the two σ terms hold 62–68%. b's share falls from 27% to 11%, and K_m's rises from 9% to 27%.

**The 0.39 against the number of draws, k_e exact, at 5%** (FIM; therapy continued past the repository's 200-dose cap):

| scenario | 25 monthly | 100 monthly | 600 monthly | 365 every 2 d | 730 daily | draws for 0.20 |
|---|---|---|---|---|---|---|
| k_e exact | 0.388 | 0.277 | 0.162 | 0.195 | 0.151 | 324 in 2 y, or 328 monthly (27 y) |
| + load bias b exact | 0.298 | 0.147 | 0.060 | 0.084 | 0.056 | 58 monthly |
| + load measured exactly (b and scatter) | 0.204 | 0.101 | 0.041 | 0.057 | 0.038 | 29 monthly |
| + K_m and b exact | 0.194 | 0.095 | 0.039 | 0.054 | 0.036 | 24 monthly |

- **No n-independent floor.** A fit of √(floor² + c²/n) over n = 6–600 gives a floor of 0.189, but it fails out of sample: it predicts 0.190 at n = 6000, where the FIM gives 0.059.
- **The load bias is what makes the fall slow.** With b at its prior, the local slope of ln sd against ln n is −0.22 (25–50 draws), −0.33 (400–600) and −0.48 (2400–6000): the prior on b is slowly overtaken by the data. The record barely updates b: its posterior/prior sd is 0.96 at 25 draws and 0.64 at 365. With b exact the fall is pure 1/√n.
- **So the answer to "mandatory or optional" is in between.**
  - At monthly or weekly sampling an independent L measurement is **necessary but not sufficient**. At 25 monthly draws, k_e exact plus the load exact gives 0.204 (Monte Carlo 0.215–0.218), and K_m must also be known: 0.194 with K_m and the load bias exact, or 0.133 with the load fully exact.
  - With near-daily sampling it is optional: 730 daily draws with k_e to 2% give 0.176 (FIM).
  - k_e to 5% fails everywhere in the grid.
- **Requirements for 0.20 (FIM)**, with the load as recorded unless stated:

  | draws | k_e needed | other conditions |
  |---|---|---|
  | 730 daily | 2.9% | — |
  | 365 every 2 d | 1.0% | — |
  | 100 every 7.3 d | 2% | load bias and scatter both to 2% |
  | 25 monthly | exact, 1% or 2% | K_m exact **and** the load bias b exact (scatter kept: 0.194 with k_e exact, 0.199 with k_e to 1%); or k_e 2% with b and scatter both to 1% (0.196) |
- **Caveat from the recheck: the CRB is not attained at the deciding designs.** Monte Carlo with k_e fixed gives:
  - 25 monthly draws: RMSE 0.643 against a CRB of 0.388 (9% of fits off by more than 2×).
  - 365 every 2 d: 0.219 against 0.195.
  - 730 daily: 0.151 against 0.151.

  So 365-every-2-day designs are at best on the edge of 0.20. Per-record CRBs there span 0.185–0.206, and 17% of records miss.
- **The k_e variance share, as a curve** (25 monthly draws, median over 24 records):

  | capacity, % of normal | 1–10 | 20 | 50 | 63 | 82.4 (band edge) | 100 | 158 | 200 | 316 | 627 | 1000 |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | k_e share of Var(ln θ̂) | 98.4–98.2% | 97.7% | 93.5% | 90.4% | 84.3% | 77.2% | 50.2% | 34.2% | 11.1% | 1.2% | 0.2% |

  Your 1 − (0.39/3.07)² = 0.9839 at 5% and the ~84% at the band edge are two points on this curve (0.9838 and 0.8425).

### F2. N = 24, and ω (`t9_omega.py`, rechecked by `t9_omega_recheck.py`: headline confirmed; the every-2-day and sample-size claims corrected)

![follow-up 2](fig/t9_omega.png)

**Task 2 read at N = 24** (FIM; target patient at 5%; the patient's marginal prior on ln k_e stays 0.8):

| ω | 24 patients at 5%: sd ln θ̂ | population sd ln k_e | patient's own sd ln k_e | mixed 5–20% cohort: sd ln θ̂ |
|---|---|---|---|---|
| 0 | **1.17** | 0.25 | 0.25 | 1.44 |
| 0.1 | 1.25 | 0.25 | 0.27 | 1.50 |
| 0.2 | 1.46 | 0.25 | 0.32 | 1.66 |
| 0.33 | **1.82** | 0.26 | 0.40 | 1.97 |
| 0.5 | 2.31 | 0.26 | 0.51 | 2.40 |

- **Monte Carlo joint MAP at N = 24** (60 replicates):

  | ω | RMSE | robust | fits more than 4.5× off |
  |---|---|---|---|
  | 0 | 2.34 | 0.75 | 8% |
  | 0.33 | 3.90 | 1.39 | 25% |

  The FIM understates the ln θ RMSE here and overstates the robust spread, as F6 explains.
- **Your plateau approximation**, √(0.39² + ((1−s)/s·ω)²), gives 1.10, 1.74 and 2.60 at ω = 0.2, 0.33 and 0.5. The computed plateau is 0.96, 1.47 and 2.07, so the approximation is 14–26% high. At ω = 0 the computed plateau is 0.38 (N → ∞ with K_m known) and 0.43 at N = 1024, against your ~0.42. The amplification actually applied is |β| = 4.50 rather than (1−s)/s = 5.13, and the patient's own record shrinks their k_e uncertainty to 0.90–0.98 ω.
- **ω cannot be estimated from the repository.** Its 24 records have one k_e value (5.0) by construction.
- **Would 24 real monthly records identify ω? No.** Simulated cohorts of 24 at 5%, each patient's data-only ln k_e summarised and ω profiled (4000 cohorts per cell):

  | design | true ω | typical 95% interval for ω | share excluding ω = 0 |
  |---|---|---|---|
  | 24 at 5%, monthly | 0 | **[0, 0.89]** | 1% |
  | 24 at 5%, monthly | 0.33 | **[0, 1.00]** | 4% |
  | mixed 5–20%, monthly | 0.33 | [0, 1.22] | 4% |
  | 24 at 5%, every 2 d | 0.33 | [0, 0.54] | 49% |
  | 24 at 5% monthly, + each patient's k_e measured to sd 0.2 | 0.33 | [0.21, 0.47] | 99% |

  The reason is that one monthly record gives each patient's ln k_e only to sd ≈ 1.06 with no k_e prior, three times larger than the ω being sought. The independent recheck simulated **raw** records rather than per-patient summaries. It used three likelihoods: θ ≥ 0 profiled, θ integrated over [0, ∞), and the exact restricted likelihood with θ unconstrained. It also ran a Bayesian posterior with a half-normal(1) prior. Its findings:

  - **Monthly: confirmed under every likelihood.** At a true ω = 0.33 the test rejects ω = 0 in 0–5% of cohorts. The posterior P(ω > 0.2) is 0.43 at a true 0 and 0.59 at 0.33, against a prior of 0.84, and never exceeds 0.95. Letting capacities vary between patients does not change this.
  - **The interval width depends on the likelihood; the conclusion does not.** The exact restricted likelihood gives upper ends of 1.05 and 1.17, so the Gaussian summary is about 15% too narrow. The θ ≥ 0 likelihoods give narrower upper ends (about 0.46 and 0.61), but those come from bias towards ω = 0, not information: every 2 days they cover the true 0.33 in only 63–70% of cohorts.
  - **Every-2-day sampling: refuted at 24 patients.** From raw records the test rejects ω = 0 in 32–53% of cohorts at a true 0.33, not ~50–60% as the Gaussian summary suggests. 80% power needs roughly 50–100 patients: 50–60% at 48, 90–100% at 96. t9's asymptotic "33 patients" is optimistic even inside its own model (72% at 33).
  - **Monthly power.** t9's "about 1200 patients" gives 77% power inside the Gaussian model; with θ ≥ 0 imposed, power is lower.
  - **A per-patient k_e measurement to sd 0.2: confirmed.** At a true 0.33 the interval excludes 0 in 100% of 40 cohorts, with coverage of 0.33 at 90%.
- **The repository's only information on ω** is brief R1's "probably a third or more", for "the nearest comparable disposal rate" in healthy adults; there is no 3-HP figure (`q14_hierarchical.py:19-30`, `PHASE8_SIMULATIONS.md:721-723`). `PHASE10_SENSORS.md:486-494` adds that R1's figure is for the marker's clearance, not the model's k_e/α, so the true spread could be larger. The only other figure is `params.py:363-367`, a 51% between-subject spread in the 3-HP *baseline* between mice of one inbred strain. The repository reads that as a lower bound on between-patient variability, and it is not a k_e figure. No number in the repository measures ω.

### F3. TB-A against TB-B, and whether 5% can reach 0.20 at all (`t10_factorial.py`, rechecked)

**Factorial attribution.** The 2⁴ FIM factorial reproduces both benches exactly (383.3 and 0.089). The change in log₁₀ sd from TB-A to TB-B is −3.63, attributed by Shapley share:

| factor | share |
|---|---|
| excursion (diet out → ×6 load) | **63.4%** |
| assay CV (5% → 1%) | 19.2% |
| capacity (5% → 20%) | 9.5% |
| sampling (1 h → 15 min) | 7.8% |

Main effects carry 99.3% of the variance across the 16 cells; the only material interaction is capacity × excursion. On the "k_e 0.8 + K_m from the record" rung the shares are excursion 48%, capacity 22%, assay CV 22% and sampling 8%.

**For NLS the excursion is the one indispensable factor.** Starting from TB-B and putting back a single TB-A setting:

| TB-B with … at its TB-A value | NLS RMSE | FIM |
|---|---|---|
| nothing (TB-B itself) | 0.100 (recheck, 200 replicates: 0.088) | 0.089 |
| capacity | 0.38 | 0.32 |
| sampling | 0.21 (recheck, 200 replicates: 0.171) | 0.17 |
| assay CV | 0.51 | 0.45 |
| excursion | **1.30, and every fit fails** | 29.7 |

The other three degrade NLS roughly in line with the FIM. Starting from TB-A, no single TB-B factor rescues it: RMSE stays at 3.3–4.1.

**A 5% patient on TB-B's other settings** (15-min sampling, 1% CV, ×6) has FIM sd 0.32 with no priors and 0.07 on the rung. **So capacity by itself does not block 0.20.** What blocks it is that the excursions able to reach 0.20 add substrate.

**The 5% design frontier.** I scanned 5280 designs:
- sampling every 5–60 min, for 24–72 h;
- assay CV 0.5–5%;
- excursions: diet withdrawal of 25% or 50%, load ×1.5 to ×8, or a drug burst at peak 450, 1500 or 2850;
- the full knowledge ladder, alone or joined to the record.

1422 pass by FIM. The cheapest passes per class, with an explicit cost order (safety class first, then knowledge that must be bought, then assay CV, then samples):

| class | best route to 0.20 (FIM) | Monte Carlo |
|---|---|---|
| **no added substrate** (drug burst, diet withdrawal) | **none at 5%, 2% or 1% CV**: 0 of the grid. At 0.5% CV: drug 2850 at N = 289 gives 0.194; drug 1500 at N = 865 gives 0.124; withdrawal 50% with L₀ to 1% at N = 577 gives 0.195. | drug 2850 0.204 [0.185, 0.222], straddles; drug 1500 **0.131 [0.114, 0.149], passes**; withdrawal 0.199 [0.179, 0.219], straddles (t10's 40 replicates gave 0.237) |
| load ×1.5 | not at 5% CV; at 2% CV with N = 865, joined, on intake seed 0 (passes in only 7 of 24 seeds) | — |
| load ×2, 5% CV | L₀ to 1%, every 15 min for 72 h, joined: 0.199 | 0.217 [0.201, 0.232] (300 replicates), **fails** |
| load ×3, 5% CV, no bought knowledge | every 5 min for 72 h (N = 865): 0.193, with 36 h of load; 0.205 if the load is held only 24 h | 0.211 [0.191, 0.230] (150 replicates), straddles |
| load ×4, 5% CV, no bought knowledge, joined | every 5 min for 48 h: 0.199 | — |

- **The earlier "×4 misses (0.205), ×6 passes (0.176)" is a sampling-and-duration edge** at hourly sampling over 48 h, not an excursion threshold. ×4 passes hourly over 72 h (0.182), and every 30 min over 48 h (0.174).
- **κ matters at small N.** At N = 289 the drug burst needs κ known exactly: with κ to 12% the 2850 burst gives 0.345. At N = 865 it passes by FIM with κ to 12% (0.182; no Monte Carlo).
- **Load excursions sit past the repository's own crisis threshold.** The parent model's fold is at 1.18× baseline load (`results/numbers/boundary.tex: fold_load_ratio = 1.177`, printed by `t0`, for a 10%-capacity patient). Every load-up excursion in this table is therefore beyond the point where the parent model has no steady state.
- **Answer.** At the clinical 5% assay CV, no design that adds no substrate reaches 0.20, even by FIM. At a research-grade 0.5% CV, the drug-1500 burst (865 samples, κ exact) passes by Monte Carlo. The drug-2850 burst and a 50% withdrawal with L₀ to 1% straddle 0.20. Load excursions are not admissible. Only five designs were checked by Monte Carlo, so the other FIM passes are unconfirmed. Two corrections to the earlier Task 7 text follow below.

### F4. One fix for the FD estimator (`t6c_iv.py`, rechecked)

![follow-up 4](fig/t6c_iv.png)

- **The noiseless FD design is nearly collinear**: ψ₀ is almost a quadratic in m along the path. In midpoint form, cond(X₀) is 4.8e4 on TB-A, 353 on TB-B and 1234 on bench C (5%, ×6, 15 min, 1% CV). The noisy cond(X) that t6 reported (about 25 on TB-B, 93–104 on TB-A) is inflated by the noise in ψ.
- **The psi column carries almost no independent signal.** Its partial signal-to-noise ratio, off [1, m, m²], is 0.0016 on TB-B and 7e-5 on C. So OLS's K_m coefficient goes to its noise-only limit: −266 on TB-B (exact limit −265.5, Monte Carlo −265.1).
- **Instrumental variables do not fix it.** I tried instruments from samples k−1 and k+2 (IV-JI, 2SLS, GMM with an MA-robust weight) and a heteroscedasticity-weighted IV. In the midpoint form they still return θ̂ ≤ 0 in 33–61% of records on TB-B and on C (up to 85% in the left-point form). The instruments look strong overall (first-stage R² 0.89) but are weak in the direction that identifies K_m (conditional R² 0.019).
- **Bias-compensated LS does not fix it either**, even with the exact lognormal moment correction, which is verified unbiased. The corrected X′X is positive definite in only about half the records.
- The unweighted linear corrections become consistent only at an assay CV of 0.01% or less. The weighted IV and weighted BCLS approach consistency at 0.1% on TB-B; on bench C everything needs 0.01% or less.
- **Generalised TLS with the row-wise MA(1) (tridiagonal) noise covariance does fix it** (GTLS-t, midpoint form).

  | | TB-B: bias / RMSE / robust | bench C: RMSE / robust |
  |---|---|---|
  | GTLS-t | −0.016 / 0.090 / 0.084 | 0.43 / 0.34 (recheck: 0.57 / 0.37) |
  | NLS | −0.011 / 0.079 / 0.085 | 0.37 / 0.31 (100 records) |
  | FIM sd | 0.089 | 0.32 |

  The left-point form is discretisation-biased. On TB-A nothing works, because the bench carries no information.
- **So:** OLS-FD and every linear correction of it is a clean negative result, with the mechanism quantified above. The FD equation itself is usable if it is estimated with its correct noise covariance. At that point it is a reparameterised NLS, and offers no advantage over fitting the ODE.

### F5. The one-pass EKF: one initialisation flaw, but linearisation is the cause; dropped (`t6d_ekf.py`, rechecked)

![follow-up 5](fig/t6d_ekf.png)

- **No algebraic error.** A copy of t6's filter reproduces t6's numbers exactly. On noiseless data at the truth, its posterior sd equals the FIM+prior bound to five digits (0.08922). The log-state Jacobian matches finite differences to 3e-4. Process noise is irrelevant: q = 0 gives the same result. And an iterated measurement update is a no-op, because the measurement is linear in the log-state.
- **One real inconsistency, which is not the cause.** The diagonal P₀ ignores that the initial marker is a deterministic function of the parameters. Tying it improves RMSE from 3.25 to 1.82 on TB-B.
- **The cause is one-pass linearisation error, which grows with the prior width.**
  - Started *at the truth* with t6's prior width, it still gives RMSE 0.50 while reporting sd 0.05.
  - Narrowing the prior to ×0.1, ×0.3, ×0.5 and ×1 gives RMSE 0.028, 0.093, 0.261 and 0.503.
  - Tied P₀, annealed parameter noise, lag-1 re-linearisation and a cubature filter all fail. The best robust spread is 0.66, against a bound of 0.088.
- **"Wider on TB-B" is mostly one record.** 83% of the EKF's sum of squared errors on TB-B comes from a single record. By robust spread it is *narrower* on TB-B (1.07) than on TB-A (1.63), where it simply returns its prior.
- **Decision: the one-pass EKF is dropped from the comparison.** It stays in t6's printed table, labelled as not a fair comparator, and is removed from t6's figure. The iterated EKF smoother is the Kalman entry, with a caveat: only its robust core matches NLS (0.082 against 0.085). Its TB-B RMSE of 0.328 is 4.2× NLS's because a few fits fail, so it too fails t6d's fairness criterion on RMSE.
- **The verifier's corrections.**
  - NLS started only from the drawn prior centre gives RMSE 0.085, identical to the truth-start in 40 of 40 records, so NLS is a like-for-like benchmark.
  - A one-pass Liu–West particle filter reaches a robust spread of 0.15–0.18, so the right statement is "no *Gaussian* one-pass filter is a fair comparator".

### F6. Where the FIM misleads, and the Task 5 spot-checks (`t11_fim_check.py`, rechecked)

![follow-up 6](fig/t11_fim_check.png)

- **Given (K_m, k_e), each θᵢ enters linearly.** So the FIM describes a Gaussian likelihood in **linear θ**, and with an sd of 3.0θ at N = 1 and 0.78θ at N = 64, θ̂ ≤ 0 is common. Every ln θ summary then depends on how that is handled.
- **N = 1.** 33.5% of fits sit on the θ ≥ 1e-3 bound. The robust spread is 5.4, 6.6 and 8.9 for bounds of 1e-2, 1e-3 and 1e-5, so t2's "6.6" measured the bound, not the data. Unconstrained, the linear-θ spread is 3.27 against the FIM's 3.04.
- **N = 64.** t2's RMSE of 0.49 was a 60-replicate artefact. With 200 replicates, patient 0's RMSE is 1.33 [0.55, 1.98], above the FIM's 0.78, and pooled over all 64 patients it is 1.47.
- **Why only the core is narrower than the FIM** (robust 0.52, about 0.67× the FIM from N = 8 to 64): positivity of the **other** patients' θ, acting through the shared k_e. The verifier tested this directly. Constraining only the other patients reproduces the full narrowing (0.556); constraining only patient 0 does nothing (0.704, equal to unconstrained). It also makes the MAP biased upward, by +0.04 to +0.12 in median ln error depending on the draws. The unbiased-estimator CRB does not bound a biased estimator, so there is no contradiction.
- **The Monte Carlo machinery was checked** on an exactly linear-Gaussian surrogate. It passes at N = 1 and 16. At N = 64 it falls just outside the 95% band on t11's draws and inside it on the verifier's independent draws.
- **Task 5 spot-checks**, 240 replicates each:

  | scenario | FIM | MC RMSE | MC robust | Monte Carlo verdict |
  |---|---|---|---|---|
  | flux reading to 10% | 0.112 | 0.118 | 0.113 | passes |
  | flux reading to 20% | 0.206 | 0.227 | 0.218 | fails |
  | k_e exact + load exact | 0.204 | 0.215 | 0.197 | straddles (95% [0.192, 0.240]) |
  | flux 20% + k_e 1.2 + load 2% | 0.199 | 0.221 | 0.202 | straddles (95% [0.197, 0.244]) |
  | k_e exact | 0.388 | 1.49 (bound hits) | 0.45 | fails |
  | k_e to 5% | 0.445 | 1.67 | 0.49 | fails |

  The two flux rows were re-run with 2000 replicates on fresh draws (`t12`): RMSE 0.113 at 10% and 0.218 at 20%, robust 0.111 and 0.206. The verdicts are unchanged.

  The ranking survives except for a tie at ranks 2–3. In linear θ the FIM is attained in every scenario (linear RMSE/FIM 0.99–1.06). With k_e exact or to 5%, the ln θ profile reaches θ near 0 in 28–38% of records on the verifier's independent draws (t11's own 24-draw profile set is skewed, P = 0.0014). So the FIM's symmetric ln θ interval is optimistic on the low side. Only the flux reading clears 0.20 by Monte Carlo. `t12_flux_mc.py` locates the crossing with 2000 replicates per point on fresh draws, checked to be representative: **18.4% [17.9, 19.0]** by RMSE and 19.35% by robust spread, against 19.33% by FIM. The FIM is attained in linear θ (linear RMSE/FIM 1.003–1.006). If the reading's recovery factor is known only to 10%, the crossing is 16.0% [15.4, 16.6].

---

## Round 3: conditioning, specification, and feasibility

### R9. What the likelihood conditions on (`t13_conditioning.py`)

The likelihood **conditions on the per-draw recorded load L(t_k) and the per-draw drug increment Δ(t_k)**, not on one constant load per patient. Each marker m(t_k) solves the balance at its own recorded load and dose (`m` recomputed from (load_k, inc_k) matches the record exactly); doubling the booked dose doubles Δ exactly. The only load quantity treated as unknown is a single systematic bias b (prior 5% of L, the one constant column of `flux_rows`).

Consequence: the L–m covariation is **already inside the 3.07**. Holding the load flat makes the record slightly *worse* (3.06 → 3.14), so the day-to-day variation is mild signal, not discarded noise; pinning the bias (measuring L exactly) moves 3.06 → 3.03, almost nothing at 5% because k_e dominates. Ideas that relied on unused L–m covariation are already counted; the only load lever is the bias (R13).

### R10. Is the model right? Residual diagnostics (`t19_residuals.py`)

**Caveat first:** all data is simulated from the one-compartment Michaelis–Menten model, so MM residuals are white by construction. This is a **power analysis of the diagnostic**, not a test of the model — a real test needs real data, which the repository does not have.

Residuals from the MM fit, regressed on marker / time / load, with power at the repository's noise (25 monthly draws, 5% CV) against three misspecified data-generating processes:

| misspecification | monthly record | load excursion |
|---|---|---|
| second compartment (slow exchange pool) | **blind** (identical steady-state balance; shows only in the transient decay, 5–15% of the rate) | — |
| saturable non-enzymatic clearance (k_e·m → Vₙₑ·m/(Kₙₑ+m)) | **nearly blind** (power ~5% at n=25; narrow on-therapy m-range looks linear) | carries the signal, but see below |
| time-varying load drift (10–20%, zero-mean) | **detected** (residual-vs-time power 88–100%) | — |

Two cautions the verifier earned. (i) The residual-vs-marker test is degenerate (the MM fit is itself a marker model). (ii) **A dense excursion's residual-slope test is size-invalid**: a noisy nonlinear regressor (errors-in-variables) inflates the false-positive rate, so at 5-min sampling (N≈577) the test rejects even the *true* MM model at ~100%. At a size-valid sampling (hourly, N≈49, true-MM rejection ~0%) only **strong** saturable clearance is genuinely detected (Kₙₑ=60: 64%); **mild** saturation (Kₙₑ=150) is not. So the load-record assumption is checkable from the monthly record, but a second compartment and mild saturable clearance are not — and even a rich excursion detects only strong misspecification, with a size-corrected test. **Model criticism and θ-identification need the same richer data, and neither is free.**

### R11. The flux observable as a specification (`t16_flux_observable.py`)

![round 3 flux contour](fig/r3_flux_contour.png)

Treating the observable as a specification, not a biological claim:

- **Generic observable O = c·θ·g(m), calibration c unknown to Z.** The 5×5 FIM collapses *exactly* (gap 3e-16) to a single calibrated reading at Y_eff = √(Y² + Z²). So the 0.20 contour is **√(Y² + Z²) ≤ 19.3% (FIM) / 18.4% (Monte Carlo)**. **Requirement sentence:** *an observable proportional to the enzymatic flux, read to Y% with calibration known to Z%, reaches sd(ln θ̂) = 0.20 at 5% of normal iff √(Y² + Z²) ≤ ~18–19%.* Averaging n readings divides Y by √n but **not** Z: a single shared calibration is the floor, so three readings still need Z ≤ 19% alone.
- **Reparameterisation (s, L).** If the reading is the fraction oxidised s = θg/L, then (s to Y%, L to X%) reaches target on the **same locus** (crossing Y ≈ 19.7%, weakly dependent on X): s and L enter as ln(flux) = ln s + ln L, so measuring L plays essentially the role of calibrating c. (Not *exactly*: X tightens the bias prior, which the 26 marker draws already partly constrain, so hypot(Y*, X) drifts 19.7 → 20.2% as X goes 0 → 5% rather than staying flat — the verifier's nuance.) Switching to the fraction oxidised buys nothing; it only moves the calibration burden onto an independent L.
- **The breath test.** It reads a fraction oxidised y = c·F(m) with recovery c unknown (round 2, `u2_breath_functional.py`). It can deliver the relative reading (a Y) but **not** the absolute calibration (its Z). With c entirely unknown, one reading is off the contour. It reaches target only if **both** its recovery is calibrated to Z ≤ ~18% **and** the load is known, with √(Y²+Z²+X²) ≤ ~18–19%. Whether such a calibration exists is a domain question, not answered here.

### R12. The k_e-precision / sampling trade-off surface (`t14_ke_sampling.py`)

![round 3 k_e × sampling](fig/r3_ke_sampling.png)

sd(ln θ̂) at 5% over (external precision sd(ln k_e), draws per year over a two-year programme), replacing the unobtainable "k_e exact" premise:

- The **0.20 contour** is reached only at k_e essentially exact: 309 draws/2 yr (≈155/yr) with k_e exact; **sd(ln k_e) = 0.05 does not reach 0.20 by 730/yr**, and ≥ 0.10 does not reach 0.30.
- At the **achievable** precision sd(ln k_e) = 0.20, sd(ln θ̂) stays **0.81–0.96** from 24 draws up to daily for two years.
- **There is no feasible point on the 0.20 contour.** Reaching it needs k_e tighter than ~0.05 log units, which the project's own estimate (~0.2) rules out. The "0.195 at 365 draws" of round 2 was conditional on k_e *exactly* known.

### R13. Does the fall asymptote? Load bias vs scatter (`t15_load_asymptote.py`)

![round 3 asymptote](fig/r3_asymptote.png)

- **n → ∞ with k_e exact: the limit is 0, not a floor.** sd(ln θ̂) falls as ~1/√n (sd·√n drifts slowly up toward ~4.6 as the priors wash out, log-log slope −0.497); 0.38 at n=25, 0.15 at 730, 0.03 at 2×10⁴. Neither the load bias nor the scatter creates a positive floor — both are learned from the dose-train modulation of the operating point. So 0.151 at 730 draws is **not** flattening toward a floor; **Task 2 is not answered by an asymptote, it is answered by feasibility** (R12).
- **The bias is modestly the more binding of the two load terms** (your premise holds). At n=25, pinning the *bias* (measuring L exactly) buys +0.094, against +0.073 for removing the per-draw *scatter* — the bias wins in 8/8 seeds and by Shapley (0.105 vs 0.084). _[An earlier draft reported the opposite; that was a bug — `core.fim` ties the default bias prior to `load_frac`, so a "scatter-off" run silently pinned the bias too. Caught by the item-13 verifier and corrected.]_ Neither is a floor.
- **Costing an L measurement.** Measuring L to 1–2% takes the n=25, k_e-exact value from 0.38 to ~0.30, but does nothing to the n → ∞ limit and almost nothing when k_e is the binding nuisance. L binds in three places — the 0.39 partition (small), the washout intercept L₀ (the biggest washout-ladder rung, but the washout is inadmissible on safety), and the flux route's calibration (R11) — and is only decisive in the third.

### R14. Task 2 is closed as a negative result

At 24 patients, all at 5%, with k_e genuinely shared, the per-patient sd is **1.17** (a factor of 3.2 either way) — still uninformative. ω is **unidentifiable** from the available records under every likelihood tried, and the repository holds no real records to estimate it from. The mechanism (the shared-k_e core is narrowed only through the other patients' positivity, which borrows little at ω ≥ 0; the between-record k_e spread is 0 by construction) is in F2 and Task 2. No further work is warranted; the recheck already committed stands.

### R15. The safety threshold as a hard constraint, and the two-stage design (`t17_two_stage.py`)

![round 3 excursion vs safety](fig/r3_excursion_safety.png)

- **The threshold, sourced and flagged.** `results/numbers/boundary.tex`: the parent model (Model C / ParamsC) loses its steady state when the load exceeds **1.177× baseline** (= 380). This is a **modelling construct** — a saddle-node of ParamsC's R0 setting — not a measured clinical crisis level; `src/model_A/run.py:35–40` states only the *ratio* (~1.5×) as "the physiological claim". The load ratio is model-independent (the fold *marker* 276 is Model C's, not Model A's — at the fold load Model A's marker is ~64). The repository's one clinical-labelled figure, M_TOL = 40, is a target, not a ceiling (the untreated 5% marker is already 54). The hard bound used here: an added-substrate excursion must keep the load below **1.177×**; withdrawals and the drug channel add no substrate and are admissible by construction.
- **The gap, as a figure.** Every sustained load excursion that drives sd(ln θ̂) below 0.20 crosses at **f ≈ 2.3× baseline** — about twice the admissible edge. The admissible band (shaded) never reaches 0.20.
- **The smallest excursion that pins k_e.** The relaxation rate λ = θg′(m) + k_e is k_e-*dominated* (≈1.13 k_e at 5%, the θ-term ~11%; R-A), not equal to k_e — so a *mild, safe* transient (even withdrawing 25% of the diet, marker dropping 54 → 47) still pins **k_e to ~0.2 log units** while leaving sd(ln θ̂) ≈ 1.2, because pinning k_e tolerates the small θ-contamination that pinning θ cannot. k_e's precision is capped by the unrecorded load L₀: ~0.2 at L₀ to 5%, ~0.06 at L₀ to 1%. Tightening k_e past 0.06 needs L₀ measured (R13).
- **The two-stage cost.** Safe k_e transient, then n steady-state draws for θ:

  | draws (2 yr) | k_e ≈ 0.20 (L₀ 5%) | k_e ≈ 0.06 (L₀ 1%) | k_e exact (ref) |
  |---|---|---|---|
  | 24 | 0.96 | 0.49 | 0.41 |
  | 52 | 0.93 | 0.42 | 0.32 |
  | 104 | 0.89 | 0.38 | 0.27 |
  | 365 | 0.81 | 0.33 | 0.19 |

  **The two-stage is safe and feasible but does not reach 0.20.** With a realistic L₀ it lands near 0.9; with L₀ measured to 1% it reaches ~0.25–0.30 only with hundreds of draws. k_e to ~0.2 is simply not tight enough (R12), and the safe transient cannot do better without also solving the load-measurement problem. The user's intuition — that pinning k_e needs a far smaller excursion than pinning θ — is **correct**, but it does not rescue θ.

### R16. Set identification elevated: what is recoverable (`t18_set_id.py`)

If no feasible design reaches 0.20, the deliverable is what the record *does* fix.

- **The one-sided bound on θ** (balance θg(m) + k_e m = L at the untreated state, L/m = 5.97/d): with K_m from the record, θ ≤ **30.7%** of normal under non-negativity alone, **25.3%** at the k_e prior's 95% low end. The lower bound is 0 — no record rules out a near-complete deficiency. The deliverable is a **ceiling on capacity**, not a point.
- **What the record pins, and what it does not** (median over 24 seeds):

  | quantity | sd(ln ·) | recoverable? |
  |---|---|---|
  | θ (capacity) | 3.05 | no |
  | a = θ·g(m)² | 3.13 | no |
  | untreated flux θ·g(m₀) | 3.09 | no |
  | fraction oxidised s | 3.08 | no |
  | K_m | 0.062 | **yes (~6%)** |
  | clearance slope b = θg′(m) + k_e (≈ λ) | 0.28 | **yes (~28%)** |
  | best (ln θ, ln k_e) combination (k_e-dominated) | 0.084 | **yes (~8%)** |

  The record identifies the **local clearance curve** (its slope and K_m), hence k_e, but not the capacity: the enzyme's own curvature a = θg² is swamped at this operating point. Neither the flux nor the fraction oxidised is a reparameterisation that escapes θ — they *are* θ, up to the pinned g and L.
- **Complete-deficiency check** (asked in Task 5, answered here): all 24 records have θ = 150; none is a θ ≈ 0 patient. If one existed at an untreated steady state, θg ≈ 0 would give k_e = L/m directly with no confounding — a tight prior at no cost. But in the parent model a θ = 0 patient has **no steady state**, so the counterfactual cannot be realised within the model.
- **Validation against the parent model.** The θ ≤ 22–25% ceiling here reproduces the parent model's published statement that the data "rule out capacity above about a quarter of normal" (`n6_profile_vs_crb.py`'s 708 = 24%, `out/n6_profile_vs_crb.txt`). An independent arrival at the same bound is a check that this pipeline matches its source, not just an internal result.

---

## Round-3 close-out

### R-A. The saturation regime, the integral collapse, and why the dynamic route fails (`t20_saturation.py`)

**A correction the brief's author flagged, and the explanation it unlocks.** The enzyme is **not** saturated in the disease range: at 5%, m/K_m ≈ 0.54 and g ≈ 0.35, so λ = θg′(m) + k_e = **1.13 k_e**, with the θ-term carrying **11%** of λ — of the same order as the precision we want, not negligible. The regime across capacity (K_m = 100):

| capacity | untreated m | m/K_m | g | on-therapy g (median) |
|---|---|---|---|---|
| 1% | 62.2 | 0.62 | 0.384 | 0.127 |
| 5% | 54.0 | 0.54 | 0.351 | 0.121 |
| 20% | 34.1 | 0.34 | 0.254 | 0.104 |
| 100% | 10.0 | 0.10 | 0.091 | 0.059 |

**This explains the TB-A integral collapse — and it is not quite the g = 1 story.** Only the *saturated* form (b) assumes g = 1 (a ~1/0.35 ≈ 3× misspecification at 5%), and it cannot identify θ structurally regardless. The general forms (b′)/(b″) use the proper g and in fact **recover θ noiselessly even on TB-A** (153.5 vs 150). Their TB-A collapse is **ill-conditioning**: θ's regressor is ∫g, and when g barely moves (0.21–0.35) the ∫g column is near-collinear with the constant and ∫m columns (design cond ≈ 2×10⁴), so the 5% assay noise is amplified ~10⁴-fold. On TB-B the ×6 excursion sweeps g to 0.75, the condition number falls ~7×, and θ survives the noise. So *sub-saturation via ill-conditioning*, not a misspecified flux, is what sinks the integral estimators.

**"Needs an excursion" and "needs to drive the enzyme toward saturation" are one constraint — and it is the one safety forbids.** As the load excursion f grows, the peak g rises from its baseline ~0.35 and θ becomes identifiable (sd and the ∫g-column condition both fall):

| f × load | peak marker | peak g | sd(ln θ̂) | admissible (≤1.177×) |
|---|---|---|---|---|
| 1.18 | 64 | 0.39 | 0.68 | yes |
| 1.5 | 83 | 0.45 | 0.32 | no |
| 2.3 | 131 | 0.57 | 0.20 | no |
| 6.0 | 364 | 0.78 | 0.14 | no |

θ first clears 0.20 at **f ≈ 2.3** (peak g ≈ 0.57, marker ≈ 131, well above K_m) — the excursion that drives the enzyme substantially toward saturation is the same one that reaches target, and it sits ~2× past the admissible edge. **The dynamic route fails for the same reason it would be unsafe: θ's information lives in the curvature of g, which only appears when substrate is pushed past the crisis fold.** This is the cleanest statement of the negative result.

**The transient route is closed in both directions, for one structural reason.** θ needs the enzyme driven *up* toward saturation (a large g-span), and g can only be raised by adding substrate — upward excursions large enough are inadmissible. The only safe excursion is *downward* (remove load), but the removable load is bounded by the gut-bacterial share (~20–50% of L), so a downward step caps g's span at ~0.13 and never approaches the ~0.22 θ needs (R-E). So no safe load change — up or down — reaches the g-range that identifies θ. That closes the transient route cleanly: it is not a sampling or estimator shortfall, it is that the identifying excitation and the safety constraint point in opposite directions.

### R-B. Sensitivity of the conclusions, and the flux contour across capacity (`t21_sensitivity.py`)

Because the data is simulated, the two load-bearing conclusions were stress-tested.

- **"No feasible point reaches 0.20" survives.** The anchor (sd(ln θ̂) at 5%, k_e to 0.2, daily draws for two years) is 0.72 at baseline and stays **far above 0.20** under every variation: K_m = 50 → 0.34, K_m = 200 → 1.46, assay CV 3–8% → 0.66–0.79, repo vs per-draw noise → 0.72, load scatter 2–10% → 0.59–0.86. Even shrinking the *achievable* k_e precision to 0.05 gives 0.27 — still above. The conclusion is governed by k_e precision, not by K_m, L, the noise model or the assay.
- **The ~19% flux contour survives.** The single-reading crossing Y* is 18.8–20.2% across K_m ∈ {50, 100, 200} and assay CV ∈ {3, 5, 8}% — a robust "~1 part in 5".
- **The flux contour across capacity** (the deliverable, extended): a single calibrated reading to **~19% recovers θ to 0.20 across the clinically relevant low-capacity range (~1–20%)**, where the requirement is essentially flat. It loosens only above ~50% (Y* 24% at 50%, >60% at 100%), where the monthly record alone already approaches 0.20 — so the reading matters most exactly where capacity is lowest. Three averaged readings relax the requirement to ~34% over the same range.
- **Provenance of "k_e achievable to ~0.2".** It is not an external measurement capability; it is what the safe transient (R15) delivers, and it is set by the unrecorded load L₀: L₀ to 5% (the repo's load prior) → sd(ln k_e) 0.21; L₀ to 2% → 0.10; L₀ to 1% → 0.06; L₀ exact → 0.04. The θ-target needs k_e ≈ 0.05, i.e. **L₀ to ~1%** — so the whole gap between the feasible (~0.2) and the required (~0.05) is the gap between knowing the unrecorded load to 5% versus 1%.

### R-C. The safety threshold cuts both ways

The 1.177× fold is a modelling construct (R15), so it must not carry a clinical conclusion in either direction. Stated cleanly:

- **Design requirement (a fact of this model):** the transient/excursion route needs a *sustained load excursion of ≈ 2.3× baseline* to reach 0.20 (R-A).
- **Admissibility (an open question for a clinician):** whether 2.3× is tolerable is not settled here; the repository's own crisis fold is at 1.177×, but that is a saddle-node of the parent model, not a measured safe limit.
- The gap between the requirement (2.3×) and the modelled threshold (1.177×) is about a factor of two. That is wide enough that the conclusion *probably* survives a corrected threshold — but it is stated as "probably", not assumed. A clinician's admissible limit anywhere below ~2.3× leaves the dynamic route infeasible.

### R-D. Realistic sampling, and longitudinal (monitoring) precision (`t22_longitudinal.py`)

**Sampling at an achievable k_e.** Holding sd(ln k_e) = 0.06 (what the safe transient delivers with L₀ to 1%, R-B), sd(ln θ̂) at 5% over a two-year programme is **0.47 monthly / 0.38 weekly / 0.31 daily**. The k_e-floor (|β|·0.06 ≈ 0.30) is approached only at daily sampling; a feasible monthly protocol sits at ~0.47. None reaches 0.20 — sampling buys about a third, not the factor of two needed.

**Longitudinal precision is the one regime that reaches target.** Every cross-sectional bar is an error on *absolute* θ. For monitoring, the quantity is the within-patient contrast sd(ln θ₂ − ln θ₁), where the shared k_e, load-bias b and K_m cancel, leaving only the independent per-visit noise:

> **One mechanism, not three escapes.** Every nuisance that defeats *classification* is a patient-level constant — the clearance k_e, the Michaelis constant K_m, the additive load-record bias b, and the two multiplicative calibration factors (ex vivo κ, breath recovery c). A within-patient contrast differences out *every* patient-level constant at once, additive and multiplicative alike. That single fact — not a separate piece of luck for each — is why monitoring survives the κ/c calibration failure that kills the ex vivo and flux routes, the k_e confounding that kills the record, and the load-bias that dominates the finite-n budget: the contrast never has to *know* any of them, only that they are stable between the two windows.

| draws per visit | sd(ln θ₂ − ln θ₁) | cross-sectional sd(ln θ) |
|---|---|---|
| 25 (monthly) | **0.26** | 2.1 |
| 104 (weekly) | **0.13** | 1.5 |
| 365 (~daily) | 0.07 | 1.0 |

The weekly contrast (0.13) is below 0.20 — **the only design in this whole study that reaches target on the marker alone.** Nobody had looked at this quantity before; it reframes the marker as a monitoring tool rather than a classifier. The stress tests (`t24_contrast_stress.py`) pin down when it actually delivers:

**Minimum programme — draws per window, not calendar span.** The contrast is set by the number of draws in each estimation window, nearly independent of how long the window spans. So two *short dense* windows work: two 6-month daily windows give **0.10**, and ~45 draws per window (≈ a 6-week daily burst) already reaches 0.20. Years of weekly draws are not needed; a longer gap between windows only allows more drift.

**Detectable change at 80% power** (two-sided 5%; minimum detectable |Δ ln θ| = 2.80·sd):

| programme | contrast sd | detectable rise / fall |
|---|---|---|
| monthly, 2-yr windows | 0.26 | +109% / −52% |
| weekly, 1-yr windows | 0.19 | +69% / −41% |
| weekly, 2-yr windows | 0.13 | **+45% / −30%** |
| daily, 6-month windows | 0.10 | +31% / −24% |

A 50% capacity rise (ln 1.5 = 0.405) needs a contrast ≤ 0.145 — weekly windows or denser.

**Breaking conditions.** The cancellation is exact only while the shared nuisances are stable between windows. The weekly contrast stays below 0.20 up to ~0.06 log-units of k_e drift and ~0.015 L of load-calibration drift; beyond that it degrades toward the cross-sectional bar (k_e drift 0.05 → 0.25, 0.10 → 0.43; b drift 0.02 L → 0.23). Physiology: k_e is the marker's clearance (renal + non-enzymatic) — over a week it is effectively fixed (< 0.01), but over a year in a growing child, where renal clearance scales with size, 0.05–0.1 is plausible; the load calibration b drifts with any diet, formula or adherence change. **So the contrast is reliable for short intervals (weeks to a few months) on a stable regimen, and degrades over a year or across a regimen change.** And it is destroyed by an antibiotic change between windows (R-E).

**The repository's own drift analysis (`PHASE11_DRIFT.md`, `PHASE11_PROPOSED_CONSTANCY.md`) supports the monitoring claim and undermines the long cross-sectional record.** These are the one inventory item that bears directly on a published conclusion, so they were read. They price the *constant-θ* assumption: a steady θ-drift of rate r over a record of span T injects a bias of order **−r·T/2 into ln θ̂**, which scales with the record's **calendar span, not its draw count**. The two implications cut in opposite directions, and both favour this report:

> - **It supports monitoring.** The monitoring result already insists on *short, dense windows* (draws-per-window, not span — the minimum-programme result above). Short windows have small T, so the drift bias is negligible precisely where the contrast is used; a decade of true θ-change is exactly what monitoring is *for*, measured as a sequence of short-window contrasts rather than assumed away.
> - **It undermines a long cross-sectional record.** The same bias makes a *single* θ estimate from a long record (e.g. a decade-spanning Tier-1 cross-sectional record) systematically wrong if θ is drifting at even ~0.5–1%/yr — the record's own length becomes a liability. The drift files price a **drift column** (letting θ vary linearly in the fit) as near-free insurance (≈ 0.3–1.6% precision premium). This is consistent with, and sharpens, R10 (the monthly record *can* detect a time-varying drift) and R-D: estimate θ over short windows, do not assume it constant over years.

This does not change the headline contrast numbers (short windows carry negligible drift bias); it confirms the design choice behind them.

**The empirical literature is hostile, and the 0.13 must be stated as an untested floor (round-6 research, `t30_monitoring_confounders.py`).** The 0.13 assumes the *only* non-shared variance between windows is the per-draw assay+load noise, and that k_e, b and the calibrations are shared exactly. The literature challenges both halves. **No published within-patient CV exists** for 3-HP, C3 or C3/C2 in stable PA. The one cohort that collected these markers serially — Moderna's n = 50 natural history study, every 6 months for 36 months — **restricted its analysis to baseline "due to high variability during follow-up,"** the closest empirical analogue to the monitoring claim, and it did not yield a usable longitudinal signal. And the named confounders map straight onto the drift axes already stress-tested: **declining eGFR drifts the marker's renal clearance k_e** (Shchelochkov 2021: C3 and 2-MC rise as eGFR falls), and **levocarnitine dose changes drift the C3 calibration** (C3 ∝ carnitine dose, r = 0.41). Re-running the contrast with these as between-window drift (`t30`): a k_e drift of 0.03 still passes (0.18) but 0.06 fails (0.28); a calibration drift of just 0.015 L already fails (0.20). Over a year, in a growing child on an adjustable carnitine regimen with unstandardised sampling, both are plausible. **So R-D's 0.13 is a best-case model floor, empirically untested; the monitoring claim holds only for short intervals on a stable regimen with standardised sampling** (time of day, fed state, timing vs dose and carnitine — none of which any PA cohort has reported). The most valuable outstanding item in the whole project is to **extract the raw within-patient CV** from Haijes 2020 (11 PA, longitudinal) and Cappuccio (serial C3): serial marker values from stable periods only, with sampling conditions, carnitine dose and renal function recorded, partitioned into within-window (the per-draw noise the 0.13 rests on) vs across-window (the non-shared drift that inflates the contrast). That calculation has never been done and the data exist.

**R-D against the only empirical serial-CV anchor — give it as a range, not 0.13 (round-10, `t42`).** The 0.13 assumes a ~5% per-draw marker scatter. The one measured serial within-patient CV for C3 is **26.4% (range 8.2–61.1%), in healthy adults** (Thompson 2012, same-day, unstandardised; 7.9% analytical, 39.2% between-subject) — there is still no PA figure. Recomputing the weekly (104-draw) contrast with the per-draw marker scatter set to the serial CV (the load bias stays shared and cancels):

| per-draw marker CV | weekly contrast | draws/window for 0.20 |
|---|---|---|
| 5% (repo assay) | 0.13 | ~46 |
| 8.2% (serial low) | 0.18 | ~81 |
| **26.4% (serial mean)** | **0.47** | **~590** |
| 61.1% (serial high) | 1.08 | ~3000 |

At the measured mean serial CV the contrast is ~0.47 and reaching 0.20 needs **~590 draws per window — not deliverable** (matching a variance ∝ CV²/n estimate). **Stated two-sided, honestly:** the anchor is healthy-adult, same-day, unstandardised data on a small labile plasma C3 pool, so it likely *overstates* the standardised stable-state PA CV; but PA's grossly elevated, larger C3 pool could also vary *less* in relative terms — the true PA serial CV is simply unknown. So the monitoring result is deliverable (≤0.20 at ~100–200 draws) **only if the true PA serial CV is near the low end (<~10%)**; at the one measured value it is not. A standardised within-patient CV in stable PA is the single number the monitoring claim turns on, and it has never been measured.

**The breath test is probably the right monitoring instrument, not plasma 3-HP (round-6 research, `t31_breath_contrast.py`).** Two facts converge. First, the one reproducibility claim anywhere in this literature attaches to the **¹³C-propionate breath test** (Shchelochkov 2021 "good intrasubject reproducibility"; Thompson 1990 "reproducible and prognostically useful") — *not* to plasma C3, whose within-patient CV is unmeasured and, by the cohort evidence, large. Second, the breath recovery factor c is subject-specific — fatal for cross-sectional classification, but it **cancels in a within-patient contrast**, like every other patient constant. The breath test reads the fraction oxidised S = θ·g(m)/L, whose elasticity d ln S/d ln θ at 5% capacity is **0.89** (repo model). So the contrast is √2·CV_S/(0.89·√n):

| CV_S per reading | n = 1 | n = 3 | n = 6 |
|---|---|---|---|
| 5% | 0.080 | 0.046 | 0.033 |
| 8% | **0.127** | 0.074 | 0.052 |
| 10% | 0.159 | 0.092 | 0.065 |

At a plausible ~8% within-subject repeatability a **single** breath reading per window (~0.13) matches 104 weekly blood draws, and 3–6 readings beat it (~0.05–0.07, detecting a ±15–20% change). Unlike the ex vivo assay (which fails cross-sectionally on the in-vitro mapping, R-G/`t29`) the breath test reads in-vivo capacity *directly*, so it is clean on both fronts. **This is a material revision: if a within-subject CV_S is measured, the monitoring instrument should probably be the breath test.**

**Stress-testing the breath result the way R-D stress-tested the marker (`t34_breath_stress.py`).** Because this is now the project's strongest result, it gets the same scrutiny:

- **Provenance of the 8%.** It is *not* a measured within-subject figure — the literature gives only between-subject spreads (Manoli controls CV 8.9%, heterozygotes 12.1%; Thompson n=5) and documents within-subject reproducibility only *qualitatively* ("good intrasubject reproducibility," Shchelochkov; "reproducible," Thompson). So CV_S is a between-subject stand-in, swept 5–10%; the headline needs a measured within-subject CV_S to stand — the same evidence gap as the marker's missing CV, but here the qualitative claim is at least on record for the right quantity.
- **c-drift is the breath analogue of the k_e/b drift.** c cancels only if *temporally* stable, and subject-specific ≠ stable: c is the ¹³C recovery (gut absorption + oxidation), so it can drift with diet, gut flora, or an **antibiotic course** — the same events that corrupt the load record (R-E). The contrast tolerates a between-window c-drift of **d(ln c) ≈ 0.14–0.19 (~15–19%)** before the bias reaches 0.20 — wider than the marker's k_e/b tolerances, but not unlimited; c must be held as fixed as k_e and the diet.
- **Control spread (Thompson n=5) matters only cross-sectionally.** For monitoring it cancels; for an absolute reading the control-mean SE is spread/√n ≈ 0.04 at n=5, 0.02 at n=30 — a 5-control reference is thin but ~15–30 controls make it negligible.
- **It dodges the C3 confounders.** The oxidised fraction, normalised to measured CO₂ production, does **not** inherit eGFR or carnitine-dose drift (those move the accumulated marker pool, not the oxidation fraction). It *is* sensitive to fed/fasted state and activity — but those are controllable within a standardised test visit, unlike the uncontrolled month-scale drift of eGFR/carnitine.
- **Cross-sectionally it has the same calibration problem as every direct route.** S = c·fraction, so with c known sd(ln θ) ≈ 0.06–0.11, and with c uncalibrated it is worthless — the R11 contour hypot(CV_S, Z)/0.89. Its advantage is *entirely longitudinal*, where c cancels; the cross-sectional verdict (only a *calibrated* direct in-vivo read reaches target) is unchanged — but unlike the ex vivo route it carries **no in-vitro/in-vivo mapping gap** (`t29`): it reads capacity in vivo by construction.

(Overall caveat: CV_S is a literature-informed sweep and the elasticity 0.89 is from the repo model; the headline holds unless the true within-subject CV_S is well above ~10%.)

### R-E. Antibiotics as a load term, and as a safe excursion (`t23_antibiotic.py`)

Gut-bacterial propionate is a substantial, antibiotic-suppressible part of the load, so **L is not a patient constant**. The **core identifiability model does not represent antibiotic status** (the load is dietary-recorded with a constant bias; a diet prescription does not record antibiotic weeks) — stated as explicitly as R9. The repository's *audit* layer does characterise it: gut flora ≥ 22.2 ± 8.5% of production, metabolites falling ~41% (12–76%) under metronidazole on a 21–35 day cycle (`params.py`, `m2_antibiotic_cycle.py`, Thompson 1990), and `u15` even treats suppression as a k_e instrument. This is a **model extension**, built on those repo numbers. **The round-6 literature confirms these are literature-grounded, not assumed:** the 22.2 ± 8.5% and the 41% (range 12–76%) metronidazole fall are both Thompson 1990's own primary figures — the 22.2% being explicitly a *minimum* (gut-bacterial production by metronidazole suppression, n = 6), not a point estimate. The commonly cited "one-third to 50%" gut-share headlines are extrapolations; the repo correctly uses the primary minimum. (A caveat worth recording: propionyl**carnitine** (C3) excretion is largely *unchanged* by metronidazole — Burns 1996 — so the ~41% suppression applies to production/turnover and the organic-acid markers, 3-HP among them, not to C3; the antibiotic effect on the load record depends on which marker the record is.)

- **An unmodelled antibiotic swing corrupts the cross-sectional estimate.** Because θ is weakly identified, the near-singular direction amplifies the unmodelled load variation without bound. A within-record cycle biases θ **downward** (suppressed troughs read as lower capacity): even a mild 30% cycle pushes ln θ down by 1–2 and collapses a fifth-to-a-half of fits to θ ≈ 0; at the realistic ~41% depth most fits collapse. A *uniform* whole-record suppression is absorbed by the per-record bias b — so the harm comes from *variation*, as with the load drift in R10. The longitudinal contrast is corrupted too if antibiotic status differs between windows — and here it is catastrophic. A load drop the record misses makes the low markers read, at the nominal recorded load, as *higher* capacity; the linear-response fake shift (`t24`) is **+1.1 ln θ at 10% depth, +2.0 at 20%, +2.9 at the realistic 41%** — a 3× to 20× apparent capacity *rise*, 8–20× the 0.13 weekly precision. **The sign matters for safety:** the fake is always toward apparent *improvement*, and antibiotics are started around clinical concern, so the failure mode is **false reassurance exactly during deterioration** — a declining patient put on metronidazole would read as improving. It would mask a real decline and mimic a large improvement. **Antibiotic status must be recorded for either route**, and held fixed across any monitoring pair.
- **A known antibiotic step is a safe (downward) excursion — good for k_e, not for θ.** It adds no substrate, so it escapes the 2.3× safety wall. For **k_e** it is a moderate instrument if sampled fast enough to catch the ~4 h relaxation: hourly gives sd(ln k_e) ≈ 0.26–0.42 (deeper = better); daily sampling misses the relaxation entirely. For **θ** it fails in *either* direction (sd ≈ 1.5–2.9): the limit is amplitude, not direction — a downward step can remove at most the gut share (~20–50% of L), moving g only 0.35 → 0.22–0.28, far short of the g ≈ 0.57 that θ-information needs (that needs a large *upward* substrate load, g-span 0.22 at +100%). The low-m collinearity of θ·g with k_e·m is a secondary effect under that amplitude cap. So a documented antibiotic step is a **safe way to pin k_e** — the stage 1 of the two-stage design (R15), without the inadmissible substrate loading — but it cannot recover θ.

### R-F. Drug response, assay precision, and a falsifiable prediction (`t25_drug_cv_marker.py`)

**Drug response.** The drug is a *known* increment Δ on the enzyme channel (at the booked dose θ+Δ ≈ 14× θ — the therapy dominates clearance, which is part of why θ is so hard to see). A **recorded** dose change therefore does *not* fake a capacity change — it is the good dual of the antibiotic: the contrast stays unbiased, at a precision cost (0.13 at a fixed dose → 0.18 if the dose doubles between windows). So the contrast is specific to the patient's *own* θ, and its drug-response use is detecting a disease-modifying change in θ itself — resolvable at ≥ +45% / −30% (weekly) over a weeks-scale window that matches the contrast's validity. (If the delivered dose/κ were *uncertain* rather than documented, it would leak into θ like the antibiotic, so dosing must be recorded too.) **Protocol:** two windows of ~45–100 draws (a ~6-week daily burst, or weekly for longer), same recorded dose and diet, antibiotic status fixed, gap set to the expected response interval.

**Assay CV × draws — a noisier marker does not substitute.** The contrast variance goes as CV²/n, so the draws per window needed to reach 0.20 scale with the square of the assay CV: **~45 at 5%, ~100 at 10%, ~190 at 15%, ~340 at 20%, ~760 at 30%.** A less-invasive but noisier collection (dried spot, saliva) at 15–30% CV is not deliverable — the protocol needs a ~5% assay. The monitoring route trades on precision per draw, not on more draws.

**A falsifiable prediction (the most valuable thing left).** The 0.13 contrast rests on the shared nuisances cancelling, which is a property of the simulation. Its checkable consequence is how much the marker *itself* fluctuates within a patient on serial sampling. The model predicts serial within-patient **sd(ln marker) ≈ 0.12 (CV ~12%) when sampling is synchronised to the dose cycle** (load 10% and assay 5% combined), rising to ~0.64 (CV ~64%) if sampling is dose-phase-random (then dominated by the 14-day dose cycle). **The test:** on real serial-monitoring data sampled at a consistent time relative to dosing, within-patient marker CV should be ≈ 12%. If real patients fluctuate materially more than that (after removing dose-phase), the extra variance is either real θ movement — the "stable θ between visits" assumption fails — or uncancelled nuisance; either way the cancellation is optimistic and the main result with it. This converts the study from simulated to testable.

### R-G. The joint multi-observable FIM (`t26_multiobservable.py`; repo sensors `q15_sensors.py`, `u3_joint_fisher.py`, `PHASE10_SENSORS.md`)

The whole report had treated this as a one-observable problem. The repository defines a sensor panel (Phase 10), and the joint FIM over it had not been placed in this report — a real gap. Built here, **repo observables and repo precisions only**.

**Enumeration (ten candidates; five dropped).** The model carries only the marker state — it has no propionyl-CoA pool — so every observable that reads that upstream pool is *inexpressible* and dropped at step 1:

| # | observable | what it reads | repo precision (source) | type | status |
|---|---|---|---|---|---|
| 1 | 1-¹³C-propionate breath test | fraction oxidised S = θ/(θ+k_e(K_m+m)), ×unknown recovery c | CV swept 5–10% (not a repo-measured figure) | repeated single-timepoint | **dropped** (no propionyl-CoA state; S not a physical flux) — modelled anyway in `u3` |
| 2 | **ex vivo PCC activity** | A = κ·θ·exp(ε) → **θ directly** | **20%** (assay repeatability; κ, τ_map unreported in any human study) | tissue/cell procedure | **kept** |
| 3 | plasma 3-HP | the balance itself (= the marker m) | 5% assay + 5% load-record | steady-state record | kept (= the record) |
| 4 | 2-methylcitrate | propionyl-CoA × oxaloacetate | — | — | dropped (needs oxaloacetate) |
| 5 | **C3 / propionylcarnitine** | near-equilibrium on propionyl-CoA, scaled by free carnitine | — | — | **dropped** (needs propionyl-CoA) |
| 6 | C3/C2 ratio | propionyl-CoA vs acetyl-CoA | — | — | dropped (needs acetyl-CoA) |
| 7 | propionylglycine | GLYAT, glycine-dependent | — | — | dropped (needs propionyl-CoA) |
| 8 | **IV ¹³C₃-3-HP micro-dose** (the "IV microbolus") | CL_total = k_e + θ/(K_m+m₀) → **k_e** | **5.4%** (P.8's four micro-dose studies; α closure assumed free) | tracer procedure | **kept** |
| 9 | **whole-body propionate turnover** (isotope dilution) | R_a = L + b → **the load** | **10%** (a turnover study at the load record's CV) | tracer procedure | **kept** |
| 10 | timed urine + paired plasma | CL_renal = f_renal·k_e → k_e, **bound only** | f_renal "not found in any population" | timed collection | kept, as a lower bound |

**Classification (the three buckets you asked for).** (a) *Known inputs already conditioned on*: dose Δ and recorded diet. (b) *Re-reads of m*: **none** — the repo has no clean re-read. Your C3 example is **not** a re-read: C3 reads the upstream propionyl-CoA pool, which would be *independent* information, but the model cannot express it without a propionyl-CoA state. (c) *Independent, parameter-specific*: ex vivo PCC → θ, IV tracer → k_e, turnover → L+b, timed urine → k_e (bound). Each buys precision on **one** existing parameter; none is a novel combination, because the only observables that would read a different mechanism (the propionyl-CoA markers) are inexpressible here.

**The joint FIM at 5% of normal, monthly record** (`t26`, reproducing `q15`):

| observable set (repo precisions) | sd(ln θ̂) |
|---|---|
| record alone (plasma 3-HP) | 3.05 |
| + ex vivo PCC (θ, 20%, **κ assumed known**) | **0.20** |
| + IV tracer (k_e, 5.4%) | 0.45 |
| + turnover (L+b, 10%) | 3.04 |
| + IV tracer + turnover (**nuisance sensors only**) | **0.44** |
| + ex vivo PCC + IV tracer + turnover (all three) | **0.18** |

- **Only the direct enzyme assay reaches target — and only if calibrated.** Ex vivo PCC at 20% gives 0.20 *when κ is assumed known*. The **nuisance sensors do not**, even together: an IV k_e tracer (5.4%) + a turnover load measurement (10%) give 0.44, and even k_e-exact + L-exact stays ~0.29 at monthly sampling. **Turnover (the load) adds almost nothing** to cross-sectional θ (3.05 → 3.04): at 10% it is looser than the record's own 5% load-bias prior, and k_e, not the load, is the binding nuisance.

- **The κ problem (item 1), which the v6 table glossed.** Ex vivo reads A = κ·θ·exp(ε): the repo's 20% is the *repeatability* (ε scatter, `q15`), while κ (cell-to-whole-body mapping) is "neither … reported in any human study" (`PHASE10 1.2`, R4) — i.e. uncalibrated, structurally the breath's recovery factor. So the 0.200 assumed κ known, and it is knife-edge even then. Carrying κ free with a calibration prior Z = sd(ln κ), the contour is the same hypot as R11, **√(Y² + Z²) ≤ 0.20** (Y the repeatability; `t27`):

  | Z = sd(ln κ) | Y = 10% | 20% | 30% |
  |---|---|---|---|
  | κ known | 0.100 | 0.200 | 0.299 |
  | 0.10 | 0.141 | 0.223 | 0.315 |
  | 0.20 | 0.223 | 0.282 | 0.358 |
  | κ free (uncalibrated) | 3.05 | 3.05 | 3.05 |

  At the repo's 20% repeatability the maximum κ uncertainty compatible with target is **≈ 0** (repeatability alone already = 0.20); tolerating any κ uncertainty needs the repeatability pushed below ~15% first (Y = 10% → sd(ln κ) ≤ 0.173). **Uncalibrated, the ex vivo route reverts to 3.0** — the same unresolved calibration problem as the breath test. The κ problem is a *classification* problem only: κ is a patient constant and cancels in the longitudinal contrast, so monitoring is immune.

- **Calibrating κ is a short, costable study — because κ is a *population* constant.** It is pinned once, not per subject. From N patients whose θ is independently known to sd(ln θ) = σ_θ, each read to repeatability Y_cal, sd(ln κ) = √(Y_cal² + σ_θ²)/√N. To reach the sd(ln κ) ≤ 0.173 the 10%-repeatability deployment needs: **N ≈ 2 patients at σ_θ = 20%, ≈ 4 at 30%, ≈ 9 at 50%** (10% calibration reads). **Deliverable:** *κ calibrated to ≤ 0.173 from ~2 patients with θ independently known to 20% makes the ex vivo route reach 0.20 at ≤ 10% repeatability.* Two caveats sharpen it: (i) **complete-deficiency patients (θ ≈ 0) do not serve** — A = κ·θ → 0 for every κ, so a zero anchor carries no κ information; the subgroup needs *known, nonzero* θ, and the repository identifies no source for it (ex vivo itself is circular), making that source the real open requirement, not N. (ii) **The same study does not fix the breath test**, whose recovery c is *subject-specific* (`u3`), not a population constant — so κ is calibratable once and c is not, the one structural advantage ex vivo holds over the breath test.
- **The ratio route, step 1 — κ need not be estimated *if* the estimand is in vitro activity (`t28_ratio_route.py`).** The previous two bullets answered *how to calibrate κ* — but the prior question is whether κ needs calibrating. For an in-vitro-activity estimand it does not. Capacity is reported **as a fraction of normal**: the figure's x-axis and the gate sd(ln θ) ≤ 0.20 are both scale-free (θ = 3000 is V_max in absolute units, but every result is reported as % of normal and sd(ln θ) is scale-invariant). If θ is *defined* as the ratio A_patient / mean(A_healthy), then A = κ·θ gives A_patient/mean(A_healthy) = θ_patient/θ_normal and **κ cancels by construction**. No known-θ anchor is needed anywhere — only a panel of healthy controls, whose mean carries an error σ_healthy/√N that falls with N. The reading's precision is then **hypot(Y/√reads, σ_healthy/√N)** with *no* κ term (the contrast is formed in log space — ln A_patient − mean_j ln A_j, i.e. the ratio to the *geometric* control mean; first-order identical to the arithmetic mean):

  | | N=10 | N=30 | N=100 |
  |---|---|---|---|
  | patient read 20%, healthy spread 30% | 0.221 | 0.207 | 0.202 |
  | patient read 20%, healthy spread 50% | 0.254 | 0.219 | 0.206 |
  | patient read 10% (or 20%×4 reads), spread 30% | **0.138** | **0.114** | **0.104** |
  | patient read 10%, healthy spread 50% | **0.187** | **0.135** | **0.112** |
  | patient read 10%, healthy spread 70% | 0.242 | **0.162** | **0.122** |

  With the patient read to ~10% (a few averaged assays) and ~30–100 obtainable healthy controls, the ratio reaches **≈ 0.11–0.19** at 5% capacity **with no known-θ anchor**. The "calibrate κ from known-θ patients" study (previous bullets, `t27`) answered a question that **need not be asked**: κ cancels. (σ_healthy — the between-healthy-subject spread of ln A — has **no repo figure**; the 51% the repo carries is the 3-HP marker baseline between mice, not PCC activity, so σ_healthy is swept and flagged.) **But this is the answer for a reading of *in vitro activity*, and the next bullet shows that is not the quantity the model's θ denotes.**

- **The ratio route, step 2 — it fails for *whole-body capacity*, which is what the repo's θ is (`t29_mapping_scatter.py`).** The round-6 literature catches a conflation step 1 glossed. The ex vivo assay reads *in vitro residual PCC activity*; the parent model's θ is *whole-body* V_max in the mass balance `dm/dt = s − θ·g(m) − k_e·m`. These are **not the same quantity**, and the literature says in vitro activity predicts in vivo capacity only weakly: Lehnert/Suormala/Baumgartner call residual activity "neither a predictive marker for severity nor for outcome," and the **only same-subject in-vitro/in-vivo study (Manoli 2021, in MMA, ¹⁴C-incorporation not a direct PCC assay) gives R² = 0.327** (n = 11) — the field's view is that the relationship is "difficult or impossible to reduce to a single scaling factor" (tissue-specific expression, hepatic vs fibroblast flux, alternative disposal). Write ln(in-vitro) = a + b·ln θ + η; the per-patient departure η is what makes the in-vitro reading an imperfect proxy, and **η does not cancel in a ratio** (unlike the population constant κ). The residual floor is σ_η = √(1−R²)·sd(ln θ_between):

  | | between-patient sd 0.6 | 0.8 | 1.0 |
  |---|---|---|---|
  | R² = 0.327 (Manoli, measured) | 0.503 | **0.654** | 0.802 |
  | R² = 0.50 (optimistic) | 0.441 | 0.571 | 0.701 |
  | R² = 0.80 (very optimistic) | 0.299 | 0.380 | 0.462 |

  (sd computed at patient read 10%, N = 30 controls, healthy spread 50%; the record contributes negligibly.) At the measured R² = 0.327 and a between-patient spread of ~0.8 log-units (Thompson in-vivo oxidation mean 36.3%, range 3.0–91.1%), σ_η ≈ **0.66 — about 5× the 0.135 of step 1, and far above the 0.20 gate.** Even R² = 0.80 with a tight 0.6 spread gives ~0.27, still over. A **tighter healthy-control spread does not rescue it** — once σ_η dominates, more or better controls barely move the answer (0.16 spread → 0.649 vs 0.50 spread → 0.654). **So the ex vivo ratio route reaches target only if the clinical question is in-vitro residual activity; for the whole-body capacity the model carries, it does not.** The one measurement that escapes this is a *direct in-vivo* read — the ¹³C-propionate breath oxidation test — which reads capacity itself, not an in-vitro proxy.

  **The floor is borrowed, not measured — and a small study would settle it.** The honest conclusion is *not* "the floor is 0.66"; it is that **there is no evidence an in-vitro assay reaches target for whole-body capacity, and the study that would decide it has never been done in PA.** The 0.66 is indicative, from the one same-subject datum (MMA). For the ratio route to reach the gate, the mapping would have to explain **R² ≥ ~0.89–0.96** (σ_η = √(1−R²)·sd_between ≤ 0.20 at sd_between 0.6–1.0) — near-perfect prediction, versus the ~⅓ of variance actually measured. The decisive experiment is **paired direct ex vivo PCC activity and ¹³C-propionate oxidation in the same PA patients** (Manoli's MMA design, with a direct PCC assay added); **~20–30 patients would reject the route** if PA's mapping is as loose as MMA's (at N=20 the 95% upper bound on R² is ~0.65, below the ~0.9 needed), and it is the single small study that decides whether the ex vivo route is viable (`t29`).

- **Say the circularity out loud (item 2).** The κ-calibration study of the previous bullet is **practically unavailable, and the ratio route is what makes it unnecessary.** Its anchor — ~2 patients with θ *independently known* to 20% — requires a source of independently-known nonzero θ, which is exactly what does not exist: **this entire project exists because θ cannot be measured to 20% by any standing method.** Complete-deficiency patients give zero κ information (A = κ·θ → 0 for every κ). Genotype (PHASE12) is the obvious candidate anchor, and it does **not** serve: PHASE12_GENOTYPE_REVIEW places genotype as giving a *class / predicted range* (a plausibility band), **not a θ identity** to 20% — and the repo has no carrier data to calibrate even that. So the absolute (κ-anchor) route is blocked on a non-existent input; the ratio route sidesteps it by never needing an absolute θ at all.

- **The assumption, and why the scatter is the worse half (item 3).** Even granting the in-vitro-activity estimand, the ratio is unbiased only if the mapping κ is the same in chronically-overloaded patients as in controls. A *systematic* difference b = ln(κ_patient/κ_control) adds in quadrature, rmse = hypot(sd, b); at an achievable sd ≈ 0.11–0.14 a difference of only **ln-ratio ≈ 0.15–0.17 (κ_patient ≈ +16–18%)** reaches 0.20 (`t28`). But this was the *lesser* problem: a systematic offset is in principle correctable, and bounded. The **per-patient scatter** of the previous bullet is the binding one — it is irreducible (it does not cancel, and no amount of averaging or control-matching removes it) and, at ~0.66, an order of magnitude larger than the 0.15–0.17 systematic tolerance. The round-6 literature is what sharpened this: the real obstruction to the ratio route is not a constant offset between groups but the fact that an in-vitro number does not pin an in-vivo one within-subject.

- **Carriers: the line closes on biology, not just missing repo data.** A 50%-of-normal carrier anchor would, in principle, give a second calibration point and test the linearity of the ex vivo reading in θ. The repository has **no carrier data**, but the round-6 literature shows the idea fails for a deeper reason even with data: carrier PCC activity is **not uniformly ~50%** — it is **subunit-specific**. Wolf & Rosenberg 1978: *PCCA* (α-subunit) carriers sit near half-normal (mean 52%, individual range ~43–64%), but *PCCB* (β-subunit) carriers are **indistinguishable from controls** (β is made in functional excess, so halving its gene dose changes nothing; Ohura 1989). About half of all PA carriers are therefore enzymatically invisible, and even the α-group overlaps the control range. So carriers would serve as a calibration anchor **only if genotype-stratified to α-subunit families** — and carrier status itself is now established by molecular, not enzyme, testing precisely because of this. The line closes on the biology, independent of the repo gap.

- **The flux pair (breath + turnover).** Turnover cleanly supplies the L+b the flux route needs (θ = S·(L+b)/g); at 10% it is ample. The unresolved half is the breath recovery: with c carried free, one breath reading is worth nothing and `u3` needs ~70 readings at 10% capacity. So the pair reaches 0.20 only if the fraction-oxidised reading is itself calibrated (S to ~19%) — exactly R11's requirement, now with the load half supplied by a named assay.
- **Item 3 with turnover's real precision.** Turnover at 10% leaves cross-sectional θ at 3.04 (negligible), and the safe-transient k_e it supports is only **sd(ln k_e) ≈ 0.29** (not the 0.06 of R-B, which needed the load to 1% — turnover does not deliver that). So a real load assay does not rescue the two-stage design.
- **Monitoring.** No sensor improves the longitudinal contrast: it already cancels the shared k_e, b and K_m, so the tracer and turnover tighten what already cancels, and a per-visit ex vivo PCC at 20% gives √2·0.20 ≈ 0.28, worse than the record's 0.13. Sensors help classification, not monitoring.

**Does it overturn the conclusion? (item 6)** After the round-6 literature, the cross-sectional picture is **more negative, not less**. The ex vivo PCC assay does *not* give a clean second route after all: granting κ cancels (the ratio route), it still reads *in vitro* activity, and that predicts the model's *whole-body* θ only to a per-patient floor of ~0.66 (R-G step 2, `t29`) — far above 0.20. So for the whole-body capacity the model carries, the only measurement that reaches target cross-sectionally is one that reads **in-vivo capacity directly**: the enzymatic-flux / ¹³C-propionate breath reading, O = c·θ·g, *if* its recovery c is calibrated (R11). The marker-only and nuisance-sensor conclusions are unchanged (k_e tracer + load turnover = 0.44). The net: **no inference-through-a-proxy route reaches 0.20 — neither the mass-balance marker, nor the nuisance sensors, nor an in-vitro enzyme assay through the weak in-vitro/in-vivo mapping; only a calibrated direct in-vivo flux read does.**

**The model extension (item 2), now the main open item.** Five of the ten candidates are dropped for want of a propionyl-CoA state, and `PHASE10 1.1(b)` flags adding one as the structural change its sensor half is waiting on. Scoping it, repo-grounded where possible:

- **The repo does not price C3** (correcting the brief's premise): C3 is dropped as inexpressible, so "compute sd(ln θ) with C3" cannot be done on repo values — it requires building route (b).
- **New parameters** beyond (θ, K_m, k_e, L, b): a non-enzymatic disposal rate k_P for the pool P, the overflow kinetics P → m (form and rate), and a C3 mapping constant κ_c — a *third* uncalibrated multiplier. None is in any repo file.
- **Are they identifiable from existing data? No.** The record observes only m, downstream of P; P and its disposal split are unobserved. Adding k_P and the overflow kinetics with only m observed adds *unidentified* parameters to an already ill-posed model — it can only make θ worse. So the extension is worth running **only together with a P-observable (C3)**.
- **Does C3 then break the degeneracy, or move it?** Structurally it **moves** it: measuring P pins the pool (the genuine promise — C3 reads the enzyme's actual substrate), but the θ-vs-k_P confound on the P balance is the *same structure* as θ-vs-k_e on m, and κ_c is a new calibration unknown. It **breaks** the degeneracy only if (a) non-enzymatic propionyl-CoA disposal k_P is negligible or measured and (b) κ_c is calibrated — neither of which the repository supplies. A definitive number needs the structural model built with those priors; that is out of the current model and is the repo's own stated main open item. Not fabricated here.

---

### R-H. Model A vs Model B — the validity question under the whole brief (`t32_massbalance.py`, `t33_model_b.py`)

**This is not a side finding. Every number in this brief — 3.07, the 0.39 partition, s = 0.16, every contour, the monitoring floor — was computed on Model A** (`dm/dt = s − θ·g(m) − k_e·m`), whose alternative clearance k_e·m is **linear and unbounded**. The round-6 literature says that is the wrong structure for the disease range, and the repository itself carries the alternative (Model B, saturating: `k_e·h(m)`, `h = m/(K_e+m)`). So Model B is a **threat to the whole analysis**, not an appendix item, and this round recomputes the headlines under it, A and B side by side.

**Which model the parent paper used.** The paper's own identifiability configuration (`sections/03_model_A.tex`, `eq:local`; `symbolic/.../observational_equivalence.py`) **is Model A**: `m′ = L + L₀ − (θ+Δ)g(m) − k_e·m`. So every published degeneracy and precision figure — and every number here — is a Model-A result. The paper does **not** hide this: claim C23 states the degeneracy is "conditioned on clearance being affine over the observed range," and Model C (`sections/05`) is where it treats a bounded alternative route and the decompensation a null genotype cannot escape. The assumption is explicit, not silent — but the round-6 literature is what shows it bites inside the disease range.

**The first question is existence, not precision.** Under a bounded alternative route, total disposal θ·g(m) + k_e·h(m) rises to a **ceiling θ + k_e** as m → ∞. A steady state exists only if the load fits under it: L < θ + k_e, i.e. θ > **θ_crit = L − k_e**. Below θ_crit the marker **diverges — there is no steady state, which is what decompensation is.** Where θ_crit lands decides whether the steady-state FIM analysis even applies:

| alt-route ceiling k_e | anchored load L | θ_crit (% of normal) | steady state at 5%? |
|---|---|---|---|
| 400 (**ParamsB**) | 285.6 | −3.8% | yes (even at θ=0) |
| 200 | 279.2 | 2.6% | yes |
| 100 | 276.0 | 5.9% | **no** |
| 50 | 274.3 | 7.5% | **no** |
| 25 | 273.5 | 8.3% | **no** |

**ParamsB sets the ceiling at 400 — above the load (~286) — so a steady state exists even at θ = 0; ParamsB does not reproduce decompensation at all**, which already fails the severe-PA phenotype (flagged: ParamsB is **miscalibrated** for severe PA, contradicting the repo's own Model C, whose comment states "the alternative route alone peaks at 500 against a baseline load of 600, so a patient with θ = 0 cannot be maintained — that is the severity of null genotypes"). With a literature-consistent ceiling *below* the load (k_e < L − θ₅%), θ_crit lands above 5%.

**But θ_crit moves with the managed load, and that is the clinically right statement (round-8 correction, `t35`).** θ_crit = L − k_e depends on L, and **L is actively managed** — protein restriction and gut suppression cut the propionate load 20–40% (Thompson: gut share ~22%, ~41% metronidazole-suppressible; protein restriction cuts the ~52% amino-acid share). So a real 5% patient is alive and stable *because their managed load sits below the ceiling*. Concretely, at a low ceiling k_e = 100 (untreated anchor L₀ ≈ 276):

| managed load L | reduction vs L₀ | θ_crit (% normal) | margin μ (% normal) | 5% patient |
|---|---|---|---|---|
| 276 (untreated) | 0% | 5.9% | −0.9% | **diverges** |
| 250 | 9% | 5.0% | 0.0% | at the edge |
| 240 | 13% | 4.7% | +0.3% | stable |
| 220 | 20% | 4.0% | +1.0% | stable |
| 200 | 28% | 3.3% | +1.7% | stable |

So the fundamental finding, correctly stated: **the steady state the brief analyses is the *managed* inter-crisis state. Untreated, a 5% patient is below θ_crit and diverging; a ~10% load reduction (well within what management delivers) brings them above it.** The analysis does not describe the untreated or decompensating patient, and the margin μ = θ − θ_crit is **thin** — the record's own antibiotic and protein-intake excursions (paper §3.5, R-E) move L and can erase it, which is exactly a θ_crit excursion, i.e. a decompensation. The benign steady state is real for managed patients; it is not a Model-A artifact. (The genuine saddle-node is the paper's own, from its non-monotone carnitine α (`prop:intrinsic`); the ceiling hypothesis loses the equilibrium to infinity rather than by annihilation, but the θ_crit = L − k_e boundary is the same idea.)

> **Provenance health warning on every Model-B number (round-9, item 4).** The verifier found that `params.py` — the paper's parameter-provenance registry, which has a row for every constant with its source — imports only `ParamsA` and `ParamsC`. **The entire Model-B parameter set is unregistered:** neither its alt-route V_max (k_e = 400) nor its half-saturation (K_e = 300) has a provenance row, and both are bare defaults with no citation. So Model B is the *least*-grounded of the three structures — weaker than A (whose k_e has a row stating "no independent human measurement," at least acknowledged) and than C (whose constants are registered, if marked arbitrary/assumed/reverse-engineered). Every B column below, and the "B is better cross-sectionally" result, inherits this: it is a statement about a parameter set the paper never sourced. Read the B numbers as *what would hold if ParamsB were right*, not as grounded values.

**Where a steady state does exist, the headlines move in magnitude but not in sign, and B is *better* on both axes.** A and B side by side at 5% of normal — at ParamsB's (miscalibrated, unregistered) anchor and at the valid managed point (k_e = 100, managed L = 220):

| quantity | Model A | Model B (ParamsB) | Model B (valid: low ceiling + managed load) |
|---|---|---|---|
| enzyme share s at 5% | 0.16 | 0.37 | 0.63 (untreated balance) |
| cross-sectional sd(ln θ), monthly | 3.05 | 0.92 | **0.29** |
| longitudinal contrast, weekly | 0.133 | 0.088 | **0.078** |

Under B the saturating alt route cannot soak up the rising marker, so the enzyme carries a **larger** share and θ is better identified — **both cross-sectionally (0.29–0.92 vs 3.05, still > 0.20) and in the contrast (0.078–0.088 vs 0.133).** *(v11 reported the contrast **worsening** to 0.167; that was a bug — the Model-B contrast failed to share the load-bias b between windows. b is a patient constant and cancels, exactly as in R-D; sharing it, B is better, not worse. There is no cross-sectional/longitudinal asymmetry — the single share mechanism moves both the same way, R-H reconciliation and `t36`.)* *Why B's sd is finite at all: over the **dosed** record m stays low (≈[3, 52]), where both channels are still rising — the paper's C23 "affine over the observed range" condition holds for the data; saturation bites only at the untreated balance (m ≈ 1273 near the managed ceiling) and below θ_crit, not in the dosed record.* (`t33`, `t35`, independently verified.)

**Reconciling "B is better for θ" with the paper's Remark 15 ("if identification fails in A it fails *a fortiori* in B") and Figure 13 (`t36_reconcile_remark15.py`).** The paper's "ceiling hypothesis" (Remark `rem:ceiling`) *is* Model B, and it prices B's cost two ways: the condition number of the column-normalised sensitivity matrix rises ~450× at a 10% excursion, and (Prop `exchange`) B carries an **exact two-point degeneracy** in which (θ,K_m) and (k_e,K_e) swap. My "B is 3.3× better" is not a contradiction — it is a *different object*, and the reconciliation refines what Remark 3 called an unestablished relation:

- **The mechanism of the marginal gain is the share.** sd(ln θ) tracks (1−s)/s (the clearance the nuisances carry relative to the enzyme's); 0.16 → 0.37 gives (1−s)/s falling 5.14 → 1.69, a factor **3.04 — 92% of the observed 3.30 sd ratio** (the ~8% residual is the extra nuisance K_e and the shifted operating point).
- **The condition number does degrade — the paper's direction.** In my monthly-record FIM, cond(B) / cond(A) = **18×** (same sign as the paper's 450×; the magnitude differs because the paper's figure is the *dynamic* sensitivity matrix at a 10% excursion). With K_e pinned, B's conditioning is actually ≤ A's — **the degradation lives entirely in the K_e column**, i.e. the k_e↔K_e exchange direction.
- **The marginal θ improves anyway because it is prior-regularised.** The ill-conditioning sits in the *nuisance* subspace, which the k_e/K_e priors absorb, protecting θ's marginal variance. **Strip the priors and B is catastrophic, exactly as the paper says:** with no nuisance priors and K_e free, B's sd(ln θ) = **55** vs A's 5.9 — the exchange degeneracy in full.

So the two objects move in **opposite directions**: the *joint* conditioning degrades under B (confirming the paper's a-fortiori for identifiability and the exchange degeneracy), while the *marginal, prior-regularised* precision on the one parameter of interest *improves*. The paper's "A is the most favourable case" holds for conditioning; it does **not** hold for the marginal CRB on θ once the nuisances carry reasonable priors — a narrow but real correction to Remark 15, bearing on the paper's stated reason for working in A (changed-claims).

**The provenance of the valid Model-B point, and its range (`t38`).** The valid column is a *construction*, not a fit: the ceiling k_e = 100 is chosen *below* the untreated load (so θ = 0 is not maintainable — the severe phenotype), literature-consistent but not a measured value (the repo gives no k_e in absolute units); the managed load L = 220 is chosen to admit a steady state with a ~1%-of-normal margin (a ~20% reduction). Because both are choices, the honest output is a **range**, not the point 0.29. Swept over the plausible grid (ceiling k_e ∈ 50–200, managed L ∈ 180–240, wherever a steady state exists), the cross-sectional sd runs **0.21 → 0.49** (median 0.38) — every valid cell above 0.20 and far below Model A's 3.05. The magnitude is model- and operating-point-dependent; the **sign** (misses 0.20 cross-sectionally, before any k_e sensor) is robust across the whole grid. The ~10–20% load reduction the margin needs is the same order as the relevant amplitudes — management supplies 20–40% (literature: Thompson 1990 et al., *not* the parent paper's own text), and the smallest illness excursion the parent paper says strands a patient is +18% (`§3.5`) — so the patient lives in a ~10–40% band where routine management keeps them just inside and a single ~18% illness pushes them out; margin and trigger are comparable, which is why decompensation recurs.

**Does the "reaches ~0.19" result survive the sweep and the K_e assumption? It only hugs the gate, and how robustly depends on the unmeasured k_e prior (`t40`).** The 0.21–0.49 range above was at the *base* priors; the ~0.19 was one constructed point with the k_e prior tightened and K_e known. Re-running the ceiling × load sweep *under those tightened conditions*: at the tracer's k_e precision (0.054) with K_e known, the **whole ceiling sweep (k_e 50–200) sits at the gate — 0.18–0.20, marginal everywhere, never comfortably below**; at the looser 0.2 prior it slips off the gate at higher ceilings (low ceilings pass, k_e 150–200 reach 0.21–0.23, not beyond). Within a ceiling the sd is nearly flat in the managed load, so **it tracks the ceiling — the enzyme share — not the margin μ** (the sorted sd-vs-μ list shows no monotone trend), and the cells nearest the gate are the *low-ceiling, smallest-margin* (most fragile) ones. So it does not survive at a less fragile operating point *by making it less fragile* — it survives only by lowering an unmeasured ceiling and tightening an unmeasured prior; a gate-hugging point held up by assumptions, not a comfortable pass.

**K_e provenance and tolerance (`t40`).** The 0.19 also requires **K_e known**, and K_e is the half-saturation of the alternative route — a **bare default** in the repo (`ParamsB.Ke = 300`, `models.py`) with **no entry in the provenance registry** (`params.py`) and **no human measurement** (the same gap as k_e). Worse, the exchange degeneracy the paper prices is *exactly* (θ, K_m) ↔ (k_e, K_e): assuming K_e known **collapses that degeneracy by fiat**, pinning the very constant whose unknownness *is* the degeneracy. Reporting sd(ln θ) as a function of sd(ln K_e) (at k_e = 100, L = 220, k_e prior 0.054): it stays ≤ 0.20 only while **K_e is known to ~30%** (sd(ln K_e) ≲ 0.30 → 0.19–0.20; 0.50 → 0.22; free → ~9), and no data fixes K_e to any precision.

*Stated as a specification* — the same shape as the flux reading's ≤19% and κ's tolerance: **an independent determination of K_e to ~30% brings the Model-B inference route to the gate.** What would that determination be? K_e is the half-saturation of the *lumped* non-PCC clearance (renal conjugate excretion, odd-chain lipid incorporation, protein propionylation, the methylcitrate shunt, carnitine buffering — a composite, not one enzyme). To fix it one would have to titrate that aggregate clearance against metabolite level across the operating range in vivo, in humans — a measurement the round-6 literature shows **has never been attempted** (there is no human k_e, let alone a K_e) and arguably *cannot* be made as a single constant, since the composite's half-saturation is not a property of any one transporter or pathway. So this is not a gap a feasible study closes, as κ's calibration or the flux reading's recovery might be; it is a request for a shape parameter of a lumped channel that may not have one. It is the weakest-grounded of the three constants, and the route that needs it is correspondingly the least promising.

**K_e is ill-posed, not merely unmeasured — which closes the route rather than leaving it pending (round-10, research doc 3).** The round-10 kinetics review is explicit that "no single half-saturation constant exists for a lumped non-PCC clearance route, and the concept is not well defined as one constant." The composite aggregates processes with *qualitatively different* kinetics: CrAT carnitine conjugation (propionyl-CoA Km ~2–50 µM, in vitro), GLYAT glycine conjugation (**sigmoidal, not hyperbolic** — so no single Km, K₀.₅ ~10–200 µM acyl-CoA / ~6 mM glycine), OCTN2-governed renal handling (a **threshold** at 40–60 µmol/L, not a Km), and **first-order, non-saturable** glomerular filtration of the free acids; and the likely rate-limiter is ACSS3 activation of free propionate, for which **no human Km is published at all**. A first-order term lumped with saturable ones has no fixed half-saturation — its apparent "K_e" drifts with concentration. So the Model-B inference route does not merely await a measurement of K_e; **it rests on a constant that is not a single constant.** That closes the route on structural grounds, independent of whether anyone measures it.

**The mass-balance conflict in physical units (round-10 bounds).** R-H's "Model A's k_e·m clears 84% of the load at 5%" can now be put against human numbers: healthy propionate turnover is ~17.6 µmol/kg/h, and *measured* non-PCC urinary disposal in PA is 1.5–3 µmol/kg/h (vs 38–55 in MMA) — i.e. the real non-PCC routes clear only **~10–20% of normal turnover**, "low-capacity overflow pathways, not a substitute for PCC" (research doc 3). Model A's linear k_e·m, carrying 84% at 5% capacity, overstates non-PCC disposal by ~4–8×; the literature puts the true figure near the enzyme-share end, confirming (in physical units now, not just as a fraction) that the bounded, low-capacity clearance of Model B/C is the physiological one and Model A's is not.

### R-I. Model C, and the pattern as a structural law (`t41_model_c.py`)

The pattern — each cross-sectional route reaches target only by importing a constant unmeasured in humans — reads as a coincidence while it is a list (κ ex vivo, c breath, K_e Model-B inference). Model C decides whether it is a law, because **C is the only structural hypothesis in which the competing channel has its own routine observable: free carnitine.** At steady state the carnitine balance `R₀ = c·k·(u(m) + ρ₀(1+m/M_i))` carries **no θ** (`rhs_C`), so free carnitine reads the alternative channel *directly*, independent of the enzyme — exactly the thing that, if measured, could break the mass-balance route's dependence on the competing channel.

Built the two-observable (m, c) FIM under C at 5% of normal, managed load 500 (a ~17% reduction below the paper's baseline 600, since at θ = 150 the load must sit below C's fold at 599.5 for a steady state — the same managed-load footing as B):

| observables / knowledge | sd(ln θ) at 5% |
|---|---|
| A — marker only, linear clearance | 3.05 |
| B — marker + K_e, managed, K_e free | 0.29 |
| C — marker only, alt channel free | ~20 |
| **C — marker + free carnitine, alt channel FREE** (K_c, R₀, M_i, ρ₀, k) | **~16** |
| **C — marker + free carnitine, alt channel KNOWN** | **0.34** |

So the carnitine reading **does** carry real information — with the channel's constants known it pins θ to 0.34, comparable to B — **but only then.** With the channel free, θ reverts to ~16: observing the competing channel is not enough; you must also know its *shape*. And the constants that must be known are exactly the ones the repo flags as unmeasured — **K_c "arbitrary," R₀ reverse-engineered** (so that the alt route peaks at 500 < the 600 baseline, making θ = 0 unmaintainable), **M_i arbitrary, ρ₀ a one-sided bound** — with K_c (free → 0.44) and M_i (free → 0.49) the most binding. (R₀ and k are additionally degenerate in c, which depends only on R₀/k, so one of the two must be known even to use the carnitine reading.) Swept over the admissible ranges with the channel known, sd stays 0.25–0.70 — robust *when known*, which is the point: the dependence is on knowing the channel, not on the operating point.

**So the pattern is a law, not a coincidence (the four instances, and the argument).** Every route that reaches the cross-sectional target imports a constant unmeasured in humans — **κ** (ex vivo transduction, R-G), **c** (breath/flux transduction, R11), **K_e** (Model-B competing-channel shape, R-H/`t40`), and now **K_c/R₀/M_i/ρ₀** (Model-C competing-channel shape, even with the channel *observed*). The decision states the structural reason: an observable reads θ either through the mass balance (needing the competing channel's shape) or by direct transduction (needing a calibration), there is no third way, and Model C shows that even *observing* the competing channel does not remove the need for its shape. The accumulation marker cannot pin absolute capacity cross-sectionally without an externally supplied constant, and across A/B/C and the full sensor panel the repository and the human literature supply none.

**Round-10 correction and a fifth instance (`t43`).** Round 9's framing slightly over-credited the carnitine observable — it said "c identifies θ (0.34) when the channel is known." On inspection, **m alone with the channel known already gives 0.34**: the carnitine reading is *redundant* for θ there, and the 0.34 came from the channel-known *assumption*, not from c. In the channel-*unknown* regime (the one the disease is in), c helps only marginally — m-only ~20, m+c ~16 — and that marginal help **erodes under the realistic plasma compartment**: Model C's c is a single well-mixed pool, but the observable is *plasma* free carnitine, which is <1% of body carnitine (>95% is in muscle, on a ~230 h vs <1 h two-timescale split) and whose partition to the pool is unknown and *drifts with the acyl load* (Stanley 1993: rising acylcarnitine displaces free carnitine). Carrying that plasma→pool partition as a free, load-dependent constant pushes m+c back toward m-only (~16 → ~17.7, `t43`). So the one model in which the competing channel looked *observable* turns out not to deliver it — plasma carnitine reports the modelled pool only up to a **fifth** unmeasured constant (the partition), the same measured-≠-modelled-compartment problem the parent paper flags for the marker m. The pattern does not merely survive the best-case counterexample; the counterexample dissolves into a fifth instance of it.

**And a sixth, of a new kind (round 11, R-J) — later corrected to the pattern's boundary (round 15, R-N).** The first five are constants the *cross-sectional estimate* of θ must import. Round 11 showed the pattern reaches the *decision* too: the dosing-interval choice is determined by data in hand (the drug trough dominates θ_crit), but only conditional on **λ**, the increment's decay rate, which swings Δ_trough(14) by a factor of 42 across a 2–7-day half-life. Round 11 called λ "the sixth instance." **Round 15 corrects that: λ is measured in humans (3.6 d, Chandler 1985) and transfers across species because it is an intensive tissue-level rate — unlike κ's extensive cell-to-whole-body map, which does not transfer. Same cell type, opposite verdict: λ is not a sixth *unmeasured* instance but the pattern's boundary, the one decision-critical constant that is both measurable and transferable (R-N item 6).** What survives is that λ gates the decision and its measured value sits on the thresholds.

**The classification, checked against the repo's full ten-sensor panel (`t26`/PHASE10).** Every observable falls on one side of the dichotomy (or, for upstream-pool reads, needs *both* kinds); the nuisance sensors read an input, not θ, so they feed the through-balance route rather than being a route themselves:

| # | observable | reads θ via | unmeasured constant it needs |
|---|---|---|---|
| 2 | ex vivo PCC activity (A = κ·θ) | **direct transduction** | κ (cell→whole-body mapping) |
| 1 | ¹³C breath / enzymatic-flux (O = c·θ·g) | **direct transduction** | c (recovery) |
| 3 | plasma 3-HP (the marker m) | **through the balance** | competing-channel shape (k_e,K_e / K_c,R₀,M_i,ρ₀) |
| 4 | 2-methylcitrate | **both** — upstream (propionyl-CoA) pool | pool calibration + the pool's own balance |
| 5 | C3 / propionylcarnitine | **both** — pool × carnitine | carnitine transduction + the pool's balance |
| 6 | C3/C2 ratio | **both** — pool read | as #5 |
| 7 | propionylglycine | **both** — pool × glycine (GLYAT) | glycine transduction + the pool's balance |
| 8 | IV ¹³C₃ micro-dose (CL = k_e + θ/(K_m+m₀)) | nuisance sensor → **feeds the balance** | reads k_e; θ still needs the balance |
| 9 | whole-body turnover (R_a = L + b) | nuisance sensor (reads the **load**, an input) | not a θ-route; supplies s |
| 10 | timed urine (CL_renal = f_renal·k_e) | nuisance sensor → feeds the balance | reads k_e (bound) |

No observable reads θ free of a constant: the direct reads (1,2) need a calibration, the state read (3) and the upstream-pool reads (4–7) need the competing channel's shape (and the pool reads a transduction on top), and the "nuisance" sensors (8–10) do not read θ at all — they supply k_e or the load into the through-balance route, which still needs the competing channel's shape to finish. The dichotomy is exhaustive, as the argument requires.

**The affine premise is checkable, and it has now been checked on the repo's data (`t39`).** The whole degeneracy rests on the paper's condition `eq:affine`: over the observed range, sup|C(m) − (a+bm)| ≤ ε (ε the measurement error). Nobody had evaluated it on the repo's records. Taking ε as the per-draw flux noise the FIM already uses: over the **dosed record** (m ≈ [3, 52]) the gap is **0.3 ε — it PASSES**, so the affine degeneracy genuinely holds on the data the model is fit to. Over the full span up to the **resting balance** (m ≈ [3, 1273], crossing the saturation knee) the gap is **~7 ε — it FAILS**, so a record that wandered up there would carry curvature information. The dose train holds the operating point low, which is exactly why the factor-75 gap between the record (m ≈ 17) and the resting balance (m ≈ 1273) does not break the analysis: the affine premise is safe where the data live and fails only where they don't go. (The far saturated tail is itself flat, hence affine; the non-affinity is the knee the operating point never reaches.)

**The sharing contract is now in one place, with a test (`contrast.py`, `t37`).** The round-7/8 contrast bug (a patient constant left unshared between windows) occurred three times — k_e once, b twice. All FIM contrasts now route through a single `window_contrast` helper that shares *every* patient constant (k_e, K_m, b, K_e) unless given an explicit drift, and `check_contract` fails if any is left independent (verified against a negative control). Re-auditing **every** FIM contrast number in the report through the helper: **R-D is clean** — the marker contrasts (0.26/0.13/0.07 and the drift breakpoints) reproduce exactly, because `t24` shared b and k_e from the start. The only number the bug ever touched is the Model-B contrast, 0.167 → 0.088.

**So does Model B change the framing?** The precision numbers shift (s = 0.16 → 0.37; sd 3.05 → 0.29–0.92; contrast 0.13 → 0.08) but **every verdict keeps its sign** — cross-sectional still fails, monitoring still passes — and B is the *easier* model for the marginal θ, not the harder one. What changes is **scope**: the framework describes the *managed, stable inter-crisis* state. Under a bounded clearance a 5% patient has a steady state only once the load is managed below θ + k_e, with a thin margin μ that antibiotic/protein excursions can erase (a θ_crit excursion is a decompensation). The report now says this at the top (see the decision's scope note). The obvious next model-side step — redoing the identifiability under Model B/C *with a bounded, literature-calibrated ceiling* and the crisis boundary treated explicitly — is outside this identifiability brief; it is the parent paper's own `PHASE5`/Model-C question, and the literature says it is the one that matters.

---

The mass-balance conflict that motivates all this, in detail. In Model A the alternative clearance k_e·m is **linear and unbounded**, so as capacity falls and the marker m rises, k_e·m rises without limit. The disposal split inverts with capacity:

| capacity | enzymatic θ·g(m) | non-PCC k_e·m | non-PCC share |
|---|---|---|---|
| 100% (normal) | 272.7 | 50.0 | 15.5% |
| 10% | 94.1 | 228.6 | 70.8% |
| 5% (disease) | 52.6 | 270.1 | **83.7%** |
| 1% | 11.5 | 311.2 | 96.4% |

So at 5% capacity the model has the alternative route silently clearing **84%** of the load — which is how a 5% patient reaches a benign steady state in the model at all.

**This conflicts with the PA literature.** Measured urinary propionate-metabolite output in PA is **1.5–3 µmol/kg/h** (Thompson & Chalmers 1990), one-to-two orders of magnitude below estimated production (~1–2 mmol/kg/day) — i.e. real non-PCC disposal clears only a *few percent* of the load, which is precisely **why PA patients accumulate and decompensate**. A route that absorbs 84% is not the PA phenotype. *(The user's own caveat is correct and I keep it: urinary output is not all non-PCC disposal — odd-chain lipid incorporation, protein propionylation, the methylcitrate shunt and acyl-CoA/carnitine buffering are not excreted — so 1.5–3 is a floor. But research doc 4 is explicit that the remainder is **buffered or accumulated, not cleared at high capacity**, so the conflict is real at the level of disposal *capacity*, not just excretion.)*

- **No absolute-unit comparison is possible from the repo.** θ = 3000 and L ≈ 323 are in abstract units (calibrated to a mouse 3-HP baseline, not an absolute propionate flux); the repository supplies **no bridge to µmol/kg/h**, so only the *fraction* (84% non-PCC) is comparable to the literature — and that is what conflicts. An absolute check needs a unit anchor the repo does not contain (flagged, not fabricated).
- **k_e provenance confirmed.** The research confirms **no human measurement of k_e (non-PCC disposal / propionate half-life) exists** in the literature — supporting the parent paper's claim that no published source supplies k_e, and the brief's choice to carry it with a wide prior.

---

### R-J. The dosing-interval decision is a trough condition, and what it needs is λ (`t44_interval_decision.py`)

Every result above prices the *measurement* of θ. The therapy forces a different, coarser act: **booking a dosing interval.** The therapeutic increment Δ(t) enters capacity as (θ+Δ)·g(m), sags to a trough just before the next dose, and the safety constraint binds at that trough (paper §7: *"a dosing interval is chosen to protect the trough, the lowest enzymatic increment the patient sees before the next dose"*): **θ_base + Δ_trough(τ) > θ_crit.** On a steady-state train of decaying doses Δ_trough(τ) = D·e^{−λτ}/(1−e^{−λτ}).

**The formula is the repository's own dose train.** Δ_trough(14) = **400.0 µM/d** and Δ_peak(14) = **6400.0** reproduce `delta()` numerically to the decimal, and match the repo's reported `train_trough_frac` 0.134·θ=3000 = 402 and `train_peak_frac` 2.13·θ = 6390. Two dosing regimes bound the trough: **(a) fixed per-dose amplitude** D (the literal form of the user's equation), and **(b) same delivered rate**, D ∝ τ (the paper's §7 convention — the mean increment is D/(λτ), so holding it fixed forces D ∝ τ). Regime (b) is the repo's: it reproduces §7's *"87 % cut in the trough at a 28-day interval"* exactly (**88 %**) and its *"doubles at seven days"* (**2.5×**). (b) gives the narrower gaps, so it is used as the conservative primary.

**1. Tolerance per interval choice (same-rate regime).** The trough increments are 1000 / 686 / 400 / 347 µM/d at τ = 7 / 10 / 14 / 15, so adjacent choices are separated in *capacity* by:

| adjacent pair | capacity gap (same-rate) | req sd(ln θ), routine load (fold 1.0) | req sd(ln θ), illness load (fold 1.5) |
|---|---|---|---|
| 14 vs 15 (fine) | 53 µM/d | 0.64 | 0.08 |
| 14 vs 10 (coarse) | 286 µM/d | **9.7** | 0.48 |
| 14 vs 7 (coarse) | 600 µM/d | **20.3** | 1.01 |

The required relative tolerance is gap/θ*, θ* the residual capacity of the boundary patient (θ_crit − Δ_trough of the longer interval). The user's linearisation δτ = δθ/(λ·Δ_trough) is a tangent approximation good to tens of per cent: against the exact capacity gaps it *overstates* the tolerance at the fine and mid grids (deriv/req ≈ 1.31 at 14-vs-15, 1.11 at 14-vs-10) and *understates* at the coarse 14-vs-7 (≈ 0.92), because Δ_trough is convex in τ so the secant departs from the tangent in both directions as the step widens — the table reports the exact gaps, not the tangent. θ_crit rests on the unmeasured s_pk/m_tol convention: at the routine-load reading (fold 1.0, θ_crit = 429.5; the regime where θ+Δ_trough ≈ 14·θ* at the booked dose, matching the user's "drug dominates" intuition) the tolerances are enormous; at the conservative illness reading (fold 1.5, θ_crit = 994.3) the drug no longer dominates (θ+Δ_trough ≈ 1.7·θ*) and the fine grid tightens to ~0.08.

**2. Against what the record delivers (grid-resolution).** Achievable: record 3.05, managed 0.29, contrast 0.13–0.48. Under the routine-load reading the coarse grid gaps (req sd 9.7 and 20.3) exceed the record 3.05, so the *intervals are resolvable at the record precision* — but this is a statement about how far apart the grid points sit, not about a given patient (item 3 is the operative per-patient version, and it is weaker). Only the fine 14-vs-15 split needs the managed/contrast precision even in the spacing sense. Under the conservative fold-1.5 reading the coarse grid moves to needing the managed 0.29, and the fine split (0.08) passes nothing.

**3. The whole posterior, not just the fibre.** The sd's used are the brief's full-uncertainty numbers (3.05 already bundles the K_m fibre, the k_e prior, the load bias b and the assay), so simulating the interval choice at them *is* the proper version. Over a PA population (θ_base log-uniform, 0.17–5 % of normal; grid {7,10,14}):

| model / convention | sd = 3.05 (record) | sd = 0.29 (managed) |
|---|---|---|
| Model A, routine load (fold 1.0) | 61 % correct, 0 % unambiguous | 93 %, 72 % |
| Model A, illness load (fold 1.5) | 78 %, 0 % | 100 %, 100 % |
| **Model C, intrinsic margin (no m_tol)** | **100 %, 100 %** | **100 %, 100 %** |

The operative per-patient quantity is the **determined subgroup** — the fraction whose booking no measurement could change (true θ more than Z·sd from the nearest boundary), the last column above. At the record precision it is **0 %** under Model A: no patient is far enough from the single boundary θ_crit − Δ_trough(14), and correct-choice (61 %/78 %) is barely above the ~52 % base rate. It reaches a usable subgroup (72 % at fold 1.0, 100 % at fold 1.5) only at the managed 0.29. So the spacing tolerance (item 2: intervals far apart) and the per-patient fraction (item 3: 0 % determined at the record) answer different questions, and the per-patient one is operative: an individual booking at the record precision is *not* settled without measurement.

**Model C reads 100 % by removing the decision, not by resolving it.** Its intrinsic critical capacity θ_crit^int(s_base = 600) = **150.8 µM/d** — from the model's own saddle-node, no m_tol convention — is *below* Δ_trough(14) = 400, so θ_base + 400 > 150.8 for every θ_base ≥ 0: **no patient is near the boundary at all.** The 100 % is "the drug clears θ_crit for everyone," not "the boundary is well resolved" — a different and weaker claim. It carries a testable consequence: **Model C predicts the uniform fortnightly regimen is adequate for every PA patient**, checkable against the trial's own ~70 % C3 reduction; if any patient is inadequately controlled at the booked interval, the intrinsic margin is too low or λ is shorter than assumed. The decision survives the full uncertainty where the parameter never does — but under Model A only at the managed precision, and under Model C only because the drug removes the boundary.

**4. λ is the sixth unmeasured constant, and its tolerance is the specification.** *(⚠ Revised by R-N, round 15: a measured human value exists (3.6 d, Chandler 1985) and λ transfers as an intensive rate, so it is not *unmeasured* — but its measured value sits on the thresholds below, which sharpens the specification rather than softening it.)* `design/run.py` sets λ = ln2/3.5 /d and labels it **ASSUMED** ("nothing measures it"); `params.py` calls it "partly corroborated, and only early" — the sole support a translational PK model's k_deg = 0.0075/h for PCC protein, fitted in mouse/rat/monkey and allometrically scaled, inside a two-compartment model whose own protein half-life is ~7 d. λ enters the trough directly, and **Δ_trough(14) swings a factor of 42 across a 2–7-day half-life** (47 → 2000 µM/d), so "the drug term dominates θ_crit" is true at one end and false at the other. The round-11 conclusion is therefore conditional on an unmeasured human constant — **the same pattern as κ (ex vivo), c (breath), K_e (Model-B shape), the Model-C shape constants, and φ (plasma/pool partition) in R-I. λ is the sixth instance, and the first that gates the *decision* rather than the *estimate* of θ.** Stated as a specification, like the 19 % flux contour and the 30 % K_e:

| headline | holds while | half-life tolerance | λ tolerance |
|---|---|---|---|
| Model C: every patient clears τ = 14 (100 % determined) | Δ_trough(14) > θ_crit^int = 151 | half-life > **2.62 d** | λ ≲ +34 % of assumed (~25 % margin) |
| routine-load coarse grid determined at the record 3.05 | req sd(14 vs 10) ≥ 3.05 | half-life > **3.32 d** | λ within **~5 %** of assumed |

*(⚠ Round 19, R-S: 2.62 d and the ~25 % margin are at Model C's K_c = 100. At the worst corner of its admissible shape box the 14-day floor is 3.06 d, i.e. λ ≲ +14 % of the assumed 3.5 d, a ~13 % margin, and the 15-day floor is 3.21 d, ~8 %; `t53`.)*

The fold-1.0 coarse-grid "no measurement" claim is **fragile** to λ — θ* ≈ 30 sits on the boundary, so req sd(14 vs 10) collapses from 9.7 at the assumed 3.5 d to below 3.05 by 3.32 d — while the Model-C version has a ~25 % margin (~13 % at the worst corner of Model C's shape box, R-S). **What would measure it:** the mRNA-3927 programme did dose selection by population PK/PD modelling, so the trial's own PK data should contain λ (the decay of PCC activity / the marker response between doses). It is the one number the whole decision turns on, it has been measured, and it is simply absent from this repository.

The through-line: the paper spends itself on whether θ is identifiable, and θ is the wrong estimand for the decision the therapy forces. The interval choice is a sign condition on θ_base + Δ_trough − θ_crit in which the drug term dominates, so where the decision is determined it is determined by data already in hand (Model C, or the coarse grid under the routine-load reading), and where it is not (a patient near the boundary, or a shorter-than-assumed λ) the binding unknown is a pharmacokinetic rate, not θ. This is the decision-estimand lesson of the parent paper (the margin is better conditioned than the parameter) carried one step further: for the *interval* the margin is dominated by a term that needs no identification of θ at all — so **the report ends not on a tolerance for θ but on a specification for λ**, a measurement that exists outside this repository and would settle the decision the whole study was about. *(Round 12, R-K, then asks whether the on-therapy record supplies λ itself — it does.)*

---

### R-K. Two years on therapy: retrospective adjustment, and the record self-calibrates λ (`t45_on_therapy.py`)

A follow-on question the study had not asked: not *choosing* the interval at the start (R-J), but *adjusting* it after two years on the booked schedule. Every patient starts at the protocol's 14 d, which Model C says is safe for all (R-J), so nothing is risked during the accumulation — this is **not** the dual-control problem §7 flags, just reading an accumulated record. The therapy's schedule is **exogenous and dated** (fixed by protocol, not chosen by the estimator), so by Prop 2 (`tworoutes`) its increment is pure excitation: Δ(t) drives the marker over an ~8–10× range between dose and trough, and dated events (interval jitter, missed doses, a supply gap) widen it, sampling g(m) across its curvature — exactly what separates θ from k_e.

**1. What two years buys — with a premise correction.** The "3.05 record" is **not** pre-therapy; it is *already* the regular on-therapy monthly record. The genuinely flat (no-increment) record is **8.56**; the 14-day dose train sampled monthly is **3.05** — so the dose train itself does most of the lift (8.6 → 3.0), the increment's phase-drift against monthly sampling already being Prop-2 excitation. Dated events add 3.05 → **2.28**, and θ precision is **sampling-limited** (≈ 1/√n):

| sampling (draws / 2 yr) | sd(ln θ) | correct-choice | determined subgroup |
|---|---|---|---|
| monthly (25) | 2.28 | 64 % | 0 % |
| fortnightly (52) | 1.69 | 69 % | 0 % |
| weekly (104) | 1.27 | 73 % | 0 % |
| alternate-day (365) | 0.75 | 83 % | **27 %** |

So the honest answer to "can any patient's interval be *re-chosen* on evidence after two years" is **mostly no at realistic cadence**: sd(ln θ) never reaches the ~0.29 that would individualise the {7,10,14} grid under the routine-load convention, and the determined subgroup is 0 % through weekly sampling, reaching ~27 % only at alternate-day. Passive observation tightens θ a lot (8.6 → 2.3) but does not make the symmetric interval choice individually determined for most patients.

**2. λ is recoverable from the same record — the clear win.** The marker's response to a longer-than-usual gap depends on λ directly, so the dated schedule identifies it. sd(ln λ) is **0.055 at monthly** sampling, tightening to 0.023 at alternate-day — comfortably inside last round's ~25 % (Model-C) specification and essentially at the ~5 % (routine-load) one. **So the accumulated dosing record self-calibrates the one constant the decision turns on, and the trial's population-PK λ becomes a cross-check rather than a prerequisite.** The honest cost: freeing λ degrades θ (monthly 2.28 → 3.45), because θ and λ share part of the dose-response direction — but θ was not individualising the decision anyway, whereas λ *is* recovered to spec. The record gives up little that was usable to pin the thing that was binding.

**3. Lengthening is one-sided, and here the record delivers.** Shortening is safe and only costs money; lengthening risks a crisis, so it is a one-sided test: lengthen only when the lower 95 % bound on μ_long = θ_base + Δ_trough(τ_long) − θ_crit exceeds 0. Under Model C (θ_crit^int = 150.8), the 15-day trough (347 µM/d) clears θ_crit for everyone — **a free burden reduction, no evidence needed** — but the 21-day trough (143) does **not** (143 < 151), so lengthening to 21 d is the genuine evidence-gated decision. The lengthen-to-21 eligible fraction (lower-bound clears θ_crit) here reads 0 % monthly → 51 % alternate-day → 73 % at the managed 0.29. **⚠ These 21-day fractions are superseded by R-L and should not be used:** they applied a *representative* (θ = 150) and *λ-fixed* FIM sd to the whole cohort. R-L redoes it per-patient, with λ estimated from the same record, by NLS simulation — and finds the 21-day margin estimate does **not** tighten with record length (sd(μ̂) ≈ 140 µM/d flat from 6 months to 5 years), so 21-day lengthening is **not** supportable from the passive record at any horizon; it needs the active managed-load protocol. The 15-day step (safe for all, no evidence) and the λ self-calibration both stand.

The round-12 upshot refines R-J's closing: the decision needs λ, not θ — and λ is **in the record**, so the specification the report ended on is met by the accumulated data itself, with the trial PK as corroboration rather than a gate. What the record cannot do is individualise θ enough to re-choose the interval symmetrically; what it can do is self-calibrate λ and support the safe, one-sided lever. *(Round 13, R-L, asks WHEN each of these arrives — and finds the 21-day lever never does.)*

---

### R-L. WHEN: the timing of each step, and why 21-day lengthening never arrives (`t46_timing.py`)

Every round-11/12 result was stated at "two years" — the record length we inherited, not a derived quantity. Round 13 asks the timing, one-sided throughout (lengthening on evidence; a wrong call costs burden, never safety).

**A method point first, because it decides the answer.** The 21-day lengthening margin μ_long = θ_base + Δ_trough(21;λ) − θ_crit^int is a *one-sided* test on a quantity whose posterior is strongly non-Gaussian (sd(ln θ) runs 1.8–36 across the disease range). FIM propagation is unreliable here: the delta-method lower bound says "never eligible," a log-normal Monte Carlo says "eligible by ~2 years" — a ~400 µM/d disagreement. This is exactly the pitfall the parent paper flags in §7b (*"a delta-method calculation … would throw away the one favourable result … μ̂ is skewed and does not cross zero as often as its sd suggests"*). So, as the paper does for its sign-of-μ decision (fig D3), round 13 **simulates**: NLS fits of (θ, K_m, b, k_e, λ) to noisy on-therapy records, measuring the decision directly. λ itself is in the Gaussian regime (small, θ-independent sd), so its timing uses the FIM, which is trustworthy there.

**1. 21-day eligibility does not accrue — at any horizon.** The margin needs θ_base resolved to sd(ln θ) < ln(θ_base/7.9)/1.645 ≈ **0.81 at θ_base = 30, 1.28 at 65, 1.79 at 150** (top of the disease range) — i.e. near the *managed-load* precision (0.29), an order below what a passive record gives. The simulation confirms the margin estimate never gets there and, decisively, **does not tighten with record length**:

| θ_base (μ_true) | 6 mo | 24 mo | 60 mo |
|---|---|---|---|
| 150 (+142) | sd(μ̂) 141, P(μ̂>0) 86% | 162, 86% | 162, 86% |
| 65 (+57) | 121, 74% | 139, 76% | 135, 74% |
| 20 (+12) | 109, 70% | 126, 66% | 119, 62% |

sd(μ̂) stays ~105–165 µM/d — comparable to θ_crit^int (151) itself — from 6 months out to 5 years, and the lower tail stays below zero at every capacity, so the one-sided bound never cleanly clears 0. **This resolves the FIM ambiguity: the MC optimism was an artefact of ignoring NLS estimator bias/heavy tails; the pessimistic reading is right.** 21-day lengthening is **not supportable from the passive on-therapy record at any horizon** — there is no crossing time to report. (This supersedes R-K item 3, whose 51–73 % lengthen fractions used a representative, λ-fixed FIM sd and so were badly optimistic — a changed claim.)

**2. λ crosses early, and uniformly.** By contrast, λ — the binding constant, which governs the trough — is pinned long before anything else, and θ_base-independently:

| spec | monthly | fortnightly | weekly |
|---|---|---|---|
| λ to 25 % (Model-C) | 6 mo | 6 mo | 6 mo |
| λ to 5 % (routine-load) | 48 mo | 24 mo | 12 mo |

So the constant the whole decision turns on arrives first — within ~6 months to the Model-C spec at any cadence, within a year to the tight spec at weekly — while θ (per-patient, for eligibility) never arrives.

**3. The 15-day step needs no record and no waiting.** Δ_trough(15) = 347 > θ_crit^int = 151, so even θ_base = 0 clears the bar at 15 d (margin +197 µM/d). The 15-day step is safe at **initiation**, for every patient — the only recommendation in the report available on day one, with no measurement and no observation period.

**4. There is nothing to review for, on the 21-day step.** The review-schedule question presumes eligibility accrues; item 1 shows it does not, so no review interval unlocks 21 d — reviewing a passive record for it is waste at any time. The time-structured recommendations that *do* hold, in order of availability: **day 0** — step to 15 d for everyone; **~6–12 months** — λ confirmed from the dosing record, corroborating the trial PK; **21 d** — only via the *active managed-load protocol* (which reaches sd ≈ 0.29), a prospective protocol change, never an eligibility to wait for under passive observation.

The round-13 upshot closes the timing arc: of the three levers, the free one (15 d) is available immediately, the constant the decision needs (λ) is confirmed within months, and the aggressive one (21 d) is unreachable from passive data at any horizon and belongs to an active protocol if it is wanted at all. What the record can give — λ, and the safe uniform step — it gives early; what it cannot — per-patient θ to individualise the aggressive step — it never gives, however long it runs.

---

### R-M. The FIM-vs-simulation audit (round 14): the lesson, and a check that it broke nothing else

Round 13 reversed a positive claim because a FIM-propagated natural-scale quantity failed under simulation. That makes the failure mode a method lesson, and obliges an audit of everywhere else the report takes the same step.

**The lesson, stated generally.** sd(ln θ) is a *log-scale* quantity: the FIM/Gaussian approximation gives it reliably, and it tightens ~1/√n with the record. A *natural-scale* or *one-sided* decision built on the same posterior — a threshold crossing, the lower bound of a margin — is a different object. Because the posterior is log-normal (a heavy natural-scale right tail, a left tail reaching 0), the one-sided bound is **tail-dominated and does not tighten even as sd(ln θ) halves.** Round 13 is the demonstration: over 6→60 months the FIM sd(ln θ) fell ~2× (3.3→1.8 at the top of the disease range) while the natural-scale margin μ̂_long did **not** tighten at all (sd ≈ 140 µM/d, flat) and its one-sided lower bound never cleared 0. The two propagations disagree by construction: a delta-method lower bound on μ_long gave ≈ **−440** while a log-normal Monte Carlo gave ≈ **+19** — a **~450 µM/d** gap, the delta method discarding the favourable skew exactly as the parent paper warns (§7b). **The rule: read a log-scale sd from the FIM; read a natural-scale or one-sided decision by simulation.**

**What that puts at risk, and its status.** Every conclusion that propagates a FIM through the (strongly non-Gaussian at 5 %, sd(ln θ) ≈ 3) posterior to a natural-scale or one-sided quantity:

| conclusion | the fragile step | status |
|---|---|---|
| R-K/R-L: 21-day lengthening eligibility (t45 "51–73 %") | lognormal-MC / representative FIM sd → one-sided margin | **broke → corrected** by NLS simulation (R-L): does not accrue; superseded |
| R-J: interval determined-subgroup (61 % correct, 0 % unambiguous at the record) | lognormal-MC(FIM sd) → interval threshold crossing | **spot-checked (t47): survives** — NLS gives 54 % / 0 %, no better (pessimistic direction) |
| cross-sectional "cannot place θ at 5 %" | FIM sd → natural-scale placement | **spot-checked (t47) + already MC-backed**: NLS sd(ln θ̂) ≈ 8.3 (λ-free), only ~9 % of fits within ±20 % — a NEGATIVE, reinforced |
| R16: "recoverable" flux / share / a (6 %/28 %/8 %) | first-order FIM propagation of θ-proportional quantities | already flagged "not Gaussian there"; θ-proportional, so inherit θ's spread — **at risk, use as order-of-magnitude only** |
| sign-of-μ (treat / not-treat) at 5 % | — | the brief makes **no FIM claim**: it defers to the parent paper's *simulated* fig D3 |

Everything else the report reports is either a **log-scale sd** (FIM-reliable: every sd(ln θ), sd(ln λ), sd(ln k_e), the contrast, the interval *grid-resolution* tolerances of R-J item 1) or **already simulated** (the scenario-table Monte Carlo column, the flux crossing 18.4 %, the estimator comparisons, round-13's eligibility, round-14's spot-checks).

**Spot-check A — round-11's determined-subgroup survives.** Redoing the interval classification with NLS θ̂ on simulated on-therapy records (not θ̂ drawn from a log-normal with the FIM sd) gives **54 % correct-choice, 0 % determined subgroup** against the reported 61 % / 0 % — no more favourable, slightly less. The claim was already pessimistic ("per-patient not settled at the record"), and NLS bias/heavy tails only reinforce it. Survives.

**Spot-check B — the cross-sectional headline survives.** NLS θ̂ at θ = 5 % (on-therapy, 24 mo monthly, λ-free) has a realised sd(ln θ̂) ≈ **8.3** (heavier than the FIM's 3.05, inflated partly by freeing λ), median θ̂/θ ≈ 1.08 but a 5–95 % band of [0.00, 2.97] and only **8.7 %** of fits within ±20 % of truth. The record cannot place θ on the natural scale at 5 %. The headline is a NEGATIVE, so the tail-dominance that broke the round-13 *optimistic* claim only strengthens it (consistent with F6's "treat FIM absolutes at 5 % as optimistic").

**Correlation scan — are there other favourable cancellations linearisation discards?** The θ–λ correlation (+0.87) that made μ_long's skew favourable could have counterparts. The on-therapy 5-parameter posterior correlation matrix (θ = 150, 24 mo weekly) has exactly two strong entries: **θ–k_e = −0.99** and **θ–λ = +0.87** (with k_e–λ = −0.85 induced). Both already feed one-sided margins that are **handled by simulation, not linearisation**: θ–k_e is the paper's own fibre result for μ = θ − θ_crit (eq:mufibre — μ linear where θ is quadratic — simulated in fig D3), and θ–λ is μ_long (R-L, t46). The remaining correlations (K_m–b = +0.88 among the nuisances; everything with λ besides θ, k_e) do not feed a reported natural-scale θ-decision. **So there is no third strong correlation feeding an only-FIM-propagated natural-scale claim** — the two favourable cancellations the posterior hides are both already captured by simulation, and they run, as the parent paper predicts, in the optimistic direction (which is why only the optimistic claim, round 13, needed the correction).

The audit's conclusion: the method failure cost exactly one claim (round-13's, already corrected). The report's load-bearing conclusions are negatives or log-scale sds, which this failure mode cannot inflate; the two natural-scale decisions that carried weight are now simulation-checked and survive; and the correlation structure hides no further optimistic cancellation that the report propagated through linearisation. The study has been wrong in both directions and is now checked in both.

---

### R-N. λ re-grounded on the literature (round 15): measured, not assumed — and sitting on the thresholds (`t48_lambda_literature.py`)

Rounds 11–14 treated λ, the therapeutic increment's decay rate, as the binding **unmeasured** constant — "the sixth instance" of the pattern, ASSUMED from rodent PK. Two literature searches overturn that provenance twice over: a **measured human value exists**, and allometric scaling was never the right operation for it. The repo's entry (`design/run.py`: λ = ln2/3.5 /d, ASSUMED, "nothing measures it") is wrong on both counts — and the measured value lands right on the decision thresholds, so the correction runs *toward* them.

**The literature (full provenance in the replacement entry below).** Chandler & Ballard 1985 measured PCC turnover in **human** fibroblasts by [³H]biotin pulse-chase: k_deg = 0.19/d, **t½ ≈ 3.6 d** — a single-exponential rate, and the *only* published human PCC turnover number. Attarwala 2023 (the mRNA-3927 translational model) fitted k_deg = 0.0075/h in **mouse** (ln2/k_deg = 3.8 d intrinsic) inside a two-compartment protein layer with an effective terminal ~7 d; allometry there was applied only to the mRNA clearance, never to k_deg. Jiang 2020 corroborates at 5.3–6.7 d (mouse-liver exogenous human PCC). The repo's 3.5 d guess was, by luck, almost exactly the human value — but its *justification* was wrong and its *comfort* was illusory.

**1. Which λ does the trough see?** The human value is single-exponential (3.6 d) — no two-compartment subtlety. For the mouse model the trough sits at 14 d, past the fast phase, so it samples the slow decline; but mapping the two-compartment decay into the repo's single-exponential increment via the **steady-state trough** (the decision-relevant quantity) gives an effective **t½ ≈ 3–5 d under both reconstructions of the (unpublished) 2-comp structure — below the paper's 7 d terminal**, because 14-day spacing samples before terminal dominance. So both sources converge on a decision-relevant half-life of **~3–5 d** (4.8 d and 3.1 d for the two reconstructions), not 7 d. Reporting "at 7 d" uses the naive terminal mapping and *overstates* the trough; the honest mapping lands at the thresholds.

**2. The decision at the measured values, and where each headline flips.** Δ_trough(14/15/21 d) and the Model-C "safe for all" status (Δ_trough > θ_crit^int = 150.8):

| t½ | Δ_tr(14) | Δ_tr(21) | 14 d | 15 d | 21 d | |
|---|---|---|---|---|---|---|
| 3.5 d | 400 | 143 | ok | ok | **FAIL** | repo's old assumption |
| 3.6 d | 434 | 161 | ok | ok | ok | **human fibroblast (primary)** |
| 3.8 d | 506 | 200 | ok | ok | ok | ln2/k_deg (mouse intrinsic) |
| 7.0 d | 2000 | 1286 | ok | ok | ok | mouse effective terminal |

The half-life at which each headline **fails**, and its margin at the human 3.6 d: Model-C **14-day safe-for-all** fails below **2.62 d** (+38 % at 3.6); **15-day** below **2.75 d** (+31 %); **coarse-grid 14-vs-10 resolvable at the record** below **3.32 d** (+8 %); and the **21-day safe-for-all floor** below **3.55 d** (+1.5 %) *(the three Model C floors, 14-, 15- and 21-day, are at K_c = 100; across its unmeasured shape constants they rise to at most 3.06, 3.21 and 4.08 d, which cuts the 14- and 15-day margins at 3.6 d to +18 % and +12 % and puts the 21-day floor above 3.6 d at some settings, R-S)*. **The correction runs the wrong way, and across a knife-edge.** The repo's old 3.5 d sat *just below* the 21-day floor, which is why round 13 found 21-day needs θ the record cannot supply; the human-measured 3.6 d sits *just above* it, flipping 21-day to "free for all" — by ~1.5 %. And human **hepatic** PCC very likely turns over **faster** than fibroblast (liver is the fastest-turnover tissue, Rolfs 2021), so the true hepatic half-life may sit *at or below* 3.6 d — toward every threshold. *(⚠ Reversed by R-P, round 17: that used the whole-liver median, but PCC is a mitochondrial-**matrix** enzyme, and the liver-mitochondrial median (4.26 d) is slower than whole-liver and slightly slower than the 3.6 d fibroblast value — so hepatic PCC likely sits at or **above** 3.6 d, central ~4.3 d, and the correction runs **away** from the thresholds. The band still straddles.)* The honest statement is a **margin, not a pass**: at 3.6 d the 14-/15-day steps hold comfortably, the coarse grid barely holds, and the 21-day step is on its boundary.

**3. λ is now the dominant model uncertainty.** The source recommends a 3-fold base-case band (t½ 1.2–11.5 d) and 10-fold worst case. **Every interval headline flips inside the 3-fold band** (14-day at 2.62, 15-day at 2.75, coarse at 3.32, 21-day at 3.55 d). The within-species-and-method spread alone is ~3-fold (2.7–8.5 d), so this is *measurement spread*, not cross-species conservatism. By the source's own rule, a conclusion that flips inside the 3-fold band makes λ the dominant model uncertainty — which it now is for every interval conclusion in the report. Even the 15-day step (round 13's "free for all") is λ-conditional: it fails below 2.75 d, inside the band, though it holds at all three literature points.

**4. The provenance entry, replaced.** λ is no longer ASSUMED. The corrected entry:
- **primary:** Chandler & Ballard 1985, human fibroblast, 0.19/d (~3.6 d) — the only human PCC turnover number.
- **corroboration:** Jiang 2020 mouse-liver 5.3–6.7 d; Attarwala 2023 mouse fit ~7 d effective (3.8 d intrinsic).
- **limitation:** no human *hepatic* PCC half-life has ever been published; the human value is skin fibroblast. *(⚠ round 17, R-P: hepatic is likely **equal or slightly slower**, not faster — PCC is mitochondrial-matrix and the liver-mito median (4.26 d) exceeds the fibroblast value.)*
- **species-transfer argument:** protein half-lives track tissue proteolytic machinery and sequence, not body mass (Swovick 2018: mouse–rat k_deg ρ ≈ 0.80; the strongest phenotypic correlate of degradation rate is maximal lifespan, not mass). So a tissue-level k_deg transfers *without* allometric scaling — allometry (BW^0.75) applies to systemic clearances (CL, hepatic blood flow, GFR), not to intensive degradation constants.
- **residual risk, named:** biotin status, not body size. Apo-carboxylase degrades far faster than holo (biotin deficiency raised PC degradation 5.1-fold in rat hepatocytes); rodent chow is biotin-replete, human intake varies. Under biotin-replete conditions essentially all PCC is holo, so the risk is small, but it is the mechanism to interrogate if a conclusion turns out λ-sensitive — which, per item 3, every interval conclusion now is.

**5. The record's λ is a test, not a self-calibration.** Round 12 pins λ to sd(ln λ) ≈ 5–9 % within the first year from the on-therapy record. With the human (3.6 d) and mouse-terminal (7 d) values ln-separated by 0.665 and the record's precision ~0.04–0.09, the two hypotheses are separated by **d′ ≈ 7 σ at 6 months, rising to ~15–18 σ by 1–2 years** — the record can decisively tell which regime a real patient is in. So round 12's self-calibration is better stated as a **prediction the trial's own marker record could settle**: read λ from the dose-timing response and it says whether human hepatic PCC turns over near 3.6 d (headlines on the edge) or near 7 d (headlines comfortable). The 3.6-vs-3.8 d distinction is *not* resolvable (< 1.5 σ) and need not be — both sit on the same side of the knife-edge. *(⚠ Round 19: they do at K_c = 100, but not across Model C's shape constants: one constant at a time, at K_c ≤ 50 or M_i = 150 the crossing is 3.64–3.74 d, so 3.6 d fails while 3.8 d clears (R-S). And (R-R) the decision turns on which side of 3.5452 d the true value is, and near 3.6 d that is exactly what one patient's record cannot resolve, at any length in `t46`. The "< 1.5 σ" also needs a qualifier: with `t52`'s sds at a 3.6 d truth, 3.6 vs 3.8 d is 0.6–1.2 σ within the first year and reaches about 2 σ only after five years of weekly draws.)* The decision does not need λ assumed; it needs λ **read**, and the record reads it. This is a stronger use of the self-calibration result than self-calibration: it turns an assumption into a testable quantity.

**6. λ is where the unmeasured-constant pattern stops — and why.** This search settles a question the report never asked: does λ have κ's problem? **It does not, and the reason is clean.** k_deg transfers across species because protein degradation is an **intensive** tissue-level rate (set by proteolytic machinery and sequence, conserved between related mammals). κ fails because the fibroblast-to-whole-body capacity map is **extensive** and tissue-weighted — a per-organism aggregate, not a rate. Same cell type (fibroblast), same literature tradition, **opposite verdicts**: κ cannot be carried from cells to whole body; k_deg can be carried from one species to another. So λ is **not** a sixth instance of the unmeasured-constant pattern (R-I) — it is the pattern's **boundary**: the one decision-critical constant that is both measurable and transferable. Rounds 11–14 mislabelled it "the sixth unmeasured constant"; round 15 corrects that to "measured, transferable, and decision-dominating." The remaining problem is not that λ is unknowable but that its known value sits on the thresholds, and its *hepatic* value — the decision-relevant one — is the single measurement still worth making.

The round-15 upshot: the study's binding constant is knowable after all (measured in humans, transferable as an intensive rate, and readable from the trial record), which is better than rounds 11–14 feared — but the measured value lands on the decision thresholds and the relevant hepatic value is likely faster still, which is worse than the repo's lucky 3.5-d guess implied. *(⚠ The "hepatic likely faster" clause is reversed by R-P, round 17: the mitochondrial-matrix comparator puts hepatic at or above the fibroblast value, central ~4.3 d; the measured value is near but on the comfortable side of the thresholds, not below them.)* Both corrections are logged; λ is promoted from "assumed" to "the report's dominant, measurable, and testable uncertainty."

---

### R-O. Can human hepatic PCC turnover be bridged? No — the trial record is the only route (`t49_hepatic_bridge.py`)

R-N left one question open: the decision-relevant value is human **hepatic** PCC turnover, and the only human datum is **fibroblast** (3.6 d). Can we bridge fibroblast → liver using the mouse tissue ratio, and does the bridged estimate fall below the 3.5452 d crossing? The answer is a clean negative, which is itself informative.

**There is no clean paired dataset.** The fibroblast numbers are cultured ([³H]biotin, Chandler; dynamic-SILAC, Swovick); the liver numbers are in-vivo (SILAC/heavy-water/¹⁵N: Rolfs 3.6 d, Price 3.5 d, Kim 4.26 d) or exogenous-mRNA (Jiang 5.3–6.7 d, Attarwala 3.8 d intrinsic). **No single study/method measures both tissues, and no mouse fibroblast PCC exists** — so a PCC-specific tissue ratio is unavailable. The only ratio is a cross-study ratio of tissue *medians*, and it is confounded: the cultured-fibroblast median (2.5 d) is inflated by cell division diluting protein, which the in-vivo liver median is not. So the three defensible readings **disagree on the direction** of the ratio:

| reading of liver/fibroblast ratio r | bridged median | P(hepatic < 3.5452 d) |
|---|---|---|
| A — cross-study medians at face value (r = 1.44, liver slower) | 5.2 d | **24 %** |
| B — like-for-like degradation measures (r ≈ 1.0) | 3.6 d | **49 %** |
| C — "liver is fastest in vivo" (r ≈ 0.8, liver faster) | 2.9 d | **65 %** |

Carrying the uncertainty properly (a ratio of medians applied to a single protein — anchor sd(ln) 0.20, ratio 0.35, per-protein 0.35; combined ~0.53, a ~1.7× one-σ band), **P(hepatic below the crossing) spans ~24–65 %, centred near 50 %.** The bridge cannot tell which side of the crossing human hepatic PCC falls on; it is consistent with either. The one suggestive anchor is a coincidence of like-for-like degradation measures — human fibroblast PCC (3.6 d) ≈ mouse liver median (3.6 d) ≈ mouse liver PCC intrinsic (3.8 d) — which argues the tissue effect is small and lands the estimate ~3.6 d, almost exactly *on* the crossing, where the probability is ~50/50 by construction.

**So the bridge is uninformative about the side of the crossing — and that is the answer: the trial record is the only route.** The decision-relevant quantity (human hepatic PCC turnover) cannot be carried from fibroblast data; it must be **read**. Round 12 shows the on-therapy record reads λ at 7–18 σ — it separates 3.6 from 7 d decisively — so the measurement already exists in principle wherever the trial's accumulated marker data lives. *(⚠ Corrected in round 19, R-R: this sentence used to add "and the two sides of a 3.5452 d crossing a fortiori". That runs the wrong way. The crossing is 43× closer to 3.6 d than 7 d is (ln gaps 0.015 vs 0.665), so per patient the record places a value near 3.6 d on a side of the crossing at only 0.2–0.6 σ, even after five years of weekly draws. It settles the side only when the true value is clear of the crossing: near R-P's 4.3 d central estimate, with an 80 % chance within a year of monthly draws, against the line at K_c = 100. Across Model C's shape constants the line can sit at 4.08 d, and there one record places a 4.3 d truth above it at 2σ with at most about a 47 % chance even after five years of weekly draws (R-S). `t52`.)*

**The two data requests the study ends on.** Everything the report can establish from the repository and the literature is now established; the two gaps that remain are not analyses but *data*:
1. **NCT04159103's on-therapy marker record** — the dated dose/response series that round 12 shows reads λ at 7–18 σ. This settles whether human hepatic PCC turns over near 3.6 d (the decision on its knife-edge) or near 7 d (comfortable), and does so per-patient. *(It does not settle, per patient, which side of the 3.5452 d crossing a value near 3.6 d is on; see R-R.)*
2. **The Gao et al. (ACoP 2025) population-model parameter table** — its indirect Kin-inhibition structure re-estimated from human data contains a human protein-decay term that is not public; it would supply λ directly, as a cross-check on (1).

With either, λ moves from "measured in the wrong tissue, sitting on the threshold" to "known where it matters." Without them, the honest statement stands: the dosing decision is on a λ knife-edge that the available data cannot resolve, and that the trial's own record resolves only if the true value is not on the knife-edge (R-R).

---

### R-P. The hepatic correction runs the other way, and the human λ is not public (round 17, `t50_hepatic_mito.py`)

Two final searches close the λ arc. The first **reverses a directional claim** carried since round 15; the second **records why the on-therapy record is the only accessible route**.

**1. The compartment correction: hepatic PCC is likely equal-or-slower than fibroblast, not faster.** Rounds 15–16 (R-N, R-O) carried the steer that human hepatic PCC turns over *faster* than fibroblast, because liver is the fastest-turnover *tissue* (Rolfs). That used the **whole-liver** median (2.4–3.6 d). But PCC is a mitochondrial-**matrix** enzyme, and the apt comparator is the **liver-mitochondrial** median — **4.26 d** (Kim 2012) — which is *slower* than whole-liver and slightly slower than the 3.6 d fibroblast value. The fibroblast proteome median is ~1.7 d, so Chandler's 3.6 d PCC already sits in the *slow tail* of its own distribution. Four anchors converge — fibroblast PCC 3.6 d, whole-liver 2.4–3.6 d, liver-mito 4.26 d, human bulk-liver FSR ~4 d — on a compartment-anchored estimate of **~3–6 d, central ~4.3 d, above the 3.5452 d crossing.** Re-running the decision (same wide band as R-O, sd(ln) ≈ 0.53, only the centre moved from 3.6 to 4.3 d):

| estimate | P(hepatic < 3.5452 d) |
|---|---|
| central 4.3 d, wide band (sd 0.53) | **36 %** |
| central 4.3 d, tight band (sd 0.35) | **29 %** |
| across the 3–6 d endpoints | 62 % (3.0 d) → 36 % (4.3 d) → 16 % (6.0 d) |

So **P(hepatic below the crossing) is ~30 %, not ~50 %** — the 21-day safe-for-all step is now **~65–70 % likely to hold** *(64 % on the wide band and 71 % on the tight, both at K_c = 100; at the least favourable Model C shape constants 54–56 % on the wide band or 56–59 % on the tight, R-S)*, the opposite of round 15's "toward every threshold." The correction runs *away* from the thresholds. **But the band still straddles:** a hepatic half-life below the ~3.55 d crossing (the lower part of the plausible 3–6 d band) still puts the 21-day step in the needs-θ regime, and P(hepatic < crossing) is ~30 % at the central 4.3 d, rising toward ~50 % only if the true value is near the 3.6 d fibroblast anchor. The direction reverses; the side is still not settled.

**2. How much the bridge can carry — it narrows without settling.** This is the limit of the fibroblast→liver inference, and it is a real limit. No single study measures both tissues; the closest same-species/same-method pair is mouse (Swovick dermal fibroblast vs Rolfs/Hasper whole liver), and for human only a cross-study pair (Mathieson hepatocytes vs human dermal fibroblast). **No numeric PCCA/PCCB half-life appears in any retrievable supplementary table**, so the ratio cannot be anchored on the protein itself. Within-tissue turnover spans 2–4 orders of magnitude and cross-tissue correlations are only R² ≈ 0.3–0.65, so a *median*-based ratio is only weakly transferable to one protein. What the compartment context *does* buy is the 3–6 d band and the reversed direction; what it cannot buy is the side of the crossing. The bridge narrows the question; it does not settle it.

**3. The human λ is not public — the route, recorded.** The one place a human-fitted PCC turnover constant exists is the Gao et al. population model (ACoP 2025 poster S-104, NCT04159103): an indirect Kin-inhibition structure re-estimated from human data. Its public abstract reports only CL = 0.266 L/h and V = 13.9 L; **the turnover term is not public** — the poster is gated on Eventscribe with no mirror, there is no full publication as of October 2026, the Nature interim paper reports only non-compartmental PK, and a FOIA request would fail (unapproved product, Exemption 4). Corresponding authors: **Wei Gao and Linh Van, Moderna.** This makes round 12's result decisive in practice: the **on-therapy marker record is the only accessible route to a human λ** (it reads λ at 7–18 σ against the 3.6-vs-7 d question), which strengthens the case for requesting NCT04159103's data over waiting on the model table. *(⚠ Round 19, R-R: this used to say the record reads λ "settling the side of the crossing". Per patient it does so only when the true value is clear of the crossing; near 3.6 d it cannot, at any record length in `t46`.)*

**The study's two remaining gaps are two emails,** not two analyses:
1. **NCT04159103's on-therapy marker record** (the data owner, via the trial sponsor) — reads λ per-patient at 7–18 σ for 3.6 vs 7 d, and settles the knife-edge directly if the true value is clear of it (near 4.3 d: an 80 % chance within a year of monthly draws, or 6 months of weekly ones, against the 3.5452 d line at K_c = 100). Across Model C's shape constants the line can sit at 4.08 d, and then one record places a 4.3 d truth above it at 2σ with at most about a 47 % chance even after five years of weekly draws, 11 % after a year of monthly ones (R-S). If the true value is near 3.6 d, no single patient's record can place it; a pooled estimate could, with 23–77 patients at 50 % power or 46–155 at 80 % (the larger counts at monthly draws), but only if λ does not vary between patients, which is unmeasured (R-R).
2. **The Gao et al. parameter table** (Wei Gao / Linh Van, Moderna) — supplies the human-fitted λ as a cross-check.

Everything the repository and the literature can establish is established and independently checked. The dosing decision reduces to one measurement — human hepatic PCC turnover, now best estimated at ~4.3 d (above the crossing, 21-day step ~65–70 % likely free at K_c = 100, about 55–60 % at the least favourable shape constants (R-S), band straddling) — that exists in two places neither of which is public, and that the trial's own accumulated record can read. That is where the λ arc stops. The round that follows re-reads the whole identifiability result through a geometry a reviewer proposed; it changes no number.

### R-Q. The line-pencil geometry (the PI's formulation): an independent check and a better figure (round 18, `t51_pencil_geometry.py`, `r18_figs.py`)

The PI reformulated the steady-state reading as geometry. From the Model-A balance (θ+Δᵢ)g(mᵢ)+k_e·mᵢ=Lᵢ, moving the **known** therapy increment to the right makes each draw a **line in the (θ, k_e) plane**,

  g(mᵢ)·θ + mᵢ·k_e = Lᵢ − Δᵢ g(mᵢ),  g(m)=m/(K_m+m).

Noise-free and at the true K_m, every line passes through the true (θ, k_e): the pencil **concurs** (checked to 10⁻¹³). With load-record and assay error the lines form a near-concurrent pencil, and the size of the near-concurrence region is the identifiability of (θ, k_e). This does not change any number in the report; it checks several from an independent direction and produces the figure the condition numbers only described.

**The mapping, stated explicitly (three confirmations, one refinement).** The PI's objects are the report's objects:

- **The −0.99 θ–k_e correlation *is* the near-parallelism of the pencil. Confirmed.** Every line has slope dk_e/dθ = −g(mᵢ)/mᵢ = −1/(K_m+mᵢ); over the dosed record's narrow band (m ≈ [4, 38]) the slopes span only [−0.0096, −0.0074], so the lines are nearly parallel. The pencil's own (weighted) concurrence region has corr(ln θ, ln k_e) = **−0.998** (the report's −0.99), an ellipse of aspect ratio **29 : 1**.
- **Proposition 1's affine degeneracy *is* the narrow-band slope collapse. Confirmed.** As the band narrows, −1/(K_m+mᵢ) becomes constant, the slopes collapse to one direction, and the pencil degenerates from "concurrent at a point" to "a sheaf of parallel lines" whose intersection smears along the fibre — exactly the affine-clearance degeneracy, now visible.
- **Rounds 13–14's NLS simulation *is* the Ωₚ sampling.** Confirmed — drawing parameter vectors, simulating, re-estimating, and scattering the estimates is what those rounds did.
- **The inverse-FIM ellipse *is* the concurrence region — with one refinement.** It is the concurrence region **at fixed K_m** (the conditional ellipse). Reading the same lines at a different K_m slides the concurrence point along a curve — **Remark 2's fibre** — so the full marginal uncertainty (the headline 3.05) is the fixed-K_m ellipse *plus* the fibre. And the pencil lives in **linear** (θ, k_e) while the FIM ellipse is Gaussian in **log** coordinates; the gap between them is precisely where round 14's heavy tail lives (below).

**1. The pencil, drawn (the figure the result deserves).** The honest content is that the lines do *not* meet — so the figure shows the lines, not a point, and why the region came out long.

![round 18 pencil: each steady-state draw is a bearing to the unknown (θ, k_e)](fig/r18_pencil.png)

Read it as **triangulation**: each draw is a bearing to an unknown position. Bearings taken from widely separated points cross sharply and give a tight fix; bearings all taken from nearly the same direction are near-parallel, so the position smears along the line of sight. *Panel A* is the actual dosed record — twenty-five bearings all from nearly one direction, because over the dosed band (m ≈ [4, 38] against K_m = 100) every slope is forced to ≈ −1/100 (the identity below); the 2σ region is a long sliver along the fibre, sd(ln θ) = 2.7, the identifiability failure as a picture. *Panel B* is a **second bearing from a different angle** — a wide excursion (the steady-state reading of the Task-7 drug-forcing design, loads swept 0.45–3.2×): the marker ranges widely, the slopes fan apart (the lines are coloured by m_i, so the high-marker draws that carry the different slopes are visible), and the pencil concurs — sd(ln θ) falls to **0.40**, a 7-fold tightening. The analogy carries the safety point without extra words: *the place you would need to stand to get the good angle is the place the patient is not allowed to go* (the excursion is past the crisis threshold, R15). Axes are shared and the positive quadrant only: the sliver continues off-panel into k_e < 0, which is not a parameter set — at the paper's σ_m/Δ = 20 %, 96 % of unconstrained fits land there (Remark `rem:eiv`; at this figure's 5 % assay, 26 %), and that inadmissibility is checkable from the fit alone.

The slope identity is the mechanism behind three of the report's results at once: every line has slope −1/(K_m + m_i), so the narrow band pins every slope near −1/100 — which is the −0.99 θ–k_e correlation, Proposition 1's affine degeneracy, and the fact that **band width is the only lever**.

A second figure makes the third parameter explicit. The pencil is really a **family indexed by K_m**: reading the same noise-free lines at K_m ∈ {50, 100, 200} places the concurrence point at three points along a curve — Remark `rem:fibre`'s fibre (the "Remark 2" of the proposal) — with the truth recovered only at the true K_m = 100.

![round 18 fibre: reading the pencil at different K_m slides the concurrence point along the fibre](fig/r18_fibre.png)

This is the refinement to the mapping above: the inverse-FIM ellipse is the concurrence region *at fixed K_m*; letting K_m vary adds this fibre, and the headline sd(ln θ) = 3.05 is the two together.

**2. Pairwise intersections — a model-free uncertainty that confirms the heavy tail from a third route.** The 25 lines give C(25,2) = 300 pairwise intersections. Their scatter estimates the uncertainty with **no noise model, no linearisation, no Gaussianity** — so, given round 14's lesson that FIM propagation to natural-scale quantities is unreliable here, it is a genuinely independent check. The PI predicted a heavy tail (near-parallel pairs meet near infinity); the record **confirms and sharpens** it:

- On the narrow dosed record the cloud is **Cauchy-like**: θ-IQR ≈ 3 800, full range ≈ 8×10⁵, and only **5 %** of intersections land within a factor 3 of the truth. The raw mean is meaningless. The tail is so heavy that even the **median is not a usable point estimate**: over 400 records the median-of-medians is 260 (a +0.55 bias in ln θ) and **35 % of records give a *negative* cloud median** — the robust centre is often nonsensical. The crossing-angle (det²) weighting the PI suggests pulls the centre back toward the MLE (θ ≈ 290 vs the weighted fit's 242), but the natural-scale spread stays enormous. **A near-parallel pencil has no natural-scale point estimate, only a log-scale width** — the same pathology round 13 found by simulation and round 14 made a method lesson, now from a third, assumption-free route.
- On the wide excursion the same construction **concurs**: median θ ≈ 117, 41 % of pairs within a factor 3, IQR ≈ 460. Width buys a usable estimate.

The three routes to the narrow-record uncertainty — all prior-free pencil geometry at fixed K_m — agree where they can and differ where they must: the FIM ellipse gives sd(θ) ≈ 444 and the weighted-NLS posterior an IQR ≈ 510 (a Gaussian of that sd has IQR ≈ 600) — **they agree in the core** — while the symmetric FIM ellipse cannot show the tail that both sampling routes do (intersection IQR ≈ 3 800). All three say θ is not placeable on the natural scale at the record precision. (The report's k_e prior tightens the fixed-K_m FIM to sd(θ) ≈ 332; the headline sd(ln θ) = 3.05 adds the K_m fibre on top — a looser, different object.)

**3. Concurrency as a model-adequacy test — and an honest account of its power.** If the lines fail to concur by more than assay noise explains, the model is misspecified over that band. The statistic is χ² = Σ(dᵢ/σ_{d,i})², the residual orthogonal distances of the best-fit point to each line in units of the distance the assay/load noise predicts, on N−2 dof. It is **calibrated** (rejects 4 % at α = 5 % over 300 records) and the **repo's simulated records PASS** (seed 0: χ² = 15.4/23, p = 0.88), as they must — they are generated from the one-channel model. Its power is asymmetric, and the asymmetry is the point:

- **A wrong / biased recorded load (an ordinate error) is detected**, with power growing in the bias: 8 / 15 / 22 / 58 % rejection for a 5 / 10 / 15 / 25 % load mis-record on half the draws.
- **A second saturable channel (a structural clearance error) is *not* detected** — even a strong one (60 % of the linear channel) on a wide band gives ≈ 0 % rejection. A second saturable term h₂(m) is nearly **collinear with {g(m), m}** over any realistic band, so the two-parameter line absorbs it: the misspecification lands **on the fibre** and concurrence is preserved.

So the test complements `t39` (which asks whether clearance is locally *linear*; this asks whether the two-channel point-model fits the pencil *at all*), and it guards against **data errors** — a mis-recorded load — not against the structural degeneracy the rest of the report is about, which no steady-state concurrency test can see.

**4. The weighting (Remark 6) — and it reproduces the FIM.** Remark 6 says the noise is in the **abscissa** (the assayed marker m), not the ordinate (the recorded load L), so ordinary least squares on the vertical residual is the wrong variance. The per-line orthogonal-distance sd is heteroscedastic (σ_d ∈ [0.6, 4.4] over the record — it grows with m and shrinks as the line flattens), so the unweighted Euclidean (P = I) fit mis-weights by the crossing angle and gives a different, worse estimate: (θ, k_e) = (725, 0.04) against the weighted (242, 3.99). The PI's weighted x′Px fit, with P the propagated noise metric, **is the report's MLE**: its normal matrix equals the **prior-free** fixed-K_m FIM of the data (corr −0.998 vs −0.998; sd(θ) 410 vs 444, an ~8 % realised-vs-expected sampling gap), and adding the report's k_e prior tightens it to sd(θ) ≈ 332 — closing the loop, the geometry and the information matrix are one object, and the det²-weighted intersection centre of item 2 lands near the same place (θ ≈ 290 vs 242). Remark 6 is that weighting made operational.

**What the geometry adds.** It checks the numbers; it does not change them. It gives (a) the pencil figure — *why* the fibre exists, not just that it does; (b) a model-free uncertainty whose heavy tail confirms rounds 13–14 from a third route and shows the pencil has no natural-scale centre; (c) a calibrated concurrency adequacy test, with power against a biased load but blind to a second channel (t39 seen from the pencil); and (d) Remark 6's weighting shown to reproduce the FIM/MLE. The one place it corrects an expectation: the robust intersection centre does **not** recover the truth on the degenerate record — the near-parallelism defeats even the median — which is itself the round-13/14 lesson, independently rediscovered.

### R-R. Can the trial record settle the knife-edge? Only if the answer is not on it (round 19, `t52_knife_resolution.py`)

Round 19 wrote a plain-language account of all eighteen rounds (`story_so_far.html`). An independent completeness critic, run on that account, found a claim in R-O and R-P that the numbers do not support, and an adversarial check confirmed it. Those rounds said the on-therapy record "reads λ at 7–18 σ" and therefore settles the two sides of the 3.5452 d crossing *a fortiori*. The 7–18 σ answers a different question: is the half-life near 3.6 d or near 7 d? Those two are ln(7/3.6) = 0.665 apart. The measured 3.6 d is ln(3.6/3.5452) = 0.015 from the crossing, **43× closer**. So placing a value near 3.6 d on one side of the crossing is much harder, not easier.

`t52` takes the per-patient FIM sd(ln λ) that `t46` computes (by record length and cadence) and re-evaluates it at each true half-life, because `t46` evaluates it at the repository's 3.5 d (re-implemented, it reproduces `t46` exactly). z is the *expected* separation, the log distance over the sd. One noisy record lands at z plus a unit Gaussian, so the chance that it actually clears 2σ on the correct side is Φ(z − 2), about 50 % when z = 2. Dead zones use the sd at the crossing; z uses the sd at each true half-life. That is why the 6-month monthly row can leave 4.3 d just outside its dead zone (sd 0.095 at the crossing) while z at 4.3 d is just under 2 (sd 0.100).

| record (one patient) | z at a true 3.6 d | z at a true 4.3 d | chance one record places 4.3 d at ≥ 2σ | 2σ dead zone (true half-life) |
|---|---|---|---|---|
| 6 months, monthly | 0.16 | 1.93 | 47 % | 2.93–4.28 d |
| 12 months, monthly | 0.23 | 2.97 | 83 % | 3.10–4.05 d |
| 6 months, weekly | 0.29 | 3.71 | 96 % | 3.19–3.94 d |
| 60 months, weekly (best in `t46`) | 0.58 | 7.22 | 100 % | 3.36–3.74 d |

**Per patient, the record never places a 3.6 d truth on a side of the crossing**, at any length or cadence in `t46` (best 0.58 σ, expected). It does settle the side when the truth is clear of the crossing. For a 4.3 d truth (R-P's central estimate), 6 months of monthly draws leave it just short of 2 σ (1.93, a 47 % chance for one record); an **80 % chance is reached within a year of monthly draws, or within 6 months of fortnightly or weekly ones** (at 12 months of monthly draws z = 2.97, an 83 % chance; below two years `t46`'s grid has only 6 and 12 months, and between them the monthly FIM sd does not fall steadily with record length, so no earlier month is quoted). Those figures are against the 3.5452 d line at K_c = 100. Across Model C's shape constants the line can sit at 4.08 d (R-S), only ln 0.052 below 4.3 d, and against that line one patient's record never reaches an expected 2 σ: at best z = 1.94 after five years of weekly draws (a 47 % chance), and z = 0.80 after a year of monthly ones (11 %). Against the largest sampled crossing, 3.97 d, a year of monthly draws gives 22 %. A pooled estimate across patients could place a 3.6 d truth at 2 σ with **23–77 patients at 50 % power, or 46–155 at 80 %** (the larger counts at monthly draws), but only if λ is the same in every patient, and nothing measured says it is.

What changes: the two emails are still the right next step, and the on-therapy record is still the best route to a human λ. But the claim that it "settles the knife-edge directly" is now conditional. **It settles the knife-edge only if the true liver half-life is not on it**, and because the knife-edge itself moves with Model C's shape constants (R-S), "not on it" has to hold wherever the line actually is: a 4.3 d truth is clear of the K_c = 100 line but close to the worst-corner one. If it is near 3.6 d, the honest result is the one the knife-edge figure already shows: under Model C at K_c = 100 the 21-day step clears by about 1.5 % in half-life (across its unmeasured shape constants it may not clear at all, R-S), and no single patient's record will confirm or refute that. R-O and R-P are flagged where the claim appeared, and the changed-claims table carries the reversal.


### R-S. The knife-edge has a second blade: Model C's unmeasured shape constants (round 19, `t53_crossing_shape.py`)

The 3.5452 d crossing is where the 21-day trough equals Model C's intrinsic critical capacity, θ_crit^int = 151 µM/d. That threshold is computed at K_c = 100, M_i = 300, ρ₀ = 0.5, and `models.py` calls K_c "ARBITRARY. There is no measurement behind it." The repository already sweeps these constants (`model_C/run.py`), re-deriving R₀ each time so the alternative route still peaks at 500. One at a time, θ_crit moves by a factor of 1.56 for K_c, 1.37 for M_i and 1.02 for ρ₀. Its joint check, `collapse_check`, draws 12 random sets and reports "a factor of about 1.4" (113.5–158.1). Twelve draws understate the joint spread: `t53` samples the same box 300 times (θ_crit^int 111–236 µM/d, median 136, a factor of 2.1) and evaluates its eight corners (110–262 µM/d, a factor of 2.4). No round carried any of this into the dosing crossing. `t53` does:

| | θ_crit^int (µM/d) | 14-day crossing | 15-day crossing | 21-day crossing | 3.6 d clears 21 d | 4.3 d clears 21 d |
|---|---|---|---|---|---|---|
| repo constants | 151 | 2.62 d | 2.75 d | 3.545 d | yes | yes |
| K_c 25 → 800 (one at a time) | 186.5 → 119.5 | | | 3.74 → 3.36 d | no at K_c ≤ 50 | yes |
| M_i 150 → 1200 (one at a time) | 171.4 → 125.5 | | | 3.66 → 3.40 d | no at M_i = 150 | yes |
| joint, 300 sampled sets | 111–236 | ≤ 2.96 d | ≤ 3.11 d | 3.30–3.97 d (median 3.46) | in 80 % of sets | in all |
| the 8 corners of the box | 110–262 | ≤ 3.06 d | ≤ 3.21 d | 3.29–4.08 d | at 6 of 8 (fails at K_c 25, M_i 100) | at all 8 |

(The 80 % is a share of sets sampled over the admissible ranges, not a probability; nothing measured puts a distribution on these constants.)

**What holds and what does not.** The 14- and 15-day conclusions do not depend on Model C's shape: even at the worst corner their crossings are 3.06 and 3.21 d, under the measured 3.6 d. They still need the liver half-life to be above those values, which sit at the bottom of R-P's plausible 3–6 d band. The 21-day conclusion at the measured 3.6 d does depend on the shape. The crossing moves by about −7 % to +15 % across the box (3.29–4.08 d against 3.545 d), wider than the 1.5 % by which 3.6 d clears it at K_c = 100. R-P's "about 65–70 % likely to hold" (64 % on its wide band, centre 4.3 d and sd(ln) 0.53; 71 % on its tight band, sd(ln) 0.35; both already at the repo constants) falls to 54–56 % on the wide band, or 56–59 % on the tight, at the least favourable shapes (the worst corner and the largest sampled crossing). **So the knife-edge has two blades:** the true liver half-life, which the two data requests address, and the carnitine channel's shape, which no data in the repository addresses. At the round-17 liver estimate of ~4.3 d the 21-day step clears the line everywhere in the box (worst case 4.08 d), so if that estimate is right the shape constants stop mattering for the answer. They still matter for confirming it: against the worst-corner line one patient's record places a 4.3 d truth above it at 2σ with at most about a 47 % chance, even after five years of weekly draws (R-R). The lead figure is drawn at K_c = 100 and its title says so.

**Addendum: the paper's own ranges, and the admissible increment (`t53`).** The box above is `collapse_check`'s, which reaches M_i = 100. Over the paper's own ranges (`tab_C_shape`: K_c 25–800, M_i 150–1200, ρ₀ 0.35–1) the corners move θ_crit^int from 110 to 232 µM/d, a factor of **2.1**, against the factor of 1.4 the paper's remark on the shape constants reports for the same joint sweep. The 21-day crossing there is 3.29–3.95 d (3.6 d clears it at 6 of 8 corners, 4.3 d at all 8), and the 14- and 15-day crossings at most 2.95 and 3.10 d. The paper's scope section also finds the largest admissible increment at 5 % of normal to be D = 2672, not the booked 6000; since Δ_trough ∝ D, the crossings there rise to 3.31 / 3.47 / 4.39 d at the nominal constants and 3.00–3.84 / 3.15–4.01 / 4.02–5.01 d over the paper's ranges. At the 4.3 d liver estimate all three intervals clear at D = 6000 for every corner; at D = 2672 the 14- and 15-day intervals still clear (narrowly) and the 21-day interval does not. P(liver above the 21-day crossing) is 64–71 % at the nominal constants, 56–60 % at the worst corner of the paper's ranges, and 48 % at D = 2672 (nominal). The technical note for the PI (`pi_note/`) states the dosing result as this specification. *(Round 20: same-rate scaling of D = 2672 to 15 or 21 days exceeds the paper's ceiling on peak enzyme; with the ceiling applied at each interval the crossings are 3.30 / 3.54 / 4.95 d at the nominal constants and 2.98–3.87 / 3.19–4.15 / 4.47–5.80 d over the paper's corners, `t53`; the PI note uses these.)*

---

## The original tasks

### 1. Identifiability (`t1_identifiability.py`)

![task 1](fig/t1_identifiability.png)

- **U-curve, 1–1000%.** 14.75 at 1%, 3.07 at 5%, 0.17 at 100%, minimum **0.1226 at 212%**, 0.29 at 1000%. The band is 82.4–627.1%.
- **Mechanism.** corr(ln θ, ln k_e) is −0.991 to −0.992 at 1–10%, −0.967 at 50%, −0.879 at 100%, −0.583 at 200% and −0.047 at 1000%. corr(ln θ, ln K_m) rises from +0.21 at 50% to +0.98 at 1000%, where m ≪ K_m and the record sees only θ/K_m. corr(ln θ, b) stays between −0.13 and +0.02.
- **Variance attribution at 5%, with each nuisance known:** none 3.07; k_e 0.388; K_m 2.885; b 3.046; k_e + b 0.298; all three 0.194. The full curve of k_e's share is in F1.
- **Elasticity.** The analytic form matches numerical differentiation to 3.3e-10 untreated, and to 9e-10 per draw on therapy.

  | | s | m/K_m | ε | d ln m/d ln θ |
  |---|---|---|---|---|
  | normal | 0.845 | 0.100 | 0.909 | −0.9154 |
  | 5% | 0.163 | 0.540 | 0.649 | −0.1729 |
- **Profile likelihood at s = 0.16**, 24 records.
  - **Lower side: no lower bound.** All 24 intervals run to the grid floor (θ = 0.05), and the MLE has median 3.8.
  - **Upper side: a wall.** The 95% upper end is at median **672 (22.4% of normal)**, where k_e has been pushed about 1.5 prior sd below its centre (1.6 on the noiseless profile; `t4`).
  - **The CRB is the wrong shape.** It implies a symmetric factor of 166,704; the likelihood is instead a one-sided plateau with a cliff about 4.5× above the truth.

### 2. Joint fit (`t2_cohort.py`; corrected by F2 and F6) — **closed as a negative result (R14)**

![task 2](fig/t2_cohort.png)

The 24 records are one diet history each of the same patient. Pooled as one patient with 600 draws they give 1.09 at 5%.

**A real cohort, FIM.** Each patient has their own θᵢ and load-diary bias bᵢ. K_m is shared, and ln k_e,ᵢ = ln k_e,pop + ηᵢ with ηᵢ ~ N(0, ω²). The population prior is σ_pop² = 0.8² − ω², and N = 1 reproduces n4 exactly. On q14's composition this agrees with `q14_hierarchical.py` within 4%.

| patients at 5% | ω = 0 | ω = 0.2 | ω = 0.33 | ω = 0.5 |
|---|---|---|---|---|
| 1 | 3.08 | 3.08 | 3.08 | 3.08 |
| 16 | 1.33 | 1.59 | 1.92 | 2.38 |
| **24** | **1.17** | 1.46 | 1.82 | 2.31 |
| 64 | 0.81 | 1.19 | 1.62 | 2.18 |
| 1024 | 0.43 | 0.98 | 1.48 | 2.08 |

**Correction.** This section first called the FIM pessimistic on the strength of `t2`'s Monte Carlo. The robust spreads stand: about 0.88 at N = 16 and 0.52 at N = 64, confirmed by 100 and 200 replicates in `t11`. What was wrong was the RMSE. `t2`'s N = 64 RMSE of 0.49 came from 60 replicates; with 200 it is 1.33, above the FIM. The robust core is narrower than the FIM for the reason tested in F6 (the other patients' θ ≥ 0, acting through the shared k_e). The FIM is 0.81 on `t2`'s seeds 0–63 and 0.78 on the Monte Carlo's own records.

**What does not change:** even exact k_e leaves 0.39 on 25 monthly draws, so borrowing cannot reach 0.20 without more draws per patient (F1). At ω = 0.33 borrowing still cuts the sd from 3.08 to 1.82 at N = 24, against 1.17 at ω = 0, but it stays far from 0.20.

### 3. Variance decomposition (`t3_variance.py`)

![task 3](fig/t3_variance.png)

Exact split: Var(ln θ̂) = Var(ln θ̂ | k_e) + β²·Var(ln k_e).

- **At the band edge** (θ = 2472) the k_e term is **84.3%** of the variance. The sd(n) fit puts the floor at 83.9% at n = 25. The implied Var(ln k_e) is 0.638, the prior itself.
- **15× sampling** (every 2 days, 365 draws) moves the edge from 82.4% to 75.1%. That is beyond the k_e-held limit of 79.9%, because dense sampling teaches k_e a little: Var(ln k_e) falls 4.3%.
- **Correction.** I wrote that n·Var(ln θ̂ | k_e) stays between 0.115 and 0.168 for n = 6–400. That was computed on n4's GRID prefixes, which go untreated after the repository's 200-dose cap. With therapy continued, the values at n = 100 and 400 are 0.178 and 0.197 (`t8_floor_recheck.py`). The floor fit and the share are essentially unaffected (0.183, 83.8%).
- The edge moves with the k_e prior: 43.7% of normal at a prior sd of 0.4, 24.6% at 0.2 and 15.6% at 0.1.

### 4. Set identification (`t4_set.py`) — **elevated to a main result (R16)**

![task 4](fig/t4_set.png)

**θ_max = (L − k_e,lo·m)(K_m,hi + m)/m**, which increases with K_m.

| k_e floor | K_m ≤ 100 | 107 | 150 | 200 | 500 | unbounded |
|---|---|---|---|---|---|---|
| 0 (non-negativity only) | 30.7% | 32.1% | 40.6% | 50.6% | 110% | ∞ |
| 1.04 (the 0.8 prior's 95% low end) | 25.3% | 26.5% | 33.5% | 41.8% | 91% | ∞ |

- On the noisy dose-train records, a constrained profile with k_e > 0 and no k_e prior gives a median of **907 = 30.2%**, finite in all 24 records. It is finite only because the dose train pins K_m to ~7% (sd(ln K_m) 0.074 with no k_e prior).
- **The reported ~25% is the k_e prior's lower tail**: 672 (22.4%) with the prior.

### 5. Value of additional measurements (`t5_value.py`; spot-checked in F6)

![task 5](fig/t5_value.png)

The full printed output of `t5_value.py` (`out/t5_value.txt`), at θ = 5% of normal on the current monthly record. These are FIM values; F6 has the Monte Carlo.

**The absolute values are lower bounds.** At 5% the FIM in ln θ is optimistic wherever θ̂ can reach 0 (F6), so every sd below is a lower bound on the achievable RMSE. Use them for ranking; read verdicts near 0.20 off the Monte Carlo column of the scenario table.

| measurement | sd(ln θ̂) |
|---|---|
| none (current design) | 3.07 |
| + one untreated baseline marker draw | 2.31 |
| k_e to a factor 1.05 / 1.1 / 1.2 / 1.3 / 1.5 / 1.75 / 2.0 (1 sd) | 0.445 / 0.576 / 0.895 / 1.21 / 1.77 / 2.32 / 2.76 |
| the same factors read as 95% (1.96 sd) | 0.403 / 0.444 / 0.568 / 0.710 / 0.991 / 1.30 / 1.57 |
| k_e exact | 0.388 |
| load measured to 1% / 2% / 5% / 10% / 20% | 2.70 / 2.77 / 3.07 / 3.41 / 3.75 |
| load measured exactly | 2.67 |
| one untreated flux reading to 5% / 10% / 15% / 20% / 25% / 30% | 0.071 / 0.112 / 0.159 / 0.206 / 0.255 / 0.303 |
| three flux readings averaged, same Y | 0.051 / 0.072 / 0.097 / 0.123 / 0.151 / 0.179 |
| k_e 1.2 + load 2% | 0.832 |
| k_e 1.05 + load 1% | 0.305 |
| k_e exact + load exact | 0.204 |
| flux 20% + k_e 2.0 / 1.5 | 0.206 / 0.205 |
| flux 30% + k_e 1.5 | 0.300 |
| flux 30% + load 2% | 0.301 |
| flux 20% + k_e 1.2 + load 2% | 0.199 |
| flux 30% + k_e 1.2 + load 2% | 0.286 |
| flux 30% + k_e exact | 0.226 |
| complete-deficiency anchor (k_e from (L − b)/m), ω = 0 / 0.1 / 0.2 / 0.33 | 0.45–0.50 / 0.63–0.67 / 0.99–1.01 / 1.49–1.51 |

- **Complete deficiency.** The repository holds no θ ≈ 0 record to read k_e from, and in the parent model (`models.ParamsC`) a θ = 0 patient has no steady state at the baseline load. The anchor row is what such a record would give if it existed.
- **Ranking:** flux ≫ k_e ≫ load. The elasticities of sd(ln θ̂) to each measurement's own precision are 0.90, 0.79 and 0.15.
- Only the flux reading reaches 0.20 alone: at Y ≤ 19.3% by FIM and Y ≤ 18.4% by Monte Carlo (`t12`), or ≤ 16.0% if its recovery factor is known only to 10%.
- **Conditions on the flux reading:** it must be taken off therapy, since on therapy any PCC flux measurement sees (θ + Δ)g, and it must be absolutely calibrated.

### 6. Dynamic estimators (`t6_estimators.py`; extended by F4 and F5)

![task 6](fig/t6_estimators.png)

- **λ = 5.6–5.9 /d at 5%** (1/λ = 4.0–4.3 h), with θg′ 11–16% of it.
- **The brief's integral form (b) cannot identify θ**: θ enters only through L₀ − θ.

| estimator | ln θ̂ on TB-A: median bias, RMSE | ln θ̂ on TB-B: median bias, RMSE (robust) | ln k̂_e median bias, A / B |
|---|---|---|---|
| (a) FD, left point | −4.6, 4.9; θ̂ ≤ 0 in 56% | θ̂ ≤ 0 in **100%** | +0.21 / +0.33 |
| (a) FD, midpoint | −4.0, 4.4; θ̂ ≤ 0 in 44% | θ̂ ≤ 0 in **100%** | +0.03 / +0.33 |
| (a) FD estimated by GTLS with the MA(1) covariance (F4) | fails (no information) | −0.016, 0.090 (0.084) | — / +0.004 |
| (b) integral, saturated (brief's) | θ not identified | θ not identified | **+0.11 / +0.23** |
| (b′) integral, general, K_m profiled | +9.8, 9.5 | +0.03, 0.32 (0.24) | +1.7 / −0.01 |
| (b″) integral, cumulative | +9.3, 8.5 | +0.04, 0.21 (0.19) | +0.7 / −0.01 |
| (c) NLS through the ODE | +1.1, 4.3 | **−0.01, 0.079 (0.085)** | 0.00 / 0.00 |
| (c) iterated EKF smoother | +0.11, 1.48 (its prior) | −0.01, 0.33 (**0.082**) | 0.00 / 0.00 |

The one-pass EKF is dropped (F5). TB-B's FIM sd is 0.089.

- **The FD bias is endogeneity of ψ**, not noise in m. The coefficients relative to truth are (−1.59, −2.06, 1.39, −2.64) with noiseless m and (−1.60, −2.07, 1.39, −2.64) with noisy m.
- **Discretisation bias:** on noiseless TB-B samples the left point unlumps to θ = 56 against 600; the midpoint gives 600.3.
- **MA(1):** −0.47 to −0.49 in the equation error at the true coefficients, and −0.30 to −0.33 in fitted residuals (`t6b_ma1.py`).
- **Correction.** I wrote "conditioning is not the issue". It is: the noiseless FD design is nearly collinear (F4).

### 7. Experiment design (`t7_design.py`; extended by F3)

![task 7](fig/t7_design.png)

The washout at 5%, hourly for 48 h, down the knowledge ladder (FIM):

| design | nothing | k_e 0.8 | + K_m from the record | + L₀ to 5% of L | + L₀ to 1% | k_e exact |
|---|---|---|---|---|---|---|
| withdraw half the load (the diet) | 383 | 9.4 | 4.1 | 1.5 | 0.63 | 0.10 |
| load ×3 for 24 h, then out | 9.4 | 7.2 | 0.86 | 0.58 | 0.24 | 0.08 |
| load ×6 for 24 h, then out | 3.1 | 3.0 | 0.57 | 0.40 | 0.18 | 0.07 |

- **Corrections from F3.**
  - I wrote that a washout needs *all of* an excursion past K_m, L₀ known to 1%, K_m from the record and the k_e prior. Too strong: with dense enough sampling, ×4 passes joined to the record with no bought knowledge (FIM).
  - I wrote that "a pure withdrawal never reaches 0.20 short of exact k_e". Wrong in the wider grid: a 50% withdrawal reaches 0.195 by FIM at 0.5% CV with L₀ to 1%. The Monte Carlo gives 0.237 on 40 replicates and 0.199 [0.179, 0.219] on the recheck's 150: it straddles.
  - The binding unknowns in a washout are still L₀, then k_e.
- **The τ tension.**
  - For the efficient estimator at fixed duration, denser is never worse.
  - At a fixed N = 24 the optimum is τ ≈ 0.3/λ.
  - Derivative estimators are best at 0.88/λ (forward) or 1.57/λ (central), with the derivative still 43% or 19% wrong.
- **Excursion.** As the peak m goes from 0.68·K_m to 7.4·K_m, cond[1, m, m²] falls from 76 to 18, and sd from 2.94 to 0.48 (k_e prior + K_m from the record).
- **Drug forcing** (13 h rise, 3.5 d decay, no added substrate, k_e untouched). Results are FIM; Monte Carlo is in F3.
  - The burst plus record gives 0.86–1.26 at admissible peaks with L₀ to 1%.
  - It reaches 0.20 if k_e is known to 2.5% with κ known, or 2.4% with κ to 12%.
  - The booked dose (peak 6400) is inadmissible at 5% and worse (1.69).

---

## Unverified or assumed

**Repository footing that this work inherits and cannot check**

- There is no patient data. Every "patient" and "record" is simulated on illustrative parameters (README).
- The marker is treated as the enzyme's substrate. The paper names this its weakest assumption.
- The noise model is the repository's: one flux-space σ per record, evaluated at the untreated level. Under per-draw noise, verdicts within a few per cent of 0.20 can flip. For example, k_e exact plus load exact at 25 monthly draws gives 0.2044 on the repository's noise and 0.1933 per draw.
- Verdicts within a few per cent of 0.20 also depend on the intake seeds and on sample timing against the 14-day dose phase. With k_e and the load exact, 25 draws at 29.2 d give 0.198, at 30.4 d 0.204, and at 28 d (2 × 14 d) 32.6.

**Assumptions I made, with their sensitivity where checked**

| assumption | detail |
|---|---|
| what the 0.8 k_e prior is | Shared-but-unknown, or between-patient spread ω: the repository does not say. Decisive for Task 2 (1.17 at ω = 0 against 1.82 at ω = 0.33, N = 24). R1's "a third or more" is qualitative, for a comparable disposal rate in healthy adults, and for the marker's clearance rather than k_e/α. |
| "k_e known to a factor F" | Read as a 1-sd factor, so sd(ln k_e) = ln F; "k_e to p%" as ln(1 + p). The 95% reading gives 0.40 rather than 0.44 at F = 1.05. |
| unrecorded load in transients | L₀ = 0.5 L and constant through the experiment. ±0.2 L moves the best washout by under 4%. Constancy under a glucose-protected withdrawal is unverified, and the repository withdrew the fast because catabolism raises the load (`u15`, `u17`). |
| load-excursion exposure | Washout and load designs hold the changed load for half the duration, so 72 h designs carry 36 h of load and the safety class counts only the multiple. With the hold capped at 24 h, a few edge passes become misses (`t10_factorial_recheck`). |
| clinical admissibility | Ranking load increases as unsafe, and placing them beyond the parent model's fold (1.18× at 10% capacity, `boundary.tex`), is a judgement drawn from the repository, not a clinical assessment. |
| assay CV below 5% | 0.5–2% assay CVs have no source in the repository. The "5%-assay equivalent" sample counts assume purely analytical, independent noise. |
| drug kernel | A difference of exponentials with a 13 h peak and a 3.5 d decay, my construction around the peak time from `q1`. The admissible ceiling (peak 2850) is n4's R2 reading, one of three that disagree. κ exact unless stated. |
| the direct flux reading | One untreated visit, off therapy, absolutely calibrated, Gaussian with sd F·√(Y² + (ε·5%)²), with m₀ known as in n6. The repository's breath test has an unknown recovery factor, under which one reading is worth nothing. |
| long records | The 27-year (328 monthly) and 500-year (6000) records assume one patient with one constant load bias for decades. They show there is no floor mathematically; they are not designs. |
| the dose train | `design.run.delta` stops after 200 doses. `t8` continues it in closed form; that continuation is mine. |
| the "K_m from the record" rung | Applied on its own with the k_e prior, it counts that prior twice. That is conservative: it created no false passes in the cases checked. |
| cohort composition | Task 2 and F2 put every patient at 5% (and q14's 5–20% mix), with each patient's load bias independent. |

**Explanations that are supported but not proven**

- **The EKF failure** is attributed to one-pass linearisation error growing with the prior width. That is supported by the exact information identity at the truth and by the error's growth with prior width; there is no second-order analysis. The drift of outliers towards the Michaelis–Menten limits is common to all records, so it does not explain the outliers.
- **GTLS-t's success** is read as MA(1) whitening restoring the level information that differencing removes. That is an interpretation, consistent with the whitened partial SNR (3.1 against 0.0016).
- **Why b is so weakly separated from θ** (its constant column against θg(mᵢ), and the limited spread of mᵢ) is plausible, not tested.
- **The linear-Gaussian surrogate at N = 64** has a whitened sd of 0.984 over 8000 replicates, just outside its 95% band. It is provisionally sampling noise.
- **That the iterated smoother converges to the MAP** is not shown. No out/ file compares its final cost with a multi-start least-squares fit, so its bench-B RMSE (0.33, against NLS's 0.079) may include convergence failures as well as the smoother's own error.
- **The ω profile-likelihood hierarchy** in `t9` uses each patient's profile, not integrated, likelihood, and holds it constant below the grid.

**Comparisons with limited support**

- The flux crossing (18.4% by Monte Carlo) is for one simulated patient on the repository's footing. The reading is Gaussian and taken with an untreated marker draw, and either perfectly calibrated or with a recovery factor known to 10% (16.0%). Whether a measurement of this kind and precision exists for PCC flux in vivo is not established here.
- The repository's only breath-test repeatability figure, a between-visit 14% (`PHASE4_SUMMARY.md:208`), is for an uncalibrated relative reading. It does not show that a 10–19% *calibrated* flux measurement exists.
- **The repository's breath test is not the flux reading modelled here.** `u2_breath_functional.py` models the 1-¹³C-propionate breath test as y = c·F(m), where F = θg/(θg + k_e·m) is the *fraction oxidised*, not the enzymatic flux θg(m), and the recovery c is unknown. So neither the 18.4% crossing nor the 16.0% recovery variant transfers to that test without further work. What Task 5 needs is an absolute in vivo flux measurement, and none is identified in the repository.
- Monte Carlo sizes: 20–40 replicates per cell in `t10`. The verifier re-ran all five designs on fresh draws (300 replicates for drug 2850 and load ×2, 150 for the others), and several verdicts moved: withdrawal 0.237 → 0.199, load ×2 0.237 → 0.217, drug 1500 0.113 → 0.131, load ×3 0.203 → 0.211. The report uses the recheck's figures. `t9` uses 60 replicates per ω. Borderline verdicts carry their bootstrap intervals in the table.
- **`t12_flux_mc.py`, the source of the 18.4% headline crossing, has no independent recheck.** It reruns `t11`'s estimator on `t11`'s own draws at Y = 10% and 20% and reproduces them exactly, and its FIM matches `t5` at every Y. The 2000-replicate fresh-draw crossing itself has not been reproduced by a second script.
- Repository-sourced numbers (the dose train, LOAD_CV, the 1.177 fold, the breath-test 14%, the ln-m/flux FIM agreement, and the saturation range under therapy) are printed by `t0_footing.py`. n6's 708 = 24% is in `out/n6_profile_vs_crb.txt`, from an unchanged run of the repository's script.
- Assumption-free conclusions are limited to the structure: the confounding mechanisms, the lack of an n-floor, the non-identifiability of ω from monthly records, and the FD mechanism. Every absolute sd depends on the repository's illustrative constants.

**Round 3 (`t13`–`t19`), assumed or of limited support**

- **Round-3 verification.** Each round-3 script was re-derived by an independent verifier (`*_recheck3.py`), by a route sharing no code with it. Items 9, 11, 12, 15, 16 came back **confirmed** to the printed precision (the item-11 hypot collapse was additionally shown to be an exact marginalisation identity; the item-16 numbers reproduced with an independent quadratic steady-state solve). Two came back **partial and were corrected above**: item 13 (the scatter-vs-bias ordering — a `core.fim` default-prior coupling that pinned the bias in the scatter-off run; fixed, the bias is now correctly the more binding term) and item 10 (the excursion's "~100% power" is largely an errors-in-variables artifact; the size-corrected result is in the R10 table). The verifier also noted that item 13's sd·√n drifts slowly to ~4.6 rather than sitting flat (the 1/√n rate and the no-floor conclusion are unaffected), and that the item-11 Monte-Carlo 18.4% is a `t12` figure not reconstructible from FIM primitives.
- **R10 is a power analysis, not a model test.** The misspecified data-generating processes (second compartment, saturable non-enzymatic clearance, time-varying load drift) are my constructions; the reported rejection rates depend on their parameterisation and on the particular residual-slope test. The claim is the **direction and ranking** (monthly record blind to the first two, sensitive to a load drift; an excursion restores power), not the exact percentages. With no real data the specification cannot actually be tested.
- **R11's (s, L) reparameterisation treats the fraction oxidised as absolutely measurable to Y%.** If the available reading carries an unknown recovery (as the breath test does), that route is back to the calibration problem and needs c calibrated to Z. The √(Y²+Z²) contour is a first-order FIM result in linear θ; near 0.20 it inherits the same ln-θ optimism as the rest (F6), so read it as a lower bound.
- **R15's safety threshold (1.177×) is a modelling construct**, the saddle-node of the parent model's R0 setting, not a measured clinical crisis level (`src/model_A/run.py`). The excursion-sweep figure uses a generous knowledge rung (k_e prior + K_m from the record + L₀ to 1%) and FIM; the authoritative excursion verdicts remain the `t10` Monte Carlo ones. The two-stage k_e precision (~0.2 at L₀ 5%) assumes the unrecorded load L₀ is known to the stated precision and constant through the transient (`dyn.py`, as in Task 7).

**Round 6 (`t29`–`t32`): literature cross-checks, borrowed figures flagged**

- **The literature figures are from different populations and assays than the repo model, and are used to test assumptions and sweep, never to replace a repo value.** Specifically: the in-vitro/in-vivo mapping R² = 0.327 is **MMA, not PA**, and a ¹⁴C-incorporation functional assay, not a direct PCC assay (Manoli 2021, n = 11) — *no PA/PCC same-subject study exists at all*; the between-patient capacity spread (~0.8 log-units) is from Thompson 1990 in-vivo oxidation (n = 6, PA+MMA); the healthy ex-vivo spread bound (~16%) is Gravel 1977 fibroblast PCC, possibly a min–max range; the breath-test within-subject CV is **not a published number** (Shchelochkov/Thompson document reproducibility only qualitatively), so it is swept 5–10%. The conclusions that depend on these (ratio route fails for whole-body capacity; breath test wins monitoring) are therefore **directional and conditional**, flagged in-text, and the honest headline is that the in-vitro/in-vivo mapping is *unmeasured in PA* — which is itself the point against the ratio route.
- **t29, t31 were re-derived by an independent verifier** (sharing no code): the mapping-scatter conditional-variance formula, the FIM injection, the breath elasticity (0.888, two independent derivations) and every table cell came back **confirmed exactly**. t30 reuses the verified `t24` drift machinery; t32 is a steady-state partition and a provenance check.
- **t33/t35/t36 (Model B) were re-derived by an independent verifier** (sharing no code): the monotone-disposal/θ_crit existence argument, the five Model-B FIM sensitivity columns (vs finite differences), `level_B`'s root-solve, the A-vs-B share/sd tables, the managed-load θ_crit margin, the Remark-15 reconciliation (share ratio 3.04 ≈ 92% of 3.30; cond ratio 18×; no-prior exchange blow-up to 55), and the b-sharing contrast fix all came back **confirmed exactly**. Two of my own errors were caught and corrected along the way: the **v11 Model-B contrast (0.167) was buggy** — it failed to share the load-bias b between windows; the correct value is 0.088 (B better than A), and there is no cross-sectional/longitudinal asymmetry. And the **"no steady state at 5%" claim overstated it** — θ_crit moves with the managed load, so a 5% patient is stable under management.
- **Round 7/8 is a model-validity round, and its sharpest claim is conditional on parameters the repo does not fix in absolute units.** The decompensation boundary θ_crit = L − k_e depends on both the alternative-route ceiling k_e (ParamsB's 400 is miscalibrated — above the untreated load; the repo's Model C and the literature put it below) and the *managed* load L, which is a clinical variable. So the robust, parameter-free part is structural: *a bounded clearance has a finite disposal ceiling, so an untreated severe patient has no steady state and needs the load managed below it, with a thin margin*. The exact θ_crit and μ need the ceiling in absolute units, which the repo lacks. Redoing the full identifiability under a calibrated bounded clearance is the open model-side step.
- **Round 8 hardening (contract, audit, affine).** The sharing contract is now centralised (`contrast.py`) with a test (`check_contract`) and a negative control, and every FIM contrast in the report was re-audited through it (`t37`): R-D unchanged, only the Model-B contrast moved (0.167→0.088). Residual conditionals worth stating: (i) the "under B a k_e sensor reaches ~0.19" result (`t38`) is a *knife-edge* at the gate and needs K_e essentially known — it is prior-bought, and with the priors stripped B is ~53, worse than A; the model-dependence is real but the cross-sectional sign (misses 0.20 before a k_e sensor) is robust across the ceiling×load grid (0.21–0.49). (ii) The affine check (`t39`) uses ε = the per-draw flux noise as the "measurement error"; a different reading of ε would shift the pass/fail margin, though the record (0.3 ε) vs excursion (7 ε) gap is wide enough that the qualitative verdict is insensitive to the exact ε. (iii) The valid Model-B point is a construction (chosen ceiling and managed load), reported as a range, not a patient fit.
- **Round 8 close (`t40`): the Model-B cross-sectional escape is a point, not a route.** Swept under the tightened conditions that produced ~0.19, it only hugs the gate — at the tracer's k_e prior the whole ceiling sweep sits at 0.18–0.20 (marginal, never comfortably below; sd tracks the enzyme share, flat in the managed load, not the margin), slipping to 0.21–0.23 at higher ceilings under a looser prior — and only while K_e is known to ~30%. K_e is a bare repo default with no provenance entry and no human measurement, and assuming it known collapses the paper's (θ,K_m)↔(k_e,K_e) exchange degeneracy by fiat. The surviving, model-independent statement is the three-constant pattern (κ, c, K_e): the cross-sectional target is reachable only by importing a constant unmeasured in humans.
- **Round 9 (`t41`): Model C, the structural-law test, is independently verified.** The two-observable (marker + free carnitine) FIM under Model C was re-derived by an independent verifier (own two-state solve, own sensitivities) and every number reproduced: m-only alt-free ~20, m+c alt-free ~16, m+c alt-known 0.34; the managed load 500 is conservative (loads nearer the fold give *smaller* sd, so it does not flatter the result); and freeing R₀ and k *together* jumps sd to 2.48 (the R₀/k degeneracy in c), confirming one of them must be known. Residual conditionals: the 0.34 known-channel figure still sits just above 0.20 and rests on K_c/R₀/M_i/ρ₀ being known, which they are not (arbitrary / reverse-engineered / one-sided-bounded); and the structural claim (no third kind of θ-route) is an argument checked exhaustively against the ten-sensor panel, not a theorem.
- **Round 10 (`t42`, `t43`): the empirical anchors and the compartment defect.** (i) The serial-CV recompute (`t42`) rescales the contrast by the per-draw marker CV; the 26.4% anchor is healthy-adult, same-day, unstandardised data (Thompson 2012) used as an explicit proxy — no PA serial CV exists, so the ~0.47 / ~590-draws result is "what R-D becomes *if* PA's serial CV matches the only measurement," stated as a range, not a point. (ii) The plasma-pool defect (`t43`) models the plasma→pool partition φ as a free (and m-dependent) scale on the carnitine observable; the finding that φ barely moves the *channel-known* case (c redundant there) and erodes the *channel-free* case (16→17.7) is robust, but the m-dependent-φ form is illustrative (a linear load-drift), not the measured partition function, which is unknown. The literature constants for Model C (biosynthesis, renal threshold, CrAT Km) are in physical units with no bridge to the model's abstract scale, so they were not literally substituted; the t41 sweep already spans the plausible ranges and the conclusion rests on known-vs-free, which the literature settles against "known."
- **Round 11 (`t44`): the interval-decision result rests on three repo conventions, each flagged.** (i) **λ is ASSUMED** — *(round-15 update, R-N: λ is now measured in humans at 3.6 d (Chandler 1985) and transferable; the sensitivity sweep is no longer hypothetical but a map of where the measured value and its 3-fold band fall — and they fall on the thresholds)*; the headline's robustness is a margin across the 1.2–11.5-day band, not a pass at a single λ. (ii) **θ_crit carries the unmeasured s_pk/m_tol convention** (the parent paper's own complaint), so Model A is reported at two folds (1.0 routine / 1.5 illness) that bracket it; the "needs no new measurement" headline is strongest at the routine-load reading and under Model C's convention-free intrinsic margin, and softens to a boundary problem at fold 1.5. (iii) The **same-rate vs fixed-amplitude** dosing regime is a modelling choice; same-rate is the repo's §7 convention (verified by the 88 % trough cut at 28 days) and is used as the conservative primary, but neither is a measured protocol. The population fractions (item 3) are over a log-uniform θ_base prior on 0.17–5 % of normal, an illustrative population, not a cohort; and the Model-C and Model-A troughs share the one µM/d flux scale, which (as throughout) has no bridge to absolute units.
- **Round 17 (`t50`): the compartment correction is an anchor shift, not a new measurement; the direction is better-grounded than round 16's, the band is not.** The reversal rests on one compartment fact (PCC is mitochondrial-matrix) plus the liver-mito median (4.26 d, Kim 2012, mouse) exceeding the fibroblast value — a qualitatively robust direction, but the central 4.3 d is still a cross-species, cross-method median applied to one protein, with no PCCA/PCCB number in any supplementary table. The P(below) ~30 % uses the same wide band as round 16 with the centre moved to 4.3 d; it should be read as "the 21-day step is more likely free than round 15 said, band straddling," not as a settled probability. The "human λ is not public" finding is a literature-access statement (Gao poster gated, no full publication as of Oct 2026, FOIA Exemption 4), not a claim about the model's internals, which remain unseen.
- **Round 16 (`t49`): the hepatic bridge is a transparent log-normal propagation of literature medians, not a measurement.** Its inputs (human fibroblast 3.6 d; mouse liver/fibroblast medians) are cross-study, cross-method, and culture-vs-in-vivo confounded, and the per-protein-to-tissue-median scatter is set by assumption (sd(ln) 0.35 each on ratio and per-protein), deliberately wide. The headline is not the point estimate but the *spread* — P(below the crossing) ranges 24–65 % across the three ratio readings, which is the finding (the bridge can't settle the side). The ratio-reading centres (1.44 / 1.0 / 0.8) are judgement calls from the cited medians, not fitted; a reader who weights the culture-dilution confound differently would shift the centre but not the conclusion that the band straddles 50 %. The crossing (3.5452 d) and the discrimination σ are carried from rounds 15/12 unchanged.
- **Round 15 (`t48`): the λ values are literature from other populations/assays, and the 2-comp mapping is reconstructed.** The human 3.6 d is *skin fibroblast* (not hepatic, the decision-relevant tissue, which is unmeasured and — per R-P, round 17 — likely *equal or slightly slower*, not faster); the mouse 7 d is a two-compartment *effective terminal* whose exact equations are unpublished (Suppl Table S1 unavailable), so the "trough sees ~3–5 d" effective-λ is computed on two bracketing symmetric reconstructions and should be read as a range, not a point. The decision table and fail half-lives (2.62/2.75/3.32/3.55 d) are exact given the same-rate trough and θ_crit conventions (which carry their own caveats, R-J) and, for the Model-C crossings, Model C's unmeasured shape constants K_c, M_i, ρ₀ (R-S). The discrimination d′ uses the round-12 on-therapy sd(ln λ); it says the record *can* distinguish 3.6 from 7 d, not that the trial has done so. These are used to replace the ASSUMED provenance and to test conclusions, not substituted silently into the model's abstract scale — the µM/d flux units still have no bridge to physical PCC concentration.
- **Round 14 (`t47`): the audit is a spot-check, not a full re-simulation of every number.** Two natural-scale claims (round-11's determined-subgroup, the cross-sectional placement at 5%) were re-derived by NLS and survive; R16's recoverable quantities are flagged at-risk but *not* re-simulated (they are explicitly order-of-magnitude). The correlation scan is at one operating point (θ=150, 24 mo weekly); the two strong correlations (θ–k_e, θ–λ) are structural and hold across the range, but the scan does not prove the *absence* of a weak correlation feeding some minor natural-scale quantity — it establishes there is no *strong* one. Spot-check B's sd(ln θ̂) ≈ 8.3 is the λ-free 5-parameter NLS and so overstates the λ-fixed headline sd (3.05); the point it makes — that θ is unplaceable on the natural scale at 5% — is direction-robust, not a revision of the 3.05 figure (which the scenario table already backs by Monte Carlo).
- **Round 13 (`t46`): the 21-day negative is simulation-based and robust; the λ timing is FIM.** The headline — that the 21-day margin estimate does not tighten with record length (sd(μ̂) ≈ 140 µM/d flat from 6 months to 5 years) — comes from NLS simulation (200 records/point), which is the right tool here precisely because FIM propagation is unreliable for this tail-dominated one-sided quantity (delta-method and log-normal MC disagree by ~400 µM/d; the parent paper's §7b warns of exactly this). The NLS fits can land in heavy-tailed minima (as F6 documents for θ), which is part of *why* μ̂ does not tighten — a real property of the estimator on this model, not a fit artefact (the same qualitative failure appears across seeds and capacities). The λ-timing table is FIM-based and trustworthy (λ is in the Gaussian regime, θ-independent sd, confirmed by the round-12 verifier's eigenvalue check). The event schedule, θ_base population, and Model-C/routine-load conventions carry the same caveats as R-J/R-K. The "needs the managed-load protocol (0.29)" alternative is asserted from the brief's own managed-load number, not re-simulated here.
- **Round 12 (`t45`): the on-therapy results are FIM lower bounds on a correctly-specified record, and the event schedule is illustrative.** (i) sd(ln θ) and sd(ln λ) are Cramér–Rao bounds (as elsewhere in the brief); the real NLS would be no better and likely worse on the low side, as F6 showed for θ. The λ self-calibration (sd(ln λ) ≈ 0.055 monthly) is the strongest new claim and assumes the dose schedule is *known and dated* — which it is for an exogenous protocol, but it also assumes the marker is generated by this one-compartment model with the increment entering exactly as (θ+Δ)g(m); a second compartment or a time-varying λ would degrade it. (ii) The event schedule (2 d jitter, 8 % missed doses, one 28 d supply gap) is a plausible but invented perturbation set; the gains scale with how much the operating point actually moves, so a more regular real schedule buys less and a more disrupted one more. (iii) "Alternate-day sampling" is a design point, not a clinical recommendation — it is included to show θ is sampling-limited, not information-limited. (iv) The determined-subgroup and lengthen-eligible fractions use the same illustrative θ_base population and the Model-C intrinsic θ_crit (routine-load for the subgroup), inheriting the same convention caveats as R-J.
- **R-H (the clearance-structure conflict) is a model-validity flag, and the re-fit it calls for is only partly done.** It shows Model A's 84 %-non-PCC split at 5 % conflicts with the PA literature and that a bounded route (Model B/C) is the literature-consistent form. Later rounds re-ran the cross-sectional identifiability under Model B at its repository constants, at a constructed literature-consistent point and over a grid of ceilings and loads (`t33`, `t35`, `t36`, `t38`, `t40`), and under Model C at its repository constants plus a one-at-a-time sweep of the channel constants (`t41`), and the cross-sectional verdict kept its sign. The monitoring contrast was recomputed under Model B (0.088) but not under Model C. What is not done is the step R-H names as the one that matters: the identifiability under a bounded clearance with a literature-calibrated ceiling and the crisis boundary treated explicitly. No absolute µmol/kg/h comparison is possible without a unit bridge the repo lacks.
- **R12 and R15's "draws per year" use a fixed two-year horizon** (n draws over 730 days), so "365 draws" is every-other-day for two years. The per-draw noise is the repository's; verdicts within a few per cent of 0.20 shift under per-draw noise and sample timing, as elsewhere.
- **R16's "recoverable" sds are first-order FIM propagation** on the monthly record's covariance; the θ-proportional quantities (flux, share, a) inherit θ's ill-conditioning and are not Gaussian there. The one-sided θ bound is exact algebra; the ~6%/28%/8% figures are local.
- **Round 19 (`t52`): the record-resolution numbers are per-patient FIM sds (`t46`'s, re-evaluated at each true half-life), not simulated.** The sds themselves rest on `t46`'s statement that ln λ is in the Gaussian regime (confirmed, per the round-13 bullet above, by the round-12 verifier's eigenvalue check), which was made at the repository's 3.5 d only; the 6-month monthly record is just 6 draws. The side-of-crossing chances Φ(z − 2) are one-sided threshold decisions, which R-M's rule says to read by simulation; these are FIM-based and have not been simulated. Two things are assumed. (i) z is the expected separation; the chance a single record clears 2σ is Φ(z − 2), which the R-R table reports. The sds are `t46`'s FIM re-evaluated at each true half-life (at 4.3 d this moves the 6-month monthly sd from 0.094 to 0.100). The pooled-patient counts (23–77 at 50 % power, 46–155 at 80 %) assume λ is identical across patients; between-patient spread in hepatic PCC turnover is unmeasured, and any real spread raises those counts or makes pooling answer a population mean that no single patient has. (ii) The dead zone is drawn at 2σ; a stricter bar widens it.
- **Round 19 (`t53`): the shape sweep is a sensitivity over the repository's admissible box, not a distribution.** The box is `collapse_check`'s (K_c 25–794, M_i 100–1259, ρ₀ 0.35–1.0), with R₀ re-derived so the alternative route still peaks at 500 against a load of 600. The extremes are taken at the eight corners, on the evidence that θ_crit^int moves monotonically in each constant in the one-at-a-time sweep; an interior extreme is not excluded. The 80 % is the share of 300 draws (log-uniform K_c and M_i, uniform ρ₀), which depends on how the box is sampled. The load (600) and the alternative-route peak (500) are held fixed; moving them moves θ_crit too.

## Appendix: repository material this report has not used

`q15_sensors.py`, `u3_joint_fisher.py` and `PHASE10_SENSORS.md` existed through seven rounds and entered the report only when asked for. This surfaces the rest. The repository holds a large pre-existing audit of the *original paper* (Phases 1–12, ~160 scripts in series a–u plus `PHASE*.md`), of which this report cites ~10. Most of the series are the paper's own checks and rechecks; listed below are the uncited items that bear on **identifiability, sensors, calibration or study design** — the category where `q15`/`u3` were found. Surfaced, not analysed; the "would it change anything" column is a first read of the title/banner only.

**Phase documents (narrative audits, none cited):**

| document | what it contains | bears on |
|---|---|---|
| `PHASE11_KE_MEASUREMENT.md`, `PHASE11_KE_ROUTE.md` | the k_e requirement frontier and what could measure it (the carnitine lever) | the binding nuisance (R-G); may refine the k_e-sensor pricing |
| `PHASE11_SAMPLING.md` | sampling modality: assay CV, home microsampling, what binds | R-F (assay-CV × draws); possibly a less-invasive option I priced as failing |
| `PHASE11_TRADE_AND_PRIOR.md` | the (n, k) trade curve, and "a caveat that attaches to every Fisher figure" | R-D/R-F; the caveat may qualify the FIM numbers |
| `PHASE11_DRIFT.md`, `PHASE11_PROPOSED_CONSTANCY.md` | the constant-θ assumption priced against record length, with a failure threshold | **read and integrated (R-D):** bias ≈ −r·T/2 scales with span not count, so it supports the short-window monitoring result and undermines a decade-long cross-sectional record; a drift column is near-free insurance |
| `PHASE11_KM_BYPASS.md` | the K_m bypass and the breath test's actual functional, with a joint Fisher | R11/R-G; a second take on the breath observable |
| `PHASE11_MARGIN_NUISANCE.md` | "the binding nuisance for the decision is not the one that binds the estimate" | the classification-vs-decision split; may reframe the verdict |
| `PHASE11_CATABOLIC_GATE.md`, `PHASE11_GLUCOSE_ARM.md` | a withdrawn result and a priced substitute; a manoeuvre "with all three properties" | the safe-excursion question (R15/R-E) |
| `PHASE12_GENOTYPE_REVIEW.md` | the genotype-prior document ("what 118 rests on") | **read (R-G, item 2):** genotype gives a *class / predicted range*, not a θ identity to 20%, and carries no carrier data — so it does **not** serve as the κ-calibration anchor; the ratio route (`t28`) makes that anchor unnecessary anyway |
| `PHASE5_ROBUSTNESS.md`, `PHASE7_TABLE2_PROVENANCE.md` | robustness to model structure; provenance of Table 2 | the specification question (R10) and the footing |

**Scripts (uncited; clearly on the brief's questions):** `q2_requirement_map.py` (the schedule/requirement map `q15` is built on), `q2b_ke_measurement.py` and `u9_ke_frontier.py` / `u11_trade_curve.py` / `u13_ke_requirement_audit.py` (the k_e sensor and its (n,k) frontier), `q5_exvivo_prior.py` and `a5_kappa.py` (the ex vivo / κ prior — directly item 1), `t6_cohort_with_kappa.py` (the breath test's per-subject repeatability range — item 1's c), `t9_measured_load.py` (how often a turnover tracer is actually run — R-G's load sensor), `t1_transient_worth.py` (the transient's worth — R15), `n3_ke_prior_worth.py` / `m1_loadbias_prior.py` (what the k_e and load-bias priors are worth), `d10_two_channels.py`, `m3_edt3_discriminator.py`, `k1_sources.py` / `k2_structures.py` / `k5_identifiability.py` (sources, structures and identifiability of the information). The remaining ~140 files (the a–u rechecks and Phases 1–9 scaffolding) are the paper's own audit and are not surfaced individually.

**Honest status:** the two that bear directly on a published conclusion — the Phase 11 **drift/constancy** documents and **Phase 12 genotype** — have now been **read and folded in** (R-D and R-G item 2 respectively): drift supports the short-window monitoring result and undermines a decade cross-sectional record; genotype gives a range not a θ anchor and so does not rescue the absolute ex vivo route, which the ratio route (`t28`) makes moot regardless. The remainder of the list is still surfaced-not-read; the next most likely to matter are the Phase 11 **k_e-measurement** and **sampling** documents (they may refine the k_e-sensor pricing and the assay-CV frontier). That is empirical/model work rather than more of this identifiability analysis.

## Reproducing

See [`README.md`](README.md) for the script list and the run order.
