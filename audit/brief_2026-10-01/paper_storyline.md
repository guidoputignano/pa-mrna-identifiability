# Storyline for a single paper

Working title: *The dose train is the identification experiment: design rules, limits of
recoverability and clinical specifications for residual enzyme capacity in propionic acidaemia.*

All numbers are from the paper (`paper/`) and the audited PI note
(`audit/brief_2026-10-01/pi_note/`), with the audit's corrections applied (REPORT.md).

## 0. The story in one sentence

The therapy is the identification experiment. With the right sampling it identifies capacity in
milder patients. In severe propionic acidaemia no sampling can, but the clinical decisions do not
need capacity: they need the enzyme's half-life, the delivered dose and a short list of
well-recorded data.

## 1. Introduction: the clinical question

- mRNA-3927 is given as an infusion every two weeks. The open question is whether the interval
  can be lengthened.
- The obvious route is to measure each patient's residual enzyme capacity θ.
- The control-theory view: each infusion is a known input acting on the enzyme, so the trial
  already runs an identification experiment.
- The questions: how should blood be sampled; can θ be identified; if not, what can be; and what
  does the clinical decision actually need?

## 2. The model: one balance, one disputed term

- The model is a metabolite balance: the load comes in, and two exits remove it: the enzyme (the
  patient's own plus the delivered one) and the other route.
- The three versions share this balance and differ only in how the other route is written. They
  are three hypotheses about one term, each chosen for a job:
  - **Parent model.** The other route runs through carnitine. It is realistic: it has a tipping
    point (a crisis) and defines the critical capacity. Used for dosing and crises.
  - **Local identification model.** The other route is a straight line. It is the most
    favourable case for telling the two exits apart. Used for identifiability and design.
  - **Saturating variant.** The other route has a maximum speed. It is a robustness check on
    which conclusions depend on the straight line.
- Structural identifiability holds, as the paper proves. The question is practical: what can be
  identified with realistic noise and sampling.

## 3. The dose train is the experiment: how to sample (clinical design rules)

- The infusions supply about 42 times the excitation of the marker's own day-to-day variation.
  The therapy does the identifying.
- **Never draw blood on dosing days.** Draws at the dose interval or a multiple see the same
  phase every time. At normal capacity this costs a factor of about 27 and never reaches the
  target.
- **The fix is to book every 14.5 days instead of 14.** That is half a day off, with the same
  number of visits, and it gives monthly-level precision. The schedule must sweep the whole dose
  cycle (the "lap rule").
- Monthly draws across dose cycles reach the target in milder patients. About a dozen visits is
  enough; after that, extra draws add little.
- **What to record:**
  - the date, dose and lot of every infusion (the dosing record is a regressor, not a nuisance);
  - intake as eaten, to better than about 10% (accuracy matters, not variability; do not smooth
    or impute);
  - antibiotic weeks.
- **What not to do:** change the dosing schedule to help identification (it costs the patient
  far more than it gains), or use draws from the first hour after an infusion in the
  steady-state fit.
- No prior on K_m is needed. K_m is recovered to 6–7% at every capacity.
- **For the clinic:** a concrete, low-cost sampling protocol.
- **Caveat:** these rules give target precision on capacity only in the 82–627% band. In severe
  patients they still serve the half-life and monitoring.

## 4. The rules hold across published models; the value of θ does not

- Records were generated with two published structures, Attarwala et al. 2023 (an algebraic
  PK/PD model of mRNA-3927) and Gao et al. 2025 (a turnover model), and analysed as ours.
- **Robust:** the sampling advice holds under every structure. The schedule ranking is
  unchanged, and the two recommended designs agree to within about 3%.
- **Structure-dependent:** the value of θ. Under the algebraic structure, θ is 1.3–2× off at the
  paper's dose and 4–20× at a trial-calibrated dose, with a fit that looks just as good.
- **Settled cheaply:** one blood draw 3 minutes after an infusion separates the two families of
  structures by about 16 standard deviations.
- **For the clinic:** trust the design rules whichever model is right, and add the 3-minute draw
  to one infusion visit.

## 5. Result 1: capacity cannot be measured in severe patients

- The paper's design reaches the target (about ±20%) only between 82% and 627% of normal. At 5%
  of normal it is 15 times worse (3.07 against 0.20).
- Sampling 15 times more often barely moves the edge. The limit is the physiology, not the
  schedule.
- What the data still say is a one-sided bound: capacity is below about 22% of normal (about 30%
  without the prior on k_e), with no lower bound.
- **For the clinic:** do not report a patient's capacity as a number; at most, report an upper
  bound.

## 6. Result 2: the paradox of direct and indirect measurement

- **Indirectly, through metabolite levels**, you see the total removed, not the split, so you
  must know the other exit.
  - **The A + B arithmetic:** the enzyme does 16% of the work and the other route 84%. A 1%
    error in the large share becomes a 5% error in the small one (the restaurant bill).
  - **The geometry:** each draw is a line, and similar levels give near-parallel lines with no
    clear crossing.
  - Measuring the other exit with tracers still leaves 0.44, against a target of 0.20.
- **Directly, through enzyme activity** (the cell assay or the breath test), you see the enzyme,
  but on an instrument with no known scale between the test tube and the whole liver. The assay
  reaches the target only if that scale is supplied, and nobody has measured it.
- There is no third way: in the model the enzyme acts on the body only through the flow it
  carries.
- **For the clinic:** no available test, alone or combined, gives absolute capacity at 5% of
  normal.

## 7. Result 3: change within one patient is measurable, as a difference

- The same unknown error appears in both periods and cancels in the difference.
  - With two periods of 104 weekly draws, the change is known to ±20 µM/d.
  - A change of about 56 µM/d is detected 80% of the time.
- The change is measured as an amount, not a percentage: the percentage needs the unknown level.
- **Conditions:**
  - the marker's scatter must be below about 10% (unmeasured in PA);
  - k_e and the calibration must be stable between the periods;
  - antibiotic status must be recorded (an unrecorded course fakes improvement).
- **For the clinic:** monitoring is feasible if the marker proves stable. Measuring its serial
  variability is the first study to run.

## 8. Result 4: the dosing interval does not need capacity

- An interval is safe for every patient if the delivered enzyme left just before the next
  infusion exceeds the critical capacity. That depends only on the half-life and the delivered
  dose.
- **Crossing half-lives:**
  - At the assumed dose: 2.6, 2.75 and 3.55 days for 14, 15 and 21 days.
  - For 21 days at a dose capped at the paper's ceiling (the enzyme never above normal): 4.5–5.8
    days.
- The likely liver half-life is about 4.3 days, from comparators; the only human measurement is
  3.6 days, in fibroblasts.
  - At the assumed dose, all three intervals are safe.
  - At the capped dose, 21 days is not guaranteed.
- **Where the half-life can be found:**
  - the sponsor's population model;
  - the trial's own dated record, pooled across patients;
  - a turnover study in human liver cells.
- **For the clinic:** lengthening the interval is a measurement problem, and it is solvable.

## 9. Result 5: how late can an infusion be?

- **Without illness:** each day late adds about 10% to the metabolite's peak, and the peak
  doubles after 6.5–9 days. A skipped infusion raises the peak 3–5 times, with no lasting effect.
- **During an illness:** the peak rises about 2.5 times even with on-time infusions, and 4–5
  times when the infusion is three days late.
- **In the least favourable settings** the peak doubles after about 2.5 days.
- **For the clinic:** a missed day is minor, but during an illness the infusion should not be
  postponed. Dated delays also help measure the half-life.

## 10. Clinical recommendations and considerations

- **Three tiers, run in order:**
  1. one blood draw 3 minutes after an infusion, which settles the model structure;
  2. a same-day pair of breath tests around one dose, whose repeatability must be better than
     about 12%;
  3. the record itself: the sampling protocol and the recording rules of section 3.
- **For the interval:** measure the hepatic half-life and the delivered dose. The specification
  table gives the thresholds.
- **For monitoring:** first measure the marker's serial variability in PA. Keep the regimen and
  antibiotics fixed between comparison periods.
- **For safety:**
  - record antibiotics, since an unrecorded course fakes improvement;
  - give infusions on time during an illness;
  - report bounds, not point estimates, of capacity.
- **Considerations:**
  - these are model-informed specifications, not clinical recommendations;
  - all results are simulated;
  - the illness and the "peak doubled" yardstick are illustrative;
  - the results apply to the stable patient between crises.

## 11. Conclusion

- The therapy is the experiment, and good sampling makes it informative.
- In severe patients, capacity stays out of reach. The decisions do not need it: the half-life,
  the delivered dose and a well-kept record do.
