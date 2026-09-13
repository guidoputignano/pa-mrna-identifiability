# Phase 2D — the three pre-prose checks. Reported before writing.

`2D.0` ("verify before writing"), `2D.3` ("CHECK FIRST") and `2D.10` ("report
before writing prose") all gate. All three are answered. **One comes back against
the handoff and changes what 2D.9 can ask for.**

Scripts: `audit/d0_spectrum.py`, `audit/d10_two_channels.py`.

---

## 2D.0 — exactly rank one. Verified in the four-column form.

Full singular spectrum of Q under constant recorded load, dense design:

| columns | singular values | s_min/s_max | rank | deficiency |
|---|---|---|---|---|
| (θ, K_m, b) | 6208.702, 884.516, **1.09e-13** | 1.76e-17 | 2 of 3 | **1** |
| (θ, K_m, k_e, b) — the footing the paper claims | 6236.548, 895.609, 2.193, **9.85e-14** | 1.58e-17 | 3 of 4 | **1** |
| (θ, K_m, k_e) — no load-bias column | 6236.518, 895.608, 2.193 | 3.52e-04 | 3 of 3 | **0** |

**Exactly rank one, in both forms.** The C2 identity vector coincides with the
smallest right singular vector at |cos| = 1.0000000000, and ‖Qc‖/‖Q‖ = 2.65e-17.

Regressing each column on the others gives residuals of 3.4e-14, 7.6e-15 and
5.1e-15 for ln θ, ln K_m and the load bias — and **3.685e-03 for ln k_e**, eleven
orders above machine zero. So there is one exact identity and no second one. What
looks like a second small direction (s = 2.193 against 895.6) is k_e being weakly
determined, not degenerate.

The third row is the theorem in one line: **delete the load-bias column and the
deficiency vanishes.** The load bias is the channel; knowing it closes the fibre.

Closure by load variation is immediate — half a per cent is enough:

| day-to-day load CV | s_min | rank |
|---|---|---|
| 0% (as coded) | 9.85e-14 | 3 |
| **0.5%** | **9.08e-02** | **4** |
| 2% | 3.26e-01 | 4 |
| 10% | 1.14e+00 | 4 |

---

## 2D.10 — COMPLEMENTS, not substitutes. And the asymmetry is the recommendation.

Gate sd(ln θ) < 0.20, two-daily × 2 years, Δ carried, 5% load prior:

| | k_e unknown | k_e to 10% | k_e exact |
|---|---|---|---|
| Δ unknown (scale flat) | 0.7364 | 0.5375 | 0.5374 |
| Δ to 20% | 0.7328 | **0.1753** | **0.1744** |
| Δ to 16% | 0.7328 | **0.1440** | **0.1429** |
| Δ exact | 0.7326 | **0.0312** | **0.0251** |

- **Δ alone, k_e genuinely unknown: never reaches the gate** at any dose precision.
- **k_e alone, dose scale unknown: never reaches the gate**, even with k_e exact (0.5374).
- Both are required. They are complements.

But they are required to wildly different standards:

| prior on ln k_e | dose scale the gate then needs |
|---|---|
| flat | unreachable |
| 80% (the module's own) | **15.9%** |
| 50% | 20.7% |
| 30% | 22.5% |
| 10% | 23.3% |
| exact | 23.4% |

**k_e is a threshold; Δ is a dial.** k_e has only to be pinned inside a factor of
about three (111% in logs) — any defensible literature prior does that, and the
module's own 0.8 already does. Tightening it further from 80% to exact buys the
dose requirement only 15.9% → 23.4%, a factor of 1.5. The dose scale is where the
money goes.

And k_e is worth nothing without the therapy channel. On the pre-2C footing, with
Δ dropped, even **exact** k_e gives sd(ln θ) = 0.4892 and fails. So the structure
of the recommendation is: the therapy channel is **necessary**; k_e is a
**permissive side-condition** the literature already satisfies.

---

## 2D.3 — THE CHECK FAILS. The 12.8% is not test-retest, and it is not the mean.

Both halves of the premise are wrong, and I take the primary source over the
repo's register.

### It is between-visit, over months to years

> *"To evaluate inter- and intraindividual variability, seven healthy volunteers
> repeated the study three times, with 1- to 3-week intervals between repeat testing"*

Patients: 15 *mut* subjects, 2–3 repeats each, at **separate NIH visits spanning
2 months to 4 years**. Not the same subject in the same state, and not the same
day. There is **no same-day paired repeatability figure anywhere in the paper.**

The paper names the causes, and they are biological, not analytical:

> *"The higher variations occurred in subjects with changing renal status between
> testing (eGFR decreased from 57 to 43 mL/kg/1.73 m², subject 15), variable VCO2
> measurements (146.0 to 159.2 ml/min, subject 17), or possibly related to
> incomplete dose administration or other factors (subject 2)."*

### And 12.8 is the standard deviation, not the mean

Results section, patients:

> *"At 120 minutes, the CV for 1-13C-propionate dose metabolized in 15 mut
> patients, who had repeated breath testing 2–3 times each, ranged from 0.39% to
> 50.84% with an average of **14.02 ± 12.86%**"*

The abstract reports *"mean coefficient of variation 6.9% and 12.8%,
respectively."* For controls that is right — results give 6.94 ± 5.3%, and 6.9 is
the mean. For patients the abstract has taken **12.86, the SD**, and reported it
as the mean. The mean is **14.02%**.

I checked the published Correction (`s41436-021-01188-9`): it restores a missing
graphical abstract and corrects no numbers. So the discrepancy stands in the
literature.

`paper/claims.csv` C48 quotes the abstract and inherits the error. This is the
third instance of the A3 pattern, and it is on the paper's own headline citation.

### What this does to the derivation

The 2C number — two paired breath tests — **does not survive.** It assumed 12.8%
was independent same-day measurement error inflating by √2 on a difference. It is
neither independent of state nor same-day, so it does not average that way, and
it is the wrong number besides.

What can honestly be said:

- The requirement is unchanged: **15.9%** on Δ (two-daily) or **13.2%** (monthly).
- The only published repeatability is between-visit: **mean 14.02%, SD 12.86%,
  range 0.39–50.84%**, in MMA not PA.
- Between-visit variation is an **upper bound** on same-day measurement error — it
  contains renal drift and VCO2 change that a same-day pair removes. So the true
  paired precision is better than 14%, possibly much better, and plausibly already
  inside 15.9%.
- But **nothing in the literature pins it**, and the range reaching 50.84% means a
  mean cannot be applied per patient.

**So the specification cannot state a number of paired tests. It has to ask for
the measurement that would fix it:** a same-day paired repeatability arm — one
patient cohort, two breath tests bracketing a single dose, on one day — which is
itself small, cheap, and the first thing the sub-study must do. That is a stronger
deliverable than a fabricated count, but it is a different deliverable from the one
2D.9 specifies, so it goes back to you before it is written.

Sources: Manoli et al., *Genet Med* 23:1522–1533 (2021),
[PMC8354855](https://pmc.ncbi.nlm.nih.gov/articles/PMC8354855/); Correction,
[PMC8553604](https://pmc.ncbi.nlm.nih.gov/articles/PMC8553604/).
