# Phase 3.0.5 and 3.0.6 — reported before any prose.

**3.0.5's gate fires, and against Option 2.** At the patient's own operating point the
model-based rate is not merely no better than reading the rate off the marker — it is
**13.5× worse**. The model's only exclusive capability is counterfactual prediction.

Reproduced by `audit/e3_baseline_and_scope.py`.

---

## 3.0.5 (i) — at the operating point, the direct read-off wins outright

Monthly design, n = 25, assay CV 5%.

| route | sd(ln m) | × β | 95% width |
|---|---|---|---|
| **direct read-off, no model, no θ** | 0.0100 | **0.0051** | **1.02×** |
| model-based, via θ̂ | 0.1354 | **0.0697** | **1.31×** |

The model reconstructs a level that was **measured directly**, and pays parameter
uncertainty to do it. Averaging 25 monthly assays gives the level to 1%; the fitted
parameters give it to 13.5%.

This is not a close call and it should not be argued away. The natural history study
already publishes the association; applying it to a measured level requires nothing this
paper contributes.

## 3.0.5 (ii) — the model's one exclusive capability is counterfactual

| counterfactual | model sd(ln rate) | direct |
|---|---|---|
| load 1.5× baseline | 0.0711 | **no answer possible** |
| load 3.0× baseline | 0.0906 | **no answer possible** |

"What would this patient's rate be if the dose were halved, or if their load rose by
half?" has no observation behind it, so the direct read-off is silent. **This is the only
column where the model has something the data do not**, and whether it is worth a paper
depends on whether that counterfactual is the clinical question. It is a different
question from the one the paper opened with.

## 3.0.5 (iii) — the map's own uncertainty dominates both routes

RR 0.70 (95% CI 0.60–0.90) per halving gives β = −log₂(RR) = **0.5146, 95% CI 0.152 to
0.737 — the exponent is known only to a factor of 4.8.** Neither route can reduce it:

| m / m_ref | contribution from the map alone |
|---|---|
| 1.25 | 0.0333 |
| 1.50 | 0.0605 |
| **2.00** | **0.1034** |
| 3.00 | 0.1639 |
| 5.00 | 0.2402 |

At a doubling of the marker the map alone contributes 0.103 — **larger than the model's
entire parameter uncertainty of 0.070.** Any argument that the model's precision matters
is an argument about the second-largest term.

## 3.0.5 (iv) — and the map has no anchor, which blocks both routes equally

`rate(m) = rate_ref · (m/m_ref)^β` needs both constants. The study reports a cohort mean
annualised MDE frequency of 1.0 (SD 1.6), so `rate_ref` exists. **It reports no
concentration range and no reference 3-HP value, so `m_ref` does not.** Without it neither
route yields a rate — only a ratio between two levels. The model adds nothing to that.

---

## 3.0.6 — where the estimand fails, in observable units

The marker rises towards a ceiling as capacity falls: as θ → 0, m → L₀/k_e = **6.45× the
managed level**. The elasticity dies with it.

| θ (% of nominal) | marker, × managed | \|elasticity\| |
|---|---|---|
| 100% | 1.00 | 0.915 |
| 50% | 1.82 | 0.807 |
| 20% | 3.41 | 0.537 |
| 10% | 4.57 | 0.321 |
| 5% | 5.40 | 0.173 |
| 2% | 6.00 | 0.072 |
| 1% | 6.22 | 0.036 |

| elasticity falls below | at θ | marker reads |
|---|---|---|
| 0.50 | 17.9% of nominal | **3.6× managed** |
| 0.25 | 7.5% | **5.0× managed** |
| 0.10 | 2.8% | **5.8× managed** |
| 0.05 | 1.4% | **6.1× managed** |

**The clinician-facing statement, and it belongs wherever the estimand is claimed:** once a
patient's managed marker sits above about three to four times the reference level, a
further fall in residual capacity barely moves it. Two patients differing tenfold in
capacity are then indistinguishable on the marker.

**But the failure is asymmetric, and that asymmetry is the one point favourable to
Option 2.** What dies is the estimate of *capacity*. The *rate* prediction does not,
because the rate is read off the marker and the marker is what is measured. So the scope
limitation is on reading θ out of the marker for severe patients — which is an argument
for the rate estimand and against the capacity one, in exactly the population the paper
opened by caring about.

---

## The rate map's provenance, to be recorded in the text as-is

Verified against Results, not the abstract:

- **Unadjusted GAMM.** Three terms: 3-HP concentration, a smooth in days since the previous
  MDE, a per-participant random-effect smooth. No adjustment for age, disease onset, or
  treatment.
- **Unresolved whether 3-HP entered as baseline-only or time-varying.** Repeated measures
  exist (6-monthly to month 36, plus during MDEs), and a figure caption says "only baseline
  biomarkers are presented", but the source does not say which governs the GAMM.
- **No concentration range**, no explicit log-linear specification, so the per-halving
  coefficient cannot be extrapolated beyond a range that is itself unstated.
- **No N** for the regression models.
- **The authors' own limitation:** the biomarkers "may not represent their overall severity",
  and "further investigation is required to determine whether they are sufficient to
  reflect the risk of MDE in future clinical trials."

Direction is as expected for a toxic metabolite, so this is not a paradoxical protective
RR. But it is unadjusted, and severity and treatment intensity both act in the same
direction, so the association is an upper bound on any causal reading.

**This is not an established dose–response and must not be presented as one.**

---

## Recommendation

**Return to Option 1**, on the evidence above, and narrow the paper's claim.

Option 2's warrant was that a rate uses the precision on θ̂ that a sign discards. 3.0.5
shows it does not use it: the rate is better estimated **without** θ̂ at the operating
point, and the map's own exponent dominates anyway. The reframing argument against the
threshold survives intact — the evidence is a smooth dose-response with no threshold — but
that argument does not require the paper to adopt the rate as its estimand.

What the paper can defensibly claim is narrower than either option promised:

1. **The identification result**, which stands on its own: the rank-one identity and its
   closure condition, the dose train supplying 42× the marker band, and sampling on dosing
   days costing a factor of 22 while half a day off recovers it.
2. **That identification is not sufficient for the decision**, demonstrated rather than
   asserted — with the 24× threshold sensitivity and the sign-flip contour as the evidence,
   and `m_tol` attached to every decision number.
3. **That the threshold framing is itself questionable**, warranted by the cited smooth
   dose-response, stated as a finding rather than adopted as a replacement.

That is a smaller paper than Phase 3 envisaged and it is the one the evidence supports.
The decision is yours; nothing is written until you take it.
