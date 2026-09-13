# Phase 2E.1 — what binds the decision. Reported before any prose.

**Gate verdict: the specification is aimed at the wrong quantity — but not for the reason
the hypothesis proposed, and the real reason is worse.**

Neither (a) nor (b) dominates. They are strongly anti-correlated, which is the margin
theorem working. What binds is that **the boundary is a convention, and the paper's chosen
convention sits within a quarter of a per cent of its own sign-flip point.**

Reproduced by `audit/e1_decomposition.py`.

---

## First, a correction to the premise

The hypothesis reads the 511–3285 band (factor 6.4) as a band on θ_crit. It is not. That is
the band on the **margin** — separatrix minus attractor, in metabolite units. Over the same
M_i sweep, `theta_crit_intrinsic` moves 125.5 → 171.4, a factor of **1.37**; over K_c,
119.5 → 186.5, a factor of **1.56**; over ρ₀, 1.02.

So "θ_crit is known only to an order of magnitude" is not what 2B.3 reports. θ_crit is
known to about a factor of 1.5. The hypothesis is tested below at that magnitude, and
separately at the magnitude that turns out to matter.

---

## (a) against (b): anti-correlated, and the full error is smaller than either

One fitted ensemble, μ̂ recomputed four ways. Monthly design, 8 seeds × 400 replicates.

| amplitude | μ_true | full | (a) θ̂ only | (b) threshold only |
|---|---|---|---|---|
| 2.98× (safe side) | +334 | **18.8% ± 1.5** | 47.9% ± 2.5 | 45.9% ± 2.0 |
| 3.58× (unsafe side) | −344 | **40.7% ± 2.4** | 40.0% ± 2.2 | 48.2% ± 2.6 |

On the safe side the full error is **less than half** either component alone. Estimating
θ̂ and the threshold *together* is far better than estimating either with the other held at
truth, because the same K̂_m enters both and the errors cancel. **That is the margin theorem
of eq:mufibre, measured.** Neither (a) nor (b) is the answer to "what binds".

## (c) the distance from the boundary — this is what binds inside the model

| illness | μ_true | \|μ\|/θ | full | (a) | (b) |
|---|---|---|---|---|---|
| 2.00× | +1441 | 0.480 | **0.0%** | 44.3% | 0.6% |
| 2.60× | +763 | 0.254 | 4.1% | 46.1% | 22.5% |
| 2.98× | +334 | 0.111 | 18.1% | 47.1% | 46.4% |
| 3.28× | −5 | 0.002 | **60.5%** | 51.2% | 48.3% |
| 3.58× | −344 | 0.115 | 40.4% | 41.0% | 47.2% |
| 4.07× | −897 | 0.299 | 17.9% | 6.2% | 46.2% |
| 5.00× | −1948 | 0.649 | 0.8% | 0.0% | 43.3% |
| 6.50× | −3642 | 1.214 | 0.0% | 0.0% | 23.1% |

**The full error goes to zero at both ends.** Beyond |μ|/θ ≈ 0.5 the decision is essentially
always right; at the boundary it is a coin toss. So within the model the design is not the
limiting factor — the patient's own distance from the boundary is, exactly as §7 says.

Note the (b) column: **flat at 43–48% across the whole sweep.** If you knew θ exactly and
estimated only the threshold, you would be wrong about half the time *at every distance*.
That is the clue to what follows.

---

## (b′) The boundary is a convention, and the paper sits on the flip

The threshold the decision actually uses is
`θ_crit_local = (s_peak − k_e·m_tol)(K_m + m_tol)/m_tol`. It depends on two quantities
**nobody measures and nobody estimates**: the illness amplitude `s_peak` and the tolerable
metabolite level `m_tol`. Both are chosen. Across defensible choices the threshold runs
**263 to 6332, a factor of 24**, against θ = 3000.

Sign of μ for **the same patient with perfect data**, at every cell:

| m_tol | 2.50× | 2.98× | 3.28× | 3.58× | 4.00× | 5.00× |
|---|---|---|---|---|---|---|
| 10 | − | − | − | − | − | − |
| 20 | − | − | − | − | − | − |
| 30 | **+** | − | − | − | − | − |
| **40** | + | + | **−5** | − | − | − |
| 60 | + | + | + | + | + | − |
| 80 | + | + | + | + | + | + |
| 120 | + | + | + | + | + | + |
| 160 | + | + | + | + | + | + |

**22 of 48 cells call this patient unsafe. The other 26 call the same patient safe.**

And the m_tol at which the sign flips, per amplitude:

| s_peak | flips at m_tol | in units of m̄ |
|---|---|---|
| 2.50× | 28.5 | 2.8 |
| 2.98× | 35.4 | 3.5 |
| **3.28×** | **40.1** | **4.0** |
| 3.58× | 45.0 | 4.5 |
| 4.00× | 52.3 | 5.2 |
| 5.00× | 71.9 | 7.2 |

**The paper uses m_tol = 40, and at the paper's own decision boundary of 3.28× the sign
flips at m_tol = 40.1.** The chosen constant sits 0.25% from the flip.

That is not a coincidence — the boundary was *defined* through that constant, so the two
move together by construction. But the consequence stands regardless: **the decision
boundary at "3.28× baseline" is a restatement of the choice m_tol = 40, and `tab_params`
marks m_tol ASSUMED with no measurement behind it.** The entire decision characteristic of
§7 is computed at a point where the answer is set by an assumption, not by the patient.

---

## What this does to the gate

The gate offered two branches. Neither is right.

- **Not (a).** Improving θ̂ does not help near the boundary and is not needed away from it.
- **Not (b) as posed.** θ_crit is known to a factor of 1.5, not an order of magnitude, and
  estimating it jointly with θ is *better* than knowing it, not worse.
- **(b′).** The boundary is not imprecisely located. It is **conventional**: it is defined by
  two chosen constants, and across their defensible range the sign of μ flips for the same
  patient with perfect data.

So the specification is indeed aimed at the wrong quantity, and rewriting it to aim at a
better-estimated θ_crit would not help, because θ_crit is not a thing the data are
uncertain about. **It is a thing the model asserts.**

---

## 2E.2 — what would locate the boundary?

2E.2 is conditional on (b) dominating. (b′) dominates instead, so the question becomes:
what observable would pin `m_tol` and `s_peak`?

**There is a direct answer in the literature this paper already cites, and it is
unfavourable.** The natural history study reports a *continuous* dose–response: a 50% lower
plasma 3-HP is associated with a decompensation rate ratio of **0.7 (95% CI 0.6–0.9,
p < 0.01)**, with C3 at 0.6 (0.4–0.9) and 2-MC at 1.0 (0.7–1.3). That is a smooth rate
ratio per halving, over the observed range. **It reports no threshold, and a smooth rate
ratio is evidence against there being one.**

So the candidate route — calibrate `m_tol` against observed event rates — does not locate a
threshold. It suggests the threshold does not exist as a feature of the data, and that
`μ = θ − θ_crit` imposes a dichotomy on a relation the evidence says is continuous.

**Two honest options follow, and they are for you rather than for me:**

1. **Keep the threshold and state it as a stipulation.** The paper then reports: given a
   stipulated tolerable level, here is the sign and here is how it moves with the
   stipulation. The sensitivity table above becomes a required part of every claim, and the
   headline decision numbers must be quoted as conditional on m_tol.
2. **Replace the estimand.** Drop the sign of μ and report the quantity the data actually
   support: θ̂ with its interval, mapped onto an expected event rate through the natural
   history study's own rate ratio. The estimand becomes a rate, not a sign, and the
   identification result then carries directly into something clinically meaningful without
   passing through an invented boundary.

Option 2 is the one the evidence points to, and it is a larger rewrite than 2D was. It
would also make the identification result *more* useful rather than less, because a rate
ratio uses the precision on θ̂ that the therapy channel buys, where a sign discards it.

**Nothing further is written until you choose.** Per 2E.4, no prose, no figures, no length
pass.
