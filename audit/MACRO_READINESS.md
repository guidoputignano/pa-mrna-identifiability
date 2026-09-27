# Macro readiness: which of the 29 can be wired, checked against its generator's current state

*Produced before any wiring, per the author's ordering. Nothing is wired by this file and nothing
under `paper/` is changed.*

**The test is the one that caught `d0c` at R141**: does a generator compute the value **today**,
and does the value agree with what the repository already carries? All six owning scripts were run
for this: `u2`, `u7`, `u9`, `u11`, `u12`, `u13`, **all exit 0**.

---

## READY: 17 macros, generator runs today and the value is confirmed

| macro | value | generator | confirmed against |
|---|---|---|---|
| `ke_frontier_3000_pct` | **80.01** | `u9` | `u9` PART A, `u13` PART B |
| `ke_frontier_600_pct` | **15.67** | `u9` | both |
| `ke_frontier_300_pct` | **6.33** | `u9` | both |
| `ke_frontier_600_factor` | **1.170** | `u13` | `u13` PART B |
| `ke_frontier_300_factor` | **1.065** | `u13` | `u13` PART B |
| `ke_prior_pct` | **80** | register, `KE_PRIOR_LOG` | `params.py` |
| `profile_nominal` | **0.2002** | `u12` | `u12` PART A, recommended design |
| `profile_lo_nominal` | **0.2679** | `u12` | `u12` PART B |
| `profile_hi_nominal` | **0.1327** | `u12` | `u12` PART B, ratio 0.495 |
| `profile_gap_tier1_pct` | **13.5** | `u12` | `u12` PART A, Tier 1 ratio 1.135 |
| `quad_width` | **0.10** | `u12` | `u12` PART B, ratio 1.000 there |
| `drift_cross_pct` | **95.4** | `u7` | `u7` PART C, `r = 0.67`/yr |
| `drift_premium_lo` | **0.3** | `u7` | `u7` PART D, ratio 1.003 |
| `drift_premium_hi` | **1.6** | `u7` | `u7` PART D, ratio 1.016 |
| `margin_target` | **0.17** | author's (g4) decision | not computed, and should not be |
| `rho_pole_q` | **0.800** | `u18` | `u18` PART B |
| `sd_monthly` *(exists)* | **0.1699** | `src/design/run.py` | already generated |

**One thing to state about `profile_nominal` before it is wired.** `0.2002` is the **recommended
25-monthly design's** profile. **Tier 1's profile is `0.2270`**, a different number on a different
design, and `profile_gap_tier1_pct` is Tier 1's ratio while `object_spread_pct` is the recommended
design's. **The wiring must not conflate them**, and a single key called `profile_...` invites
exactly that.

---

## BLOCKED: 2 macros, on a footing conflict rather than a wrong number

| macro | value | the conflict |
|---|---|---|
| `object_spread_pct` | **18.2** or **17.8** | profile against `u12`'s marginal 0.1695, or against the paper's committed `sd_monthly` 0.1699 |
| `sim_spread_nominal` | 0.1777 | sits on the same choice of footing |

`PHASE8_SIMULATIONS.md:623` documents the difference exactly: the quantity is `0.1699` **as
committed, with the `t = 0` draw read post-dose** inside the settling window, `0.1695` **with that
draw read pre-dose**, `0.1721` with it dropped. **Neither is wrong.** R142.

> **This is the case the ordering exists for.** Wire either and `check_artefacts` certifies a
> footing nobody chose. It needs one line from the author, and (g4)'s target of 0.17 stands either
> way.

---

## DEGENERATE: 3 macros where macro-ising the sentence made it worse

| macro | value |
|---|---|
| `rho_above_fail_pct` | 100.0 |
| `rho_above_fail_lo` | **100.0** |
| `rho_above_fail_hi` | **100.0** |

The proposed sentence reads *"at $\num{rho_above_fail_pct}\%$, with a seed-to-seed range of
$\num{rho_above_fail_lo}$ to $\num{rho_above_fail_hi}$"*. The spread is **zero across all eight
seeds**, so that renders as **"a range of 100.0 to 100.0"**, which is absurd.

**The original prose said "with no seed-to-seed variation", and it was better.** I converted it to
macros at R134 to fix the single-seed defect, and in doing so replaced a correct English statement
with a degenerate numeric one. **Recommendation: keep `rho_above_fail_pct`, drop the other two, and
restore the prose clause.** A macro is the right answer for a number that could move and the wrong
one for a quantity whose content is that it does not.

---

## NOT YET EXTRACTED: 7 breath macros

`breath_el_nominal`, `breath_el_disease`, `marker_el_nominal`, `marker_el_disease`,
`breath_k_600`, `breath_k_300`, `breath_k_300_pub`. Owners `u2` and `u11`, both exit 0. Their
values are row-selections from tables indexed by breath CV, so each needs its row named before it
can be promoted, and naming the row is a choice rather than a lookup. **Deferred rather than
guessed.**

---

## What this pass established

- **17 ready, 2 blocked on one author decision, 3 that should not be macros at all, 7 deferred.**
- **No repeat of the `d0c` failure**: every value in the READY table is produced by a script that
  runs today, not asserted in a docstring or a findings file.
- **The blocking item is a footing conflict, which is a class the `d0c` test was not designed to
  find.** It was found by checking the four macros that *already exist*, which was not part of the
  plan and turned out to be the highest-yield step.
