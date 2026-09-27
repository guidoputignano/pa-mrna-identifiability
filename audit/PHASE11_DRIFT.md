# Phase 11, item 3. The constant-theta assumption, priced against record length

*Proposals only. Nothing under `paper/` is changed by this file.*

Scripts: `audit/u7_drift.py`. Independent re-check: `audit/u8_recheck_drift.py`.

---

## 0. The question, and why it is not idle

`eq:input` fixes `theta_base` over the whole record. `04_local.tex` warrants that on the
kinetics: most pathogenic variants reduce the *quantity* of functional enzyme rather than
the turnover number, so `theta` is genotype-determined and constant **over the afebrile
range**. Tier 1 needs 118 monthly draws, which is **9.82 years**. "Constant over the
afebrile range" and "constant over a decade" are different assumptions, and growth, disease
progression and any change in endogenous expression all act on the second and none on the
first.

Let `ln theta(t) = ln theta_ref + r (t - t_ref)`. Then variance falls with the number of
visits and bias grows with the span those visits cover, a monthly schedule cannot buy one
without the other, so there is an interior optimum and a floor. The question is whether that
floor sits above or below the 0.20 gate.

**Footing** as `audit/PHASE11_SAMPLING.md` section 0. The estimand is capacity at the
**last visit**, because that is what a dosing interval needs; the midpoint convention is
reported alongside to show the finding is about which `theta` is wanted and not about a
badly chosen origin.

---

## 1. The bias has to be computed exactly, and the obvious way is wrong

The first-order projection of the drift through the estimator,

```
E[beta_hat] - beta_ref = F^{-1} (1/sigma^2) SUM_i q_i q_{i,theta} delta_i
```

is exact **as a projection** and is exactly the weighted least squares it claims to be. It
is also a projection of a perturbation assumed small, and `r` times the record's span is
not small on the records at issue:

| n | r /yr | span yr | `r x span` | linearised bias | exact bias | gap |
|---|---|---|---|---|---|---|
| 13 | 0.01 | 1.00 | 0.010 | -0.00196 | -0.00211 | 0.0002 |
| 13 | 0.10 | 1.00 | 0.100 | -0.00726 | -0.02258 | 0.0153 |
| 25 | 0.05 | 2.00 | 0.100 | 0.00546 | -0.02980 | 0.0353 |
| 60 | 0.05 | 4.91 | 0.246 | 0.25106 | -0.10900 | 0.3601 |
| **118** | **0.05** | **9.74** | **0.487** | **0.59635** | **-0.25166** | **0.8480** |
| 240 | 0.10 | 19.89 | 1.989 | 0.14722 | -1.54787 | 1.6951 |

The linearised column is not merely inaccurate: it has the **wrong sign** from `n = 60`
onward, and it is not monotone in the record length. The exact column is monotone,
correctly signed (a constant fit lags a rising capacity) and tracks `-r x span / 2`, which
is where an information-weighted mean time puts it.

**Every bias below is the exact column**, from a nonlinear weighted least squares on
noiseless drifting data. One bug in that computation is recorded rather than tidied away:
comparing the model's flux at the model's **own root** against the data's flux at the
data's own root makes the residual identically the load, so it loses every dependence on
`theta` and the fit returns its start point. The residual must evaluate the model clearance
**at the observed level**.

---

## 2. The interior optimum exists, and at nominal capacity it is harmless

Monthly draws, capacity wanted at the last visit, medians over 24 intake records,
`rmse = hypot(sd, bias)`:

| n | years | r=0 | r=0.02 | r=0.05 | r=0.10 | r=0.20 |
|---|---|---|---|---|---|---|
| 13 | 1.08 | 0.1842 | 0.1846 | 0.1852 | 0.1861 | **0.1902** |
| 25 | 2.08 | 0.1696 | 0.1712 | **0.1733** | **0.1782** | 0.2011 |
| 50 | 4.16 | 0.1604 | 0.1674 | 0.1838 | 0.2388 | 0.4393 |
| 100 | 8.32 | 0.1553 | **0.1566** | 0.2064 | 0.4382 | 1.0493 |
| 118 | 9.82 | 0.1543 | 0.1597 | 0.2557 | 0.5762 | 1.3717 |
| 240 | 19.98 | **0.1518** | 0.2013 | 0.6441 | 1.6637 | 11.1600 |

The optimum moves in from 240 draws at `r = 0` to 13 at `r = 0.20`, exactly as the
arithmetic predicts.

Swept over the drift rate, at nominal capacity with the increment known exactly:

| r /yr | halving time (yr) | n* | years* | floor | vs gate |
|---|---|---|---|---|---|
| 0.005 | 138.6 | 240 | 19.98 | 0.1444 | PASSES |
| 0.020 | 34.7 | 135 | 11.24 | 0.1534 | PASSES |
| 0.050 | 13.9 | 30 | 2.50 | 0.1687 | PASSES |
| 0.100 | 6.9 | 30 | 2.50 | 0.1687 | PASSES |
| 0.200 | 3.5 | 12 | 1.00 | 0.1828 | PASSES |
| 0.500 | 1.4 | 12 | 1.00 | 0.1709 | PASSES |
| 0.690 | 1.0 | 7 | 0.58 | 0.2030 | FAILS |

**The floor crosses the gate only at `r = 0.67` per year**, a halving time of one year, 95%
a year. On this arm the answer to the question as posed is that the floor sits comfortably
**below** the target for any drift rate one could take seriously.

The reason is not that drift is harmless. It is that **this arm does not need a long
record**: the gate is met at 13 to 25 draws, one to two years, and the drift has no span to
act over. The optimiser simply shortens the record and pays a little variance.

---

## 3. On the arm that actually needs the decade, the answer inverts

Tier 1 is the arm that needs 118 draws, and it needs them because the increment is
characterised only to the published 14.02%. Carrying that throughout:

| r /yr | n* | years* | sd at n* | bias at n* | floor | vs gate |
|---|---|---|---|---|---|---|
| 0.0000 | 240 | 19.98 | 0.1980 | 0.0000 | 0.1980 | PASSES |
| 0.0020 | 240 | 19.98 | 0.1975 | -0.0198 | 0.1985 | PASSES |
| 0.0050 | 240 | 19.98 | 0.1890 | -0.0527 | 0.1962 | PASSES |
| 0.0075 | 222 | 18.48 | 0.1830 | -0.0747 | 0.1976 | PASSES |
| **0.0100** | **181** | **15.06** | 0.1843 | -0.0811 | **0.2013** | **FAILS** |
| 0.0150 | 105 | 8.74 | 0.1952 | -0.0624 | 0.2049 | FAILS |
| 0.0200 | 100 | 8.32 | 0.1917 | -0.0803 | 0.2083 | FAILS |
| 0.0500 | 30 | 2.50 | 0.2101 | -0.0313 | 0.2124 | FAILS |

**The floor crosses the gate at about `r = 0.009` per year: 0.92% a year, a halving time of
76 years.**

And at the requirement's own 118 draws, the arithmetic is even tighter, because 118 was
*defined* as the count that lands on the gate:

| r /yr | sd at n=118 | bias at n=118 | rmse | verdict |
|---|---|---|---|---|
| 0.000 | 0.2000 | 0.0000 | **0.2000** | PASSES (exactly) |
| 0.005 | 0.2006 | -0.0241 | 0.2020 | FAILS |
| 0.010 | 0.1976 | -0.0501 | 0.2037 | FAILS |
| 0.020 | 0.1857 | -0.1063 | 0.2138 | FAILS |
| 0.050 | 0.1590 | -0.3118 | 0.3502 | FAILS |

**The Tier 1 requirement has no margin against drift at all**, because it sits on the gate
by construction. The variance floor the increment's prior imposes is 0.174 at a 14% CV
(`tab:incfloor`), which leaves `sqrt(0.20^2 - 0.174^2) = 0.098` of total budget for bias,
and 9.82 years of any drift exhausts it.

For the midpoint estimand, which is the most favourable choice, the bias at `n = 118` falls
from -0.0501 to a much smaller number at `r = 0.01`, and that is the honest statement of
what the assumption costs when the quantity wanted is the record's **average** capacity
rather than the patient's capacity **now**. A dosing interval needs the second.

---

## 4. The fix is to fit the drift, and it is nearly free

A sixth column, `d/dr`, with a flat prior: the record is asked to find the drift itself.

**The honest price is what it costs when there is no drift to find**, since that is the
premium paid for the insurance. At `r = 0`:

| n | years | constant | drift column | premium |
|---|---|---|---|---|
| 13 | 1.08 | 0.1842 | 0.1859 | **1.009** |
| 25 | 2.08 | 0.1696 | 0.1714 | **1.010** |
| 50 | 4.16 | 0.1604 | 0.1630 | 1.016 |
| 118 | 9.82 | 0.1543 | 0.1553 | 1.007 |
| 240 | 19.98 | 0.1518 | 0.1523 | 1.003 |

**Between 0.3% and 1.6%.** And on the Tier 1 arm it restores the gate at every rate tested:

| r /yr | constant: floor | drift fitted, n=240 | vs gate |
|---|---|---|---|
| 0.005 | 0.1962 | 0.1985 | PASSES |
| 0.010 | 0.2013 | **0.1984** | **PASSES** |
| 0.020 | 0.2083 | **0.1965** | **PASSES** |
| 0.050 | 0.2124 | **0.1636** | **PASSES** |

At larger drift rates it is dramatically better than the constant fit (0.0426 against
1.3717 at `r = 0.20`, `n = 118`), but that comparison flatters it: a drift is itself
excitation in `theta`, so the arm that models it gains information the constant arm cannot
use. The premium table above is the number to quote, not the ratio.

---

## 5. Independent re-check

`u8` adds noise, writes the likelihood in **level** space rather than flux space, fits every
replicate by nonlinear least squares over 240 replicates and measures the mean offset and
robust spread directly. Neither the bias nor the variance comes from anything `u7` computed.

| n | r /yr | kappa | u7 bias | sim bias | u7 rmse | sim rmse | ratio | collapse (Wilson 95%) |
|---|---|---|---|---|---|---|---|---|
| 25 | 0.000 | - | 0.0000 | 0.0174 | 0.1696 | 0.1879 | 1.108 | 0/240 [0.000,0.016] |
| 118 | 0.000 | - | 0.0000 | 0.0033 | 0.1543 | 0.1568 | 1.016 | 0/240 [0.000,0.016] |
| 118 | 0.010 | - | -0.0458 | -0.0303 | 0.1569 | 0.1506 | 0.960 | 0/240 [0.000,0.016] |
| 118 | 0.050 | - | -0.2422 | -0.2375 | 0.2557 | 0.3193 | 1.249 | 3/240 [0.004,0.036] |
| 118 | 0.000 | 0.140 | 0.0000 | 0.0350 | **0.2000** | **0.2099** | 1.050 | 0/240 [0.000,0.016] |
| 118 | 0.010 | 0.140 | -0.0501 | -0.0178 | **0.2037** | **0.2109** | 1.035 | 0/240 [0.000,0.016] |
| 118 | 0.020 | 0.140 | -0.1063 | -0.0647 | **0.2138** | **0.2209** | 1.033 | 0/240 [0.000,0.016] |
| 118 | 0.050 | 0.140 | -0.3118 | -0.1925 | 0.3502 | 0.3188 | 0.910 | 4/240 [0.006,0.042] |

rmse ratios span 0.804 to 1.249. **Every Tier 1 row fails the gate under both methods**, so
the verdict of section 3 is unchanged by the re-check.

The bias comparison is the load-bearing one, since `u7`'s bias is a noiseless flux-space fit
and the simulated one is the mean of a noisy level-space fit. They agree in sign and order
everywhere and closely at the larger rates. At the smaller rates the simulated bias is
smaller than `u7`'s (-0.0303 against -0.0458 at `r = 0.01`), which would move the crossing
of section 3 up somewhat. **Read the crossing as of order one per cent a year, not as a
sharp threshold at 0.0091.**

---

## 6. What this does and does not establish

**It does not establish that capacity drifts.** Nothing here measures `r`, and the paper's
own kinetic argument for constancy is untouched. What it establishes is the **conditional**:
*if* capacity drifts at a rate of order one per cent a year or more, the Tier 1 requirement
as written does not deliver its stated precision, and no schedule of monthly draws recovers
it.

Three things follow, in order of how cheaply they can be had.

1. **The drift column is nearly free insurance** (0.3% to 1.6% when there is no drift) and
   removes the exposure entirely. If a decade-long record is going to be fitted, it should
   be fitted with a drift term whether or not anyone believes in the drift.
2. **A record's span is an assumption-carrying quantity and the paper does not treat it as
   one.** `04_local.tex` warrants constancy over the afebrile range and `08_specification`
   then spends a decade of it. Those two are in different sections and nothing currently
   connects them.
3. **`r` is measurable from the data the specification already asks for**, since fitting
   the drift column returns an estimate of it. A study that runs Tier 1 at all gets the
   drift rate for the price of one parameter.

**Unsettled.** Everything about the plausible magnitude of `r`. Growth in a paediatric
cohort, hepatic disease progression and the therapy's own effect on endogenous expression
all act in this window and none is quantified here or anywhere in the repository. The sweep
is a sweep precisely because the container has no egress and no source could be checked.
