# Phase 4, the stress test: eight attempts to refute P5 before a word was retracted

`audit/PHASE4_PI.md` P5 claims the constant-load degeneracy is an exact fibre and that the
paper's "recovered through curvature" and `sqrt(noise)` sentences are false. That is a
retraction of published text, and this repository's review log already records two
retractions issued without checking (`R5`, `R7`) and one number that was right about its own
computation and wrong about what it computed (`R8`). So the finding was attacked before it
was used.

Eight independent lenses, each prompted to **refute** P5 and to default to "refuted" on
anything material, each then re-run from scratch by a second agent told to check the first
rather than to check P5. Thirteen of twenty agents completed on the first pass; the rest
were re-run.

**Outcome: sub-claims (1) to (4) survived all eight lenses. Sub-claim (5) was refuted, and
it was the audit's, not the paper's.** Four further errors in the audit were found, one of
them a bug in `audit/f5_fibre.py`. The corrections are folded into `PHASE4_PI.md` as C1-C8.

---

## The lenses, and what each returned

| lens | verdict | what it did |
|---|---|---|
| re-derivation | PARTIAL, (5) fails | eliminated `m` by resultant and Groebner rather than solving for `Delta` |
| exact arithmetic | survives | sympy identity, exact rationals, and a precision-tracking test at 15 to 120 digits |
| basis exhaustion | survives | attacked the step from "independent as functions" to "exhausts a finite record" |
| ODE transient | survives | forward-sensitivity profile Fisher information about `b` from the transient |
| noise model | survives | multiplicative-on-`m`, heteroscedastic, and full errors-in-variables likelihoods |
| reconstruct `36.9` | survives | tried to vindicate the paper: 88 designs x 7 sigma x 3 scales, plus every profile variant |
| already known | survives | tried to show the paper already says it, so the fix would be wording |
| scope overreach | PARTIAL, (4) in part | granted the fibre and attacked the consequences instead |

### What could not be broken

**The fibre is an identity, not a numerical coincidence.**
`R_fibre(m) = R_true(m) * [1 + K_0 b/(K_0 L + (L+b) m)]`, so the fibre's residual vanishes
wherever the truth's does, for every `b`, `Delta` and `m`. Any true profile leaving `theta`
and `K_m` free is therefore exactly flat, and no profile construction can produce a finite
interval.

**There is no fourth invariant.** By the opposite elimination, the resultant in `m` is
quadratic in the known `Delta`, and requiring every coefficient to vanish gives sympy the
single solution `{A2: A1, B2: B1, k_e2: k_e}`. Its Groebner basis has seven generators,
among them `(A1-A2)(k_e-k_e2)` and `(B1-B2)(2 A1 k_e - A1 k_e2 - A2 k_e)`. The map
`(theta, K_m, k_e, b) -> (A, B, C)` has Jacobian rank 3 of 4 symbolically, and the analytic
level-map Jacobian has fourth singular value `9.0e-17` to `9.4e-16` against a first of
`~1.99`.

**The `1e-25` is not a cancellation artefact.** Re-run in mpmath at 15 / 30 / 60 / 120
digits, `chi^2` at `b = 30` falls `2.20e-25 -> 1.91e-55 -> 2.14e-115 -> 8.22e-235`, and the
same at `b = 0` and `b = 200`. It tracks working precision across 210 orders of magnitude,
which a manufactured zero cannot do. In exact rationals, at
`K_m = 8875/97`, `theta = 298035/97` (`b = 30`), the residual is the integer `0` at every
increment tested.

**The observable itself is invariant**, which is stronger than the argument the audit made.
`max|m(b) - m(0)| = 0.000000e+00` in floating point at `b = 1, 10, 30, 100, 200, 1000, -100,
-200`, and a marker-scale profile under multiplicative noise moves `chi^2` by at most
`4.3e-4` out to `b = 5000` against `chi^2(0) = 20.73`. So every likelihood whose
noise acts on `m` with parameter-free constants is exactly flat, whatever its shape.

**The `sqrt(noise)` price is refuted quantitatively, not just structurally.** Halving the
noise and reading the 95% true-profile half-width in `b` (`audit/f6_fibre_adversarial.py` D):

| load CV | `sigma` | `sigma/2` | `sigma/4` | `sigma/8` | ratio per halving |
|---|---|---|---|---|---|
| 0% | unbounded | unbounded | unbounded | unbounded | halving buys nothing |
| 0.5% | 2108.5 | 905.8 | 422.4 | 204.3 | 0.430, 0.466, 0.484 |
| 2% | 423.5 | 205.3 | 101.6 | 50.9 | 0.485, 0.495, 0.501 |
| 10% | 80.3 | 40.9 | 21.4 | 11.6 | 0.509, 0.523, 0.544 |

`sqrt(noise)` predicts `0.7071` per halving and a first-order direction predicts `0.500`. The
ratios sit on `0.5`. At a constant load, where the paper's sentence is aimed, halving the
assay CV buys nothing at all.

**The transient does not rescue it.** Re-optimising the other three parameters absorbs
`99.9934%` of the ODE's departure from the locus; `sd(b_hat)` from the transient alone is
`3.90e7`, 2.4 million times the load prior.

**The consequence claim holds.** The load-bias prior is worth 1.4% of `sd(ln theta)` on the
monthly design at load CV 10% and 0.6-3.5% elsewhere. No printed number in the paper changes.

---

## What broke, and it was this audit's

### The refutation: claim (5)

`chi^2/eps^4` is a **second-order** quantity, so it is parametrisation dependent. `prop:rank`
is stated in log parameters (`06_excitation.tex:332`) and the information filter's columns
are derivatives in those coordinates. The audit drew the straight line in natural
`(theta, K_m)` coordinates and then went looking for a design that matched, and found a
plausible wrong answer.

In the paper's own coordinates the two-daily x 2 y design (n = 365) gives
`chi^2/eps^4 = 2.0900e-06` against the printed `2.08e-6`, and a quartic half-width of `36.82`
against the printed `36.9`. The provenance is documented rather than inferred:
`audit/c2b_degeneracy.py:95` runs `np.arange(0.0, 730.0, 2.0)`; `audit/PHASE2C.md:92-94`
attaches the `2.08e-6` to `sd = 0.0251`; `audit/PHASE2A.md:35` shows `0.0251` is the
two-daily row.

So "a design of roughly weekly size", "not the monthly one", "`+/-54.2`" and "a factor of
1.72, not within a fifth" are withdrawn. On the design the number was computed for, the
paper's "the two agree to within a fifth" is right as computed.

What survives is the first half, and it survives better than before: the `+/-36.9` is a ray
and not a profile. The printed pair is internally the quartic extrapolation,
`(3.84/2.08e-6)^(1/4) = 36.86` against a printed `36.9`, and a true profile does not yield a
constant `chi^2/eps^4` to read a half-width off. A weakened design charge also survives,
because the paper points `rem:second_order` at a monthly covariance where the same
construction gives `69.4` against `31.6`.

### Four more

- **A bug.** `audit/f5_fibre.py`'s bisection for `theta(b) = 0` inverted its bracket update
  and printed `-215.1515`, which is where `K_m` reaches 300 and `theta` is still 1785. The
  root is `-272.7273`. Fixed with `brentq`.
- **A single seed.** S4's `k_e` median of `-0.879` is the minimum over eight seeds; the
  median across them is `+4.67` and the sign is a coin toss (negative in 42-51% of
  replicates). The `theta_hat <= 0` rate survives the sweep at 0.22-0.29. S4 now quotes
  spreads, which is this repository's own rule and which the first version broke.
- **Scope, too wide.** "Flat only to first order" is a true statement about `prop:rank`'s
  null vector. Only two clauses are false.
- **Scope, too narrow.** `rem:rankvsprec` was called "correct as it stands" in one place and
  "generous" in another. One of its sentences is false at two of the four rows it prints.

---

## Two things the lenses found that the audit had missed

**The paper's proof drops the bias.** `06_excitation.tex:371` writes
`(theta+Delta) u = L - k_e m` where the balance gives `L + b - k_e m`, so `eq:nullvec` as
printed is the `b = 0` value of a field that depends on `b`. Carrying `b` through gives
`n = (L + b + k_e K_m, -K_m, 0, L + b)`, and integrating that reproduces the fibre with
symbolic difference exactly `0`; integrating the printed frozen-`L` vector gives
`K_m = K_0 exp(-b/L)`, which leaves the fibre at second order. This is why the exact fibre is
not derivable from the paper as written, and it should be fixed whatever is decided about
the rest.

**"Carries no information" overstates, by two channels.** Neither rescues the paper. Under
errors-in-variables with an unrecorded per-visit load fluctuation the lever arm
`d(load)/dm` varies along the fibre even though `m` does not, giving a genuine 95% interval
in `b` of about `[-120, +383]`, four to twelve times looser than the prior, saturating rather
than following `sqrt(noise)`. And a record sampling inside the post-dose settling window does
break the fibre: 225 such draws would reach `+/-24.2`. That design is excluded, but
replacement text should not say the combination is never recoverable.

---

## One thing a lens got wrong, and how it was caught

The re-derivation lens reported that the synchronous 14 d and 28 d designs drop to rank 2, so
the fibre there would be two-dimensional. Its Jacobian was built by central finite differences
at `h = 1e-6`, and its small singular values floored at the differencing noise. Rebuilt
analytically by implicit differentiation and evaluated at 60 digits, those designs are rank 3
with a 24-order gap, exactly as the lens's own resultant argument predicts. The claim was
dropped.

That is the reason each lens was checked by a second agent rather than taken at its word, and
it is the same discipline this document is asking of the audit it is checking.

---

## Reproducing it

The lenses ran in throwaway scripts. The four checks the retraction actually rests on are
re-implemented in `audit/f6_fibre_adversarial.py`, which is committed and which reproduces
every number quoted above under those four headings:

```
python audit/f6_fibre_adversarial.py
```

  A  the resultant and Groebner elimination, and the rank of `(theta, K_m, k_e, b) -> (A, B, C)`
  B  the precision-tracking test at 15, 30, 60 and 120 digits, and the exact-rational residual
  C  the bit-identity of the marker along the fibre, and a marker-scale profile
  D  the half-width against the noise level, which is the `sqrt(noise)` test

The coordinate-system finding that refuted the audit's claim (5) is in
`audit/f5_fibre.py` S3b, which now prints both parametrisations side by side.

---

## Method note

The stress test is not part of the finding and does not make it true. What it does is bound
the ways it could have been false that anyone bothered to try. Three of the eight lenses
returned a correction that changed the text; one returned a correction that was itself wrong;
four returned nothing but confirmation. The single most useful design choice was prompting
each lens to refute rather than to evaluate, because the two lenses that produced the
refutation of claim (5) both began by trying to reproduce the paper's number rather than the
audit's.
