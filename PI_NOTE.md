# Note for the PI

Three months on from your feedback, this is where the work stands. Nothing here is new
computation. It is assembled from the paper, the supplement and the audit records, so every
number in it is one you can find in those documents.

---

## Page 1. What your feedback produced

You made two comments. Each one changed the work, and the second one produced the result the
paper is now about.

### "What is the model"

At the time there were three models, called A, B and C, and the note did not say how they
were related. There are now three objects and they are not three models. They are one model
read at three resolutions, and the relations between them were checked rather than asserted.

The dynamic model is the parent ODE, with a two-state carnitine system behind it. Solving
the fast state for its equilibrium and substituting returns the reduced clearance exactly
and symbolically, with the conjugation rate constant entering only through a ratio. Setting
the marker's derivative to zero gives the equilibrium reduction, whose solutions are the
attractor and the separatrix and whose fold defines the critical capacity. Integrating the
two-state system from below the separatrix settles on the reduced attractor to fourteen
digits, and from above it diverges.

The third object is the local identification model, and it is the one the check did not
pass. It is not a reduction of the parent. The tangent to the parent's clearance curve
carries an extra bracket in its intercept, and at the illustrative constants that bracket is
374 against the enzymatic term's 61. The alternative channel owns six-sevenths of the
intercept that the residual capacity was supposed to be carried by. The local model is
retained, but it is now stated for what it is: a different hypothesis about one term, kept
because it is the most favourable case for identification rather than because it
approximates anything. That matters for the conclusion. The degeneracy the paper is about
appears in the favourable structure and disappears in the awkward one, so it is not an
artefact of a bad approximation.

The letters A, B and C are retired.

### "Time is a constraint"

This produced the paper's central result. The reason is one sentence. A static model has no
coordinate in which the therapeutic increment can vary, so the dose train is invisible until
the model is dynamic.

Once the model is dynamic the increment enters as a known, exogenous input acting selectively
on the parameter of interest, and the fixed biweekly schedule the trial already runs turns
out to be the identification experiment. Everything the paper now recommends is downstream
of that. It is also why the recommendation is a sampling schedule and not a new assay.

### Where your other suggestions landed

The filter became augmented-state estimation on the coordinates (ln theta, ln K_m, L_0,
ln k_e), with no process noise on theta. It earned its place by displaying the degeneracy
rather than concealing it. A random walk on theta would manufacture a finite stationary
covariance in a direction the data may never see, which is precisely the error the paper
exists to warn about. As constructed, what does not collapse stays visibly uncollapsed.

The multiple-model idea became the nested reduction ladder above, and then a mechanical
check. StructuralIdentifiability.jl returns globally identifiable for eight configurations
of the model: time-varying load, constant load, two patients sharing population constants, a
null-genotype anchor, carnitine as a known input, therapy as a known input, and a two-state
carnitine model in which the second state is never observed. The single exception is the
variant in which both clearance channels saturate, and that exception is proved rather than
reported.

Euler stayed in the exposition, where it is clearer. Every computed trajectory uses an
adaptive stiff-capable solver at relative tolerance 1e-8 or tighter, and two tests guard it:
one halves the step and asserts the answer does not move, the other asserts that an
integrated steady state matches the closed-form equilibrium reduction.

---

## Page 2. What the paper now claims

Four claims, in the order the paper states them.

1. The dose train is the identification experiment, and it is already being run. The
therapy is administered on a fixed biweekly schedule and read as treatment. It is also a
known, exogenous, selective input on exactly the parameter in question, and it supplies
forty-two times the excitation the marker's own excursions do. That is what makes theta
identifiable. Monthly outpatient draws across dose cycles reach 17 per cent over two years,
twenty-five draws, no excursion and no admission. Thirteen two-monthly draws reach 19 per
cent. Both are conditional on the non-enzymatic clearance being known within a factor of
three and the delivered increment characterised to about 12 per cent. Nothing is assumed
about the half-saturation constant.

One scheduling choice destroys it, and it is the choice a protocol writer would make for
convenience. Sampling on dosing days, or at any multiple of the dose interval, puts every
draw at the same phase of every cycle and costs a factor of twenty-two in the precision of
the answer. Moving the schedule half a day recovers it.

The corollary is a limit on reach rather than on precision. Nothing else supplies enough
excitation. Spontaneous load variability alone gives 1.77 against a 0.20 target, and no level
of dietary variability reaches the target at all: the precision floors near 0.376 and rises
again. Recorded carnitine steps make it worse, not better, because carnitine acts on the
wrong channel. So the design applies to patients receiving the therapy and to no one else.

2. Identifying theta does not decide the dosing interval. The decision is usually posed
as the sign of a margin against a critical capacity. That threshold is not a measured
quantity. It is defined through two chosen constants, and across their defensible range it
moves by a factor of 24 and the sign flips for the same patient with perfect data.

3. The threshold framing is itself doubtful. The only endpoint evidence available
reports a smooth dose-response and no threshold at all. The paper records that as a finding
and does not adopt a replacement estimand.

4. Capacity stops being recoverable exactly where it matters most. Above roughly three to
four times the managed marker level, two patients differing tenfold in residual capacity are
indistinguishable. That is the severely affected population.

### The qualification that belongs next to the headline

The 17 per cent is one patient, the illustrative parameter set. Across four thousand
simulated patients the achieved precision has median 0.192, quartiles 0.132 and 0.329, and
just over half meet the 0.20 target. The shortfall is systematic rather than random. The
dose is fixed in absolute terms and the capacity it is measured against is not, so patients
with more residual capacity receive less relative excitation and are identified less well.
The design is specified as a cohort-level demonstration with a per-patient interval
reported, never as a guarantee attached to an individual record.

---

## Page 3. What is open, and what is asked of you

### The missing measurement

Everything the therapy channel buys is stated in units of theta, so the ratio of increment
to capacity is only as well known as the increment is. The design needs it characterised to
about 12 per cent, and it has never been measured. What is published is an upper bound and
it is not the quantity needed.

The first tier of the data specification is therefore one study, and it is small. One
cohort, two 1-13C-propionate breath tests bracketing a single dose, on one day. What it
measures is the repeatability of the assay within a day and within a state, which is what a
pre-and-post pair around a dose actually depends on. The paper ends in that specification
rather than in a dosing rule.

### The two questions in your field

The first is what estimator, and what uncertainty statement, is legitimate here. The information
matrix is singular at the truth, and the degenerate combination is flat only to first order,
so it is recovered through curvature rather than through slope. The consequence is a price
rather than a prohibition: the recoverable displacement scales as the square root of the
noise level where a first-order direction would scale linearly. A singular Fisher matrix at
a point does not contradict global injectivity of the map; it sets the price. The paper
prices that and reports intervals from the inverse information matrix. What it does not have
is an argument that this is the right object to report when the local quadratic model is the
thing that fails.

The second is whether the reduction error can be bounded a priori. The 374 against 61 on page 1 was
measured after the fact, by writing down both intercepts and comparing them. A model
reduction that preserves the input-output response to periodic impulsive inputs does not
obviously yield a reduced model whose parameters are identifiable from sparse sampling.
These are distinct objectives and they may conflict. Whether the gap can be bounded before
choosing the reduction, rather than discovered afterwards, is the question the paper would
most like answered.

### The state of the documents

The paper is 26 pages and the supplement is 52. Both are frozen. The supplement carries
every section that was moved out of the main text, with the superseded positions kept and
their resolutions attached rather than deleted. Thirteen retractions and method failures are
recorded in the review log with their causes, four of them found in the final passes. One of
those four had reached the abstract.
