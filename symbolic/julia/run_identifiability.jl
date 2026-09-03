# =====================================================================
#  Structural identifiability of the propionate models
#  Requires: julia >= 1.9,  ] add StructuralIdentifiability
#  Run with:  julia run_identifiability.jl   (expect several minutes)
# =====================================================================

using StructuralIdentifiability
using Logging

banner(t) = (println("\n" * "="^70); println("  ", t); println("="^70))

results = Dict{String,Any}()

function report(name, ode)
    banner(name)
    println("\nMODEL:")
    println(ode)

    println("\n--- LOCAL identifiability ---")
    loc = assess_local_identifiability(ode)
    println(loc)

    println("\n--- GLOBAL identifiability ---")
    glob = assess_identifiability(ode)
    println(glob)

    println("\n--- IDENTIFIABLE FUNCTIONS of the parameters ---")
    println("(if a parameter is not identifiable alone, this shows which")
    println(" combination of parameters is)")
    funcs = find_identifiable_functions(ode)
    for f in funcs
        println("   ", f)
    end

    results[name] = (local_id = loc, global_id = glob, functions = funcs)
    return nothing
end

# ---------------------------------------------------------------------
# A1.  Baseline. Generic time-varying input s(t).
#      EXPECTED: theta, Km, ke all globally identifiable.
#      This is Proposition 1 of the note.
# ---------------------------------------------------------------------
A1 = @ODEmodel(
    m'(t) = s(t) - theta * m(t) / (Km + m(t)) - ke * m(t),
    y(t)  = m(t)
)
report("A1  baseline, time-varying input s(t)", A1)

# ---------------------------------------------------------------------
# A2.  Constant input. s is promoted to an unknown constant parameter
#      with no dynamics, which encodes a patient held at steady state.
#      EXPECTED: individual parameters NOT identifiable; the identifiable
#      function should be the degenerate combination theta/Km + ke
#      (possibly returned in an equivalent algebraic form).
#      This is the formal version of the clinical situation.
# ---------------------------------------------------------------------
A2 = @ODEmodel(
    m'(t)  = s0 - theta * m(t) / (Km + m(t)) - ke * m(t),
    s0'(t) = 0,
    y(t)   = m(t)
)
report("A2  constant input (steady-state patient)", A2)

# ---------------------------------------------------------------------
# A3.  Two patients sharing the population parameters (Km, ke),
#      each with an individual theta_i, both at constant input.
#      EXPECTED: cohort structure ALONE does not resolve the degeneracy.
#      Section 3 of the note predicts deficiency 2 regardless of N.
# ---------------------------------------------------------------------
A3 = @ODEmodel(
    m1'(t) = s1 - th1 * m1(t) / (Km + m1(t)) - ke * m1(t),
    m2'(t) = s2 - th2 * m2(t) / (Km + m2(t)) - ke * m2(t),
    s1'(t) = 0,
    s2'(t) = 0,
    y1(t)  = m1(t),
    y2(t)  = m2(t)
)
report("A3  two patients, shared (Km, ke), constant inputs", A3)

# ---------------------------------------------------------------------
# A4.  As A3 but patient 2 carries a NULL genotype: theta_2 = 0.
#      Its equation reduces to  m2' = s2 - ke*m2, which involves ke
#      alone. EXPECTED: ke becomes identifiable.
#      This is condition (i) of the anchor proposition.
# ---------------------------------------------------------------------
A4 = @ODEmodel(
    m1'(t) = s1 - th1 * m1(t) / (Km + m1(t)) - ke * m1(t),
    m2'(t) = s2 - ke * m2(t),
    s1'(t) = 0,
    s2'(t) = 0,
    y1(t)  = m1(t),
    y2(t)  = m2(t)
)
report("A4  null-genotype anchor (theta_2 = 0)", A4)

# ---------------------------------------------------------------------
# A5.  Carnitine as a known input acting on the ALTERNATIVE channel.
#      ke -> ke + kappa * uc(t).
#      Caveat outside the mathematics: carnitine dose is adjusted in
#      response to patient state, so uc is endogenous in real data.
# ---------------------------------------------------------------------
A5 = @ODEmodel(
    m'(t) = s0 - theta * m(t) / (Km + m(t)) - (ke + kappa * uc(t)) * m(t),
    s0'(t) = 0,
    y(t)  = m(t)
)
report("A5  carnitine input on the alternative channel", A5)

# ---------------------------------------------------------------------
# A6.  Therapy as a known input acting on the ENZYMATIC channel.
#      Capacity becomes theta + d(t), with d(t) the known increment.
#      This is the selective input on the other channel.
# ---------------------------------------------------------------------
A6 = @ODEmodel(
    m'(t) = s0 - (theta + d(t)) * m(t) / (Km + m(t)) - ke * m(t),
    s0'(t) = 0,
    y(t)  = m(t)
)
report("A6  therapy input on the enzymatic channel", A6)

# ---------------------------------------------------------------------
# B.   Both clearance terms saturating.
#      EXPECTED: weaker than A, since the two channels now share a
#      functional form and are distinguished only by Km != Ke.
# ---------------------------------------------------------------------
B = @ODEmodel(
    m'(t) = s(t) - theta * m(t) / (Km + m(t)) - ke * m(t) / (Ke + m(t)),
    y(t)  = m(t)
)
report("B   both channels saturating", B)

# ---------------------------------------------------------------------
# C.   Two-state carnitine depletion model.
#      c(t) is NOT measured, which is the point: only m is observed.
# ---------------------------------------------------------------------
C = @ODEmodel(
    m'(t) = s(t) - theta * m(t) / (Km + m(t)) - k * c(t) * m(t),
    c'(t) = R - k * c(t) * m(t),
    y(t)  = m(t)
)
report("C   carnitine depletion, c(t) unobserved", C)

# ---------------------------------------------------------------------
banner("SUMMARY")
for (k, v) in results
    println("\n", k)
    println("  global: ", v.global_id)
end
println("\nDone.")
