# =====================================================================
#  12.1 -- the channel-exchange theorem, checked against a general tool.
#
#  CLAIM UNDER TEST (derived in audit/e6_channel_exchange.py): for an
#  ordinate that is a sum of n channels sharing one functional form,
#  a permutation of the (amplitude, shape) pairs is an observational
#  symmetry exactly when it permutes channels carrying the SAME known
#  input. So a single selective input resolves the exchange only at
#  n = 2; at n = 3 it should leave a two-point fibre.
#
#  Expectations are NOT written here. Each model prints what the tool
#  returned; the comparison against the claim is made in the audit
#  record, from this output.
#
#  Run with:  julia symbolic/julia/run_channel_exchange.jl
# =====================================================================

using StructuralIdentifiability

banner(t) = (println("\n" * "="^70); println("  ", t); println("="^70))

function report(name, ode)
    banner(name)
    println("\nMODEL:"); println(ode)
    println("\n--- LOCAL identifiability ---")
    println(assess_local_identifiability(ode))
    println("\n--- GLOBAL identifiability ---")
    println(assess_identifiability(ode))
    println("\n--- IDENTIFIABLE FUNCTIONS ---")
    for f in find_identifiable_functions(ode)
        println("   ", f)
    end
    return nothing
end

# ---------------------------------------------------------------------
# N2.  Two saturating channels, NO input. The baseline degeneracy.
# ---------------------------------------------------------------------
N2_noinput = @ODEmodel(
    x'(t) = s(t) - a1 * x(t)/(K1 + x(t)) - a2 * x(t)/(K2 + x(t)),
    y(t)  = x(t)
)
report("N2  two saturating channels, no input", N2_noinput)

# ---------------------------------------------------------------------
# N2d. Two saturating channels, a known input d(t) on channel 1 ONLY.
# ---------------------------------------------------------------------
N2_input = @ODEmodel(
    x'(t) = s(t) - (a1 + d(t)) * x(t)/(K1 + x(t)) - a2 * x(t)/(K2 + x(t)),
    y(t)  = x(t)
)
report("N2d two saturating channels, known input on channel 1", N2_input)

# ---------------------------------------------------------------------
# N3d. THREE saturating channels, a known input on channel 1 ONLY.
#      This is the discriminating case.
# ---------------------------------------------------------------------
N3_input = @ODEmodel(
    x'(t) = s(t) - (a1 + d(t)) * x(t)/(K1 + x(t))
                 - a2 * x(t)/(K2 + x(t))
                 - a3 * x(t)/(K3 + x(t)),
    y(t)  = x(t)
)
report("N3d three saturating channels, known input on channel 1 only", N3_input)

# ---------------------------------------------------------------------
# N3dd. THREE saturating channels, DISTINCT known inputs on 1 and 2.
# ---------------------------------------------------------------------
N3_two_inputs = @ODEmodel(
    x'(t) = s(t) - (a1 + d(t)) * x(t)/(K1 + x(t))
                 - (a2 + e(t)) * x(t)/(K2 + x(t))
                 - a3 * x(t)/(K3 + x(t)),
    y(t)  = x(t)
)
report("N3dd three saturating channels, distinct inputs on 1 and 2", N3_two_inputs)

# ---------------------------------------------------------------------
# MIX. One saturating and one LINEAR channel, NO input. The note's own
#      configuration. A linear channel is the K -> infinity limit, so
#      the two shapes differ and no exchange should be available.
# ---------------------------------------------------------------------
MIX = @ODEmodel(
    x'(t) = s(t) - a1 * x(t)/(K1 + x(t)) - ke * x(t),
    y(t)  = x(t)
)
report("MIX one saturating + one linear channel, no input", MIX)

banner("Done.")
