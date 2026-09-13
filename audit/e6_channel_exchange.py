"""
12.1 -- does the Model B result generalise?

Model B gave: two saturating channels are exchangeable, and a known input on one
of them breaks the exchange, the difference being

    Delta * m * (K_e - K_m) / ((K_m + m)(K_e + m)).

This asks whether that is an accident of two Michaelis-Menten channels or an
instance of something general. Nothing below is asserted; each claim is derived
and then checked.

SETTING. Let the observed ordinate be a sum of n channels sharing one functional
form f, each with its own amplitude a_i, its own shape constant K_i, and its own
KNOWN exogenous input Delta_i(t) riding on the amplitude:

    C(m, t) = sum_i ( a_i + Delta_i(t) ) f(m; K_i).

The unknown parameter is the list of pairs p = ((a_1,K_1), ..., (a_n,K_n)). The
inputs are indexed by CHANNEL and are not part of the parameter: they stay where
they are when the parameter is permuted. That is the whole content of "an input
on one channel".
"""
import itertools
import sympy as sp

m, t = sp.symbols('m t', positive=True)


def banner(s):
    print("\n" + "=" * 78)
    print("  " + s)
    print("=" * 78)


# ---------------------------------------------------------------------------
def part1_the_general_difference(n=3):
    banner(f"PART 1  --  the difference under a permutation, for general f, n = {n}")

    a = sp.symbols(f'a1:{n+1}', positive=True)
    K = sp.symbols(f'K1:{n+1}', positive=True)
    D = sp.symbols(f'Delta1:{n+1}')            # known inputs, one per CHANNEL
    f = sp.Function('f')                       # left completely general

    def C(perm):
        # p^pi has (a^pi_i, K^pi_i) = (a_{pi(i)}, K_{pi(i)}); inputs stay put
        return sum((a[perm[i]] + D[i]) * f(m, K[perm[i]]) for i in range(n))

    ident = tuple(range(n))
    print("  C(p)      = sum_i (a_i + Delta_i) f(m; K_i)")
    print("  C(p^pi)   = sum_i (a_{pi(i)} + Delta_i) f(m; K_{pi(i)})")
    print("            = sum_j (a_j + Delta_{pi^-1(j)}) f(m; K_j)   [reindex j = pi(i)]")
    print()
    print("  so the difference collects on the SHAPE functions, and the amplitudes")
    print("  cancel identically:")
    print()
    for perm in itertools.permutations(range(n)):
        if perm == ident:
            continue
        diff = sp.simplify(sp.expand(C(perm) - C(ident)))
        inv = [perm.index(j) for j in range(n)]
        pred = sp.simplify(sp.expand(
            sum((D[inv[j]] - D[j]) * f(m, K[j]) for j in range(n))))
        ok = sp.simplify(diff - pred) == 0
        print(f"    pi = {perm}   C(p^pi) - C(p) == sum_j (D_[pi^-1(j)] - D_j) f(m;K_j) ?  {ok}")
    print()
    print("  Every permutation, every n: the difference is a linear combination of")
    print("  the shape functions with coefficients (Delta_{pi^-1(j)} - Delta_j).")
    print("  f was never specialised, so this step is independent of the channel's")
    print("  functional form.")


# ---------------------------------------------------------------------------
def part2_linear_independence():
    banner("PART 2  --  when do those coefficients have to vanish separately?")

    print("  Only if the shape functions are linearly independent. For")
    print("  Michaelis-Menten f(m;K) = m/(K+m) with DISTINCT positive K_i:")
    n = 3
    K = sp.symbols(f'K1:{n+1}', positive=True)
    c = sp.symbols(f'c1:{n+1}')
    expr = sum(c[i] * m / (K[i] + m) for i in range(n))
    num = sp.Poly(sp.numer(sp.cancel(sp.together(expr))), m)
    sols = sp.solve([sp.Eq(co, 0) for co in num.all_coeffs()], list(c), dict=True)
    print(f"    sum_i c_i m/(K_i+m) == 0 for all m  =>  {sols}")
    print()
    print("  and the reason, which does not need a solver: multiplying through by")
    print("  prod_j (K_j + m) gives a polynomial identity; evaluating it at")
    print("  m = -K_i kills every term but the i-th and leaves")
    print("      c_i * (-K_i) * prod_{j != i} (K_j - K_i) = 0,")
    print("  non-zero factors throughout when the K_i are distinct and positive.")
    print("  So c_i = 0. Distinctness of the shape constants is exactly what is")
    print("  needed, and equal K_i is the degenerate case where two channels are")
    print("  one channel and only the SUM of their amplitudes exists.")


# ---------------------------------------------------------------------------
def part3_the_theorem(n=3):
    banner(f"PART 3  --  the theorem, and the residual group, n = {n}")

    print("  Combining: pi is an observational symmetry")
    print("      <=>  Delta_{pi^-1(j)} = Delta_j  for every j")
    print("      <=>  pi permutes only channels carrying the SAME input.")
    print()
    print("  So the residual symmetry group is the Young subgroup prod_B S_|B|,")
    print("  B running over the blocks of channels sharing an input, and the fibre")
    print("  has prod_B |B|! points.")
    print()
    print("  Checked exhaustively against the symbolic difference, on Michaelis-")
    print("  Menten channels with distinct K, for each input pattern:")

    a = sp.symbols(f'a1:{n+1}', positive=True)
    K = sp.symbols(f'K1:{n+1}', positive=True)
    g = lambda KK: m / (KK + m)

    # input patterns given as a signature: equal labels = equal input functions
    d = sp.Symbol('Delta', positive=True)
    e = sp.Symbol('Epsilon', positive=True)
    patterns = {
        "no input anywhere":            [0, 0, 0],
        "one input, channel 1":         [d, 0, 0],
        "the SAME input on 1 and 2":    [d, d, 0],
        "DIFFERENT inputs on 1 and 2":  [d, e, 0],
        "a different input on each":    [d, e, d + e],
    }
    print(f"\n    {'input pattern':<30} {'predicted |fibre|':>18} {'found':>8}   match")
    for name, D in patterns.items():
        blocks = {}
        for i, Di in enumerate(D):
            blocks.setdefault(sp.srepr(sp.sympify(Di)), []).append(i)
        pred = 1
        for b in blocks.values():
            pred *= sp.factorial(len(b))
        found = 0
        for perm in itertools.permutations(range(n)):
            diff = sp.simplify(
                sum((a[perm[i]] + D[i]) * g(K[perm[i]]) for i in range(n))
                - sum((a[i] + D[i]) * g(K[i]) for i in range(n)))
            if diff == 0:
                found += 1
        print(f"    {name:<30} {int(pred):>18} {found:>8}   {'ok' if int(pred)==found else 'MISMATCH'}")


# ---------------------------------------------------------------------------
def part3b_amplitudes_re_absorb():
    """The refinement a pure-permutation test cannot see.

    Parts 1 and 3 ask whether permuting the pairs leaves C unchanged. The
    operative question is whether some OTHER parameter vector gives the same C,
    and the amplitudes are free. Since a_i is an additive constant on channel i
    and Delta_i adds to it, any CONSTANT part of Delta_i is indistinguishable
    from a_i, so a permutation can be repaired by shifting the amplitudes.
    """
    banner("PART 3b  --  amplitudes re-absorb a constant input, so PRESENCE of "
           "an input is not enough")

    a1, a2, K1, K2 = sp.symbols('a1 a2 K1 K2', positive=True)
    b1, b2 = sp.symbols('b1 b2', real=True)
    c = sp.Symbol('c', positive=True)
    f = lambda K: m / (K + m)

    print("  A KNOWN CONSTANT input c on channel 1. Swap the shapes, leave the")
    print("  amplitudes free:")
    resid = ((a1 + c) * f(K1) + a2 * f(K2)) - ((b1 + c) * f(K2) + b2 * f(K1))
    num = sp.Poly(sp.numer(sp.cancel(sp.together(resid))), m)
    sol = sp.solve([sp.Eq(co, 0) for co in num.all_coeffs()], [b1, b2], dict=True)
    print(f"    solution: {sol}")
    if sol:
        print(f"    residual there: {sp.simplify(resid.subs(sol[0]))}")
    print("  => a constant known input breaks NOTHING. The swap is absorbed by")
    print("     b_1 = a_2 - c, b_2 = a_1 + c, and the fibre is still two points.")

    print()
    print("  Now an input taking TWO DISTINCT VALUES. One (b_1,b_2) must serve")
    print("  both, so collect on powers of d as well as of m:")
    d = sp.Symbol('d', real=True)
    resid = ((a1 + d) * f(K1) + a2 * f(K2)) - ((b1 + d) * f(K2) + b2 * f(K1))
    num = sp.Poly(sp.numer(sp.cancel(sp.together(resid))), m)
    eqs = []
    for co in num.all_coeffs():
        for dco in sp.Poly(sp.expand(co), d).all_coeffs():
            if sp.expand(dco) != 0:
                eqs.append(sp.expand(dco))
    for e in dict.fromkeys(eqs):
        print(f"      {sp.factor(e)} = 0")
    print("  The d-linear equation is K_2 - K_1 = 0, so with distinct shapes")
    print("  there is no solution. Demanding one (b_1,b_2) at d = 100 and 400:")
    r = []
    for dv in (100, 400):
        r += sp.Poly(sp.numer(sp.cancel(sp.together(resid.subs(d, dv)))), m).all_coeffs()
    got = sp.solve([sp.Eq(x, 0) for x in r], [b1, b2], dict=True)
    print(f"      {got if got else 'no solution  =>  the exchange IS broken'}")

    print()
    print("  SO THE CONDITION IS SHARPER THAN PART 3 STATES:")
    print()
    print("      the exchange (i j) is broken")
    print("        <=>  Delta_i - Delta_j is NON-CONSTANT in t")
    print("        AND  f(.;K_i) is not identically f(.;K_j).")
    print()
    print("  Part 3's 'same input' must be read as 'same input up to an additive")
    print("  constant'. Two channels are separated only by what VARIES between")
    print("  them. That is the same statement the note already makes about the")
    print("  load record -- day-to-day variability, not precision, closes the")
    print("  load-bias degeneracy -- reached here by a completely different route.")


# ---------------------------------------------------------------------------
def part4_the_two_channel_instance():
    banner("PART 4  --  the Model B result as the n = 2 instance")

    th, Km, ke, Ke, d = sp.symbols('theta K_m k_e K_e Delta', positive=True)
    g = lambda KK: m / (KK + m)
    diff = sp.simplify(sp.together(
        ((th + d) * g(Km) + ke * g(Ke)) - ((ke + d) * g(Ke) + th * g(Km))))
    print("  C(p) - C(p^swap) with the input on channel 1 only:")
    sp.pprint(diff)
    print()
    general = sp.simplify(d * (g(Km) - g(Ke)) - diff)
    print(f"  equals  Delta * ( f(m;K_m) - f(m;K_e) )  ?   {general == 0}")
    print()
    print("  which is the general formula with n = 2, Delta_1 = Delta, Delta_2 = 0.")
    print("  The (K_e - K_m) numerator reported for Model B is what")
    print("  f(m;K_m) - f(m;K_e) happens to look like for Michaelis-Menten; it")
    print("  carries no information the general statement does not.")


# ---------------------------------------------------------------------------
def part5_where_it_stops(n=3):
    banner("PART 5  --  the scope limit, exhibited rather than described")

    print("  For n >= 3 a single selective input does NOT resolve the exchange.")
    print("  An input on channel 1 breaks every transposition involving 1 and")
    print("  leaves S_{n-1} on the rest. Explicitly, at n = 3 with an input on")
    print("  channel 1, these two parameter vectors are observationally identical:")

    a = [sp.Integer(300), sp.Integer(50), sp.Integer(8)]
    K = [sp.Integer(100), sp.Integer(20), sp.Integer(700)]
    d = sp.Symbol('Delta', positive=True)
    D = [d, 0, 0]
    g = lambda KK: m / (KK + m)
    C0 = sum((a[i] + D[i]) * g(K[i]) for i in range(n))
    swap23 = (0, 2, 1)
    C1 = sum((a[swap23[i]] + D[i]) * g(K[swap23[i]]) for i in range(n))
    print(f"\n    p  = {[(int(a[i]), int(K[i])) for i in range(n)]}")
    print(f"    p' = {[(int(a[swap23[i]]), int(K[swap23[i]])) for i in range(n)]}"
          "   (channels 2 and 3 exchanged)")
    print(f"\n    C(p) - C(p') simplifies to: {sp.simplify(C0 - C1)}")
    for mv in (1, 10, 100, 1000):
        v0 = C0.subs({m: mv, d: 400}); v1 = C1.subs({m: mv, d: 400})
        print(f"      m = {mv:5d}   C(p) = {sp.nsimplify(v0).evalf():.12f}   "
              f"C(p') = {sp.nsimplify(v1).evalf():.12f}")
    print("\n  Identical for every m and every Delta. One input, three channels:")
    print("  the fibre drops from 3! = 6 points to 2, not to 1.")
    print()
    print("  To break S_n completely, the inputs must separate every pair, which")
    print("  needs n-1 channels carrying pairwise-distinct inputs. An input")
    print("  duplicated across two channels breaks neither of them against the")
    print("  other -- row 'the SAME input on 1 and 2' in Part 3.")


# ---------------------------------------------------------------------------
def part6_scope_for_this_note():
    banner("PART 6  --  what this says about the note's own model. A SCOPE LIMIT.")

    print("  The note's alternative channel is LINEAR, so its two channels are not")
    print("  two members of one parameterised family: the enzymatic term carries")
    print("  (theta, K_m) and the alternative term carries k_e alone. There is no")
    print("  permutation of (theta, K_m, k_e) that exchanges them, so the note's")
    print("  model has no channel-exchange degeneracy to break. The tool agrees:")
    print("  saturating + linear with NO input returns globally identifiable.")
    print()
    print("  An earlier version of this script argued that from the limit")
    print("      lim_(K_e -> oo) k_e K_e m/(K_e + m) = k_e m,")
    print("  i.e. that the shapes differ. That reasoning is poor and is withdrawn.")
    print("  The limit is singular -- it needs k_e K_e held fixed while K_e grows,")
    print("  and K_e = oo is not a value of a positive parameter -- so it does not")
    print("  put the two channels in one family, and 'the shapes differ' is not")
    print("  what the argument turns on. What it turns on is that there is no")
    print("  permutation available at all.")
    print()
    print("  AND THE LIMIT CUTS THE OTHER WAY, WHICH IS THE POINT WORTH KEEPING.")
    print("  Model B contains the note's model only as a boundary case. At any")
    print("  FINITE K_e != K_m with no input, the exchange is present and exact.")
    print("  So the note's freedom from it is a property OF THE LINEARISATION,")
    print("  not evidence that the linearisation was harmless. That is a scope")
    print("  limit on the note and should be reported as one.")
    print()
    print("  It is not fatal, for a reason that is itself a result: the note's")
    print("  therapy is a decaying dose train, so Delta VARIES, and by Part 3b a")
    print("  varying selective input breaks the exchange whenever the shapes")
    print("  differ. A note built on a saturating alternative channel would")
    print("  therefore still identify -- but it would depend on the therapy for")
    print("  something this note gets for free, and a CONSTANT infusion would")
    print("  not do it.")


if __name__ == "__main__":
    print(__doc__)
    part1_the_general_difference()
    part2_linear_independence()
    part3_the_theorem()
    part3b_amplitudes_re_absorb()
    part4_the_two_channel_instance()
    part5_where_it_stops()
    part6_scope_for_this_note()
