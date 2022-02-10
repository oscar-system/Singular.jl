include("GroebnerWalkFinal.jl")

#function test(case::Int)
test_successfull = true
#if case == 1 || case == 99

dim = 4
ve = [1, 1, 1, 1]
StartOrd = ordering_as_matrix(:degrevlex, dim)
TarOrd = ordering_as_matrix(:lex, dim)
R, (a, b, c, d) = Singular.PolynomialRing(
    Singular.N_ZpField(32003),
    ["a", "b", "c", "d"],
    ordering = Singular.ordering_M(StartOrd),
)
S = change_order(R, TarOrd)

I = Singular.std(
    Singular.Ideal(
        R,
        [
            5 + a^2 * b + 2 * a * c^2 + a^4 + 3 * a^2 * b * c,
            c +
            5 * a * b +
            4 * c^2 +
            2 * a * c * d +
            2 * b^4 +
            3 * c^4 +
            2 * a * d^3,
        ],
    ),
    complete_reduction = true,
)
println(I)
println(groebnerwalk(I, StartOrd, TarOrd, :pertubed, 4))
# end
#end
