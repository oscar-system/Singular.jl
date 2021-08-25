using Pkg
include("GroebnerWalk.jl")

using Oscar


case = 99
success = true

if case == 1 || case == 99
    R, (x, y) = Oscar.Singular.PolynomialRing(
        QQ,
        ["x", "y"],
        ordering = Singular.ordering_M(ordering_as_matrix(:degrevlex, 2)),
    )

    f1 = x^2 - y^3
    f2 = x^3 - y^2 - x
    I = Singular.Ideal(R, [f1, f2])

    I = Oscar.Singular.std(I, complete_reduction = true)

    S, V = Oscar.Singular.PolynomialRing(QQ, ["x", "y"], ordering = :lex)

    @time J = fractal_walk(
        I,
        MonomialOrder(ordering_as_matrix(:degrevlex, 2), [1, 1], [0]),
        MonomialOrder(ordering_as_matrix(:lex, 2), [1, 0], [1, 0]),
    )

    @time L = groebnerwalk(
        I,
        ordering_as_matrix(:degrevlex, 2),
        ordering_as_matrix(:lex, 2),
        :standard
    )
    @time K = groebnerwalk(
        I,
        ordering_as_matrix(:degrevlex, 2),
        ordering_as_matrix(:lex, 2),
        :pertubed,
        2
    )

    @time M = groebnerwalk(
        I,
        ordering_as_matrix(:degrevlex, 2),
        ordering_as_matrix(:lex, 2),
        :generic
    )
    T1 = Oscar.Singular.Ideal(S, [change_ring(x, S) for x in gens(J)])
    T2 = Oscar.Singular.Ideal(S, [change_ring(x, S) for x in gens(K)])
    T3 = Oscar.Singular.Ideal(S, [change_ring(x, S) for x in gens(L)])
    T4 = Oscar.Singular.Ideal(S, [change_ring(x, S) for x in gens(M)])

    T0 = Singular.std(
        Oscar.Singular.Ideal(S, [change_ring(x, S) for x in gens(I)]),
        complete_reduction = true,
    )

    println("test fractal: ", equalitytest(T0, T1))
    println("test pertubed: ", equalitytest(T2, T0))
    println("test standard: ", equalitytest(T3, T0))
    println("test generic: ", equalitytest(T4, T0))


    if !(equalitytest(T2, T1) && equalitytest(T3, T4))
       success = false
    end
end

if case == 2 || case == 99

    R, (x, y, z) = Oscar.Singular.PolynomialRing(
        QQ,
        ["x", "y", "z"],
        ordering = Singular.ordering_M(ordering_as_matrix([5, 4, 1], :deglex)),
    )

    f1 = x^2 - y
    f2 = y^2 - x * z - y * z
    I = Singular.Ideal(R, [f1, f2])

    I = Oscar.Singular.std(I, complete_reduction = true)


    @time J = fractal_walk(
        I,
        MonomialOrder(ordering_as_matrix(:deglex, 3), [5, 4, 1], [0]),
        MonomialOrder(ordering_as_matrix(:lex, 3), [6, 1, 3], [6, 1, 3]),
    )
    @time K = groebnerwalk(
        I,
        ordering_as_matrix([5, 4, 1], :deglex),
        ordering_as_matrix([6, 1, 3], :lex),
        :pertubed,
        2
    )
    @time L = groebnerwalk(
        I,
        ordering_as_matrix([5, 4, 1], :deglex),
        ordering_as_matrix([6, 1, 3], :lex),
        :standard,
    )
    @time M = groebnerwalk(
        I,
        ordering_as_matrix([5, 4, 1], :deglex),
        ordering_as_matrix([6, 1, 3], :lex),
        :generic,
    )

    S, V = Oscar.Singular.PolynomialRing(
        QQ,
        ["x", "y", "z"],
        ordering = Singular.ordering_M(ordering_as_matrix([6, 1, 3], :lex)),
    )


    T0 = Singular.std(
        Oscar.Singular.Ideal(S, [change_ring(x, S) for x in gens(I)]),
        complete_reduction = true,
    )
    T1 = Oscar.Singular.Ideal(S, [change_ring(x, S) for x in gens(J)])
    T2 = Oscar.Singular.Ideal(S, [change_ring(x, S) for x in gens(K)])
    T3 = Oscar.Singular.Ideal(S, [change_ring(x, S) for x in gens(L)])
    T4 = Oscar.Singular.Ideal(S, [change_ring(x, S) for x in gens(M)])

    println("test fractal: ", equalitytest(T0, T1))
    println("test pertubed: ", equalitytest(T2, T0))
    println("test standard: ", equalitytest(T3, T0))
    println("test generic: ", equalitytest(T4, T0))


    if !(equalitytest(T2, T1) && equalitytest(T3, T4))
        success = false
    end
end


if case == 3 || case == 99
    R, (x, y, z, u, v, w) = Oscar.Singular.PolynomialRing(
        QQ,
        ["x", "y", "z", "u", "v", "w"],
        ordering = Singular.ordering_M(ordering_as_matrix(:degrevlex, 6)),
    )


    f1 = y^4 + y * z^2 - u^2 * w
    f2 = 2 * x^2 * y + x^3 * w * u^2 + x
    f3 = 2 - 3 * x^2 * v * z^4 * w
    I = Singular.Ideal(R, [f1, f2, f3])
    I = Oscar.Singular.std(I, complete_reduction = true)


    @time J = fractal_walk(
        I,
        MonomialOrder(
            ordering_as_matrix(:degrevlex, 6),
            [1, 1, 1, 1, 1, 1],
            [0],
        ),
        MonomialOrder(ordering_as_matrix(:lex, 6), [1, 0, 0, 0, 0, 0], [0]),
    )
        @time K = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 6),
            ordering_as_matrix(:lex, 6),
            :pertubed,
            3
        )
    @time L = groebnerwalk(
        I,
        ordering_as_matrix(:degrevlex, 6),
        ordering_as_matrix(:lex, 6),
        :standard,
    )
    @time M = groebnerwalk(
        I,
        ordering_as_matrix(:degrevlex, 6),
        ordering_as_matrix(:lex, 6),
        :generic,
    )


    S, V = Oscar.Singular.PolynomialRing(
        QQ,
        ["x", "y", "z", "u", "v", "w"],
        ordering = Oscar.Singular.ordering_M(ordering_as_matrix(:lex, 6)),
    )
    T0 = Singular.std(
        Oscar.Singular.Ideal(S, [change_ring(x, S) for x in gens(I)]),
        complete_reduction = true,
    )
    T1 = Oscar.Singular.Ideal(S, [change_ring(x, S) for x in gens(J)])
    T2 = Oscar.Singular.Ideal(S, [change_ring(x, S) for x in gens(K)])
    T3 = Oscar.Singular.Ideal(S, [change_ring(x, S) for x in gens(L)])
    T4 = Oscar.Singular.Ideal(S, [change_ring(x, S) for x in gens(M)])

    println("test fractal: ", equalitytest(T0, T1))
    println("test pertubed: ", equalitytest(T2, T0))
    println("test standard: ", equalitytest(T3, T0))
    println("test generic: ", equalitytest(T4, T0))


    if !(equalitytest(T2, T1) && equalitytest(T3, T4))
        success = false
    end
end
println("All tests were: ", success)
