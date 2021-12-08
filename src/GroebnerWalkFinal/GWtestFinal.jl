include("GroebnerWalkFinal.jl")

function test(case::Int)

    test_successfull = true

    if case == 1 || case == 99
        R, (x, y) = Singular.PolynomialRing(
            Singular.QQ,
            ["x", "y"],
            ordering = Singular.ordering_M(ordering_as_matrix(:degrevlex, 2)),
        )

        f1 = x^2 - y^3
        f2 = x^3 - y^2 - x
        I = Singular.Ideal(R, [f1, f2])

        I = Singular.std(I, complete_reduction = true)

        S, V = Singular.PolynomialRing(Singular.QQ, ["x", "y"], ordering = :lex)

        @time H = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 2),
            ordering_as_matrix(:lex, 2),
            :tran,
        )
        @time J = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 2),
            ordering_as_matrix(:lex, 2),
            :fractal,
        )
        @time JJ = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 2),
            ordering_as_matrix(:lex, 2),
            :fractal_look_ahead,
        )
        @time L = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 2),
            ordering_as_matrix(:lex, 2),
            :standard,
        )
        @time K = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 2),
            ordering_as_matrix(:lex, 2),
            :pertubed,
            2,
        )
        @time M = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 2),
            ordering_as_matrix(:lex, 2),
            :generic,
        )
        @time T0 = Singular.std(
            Singular.Ideal(S, [change_ring(x, S) for x in gens(I)]),
            complete_reduction = true,
        )
        T1 = Singular.Ideal(S, [change_ring(x, S) for x in gens(J)])
        T2 = Singular.Ideal(S, [change_ring(x, S) for x in gens(K)])
        T3 = Singular.Ideal(S, [change_ring(x, S) for x in gens(L)])
        T4 = Singular.Ideal(S, [change_ring(x, S) for x in gens(M)])
        T5 = Singular.Ideal(S, [change_ring(x, S) for x in gens(JJ)])
        T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(H)])

        println("test tran: ", equalitytest(T0, T6))
        println("test fractal: ", equalitytest(T0, T5))
        println("test fractal: ", equalitytest(T0, T1))
        println("test pertubed: ", equalitytest(T2, T0))
        println("test standard: ", equalitytest(T3, T0))
        println("test generic: ", equalitytest(T4, T0))

        if !(
            equalitytest(T2, T1) &&
            equalitytest(T3, T4) &&
            equalitytest(T6, T5)
        )
            test_successfull = false
        end
    end
    if case == 2 || case == 99

        R, (x, y, z) = Singular.PolynomialRing(
            Singular.QQ,
            ["x", "y", "z"],
            ordering = Singular.ordering_M(
                ordering_as_matrix([5, 4, 1], :deglex),
            ),
        )

        f1 = x^2 - y
        f2 = y^2 - x * z - y * z
        I = Singular.Ideal(R, [f1, f2])

        I = Singular.std(I, complete_reduction = true)
        @time H = groebnerwalk(
            I,
            ordering_as_matrix([5, 4, 1], :deglex),
            ordering_as_matrix([6, 1, 3], :lex),
            :tran,
        )
        @time J = groebnerwalk(
            I,
            ordering_as_matrix([5, 4, 1], :deglex),
            ordering_as_matrix([6, 1, 3], :lex),
            :fractal,
        )
        @time JJ = groebnerwalk(
            I,
            ordering_as_matrix([5, 4, 1], :deglex),
            ordering_as_matrix([6, 1, 3], :lex),
            :fractal_look_ahead,
        )
        @time K = groebnerwalk(
            I,
            ordering_as_matrix([5, 4, 1], :deglex),
            ordering_as_matrix([6, 1, 3], :lex),
            :pertubed,
            3,
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

        S, V = Singular.PolynomialRing(
            Singular.QQ,
            ["x", "y", "z"],
            ordering = Singular.ordering_M(ordering_as_matrix([6, 1, 3], :lex)),
        )

        T0 = Singular.std(
            Singular.Ideal(S, [change_ring(x, S) for x in gens(I)]),
            complete_reduction = true,
        )
        T1 = Singular.Ideal(S, [change_ring(x, S) for x in gens(J)])
        T2 = Singular.Ideal(S, [change_ring(x, S) for x in gens(K)])
        T3 = Singular.Ideal(S, [change_ring(x, S) for x in gens(L)])
        T4 = Singular.Ideal(S, [change_ring(x, S) for x in gens(M)])
        T5 = Singular.Ideal(S, [change_ring(x, S) for x in gens(JJ)])
        T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(H)])

        println("test tran: ", equalitytest(T0, T6))
        println("test fractal2: ", equalitytest(T0, T5))
        println("test fractal: ", equalitytest(T0, T1))
        println("test pertubed: ", equalitytest(T2, T0))
        println("test standard: ", equalitytest(T3, T0))
        println("test generic: ", equalitytest(T4, T0))

        if !(
            equalitytest(T2, T1) &&
            equalitytest(T3, T4) &&
            equalitytest(T6, T5)
        )
            test_successfull = false
        end
    end

    if case == 3 || case == 99
        R, (x, y, z, u, v, w) = Singular.PolynomialRing(
            Singular.QQ,
            ["x", "y", "z", "u", "v", "w"],
            ordering = Singular.ordering_M(ordering_as_matrix(:degrevlex, 6)),
        )

        f1 = y^4 + y * z^2 - u^2 * w
        f2 = 2 * x^2 * y + x^3 * w * u^2 + x
        f3 = 2 - 3 * x^2 * v * z^4 * w
        I = Singular.Ideal(R, [f1, f2, f3])
        I = Singular.std(I, complete_reduction = true)

        @time H = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 6),
            ordering_as_matrix(:lex, 6),
            :standard,
        )
        #=    @time J = fractal_walk(
                I,
                MonomialOrder(
                    ordering_as_matrix(:degrevlex, 6),
                    [1, 1, 1, 1, 1, 1],
                    [0],
                ),
                MonomialOrder(
                    ordering_as_matrix(:lex, 6),
                    [1, 0, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                ),
            ) =#


        @time JJ = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 6),
            ordering_as_matrix(:lex, 6),
            :fractal_look_ahead,
        )
        @time K = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 6),
            ordering_as_matrix(:lex, 6),
            :pertubed,
            4,
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

        S, V = Singular.PolynomialRing(
            Singular.QQ,
            ["x", "y", "z", "u", "v", "w"],
            ordering = Singular.ordering_M(ordering_as_matrix(:lex, 6)),
        )
        @time T0 = Singular.std(
            Singular.Ideal(S, [change_ring(x, S) for x in gens(I)]),
            complete_reduction = true,
        )

        #    T1test(4) = Singular.Ideal(S, [change_ring(x, S) for x in gens(J)])
        T2 = Singular.Ideal(S, [change_ring(x, S) for x in gens(K)])
        T3 = Singular.Ideal(S, [change_ring(x, S) for x in gens(L)])
        T4 = Singular.Ideal(S, [change_ring(x, S) for x in gens(M)])
        T5 = Singular.Ideal(S, [change_ring(x, S) for x in gens(JJ)])
        T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(H)])

        println("test tran: ", equalitytest(T0, T6))
        println("test fractal2: ", equalitytest(T0, T5))
        #println("test fractal: ", equalitytest(T0, T1))
        println("test pertubed: ", equalitytest(T2, T0))
        println("test standard: ", equalitytest(T3, T0))
        println("test generic: ", equal(T4, T0))

        if !(
            equalitytest(T2, T0) &&
            equalitytest(T3, T6) &&
            equalitytest(T0, T5) &&
            Singular.equal(T0, T4)
        )
            test_successfull = false
        end
    end
    if case == 4 || case == 99
        R, (q, c, p, d) = Singular.PolynomialRing(
            Singular.QQ,
            ["q", "c", "p", "d"],
            ordering = Singular.ordering_M(ordering_as_matrix(:degrevlex, 4)),
        )


        f1 =
            2 * c * d * p * q - 4 * c * d + 2 * p * q + 4 * p * q * c -
            4 * d^2 * q - 2 * d^2 * p * q + p^2 * d^2 - 2 * c^2 * q +
            c^2 * q^2 +
            2 * c * d * q - 2 * c * d * q^2 + d^2 * q^2 - 2 * c * d * p -
            8 * p +
            c^2 +
            4 * d^2 - 2 * q +
            10 * p^2 +
            2
        f2 =
            2 * d * p * q + 4 * d * p^2 + d * p - 7 * d * p + c * p -
            3 * p * q * c + 4 * c
        f3 = -2 * p^2 + 8 * p - 2 - 2 * p * q - 2 * q
        f4 =
            4 - 4 * p - 4 * q^2 + 3 * c^2 * q^2 - 6 * c^2 * q +
            3 * c^2 +
            9 * p^2 * d^2 +
            6 * d^2 * p * q - 3d^2 * q^2 - 24 * p * d^2 +
            12 * d^2 +
            4 * p^2 +
            12 * c * d * p +
            12 * c * d * q +
            12 * c * d * p * q - 12 * c * d
        I = Singular.Ideal(R, [f1, f2, f3, f4])
        I = Singular.std(I, complete_reduction = true)

        @time H = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 4),
            ordering_as_matrix(:lex, 4),
            :tran,
        )
        @time J = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 4),
            ordering_as_matrix(:lex, 4),
            :fractal,
        )
        @time JJ = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 4),
            ordering_as_matrix(:lex, 4),
            :fractal_look_ahead,
        )

        @time K = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 4),
            ordering_as_matrix(:lex, 4),
            :pertubed,
            4,
        )
        @time L = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 4),
            ordering_as_matrix(:lex, 4),
            :fractal_lex,
        )
        @time M = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 4),
            ordering_as_matrix(:lex, 4),
            :generic,
        )


        S, V = Singular.PolynomialRing(
            Singular.QQ,
            ["q", "c", "p", "d"],
            ordering = Singular.ordering_M(ordering_as_matrix(:lex, 4)),
        )
        @time T0 = Singular.std(
            Singular.Ideal(S, [change_ring(x, S) for x in gens(I)]),
            complete_reduction = true,
        )
        T1 = Singular.Ideal(S, [change_ring(x, S) for x in gens(J)])
        T2 = Singular.Ideal(S, [change_ring(x, S) for x in gens(K)])
        T3 = Singular.Ideal(S, [change_ring(x, S) for x in gens(L)])
        T4 = Singular.Ideal(S, [change_ring(x, S) for x in gens(M)])
        T5 = Singular.Ideal(S, [change_ring(x, S) for x in gens(JJ)])
        T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(H)])

        println("test tran: ", equalitytest(T0, T6))
        println("test fractal: ", equalitytest(T0, T1))
        println("test fractal2: ", equalitytest(T0, T5))
        println("test pertubed: ", equalitytest(T2, T0))
        println("test standard: ", equalitytest(T3, T0))
        println("test generic: ", Singular.equal(T4, T0))


        if !(
            equalitytest(T2, T1) &&
            equalitytest(T3, T6) &&
            equalitytest(T0, T5) &&
            Singular.equal(T0, T4)
        )
            test_successfull = false
        end
    end
    if case == 5 || case == 99
        R, (x, y, z) = Singular.PolynomialRing(
            Singular.QQ,
            ["x", "y", "z"],
            ordering = Singular.ordering_M(ordering_as_matrix(:degrevlex, 3)),
        )

        f1 = x^2 + x * y^2 * z - 2x * y + y^4 + y^2 + z^2
        f2 = -x^3 * y^2 + x * y^2 * z + x * y * z^3 - 2x * y + y^4
        f3 = -2x^2 * y + x * y^4 + y * z^4 - 3
        I = Singular.Ideal(R, [f1, f2, f3])

        I = Singular.std(I, complete_reduction = true)

        S, V = Singular.PolynomialRing(
            Singular.QQ,
            ["x", "y", "z"],
            ordering = :lex,
        )

        @time H = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 3),
            ordering_as_matrix(:lex, 3),
            :tran,
        )
        @time J = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 3),
            ordering_as_matrix(:lex, 3),
            :fractal,
        )
        @time JJ = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 3),
            ordering_as_matrix(:lex, 3),
            :fractal_look_ahead,
        )

        @time L = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 3),
            ordering_as_matrix(:lex, 3),
            :fractal_lex,
        )
        @time K = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 3),
            ordering_as_matrix(:lex, 3),
            :pertubed,
            3,
        )
        @time M = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 3),
            ordering_as_matrix(:lex, 3),
            :generic,
        )
        T1 = Singular.Ideal(S, [change_ring(x, S) for x in gens(J)])
        T2 = Singular.Ideal(S, [change_ring(x, S) for x in gens(K)])
        T3 = Singular.Ideal(S, [change_ring(x, S) for x in gens(L)])
        T4 = Singular.Ideal(S, [change_ring(x, S) for x in gens(M)])
        T5 = Singular.Ideal(S, [change_ring(x, S) for x in gens(JJ)])

        @time T0 = Singular.std(
            Singular.Ideal(S, [change_ring(x, S) for x in gens(I)]),
            complete_reduction = true,
        )
        T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(H)])

        println("test tran: ", equalitytest(T0, T6))
        println("test fractal2: ", equalitytest(T0, T5))
        println("test fractal: ", equalitytest(T0, T1))
        println("test pertubed: ", equalitytest(T2, T0))
        println("test standard: ", equalitytest(T3, T0))
        println("test generic: ", Singular.equal(T4, T0))

        if !(
            equalitytest(T2, T1) &&
            equalitytest(T3, T6) &&
            equalitytest(T0, T5) &&
            Singular.equal(T0, T4)
        )
            test_successfull = false
        end

    end



    if case == 6 || case == 99
        R, (x, y, z, u, v, w) = Singular.PolynomialRing(
            Singular.QQ,
            ["x", "y", "z", "u", "v", "w"],
            ordering = Singular.ordering_M(ordering_as_matrix(:degrevlex, 6)),
        )


        f1 = 2 * w * y + 2 * v * z + u^2 + x^2 + x
        f2 = 2 * w * z + 2 * v * u + 2 * y * x + y
        f3 = 2 * w * u + v^2 + 2 * z * x + z + y^2
        f4 = 2 * w * v + 2 * u * x + u + 2 * z * y
        f5 = w^2 + 2 * v * x + v + 2 * u * y + z^2
        f6 = 2 * w * x + w + 2 * v * y + 2 * u * z
        I = Singular.Ideal(R, [f1, f2, f3, f4, f5, f6])
        I = Singular.std(I, complete_reduction = true)

        @time H = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 6),
            ordering_as_matrix(:lex, 6),
            :standard,
        )
        #=    @time J = fractal_walk(
                I,
                MonomialOrder(
                    ordering_as_matrix(:degrevlex, 6),
                    [1, 1, 1, 1, 1, 1],
                    [0],
                ),
                MonomialOrder(
                    ordering_as_matrix(:lex, 6),
                    [1, 0, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                ),
            ) =#
        @time JJ = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 6),
            ordering_as_matrix(:lex, 6),
            :fractal_look_ahead,
        )

        @time K = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 6),
            ordering_as_matrix(:lex, 6),
            :pertubed,
            6,
        )
        @time L = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 6),
            ordering_as_matrix(:lex, 6),
            :fractal,
        )
        @time M = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 6),
            ordering_as_matrix(:lex, 6),
            :generic,
        )
        S, V = Singular.PolynomialRing(
            Singular.QQ,
            ["x", "y", "z", "u", "v", "w"],
            ordering = Singular.ordering_M(ordering_as_matrix(:lex, 6)),
        )
        @time T0 = Singular.std(
            Singular.Ideal(S, [change_ring(x, S) for x in gens(I)]),
            complete_reduction = true,
        )
        #    T1test(4) = Singular.Ideal(S, [change_ring(x, S) for x in gens(J)])
        T2 = Singular.Ideal(S, [change_ring(x, S) for x in gens(K)])
        T3 = Singular.Ideal(S, [change_ring(x, S) for x in gens(L)])
        T4 = Singular.Ideal(S, [change_ring(x, S) for x in gens(M)])
        T5 = Singular.Ideal(S, [change_ring(x, S) for x in gens(JJ)])
        T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(H)])

        println("test standard: ", equalitytest(T0, T6))
        println("test fractal2: ", equalitytest(T0, T5))
        #println("test fractal: ", equalitytest(T0, T1))
        println("test pertubed: ", equalitytest(T2, T0))
        println("test standard: ", equalitytest(T3, T0))
        println("test generic: ", equal(T4, T0))

        if !(
            equalitytest(T2, T0) &&
            equalitytest(T3, T6) &&
            equalitytest(T0, T5) &&
            Singular.equal(T0, T4)
        )
            test_successfull = false
        end

    end
    println("All tests were: ", test_successfull)
end
