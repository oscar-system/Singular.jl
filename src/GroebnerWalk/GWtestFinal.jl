include("GroebnerWalk.jl")
include("Examples")
function test(case::Int)

    test_successfull = true
    if case == 1 || case == 99
        id = katsura5()
        R = base_ring(id)
        dim = nvars(R)
        ve = ones(Int, dim)
        StartOrd = ordering_as_matrix(:degrevlex, dim)
        TarOrd = ordering_as_matrix(:lex, dim)
        R2 = change_order(R, StartOrd)
        S = change_order(R, TarOrd)
        I = Singular.std(
            Singular.Ideal(R2, [change_ring(x, R2) for x in gens(id)]),
            complete_reduction = true,
        )
        #I = Singular.std(id, complete_reduction = true)
        @time T = groebnerwalk(
            deepcopy(I),
            ordering_as_matrix(:degrevlex, dim),
            ordering_as_matrix(:lex, dim),
            :generic,
            2,
        )
        @time F = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, dim),
            ordering_as_matrix(:lex, dim),
            :fractal,
            2,
        )



        @time FA = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, dim),
            ordering_as_matrix(:lex, dim),
            :fractal_combined,
        )

        @time St = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, dim),
            ordering_as_matrix(:lex, dim),
            :pertubed,
            3,
        )

        @time Pe = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, dim),
            ordering_as_matrix(:lex, dim),
            :pertubed,
            4,
        )
        @time Ge = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, dim),
            ordering_as_matrix(:lex, dim),
            :generic,
            4,
        )

        #@time T0 = Singular.std(
        #    Singular.Ideal(S, [change_ring(x, S) for x in gens(I)]),
        #    complete_reduction = true,
        #)

        T1 = Singular.Ideal(S, [change_ring(x, S) for x in gens(T)])
        T2 = Singular.Ideal(S, [change_ring(x, S) for x in gens(F)])
        T3 = Singular.Ideal(S, [change_ring(x, S) for x in gens(FA)])
        T4 = Singular.Ideal(S, [change_ring(x, S) for x in gens(St)])
        T5 = Singular.Ideal(S, [change_ring(x, S) for x in gens(Pe)])
        T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(Ge)])




        println("test tran: ", equalitytest(T6, T1))
        println("test fractal: ", equalitytest(T6, T2))
        println("test fractal: ", equalitytest(T6, T3))
        println("test pertubed: ", equalitytest(T5, T6))
        println("test standard: ", equalitytest(T4, T6))
        println("test generic: ", equalitytest(T6, T6))
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
            ordering = Singular.ordering_M(ordering_as_matrix(:degrevlex, 3)),
        )


        f1 = 16 + 3 * x^2 + 16 * x^2 * z + 14 * x^2 * y^3
        f2 = 6 + y^3 * z + 17 * x^2 * z^2 + 7 * x * y^2 * z^2 + 13 * x^3 * z^2
        J = Singular.Ideal(R, [f1, f2])



        I = Singular.std(J, complete_reduction = true)
        @time T = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 3),
            ordering_as_matrix(:lex, 3),
            :tran,
            3,
        )

        @time F = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 3),
            ordering_as_matrix(:lex, 3),
            :fractal_start_order,
        )


        @time FA = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 3),
            ordering_as_matrix(:lex, 3),
            :fractal_combined,
        )
        @time Pe = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 3),
            ordering_as_matrix(:lex, 3),
            :pertubed,
            2,
        )
        @time St = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 3),
            ordering_as_matrix(:lex, 3),
            :standard,
        )
        @time Ge = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 3),
            ordering_as_matrix(:lex, 3),
            :generic,
        )


        S, V = Singular.PolynomialRing(
            Singular.QQ,
            ["x", "y", "z"],
            ordering = Singular.ordering_M(ordering_as_matrix(:lex, 3)),
        )

        T0 = Singular.std(
            Singular.Ideal(S, [change_ring(x, S) for x in gens(I)]),
            complete_reduction = true,
        )

        T1 = Singular.Ideal(S, [change_ring(x, S) for x in gens(T)])
        T2 = Singular.Ideal(S, [change_ring(x, S) for x in gens(F)])
        T3 = Singular.Ideal(S, [change_ring(x, S) for x in gens(FA)])
        T4 = Singular.Ideal(S, [change_ring(x, S) for x in gens(St)])
        T5 = Singular.Ideal(S, [change_ring(x, S) for x in gens(Pe)])
        T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(Ge)])
        println("test tran: ", equalitytest(T0, T1))
        println("test fractal: ", equalitytest(T0, T2))
        println("test fractal: ", equalitytest(T0, T3))
        println("test pertubed: ", equalitytest(T5, T0))
        println("test standard: ", equalitytest(T4, T0))
        println("test generic: ", equalitytest(T6, T0))



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
            :pertubed,
            4,
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
            :fractal_combined,
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
            equalitytest(T0, T4)
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



        @time T = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 4),
            ordering_as_matrix(:lex, 4),
            :tran,
        )
        @time F = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 4),
            ordering_as_matrix(:lex, 4),
            :pertubed,
            3
        )

        @time Fa = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 4),
            ordering_as_matrix(:lex, 4),
            :fractal_combined,
        )

        @time Pe = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 4),
            ordering_as_matrix(:lex, 4),
            :pertubed,
            4,
        )
        @time St = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 4),
            ordering_as_matrix(:lex, 4),
            :standard,
        )
        @time Ge = groebnerwalk(
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
        T1 = Singular.Ideal(S, [change_ring(x, S) for x in gens(T)])
        T2 = Singular.Ideal(S, [change_ring(x, S) for x in gens(F)])
        T3 = Singular.Ideal(S, [change_ring(x, S) for x in gens(Fa)])
        T4 = Singular.Ideal(S, [change_ring(x, S) for x in gens(St)])
        T5 = Singular.Ideal(S, [change_ring(x, S) for x in gens(Pe)])
        T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(Ge)])

        println("test tran: ", equalitytest(T0, T1))
        println("test fractal: ", equalitytest(T0, T2))
        println("test fractal: ", equalitytest(T0, T3))
        println("test pertubed: ", equalitytest(T5, T0))
        println("test standard: ", equalitytest(T4, T0))
        println("test generic: ", equalitytest(T6, T0))

        if !(
            equalitytest(T2, T1) &&
            equalitytest(T3, T6) &&
            equalitytest(T0, T5) &&
            equalitytest(T0, T4)
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

        @time T = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 3),
            ordering_as_matrix(:lex, 3),
            :tran,
        )
        @time F = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 3),
            ordering_as_matrix(:lex, 3),
            :fractal_combined,
        )
        @time FA = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 3),
            ordering_as_matrix(:lex, 3),
            :pertubed,
            3
        )

        @time St = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 3),
            ordering_as_matrix(:lex, 3),
            :standard,
        )
        @time Pe = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 3),
            ordering_as_matrix(:lex, 3),
            :pertubed,
            3,
        )



        @time Ge = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 3),
            ordering_as_matrix(:lex, 3),
            :generic,
        )

        @time T0 = Singular.std(
            Singular.Ideal(S, [change_ring(x, S) for x in gens(I)]),
            complete_reduction = true,
        )

        T1 = Singular.Ideal(S, [change_ring(x, S) for x in gens(T)])
        T2 = Singular.Ideal(S, [change_ring(x, S) for x in gens(F)])
        T3 = Singular.Ideal(S, [change_ring(x, S) for x in gens(FA)])
        T4 = Singular.Ideal(S, [change_ring(x, S) for x in gens(St)])
        T5 = Singular.Ideal(S, [change_ring(x, S) for x in gens(Pe)])
        T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(Ge)])

        println("test tran: ", equalitytest(T0, T1))
        println("test fractal: ", equalitytest(T0, T2))
        println("test fractal: ", equalitytest(T0, T3))
        println("test pertubed: ", equalitytest(T5, T0))
        println("test standard: ", equalitytest(T4, T0))
        println("test generic: ", equalitytest(T6, T0))
        if !(
            equalitytest(T2, T1) &&
            equalitytest(T3, T6) &&
            equalitytest(T0, T5) &&
            equalitytest(T0, T4)
        )
            test_successfull = false
        end

    end



    if case == 6 || case == 99
        R, (x, y, z, u, v, w) = Singular.PolynomialRing(
            Singular.QQ,
            ["x", "y", "z", "u", "v", "w"],
            ordering = :degrevlex,
        )


        f1 = 2 * w * y + 2 * v * z + u^2 + x^2 + x
        f2 = 2 * w * z + 2 * v * u + 2 * y * x + y
        f3 = 2 * w * u + v^2 + 2 * z * x + z + y^2
        f4 = 2 * w * v + 2 * u * x + u + 2 * z * y
        f5 = w^2 + 2 * v * x + v + 2 * u * y + z^2
        f6 = 2 * w * x + w + 2 * v * y + 2 * u * z
        J = Singular.Ideal(R, [f1, f2, f3, f4, f5, f6])
        I = Singular.std(J, complete_reduction = true)

        @time H = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 6),
            ordering_as_matrix(:lex, 6),
            :standard,
        )
        #=        @time J = fractal_walk_combined(
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
                )=#

        @time JJ = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, 6),
            ordering_as_matrix(:lex, 6),
            :fractal_combined,
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
            :fractal_combined,
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
            equalitytest(T0, T4)
        )
            test_successfull = false
        end

    end

    if case == 7 || case == 99
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
        I = Singular.Ideal(
            R,
            [
                2 * a^2 * b +
                3 * a * b^2 +
                3 * b^3 +
                4 * c^3 +
                4 * a * b * d +
                c^2 * d +
                2 * b * d^2 +
                2 * d^3 +
                4 * c^2 +
                2 * c * d +
                2 * c,
                2 * a^2 * b +
                5 * a^2 * c +
                2 * b * c^2 +
                c^2 * d +
                a * c +
                2 * b * d +
                2 * c * d +
                d^2 +
                a +
                4 * d,
            ],
        )


        I = Singular.std(I, complete_reduction = true)


        @time T = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, dim),
            ordering_as_matrix(:lex, dim),
            :fractal_combined,
            4,
        )
        @time F = groebnerwalk(
            deepcopy(I),
            ordering_as_matrix(:degrevlex, dim),
            ordering_as_matrix(:lex, dim),
            :fractal_start_order,
        )
        @time FA = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, dim),
            ordering_as_matrix(:lex, dim),
            :pertubed,
            3
        )
        @time St = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, dim),
            ordering_as_matrix(:lex, dim),
            :pertubed,
            4,
        )
        @time Pe = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, dim),
            ordering_as_matrix(:lex, dim),
            :fractal_combined,
            3,
        )
        @time Ge = groebnerwalk(
            I,
            ordering_as_matrix(:degrevlex, dim),
            ordering_as_matrix(:lex, dim),
            :generic,
        )


    #@time T6 = Singular.std(
    #        Singular.Ideal(S, [change_ring(x, S) for x in gens(I)]),
    #        complete_reduction = true,
    #    )

        T1 = Singular.Ideal(S, [change_ring(x, S) for x in gens(T)])
        T2 = Singular.Ideal(S, [change_ring(x, S) for x in gens(F)])
        T3 = Singular.Ideal(S, [change_ring(x, S) for x in gens(FA)])
        T4 = Singular.Ideal(S, [change_ring(x, S) for x in gens(St)])
        T5 = Singular.Ideal(S, [change_ring(x, S) for x in gens(Pe)])
        T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(Ge)])

        println("test tran: ", equalitytest(T6, T1))
        println("test fractal: ", equalitytest(T6, T2))
        println("test fractal: ", equalitytest(T6, T3))
        println("test pertubed: ", equalitytest(T5, T6))
        println("test standard: ", equalitytest(T4, T6))
        println("test generic: ", equalitytest(T6, T6))
        if !(
            equalitytest(T2, T1) &&
            equalitytest(T3, T4) &&
            equalitytest(T6, T5)
        )
            test_successfull = false
        end
    end
    println("All tests were: ", test_successfull)
end
