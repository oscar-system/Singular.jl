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


        @time T0 = Singular.std(
            Singular.Ideal(S, [change_ring(x, S) for x in gens(I)]),
            complete_reduction = true,
        )
        T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(H)])
        println("test tran: ", equalitytest(T0, T6))

        if !(equalitytest(T6, T0))
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
            :fractal,
            3,
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

        T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(H)])

        println("test tran: ", equalitytest(T0, T6))

        if !(equalitytest(T0, T6))
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
            :fractal,
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

        T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(H)])

        println("test tran: ", Singular.equal(T0, T6))
        if !(Singular.equal(T0, T6))
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
            :generic,
            3
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

        T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(H)])

        println("test tran: ", equalitytest(T0, T6))

        if !(equalitytest(T6, T0))
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
            :fractal,
            2,
        )

        T1 = Singular.Ideal(S, [change_ring(x, S) for x in gens(H)])

        @time T0 = Singular.std(
            Singular.Ideal(S, [change_ring(x, S) for x in gens(I)]),
            complete_reduction = true,
        )
        T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(H)])

        println("test tran: ", equalitytest(T0, T6))

        if !(Singular.equal(T0, T6))
            test_successfull = false
        end

    end



    if case == 6 || case == 99
        R, (v, u, t, z, y, x) = Singular.PolynomialRing(
            Singular.QQ,
            ["v", "u", "t", "z", "y", "x"],
            ordering = Singular.ordering_M(
                ordering_as_matrix([1, 1, 1, 1, 1, 1], :lex),
            ),
        )

        f1 = 2 * x^2 + 2 * y^2 + 2 * z^2 + 2 * t^2 + 2 * u^2 + v^2 - v
        f2 = 2*x * y +2* y * z + 2 * z * t + 2 * t * u + 2 * u * v - u
        f3 = 2 * x * z + 2 * y * t + 2 * z * u +2* u^2 + 2 * t * v - t
        f4 = 2 * x * t + 2 * y * u + 2 * t * u + 2 * z * v - z
        f5 = t^2 + 2 * x * v + 2 * y * v + 2 * z * v - y
        f6 = 2 * x + 2 * y + 2 * z + 2 * t + 2 * u + v - 1

        I = Singular.Ideal(R, [f1, f2, f3, f4, f5, f6])
        I = Singular.std(I, complete_reduction = true)

        @time H = groebnerwalk(
            I,
            ordering_as_matrix([1, 1, 1, 1, 1, 1], :lex),
            ordering_as_matrix(:lex, 6),
            :pertubed,
            5,
        )
        @time J = groebnerwalk(
            I,
            ordering_as_matrix([1, 1, 1, 1, 1, 1], :lex),
            ordering_as_matrix(:lex, 6),
            :generic,
            5,
        )

        S, V = Singular.PolynomialRing(
            Singular.QQ,
            ["v", "u", "t", "z", "y", "x"],
            ordering = Singular.ordering_M(ordering_as_matrix(:lex, 6)),
        )
        T0 = Singular.Ideal(S, [change_ring(x, S) for x in gens(J)])
        T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(H)])

        println("test standard: ", equalitytest(T0, T6))
        if !(isequal(T6, T0))
            test_successfull = false
        end
    end


if case == 7 || case == 99
    #Katsura 5
    R, (v, u, t, z, y) = Singular.PolynomialRing(
        Singular.QQ,
        ["v", "u", "t", "z", "y"],
        ordering = Singular.ordering_M(
            ordering_as_matrix([1, 1, 1, 1, 1], :lex),
        ),
    )

    f1 = 2 * y^2 + 2 * z^2 + 2 * t^2 + 2 * u^2 + v^2 - v
    f2 =  2*z*y + 2 * z * t + 2 * t * u + 2 * u * v - u
    f3 =  2 * y * t + 2 * z * u + 2*u^2 + 2 * t * v - t
    f4 =  2 * y * u + 2 * t * u + 2 * z * v - z
    f5 =  2 * y + 2 * z + 2 * t + 2 * u + v - 1

    I = Singular.Ideal(R, [f1, f2, f3, f4, f5])
    I = Singular.std(I, complete_reduction = true)

    @time H = groebnerwalk(
        I,
        ordering_as_matrix([1, 1, 1, 1, 1], :lex),
        ordering_as_matrix(:lex, 5),
        :fractal_lex,
        5,
    )
    @time J = groebnerwalk(
        I,
        ordering_as_matrix([1, 1, 1, 1, 1], :lex),
        ordering_as_matrix(:lex, 5),
        :fractal_look_ahead,
        5,
    )


    S, V = Singular.PolynomialRing(
        Singular.QQ,
        ["v", "u", "t", "z", "y"],
        ordering = Singular.ordering_M(ordering_as_matrix(:lex, 5)),
    )
    T0 = Singular.Ideal(S, [change_ring(x, S) for x in gens(J)])
    T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(H)])

    println("test standard: ", equalitytest(T0, T6))
    if !(isequal(T6, T0))
        test_successfull = false
    end
end

if case == 8 || case == 99

    #Katsura 7
    R, (x1,x2, x3, x4, x5, x6, x7, x8) = Singular.PolynomialRing(
        Singular.QQ,
        ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8"],
        ordering = Singular.ordering_M(
            ordering_as_matrix([1, 1, 1, 1, 1,1,1,1], :lex),
        ),
    )

    f1=-x1+2*x8^2+2*x7^2+2*x6^2+2*x5^2+2*x4^2+2*x3^2+2*x2^2+x1^2
     f2=-x2+2*x8*x7+2*x7*x6+2*x6*x5+2*x5*x4+2*x4*x3+2*x3*x2+2*x2*x1
    f3= -x3+2*x8*x6+2*x7*x5+2*x6*x4+2*x5*x3+2*x4*x2+2*x3*x1+x2^2
    f4= -x4+2*x8*x5+2*x7*x4+2*x6*x3+2*x5*x2+2*x4*x1+2*x3*x2
    f5= -x5+2*x8*x4+2*x7*x3+2*x6*x2+2*x5*x1+2*x4*x2+x3^2
     f6=-x6+2*x8*x3+2*x7*x2+2*x6*x1+2*x5*x2+2*x4*x3
     f7=-x7+2*x8*x2+2*x7*x1+2*x6*x2+2*x5*x3+x4^2
     f8=-1+2*x8+2*x7+2*x6+2*x5+2*x4+2*x3+2*x2+x1
    I = Singular.Ideal(R, [f1, f2, f3, f4, f5,f6,f7,f8])
    I = Singular.std(I, complete_reduction = true)

    @time H = groebnerwalk(
        I,
        ordering_as_matrix([1, 1, 1, 1, 1,1,1,1], :lex),
        ordering_as_matrix(:lex, 8),
        :generic,
        5,
    )
    @time J = groebnerwalk(
        I,
        ordering_as_matrix([1, 1, 1, 1, 1,1,1,1], :lex),
        ordering_as_matrix(:lex, 8),
        :generic,
        5,
    )


    S, V = Singular.PolynomialRing(
        Singular.QQ,
        ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8"],
        ordering = Singular.ordering_M(ordering_as_matrix(:lex, 8)),
    )
    T0 = Singular.Ideal(S, [change_ring(x, S) for x in gens(J)])
    T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(H)])

    println("test standard: ", equalitytest(T0, T6))
    if !(isequal(T6, T0))
        test_successfull = false
    end
end

if case == 9 || case == 99

    #eco10
    R, (x0,x1,x2, x3, x4, x5, x6, x7, x8, x9) = Singular.PolynomialRing(
        Singular.QQ,
        ["x0","x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9"],
        ordering = Singular.ordering_M(
            ordering_as_matrix([1, 1, 1, 1, 1,1,1,1,1,1], :lex),
        ),
    )

    f1=x0*x1*x9+x1*x2*x9+x2*x3*x9+x3*x4*x9+x4*x5*x9+x5*x6*x9+x6*x7*x9+x7*x8*x9+x0*x9-1
    f2=x0*x2*x9+x1*x3*x9+x2*x4*x9+x3*x5*x9+x4*x6*x9+x5*x7*x9+x6*x8*x9+x1*x9-2
    f3=x0*x3*x9+x1*x4*x9+x2*x5*x9+x3*x6*x9+x4*x7*x9+x5*x8*x9+x2*x9-3
    f4=x0*x4*x9+x1*x5*x9+x2*x6*x9+x3*x7*x9+x4*x8*x9+x3*x9-4
    f5=x0*x5*x9+x1*x6*x9+x2*x7*x9+x3*x8*x9+x4*x9-5
    f6=x0*x6*x9+x1*x7*x9+x2*x8*x9+x5*x9-6
    f7=x0*x7*x9+x1*x8*x9+x6*x9-7
    f8=x0*x8*x9+x7*x9-8
    f9=x8*x9-9
    f10=x0+x1+x2+x3+x4+x5+x6+x7+x8+1
    I = Singular.Ideal(R, [f1, f2, f3, f4, f5,f6,f7,f8,f9,f10])
    I = Singular.std(I, complete_reduction = true)

    @time H = groebnerwalk(
        I,
        ordering_as_matrix([1, 1, 1, 1, 1,1,1,1,1,1], :lex),
        ordering_as_matrix(:lex, 10),
        :generic,
        5,
    )
    @time J = groebnerwalk(
        I,
        ordering_as_matrix([1, 1, 1, 1, 1,1,1,1,1,1], :lex),
        ordering_as_matrix(:lex, 10),
        :generic,
        5,
    )


    S, V = Singular.PolynomialRing(
        Singular.QQ,
        ["x0","x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9"],
        ordering = Singular.ordering_M(ordering_as_matrix(:lex, 10)),
    )
    T0 = Singular.Ideal(S, [change_ring(x, S) for x in gens(J)])
    T6 = Singular.Ideal(S, [change_ring(x, S) for x in gens(H)])

    println("test standard: ", equalitytest(T0, T6))
    if !(isequal(T6, T0))
        test_successfull = false
    end
end
println("All tests were: ", test_successfull)
end
