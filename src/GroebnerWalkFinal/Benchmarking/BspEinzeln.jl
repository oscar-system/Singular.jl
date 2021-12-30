include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/Benchmarking/bechmarkingEveryProcedure/GroebnerWalkFinalBenchmarkProcedures.jl")
include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/Benchmarking/bechmarkingEveryProcedure/runbenchmark.jl")
include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/Benchmarking/BenchmarkingAlg/GroebnerWalkFinalBenchmark.jl")
include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/Benchmarking/BenchmarkingAlg/runbenchmark2.jl")
include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/Benchmarking/BenchmarkHelper")
include("readWriteHelper.jl")



using DataFrames
using CSV
function runAllSingleExample()
    cd("/Users/JordiWelp/Results")
    prepare()
    prepare2()
    prepareAlloc()

    #Katsura 5
    dim = 5
    ve = [1, 1, 1, 1, 1]
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (v, u, t, z, y) = Singular.PolynomialRing(
        Singular.QQ,
        ["v", "u", "t", "z", "y"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)

    f1 = 2 * y^2 + 2 * z^2 + 2 * t^2 + 2 * u^2 + v^2 - v
    f2 = 2 * z * y + 2 * z * t + 2 * t * u + 2 * u * v - u
    f3 = 2 * y * t + 2 * z * u + 2 * u^2 + 2 * t * v - t
    f4 = 2 * y * u + 2 * t * u + 2 * z * v - z
    f5 = 2 * y + 2 * z + 2 * t + 2 * u + v - 1
    I = Singular.Ideal(R, [f1, f2, f3, f4, f5])
    runAll("Katsura5", I, S, StartOrd, TarOrd)


    #Katsura6

    dim = 6
    ve = [1, 1, 1, 1, 1, 1]
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (v, u, t, z, y, x) = Singular.PolynomialRing(
        Singular.QQ,
        ["v", "u", "t", "z", "y", "x"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)
    f1 = 2 * x^2 + 2 * y^2 + 2 * z^2 + 2 * t^2 + 2 * u^2 + v^2 - v
    f2 = 2 * x * y + 2 * y * z + 2 * z * t + 2 * t * u + 2 * u * v - u
    f3 = 2 * x * z + 2 * y * t + 2 * z * u + 2 * u^2 + 2 * t * v - t
    f4 = 2 * x * t + 2 * y * u + 2 * t * u + 2 * z * v - z
    f5 = t^2 + 2 * x * v + 2 * y * v + 2 * z * v - y
    f6 = 2 * x + 2 * y + 2 * z + 2 * t + 2 * u + v - 1

    I = Singular.Ideal(R, [f1, f2, f3, f4, f5, f6])
    runAll("Katsura6", I, S, StartOrd, TarOrd)


    dim = 7
    ve = [1, 1, 1, 1, 1, 1, 1]
    example = "Cyclic7"
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (x1, x2, x3, x4, x5, x6, x7) = Singular.PolynomialRing(
        Singular.QQ,
        ["x1", "x2", "x3", "x4", "x5", "x6", "x7"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)
    f1 = x1 + x2 + x3 + x4 + x5 + x6 + x7
    f2 = x1 * x2 + x2 * x3 + x3 * x4 + x4 * x5 + x5 * x6 + x1 * x7 + x6 * x7
    f3 =
        x1 * x2 * x3 +
        x2 * x3 * x4 +
        x3 * x4 * x5 +
        x4 * x5 * x6 +
        x1 * x2 * x7 +
        x1 * x6 * x7 +
        x5 * x6 * x7
    f4 =
        x1 * x2 * x3 * x4 +
        x2 * x3 * x4 * x5 +
        x3 * x4 * x5 * x6 +
        x1 * x2 * x3 * x7 +
        x1 * x2 * x6 * x7 +
        x1 * x5 * x6 * x7 +
        x4 * x5 * x6 * x7
    f5 =
        x1 * x2 * x3 * x4 * x5 +
        x2 * x3 * x4 * x5 * x6 +
        x1 * x2 * x3 * x4 * x7 +
        x1 * x2 * x3 * x6 * x7 +
        x1 * x2 * x5 * x6 * x7 +
        x1 * x4 * x5 * x6 * x7 +
        x3 * x4 * x5 * x6 * x7
    f6 =
        x1 * x2 * x3 * x4 * x5 * x6 +
        x1 * x2 * x3 * x4 * x5 * x7 +
        x1 * x2 * x3 * x4 * x6 * x7 +
        x1 * x2 * x3 * x5 * x6 * x7 +
        x1 * x2 * x4 * x5 * x6 * x7 +
        x1 * x3 * x4 * x5 * x6 * x7 +
        x2 * x3 * x4 * x5 * x6 * x7
    f7 = x1 * x2 * x3 * x4 * x5 * x6 * x7 - 1

    I = Singular.Ideal(R, [f1, f2, f3, f4, f5, f6, f7])
    runAll(example, I, S, StartOrd, TarOrd)




    dim = 6
    ve = [1, 1, 1, 1, 1, 1]
    example = "Cyclic6"
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (x1, x2, x3, x4, x5, x6) = Singular.PolynomialRing(
        Singular.QQ,
        ["x1", "x2", "x3", "x4", "x5", "x6"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)


    f1 = x1 + x2 + x3 + x4 + x5 + x6
    f2 = x1 * x2 + x2 * x3 + x3 * x4 + x4 * x5 + x1 * x6 + x5 * x6
    f3 =
        x1 * x2 * x3 +
        x2 * x3 * x4 +
        x3 * x4 * x5 +
        x1 * x2 * x6 +
        x1 * x5 * x6 +
        x4 * x5 * x6
    f4 =
        x1 * x2 * x3 * x4 +
        x2 * x3 * x4 * x5 +
        x1 * x2 * x3 * x6 +
        x1 * x2 * x5 * x6 +
        x1 * x4 * x5 * x6 +
        x3 * x4 * x5 * x6
    f5 =
        x1 * x2 * x3 * x4 * x5 +
        x1 * x2 * x3 * x4 * x6 +
        x1 * x2 * x3 * x5 * x6 +
        x1 * x2 * x4 * x5 * x6 +
        x1 * x3 * x4 * x5 * x6 +
        x2 * x3 * x4 * x5 * x6
    f6 = x1 * x2 * x3 * x4 * x5 * x6 - 1

    I = Singular.Ideal(R, [f1, f2, f3, f4, f5, f6])
    runAll(example, I, S, StartOrd, TarOrd)

    dim = 5
    ve = [1, 1, 1, 1, 1]
    example = "Cyclic5"
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (v, w, x, y, z) = Singular.PolynomialRing(
        Singular.QQ,
        ["v", "w", "x", "y", "z"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)

    f1 = v + w + x + y + z
    f2 = v * w + w * x + x * y + v * z + y * z
    f3 = v * w * x + w * x * y + v * w * z + v * y * z + x * y * z
    f4 =
        v * w * x * y +
        v * w * x * z +
        v * w * y * z +
        v * x * y * z +
        w * x * y * z
    f5 = v * w * x * y * z - 1


    I = Singular.Ideal(R, [f1, f2, f3, f4, f5])
    runAll(example, I, S, StartOrd, TarOrd)


    dim = 6
    ve = [1, 1, 1, 1, 1, 1]
    example = "eco6"
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (x1, x2, x3, x4, x5, x6) = Singular.PolynomialRing(
        Singular.QQ,
        ["x1", "x2", "x3", "x4", "x5", "x6"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)


    f1 = x1 + x2 + x3 + x4 + x5 + 1
    f2 = x5 * x6 - 5
    f3 = x1 * x5 * x6 + x4 * x6 - 4
    f4 = x1 * x4 * x6 + x2 * x5 * x6 + x3 * x6 - 3
    f5 = x1 * x3 * x6 + x2 * x4 * x6 + x3 * x5 * x6 + x2 * x6 - 2
    f6 = x1 * x2 * x6 + x2 * x3 * x6 + x3 * x4 * x6 + x4 * x5 * x6 + x1 * x6 - 1

    I = Singular.Ideal(R, [f1, f2, f3, f4, f5, f6])
    runAll(example, I, S, StartOrd, TarOrd)


    dim = 6
    ve = [1, 1, 1, 1, 1, 1]
    example = "eco7"
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (x1, x2, x3, x4, x5, x6) = Singular.PolynomialRing(
        Singular.QQ,
        ["x1", "x2", "x3", "x4", "x5", "x6"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)

    f1 = x1 + x2 + x3 + x4 + x5 + x6 + 1
    f2 = x6 * x7 - 6
    f3 = x1 * x6 * x7 + x5 * x7 - 5
    f4 = x1 * x5 * x7 + x2 * x6 * x7 + x4 * x7 - 4
    f5 = x1 * x4 * x7 + x2 * x5 * x7 + x3 * x6 * x7 + x3 * x7 - 3
    f6 = x1 * x3 * x7 + x2 * x4 * x7 + x3 * x5 * x7 + x4 * x6 * x7 + x2 * x7 - 2
    f7 =
        x1 * x2 * x7 +
        x2 * x3 * x7 +
        x3 * x4 * x7 +
        x4 * x5 * x7 +
        x5 * x6 * x7 +
        x1 * x7 - 1


    I = Singular.Ideal(R, [f1, f2, f3, f4, f5, f6, f7])
    runAll(example, I, S, StartOrd, TarOrd)



    dim = 5
    ve = [1, 1, 1, 1, 1]
    example = "noon5"
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (x1, x2, x3, x4, x5) = Singular.PolynomialRing(
        Singular.QQ,
        ["x1", "x2", "x3", "x4", "x5"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)

    f1 =
        10 * x1^2 * x5 + 10 * x2^2 * x5 + 10 * x3^2 * x5 + 10 * x4^2 * x5 -
        11 * x5 + 10
    f2 =
        10 * x1^2 * x4 + 10 * x2^2 * x4 + 10 * x3^2 * x4 + 10 * x4 * x5^2 -
        11 * x4 + 10
    f3 =
        10 * x1^2 * x3 + 10 * x2^2 * x3 + 10 * x3 * x4^2 + 10 * x3 * x5^2 -
        11 * x3 + 10
    f4 =
        10 * x1 * x2^2 + 10 * x1 * x3^2 + 10 * x1 * x4^2 + 10 * x1 * x5^2 -
        11 * x1 + 10
    f5 =
        10 * x1^2 * x2 + 10 * x2 * x3^2 + 10 * x2 * x4^2 + 10 * x2 * x5^2 -
        11 * x2 + 10

    I = Singular.Ideal(R, [f1, f2, f3, f4, f5])
    runAll(example, I, S, StartOrd, TarOrd)


    dim = 6
    ve = [1, 1, 1, 1, 1, 1]
    example = "noon6"
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (x1, x2, x3, x4, x5, x6) = Singular.PolynomialRing(
        Singular.QQ,
        ["x1", "x2", "x3", "x4", "x5", "x6"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)

    f1 =
        10 * x1^2 * x6 +
        10 * x2^2 * x6 +
        10 * x3^2 * x6 +
        10 * x4^2 * x6 +
        10 * x5^2 * x6 - 11 * x6 + 10
    f2 =
        10 * x1^2 * x5 +
        10 * x2^2 * x5 +
        10 * x3^2 * x5 +
        10 * x4^2 * x5 +
        10 * x5 * x6^2 - 11 * x5 + 10
    f3 =
        10 * x1^2 * x4 +
        10 * x2^2 * x4 +
        10 * x3^2 * x4 +
        10 * x4 * x5^2 +
        10 * x4 * x6^2 - 11 * x4 + 10
    f4 =
        10 * x1^2 * x3 +
        10 * x2^2 * x3 +
        10 * x3 * x4^2 +
        10 * x3 * x5^2 +
        10 * x3 * x6^2 - 11 * x3 + 10
    f5 =
        10 * x1 * x2^2 +
        10 * x1 * x3^2 +
        10 * x1 * x4^2 +
        10 * x1 * x5^2 +
        10 * x1 * x6^2 - 11 * x1 + 10
    f6 =
        10 * x1^2 * x2 +
        10 * x2 * x3^2 +
        10 * x2 * x4^2 +
        10 * x2 * x5^2 +
        10 * x2 * x6^2 - 11 * x2 + 10

    I = Singular.Ideal(R, [f1, f2, f3, f4, f5, f6])
    runAll(example, I, S, StartOrd, TarOrd)


    dim = 7
    ve = [1, 1, 1, 1, 1, 1, 1]
    example = "noon7"
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (x1, x2, x3, x4, x5, x6, x7) = Singular.PolynomialRing(
        Singular.QQ,
        ["x1", "x2", "x3", "x4", "x5", "x6", "x7"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)

    f1 =
        10 * x1^2 * x7 +
        10 * x2^2 * x7 +
        10 * x3^2 * x7 +
        10 * x4^2 * x7 +
        10 * x5^2 * x7 +
        10 * x6^2 * x7 - 11 * x7 + 10
    f2 =
        10 * x1^2 * x6 +
        10 * x2^2 * x6 +
        10 * x3^2 * x6 +
        10 * x4^2 * x6 +
        10 * x5^2 * x6 +
        10 * x6 * x7^2 - 11 * x6 + 10
    f3 =
        10 * x1^2 * x5 +
        10 * x2^2 * x5 +
        10 * x3^2 * x5 +
        10 * x4^2 * x5 +
        10 * x5 * x6^2 +
        10 * x5 * x7^2 - 11 * x5 + 10
    f4 =
        10 * x1^2 * x4 +
        10 * x2^2 * x4 +
        10 * x3^2 * x4 +
        10 * x4 * x5^2 +
        10 * x4 * x6^2 +
        10 * x4 * x7^2 - 11 * x4 + 10
    f5 =
        10 * x1^2 * x3 +
        10 * x2^2 * x3 +
        10 * x3 * x4^2 +
        10 * x3 * x5^2 +
        10 * x3 * x6^2 +
        10 * x3 * x7^2 - 11 * x3 + 10
    f6 =
        10 * x1 * x2^2 +
        10 * x1 * x3^2 +
        10 * x1 * x4^2 +
        10 * x1 * x5^2 +
        10 * x1 * x6^2 +
        10 * x1 * x7^2 - 11 * x1 + 10
    f7 =
        10 * x1^2 * x2 +
        10 * x2 * x3^2 +
        10 * x2 * x4^2 +
        10 * x2 * x5^2 +
        10 * x2 * x6^2 +
        10 * x2 * x7^2 - 11 * x2 + 10

    I = Singular.Ideal(R, [f1, f2, f3, f4, f5, f6, f7])
    runAll(example, I, S, StartOrd, TarOrd)



    dim = 8
    ve = [1, 1, 1, 1, 1, 1, 1, 1]
    example = "noon8"
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (x1, x2, x3, x4, x5, x6, x7, x8) = Singular.PolynomialRing(
        Singular.QQ,
        ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)

    f1 =
        10 * x1^2 * x8 +
        10 * x2^2 * x8 +
        10 * x3^2 * x8 +
        10 * x4^2 * x8 +
        10 * x5^2 * x8 +
        10 * x6^2 * x8 +
        10 * x7^2 * x8 - 11 * x8 + 10
    f2 =
        10 * x1^2 * x7 +
        10 * x2^2 * x7 +
        10 * x3^2 * x7 +
        10 * x4^2 * x7 +
        10 * x5^2 * x7 +
        10 * x6^2 * x7 +
        10 * x7 * x8^2 - 11 * x7 + 10
    f3 =
        10 * x1^2 * x6 +
        10 * x2^2 * x6 +
        10 * x3^2 * x6 +
        10 * x4^2 * x6 +
        10 * x5^2 * x6 +
        10 * x6 * x7^2 +
        10 * x6 * x8^2 - 11 * x6 + 10
    f4 =
        10 * x1^2 * x5 +
        10 * x2^2 * x5 +
        10 * x3^2 * x5 +
        10 * x4^2 * x5 +
        10 * x5 * x6^2 +
        10 * x5 * x7^2 +
        10 * x5 * x8^2 - 11 * x5 + 10
    f5 =
        10 * x1^2 * x4 +
        10 * x2^2 * x4 +
        10 * x3^2 * x4 +
        10 * x4 * x5^2 +
        10 * x4 * x6^2 +
        10 * x4 * x7^2 +
        10 * x4 * x8^2 - 11 * x4 + 10
    f6 =
        10 * x1^2 * x3 +
        10 * x2^2 * x3 +
        10 * x3 * x4^2 +
        10 * x3 * x5^2 +
        10 * x3 * x6^2 +
        10 * x3 * x7^2 +
        10 * x3 * x8^2 - 11 * x3 + 10
    f7 =
        10 * x1 * x2^2 +
        10 * x1 * x3^2 +
        10 * x1 * x4^2 +
        10 * x1 * x5^2 +
        10 * x1 * x6^2 +
        10 * x1 * x7^2 +
        10 * x1 * x8^2 - 11 * x1 + 10
    f8 =
        10 * x1^2 * x2 +
        10 * x2 * x3^2 +
        10 * x2 * x4^2 +
        10 * x2 * x5^2 +
        10 * x2 * x6^2 +
        10 * x2 * x7^2 +
        10 * x2 * x8^2 - 11 * x2 + 10

    I = Singular.Ideal(R, [f1, f2, f3, f4, f5, f6, f7, f8])
    runAll(example, I, S, StartOrd, TarOrd)


    dim = 7
    ve = [1, 1, 1, 1, 1, 1, 1]
    example = "redeco7"
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (x1, x2, x3, x4, u7, x5, x6) = Singular.PolynomialRing(
        Singular.QQ,
        ["x1", "x2", "x3", "x4", "u7", "x5", "x6"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)

    f1 = -6 * u7 + x6
    f2 = x1 + x2 + x3 + x4 + x5 + x6 + 1
    f3 = x1 * x6 - 5 * u7 + x5
    f4 = x1 * x5 + x2 * x6 + x4 - 4 * u7
    f5 = x1 * x4 + x2 * x5 + x3 * x6 + x3 - 3 * u7
    f6 = x1 * x3 + x2 * x4 + x3 * x5 + x4 * x6 + x2 - 2 * u7
    f7 = x1 * x2 + x2 * x3 + x3 * x4 + x4 * x5 + x5 * x6 + x1 - u7

    I = Singular.Ideal(R, [f1, f2, f3, f4, f5, f6, f7])
    runAll(example, I, S, StartOrd, TarOrd)

    dim = 8
    ve = [1, 1, 1, 1, 1,1,1,1]
    example ="redeco8"
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    R, (x1,x2,x3,x4,u8,x5,x6,x7) = Singular.PolynomialRing(
        Singular.QQ,
        ["x1","x2","x3","x4","u8","x5","x6","x7"],
        ordering = Singular.ordering_M(StartOrd),
    )
    S = change_order(R, TarOrd)

      f1 =-7*u8+x7
      f2 =x1+x2+x3+x4+x5+x6+x7+1
      f3 =x1*x7-6*u8+x6
      f4 =x1*x6+x2*x7+x5-5*u8
      f5 =x1*x5+x2*x6+x3*x7+x4-4*u8
      f6 =x1*x4+x2*x5+x3*x6+x4*x7+x3-3*u8
      f7 =x1*x3+x2*x4+x3*x5+x4*x6+x5*x7+x2-2*u8
      f8 =x1*x2+x2*x3+x3*x4+x4*x5+x5*x6+x6*x7+x1-u8

      I = Singular.Ideal(R, [f1, f2,f3,f4,f5,f6,f7, f8])
      runAll(example, I, S, StartOrd, TarOrd)


      dim = 6
      ve = [1, 1, 1, 1, 1,1]
      example ="Wang-91"
      StartOrd = ordering_as_matrix(ve, :lex)
      TarOrd = ordering_as_matrix(:lex, dim)
      R, (x3,x2,x1,x0,b,a) = Singular.PolynomialRing(
          Singular.QQ,
          ["x3","x2","x1","x0","b","a"],
          ordering = Singular.ordering_M(StartOrd),
      )
      S = change_order(R, TarOrd)

  f1 =3*x2*x1*a+3*x0^2
  f2 =3*x2*x1*b+3*x3^2
  f3 =3*x3*x1*b+3*x1*x0*a+3*x2^2
  f4 =3*x3*x2*b+3*x2*x0*a+3*x1^2

  I = Singular.Ideal(R, [f1, f2,f3,f4])
  runAll(example, I, S, StartOrd, TarOrd)



  dim = 4
  ve = [1, 1, 1, 1]
  example ="cohn4"
  StartOrd = ordering_as_matrix(ve, :lex)
  TarOrd = ordering_as_matrix(:lex, dim)
  R, (x,y,z,t) = Singular.PolynomialRing(
      Singular.QQ,
      ["x","y","z","t"],
      ordering = Singular.ordering_M(StartOrd),
  )
  S = change_order(R, TarOrd)

  f1 =-x^3*y^2+2*x^2*y^2*z-x^2*y*z^2-144*x^2*y^2-207*x^2*y*z+288*x*y^2*z+78*x*y*z^2+x*z^3-3456*x^2*y-5184*x*y^2-9504*x*y*z-432*x*z^2-248832*x*y+62208*x*z-2985984*x
  f2 =y^3*t^3-y^2*z*t^3+4*y^3*t^2-2*y^2*z*t^2+72*y^2*t^3+71*y*z*t^3+288*y^2*t^2+360*y*z*t^2+6*z^2*t^2+1728*y*t^3-464*z*t^3+432*y*z*t+8*z^2*t+6912*y*t^2-4320*z*t^2+13824*t^3+z^2-13824*z*t+55296*t^2-13824*z
  f3 =x^2*y*t^3-2*x*y^2*t^3+y^3*t^3+8*x^2*y*t^2-12*x*y^2*t^2+4*y^3*t^2-24*x*y*t^3+24*y^2*t^3+20*x^2*y*t-20*x*y^2*t-160*x*y*t^2+96*y^2*t^2+128*x*t^3+16*x^2*y+96*x*y*t+2304*x*t^2+1152*x*y+13824*x*t+27648*x
  f4 =-x^3*z*t^2+x^2*z^2*t^2-6*x^3*z*t+4*x^2*z^2*t+32*x^3*t^2-72*x^2*z*t^2-87*x*z^2*t^2-z^3*t^2-8*x^3*z-432*x^2*z*t-414*x*z^2*t+2592*x*z*t^2+864*z^2*t^2-1728*x^2*z-20736*x*z*t+3456*z^2*t-186624*z*t^2-124416*x*z-1492992*z*t-2985984*z
  I = Singular.Ideal(R, [f1, f2,f3,f4])
  runAll(example, I, S, StartOrd, TarOrd)



  dim = 3
  ve = [9, 9, 8]
  example ="Oberfranz"
  StartOrd = ordering_as_matrix(ve, :lex)
  TarOrd = ordering_as_matrix(:lex, dim)
  R, (x1,x2,x3) = Singular.PolynomialRing(
      Singular.QQ,
      ["x1","x2","x3"],
      ordering = Singular.ordering_M(StartOrd),
  )
  S = change_order(R, TarOrd)
  f1=x2^3+x1*x2*x3+x2^2*x3+x1*x3^3
  f2=3+x1*x2+x1^2*x2+x2^2*x3

  I = Singular.Ideal(R, [f1, f2])
  runAll(example, I, S, StartOrd, TarOrd)


end
#=

dim = 7
ve = [1, 1, 1, 1, 1,1,1]
example ="redeco7"
StartOrd = ordering_as_matrix(ve, :lex)
TarOrd = ordering_as_matrix(:lex, dim)
R, (x1,x2,x3,x4,u7,x5,x6) = Singular.PolynomialRing(
    Singular.QQ,
    ["x1","x2","x3","x4","u7","x5","x6"],
    ordering = Singular.ordering_M(StartOrd),
)
S = change_order(R, TarOrd)

  f1 =-6*u7+x6
  f2 =x1+x2+x3+x4+x5+x6+1
  f3 =x1*x6-5*u7+x5
  f4 =x1*x5+x2*x6+x4-4*u7
  f5 =x1*x4+x2*x5+x3*x6+x3-3*u7
  f6 =x1*x3+x2*x4+x3*x5+x4*x6+x2-2*u7
  f7 =x1*x2+x2*x3+x3*x4+x4*x5+x5*x6+x1-u7

I = Singular.Ideal(R, [f1, f2,f3,f4,f5,f6,f7])
runb(example, I, S, StartOrd, TarOrd)

s = Singular.std(
    Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(I)]),
    complete_reduction = true,
)
df = DataFrame(a = ["test1"], b = ["test2"],c=["example"])
savew(df, "correct.txt")
for id in ideals
    a = isequal(
        Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
        s,
    )
    b = equalitytest(s, id)
    df = DataFrame(a = [a], b = [b], c=[example])
    savea(df, "correct.txt")
end
=#
