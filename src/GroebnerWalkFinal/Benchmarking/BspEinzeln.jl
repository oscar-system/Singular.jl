include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/Benchmarking/bechmarkingEveryProcedure/GroebnerWalkFinalBenchmarkProcedures.jl")
include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/Benchmarking/bechmarkingEveryProcedure/runbenchmark.jl")
include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/Benchmarking/BenchmarkingAlg/GroebnerWalkFinalBenchmark.jl")
include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/Benchmarking/BenchmarkingAlg/runbenchmark2.jl")
include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/Benchmarking/BenchmarkHelper")
include("readWriteHelper.jl")
include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/Examples")



using DataFrames
using CSV
function runAllSingleExample()
    cd("/Users/JordiWelp/Results")
    prepare()
    prepare2()
    prepareAlloc()

    id = katsura5()
    dim = nvars(base_ring(id))
    ve=ones(Int64,dim)
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    S = change_order(R, TarOrd)
    runAll("Katsura5", I, S, StartOrd, TarOrd)


    #Katsura6

    id = katsura6()
    dim = nvars(base_ring(id))
    ve=ones(Int64,dim)
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    S = change_order(R, TarOrd)
    runAll("Katsura6", I, S, StartOrd, TarOrd)



    example = "Cyclic7"
    id = cyclic7()
    dim = nvars(base_ring(id))
    ve=ones(Int64,dim)
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    S = change_order(R, TarOrd)
    runAll(example, I, S, StartOrd, TarOrd)




    example = "Cyclic6"
    id = cyclic6()
    dim = nvars(base_ring(id))
    ve=ones(Int64,dim)
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    S = change_order(R, TarOrd)
    runAll(example, I, S, StartOrd, TarOrd)


    example = "Cyclic5"
    id = cyclic5()
    dim = nvars(base_ring(id))
    ve=ones(Int64,dim)
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    S = change_order(R, TarOrd)
    runAll(example, I, S, StartOrd, TarOrd)


    example = "eco6"
    id = eco6()
    dim = nvars(base_ring(id))
    ve=ones(Int64,dim)
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    S = change_order(R, TarOrd)
    runAll(example, I, S, StartOrd, TarOrd)

    example = "eco7"
    id = eco7()
    dim = nvars(base_ring(id))
    ve=ones(Int64,dim)
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    S = change_order(R, TarOrd)
    runAll(example, I, S, StartOrd, TarOrd)

    example = "noon5"
    id = noon5()
    dim = nvars(base_ring(id))
    ve=ones(Int64,dim)
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    S = change_order(R, TarOrd)
    runAll(example, I, S, StartOrd, TarOrd)



    example = "noon6"
    id = noon6()
    dim = nvars(base_ring(id))
    ve=ones(Int64,dim)
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    S = change_order(R, TarOrd)
    runAll(example, I, S, StartOrd, TarOrd)

    example = "noon7"
    id = noon7()
    dim = nvars(base_ring(id))
    ve=ones(Int64,dim)
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    S = change_order(R, TarOrd)
    runAll(example, I, S, StartOrd, TarOrd)

    example = "noon8"
    id = noon8()
    dim = nvars(base_ring(id))
    ve=ones(Int64,dim)
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    S = change_order(R, TarOrd)
    runAll(example, I, S, StartOrd, TarOrd)

    example = "redeco7"
    id = redeco7()
    dim = nvars(base_ring(id))
    ve=ones(Int64,dim)
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    S = change_order(R, TarOrd)
    runAll(example, I, S, StartOrd, TarOrd)

    example = "redeco8"
    id = redeco8()
    dim = nvars(base_ring(id))
    ve=ones(Int64,dim)
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    S = change_order(R, TarOrd)
    runAll(example, I, S, StartOrd, TarOrd)

    example = "wang91"
    id = wang91()
    dim = nvars(base_ring(id))
    ve=ones(Int64,dim)
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    S = change_order(R, TarOrd)
    runAll(example, I, S, StartOrd, TarOrd)

    example = "cohn4"
    id = cohn4()
    dim = nvars(base_ring(id))
    ve=ones(Int64,dim)
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    S = change_order(R, TarOrd)
    runAll(example, I, S, StartOrd, TarOrd)

    example = "of"
    id = oberfr()
    dim = nvars(base_ring(id))
    ve=ones(Int64,dim)
    StartOrd = ordering_as_matrix(ve, :lex)
    TarOrd = ordering_as_matrix(:lex, dim)
    S = change_order(R, TarOrd)
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
