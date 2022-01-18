using BenchmarkTools


function runb(
    v::String,
    ideal::Singular.sideal,
    S::Singular.PolyRing,
    StartOrd::Matrix{Int64},
    TarOrd::Matrix{Int64},
)
    println("starting Benchmark")
    prepareExampleAlloc(v)
    prepareExampleElapsed(v)

    println("Computing deglex-Basis")
    I = Singular.std(ideal, complete_reduction = true)

    println("Benchmarking GroebnerWalk")
    ideals = []
    for i = 2:nvars(S)
        push!(ideals, groebnerwalk(I, StartOrd, TarOrd, :pertubed, i))
    end

    push!(ideals, groebnerwalk(I, StartOrd, TarOrd, :standard))
    push!(ideals, groebnerwalk(I, StartOrd, TarOrd, :generic))
    push!(ideals, groebnerwalk(I, StartOrd, TarOrd, :fractal))
    push!(ideals, groebnerwalk(I, StartOrd, TarOrd, :fractal_look_ahead))
    push!(ideals, groebnerwalk(I, StartOrd, TarOrd, :fractal_lex))
    push!(ideals, groebnerwalk(I, StartOrd, TarOrd, :fractal_combined))
    push!(ideals, groebnerwalk(I, StartOrd, TarOrd, :tran))


    s = Singular.std(
        Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(ideal)]),
        complete_reduction = true,
    )

    df = DataFrame(test1 = ["-"], test2 = ["-"], example = ["-"])
    savea(df, "correct")

    println("Benchmarking ideals")
    for id in ideals
        a = Singular.isequal(
            Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
            s,
        )
        b = "-"
        df = DataFrame(a = [a], b = [b], c = [v])
        savea(df, "correct")
    end
end
