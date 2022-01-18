using BenchmarkTools


function runb2(
    v::String,
    ideal::Singular.sideal,
    S::Singular.PolyRing,
    StartOrd::Matrix{Int64},
    TarOrd::Matrix{Int64},
)
    println("starting Benchmark2")


    println("Computing deglex-Basis")
    stime =
        @belapsed Singular.std($ideal, complete_reduction = true) evals = 5 samples =
            1
    I = Singular.std(ideal, complete_reduction = true)
    ideals = []
    println("Benchmarking groebnerwalk2")
    pertubed2 = "-"
    pertubed3 = "-"
    pertubed4 = "-"
    pertubed5 = "-"
    pertubed6 = "-"
    pertubed7 = "-"
    pertubed8 = "-"
    pertubed9 = "-"
    pertubed10 = "-"
    for i = 2:nvars(S)
        if i == 2
            pertubed2 =
                @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =
                    1 samples = 1
        elseif i == 3
            pertubed3 =
                @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =
                    1 samples = 1
        elseif i == 4
            pertubed4 =
                @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =
                    1 samples = 1
        elseif i == 5
            pertubed5 =
                @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =
                    1 samples = 1
        elseif i == 6
            pertubed6 =
                @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =
                    1 samples = 1
        elseif i == 7
            pertubed7 =
                @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =
                    1 samples = 1
        elseif i == 8
            pertubed8 =
                @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =
                    1 samples = 1
        elseif i == 9
            pertubed9 =
                @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =
                    1 samples = 1
        elseif i == 10
            pertubed10 =
                @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =
                    1 samples = 1
        end
    end
    standard =
        @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:standard) evals = 1 samples =
            1
    generic =
        @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:generic) evals = 1 samples =
            1
    fractal =
        @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:fractal) evals = 1 samples =
            1
    fractallex =
        @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:fractal_lex) evals = 1 samples =
            1
    fractallookahead =
        @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:fractal_look_ahead) evals =
            1 samples = 1
    fractalcombined =
        @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:fractal_combined) evals =
            1 samples = 1
    tran =
        @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:tran) evals = 1 samples =
            1

    println("Computing lex-Basis")
    ttime = @belapsed Singular.std(
        Singular.Ideal($S, [change_ring(x, $S) for x in Singular.gens($ideal)]),
        complete_reduction = true,
    ) evals = 5 samples = 1

    df = DataFrame(
        example = [v],
        standard = [standard],
        pertubed2 = [pertubed2],
        pertubed3 = [pertubed3],
        pertubed4 = [pertubed4],
        pertubed5 = [pertubed5],
        pertubed6 = [pertubed6],
        pertubed7 = [pertubed7],
        pertubed8 = [pertubed8],
        pertubed9 = [pertubed9],
        pertubed10 = [pertubed10],
        fractal = [fractal],
        fractallex = [fractallex],
        fractallookahead = [fractallookahead],
        fractalcombined = [fractalcombined],
        generic = [generic],
        tran = [tran],
        stime = [stime],
        ttime = [ttime],
    )
    savea(df, "CompareAlg")

end
