using BenchmarkTools
function prepare2()

    df = DataFrame(
    example=[],
    standard=[],
    pertubed2=[],
    pertubed3=[],
    pertubed4=[],
    pertubed5=[],
    pertubed6=[],
    pertubed7=[],
    pertubed8=[],
    pertubed9=[],
    pertubed10=[],
    fractal=[],
    fractallex=[],
    fractallookahead[],
    fractalcombined=[],
    generic=[],
    tran=[],
    stime=[],
    ttime=[]
    )
    savew(df, "CompareAlg.txt")
end

function runb2(
    v::String,
    ideal::Singular.sideal,
    S::Singular.PolyRing,
    StartOrd::Matrix{Int64},
    TarOrd::Matrix{Int64},
)
    println("starting Benchmark")


    example=[]
    standard=[]
    pertubed2=[]
    pertubed3=[]
    pertubed4=[]
    pertubed5=[]
    pertubed6=[]
    pertubed7=[]
    pertubed8=[]
    pertubed9=[]
    pertubed10=[]
    fractal=[]
    fractallex=[]
    fractallookahead[]
    fractalcombined=[]
    generic=[]
    tran=[]

    @belapsed stime = Singular.std($ideal, complete_reduction = true) evals = 1 samples = 1
    I = Singular.std(ideal, complete_reduction = true)
    ideals = []

    for i = 2:nvars(S)
        if i ==2
        pertubed2 = @belapsed groebnerwalk($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =1 samples =1
    elseif i==3
        pertubed3 = @belapsed groebnerwalk($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =1 samples =1
    elseif i==4
        pertubed4 = @belapsed groebnerwalk($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =1 samples =1
    elseif i==5
        pertubed5 = @belapsed groebnerwalk($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =1 samples =1
    elseif i==6
        pertubed6 = @belapsed groebnerwalk($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =1 samples =1
    elseif i==7
        pertubed7 = @belapsed groebnerwalk($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =1 samples =1
    elseif i==8
        pertubed8 = @belapsed groebnerwalk($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =1 samples =1
    elseif i==9
        pertubed9 = @belapsed groebnerwalk($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =1 samples =1
    elseif i==10
        pertubed10 = @belapsed groebnerwalk($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =1 samples =1
    end
        gb = groebnerwalk(I, StartOrd, TarOrd, :pertubed, i)
        push!(ideals, gb)
    end
    standard = @belapsed groebnerwalk($I, $StartOrd, $TarOrd, $:standard, $i) evals =1 samples =1
    push!(ideals, groebnerwalk(I, StartOrd, TarOrd, :standard))
    generic = @belapsed groebnerwalk($I, $StartOrd, $TarOrd, $:generic, $i) evals =1 samples =1
    push!(ideals, groebnerwalk(I, StartOrd, TarOrd, :generic))
    fractal = @belapsed groebnerwalk($I, $StartOrd, $TarOrd, $:fractal, $i) evals =1 samples =1
    push!(ideals, groebnerwalk(I, StartOrd, TarOrd, :fractal))
    fractallex = @belapsed groebnerwalk($I, $StartOrd, $TarOrd, $:fractallex, $i) evals =1 samples =1
    push!(ideals, groebnerwalk(I, StartOrd, TarOrd, :fractal_look_ahead))
    fractallookahead = @belapsed groebnerwalk($I, $StartOrd, $TarOrd, $:fractallookahead, $i) evals =1 samples =1
    push!(ideals, groebnerwalk(I, StartOrd, TarOrd, :fractal_lex))
    fractalcombined = @belapsed groebnerwalk($I, $StartOrd, $TarOrd, $:fractalcombined, $i) evals =1 samples =1
    push!(ideals, groebnerwalk(I, StartOrd, TarOrd, :fractal_combined))
    tran = @belapsed groebnerwalk($I, $StartOrd, $TarOrd, $:tran, $i) evals =1 samples =1
    push!(ideals, groebnerwalk(I, StartOrd, TarOrd, :tran))

    println("Computing GB")
    ttime = @belapsed Singular.std(
        Singular.Ideal($S, [change_ring(x, $S) for x in Singular.gens($ideal)]),
        complete_reduction = true,
    ) evals=1 samples =1

    df = DataFrame(startTime = [stime], targetTime =[ttime], example=[v])
    savea(df, "SingularComputationTimings")
    s = Singular.std(
        Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(ideal)]),
        complete_reduction = true,
    )

    df = DataFrame(test1 = ["-"], test2 = ["-"], example = ["-"])
    savea(df, "correct.txt")
    df = DataFrame(
    example=[v],
    standard=[standard],
    pertubed2=[pertubed2],
    pertubed3=[pertubed3],
    pertubed4=[pertubed4],
    pertubed5=[pertubed5],
    pertubed6=[pertubed6],
    pertubed7=[pertubed7],
    pertubed8=[pertubed8],
    pertubed9=[pertubed9],
    pertubed10=[pertubed10],
    fractal=[fractal],
    fractallex=[fractallex],
    fractallookahead[fractallookahead],
    fractalcombined=[fractalcombined],
    generic=[generic],
    tran=[tran],
    stime=[stime],
    ttime=[ttime]
    )
    savea(df, "CompareAlg.txt")

    for id in ideals
        a = isequal(
            Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
            s,
        )
        b = equalitytest(s, id)
        df = DataFrame(a = [a], b = [b],c=[v])
        savea(df, "correct.txt")
    end
end
