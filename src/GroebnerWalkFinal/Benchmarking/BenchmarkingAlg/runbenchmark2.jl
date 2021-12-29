using BenchmarkTools
function prepare2()

    df = DataFrame(
    example=["-"],
    standard=["-"],
    pertubed2=["-"],
    pertubed3=["-"],
    pertubed4=["-"],
    pertubed5=["-"],
    pertubed6=["-"],
    pertubed7=["-"],
    pertubed8=["-"],
    pertubed9=["-"],
    pertubed10=["-"],
    fractal=["-"],
    fractallex=["-"],
    fractallookahead=["-"],
    fractalcombined=["-"],
    generic=["-"],
    tran=["-"],
    stime=["-"],
    ttime=["-"]
    )
    savew(df, "CompareAlg")
end

function runb2(
    v::String,
    ideal::Singular.sideal,
    S::Singular.PolyRing,
    StartOrd::Matrix{Int64},
    TarOrd::Matrix{Int64},
)
    println("starting Benchmark2")


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

    println("Computing deglex-Basis")
    @belapsed stime = Singular.std($ideal, complete_reduction = true) evals = 1 samples = 1
    I = Singular.std(ideal, complete_reduction = true)
    ideals = []
    println("Benchmarking groebnerwalk2")

    for i = 2:nvars(S)
        if i ==2
        pertubed2 = @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =1 samples =1
    elseif i==3
        pertubed3 = @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =1 samples =1
    elseif i==4
        pertubed4 = @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =1 samples =1
    elseif i==5
        pertubed5 = @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =1 samples =1
    elseif i==6
        pertubed6 = @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =1 samples =1
    elseif i==7
        pertubed7 = @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =1 samples =1
    elseif i==8
        pertubed8 = @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =1 samples =1
    elseif i==9
        pertubed9 = @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =1 samples =1
    elseif i==10
        pertubed10 = @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:pertubed, $i) evals =1 samples =1
    end
        gb = groebnerwalk2(I, StartOrd, TarOrd, :pertubed, i)
        push!(ideals, gb)
    end
    standard = @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:standard, $i) evals =1 samples =1
    push!(ideals, groebnerwalk2(I, StartOrd, TarOrd, :standard))
    generic = @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:generic, $i) evals =1 samples =1
    push!(ideals, groebnerwalk2(I, StartOrd, TarOrd, :generic))
    fractal = @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:fractal, $i) evals =1 samples =1
    push!(ideals, groebnerwalk2(I, StartOrd, TarOrd, :fractal))
    fractallex = @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:fractallex, $i) evals =1 samples =1
    push!(ideals, groebnerwalk2(I, StartOrd, TarOrd, :fractal_look_ahead))
    fractallookahead = @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:fractallookahead, $i) evals =1 samples =1
    push!(ideals, groebnerwalk2(I, StartOrd, TarOrd, :fractal_lex))
    fractalcombined = @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:fractalcombined, $i) evals =1 samples =1
    push!(ideals, groebnerwalk2(I, StartOrd, TarOrd, :fractal_combined))
    tran = @belapsed groebnerwalk2($I, $StartOrd, $TarOrd, $:tran, $i) evals =1 samples =1
    push!(ideals, groebnerwalk2(I, StartOrd, TarOrd, :tran))


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
    savea(df, "CompareAlg")

end
