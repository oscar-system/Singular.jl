using BenchmarkTools

function prepare3()
df = DataFrame(
    weights = ["weights"],
    weight = ["weight"],
    nGens = ["nGens"],
    initials = ["initials"],
    stdh = ["stdh"],
    liftGW2 = ["liftGW2"],
    lift2 = ["lift2"],
    interred = ["interred"],
    example = ["example"],
)
savew(df, "allocsStandardWalk")
df = DataFrame(
    weights = ["weights"],
    weight = ["weight"],
    nGens = ["nGens"],
    initials = ["initials"],
    stdh = ["stdh"],
    liftGW2 = ["liftGW2"],
    lift2 = ["lift2"],
    interred = ["interred"],
    pert = ["pert"],
    inCone = ["inCone"],
    laststd = ["laststd"],
    degree = ["degree"],
    example = ["example"],
)
for i in 2:10
savew(df, "allocsPertubedWalk",i)
end

df = DataFrame(
    weights = ["weights"],
    facetnormal = ["-"],
    initials = ["initials"],
    stdh = ["stdh"],
    liftgeneric = ["liftgeneric"],
    interred = ["interred"],
    example = ["example"],
)
savew(df, "allocsGenericWalk")

df = DataFrame(
    weights = ["weights"],
    weight = ["weight"],
    nGens = ["nGens"],
    initials = ["initials"],
    stdh = ["stdh"],
    lift2 = ["lift2"],
    lift = ["lift"],
    interred = ["interred"],
    inCone = ["inCone"],
    pertvec = ["pertvec"],
    depth = ["depth"],
    example = ["example"],
)
savew(df, "allocsFractalWalk")
df = DataFrame(
    weights = ["weights"],
    weight = ["weight"],
    nGens = ["nGens"],
    initials = ["initials"],
    stdh = ["stdh"],
    lift2 = ["lift2"],
    lift = ["lift"],
    interred = ["interred"],
    inCone = ["inCone"],
    pertvec = ["pertvec"],
    depth = ["depth"],
    example = ["example"],
)
savew(df, "allocsFractalWalklex")
df = DataFrame(
    weights = ["weights"],
    weight = ["weight"],
    nGens = ["nGens"],
    initials = ["initials"],
    stdh = ["stdh"],
    lift2 = ["lift2"],
    lift = ["lift"],
    interred = ["interred"],
    inCone = ["inCone"],
    pertvec = ["pertvec"],
    depth = ["depth"],
    example = ["example"],
)
savew(df, "allocsFractalWalklookahead")
df = DataFrame(
    weights = ["weights"],
    weight = ["weight"],
    nGens = ["nGens"],
    initials = ["initials"],
    stdh = ["stdh"],
    lift2 = ["lift2"],
    lift = ["lift"],
    interred = ["interred"],
    inCone = ["inCone"],
    pertvec = ["pertvec"],
    depth = ["depth"],
    example = ["example"],
)
savew(df, "allocsFractalWalkcombined")
df = DataFrame(
    weights = "weights",
    weight = "weight",
    nGens = ["nGens"],
    initials = ["initials"],
    stdh = ["stdh"],
    liftGW2 = ["liftGW2"],
    lift2 = ["lift2"],
    interred = ["interred"],
    rep = ["rep"],
    inCone = ["inCone"],
    inseveral = ["inseveral"],
    example = ["example"],
)
savew(df, "allocsTranWalk")
end
function runb3(
    v::String,
    ideal::Singular.sideal,
    S::Singular.PolyRing,
    StartOrd::Matrix{Int64},
    TarOrd::Matrix{Int64},
)
    println("starting Benchmark")
    df = DataFrame(
        weights = ["-"],
        weight = ["-"],
        nGens = ["-"],
        initials = ["-"],
        stdh = ["-"],
        liftGW2 = ["-"],
        lift2 = ["-"],
        interred = ["-"],
        example = [v],
    )
    savea(df, "allocsStandardWalk")
    df = DataFrame(
        weights = ["-"],
        weight = ["-"],
        nGens = ["-"],
        initials = ["-"],
        stdh = ["-"],
        liftGW2 = ["-"],
        lift2 = ["-"],
        interred = ["-"],
        pert = ["-"],
        inCone = ["-"],
        laststd = ["-"],
        degree = ["-"],
        example = [v],
    )
    for i in 2:10
    savea(df, "allocsPertubedWalk",i)
end

    df = DataFrame(
        weights = ["-"],
        facetnormal = ["-"],
        initials = ["-"],
        stdh = ["-"],
        liftgeneric = ["-"],
        interred = ["-"],
        example = [v],
    )
    savea(df, "allocsGenericWalk")

    df = DataFrame(
        weights = ["-"],
        weight = ["-"],
        nGens = ["-"],
        initials = ["-"],
        stdh = ["-"],
        lift2 = ["-"],
        lift = ["-"],
        interred = ["-"],
        inCone = ["-"],
        pertvec = ["-"],
        depth = ["-"],
        example = [v],
    )
    savea(df, "allocsFractalWalk")
    df = DataFrame(
        weights = ["-"],
        weight = ["-"],
        nGens = ["-"],
        initials = ["-"],
        stdh = ["-"],
        lift2 = ["-"],
        lift = ["-"],
        interred = ["-"],
        inCone = ["-"],
        pertvec = ["-"],
        depth = ["-"],
        example = [v],
    )
    savea(df, "allocsFractalWalklex")
    df = DataFrame(
        weights = ["-"],
        weight = ["-"],
        nGens = ["-"],
        initials = ["-"],
        stdh = ["-"],
        lift2 = ["-"],
        lift = ["-"],
        interred = ["-"],
        inCone = ["-"],
        pertvec = ["-"],
        depth = ["-"],
        example = [v],
    )
    savea(df, "allocsFractalWalklookahead")
    df = DataFrame(
        weights = ["-"],
        weight = ["-"],
        nGens = ["-"],
        initials = ["-"],
        stdh = ["-"],
        lift2 = ["-"],
        lift = ["-"],
        interred = ["-"],
        inCone = ["-"],
        pertvec = ["-"],
        depth = ["-"],
        example = [v],
    )
    savea(df, "allocsFractalWalkcombined")
    df = DataFrame(
        weights = "-",
        weight = "-",
        nGens = ["-"],
        initials = ["-"],
        stdh = ["-"],
        liftGW2 = ["-"],
        lift2 = ["-"],
        interred = ["-"],
        rep = ["-"],
        inCone = ["-"],
        inseveral = ["-"],
        example = [v],
    )
    savea(df, "allocsTranWalk")

    println("Computing deglex-Basis")
    stime =@belapsed Singular.std($ideal, complete_reduction = true) evals = 1 samples = 1
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

    println("Computing lex-Basis")
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
    savea(df, "correct")

    println("Benchmarking ideals")
    for id in ideals
        a = isequal(
            Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
            s,
        )
        b = equalitytest(s, id)
        df = DataFrame(a = [a], b = [b],c=[v])
        savea(df, "correct")
    end
end
