using BenchmarkTools

function prepare()
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
savew(df, "standardWalk.txt")
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
savew(df, "pertubedWalk2.txt")
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
savew(df, "pertubedWalk3.txt")
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
savew(df, "pertubedWalk4.txt")
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
savew(df, "pertubedWalk5.txt")
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
savew(df, "pertubedWalk6.txt")
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
savew(df, "pertubedWalk7.txt")
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
savew(df, "pertubedWalk8.txt")
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
savew(df, "pertubedWalk9.txt")
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
savew(df, "pertubedWalk10.txt")
df = DataFrame(
    weights = ["weights"],
    facetnormal = ["-"],
    initials = ["initials"],
    stdh = ["stdh"],
    liftgeneric = ["liftgeneric"],
    interred = ["interred"],
    example = ["example"],
)
savew(df, "genericWalk.txt")

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
savew(df, "fractalWalk.txt")
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
savew(df, "fractalWalklex.txt")
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
savew(df, "fractalWalklookahead.txt")
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
savew(df, "fractalWalkcombined.txt")
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
savew(df, "tranWalk.txt")
end
function runb(
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
    savea(df, "standardWalk.txt")
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
    savea(df, "pertubedWalk2.txt")
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
    savea(df, "pertubedWalk3.txt")
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
    savea(df, "pertubedWalk4.txt")
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
    savea(df, "pertubedWalk5.txt")
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
    savea(df, "pertubedWalk6.txt")
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
    savea(df, "pertubedWalk7.txt")
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
    savea(df, "pertubedWalk8.txt")
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
    savea(df, "pertubedWalk9.txt")
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
    savea(df, "pertubedWalk10.txt")
    df = DataFrame(
        weights = ["-"],
        facetnormal = ["-"],
        initials = ["-"],
        stdh = ["-"],
        liftgeneric = ["-"],
        interred = ["-"],
        example = [v],
    )
    savea(df, "genericWalk.txt")

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
    savea(df, "fractalWalk.txt")
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
    savea(df, "fractalWalklex.txt")
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
    savea(df, "fractalWalklookahead.txt")
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
    savea(df, "fractalWalkcombined.txt")
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
    savea(df, "tranWalk.txt")


    stime =@belapsed Singular.std($ideal, complete_reduction = true) evals = 1 samples = 1
    I = Singular.std(ideal, complete_reduction = true)
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

    println("Benchmark ideals")
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
