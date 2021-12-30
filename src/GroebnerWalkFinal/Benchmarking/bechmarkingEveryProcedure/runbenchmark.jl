using BenchmarkTools

function prepare()
df = DataFrame(
    NextWeight = ["NextWeight"],
    weight = ["weight"],
    nGens = ["nGens"],
    initials = ["initials"],
    stdH = ["stdH"],
    liftGW2 = ["liftGW2"],
    lift = ["lift"],
    interred = ["interred"],
    example = ["example"],
)
savew(df, "standardWalk")
df = DataFrame(
    NextWeight = ["NextWeight"],
    weight = ["weight"],
    nGens = ["nGens"],
    initials = ["initials"],
    stdh = ["stdh"],
    liftGW2 = ["liftGW2"],
    lift = ["lift"],
    interred = ["interred"],
    pert = ["pert"],
    inCone = ["inCone"],
    laststd = ["laststd"],
    degree = ["degree"],
    example = ["example"],
)
for i in 2:10
savew(df, "pertubedWalk",i)
end
df = DataFrame(
    NextWeight = ["NextWeight"],
    facetnormal = ["-"],
    initials = ["initials"],
    stdh = ["stdh"],
    liftgeneric = ["liftgeneric"],
    interred = ["interred"],
    example = ["example"],
)
savew(df, "genericWalk")

df = DataFrame(
    NextWeight = ["NextWeight"],
    weight = ["weight"],
    nGens = ["nGens"],
    initials = ["initials"],
    stdh = ["stdh"],
    liftGW2 = ["liftGW2"],
    lift = ["lift"],
    interred = ["interred"],
    inCone = ["inCone"],
    pertvec = ["pertvec"],
    depth = ["depth"],
    example = ["example"],
)
savew(df, "fractalWalk")
df = DataFrame(
    NextWeight = ["NextWeight"],
    weight = ["weight"],
    nGens = ["nGens"],
    initials = ["initials"],
    stdh = ["stdh"],
    liftGW2 = ["liftGW2"],
    lift = ["lift"],
    interred = ["interred"],
    inCone = ["inCone"],
    pertvec = ["pertvec"],
    depth = ["depth"],
    example = ["example"],
)
savew(df, "fractalWalklex")
df = DataFrame(
    NextWeight = ["NextWeight"],
    weight = ["weight"],
    nGens = ["nGens"],
    initials = ["initials"],
    stdh = ["stdh"],
    liftGW2 = ["liftGW2"],
    lift = ["lift"],
    interred = ["interred"],
    inCone = ["inCone"],
    pertvec = ["pertvec"],
    depth = ["depth"],
    example = ["example"],
)
savew(df, "fractalWalklookahead")
df = DataFrame(
    NextWeight = ["NextWeight"],
    weight = ["weight"],
    nGens = ["nGens"],
    initials = ["initials"],
    stdh = ["stdh"],
    liftGW2 = ["liftGW2"],
    lift = ["lift"],
    interred = ["interred"],
    inCone = ["inCone"],
    pertvec = ["pertvec"],
    depth = ["depth"],
    example = ["example"],
)
savew(df, "fractalWalkcombined")
df = DataFrame(
    NextWeight = "NextWeight",
    weight = "weight",
    nGens = ["nGens"],
    initials = ["initials"],
    stdh = ["stdh"],
    liftGW2 = ["liftGW2"],
    lift = ["lift"],
    interred = ["interred"],
    rep = ["rep"],
    inCone = ["inCone"],
    inseveral = ["inseveral"],
    example = ["example"],
)
savew(df, "tranWalk")
end
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
        a = isequal(
            Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
            s,
        )
        b = equalitytest(s, id)
        df = DataFrame(a = [a], b = [b],c=[v])
        savea(df, "correct")
    end
end
