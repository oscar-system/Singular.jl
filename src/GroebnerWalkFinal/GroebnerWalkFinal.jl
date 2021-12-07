include("GroebnerWalkUtilitysFinal.jl")
include("FractalWalkUtilitysFinal.jl")
include("GenericWalkUtilitysFinal.jl")
include("StandardWalkUtilitysFinal.jl")
include("TranWalkUtilitysFinal.jl")



export groebnerwalk

###############################################################
#Implementation of the gröbner walk.
###############################################################

#for counting the steps of the groebnerwalk.
counter = 0
function getCounter()
    global counter
    temp = counter
    counter = 0
    return temp
end

#Top-level function
#S,T have to be nxn-matrices with full rank
function groebnerwalk(
    G::Singular.sideal,
    S::Matrix{Int64},
    T::Matrix{Int64},
    grwalktype::Symbol = :standard,
    p::Int64 = 0,
)

    if grwalktype == :standard
        walk = (x, y, z) -> StandardWalk(x, y, z)
    elseif grwalktype == :generic
        walk = (x, y, z) -> GenericWalk(x, y, z)
    elseif grwalktype == :pertubed
        walk = (x, y, z) -> PertubedWalk(x, y, z, p)
    elseif grwalktype == :fractal
        walk =
            (x, y, z) -> FractalWalk(
                x,
                MonomialOrder(S, S[1, :], [0]),
                MonomialOrder(T, T[1, :], T[1, :]),
            )
    elseif grwalktype == :fractalStartPertubation
        walk =
            (x, y, z) -> FractalWalkStartOrder(
                x,
                MonomialOrder(S, S[1, :], [0]),
                MonomialOrder(T, T[1, :], T[1, :]),
            )
    elseif grwalktype == :fractalLex
        walk =
            (x, y, z) -> FractalWalkLex(
                x,
                MonomialOrder(S, S[1, :], [0]),
                MonomialOrder(T, T[1, :], T[1, :]),
            )
    elseif grwalktype == :fractalLookAhead
        walk =
            (x, y, z) -> FractalWalkLookAhead(
                x,
                MonomialOrder(S, S[1, :], [0]),
                MonomialOrder(T, T[1, :], T[1, :]),
            )
    elseif grwalktype == :tran
        walk = (x, y, z) -> TranWalk(x, y, z)
    elseif grwalktype == :alternative
        walk =
            (x, y, z) -> alternative_algorithm_top(
                x,
                MonomialOrder(S, S[1, :], [0]),
                MonomialOrder(T, T[1, :], T[1, :]),
            )
    end

    ######TODO:Check the parameter#####
    R = base_ring(G)
    I = Singular.Ideal(R, [R(x) for x in gens(G)])

    Gb = walk(I, S, T)
    println("Cones crossed: ", getCounter())

    S, V = change_order(Gb.base_ring, T)
    return Singular.Ideal(S, [change_ring(gen, S) for gen in gens(Gb)])
end


###############################################################
#Standard-version of the groebner walk by Collart et al.
###############################################################
function StandardWalk(G::Singular.sideal, S::Matrix{Int64}, T::Matrix{Int64})
    println("Standardwalk results")
    println("Crossed Cones in: ")
    StandardWalk(G, S, T, S[1, :], T[1, :])
end

function StandardWalk(
    G::Singular.sideal,
    S::Matrix{Int64},
    T::Matrix{Int64},
    cweight::Vector{Int64},
    tweight::Vector{Int64},
)
    #loop
    terminate = false
    while !terminate
        G = StandardStep(G, cweight, tweight, T)
        println(cweight)
        global counter = getCounter() + 1
        if cweight == tweight
            terminate = true
        else
            cweight = NextWeight(G, cweight, tweight)
        end
    end
    return G
end

function StandardStep(
    G::Singular.sideal,
    cw::Vector{Int64},
    tw::Vector{Int64},
    T::Matrix{Int64},
)
    R = base_ring(G)
    S, V = change_order(R, cw, T)

    Gw = Initials(S, gens(G), cw)

    H = Singular.std(
        Singular.Ideal(S, [S(x) for x in Gw]),
        complete_reduction = true,
    )

    G = Lift(G, H, R, S)
    return Singular.std(G, complete_reduction = true)
end

###############################################################
#Generic-version of the groebner walk by Fukuda et al. (2007)
###############################################################


function GenericWalk(G::Singular.sideal, S::Matrix{Int64}, T::Matrix{Int64})
    R, V = change_order(G.base_ring, T)

    v = NextGamma(G, [0], S, T)
    Lm = [change_ring(Singular.leading_term(g), R) for g in gens(G)]
    G = Singular.Ideal(R, [change_ring(x, R) for x in gens(G)])

    println("GenericWalk results")
    println("Crossed Cones with facetNormal: ")
    while !isempty(v)
        global counter = getCounter() + 1
        println(v)
        G, Lm = GenericStep(G, R, v, T, Lm)
        v = NextGamma(G, Lm, v, S, T)
    end
    return Singular.interreduce(G)
end

function GenericStep(
    G::Singular.sideal,
    S::Singular.PolyRing,
    v::Vector{Int64},
    T::Matrix{Int64},
    Lm::Vector{Singular.spoly{L}},
) where {L<:Nemo.RingElem}

    R = Singular.base_ring(G)

    facet_Generators = FacetInitials(S, G, v, Lm)
    H = Singular.std(
        Singular.Ideal(S, [S(x) for x in facet_Generators]),
        complete_reduction = true,
    )
    H, Lm = LiftGeneric(G, Lm, H)
    G = Singular.Ideal(S, [x for x in H])
    G = interreduce(G, Lm)
    G = Singular.Ideal(R, [change_ring(x, R) for x in gens(G)])
    G.isGB = true
    return G, Lm
end


###############################################################
#Pertubed-version of the groebner walk Amrhein et al.
###############################################################
function PertubedWalk(
    G::Singular.sideal,
    S::Matrix{Int64},
    T::Matrix{Int64},
    p::Int64,
)
    #sweight = pert_Vectors(G, S, p)
    sweight = S[1, :]
    terminate = false
    println("PertubedWalk results")
    println("Crossed Cones in: ")

    while !terminate
        tweight = PertubedVector(G, T, p)
        G = StandardWalk(G, S, T, sweight, tweight)
        if InCone(G, T, tweight)
            terminate = true
        else
            if p == 1
                R, V = change_order(G.base_ring, T)
                G =
                    Singular.Ideal(R, [change_ring(x, R) for x in gens(G)])
                G = Singular.std(G, complete_reduction = true)
                term = true
            end
            p = p - 1
            sweight = tweight
        end
    end
    return G
end

###############################################################
#fractal-walk by Amrhein et al.
#Inlcuding:
#fractal_walk -> standard-version
#fractal_walk2 -> checks if the starting weight is in the inner of a cone.
#fractal_walk3 -> fractal walk expecially for conversion to the lexikographic orderig.
#                 checks if the starting weight is in the inner of a cone.
###############################################################

########################################
#Counter for the steps in the fractalwalk
########################################
counterFr = 0
function deleteCounterFr()
    global counterFr
    temp = counterFr
    counterFr = 0
    return temp
end
function getCounterFr()
    global counterFr
    return counterFr
end
function raiseCounterFr(std::Bool = false)
        global counterFr = getCounterFr() + 1

end

PVecs = []

sigma = []

function FractalWalk(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
)
    global PVecs = [PertubedVector(G, T, i) for i = 1:nvars(base_ring(G))]
    println("FacrtalWalk_standard results")
    println("Crossed Cones in: ")
    Gb = FractalRecursive(G, S, T, S.w, PVecs, 1)
    println("Cones crossed: ", deleteCounterFr())
    return Gb
end

function FractalRecursive(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    s::Vector{Int64},
    PVecs::Vector{Vector{Int64}},
    p::Int64,
)
    R = base_ring(G)
    terminate = false
    G.isGB = true
    w = s

    while !terminate
        #w = nextW_2(G, cweight, PVecs[p])
        t = NextT(G, w, PVecs[p])
        if (t == [0])
            if InCone(G, T, PVecs[p])
                return G
            else
                global PVecs = [PertubedVector(G, T, i) for i = 1:nvars(R)]
                continue
            end
        end
        w = w + t * (PVecs[p] - w)
        w = convertBoundingVector(w)
        T.w = w
        Rn, V = change_order(R, T)
        Gw = Initials(R, gens(G), w)
        if p == nvars(R)
            H = Singular.std(
                Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]),
                complete_reduction = true,
            )
            println(w, " in depth", p)
            raiseCounterFr(true)
        else
            println("up in: ", p, " with: ", w)

            H = FractalRecursive(
                Singular.Ideal(R, [x for x in Gw]),
                S,
                T,
                s,
                PVecs,
                p + 1,
            )
        end
        H = LiftFractalWalk(G, H, R, Rn)
        G = Singular.std(H, complete_reduction = true)
        R = Rn
        S = T
    end
    return G
end

cwPert = []
firstStepMode = false
function cwpert(p::Int64)
    cwPert[p]
end

function FractalWalkStartOrder(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
)
    global PVecs =
        [PertubedVector(G, T, i) for i = 1:nvars(Singular.base_ring(G))]
    global sigma = S.w
    println("FractalWalk_withStartorder results")
    println("Crossed Cones in: ")
    Gb = FractalWalkRecursiveStartOrder(G, S, T,S.w, PVecs, 1)
    println("Cones crossed: ", deleteCounterFr())
    return Gb
end

function FractalWalkRecursiveStartOrder(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    s::Vector{Int64},
    PVecs::Vector{Vector{Int64}},
    p::Int64,
)
    R = Singular.base_ring(G)
    terminate = false
    G.isGB = true
    if (p == 1)
        if !IsMonomial(Initials(R, Singular.gens(G), S.w))
            global cwPert = [PertubedVector(G, S, S.w, i) for i = 1:nvars(R)]
            global firstStepMode = true
        end
    end
    if firstStepMode
        w = cwPert[p]
    else
        w = s
    end

    while !terminate
        t = NextT(G, w, PVecs[p])
        if t == [0]
            if InCone(G, T, PVecs[p])
                println(PVecs[p], " in depth", p)
                return G
            else
                global PVecs = [PertubedVector(G, T, i) for i = 1:nvars(R)]
                continue
            end
        end
        w = w + t * (PVecs[p] - w)
        w = convertBoundingVector(w)
        T.w = w
        Rn, V = change_order(R, T)
        Gw = Initials(R, gens(G), w)
        if p == Singular.nvars(R)
            H = Singular.std(
                Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]),
                complete_reduction = true,
            )
            println(w, " in depth", p)
            raiseCounterFr(true)
        else
            println("up in: ", p, " with: ", w)

            H = FractalWalkRecursiveStartOrder(
                Singular.Ideal(R, [x for x in Gw]),
                S,
                T,
                s,
                PVecs,
                p + 1,
            )
            global firstStepMode = false
        end
        H = LiftFractalWalk(G, H, R, Rn)
        G = Singular.std(H, complete_reduction = true)
        R = Rn
        S = T
    end
    return G
end
function FractalWalkLex(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
)
    global PVecs = [PertubedVector(G, T, i) for i = 1:nvars(base_ring(G))]
    println("FractalWalkLex results")
    println("Crossed Cones in: ")
    Gb = FractalWalkRecursiveLex(G, S, T,S.w, PVecs, 1)
    println("Cones crossed: ", deleteCounterFr())
    return Gb
end

function FractalWalkRecursiveLex(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    s::Vector{Int64},
    PVecs::Vector{Vector{Int64}},
    p::Int64,
)
    R = Singular.base_ring(G)
    terminate = false
    G.isGB = true
    w = s
    while !terminate
        t = NextT(G, w, PVecs[p])
        if t == [0]
            if InCone(G, T, PVecs[p])
                return G
            else
                global PVecs =
                    [PertubedVector(G, T, i) for i = 1:Singular.nvars(R)]
                continue
            end
        end
        w = w + t * (PVecs[p] - w)
        w = convertBoundingVector(w)
        T.w = w
        Rn, V = change_order(R, T)
        Gw = Initials(R, Singular.gens(G), w)
        if p == Singular.nvars(R)
            H = Singular.std(
                Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]),
                complete_reduction = true,
            )
            println(w, " in depth", p)
            raiseCounterFr(true)
        else
            println("up in: ", p, " with: ", w)

            H = FractalWalkRecursiveLex(
                Singular.Ideal(R, [x for x in Gw]),
                S,
                T,
                s,
                PVecs,
                p + 1,
            )
            global firstStepMode = false
        end
        H = LiftFractalWalk(G, H, R, Rn)
        G = Singular.std(H, complete_reduction = true)
        R = Rn
        S = T
    end
    return G
end
function FractalWalkLookAhead(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
)
    println("FractalWalkLookAhead results")
    println("Crossed Cones in: ")
    global PVecs = [PertubedVector(G, T, i) for i = 1:nvars(base_ring(G))]
    Gb = FractalWalkLookAheadRecursive(G, S, T, PVecs, 1)
    println("Cones crossed: ", deleteCounterFr())
    return Gb
end

function FractalWalkLookAheadRecursive(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    s::Vector{Int64},
    PVecs::Vector{Vector{Int64}},
    p::Int64,
)
    R = Singular.base_ring(G)
    terminate = false
    G.isGB = true
    w = s

    while !terminate
        t = NextT(G, cweight, PVecs[p])
        if t == [0]
            if inCone(G, T, PVecs[p])
                return G
            else
                global PVecs = [PertubedVector(G, T, i) for i = 1:nvars(R)]
                continue
            end
        end
        w = w + t * (PVecs[p] - w)
        w = convertBoundingVector(w)
        T.w = w
        Rn, V = change_order(R, T)
        Gw = Initials(R, Singular.gens(G), w)
        if (p == Singular.nvars(R) || IsBinomial(Gw))
            H = Singular.std(
                Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]),
                complete_reduction = true,
            )
            println(w, " in depth", p)
            raiseCounterFr(true)
        else
            println("up in: ", p, " with: ", w)
            H = FractalWalkLookAheadRecursive(
                Singular.Ideal(R, [x for x in Gw]),
                S,
                T,
                s,
                PVecs,
                p + 1,
            )
        end
        H = LiftFractalWalk(G, H, R, Rn)
        G = Singular.std(H, complete_reduction = true)
        R = Rn
        S = T
    end
    return G
end

###############################################################
#Tran-version of the groebner walk by Tran (2002)
###############################################################

function TranWalk(G::Singular.sideal, S::Matrix{Int64}, T::Matrix{Int64})
    cweight = S[1, :]
    tweight = T[1, :]
    println("Tranwalk results")
    println("Crossed Cones in: ")

    #loop
    term = false
    while !term
        w = NextWeight(G, cweight, tweight)
        if w == tweight
            if InCone(G, T, cweight)
                return G
            else
                if InSeveralCones(Initials(base_ring(G), gens(G), w))
                    tweight = RepresentationVector(G, T)
                    continue
                end
            end
        end
        G = StandardStep(G, w, tweight, T)
        global counter = getCounter() + 1
        println(w)
        cweight = w
    end
end
