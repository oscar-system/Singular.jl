include("GroebnerWalkUtilitys.jl")
include("FractalWalkUtilitys.jl")
include("GenericWalkUtilitys.jl")
include("StandardWalkUtilitys.jl")
include("TranWalkUtilitys.jl")



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
function groebnerwalk(
    G::Singular.sideal,
    S::Matrix{Int64},
    T::Matrix{Int64},
    grwalktype::Symbol = :standard,
    k::Int64 = 0,
)

    if grwalktype == :standard
        walk = (x, y, z) -> standard_walk(x, y, z)
    elseif grwalktype == :generic
        walk = (x, y, z) -> generic_walk(x, y, z)
    elseif grwalktype == :pertubed
        walk = (x, y, z) -> pertubed_walk(x, y, z, k)
    elseif grwalktype == :fractal
        walk =
            (x, y, z) -> fractal_walk(
                x,
                MonomialOrder(S, S[1, :], [0]),
                MonomialOrder(T, T[1, :], T[1, :]),
            )
    elseif grwalktype == :fractal2
        walk =
            (x, y, z) -> fractal_walk2(
                x,
                MonomialOrder(S, S[1, :], [0]),
                MonomialOrder(T, T[1, :], T[1, :]),
            )
    elseif grwalktype == :fractal3
        walk =
            (x, y, z) -> fractal_walk3(
                x,
                MonomialOrder(S, S[1, :], [0]),
                MonomialOrder(T, T[1, :], T[1, :]),
            )
    elseif grwalktype == :fractal_walk_lookahead
        walk =
            (x, y, z) -> fractal_walk_lookahead(
                x,
                MonomialOrder(S, S[1, :], [0]),
                MonomialOrder(T, T[1, :], T[1, :]),
            )
    elseif grwalktype == :tran
        walk = (x, y, z) -> tran_walk(x, y, z)
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
#Standard-version of the groebner walk by Collart et al
###############################################################
function standard_walk(G::Singular.sideal, S::Matrix{Int64}, T::Matrix{Int64})
    cweight = S[1, :]
    tweight = T[1, :]
    println("Standardwalk results")
    println("Crossed Cones in: ")

    standard_walk(G, S, T, cweight, tweight)
end

function standard_walk(
    G::Singular.sideal,
    S::Matrix{Int64},
    T::Matrix{Int64},
    cweight::Vector{Int64},
    tweight::Vector{Int64},
)

    #loop
    term = false
    while !term
        G = standard_step(G, cweight, tweight, T)
        global counter = getCounter() + 1
        println(cweight)
        if cweight == tweight
            term = true
        else
            cweight = nextw(G, cweight, tweight)
        end
    end
    return G
end

function standard_step(
    G::Singular.sideal,
    cw::Array{Int64,1},
    tw::Array{Int64,1},
    T::Matrix{Int64},
)
    R = base_ring(G)
    S, V = change_order(R, cw, T)

    #taking initials
    inwG = initials(S, gens(G), cw)
    #Initial Ideal
    IinwG = Singular.std(
        Singular.Ideal(S, [S(x) for x in inwG]),
        complete_reduction = true,
    )

    #Lifting to GB of new cone
    Gnew = liftGW(G, IinwG, R, S)
    #Interreduce
    return Singular.std(Gnew, complete_reduction = true)
end

###############################################################
#Generic-version of the groebner walk by Fukuda et al. (2007)
###############################################################


function generic_walk(G::Singular.sideal, S::Matrix{Int64}, T::Matrix{Int64})
    R, V = change_order(G.base_ring, T)

    v = nextV(G, [0], S, T)
    lm = [change_ring(Singular.leading_term(g), R) for g in gens(G)]
    G = Singular.Ideal(R, [change_ring(x, R) for x in gens(G)])

    println("GenericWalk results")
    println("Crossed Cones with facetNormal: ")
    while !isempty(v)
        global counter = getCounter() + 1
        println(v)
        G, lm = generic_step(G, R, v, T, lm)
        #@info "Current Gröbnerbase: " G
        #println("nextv")
         v = nextV(G, lm, v, S, T)
    end
    return G
end

function generic_step(
    G::Singular.sideal,
    S::Singular.PolyRing,
    v::Vector{Int64},
    T::Matrix{Int64},
    lm::Vector{Singular.spoly{L}},
) where {L<:Nemo.RingElem}

    R = Singular.base_ring(G)

    #println("initials")
    facet_Generators = facet_initials(S, G, v, lm)
    #println("std")
    facet_Ideal = Singular.std(
        Singular.Ideal(S, [S(x) for x in facet_Generators]),
        complete_reduction = true,
    )

    #println("liftgeneric")
    liftArray, Newlm = liftgeneric(G, lm, facet_Ideal)

    Gnew = Singular.Ideal(S, [x for x in liftArray])

    #println("interred")
    Gnew = interreduce_generic(Gnew, Newlm)
    Gnew = Singular.Ideal(R, [change_ring(x, R) for x in gens(Gnew)])
    Gnew.isGB = true
    return Gnew, Newlm
end


###############################################################
#Pertubed-version of the groebner walk Amrhein et al.
###############################################################
function pertubed_walk(
    G::Singular.sideal,
    S::Matrix{Int64},
    T::Matrix{Int64},
    p::Int64,
)
    #sweight = pert_Vectors(G, S, p)
    sweight = S[1,:]
    Gnew = G
    #loop
    term = false
    println("PertubedWalk results")
    println("Crossed Cones in: ")

    while !term
        #test = [pert_Vectors(Gnew, T, i) for i = 1:p]
        tweight = pert_Vectors(Gnew, T, p)
        Gnew = standard_walk(Gnew, S, T, sweight, tweight)
        if inCone(Gnew, T, tweight)
            term = true
        else
            if p == 1
                R, V = change_order(Gnew.base_ring, T)
                Gnew =
                    Singular.Ideal(R, [change_ring(x, R) for x in gens(Gnew)])
                Gnew = Singular.std(Gnew, complete_reduction = true)
                term = true
            end
            p = p - 1
            sweight = tweight
        end
    end
    return Gnew
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
#Counter for the steps in a fractalwalk
########################################
counterFr = 0
counterFrFin = false
function deleteCounterFr()
    global counterFr
    temp = counterFr
    counterFr = 0
    global counterFrFin = false
    return temp
end
function getCounterFr()
    global counterFr
    return counterFr
end
function raiseCounterFr(std::Bool = false)
    global counterFrFin
    if !counterFrFin || std
        counterFrFin = true
        global counterFr = getCounterFr() + 1
        return true
    end
    return false
end
function resetCounterFr()
    global counterFrFin = false
end

PVecs = []
function checkPvecs()
    global PVecs
    println(PVecs)
end

sigma = []
#=
function getSigma()
    global sigma
    return sigma
end =#


function fractal_walk(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
)
R= base_ring(G)
    global PVecs = [pert_Vectors(G, T, i) for i = 1:nvars(R)]
    println("FacrtalWalk_standard results")
    println("Crossed Cones in: ")
    Gb = fractal_recursiv(G, S, T, S.w, PVecs, 1)
    println("Cones crossed: ", deleteCounterFr())
    return Gb
end

function fractal_recursiv(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    s::Vector{Int64},
    PVecs::Vector{Vector{Int64}},
    p::Int64,
)
    R = base_ring(G)
    term = false
    G.isGB = true
    cweight = s

    while !term
        w = nextW_2(G, cweight, PVecs[p])
        if (w == [0])
            #Look up if inCone is special for a deeper step
            if inCone(G, T, PVecs[p])
                if raiseCounterFr()
                    println(PVecs[p], " in depth", p)
                end
                return G
            else
                global PVecs = [pert_Vectors(G, T, 2, i) for i = 1:nvars(R)]
                println(PVecs)
                continue
            end
        end
        T.w = w
        Rn, V = change_order(R, T)
        Gw = initials(R, gens(G), w)

        if p == nvars(R)
             Gnew = Singular.std(
                Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]),
                complete_reduction = true,
            )
            println(w, " in depth", p)
            raiseCounterFr(true)
        else
            resetCounterFr()
            println("up in: ", p, " with: ", w)

             Gnew = fractal_recursiv(
                Singular.Ideal(R, [x for x in Gw]),
                S,
                T,
                s,
                PVecs,
                p + 1,
            )
            resetCounterFr()

        end

        G = liftGWfr(G, Gnew, R, Rn)

        G = Singular.std(G, complete_reduction = true)
        R = Rn
        S = T
        cweight = w
    end

    return G
end

cwPert = []
firstStepMode = false
function cwpert(p::Int64)
    cwPert[p]
end

function fractal_walk2(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
)
    global PVecs = [pert_Vectors(G, T, i) for i = 1:nvars(R)]
    global sigma = S.w
    println("FractalWalk_withStartorder results")
    println("Crossed Cones in: ")
    Gb = fractal_recursiv2(G, S, T, PVecs, 1)
    println("Cones crossed: ", deleteCounterFr())
    return Gb
end

function fractal_recursiv2(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    PVecs::Vector{Vector{Int64}},
    p::Int64,
)
    R = Singular.base_ring(G)
    term = false
    G.isGB = true
    if (p == 1)
        if !ismonomial(initials(R, Singular.gens(G), S.w))
            println("nomonomial")
            global cwPert = [pert_Vectors(G, S, S.w, i) for i = 1:nvars(R)]
            global firstStepMode = true
        end
    end
    if firstStepMode
        cweight = cwPert[p]
    else
        cweight = getSigma()
    end

    while !term
        w = nextW_2(G, cweight, PVecs[p])
        if w == [0]

            if inCone(G, T, PVecs[p])
                if raiseCounterFr()
                    println(PVecs[p], " in depth", p)
                end
                return G
            else
                global PVecs = [pert_Vectors(G, T, 2, i) for i = 1:nvars(R)]
                continue
            end
        end
        T.w = w
        Rn, V = change_order(R, T)
        Gw = initials(R, gens(G), w)
        if p == Singular.nvars(R)
            Gnew = Singular.std(
                Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]),
                complete_reduction = true,
            )
            println(w, " in depth", p)
            raiseCounterFr(true)
            #println("Computed Gröbnerbase: ", w, "in depth: ", p)

        else
            Gnew = fractal_recursiv2(
                Singular.Ideal(R, [x for x in Gw]),
                S,
                T,
                PVecs,
                p + 1,
            )
            resetCounterFr()

            global firstStepMode = false
        end

        G = liftGWfr(G, Gnew, R, Rn)
        G = Singular.std(G, complete_reduction = true)
        R = Rn
        S = T
        cweight = w
    end
    return G
end
function fractal_walk3(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
)
    global PVecs = [pert_Vectors(G, T, i) for i = 1:nvars(R)]
    println("FractalWalk_withStartorderLex results")
    println("Crossed Cones in: ")
    Gb = fractal_recursiv3(G, S, T, PVecs, 1)
    println("Cones crossed: ", deleteCounterFr())
    return Gb
end

function fractal_recursiv3(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    PVecs::Vector{Vector{Int64}},
    p::Int64,
)
    R = Singular.base_ring(G)
    term = false
    G.isGB = true
    if (p == 1)
        if !ismonomial(initials(R, Singular.gens(G), S.w))
            global cwPert =
                [pert_Vectors(G, S, S.w, i) for i = 1:Singular.nvars(R)]
            global firstStepMode = true
        end
    end
    if firstStepMode
        cweight = cwPert[p]
    else
        cweight = S.w
    end

    #println(cwPert)
    while !term
        w = nextw_fr(G, cweight, PVecs[p])
        #println(cweight, w)

        if (p == 1 && w == PVecs[p])
            p = p + 1
        end
        if w == [0]

            if inCone(G, T, PVecs[p])
                if raiseCounterFr()
                    println(PVecs[p], " in depth", p)
                end
                return G
            else
                global PVecs =
                    [pert_Vectors(G, T, 2, i) for i = 1:Singular.nvars(R)]
                continue
            end
        end
        T.w = w
        Rn, V = change_order(R, T)
        Gw = initials(R, Singular.gens(G), w)
        if p == Singular.nvars(R)
            Gnew = Singular.std(
                Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]),
                complete_reduction = true,
            )
            println(w, " in depth", p)
            raiseCounterFr(true)
            #println("Computed Gröbnerbase: ", w, "in depth: ", p)

        else
            resetCounterFr
            Gnew = fractal_recursiv3(
                Singular.Ideal(R, [x for x in Gw]),
                S,
                T,
                PVecs,
                p + 1,
            )
            global firstStepMode = false
        end

        G = liftGWfr(G, Gnew, R, Rn)
        G = Singular.std(G, complete_reduction = true)
        R = Rn
        S = T
        cweight = w
    end
    return G
end
function fractal_walk_lookahead(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
)
    println("FractalWalk_lookahead results")
    println("Crossed Cones in: ")
    R = base_ring(G)
    global PVecs = [pert_Vectors(G, T, i) for i = 1:nvars(R)]

    Gb = fractal_recursiv_lookahead(G, S, T, PVecs, 1)
    println("Cones crossed: ", deleteCounterFr())
    return Gb
end

function fractal_recursiv_lookahead(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    PVecs::Vector{Vector{Int64}},
    p::Int64,
)
    R = Singular.base_ring(G)
    term = false
    G.isGB = true
    cweight = S.w

    while !term
        w = nextW_2(G, cweight, PVecs[p])
        if (w == [0])
            if inCone(G, T, PVecs[p])
                if raiseCounterFr()
                    println(PVecs[p], " in depth", p)
                end
                return G
            else
                global PVecs = [pert_Vectors(G, T, 2, i) for i = 1:nvars(R)]
                #checkPvecs()
                continue
            end
        end
        T.w = w
        Rn, V = change_order(R, T)
        Gw = initials(R, Singular.gens(G), w)
        if (p == Singular.nvars(R) || isBinomial(Gw))
            Gnew = Singular.std(
                Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]),
                complete_reduction = true,
            )
            println(w, " in depth", p)
            raiseCounterFr(true)
        else
            println("up in: ", p, " with: ", w)

            resetCounterFr()
            Gnew = fractal_recursiv_lookahead(
                Singular.Ideal(R, [x for x in Gw]),
                S,
                T,
                PVecs,
                p + 1,
            )
        end

        G = liftGWfr(G, Gnew, R, Rn)
        G = Singular.std(G, complete_reduction = true)
        R = Rn
        S = T
        cweight = w
    end
    return G
end

###############################################################
#Tran-version of the groebner walk by Tran (2002)
###############################################################

function tran_walk(G::Singular.sideal, S::Matrix{Int64}, T::Matrix{Int64})
    cweight = S[1, :]
    tweight = T[1, :]
    println("Tranwalk results")
    println("Crossed Cones in: ")

    #loop
    term = false
    while !term
        w = nextw(G, cweight, tweight)
        if w == tweight
            if inCone(G, T, cweight)
                return G
            else
                #TODO: Implement: if in several cones
                tweight = representationVector(G, T)
                continue
            end
        end
        G = standard_step(G, w, tweight, T)
        global counter = getCounter() + 1
        println(w)
        cweight = w
    end
end
