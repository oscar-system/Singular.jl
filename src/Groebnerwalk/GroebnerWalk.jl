using Oscar
include("GroebnerWalkUtilitys.jl")
include("FractalWalkUtilitys.jl")
include("GenericWalkUtilitys.jl")
include("StandardWalkUtilitys.jl")


export groebnerwalk

###############################################################
#Implementation of the gröbner walk.
#TODO: Improve pertubed version.
#TODO: Improve GenWalk algorithm
#TODO: Improve fractal version.
###############################################################

###############################################################
#Top-Level
#TODO: Implement input checks
#Fragen:
#Singular.reduce für f % G selbst programmieren?
#Für divrem in Oscar-Ideal: sinnvoller möglich? reduce liefert nur den rest
#Arithmetik von Polynomen -> Immer MPolyBuildCTX?
###############################################################

counter = 0
function getCounter()
    global counter
    temp = counter
    counter = 0
    return temp
end
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
        walk = (x, y, z) -> fractal_walk(x, y, z)
    elseif grwalktype == :fractal2
        walk = (x, y, z) -> fractal_walk2(x, y, z)
    elseif grwalktype == :fractal3
        walk = (x, y, z) -> fractal_walk3(x, y, z)
    elseif grwalktype == :fractal4
        walk = (x, y, z) -> fractal_walk_lookahead(x, y, z)
        #error("Choose a strategy from: :generic, :standard ...")
        # fill the list or implement a functionality to choose the best
        # way w.r.t. the type of the ideal
    end

    ######TODO:Check the parameter#####
    R = singular_ring(base_ring(G))
    I = Singular.Ideal(R, [R(x) for x in gens(G)])

    Gb = walk(I, S, T)
    println("Cones crossed: ", getCounter())

    S, V = change_order(Gb, T)
    return Oscar.Singular.Ideal(S, [change_ring(gen, S) for gen in gens(Gb)])
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
    S, V = change_order(G, cw, T)

    #taking initials
    inwG = initials(S, gens(G), cw)
    #Initial Ideal
    IinwG = Singular.std(
        Oscar.Singular.Ideal(S, [S(x) for x in inwG]),
        complete_reduction = true,
    )

    #Lifting to GB of new cone
    Gnew = liftGW(G, IinwG, R, S)

    #Interreduce
    return Oscar.Singular.std(Gnew, complete_reduction = true)
end

###############################################################
#Generic-version of the groebner walk by Fukuda et al. (2007)
###############################################################


function generic_walk(G::Singular.sideal, S::Matrix{Int64}, T::Matrix{Int64})
    R, V = change_order(G, T)

    v = nextV(G, [0], S, T)
    lm = [Oscar.Singular.leading_monomial(g) for g in gens(G)]
    println("GenericWalk results")
    println("Crossed Cones with facetNormal: ")
    while !isempty(v)
        global counter = getCounter() + 1
        println(v)
        G, lm = generic_step(G, R, v, T, lm)
        #@info "Current Gröbnerbase: " G
        v = nextV(G, lm, v, S, T)
    end
    return G
end

function generic_step(
    G::Singular.sideal,
    S::MPolyRing,
    v::Vector{Int64},
    T::Matrix{Int64},
    lm::Vector{Singular.spoly{Singular.n_FieldElem{fmpq}}},
)
    R = base_ring(G)
    facet_Generators = facet_initials(S, G, v, lm)

    facet_Ideal = Singular.std(
        Oscar.Singular.Ideal(S, [S(x) for x in facet_Generators]),
        complete_reduction = true,
    )

    liftArray, Newlm = liftgeneric(G, lm, facet_Ideal)
    Gnew = Oscar.Singular.Ideal(S, [x for x in liftArray])

    Gnew = interreduce_new(Gnew, Newlm)
    Gnew = Oscar.Singular.Ideal(R, [change_ring(x, R) for x in gens(Gnew)])
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
    k::Int64,
)
    p = k
    sweight = pert_Vectors(G, S, p)
    Gnew = G
    #loop
    term = false
    println("PertubedWalk results")
    println("Crossed Cones in: ")

    while !term
        tweight = pert_Vectors(Gnew, T, p)
        Gnew = standard_walk(Gnew, S, T, sweight, tweight)

        if inCone(Gnew, T, tweight)
            term = true
        else
            if k == 1
                R, V = change_order(Gnew, T)
                Gnew = Oscar.Singular.Ideal(
                    R,
                    [change_ring(x, R) for x in gens(Gnew)],
                )
                Gnew = Singular.std(Gnew, complete_reduction = true)
                term = true
            end

            p = p - 1
            sweight = tweight
        end
    end
    return Gnew
end
#=
function fractal_walk(G::Singular.sideal, S::Matrix{Int64}, T::Matrix{Int64})
    PVecs = [pert_Vectors(G, T, i) for i = 1:nvars(R)]

    return fractal_recursiv(G, S, T,S[1,:], PVecs, 1)
end

function fractal_recursiv(
    G::Singular.sideal,
    S::Matrix{Int64},
    T::Matrix{Int64},
    cweight::Vector{Int64},
    PVecs::Vector{Vector{Int64}},
    p::Int64,
)
    R = base_ring(G)
    term = false
    G.isGB = true

    while !term
        w = nextw(G, cweight, PVecs[p])
        if w == PVecs[p]
            if inCone(G, T, w)
                println(G, T, w)
                println(inCone(G,T,w))
                term = true
                break
            else
                PVecs = [pert_Vectors(G, T, i) for i = 1:nvars(R)]
            end
        end
        Rn, V = change_order(G, w, T)
        Gw = initials(R, gens(G), w)
        if p == nvars(R)
            Gnew = Singular.std(Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]), complete_reduction = true)
        else
            Gnew = fractal_recursiv(Singular.Ideal(R, [x for x in Gw]), S, T, PVecs, p + 1)
        end

        rest = [
            change_ring(gen, Rn) - change_ring(Oscar.Singular.reduce(change_ring(gen, R), G), Rn)
            for gen in gens(Gnew)
        ]
        G = Oscar.Singular.Ideal(Rn, [Rn(x) for x in rest])
        G.isGB = true
        G = Singular.std(G, complete_reduction = true)
        R = Rn
    end
    @info "Current Gröbnerbase: " G

    return G
end
=#

PVecs = []
function checkPvecs()
    global PVecs
    println(PVecs)
end
function fractal_walk(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
)
global PVecs = [pert_Vectors(G, T, i) for i = 1:nvars(R)]

println("FacrtalWalk_standard results")
println("Crossed Cones in: ")
Gb = fractal_recursiv(G, S, T, PVecs, 1)
println("Cones crossed: ", getCounter())
return Gb
end

function fractal_recursiv(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    PVecs::Vector{Vector{Int64}},
    p::Int64,
)
    R = base_ring(G)
    term = false
    G.isGB = true
    cweight = S.w


    while !term
        w = nextw_fr(G, cweight, PVecs[p])
        if (w == [0])
            if inCone(G, T, PVecs[p])
                #@info "Incone Gröbnerbase: " G "in depth: " p
                #println(inCone)
                return G
            else
                #println("test")
                global PVecs = [pert_Vectors(G, T, i) for i = 1:nvars(R)]
                #checkPvecs()
                continue
            end
        end
        T.w = w
        Rn, V = change_order(G, T)
        Gw = initials(R, gens(G), w)
        if p == nvars(R)
            Gnew = Singular.std(
                Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]),
                complete_reduction = true,
            )
            println(w, " in depth", p)
            #@info "Computed Gröbnerbase: " Gnew "in depth: " p
            #println("Computed Gröbnerbase: ", w, "in depth: ", p)
            global counter = getCounter() + 1
        else
            Gnew = fractal_recursiv(
                Singular.Ideal(R, [x for x in Gw]),
                S,
                T,
                PVecs,
                p + 1,
            )
            #println("Computed Gröbnerbase: ", w, "in depth: ", p)

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

println("FractalWalk_withStartorder results")
println("Crossed Cones in: ")
Gb = fractal_recursiv2(G, S, T, PVecs, 1)
println("Cones crossed: ", getCounter())
return Gb
end

function fractal_recursiv2(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    PVecs::Vector{Vector{Int64}},
    p::Int64,
)
    R = base_ring(G)
    term = false
    G.isGB = true
    if (p == 1)
        if !ismonomial(initials(R, gens(G), S.w))
            println("nomonomial")
            global cwPert = [pert_Vectors(G, S, S.w, i) for i = 1:nvars(R)]
            global firstStepMode = true
        end
    end
    if firstStepMode
        cweight = cwPert[p]
    else
        cweight = S.w
    end

    while !term
        w = nextw_fr(G, cweight, PVecs[p])
        if w == [0]

            if inCone(G, T, PVecs[p])
                #@info "Incone Gröbnerbase: " G "in depth: " p
                #println(inCone)
                return G
            else
                global PVecs = [pert_Vectors(G, T, i) for i = 1:nvars(R)]
                continue
            end
        end
        T.w = w
        Rn, V = change_order(G, T)
        Gw = initials(R, gens(G), w)
        if p == nvars(R)
            Gnew = Singular.std(
                Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]),
                complete_reduction = true,
            )
            println(w, " in depth: ", p)
            global counter = getCounter() + 1
            #println("Computed Gröbnerbase: ", w, "in depth: ", p)

        else
            Gnew = fractal_recursiv2(
                Singular.Ideal(R, [x for x in Gw]),
                S,
                T,
                PVecs,
                p + 1,
            )
            global firstStepMode = false
            #println("Computed Gröbnerbase: ", w, "in depth: ", p)
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
println("Cones crossed: ", getCounter())
return Gb
end

function fractal_recursiv3(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    PVecs::Vector{Vector{Int64}},
    p::Int64,
)
    R = base_ring(G)
    term = false
    G.isGB = true
    if (p == 1)
        if !ismonomial(initials(R, gens(G), S.w))
            global cwPert = [pert_Vectors(G, S, S.w, i) for i = 1:nvars(R)]
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
                #@info "Incone Gröbnerbase: " G "in depth: " p
                #println(inCone)
                return G
            else
                global PVecs = [pert_Vectors(G, T, i) for i = 1:nvars(R)]
                continue
            end
        end
        T.w = w
        Rn, V = change_order(G, T)
        Gw = initials(R, gens(G), w)
        if p == nvars(R)
            Gnew = Singular.std(
                Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]),
                complete_reduction = true,
            )
            global counter = getCounter() + 1
            println(w, " in depth: ", p)
            #println("Computed Gröbnerbase: ", w, "in depth: ", p)

        else
            Gnew = fractal_recursiv3(
                Singular.Ideal(R, [x for x in Gw]),
                S,
                T,
                PVecs,
                p + 1,
            )
            global firstStepMode = false
            #println("Computed Gröbnerbase: ", w, "in depth: ", p)
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
    global PVecs = [pert_Vectors(G, T, i) for i = 1:nvars(R)]

    Gb = fractal_recursiv_lookahead(G, S, T, PVecs, 1)
    println("Cones crossed: ", getCounter())
    return Gb
end

function fractal_recursiv_lookahead(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    PVecs::Vector{Vector{Int64}},
    p::Int64,
)
    R = base_ring(G)
    term = false
    G.isGB = true
    cweight = S.w


    while !term
        w = nextw_fr(G, cweight, PVecs[p])
        if (w == [0])
            if inCone(G, T, PVecs[p])
                #@info "Incone Gröbnerbase: " G "in depth: " p
                return G
            else
                global PVecs = [pert_Vectors(G, T, i) for i = 1:nvars(R)]
                #checkPvecs()
                continue
            end
        end
        T.w = w
        Rn, V = change_order(G, T)
        Gw = initials(R, gens(G), w)
        if (p == nvars(R) || isBinomial(Gw))
            Gnew = Singular.std(
                Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]),
                complete_reduction = true,
            )
            println(w, " in depth: ", p)
            global counter = getCounter() + 1
            #@info "Computed Gröbnerbase: " Gnew "in depth: " p
            #println("Computed Gröbnerbase: ", w, "in depth: ", p)

        else
            Gnew = fractal_recursiv_lookahead(
                Singular.Ideal(R, [x for x in Gw]),
                S,
                T,
                PVecs,
                p + 1,
            )
            #println("Computed Gröbnerbase: ", w, "in depth: ", p)

        end

        G = liftGWfr(G, Gnew, R, Rn)
        G = Singular.std(G, complete_reduction = true)
        R = Rn
        S = T
        cweight = w
    end
    return G
end
