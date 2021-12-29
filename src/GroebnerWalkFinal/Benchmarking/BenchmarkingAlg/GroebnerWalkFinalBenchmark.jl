include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/GroebnerWalkUtilitysFinal.jl")
include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/FractalWalkUtilitysFinal.jl")
include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/GenericWalkUtilitysFinal.jl")
include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/StandardWalkUtilitysFinal.jl")
include("/Users/JordiWelp/github/Singular.jl/src/GroebnerWalkFinal/TranWalkUtilitysFinal.jl")
using BenchmarkTools

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
#=
@doc Markdown.doc"""
function groebnerwalk(
    G::Singular.sideal,
    S::Matrix{Int64},
    T::Matrix{Int64},
    grwalktype::Symbol = :standard,
    p::Int64 = 0,
)
Given an Ideal $G$ generated by a reduced Groebner Basis w.r.t. the monomial ordering $S$ this function
returns a reduced Groebner Basis w.r.t. the monomial ordering $T$ by converting it using the Groebner Walk.
The Groebner Walk is proposed by Collart et al. (1993)
One can choose a strategy of:
Standard Walk (:standard) computes the Walk like it´s presented in Cox et al. (2005).
Generic Walk (:generic) computes the Walk like it´s presented in Fukuda et al. (2006).
Pertubed Walk (:pertubed, with $p$ = Pertubation degree) computes the Walk like it´s presented in Amrhein et al. (1997).
Tran´s Walk (:tran) computes the Walk like it´s presented in Tran (2000).
Fractal Walk (:fractal) computes the Walk like it´s presented in Amrhein & Gloor (1998). Pertubes only the target vector.
Fractal Walk (:fractal_start_order) computes the Walk like it´s presented in Amrhein & Gloor (1998). Pertubes oth, the start und the target vector.
Fractal Walk (:fractal_lex) computes the Walk like it´s presented in Amrhein & Gloor (1998) in the special case that $T$ represents the lex ordering. Pertubes only the target vector.
Fractal Walk (:factal_look_ahead) computes the Walk like it´spresented in Amrhein & Gloor (1998). This Version uses the buchberger algorithm under certain circumstances before reaching the maximal pertubation depth.
"""=#
function groebnerwalk2(
    G::Singular.sideal,
    S::Matrix{Int64},
    T::Matrix{Int64},
    grwalktype::Symbol = :standard,
    p::Int64 = 0,
)
    if grwalktype == :standard
        walk = (x, y, z) -> standard_walk2(x, y, z)
    elseif grwalktype == :generic
        walk = (x, y, z) -> generic_walk2(x, y, z)
    elseif grwalktype == :pertubed
        walk = (x, y, z) -> pertubed_walk2(x, y, z, p)
    elseif grwalktype == :fractal
        walk =
            (x, y, z) -> fractal_walk2(
                x,
                MonomialOrder(S, S[1, :], [0]),
                MonomialOrder(T, T[1, :], T[1, :]),
            )
    elseif grwalktype == :fractal_start_order
        walk =
            (x, y, z) -> fractal_walk2_start_order(
                x,
                MonomialOrder(S, S[1, :], [0]),
                MonomialOrder(T, T[1, :], T[1, :]),
            )
    elseif grwalktype == :fractal_lex
        walk =
            (x, y, z) -> fractal_walk2_lex(
                x,
                MonomialOrder(S, S[1, :], [0]),
                MonomialOrder(T, T[1, :], T[1, :]),
            )
    elseif grwalktype == :fractal_look_ahead
        walk =
            (x, y, z) -> fractal_walk2_look_ahead(
                x,
                MonomialOrder(S, S[1, :], [0]),
                MonomialOrder(T, T[1, :], T[1, :]),
            )
    elseif grwalktype == :tran
        walk = (x, y, z) -> tran_walk2(x, y, z)
    elseif grwalktype == :fractal_combined2
        walk =
            (x, y, z) -> fractal_combined2(
                x,
                MonomialOrder(S, S[1, :], [0]),
                MonomialOrder(T, T[1, :], T[1, :]),
            )
    end

    ######TODO:Check the parameter#####
    R = base_ring(G)
    I = Singular.Ideal(R, [R(x) for x in gens(G)])

    Gb = walk(I, S, T)
    #println("Cones crossed: ", getCounter())

    S = change_order(Gb.base_ring, T)
    return Singular.Ideal(S, [change_ring(gen, S) for gen in gens(Gb)])
end


function standard_walk22(G::Singular.sideal, S::Matrix{Int64}, T::Matrix{Int64})
    #println("standard_walk22 results")
    #println("Crossed Cones in: ")
    standard_walk22(G, S, T, S[1, :], T[1, :])
end

function standard_walk22(
    G::Singular.sideal,
    S::Matrix{Int64},
    T::Matrix{Int64},
    cweight::Vector{Int64},
    tweight::Vector{Int64},
)
    R = base_ring(G)
    Rn = change_order(R, cweight, T)
    terminate = false
    while !terminate
        G = standard_step(G, R, cweight, Rn)
        #println(cweight)
        #global counter = getCounter() + 1
        if cweight == tweight
            terminate = true
        else
            cweight = next_weight(G, cweight, tweight)
            R = Rn
            Rn = change_order(Rn, cweight, T)
        end
    end
    return G
end

function standard_step(
    G::Singular.sideal,
    R::Singular.PolyRing,
    cw::Vector{Int64},
    Rn::Singular.PolyRing
)
    Gw = initials(Rn, gens(G), cw)
    H = Singular.std(
        Singular.Ideal(Rn, Gw),
        complete_reduction = true,
    )
    #H = liftGW2(G, R, Gw, H, Rn)
    H = lift(G, R, H, Rn)
    return Singular.std(H, complete_reduction = true)
end

###############################################################
#Generic-version of the groebner walk by Fukuda et al. (2007)
###############################################################

function generic_walk2(G::Singular.sideal, S::Matrix{Int64}, T::Matrix{Int64})
    R = base_ring(G)
    Rn = change_order(G.base_ring, T)
    v = next_gamma(G, [0], S, T)
    Lm = [change_ring(Singular.leading_term(g), Rn) for g in gens(G)]
    G = Singular.Ideal(Rn, [change_ring(x, Rn) for x in gens(G)])

    #println("generic_walk2 results")
    #println("Crossed Cones with facetNormal: ")
    while !isempty(v)
        #global counter = getCounter() + 1
        #println(v)
        G, Lm = generic_step(G, Lm, v, T,R)
        v = next_gamma(G, Lm, v, S, T)
    end
    return Singular.interreduce(G)
end

function generic_step(
    G::Singular.sideal,
    Lm::Vector{Singular.spoly{L}},
    v::Vector{Int64},
    T::Matrix{Int64},
    R::Singular.PolyRing
) where {L<:Nemo.RingElem}

    Rn = Singular.base_ring(G)

    facet_Generators = facet_initials(G,Lm, v)
    H = Singular.std(
        Singular.Ideal(Rn, facet_Generators),
        complete_reduction = true,
    )
    H, Lm = lift_generic(G, Lm, H)
    G = interreduce(H, Lm)
    G = Singular.Ideal(Rn, G)
    G.isGB = true
    return G, Lm
end


###############################################################
#Pertubed-version of the groebner walk Amrhein et al.
###############################################################
function pertubed_walk2(
    G::Singular.sideal,
    S::Matrix{Int64},
    T::Matrix{Int64},
    p::Int64,
)
    #cweight = pertubed_vector(G, S, p)
    cweight = S[1, :]
    terminate = false
    #println("pertubed_walk2 results")
    #println("Crossed Cones in: ")

    while !terminate
        tweight = pertubed_vector(G, T, p)
        G = standard_walk22(G, S, T, cweight, tweight)
        if inCone(G, T, tweight)
            terminate = true
        else
            if p == 1
                R = change_order(G.base_ring, T)
                G = Singular.Ideal(R, [change_ring(x, R) for x in gens(G)])
                G = Singular.std(G, complete_reduction = true)
                terminate = true
            end
            p = p - 1
            cweight = tweight
        end
    end
    return G
end

###############################################################
#fractal-walk by Amrhein et al.
#Inlcuding:
#fractal_walk2 -> standard-version
#fractal_walk22 -> checks if the starting weight is in the inner of a cone.
#fractal_walk23 -> fractal walk expecially for conversion to the lexikographic orderig.
#                 checks if the starting weight is in the inner of a cone.
###############################################################

########################################
#Counter for the steps in the fractal_walk2
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
function raiseCounterFr()
    global counterFr = getCounterFr() + 1
end
PertVecs = []
sigma = []

function fractal_walk2(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
)
    global PertVecs = [pertubed_vector(G, T, i) for i = 1:nvars(base_ring(G))]
    #println(PertVecs)
    #println("FacrtalWalk_standard results")
    #println("Crossed Cones in: ")
    Gb = fractal_recursiv2(G, S,T, PertVecs, 1)
    #println("Cones crossed: ", deleteCounterFr())
    return Gb
end

function fractal_recursiv2(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    PertVecs::Vector{Vector{Int64}},
    p::Int64,
)
    R = base_ring(G)
    terminate = false
    G.isGB = true
    w = S.w

    while !terminate
        t = nextT(G, w, PertVecs[p])
        if (t == [0])
            if inCone(G, T,PertVecs[p])
                return G
            else
                global PertVecs = [pertubed_vector(G, T, i) for i = 1:nvars(R)]
                continue
            end
        end
        w = w + t * (PertVecs[p] - w)
        w = convert_bounding_vector(w)
        T.w = w
        Rn = change_order(R, T)
        Gw = initials(R, Singular.gens(G), w)
        if p == nvars(R)
            H = Singular.std(
                Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]),
                complete_reduction = true,
            )
            #println(w, " in depth", p)
            #raiseCounterFr()
        else
            #println("up in: ", p, " with: ", w)
            H = fractal_recursiv2(
                Singular.Ideal(R, [x for x in Gw]),
                S,
                T,
                PertVecs,
                p + 1,
            )
        end
        H = liftGW2(G, R, Gw, H, Rn)
        #H = lift(G, R, H, Rn)
        G = Singular.std(H, complete_reduction = true)
        R = Rn
    end
    return G
end

cwPert = []
firstStepMode = false
function cwpert(p::Int64)
    cwPert[p]
end

function fractal_walk2_start_order(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
)
    global PertVecs =
        [pertubed_vector(G, T, i) for i = 1:nvars(Singular.base_ring(G))]
    global sigma = S.w
    #println("fractal_walk2_withStartorder results")
    #println("Crossed Cones in: ")
    Gb = fractal_walk2_recursiv_startorder2(G, S, T, PertVecs, 1)
    #println("Cones crossed: ", deleteCounterFr())
    return Gb
end

function fractal_walk2_recursiv_startorder2(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    PertVecs::Vector{Vector{Int64}},
    p::Int64,
)
    R = Singular.base_ring(G)
    terminate = false
    G.isGB = true
    if (p == 1)
        if !isMonomial(initials(R, Singular.gens(G), S.w))
            global cwPert = [pertubed_vector(G, S, S.w, i) for i = 1:nvars(R)]
            global firstStepMode = true
        end
    end
    if firstStepMode
        w = cwPert[p]
    else
        w = S.w
    end

    while !terminate
        t = nextT(G, w, PertVecs[p])
        if t == [0]
            if inCone(G, T, PertVecs[p])
                #println(PertVecs[p], " in depth", p)
                return G
            else
                global PertVecs = [pertubed_vector(G, T, i) for i = 1:nvars(R)]
                continue
            end
        end
        w = w + t * (PertVecs[p] - w)
        w = convert_bounding_vector(w)
        T.w = w
        Rn = change_order(R, T)
        Gw = initials(R, gens(G), w)
        if p == Singular.nvars(R)
            H = Singular.std(
                Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]),
                complete_reduction = true,
            )
            #println(w, " in depth", p)
            #raiseCounterFr()
        else
            #println("up in: ", p, " with: ", w)

            H = fractal_walk2_recursiv_startorder2(
                Singular.Ideal(R, [x for x in Gw]),
                S,
                T,
                PertVecs,
                p + 1,
            )
            global firstStepMode = false
        end
        H = liftGW2(G, R, Gw, H, Rn)
        #H = lift(G, R, H, Rn)
        G = Singular.std(H, complete_reduction = true)
        R = Rn
    end
    return G
end
function fractal_walk2_lex(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
)
    global PertVecs = [pertubed_vector(G, T, i) for i = 1:nvars(base_ring(G))]
    #println("fractal_walk2_lex results")
    #println("Crossed Cones in: ")
    Gb = fractal_walk2_recursive_lex2(G, S, T, PertVecs, 1)
    #println("Cones crossed: ", deleteCounterFr())
    return Gb
end

function fractal_walk2_recursive_lex2(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    PertVecs::Vector{Vector{Int64}},
    p::Int64,
)
    R = Singular.base_ring(G)
    terminate = false
    G.isGB = true
    w = S.w
    while !terminate
        t = nextT(G, w, PertVecs[p])
        if t == [0]
            if inCone(G, T, PertVecs[p])
                return G
            else
                global PertVecs =
                    [pertubed_vector(G, T, i) for i = 1:Singular.nvars(R)]
                    #println(PertVecs)
                continue
            end
        end
        if t == 1 && p==1
            return fractal_walk2_recursive_lex2(
                G,
                S,
                T,
                PertVecs,
                p + 1,
            )
        else
        w = w + t * (PertVecs[p] - w)
        w = convert_bounding_vector(w)
        T.w = w
        Rn = change_order(R, T)
        Gw = initials(R, Singular.gens(G), w)
        if p == Singular.nvars(R)
            H = Singular.std(
                Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]),
                complete_reduction = true,
            )
            #println(w, " in depth", p)
            #raiseCounterFr()
        else
            #println("up in: ", p, " with: ", w)
            H = fractal_walk2_recursive_lex2(
                Singular.Ideal(R, [x for x in Gw]),
                S,
                T,
                PertVecs,
                p + 1,
            )
            global firstStepMode = false
        end
    end
        H = liftGW2(G, R, Gw, H, Rn)
        #H = lift(G, R, H, Rn)
        G = Singular.std(H, complete_reduction = true)
        R = Rn
    end
    return G
end
function fractal_walk2_look_ahead(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
)
    #println("fractal_walk2_look_ahead results")
    #println("Crossed Cones in: ")
    global PertVecs = [pertubed_vector(G, T, i) for i = 1:nvars(base_ring(G))]
    Gb = fractal_walk2_look_ahead_recursiv2(G, S, T, PertVecs, 1)
    #println("Cones crossed: ", deleteCounterFr())
    return Gb
end

function fractal_walk2_look_ahead_recursiv2(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    PertVecs::Vector{Vector{Int64}},
    p::Int64,
)
    R = Singular.base_ring(G)
    terminate = false
    G.isGB = true
    w = S.w

    while !terminate
        t = nextT(G, w, PertVecs[p])
        if t == [0]
            if inCone(G, T, PertVecs[p])
                return G
            else
                global PertVecs = [pertubed_vector(G, T, i) for i = 1:nvars(R)]
                continue
            end
        end
        w = w + t * (PertVecs[p] - w)
        w = convert_bounding_vector(w)
        T.w = w
        Rn = change_order(R, T)
        Gw = initials(R, Singular.gens(G), w)
        if (p == Singular.nvars(R) || isbinomial(Gw))
            H = Singular.std(
                Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]),
                complete_reduction = true,
            )
            #println(w, " in depth", p)
            #raiseCounterFr()
        else
            #println("up in: ", p, " with: ", w)
            H = fractal_walk2_look_ahead_recursiv2(
                Singular.Ideal(R, Gw),
                S,
                T,
                PertVecs,
                p + 1,
            )
        end

        H = liftGW2(G, R, Gw, H, Rn)
        #H = lift(G, R H, Rn)
        G = Singular.std(H, complete_reduction = true)
        R = Rn
    end
    return G
end




function fractal_walk2_combined(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
)
    global PertVecs =
        [pertubed_vector(G, T, i) for i = 1:nvars(Singular.base_ring(G))]
    #println("fractal_walk2_withStartorder results")
    #println("Crossed Cones in: ")
    Gb = fractal_walk2_combined(G, S, T, PertVecs, 1)
    #println("Cones crossed: ", deleteCounterFr())
    return Gb
end

function fractal_walk2_combined(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    PertVecs::Vector{Vector{Int64}},
    p::Int64,
)
    R = Singular.base_ring(G)
    terminate = false
    G.isGB = true
    if (p == 1)
        if !isMonomial(initials(R, Singular.gens(G), S.w))
            global cwPert = [pertubed_vector(G, S, S.w, i) for i = 1:nvars(R)]
            global firstStepMode = true
        end
    end
    if firstStepMode
        w = cwPert[p]
    else
        w = S.w
    end
    while !terminate
        t = nextT(G, w, PertVecs[p])
        if t == [0]
            if inCone(G, T, PertVecs[p])
                #println(PertVecs[p], " in depth", p)
                return G
            else
                global PertVecs = [pertubed_vector(G, T, i) for i = 1:nvars(R)]
                continue
            end
        end
        if t == 1 && p==1
            return fractal_walk2_combined(
                G,
                S,
                T,
                PertVecs,
                p + 1,
            )
        else
        w = w + t * (PertVecs[p] - w)
        w = convert_bounding_vector(w)
        T.w = w
        b = w
        Rn = change_order(R, T)
        Gw = initials(R, gens(G), w)
        if (p == Singular.nvars(R) || isbinomial(Gw))
            H = Singular.std(
                Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]),
                complete_reduction = true,
            )
            #println(w, " in depth", p)
            #raiseCounterFr()
        else
            #println("up in: ", p, " with: ", w)
            H = fractal_walk2_combined(
                Singular.Ideal(R, [x for x in Gw]),
                S,
                T,
                PertVecs,
                p + 1,
            )
            global firstStepMode = false
        end
    end
        H = liftGW2(G, R, Gw, H, Rn)
        #H = lift(G, R, H, Rn)
        G = Singular.std(H, complete_reduction = true)
        R = Rn
    end
    return G
end
###############################################################
#Tran-version of the groebner walk by Tran (2002)
###############################################################

function tran_walk2(G::Singular.sideal, S::Matrix{Int64}, T::Matrix{Int64})
    cweight = S[1, :]
    tweight = T[1, :]
    #println("tran_walk2 results")
    #println("Crossed Cones in: ")
    R = base_ring(G)
    if !isMonomial(initials(R, Singular.gens(G), cweight))
        cweight = pertubed_vector(G, S, nvars(R))
    end

    terminate = false
    while !terminate
        w = next_weight(G, cweight, tweight)
        if tryparse(string(w), Int32) == nothing
            #println("w bigger than int32")
            return G
        end
        Rn= change_order(R, w, T)
        if w == tweight
            if inCone(G, T, cweight)
                return G
            else
                if inSeveralCones(initials(base_ring(G), gens(G), w))
                    tweight = representation_vector(G, T)
                    continue
                end
            end
        end
        G = standard_step(G, R, w, Rn)
        #global counter = getCounter() + 1
        #println(w)
        R = Rn
        cweight = w
    end
end
