using Oscar

export groebner_walk


###############################################################
#Implementation of the gröbner walk.
#TODO: Implement pertubed version.
#TODO: Implement generic version.
#TODO: Implement fractal version.
###############################################################


###############################################################
#Top-Level
#TODO: Implement input checks
#Fragen:
#Singular.reduce für f % G selbst programmieren?
#Interreduction selbst endgültig programmieren programmieren?
#Für divrem in Oscar-Ideal: sinnvoller möglich?
#Arithmetik von Polynomen -> Immer MPolyBuildCTX?
###############################################################

function groebner_walk(
    G::Singular.sideal,
    S::Matrix{Int64},
    T::Matrix{Int64},
    grwalktype::Symbol = :standard,
)

    if grwalktype == :standard
        walk = (x, y, z) -> standard_walk(x, y, z)
    elseif grwalktype == :generic
        walk = (x, y, z) -> generic_walk(x, y, z)
    elseif grwalktype == :pertubed
        walk = (x, y, z) -> pertubed_walk(x, y, z)
    elseif grwalktype == :fractal
        walk = (x, y, z) -> fractal_walk(x, y, z)
        #error("Choose a strategy from: :generic, :standard ...")
        # fill the list or implement a functionality to choose the best
        # way w.r.t. the type of the ideal
    end

    ######TODO:Check the parameter#####
    R = singular_ring(base_ring(G))
    I = Singular.Ideal(R, [R(x) for x in gens(G)])

    return walk(I, S, T)
end


###############################################################
#Standard-version of the groebner walk by Collart et al
###############################################################
function standard_walk(G::Singular.sideal, S::Matrix{Int64}, T::Matrix{Int64})
    cweight = S[1, :]
    tweight = T[1, :]
    G = standard_step(I, cweight, tweight, T)

    #loop
    while cweight != tweight
        cweight = nextw(G, cweight, tweight)
        G = standard_step(G, cweight, tweight, T)
        #@info "Current Gröbnerbase: " G
    end

    #finalization
    S, V = change_order(G, T)
    return Oscar.Singular.Ideal(S, [change_ring(gen, S) for gen in gens(G)])

end

function standard_step(
    G::Singular.sideal,
    cw::Array{Int64,1},
    tw::Array{Int64,1},
    T::Matrix{Int64},
)
    R = base_ring(G)
    S, V = change_order(G, cw, T)
    inwG = initials(S, gens(G), cw)
    IinwG = Singular.std(
        Oscar.Singular.Ideal(S, [S(x) for x in inwG]),
        complete_reduction = true,
    )

    #Lifting to GB of new cone
    rest = [
        gen - change_ring(Oscar.Singular.reduce(change_ring(gen, R), G), S)
        for gen in gens(IinwG)
    ]
    Gnew = Oscar.Singular.Ideal(S, [S(x) for x in rest])
    Gnew.isGB = true

    return Oscar.Singular.std(Gnew, complete_reduction = true)
end

###############################################################
#Generic-version of the groebner walk by Fukuda et al. (2007)
###############################################################


function generic_walk(G::Singular.sideal, S::Matrix{Int64}, T::Matrix{Int64})
    R, V = change_order(G, T)

    v = nextV(G, [0], S, T)
    lm = [Singular.leading_monomial(g) for g in gens(G)]

    while !isempty(v)
        G, lm = generic_step(G, R, v, T, lm)
        @info "Current Gröbnerbase: " G
        v = nextV(G,lm, v, S, T)
    end
    G = Oscar.Singular.Ideal(R, [change_ring(x, R) for x in gens(G)])
    return G
end

function generic_step(
    G::Singular.sideal,
    S::MPolyRing,
    v::Vector{Int64},
    T::Matrix{Int64},
    lm::Vector{Singular.spoly{Singular.n_unknown{fmpq}}},
)
    R = base_ring(G)
    facet_Generators = facet_initials(S, G, v, lm)

    facet_Ideal = Singular.std(
        Oscar.Singular.Ideal(S, [S(x) for x in facet_Generators]),
        complete_reduction = true,
    )

    println("Facet Ideal: ", facet_Ideal)

    liftArray, Newlm = liftgeneric(G, lm, facet_Ideal)
    Gnew = Oscar.Singular.Ideal(S, [S(x) for x in liftArray])

    #println("New lifted GB:   ", Gnew, "  with LM:  ", Newlm)
    Gnew = interreduce_new(Gnew, Newlm)
    Gnew = Oscar.Singular.Ideal(R, [change_ring(x, R) for x in gens(Gnew)])
    Gnew.isGB = true
    return Gnew, Newlm
end


###############################################################
#Fractal-version of the groebner walk by Amrhein & Gloor (200)
###############################################################

function fractal_walk(
    G::Singular.sideal,
    cweight::Array{Int,1},
    tweight::Array{Int,1},
    tord::Symbol,
)

    T = ordering_as_Matrix(tweight, tord)
    M  = nextw(G, cweight, tweight)
    R = base_ring(G)
    S, V = change_order(G, cw, tw, ordering)

    inwG = initials(S, gens(G), cw)

    IinwG = Singular.std(
        Oscar.Singular.Ideal(S, [S(x) for x in inwG]),
        complete_reduction = true,
    )
    G.isGB = true

    #Lifting to GB of new cone
    rest = [
        gen - change_ring(Oscar.Singular.reduce(change_ring(gen, R), G), S)
        for gen in gens(IinwG)
    ]
    Gnew = Oscar.Singular.Ideal(S, [S(x) for x in rest])
    Gnew.isGB = true
    return std(Gnew, complete_reduction = true)
end
