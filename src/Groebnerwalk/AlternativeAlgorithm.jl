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
P = Singular.sideal
function checkPvecs()
    global PVecs
    println(PVecs)
end

sigma = []
function alternative_algorithm_top(
    G::Singular.sideal,
    S::MonomialOrder{Matrix{Int64},Vector{Int64}},
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
)
    println("FractalWalk_lookahead results")
    println("Crossed Cones in: ")
    R = base_ring(G)
    global PVecs = [pert_Vectors(G, T, i) for i = 1:nvars(R)]

    Gb = alternative_algorithm(G, S, T, PVecs, 1)
    println("Cones crossed: ", deleteCounterFr())
    return Gb
end
function alternative_algorithm(
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
        println("d")
        w = nextW_2(G, cweight, PVecs[p])
        println("e")
        if (w == [0])
            if inCone(G, T, PVecs[p])
                println("f")
                if raiseCounterFr()
                    println(PVecs[p], " in depth", p)
                end
                return G
            else
                println("g")
                global PVecs = [pert_Vectors(G, T, 2, i) for i = 1:nvars(R)]
                #checkPvecs()
                continue
            end
        end
        println("c")
        T.w = w
        Rn, V = change_order(R, T)
        Gw = initials(R, Singular.gens(G), w)
        bin = isBinomial(Gw)
        if p > 3 || bin
            if bin
                Gnew = Singular.std(
                    Singular.Ideal(Rn, [change_ring(x, Rn) for x in Gw]),
                    complete_reduction = true,
                )
                println(w, " in depth", p)
                raiseCounterFr(true)
            elseif p > 3
                Gnew = generic_walk(Singular.Ideal(R, [x for x in Gw]), S.m, T.m)
                println(w, " in depth", p)
                raiseCounterFr(true)
            end
            else
                println("up in: ", p, " with: ", w)
                resetCounterFr()
                Gnew = alternative_algorithm(
                    Singular.Ideal(R, [x for x in Gw]),
                    S,
                    T,
                    PVecs,
                    p + 1,
                )
        end
        G = liftGWfr(G, Gnew, R, Rn)
        println("k")
        G = Singular.std(G, complete_reduction = true)
        println("a")
        R = Rn
        println("b")
        S = T
        cweight = w
    end
    return G
end

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
