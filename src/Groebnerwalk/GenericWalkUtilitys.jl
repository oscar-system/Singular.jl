using Oscar
include("GroebnerWalkUtilitys.jl")

###############################################################
#Utilitys for generic_walk
###############################################################

function facet_initials(
    R::MPolyRing,
    G::Singular.sideal,
    v::Vector{Int64},
    lm::Vector{Singular.spoly{Singular.n_FieldElem{fmpq}}},
)
    inits = []
    count = 1
    for g in gens(G)
        inw = MPolyBuildCtx(R)
        #V = filter(x -> less_facet(S, x), diff_vectors(G)
        el = collect(Singular.exponent_vectors(lm[count]))[1]

        for i = 1:length(g)
            e = collect(Singular.exponent_vectors(g))[i]
            if el == e || isparallel(el - e, v)
                push_term!(inw, collect(Singular.coefficients(g))[i], e)
            end
        end

        #push_term!(inw, collect(Singular.coefficients(lm[count]))[1], collect(Singular.exponent_vectors(lm[count]))[1])
        h = finish(inw)
        push!(inits, h)
        count = count + 1

    end
    #@info "Facet Initials: " inits
    return inits
end

function isparallel(u::Vector{Int64}, v::Vector{Int64})
    #TODO: maybe false this way
    count = 1
    x = 0
    for i = 1:length(u)
        if u[i] == 0
            if v[count] == 0
                count = count + 1
            else
                return false
            end
        else
            x = v[count] // u[i]
            count = count + 1
            break
        end
    end
    if count > length(v)
        return true
    end
    for i = count:length(v)
        if v[i] != x * u[i]
            return false
        end
    end
    return true
end

function liftgeneric(
    G::Singular.sideal,
    Lm::Vector{Singular.spoly{Singular.n_FieldElem{fmpq}}},
    I::Singular.sideal,
)
    S = base_ring(G)
    count = 1
    Newlm = Array{Singular.spoly{Singular.n_FieldElem{fmpq}},1}(undef, 0)
    liftArray = Array{Singular.spoly{Singular.n_FieldElem{fmpq}},1}(undef, 0)
    for g in gens(I)
        diff =
            subtractTerms(change_ring(g, S),
            S(reducegeneric_recursiv2(change_ring(g, S), G, Lm)), nvars(S))
    #=    println("red 2:",
            change_ring(g, S) -
            S(reducegeneric_recursiv2(change_ring(g, S), G, Lm)),
        )
        println("red 1:",
            change_ring(g, S) -
            S(reducegeneric_recursiv(change_ring(g, S), G, Lm)),
        ) =#
        if diff != 0
            push!(Newlm, Singular.leading_monomial(g))
            push!(liftArray, diff)
        end
    end
    return check_zeros(liftArray, Newlm)
end

function check_zeros(
    G::Vector{Singular.spoly{Singular.n_FieldElem{fmpq}}},
    Lm::Vector{Singular.spoly{Singular.n_FieldElem{fmpq}}},
)
    Lm = filter(x -> x != 0, Lm)
    G = filter(x -> x != 0, G)
    return G, Lm
end

function nextV(
    G::Singular.sideal,
    Lm::Vector{Singular.spoly{Singular.n_FieldElem{fmpq}}},
    w::Vector{Int64},
    S::Matrix{Int64},
    T::Matrix{Int64},
)
    V = filter(x -> bigger_than_zero(S, x), diff_vectors(G, Lm))
    #@info "#1 Current V: " V

    filter!(x -> less_than_zero(T, x), V)
    #@info "#2 Current V: " V

    if (w != [0])
        filter!(x -> less_facet(w, x, S, T), V)
    end
    #@info "#3 Current V: " V
    if isempty(V)
        return V
    end


    minV = V[1]
    for i = 2:length(V)
        if less_facet(V[i], minV, S, T)
            minV = V[i]
        end
    end

    #@info "#4 Current V: " minV
    return minV
end
function nextV(
    G::Singular.sideal,
    w::Vector{Int64},
    S::Matrix{Int64},
    T::Matrix{Int64},
)
    V = filter(x -> bigger_than_zero(S, x), diff_vectors(G))
    #@info "#1 Current V: " V

    filter!(x -> less_than_zero(T, x), V)
    #@info "#2 Current V: " V

    if (w != [0])
        filter!(x -> less_facet(w, x, S, T), V)
    end
    #@info "#3 Current V: " V
    if isempty(V)
        return V
    end


    minV = V[1]
    for i = 2:length(V)
        if less_facet(V[i], minV, S, T)
            minV = V[i]
        end
    end

    #@info "#4 Current V: " minV
    return minV
end


function bigger_than_zero(M::Matrix{Int64}, v::Vector{Int64})
    for i = 1:size(M)[1]
        d = dot(M[i, :], v)
        if (d != 0)
            return d > 0
        end
    end
    return false
end


function less_than_zero(M::Matrix{Int64}, v::Vector{Int64})
    for i = 1:size(M)[1]
        d = dot(M[i, :], v)
        if (d != 0)
            return d < 0
        end
    end
    return false
end

function less_facet(
    u::Vector{Int64},
    v::Vector{Int64},
    S::Matrix{Int64},
    T::Matrix{Int64},
)
    for i = 1:size(T)[1]
        for j = 1:size(S)[1]
            Tuv = dot(T[i, :], u) * dot(S[j, :], v)
            Tvu = dot(T[i, :], v) * dot(S[j, :], u)
            if Tuv != Tvu
                return Tuv < Tvu
            end
        end
    end
    return false
end
function reduce_modulo(
    p::Singular.spoly{Singular.n_FieldElem{fmpq}},
    lm::Singular.spoly{Singular.n_FieldElem{fmpq}},
    S::MPolyRing,
)
    a = []
    r = 0
    div = false
    newpoly = Singular.MPolyBuildCtx(S)
    for term in Singular.terms(p)
        #println(term, "and p", lm)
        (b, c) = Singular.divides(change_ring(term, S), change_ring(lm, S))
        #println("mult ",c)
        if b
            #=    if c == 1
                    push_term!(newpoly, 1, zeros(Int64, 1, nvars(S)))
                else =#
            push_term!(
                newpoly,
                collect(Singular.coefficients(c))[1],
                collect(Singular.exponent_vectors(c))[1],
            )
            #    end
            div = true
        end
    end
    if div
        np = finish(newpoly)
        return (np, true)
    end
    return (nothing, false)
end
changedInReduction = false
function hasChangedInReduction()
    global changedInReduction
    temp = changedInReduction
    global changedInReduction = false
    return temp
end
function reducegeneric_recursiv2(
    p::Singular.spoly,
    G::Singular.sideal,
    Lm::Vector{Singular.spoly{Singular.n_FieldElem{fmpq}}},
) where {T<:RingElement}
    R = base_ring(G)

    S, V = Singular.PolynomialRing(
        base_ring(R).base_ring,
        [String(s) for s in symbols(R)],
        ordering = :lex,
    )
        Gh = []
    for i = 1:ngens(G)
        (q, b) = reduce_modulo(change_ring(p, S), change_ring(Lm[i], S), S)
        if b
            push!(Gh, (i, q))
            break
        end
    end
    if !isempty(Gh)
        mul = [x for x in gens(G)]
        global changedInReduction = true
        for (i, q) in Gh
        #=    println("pretet ", p, "and", q * change_ring(mul[i], S))
            println(
                "test ",
                subtractTerms(change_ring(p, S), (q * change_ring(mul[i],S)), nvars(S)),
            ) =#
            return reducegeneric_recursiv2(
                subtractTerms(change_ring(p, S), (q * change_ring(mul[i],S)), nvars(S)),
                G,
                Lm,
            )
        end
    else
        return p
    end
end
function subtractTerms(p::Singular.spoly, q::Singular.spoly, dim::Int64)
    sol = MPolyBuildCtx(R)
    ep = collect(Singular.exponent_vectors(p))
    eq = collect(Singular.exponent_vectors(q))
    cp = collect(Singular.coefficients(p))
    cq = collect(Singular.coefficients(q))
    dummy = zeros(Int64, dim)
    dummy[1] = -1
    for j = 1:length(ep)
        fin = false
        for k = 1:length(eq)
            equals = true
            for i = 1:dim
                if ep[j][i] != eq[k][i]
                    equals = false
                    break
                end
            end
            if equals
                diff = cp[j] - cq[k]
                if diff != 0
                push_term!(
                    sol,
                    diff,
                    ep[j],
                )
            end
                fin = true
                eq[k] = dummy #delete Vector if equal
                break
            end
        end
        if !fin
            push_term!(sol, cp[j], ep[j])
        end
    end
        for i in 1:length(eq)
            if eq[i][1] != -1
        push_term!(sol, - cq[i], eq[i])
end
    end
    h = finish(sol)
    return h
end
#=
function reducegeneric_recursiv(
    p::T,
    G::Singular.sideal,
    Lm::Vector{Singular.spoly{Singular.n_FieldElem{fmpq}}},
) where {T<:RingElement}
    R = base_ring(G)
    S, V = PolynomialRing(
        base_ring(R).base_ring,
        [String(s) for s in symbols(R)],
        ordering = :lex,
    )
    Gh = []
    for i = 1:ngens(G)
        for mon in terms(first(gens(ideal(S, [p]))))
            (q, r) = divrem(mon, gens(ideal(S, Lm))[i])
            #printn("q ",q,"r ", r, "mon", mon, "p ", p)
            if r == 0
                push!(Gh, (i, q))
                #println(q)
                break
            end
        end
    end
    mul = gens(ideal(S, gens(G)))
    if !isempty(Gh)
        for (i, q) in Gh
            return reducegeneric_recursiv(
                first(gens(ideal(S, [p]))) - (q * mul[i]),
                G,
                Lm,
            )
        end
    else
        return p
    end
end
=#
function interreduce_new(
    G::Singular.sideal,
    Lm::Vector{Singular.spoly{Singular.n_FieldElem{fmpq}}},
)
    R = base_ring(G)
    gens = collect(Singular.gens(G))
    changed = true
    while changed
        changed = false
        for i = 1:ngens(G)
            #gensrest = filter(x -> x != gens[i], gens)
            gensrest =
                Array{Singular.spoly{Singular.n_FieldElem{fmpq}},1}(undef, 0)
            Lmrest =
                Array{Singular.spoly{Singular.n_FieldElem{fmpq}},1}(undef, 0)
            ## New function
            for j = 1:ngens(G)
                if i != j
                    push!(gensrest, gens[j])
                    push!(Lmrest, Lm[j])
                end
            end
            #Lmrest = filter(x -> x != Lm[i], Lm)
            gen = reducegeneric_recursiv2(
                gens[i],
                Singular.Ideal(R, gensrest),
                Lmrest,
            )
        #=    println("interred 2 :", reducegeneric_recursiv2(
                            gens[i],
                            Singular.Ideal(R, gensrest),
                            Lmrest,
                        ))
            println("interred 1 :", reducegeneric_recursiv(
                                        gens[i],
                                        Singular.Ideal(R, gensrest),
                                        Lmrest,
                                    )) =#
            if hasChangedInReduction()
                changed = true
                gens[i] = first(Singular.gens(Singular.Ideal(R, change_ring(gen, R))))
                break
            end
        end
    end
    return Oscar.Singular.Ideal(R, [R(p) for p in gens])
end
