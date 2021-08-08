using Oscar

###############################################################
#Utilitys for generic_walk
###############################################################

function facet_initials(
    R::MPolyRing,
    G::Singular.sideal,
    v::Vector{Int64},
    lm::Vector{Singular.spoly{Singular.n_unknown{fmpq}}},
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
    for i in 1:length(u)
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
    Lm::Vector{Singular.spoly{Singular.n_unknown{fmpq}}},
    I::Singular.sideal,
)
    S = base_ring(G)
    count = 1
    Newlm = Array{Singular.spoly{Singular.n_unknown{fmpq}},1}(undef, 0)
    liftArray = Array{fmpq_mpoly,1}(undef, 0)
    for g in gens(I)
        diff =
            change_ring(g, S) - reducegeneric_recursiv(change_ring(g, S), G, Lm)
        if diff != 0
            push!(Newlm, Oscar.Singular.leading_monomial(g))
            push!(liftArray, diff)
        end
    end
    return check_zeros(liftArray, Newlm)
end

function check_zeros(
    G::Vector{fmpq_mpoly},
    Lm::Vector{Singular.spoly{Singular.n_unknown{fmpq}}},
)
    Lm = filter(x -> x != 0, Lm)
    G = filter(x -> x != 0, G)
    return G, Lm
end

function nextV(
    G::Singular.sideal,
    Lm::Vector{Singular.spoly{Singular.n_unknown{fmpq}}},
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

function reducegeneric_recursiv(
    p::T,
    G::Singular.sideal,
    Lm::Vector{Singular.spoly{Singular.n_unknown{fmpq}}},
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
            if r == 0
                push!(Gh, (i, q))
                break
            end
        end
    end
    if !isempty(Gh)
        for (i, q) in Gh
            return reducegeneric_recursiv(
                p - q * gens(ideal(S, gens(G)))[i],
                G,
                Lm,
            )
        end
    else
        return p
    end
end
