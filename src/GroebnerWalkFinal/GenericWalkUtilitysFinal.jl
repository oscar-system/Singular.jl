include("GroebnerWalkUtilitysFinal.jl")

###############################################################
#Utilitys for generic_walk
###############################################################

#Return the facet_initials of polynomials w.r.t. a weight vector.
function FacetInitials(
    R::Singular.PolyRing,
    G::Singular.sideal,
    v::Vector{Int64},
    lm::Vector{spoly{L}},
) where {L<:Nemo.RingElem}
    inits = []
    count = 1
    for g in Singular.gens(G)
        inw = Singular.MPolyBuildCtx(R)
        el = first(Singular.exponent_vectors(lm[count]))
        for (e, c) in zip(Singular.exponent_vectors(g), Singular.coefficients(g))
            if el == e || isParallel(el - e, v)
                Singular.push_term!(inw, c, e)
            end
        end
        h = finish(inw)
        push!(inits, h)
        count += 1
    end
    return inits
end

#
function isParallel(u::Vector{Int64}, v::Vector{Int64})
    count = 1
    x = 0
    for i = 1:length(u)
        if u[i] == 0
            if v[count] == 0
                count += +1
            else
                return false
            end
        else
            x = v[count] // u[i]
            count += 1
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

#lifting step of the generic_walk
function LiftGeneric(
    G::Singular.sideal,
    Lm::Vector{Singular.spoly{L}},
    H::Singular.sideal,
) where {L<:Nemo.RingElem}
    S = base_ring(G)
    Newlm = Array{Singular.elem_type(S),1}(undef, 0)
    liftPolys = Array{Singular.elem_type(S),1}(undef, 0)
    for g in Singular.gens(H)
        r, b = reduce_recursive(g, gens(G), Lm, S)
        diff = g - r
        if diff != 0
            push!(Newlm, Singular.leading_term(g))
            push!(liftPolys, diff)
        end
    end
    return liftPolys, Newlm, S
end

function filter_btz(S::Matrix{Int64}, V::Vector{Vector{Int64}})
    btz = Set{Vector{Int64}}()
    for v in V
        if bigger_than_zero(S, v)
            push!(btz, v)
        end
    end
    return btz
end

function filter_ltz(S::Matrix{Int64}, V::Set{Vector{Int64}})
    btz = Set{Vector{Int64}}()
    for v in V
        if less_than_zero(S, v)
            push!(btz, v)
        end
    end
    return btz
end
function filter_lf(
    w::Vector{Int64},
    S::Matrix{Int64},
    T::Matrix{Int64},
    V::Set{Vector{Int64}},
)
    btz = Set{Vector{Int64}}()
    for v in V
        if less_facet(w, v, S, T)
            push!(btz, v)
        end
    end
    return btz
end

#return the next facet_normal.
function NextGamma(
    G::Singular.sideal,
    Lm::Vector{spoly{L}},
    w::Vector{Int64},
    S::Matrix{Int64},
    T::Matrix{Int64},
) where {L<:Nemo.RingElem}
    V = filter_btz(S, LmMinusV(G, Lm))
    V = filter_ltz(T, V)
    if (w != [0])
        V = filter_lf(w, S, T, V)
    end
    if isempty(V)
        return V
    end
    minV = first(V)
    for v in V
        if less_facet(v, minV, S, T)
            minV = v
        end
    end
    return minV
end

#return the next facet_normal.
function NextGamma(
    G::Singular.sideal,
    w::Vector{Int64},
    S::Matrix{Int64},
    T::Matrix{Int64},
)
    V = filter_btz(S, LmMinusV(G))
    V = filter_ltz(T, V)
    if (w != [0])
        V = filter_lf(w, S, T, V)
    end
    if isempty(V)
        return V
    end
    minV = first(V)
    for v in V
        if less_facet(v, minV, S, T)
            minV = v
        end
    end
    return minV
end

function bigger_than_zero(M::Matrix{Int64}, v::Vector{Int64})
    for i = 1:size(M)[1]
        d = dot(M[i, :], v)
        if d != 0
            return d > 0
        end
    end
    return false
end

function less_than_zero(M::Matrix{Int64}, v::Vector{Int64})
    nrows, ncols = size(M)
    for i = 1:nrows
        d = 0
        for j = 1:ncols
            @inbounds d += M[i, j] * v[j]
        end
        if d != 0
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

#returns divrem()
function divrem(
    p::Singular.spoly,
    lm::Singular.spoly,
    S::Singular.PolyRing,
)
    div = false
    newpoly = Singular.MPolyBuildCtx(S)
    for term in Singular.terms(p)
        (b, c) = Singular.divides(term, lm)
        if b
            push_term!(
                newpoly,
                first(Singular.coefficients(c)),
                first(Singular.exponent_vectors(c)),
            )
            div = true
        end
    end
    return (finish(newpoly), div)
end

function reduce_recursive(
    p::Singular.spoly,
    G::Vector{spoly{L}},
    Lm::Vector{spoly{L}},
    R::Singular.PolyRing,
    c::Bool = false,
) where {L<:Nemo.RingElem}
    I = 0
    Q = zero(R)
    for i = 1:length(G)
        (q, b) = divrem(p, Lm[i], R)
        if b
            I = i
            Q = q
            break
        end
    end
    if I != 0
        r, b = reduce_recursive(p - (Q * G[I]), G, Lm, R)
        return r, true
    else
        return p, false
    end
end

function interreduce(
    G::Singular.sideal,
    Lm::Vector{spoly{L}},
) where {L<:Nemo.RingElem}
    R = Singular.base_ring(G)
    gens = collect(Singular.gens(G))
    changed = true
    while changed
        changed = false
        for i = 1:Singular.ngens(G)
            gensrest = Array{Singular.elem_type(R),1}(undef, 0)
            Lmrest = Array{Singular.elem_type(R),1}(undef, 0)
            for j = 1:Singular.ngens(G)
                if i != j
                    push!(gensrest, gens[j])
                    push!(Lmrest, Lm[j])
                end
            end
            r, b = reduce_recursive(gens[i], gensrest, Lmrest, R)
            if b
                changed = true
                gens[i] = r
                break
            end
        end
    end
    return Singular.Ideal(R, [R(p) for p in gens])
end
