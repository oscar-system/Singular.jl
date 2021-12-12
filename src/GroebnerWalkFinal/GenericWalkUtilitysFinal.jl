include("GroebnerWalkUtilitysFinal.jl")

###############################################################
#Utilitys for generic_walk
###############################################################

#Return the facet_facet_facet_initials of polynomials w.r.t. a weight vector.
function facet_initials(
    G::Singular.sideal,
    lm::Vector{spoly{L}},
    v::Vector{Int64},
) where {L<:Nemo.RingElem}
    Rn = base_ring(G)
    initials = Array{Singular.elem_type(Rn),1}(undef, 0)
    count = 1
    for g in Singular.gens(G)
        inw = Singular.MPolyBuildCtx(Rn)
        el = first(Singular.exponent_vectors(lm[count]))
        for (e, c) in
            zip(Singular.exponent_vectors(g), Singular.coefficients(g))
            if el == e || isParallel(el - e, v)
                Singular.push_term!(inw, c, e)
            end
        end
        h = finish(inw)
        push!(initials, h)
        count += 1
    end
    return initials
end

#Return the difference of the exponents of the leading terms (Lm) and the
#exponent vectors of the tail of all polynomials of the ideal.
function difference_lead_tail(
    I::Singular.sideal,
    Lm::Vector{spoly{L}},
) where {L<:Nemo.RingElem}
    v = Vector{Int}[]
    for i = 1:ngens(I)
        ltu = Singular.leading_exponent_vector(Lm[i])
        for e in Singular.exponent_vectors(gens(I)[i])
            if ltu != e
                push!(v, ltu .- e)
            end
        end
    end
    return unique!(v)
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
function lift_generic(
    G::Singular.sideal,
    Lm::Vector{Singular.spoly{L}},
    H::Singular.sideal,
) where {L<:Nemo.RingElem}
    Rn = base_ring(G)
    Newlm = Array{Singular.elem_type(Rn),1}(undef, 0)
    liftPolys = Array{Singular.elem_type(Rn),1}(undef, 0)
    for g in Singular.gens(H)
        r, b = modulo(g, gens(G), Lm)
        diff = g - r
        if diff != 0
            push!(Newlm, Singular.leading_term(g))
            push!(liftPolys, diff)
        end
    end
    return liftPolys, Newlm
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
function next_gamma(
    G::Singular.sideal,
    Lm::Vector{spoly{L}},
    w::Vector{Int64},
    S::Matrix{Int64},
    T::Matrix{Int64},
) where {L<:Nemo.RingElem}
    V = filter_btz(S, difference_lead_tail(G, Lm))
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
function next_gamma(
    G::Singular.sideal,
    w::Vector{Int64},
    S::Matrix{Int64},
    T::Matrix{Int64},
)
    V = filter_btz(S, difference_lead_tail(G))
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
function divremGW(p::Singular.spoly, lm::Singular.spoly, S::Singular.PolyRing)
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

function modulo(
    p::Singular.spoly,
    G::Vector{spoly{L}},
    Lm::Vector{spoly{L}},
    c::Bool = false,
) where {L<:Nemo.RingElem}
    I = 0
    R = parent(p)
    Q = zero(R)
    for i = 1:length(G)
        (q, b) = divremGW(p, Lm[i], R)
        if b
            I = i
            Q = q
            break
        end
    end
    if I != 0
        r, b = modulo(p - (Q * G[I]), G, Lm)
        return r, true
    else
        return p, false
    end
end

function interreduce(
    G::Vector{spoly{L}},
    Lm::Vector{spoly{L}},
) where {L<:Nemo.RingElem}
    Rn = parent(first(G))
    changed = true
    while changed
        changed = false
        for i = 1:Singular.length(G)
            gensrest = Array{Singular.elem_type(Rn),1}(undef, 0)
            Lmrest = Array{Singular.elem_type(Rn),1}(undef, 0)
            for j = 1:length(G)
                if i != j
                    push!(gensrest, G[j])
                    push!(Lmrest, Lm[j])
                end
            end
            r, b = modulo(G[i], gensrest, Lmrest)
            if b
                changed = true
                G[i] = r
                break
            end
        end
    end
    return G
end
