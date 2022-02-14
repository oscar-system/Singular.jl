include("GroebnerWalkUtilitysFinal.jl")

###############################################################
#Utilitys for generic_walk
###############################################################

#=
@doc Markdown.doc"""
function facet_initials(
    G::Singular.sideal,
    lm::Vector{spoly{L}},
    v::Vector{Int},
) where {L<:Nemo.RingElem}
Returns the facet initials of the polynomials w.r.t. the vector v.
"""=#
function facet_initials(
    G::Singular.sideal,
    lm::Vector{spoly{L}},
    v::Vector{Int},
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

#=
@doc Markdown.doc"""
function difference_lead_tail(
    I::Singular.sideal,
    Lm::Vector{spoly{L}},
) where {L<:Nemo.RingElem}
Returns the differences of the exponent vectors of the leading terms and the polynomials of the generators of I.
"""=#
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

#=
@doc Markdown.doc"""
function isParallel(u::Vector{Int}, v::Vector{Int})
Returns true if the vector u is parallel to the vector v.
"""=#
function isParallel(u::Vector{Int}, v::Vector{Int})
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

#=
@doc Markdown.doc"""
function lift_generic(
    G::Singular.sideal,
    Lm::Vector{Singular.spoly{L}},
    H::Singular.sideal,
) where {L<:Nemo.RingElem}
Performs a lifting step in the Groebner Walk proposed by Fukuda et. al.
"""=#
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

function filter_btz(S::Matrix{Int}, V::Vector{Vector{Int}})
    btz = Set{Vector{Int}}()
    for v in V
        if bigger_than_zero(S, v)
            push!(btz, v)
        end
    end
    return btz
end

function filter_ltz(S::Matrix{Int}, V::Set{Vector{Int}})
    btz = Set{Vector{Int}}()
    for v in V
        if less_than_zero(S, v)
            push!(btz, v)
        end
    end
    return btz
end
function filter_lf(
    w::Vector{Int},
    S::Matrix{Int},
    T::Matrix{Int},
    V::Set{Vector{Int}},
)
    btz = Set{Vector{Int}}()
    for v in V
        if less_facet(w, v, S, T)
            push!(btz, v)
        end
    end
    return btz
end

function next_gamma(
    G::Singular.sideal,
    Lm::Vector{spoly{L}},
    w::Vector{Int},
    S::Matrix{Int},
    T::Matrix{Int},
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

function next_gamma(
    G::Singular.sideal,
    w::Vector{Int},
    S::Matrix{Int},
    T::Matrix{Int},
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

function bigger_than_zero(M::Matrix{Int}, v::Vector{Int})
    for i = 1:size(M)[1]
        d = dot(M[i, :], v)
        if d != 0
            return d > 0
        end
    end
    return false
end

function less_than_zero(M::Matrix{Int}, v::Vector{Int})
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
    u::Vector{Int},
    v::Vector{Int},
    S::Matrix{Int},
    T::Matrix{Int},
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

#=
@doc Markdown.doc"""
function dividesGW(p::Singular.spoly, lm::Singular.spoly, S::Singular.PolyRing)
Returns the multiple m for all terms q in p with lm * m = q.
"""=#
function dividesGW(p::Singular.spoly, lm::Singular.spoly, S::Singular.PolyRing)
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

#=
@doc Markdown.doc"""
function modulo(
    p::Singular.spoly,
    G::Vector{spoly{L}},
    Lm::Vector{spoly{L}},
    c::Bool = false,
) where {L<:Nemo.RingElem}
Returns p modulo G w.r.t. the leading terms Lm and true if a division occured.
"""=#

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
        (q, b) = dividesGW(p, Lm[i], R)
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
#=
function modulo(
    p::Singular.spoly,
    G::Vector{spoly{L}},
    Lm::Vector{spoly{L}},
    c::Bool = false,
) where {L<:Nemo.RingElem}
    I = 0
    R = parent(p)
    Q = zero(R)
    result = p
    b = true
    r = false
    while b
        b = false
        for i = 1:length(G)
            (q, b) = dividesGW(result, Lm[i], R)
            if b
                result = result - (q * G[i])
                b = true
                r = true
                break
            end
        end
    end
    return result, r
end=#
#=
@doc Markdown.doc"""
function interreduce(
    G::Vector{spoly{L}},
    Lm::Vector{spoly{L}},
) where {L<:Nemo.RingElem}
G represents a Gröbnerbasis. This function interreduces G w.r.t. the leading terms Lm with tail-reduction.
"""=#
function interreduce(
    G::Vector{spoly{L}},
    Lm::Vector{spoly{L}},
) where {L<:Nemo.RingElem}
    Rn = parent(first(G))
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
            G[i] = r
        end
    end
    return G
end
