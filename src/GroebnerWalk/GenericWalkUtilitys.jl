include("GroebnerWalkUtilitys.jl")

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
    G::Vector{Singular.spoly{L}},
    lm::Vector{spoly{L}},
    v::Vector{Int},
) where {L<:Nemo.RingElem}
    Rn = parent(first(G))
    initials = Array{Singular.elem_type(Rn),1}(undef, 0)
    count = 1
    for g in G
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
    G::Vector{spoly{L}},
    Lm::Vector{spoly{L}},
) where {L<:Nemo.RingElem}
    v = Vector{Int}[]
    for i = 1:length(G)
        ltu = Singular.leading_exponent_vector(Lm[i])
        for e in Singular.exponent_vectors(G[i])
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
        @inbounds if v[i] != x * u[i]
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
    G::Vector{spoly{L}},
    Lm::Vector{Singular.spoly{L}},
    H::Singular.sideal,
) where {L<:Nemo.RingElem}
    Rn = parent(first(G))
    Newlm = Array{Singular.elem_type(Rn),1}(undef, 0)
    liftPolys = Array{Singular.elem_type(Rn),1}(undef, 0)
    for g in Singular.gens(H)
        push!(Newlm, Singular.leading_term(g))
        push!(liftPolys, g - modulo(g, G, Lm))
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
    G::Vector{spoly{L}},
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
#=
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
=#
function bigger_than_zero(M::Matrix{Int}, v::Vector{Int})
    nrows, ncols = size(M)
    for i = 1:nrows
        d = 0
        for j = 1:ncols
            @inbounds d += M[i, j] * v[j]
        end
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
    for i = 1:size(T, 1)
        for j = 1:size(S, 1)
            @inbounds Tuv = dot(T[i, :], u) * dot(S[j, :], v)
            @inbounds Tvu = dot(T[i, :], v) * dot(S[j, :], u)
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
    return finish(newpoly), div
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
) where {L<:Nemo.RingElem}
    for i = 1:length(G)
        (q, b) = dividesGW(p, Lm[i], parent(p))
        if b
            return modulo(p - (q * G[i]), G, Lm)
        end
    end
    return p
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
        G[i] = modulo(G[i], gensrest, Lmrest)
    end
    return G
end

function submult(
    p::Singular.spoly{L},
    q::Singular.spoly{L},
    c::Singular.spoly{L},
    dim::Int64,
) where {L<:Nemo.RingElem}
    sol = MPolyBuildCtx(parent(p))
    ep = collect(Singular.exponent_vectors(p))
    eq = collect(Singular.exponent_vectors(q))
    cp = collect(Singular.coefficients(p))
    cq = collect(Singular.coefficients(q))
    ec = collect(Singular.exponent_vectors(c))
    cc = collect(Singular.coefficients(c))
    multc = Array{RingElem,1}(undef, 0)
    multe = Array{Vector{Int},1}(undef, 0)

    for m = 1:length(cq)
        for n = 1:length(cc)
            skip = false
            so = eq[m] + ec[n]
            for i = 1:length(multe)
                if multe[i] == so
                    multc[i] = multc[i] + (cq[m] * cc[n])
                    skip = true
                end
            end
            if !skip
                push!(multe, so)
                push!(multc, cq[m] * cc[n])
            end
        end
    end
    for j = 1:length(ep)
        fin = true
        for k = 1:length(multe)
            equals = true
            if ep[j] != multe[k]
                equals = false
            end
            if equals
                diff = cp[j] - multc[k]
                if diff != 0
                    Singular.push_term!(sol, diff, ep[j])
                    diff = nothing
                end
                fin = false
                multe[k][1] = -1 #delete Vector if equal
                break
            end
        end
        if fin
            Singular.push_term!(sol, cp[j], ep[j])
        end
    end
    for i = 1:length(multe)
        if multe[i][1] != -1
            Singular.push_term!(sol, -multc[i], multe[i])
        end
    end
    ep = nothing
    eq = nothing
    cp = nothing
    cq = nothing
    ec = nothing
    cc = nothing
    multc = nothing
    multe = nothing
    return Singular.finish(sol)
end
