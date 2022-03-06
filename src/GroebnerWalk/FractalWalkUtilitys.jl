include("GroebnerWalkUtilitys.jl")

#=
@doc Markdown.doc"""
function lift_fractal_walk(
G::Singular.sideal,
H::Singular.sideal,
Rn::Singular.PolyRing,
)
Performs a lifting step in the Fractal Walk.
"""=#
function lift_fractal_walk(
    G::Singular.sideal,
    H::Singular.sideal,
    Rn::Singular.PolyRing,
)
    R = base_ring(G)
    G.isGB = true
    G = Singular.Ideal(
        Rn,
        [
            change_ring(gen, Rn) -
            change_ring(Singular.reduce(change_ring(gen, R), G), Rn) for
            gen in Singular.gens(H)
        ],
    )
    G.isGB = true
    return G
end

#=
@doc Markdown.doc"""
function ismonomial(
Gw::Vector{spoly{L}},
) where {L<:Nemo.RingElem}
Returns ´true´ if all polynomials of the given array are monomials.
"""=#
function ismonomial(Gw::Vector{spoly{L}}) where {L<:Nemo.RingElem}
    for g in Gw
        if length(Singular.coefficients(g)) > 1
            return false
        end
    end
    return true
end

#=
@doc Markdown.doc"""
function isbinomial(
Gw::Vector{spoly{L}},
)
Returns ´true´ if all polynomials of the given array are binomials or less.
"""=#
function isbinomial(Gw::Vector{spoly{L}}) where {L<:Nemo.RingElem}
    for g in Gw
        if length(Singular.coefficients(g)) > 2
            return false
        end
    end
    return true
end


#=
@doc Markdown.doc"""
function nextT(
    G::Singular.sideal,
    w::Array{T,1},
    tw::Array{K,1},
) where {T<:Number,K<:Number}
Returns the next t to compute the next weight vector $w(t) = w + t * (tw - w)$ like it´s done in Amrhein & Gloor (1998)
"""=#
function nextT(
    G::Singular.sideal,
    w::Array{T,1},
    tw::Array{K,1},
) where {T<:Number,K<:Number}
    if (w == tw)
        return [0]
    end
    tmin = 2
    t = 0
    for g in gens(G)
        a = Singular.leading_exponent_vector(g)
        d = Singular.exponent_vectors(tail(g))
        for v in d
            frac = (dot(w, a) - dot(w, v) + dot(tw, v) - dot(tw, a))
            if frac != 0
                t = (dot(w, a) - dot(w, v)) // frac
            end
            if t > 0 && t < tmin
                tmin = t
            end
        end
    end
    if tmin <= 1
        return tmin
    else
        return [0]
    end
end

#=
@doc Markdown.doc"""
function next_weightfr(
    G::Singular.sideal,
    cweight::Array{T,1},
    tweight::Array{K,1},
) where {T<:Number,K<:Number}
Returns the next t to compute the next weight vector $w(t) = w + t * (tw - w)$ like it´s done in Fukuda et al. (2005).
"""=#
function next_weightfr(
    G::Singular.sideal,
    cweight::Array{T,1},
    tweight::Array{K,1},
) where {T<:Number,K<:Number}
    if (cweight == tweight)
        return [0]
    end
    tmin = BigInt(1)
    for v in difference_lead_tail(G)
        cw = BigInt(dot(cweight, v))
        tw = BigInt(dot(tweight, v))
        if tw < 0
            t = cw // (cw - tw)
            if t < tmin
                tmin = t
            end
        end
    end
    return tmin
end

#=
@doc Markdown.doc"""
function inCone(
    G::Singular.sideal,
    T::Matrix{Int},
    pvecs::Vector{Vector{Int}},
    p::Int,
)
Returns 'true' if the leading tems of $G$ w.r.t the matrixordering $T$ are the same as the leading terms of $G$ w.r.t the weighted monomial ordering with weight vector of pvecs[p] (pvecs[p-1]) and the Matrixordering $T$.
"""=#
function inCone(
    G::Singular.sideal,
    T::Matrix{Int},
    pvecs::Vector{Vector{Int}},
    p::Int,
)
    R = change_order(G.base_ring, T)
    I = Singular.Ideal(R, [change_ring(x, R) for x in gens(G)])
    cvzip = zip(
        Singular.gens(I),
        initials(R, Singular.gens(I), pvecs[p-1]),
        initials(R, Singular.gens(I), pvecs[p]),
    )
    for (g, in, in2) in cvzip
        if !isequal(Singular.leading_term(g), Singular.leading_term(in)) ||
           !isequal(Singular.leading_term(g), Singular.leading_term(in2))
            return false
        end
    end
    return true
end
