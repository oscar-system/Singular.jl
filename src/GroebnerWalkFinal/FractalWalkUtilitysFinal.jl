include("GroebnerWalkUtilitysFinal.jl")

#Structure which is used to define a MonomialOrdering a(v)*a(tv)*ordering_M(T)
#Maybe not needed
mutable struct MonomialOrder{T<:Matrix{Int},v<:Vector{Int},tv<:Vector{Int}}
    m::T
    w::v
    t::tv
end
#=
#=
@doc Markdown.doc"""
   convert_bounding_vector(wtemp::Vector{T}) where {T<:Number}
Given a Vector{Number} $v$ this function computes a Vector{Int} w with w = v:gcd(v).
"""=#
function convert_bounding_vector(v::Vector{T}) where {T<:Number}
    w = Vector{Int}()
    for i = 1:length(v)
        push!(w, float(divexact(v[i], gcd(v))))
    end
    return w
end=#

#=
@doc Markdown.doc"""
function inCone(
    G::Singular.sideal,
    T::MonomialOrder{Matrix{Int},Vector{Int}},
    t::Vector{Int},
)
Returns 'true' if the leading tems of $G$ w.r.t the monomial ordering $T$ are the same as the leading terms of $G$ w.r.t the weighted monomial ordering with weight vector $t$ and the monomial ordering $T$.
"""=#
function inCone(
    G::Singular.sideal,
    T::MonomialOrder{Matrix{Int},Vector{Int}},
    t::Vector{Int},
)
    R = change_order(G.base_ring, T.m)
    I = Singular.Ideal(R, [change_ring(x, R) for x in Singular.gens(G)])
    cvzip = zip(Singular.gens(I), initials(R, Singular.gens(I), t))
    for (g, ing) in cvzip
        if !isequal(Singular.leading_term(g), Singular.leading_term(ing))
            return false
        end
    end
    return true
end

#=
@doc Markdown.doc"""
function lift_fractal_walk(
    G::Singular.sideal,
    Gw::Singular.sideal,
    R::Singular.PolyRing,
    S::Singular.PolyRing,
)
Performs a lifting step in the Groebner Walk proposed by Amrhein et. al. and Cox Little Oshea
"""=#
function lift_fractal_walk(
    G::Singular.sideal,
    R::Singular.PolyRing,
    H::Singular.sideal,
    Rn::Singular.PolyRing,
)
    G.isGB = true
    G = Singular.Ideal(Rn, [change_ring(gen, Rn) -
    change_ring(Singular.reduce(change_ring(gen, R), G), Rn) for
    gen in Singular.gens(H)])
    G.isGB = true
    return G
end

#=
@doc Markdown.doc"""
function isMonomial(
Gw::Vector{spoly{L}},
) where {L<:Nemo.RingElem}
Returns ´true´ if all polynomials of the given array are monomials.
"""=#
function isMonomial(Gw::Vector{spoly{L}}) where {L<:Nemo.RingElem}
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
Returns the next t to compute the next weight vector w(t) = w + t * (tw - w).
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
function change_order(
    R::Singular.PolyRing,
    T::MonomialOrder{Matrix{Int},Vector{Int}},
) where {L<:Number,K<:Number}
Returns a new PolynomialRing w.r.t. the monomial ordering T.
"""=#
function change_order(
    R::Singular.PolyRing,
    T::MonomialOrder{Matrix{Int},Vector{Int}},
) where {L<:Number,K<:Number}
    G = Singular.gens(R)
    Gstrich = string.(G)
    S, H = Singular.PolynomialRing(
        R.base_ring,
        Gstrich,
        ordering = Singular.ordering_a(T.w) *
                   Singular.ordering_a(T.t) *
                   Singular.ordering_M(T.m),
    )
    return S
end

function next_weightfr(
    G::Singular.sideal,
    cweight::Array{T,1},
    tweight::Array{K,1},
) where {T<:Number,K<:Number}
    if (cweight == tweight)
        return [0]
    end
    tmin = 1
    for v in difference_lead_tail(G)
        cw = dot(cweight, v)
        tw = dot(tweight, v)
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
    t::Vector{Int},
)
Returns 'true' if the leading tems of $G$ w.r.t the matrixordering $T$ are the same as the leading terms of $G$ w.r.t the weighted monomial ordering with weight vector $t$ and the Matrixordering $T$.
"""=#
function inCone(G::Singular.sideal, T::Matrix{Int}, pvecs::Vector{Vector{Int}}, p::Int)
    R = change_order(G.base_ring, T)
    I = Singular.Ideal(R, [change_ring(x, R) for x in gens(G)])
    cvzip = zip(Singular.gens(I), initials(R, Singular.gens(I), pvecs[p-1]), initials(R, Singular.gens(I), pvecs[p]))
    for (g, in, in2) in cvzip
        if !isequal(Singular.leading_term(g), Singular.leading_term(in)) || !isequal(Singular.leading_term(g), Singular.leading_term(in2))
            return false
        end
    end
    return true
end
