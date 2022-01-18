include("GroebnerWalkUtilitysFinal.jl")

#Structure which is used to define a MonomialOrdering a(v)*a(tv)*ordering_M(T)
#Maybe not needed
mutable struct MonomialOrder{
    T<:Matrix{Int64},
    v<:Vector{Int64},
    tv<:Vector{Int64},
}
    m::T
    w::v
    t::tv
end
#=
@doc Markdown.doc"""
   convert_bounding_vector(wtemp::Vector{T}) where {T<:Number}
Given a Vector{Number} $v$ this function computes a Vector{Int64} w with w = v*gcd(v).
"""=#
function convert_bounding_vector(v::Vector{T}) where {T<:Number}
    w = Vector{Int64}()
    for i = 1:length(v)
        push!(w, float(divexact(v[i], gcd(v))))
    end
    return w
end

#=
@doc Markdown.doc"""
function inCone(
    G::Singular.sideal,
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    t::Vector{Int64},
)
Returns 'true' if the leading tems of $G$ w.r.t the monomial ordering $T$ are the same as the leading terms of $G$ w.r.t the weighted monomial ordering with weight vector $t$ and the monomial ordering $T$.
"""=#
function inCone(
    G::Singular.sideal,
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    t::Vector{Int64},
)
    R = change_order(G.base_ring, T.t, T.m)
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
Performs a lifting step in the Groebner Walk proposed by Amrhein et. al.
"""=#
function lift_fractal_walk(
    G::Singular.sideal,
    R::Singular.PolyRing,
    H::Singular.sideal,
    Rn::Singular.PolyRing,
)
    G.isGB = true
    rest = [
        change_ring(gen, Rn) -
        change_ring(Singular.reduce(change_ring(gen, R), G), Rn) for
        gen in Singular.gens(H)
    ]
    G = Singular.Ideal(Rn, [Rn(x) for x in rest])
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
        if size(collect(Singular.coefficients(g)))[1] > 1
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
        if size(collect(Singular.coefficients(g)))[1] > 2
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
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
) where {L<:Number,K<:Number}
Returns a new PolynomialRing w.r.t. the monomial ordering T.
"""=#
function change_order(
    R::Singular.PolyRing,
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
) where {L<:Number,K<:Number}
    G = Singular.gens(R)
    Gstrich = string.(G)
    S, H = Singular.PolynomialRing(
        R.base_ring,
        Gstrich,
        ordering = Singular.ordering_a(T.w) * Singular.ordering_M(T.m),
    )
    return S
end

#=
@doc Markdown.doc"""
function pertubed_vector(
G::Singular.sideal,
Mo::MonomialOrder{Matrix{Int64}},
t::Vector{Int64},
p::Integer,
)
Computes a p-pertubed weight vector of the current weight vector t by using the monomial ordering Mo.
"""=#
function pertubed_vector(
    G::Singular.sideal,
    Mo::MonomialOrder{Matrix{Int64}},
    t::Vector{Int64},
    p::Integer,
)
    if t == Mo.m[1, :]
        M = Mo.m
    else
        M = insert_weight_vector(t, Mo.m)
    end
    m = []
    n = size(M)[1]
    for i = 1:p
        max = M[i, 1]
        for j = 2:n
            temp = abs(M[i, j])
            if temp > max
                max = temp
            end
        end
        push!(m, max)
    end
    msum = 0
    for i = 2:p
        msum += m[i]
    end
    maxdeg = 0
    for g in gens(G)
        td = deg(g, n)
        if (td > maxdeg)
            maxdeg = td
        end
    end
    e = maxdeg * msum + 1
    w = view(M,1, :) * e^(p - 1)
    for i = 2:p
        w += e^(p - i) * view(M,i, :)
    end
    return w
end

#=
@doc Markdown.doc"""
function pertubed_vector(
G::Singular.sideal,
Mo::MonomialOrder{Matrix{Int64}},
t::Vector{Int64},
p::Integer,
)
Computes a p-pertubed weight vector of the current weight vector given by the first row of the matrix corresponding the the monomial ordering T.
"""=#
function pertubed_vector(
    G::Singular.sideal,
    T::MonomialOrder{Matrix{Int64},Vector{Int64}},
    p::Integer,
)
    m = []
    if T.t == T.m[1, :]
        M = T.m
    else
        M = insert_weight_vector(T.t, T.m)
    end
    n = size(M)[1]
    for i = 1:p
        max = M[i, 1]
        for j = 2:n
            temp = abs(M[i, j])
            if temp > max
                max = temp
            end
        end
        push!(m, max)
    end
    msum = 0
    for i = 2:p
        msum += m[i]
    end
    maxdeg = 0
    for g in gens(G)
        td = deg(g, n)
        if (td > maxdeg)
            maxdeg = td
        end
    end
    e = maxdeg * msum + 1
    w = view(M, 1, :) * e^(p - 1)
    for i = 2:p
        w += e^(p - i) * view(M, i, :)
    end
    return w
end
