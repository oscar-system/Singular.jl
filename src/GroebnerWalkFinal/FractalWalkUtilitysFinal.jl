include("GroebnerWalkUtilitysFinal.jl")
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
    R, V = change_order(G.base_ring, T.t, T.m)
    I = Singular.Ideal(R, [change_ring(x, R) for x in Singular.gens(G)])
    cvzip = zip(Singular.gens(I), initials(R, Singular.gens(I), t))
    for (g, ing) in cvzip
        if !isequal(Singular.leading_term(g), Singular.leading_term(ing))
            return false
        end
    end
    return true
end
#return a copy of the PolynomialRing I, equipped with the ordering represented by T.
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
end

function lift_fractal_walk(
    G::Singular.sideal,
    Gw::Singular.sideal,
    R::Singular.PolyRing,
    S::Singular.PolyRing,
)
    G.isGB = true
    rest = [
        change_ring(gen, S) -
        change_ring(Singular.reduce(change_ring(gen, R), G), S) for
        gen in Singular.gens(Gw)
    ]
    G = Singular.Ideal(S, [S(x) for x in rest])
    G.isGB = true
    return G
end

#returns ´true´ if all polynomials of the array are monomials.
function isMonomial(Gw::Vector{spoly{L}}) where {L<:Nemo.RingElem}
    for g in Gw
        if size(collect(Singular.coefficients(g)))[1] > 1
            return false
        end
    end
    return true
end

#returns ´true´ if all polynomials of the array are binomial or less.
function isbinomial(Gw::Vector{spoly{L}}) where {L<:Nemo.RingElem}
    for g in Gw
        if size(collect(Singular.coefficients(g)))[1] > 2
            return false
        end
    end
    return true
end

#Computes next weight vector. Version used in Amrhein Gloor.
function next_weight(
    G::Singular.sideal,
    w::Array{T,1},
    tw::Array{K,1},
) where {T<:Number,K<:Number}
    if (w == t)
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
        w = w + tmin * (t - w)
    else
        return [0]
    end
    return convert_bounding_vector(w)
end

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
