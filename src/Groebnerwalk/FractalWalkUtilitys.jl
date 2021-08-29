include("GroebnerWalkUtilitys.jl")

#Solves problems with weight vectors of floats.
function convertBoundingVector(wtemp::Vector{T}) where {T<:Number}
    w = Vector{Int64}()

    for i = 1:length(wtemp)
        push!(w, float(divexact(wtemp[i], gcd(wtemp))))
    end
    return w
end


function inCone(G::Singular.sideal, T::MonomialOrder{Matrix{Int64}, Vector{Int64}},t::Vector{Int64})
    R, V = change_order(G, T.t, T.m)
    I = Oscar.Singular.Ideal(R, [change_ring(x, R) for x in gens(G)])
    cvzip = zip(gens(I), initials(R, gens(I), t))
    for (g, ing) in cvzip
        if !isequal(Singular.leading_term(g), Singular.leading_term(ing))
            return false
        end
    end
    return true
end

function change_order(
    I::Singular.sideal,
    T::MonomialOrder{Matrix{Int64}, Vector{Int64}},
) where {L<:Number,K<:Number}
    R = I.base_ring
    G = Singular.gens(I.base_ring)
    Gstrich = string.(G)
        S, H = Oscar.Singular.PolynomialRing(
            R.base_ring,
            Gstrich,
            ordering = Oscar.Singular.ordering_a(T.w)*
                       Oscar.Singular.ordering_M(T.m),
        )
#=    S, H = Oscar.Singular.PolynomialRing(
        R.base_ring,
        Gstrich,
        ordering = Oscar.Singular.ordering_a(T.t) *Oscar.Singular.ordering_a(T.w)*
                   Oscar.Singular.ordering_M(T.m)
    )
    return S, H
    =#
end

function liftGWfr(G::Singular.sideal, Gw::Singular.sideal, R::Singular.PolyRing, S::Singular.PolyRing)
rest = [
    change_ring(gen, S) - change_ring(Oscar.Singular.reduce(change_ring(gen, R), G), S)
    for gen in gens(Gw)
]
G = Oscar.Singular.Ideal(S, [S(x) for x in rest])
G.isGB = true
return G
end

function ismonomial(Gw::Vector{Any})
    for g in Gw
        if size(collect(Singular.coefficients(g)))[1] > 1
            return false
        end
    end
    return true
end
function isBinomial(Gw::Vector{Any})
    for g in Gw
        if size(collect(Singular.coefficients(g)))[1] > 2
            return false
        end
    end
    return true
end
function nextw_fr(
    G::Singular.sideal,
    cweight::Array{T,1},
    tweight::Array{K,1},
) where {T<:Number,K<:Number}
if (cweight == tweight)
    return [0]
end
    tv = []
    for v in diff_vectors(G)
        tw = dot(tweight, v)
        if tw < 0
            push!(tv, dot(cweight, v) // (dot(cweight, v) - tw))
        end
    end
    push!(tv, 1)
    t = minimum(tv)
    if (0 == numerator(t))
        filter!(x -> (numerator(x) > 0 && x!=1), tv)
        if isempty(tv)
            return[0]
        else
        t = minimum(tv)
    end
    end
    if (t == 1)
        return tweight
    end
    w = (1 - t) * cweight + t * tweight

    return convertBoundingVector(w)
end

function diff_vectors_lt(I::Singular.sideal)
    zip = []
    v = []
    for g in gens(I)
        ltu = Singular.leading_exponent_vector(g)
        for e in Singular.exponent_vectors(tail(g))
            push!(v, (ltu .- e))
        end
        push!(zip, (unique(v), ltu))
        v = []
    end
    return zip
end

function nextW_2(G::Singular.sideal,
    w::Array{T,1},
    sigma::Array{K,1},
    ) where {T<:Number,K<:Number}
    if (w == sigma)
        return [0]
    end

tmin = 2
t = 0
    for g in gens(G)
        a = Singular.leading_exponent_vector(g)
        d = Singular.exponent_vectors(tail(g))
        for v in d
            frac = (dot(w,a)- dot(w,v) + dot(sigma,v) - dot(sigma, a))
            if frac != 0
            t = (dot(w,a) - dot(w,v)) // frac
        end
            if t > 0 && t < tmin
                tmin = t
            end
        end
    end

    if tmin <= 1
    w = w + tmin * (sigma - w)
else
    return [0]
end
    return convertBoundingVector(w)
        end
