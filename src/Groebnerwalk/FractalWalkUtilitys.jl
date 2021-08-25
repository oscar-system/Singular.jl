include("GroebnerWalkUtilitys.jl")

#Solves problems with weight vectors of floats.
function convertBoundingVector(wtemp::Vector{T}) where {T<:Number}
    w = Vector{Int64}()
    for i = 1:length(wtemp)
        push!(w, float(divexact(wtemp[i], gcd(wtemp))))
    end
    return w
end

function nextw_fr(
    G::Singular.sideal,
    cweight::Array{T,1},
    tweight::Array{K,1},
) where {T<:Number,K<:Number}
    tv = []
    for v in diff_vectors(G)
        cw = dot(cweight, v)
        tw = dot(tweight, v)
        ctw = cw - tw
        if tw < 0
            push!(tv, cw // ctw)
        end
    end
    #tv = [dot(cweight, v) < 0 ? dot(cweight, v) / (dot(cweight, v) - dot(tweight, v)) : nothing for v = V ]
    if isempty(tv)
        return [0]
    end
    t = minimum(tv)
    if (0 == float(t))
        return [0]
    end
    w = (1 - t) * cweight + t * tweight
    return convertBoundingVector(w)
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
