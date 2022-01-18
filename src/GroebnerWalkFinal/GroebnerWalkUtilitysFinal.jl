###############################################################
#Utilitys for Groebnerwalks
###############################################################

#Computes next weight vector. Version used in Cox O´shea and Fukuda et al.
#This Version is used by the standard_walk, pertubed_walk and tran_walk.
function next_weight(
    G::Singular.sideal,
    cweight::Array{T,1},
    tweight::Array{K,1},
) where {T<:Number,K<:Number}
    tv = []
    for v in difference_lead_tail(G)
        cw = dot(cweight, v)
        tw = dot(tweight, v)
        if tw < 0
            push!(tv, cw // (cw - tw))
        end
    end
    push!(tv, 1)
    t = minimum(tv)
    w = (1 - t) * cweight + t * tweight
    return convert_bounding_vector(w)
end

function checkInt32(w::Vector{Int64})
for i = 1:length(w)
    if tryparse(Int32, string(w[i])) == nothing
        println("w bigger than int32")
        return false
    end
end
return true
end

#Return the initials of polynomials w.r.t. a weight vector.
function initials(R::Singular.PolyRing, G::Vector{spoly{n_Q}}, w::Vector{Int64})
    inits = spoly{elem_type(base_ring(R))}[]
    indexw = Tuple{Vector{Int},elem_type(base_ring(R))}[]
    for g in G
        empty!(indexw)
        maxw = 0
        eczip = zip(Singular.exponent_vectors(g), Singular.coefficients(g))
        for (e, c) in eczip
            tmpw = dot(w, e)
            if maxw == tmpw
                push!(indexw, (e, c))
            elseif maxw < tmpw
                empty!(indexw)
                push!(indexw, (e, c))
                maxw = tmpw
            end
        end
        inw = MPolyBuildCtx(R)
        for (e, c) in indexw
            Singular.push_term!(inw, c, e)
        end
        h = finish(inw)
        push!(inits, h)
    end
    return inits
end

#Return the difference of the exponents of the leading terms (Lm) and the
#exponent vectors of the tail of all polynomials of the ideal.
function difference_lead_tail(I::Singular.sideal)
    v = Vector{Int}[]
    for g in gens(I)
        ltu = Singular.leading_exponent_vector(g)
        for e in Singular.exponent_vectors(tail(g))
            push!(v, ltu .- e)
        end
    end
    return unique!(v)
end

function pertubed_vector(G::Singular.sideal, M::Matrix{Int64}, p::Integer)
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
    w = M[1, :] * e^(p - 1)
    for i = 2:p
        w += e^(p - i) * M[i, :]
    end
    return w
end

function inCone(G::Singular.sideal, T::Matrix{Int64}, t::Vector{Int64})
    R = change_order(G.base_ring, T)
    I = Singular.Ideal(R, [change_ring(x, R) for x in gens(G)])
    cvzip = zip(Singular.gens(I), initials(R, Singular.gens(I), t))
    for (g, ing) in cvzip
        if !isequal(Singular.leading_term(g), Singular.leading_term(ing))
            return false
        end
    end
    return true
end
#Fukuda et al
function lift(
    G::Singular.sideal,
    R::Singular.PolyRing,
    H::Singular.sideal,
    Rn::Singular.PolyRing
)
    G.isGB = true
    rest = [
        gen - change_ring(Singular.reduce(change_ring(gen, R), G), Rn) for
        gen in gens(H)
    ]
    G = Singular.Ideal(Rn, [Rn(x) for x in rest])
    G.isGB = true
    return G
end
#Amrhein  & Gloor
function liftGW2(
    G::Singular.sideal,
    R::Singular.PolyRing,
    inG::Vector{spoly{L}},
    H::Singular.sideal,
    Rn::Singular.PolyRing
) where {L<:Nemo.RingElem}

    gH = collect(gens(H))
    gG = collect(gens(G))
    s = length(inG)
    for i = 1:length(gH)
        q = divalg(change_ring(gH[i], R), [change_ring(x, R) for x in inG], R)
        gH[i] = Rn(0)
        for j = 1:s
            gH[i] = change_ring(gH[i], Rn) + change_ring(q[j],Rn) * change_ring(gG[j],Rn)
        end
    end
    G = Singular.Ideal(Rn, [x for x in gH])
    G.isGB = true
    return G
end

function divalg(
    p::spoly{L},
    f::Vector{spoly{L}},
    R::Singular.PolyRing,
) where {L<:Nemo.RingElem}
    s = length(f)
    q = Array{Singular.elem_type(R),1}(undef, s)
    for i = 1:s
        q[i] = R(0)
    end
    while !isequal(p, R(0))
        i = 1
        div = false
        while (div == false && i <= s)
            b, m = divides(leading_term(p), leading_term(f[i]))
            if b
                q[i] = q[i] + m
                p = p - (m * f[i])
                div = true
            else
                i = i + 1
            end
        end
            if div == false
                p = p - leading_term(p)
        end
    end
    return q
end

#Solves problems with weight vectors of floats.
function convert_bounding_vector(wtemp::Vector{T}) where {T<:Number}
    w = Vector{Int64}()
    for i = 1:length(wtemp)
        push!(w, float(divexact(wtemp[i], gcd(wtemp))))
    end
    return w
end

#return a copy of the PolynomialRing I, equipped with the ordering a(cweight)*ordering_M(T)
function change_order(
    R::Singular.PolyRing,
    cweight::Array{L,1},
    T::Matrix{Int64},
) where {L<:Number,K<:Number}
    G = Singular.gens(R)
    Gstrich = string.(G)
    S, H = Singular.PolynomialRing(
        R.base_ring,
        Gstrich,
        ordering = Singular.ordering_a(cweight) * Singular.ordering_M(T),
        cached = false,
    )
    return S
end

#return a copy of the PolynomialRing I, equipped with the ordering ordering_M(T)
function change_order(
    R::Singular.PolyRing,
    M::Matrix{Int64},
) where {T<:Number,K<:Number}
    G = Singular.gens(R)
    Gstrich = string.(G)
    S, H = Singular.PolynomialRing(
        R.base_ring,
        Gstrich,
        ordering = Singular.ordering_M(M),
        cached = false,
    )
    #@error("Not implemented yet")
    return S
end

function change_ring(p::Singular.spoly, R::Singular.PolyRing)
    cvzip = zip(Singular.coefficients(p), Singular.exponent_vectors(p))
    M = MPolyBuildCtx(R)
    for (c, v) in cvzip
        Singular.push_term!(M, c, v)
    end
    return finish(M)
end
function change_ring(p::Singular.spoly, R::Singular.PolyRing)
    cvzip = zip(Singular.coefficients(p), Singular.exponent_vectors(p))
    M = MPolyBuildCtx(R)
    for (c, v) in cvzip
        Singular.push_term!(M, c, v)
    end
    return finish(M)
end


#############################################
# unspecific help functions
#############################################

function ident_matrix(n::Int64)
    M = zeros(Int64, n, n)
    for i = 1:n
        M[i, i] = 1
    end
    return M
end

function anti_diagonal_matrix(n::Int64)
    M = zeros(Int64, n, n)
    for i = 1:n
        M[i, n+1-i] = -1
    end
    return M
end

# Singular.isequal depends on order of generators
function equalitytest(G::Singular.sideal, K::Singular.sideal)
    generators = Singular.gens(G)
    count = 0
    for gen in generators
        for r in Singular.gens(K)
            if gen - r == 0
                count += 1
            end
        end
    end
    if count == Singular.ngens(G)
        return true
    end
    return false
end

function dot(v::Vector{Int64}, w::Vector{Int64})
    n = length(v)
    sum = 0
    for i = 1:n
        sum += v[i] * w[i]
    end
    return sum
end

function ordering_as_matrix(w::Vector{Int64}, ord::Symbol)
    if length(w) > 2
        if ord == :lex
            return [
                w'
                ident_matrix(length(w))[1:length(w)-1, :]
            ]
        end
        if ord == :deglex
            return [
                w'
                ones(Int64, length(w))'
                ident_matrix(length(w))[1:length(w)-2, :]
            ]
        end
        if ord == :degrevlex
            return [
                w'
                ones(Int64, length(w))'
                anti_diagonal_matrix(length(w))[1:length(w)-2, :]
            ]
        end
        if ord == :revlex
            return [
                w'
                anti_diagonal_matrix(length(w))[1:length(w)-1, :]
            ]
        end
    else
        error("not implemented")
    end
end

function change_weight_vector(w::Vector{Int64}, M::Matrix{Int64})
    return [
        w'
        M[2:length(w), :]
    ]
end
function insert_weight_vector(w::Vector{Int64}, M::Matrix{Int64})
    return [
        w'
        M[1:length(w)-1, :]
    ]
end


function ordering_as_matrix(ord::Symbol, nvars::Int64)
    if ord == :lex
        return ident_matrix(nvars)
    end
    if ord == :deglex
        return [
            ones(Int64, nvars)'
            ident_matrix(nvars)[1:nvars-1, :]
        ]
    end
    if ord == :degrevlex
        return [
            ones(Int64, nvars)'
            anti_diagonal_matrix(nvars)[1:nvars-1, :]
        ]
    end
    if ord == :revlex
        return [
            w'
            anti_diagonal_matrix(length(w))[1:length(w)-1, :]
        ]
    else
        error("not implemented")
    end
end


function deg(p::Singular.spoly, n::Int64)
    max = 0
    for mon in Singular.monomials(p)
        ev = Singular.exponent_vectors(mon)
        sum = 0
        for e in ev
            for i = 1:n
                sum += e[i]
            end
        end
        if (max < sum)
            max = sum
        end
    end
    return max
end
