###############################################################
#Utilitys for Groebnerwalks
###############################################################

#Computes next weight vector. Version used in Cox O´shea and Fukuda et al.
#This Version is used by the standard_walk, pertubed_walk and tran_walk.
function next_weight(
    G::Singular.sideal,
    cw::Vector{Int},
    tw::Vector{Int},
) where {K<:Number}
    tmin = BigInt(1)
    for v in difference_lead_tail(G)
        cdotw = BigInt(dot(cw, v))
        tdotw = BigInt(dot(tw, v))
        if tdotw < 0
            t = cdotw // (cdotw - tdotw)
            if t < tmin
                tmin = t
            end
        end
    end
    return convert_bounding_vector(cw + tmin * (tw - cw))
end

function checkInt32(w::Vector{Int})
    for i = 1:length(w)
        if tryparse(Int32, string(w[i])) == nothing
            println("int32")
            return false
        end
    end
    return true
end
function truncw(
    G::Singular.sideal,
    w::Vector{Int},
    inw::Vector{L},
) where {L<:Nemo.RingElem}
    while !checkInt32(w)
        println("$w > Int32")
        R = base_ring(G)
        for i = 1:length(w)
            w[i] = round(w[i] * 0.10)
        end
        w = convert_bounding_vector(w)
        if inw != initials(R, gens(G), w)
            println("Initials are different - return not rounded weight")
            return w, false
        end
    end
    println("Converted to Vector $w of the same face.")
    return w, true
end
#Return the initials of polynomials w.r.t. a weight vector.
function initials(
    R::Singular.PolyRing,
    G::Vector{L},
    w::Vector{Int},
) where {L<:Nemo.RingElem}
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

#=
@doc Markdown.doc"""
function pertubed_vector(
G::Singular.sideal,
M::Matrix{Int},
p::Integer
)
Computes a p-pertubed weight vector of M.
"""=#
function pertubed_vector(G::Singular.sideal, M::Matrix{Int}, p::Integer)
    m = Int[]
    n = size(M, 1)
    for i = 1:p
        max = M[i, 1]
        for j = 1:n
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
    return convert_bounding_vector(w)
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
function inCone(G::Singular.sideal, T::Matrix{Int}, t::Vector{Int})
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

#=
@doc Markdown.doc"""
function isGb(
    G::Singular.sideal,
    T::Matrix{Int},
)
"""=#
function isGb(G::Singular.sideal, T::Matrix{Int})
    R = change_order(G.base_ring, T)
    for g in Singular.gens(G)
        if !isequal(
            Singular.leading_term(g),
            change_ring(Singular.leading_term(change_ring(g, R)), G.base_ring),
        )
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
    Rn::Singular.PolyRing,
)
    G.isGB = true
    G = Singular.Ideal(
        Rn,
        [
            gen - change_ring(Singular.reduce(change_ring(gen, R), G), Rn)
            for gen in gens(H)
        ],
    )
    G.isGB = true
    return G
end

#=
@doc Markdown.doc"""
function lift_fractal_walk(
G::Singular.sideal,
R::Singular.PolyRing,
inG::Vector{spoly{L}},
H::Singular.sideal,
Rn::Singular.PolyRing,
)
Performs a lifting step in the Groebner Walk proposed by Amrhein et. al. and Cox Little Oshea
"""=#
function liftGW2(
    G::Singular.sideal,
    R::Singular.PolyRing,
    inG::Vector{spoly{L}},
    H::Singular.sideal,
    Rn::Singular.PolyRing,
) where {L<:Nemo.RingElem}

    gH = collect(gens(H))
    gG = collect(gens(G))
    inG = [change_ring(x, R) for x in inG]
    for i = 1:length(gH)
        q = division_algorithm(change_ring(gH[i], R), inG, R)
        gH[i] = Rn(0)
        for j = 1:length(inG)
            gH[i] =
                change_ring(gH[i], Rn) +
                change_ring(q[j], Rn) * change_ring(gG[j], Rn)
        end
    end
    G = Singular.Ideal(Rn, [x for x in gH])
    G.isGB = true
    return G
end

function division_algorithm(
    p::spoly{L},
    f::Vector{spoly{L}},
    R::Singular.PolyRing,
) where {L<:Nemo.RingElem}
    q = Array{Singular.elem_type(R),1}(undef, length(f))
    for i = 1:length(f)
        q[i] = R(0)
    end
    while !isequal(p, R(0))
        i = 1
        div = false
        while (div == false && i <= length(f))
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

#=
@doc Markdown.doc"""
   convert_bounding_vector(wtemp::Vector{T}) where {T<:Number}
Given a Vector{Number} $v$ this function computes a Vector{Int} w with w = v:gcd(v).
"""=#
function convert_bounding_vector(wtemp::Vector{T}) where {T<:Rational{BigInt}}
    w = Vector{Int}()
    g = gcd(wtemp)
    for i = 1:length(wtemp)
        push!(w, float(divexact(wtemp[i], g)))
    end
    return w
end
#=
@doc Markdown.doc"""
   convert_bounding_vector(wtemp::Vector{T}) where {T<:Number}
Given a Vector{Number} $v$ this function computes a Vector{Int} w with w = v:gcd(v).
"""=#
function convert_bounding_vector(wtemp::Vector{T}) where {T<:Number}
    w = Vector{Int}()
    g = gcd(wtemp)
    for i = 1:length(wtemp)
        push!(w, float(divexact(wtemp[i], g)))
    end
    return w
end


#return a copy of the PolynomialRing I, equipped with the ordering a(cw)*ordering_M(T)
function change_order(
    R::Singular.PolyRing,
    cw::Array{L,1},
    T::Matrix{Int},
) where {L<:Number,K<:Number}
    G = Singular.gens(R)
    Gstrich = string.(G)
    s = size(T)
    if s[1] == s[2]
        S, H = Singular.PolynomialRing(
            R.base_ring,
            Gstrich,
            ordering = Singular.ordering_a(cw) * Singular.ordering_M(T),
            cached = false,
        )
    elseif s[1] - s[2] == 1
        S, H = Singular.PolynomialRing(
            R.base_ring,
            Gstrich,
            ordering = Singular.ordering_a(cw) *
                       Singular.ordering_a(T[1, :]) *
                       Singular.ordering_M(T[2:end, :]),
            cached = false,
        )
    else
        S, H = Singular.PolynomialRing(
            R.base_ring,
            Gstrich,
            ordering = Singular.ordering_a(cw) *
                       Singular.ordering_a(T[1, :]) *
                       Singular.ordering_a(T[2, :]) *
                       Singular.ordering_M(T[3:end, :]),
            cached = false,
        )
    end
    return S
end

#return a copy of the PolynomialRing I, equipped with the ordering a(cw)*ordering_M(T)
function change_order(
    R::Singular.PolyRing,
    w::Vector{Int},
    t::Vector{Int},
    T::Matrix{Int},
) where {}
    G = Singular.gens(R)
    Gstrich = string.(G)
    s = size(T)
    S, H = Singular.PolynomialRing(
        R.base_ring,
        Gstrich,
        ordering = Singular.ordering_a(w) *
                   Singular.ordering_a(t) *
                   Singular.ordering_M(T),
        cached = false,
    )
    return S
end

#return a copy of the PolynomialRing I, equipped with the ordering ordering_M(T)
function change_order(
    R::Singular.PolyRing,
    M::Matrix{Int},
) where {T<:Number,K<:Number}
    G = Singular.gens(R)
    Gstrich = string.(G)
    S, H = Singular.PolynomialRing(
        R.base_ring,
        Gstrich,
        ordering = Singular.ordering_M(M),
        cached = false,
    )
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


#=
@doc Markdown.doc"""
function interreduce(
    G::Vector{spoly{L}},
    Lm::Vector{spoly{L}},
) where {L<:Nemo.RingElem}
G represents a Gröbnerbasis. This function interreduces G w.r.t. the leading terms Lm with tail-reduction.
"""=#
function interreduce_walk(G::Singular.sideal) where {L<:Nemo.RingElem}
    Rn = base_ring(G)
    Generator = collect(gens(G))
    I = 0
    for i = 1:ngens(G)
        I = Singular.Ideal(Rn, Generator[1:end.!=i])
        I.isGB = true
        Generator[i] = reduce(Generator[i], I)
    end
    G = Singular.Ideal(Rn, Generator)
    return G
end
