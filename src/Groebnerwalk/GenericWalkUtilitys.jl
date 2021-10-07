include("GroebnerWalkUtilitys.jl")

###############################################################
#Utilitys for generic_walk
###############################################################

#Return the facet_initials of polynomials w.r.t. a weight vector.
function facet_initials(
    R::Singular.PolyRing,
    G::Singular.sideal,
    v::Vector{Int64},
    lm::Vector{spoly{n_Q}},
)
    inits = []
    count = 1
    for g in Singular.gens(G)
        inw = Singular.MPolyBuildCtx(R)
        el = first(Singular.exponent_vectors(lm[count]))
        mzip = zip(Singular.exponent_vectors(g), Singular.coefficients(g))
        for (e, c) in mzip
            if el == e || isparallel(el - e, v)
                Singular.push_term!(inw, c, e)
            end
        end
        h = finish(inw)
        push!(inits, h)
        count += 1

    end
    #@info "Facet Initials: " inits
    return inits
end

#
function isparallel(u::Vector{Int64}, v::Vector{Int64})
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

#lifting step of the generic_walk
function liftgeneric(
    G::Singular.sideal,
    Lm::Vector{Singular.spoly{L}},
    H::Singular.sideal,
) where {L<:Nemo.RingElem}
    S = base_ring(G)
    Newlm = Array{Singular.elem_type(S),1}(undef, 0)
    liftArray = Array{Singular.elem_type(S),1}(undef, 0)
    for g in Singular.gens(H)
        r, b = reducegeneric_recursiv(g, gens(G), Lm, S)
        diff = g - r
        if diff != 0
            push!(Newlm, Singular.leading_term(g))
            push!(liftArray, diff)
        end
    end
    return liftArray, Newlm, S
end

#removes polynomials with zero coefficients.
function check_zeros(
    G::Vector{spoly{L}},
    Lm::Vector{spoly{L}},
    S::Singular.PolyRing,
) where {L<:Nemo.RingElem}

    lmfin = Array{Singular.elem_type(S),1}(undef, 0)
    Gfin = Array{Singular.elem_type(S),1}(undef, 0)
    for lm in Lm
        if lm != 0
            push!(lmfin, lm)
        end
    end
    for g in G
        if g != 0
            push!(Gfin, g)
        end
    end
    return Gfin, lmfin
end


function filter_btz(S::Matrix{Int64}, V::Vector{Any})
    btz = Set()
    for v in V
        if bigger_than_zero(S, v)
            push!(btz, v)
        end
    end
    return btz
end

function filter_ltz(S::Matrix{Int64}, V::Set{Any})
    btz = Set()
    for v in V
        if less_than_zero(S, v)
            push!(btz, v)
        end
    end
    return btz
end
function filter_lf(
    w::Vector{Int64},
    S::Matrix{Int64},
    T::Matrix{Int64},
    V::Set{Any},
)
    btz = Set()
    for v in V
        if less_facet(w, v, S, T)
            push!(btz, v)
        end
    end
    return btz
end

#return the next facet_normal.
function nextV(
    G::Singular.sideal,
    Lm::Vector{spoly{L}},
    w::Vector{Int64},
    S::Matrix{Int64},
    T::Matrix{Int64},
) where {L<:Nemo.RingElem}
    V = filter_btz(S, diff_vectors(G, Lm))
    #@info "#1 Current V: " V

    V = filter_ltz(T, V)
    #@info "#2 Current V: " V

    if (w != [0])
        V = filter_lf(w, S, T, V)
    end
    #@info "#3 Current V: " V
    if isempty(V)
        return V
    end

    minV = first(V)
    for v in V
        if less_facet(v, minV, S, T)
            minV = v
        end
    end

    #@info "#4 Current V: " minV
    return minV
end

#return the next facet_normal.
function nextV(
    G::Singular.sideal,
    w::Vector{Int64},
    S::Matrix{Int64},
    T::Matrix{Int64},
)
    V = filter_btz(S, diff_vectors(G))
    #@info "#1 Current V: " V

    V = filter_ltz(T, V)
    #@info "#2 Current V: " V

    if (w != [0])
        V = filter_lf(w, S, T, V)
    end
    #@info "#3 Current V: " V
    if isempty(V)
        return V
    end

    minV = first(V)
    for v in V
        if less_facet(v, minV, S, T)
            minV = v
        end
    end

    #@info "#4 Current V: " minV
    return minV
end

function bigger_than_zero(M::Matrix{Int64}, v::Vector{Int64})
    for i = 1:size(M)[1]
        d = dot(M[i, :], v)
        if d != 0
            return d > 0
        end
    end
    return false
end

function less_than_zero(M::Matrix{Int64}, v::Vector{Int64})
    nrows, ncols = size(M)
    for i = 1:nrows
        d = 0
        for j = 1:ncols
            @inbounds d += M[i, j] * v[j]
            #d = dot(M[i, :], v)
        end
        if d != 0
            return d < 0
        end
    end
    return false
end

function less_facet(
    u::Vector{Int64},
    v::Vector{Int64},
    S::Matrix{Int64},
    T::Matrix{Int64},
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

#returns p mod lm
function reduce_modulo(
    p::Singular.spoly,
    lm::Singular.spoly,
    S::Singular.PolyRing,
)
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

function reducegeneric_recursiv(
    p::Singular.spoly,
    G::Vector{spoly{L}},
    Lm::Vector{spoly{L}},
    R::Singular.PolyRing,
    c::Bool = false,
) where {L<:Nemo.RingElem}
    I = 0
    Q = zero(R)
    for i = 1:length(G)
        (q, b) = reduce_modulo(p, Lm[i], R)
        if b
            I = i
            Q = q
            break
        end
    end
    if I != 0
        r, b = reducegeneric_recursiv(p - (Q * G[I]), G, Lm, R)
        return r, true
    else
        return p, false
    end
end
function subtractTerms(
    p::Singular.spoly{L},
    q::Singular.spoly{L},
    dim::Int64,
) where {L<:Nemo.RingElem}
    sol = MPolyBuildCtx(parent(p))
    ep = collect(Singular.exponent_vectors(p))
    eq = collect(Singular.exponent_vectors(q))
    cp = collect(Singular.coefficients(p))
    cq = collect(Singular.coefficients(q))
    for j = 1:length(ep)
        fin = true
        for k = 1:length(eq)
            equals = true
            if ep[j] != eq[k]
                equals = false
            end
            if equals
                diff = cp[j] - cq[k]
                if diff != 0
                    Singular.push_term!(sol, diff, ep[j])
                end
                fin = false
                eq[k][1] = -1 #delete Vector if equal
                break
            end
        end
        if fin
            Singular.push_term!(sol, cp[j], ep[j])
        end
    end
    for i = 1:length(eq)
        if eq[i][1] != -1
            Singular.push_term!(sol, -cq[i], eq[i])
        end
    end
    h = Singular.finish(sol)
    return h
end

function interreduce_generic(
    G::Singular.sideal,
    Lm::Vector{spoly{L}},
) where {L<:Nemo.RingElem}
    R = Singular.base_ring(G)
    gens = collect(Singular.gens(G))
    changed = true
    while changed
        changed = false
        for i = 1:Singular.ngens(G)
            gensrest = Array{Singular.elem_type(R),1}(undef, 0)
            Lmrest = Array{Singular.elem_type(R),1}(undef, 0)
            for j = 1:Singular.ngens(G)
                if i != j
                    push!(gensrest, gens[j])
                    push!(Lmrest, Lm[j])
                end
            end
            r, b = reducegeneric_recursiv(gens[i], gensrest, Lmrest, R)
            if b
                changed = true
                gens[i] = r
                break
            end
        end
    end

    return Singular.Ideal(R, [R(p) for p in gens])
end
