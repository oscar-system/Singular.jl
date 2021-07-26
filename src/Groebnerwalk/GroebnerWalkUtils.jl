
using Oscar





#Returns the initials of polynomials w.r.t. a weight vector.
#The ordering doesn´t affect this.
function initials(
    R::MPolyRing,
    G::Array{T,1},
    w::Array{K,1},
) where {T<:RingElement,K<:Number}
    inits = []
    for g in G
        maxw = 0
        indexw = []
        e = collect(Singular.exponent_vectors(g))
        for i = 1:length(g)
            tmpw = dot(w, e[i])
            if maxw == tmpw
                push!(indexw, (i, e[i]))
                #rethink this. gens are preordered
            elseif maxw < tmpw
                indexw = []
                push!(indexw, (i, e[i]))
                maxw = tmpw
            end
        end
        inw = MPolyBuildCtx(R)
        for (k, j) in indexw
            push_term!(inw, collect(Singular.coefficients(g))[k], j)
        end
        h = finish(inw)
        push!(inits, h)
    end
    return inits
end

function diff_vectors(
    I::Singular.sideal,
    Lm::Vector{Singular.spoly{Singular.n_unknown{fmpq}}},
    )
    v = []
    for i in 1:ngens(I)
        ltu = Singular.leading_exponent_vector(Lm[i])
        for e in filter(x -> ltu != x, collect(Singular.exponent_vectors(gens(I)[i])))
            push!(v, ltu .- e)
        end
    end
    return unique!(v)
end

function diff_vectors(I::Singular.sideal)
    v = []
    for g in gens(I)
        ltu = Singular.leading_exponent_vector(g)
        for e in Singular.exponent_vectors(tail(g))
            push!(v, ltu .- e)
        end
    end
    return unique!(v)
end

###############################################################
#TODO: Change T instead of using a()
###############################################################
function change_order(
    I::Singular.sideal,
    cweight::Array{L,1},
    T::Matrix{Int64},
) where {L<:Number,K<:Number}
    R = I.base_ring
    G = Singular.gens(I.base_ring)
    Gstrich = string.(G)
    S, H = Oscar.Singular.PolynomialRing(
        R.base_ring,
        Gstrich,
        ordering = Oscar.Singular.ordering_a(cweight) *
                   Oscar.Singular.ordering_M(T),
    )
    return S, H
end

function change_order(
    I::Singular.sideal,
    M::Matrix{Int64},
) where {T<:Number,K<:Number}
    R = I.base_ring
    G = Singular.gens(I.base_ring)
    Gstrich = string.(G)
    S, H = Oscar.Singular.PolynomialRing(
        R.base_ring,
        Gstrich,
        ordering = Oscar.Singular.ordering_M(M),
    )
    #@error("Not implemented yet")
    return S, H
end

function change_ring(p::Singular.spoly, R::Singular.PolyRing)
    cvzip = zip(coefficients(p), exponent_vectors(p))
    M = MPolyBuildCtx(R)
    for (c, v) in cvzip
        push_term!(M, c, v)
    end
    return finish(M)
end
function change_ring(p::Singular.spoly, R::MPolyRing)
    cvzip = zip(coefficients(p), exponent_vectors(p))
    M = MPolyBuildCtx(R)
    for (c, v) in cvzip
        push_term!(M, c, v)
    end
    return finish(M)
end

function ordering_as_matrix(w::Vector{Int64}, ord::Symbol)
    if length(w) > 2
        if ord == :lex || ord == Symbol("Singular(lp)")
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
        if ord == :degrevlex || a.ord == Symbol("Singular(dp)")
            return [
                w'
                ones(Int64, length(w))'
                anti_diagonal_matrix(length(w))[1:length(w)-2, :]
            ]
        end
    else
        error("not implemented yet")
    end
end

function ordering_as_matrix(ord::Symbol, nvars::Int64)
    if ord == :lex || ord == Symbol("Singular(lp)")
        #return [w'; ident_matrix(length(w))[1:length(w)-1,:]]
        return ident_matrix(nvars)
    end
    if ord == :deglex
        return [
            ones(Int64, nvars)'
            ident_matrix(nvars)[1:nvars-1, :]
        ]
    end
    if ord == :degrevlex || a.ord == Symbol("Singular(dp)")
        return [
            ones(Int64, nvars)'
            anti_diagonal_matrix(nvars)[1:nvars-1, :]
        ]
    end
end



#Not needed at the moment. Use this instead of a(),ordering_M
function current_Order_of_Cone(v::Vector{Int64}, T::Matrix{Int64})
    M = [
        v'
        T[1:length(v)-1, :]
    ]
    return M
end

function pert_Vectors(M::Matrix{Int64}, p::Integer)

end

function tdeg(p::Singular.spoly)
    max = 0
    for mon in monomials
        ev = collect(exponent_vectors(mon))
        sum = 0
        for e in ev
            sum = e + sum
        end
        if (max < sum)
            max = sum
        end
    end
    return max
end

function interreduce_new(
    G::Singular.sideal,
    Lm::Vector{Singular.spoly{Singular.n_unknown{fmpq}}})
    R=base_ring(G)
    gens = collect(Singular.gens(G))
    changed = true
    while changed
        changed = false
    for i = 1:ngens(G)
        #gensrest = filter(x -> x != gens[i], gens)
        gensrest = Array{Singular.spoly{Singular.n_unknown{fmpq}},1}(undef, 0)
        Lmrest = Array{Singular.spoly{Singular.n_unknown{fmpq}},1}(undef, 0)
        ## New function
        for j = 1:ngens(G)
            if i != j
        push!(gensrest, gens[j])
        push!(Lmrest, Lm[j])
    end
    end
        #Lmrest = filter(x -> x != Lm[i], Lm)
        gen =reducegeneric_recursiv(gens[i], Singular.Ideal(R,gensrest), Lmrest)
        if gens[i] != gen
            changed = true
            gens[i]= first(Singular.gens(Singular.Ideal(R, R(gen))))
            break
        end
    end
end
    return Oscar.Singular.Ideal(R, [R(p) for p in gens])
end

#Use MPolybuildCTX
function vec_sum(p::Vector{fmpq_mpoly}, q::Vector{fmpq_mpoly})
    poly = 0
    for i = 1:length(p)
        poly = poly + p[i] * q[i]
    end
    return poly
end

# Singular.isequal depends on order of generators
function equalitytest(G::Singular.sideal, K::Singular.sideal)
    generators = gens(G)
    count = 0
    for gen in generators
        for r in gens(K)
            if gen - r == 0
                count = count + 1
            end
        end
    end
    if count == ngens(G)
        return true
    end
    return false
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
