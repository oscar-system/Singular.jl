include("GroebnerWalkUtilitysFinal.jl")


function RepresentationVector(G::Singular.sideal, T::Matrix{Int64})
    n = size(T)[1]
    M = 0
    for i = 1:n
        for j = 1:n
            temp = T[i, j]
            if M < temp
                M = temp
            end
        end
    end

    d0 = 0
    for g in Singular.gens(G)
        #        println(g, " in tedeg", tdeg(g))
        temp = tdeg(g, n)
        if d0 < temp
            d0 = temp
        end
    end
    d = M * (2 * d0^2 + (n + 1) * d0)
    w = d^(n - 1) * T[1, :]
    for i = 2:n
        w = w + d^(n - i) * T[i, :]
    end
    return w
end

function InSeveralCones(Gw::Vector{spoly{L}}) where {L<:Nemo.RingElem}
    counter = 0
    for g in Gw
        if size(collect(Singular.coefficients(g)))[1] > 2
            return true
        end
        if size(collect(Singular.coefficients(g)))[1] == 2
            counter = counter +1
        end
end
    if counter > 1
        return true
    end
    return false
end
