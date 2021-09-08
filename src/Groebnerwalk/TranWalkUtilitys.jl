include("GroebnerWalkUtilitys.jl")


function representationVector(G::Singular.sideal, T::Matrix{Int64})
    n = size(T)[1]
    M = 0
    for i in 1:n
        for j in 1:n
    temp = T[i,j]
    if M < temp
        M = temp
    end
end
end

    d0 = 0
    for g in Singular.gens(G)
#        println(g, " in tedeg", tdeg(g))
    temp = tdeg(g,n)
    if d0 < temp
        d0 = temp
    end
end
d = M * (2*d0^2 + (n + 1)* d0)
w = d^(n-1) * T[1,:]
for i in 2:n
    w = w + d^(n-i) * T[i,:]
end
return w
end
