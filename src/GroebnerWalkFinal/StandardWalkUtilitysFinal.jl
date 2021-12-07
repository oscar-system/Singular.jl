include("GroebnerWalkUtilitysFinal.jl")

###############################################################
#Utilitys for standard_walk
###############################################################


#Solves problems with weight vectors of floats.
function convertBoundingVector(wtemp::Vector{T}) where {T<:Number}
    w = Vector{Int64}()
    for i = 1:length(wtemp)
        push!(w, float(divexact(wtemp[i], gcd(wtemp))))
    end
    return w
end

#Computes next weight vector. Version used in Cox O´shea and Fukuda et al.
#This Version is used by the standard_walk, pertubed_walk and tran_walk.
function NextWeight(
    G::Singular.sideal,
    cweight::Array{T,1},
    tweight::Array{K,1},
) where {T<:Number,K<:Number}
    tv = []
    for v in LmMinusV(G)
        cw = dot(cweight, v)
        tw = dot(tweight, v)
        if tw < 0
                push!(tv, cw // (cw - tw))
        end
    end
    push!(tv, 1)
    #filter!(x -> numerator(x) != 0, tv)
    t = minimum(tv)
    w = (1 - t) * cweight + t * tweight
    return convertBoundingVector(w)
end
