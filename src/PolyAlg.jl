include("poly/PolyTypes.jl")
include("poly/GPolyTypes.jl")
include("weyl/WeylTypes.jl")
include("exterior/ExteriorTypes.jl")

const polyalg{T} = Union{spoly{T}, pweyl{T}, pexterior{T}} where T <: Nemo.RingElem

include("poly/poly.jl")
include("poly/gpoly.jl")
include("weyl/weyl.jl")
include("exterior/exterior.jl")

