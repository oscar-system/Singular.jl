module libSingular

import Libdl
using CxxWrap

import ..Singular: libflint, libantic, libsingular_julia, AbstractAlgebra

@wrapmodule(libsingular_julia)

function __init__()
   @initcxx
   initialize_jl_c_types(@__MODULE__)
end

include("libsingular/LibSingularTypes.jl")

include("libsingular/coeffs.jl")

include("libsingular/rings.jl")

include("libsingular/matrices.jl")

include("libsingular/ideals.jl")

include("libsingular/resolutions.jl")

end # module
