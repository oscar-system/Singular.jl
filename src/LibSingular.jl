module libSingular

import Libdl
using CxxWrap
@wrapmodule(realpath(joinpath(@__DIR__, "..", "local", "lib", "libsingularwrap." * Libdl.dlext)))

function __init__()
   @initcxx
end

include("libsingular/LibSingularTypes.jl")

include("libsingular/coeffs.jl")

include("libsingular/rings.jl")

include("libsingular/matrices.jl")

include("libsingular/ideals.jl")

include("libsingular/resolutions.jl")

end # module
