module libSingular

using CxxWrap
@wrapmodule(realpath(joinpath(@__DIR__, "..", "local", "lib", "libsingularwrap.so")))

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
