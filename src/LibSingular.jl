module libSingular

using CxxWrap

wrap_module(Pkg.dir("Singular", "local", "lib", "libsingularwrap.so"), libSingular)

include("libsingular/LibSingularTypes.jl")

include("libsingular/coeffs.jl")

include("libsingular/rings.jl")

include("libsingular/matrices.jl")

include("libsingular/ideals.jl")

include("libsingular/resolutions.jl")

end # module
