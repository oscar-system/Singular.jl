module libSingular

import Libdl
using CxxWrap

import ..Singular: libflint, libantic

const libsingularwrap_path = joinpath(@__DIR__, "..", "deps", "usr",
                "lib", "libsingularwrap." * Libdl.dlext)
if !isfile(libsingularwrap_path)
    error("""Singular.jl needs to be compiled; please run `using Pkg; Pkg.build("Singular")`""")
end
@wrapmodule(realpath(libsingularwrap_path))

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
