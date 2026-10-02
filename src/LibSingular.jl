module libSingular

using CxxWrap

import ..Singular: libflint, libsingular_julia, AbstractAlgebra, Setup

@wrapmodule(libsingular_julia)

function _check_omalloc_page_size(wrapper_page_size::Integer, singular_page_size::Integer)
   wrapper_page_size == singular_page_size && return nothing
   error(
      "Incompatible omalloc page sizes: Singular.jl requires $(wrapper_page_size) " *
      "bytes, but the selected Singular_jll was built for $(singular_page_size) bytes. " *
      "Rebuild libsingular_julia against the selected Singular_jll."
   )
end

_check_omalloc_page_size(config_path::AbstractString) =
   _check_omalloc_page_size(omalloc_page_size(), Setup.omalloc_page_size(config_path))

_check_omalloc_page_size() =
   _check_omalloc_page_size(omalloc_page_size(), Setup.omalloc_page_size())

function __init__()
   @initcxx
   _check_omalloc_page_size()
   initialize_jl_c_types(@__MODULE__)
end

include("libsingular/LibSingularTypes.jl")

include("libsingular/errors.jl")

include("libsingular/coeffs.jl")

include("libsingular/rings.jl")

include("libsingular/matrices.jl")

include("libsingular/ideals.jl")

include("libsingular/resolutions.jl")

end # module
