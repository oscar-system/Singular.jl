module libSingular

using CxxWrap

using Nemo

import Base:setindex!, getindex

import Singular: init_wrap
    
function __init__()
   init_wrap()
end

# include("libsingular/LibSingularTypes.jl")

# include("libsingular/coeffs.jl")

# include("libsingular/rings.jl")

# include("libsingular/matrices.jl")

# include("libsingular/ideals.jl")

# include("libsingular/resolutions.jl")

end # module
