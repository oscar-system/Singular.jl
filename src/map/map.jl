export get_field, set_field!

################################################################################
#
#  All maps
#
################################################################################

get_field(M::AbstractAlgebra.Map, f) = getfield(M, f) # fall back to Julia builtin
set_field!(M::AbstractAlgebra.Map, f) = setfield!(M, f) # fall back to Julia builtin
set_field!(M::AbstractAlgebra.Map, f, g) = setfield!(M, f, g)

@doc Markdown.doc"""
   check_composable(a::AbstractAlgebra.Map, b::AbstractAlgebra.Map)
> Returns true if the maps $a$ and $b$ are composable, otherwise false.
"""
function check_composable(a::AbstractAlgebra.Map, b::AbstractAlgebra.Map)
   codomain(a) != domain(b) && error("Incompatible maps")
end

