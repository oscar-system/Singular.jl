###############################################################################
#
#   AlgebraHomomorphism
#
###############################################################################

mutable struct salghom{T <: Union{Singular.Ring, Singular.Field}} <:
Singular.Map{Singular.PolyRing, Singular.PolyRing, Singular.Map, salghom}

   domain::PolyRing
   codomain::PolyRing
   image::sideal

   function salghom{T}(D::PolyRing, C::PolyRing, m::sideal) where T <:
   Union{Singular.Ring, Singular.Field}
      z = new(D, C, m)
      return z
   end
end

