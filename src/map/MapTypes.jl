export Map

###############################################################################
#
#   GeneralType
#
###############################################################################

Map(::Type{T}) where T <: AbstractAlgebra.Map = AbstractAlgebra.Map(T)

###############################################################################
#
#   IdentityAlgebraHomomorphism
#
###############################################################################

mutable struct SIdAlgHom{T} <: AbstractAlgebra.Map{PolyRing, PolyRing,
        Singular.AbstractAlgebraHomomorphism, SIdAlgHom} where T <: Union{Ring, Field}
   domain::PolyRing
   image::sideal
   function SIdAlgHom{T}(R::PolyRing) where T <: Union{Ring, Field}
   z = new(R, Singular.MaximalIdeal(R, 1))
   end
end

###############################################################################
#
#   AlgebraHomomorphism
#
###############################################################################

mutable struct SAlgHom{T} <: AbstractAlgebra.Map{PolyRing, PolyRing,
          Singular.AbstractAlgebraHomomorphism, SAlgHom} where T <: Union{Ring, Field}

   domain::PolyRing
   codomain::PolyRing
   image::sideal

   function SAlgHom{T}(domain::PolyRing, codomain::PolyRing,
                                   m::sideal) where T <: Union{Ring, Field}
      z = new(domain, codomain, m)
      return z
   end
end
