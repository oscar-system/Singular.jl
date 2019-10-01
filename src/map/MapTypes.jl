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
   image::Vector
   ptr::libSingular.ideal

   function SIdAlgHom{T}(R::PolyRing) where T <: Union{Ring, Field}
      M = Singular.MaximalIdeal(R, 1)
      z = new(R, gens(M), M.ptr)
      return z
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
   image::Vector
   ptr::libSingular.ideal

   function SAlgHom{T}(domain::PolyRing, codomain::PolyRing,
             V::Vector) where T <: Union{Ring, Field}

      I = Ideal(codomain, V)
      z = new(domain, codomain, V, I.ptr)
      return z
   end

   function SAlgHom{T}(domain::PolyRing, codomain::PolyRing,
             V::Vector, ptr::libSingular.ideal) where T <: Union{Ring, Field}

      z = new(domain, codomain, V, ptr)
      return z
   end
end
