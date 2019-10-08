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
      V = gens(R)
      n = nvars(R)
      ptr = libSingular.idInit(Cint(n), 1)
      z = new(R, V, ptr)
      R.refcount += 1
      finalizer(_SIdAlgHom_clear_fn, z)
      for i = 1:n
         p = libSingular.p_Copy(V[i].ptr, R.ptr)
         libSingular.setindex_internal(ptr, p, Cint(i - 1))
      end
      return z
   end
end

function _SIdAlgHom_clear_fn(f::SIdAlgHom)
   R = f.domain
   libSingular.id_Delete(f.ptr, R.ptr)
   _PolyRing_clear_fn(R)
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

      n = length(V)
      ptr = libSingular.idInit(Cint(n), 1)
      z = new(domain, codomain, V, ptr)
      codomain.refcount += 1
      for i = 1:n
         p = libSingular.p_Copy(V[i].ptr, codomain.ptr)
         libSingular.setindex_internal(ptr, p, Cint(i - 1))
      end
      return z
   end
end

function _SAlgHom_clear_fn(f::SAlgHom)
   R = f.codomain
   libSingular.id_Delete(f.ptr, R.ptr)
   _PolyRing_clear_fn(R)
end

