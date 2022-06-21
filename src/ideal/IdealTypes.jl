###############################################################################
#
#   IdealSet/sideal
#
###############################################################################

const IdealSetID = Dict{PolyRing, Set}()

mutable struct IdealSet{T <: AbstractAlgebra.NCRingElem} <: Set
   base_ring::PolyRingUnion

   function IdealSet{T}(R::PolyRing) where T
      return get!(IdealSetID, R) do
         new(R)
      end
   end
end

mutable struct sideal{T <: AbstractAlgebra.NCRingElem} <: Set
   ptr::libSingular.ideal_ptr
   base_ring::PolyRingUnion
   isGB::Bool
   # Singular supports left ideals and two-sided ideals.
   #    spoly: isTwoSided doesn't matter because is it commutative
   # spluralg: isTwoSided = false is supported and means left ideal
   #   slpalg: isTwoSided = false is currently NOT supported, that is all
   #           letterplace ideals must be two-sided at this point in time.
   isTwoSided::Bool

   # take ownership of the pointer - not for general users
   function sideal{T}(R::PolyRingUnion, id::libSingular.ideal_ptr, isGB::Bool, isTwoSided::Bool) where T
      z = new{T}(id, R, isGB, isTwoSided)
      R.refcount += 1
      finalizer(_sideal_clear_fn, z)
      return z
   end
end

function _sideal_clear_fn(I::sideal{T}) where T <: SPolyUnion
   R = I.base_ring
   libSingular.id_Delete(I.ptr, R.ptr)
   _PolyRing_clear_fn(R)
end

is_default_twosided_ideal(R::PolyRing) = true
is_default_twosided_ideal(R::LPRing) = true
is_default_twosided_ideal(R::PluralRing) = false

is_default_twosided_ideal(::Type{<:spoly}) = true
is_default_twosided_ideal(::Type{<:slpalg}) = true
is_default_twosided_ideal(::Type{<:spluralg}) = false

# take ownership of the pointer - not for general users
function sideal{S}(R::PolyRingUnion, id::libSingular.ideal_ptr, isGB::Bool) where S
   return sideal{S}(R, id, isGB, is_default_twosided_ideal(S))
end

# take ownership of the pointer - not for general users
function sideal{S}(R::PolyRingUnion, id::libSingular.ideal_ptr) where S
   return sideal{S}(R, id, false, is_default_twosided_ideal(S))
end

function sideal{S}(R::PolyRingUnion, ids::Vector{<:SPolyUnion}, isTwoSided::Bool) where S
   n = length(ids)
   id = libSingular.idInit(Cint(n), 1)
   for i = 1:n
      p = libSingular.p_Copy(ids[i].ptr, R.ptr)
      libSingular.setindex_internal(id, p, Cint(i - 1))
   end
   # n < 2 && !is_quotient_ring => isGB
   return sideal{S}(R, id, (n < 2) && !is_quotient_ring(R), isTwoSided)
end

function sideal{T}(R::PolyRingUnion, ids::SPolyUnion...) where T
   n = length(ids)
   id = libSingular.idInit(Cint(n), 1)
   for i = 1:n
      p = libSingular.p_Copy(ids[i].ptr, R.ptr)
      libSingular.setindex_internal(id, p, Cint(i - 1))
   end
   # n < 2 && !is_quotient_ring => isGB
   return sideal{T}(R, id, (n < 2) && !is_quotient_ring(R), is_default_twosided_ideal(R))
end

