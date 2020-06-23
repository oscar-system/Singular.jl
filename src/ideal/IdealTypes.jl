###############################################################################
#
#   IdealSet/sideal
#
###############################################################################

const IdealSetID = Dict{PolyRing, Set}()

mutable struct IdealSet{T <: AbstractAlgebra.NCRingElem} <: Set
   base_ring::Union{PolyRing, WeylAlgebra, ExteriorAlgebra}

   function IdealSet{T}(R::PolyRing) where T
      if haskey(IdealSetID, R)
         return IdealSetID[R]
      else
         return IdealSetID[R] = new(R)
      end
   end
end

mutable struct sideal{T <: AbstractAlgebra.NCRingElem} <: Set
   ptr::libSingular.ideal_ptr
   base_ring::Union{PolyRing, WeylAlgebra, ExteriorAlgebra}
   isGB::Bool
   isTwoSided::Bool

   function sideal{T}(R::PolyRing, ids::spoly...) where T
      n = length(ids)
      id = libSingular.idInit(Cint(n),1)
      z = new(id, R, false, true)
      R.refcount += 1
      finalizer(_sideal_clear_fn, z)
      for i = 1:n
         p = libSingular.p_Copy(ids[i].ptr, R.ptr)
         libSingular.setindex_internal(id, p, Cint(i - 1))
      end
      return z
   end

   function sideal{T}(R::Union{WeylAlgebra, ExteriorAlgebra}, ids::Union{pweyl, pexterior}...) where T
      n = length(ids)
      id = libSingular.idInit(Cint(n),1)
      z = new(id, R, false, false)
      R.refcount += 1
      finalizer(_sideal_clear_fn, z)
      for i = 1:n
         p = libSingular.p_Copy(ids[i].ptr, R.ptr)
         libSingular.setindex_internal(id, p, Cint(i - 1))
      end
      return z
   end

   function sideal{T}(R::PolyRing, id::libSingular.ideal_ptr) where T
      z = new(id, R, false, true)
      R.refcount += 1
      finalizer(_sideal_clear_fn, z)
      return z
   end

   function sideal{T}(R::Union{WeylAlgebra, ExteriorAlgebra}, id::libSingular.ideal_ptr) where T
      z = new(id, R, false, false)
      R.refcount += 1
      finalizer(_sideal_clear_fn, z)
      return z
   end
end

function _sideal_clear_fn(I::sideal{spoly{T}}) where T
   R = I.base_ring
   libSingular.id_Delete(I.ptr, R.ptr)
   _PolyRing_clear_fn(R)
end

function _sideal_clear_fn(I::sideal{pweyl{T}}) where T
   R = I.base_ring
   libSingular.id_Delete(I.ptr, R.ptr)
   _WeylAlgebra_clear_fn(R)
end

function _sideal_clear_fn(I::sideal{pexterior{T}}) where T
   R = I.base_ring
   libSingular.id_Delete(I.ptr, R.ptr)
   _ExteriorAlgebra_clear_fn(R)
end

