###############################################################################
#
#   IdealSet/sideal
#
###############################################################################

const IdealSetID = Dict{PolyRing, Set}()

mutable struct IdealSet{T <: Nemo.RingElem} <: Set
   base_ring::PolyRing

   function IdealSet{T}(R::PolyRing) where T
      return get!(IdealSetID, R) do
         new(R)
      end
   end
end

mutable struct sideal{T <: Nemo.RingElem} <: Module{T}
   base_ring::PolyRing
   ptr::libSingular.ideal_ptr
   isGB::Bool

   function sideal{T}(R::PolyRing, id::libSingular.ideal_ptr) where T
      T === elem_type(R) || error("type mismatch")
      z = new(R, id, false)
      R.refcount += 1
      finalizer(_sideal_clear_fn, z)
      return z
   end
end

function sideal{T}(R::PolyRing, ids::spoly...) where T
   n = length(ids)
   id = libSingular.idInit(Cint(n),1)
   z = sideal{T}(R, id)
   for i = 1:n
      p = libSingular.p_Copy(ids[i].ptr, R.ptr)
      libSingular.setindex_internal(id, p, Cint(i - 1))
   end
   return z
end

function _sideal_clear_fn(I::sideal)
   R = I.base_ring
   libSingular.id_Delete(I.ptr, R.ptr)
   _PolyRing_clear_fn(R)
end
