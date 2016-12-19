###############################################################################
#
#   SingularIdealSet/sideal 
#
###############################################################################

const SingularIdealSetID = ObjectIdDict()

type SingularIdealSet{T <: Nemo.RingElem} <: Nemo.Set
   base_ring::SingularPolyRing

   function SingularIdealSet(R::SingularPolyRing)
      if haskey(SingularIdealSetID, R)
         return SingularIdealSetID[R]
      else
         return SingularIdealSetID[R] = new(R)
      end
   end
end

type sideal{T <: Nemo.RingElem} <: Nemo.Module{T}
   ptr::libSingular.ideal
   base_ring::SingularPolyRing
   isGB::Bool

   function sideal(R::SingularPolyRing, ids::spoly...)
      n = length(ids)
      id = libSingular.idInit(Cint(n))
      z = new(id, R, false)
      finalizer(z, _sideal_clear_fn)
      for i = 1:n
         p = libSingular.p_Copy(ids[i].ptr, R.ptr)
         libSingular.setindex!(id, p, Cint(i - 1))
      end
      return z
   end

   function sideal(R::SingularPolyRing, id::libSingular.ideal)
      z = new(id, R, false)
      finalizer(z, _sideal_clear_fn)
      return z
   end
end

function _sideal_clear_fn(I::sideal)
   libSingular.id_Delete(I.ptr, I.base_ring.ptr)
end
