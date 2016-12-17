###############################################################################
#
#   SingularModuleClass/smodule 
#
###############################################################################

const SingularModuleClassID = ObjectIdDict()

type SingularModuleClass{T <: Nemo.RingElem} <: Nemo.Set
   base_ring::SingularPolyRing

   function SingularModuleClass(R::SingularPolyRing)
      if haskey(SingularModuleClassID, R)
         return SingularModuleClassID[R]
      else
         return SingularModuleClassID[R] = new(R)
      end
   end
end

type smodule{T <: Nemo.RingElem} <: Nemo.Module{T}
   ptr::libSingular.ideal # ideal and module types are the same in Singular
   base_ring::SingularPolyRing

   function smodule(R::SingularPolyRing, m::libSingular.ideal)
      z = new(m, R)
      finalizer(z, _smodule_clear_fn)
      return z
   end
end

function _smodule_clear_fn(I::smodule)
   libSingular.id_Delete(I.ptr, I.base_ring.ptr)
end
