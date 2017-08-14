###############################################################################
#
#   SingularModuleClass/smodule 
#
###############################################################################

const SingularModuleClassID = ObjectIdDict()

type SingularModuleClass{T <: Nemo.RingElem} <: Nemo.Set
   base_ring::SingularPolyRing

   function SingularModuleClass{T}(R::SingularPolyRing) where T
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
   isGB::Bool

   function smodule{T}(R::SingularPolyRing, m::libSingular.ideal) where T
      z = new(m, R, false)
      R.refcount += 1
      finalizer(z, _smodule_clear_fn)
      return z
   end
end

function _smodule_clear_fn(I::smodule)
   R = I.base_ring
   libSingular.id_Delete(I.ptr, R.ptr)
   _SingularPolyRing_clear_fn(R)
end

###############################################################################
#
#   SingularPolyRing/spoly 
#
###############################################################################

const SingularFreeModID = ObjectIdDict()

type SingularFreeMod{T <: Nemo.RingElem} <: Nemo.Module{T}
   base_ring::SingularPolyRing
   rank::Int

   function SingularFreeMod{T}(R::SingularPolyRing, r::Int) where T
      if haskey(SingularFreeModID, (R, r))
         return SingularFreeModID[R, r]::SingularFreeMod{T}
      else
         return SingularFreeModID[R, r] = new(R, r)
      end
   end
end

type svector{T <: Nemo.RingElem} <: Nemo.ModuleElem{T}
   ptr::libSingular.poly # not really a polynomial
   rank::Int
   base_ring::SingularPolyRing

   function svector{T}(R::SingularPolyRing, r::Int, p::libSingular.poly) where T
      z = new(p, r, R)
      R.refcount += 1
      finalizer(z, _svector_clear_fn)
      return z
   end
end

function _svector_clear_fn(p::svector)
   R = p.base_ring
   libSingular.p_Delete(p.ptr, R.ptr)
   _SingularPolyRing_clear_fn(R)
end
