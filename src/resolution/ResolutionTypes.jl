###############################################################################
#
#   ResolutionSet/sresolution 
#
###############################################################################

type ResolutionSet{T <: Nemo.RingElem} <: Nemo.Set
   base_ring::PolyRing

   function ResolutionSet{T}(R::PolyRing) where T
      return new(R)
   end
end

type sresolution{T <: Nemo.RingElem} <: Nemo.SetElem
   ptr::libSingular.resolvente
   len::Int
   base_ring::PolyRing

   # really takes a Singular module, which has type ideal
   function sresolution{T}(R::PolyRing, n::Int, ptr::libSingular.resolvente) where T
      z = new(ptr, n, R)
      finalizer(z, _sresolution_clear_fn)
      return z
   end
end

function _sresolution_clear_fn(r::sresolution)
   for i = 1:r.len
      libSingular.sy_Delete(r.ptr, Cint(i), r.base_ring.ptr)
   end
end