###############################################################################
#
#   ResolutionSet/sresolution 
#
###############################################################################

const ResolutionSetID = Dict{Ring, Set}()

type ResolutionSet{T <: Nemo.RingElem} <: Set
   base_ring::PolyRing

   function ResolutionSet{T}(R::PolyRing) where T
      if haskey(ResolutionSetID, R)
         return ResolutionSetID[R]
      else
         return ResolutionSetID[R] = new(R)
      end
   end
end

type sresolution{T <: Nemo.RingElem} <: Nemo.SetElem
   ptr::libSingular.resolvente
   len::Int
   minimal::Bool
   base_ring::PolyRing

   # really takes a Singular module, which has type ideal
   function sresolution{T}(R::PolyRing, n::Int, ptr::libSingular.resolvente, minimal::Bool=false) where T
      z = new(ptr, n, minimal, R)
      finalizer(z, _sresolution_clear_fn)
      return z
   end
end

function _sresolution_clear_fn(r::sresolution)
   libSingular.res_Delete(r.ptr, Cint(r.len), r.base_ring.ptr)
end
