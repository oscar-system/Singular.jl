###############################################################################
#
#   ResolutionSet/sresolution 
#
###############################################################################

const ResolutionSetID = Dict{Ring, Set}()

mutable struct ResolutionSet{T <: Nemo.RingElem} <: Set
   base_ring::PolyRing

   function ResolutionSet{T}(R::PolyRing) where T
      if haskey(ResolutionSetID, R)
         return ResolutionSetID[R]
      else
         return ResolutionSetID[R] = new(R)
      end
   end
end

mutable struct sresolution{T <: Nemo.RingElem} <: Nemo.SetElem
   ptr::Ptr{Nothing}
   len::Int
   minimal::Bool
   base_ring::PolyRing

   # really takes a Singular module, which has type ideal
   function sresolution{T}(R::PolyRing, n::Int, ptr::Ptr{Nothing}, minimal::Bool=false) where T
      R.refcount += 1
      z = new(ptr, n, minimal, R)
      finalizer(_sresolution_clear_fn, z)
      return z
   end
end

function _sresolution_clear_fn(r::sresolution)
   R = base_ring(r)
   libSingular.res_Delete_helper(r.ptr, Cint(r.len), R.ptr)
    _PolyRing_clear_fn(R)
end
