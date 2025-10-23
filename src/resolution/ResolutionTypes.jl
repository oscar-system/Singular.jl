###############################################################################
#
#   ResolutionSet/sresolution
#
###############################################################################

const ResolutionSetID = Dict{PolyRing, Set}()

mutable struct ResolutionSet{T <: Nemo.RingElem} <: Set
   base_ring::PolyRing

   function ResolutionSet{T}(R::PolyRing) where T
      @assert isconcretetype(T)
      return get!(ResolutionSetID, R) do
         new(R)
      end
   end
end

mutable struct sresolution{T <: Nemo.RingElem} <: Nemo.SetElem
   base_ring::PolyRing
   ptr::libSingular.syStrategy_ptr
   minimal::Bool
   index_1_is_an_ideal::Bool

   # really takes a Singular module, which has type ideal
   # take ownership of the pointer - not for general users
   function sresolution{T}(R::PolyRing, ptr::libSingular.syStrategy_ptr,
                           minimal::Bool=false,
                           index_1_is_an_ideal::Bool=false) where T
      @assert isconcretetype(T)
      T === elem_type(R) || error("type mismatch")
      R.refcount += 1
      z = new(R, ptr, minimal, index_1_is_an_ideal)
      finalizer(_sresolution_clear_fn, z)
      return z
   end
end

function _sresolution_clear_fn(r::sresolution)
    R = base_ring(r)
    libSingular.res_Delete_helper(r.ptr, R.ptr)
    _PolyRing_clear_fn(R)
end
