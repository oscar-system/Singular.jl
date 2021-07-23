###############################################################################
#
#   MatrixSpace/smatrix
#
###############################################################################

mutable struct MatrixSpace{T <: Nemo.RingElem} <: AbstractAlgebra.MatSpace{T}
   base_ring::PolyRing
   nrows::Int
   ncols::Int

   function MatrixSpace{T}(R::PolyRing, r::Int, c::Int) where T
      return get!(MatrixSpaceID, (R, r, c)) do
         new{T}(R, r, c)
      end
   end
end

const MatrixSpaceID = Dict{Tuple{PolyRing, Int, Int}, MatrixSpace}()

mutable struct smatrix{T <: Nemo.RingElem} <: AbstractAlgebra.MatElem{T}
   ptr::libSingular.matrix_ptr
   base_ring::PolyRing

   function smatrix{T}(R::PolyRing, ptr::libSingular.matrix_ptr) where {T}
      T === elem_type(R) || error("type mismatch")
      z = new(ptr, R)
      finalizer(_smatrix_clear_fn, z)
      return z
   end
end

# really takes a Singular module, which has type ideal
function smatrix{T}(R::PolyRing, m::libSingular.ideal_ptr) where T
   ptr = libSingular.id_Copy(m, R.ptr)
   ptr = libSingular.id_Module2Matrix(ptr, R.ptr)
   return smatrix{T}(R, ptr)
end

function _smatrix_clear_fn(I::smatrix)
   libSingular.mp_Delete(I.ptr, I.base_ring.ptr)
end
