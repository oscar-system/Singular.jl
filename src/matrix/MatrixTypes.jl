###############################################################################
#
#   MatrixSpace/smatrix 
#
###############################################################################

const MatrixSpaceID = Dict{Ring, Set}()

type MatrixSpace{T <: Nemo.RingElem} <: Set
   base_ring::PolyRing
   nrows::Int
   ncols::Int

   function MatrixSpace{T}(R::PolyRing, r::Int, c::Int) where T
      if haskey(MatrixSpaceID, (R, r, c))
         return MatrixSpaceID[R, r, c]
      else
         return MatrixSpaceID[R, r, c] = new(R, r, c)
      end
   end
end

type smatrix{T <: Nemo.RingElem} <: Nemo.SetElem
   ptr::libSingular.matrix
   base_ring::PolyRing

   # really takes a Singular module, which has type ideal
   function smatrix{T}(R::PolyRing, m::libSingular.ideal) where T
      ptr = libSingular.id_Copy(m, R.ptr)
      ptr = libSingular.id_Module2Matrix(ptr, R.ptr)
      z = new(ptr, R)
      finalizer(z, _smatrix_clear_fn)
      return z
   end

   function smatrix{T}(R::PolyRing, ptr::libSingular.matrix) where {T}
      z = new(ptr, R)
      finalizer(z, _smatrix_clear_fn)
      return z
   end
end

function _smatrix_clear_fn(I::smatrix)
   libSingular.mp_Delete(I.ptr, I.base_ring.ptr)
end

