###############################################################################
#
#   SingularMatrixSpace/smatrix 
#
###############################################################################

const SingularMatrixSpaceID = ObjectIdDict()

type SingularMatrixSpace{T <: Nemo.RingElem} <: Nemo.Set
   base_ring::SingularPolyRing
   nrows::Int
   ncols::Int

   function SingularMatrixSpace(R::SingularPolyRing, r::Int, c::Int)
      if haskey(SingularMatrixSpaceID, (R, r, c))
         return SingularMatrixSpaceID[R, r, c]
      else
         return SingularMatrixSpaceID[R, r, c] = new(R, r, c)
      end
   end
end

type smatrix <: Nemo.SetElem
   ptr::libSingular.matrix
   base_ring::SingularPolyRing

   # really takes a Singular module, which has type ideal
   function smatrix(R::SingularPolyRing, m::libSingular.ideal)
      ptr = libSingular.id_Copy(m, R.ptr)
      ptr = libSingular.id_Module2Matrix(ptr, R.ptr)
      z = new(ptr, R)
      finalizer(z, _smatrix_clear_fn)
      return z
   end
end

function _smatrix_clear_fn(I::smatrix)
   libSingular.mp_Delete(I.ptr, I.base_ring.ptr)
end

