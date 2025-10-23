###############################################################################
#
#   matrix_space/smatrix
#
###############################################################################

const MatrixSpaceID = Dict{Tuple{PolyRing, Int, Int}, Set}()

mutable struct matrix_space{T <: Nemo.RingElem} <: Set
   base_ring::PolyRing
   nrows::Int
   ncols::Int

   function matrix_space{T}(R::PolyRing, r::Int, c::Int) where T
      @assert isconcretetype(T)
      return get!(MatrixSpaceID, (R, r, c)) do
         new{T}(R, r, c)
      end::matrix_space{T}
   end
end

mutable struct smatrix{T <: Nemo.RingElem} <: Nemo.SetElem
   ptr::libSingular.matrix_ptr
   base_ring::PolyRing

   # take ownership of the pointer - not for general users
   function smatrix{T}(R::PolyRing, ptr::libSingular.matrix_ptr) where {T}
      @assert isconcretetype(T)
      T === elem_type(R) || error("type mismatch")
      z = new(ptr, R)
      finalizer(_smatrix_clear_fn, z)
      return z
   end
end

# really takes a Singular module, which has type ideal
# ownership of the pointer is NOT taken - not for general users
function smatrix{T}(R::PolyRing, m::libSingular.ideal_ptr) where T
   ptr = libSingular.id_Copy(m, R.ptr)
   ptr = libSingular.id_Module2Matrix(ptr, R.ptr)
   return smatrix{T}(R, ptr)
end

function _smatrix_clear_fn(I::smatrix)
   libSingular.mp_Delete(I.ptr, I.base_ring.ptr)
end


###############################################################################
#
#   sbigintmat: just a plain array of bigint's (coeffs_BIGINT)
#
###############################################################################

mutable struct sbigintmat
    ptr::libSingular.bigintmat_ptr

   # take ownership of the pointer - not for general users
   function sbigintmat(ptr::libSingular.bigintmat_ptr)
      z = new(ptr)
      finalizer(_sbigintmat_clear_fn, z)
      return z
   end
end

function _sbigintmat_clear_fn(m::sbigintmat)
   libSingular.bigintmat_clear(m.ptr)
end
