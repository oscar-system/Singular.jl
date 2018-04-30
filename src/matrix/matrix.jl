export nrows, ncols

###############################################################################
#
#   Basic manipulation
#
###############################################################################

nrows(M::smatrix) = Int(libSingular.nrows(M.ptr))

ncols(M::smatrix) = Int(libSingular.ncols(M.ptr))

function parent(M::smatrix)
   return MatrixSpace(M.base_ring, nrows(M), ncols(M))
end

base_ring(S::MatrixSpace) = S.base_ring

base_ring(M::smatrix) = M.base_ring

elem_type(S::MatrixSpace) = smatrix

parent_type(M::smatrix) = MatrixSpace

function getindex(M::smatrix, i::Int, j::Int)
   R = base_ring(M)
   ptr = libSingular.getindex(M.ptr, Cint(i), Cint(j))
   return R(libSingular.p_Copy(ptr, R.ptr))
end

function iszero(M::smatrix)
   for i = 1:nrows(M)
      for j = 1:ncols(M)
         if !iszero(M[i, j])
            return false
         end
      end
   end
   return true
end

###############################################################################
#
#   String I/O 
#
###############################################################################

function show(io::IO, S::MatrixSpace)
   print(io, "Space of ", S.nrows, "x", S.ncols, " Singular Matrices over ")
   show(io, base_ring(S))
end

function show(io::IO, M::smatrix)
   print(io, "[")
   m = nrows(M)
   n = ncols(M)
   for i = 1:m
      for j = 1:n
         show(io, M[i, j])
         if j != n
            print(io, ", ")
         elseif i != m
            println(io, "")
         end
      end
   end
   print(io, "]")
end

###############################################################################
#
#   Binary operations
#
###############################################################################

function +(M::smatrix{T}, N::smatrix{T}) where T <: AbstractAlgebra.RingElem
   R = base_ring(M)
   R != base_ring(N) && error("Matrices are not over the same ring")
   (nrows(M) != nrows(N)) && error("Incompatible dimensions")
   (ncols(M) != ncols(N)) && error("Incompatible dimensions")
   ptr = libSingular.mp_Add(M.ptr, N.ptr, R.ptr)
   return smatrix{T}(R, ptr)
end

function -(M::smatrix{T}, N::smatrix{T}) where T <: AbstractAlgebra.RingElem
   R = base_ring(M)
   R != base_ring(N) && error("Matrices are not over the same ring")
   (nrows(M) != nrows(N)) && error("Incompatible dimensions")
   (ncols(M) != ncols(N)) && error("Incompatible dimensions")
   ptr = libSingular.mp_Sub(M.ptr, N.ptr, R.ptr)
   return smatrix{T}(R, ptr)
end

function *(M::smatrix{T}, N::smatrix{T}) where T <: AbstractAlgebra.RingElem
   R = base_ring(M)
   R != base_ring(N) && error("Matrices are not over the same ring")
   (ncols(M) != nrows(N)) && error("Incompatible dimensions")
   ptr = libSingular.mp_Mult(M.ptr, N.ptr, R.ptr)
   return smatrix{T}(R, ptr)
end

###############################################################################
#
#   Matrix constructors
#
###############################################################################

function Matrix{T <: Nemo.RingElem}(I::smodule{T})
   return smatrix{T}(base_ring(I), I.ptr)
end

