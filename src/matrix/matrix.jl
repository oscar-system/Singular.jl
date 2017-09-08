###############################################################################
#
#   Basic manipulation
#
###############################################################################

nrows(M::smatrix) = Int(libSingular.nrows(M.ptr))

ncols(M::smatrix) = Int(libSingular.nrows(M.ptr))

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
#   Matrix constructors
#
###############################################################################

function Matrix{T <: Nemo.RingElem}(I::smodule{T})
   return smatrix{T}(base_ring(I), I.ptr)
end
