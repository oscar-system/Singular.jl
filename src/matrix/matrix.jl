export identity_matrix, MatrixSpace, nrows, ncols, smatrix, zero_matrix

###############################################################################
#
#   Basic manipulation
#
###############################################################################

@doc Markdown.doc"""
   nrows(M::smatrix)
> Return the number of rows of $M$.
"""
nrows(M::smatrix) = Int(libSingular.nrows(M.ptr))

@doc Markdown.doc"""
   nrows(M::smatrix)
> Return the number of colums of $M$.
"""
ncols(M::smatrix) = Int(libSingular.ncols(M.ptr))

function parent(M::smatrix{T}) where T <: AbstractAlgebra.RingElem
   return MatrixSpace{T}(M.base_ring, nrows(M), ncols(M))
end

base_ring(S::MatrixSpace) = S.base_ring

base_ring(M::smatrix) = M.base_ring

elem_type(::Type{MatrixSpace{T}}) where T <: AbstractAlgebra.RingElem = smatrix{T}

elem_type(::MatrixSpace{T}) where T <: AbstractAlgebra.RingElem = smatrix{T}

parent_type(::Type{smatrix{T}}) where T <: AbstractAlgebra.RingElem = MatrixSpace{T}

@doc Markdown.doc"""
   getindex(M::smatrix{T}, i::Int, j::Int) where T <: AbstractAlgebra.RingElem
> Given a matrix $M=(m_{ij})_{i, j}$, return the entry $m_{ij}$.
"""
function getindex(M::smatrix{T}, i::Int, j::Int) where T <: AbstractAlgebra.RingElem
   (i > nrows(M) || j > ncols(M)) && error("Incompatible dimensions")
   R = base_ring(M)
   ptr = libSingular.getindex(M.ptr, Cint(i), Cint(j))
   return R(libSingular.p_Copy(ptr, R.ptr))
end

@doc Markdown.doc"""
   setindex!(M::smatrix, p::spoly i::Int, j::Int)
> Given a matrix $M=(m_{ij})_{i, j}$, set the entry $m_{ij}$ to the value $p$.
"""
function setindex!(M::smatrix, p::spoly, i::Int, j::Int)
   (i > nrows(M) || j > ncols(M)) && error("Incompatible dimensions")
   R = base_ring(M)
   R != parent(p) && error(" Base rings do not match.")
   libSingular.setindex(M.ptr, p.ptr, Cint(i), Cint(j), R.ptr)
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

@doc Markdown.doc"""
   transpose(M::smatrix{T}) where T <: AbstractAlgebra.RingElem
> Given a matrix $M=(m_{ij})_{i, j}$, return the matrix $M^T=(m_{ji})_{j, i}$.
"""
function transpose(M::smatrix{T}) where T <: AbstractAlgebra.RingElem
   R = base_ring(M)
   ptr = libSingular.mp_Transp(M.ptr, R.ptr)
   return smatrix{T}(R, ptr)
end

function deepcopy_internal(M::smatrix, dict::IdDict)
   R = base_ring(M)
   ptr = libSingular.mp_Copy(M.ptr, R.ptr)
   return parent(M)(ptr)
end

function hash(M::smatrix, h::UInt)
   v = 0x68bf046fc9d6afbf%UInt
   for i in 1:nrows(M)
      for j in 1:ncols(M)
         v = xor(hash(M[i, j], h), v)
      end
   end
   return v
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

function *(p::spoly, M::smatrix{T}) where T <: AbstractAlgebra.RingElem
   R = base_ring(M)
   R != parent(p) && error("Base rings do not match.")
   x = libSingular.mp_Copy(M.ptr, R.ptr)
   y = libSingular.p_Copy(p.ptr, R.ptr)
   ptr = libSingular.mp_MultP(x, y, R.ptr)
   return smatrix{T}(R, ptr)
end

function *(i::Int, M::smatrix{T}) where T <: AbstractAlgebra.RingElem
   R = base_ring(M)
   return R(i) * M
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(M::smatrix{T}, N::smatrix{T}) where T <: AbstractAlgebra.RingElem
   R = base_ring(M)
   R != base_ring(N) && error("Matrices are not over the same ring")
   (nrows(M) != nrows(N)) && error("Incompatible dimensions")
   (ncols(M) != ncols(N)) && error("Incompatible dimensions")
   return Bool(libSingular.mp_Equal(M.ptr, N.ptr, R.ptr))
end

###############################################################################
#
#   Concatenation
#
###############################################################################

@doc Markdown.doc"""
   hcat(A::smatrix{T}, B::smatrix{T}) where T <: AbstractAlgebra.RingElem
> Return the horizontal concatenation of $A$ and $B$.
> Assumes that the number of rows is the same in $A$ and $B$.
"""
function hcat(A::smatrix{T}, B::smatrix{T}) where T <: AbstractAlgebra.RingElem
   nr = nrows(A)
   R = base_ring(A)
   (nrows(B) != nr) && error("Matrices must have same number of rows.")
   (R != base_ring(B)) && error("Matrices are not over the same ring")

   nca = ncols(A)
   ncb = ncols(B)

   Z = zero_matrix(R, nr, nca + ncb)
   for i in 1:nr
      for j in 1:nca
         Z[i, j] = A[i, j]
      end

      for j in 1:ncb
         Z[i, nca + j] = B[i, j]
      end
   end
   return Z
end

@doc Markdown.doc"""
   hcat(A::Vector{ <: smatrix{T}) where T <: AbstractAlgebra.RingElem
> Return the horizontal concatenation of the matrices $A$.
> All component matrices need to have the same base ring and number of rows.
"""
function hcat(A::Vector{ <: smatrix{T} }) where T <: AbstractAlgebra.RingElem
  R = base_ring(A[1])
  nr = nrows(A[1])
  !all( M -> base_ring(M) == R, A) && error("Matrices are not over the same ring")
  !all( M -> nrows(M) == nr, A) && error("Matrices must have same number of rows.")

  l = length(A)
  nc = sum(ncols.(A))
  m = 0
  Z = zero_matrix(R, nr, nc)
  for i in 1:l
     nca = ncols(A[i])
     for j in m + 1:m + nca
        for k in 1:nr
           Z[k, j] = A[i][k, j - m]
        end
     end
     m += nca
  end
  return Z
end

@doc Markdown.doc"""
   vcat(A::smatrix{T}, B::smatrix{T}) where T <: AbstractAlgebra.RingElem
> Return the vertical concatenation of $A$ and $B$.
> Assumes that the number of columns is the same in $A$ and $B$.
"""
function vcat(A::smatrix{T}, B::smatrix{T}) where T <: AbstractAlgebra.RingElem
   nc = ncols(A)
   R = base_ring(A)
   (ncols(B) != nc) && error("Matrices must have same number of columns.")
   (R != base_ring(B)) && error("Matrices are not over the same ring")

   nra = nrows(A)
   nrb = nrows(B)

   Z = zero_matrix(R, nra + nrb, nc)
   for i in 1:nc
      for j in 1:nra
         Z[j ,i] = A[j, i]
      end

      for j in 1:nrb
         Z[nra + j, i] = B[j, i]
      end
   end
   return Z
end

@doc Markdown.doc"""
   vcat(A::Vector{ <: smatrix{T} })
> Return the vertical concatenation of the matrices $A$.
> All component matrices need to have the same base ring and number of columns.
"""
function vcat(A::Vector{ <: smatrix{T} }) where T <: AbstractAlgebra.RingElem
  R = base_ring(A[1])
  nc = ncols(A[1])
  !all( M -> base_ring(M) == R, A) && error("Matrices are not over the same ring")
  !all( M -> ncols(M) == nc, A) && error("Matrices must have same number of columns.")

  l = length(A)
  nr = sum(nrows.(A))
  m = 0
  Z = zero_matrix(R, nr, nc)
  for i in 1:l
     nra = nrows(A[i])
     for j in m + 1:m + nra
        for k in 1:nc
           Z[j, k] = A[i][j - m, k]
        end
     end
     m += nra
  end
  return Z
end

###############################################################################
#
#   Parent object call
#
###############################################################################

function (S::MatrixSpace{T})(ptr::libSingular.matrix_ref) where T <: AbstractAlgebra.RingElem
   M = smatrix{T}(base_ring(S), ptr)
   (S.ncols != ncols(M) || S.nrows != nrows(M)) && error("Incompatible dimensions")
   return M
end

###############################################################################
#
#   Matrix constructors
#
###############################################################################

function Matrix(I::smodule{T}) where T <: Nemo.RingElem
   return smatrix{T}(base_ring(I), I.ptr)
end

function Matrix(I::sideal{T}) where T <: Nemo.RingElem
   return smatrix{T}(base_ring(I), I.ptr)
end

@doc Markdown.doc"""
   identity_matrix(R::PolyRing, n::Int)
> Returns the $n \times n$ identity matrix over $R.$
"""
function identity_matrix(R::PolyRing, n::Int)
   p = R(1)
   return smatrix{elem_type(R)}(R, libSingular.mp_InitP(n, p.ptr, R.ptr))
end

@doc Markdown.doc"""
   zero_matrix(R::PolyRing, r::Int, c::Int)
> Returns the $r \times c$ zero matrix over $R.$
"""
function zero_matrix(R::PolyRing, r::Int, c::Int)
   return smatrix{elem_type(R)}(R, libSingular.mpNew(R.ptr, r, c))
end
