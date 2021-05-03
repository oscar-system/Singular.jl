export identity_matrix, MatrixSpace, nrows, ncols, smatrix, zero_matrix

###############################################################################
#
#   Basic manipulation
#
###############################################################################

@doc Markdown.doc"""
    nrows(M::smatrix)

Return the number of rows of $M$.
"""
nrows(M::smatrix) = Int(libSingular.nrows(M.ptr))

@doc Markdown.doc"""
   nrows(M::smatrix)
Return the number of colums of $M$.
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
Given a matrix $M = (m_{ij})_{i, j}$, return the entry $m_{ij}$.
"""
function getindex(M::smatrix{T}, i::Int, j::Int) where T <: AbstractAlgebra.RingElem
   (i > nrows(M) || j > ncols(M)) && error("Incompatible dimensions")
   R = base_ring(M)
   GC.@preserve M R begin
      ptr = libSingular.getindex(M.ptr, Cint(i), Cint(j))
      return R(libSingular.p_Copy(ptr, R.ptr))
   end
end

@doc Markdown.doc"""
   setindex!(M::smatrix, p::spoly i::Int, j::Int)
Given a matrix $M = (m_{ij})_{i, j}$, set the entry $m_{ij}$ to the value $p$.
"""
function setindex!(M::smatrix, p::spoly, i::Int, j::Int)
   (i > nrows(M) || j > ncols(M)) && error("Incompatible dimensions")
   R = base_ring(M)
   R != parent(p) && error(" Base rings do not match.")
   GC.@preserve M p R libSingular.setindex(M.ptr, p.ptr, Cint(i), Cint(j), R.ptr)
end

"""
    iszero(M::smatrix)

Return whether the supplied matrix `M` is the zero matrix.
"""
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
Given a matrix $M=(m_{ij})_{i, j}$, return the matrix $M^T=(m_{ji})_{j, i}$.
"""
function transpose(M::smatrix{T}) where T <: AbstractAlgebra.RingElem
   R = base_ring(M)
   ptr = GC.@preserve M R libSingular.mp_Transp(M.ptr, R.ptr)
   return smatrix{T}(R, ptr)
end

function deepcopy_internal(M::smatrix, dict::IdDict)
   R = base_ring(M)
   ptr = GC.@preserve M R libSingular.mp_Copy(M.ptr, R.ptr)
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
   ptr = GC.@preserve M N R libSingular.mp_Add(M.ptr, N.ptr, R.ptr)
   return smatrix{T}(R, ptr)
end

function -(M::smatrix{T}, N::smatrix{T}) where T <: AbstractAlgebra.RingElem
   R = base_ring(M)
   R != base_ring(N) && error("Matrices are not over the same ring")
   (nrows(M) != nrows(N)) && error("Incompatible dimensions")
   (ncols(M) != ncols(N)) && error("Incompatible dimensions")
   ptr = GC.@preserve M N R libSingular.mp_Sub(M.ptr, N.ptr, R.ptr)
   return smatrix{T}(R, ptr)
end

function *(M::smatrix{T}, N::smatrix{T}) where T <: AbstractAlgebra.RingElem
   R = base_ring(M)
   R != base_ring(N) && error("Matrices are not over the same ring")
   (ncols(M) != nrows(N)) && error("Incompatible dimensions")
   ptr = GC.@preserve M N R libSingular.mp_Mult(M.ptr, N.ptr, R.ptr)
   return smatrix{T}(R, ptr)
end

function *(p::spoly{T}, M::smatrix{spoly{T}}) where T <: AbstractAlgebra.RingElem
   R = base_ring(M)
   R != parent(p) && error("Base rings do not match.")
   GC.@preserve M R begin
      x = libSingular.mp_Copy(M.ptr, R.ptr)
      y = libSingular.p_Copy(p.ptr, R.ptr)
      ptr = libSingular.mp_MultP(x, y, R.ptr)
      return smatrix{spoly{T}}(R, ptr)
   end
end

function *(M::smatrix{spoly{T}}, p::spoly{T}) where T <: AbstractAlgebra.RingElem
   return p*M
end

function *(i::Int, M::smatrix{T}) where T <: AbstractAlgebra.RingElem
   R = base_ring(M)
   return R(i) * M
end

function *(M::smatrix{T}, i::Int) where T <: AbstractAlgebra.RingElem
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
   GC.@preserve M N R return Bool(libSingular.mp_Equal(M.ptr, N.ptr, R.ptr))
end

###############################################################################
#
#   Parent object call
#
###############################################################################

function (S::MatrixSpace{T})(ptr::libSingular.matrix_ptr) where T <: AbstractAlgebra.RingElem
   M = smatrix{T}(base_ring(S), ptr)
   (S.ncols != ncols(M) || S.nrows != nrows(M)) && error("Incompatible dimensions")
   return M
end

(S::MatrixSpace)() = zero_matrix(S.base_ring, S.nrows, S.ncols)

(S::MatrixSpace)(a::smatrix) =
   S.base_ring == a.base_ring && S.nrows == nrows(a) && S.ncols == ncols(a) ?
      a :
      throw(ArgumentError("unable to coerce matrix"))

(S::MatrixSpace)(x::RingElement) = diagonal_matrix(S.base_ring(x), S.nrows, S.ncols)

###############################################################################
#
#   Matrix constructors
#
###############################################################################

function Matrix(I::smodule{T}) where T <: Nemo.RingElem
   GC.@preserve I return smatrix{T}(base_ring(I), I.ptr)
end

function Matrix(I::sideal{T}) where T <: Nemo.RingElem
   GC.@preserve I return smatrix{T}(base_ring(I), I.ptr)
end

function Module(vecs::smatrix{spoly{T}}) where T <: Nemo.RingElem
   R = base_ring(vecs)
   S = elem_type(R)
   GC.@preserve vecs return smodule{S}(R, vecs.ptr)
end

@doc Markdown.doc"""
   identity_matrix(R::PolyRing, n::Int)
Returns the $n \times n$ identity matrix over $R.$
"""
function identity_matrix(R::PolyRing, n::Int)
   p = R(1)
   GC.@preserve p R return smatrix{elem_type(R)}(R, libSingular.mp_InitP(n, p.ptr, R.ptr))
end

@doc Markdown.doc"""
   zero_matrix(R::PolyRing, r::Int, c::Int)
Returns the $r \times c$ zero matrix over $R.$
"""
function zero_matrix(R::PolyRing, r::Int, c::Int)
   return smatrix{elem_type(R)}(R, libSingular.mpNew(r, c))
end
