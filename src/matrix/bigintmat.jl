# The Singular bigintmat is in an unfortunate place. The entries belong
# to coeffs_BIGINT and not to libSingular's n_Z. Therefore, it does not make
# much sense to have a constructor
#   matrix(n_Z, [1 2; 3 4])
# return a matrix over the wrong integers. Nor does it make sense to make a new
# ring of integers for use in Singular.jl. Therefore, we keep the api very
# simple and allow conversions between fmpz_mat and Matrix{BigInt}.

function sbigintmat(r::Int, c::Int)
   return sbigintmat(libSingular.bigintmat_init(r, c))
end

function ncols(m::sbigintmat)
   return Int(libSingular.bigintmat_ncols(m.ptr))
end

function nrows(m::sbigintmat)
   return Int(libSingular.bigintmat_nrows(m.ptr))
end

function size(m::sbigintmat)
   return (ncols(m), nrows(m))
end

function getindex(m::sbigintmat, i::Int, j::Int)
   (0 < i <= nrows(m) && 0 < j <= ncols(m)) || error("index out of range")
   GC.@preserve m begin
      # not our own ptr
      ptr = libSingular.bigintmat_viewindex(m.ptr, Cint(i), Cint(j))
      return libSingular.n_GetMPZ(ptr, libSingular.get_coeffs_BIGINT())
   end
end

function setindex!(m::sbigintmat, a::Union{Integer, fmpz}, i::Int, j::Int)
   (0 < i <= nrows(m) && 0 < j <= ncols(m)) || error("index out of range")
   # rawset consumes our ptr
   ptr = libSingular.number_ptr(a, libSingular.get_coeffs_BIGINT())
   GC.@preserve m libSingular.bigintmat_rawset(m.ptr, ptr, Cint(i), Cint(j))
end

function show(io::IO, M::sbigintmat)
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

# conversion to bigintmat

function sbigintmat(a::Nemo.MatElem{ <: Union{Nemo.Integer, Nemo.fmpz}})
   (r, c) = (nrows(a), ncols(a))
   z = sbigintmat(r, c)
   for i in 1:r, j in 1:c
      z[i, j] = a[i, j]
   end
   return z
end

function sbigintmat(a::Matrix{ <: Integer})
   (r, c) = size(a)
   z = sbigintmat(r, c)
   for i in 1:r, j in 1:c
      z[i, j] = a[i, j]
   end
   return z
end


# conversion out of bigintmat

function Base.Array(a::sbigintmat)
   (r, c) = size(a)
   z = zeros(BigInt, r, c)
   for i in 1:r, j in 1:c
      z[i, j] = a[i, j]
   end
   return z
end

function Nemo.matrix(R::Nemo.Ring, a::sbigintmat)
   (r, c) = size(a)
   z = zero_matrix(R, r, c)
   for i in 1:r, j in 1:c
      z[i, j] = a[i, j]
   end
   return z
end

