export FreeMod, svector, gens, rank, vector

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent{T <: Nemo.RingElem}(v::svector{T}) = FreeMod{T}(v.base_ring, v.rank)

base_ring(R::FreeMod) = R.base_ring

base_ring(v::svector) = v.base_ring

elem_type{T <: Nemo.RingElem}(::FreeMod{T}) = svector{T}

elem_type{T <: Nemo.RingElem}(::Type{FreeMod{T}}) = svector{T}

parent_type{T <: Nemo.RingElem}(::Type{svector{T}}) = FreeMod{T}

doc"""
    rank(M::FreeMod)
> Return the rank of the given free module.
"""
rank(M::FreeMod) = M.rank

doc"""
    gens{T <: AbstractAlgebra.RingElem}(M::FreeMod{T})
> Return a Julia array whose entries are the generators of the given free module.
"""
function gens(M::FreeMod{T}) where T <: AbstractAlgebra.RingElem
   R = base_ring(M)
   ptr = libSingular.id_FreeModule(Cint(M.rank), R.ptr)
   return [svector{T}(R, M.rank, libSingular.getindex(ptr, Cint(i - 1))) for i in 1:M.rank]
end

function deepcopy_internal(p::svector{T}, dict::ObjectIdDict) where T <: AbstractAlgebra.RingElem
   p2 = libSingular.p_Copy(p.ptr, base_ring(p).ptr)
   return svector{T}(p.base_ring, p.rank, p2)
end

function check_parent{T <: Nemo.RingElem}(a::svector{T}, b::svector{T})
   base_ring(a) != base_ring(b) && error("Incompatible base rings")
   a.rank != b.rank && error("Vectors of incompatible rank")
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::FreeMod)
   print(io, "Free Module of rank ", R.rank, " over ")
   show(io, R.base_ring)
end

function show(io::IO, a::svector)
   m = libSingular.p_String(a.ptr, base_ring(a).ptr)
   s = unsafe_string(m)
   libSingular.omFree(Ptr{Void}(m))
   print(io, s)
end


###############################################################################
#
#   Unary functions
#
###############################################################################

function -(a::svector{T}) where T <: AbstractAlgebra.RingElem
   R = base_ring(a)
   a1 = libSingular.p_Copy(a.ptr, R.ptr)
   s = libSingular.p_Neg(a1, R.ptr)
   return svector{T}(R, a.rank, s) 
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(a::svector{T}, b::svector{T}) where T <: AbstractAlgebra.RingElem
   check_parent(a, b)
   R = base_ring(a)
   a1 = libSingular.p_Copy(a.ptr, R.ptr)
   b1 = libSingular.p_Copy(b.ptr, R.ptr)
   s = libSingular.p_Add_q(a1, b1, R.ptr)
   return svector{T}(R, a.rank, s) 
end

function -(a::svector{T}, b::svector{T}) where T <: AbstractAlgebra.RingElem
   check_parent(a, b)
   R = base_ring(a)
   a1 = libSingular.p_Copy(a.ptr, R.ptr)
   b1 = libSingular.p_Copy(b.ptr, R.ptr)
   s = libSingular.p_Sub(a1, b1, R.ptr)
   return svector{T}(R, a.rank, s) 
end

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

function *(a::svector{spoly{T}}, b::spoly{T}) where T <: AbstractAlgebra.RingElem
   base_ring(a) != parent(b) && error("Incompatible base rings")
   R = base_ring(a)
   a1 = libSingular.p_Copy(a.ptr, R.ptr)
   b1 = libSingular.p_Copy(b.ptr, R.ptr)
   s = libSingular.p_Mult_q(a1, b1, R.ptr)
   return svector{spoly{T}}(R, a.rank, s)
end

function *{T <: Nemo.RingElem}(a::spoly{T}, b::svector{spoly{T}})
   base_ring(b) != parent(a) && error("Incompatible base rings")
   R = base_ring(b)
   a1 = libSingular.p_Copy(a.ptr, R.ptr)
   b1 = libSingular.p_Copy(b.ptr, R.ptr)
   s = libSingular.p_Mult_q(a1, b1, R.ptr)
   return svector{spoly{T}}(R, b.rank, s)
end

*{T <: Nemo.RingElem}(a::svector{spoly{T}}, b::T) = a*base_ring(a)(b)

*{T <: Nemo.RingElem}(a::T, b::svector{spoly{T}}) = base_ring(b)(a)*b

*(a::svector, b::Integer) = a*base_ring(a)(b)

*(a::Integer, b::svector) = base_ring(b)(a)*b

###############################################################################
#
#   Comparison
#
###############################################################################

function =={T <: Nemo.RingElem}(x::svector{T}, y::svector{T})
   check_parent(x, y)
   return Bool(libSingular.p_EqualPolys(x.ptr, y.ptr, base_ring(x).ptr))
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (S::FreeMod{T})(a::Array{T, 1}) where T <: AbstractAlgebra.RingElem
   R = base_ring(S) # polynomial ring
   n = size(a)[1]
   aa = [p.ptr.cpp_object for p in a]
   v = libSingular.id_Array2Vector(reinterpret(Ptr{Void},pointer(aa)), n, base_ring(S).ptr)
   return svector{T}(R, n, v)
end

###############################################################################
#
#   Conversions
#
###############################################################################

function Array{T <: Nemo.RingElem}(v::svector{spoly{T}})
   n = v.rank
   aa_val = Array{Ptr{Void},1}(n)
   R = base_ring(v)
   libSingular.p_Vector2Array(v.ptr, reinterpret(Ptr{Void},pointer(aa_val)), n, R.ptr)
   aa = [libSingular.internal_void_to_poly_helper(p) for p in aa_val]
   return [spoly{T}(R, p) for p in aa]
end

###############################################################################
#
#   Vector constructors
#
###############################################################################

function vector(R::PolyRing{T}, coords::spoly{T}...) where T <: AbstractAlgebra.RingElem
   n = length(coords)
   aa = [p.ptr.cpp_object for p in coords]
   v = libSingular.id_Array2Vector(reinterpret(Ptr{Void},pointer(aa)), n, R.ptr)
   return svector{spoly{T}}(R, n, v)
end

###############################################################################
#
#   Free module constructor
#
###############################################################################

# free module of rank n
function FreeModule{T <: Nemo.RingElem}(R::PolyRing{T}, n::Int)
   (n > typemax(Cint) || n < 0) && throw(DomainError())
   S = elem_type(R)
   return FreeMod{S}(R, n)
end

