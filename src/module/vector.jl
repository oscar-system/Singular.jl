export FreeMod, svector, gens, rank, vector, jet

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent(v::svector{T}) where {T <: Nemo.RingElem} = FreeMod{T}(v.base_ring, v.rank)

base_ring(R::FreeMod) = R.base_ring

base_ring(v::svector) = v.base_ring

elem_type(::FreeMod{T}) where {T <: Nemo.RingElem} = svector{T}

elem_type(::Type{FreeMod{T}}) where {T <: Nemo.RingElem} = svector{T}

parent_type(::Type{svector{T}}) where {T <: Nemo.RingElem} = FreeMod{T}

@doc Markdown.doc"""
    rank(M::FreeMod)

> Return the rank of the given free module.
"""
rank(M::FreeMod) = M.rank

@doc Markdown.doc"""
    gens{T <: AbstractAlgebra.RingElem}(M::FreeMod{T})

> Return a Julia array whose entries are the generators of the given free module.
"""
function gens(M::FreeMod{T}) where T <: AbstractAlgebra.RingElem
   R = base_ring(M)
   ptr = libSingular.id_FreeModule(Cint(M.rank), R.ptr)
   return [svector{T}(R, M.rank, libSingular.getindex(ptr, Cint(i - 1))) for i in 1:M.rank]
end

function deepcopy_internal(p::svector{T}, dict::IdDict) where T <: AbstractAlgebra.RingElem
   p2 = libSingular.p_Copy(p.ptr, base_ring(p).ptr)
   return svector{T}(p.base_ring, p.rank, p2)
end

function check_parent(a::svector{T}, b::svector{T}) where T <: Nemo.RingElem
   base_ring(a) != base_ring(b) && error("Incompatible base rings")
   a.rank != b.rank && error("Vectors of incompatible rank")
end

function hash(V::svector, h::UInt)
   v = 0xd2fd44bc67ee655e%UInt
   for p in Array(V)
      v = xor(hash(p, h), v)
   end
   return v
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
   print(io, m)
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

iszero(p::svector) = p.ptr.cpp_object == C_NULL

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
   s = libSingular.pp_Mult_qq(a.ptr, b.ptr, R.ptr)
   return svector{spoly{T}}(R, a.rank, s)
end

function (a::spoly{T} * b::svector{spoly{T}}) where T <: Nemo.RingElem
   base_ring(b) != parent(a) && error("Incompatible base rings")
   R = base_ring(b)
   s = libSingular.pp_Mult_qq(a.ptr, b.ptr, R.ptr)
   return svector{spoly{T}}(R, b.rank, s)
end

(a::svector{spoly{T}} * b::T) where T <: Nemo.RingElem = a*base_ring(a)(b)

(a::T * b::svector{spoly{T}}) where T <: Nemo.RingElem = base_ring(b)(a)*b

*(a::svector, b::Integer) = a*base_ring(a)(b)

*(a::Integer, b::svector) = base_ring(b)(a)*b

###############################################################################
#
#   Comparison
#
###############################################################################

function (x::svector{T} == y::svector{T}) where T <: Nemo.RingElem
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
   v = libSingular.id_Array2Vector(reinterpret(Ptr{Nothing},pointer(aa)), n, base_ring(S).ptr)
   return svector{T}(R, n, v)
end

function (S::FreeMod{T})() where T <: AbstractAlgebra.RingElem
   R = base_ring(S) # polynomial ring
   n = S.rank
   z = zero(R)
   return svector{T}(R, n, z.ptr)
end


###############################################################################
#
#   Conversions
#
###############################################################################

function Array(v::svector{spoly{T}}) where T <: Nemo.RingElem
   n = v.rank
   aa_val = Array{Ptr{Nothing},1}(undef, n)
   R = base_ring(v)
   libSingular.p_Vector2Array(v.ptr, reinterpret(Ptr{Nothing},pointer(aa_val)), n, R.ptr)
   aa = [libSingular.internal_void_to_poly_helper(p) for p in aa_val]
   return [spoly{T}(R, p) for p in aa]
end

###############################################################################
#
#   iterators
#
###############################################################################

function Base.iterate(p::svector{spoly{T}}) where T <: Nemo.RingElem
   ptr = p.ptr
   if ptr.cpp_object == C_NULL
      return nothing
   else
      R = base_ring(p)
      A = Array{Int}(undef, nvars(R))
      c = libSingular.p_GetExpVLV(ptr, A, R.ptr)
      S = base_ring(R)
      a = S(libSingular.n_Copy(libSingular.pGetCoeff(ptr), S.ptr))
      return (c, A, a), ptr
   end
end

function Base.iterate(p::svector{spoly{T}}, state) where T <: Nemo.RingElem
   state = libSingular.pNext(state)
   if state.cpp_object == C_NULL
      return nothing
   else
      R = base_ring(p)
      A = Array{Int}(undef, nvars(R))
      c = libSingular.p_GetExpVLV(state, A, R.ptr)
      S = base_ring(R)
      a = S(libSingular.n_Copy(libSingular.pGetCoeff(state), S.ptr))
      return (c, A, a), state
   end
end

Base.IteratorSize(::svector{spoly{T}}) where {T} = Base.SizeUnknown()
###############################################################################
#
#   Vector constructors
#
###############################################################################

function vector(R::PolyRing{T}, coords::spoly{T}...) where T <: AbstractAlgebra.RingElem
   n = length(coords)
   aa = [p.ptr.cpp_object for p in coords]
   v = libSingular.id_Array2Vector(reinterpret(Ptr{Nothing},pointer(aa)), n, R.ptr)

   return svector{spoly{T}}(R, n, v)
end

###############################################################################
#
#   Free module constructor
#
###############################################################################

# free module of rank n
function FreeModule(R::PolyRing{T}, n::Int) where T <: Nemo.RingElem
   (n > typemax(Cint) || n < 0) && throw(DomainError())
   S = elem_type(R)
   return FreeMod{S}(R, n)
end

###############################################################################
#
#   Differential functions
#
###############################################################################

@doc Markdown.doc"""
    jet(x::svector{spoly{T}}, n::Int)

> Given a vector $x$ this function truncates each entry of $x$ up to degree $n$.
"""
function jet(x::svector{spoly{T}}, n::Int) where T <: AbstractAlgebra.RingElem
   R = base_ring(x)
   p = libSingular.p_Copy(x.ptr, R.ptr)
   s = libSingular.p_Jet(p, Cint(n), R.ptr)
   return svector{spoly{T}}(R, x.rank, s)
end
