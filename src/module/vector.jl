export FreeMod, svector, gens, rank, vector, jet

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent(v::svector{T}) where {T <: Nemo.NCRingElem} = FreeMod{T}(v.base_ring, v.rank)

base_ring(R::FreeMod) = R.base_ring

base_ring(v::svector) = v.base_ring

elem_type(::Type{FreeMod{T}}) where {T <: Nemo.NCRingElem} = svector{T}

parent_type(::Type{svector{T}}) where {T <: Nemo.NCRingElem} = FreeMod{T}

@doc raw"""
    rank(M::FreeMod)

Return the rank of the given free module.
"""
rank(M::FreeMod) = M.rank

@doc raw"""
    gens{T <: AbstractAlgebra.RingElem}(M::FreeMod{T})

Return a Julia array whose entries are the generators of the given free module.
"""
function gens(M::FreeMod{T}) where T <: Nemo.NCRingElem
   R = base_ring(M)
   ptr = GC.@preserve R libSingular.id_FreeModule(Cint(M.rank), R.ptr)
   return [svector{T}(R, M.rank, libSingular.getindex(ptr, Cint(i - 1))) for i in 1:M.rank]
end

function gen(M::FreeMod{T}, i::Int) where T <: Nemo.NCRingElem
   1 <= i <= M.rank || error("index out of range")
   R = base_ring(M)
   ptr = GC.@preserve R libSingular.id_FreeModule(Cint(M.rank), R.ptr)
   return svector{T}(R, M.rank, libSingular.getindex(ptr, Cint(i - 1)))
end

number_of_generators(M::FreeMod{T}) where T <: Nemo.NCRingElem = rank(M)

function deepcopy_internal(p::svector{T}, dict::IdDict) where T <: Nemo.NCRingElem
   R = base_ring(p)
   p2 = GC.@preserve p R libSingular.p_Copy(p.ptr, R.ptr)
   return svector{T}(p.base_ring, p.rank, p2)
end

function check_parent(a::svector{T}, b::svector{T}) where T <: Nemo.NCRingElem
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
   print(io, "Free module of rank ", R.rank, " over ")
   show(terse(io), R.base_ring)
end

function show(io::IO, a::svector)
   R = base_ring(a)
   m = GC.@preserve a R libSingular.p_String(a.ptr, base_ring(a).ptr)
   print(io, m)
end


###############################################################################
#
#   Unary functions
#
###############################################################################

function -(a::svector{T}) where T <: Nemo.NCRingElem
   R = base_ring(a)
   GC.@preserve a R begin
      a1 = libSingular.p_Copy(a.ptr, R.ptr)
      s = libSingular.p_Neg(a1, R.ptr)
      return svector{T}(R, a.rank, s)
   end
end

iszero(p::svector) = p.ptr.cpp_object == C_NULL

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(a::svector{T}, b::svector{T}) where T <: Nemo.NCRingElem
   check_parent(a, b)
   R = base_ring(a)
   GC.@preserve a b R begin
      a1 = libSingular.p_Copy(a.ptr, R.ptr)
      b1 = libSingular.p_Copy(b.ptr, R.ptr)
      s = libSingular.p_Add_q(a1, b1, R.ptr)
      return svector{T}(R, a.rank, s)
   end
end

function -(a::svector{T}, b::svector{T}) where T <: Nemo.NCRingElem
   check_parent(a, b)
   R = base_ring(a)
   GC.@preserve a b R begin
      a1 = libSingular.p_Copy(a.ptr, R.ptr)
      b1 = libSingular.p_Copy(b.ptr, R.ptr)
      s = libSingular.p_Sub(a1, b1, R.ptr)
      return svector{T}(R, a.rank, s)
   end
end

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

function *(a::svector{spoly{T}}, b::spoly{T}) where T <: AbstractAlgebra.RingElem
   base_ring(a) != parent(b) && error("Incompatible base rings")
   R = base_ring(a)
   s = GC.@preserve a b R libSingular.pp_Mult_qq(a.ptr, b.ptr, R.ptr)
   return svector{spoly{T}}(R, a.rank, s)
end

function *(a::svector{Singular.spluralg{T}}, b::Singular.spluralg{T}) where T <: AbstractAlgebra.RingElem
   base_ring(a) != parent(b) && error("Incompatible base rings")
   R = base_ring(a)
   s = GC.@preserve a b R libSingular.pp_Mult_qq(a.ptr, b.ptr, R.ptr)
   return svector{spoly{T}}(R, a.rank, s)
end

function (a::spoly{T} * b::svector{spoly{T}}) where T <: Nemo.RingElem
   base_ring(b) != parent(a) && error("Incompatible base rings")
   R = base_ring(b)
   s = GC.@preserve a b R libSingular.pp_Mult_qq(a.ptr, b.ptr, R.ptr)
   return svector{spoly{T}}(R, b.rank, s)
end

function *(b::Singular.spluralg{T}, a::svector{Singular.spluralg{T}}) where T
   base_ring(a) != parent(b) && error("Incompatible base rings")
   R = base_ring(a)
   s = GC.@preserve a b R libSingular.pp_Mult_qq(a.ptr, b.ptr, R.ptr)
   return svector{spluralg{T}}(R, a.rank, s)
end

(a::svector{spoly{T}} * b::T) where T <: Nemo.RingElem = a*base_ring(a)(b)

(a::T * b::svector{spoly{T}}) where T <: Nemo.RingElem = base_ring(b)(a)*b

#(a::Singular.spluralg * b::svector{Singular.spluralg}) = base_ring(b)(a)*b

*(a::svector, b::Integer) = a*base_ring(a)(b)

*(a::Integer, b::svector) = base_ring(b)(a)*b

###############################################################################
#
#   Comparison
#
###############################################################################

function (x::svector{T} == y::svector{T}) where T <: Union{Nemo.RingElem, Singular.spluralg}
   check_parent(x, y)
   R = base_ring(x)
   GC.@preserve x y R return Bool(libSingular.p_EqualPolys(x.ptr, y.ptr, R.ptr))
end

###############################################################################
#
#   Leading terms
#
###############################################################################

@doc raw"""
    lead(a::svector{T}) where T <: Nemo.RingElem

Return the initial terms of $a$.
"""
function lead(a::svector{T}) where T <: Nemo.RingElem
   R = base_ring(a)
   ptr = GC.@preserve a R libSingular.p_Head(a.ptr, R.ptr)
   return svector{T}(R, a.rank, ptr)
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (S::FreeMod{T})(a::Vector{T}) where T <: Union{AbstractAlgebra.RingElem, Singular.spluralg}
   R = base_ring(S) # polynomial ring
   GC.@preserve a R begin
      n = size(a)[1]
      aa = [p.ptr.cpp_object for p in a]
      v = libSingular.id_Array2Vector(reinterpret(Ptr{Nothing},pointer(aa)), n, R.ptr)
      return svector{T}(R, n, v)
   end
end

function (S::FreeMod{T})() where T <: Union{AbstractAlgebra.RingElem, Singular.spluralg}
   R = base_ring(S) # polynomial ring
   n = S.rank
   z = zero(R)
   GC.@preserve z return svector{T}(R, n, z.ptr)
end


###############################################################################
#
#   Conversions
#
###############################################################################

function Base.Array(v::svector{spoly{T}}) where T <: Nemo.RingElem
   R = base_ring(v)
   GC.@preserve v R begin
      n = v.rank
      aa_val = Vector{Ptr{Nothing}}(undef, n)
      libSingular.p_Vector2Array(v.ptr, reinterpret(Ptr{Nothing},pointer(aa_val)), n, R.ptr)
      aa = [libSingular.internal_void_to_poly_helper(p) for p in aa_val]
      return [spoly{T}(R, p) for p in aa]
   end
end

###############################################################################
#
#   iterators
#
###############################################################################

function Base.iterate(p::svector{spoly{T}}) where T <: Nemo.RingElem
   GC.@preserve p begin
      ptr = p.ptr
      if ptr.cpp_object == C_NULL
         return nothing
      else
         R = base_ring(p)
         A = Array{Int}(undef, nvars(R))
         c = GC.@preserve R libSingular.p_GetExpVLV(ptr, A, R.ptr)
         S = base_ring(R)
         a = GC.@preserve S S(libSingular.n_Copy(libSingular.pGetCoeff(ptr), S.ptr))
         return (c, A, a), ptr
      end
   end
end

function Base.iterate(p::Singular.svector{Singular.spluralg{T}}) where T 
   GC.@preserve p begin
      ptr = p.ptr
      if ptr.cpp_object == C_NULL
         return nothing
      else
         R = base_ring(p)
         A = Array{Int}(undef, nvars(R))
         c = GC.@preserve R libSingular.p_GetExpVLV(ptr, A, R.ptr)
         S = base_ring(R)
         a = GC.@preserve S S(libSingular.n_Copy(libSingular.pGetCoeff(ptr), S.ptr))
         return (c, A, a), ptr
      end
   end
end

function Base.iterate(p::svector{spoly{T}}, state) where T <: Nemo.RingElem
   GC.@preserve p begin
      state = libSingular.pNext(state)
      if state.cpp_object == C_NULL
         return nothing
      else
         R = base_ring(p)
         A = Array{Int}(undef, nvars(R))
         c = GC.@preserve R libSingular.p_GetExpVLV(state, A, R.ptr)
         S = base_ring(R)
         a = GC.@preserve S S(libSingular.n_Copy(libSingular.pGetCoeff(state), S.ptr))
         return (c, A, a), state
      end
   end
end

function Base.iterate(p::svector{spluralg{T}}, state) where T
   GC.@preserve p begin
      state = libSingular.pNext(state)
      if state.cpp_object == C_NULL
         return nothing
      else
         R = base_ring(p)
         A = Array{Int}(undef, nvars(R))
         c = GC.@preserve R libSingular.p_GetExpVLV(state, A, R.ptr)
         S = base_ring(R)
         a = GC.@preserve S S(libSingular.n_Copy(libSingular.pGetCoeff(state), S.ptr))
         return (c, A, a), state
      end
   end
end

Base.IteratorSize(::svector{spoly{T}}) where {T} = Base.SizeUnknown()
###############################################################################
#
#   Vector constructors
#
###############################################################################

function vector(R::PolyRing{T}, coords::spoly{T}...) where T <: AbstractAlgebra.RingElem
   GC.@preserve R coords begin
      n = length(coords)
      aa = [p.ptr.cpp_object for p in coords]
      v = libSingular.id_Array2Vector(reinterpret(Ptr{Nothing},pointer(aa)), n, R.ptr)
      return svector{spoly{T}}(R, n, v)
   end
end

###############################################################################
#
#   Free module constructor
#
###############################################################################

# free module of rank n
function FreeModule(R::PolyRing{T}, n::Int) where T <: Nemo.RingElem
   (n > typemax(Cint) || n < 0) &&
      throw(DomainError(n, "rank must be non-negative and <= $(typemax(Cint))"))
   S = elem_type(R)
   return FreeMod{S}(R, n)
end

function FreeModule(R::Singular.PluralRing{T}, n::Int) where T <: Singular.n_Q
   (n > typemax(Cint) || n < 0) &&
      throw(DomainError(n, "rank must be non-negative and <= $(typemax(Cint))"))
   S = elem_type(R)
   return FreeMod{S}(R, n)
end

###############################################################################
#
#   Differential functions
#
###############################################################################

@doc raw"""
    jet(x::svector{spoly{T}}, n::Int)

Given a vector $x$ this function truncates each entry of $x$ up to degree $n$.
"""
function jet(x::svector{spoly{T}}, n::Int) where T <: AbstractAlgebra.RingElem
   R = base_ring(x)
   GC.@preserve x R begin
      p = libSingular.p_Copy(x.ptr, R.ptr)
      s = libSingular.p_Jet(p, Cint(n), R.ptr)
   end
   return svector{spoly{T}}(R, x.rank, s)
end
