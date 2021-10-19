export slppoly, LPPolyRing, FreeAlgebra

###############################################################################
#
# String I/O
#
###############################################################################

function show(io::IO, R::LPPolyRing)
   GC.@preserve R s = libSingular.rString(R.ptr)
   if libSingular.rIsQuotientRing(R.ptr)
      print(io, "Singular letterplace Quotient Ring ", s)
   else
      print(io, "Singular letterplace Ring ", s)
   end
end

function show(io::IO, a::slppoly)
   R = parent(a)
   GC.@preserve a R s = libSingular.p_String(a.ptr, R.ptr)
   print(io, s)
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent(p::slppoly) = p.parent

base_ring(R::LPPolyRing{T}) where T <: Nemo.RingElem = R.base_ring::parent_type(T)

base_ring(p::slppoly) = base_ring(parent(p))

elem_type(::Type{LPPolyRing{T}}) where T <: Nemo.RingElem = slppoly{T}

parent_type(::Type{slppoly{T}}) where T <: Nemo.RingElem = LPPolyRing{T}

function degree_bound(R::LPPolyRing)
   return R.deg_bound
end

function nvars(R::LPPolyRing)
   # note libSingular.rVar returns nvars*deg_bound
   return length(R.S)
end

has_global_ordering(R::LPPolyRing) = Bool(libSingular.rHasGlobalOrdering(R.ptr))

has_mixed_ordering(R::LPPolyRing) = Bool(libSingular.rHasMixedOrdering(R.ptr))

function has_local_ordering(R::LPPolyRing)
   return !has_global_ordering(R) && !has_mixed_ordering(R)
end

function isgen(p::slppoly)
   R = parent(p)
   GC.@preserve R p begin
      if p.ptr.cpp_object == C_NULL || libSingular.pNext(p.ptr).cpp_object != C_NULL ||
        !Bool(libSingular.n_IsOne(libSingular.pGetCoeff(p.ptr), base_ring(p).ptr))
         return false
      end
      n = 0
      for i = 1:nvars(R)
         d = libSingular.p_GetExp(p.ptr, Cint(i), R.ptr)
         if d > 1
            return false
         elseif d == 1
            n == 0 || return false
            n = 1
         end
      end
      n == 1 || return false
      for i = nvars(R)+1:nvars(R)*degree_bound(R)
         libSingular.p_GetExp(p.ptr, Cint(i), R.ptr) == 0 || return false
      end
      return true
   end
end

function isquotient_ring(R::LPPolyRing)
   GC.@preserve R return Bool(Singular.libSingular.rIsQuotientRing(R.ptr))
end

function characteristic(R::LPPolyRing)
   GC.@preserve R return Int(libSingular.rChar(R.ptr))
end

function gens(R::LPPolyRing)
   n = nvars(R)
   GC.@preserve R return [R(libSingular.rGetVar(Cint(i), R.ptr)) for i = 1:n]
end

function gen(R::LPPolyRing, i::Int)
   GC.@preserve R return R(libSingular.rGetVar(Cint(i), R.ptr))
end

function symbols(R::LPPolyRing)
   return R.S
end

function singular_symbols(R::LPPolyRing)
   GC.@preserve R return singular_symbols(R.ptr)
end

ordering(R::LPPolyRing) = R.ord

zero(R::LPPolyRing) = R()

one(R::LPPolyRing) = R(1)

function iszero(p::slppoly)
   GC.@preserve p p.ptr.cpp_object == C_NULL
end

function isone(p::slppoly)
   R = parent(p)
   GC.@preserve R p return Bool(libSingular.p_IsOne(p.ptr, R.ptr))
end

###############################################################################
#
# Exponents
#
###############################################################################

function _lpvec_to_exponent_word(v::Vector{Int}, nvars::Int, deg_bound::Int)
   z = Int[]
   for d in 0:deg_bound-1
      for i in 1:nvars
         if v[i + d*nvars] != 0
            push!(z, i)
            break
         end
      end
   end
   return z
end

function _exponent_word(p::libSingular.poly_ptr, R::LPPolyRing, tmp::Vector{Int})
   @assert length(tmp) >= nvars(R)*degree_bound(R)
   GC.@preserve R libSingular.p_GetExpVL(p, tmp, R.ptr)
   return _lpvec_to_exponent_word(tmp, nvars(R), degree_bound(R))
end

function exponent_words(x::slppoly)
   R = parent(x)
   return SPolyExponentWords(x, zeros(Int, nvars(R)*degree_bound(R)))
end

function Base.iterate(x::SPolyExponentWords{<:slppoly{T}}) where T <: Nemo.RingElem
   p = x.poly
   ptr = p.ptr
   if ptr.cpp_object == C_NULL
      return nothing
   else
      return _exponent_word(ptr, parent(p), x.tmp), ptr
   end
end

function Base.iterate(x::SPolyExponentWords{<:slppoly{T}}, state) where T <: Nemo.RingElem
   state = libSingular.pNext(state)
   if state.cpp_object == C_NULL
      return nothing
   else
      return _exponent_word(state, parent(x.poly), x.tmp), state
   end
end

function Base.hash(p::slppoly{T}, h::UInt) where T <: Nemo.RingElem
   v = 0xaf708b07f940b4d2%UInt
   for (c, e) in zip(coefficients(p), exponent_words(p))
      v = xor(hash(e, h), v)
      v = xor(hash(c, h), v)
      v = (v << 1) | (v >> (sizeof(Int)*8 - 1))
   end
   return v
end

function _exponent_word_to_lp_vec!(v::Vector{Int}, w::Vector{Int},
                                   nvars::Int, deg_bound::Int)
   length(w) <= deg_bound || error("monomial length exceeds degree bound $deg_bound")
   @assert length(v) >= nvars*deg_bound
   for i in 1:nvars*deg_bound
      v[i] = 0
   end
   for d in 0:length(w)-1
      j = w[1 + d]
      1 <= j <= nvars || error("variable index out of range")
      v[j + d*nvars] = 1
   end
end

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function (R::LPPolyRing{T})() where T <: Nemo.RingElem
   return slppoly{T}(R)
end

function (R::LPPolyRing{T})(n::Int) where T <: Nemo.RingElem
   return slppoly{T}(R, n)
end

function (R::LPPolyRing{T})(n::Integer) where T <: Nemo.RingElem
   return slppoly{T}(R, BigInt(n))
end

function (R::LPPolyRing{T})(n::n_Z) where T <: Nemo.RingElem
   n = base_ring(R)(n)
   ptr = libSingular.n_Copy(n.ptr, parent(n).ptr)
   return slppoly{T}(R, ptr)
end

function (R::LPPolyRing)(n::Rational)
   return R(base_ring(R)(n))
end

function (R::LPPolyRing{T})(n::libSingular.poly_ptr) where T <: Nemo.RingElem
   return slppoly{T}(R, n)
end

function (R::LPPolyRing{T})(n::T) where T <: Nemo.RingElem
   parent(n) != base_ring(R) && error("Unable to coerce into Exterior algebra")
   return slppoly{T}(R, n.ptr)
end

function (R::LPPolyRing{S})(n::T) where {S <: Nemo.RingElem, T <: Nemo.RingElem}
   return R(base_ring(R)(n))
end

function (R::LPPolyRing)(p::slppoly)
   parent(p) != R && error("Unable to coerce")
   return p
end

function (R::LPPolyRing)(n::libSingular.number_ptr)
    return R.base_ring(n)
end

###############################################################################
#
#   FreeAlgebra constructor
#
###############################################################################

function _FreeAlgebra(R, s::Vector{String}, degree_bound, ordering, ordering2, cached)
   T = elem_type(R)
   if isa(ordering, Symbol)
      ord1 = sym2ringorder[ordering]
      ord2 = sym2ringorder[ordering2]
      if (Int(ord1 == ringorder_c || ord1 == ringorder_C) +
          Int(ord2 == ringorder_c || ord2 == ringorder_C) != 1)
         error("ordering from symbol requires exactly one module ordering")
      end
      fancy_ordering = sordering([sorder_block(ord1, 0, Int[])])
      if ord2 != ringorder_C
         push!(fancy_ordering.data, sorder_block(ord2, 0, Int[]))
      end
   elseif isa(ordering, sordering)
      fancy_ordering = ordering
   else
      error("ordering must be a Symbol or an sordering")
   end
   parent_obj = LPPolyRing{T}(R, Symbol.(s), degree_bound, fancy_ordering, cached)
   return (parent_obj, gens(parent_obj))
end

function FreeAlgebra(R::Field, s::Vector{String}, degree_bound::Int;
                     ordering = :degrevlex, ordering2::Symbol = :comp1min,
                     cached::Bool = true)
   return _FreeAlgebra(R, s, degree_bound, ordering, ordering2, cached)
end

