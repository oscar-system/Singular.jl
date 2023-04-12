export slpalg, LPRing, FreeAlgebra

###############################################################################
#
# String I/O
#
###############################################################################

function show(io::IO, R::LPRing)
   GC.@preserve R s = libSingular.rString(R.ptr)
   if libSingular.rIsQuotientRing(R.ptr)
      print(io, "Singular letterplace Quotient Ring ", s)
   else
      print(io, "Singular letterplace Ring ", s)
   end
end

function expressify(a::slpalg, x = symbols(parent(a)); context = nothing)
   sum = Expr(:call, :+)
   for (c, v) in zip(coefficients(a), exponent_words(a))
      prod = Expr(:call, :*)
      if !isone(c)
         push!(prod.args, expressify(c, context = context))
      end
      j = -1
      e = 0
      for i in v
         if j != i
            if j > 0 && !iszero(e)
               push!(prod.args, e == 1 ? x[j] : Expr(:call, :^, x[j], e))
            end
            e = 0
         end
         j = i
         e += 1
      end
      if j > 0 && !iszero(e)
         push!(prod.args, e == 1 ? x[j] : Expr(:call, :^, x[j], e))
      end
      push!(sum.args, prod)
   end
   return sum
end

AbstractAlgebra.@enable_all_show_via_expressify slpalg

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent(p::slpalg) = p.parent

base_ring(R::LPRing{T}) where T <: Nemo.RingElem = R.base_ring::parent_type(T)

base_ring(p::slpalg) = base_ring(parent(p))

elem_type(::Type{LPRing{T}}) where T <: Nemo.RingElem = slpalg{T}

parent_type(::Type{slpalg{T}}) where T <: Nemo.RingElem = LPRing{T}

@doc raw"""
    degree_bound(R::LPRing)

Return the maximum length of a monomial word as requested by the user via the
`degree_bound` parameter of the `FreeAlgebra` constructor.
"""
function degree_bound(R::LPRing)
   return R.deg_bound
end

function is_gen(p::slpalg)
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

# ownership of the pointer is NOT taken - not for general users
function _exponent_word(p::libSingular.poly_ptr, R::LPRing, tmp::Vector{Int})
   @assert length(tmp) >= nvars(R)*degree_bound(R)
   GC.@preserve R libSingular.p_GetExpVL(p, tmp, R.ptr)
   return _lpvec_to_exponent_word(tmp, nvars(R), degree_bound(R))
end

function exponent_words(x::slpalg)
   R = parent(x)
   return SPolyExponentWords(x, zeros(Int, nvars(R)*degree_bound(R)))
end

function Base.iterate(x::SPolyExponentWords{<:slpalg{T}}) where T <: Nemo.RingElem
   p = x.poly
   ptr = p.ptr
   if ptr.cpp_object == C_NULL
      return nothing
   else
      return _exponent_word(ptr, parent(p), x.tmp), ptr
   end
end

function Base.iterate(x::SPolyExponentWords{<:slpalg{T}}, state) where T <: Nemo.RingElem
   state = libSingular.pNext(state)
   if state.cpp_object == C_NULL
      return nothing
   else
      return _exponent_word(state, parent(x.poly), x.tmp), state
   end
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
#   Promote rules
#
###############################################################################

promote_rule(::Type{slpalg{T}}, ::Type{slpalg{T}}) where T <: Nemo.RingElem = slpalg{T}

function promote_rule(::Type{slpalg{T}}, ::Type{U}) where {T <: Nemo.RingElem, U <: Nemo.RingElem}
   promote_rule(T, U) == T ? slpalg{T} : Union{}
end

function promote_rule(::Type{slpalg{T}}, ::Type{T}) where {T <: Nemo.RingElem}
   return slpalg{T}
end

function promote_rule(::Type{slpalg{n_RingElem{RingElemWrapper{S, T}}}}, ::Type{U}) where {T <: Nemo.RingElem, U <: Nemo.RingElem, S}
   return slpalg{n_RingElem{RingElemWrapper{S, T}}}
end

function promote_rule(::Type{slpalg{n_FieldElem{FieldElemWrapper{S, T}}}}, ::Type{U}) where {T <: Nemo.FieldElem, U <: Nemo.FieldElem, S}
   return slpalg{n_FieldElem{FieldElemWrapper{S, T}}}
end

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function (R::LPRing{T})() where T <: Nemo.RingElem
   return slpalg{T}(R)
end

function (R::LPRing{T})(n::Int) where T <: Nemo.RingElem
   return slpalg{T}(R, n)
end

function (R::LPRing{T})(n::Integer) where T <: Nemo.RingElem
   return slpalg{T}(R, BigInt(n))
end

function (R::LPRing{T})(n::n_Z) where T <: Nemo.RingElem
   return slpalg{T}(R, base_ring(R)(n))
end

function (R::LPRing{T})(n::Rational) where T <: Nemo.RingElem
   return slpalg{T}(R, base_ring(R)(n))
end

# take ownership of the pointer - not for general users
function (R::LPRing{T})(ptr::libSingular.poly_ptr) where T <: Nemo.RingElem
   return slpalg{T}(R, ptr)
end

function (R::LPRing{T})(n::T) where T <: Nemo.RingElem
   parent(n) != base_ring(R) && error("Unable to coerce into free algebra")
   return slpalg{T}(R, n)
end

function (R::LPRing{S})(n::T) where {S <: Nemo.RingElem, T <: Nemo.RingElem}
   return R(base_ring(R)(n))
end

function (R::LPRing)(p::slpalg)
   parent(p) !== R && error("Unable to coerce into free algebra")
   return p
end

###############################################################################
#
#   FreeAlgebra constructor
#
###############################################################################

function _FreeAlgebra(R, s::Union{Vector{String},Vector{Symbol}}, degree_bound, ordering, ordering2, cached)
   s = map(Symbol, s)
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
   parent_obj = LPRing{T}(R, s, degree_bound, fancy_ordering, cached)
   return (parent_obj, gens(parent_obj))
end

function FreeAlgebra(R::Union{Ring, Field}, s::Union{Vector{String},Vector{Symbol}}, degree_bound::Int;
                     ordering = :degrevlex, ordering2::Symbol = :comp1min,
                     cached::Bool = true)
   return _FreeAlgebra(R, s, degree_bound, ordering, ordering2, cached)
end

function FreeAlgebra(R::Nemo.Ring, s::Union{Vector{String},Vector{Symbol}}, degree_bound::Int;
                     ordering = :degrevlex, ordering2::Symbol = :comp1min,
                     cached::Bool = true)
   R = CoefficientRing(R)
   return _FreeAlgebra(R, s, degree_bound, ordering, ordering2, cached)
end

