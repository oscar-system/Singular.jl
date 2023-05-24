export spoly, PolyRing, change_base_ring, coeff, coefficients,
       constant_coefficient, content, deflation, deflate,
       degree, degrees, degree_bound,
       derivative, div, divides, evaluate, exponent,
       exponent_vectors, exponent_words, factor, factor_squarefree, finish, gen,
       has_global_ordering, has_mixed_ordering, has_local_ordering,
       homogenize,
       inflate, is_gen, is_monomial, is_ordering_symbolic, is_term,
       jacobian_ideal, jacobian_matrix, jet,
       leading_coefficient, leading_exponent_vector, leading_term,
       leading_monomial, lead_exponent,
       monomials, MPolyBuildCtx,
       nvars, order, ordering, ordering_as_symbol, ordering_size, ordering_weights,
       @polynomial_ring, primpart, push_term!,
       remove, sort_terms!, symbols,
       tail, terms, total_degree, trailing_coefficient,
       valuation, var_index, vars

export ordering_lp, ordering_rp, ordering_dp, ordering_Dp, ordering_wp, ordering_Wp,
       ordering_ls, ordering_rs, ordering_ds, ordering_Ds, ordering_ws, ordering_Ws,
       ordering_a, ordering_M, ordering_c, ordering_C, ordering_s, ordering_S


###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent(p::SPolyUnion) = p.parent

base_ring(R::PolyRing{T}) where T <: Nemo.RingElem = R.base_ring::parent_type(T)

base_ring(p::spoly) = base_ring(parent(p))

elem_type(::Type{PolyRing{T}}) where T <: Nemo.RingElem = spoly{T}

parent_type(::Type{spoly{T}}) where T <: Nemo.RingElem = PolyRing{T}

@doc raw"""
    has_global_ordering(R::PolyRingUnion)

Return `true` if the given ring has a global ordering, i.e. if $1 < x$ for
each variable $x$ in the ring. This include `:lex`, `:deglex` and `:degrevlex`
orderings.
"""
function has_global_ordering(R::PolyRingUnion)
   GC.@preserve R return Bool(libSingular.rHasGlobalOrdering(R.ptr))
end

@doc raw"""
    has_mixed_ordering(R::PolyRingUnion)

Return `true` if the given ring has a mixed ordering, i.e. if $1 < x_i$ for
a variable $x_i$ and $1>x_j$ for another variable $x_j$.
"""
function has_mixed_ordering(R::PolyRingUnion)
   GC.@preserve R return Bool(libSingular.rHasMixedOrdering(R.ptr))
end

@doc raw"""
    has_local_ordering(R::PolyRingUnion)

Return `true` if the given ring has a local ordering, i.e. if $1 > x$ for
all variables $x$.
"""
function has_local_ordering(R::PolyRingUnion)
   return !has_global_ordering(R) && !has_mixed_ordering(R)
end

function characteristic(R::PolyRingUnion)
   GC.@preserve R return Int(libSingular.rChar(R.ptr))
end

function nvars(R::PolyRingUnion)
   return length(R.S)
end

function gens(R::PolyRingUnion)
   n = nvars(R)
   GC.@preserve R return [R(libSingular.rGetVar(Cint(i), R.ptr)) for i = 1:n]
end

function gen(R::PolyRingUnion, i::Int)
   1 <= i <= nvars(R) || error("index out of range")
   GC.@preserve R return R(libSingular.rGetVar(Cint(i), R.ptr))
end

function symbols(R::PolyRingUnion)
   return R.S
end

# ownership of the pointer is NOT taken - not for general users
function singular_symbols(r::libSingular.ring_ptr)
   n = Int(libSingular.rIsLPRing(r))
   n = n > 0 ? n : Int(libSingular.rVar(r))
   return [Symbol(libSingular.rRingVar(Cshort(i - 1), r)) for i in 1:n]
end

function singular_symbols(R::PolyRingUnion)
   GC.@preserve R return singular_symbols(R.ptr)
end

function ordering(R::PolyRingUnion)
   return R.ord
end

@doc raw"""
    is_ordering_symbolic(R::PolyRing)

Return `true` if the ordering of `R` can be represented as a symbol.
"""
is_ordering_symbolic(R::PolyRing) = is_ordering_symbolic(R.ord)

@doc raw"""
    ordering_as_symbol(R::PolyRing)

Assuming the ordering of `R` can be represented as a symbol, return that symbol.
"""
ordering_as_symbol(R::PolyRing) = ordering_as_symbol(R.ord)

@doc raw"""
    degree_bound(R::PolyRing)

Return the internal degree bound in each variable, enforced by Singular. This is the
largest positive value any degree can have before an overflow will occur. This
internal bound may be higher than the bound requested by the user via the
`degree_bound` parameter of the `polynomial_ring` constructor.
"""
function degree_bound(R::PolyRing)
   GC.@preserve R return Int(libSingular.rBitmask(R.ptr))
end

function zero(R::PolyRingUnion)
   return R()
end

function one(R::PolyRingUnion)
   return R(1)
end

function iszero(p::SPolyUnion)
   GC.@preserve p p.ptr.cpp_object == C_NULL
end

function isone(p::SPolyUnion)
   R = parent(p)
   GC.@preserve R p return Bool(libSingular.p_IsOne(p.ptr, R.ptr))
end

# return 0 for non-generators, otherwise the index of the variable
function _var_index(p::Union{spoly, spluralg})
   R = parent(p)
   GC.@preserve R p begin
      if p.ptr.cpp_object == C_NULL || libSingular.pNext(p.ptr).cpp_object != C_NULL ||
        !Bool(libSingular.n_IsOne(libSingular.pGetCoeff(p.ptr), base_ring(p).ptr))
         return 0
      end
      idx = 0
      for i = 1:nvars(R)
         d = libSingular.p_GetExp(p.ptr, Cint(i), R.ptr)
         if d > 1
            return 0
         elseif d == 1
            if idx > 0
               return 0
            end
            idx = i
         end
      end
      return idx
   end
end

function is_gen(p::Union{spoly, spluralg})
   i = _var_index(p)
   return i > 0
end

function var_index(p::Union{spoly, spluralg})
   i = _var_index(p)
   i > 0 || error("$p is not a variable")
   return i
end

function is_constant(p::SPolyUnion)
   R = parent(p)
   GC.@preserve R p begin
      if p.ptr.cpp_object == C_NULL
         return true
      end
      if libSingular.pNext(p.ptr).cpp_object != C_NULL
         return false
      end
      for i = 1:nvars(R)
         if libSingular.p_GetExp(p.ptr, Cint(i), R.ptr) != 0
            return false
         end
      end
      return true
   end
end

function is_term(p::SPolyUnion)
   GC.@preserve p begin
      return p.ptr.cpp_object != C_NULL &&
             libSingular.pNext(p.ptr).cpp_object == C_NULL
   end
end

function is_unit(p::SPolyUnion)
   GC.@preserve p begin
      return p.ptr.cpp_object != C_NULL &&
             libSingular.pNext(p.ptr).cpp_object == C_NULL &&
             Bool(libSingular.p_IsUnit(p.ptr, parent(p).ptr))
   end
end


length(p::SPolyUnion) = Int(libSingular.pLength(p.ptr))

@doc raw"""
    total_degree(p::spoly)

Return the total degree (largest sum of exponents of any monomial) of $p$.
"""
function total_degree(p::spoly)
   R = parent(p)
   GC.@preserve p R return libSingular.pLDeg(p.ptr, R.ptr)
end

@doc raw"""
    order(p::spoly)

Returns the order of $p$.
"""
function order(p::spoly)
   if p.ptr.cpp_object == C_NULL
      return -1
   end

   R = parent(p)
   x = deepcopy(p)
   GC.@preserve R x begin
      hptr = Singular.libSingular.p_Head(x.ptr, R.ptr)
      ord = Singular.libSingular.pLDeg(hptr, R.ptr)
      xptr = libSingular.pNext(p.ptr)
      while xptr.cpp_object != C_NULL
         hptr = Singular.libSingular.p_Head(xptr, R.ptr)
         ord = min(ord, Singular.libSingular.pLDeg(hptr, R.ptr))
         xptr = libSingular.pNext(xptr)
      end
      return ord
   end
end

function leading_exponent_vector(p::Union{spoly, spluralg})
   R = parent(p)
   n = nvars(R)
   if iszero(p)
      throw(ArgumentError("Zero polynomial does not have a leading exponent vector"))
   end
   A = Array{Int}(undef, n)
   GC.@preserve p R libSingular.p_GetExpVL(p.ptr, A, R.ptr)
   return A
end

@deprecate lead_exponent(p::spoly) leading_exponent_vector(p)

function leading_coefficient(p::SPolyUnion)
   R = base_ring(p)
   if p.ptr.cpp_object == C_NULL
      return zero(R)
   else
      return R(libSingular.n_Copy(libSingular.pGetCoeff(p.ptr), R.ptr))
   end
end

function trailing_coefficient(p::SPolyUnion)
   R = base_ring(p)
   GC.@preserve p R begin
      P = p.ptr
      while P.cpp_object != C_NULL
         Q = libSingular.pNext(P)
         if Q.cpp_object == C_NULL
            return R(libSingular.n_Copy(libSingular.pGetCoeff(P), R.ptr))
         end
         P = Q
      end
      return zero(R)
   end
end

function constant_coefficient(p::SPolyUnion)
   R = base_ring(p)
   S = parent(p)
   GC.@preserve p R S begin
      P = p.ptr
      while P.cpp_object != C_NULL
         if libSingular.p_LmIsConstant(P, S.ptr)
            return R(libSingular.n_Copy(libSingular.pGetCoeff(P), R.ptr))
         end
         P = libSingular.pNext(P)
      end
      return zero(R)
   end
end

function leading_term(p::SPolyUnion)
   R = parent(p)
   GC.@preserve p R begin
      P = p.ptr
      if P.cpp_object == C_NULL
         throw(ArgumentError("Zero polynomial does not have a leading term"))
      end
      return  R(libSingular.p_Head(P, R.ptr))
   end
end

function leading_monomial(p::SPolyUnion)
   R = parent(p)
   GC.@preserve p R begin
      P = p.ptr
      if P.cpp_object == C_NULL
         throw(ArgumentError("Zero polynomial does not have a leading exponent vector"))
      end
      B = base_ring(R)
      t = one(B)
      GC.@preserve B t begin
         mptr = libSingular.p_Head(P, R.ptr)
         n = libSingular.n_Copy(t.ptr, B.ptr)
         libSingular.p_SetCoeff0(mptr, n, R.ptr)
         return R(mptr)
      end
   end
end

function tail(p::SPolyUnion)
   R = parent(p)
   GC.@preserve p R begin
      P = p.ptr
      if P.cpp_object == C_NULL
         return zero(R)
      else
         return R(libSingular.p_Copy(libSingular.pNext(P), R.ptr))
      end
   end
end

function deepcopy_internal(p::SPolyUnion, dict::IdDict)
   R = parent(p)
   p2 = GC.@preserve p R libSingular.p_Copy(p.ptr, R.ptr)
   return R(p2)
end

function check_parent(a::SPolyUnion{T}, b::SPolyUnion{T}) where T <: Nemo.RingElem
   parent(a) != parent(b) && error("Incompatible parent objects")
end

function canonical_unit(a::SPolyUnion{T}) where T <: Nemo.RingElem
  return iszero(a) ? one(base_ring(a)) : canonical_unit(leading_coefficient(a))
end

function Base.hash(p::SPolyUnion, h::UInt)
   z = 0xaf708b07f940b4d2%UInt
   z = bitrotate(xor(z, hash(parent(p), h)), 1)
   es = p isa slpalg ? exponent_words(p) : exponent_vectors(p)
   for (c, e) in zip(coefficients(p), es)
      z = bitrotate(xor(z, hash(e, h)), 1)
      z = bitrotate(xor(z, hash(c, h)), 1)
   end
   return z
end

###############################################################################
#
#   Iterators
#
###############################################################################

function Base.length(x::Union{SPolyCoeffs, SPolyExponentVectors, SPolyExponentWords, SPolyTerms, SPolyMonomials})
   return length(x.poly)
end

function Base.eltype(x::SPolyCoeffs{T}) where T <: SPolyUnion{S} where S
   return S
end

function Base.eltype(x::SPolyExponentVectors)
   return Vector{Int}
end

function Base.eltype(x::SPolyExponentWords)
   return Vector{Int}
end

function Base.eltype(x::SPolyMonomials{T}) where T
   return T
end

function Base.eltype(x::SPolyTerms{T}) where T
   return T
end

function coefficients(x::SPolyUnion)
   return SPolyCoeffs(x)
end

function exponent_vectors(x::SPolyUnion)
   return SPolyExponentVectors(x)
end

function terms(x::SPolyUnion)
   return SPolyTerms(x)
end

function monomials(x::SPolyUnion)
   return SPolyMonomials(x)
end

function Base.iterate(x::SPolyCoeffs{<:SPolyUnion{T}}) where T <: Nemo.RingElem
   p = x.poly
   ptr = p.ptr
   if ptr.cpp_object == C_NULL
      return nothing
   else
      R = base_ring(p)
      GC.@preserve R return R(libSingular.n_Copy(libSingular.pGetCoeff(ptr), R.ptr)), ptr
   end
end

function Base.iterate(x::SPolyCoeffs{<:SPolyUnion{T}}, state) where T <: Nemo.RingElem
   state = libSingular.pNext(state)
   if state.cpp_object == C_NULL
      return nothing
   else
      R = base_ring(x.poly)
      GC.@preserve R return R(libSingular.n_Copy(libSingular.pGetCoeff(state), R.ptr)), state
   end
end

function Base.iterate(x::SPolyExponentVectors{<:SPolyUnion{T}}) where T <: Nemo.RingElem
   p = x.poly
   ptr = p.ptr
   if ptr.cpp_object == C_NULL
      return nothing
   else
      R = parent(p)
      A = Array{Int}(undef, nvars(R))
      GC.@preserve R libSingular.p_GetExpVL(ptr, A, R.ptr)
      return A, ptr
   end
end

function Base.iterate(x::SPolyExponentVectors{<:SPolyUnion{T}}, state) where T <: Nemo.RingElem
   state = libSingular.pNext(state)
   if state.cpp_object == C_NULL
      return nothing
   else
      R = parent(x.poly)
      A = Array{Int}(undef, nvars(R))
      GC.@preserve R libSingular.p_GetExpVL(state, A, R.ptr)
      return A, state
   end
end

function Base.iterate(x::SPolyTerms{<:SPolyUnion{T}}) where T <: Nemo.RingElem
   p = x.poly
   ptr = p.ptr
   if ptr.cpp_object == C_NULL
      return nothing
   else
      R = parent(p)
      GC.@preserve R return R(libSingular.p_Head(ptr, R.ptr)), ptr
   end
end

function Base.iterate(x::SPolyTerms{<:SPolyUnion{T}}, state) where T <: Nemo.RingElem
   state = libSingular.pNext(state)
   if state.cpp_object == C_NULL
      return nothing
   else
      R = parent(x.poly)
      GC.@preserve R return R(libSingular.p_Head(state, R.ptr)), state
   end
end

function Base.iterate(x::SPolyMonomials{<:SPolyUnion{T}}) where T <: Nemo.RingElem
   p = x.poly
   ptr = p.ptr
   if ptr.cpp_object == C_NULL
      return nothing
   else
      S = parent(p)
      R = base_ring(p)
      t = one(R)
      GC.@preserve t R S begin
         mptr = libSingular.p_Head(ptr, S.ptr)
         n = libSingular.n_Copy(t.ptr, R.ptr)
         libSingular.p_SetCoeff0(mptr, n, S.ptr)
         return S(mptr), ptr
      end
   end
end

function Base.iterate(x::SPolyMonomials{<:SPolyUnion{T}}, state) where T <: Nemo.RingElem
   state = libSingular.pNext(state)
   if state.cpp_object == C_NULL
      return nothing
   else
      p = x.poly
      S = parent(p)
      R = base_ring(p)
      t = one(R)
      GC.@preserve t R S begin
         mptr = libSingular.p_Head(state, S.ptr)
         n = libSingular.n_Copy(t.ptr, R.ptr)
         libSingular.p_SetCoeff0(mptr, n, S.ptr)
         return S(mptr), state
      end
   end
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::PolyRing)
   s = libSingular.rString(R.ptr)
   if libSingular.rIsQuotientRing(R.ptr)
      print(io, "Singular Polynomial Quotient Ring ", s)
   else
      print(io, "Singular Polynomial Ring ", s)
   end
end

function show(io::IO, ::MIME"text/plain", R::PolyRing)
   s = libSingular.rString(R.ptr)
   if libSingular.rIsQuotientRing(R.ptr)
      print(io, "Singular Polynomial Quotient Ring ", s)
   else
      print(io, "Singular Polynomial Ring ", s)
   end
end

function expressify(a::Union{spoly, spluralg}, x = symbols(parent(a)); context = nothing)
   sum = Expr(:call, :+)
   for (c, v) in zip(coefficients(a), exponent_vectors(a))
      prod = Expr(:call, :*)
      if !isone(c)
         push!(prod.args, expressify(c, context = context))
      end
      for i in 1:length(v)
         if v[i] > 1
            push!(prod.args, Expr(:call, :^, x[i], v[i]))
         elseif v[i] == 1
            push!(prod.args, x[i])
         end
      end
      push!(sum.args, prod)
   end
   return sum
end

AbstractAlgebra.@enable_all_show_via_expressify spoly
AbstractAlgebra.@enable_all_show_via_expressify spluralg

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(a::SPolyUnion)
   R = parent(a)
   GC.@preserve a R begin
      a1 = libSingular.p_Copy(a.ptr, R.ptr)
      s = libSingular.p_Neg(a1, R.ptr)
      return R(s)
   end
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(a::SPolyUnion{T}, b::SPolyUnion{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   GC.@preserve a b R begin
      a1 = libSingular.p_Copy(a.ptr, R.ptr)
      b1 = libSingular.p_Copy(b.ptr, R.ptr)
      s = libSingular.p_Add_q(a1, b1, R.ptr)
      return R(s)
   end
end

function +(a::SPolyUnion{T}, b::T) where T <: Nemo.RingElem
   R = parent(a)
   S = base_ring(a)
   S != parent(b) && error("Incompatible rings")
   GC.@preserve a b R S begin
      a1 = libSingular.p_Copy(a.ptr, R.ptr)
      b1 = libSingular.p_NSet(libSingular.n_Copy(b.ptr, S.ptr), R.ptr)
      z = libSingular.p_Add_q(a1, b1, R.ptr)
      return R(z)
   end
end

function +(a::T, b::SPolyUnion{T}) where T <: Nemo.RingElem
   return b + a
end

function -(a::SPolyUnion{T}, b::SPolyUnion{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   GC.@preserve a b R begin
      a1 = libSingular.p_Copy(a.ptr, R.ptr)
      b1 = libSingular.p_Copy(b.ptr, R.ptr)
      s = libSingular.p_Sub(a1, b1, R.ptr)
      return R(s)
   end
end

function -(a::SPolyUnion{T}, b::T) where T <: Nemo.RingElem
   R = parent(a)
   S = base_ring(a)
   S != parent(b) && error("Incompatible rings")
   GC.@preserve a b R S begin
      a1 = libSingular.p_Copy(a.ptr, R.ptr)
      b1 = libSingular.p_NSet(libSingular.n_Copy(b.ptr, S.ptr), R.ptr)
      z = libSingular.p_Sub(a1, b1, R.ptr)
      return R(z)
   end
end

function -(a::T, b::SPolyUnion{T}) where T <: Nemo.RingElem
   R = parent(b)
   S = base_ring(b)
   S != parent(a) && error("Incompatible rings")
   GC.@preserve a b R S begin
      b1 = libSingular.p_Copy(b.ptr, R.ptr)
      a1 = libSingular.p_NSet(libSingular.n_Copy(a.ptr, S.ptr), R.ptr)
      z = libSingular.p_Sub(a1, b1, R.ptr)
      return R(z)
   end
end

function *(a::SPolyUnion{T}, b::SPolyUnion{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   s = GC.@preserve a b R libSingular.pp_Mult_qq(a.ptr, b.ptr, R.ptr)
   z = R(s)
   if a isa slpalg
      libSingular.check_error()
   end
   return z
end

function *(a::SPolyUnion{T}, b::T) where T <: Nemo.RingElem
   R = parent(a)
   base_ring(a) != parent(b) && error("Incompatible rings")
   GC.@preserve a b R begin
      a1 = libSingular.p_Copy(a.ptr, R.ptr)
      p = libSingular.p_Mult_nn(a1, b.ptr, R.ptr)
      return R(p)
   end
end

function *(a::T, b::SPolyUnion{T}) where T <: Nemo.RingElem
   return b*a  # multiplication by coeffs should be commutative
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::SPolyUnion, y::Int)
   0 <= y <= typemax(Cint) ||
      throw(DomainError(y, "exponent must be non-negative and <= $(typemax(Cint))"))
   R = parent(x)
   if isone(x)
      return deepcopy(x)
   elseif y == 0
      return one(R)
   elseif y == 1
      return deepcopy(x)
   end
   GC.@preserve x R begin
      x1 = libSingular.p_Copy(x.ptr, R.ptr)
      p = libSingular.p_Power(x1, Cint(y), R.ptr)
      z = R(p)
      if x isa slpalg
         libSingular.check_error()
      end
      return z
   end
end

###############################################################################
#
#   Comparison
#
###############################################################################

function (x::SPolyUnion{T} == y::SPolyUnion{T}) where T <: Nemo.RingElem
   check_parent(x, y)
   GC.@preserve x y return Bool(libSingular.p_EqualPolys(x.ptr, y.ptr, parent(x).ptr))
end

function Base.isless(x::SPolyUnion{T}, y::SPolyUnion{T}) where T <: Nemo.RingElem
   check_parent(x, y)
   R = parent(x)
   GC.@preserve x y R return libSingular.p_Compare(x.ptr, y.ptr, R.ptr) < 0
end

function ==(a::SPolyUnion{T}, b::T) where T <: Nemo.RingElem
   return a == parent(a)(b)
end

function ==(a::T, b::SPolyUnion{T}) where T <: Nemo.RingElem
   return parent(b)(a) == a
end

function ==(a::SPolyUnion, b::Union{Integer, ZZRingElem, Rational, QQFieldElem})
   return a == parent(a)(b)
end

function ==(a::Union{Integer, ZZRingElem, Rational, QQFieldElem}, b::SPolyUnion)
   return parent(b)(a) == a
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::SPolyUnion, y::SPolyUnion; check::Bool=true)
   check_parent(x, y)
   R = parent(x)
   GC.@preserve x y R begin
      x1 = libSingular.p_Copy(x.ptr, R.ptr)
      y1 = libSingular.p_Copy(y.ptr, R.ptr)
      p = libSingular.p_Divide(x1, y1, R.ptr)
      return R(p)
   end
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(x::SPolyUnion{T}, y::T; check::Bool=true) where T <: Nemo.RingElem
   R = parent(x)
   base_ring(x) != parent(y) && error("Incompatible rings")
   GC.@preserve x y R begin
      x1 = libSingular.p_Copy(x.ptr, R.ptr)
      p = libSingular.p_Div_nn(x1, y.ptr, R.ptr)
      return R(p)
   end
end

function divexact(x::SPolyUnion, y::n_Z; check::Bool=true)
   y1 = base_ring(x)(y)
   R = parent(x)
   GC.@preserve x y1 R begin
      x1 = libSingular.p_Copy(x.ptr, R.ptr)
      p = libSingular.p_Div_nn(x1, y1.ptr, R.ptr)
      return R(p)
   end
end

function divexact(x::SPolyUnion, y::n_Q; check::Bool=true)
   y1 = base_ring(x)(y)
   R = parent(x)
   GC.@preserve x y1 R begin
      x1 = libSingular.p_Copy(x.ptr, R.ptr)
      p = libSingular.p_Div_nn(x1, y1.ptr, R.ptr)
      return R(p)
   end
end

function divexact(x::SPolyUnion, y::Int; check::Bool=true)
   R = base_ring(x)
   S = parent(x)
   GC.@preserve x R S begin
      y1 = libSingular.n_Init(y, R.ptr)
      x1 = libSingular.p_Copy(x.ptr, S.ptr)
      p = libSingular.p_Div_nn(x1, y1, S.ptr)
      libSingular.n_Delete(y1, R.ptr)
      return S(p)
   end
end

################################################################################
#
#   Ad hoc binary
#
################################################################################

# We cannot use the promote_rule mechanism, since n_Q and QQFieldElem have no and
# should not have a promote_rule

# cast the integers/rationals to just the coefficient ring, as more efficient
# scalar binary operators are available
function +(x::SPolyUnion, y::Union{Integer, ZZRingElem, Rational, QQFieldElem})
  return x + base_ring(x)(y)
end

function +(x::Union{Integer, ZZRingElem, Rational, QQFieldElem}, y::SPolyUnion)
  return base_ring(y)(x) + y
end

function -(x::SPolyUnion, y::Union{Integer, ZZRingElem, Rational, QQFieldElem})
  return x - base_ring(x)(y)
end

function -(x::Union{Integer, ZZRingElem, Rational, QQFieldElem}, y::SPolyUnion)
  return base_ring(y)(x) - y
end

function *(x::SPolyUnion, y::Union{Integer, ZZRingElem, Rational, QQFieldElem})
  return x * base_ring(x)(y)
end

function *(x::Union{Integer, ZZRingElem, Rational, QQFieldElem}, y::SPolyUnion)
  return base_ring(y)(x) * y
end

function divexact(x::SPolyUnion, y::Union{Integer, ZZRingElem, Rational, QQFieldElem}; check::Bool=true)
  return divexact(x, base_ring(x)(y), check=check)
end

###############################################################################
#
#   Divisibility testing
#
###############################################################################

function divides(x::SPolyUnion{T}, y::SPolyUnion{T}) where T <: Nemo.FieldElem
   check_parent(x, y)
   R = parent(x)
   GC.@preserve x y R begin
      # First check divisibility in a cheap way
      x2 = libSingular.p_Copy(x.ptr, R.ptr)
      y2 = libSingular.p_Copy(y.ptr, R.ptr)
      flag = libSingular.p_IsDivisibleBy(x2, y2, R.ptr)
      if flag
         # now compute exact quotient
         x2 = libSingular.p_Copy(x.ptr, R.ptr)
         y2 = libSingular.p_Copy(y.ptr, R.ptr)
         q = libSingular.p_Divide(x2, y2, R.ptr)
         return true, R(q)
      else
         return false, R()
      end
   end
end

###############################################################################
#
#   Euclidean division
#
###############################################################################

function divrem(x::SPolyUnion{T}, y::SPolyUnion{T}) where T <: Nemo.FieldElem
   check_parent(x, y)
   R = parent(x)
   GC.@preserve x y R begin
      px = libSingular.p_Copy(x.ptr, R.ptr)
      py = libSingular.p_Copy(y.ptr, R.ptr)
      q, r = libSingular.p_DivRem(px, py, R.ptr)
      qref = libSingular.toPolyRef(q)
      rref = libSingular.toPolyRef(r)
      return R(qref), R(rref)
   end
end

function div(x::SPolyUnion{T}, y::SPolyUnion{T}) where T <: Nemo.FieldElem
   check_parent(x, y)
   R = parent(x)
   GC.@preserve x y R begin
      px = libSingular.p_Copy(x.ptr, R.ptr)
      py = libSingular.p_Copy(y.ptr, R.ptr)
      q = R(libSingular.p_Divide(px, py, R.ptr))
      libSingular.check_error()
      return q
   end
end

function divrem(a::spoly{T}, b::spoly{T}) where T <: Nemo.FieldElem
    check_parent(a, b)
    iszero(b) && throw(DivideError())
    R = parent(a)
    q, r, _ = lift(Module(R, vector(R, b)), Module(R, vector(R, a)),
                   false, false, true)
    return (Array(q[1])[1], Array(r[1])[1])
end

function div(a::spoly{T}, b::spoly{T}) where T <: Nemo.FieldElem
    check_parent(a, b)
    iszero(b) && throw(DivideError())
    R = parent(a)
    q, _, _ = lift(Module(R, vector(R, b)), Module(R, vector(R, a)),
                   false, false, true)
    return Array(q[1])[1]
end

###############################################################################
#
#   GCD
#
###############################################################################

function gcd(x::spoly{T}, y::spoly{T}) where T <: Nemo.RingElem
   check_parent(x, y)
   R = parent(x)
   GC.@preserve x y R begin
      x1 = libSingular.p_Copy(x.ptr, R.ptr)
      y1 = libSingular.p_Copy(y.ptr, R.ptr)
      p = libSingular.singclap_gcd(x1, y1, R.ptr)
      return R(p)
   end
end

function gcdx(x::spoly{T}, y::spoly{T}) where T <: Nemo.FieldElem
   check_parent(x, y)
   R = parent(x)
   GC.@preserve x y R begin
      x1 = libSingular.p_Copy(x.ptr, R.ptr)
      y1 = libSingular.p_Copy(y.ptr, R.ptr)
      p, s, t = libSingular.singclap_extgcd(x1, y1, R.ptr)
      res = (R(p), R(s), R(t))
      libSingular.check_error()
      return res
   end
end

function lcm(x::spoly{T}, y::spoly{T}) where T <: Nemo.RingElem
   if iszero(x) && iszero(y)
      return parent(x)()
   end
   return divexact(x*y, gcd(x, y))
end

@doc raw"""
    primpart(x::SPolyUnion)

Return the primitive part of the polynomial, i.e. the polynomial divided by the GCD
of its coefficients.
"""
function primpart(x::SPolyUnion)
   R = parent(x)
   p = deepcopy(x)
   libSingular.p_Content(p.ptr, R.ptr)
   return p
end

@doc raw"""
    content(x::SPolyUnion)

Return the content of the polynomial, i.e. the GCD of its coefficients.
"""
function content(x::SPolyUnion)
   R = base_ring(x)
   d = R()
   for c in coefficients(x)
      d = gcd(d, c)
      if isone(d)
         break
      end
   end
   return d
end

###############################################################################
#
#   Evaluation
#
###############################################################################

function evaluate(a::SPolyUnion{T}, C::Vector{T}) where T <: Nemo.RingElem
   S = parent(a)
   GC.@preserve a C S begin
      carr = [c.ptr.cpp_object for c in C]
      n = libSingular.maEvalAt(a.ptr, carr, S.ptr)
      return base_ring(a)(n)
   end
end

function evaluate(a::SPolyUnion{T}, C::Vector{U}) where {T <: Nemo.RingElem, U <: Union{Integer, Rational}}
   C2 = [base_ring(a)(c) for c in C]
   return evaluate(a, C2)
end

function evaluate(a::SPolyUnion{T}, C::Vector{n_Z}) where T <: Nemo.RingElem
   C2 = [base_ring(a)(c) for c in C]
   return evaluate(a, C2)
end

function (a::SPolyUnion{T})() where T <: RingElem
   0 != nvars(parent(a)) && error("Number of variables does not match number of values")
   return evaluate(a, T[])
end

function (a::SPolyUnion{T})(val::T, vals::T...) where T <: RingElem
   1 + length(vals) != nvars(parent(a)) && error("Number of variables does not match number of values")
   return evaluate(a, T[val, vals...])
end

function (a::SPolyUnion{T})(val::U, vals::U...) where {T <: RingElem, U <: Union{Integer, Rational, AbstractFloat}}
   1 + length(vals) != nvars(parent(a)) && error("Number of variables does not match number of values")
   return evaluate(a, U[val, vals...])
end

function (a::SPolyUnion{T})(vals::Union{Nemo.NCRingElem, Nemo.RingElement}...) where T <: RingElem
   length(vals) != nvars(parent(a)) && error("Number of variables does not match number of values")
   R = base_ring(a)
   # The best we can do here is to cache previously used powers of the values
   # being substituted, as we cannot assume anything about the relative
   # performance of powering vs multiplication. The function should not try
   # to optimise computing new powers in any way.
   # Note that this function accepts values in a non-commutative ring, so operations
   # must be done in a certain order.
   powers = [Dict{Int, Any}() for i in 1:length(vals)]
   # First work out types of products
   r = R()
   c = zero(R)
   U = Vector{Any}(undef, length(vals))
   for j = 1:length(vals)
      W = typeof(vals[j])
      if ((W <: Integer && W != BigInt) ||
          (W <: Rational && W != Rational{BigInt}))
         c = c*zero(W)
         U[j] = parent(c)
      else
         U[j] = parent(vals[j])
         c = c*zero(parent(vals[j]))
      end
   end
   cvzip = zip(coefficients(a), exponent_vectors(a))
   for (c, v) in cvzip
      t = c
      for j = 1:length(vals)
         exp = v[j]
         if !haskey(powers[j], exp)
            powers[j][exp] = (U[j](vals[j]))^exp
         end
         t = t*powers[j][exp]
      end
      r += t
   end
   return r
end

###############################################################################
#
#   Factorization
#
###############################################################################

function factor_squarefree(x::spoly)
  R = parent(x)
  br_type = typeof(base_ring(R))
  bool = (br_type <: Union{Rationals, N_FField, N_ZpField, N_AlgExtField})

  # Check if base ring is valid
  if !bool
    error("Base ring not supported.")
  end

  a = Vector{Int32}()
  I = GC.@preserve x R Ideal(R, libSingular.singclap_sqrfree(x.ptr, a, R.ptr))
  D = Dict{typeof(I[1]), Int64}()
  n = ngens(I)
  if n == 1
    return Fac(I[1], D)
  else
    for i in 2:n
      push!(D, I[i] => Int64(a[i]))
    end
  end
  return Fac(I[1], D)
end

function factor(x::spoly)
  R = parent(x)
  br_type = typeof(base_ring(R))
  bool = (br_type <: Union{Rationals, Integers, N_FField, N_ZpField, N_AlgExtField})

  # Check if base ring is valid
  if !bool
    error("Base ring not supported.")
  end

  a = Vector{Int32}()
  I = GC.@preserve x R Ideal(R, libSingular.singclap_factorize(x.ptr, a, R.ptr))
  D = Dict{typeof(I[1]), Int64}()
  n = ngens(I)
  if n == 1
    return Fac(I[1], D)
  else
    for i in 2:n
      push!(D, I[i] => Int64(a[i]))
    end
  return Fac(I[1], D)
  end
end

###############################################################################
#
#   Variable substitution
#
###############################################################################

@doc raw"""
    substitute_variable(p::SPolyUnion,i::Int64,q::spoly)

Substitutes the `i`-th variable of the polynomial `p` with the polynomial `q`.
Returns a new polynomial.
"""
function substitute_variable(p::SPolyUnion, i::Int64, q::SPolyUnion)
    R = parent(p)
    check_parent(p, q)
    GC.@preserve p q R return R(libSingular.p_Subst(p.ptr,i, q.ptr, R.ptr))
end

@doc raw"""
    permute_variables(p::SPolyUnion, perm::Vector{Int64}, new_ring::PolyRing)

Permutes the indeterminates of `p` according to `perm` to the indeterminates
of the ring `new_ring`.
"""
function permute_variables(p::SPolyUnion, perm::Vector{Int64}, new_ring::PolyRing)
   old_ring = parent(p)
   old_base = base_ring(old_ring)
   new_base = base_ring(new_ring)
   GC.@preserve p new_ring old_ring old_base new_base begin
       perm_64 = [0]
       append!(perm_64,perm)
       perm_32 = convert(Vector{Int32},perm_64)
       map_ptr = libSingular.n_SetMap(old_base.ptr, new_base.ptr)
       poly_ptr = libSingular.p_PermPoly(p.ptr, perm_32, old_ring.ptr,
                                     new_ring.ptr, map_ptr, Ptr{Int32}(C_NULL))
       poly = new_ring(poly_ptr)
       return poly
   end
end

@doc raw"""
    homogenize(p::spoly{T}, v::spoly{T}) where T <: Nemo.RingElem

Multiply each monomial in p by a suitable power of the
variable `v` and return the corresponding homogeneous polynomial.
The variable `v` must have weight `1`.
"""
function homogenize(p::spoly{T}, v::spoly{T}) where T <: Nemo.RingElem
   R = parent(p)
   check_parent(p, v)
   i = var_index(v)
   GC.@preserve p v R begin
      isone(libSingular.p_WTotaldegree(v.ptr, R.ptr)) ||
               error("variable must have weight 1")
      z = R(libSingular.p_Homogen(p.ptr, i, R.ptr))
      return z
   end
end

###############################################################################
#
#   Conversion functions
#
###############################################################################

@doc raw"""
    (R::PolyRing){T <: RingElem}(p::AbstractAlgebra.Generic.MPoly{T})

Return a Singular polynomial in $R$ with the same coefficients and exponents as $p$.
"""
function (R::PolyRing)(p::AbstractAlgebra.Generic.MPoly{T}) where T <: Nemo.RingElem
   S = base_ring(R)
   B = MPolyBuildCtx(R)
   cvzip = zip(coefficients(p), exponent_vectors(p))
   for (c, v) in cvzip
      push_term!(B, S(c), v)
   end
   return finish(B)
end

@doc raw"""
    (R::AbstractAlgebra.Generic.MPolyRing{T}) where T <: Nemo.RingElem

Return an AbstractAlgebra polynomial in the ring $R$ with the same
coefficients and exponents as $p$.
"""
function (R::AbstractAlgebra.Generic.MPolyRing{T})(p::Singular.spoly{Singular.n_RingElem{T}}) where T <: Nemo.RingElem
   B = MPolyBuildCtx(R)
   cvzip = zip(coefficients(p), exponent_vectors(p))
   for (c, v) in cvzip
      GC.@preserve c push_term!(B, libSingular.julia(libSingular.cast_number_to_void(c.ptr)), v)
   end
   return finish(B)
end

function (R::AbstractAlgebra.Generic.MPolyRing{T})(p::Singular.spoly{Singular.n_FieldElem{T}}) where T <: Nemo.FieldElem
   B = MPolyBuildCtx(R)
   cvzip = zip(coefficients(p), exponent_vectors(p))
   for (c, v) in cvzip
      GC.@preserve c push_term!(B, libSingular.julia(libSingular.cast_number_to_void(c.ptr)), v)
   end
   return finish(B)
end


function (R::AbstractAlgebra.Generic.MPolyRing{T})(
   p::Singular.spoly{Singular.n_RingElem{Singular.RingElemWrapper{S, T}}}
) where {S, T <: Nemo.RingElem}

   B = MPolyBuildCtx(R)
   cvzip = zip(coefficients(p), exponent_vectors(p))
   for (c, v) in cvzip
      cc = GC.@preserve c libSingular.julia(libSingular.cast_number_to_void(c.ptr))
      push_term!(B, cc.data, v)
   end
   return finish(B)
end

function (R::AbstractAlgebra.Generic.MPolyRing{T})(
   p::Singular.spoly{Singular.n_FieldElem{Singular.FieldElemWrapper{S, T}}}
) where {S, T <: Nemo.FieldElem}

   B = MPolyBuildCtx(R)
   cvzip = zip(coefficients(p), exponent_vectors(p))
   for (c, v) in cvzip
      cc = GC.@preserve c libSingular.julia(libSingular.cast_number_to_void(c.ptr))
      push_term!(B, cc.data, v)
   end
   return finish(B)
end

###############################################################################
#
#   Differential functions
#
###############################################################################

@doc raw"""
    jet(x::spoly{T}, n::Int) where T <: Nemo.RingElem

Return the truncation of $x$ up to degree $n$.
"""
function jet(x::spoly{T}, n::Int) where T <: Nemo.RingElem
   R = parent(x)
   GC.@preserve x R return R(libSingular.p_Jet(x.ptr, Cint(n), R.ptr))
end

@doc raw"""
    derivative(x::spoly{T}, n::Int) where T <: Nemo.RingElem

Return the derivative of $x$ with respect to the variable of index $n$.
"""
function derivative(x::spoly{T}, n::Int) where T <: Nemo.RingElem
   R = parent(x)
   1 <= n <= nvars(R) || error("Variable does not exist")
   GC.@preserve x R return R(libSingular.p_Diff(x.ptr, Cint(n), R.ptr))
end

@doc raw"""
    derivative(x::spoly{T}, v::spoly{T}) where T <: Nemo.RingElem

Return the derivative of $x$ with respect to the variable $v$.
"""
function derivative(x::spoly{T}, v::spoly{T}) where T <: Nemo.RingElem
   R = parent(x)
   R == parent(v) && is_gen(v) || error("Second argument is not a variable")
   return derivative(x, var_index(v))
end

@doc raw"""
    jacobian_ideal(p::spoly{T}) where T <: Nemo.RingElem

Returns the ideal generated by all partial derivatives of $x$.
"""
function jacobian_ideal(p::spoly{T}) where T <: Nemo.RingElem
   R = parent(p)
   n = nvars(R)
   return Ideal(R, spoly{T}[derivative(p, i) for i in 1:n])
end

@doc raw"""
    jacobian_matrix(p::spoly{T}) where T <: Nemo.RingElem

Returns the column matrix $\{\frac{\partial p}{\partial x_i}\}_i$ of partial derivatives.
"""
function jacobian_matrix(p::spoly{T}) where T <: Nemo.RingElem
   R = parent(p)
   n = nvars(R)
   J = zero_matrix(R, n, 1)
   for i in 1:n
      J[i, 1] = derivative(p, i)
   end
   return J
end

@doc raw"""
    jacobian_matrix(a::Vector{spoly{T}}) where T <: Nemo.RingElem

Returns the matrix $\{\frac{\partial a_i}{\partial x_j}\}_{ij}$ of partial derivatives.
"""
function jacobian_matrix(A::Vector{spoly{T}}) where T <: Nemo.RingElem
   m = length(A)
   m == 0 && error("Array has to be non-empty.")
   R = parent(A[1])
   n = nvars(R)
   J = zero_matrix(R, m, n)
   for i in 1:n
      for j in 1:m
         J[j, i] = derivative(A[j], i)
      end
   end
   return J
end

###############################################################################
#
#   Unsafe operations
#
###############################################################################

function sort_terms!(x::SPolyUnion)
   S = parent(x)
   x.ptr = GC.@preserve x S libSingular.p_SortMerge(x.ptr, S.ptr)
   return x
end

function addeq!(x::SPolyUnion, y::SPolyUnion)
   R = parent(x)
   GC.@preserve x y R begin
       if y.ptr == C_NULL
       elseif x.ptr == C_NULL
          x.ptr = libSingular.p_Copy(y.ptr, R.ptr)
       else
          x.ptr = libSingular.p_Add_q(x.ptr, libSingular.p_Copy(y.ptr, R.ptr), R.ptr)
       end
       return x
   end
end

function mul!(c::SPolyUnion, x::SPolyUnion, y::SPolyUnion)
   R = parent(x)
   GC.@preserve c x y R begin
      ptr = libSingular.pp_Mult_qq(x.ptr, y.ptr, R.ptr)
      if c.ptr != C_NULL
         libSingular.p_Delete(c.ptr, R.ptr)
      end
      c.ptr = ptr
      return c
   end
end

function add!(c::SPolyUnion, x::SPolyUnion, y::SPolyUnion)
   R = parent(x)
   GC.@preserve c x y R begin
      x1 = libSingular.p_Copy(x.ptr, R.ptr)
      y1 = libSingular.p_Copy(y.ptr, R.ptr)
      ptr = libSingular.p_Add_q(x1, y1, R.ptr)
      if c.ptr != C_NULL
         libSingular.p_Delete(c.ptr, R.ptr)
      end
      c.ptr = ptr
      return c
   end
end

function zero!(x::SPolyUnion)
   GC.@preserve x begin
      if x.ptr != C_NULL
         libSingular.p_Delete(x.ptr, parent(x).ptr)
         x.ptr = libSingular.p_ISet(0, parent(x).ptr)
      end
   end
   return x
end

###############################################################################
#
#   Changing base ring
#
###############################################################################

function change_base_ring(C::T, p::spoly) where T <: Union{Ring, Field}
   S, = Singular.polynomial_ring(C, symbols(parent(p)), ordering = parent(p).ord)
   return change_base_ring(C, p, parent = S)
end

###############################################################################
#
#   Promote rules
#
###############################################################################

promote_rule(::Type{spoly{T}}, ::Type{spoly{T}}) where T <: Nemo.RingElem = spoly{T}

function promote_rule(::Type{spoly{T}}, ::Type{U}) where {T <: Nemo.RingElem, U <: Nemo.RingElem}
   promote_rule(T, U) == T ? spoly{T} : Union{}
end

function promote_rule(::Type{spoly{T}}, ::Type{T}) where {T <: Nemo.RingElem}
   return spoly{T}
end

function promote_rule(::Type{spoly{n_RingElem{RingElemWrapper{S, T}}}}, ::Type{U}) where {T <: Nemo.RingElem, U <: Nemo.RingElem, S}
   return spoly{n_RingElem{RingElemWrapper{S, T}}}
end

function promote_rule(::Type{spoly{n_FieldElem{FieldElemWrapper{S, T}}}}, ::Type{U}) where {T <: Nemo.FieldElem, U <: Nemo.FieldElem, S}
   return spoly{n_FieldElem{FieldElemWrapper{S, T}}}
end

###############################################################################
#
#   Build context
#
###############################################################################

function MPolyBuildCtx(R::Union{PolyRing, PluralRing, LPRing})
   # MPolyBuildCtx's inner constructor is too uncomfortable to use
   t = zero(R)
   M = Nemo.@new_struct(MPolyBuildCtx{typeof(t), typeof(t.ptr)}, t, t.ptr)
   return M
end

function push_term!(M::MPolyBuildCtx{S, U}, c::T, e::Vector{Int}) where {U,
                                        T <: Nemo.RingElem, S <: SPolyUnion{T}}
   R = parent(M.poly)
   RR = base_ring(R)
   GC.@preserve c R RR begin
      if S <: slpalg
         nv = nvars(R)*degree_bound(R)
         v = zeros(Int, nv)
         _exponent_word_to_lp_vec!(v, e, nvars(R), degree_bound(R))
      else
         nv = nvars(R)
         v = e
         nv == length(e) || error("Incorrect number of exponents in push_term!")
      end
      if iszero(c)
         return M.poly
      end
      ptr = libSingular.p_Init(R.ptr)
      libSingular.p_SetCoeff0(ptr, libSingular.n_Copy(c.ptr, RR.ptr), R.ptr)
      for i = 1:nv
         libSingular.p_SetExp(ptr, i, v[i], R.ptr)
      end
      libSingular.p_Setm(ptr, R.ptr)
      if M.state.cpp_object == C_NULL
         @assert M.poly.ptr.cpp_object == C_NULL
         M.poly.ptr = ptr
      else
         libSingular.p_SetNext(M.state, ptr)
      end
      M.state = ptr
      return M.poly
   end
end

function finish(M::MPolyBuildCtx{<:SPolyUnion{T}, U}) where {U, T <: Nemo.RingElem}
   p = sort_terms!(M.poly)
   M.poly = zero(parent(p))
   M.state = M.poly.ptr
   return p
end

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function (R::PolyRing)()
   T = elem_type(base_ring(R))
   return spoly{T}(R)
end

function (R::PolyRing)(n::Int)
   T = elem_type(base_ring(R))
   return spoly{T}(R, n)
end

function (R::PolyRing)(n::Integer)
   T = elem_type(base_ring(R))
   return spoly{T}(R, BigInt(n))
end

function (R::PolyRing)(n::n_Z)
   T = elem_type(base_ring(R))
   return spoly{T}(R, base_ring(R)(n))
end

function (R::PolyRing)(n::Rational)
   T = elem_type(base_ring(R))
   return spoly{T}(R, base_ring(R)(n))
end

# take ownership of the pointer - not for general users
function (R::PolyRing)(ptr::libSingular.poly_ptr)
   T = elem_type(base_ring(R))
   return spoly{T}(R, ptr)
end

function (R::PolyRing{T})(n::T) where T <: Nemo.RingElem
   parent(n) != base_ring(R) && error("Unable to coerce into polynomial ring")
   return spoly{T}(R, n)
end

# needed only to prevent ambiguity with Oscar/src/Rings/mpoly.jl
function (R::PolyRing{T})(n::T) where T <: n_unknown
   parent(n) != base_ring(R) && error("Unable to coerce into polynomial ring")
   return spoly{T}(R, n)
end

function (R::PolyRing)(f::T) where T <: Nemo.MPolyRingElem
  parent(f) == R && return f
  B = base_ring(R)
  g = MPolyBuildCtx(R)
  for (c, e) = zip(Nemo.coefficients(f), Nemo.exponent_vectors(f))
    push_term!(g, B(c), e)
  end
  return finish(g)
end

function (R::PolyRing{S})(n::T) where {S <: Nemo.RingElem, T <: Nemo.RingElem}
   return spoly{S}(R, base_ring(R)(n))
end

function (R::PolyRing)(p::spoly)
   parent(p) != R && error("Unable to coerce polynomial")
   return p
end

###############################################################################
#
#   polynomial_ring constructor
#
###############################################################################

# turn user input into an sordering
function get_fancy_ordering(ordering, ordering2)
   if ordering isa Symbol
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
   elseif ordering isa sordering
      fancy_ordering = ordering
   else
      error("ordering must be a Symbol or an sordering")
   end
   return fancy_ordering::sordering
end

function _PolynomialRing(R, s::AbstractVector{<:VarName}, ordering, ordering2, cached, degree_bound)
   return _PolynomialRing(R, map(Symbol, s), ordering, ordering2, cached, degree_bound)
end

function _PolynomialRing(R, s::Vector{Symbol}, ordering, ordering2, cached, degree_bound)
   sord = get_fancy_ordering(ordering, ordering2)
   z = PolyRing{elem_type(R)}(R, s, sord, cached, degree_bound)
   return (z, gens(z))
end

# keyword arguments do not participate in dispatch
function polynomial_ring(R::Union{Ring, Field}, s::AbstractVector{<:VarName};
                        ordering = :degrevlex, ordering2::Symbol = :comp1min,
                        cached::Bool = true, degree_bound::Int = 0)
   return _PolynomialRing(R, s, ordering, ordering2, cached, degree_bound)
end

function polynomial_ring(R::Nemo.Ring, s::AbstractVector{<:VarName};
                        ordering = :degrevlex, ordering2::Symbol = :comp1min,
                        cached::Bool = true, degree_bound::Int = 0)
   R = CoefficientRing(R)
   return _PolynomialRing(R, s, ordering, ordering2, cached, degree_bound)
end

macro polynomial_ring(R, s, n, o)
   S = gensym()
   y = gensym()
   v0 = [s*string(i) for i in 1:n]
   exp1 = :(($S, $y) = polynomial_ring($R, $v0; ordering=$o))
   v = [:($(Symbol(s, i)) = $y[$i]) for i in 1:n]
   v1 = Expr(:block, exp1, v..., S)
   return esc(v1)
end

macro polynomial_ring(R, s, n)
   S = gensym()
   y = gensym()
   v0 = [s*string(i) for i in 1:n]
   exp1 = :(($S, $y) = polynomial_ring($R, $v0))
   v = [:($(Symbol(s, i)) = $y[$i]) for i in 1:n]
   v1 = Expr(:block, exp1, v..., S)
   return esc(v1)
end
