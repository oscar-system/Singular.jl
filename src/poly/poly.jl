export spoly, PolyRing, change_base_ring, coeff, coefficients,
       constant_coefficient, content, deflation, deflate,
       degree, degrees, degree_bound,
       derivative, div, divides, evaluate, exponent,
       exponent_vectors, factor, factor_squarefree, finish, gen,
       has_global_ordering, has_mixed_ordering, has_local_ordering,
       inflate, isgen,
       ismonomial, isordering_symbolic, isquotient_ring, isterm,
       jacobian_ideal, jacobian_matrix, jet,
       leading_coefficient, leading_exponent_vector, leading_term,
       leading_monomial, lead_exponent,
       monomials, MPolyBuildCtx,
       nvars, order, ordering, ordering_as_symbol, ordering_size, ordering_weights,
       @PolynomialRing, primpart, push_term!,
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

parent(p::polyalg) = p.parent

base_ring(R::PolyRing{T}) where T <: Nemo.RingElem = R.base_ring::parent_type(T)

base_ring(p::spoly) = base_ring(parent(p))

elem_type(::Type{PolyRing{T}}) where T <: Nemo.RingElem = spoly{T}

parent_type(::Type{spoly{T}}) where T <: Nemo.RingElem = PolyRing{T}

@doc Markdown.doc"""
    nvars(R::PolyRing)

Return the number of variables in the given polynomial ring.
"""
function nvars(R::PolyRing)
   GC.@preserve R return Int(libSingular.rVar(R.ptr))
end

@doc Markdown.doc"""
    has_global_ordering(R::PolyRing)

Return `true` if the given ring has a global ordering, i.e. if $1 < x$ for
each variable $x$ in the ring. This include `:lex`, `:deglex` and `:degrevlex`
orderings.
"""
function has_global_ordering(R::PolyRing)
   GC.@preserve R return Bool(libSingular.rHasGlobalOrdering(R.ptr))
end

@doc Markdown.doc"""
    has_mixed_ordering(R::PolyRing)

Return `true` if the given ring has a mixed ordering, i.e. if $1 < x_i$ for
a variable $x_i$ and $1>x_j$ for another variable $x_j$.
"""
function has_mixed_ordering(R::PolyRing)
   GC.@preserve R return Bool(libSingular.rHasMixedOrdering(R.ptr))
end

@doc Markdown.doc"""
    has_local_ordering(R::PolyRing)

Return `true` if the given ring has a local ordering, i.e. if $1 > x$ for
all variables $x$.
"""
function has_local_ordering(R::PolyRing)
   return !has_global_ordering(R) && !has_mixed_ordering(R)
end

@doc Markdown.doc"""
    isquotient_ring(R::PolyRing)

Return `true` if the given ring is the quotient of a polynomial ring with
a non - zero ideal.
"""
function isquotient_ring(R::PolyRing)
   GC.@preserve R return Bool(Singular.libSingular.rIsQuotientRing(R.ptr))
end

function characteristic(R::PolyRing)
   GC.@preserve R return Int(libSingular.rChar(R.ptr))
end

function gens(R::PolyRing)
   n = nvars(R)
   GC.@preserve R return [R(libSingular.rGetVar(Cint(i), R.ptr)) for i = 1:n]
end

function gen(R::PolyRing, i::Int)
   GC.@preserve R return R(libSingular.rGetVar(Cint(i), R.ptr))
end

@doc Markdown.doc"""
    symbols(R::PolyRing)

Return symbols for the generators of the polynomial ring $R$.
"""
function symbols(R::PolyRing)
   GC.@preserve R return [Symbol(libSingular.rRingVar(Cshort(i - 1), R.ptr)) for i in 1:nvars(R)]
end

ordering(R::PolyRing) = R.ord

@doc Markdown.doc"""
    isordering_symbolic(R::PolyRing)

Return `true` if the ordering of `R` can be represented as a symbol.
"""
isordering_symbolic(R::PolyRing) = isordering_symbolic(R.ord)

@doc Markdown.doc"""
    ordering_as_symbol(R::PolyRing)

Assuming the ordering of `R` can be represented as a symbol, return that symbol.
"""
ordering_as_symbol(R::PolyRing) = ordering_as_symbol(R.ord)

@doc Markdown.doc"""
    degree_bound(R::PolyRing)

Return the internal degree bound in each variable, enforced by Singular. This is the
largest positive value any degree can have before an overflow will occur. This
internal bound may be higher than the bound requested by the user via the
`degree_bound` parameter of the `PolynomialRing` constructor.
"""
function degree_bound(R::PolyRing)
   GC.@preserve R return Int(libSingular.rBitmask(R.ptr))
end

zero(R::PolyRing) = R()

one(R::PolyRing) = R(1)

iszero(p::spoly) = p.ptr.cpp_object == C_NULL

function isone(p::spoly)
   GC.@preserve p return Bool(libSingular.p_IsOne(p.ptr, parent(p).ptr))
end

function isgen(p::spoly)
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
            if n == 1
               return false
            end
            n = 1
         end
      end
      return n == 1
   end
end

function isconstant(p::spoly)
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

function isunit(p::spoly)
   GC.@preserve p begin
      return p.ptr.cpp_object != C_NULL &&
             libSingular.pNext(p.ptr).cpp_object == C_NULL &&
             Bool(libSingular.p_IsUnit(p.ptr, parent(p).ptr))
   end
end

length(p::spoly) = Int(libSingular.pLength(p.ptr))

@doc Markdown.doc"""
    total_degree(p::spoly)

Return the total degree (largest sum of exponents of any monomial) of $p$.
"""
function total_degree(p::spoly)
   R = parent(p)
   GC.@preserve p R return libSingular.pLDeg(p.ptr, R.ptr)
end

@doc Markdown.doc"""
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

@doc Markdown.doc"""
    leading_exponent_vector(p::spoly)

Return the exponent vector of the leading term of the given polynomial. The return
value is a Julia 1-dimensional array giving the exponent for each variable of the
leading term.
"""
function leading_exponent_vector(p::spoly)
   R = parent(p)
   n = nvars(R)
   A = Array{Int}(undef, n)
   GC.@preserve p R libSingular.p_GetExpVL(p.ptr, A, R.ptr)
   return A
end

@deprecate lead_exponent(p::spoly) leading_exponent_vector(p)

function leading_coefficient(p::spoly)
   R = base_ring(p)
   if p.ptr.cpp_object == C_NULL
      return zero(R)
   else
      return R(libSingular.n_Copy(libSingular.pGetCoeff(p.ptr), R.ptr))
   end
end

function deepcopy_internal(p::spoly, dict::IdDict)
   R = parent(p)
   p2 = GC.@preserve p R libSingular.p_Copy(p.ptr, R.ptr)
   return R(p2)
end

function check_parent(a::spoly{T}, b::spoly{T}) where T <: Nemo.RingElem
   parent(a) != parent(b) && error("Incompatible parent objects")
end

function canonical_unit(a::spoly{T}) where T <: Nemo.RingElem
  return iszero(a) ? one(base_ring(a)) : canonical_unit(leading_coefficient(a))
end

function Base.hash(p::spoly{T}, h::UInt) where T <: Nemo.RingElem
   v = 0x37eec82e994ab710%UInt
   v = xor(hash(collect(exponent_vectors(p)), h), v)
   for c in coefficients(p)
      v = xor(hash(c, h), v)
      v = (v << 1) | (v >> (sizeof(Int)*8 - 1))
   end
   return v
end

###############################################################################
#
#   Iterators
#
###############################################################################

function Base.iterate(x::Nemo.Generic.MPolyCoeffs{polyalg{T}}) where T <: Nemo.RingElem
   p = x.poly
   ptr = p.ptr
   if ptr.cpp_object == C_NULL
      return nothing
   else
      R = base_ring(p)
      GC.@preserve R return R(libSingular.n_Copy(libSingular.pGetCoeff(ptr), R.ptr)), ptr
   end
end

function Base.iterate(x::Nemo.Generic.MPolyCoeffs{polyalg{T}}, state) where T <: Nemo.RingElem
   state = libSingular.pNext(state)
   if state.cpp_object == C_NULL
      return nothing
   else
      R = base_ring(x.poly)
      GC.@preserve R return R(libSingular.n_Copy(libSingular.pGetCoeff(state), R.ptr)), state
   end
end

function Base.iterate(x::Nemo.Generic.MPolyExponentVectors{polyalg{T}}) where T <: Nemo.RingElem
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

function Base.iterate(x::Nemo.Generic.MPolyExponentVectors{polyalg{T}}, state) where T <: Nemo.RingElem
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

function Base.iterate(x::Nemo.Generic.MPolyTerms{polyalg{T}}) where T <: Nemo.RingElem
   p = x.poly
   ptr = p.ptr
   if ptr.cpp_object == C_NULL
      return nothing
   else
      R = parent(p)
      GC.@preserve R return R(libSingular.p_Head(ptr, R.ptr)), ptr
   end
end

function Base.iterate(x::Nemo.Generic.MPolyTerms{polyalg{T}}, state) where T <: Nemo.RingElem
   state = libSingular.pNext(state)
   if state.cpp_object == C_NULL
      return nothing
   else
      R = parent(x.poly)
      GC.@preserve R return R(libSingular.p_Head(state, R.ptr)), state
   end
end

function Base.iterate(x::Nemo.Generic.MPolyMonomials{polyalg{T}}) where T <: Nemo.RingElem
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

function Base.iterate(x::Nemo.Generic.MPolyMonomials{polyalg{T}}, state) where T <: Nemo.RingElem
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

function Base.show(io::IO, a::polyalg)
   print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

isnegative(x::spoly) = isconstant(x) && !iszero(x) && isnegative(leading_coefficient(x))

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(a::polyalg)
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

function (a::polyalg{T} + b::polyalg{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   GC.@preserve a b R begin
      a1 = libSingular.p_Copy(a.ptr, R.ptr)
      b1 = libSingular.p_Copy(b.ptr, R.ptr)
      s = libSingular.p_Add_q(a1, b1, R.ptr)
      return R(s)
   end
end

function (a::polyalg{T} - b::polyalg{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   GC.@preserve a b R begin
      a1 = libSingular.p_Copy(a.ptr, R.ptr)
      b1 = libSingular.p_Copy(b.ptr, R.ptr)
      s = libSingular.p_Sub(a1, b1, R.ptr)
      return R(s)
   end
end

function (a::polyalg{T} * b::polyalg{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   s = GC.@preserve a b R libSingular.pp_Mult_qq(a.ptr, b.ptr, R.ptr)
   return R(s)
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::polyalg, y::Int)
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
      return R(p)
   end
end

###############################################################################
#
#   Comparison
#
###############################################################################

function (x::polyalg{T} == y::polyalg{T}) where T <: Nemo.RingElem
   check_parent(x, y)
   GC.@preserve x y return Bool(libSingular.p_EqualPolys(x.ptr, y.ptr, parent(x).ptr))
end

function Base.isless(x::polyalg{T}, y::polyalg{T}) where T <: Nemo.RingElem
   check_parent(x, y)
   R = parent(x)
   GC.@preserve x y R return libSingular.p_Compare(x.ptr, y.ptr, R.ptr) < 0
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::polyalg, y::polyalg)
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

function divexact(x::polyalg{T}, y::T) where T <: Nemo.RingElem
   R = parent(x)
   base_ring(x) != parent(y) && error("Incompatible rings")
   GC.@preserve x y R begin
      x1 = libSingular.p_Copy(x.ptr, R.ptr)
      p = libSingular.p_Div_nn(x1, y.ptr, R.ptr)
      return R(p)
   end
end

function divexact(x::polyalg, y::n_Z)
   y1 = base_ring(x)(y)
   R = parent(x)
   GC.@preserve x y1 R begin
      x1 = libSingular.p_Copy(x.ptr, R.ptr)
      p = libSingular.p_Div_nn(x1, y1.ptr, R.ptr)
      return R(p)
   end
end

function divexact(x::polyalg, y::n_Q)
   y1 = base_ring(x)(y)
   R = parent(x)
   GC.@preserve x y1 R begin
      x1 = libSingular.p_Copy(x.ptr, R.ptr)
      p = libSingular.p_Div_nn(x1, y1.ptr, R.ptr)
      return R(p)
   end
end

function divexact(x::polyalg, y::Int)
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

divexact(x::polyalg, y::Integer) = divexact(x, base_ring(x)(y))

function divexact(x::polyalg, y::Rational)
   return divexact(x, base_ring(x)(y))
end

###############################################################################
#
#   Divisibility testing
#
###############################################################################

function divides(x::polyalg{T}, y::polyalg{T}) where T <: Nemo.FieldElem
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

function divrem(x::polyalg{T}, y::polyalg{T}) where T <: Nemo.FieldElem
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

function div(x::polyalg{T}, y::polyalg{T}) where T <: Nemo.FieldElem
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
      s = [libSingular.p_ISet(0,R.ptr)]
      t = [libSingular.p_ISet(0,R.ptr)]
      p = [libSingular.p_ISet(0,R.ptr)]
      libSingular.p_ExtGcd(x1, y1, pointer(p), pointer(s), pointer(t), R.ptr)
      return R(p[]), R(s[]), R(t[])
   end
end

function lcm(x::spoly{T}, y::spoly{T}) where T <: Nemo.RingElem
   if iszero(x) && iszero(y)
      return parent(x)()
   end
   return divexact(x*y, gcd(x, y))
end

@doc Markdown.doc"""
    primpart(x::polyalg)

Return the primitive part of the polynomial, i.e. the polynomial divided by the GCD
of its coefficients.
"""
function primpart(x::polyalg)
   R = parent(x)
   p = deepcopy(x)
   libSingular.p_Content(p.ptr, R.ptr)
   return p
end

@doc Markdown.doc"""
    content(x::polyalg)

Return the content of the polynomial, i.e. the GCD of its coefficients.
"""
function content(x::polyalg)
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

function evaluate(a::polyalg{T}, C::Vector{T}) where T <: Nemo.RingElem
   S = parent(a)
   R = base_ring(a)
   @GC.preserve C begin
      carr = [c.ptr.cpp_object for c in C]
      n = libSingular.maEvalAt(a.ptr, carr, S.ptr)
      return base_ring(a)(n)
   end
end

function evaluate(a::polyalg{T}, C::Vector{U}) where {T <: Nemo.RingElem, U <: Union{Integer, Rational}}
   C2 = [base_ring(a)(c) for c in C]
   return evaluate(a, C2)
end

function evaluate(a::polyalg{T}, C::Vector{n_Z}) where T <: Nemo.RingElem
   C2 = [base_ring(a)(c) for c in C]
   return evaluate(a, C2)
end

function (a::polyalg{T})(vals::T...) where T <: RingElem
   length(vals) != nvars(parent(a)) && error("Number of variables does not match number o
f values")
   return evaluate(a, [vals...])
end

function (a::polyalg{T})(vals::U...) where {T <: RingElem, U <: Union{Integer, Rational,
 AbstractFloat}}
   length(vals) != nvars(parent(a)) && error("Number of variables does not match number o
f values")
   return evaluate(a, [vals...])
end

function (a::polyalg{T})(vals::Union{Nemo.NCRingElem, Nemo.RingElement}...) where T <: RingElem
   length(vals) != nvars(parent(a)) && error("Number of variables does not match number o
f values")
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
   U = Array{Any, 1}(undef, length(vals))
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

@doc Markdown.doc"""
    factor_squarefree(x::spoly)

Returns a squarefree factorization of $x$.
"""
function factor_squarefree(x::spoly)
  R = parent(x)
  br_type = typeof(base_ring(R))
  bool = (br_type <: Union{Rationals, N_FField, N_ZpField, N_AlgExtField})

  # Check if base ring is valid
  if !bool
    error("Base ring not supported.")
  end

  a = Array{Int32, 1}()
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

@doc Markdown.doc"""
    factor(x::spoly)

Returns the factorization of $x$.
"""
function factor(x::spoly)
  R = parent(x)
  br_type = typeof(base_ring(R))
  bool = (br_type <: Union{Rationals, Integers, N_FField, N_ZpField, N_AlgExtField})

  # Check if base ring is valid
  if !bool
    error("Base ring not supported.")
  end

  a = Array{Int32, 1}()
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

@doc Markdown.doc"""
    substitute_variable(p::polyalg,i::Int64,q::spoly)

Substitutes the `i`-th variable of the polynomial `p` with the polynomial `q`.
Returns a new polynomial.
"""
function substitute_variable(p::polyalg, i::Int64, q::polyalg)
    R = parent(p)
    check_parent(p, q)
    GC.@preserve p q R return R(libSingular.p_Subst(p.ptr,i, q.ptr, R.ptr))
end

@doc Markdown.doc"""
    permute_variables(p::polyalg, perm::Array{Int64,1}, new_ring::PolyRing)

Permutes the indeterminates of `p` according to `perm` to the indeterminates
of the ring `new_ring`.
"""
function permute_variables(p::polyalg, perm::Array{Int64,1}, new_ring::PolyRing)
   old_ring = parent(p)
   old_base = base_ring(old_ring)
   new_base = base_ring(new_ring)
   GC.@preserve p new_ring old_ring old_base new_base begin
       perm_64 = [0]
       append!(perm_64,perm)
       perm_32 = convert(Array{Int32,1},perm_64)
       map_ptr = libSingular.n_SetMap(old_base.ptr, new_base.ptr)
       poly_ptr = libSingular.p_PermPoly(p.ptr, perm_32, old_ring.ptr,
                                     new_ring.ptr, map_ptr, Ptr{Int32}(C_NULL))
       poly = new_ring(poly_ptr)
       return poly
   end
end

###############################################################################
#
#   Conversion functions
#
###############################################################################

@doc Markdown.doc"""
    AsEquivalentSingularPolynomialRing(R::AbstractAlgebra.Generic.MPolyRing{T}; cached::Bool = true,
      ordering::Symbol = :degrevlex, ordering2::Symbol = :comp1min,
      degree_bound::Int = 0)  where {T <: RingElem}
Return a Singular (multivariate) polynomial ring over the base ring of $R$ in variables having the same names as those of R.
"""
function AsEquivalentSingularPolynomialRing(R::AbstractAlgebra.Generic.MPolyRing{T}; cached::Bool = true,
      ordering::Symbol = :degrevlex, ordering2::Symbol = :comp1min,
      degree_bound::Int = 0)  where {T <: RingElem}
   return PolynomialRing(AbstractAlgebra.Generic.base_ring(R), [string(v) for v in AbstractAlgebra.Generic.symbols(R)], cached=cached, ordering=ordering, ordering2=ordering2, degree_bound=degree_bound)
end

@doc Markdown.doc"""
    AsEquivalentAbstractAlgebraPolynomialRing(R::Singular.PolyRing{Singular.n_unknown{T}}; ordering::Symbol = :degrevlex)  where {T <: RingElem}

Return an AbstractAlgebra (multivariate) polynomial ring over the base ring of $R$ in variables having the same names as those of R.
"""
function AsEquivalentAbstractAlgebraPolynomialRing(R::Singular.PolyRing{Singular.n_unknown{T}}; ordering::Symbol = :degrevlex)  where {T <: RingElem}
   return AbstractAlgebra.Generic.PolynomialRing(base_ring(R).base_ring,
	       [String(s) for s in symbols(R)], ordering=ordering)
end


@doc Markdown.doc"""
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

@doc Markdown.doc"""
    (R::AbstractAlgebra.Generic.MPolyRing{T}){T <: Nemo.RingElem}(p::Singular.spoly{Singular.n_unknown{T}})

Return an AbstractAlgebra polynomial in the ring $R$ with the same
coefficients and exponents as $p$.
"""
function (R::AbstractAlgebra.Generic.MPolyRing{T})(p::Singular.spoly{Singular.n_unknown{T}}) where T <: Nemo.RingElem
   B = MPolyBuildCtx(R)
   cvzip = zip(coefficients(p), exponent_vectors(p))
   for (c, v) in cvzip
      GC.@preserve c push_term!(B, libSingular.julia(libSingular.cast_number_to_void(c.ptr)), v)
   end
   return finish(B)
end

###############################################################################
#
#   Differential functions
#
###############################################################################

@doc Markdown.doc"""
   jet(x::polyalg, n::Int)
Given a polynomial $x$ this function truncates $x$ up to degree $n$.
"""
function jet(x::polyalg, n::Int)
   p = deepcopy(x)
   R = parent(x)
   p.ptr = GC.@preserve x R libSingular.p_Jet(x.ptr, Cint(n), R.ptr)
   return p
end

@doc Markdown.doc"""
   derivative(x::polyalg, n::Int)
Given a polynomial $x$ this function returns the derivative of $x$
with respect to the variable with number $n$.
"""
function derivative(x::polyalg, n::Int)
   R = parent(x)
   if n > nvars(R) || n < 1
       error("Variable does not exist")
   else
       p = deepcopy(x)
       p.ptr = GC.@preserve p R libSingular.p_Diff(p.ptr, Cint(n), R.ptr)
       return p
   end
end

@doc Markdown.doc"""
   derivative(x::polyalg, v::polyalg)
Given a polynomial $x$ this function returns the derivative of $x$
with respect to the variable $v$.
"""
function derivative(x::polyalg, v::polyalg)
   R = parent(x)
   if R == parent(v) && isgen(v)
       p = deepcopy(x)
       return(derivative(p,var_index(v)))
   else
       error("Second argument is not a variable")
       return p
   end
end

@doc Markdown.doc"""
   jacobian_ideal(x::polyalg)
Given a polynomial $x$ this function returns the Jacobian ideal of $x$.
"""
function jacobian_ideal(p::polyalg)
   R = parent(p)
   B = base_ring(R)
   n = nvars(R)
   J = Array{spoly{elem_type(B)}, 1}()
   for i in 1:n
       push!(J, derivative(p, i))
   end
  return Ideal(R, J)
end

@doc Markdown.doc"""
   jacobian_matrix(x::polyalg)
Given a polynomial $x$ this function returns the Jacobian matrix of $x$.
"""
function jacobian_matrix(p::polyalg)
   R = parent(p)
   n = nvars(R)
   J = zero_matrix(R, n, 1)
   for i in 1:n
      J[i, 1] = derivative(p, i)
   end
   return J
end

@doc Markdown.doc"""
   jacobian_matrix(A::Vector{polyalg, 1})
Given an array $A$ of polynomials over the same base ring,
this function returns the Jacobian matrix of $A$.
"""
function jacobian_matrix(A::Vector{polyalg{T}}) where T <: Nemo.RingElem
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

function sort_terms!(x::polyalg)
   S = parent(x)
   x.ptr = GC.@preserve x S libSingular.p_SortMerge(x.ptr, S.ptr)
   return x
end

function addeq!(x::polyalg, y::polyalg)
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

function mul!(c::polyalg, x::polyalg, y::polyalg)
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

function add!(c::polyalg, x::polyalg, y::polyalg)
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

function zero!(x::polyalg)
   GC.@preserve x begin
      if x.ptr != C_NULL
         libSingular.p_Delete(x.ptr, parent(x).ptr)
         x.ptr = libSingular.p_ISet(0, parent(x).ptr)
      end
      return x
end

###############################################################################
#
#   Changing base ring
#
###############################################################################

@doc Markdown.doc"""
   change_base_ring(C::T, p::spoly) where T <: Union{Ring, Field}
   > Return a polynomial ring, whose coefficient ring is subsituted by the ring
   > $C$.
   """
function change_base_ring(C::T, p::spoly) where T <: Union{Ring, Field}
   S, = Singular.PolynomialRing(C, [String(v) for v in symbols(parent(p))],
				             ordering = parent(p).ord)
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

###############################################################################
#
#   Build context
#
###############################################################################

function MPolyBuildCtx(R::PolyRing)
   t = R()
   GC.@preserve t return MPolyBuildCtx(R, t.ptr)
end

function push_term!(M::MPolyBuildCtx{polyalg{S}, U}, c::S, expv::Vector{Int}) where {U, S <: Nemo.RingElem}
   if c == 0
      return
   end
   R = parent(M.poly)
   RR = base_ring(R)
   GC.@preserve c R begin
      nv = nvars(R)
      nv != length(expv) && error("Incorrect number of exponents in push_term!")
      p = M.poly
      ptr = libSingular.p_Init(R.ptr)
      libSingular.p_SetCoeff0(ptr, libSingular.n_Copy(c.ptr, RR.ptr), R.ptr)
      for i = 1:nv
         libSingular.p_SetExp(ptr, i, expv[i], R.ptr)
      end
      libSingular.p_Setm(ptr, R.ptr)
      if iszero(p)
         p.ptr = ptr
      else
         libSingular.p_SetNext(M.state, ptr)
      end
      M.state = ptr
      return M.poly
   end
end

function finish(M::MPolyBuildCtx{polyalg{S}, U}) where {S, U}
   p = sort_terms!(M.poly)
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
   n = base_ring(R)(n)
   ptr = GC.@preserve n libSingular.n_Copy(n.ptr, parent(n).ptr)
   T = elem_type(base_ring(R))
   return spoly{T}(R, ptr)
end

function (R::PolyRing)(n::libSingular.poly_ptr)
   T = elem_type(base_ring(R))
   return spoly{T}(R, n)
end

function (R::PolyRing{T})(n::T) where T <: Nemo.RingElem
   parent(n) != base_ring(R) && error("Unable to coerce into polynomial ring")
   GC.@preserve n return spoly{T}(R, n.ptr)
end

function (R::Singular.PolyRing{T})(n::T) where T<:Singular.n_unknown
   parent(n) != base_ring(R) && error("Unable to coerce into polynomial ring")
   GC.@preserve n return spoly{T}(R, n.ptr)
end

function (R::PolyRing)(f::T) where T <: Nemo.MPolyElem
  parent(f) == R && return f
  B = base_ring(R)
  g = MPolyBuildCtx(R)
  for (c, e) = zip(Nemo.coefficients(f), Nemo.exponent_vectors(f))
    push_term!(g, B(c), e)
  end
  return finish(g)
end

function (R::PolyRing{S})(n::T) where {S <: Nemo.RingElem, T <: Nemo.RingElem}
   return R(base_ring(R)(n))
end

function (R::PolyRing)(p::spoly)
   parent(p) != R && error("Unable to coerce polynomial")
   return p
end

function(R::PolyRing)(n::libSingular.number_ptr)
    return R.base_ring(n)
end

###############################################################################
#
#   Fancy orderings
#
###############################################################################

function AbstractAlgebra.expressify(a::sordering; context = nothing)
   prod = Expr(:call, :cdot)
   for i in a.data
      if i.order == ringorder_lp
         this = Expr(:call, :ordering_lp, i.size)
      elseif i.order == ringorder_rp
         this = Expr(:call, :ordering_rp, i.size)
      elseif i.order == ringorder_dp
         this = Expr(:call, :ordering_dp, i.size)
      elseif i.order == ringorder_Dp
         this = Expr(:call, :ordering_Dp, i.size)
      elseif i.order == ringorder_wp
         this = Expr(:call, :ordering_wp, string(i.weights))
      elseif i.order == ringorder_Wp
         this = Expr(:call, :ordering_Wp, string(i.weights))
      elseif i.order == ringorder_ls
         this = Expr(:call, :ordering_ls, i.size)
      elseif i.order == ringorder_rs
         this = Expr(:call, :ordering_rs, i.size)
      elseif i.order == ringorder_ds
         this = Expr(:call, :ordering_ds, i.size)
      elseif i.order == ringorder_Ds
         this = Expr(:call, :ordering_Ds, i.size)
      elseif i.order == ringorder_ws
         this = Expr(:call, :ordering_ws, string(i.weights))
      elseif i.order == ringorder_Ws
         this = Expr(:call, :ordering_Ws, string(i.weights))
      elseif i.order == ringorder_a
         this = Expr(:call, :ordering_a, string(i.weights))
      elseif i.order == ringorder_M         
         this = Expr(:call, :ordering_M, string(transpose(reshape(i.weights, (i.size, i.size)))))
      elseif i.order == ringorder_C
         this = Expr(:call, :ordering_C)
      elseif i.order == ringorder_c
         this = Expr(:call, :ordering_c)
      elseif i.order == ringorder_S
         this = Expr(:call, :ordering_S)
      elseif i.order == ringorder_s
         this = Expr(:call, :ordering_s, i.size)
      else
         this = Expr(:call, :ordering_unknown)
      end
      push!(prod.args, this)
   end
   return prod
end

function Base.show(io::IO, mi::MIME"text/plain", a::sordering)
   Singular.AbstractAlgebra.show_via_expressify(io, mi, a)
end

function Base.show(io::IO, a::sordering)
   Singular.AbstractAlgebra.show_via_expressify(io, a)
end

function _is_basic_ordering(t::libSingular.rRingOrder_t)
    return t == ringorder_lp || t == ringorder_ls ||
           t == ringorder_rp || t == ringorder_rs ||
           t == ringorder_dp || t == ringorder_ds ||
           t == ringorder_Dp || t == ringorder_Ds
end

function _is_weighted_ordering(t::libSingular.rRingOrder_t)
    return t == ringorder_wp || t == ringorder_ws ||
           t == ringorder_Wp || t == ringorder_Ws 
end

function _basic_ordering(t::libSingular.rRingOrder_t, size::Int)
   size >= 0 || throw(ArgumentError("block size must be nonnegative"))
   return sordering([sorder_block(t, size, Int[])])
end

function _global_weighted_ordering(t::libSingular.rRingOrder_t, v::Vector{Int})
   len = length(v)
   len > 0 || throw(ArgumentError("weight vector must be non-empty"))
   all(x->x>0, v) || throw(ArgumentError("all weights must be positive"))
   return sordering([sorder_block(t, len, v)])
end

function _local_weighted_ordering(t::libSingular.rRingOrder_t, v::Vector{Int})
   len = length(v)
   len > 0 || throw(ArgumentError("weight vector must be non-empty"))
   v[1] != 0 || throw(ArgumentError("first weight must be nonzero"))
   return sordering([sorder_block(t, len, v)])
end

@doc Markdown.doc"""
    ordering_lp(nvars::Int = 1)

Represents a block of at least `nvars` variables with the
lexicographical ordering (:lex).
"""
ordering_lp(nvars::Int = 1) = _basic_ordering(Singular.ringorder_lp, nvars)

@doc Markdown.doc"""
    ordering_rp(nvars::Int = 1)

Represents a block of at least `nvars` variables with the
reverse lexicographical ordering (:revlex).
"""
ordering_rp(nvars::Int = 1) = _basic_ordering(Singular.ringorder_rp, nvars)

@doc Markdown.doc"""
    ordering_dp(nvars::Int = 1)

Represents a block of at least `nvars` variables with the
degree reverse lexicographical ordering (:degrevlex).
"""
ordering_dp(nvars::Int = 1) = _basic_ordering(Singular.ringorder_dp, nvars)

@doc Markdown.doc"""
    ordering_Dp(nvars::Int = 1)

Represents a block of at least `nvars` variables with the
degree lexicographical ordering (:deglex).
"""
ordering_Dp(nvars::Int = 1) = _basic_ordering(Singular.ringorder_Dp, nvars)

@doc Markdown.doc"""
    ordering_wp(w::Vector{Int})

Represents a block of variables with the
weighted reverse lexicographical ordering.
The weight vector `w` is expected to consist of positive integers only.
"""
ordering_wp(w::Vector{Int}) = _global_weighted_ordering(Singular.ringorder_wp, w)

@doc Markdown.doc"""
    ordering_Wp(w::Vector{Int})

Represents a block of variables with the
weighted lexicographical ordering.
The weight vector is expected to consist of positive integers only.
"""
ordering_Wp(w::Vector{Int}) = _global_weighted_ordering(Singular.ringorder_Wp, w)

@doc Markdown.doc"""
    ordering_ls(nvars::Int = 1)

Represents a block of at least `nvars` variables with the
negative lexicographical ordering (:neglex).
"""
ordering_ls(nvars::Int = 1) = _basic_ordering(Singular.ringorder_ls, nvars)

@doc Markdown.doc"""
    ordering_rs(nvars::Int = 1)

Represents a block of at least `nvars` variables with the
negative reverse lexicographical ordering (:negrevlex).
"""
ordering_rs(nvars::Int = 1) = _basic_ordering(Singular.ringorder_rs, nvars)

@doc Markdown.doc"""
    ordering_ds(nvars::Int = 1)

Represents a block of at least `nvars` variables with the
negative degree reverse lexicographical ordering (:negdegrevlex).
"""
ordering_ds(nvars::Int = 1) = _basic_ordering(Singular.ringorder_ds, nvars)

@doc Markdown.doc"""
    ordering_Ds(nvars::Int = 1)

Represents a block of at least `nvars` variables with the
negative degree reverse lexicographical ordering (:negdeglex).
"""
ordering_Ds(nvars::Int = 1) = _basic_ordering(Singular.ringorder_Ds, nvars)

@doc Markdown.doc"""
    ordering_ws(w::Vector{Int})

Represents a block of variables with the
general weighted reverse lexicographical ordering.
The weight vector `w` is expected to have a nonzero first entry.
"""
ordering_ws(w::Vector{Int}) = _local_weighted_ordering(Singular.ringorder_ws, w)

@doc Markdown.doc"""
    ordering_Ws(w::Vector{Int})

Represents a block of variables with the
general weighted lexicographical ordering.
The weight vector `w` is expected to have a nonzero first entry.
"""
ordering_Ws(w::Vector{Int}) = _local_weighted_ordering(Singular.ringorder_Ws, w)

@doc Markdown.doc"""
    ordering_a(w::Vector{Int})

Represents an extra weight vector that may precede any monomial ordering.
An extra weight vector does not define a monomial ordering by itself: it can
only be used in combination with other orderings to insert an extra line of
weights into the ordering matrix.
"""
ordering_a(w::Vector{Int}) = sordering([sorder_block(ringorder_a, 0, w)])

@doc Markdown.doc"""
    ordering_M(m::Matrix{Int}; checked::Bool = true)

Represents a block of variables with a general matrix ordering.
The matrix `m` is expected to be invertible, and this is checked by default.
"""
function ordering_M(m::Matrix{Int}; checked::Bool = true)
   (nr, nc) = size(m)
   nr > 0 && nr == nc || throw(ArgumentError("weight matrix must be square"))
   !checked || !iszero(Nemo.det(Nemo.matrix(Nemo.ZZ, m))) || throw(ArgumentError("weight matrix must nonsingular"))
   return sordering([sorder_block(ringorder_M, nr, vec(transpose(m)))])
end

function ordering_M(m::fmpz_mat, checked::Bool = true)
   !checked || !iszero(Nemo.det(m)) || throw(ArgumentError("weight matrix must nonsingular"))
   return ordering_M(Int.(m), checked = false)
end

# C, c, and S can take a dummy int in singular, but they do nothing with it?

@doc Markdown.doc"""
    ordering_C()

Represents an ascending ordering on vector components `gen(1) < gen(2) < ...`.
All monomial block orderings preceding the component ordering have higher
precedence, and all succeeding monomial block orderings have lower precedence.
It is not necessary to specify this ordering explicitly since it appended
automatically to an ordering lacking a component specification.
"""
ordering_C(dummy::Int = 0) = _basic_ordering(Singular.ringorder_C, 0)

@doc Markdown.doc"""
    ordering_c()

Represents a decending ordering on vector components `gen(1) > gen(2) > ...`.
All monomial block orderings preceding the component ordering have higher
precedence, and all succeeding monomial block orderings have lower precedence.
"""
ordering_c(dummy::Int = 0) = _basic_ordering(Singular.ringorder_c, 0)

ordering_S(dummy::Int = 0) = _basic_ordering(Singular.ringorder_S, 0)
ordering_s(syz_comp::Int = 0) = sordering([sorder_block(Singular.ringorder_s, syz_comp, Int[])])

@doc Markdown.doc"""
    *(a::sordering, b::sordering)

Return the concatenation two orderings. Some simplification may take place,
i.e. ordering_lp(2)*ordering_lp(3) may return ordering_lp(5)
"""
function Base.:*(a::sordering, b::sordering)
   return sordering(vcat(a.data, b.data))
end

function _ispure_block(a::sordering)
   if length(a.data) == 1
      return true
   elseif length(a.data) == 2
      return a.data[2].order == ringorder_C
   else
      return false
   end
end

@doc Markdown.doc"""
    ordering_size(a::sordering)

Return the size of the block of the ordering `a`, which must be a pure block.
"""
function ordering_size(a::sordering)
   _ispure_block(a) || error("ordering must be a pure block")
   return a.data[1].size
end

@doc Markdown.doc"""
    ordering_weights(a::sordering)

Return the weights of the ordering `a`, which must be a pure block.
Note that for a block with an ordering specified by a matrix,
`ordering_as_symbol(a)` will return `:matrix` and the return of
`ordering_weights(a)` can be reshaped into a square matrix of dimension
`ordering_size(a)`.
"""
function ordering_weights(a::sordering)
   _ispure_block(a) || error("ordering must be a pure block")
   return a.data[1].weights
end

isordering_symbolic(a::sordering) = isordering_symbolic_with_symbol(a)[1]

@doc Markdown.doc"""
    ordering_as_symbol(a::sordering)

If the ordering `a` is a pure block, return a symbol representing its type.
The symbol `:unknown` is returned if `a` is not a pure block.
"""
ordering_as_symbol(a::sordering) = isordering_symbolic_with_symbol(a)[2]

function isordering_symbolic_with_symbol(a::sordering)
   _ispure_block(a) || return (false, :unknown)
   o = a.data[1].order
   if o == ringorder_lp
      return (true, :lex)
   elseif o == ringorder_rp
      return (true, :revlex)
   elseif o == ringorder_ls
      return (true, :neglex)
   elseif o == ringorder_rs
      return (true, :negrevlex)
   elseif o == ringorder_dp
      return (true, :degrevlex)
   elseif o == ringorder_Dp
      return (true, :deglex)
   elseif o == ringorder_ds
      return (true, :negdegrevlex)
   elseif o == ringorder_Ds
      return (true, :negdeglex)
   elseif o == ringorder_Ds
      return (true, :negdeglex)
   elseif o == ringorder_Ds
      return (true, :negdeglex)
   elseif o == ringorder_wp
      return (true, :weightedrevlex)
   elseif o == ringorder_Wp
      return (true, :weightedlex)
   elseif o == ringorder_ws
      return (false, :negweightedrevlex)
   elseif o == ringorder_Ws
      return (false, :negweightedlex)
   elseif o == ringorder_a
      return (false, :extraweight)
   elseif o == ringorder_M
      return (false, :matrix)
   elseif o == ringorder_c
      return (false, :comp1max)
   elseif o == ringorder_C
      return (false, :comp1min)
   else
      return (false, :unknown)
   end
end

function Base.eltype(a::sordering)
   return sordering
end

function Base.length(a::sordering)
   return length(a.data)
end

function Base.iterate(a::sordering, state = 1)
   if state > length(a.data)
      return nothing
   else
      return sordering([a.data[state]]), state + 1
   end
end

function serialize_ordering(nvars::Int, ord::sordering)
   b = Cint[length(ord.data)]
   lastvar = 0
   cC_count = 0
   for l in 1:length(ord.data)
      i = ord.data[l]
      if i.order == ringorder_c || i.order == ringorder_C
         cC_count += 1
         if cC_count == 1
            push!(b, libSingular.ringorder_to_int(i.order))
            push!(b, 0)
            push!(b, 0)
            push!(b, 0)
         else
            error("more than one ordering c/C specified")
         end
      elseif i.order == ringorder_s || i.order == ringorder_S
         push!(b, libSingular.ringorder_to_int(i.order))
         push!(b, i.size) # blk0 and
         push!(b, i.size) # blk1 set to syz_comp for ringorder_s
         push!(b, 0)
      elseif i.order == ringorder_a
         push!(b, libSingular.ringorder_to_int(i.order))
         push!(b, lastvar + 1)
         nweights = min(length(i.weights), nvars - lastvar)
         push!(b, lastvar + nweights)
         push!(b, length(i.weights))
         for j in 1:nweights
            push!(b, i.weights[j])
         end         
      else
         blksize = i.size
         if _is_weighted_ordering(i.order)
            @assert blksize > 0 && length(i.weights) == blksize
         elseif i.order == ringorder_M
            @assert blksize > 0 && length(i.weights) == blksize*blksize
         elseif _is_basic_ordering(i.order)
            @assert length(i.weights) == 0
            @assert blksize >= 0
            # consume all remaining variables when succeeded by only C,c,S,s,IS
            at_end = true
            for ll in l+1:length(ord.data)
               o = ord.data[ll].order
               if o != ringorder_C && o != ringorder_c && o != ringorder_S &&
                                      o != ringorder_s && o != ringorder_IS
                  at_end = false
                  break
               end
            end
            if at_end
               blksize = max(blksize, nvars - lastvar)
            end
         else
            error("unknown ordering $(i.order)")
         end
         push!(b, libSingular.ringorder_to_int(i.order))
         push!(b, lastvar + 1)
         lastvar += blksize
         push!(b, lastvar)
         push!(b, length(i.weights))
         for j in i.weights
            push!(b, j)
         end
      end
   end

   if nvars != lastvar
      error("mismatch of number of variables (", nvars, ") and ordering (", lastvar, ")")
   end

   # add order C if none exists
   if cC_count == 0
      b[1] += 1
      push!(b, libSingular.ringorder_to_int(ringorder_C))
      push!(b, 0)
      push!(b, 0)
      push!(b, 0)
   end

   return b
end

function deserialize_ordering(b::Vector{Cint})
   off = 0
   nblocks = b[off+=1]
   data = sorder_block[]
   for i in 1:nblocks
      o = libSingular.ringorder_from_int(b[off+=1])
      blk0 = b[off+=1]
      blk1 = b[off+=1]
      nweights = b[off+=1]
      weights = Int.(b[(off+1):(off+nweights)])
      off += nweights
      blksize = blk1 - blk0 + 1
      if _is_basic_ordering(o)
         @assert nweights == 0
         push!(data, sorder_block(o, blksize, Int[]))
      elseif _is_weighted_ordering(o) || o == ringorder_M || o == ringorder_a
         @assert nweights > 0
         @assert nweights == (o == ringorder_M ? blksize*blksize : blksize)
         push!(data, sorder_block(o, blksize, weights))
      elseif o == ringorder_C || o == ringorder_c || o == ringorder_S
         @assert nweights == 0
         push!(data, sorder_block(o, 0, Int[]))
      elseif o == ringorder_s
         push!(data, sorder_block(o, blk0, Int[]))         
      else
         error("unknown ordering $o")
      end
    end
    return sordering(data)
end


###############################################################################
#
#   PolynomialRing constructor
#
###############################################################################

function _PolynomialRing(R, s::Array{String, 1}, ordering, ordering2, cached, degree_bound)
   S = [Symbol(v) for v in s]
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
   parent_obj = PolyRing{T}(R, S, fancy_ordering, cached, degree_bound)
   return tuple(parent_obj, gens(parent_obj))
end

# keyword arguments do not participate in dispatch
function PolynomialRing(R::Union{Ring, Field}, s::Array{String, 1};
                        ordering = :degrevlex, ordering2::Symbol = :comp1min,
                        cached::Bool = true, degree_bound::Int = 0)
   return _PolynomialRing(R, s, ordering, ordering2, cached, degree_bound)
end

function PolynomialRing(R::Nemo.Ring, s::Array{String, 1}; cached::Bool = true,
      ordering = :degrevlex, ordering2::Symbol = :comp1min,
      degree_bound::Int = 0)
   R = CoefficientRing(R)
   return _PolynomialRing(R, s, ordering, ordering2, cached, degree_bound)
end

macro PolynomialRing(R, s, n, o)
   S = gensym()
   y = gensym()
   v0 = [s*string(i) for i in 1:n]
   exp1 = :(($S, $y) = PolynomialRing($R, $v0; ordering=$o))
   v = [:($(Symbol(s, i)) = $y[$i]) for i in 1:n]
   v1 = Expr(:block, exp1, v..., S)
   return esc(v1)
end

macro PolynomialRing(R, s, n)
   S = gensym()
   y = gensym()
   v0 = [s*string(i) for i in 1:n]
   exp1 = :(($S, $y) = PolynomialRing($R, $v0))
   v = [:($(Symbol(s, i)) = $y[$i]) for i in 1:n]
   v1 = Expr(:block, exp1, v..., S)
   return esc(v1)
end
