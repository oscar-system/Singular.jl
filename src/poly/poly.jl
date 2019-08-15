export spoly, PolyRing, change_base_ring, coeff, coeffs, coeffs_expos,
       content, deflation, deflate, degree, degrees, degree_bound,
       derivative, div, divides, direm, evaluate, exponent, exponent!,
       exponent_vectors, factor, factor_squarefree, finish, gen,
       has_global_ordering, inflate, isgen,
       ismonomial, isquotientring, isterm, jacobi, jet, lc, lt, lm,
       lead_exponent, monomials, MPolyBuildCtx,
       nvars, ordering, @PolynomialRing, primpart,
       push_term!, remove, sort_terms!, symbols, terms, total_degree,
       valuation, var_index, vars

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent(p::spoly) = p.parent

base_ring(R::PolyRing{T}) where T <: Nemo.RingElem = R.base_ring

base_ring(p::spoly) = base_ring(parent(p))

elem_type(::Type{PolyRing{T}}) where T <: Nemo.RingElem = spoly{T}

parent_type(::Type{spoly{T}}) where T <: Nemo.RingElem = PolyRing{T}

@doc Markdown.doc"""
    nvars(R::PolyRing)
> Return the number of variables in the given polynomial ring.
"""
nvars(R::PolyRing) = Int(libSingular.rVar(R.ptr))

@doc Markdown.doc"""
    has_global_ordering(R::PolyRing)
> Return `true` if the given ring has a global ordering, i.e. if $1 < x$ for
> each variable $x$ in the ring. This include `:lex`, `:deglex` and `:degrevlex`
> orderings.
"""
has_global_ordering(R::PolyRing) = Bool(libSingular.rHasGlobalOrdering(R.ptr))

@doc Markdown.doc"""
    isquotientring(R::PolyRing)
> Return `true` if the given ring is the quotient of a polynomial ring with
> a non - zero ideal.
"""
isquotientring(R::PolyRing) = Bool(Singular.libSingular.rIsQuotientRing(R.ptr))

@doc Markdown.doc"""
    characteristic(R::PolyRing)
> Return the characteristic of the polynomial ring, i.e. the characteristic of the
> coefficient ring.
"""
characteristic(R::PolyRing) = Int(libSingular.rChar(R.ptr))

function gens(R::PolyRing)
   n = nvars(R)
   return [R(libSingular.rGetVar(Cint(i), R.ptr)) for i = 1:n]
end

function gen(R::PolyRing, i::Int)
   return R(libSingular.rGetVar(Cint(i), R.ptr))
end

@doc Markdown.doc"""
    symbols(R::PolyRing)
> Return symbols for the generators of the polynomial ring $R$.
"""
function symbols(R::PolyRing)
   io = IOBuffer();
   symbols = Array{Symbol,1}(undef, 0)
   for g in gens(R)
      show(io,g)
      push!(symbols, Symbol(String(take!(io))))
   end
   return symbols
end

ordering(R::PolyRing) = R.ord

@doc Markdown.doc"""
    degree_bound(R::PolyRing)
> Return the internal degree bound in each variable, enforced by Singular. This is the
> largest positive value any degree can have before an overflow will occur. This
> internal bound may be higher than the bound requested by the user via the
> `degree_bound` parameter of the `PolynomialRing` constructor.
"""
function degree_bound(R::PolyRing)
   return Int(libSingular.rBitmask(R.ptr))
end

zero(R::PolyRing) = R()

one(R::PolyRing) = R(1)

iszero(p::spoly) = p.ptr.cpp_object == C_NULL

function isone(p::spoly)
   return Bool(libSingular.p_IsOne(p.ptr, parent(p).ptr))
end

function isgen(p::spoly)
   R = parent(p)
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

function isconstant(p::spoly)
   R = parent(p)
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

function isunit(p::spoly)
   return p.ptr.cpp_object != C_NULL && libSingular.pNext(p.ptr).cpp_object == C_NULL &&
          Bool(libSingular.p_IsUnit(p.ptr, parent(p).ptr))
end

length(p::spoly) = Int(libSingular.pLength(p.ptr))

function total_degree(p::spoly)
   R = parent(p)
   libSingular.pLDeg(p.ptr, R.ptr)
end

@doc Markdown.doc"""
    lead_exponent(p::spoly)
> Return the exponent vector of the leading term of the given polynomial. The return
> value is a Julia 1-dimensional array giving the exponent for each variable of the
> leading term.
"""
function lead_exponent(p::spoly)
   R = parent(p)
   n = nvars(R)
   A = Array{Int}(undef, n)
   libSingular.p_GetExpVL(p.ptr, A, R.ptr)
   return A
end

function deepcopy_internal(p::spoly, dict::IdDict)
   p2 = libSingular.p_Copy(p.ptr, parent(p).ptr)
   return parent(p)(p2)
end

function check_parent(a::spoly{T}, b::spoly{T}) where T <: Nemo.RingElem
   parent(a) != parent(b) && error("Incompatible parent objects")
end

function canonical_unit(a::spoly{T}) where T <: Nemo.RingElem
  return a == 0 ? one(base_ring(a)) : canonical_unit(coeff(a, 0))
end

###############################################################################
#
#   Iterators
#
###############################################################################

function Base.iterate(x::Nemo.Generic.MPolyCoeffs{spoly{T}}) where T <: Nemo.RingElem
   p = x.poly
   ptr = p.ptr
   if ptr.cpp_object == C_NULL
      return nothing
   else
      R = base_ring(p)
      return R(libSingular.n_Copy(libSingular.pGetCoeff(ptr), R.ptr)), ptr
   end
end

function Base.iterate(x::Nemo.Generic.MPolyCoeffs{spoly{T}}, state) where T <: Nemo.RingElem
   state = libSingular.pNext(state)
   if state.cpp_object == C_NULL
      return nothing
   else
      R = base_ring(x.poly)
      return R(libSingular.n_Copy(libSingular.pGetCoeff(state), R.ptr)), state
   end
end

function Base.iterate(x::Nemo.Generic.MPolyExponentVectors{spoly{T}}) where T <: Nemo.RingElem
   p = x.poly
   ptr = p.ptr
   if ptr.cpp_object == C_NULL
      return nothing
   else
      R = parent(p)
      A = Array{Int}(undef, nvars(R))
      libSingular.p_GetExpVL(ptr, A, R.ptr)
      return A, ptr
   end
end

function Base.iterate(x::Nemo.Generic.MPolyExponentVectors{spoly{T}}, state) where T <: Nemo.RingElem
   state = libSingular.pNext(state)
   if state.cpp_object == C_NULL
      return nothing
   else
      R = parent(x.poly)
      A = Array{Int}(undef, nvars(R))
      libSingular.p_GetExpVL(state, A, R.ptr)
      return A, state
   end
end

function Base.iterate(x::Nemo.Generic.MPolyTerms{spoly{T}}) where T <: Nemo.RingElem
   p = x.poly
   ptr = p.ptr
   if ptr.cpp_object == C_NULL
      return nothing
   else
      R = parent(p)
      return R(libSingular.p_Head(ptr, R.ptr)), ptr
   end
end

function Base.iterate(x::Nemo.Generic.MPolyTerms{spoly{T}}, state) where T <: Nemo.RingElem
   state = libSingular.pNext(state)
   if state.cpp_object == C_NULL
      return nothing
   else
      R = parent(x.poly)
      return R(libSingular.p_Head(state, R.ptr)), state
   end
end

function Base.iterate(x::Nemo.Generic.MPolyMonomials{spoly{T}}) where T <: Nemo.RingElem
   p = x.poly
   ptr = p.ptr
   if ptr.cpp_object == C_NULL
      return nothing
   else
      S = parent(p)
      R = base_ring(p)
      mptr = libSingular.p_Head(ptr, S.ptr)
      n = libSingular.n_Copy(one(R).ptr, R.ptr)
      libSingular.p_SetCoeff0(mptr, n, S.ptr)
      return S(mptr), ptr
   end
end

function Base.iterate(x::Nemo.Generic.MPolyMonomials{spoly{T}}, state) where T <: Nemo.RingElem
   state = libSingular.pNext(state)
   if state.cpp_object == C_NULL
      return nothing
   else
      p = x.poly
      S = parent(p)
      R = base_ring(p)
      mptr = libSingular.p_Head(state, S.ptr)
      n = libSingular.n_Copy(one(R).ptr, R.ptr)
      libSingular.p_SetCoeff0(mptr, n, S.ptr)
      return S(mptr), state
   end
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::PolyRing)
   s = libSingular.rString(R.ptr)
   print(io, "Singular Polynomial Ring ", s)
end

function show(io::IO, a::spoly)
   s = libSingular.p_String(a.ptr, parent(a).ptr)
   print(io, s)
end

show_minus_one(::Type{spoly{T}}) where T <: Nemo.RingElem = show_minus_one(T)

needs_parentheses(x::spoly) = length(x) > 1

isnegative(x::spoly) = isconstant(x) && !iszero(x) && isnegative(coeff(x, 0))

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(a::spoly)
   a1 = libSingular.p_Copy(a.ptr, parent(a).ptr)
   s = libSingular.p_Neg(a1, parent(a).ptr)
   return parent(a)(s)
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function (a::spoly{T} + b::spoly{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   a1 = libSingular.p_Copy(a.ptr, parent(a).ptr)
   b1 = libSingular.p_Copy(b.ptr, parent(a).ptr)
   s = libSingular.p_Add_q(a1, b1, parent(a).ptr)
   return parent(a)(s)
end

function (a::spoly{T} - b::spoly{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   a1 = libSingular.p_Copy(a.ptr, parent(a).ptr)
   b1 = libSingular.p_Copy(b.ptr, parent(a).ptr)
   s = libSingular.p_Sub(a1, b1, parent(a).ptr)
   return parent(a)(s)
end

function (a::spoly{T} * b::spoly{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   s = libSingular.pp_Mult_qq(a.ptr, b.ptr, parent(a).ptr)
   return parent(a)(s)
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::spoly, y::Int)
    y > typemax(Cint) || y < typemin(Cint) && throw(DomainError())
    R = parent(x)
    if isone(x)
       return deepcopy(x)
    elseif y == 0
       return one(R)
    elseif y == 1
       return deepcopy(x)
    end
    x1 = libSingular.p_Copy(x.ptr, R.ptr)
    p = libSingular.p_Power(x1, Cint(y), R.ptr)
    return R(p)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function (x::spoly{T} == y::spoly{T}) where T <: Nemo.RingElem
    check_parent(x, y)
    return Bool(libSingular.p_EqualPolys(x.ptr, y.ptr, parent(x).ptr))
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::spoly, y::spoly)
   check_parent(x, y)
   R = parent(x)
   x1 = libSingular.p_Copy(x.ptr, R.ptr)
   y1 = libSingular.p_Copy(y.ptr, R.ptr)
   p = libSingular.p_Divide(x1, y1, R.ptr)
   return R(p)
end

###############################################################################
#
#   Ad hoc exact division
#
###############################################################################

function divexact(x::spoly{T}, y::T) where T <: Nemo.RingElem
   R = parent(x)
   base_ring(x) != parent(y) && error("Incompatible rings")
   x1 = libSingular.p_Copy(x.ptr, R.ptr)
   p = libSingular.p_Div_nn(x1, y.ptr, R.ptr)
   return R(p)
end

function divexact(x::spoly, y::n_Z)
   y1 = base_ring(x)(y)
   R = parent(x)
   x1 = libSingular.p_Copy(x.ptr, R.ptr)
   p = libSingular.p_Div_nn(x1, y1.ptr, R.ptr)
   return R(p)
end

function divexact(x::spoly, y::n_Q)
   y1 = base_ring(x)(y)
   R = parent(x)
   x1 = libSingular.p_Copy(x.ptr, R.ptr)
   p = libSingular.p_Div_nn(x1, y1.ptr, R.ptr)
   return R(p)
end

function divexact(x::spoly, y::Int)
   R = base_ring(x)
   S = parent(x)
   y1 = libSingular.n_Init(y, R.ptr)
   x1 = libSingular.p_Copy(x.ptr, S.ptr)
   p = libSingular.p_Div_nn(x1, y1, S.ptr)
   libSingular.n_Delete(y1, R.ptr)
   return S(p)
end

divexact(x::spoly, y::Integer) = divexact(x, base_ring(x)(y))

function divexact(x::spoly, y::Rational)
   return divexact(x, base_ring(x)(y))
end

###############################################################################
#
#   Divisibility testing
#
###############################################################################

function divides(x::spoly{T}, y::spoly{T}) where T <: Nemo.FieldElem
   check_parent(x, y)
   R = parent(x)
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

###############################################################################
#
#   Euclidean division
#
###############################################################################

function divrem(x::spoly{T}, y::spoly{T}) where T <: Nemo.FieldElem
   check_parent(x, y)
   R = parent(x)
   px = libSingular.p_Copy(x.ptr, R.ptr)
   py = libSingular.p_Copy(y.ptr, R.ptr)
   q, r = libSingular.p_DivRem(px, py, R.ptr)
   qref = libSingular.toPolyRef(q)
   rref = libSingular.toPolyRef(r)
   return R(qref), R(rref)
end

function div(x::spoly{T}, y::spoly{T}) where T <: Nemo.FieldElem
   check_parent(x, y)
   R = parent(x)
   px = libSingular.p_Copy(x.ptr, R.ptr)
   py = libSingular.p_Copy(y.ptr, R.ptr)
   q = libSingular.p_Divide(px, py, R.ptr)
   return R(q)
end

###############################################################################
#
#   GCD
#
###############################################################################

function gcd(x::spoly{T}, y::spoly{T}) where T <: Nemo.RingElem
   check_parent(x, y)
   R = parent(x)
   x1 = libSingular.p_Copy(x.ptr, R.ptr)
   y1 = libSingular.p_Copy(y.ptr, R.ptr)
   p = libSingular.singclap_gcd(x1, y1, R.ptr)
   return R(p)
end

function gcdx(x::spoly{T}, y::spoly{T}) where T <: Nemo.FieldElem
   check_parent(x, y)
   R = parent(x)
   x1 = libSingular.p_Copy(x.ptr, R.ptr)
   y1 = libSingular.p_Copy(y.ptr, R.ptr)
   s = [libSingular.p_ISet(0,R.ptr)]
   t = [libSingular.p_ISet(0,R.ptr)]
   p = [libSingular.p_ISet(0,R.ptr)]
   libSingular.p_ExtGcd(x1, y1, pointer(p), pointer(s), pointer(t), R.ptr)
   return R(p[]), R(s[]), R(t[])
end

function lcm(x::spoly{T}, y::spoly{T}) where T <: Nemo.RingElem
   if iszero(x) && iszero(y)
      return parent(x)()
   end
   return divexact(x*y, gcd(x, y))
end

@doc Markdown.doc"""
    primpart(x::spoly)
> Return the primitive part of the polynomial, i.e. the polynomial divided by the GCD
> of its coefficients.
"""
function primpart(x::spoly)
   R = parent(x)
   p = deepcopy(x)
   libSingular.p_Content(p.ptr, R.ptr)
   return p
end

@doc Markdown.doc"""
    content(x::spoly)
> Return the content of the polynomial, i.e. the GCD of its coefficients.
"""
function content(x::spoly)
   R = base_ring(x)
   d = R()
   for c in coeffs(x)
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

function evaluate(a::spoly{T}, C::Vector{T}) where T <: Nemo.RingElem
   S = parent(a)
   R = base_ring(a)
   @GC.preserve C begin
      carr = [c.ptr.cpp_object for c in C]
      n = libSingular.maEvalAt(a.ptr, carr, S.ptr)
      return base_ring(a)(n)
   end
end

function evaluate(a::spoly{T}, C::Vector{U}) where {T <: Nemo.RingElem, U <: Union{Integer, Rational}}
   C2 = [base_ring(a)(c) for c in C]
   return evaluate(a, C2)
end

function evaluate(a::spoly{T}, C::Vector{n_Z}) where T <: Nemo.RingElem
   C2 = [base_ring(a)(c) for c in C]
   return evaluate(a, C2)
end

function (a::spoly{T})(vals::T...) where T <: RingElem
   length(vals) != nvars(parent(a)) && error("Number of variables does not match number o
f values")
   return evaluate(a, [vals...])
end

function (a::spoly{T})(vals::U...) where {T <: RingElem, U <: Union{Integer, Rational,
 AbstractFloat}}
   length(vals) != nvars(parent(a)) && error("Number of variables does not match number o
f values")
   return evaluate(a, [vals...])
end

function (a::spoly{T})(vals::Union{Nemo.NCRingElem, Nemo.RingElement}...) where T <: RingElem
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
   cvzip = zip(coeffs(a), exponent_vectors(a))
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
> Returns a squarefree factorization of $x$.
"""
function factor_squarefree(x::spoly)
  R = parent(x)

  # Check if base ring is valid
  if((R.base_ring != Singular.QQ) && (typeof(R.base_ring) != N_ZpField))
    error("Base ring is not supported.")
  end

  a = Array{Int32, 1}()
  I = Ideal(R, libSingular.singclap_sqrfree(x.ptr, a, R.ptr))
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
> Returns the factorization of $x$.
"""
function factor(x::spoly)
  R = parent(x)

  # Check if base ring is valid
  if((R.base_ring != Singular.QQ) && (R.base_ring != Singular.ZZ) &&
  (typeof(R.base_ring) != N_ZpField))
    error("Base ring is not supported.")
  end

  a = Array{Int32, 1}()
  I = Ideal(R, libSingular.singclap_factorize(x.ptr, a, R.ptr))
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
    substitute_variable(p::spoly,i::Int64,q::spoly)
> Substitutes the `i`-th variable of the polynomial `p` with the polynomial `q`.
> Returns a new polynomial.
"""
function substitute_variable(p::spoly, i::Int64, q::spoly)
    R = parent(p)
    check_parent(p, q)
    return R( libSingular.p_Subst(p.ptr,i, q.ptr, R.ptr) )
end

@doc Markdown.doc"""
    permute_variables(p::spoly, perm::Array{Int64,1}, new_ring::PolyRing)
> Permutes the indeterminates of `p` according to `perm` to the indeterminates
> of the ring `new_ring`.
"""
function permute_variables(p::spoly, perm::Array{Int64,1}, new_ring::PolyRing)
    old_ring = parent(p)
    perm_64 = [0]
    append!(perm_64,perm)
    perm_32 = convert(Array{Int32,1},perm_64)
    map_ptr = libSingular.n_SetMap(base_ring(old_ring).ptr, base_ring(new_ring).ptr)
    poly_ptr = libSingular.p_PermPoly(p.ptr, perm_32, old_ring.ptr, new_ring.ptr, map_ptr)
    poly = new_ring(poly_ptr)
    return poly
end

###############################################################################
#
#   Unsafe operations
#
###############################################################################

function sort_terms!(x::spoly)
   S = parent(x)
   x.ptr = libSingular.p_SortMerge(x.ptr, S.ptr)
   return x
end

function addeq!(x::spoly, y::spoly)
    R = parent(x)
    if y.ptr == C_NULL
    elseif x.ptr == C_NULL
       x.ptr = libSingular.p_Copy(y.ptr, R.ptr)
    else
       x.ptr = libSingular.p_Add_q(x.ptr, libSingular.p_Copy(y.ptr, R.ptr), R.ptr)
    end
    return x
end

function mul!(c::spoly, x::spoly, y::spoly)
   R = parent(x)
   ptr = libSingular.pp_Mult_qq(x.ptr, y.ptr, R.ptr)
   if c.ptr != C_NULL
      libSingular.p_Delete(c.ptr, R.ptr)
   end
   c.ptr = ptr
   return c
end

function add!(c::spoly, x::spoly, y::spoly)
   R = parent(x)
   x1 = libSingular.p_Copy(x.ptr, R.ptr)
   y1 = libSingular.p_Copy(y.ptr, R.ptr)
   ptr = libSingular.p_Add_q(x1, y1, R.ptr)
   if c.ptr != C_NULL
      libSingular.p_Delete(c.ptr, R.ptr)
   end
   c.ptr = ptr
   return c
end

function zero!(x::spoly)
   if x.ptr != C_NULL
      libSingular.p_Delete(x.ptr, parent(x).ptr)
      x.ptr = libSingular.p_ISet(0, parent(x).ptr)
   end
   return x
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

###############################################################################
#
#   Build context
#
###############################################################################

function MPolyBuildCtx(R::PolyRing)
   return MPolyBuildCtx(R, R().ptr)
end

function push_term!(M::MPolyBuildCtx{spoly{S}, U}, c::S, expv::Vector{Int}) where {U, S <: Nemo.RingElem}
   if c == 0
      return
   end
   R = parent(M.poly)
   nv = nvars(R)
   nv != length(expv) && error("Incorrect number of exponents in push_term!")
   p = M.poly
   ptr = libSingular.p_Init(R.ptr)
   libSingular.p_SetCoeff0(ptr, libSingular.n_Copy(c.ptr, base_ring(R).ptr), R.ptr)
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

function finish(M::MPolyBuildCtx{spoly{S}, U}) where {S, U}
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
   z = spoly{T}(R)
   z.parent = R
   return z
end

function (R::PolyRing)(n::Int)
   T = elem_type(base_ring(R))
   z = spoly{T}(R, n)
   z.parent = R
   return z
end

function (R::PolyRing)(n::Integer)
   T = elem_type(base_ring(R))
   z = spoly{T}(R, BigInt(n))
   z.parent = R
   return z
end

function (R::PolyRing)(n::n_Z)
   n = base_ring(R)(n)
   ptr = libSingular.n_Copy(n.ptr, parent(n).ptr)
   T = elem_type(base_ring(R))
   z = spoly{T}(R, ptr)
   z.parent = R
   return z
end

function (R::PolyRing)(n::libSingular.poly)
   T = elem_type(base_ring(R))
   z = spoly{T}(R, n)
   z.parent = R
   return z
end

function (R::PolyRing{T})(n::T) where T <: Nemo.RingElem
   parent(n) != base_ring(R) && error("Unable to coerce into polynomial ring")
   z = spoly{T}(R, n.ptr)
   z.parent = R
   return z
end

function (R::PolyRing{S})(n::T) where {S <: Nemo.RingElem, T <: Nemo.RingElem}
   return R(base_ring(R)(n))
end

function (R::PolyRing)(p::spoly)
   parent(p) != R && error("Unable to coerce polynomial")
   return p
end

###############################################################################
#
#   PolynomialRing constructor
#
###############################################################################

function PolynomialRing(R::Union{Ring, Field}, s::Array{String, 1};
      cached::Bool = true, ordering::Symbol = :degrevlex,
      ordering2::Symbol = :comp1min, degree_bound::Int = 0)
   U = [Symbol(v) for v in s]
   T = elem_type(R)
   parent_obj = PolyRing{T}(R, U, ordering, cached, sym2ringorder[ordering],
         sym2ringorder[ordering2], degree_bound)
   return tuple(parent_obj, gens(parent_obj))
end

function PolynomialRing(R::Nemo.Ring, s::Array{String, 1}; cached::Bool = true,
      ordering::Symbol = :degrevlex, ordering2::Symbol = :comp1min,
      degree_bound::Int = 0)
   S = CoefficientRing(R)
   U = [Symbol(v) for v in s]
   T = elem_type(S)
   parent_obj = PolyRing{T}(S, U, ordering, cached, sym2ringorder[ordering],
         sym2ringorder[ordering2], degree_bound)
   return tuple(parent_obj, gens(parent_obj))
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

###############################################################################
#
#   Conversion functions
#
###############################################################################

@doc Markdown.doc"""
    AsEquivalentSingularPolynomialRing(R::AbstractAlgebra.Generic.MPolyRing{T}; cached::Bool = true,
      ordering::Symbol = :degrevlex, ordering2::Symbol = :comp1min,
      degree_bound::Int = 0)  where {T <: RingElem}
> Return a Singular (multivariate) polynomial ring over the base ring of $R$ in variables having the same names as those of R.
"""
function AsEquivalentSingularPolynomialRing(R::AbstractAlgebra.Generic.MPolyRing{T}; cached::Bool = true,
      ordering::Symbol = :degrevlex, ordering2::Symbol = :comp1min,
      degree_bound::Int = 0)  where {T <: RingElem}
   return PolynomialRing(AbstractAlgebra.Generic.base_ring(R), [string(v) for v in AbstractAlgebra.Generic.symbols(R)], cached=cached, ordering=ordering, ordering2=ordering2, degree_bound=degree_bound)
end

@doc Markdown.doc"""
    AsEquivalentAbstractAlgebraPolynomialRing(R::Singular.PolyRing{Singular.n_unknown{T}}; ordering::Symbol = :degrevlex)  where {T <: RingElem}
> Return an AbstractAlgebra (multivariate) polynomial ring over the base ring of $R$ in variables having the same names as those of R.
"""
function AsEquivalentAbstractAlgebraPolynomialRing(R::Singular.PolyRing{Singular.n_unknown{T}}; ordering::Symbol = :degrevlex)  where {T <: RingElem}
   return AbstractAlgebra.Generic.PolynomialRing(base_ring(R).base_ring, Base.convert( Array{String,1}, symbols(R)), ordering=ordering)
end


@doc Markdown.doc"""
    (R::PolyRing){T <: RingElem}(p::AbstractAlgebra.Generic.MPoly{T})
> Return a Singular polynomial in $R$ with the same coefficients and exponents as $p$.
"""
function (R::PolyRing)(p::AbstractAlgebra.Generic.MPoly{T}) where T <: Nemo.RingElem
   parent_ring = parent(p)
   n = nvars(parent_ring)
   exps = p.exps

   size_exps = size(p.exps)
   if parent_ring.ord != :lex
      exps = exps[2:size_exps[1],:]
   end
   if parent_ring.ord == :degrevlex
      exps = exps[end:-1:1,:]
   end

   new_p = zero(R)

   for i = 1:length(p)
      prod = p.coeffs[i]
      for j = 1:n
         prod = prod * gens(R)[j]^Int(exps[j,i])
      end
      new_p = new_p + prod
   end
   return(new_p)
end

@doc Markdown.doc"""
   (R::AbstractAlgebra.Generic.MPolyRing{T}){T <: Nemo.RingElem}(p::Singular.spoly{Singular.n_unknown{T}})
> Return a AbstractAlgebra polynomial in R with the same coefficients and exponents as $p$.
"""
function (R::AbstractAlgebra.Generic.MPolyRing{T})(p::Singular.spoly{Singular.n_unknown{T}}) where T <: Nemo.RingElem
   parent_ring = parent(p)
   n = nvars(parent_ring)
   new_p = zero(R)
   gens = AbstractAlgebra.Generic.gens(R)
   iterator = coeffs_expos(p)
   t = start(iterator)
   while !done(iterator,t)
      (coeff, exps),t = next(iterator,t)
      new_p = new_p + libSingular.julia(coeff.ptr) * prod(gens.^exps)
   end
   return(new_p)
end

###############################################################################
#
#   Differential functions
#
###############################################################################

@doc Markdown.doc"""
   jet(x::spoly, n::Int)
> Given a polynomial $x$ this function truncates $x$ up to degree $n$.
"""
function jet(x::spoly, n::Int)
   p = deepcopy(x)
   p.ptr = libSingular.p_Jet(x.ptr, Cint(n), parent(x).ptr)
   return p
end

@doc Markdown.doc"""
   derivative(x::spoly, n::Int)
> Given a polynomial $x$ this function returns the derivative of $x$
> with respect to the variable with number $n$.
"""
function derivative(x::spoly, n::Int)
   R = parent(x)
   if n > nvars(R) || n < 1
       error("Variable does not exist")
   else
       p = deepcopy(x)
       p.ptr = libSingular.p_Diff(p.ptr, Cint(n), R.ptr)
       return p
   end
end

@doc Markdown.doc"""
   derivative(x::spoly, v::spoly)
> Given a polynomial $x$ this function returns the derivative of $x$
> with respect to the variable $v$.
"""
function derivative(x::spoly, v::spoly)
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
   jacobi(x::spoly)
> Given a polynomial $x$ this function the Jacobian ideal of $x$.
"""
function jacobi(p::spoly)
   R = parent(p)
   n = nvars(R)
   I = Ideal(R, derivative(p, 1))
   for i in 2:n
       I = I + Ideal(R, derivative(p, i))
   end
  return I
end
