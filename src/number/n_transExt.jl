export n_transExt, N_FField, transcendence_degree, transcendence_basis,
       n_transExt_to_spoly

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

elem_type(::Type{N_FField}) = n_transExt

parent(a::n_transExt) = a.parent

parent_type(::Type{n_transExt}) = N_FField

base_ring(a::n_transExt) = base_ring(parent(a))

base_ring(a::N_FField) = a.base_ring

@doc raw"""
    transcendence_degree(F::N_FField)

Return the transcendence degree of the given function field.
"""
transcendence_degree(F::N_FField) = Int(libSingular.rPar(F.ptr))

@doc raw"""
    basis(F::N_FField)

Return the transcendence basis of the given function field.
"""
function transcendence_basis(F::N_FField)
   n = transcendence_degree(F)
   GC.@preserve F return [F(libSingular.n_Param(Cint(i), F.ptr)) for i = 1:n]
end

function characteristic(R::N_FField)
   GC.@preserve R return ZZ(libSingular.n_GetChar(R.ptr))
end

function symbols(R::N_FField)
   return R.S
end

function singular_symbols(R::N_FField)
   n = transcendence_degree(R)
   GC.@preserve R return [Symbol(libSingular.n_ParameterName(Cint(i - 1), R.ptr)) for i = 1:n]
end

function deepcopy_internal(a::n_transExt, dict::IdDict)
   c = parent(a)
   GC.@preserve a c return c(libSingular.n_Copy(a.ptr, c.ptr))
end

function hash(a::n_transExt, h::UInt)
   n = n_transExt_to_spoly(numerator(a))
   d = n_transExt_to_spoly(denominator(a))
   nhash = hash(n, h)
   dhash = hash(d, h)
   return xor(xor(nhash, dhash), 0xf348b78c190e8fc1%UInt)
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

one(R::N_FField) = R(1)

zero(R::N_FField) = R(0)

@doc raw"""
    numerator(x::n_transExt)

Return the numerator of $x$.
"""
function numerator(x::n_transExt)
   c = parent(x)
   GC.@preserve x c begin
      xref = Ref(x.ptr)
      p = libSingular.n_GetNumerator(xref, c.ptr)
      x.ptr = xref[]
      return c(p)
   end
end

@doc raw"""
    denominator(x::n_transExt)

Return the denominator of $x$.
"""
function denominator(x::n_transExt)
   c = parent(x)
   GC.@preserve x c begin
      xref = Ref(x.ptr)
      p = libSingular.n_GetDenom(xref, c.ptr)
      x.ptr = xref[]
      return c(p)
   end
end

function isone(n::n_transExt)
   c = parent(n)
   GC.@preserve n c return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_transExt)
   c = parent(n)
   GC.@preserve n c return libSingular.n_IsZero(n.ptr, c.ptr)
end

is_unit(n::n_transExt) = !iszero(n)

function _tExt_to_poly(R::N_FField, cached)
   n = transcendence_degree(R)
   P, = polynomial_ring(base_ring(R), map(string, symbols(R)), cached = cached)
   return P
end

@doc raw"""
    n_transExt_to_spoly(x::n_transExt; parent::PolyRing)

Return the numerator of `x` as a polynomial in a polynomial ring with at least
as many variables as the transcendence degree of `parent(x)`. If a ring `parent`
is given to the function, it will be the parent ring of the output.
"""
function n_transExt_to_spoly(x::n_transExt; cached = true,
     parent_ring::AbstractAlgebra.MPolyRing = _tExt_to_poly(parent(x), cached))
   R = parent(x)
   B = base_ring(R)
   n = transcendence_degree(R)
   S = parent_ring

   if B != base_ring(S) || n > nvars(S)
      error("Base rings do not match.")
   end

   GC.@preserve x R S return S(libSingular.transExt_to_poly(x.ptr, R.ptr, S.ptr))
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(x::n_transExt) = x

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, F::N_FField)
   print(io, "Function Field over ", base_ring(F), " with transcendence basis ", transcendence_basis(F))
end

function expressify(a::n_transExt; context = nothing)::Any
   n = numerator(a)
   n = expressify(n_transExt_to_spoly(n), context = context)
   d = denominator(a)
   if isone(d)
      return n
   else
      d = expressify(n_transExt_to_spoly(d), context = context)
      return Expr(:call, ://, n, d)
   end
end

AbstractAlgebra.@enable_all_show_via_expressify n_transExt

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(x::n_transExt)
   c = parent(x)
   p = GC.@preserve x c libSingular.n_Neg(x.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(x::n_transExt, y::n_transExt)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Add(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function -(x::n_transExt, y::n_transExt)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Sub(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function *(x::n_transExt, y::n_transExt)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Mult(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(x::n_transExt, y::n_transExt)
   c = parent(x)
   GC.@preserve x y c return c == parent(y) && libSingular.n_Equal(x.ptr, y.ptr, c.ptr)
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::n_transExt, y::Int)
   y < 0 && throw(DomainError(y, "exponent must be non-negative"))
   if isone(x)
      return x
   elseif y == 0
      return one(parent(x))
   elseif y == 1
      return x
   else
      c = parent(x)
      p = GC.@preserve x y libSingular.n_Power(x.ptr, y, c.ptr)
      return c(p)
   end
end

###############################################################################
#
#   Exact division
#
###############################################################################

function inv(x::n_transExt)
   c = parent(x)
   p = GC.@preserve x c libSingular.n_Invers(x.ptr, c.ptr)
   z = c(p)
   libSingular.check_error()
   return z
end

function divexact(x::n_transExt, y::n_transExt; check::Bool=true)
   # low priority FIXME Singular still crashes upon division by zero here
   is_unit(y) || error("division by zero")
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Div(x.ptr, y.ptr, c.ptr)
   z = c(p)
   libSingular.check_error()
   return z
end

###############################################################################
#
#   GCD and LCM
#
###############################################################################

function gcd(x::n_transExt, y::n_transExt)
   if x == 0 && y == 0
      return zero(parent(x))
   end
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Gcd(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addeq!(x::n_transExt, y::n_transExt)
   c = parent(x)
   x.ptr = GC.@preserve x y c libSingular.n_InpAdd(x.ptr, y.ptr, c.ptr)
   return x
end

function mul!(x::n_transExt, y::n_transExt, z::n_transExt)
   c = parent(x)
   GC.@preserve x y z c begin
      ptr = libSingular.n_Mult(y.ptr, z.ptr, c.ptr)
      libSingular.n_Delete(x.ptr, c.ptr)
      x.ptr = ptr
      return x
   end
end

function add!(x::n_transExt, y::n_transExt, z::n_transExt)
   c = parent(x)
   GC.@preserve x y z c begin
      ptr = libSingular.n_Add(y.ptr, z.ptr, c.ptr)
      libSingular.n_Delete(x.ptr, c.ptr)
      x.ptr = ptr
      return x
   end
end

function zero!(x::n_transExt)
   c = parent(x)
   GC.@preserve x c begin
      ptr = libSingular.n_Init(0, c.ptr)
      libSingular.n_Delete(x.ptr, c.ptr)
      x.ptr = ptr
      return x
   end
end

###############################################################################
#
#   Conversions and promotions
#
###############################################################################

promote_rule(C::Type{n_transExt}, ::Type{T}) where T <: Integer = n_transExt

promote_rule(C::Type{n_transExt}, ::Type{T}) where T <: Rational = n_transExt

promote_rule(C::Type{n_transExt}, ::Type{Nemo.ZZRingElem}) = n_transExt

promote_rule(C::Type{n_transExt}, ::Type{Nemo.QQFieldElem}) = n_transExt

promote_rule(C::Type{n_transExt}, ::Type{n_Z}) = n_transExt

promote_rule(C::Type{n_transExt}, ::Type{n_Q}) = n_transExt

###############################################################################
#
#   Parent call functions
#
###############################################################################

(R::N_FField)(n::IntegerLikeTypes = 0) = n_transExt(R, n)

function (R::N_FField)(n::n_transExt)
   R != parent(n) && error("Parent ring does not match.")
   return n
end

function (F::Singular.N_FField)(x::fpFieldElem)
   if characteristic(F) != characteristic(parent(x))
      throw(ArgumentError("wrong characteristic"))
   end
   return n_transExt(F, x.data)
end

function (F::Singular.N_FField)(x::FqFieldElem)
   if characteristic(F) != characteristic(parent(x))
      throw(ArgumentError("wrong characteristic"))
   end
   return n_transExt(F, lift(Nemo.ZZ, x))
end

function (F::Singular.N_FField)(x::Union{n_Q, Nemo.QQFieldElem, Rational})
   return F(numerator(x)) // F(denominator(x))
end

# take ownership of the pointer - not for general users
(R::N_FField)(n::libSingular.number_ptr) = n_transExt(R, n)

###############################################################################
#
#   FunctionField constructor
#
###############################################################################

@doc raw"""
    FunctionField(F::Singular.Field, S::AbstractVector{<:VarName})

Return a tuple $K, a$ consisting of a function field $K$ over the field $F$
with transcendence basis stored in the array $S$.
"""
function FunctionField(F::Singular.Field, S::AbstractVector{<:VarName}; cached::Bool=true)
   isa(F, Rationals) || isa(F, N_ZpField) ||
             error("Only transcendental extensions of Q and Fp are supported.")
   isempty(S) && throw(ArgumentError("array must be non-empty"))
   R = N_FField(F, S, cached)
   return tuple(R, transcendence_basis(R))
end

@doc raw"""
    FunctionField(F::Singular.Field, n::Int)

Return a tuple $K, a$ consisting of a function field $K$ over the field $F$
with transcendence degree $n$ and transcendence basis $a1, ..., an$.
"""
function FunctionField(F::Singular.Field, n::Int; cached::Bool=true)
   S = [Symbol(:a, i) for i in 1:n]
   return FunctionField(F, S, cached = cached)
end
