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

@doc Markdown.doc"""
    transcendence_degree(F::N_FField)

Return the transcendence degree of the given function field.
"""
transcendence_degree(F::N_FField) = Int(libSingular.rPar(F.ptr))

@doc Markdown.doc"""
    basis(F::N_FField)

Return the transcendence basis of the given function field.
"""
function transcendence_basis(F::N_FField)
   n = transcendence_degree(F)
   GC.@preserve F return [F(libSingular.n_Param(Cint(i), F.ptr)) for i = 1:n]
end

@doc Markdown.doc"""
    characteristic(R::N_FField)

Return the characteristic of the field.
"""
function characteristic(R::N_FField)
   GC.@preserve R return ZZ(libSingular.n_GetChar(R.ptr))
end

function symbols(R::N_FField)
   return R.S
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

@doc Markdown.doc"""
    numerator(n::n_transExt)

Return the numerator of the given fraction.
"""
function numerator(n::n_transExt)
   c = parent(n)
   GC.@preserve n c return c(libSingular.n_GetNumerator(n.ptr, c.ptr))
end

@doc Markdown.doc"""
    denominator(n::n_transExt)

Return the denominator of the given fraction.
"""
function denominator(n::n_transExt)
   c = parent(n)
   GC.@preserve n c return c(libSingular.n_GetDenom(n.ptr, c.ptr))
end

function isone(n::n_transExt)
   c = parent(n)
   GC.@preserve n c return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_transExt)
   c = parent(n)
   GC.@preserve n c return libSingular.n_IsZero(n.ptr, c.ptr)
end

@doc Markdown.doc"""
   isunit(n::n_transExt)
Return `true` if $n$ is a unit in the field, i.e. nonzero.
"""
isunit(n::n_transExt) = !iszero(n)

function _tExt_to_poly(R::N_FField, cached)
   n = transcendence_degree(R)
   P, = PolynomialRing(base_ring(R), map(string, symbols(R)), cached = cached)
   return P
end

@doc Markdown.doc"""
   n_transExt_to_spoly(x::n_transExt; parent::PolyRing)
Returns the numerator of $x$ as a polynomial in a polynomial
ring with at least as many variables, as the
transcendence degree of $parent(x)$. If a ring $parent_ring$ is
given to the function, it will be the parent ring of the output.
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

function AbstractAlgebra.expressify(a::n_transExt; context = nothing)::Any
   n = numerator(a)
   n = AbstractAlgebra.expressify(n_transExt_to_spoly(n), context = context)
   d = denominator(a)
   if isone(d)
      return n
   else
      d = AbstractAlgebra.expressify(n_transExt_to_spoly(d), context = context)
      return Expr(:call, ://, n, d)
   end
end

function Base.show(io::IO, a::n_transExt)
   print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

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

isequal(x::n_transExt, y::n_transExt) = (x == y)

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
   return c(p)
end

function divexact(x::n_transExt, y::n_transExt)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Div(x.ptr, y.ptr, c.ptr)
   return c(p)
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

promote_rule(C::Type{n_transExt}, ::Type{Nemo.fmpz}) = n_transExt

promote_rule(C::Type{n_transExt}, ::Type{n_Z}) = n_transExt

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

(R::N_FField)(n::libSingular.number_ptr) = n_transExt(R, n)

###############################################################################
#
#   FunctionField constructor
#
###############################################################################

@doc Markdown.doc"""
    FunctionField(F::Singular.Field, S::Vector{String})

Returns a tuple $K, a$ consisting of a function field $K$ over the field $F$
with transcendence basis stored in the array $S$.
"""
function FunctionField(F::Singular.Field, S::Vector{String}; cached::Bool=true)
   isa(F, Rationals) || isa(F, N_ZpField) ||
             error("Only transcendental extensions of Q and Fp are supported.")
   isempty(S) && throw(ArgumentError("array must be non-empty"))
   any(isempty, S) && throw(ArgumentError("strings in array must be non-empty"))
   allunique(S) || throw(ArgumentError("strings in array must be pairwise different"))
   R = N_FField(F, Symbol.(S), cached)
   return tuple(R, transcendence_basis(R))
end

@doc Markdown.doc"""
    FunctionField(F::Singular.Field, n::Int)

Returns a tuple $K, a$ consisting of a function field $K$ over the field $F$
with transcendence degree $n$ and transcendence basis $a1, ..., an$.
"""
function FunctionField(F::Singular.Field, n::Int; cached::Bool=true)
   S = ["a$i" for i in 1:n]
   return FunctionField(F, S, cached)
end
