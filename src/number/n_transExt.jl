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

base_ring(a::n_transExt) = Union{}

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
   return [F(libSingular.n_Param(Cint(i), F.ptr)) for i = 1:n]
end

@doc Markdown.doc"""
    characteristic(R::N_FField)

Return the characteristic of the field.
"""
function characteristic(R::N_FField)
   return ZZ(libSingular.n_GetChar(R.ptr))
end

function deepcopy_internal(a::n_transExt, dict::IdDict)
   return parent(a)(libSingular.n_Copy(a.ptr, parent(a).ptr))
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
   F = parent(n)
   return F(libSingular.n_GetNumerator(n.ptr, F.ptr))
end

@doc Markdown.doc"""
    denominator(n::n_transExt)

Return the denominator of the given fraction.
"""
function denominator(n::n_transExt)
   F = parent(n)
   return F(libSingular.n_GetDenom(n.ptr, F.ptr))
end

function isone(n::n_transExt)
   c = parent(n)
   return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_transExt)
   c = parent(n)
   return libSingular.n_IsZero(n.ptr, c.ptr)
end

@doc Markdown.doc"""
   isunit(n::n_transExt)
Return `true` if $n$ is a unit in the field, i.e. nonzero.
"""
isunit(n::n_transExt) = !iszero(n)

function _tExt_to_poly(R, cached)
   n = transcendence_degree(R)
   P, = PolynomialRing(base_ring(R), ["a$i" for i in 1:n], cached = cached)
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
   return S(Singular.libSingular.transExt_to_poly(x.ptr, R.ptr, S.ptr))
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

function show(io::IO, n::n_transExt)
   libSingular.StringSetS("")

   libSingular.n_Write(n.ptr, parent(n).ptr, false)

   m = libSingular.StringEndS()

   print(io, m)
end

needs_parentheses(x::n_transExt) = false

show_minus_one(::Type{n_transExt}) = false

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(x::n_transExt)
    C = parent(x)
    ptr = libSingular.n_Neg(x.ptr, C.ptr)
    return C(ptr)
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(x::n_transExt, y::n_transExt)
   c = parent(x)
   p = libSingular.n_Add(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function -(x::n_transExt, y::n_transExt)
   c = parent(x)
   p = libSingular.n_Sub(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function *(x::n_transExt, y::n_transExt)
   c = parent(x)
   p = libSingular.n_Mult(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

+(x::n_transExt, y::Integer) = x + parent(x)(y)

+(x::Integer, y::n_transExt) = parent(y)(x) + y

-(x::n_transExt, y::Integer) = x - parent(x)(y)

-(x::Integer, y::n_transExt) = parent(y)(x) - y

*(x::n_transExt, y::Integer) = x*parent(x)(y)

*(x::Integer, y::n_transExt) = parent(y)(x)*y

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(x::n_transExt, y::n_transExt)
    return libSingular.n_Equal(x.ptr, y.ptr, parent(x).ptr)
end


isequal(x::n_transExt, y::n_transExt) = (x == y)

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::n_transExt, y::Int)
    y < 0 && throw(DomainError())
    if isone(x)
       return x
    elseif y == 0
       return one(parent(x))
    elseif y == 1
       return x
    else
       p = libSingular.n_Power(x.ptr, y, parent(x).ptr)
       return parent(x)(p)
    end
end

###############################################################################
#
#   Exact division
#
###############################################################################

function inv(x::n_transExt)
   c = parent(x)
   p = libSingular.n_Invers(x.ptr, c.ptr)
   return c(p)
end

function divexact(x::n_transExt, y::n_transExt)
   c = parent(x)
   p = libSingular.n_Div(x.ptr, y.ptr, c.ptr)
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
   par = parent(x)
   p = libSingular.n_Gcd(x.ptr, y.ptr, par.ptr)
   return par(p)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addeq!(x::n_transExt, y::n_transExt)
   x.ptr = libSingular.n_InpAdd(x.ptr, y.ptr, parent(x).ptr)
   return x
end

function mul!(x::n_transExt, y::n_transExt, z::n_transExt)
   ptr = libSingular.n_Mult(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   return x
end

function add!(x::n_transExt, y::n_transExt, z::n_transExt)
   ptr = libSingular.n_Add(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   return x
end

function zero!(x::n_transExt)
   ptr = libSingular.n_Init(0, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   return x
end

###############################################################################
#
#   Conversions and promotions
#
###############################################################################

promote_rule(C::Type{n_transExt}, ::Type{T}) where T <: Integer = n_transExt

promote_rule(C::Type{n_transExt}, ::Type{n_Z}) = n_transExt

###############################################################################
#
#   Parent call functions
#
###############################################################################

function (R::N_FField)()
   z = n_transExt(R)
   z.parent = R
   return z
end

function (R::N_FField)(x::Integer)
   z = R(libSingular.n_InitMPZ(BigInt(x), R.ptr))
   z.parent = R
   return z
end

function (R::N_FField)(n::Int)
   z = n_transExt(R, n)
   z.parent = R
   return z
end

function (R::N_FField)(n::n_transExt)
   R != parent(n) && error("Parent ring does not match.")
   return n
end

function (R::N_FField)(n::libSingular.number_ptr)
   z = n_transExt(R, n)
   z.parent = R
   return z
end

function (R::N_FField)(x::Nemo.fmpz)
   z = convert_from_fmpz(R, x)
   z.parent = R
   return z
end

###############################################################################
#
#   FunctionField constructor
#
###############################################################################

@doc Markdown.doc"""
    FunctionField(F::Singular.Field, S::Array{String, 1})

Returns a tuple $K, a$ consisting of a function field $K$ over the field $F$
with transcendence basis stored in the array $S$.
"""
function FunctionField(F::Singular.Field, S::Array{String, 1}; cached::Bool=true)
   typeof(F) == N_GField && error("Finite field extensions of Fp not supported.")
   l = length(S)
   isempty(l) && throw(DomainError())
   R = N_FField(F, Symbol.(S))
   return tuple(R, transcendence_basis(R))
end

@doc Markdown.doc"""
    FunctionField(F::Singular.Field, n::Int)

Returns a tuple $K, a$ consisting of a function field $K$ over the field $F$
with transcendence degree $n$ and transcendence basis $a1, ..., an$.
"""
function FunctionField(F::Singular.Field, n::Int; cached::Bool=true)
   S = ["a$i" for i in 1:n]
   return FunctionField(F, S)
end
