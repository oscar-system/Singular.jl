export n_Q, Rationals, reconstruct, isone, iszero, isunit, divexact,
       canonical_unit

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

elem_type(::Type{Rationals}) = n_Q

parent(a::n_Q) = QQ

parent_type(::Type{n_Q}) = Rationals

base_ring(a::n_Q) = ZZ

base_ring(a::Rationals) = ZZ

function deepcopy_internal(a::n_Q, dict::IdDict)
   return parent(a)(libSingular.n_Copy(a.ptr, parent(a).ptr))
end

function hash(a::n_Q, h::UInt)
   n = numerator(a)
   d = denominator(a)
   nhash = hash(n, h)
   dhash = hash(d, h)
   return xor(xor(nhash, dhash), 0xf348b78c190e8fc1%UInt)
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

one(::Rationals) = QQ(1)

zero(::Rationals) = QQ(0)

function isone(n::n_Q)
   c = parent(n)
   return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_Q)
   c = parent(n)
   return libSingular.n_IsZero(n.ptr, c.ptr)
end

isunit(n::n_Q) = !iszero(n)

@doc Markdown.doc"""
    numerator(n::n_Q)

Return the numerator of the given fraction.
"""
function numerator(n::n_Q)
   p = libSingular.n_GetNumerator(n.ptr, parent(n).ptr)
   pp = libSingular.nApplyMapFunc(n_Q_2_n_Z, p, QQ.ptr, ZZ.ptr)
   libSingular.n_Delete(p, QQ.ptr)
   return ZZ(pp)
end

@doc Markdown.doc"""
    denominator(n::n_Q)

Return the denominator of the given fraction.
"""
function denominator(n::n_Q)
   p = libSingular.n_GetDenom(n.ptr, parent(n).ptr)
   pp = libSingular.nApplyMapFunc(n_Q_2_n_Z, p, QQ.ptr, ZZ.ptr)
   libSingular.n_Delete(p, QQ.ptr)
   return ZZ(pp)
end

@doc Markdown.doc"""
    abs(n::n_Q)

Return the absolute value of the given fraction.
"""
function abs(n::n_Q)
   if libSingular.n_GreaterZero(n.ptr, parent(n).ptr) || iszero(n)
      return deepcopy(n)
   else
      return -n
   end
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(x::n_Q) = x

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, c::Rationals)
   print(io, "Rational Field")
end

function show(io::IO, n::n_Q)
    libSingular.StringSetS("")

    libSingular.n_Write(n.ptr, parent(n).ptr, false)
 
    m = libSingular.StringEndS()
 
    print(io,m)
end

needs_parentheses(x::n_Q) = false

isnegative(x::n_Q) = !libSingular.n_GreaterZero(x.ptr, parent(x).ptr) &&
                      !iszero(x)

show_minus_one(::Type{n_Q}) = false

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(x::n_Q) 
    C = parent(x)
    ptr = libSingular.n_Neg(x.ptr, C.ptr)
    return C(ptr) 
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(x::n_Q, y::n_Q)
   c = parent(x)
   p = libSingular.n_Add(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function -(x::n_Q, y::n_Q)
   c = parent(x)
   p = libSingular.n_Sub(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function *(x::n_Q, y::n_Q)
   c = parent(x)
   p = libSingular.n_Mult(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

+(x::n_Q, y::Integer) = x + parent(x)(y)

+(x::Integer, y::n_Q) = parent(y)(x) + y

-(x::n_Q, y::Integer) = x - parent(x)(y)

-(x::Integer, y::n_Q) = parent(y)(x) - y

*(x::n_Q, y::Integer) = x*parent(x)(y)

*(x::Integer, y::n_Q) = parent(y)(x)*y

###############################################################################
#
#   Comparison
#
###############################################################################

function isless(x::n_Q, y::n_Q)
    libSingular.n_Greater(y.ptr, x.ptr, parent(x).ptr)
end

function ==(x::n_Q, y::n_Q)
    return libSingular.n_Equal(x.ptr, y.ptr, parent(x).ptr)
end

isequal(x::n_Q, y::n_Q) = (x == y)

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

==(x::n_Q, y::Integer) = (x ==  parent(x)(y))

==(x::Integer, y::n_Q) = (parent(y)(x) == y)

isequal(x::n_Q, y::Integer) = (x == y)

isequal(x::Integer, y::n_Q) = (x == y)

isless(x::n_Q, y::Integer) = isless(x, parent(x)(y))

isless(x::Integer, y::n_Q) = isless(parent(y)(x), y)

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::n_Q, y::Int)
    if y < 0
       return div(parent(x)(1), x^(-y))
    elseif isone(x)
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

function inv(x::n_Q)
   c = parent(x)
   p = libSingular.n_Invers(x.ptr, c.ptr)
   return c(p)
end

function divexact(x::n_Q, y::n_Q)
   c = parent(x)
   p = libSingular.n_Div(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   GCD and LCM
#
###############################################################################

function gcd(x::n_Q, y::n_Q)
   if x == 0 && y == 0
      return zero(parent(x))
   end
   par = parent(x)
   p = libSingular.n_Gcd(x.ptr, y.ptr, par.ptr)
   return par(p)
end

function lcm(x::n_Q, y::n_Q)
   if x == 0 && y == 0
      return zero(parent(x))
   end
   return divexact(x*y, gcd(x, y))
end

###############################################################################
#
#   Rational reconstruction
#
###############################################################################

@doc Markdown.doc"""
    reconstruct(x::n_Z, y::n_Z)

Given $x$ modulo $y$, find $r/s$ such that $x \equiv r/s \pmod{y}$ for values
$r$ and $s$ satisfying the bound $y > 2(|r| + 1)(s + 1)$.
"""
function reconstruct(x::n_Z, y::n_Z)
   p = libSingular.n_Farey(x.ptr, y.ptr, QQ.ptr)
   return QQ(p)
end

reconstruct(x::n_Z, y::Integer) = reconstruct(x, ZZ(y))

reconstruct(x::Integer, y::n_Z) = reconstruct(ZZ(x), y)

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addeq!(x::n_Q, y::n_Q)
   x.ptr = libSingular.n_InpAdd(x.ptr, y.ptr, parent(x).ptr)
   return x
end

function mul!(x::n_Q, y::n_Q, z::n_Q)
   ptr = libSingular.n_Mult(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   return x
end

function add!(x::n_Q, y::n_Q, z::n_Q)
   ptr = libSingular.n_Add(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   return x
end

function zero!(x::n_Q)
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

promote_rule(C::Type{n_Q}, ::Type{T}) where {T <: Integer} = n_Q

promote_rule(C::Type{n_Q}, ::Type{Nemo.fmpz}) = n_Q

promote_rule(C::Type{n_Q}, ::Type{n_Q}) = n_Z

Rational{T}(x::n_Q) where {T} = convert(T, numerator(x)) // convert(T, denominator(x))
Rational(x::n_Q) = Integer(numerator(x)) // Integer(denominator(x))

convert(::Type{T}, x::n_Q) where {T <: Rational} = T(x)

Nemo.fmpq(x::n_Q) = fmpq(fmpz(numerator(x)), fmpz(denominator(x)))
convert(::Type{fmpq}, x::n_Q) = fmpq(x)
(::FlintRationalField)(x::n_Q) = fmpq(x)

###############################################################################
#
#   Parent call functions
#
###############################################################################

(::Rationals)() = n_Q()

(R::Rationals)(n::Union{Integer,fmpz,n_Z}, m::Union{Integer,fmpz,n_Z}) = R(n) // R(m)

(::Rationals)(n::libSingular.number_ptr) = n_Q(n)

# from integers

(::Rationals)(n::Int) = n_Q(n)

(R::Rationals)(x::Integer) = R(libSingular.n_InitMPZ(BigInt(x), R.ptr))

(::Rationals)(n::n_Z) = n_Q(n)

(R::Rationals)(x::fmpz) = convert_from_fmpz(R, x)

# from rationals

(::Rationals)(x::Rational{Int}) = n_Q(numerator(x)) // n_Q(denominator(x))

(R::Rationals)(x::Union{Rational,fmpq}) = R(numerator(x)) // R(denominator(x))

(::Rationals)(n::n_Q) = n
