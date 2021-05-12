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

characteristic(R::Rationals) = ZZ(0)

function deepcopy_internal(a::n_Q, dict::IdDict)
   c = parent(a)
   GC.@preserve a c return c(libSingular.n_Copy(a.ptr, c.ptr))
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
   GC.@preserve n c return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_Q)
   c = parent(n)
   GC.@preserve n c return libSingular.n_IsZero(n.ptr, c.ptr)
end

isunit(n::n_Q) = !iszero(n)

@doc Markdown.doc"""
    numerator(x::n_Q)

Return the numerator of the given fraction.
"""
function numerator(x::n_Q)
   c = parent(x)
   GC.@preserve x c QQ ZZ begin
      xref = Ref(x.ptr)
      p = libSingular.n_GetNumerator(xref, c.ptr)
      x.ptr = xref[]
      pp = libSingular.nApplyMapFunc(n_Q_2_n_Z, p, QQ.ptr, ZZ.ptr)
      libSingular.n_Delete(p, QQ.ptr)
      return ZZ(pp)
   end
end

@doc Markdown.doc"""
    denominator(x::n_Q)

Return the denominator of the given fraction.
"""
function denominator(x::n_Q)
   c = parent(x)
   GC.@preserve x c QQ ZZ begin
      xref = Ref(x.ptr)
      p = libSingular.n_GetDenom(xref, c.ptr)
      x.ptr = xref[]
      pp = libSingular.nApplyMapFunc(n_Q_2_n_Z, p, QQ.ptr, ZZ.ptr)
      libSingular.n_Delete(p, QQ.ptr)
      return ZZ(pp)
   end
end

@doc Markdown.doc"""
    abs(n::n_Q)

Return the absolute value of the given fraction.
"""
function abs(n::n_Q)
   c = parent(n)
   GC.@preserve n c begin
      if libSingular.n_GreaterZero(n.ptr, c.ptr) || iszero(n)
         return deepcopy(n)
      else
         return -n
      end
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

function AbstractAlgebra.expressify(n::n_Q; context = nothing)::Any
   return AbstractAlgebra.expressify(Rational{BigInt}(n), context = context)
end

function show(io::IO, n::n_Q)
   libSingular.StringSetS("")
   c = parent(n)
   GC.@preserve n c libSingular.n_Write(n.ptr, c.ptr, false)
   print(io, libSingular.StringEndS())
end

isnegative(x::n_Q) = !libSingular.n_GreaterZero(x.ptr, parent(x).ptr) &&
                      !iszero(x)

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(x::n_Q)
   c = parent(x)
   ptr = GC.@preserve x c libSingular.n_Neg(x.ptr, c.ptr)
   return c(ptr)
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(x::n_Q, y::n_Q)
   c = parent(x)
   GC.@preserve x y c begin
      p = libSingular.n_Add(x.ptr, y.ptr, c.ptr)
      p = libSingular.n_Normalize(p, c.ptr)
      return c(p)
   end
end

function -(x::n_Q, y::n_Q)
   c = parent(x)
   GC.@preserve x y c begin
      p = libSingular.n_Sub(x.ptr, y.ptr, c.ptr)
      p = libSingular.n_Normalize(p, c.ptr)
      return c(p)
   end
end

function *(x::n_Q, y::n_Q)
   c = parent(x)
   GC.@preserve x y c begin
      p = libSingular.n_Mult(x.ptr, y.ptr, c.ptr)
      p = libSingular.n_Normalize(p, c.ptr)
      return c(p)
   end
end

###############################################################################
#
#   Comparison
#
###############################################################################

function isless(x::n_Q, y::n_Q)
   c = parent(x)
   GC.@preserve x y c return libSingular.n_Greater(y.ptr, x.ptr, c.ptr)
end

function ==(x::n_Q, y::n_Q)
   c = parent(x)
   GC.@preserve x y c return libSingular.n_Equal(x.ptr, y.ptr, c.ptr)
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
      c = parent(x)
      p = GC.@preserve x c libSingular.n_Power(x.ptr, y, c.ptr)
      return c(p)
   end
end

###############################################################################
#
#   Exact division
#
###############################################################################

function inv(x::n_Q)
   c = parent(x)
   p = GC.@preserve x c libSingular.n_Invers(x.ptr, c.ptr)
   return c(p)
end

function divexact(x::n_Q, y::n_Q)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Div(x.ptr, y.ptr, c.ptr)
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
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Gcd(x.ptr, y.ptr, c.ptr)
   return c(p)
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
   p = GC.@preserve x y QQ libSingular.n_Farey(x.ptr, y.ptr, QQ.ptr)
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
   c = parent(x)
   GC.@preserve x y c begin
      x.ptr = libSingular.n_InpAdd(x.ptr, y.ptr, c.ptr)
      x.ptr = libSingular.n_Normalize(x.ptr, c.ptr)
      return x
   end
end

function mul!(x::n_Q, y::n_Q, z::n_Q)
   c = parent(x)
   GC.@preserve x y z c begin
      ptr = libSingular.n_Mult(y.ptr, z.ptr, c.ptr)
      ptr = libSingular.n_Normalize(ptr, c.ptr)
      libSingular.n_Delete(x.ptr, c.ptr)
      x.ptr = ptr
      return x
   end
end

function add!(x::n_Q, y::n_Q, z::n_Q)
   c = parent(x)
   GC.@preserve x y z c begin
      ptr = libSingular.n_Add(y.ptr, z.ptr, c.ptr)
      ptr = libSingular.n_Normalize(ptr, c.ptr)
      libSingular.n_Delete(x.ptr, c.ptr)
      x.ptr = ptr
      return x
   end
end

function zero!(x::n_Q)
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

promote_rule(C::Type{n_Q}, ::Type{T}) where {T <: Integer} = n_Q

promote_rule(C::Type{n_Q}, ::Type{Nemo.fmpz}) = n_Q

promote_rule(C::Type{n_Q}, ::Type{n_Z}) = n_Q

Rational{T}(x::n_Q) where {T} = convert(T, numerator(x)) // convert(T, denominator(x))
Rational(x::n_Q) = Integer(numerator(x)) // Integer(denominator(x))

convert(::Type{T}, x::n_Q) where {T <: Rational} = T(x)

Nemo.fmpq(x::n_Q) = fmpq(fmpz(numerator(x)), fmpz(denominator(x)))
convert(::Type{fmpq}, x::n_Q) = fmpq(x)
(::FlintRationalField)(x::n_Q) = fmpq(x)

###############################################################################
#
#   Random generation
#
###############################################################################

RandomExtensions.maketype(R::Rationals, _) = n_Q

# define rand(make(QQ, n:m))
function rand(rng::AbstractRNG,
              sp::SamplerTrivial{<:Make2{n_Q, Rationals,
                                         <:AbstractArray{<:Integer}}})
   R, v = sp[][1:end]
   Z = base_ring(R)
   n = rand(rng, v)
   local d
   while true
      d = rand(rng, v)
      iszero(d) || return R(n, d)
   end
end

rand(rng::AbstractRNG, R::Rationals, n) = R(rand(rng, n))

rand(R::Rationals, n) = rand(Random.GLOBAL_RNG, R, n)

###############################################################################
#
#   Parent call functions
#
###############################################################################

(::Rationals)() = n_Q()

(R::Rationals)(n::Union{Integer,fmpz,n_Z}, m::Union{Integer,fmpz,n_Z}) = R(n) // R(m)

(::Rationals)(n::libSingular.number_ptr) = n_Q(n)

# from integers

(::Rationals)(n::IntegerLikeTypes) = n_Q(n)

# from rationals

(::Rationals)(x::Rational{Int}) = n_Q(numerator(x)) // n_Q(denominator(x))

(R::Rationals)(x::Union{Rational,fmpq}) = R(numerator(x)) // R(denominator(x))

(::Rationals)(n::n_Q) = n
