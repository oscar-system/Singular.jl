export n_Z, Integers, crt

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

elem_type(::Type{Integers}) = n_Z

parent(a::n_Z) = ZZ

parent_type(::Type{n_Z}) = Integers

base_ring(a::n_Z) = Union{}

base_ring(a::Integers) = Union{}

characteristic(R::Integers) = ZZ(0)

function deepcopy_internal(a::n_Z, dict::IdDict)
   c = parent(a)
   GC.@preserve a c return c(libSingular.n_Copy(a.ptr, c.ptr))
end

function hash(a::n_Z, h::UInt)
   return xor(hash(convert(BigInt, a), h), 0x8ab976d24b6a0cca%UInt)
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

one(::Integers) = ZZ(1)

zero(::Integers) = ZZ(0)

function isone(n::n_Z)
   c = parent(n)
   GC.@preserve n c return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_Z)
   c = parent(n)
   GC.@preserve n c return libSingular.n_IsZero(n.ptr, c.ptr)
end

@doc Markdown.doc"""
    isunit(n::n_Z)

Return `true` if $n$ is $\pm 1$.
"""
isunit(n::n_Z) = n == 1 || n == -1

@doc Markdown.doc"""
    numerator(n::n_Z)

Return the numerator of $n$ (which is $n$ itself).
"""
function numerator(n::n_Z)
   c = parent(n)
   ptr = GC.@preserve n c libSingular.n_Copy(n.ptr, c.ptr)
   return c(ptr)
end

@doc Markdown.doc"""
    denominator(n::n_Z)

Return the denominator of $n$ (which will always be $1$).
"""
function denominator(n::n_Z)
   return one(parent(n))
end

@doc Markdown.doc"""
    abs(n::n_Z)

Return the absolute value of $n$.
"""
function abs(n::n_Z)
   c = parent(n)
   GC.@preserve n c begin
      if libSingular.n_GreaterZero(n.ptr, parent(n).ptr) || iszero(n)
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

canonical_unit(x::n_Z) = isnegative(x) ? -one(parent(x)) : one(parent(x))

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, c::Integers)
   print(io, "Integer Ring")
end

function AbstractAlgebra.expressify(n::n_Z; context = nothing)::Any
  return AbstractAlgebra.expressify(BigInt(n), context = context)
end

function show(io::IO, n::n_Z)
   libSingular.StringSetS("")
   c = parent(n)
   GC.@preserve n libSingular.n_Write(n.ptr, c.ptr, false)
   print(io, libSingular.StringEndS())
end

isnegative(x::n_Z) = !libSingular.n_GreaterZero(x.ptr, parent(x).ptr) &&
                      !iszero(x)

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(x::n_Z)
   c = parent(x)
   ptr = GC.@preserve x c libSingular.n_Neg(x.ptr, c.ptr)
   return c(ptr)
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(x::n_Z, y::n_Z)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Add(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function -(x::n_Z, y::n_Z)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Sub(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function *(x::n_Z, y::n_Z)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Mult(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function isless(x::n_Z, y::n_Z)
   c = parent(x)
   GC.@preserve x y c libSingular.n_Greater(y.ptr, x.ptr, c.ptr)
end

function ==(x::n_Z, y::n_Z)
   c = parent(x)
   GC.@preserve x y c return libSingular.n_Equal(x.ptr, y.ptr, c.ptr)
end

isequal(x::n_Z, y::n_Z) = (x == y)

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

==(x::n_Z, y::Integer) = (x ==  parent(x)(y))

==(x::Integer, y::n_Z) = (parent(y)(x) == y)

==(x::n_Z, y::Nemo.fmpz) = Nemo.fmpz(convert(BigInt, x)) == y
==(x::Nemo.fmpz, y::n_Z) = y == x

isequal(x::n_Z, y::Integer) = (x == y)

isequal(x::Integer, y::n_Z) = (x == y)

isless(x::n_Z, y::Integer) = isless(x, parent(x)(y))

isless(x::Integer, y::n_Z) = isless(parent(y)(x), y)

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::n_Z, y::Int)
   y < 0 && throw(DomainError(y, "exponent must be non-negative"))
   if isone(x)
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

function div(x::n_Z, y::n_Z)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Div(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function divexact(x::n_Z, y::n_Z; check::Bool=true)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_ExactDiv(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Euclidean division
#
###############################################################################

function divrem(x::n_Z, y::n_Z)
   c = parent(x)
   r = GC.@preserve c [libSingular.n_Init(0, c.ptr)]
   q = GC.@preserve x y r c libSingular.n_QuotRem(x.ptr, y.ptr, pointer(r), c.ptr)
   return c(q), c(r[])
end

function rem(x::n_Z, y::n_Z)
   c = parent(x)
   r = GC.@preserve x y c libSingular.n_IntMod(x.ptr, y.ptr, c.ptr)
   return c(r)
end

function mod(x::n_Z, y::n_Z)
   d = rem(x, y)
   return d < 0 ? d + abs(y) : d
end

###############################################################################
#
#   GCD and LCM
#
###############################################################################

function gcd(x::n_Z, y::n_Z)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Gcd(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function lcm(x::n_Z, y::n_Z)
   if x == 0 && y == 0
      return zero(parent(x))
   end
   return div(x*y, gcd(x, y))
end

###############################################################################
#
#   Extended GCD
#
###############################################################################

function gcdx(x::n_Z, y::n_Z)
   c = parent(x)
   GC.@preserve x y c begin
      s = Ref(libSingular.n_Init(0, c.ptr))
      t = Ref(libSingular.n_Init(0, c.ptr))
      g = libSingular.n_ExtGcd(x.ptr, y.ptr, s, t, c.ptr)
      return c(g), c(s[]), c(t[])
   end
end

###############################################################################
#
#   Chinese remainder
#
###############################################################################

function crt(r1::n_Z, m1::n_Z, r2::n_Z, m2::n_Z, signed=false)
   par = parent(r1)
   r = [r1.ptr, r2.ptr]
   m = [m1.ptr, m2.ptr]
   println("type par",typeof(par.ptr))
   println("type r",typeof(r))
   println("type m",typeof(m))
   a = libSingular.n_ChineseRemainderSym(reinterpret(Ptr{Nothing},pointer(r)), reinterpret(Ptr{Nothing},pointer(m)), Cint(2), Cint(signed), par.ptr)
   println("Got through")
   return par(a)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addeq!(x::n_Z, y::n_Z)
   c = parent(x)
   x.ptr = GC.@preserve x y c libSingular.n_InpAdd(x.ptr, y.ptr, c.ptr)
   return x
end

function mul!(x::n_Z, y::n_Z, z::n_Z)
   c = parent(x)
   GC.@preserve x y z c begin
      ptr = libSingular.n_Mult(y.ptr, z.ptr, c.ptr)
      libSingular.n_Delete(x.ptr, c.ptr)
      x.ptr = ptr
      return x
   end
end

function add!(x::n_Z, y::n_Z, z::n_Z)
   c = parent(x)
   GC.@preserve x y z c begin
      ptr = libSingular.n_Add(y.ptr, z.ptr, c.ptr)
      libSingular.n_Delete(x.ptr, c.ptr)
      x.ptr = ptr
      return x
   end
end

function zero!(x::n_Z)
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

promote_rule(C::Type{n_Z}, ::Type{T}) where {T <: Integer} = n_Z

promote_rule(C::Type{n_Z}, ::Type{Nemo.fmpz}) = n_Z

BigInt(n::n_Z) = libSingular.n_GetMPZ(n.ptr, parent(n).ptr)
Integer(n::n_Z) = BigInt(n)
(::Type{T})(n::n_Z) where {T <: Integer} = T(BigInt(n))

convert(::Type{T}, n::n_Z) where {T <: Integer} = T(n)

# fmpz
Nemo.fmpz(n::n_Z) = Nemo.fmpz(BigInt(n))
convert(::Type{Nemo.fmpz}, n::n_Z) = Nemo.fmpz(n)
(::FlintIntegerRing)(n::n_Z) = Nemo.fmpz(n)

###############################################################################
#
#   Random generation
#
###############################################################################

RandomExtensions.maketype(R::Integers, _) = n_Z

# define rand(make(ZZ, n:m))
rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make2{n_Z, Integers,
                                                  <:AbstractArray{<:Integer}}}) =
   sp[][1](rand(rng, sp[][2]))

rand(rng::AbstractRNG, R::Integers, n) = R(rand(rng, n))

rand(R::Integers, n) = rand(Random.GLOBAL_RNG, R, n)

###############################################################################
#
#   Parent call functions
#
###############################################################################

(::Integers)(n::IntegerLikeTypes = 0) = n_Z(n)

(::Integers)(n::n_Z) = n

(::Integers)(n::libSingular.number_ptr) = n_Z(n)
