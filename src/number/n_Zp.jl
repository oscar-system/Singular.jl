export n_Zp, N_ZpField

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

elem_type(::Type{N_ZpField}) = n_Zp

parent(a::n_Zp) = a.parent

parent_type(::Type{n_Zp}) = N_ZpField

base_ring(a::n_Zp) = Union{}

base_ring(a::N_ZpField) = Union{}

function characteristic(R::N_ZpField)
   return ZZ(libSingular.n_GetChar(R.ptr))
end

function deepcopy_internal(a::n_Zp, dict::IdDict)
   c = parent(a)
   GC.@preserve a c return parent(a)(libSingular.n_Copy(a.ptr, c.ptr))
end

function hash(a::n_Zp, h::UInt)
   chash = hash(characteristic(parent(a)), h)
   ahash = hash(Int(a), h)
   return xor(xor(chash, ahash), 0x77dc334c1532ce3c%UInt)
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

one(R::N_ZpField) = R(1)

zero(R::N_ZpField) = R(0)

function isone(n::n_Zp)
   c = parent(n)
   GC.@preserve n c return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_Zp)
   c = parent(n)
   GC.@preserve n return libSingular.n_IsZero(n.ptr, c.ptr)
end

is_unit(n::n_Zp) = !iszero(n)

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(x::n_Zp) = x

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, c::N_ZpField)
   print(io, "Finite Field of Characteristic ", characteristic(c))
end

function expressify(n::n_Zp; context = nothing)::Any
  nn = rem(Int(n), Int(characteristic(parent(n))), RoundNearest)
  return expressify(nn, context = context)
end

AbstractAlgebra.@enable_all_show_via_expressify n_Zp

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(x::n_Zp)
   c = parent(x)
   p = GC.@preserve x c libSingular.n_Neg(x.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(x::n_Zp, y::n_Zp)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Add(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function -(x::n_Zp, y::n_Zp)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Sub(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function *(x::n_Zp, y::n_Zp)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Mult(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(x::n_Zp, y::n_Zp)
   c = parent(x)
   GC.@preserve x y c return libSingular.n_Equal(x.ptr, y.ptr, c.ptr)
end

isequal(x::n_Zp, y::n_Zp) = (x == y)

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

==(x::n_Zp, y::Integer) = (x ==  parent(x)(y))

==(x::Integer, y::n_Zp) = (parent(y)(x) == y)

==(x::n_Zp, y::n_Z) = (x ==  parent(x)(y))

==(x::n_Z, y::n_Zp) = (parent(y)(x) == y)

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::n_Zp, y::Int)
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

function inv(x::n_Zp)
   c = parent(x)
   p = GC.@preserve x c libSingular.n_Invers(x.ptr, c.ptr)
   z = c(p)
   libSingular.check_error()
   return z
end

function divexact(x::n_Zp, y::n_Zp; check::Bool=true)
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

function gcd(x::n_Zp, y::n_Zp)
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

function addeq!(x::n_Zp, y::n_Zp)
   c = parent(x)
   x.ptr = GC.@preserve x y c libSingular.n_InpAdd(x.ptr, y.ptr, c.ptr)
   return x
end

function mul!(x::n_Zp, y::n_Zp, z::n_Zp)
   c = parent(x)
   GC.@preserve x y z c begin
      ptr = libSingular.n_Mult(y.ptr, z.ptr, c.ptr)
      libSingular.n_Delete(x.ptr, c.ptr)
      x.ptr = ptr
      return x
   end
end

function add!(x::n_Zp, y::n_Zp, z::n_Zp)
   c = parent(x)
   GC.@preserve x y z c begin
      ptr = libSingular.n_Add(y.ptr, z.ptr, c.ptr)
      libSingular.n_Delete(x.ptr, c.ptr)
      x.ptr = ptr
      return x
   end
end

function zero!(x::n_Zp)
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
#   Random functions
#
###############################################################################

# define rand(::GaloisField)

Random.gentype(::Type{N_ZpField}) = elem_type(N_ZpField)

Random.Sampler(::Type{RNG}, R::N_ZpField, n::Random.Repetition) where {RNG<:AbstractRNG} =
   Random.SamplerSimple(R, Random.Sampler(RNG, Int(0):Int(characteristic(R)) - 1, n))

rand(rng::AbstractRNG, R::Random.SamplerSimple{N_ZpField}) = R[](rand(rng, R.data))

# define rand(make(::N_ZpField, dist))

RandomExtensions.maketype(R::N_ZpField, _) = elem_type(R) # n_Zp

rand(rng::AbstractRNG, sp::SamplerTrivial{<:Make2{n_Zp,N_ZpField}}) =
   sp[][1](rand(rng, sp[][2]))

# define rand(::N_ZpField, integer_array)
# we restrict to array so that the `rand` method producing arrays (e.g. rand(R, 3)) works

rand(rng::AbstractRNG, R::N_ZpField, b::AbstractArray{<:Integer}) = rand(rng, make(R, b))

rand(R::N_ZpField, b::AbstractArray{<:Integer}) = rand(Random.GLOBAL_RNG, R, b)


###############################################################################
#
#   Conversions and promotions
#
###############################################################################

promote_rule(C::Type{n_Zp}, ::Type{T}) where {T <: Integer} = n_Zp

promote_rule(C::Type{n_Zp}, ::Type{Nemo.fmpz}) = n_Zp

promote_rule(C::Type{n_Zp}, ::Type{n_Z}) = n_Zp

###############################################################################
#
#   Parent call functions
#
###############################################################################

(R::N_ZpField)(n::IntegerLikeTypes = 0) = n_Zp(R, n)

(R::N_ZpField)(n::n_Zp) = n

# take ownership of the pointer - not for general users
(R::N_ZpField)(n::libSingular.number_ptr) = n_Zp(R, n)

###############################################################################
#
#   Fp constructor
#
###############################################################################

function Fp(a::Int; cached=true)
   a == 0 && throw(DivideError(a))
   a < 0 && throw(DomainError(a, "prime must be positive"))
   a > 2^29 && throw(DomainError(a, "prime must be <= 2^29"))
   !Nemo.isprime(Nemo.fmpz(a)) && throw(DomainError(a, "characteristic must be prime"))

   return N_ZpField(a, cached)
end

function Base.Int(a::n_Zp)
   GC.@preserve a return reinterpret(Int, a.ptr.cpp_object)
end
