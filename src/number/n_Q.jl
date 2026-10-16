export n_Q, Rationals, reconstruct, isone, iszero, is_unit, divexact,
       canonical_unit

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

elem_type(::Type{Rationals}) = n_Q

parent(a::n_Q) = QQ
parent_type(::Type{n_Q}) = Rationals

base_ring(a::Rationals) = ZZ
base_ring_type(::Type{Rationals}) = Integers

characteristic(R::Rationals) = ZZ(0)

function deepcopy_internal(a::n_Q, dict::IdDict)
   c = parent(a)
   GC.@preserve a c return c(libSingular.n_Copy(a.ptr, c.ptr))
end

# An n_Q is either an immediate integer or a pointer to Singular's
# snumber { mpz_t z; mpz_t n; int s }; see src/libsingular/rawnumbers.jl. Views
# into the latter borrow Singular's limbs, so every use must sit inside a
# GC.@preserve of the n_Q.
#
# Singular may leave a rational unnormalized (s == 0), in which case z and n
# need not be coprime. Normalizing can replace the snumber, so this returns the
# pointer to use from here on.
@inline function _normalized(x::n_Q)
   p = x.ptr
   if !libSingular.rational_is_imm(p) &&
         libSingular.rational_state(p) == libSingular.RAT_UNNORMALIZED
      c = parent(x)
      p = GC.@preserve x c libSingular.n_Normalize(p, c.ptr)
      x.ptr = p
   end
   return p
end

@inline _is_fraction(p::libSingular.number_ptr) =
   !libSingular.rational_is_imm(p) &&
      libSingular.rational_state(p) != libSingular.RAT_INTEGER

function hash(a::n_Q, h::UInt)
   GC.@preserve a begin
      p = _normalized(a)
      if libSingular.rational_is_imm(p)
         nhash = hash(libSingular.rational_imm(p), h)
         dhash = hash(1, h)
      else
         nhash = libSingular.mpz_hash(libSingular.rational_num_mpz(p), h)
         dhash = _is_fraction(p) ?
                 libSingular.mpz_hash(libSingular.rational_den_mpz(p), h) :
                 hash(1, h)
      end
   end
   return xor(xor(nhash, dhash), 0xf348b78c190e8fc1%UInt)
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

one(::Rationals) = QQ(1)

zero(::Rationals) = QQ(0)

# Singular shortens every representable integer to an immediate, so 0 and 1
# only ever occur in that form
isone(n::n_Q) = libSingular.rational_is_imm(n.ptr) &&
                isone(libSingular.rational_imm(n.ptr))

iszero(n::n_Q) = libSingular.rational_is_imm(n.ptr) &&
                 iszero(libSingular.rational_imm(n.ptr))

is_unit(n::n_Q) = !iszero(n)

@doc raw"""
    numerator(x::n_Q)

Return in `ZZ` the numerator of $x$.
"""
function numerator(x::n_Q)
   GC.@preserve x begin
      p = _normalized(x)
      libSingular.rational_is_imm(p) && return ZZ(libSingular.rational_imm(p))
      return ZZ(libSingular.n_InitMPZ(libSingular.rational_num_mpz(p), ZZ.ptr))
   end
end

@doc raw"""
    denominator(x::n_Q)

Return in `ZZ` the denominator of $x$.
"""
function denominator(x::n_Q)
   GC.@preserve x begin
      p = _normalized(x)
      _is_fraction(p) || return one(ZZ)
      return ZZ(libSingular.n_InitMPZ(libSingular.rational_den_mpz(p), ZZ.ptr))
   end
end

function is_negative(x::n_Q)
   GC.@preserve x begin
      p = x.ptr
      libSingular.rational_is_imm(p) && return libSingular.rational_imm(p) < 0
      return libSingular.mpz_sign(
                libSingular.mpz_view(libSingular.rational_num_mpz(p))) < 0
   end
end

abs(n::n_Q) = is_negative(n) ? -n : deepcopy(n)

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
   print(io, "Rational field")
end

function expressify(n::n_Q; context = nothing)::Any
   return expressify(Rational{BigInt}(n), context = context)
end

AbstractAlgebra.@enable_all_show_via_expressify n_Q

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

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

# Sign of x - y, or `nothing` when x is a proper fraction, where deciding the
# comparison needs more than a look at the numerator.
@inline function _cmp_machine_int(x::n_Q, y::Int)
   GC.@preserve x begin
      p = _normalized(x)
      libSingular.rational_is_imm(p) &&
         return cmp(libSingular.rational_imm(p), y)
      _is_fraction(p) && return nothing
      return libSingular.mpz_cmp_int(
                libSingular.mpz_view(libSingular.rational_num_mpz(p)), y)
   end
end

function _cmp(x::n_Q, y::n_Q)
   c = parent(x)
   GC.@preserve x y c begin
      libSingular.n_Equal(x.ptr, y.ptr, c.ptr) && return 0
      return libSingular.n_Greater(x.ptr, y.ptr, c.ptr) ? 1 : -1
   end
end

function _cmp(x::n_Q, y::Integer)
   yy = libSingular.machine_int(y)
   if yy !== nothing
      c = _cmp_machine_int(x, yy)
      c !== nothing && return c
   end
   return _cmp(x, QQ(y))
end

_cmp(x::n_Q, y::Nemo.ZZRingElem) =
   libSingular.fmpz_is_small(y) ? _cmp(x, y.d) : _cmp(x, QQ(y))

_cmp(x::n_Q, y::Union{Rational, Nemo.QQFieldElem}) = _cmp(x, QQ(y))

const RationalLikeTypes = Union{Integer, Rational, Nemo.ZZRingElem,
                                Nemo.QQFieldElem}

==(x::n_Q, y::RationalLikeTypes) = _cmp(x, y) == 0

==(x::RationalLikeTypes, y::n_Q) = _cmp(y, x) == 0

isless(x::n_Q, y::RationalLikeTypes) = _cmp(x, y) < 0

isless(x::RationalLikeTypes, y::n_Q) = _cmp(y, x) > 0

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
   z = c(p)
   libSingular.check_error()
   return z
end

function divexact(x::n_Q, y::n_Q; check::Bool=true)
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

@doc raw"""
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

function add!(x::n_Q, y::n_Q)
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

promote_rule(C::Type{n_Q}, ::Type{Nemo.ZZRingElem}) = n_Q

promote_rule(C::Type{n_Q}, ::Type{n_Z}) = n_Q

function Rational{T}(x::n_Q) where {T <: Integer}
   GC.@preserve x begin
      p = _normalized(x)
      libSingular.rational_is_imm(p) &&
         return Rational{T}(libSingular.rational_imm(p))
      n = libSingular.mpz_to(T, libSingular.rational_num_mpz(p))
      d = _is_fraction(p) ?
          libSingular.mpz_to(T, libSingular.rational_den_mpz(p)) : one(T)
      # Singular keeps z/n in lowest terms with n > 0
      return Base.unsafe_rational(n, d)
   end
end

Rational(x::n_Q) = Rational{BigInt}(x)

convert(::Type{T}, x::n_Q) where {T <: Rational} = T(x)

function Nemo.QQFieldElem(x::n_Q)
   GC.@preserve x begin
      p = _normalized(x)
      libSingular.rational_is_imm(p) &&
         return QQFieldElem(libSingular.rational_imm(p))
      # Singular keeps z/n in lowest terms with n > 0, which is exactly FLINT's
      # canonical form, so both halves can be filled in directly
      z = QQFieldElem()
      GC.@preserve z begin
         libSingular.fmpz_set!(libSingular.fmpq_num_ptr(z),
                               libSingular.rational_num_mpz(p))
         _is_fraction(p) && libSingular.fmpz_set!(libSingular.fmpq_den_ptr(z),
                                                  libSingular.rational_den_mpz(p))
      end
      return z
   end
end

convert(::Type{QQFieldElem}, x::n_Q) = QQFieldElem(x)
(::QQField)(x::n_Q) = QQFieldElem(x)

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

(R::Rationals)(n::Union{Integer,ZZRingElem,n_Z}, m::Union{Integer,ZZRingElem,n_Z}) = R(n) // R(m)

# take ownership of the pointer - not for general users
(::Rationals)(n::libSingular.number_ptr) = n_Q(n)

# from integers

(::Rationals)(n::IntegerLikeTypes) = n_Q(n)

# from rationals

(::Rationals)(x::Union{Rational, QQFieldElem}) = n_Q(x)

(::Rationals)(n::n_Q) = n
