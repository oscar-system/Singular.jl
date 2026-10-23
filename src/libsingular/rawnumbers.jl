###############################################################################
#
#   Decoding n_Z and n_Q coefficients in place
#
#   A coefficient is one machine word, either a pointer or a value carrying a
#   tag in its low bits.
#
#     n_Z   always a pointer to a GMP mpz_t. SI_INTEGER_VARIANT 3 would tag
#           small values the way n_Q does, but Singular is built with variant
#           2, which has no immediate form at all.
#
#     n_Q   (SI_BIGINT_VARIANT 1, longrat.cc) an immediate *integer*, or a
#           pointer to an snumber. A fraction is always an snumber.
#
#   An immediate n_Q, by bit:
#
#      63                              2   1   0
#     +---------------------------------+---+---+
#     |             value               | 0 | 1 |
#     +---------------------------------+---+---+
#
#   Only -2^28 .. 2^28-1 ever turns up there, because longrat.cc computes with
#   the decoded value as a C int (MAX_NUM_SIZE 28; its 64-bit branch is
#   compiled out).
#
#   An snumber, by byte:
#
#      0         16         32  36
#     +----------+----------+---+
#     | mpz_t z  | mpz_t n  | s |
#     +----------+----------+---+
#
#   z/n is the value and s says how to read it: 3 for an integer, where n is
#   undefined, 1 for a fraction in lowest terms, 0 for one still to be reduced.
#
#   Every call into Singular costs a ccall through CxxWrap, and one that hands
#   back a number costs a Julia allocation on top; for small coefficients that
#   is the bulk of the work. Decoding the word here costs a few loads instead.
#
###############################################################################

import Nemo

const Limb = Base.GMP.Limb

###############################################################################
#
#   GMP mpz_t
#
###############################################################################

"""
Non-owning view of the header of a GMP `mpz_t` that belongs to Singular. The
limbs live behind `d`, so a view is only valid while the Singular object
holding the `mpz_t` is reachable: every use needs a `GC.@preserve`.
"""
struct mpz_view
   alloc::Cint
   size::Cint
   d::Ptr{Limb}
end

const mpz_view_ptr = Ptr{mpz_view}

@inline mpz_view(p::mpz_view_ptr) = unsafe_load(p)

@inline mpz_is_zero(m::mpz_view) = m.size == 0

@inline mpz_sign(m::mpz_view) = m.size > 0 ? 1 : (m.size < 0 ? -1 : 0)

@inline mpz_is_abs_one(m::mpz_view) = abs(m.size) == 1 && unsafe_load(m.d) == 1

@inline mpz_is_one(m::mpz_view) = m.size == 1 && unsafe_load(m.d) == 1

# The magnitude of a one-limb mpz_t is representable as an Int unless it is
# larger than typemax(Int); -2^63 is the one negative value that still fits.
@inline function mpz_fits_int(m::mpz_view)
   m.size == 0 && return true
   abs(m.size) == 1 || return false
   limb = unsafe_load(m.d)
   return limb <= (m.size > 0 ? Limb(typemax(Int)) : Limb(typemin(Int) % UInt))
end

@inline function mpz_to_int(m::mpz_view)
   m.size == 0 && return 0
   limb = unsafe_load(m.d) % Int
   return m.size > 0 ? limb : -limb
end

"""
    mpz_cmp_int(m::mpz_view, y::Int)

Sign of `m - y`.
"""
@inline function mpz_cmp_int(m::mpz_view, y::Int)
   # anything too large for an Int has a bigger magnitude than y
   mpz_fits_int(m) || return mpz_sign(m)
   return cmp(mpz_to_int(m), y)
end

"""
    mpz_cmp(a::mpz_view_ptr, b::mpz_view_ptr)

Sign of `a - b`.
"""
@inline mpz_cmp(a::mpz_view_ptr, b::mpz_view_ptr) =
   Int(@ccall Base.GMP.libgmp.__gmpz_cmp(a::mpz_view_ptr, b::mpz_view_ptr)::Cint)

"""
    bigint_set!(z::BigInt, m::mpz_view_ptr)

Copy the integer at `m` into `z`.
"""
@inline function bigint_set!(z::BigInt, m::mpz_view_ptr)
   @ccall Base.GMP.libgmp.__gmpz_set(z::Ref{BigInt}, m::mpz_view_ptr)::Nothing
   return z
end

"""
    mpz_to(T::Type{<:Integer}, m::mpz_view_ptr)

The integer at `m` as a `T`.
"""
function mpz_to(::Type{T}, m::mpz_view_ptr) where {T <: Integer}
   v = mpz_view(m)
   mpz_fits_int(v) && return T(mpz_to_int(v))
   return T(bigint_set!(BigInt(), m))
end

"""
    mpz_hash(m::mpz_view_ptr, h::UInt)

`hash(x, h)` for the integer `x` at `m`, without copying it where possible.
"""
function mpz_hash(m::mpz_view_ptr, h::UInt)
   v = mpz_view(m)
   # Base hashes a one-limb BigInt as the corresponding Int
   mpz_fits_int(v) && return hash(mpz_to_int(v), h)
   return hash(bigint_set!(BigInt(), m), h)
end

"""
    fmpz_set!(z, m::mpz_view_ptr)

Copy the integer at `m` into the `ZZRingElem` `z`, which may also be given as a
pointer into a larger FLINT object.
"""
@inline function fmpz_set!(z::Nemo.ZZRingElem, m::mpz_view_ptr)
   @ccall libflint.fmpz_set_mpz(z::Ref{Nemo.ZZRingElem}, m::mpz_view_ptr)::Nothing
   return z
end

@inline function fmpz_set!(z::Ptr{Nemo.ZZRingElem}, m::mpz_view_ptr)
   @ccall libflint.fmpz_set_mpz(z::Ptr{Nemo.ZZRingElem}, m::mpz_view_ptr)::Nothing
   return z
end

"""
    machine_int(y::Integer)

`y` as an `Int`, or `nothing` if it does not fit. For an `Int` argument this
compiles away, so callers can offer a fast path for machine integers without
paying for it.
"""
@inline machine_int(y::Integer) =
   typemin(Int) <= y <= typemax(Int) ? Int(y) : nothing

###############################################################################
#
#   n_Z: the number is the mpz_t
#
###############################################################################

@inline raw(p::number_ptr) = reinterpret(Int, p.cpp_object)

@inline integer_mpz(p::number_ptr) = mpz_view_ptr(raw(p))

###############################################################################
#
#   n_Q: immediate integer or snumber
#
###############################################################################

# values of snumber.s
const RAT_UNNORMALIZED = Cint(0)   # z/n, gcd(z, n) not necessarily 1
const RAT_NORMALIZED   = Cint(1)   # z/n in lowest terms, n > 1
const RAT_INTEGER      = Cint(3)   # z, n undefined

const _SR_INT = 1
const _SNUMBER_N_OFFSET = sizeof(mpz_view)

# `s` follows z and n, except in a Singular built for debugging, where snumber
# carries an extra `debug` field first; `check_number_layout` finds out which.
global _snumber_s_offset::Int = 2 * sizeof(mpz_view)

# Singular computes with immediate values as C ints (MAX_NUM_SIZE 28 in
# longrat.cc), so a wider immediate would silently overflow its arithmetic.
const _SR_MAX = (1 << 28) - 1
const _SR_MIN = -(1 << 28)

@inline rational_is_imm(p::number_ptr) = raw(p) & _SR_INT != 0

@inline rational_imm(p::number_ptr) = raw(p) >> 2

@inline rational_imm_fits(y::Int) = _SR_MIN <= y <= _SR_MAX

@inline rational_from_imm(y::Int) = number_ptr(Ptr{number}((y << 2) + _SR_INT))

@inline rational_state(p::number_ptr) =
   unsafe_load(Ptr{Cint}(raw(p) + _snumber_s_offset))

@inline rational_num_mpz(p::number_ptr) = mpz_view_ptr(raw(p))

@inline rational_den_mpz(p::number_ptr) =
   mpz_view_ptr(raw(p) + _SNUMBER_N_OFFSET)

"""
    rational_to_bigint(p::number_ptr)

The value of `p` as a `BigInt`, where `p` is an `n_Q` or `coeffs_BIGINT` number
that holds an integer. Unlike `n_GetMPZ` this never normalizes `p`, so it is
also safe on a number that Singular owns.
"""
function rational_to_bigint(p::number_ptr)
   rational_is_imm(p) && return BigInt(rational_imm(p))
   return bigint_set!(BigInt(), rational_num_mpz(p))
end

###############################################################################
#
#   FLINT plays the same trick
#
#   An fmpz is one word: the value itself, or a pointer to an mpz_t shifted
#   right by two with the top two bits marking it as such. A fmpq is just two
#   fmpz in a row, kept coprime with a positive denominator.
#
###############################################################################

@inline fmpz_is_small(d::Int) = unsigned(d) >> (Sys.WORD_SIZE - 2) != 1

@inline fmpz_is_small(a::Nemo.ZZRingElem) = fmpz_is_small(a.d)

# pointer to the mpz_t backing an fmpz that is not an immediate value
@inline fmpz_mpz(d::Int) = mpz_view_ptr(d << 2)

@inline fmpz_mpz(a::Nemo.ZZRingElem) = fmpz_mpz(a.d)

# The two halves of a QQFieldElem can be written in place; only valid under
# GC.@preserve of `a`.
@inline fmpq_num_ptr(a::Nemo.QQFieldElem) =
   Ptr{Nemo.ZZRingElem}(pointer_from_objref(a))

@inline fmpq_den_ptr(a::Nemo.QQFieldElem) =
   fmpq_num_ptr(a) + sizeof(Nemo.ZZRingElem)

###############################################################################
#
#   Building numbers
#
###############################################################################

# initialise a number from an mpz owned by someone else, who keeps it alive
@inline n_InitMPZ(m::mpz_view_ptr, cf::coeffs_ptr) =
   n_InitMPZ_internal(reinterpret(Ptr{Nothing}, m), cf)

"""
    number_ptr_from_fmpz(d::Int, cf::coeffs_ptr)

A number in `cf` equal to the integer stored in the fmpz word `d`. The caller
keeps the fmpz alive.
"""
function number_ptr_from_fmpz(d::Int, cf::coeffs_ptr)
   fmpz_is_small(d) && return n_Init(d, cf)
   return n_InitMPZ(fmpz_mpz(d), cf)
end

"""
    rational_ptr_from_fmpz(d::Int, qq::coeffs_ptr)

An `n_Q` in `qq` equal to the integer stored in the fmpz word `d`, as an
immediate value where possible. The caller keeps the fmpz alive.
"""
function rational_ptr_from_fmpz(d::Int, qq::coeffs_ptr)
   fmpz_is_small(d) && rational_imm_fits(d) && return rational_from_imm(d)
   return number_ptr_from_fmpz(d, qq)
end

###############################################################################
#
#   Layout check
#
#   Everything above hardcodes how Singular represents a coefficient, so
#   `check_number_layout` confirms that at startup: a mismatched libsingular
#   must fail loudly here rather than read whatever the word happens to point
#   at. It also measures the one part of the layout that legitimately varies.
#
###############################################################################

"""
    check_number_layout(zz::coeffs_ptr, qq::coeffs_ptr)

Verify that `n_Z` and `n_Q` coefficients are laid out as the accessors in this
file assume, raising an error if they are not, and record where `snumber` keeps
its state field.
"""
function check_number_layout(zz::coeffs_ptr, qq::coeffs_ptr)
   _check_integer_layout(zz)
   _check_immediate_layout(qq)
   _check_snumber_layout(qq)
   nothing
end

function _layout_error(what)
   error("""
         Singular.jl cannot decode $what as expected, which means the
         libsingular it runs against represents coefficients differently.
         Check SI_INTEGER_VARIANT and SI_BIGINT_VARIANT in Singular's
         libpolys/misc/auxiliary.h against src/libsingular/rawnumbers.jl, and
         please report this at https://github.com/oscar-system/Singular.jl/issues
         """)
end

# everything the accessors above assume about one n_Z
function _check_integer(p::number_ptr, expected::BigInt)
   raw(p) & 3 == 0 || _layout_error("the n_Z $expected (expected a pointer)")
   m = mpz_view(integer_mpz(p))
   fits = typemin(Int) <= expected <= typemax(Int)
   (mpz_sign(m) == sign(expected) &&
    mpz_fits_int(m) == fits &&
    (!fits || mpz_to_int(m) == expected) &&
    bigint_set!(BigInt(), integer_mpz(p)) == expected) ||
      _layout_error("the n_Z $expected")
end

function _check_integer_layout(zz::coeffs_ptr)
   # the ends of the one-limb range, where mpz_fits_int and mpz_to_int are
   # delicate: -2^63 still fits an Int, 2^63 would not
   for v in (0, 1, -1, typemax(Int), typemin(Int))
      p = n_Init(v, zz)
      _check_integer(p, BigInt(v))
      n_Delete(p, zz)
   end

   # and two limbs, which must not be mistaken for a machine integer
   two = n_Init(2, zz)
   big = n_Power(two, 64, zz)
   negbig = n_Neg(big, zz)
   _check_integer(big, BigInt(2)^64)
   _check_integer(negbig, -BigInt(2)^64)
   for p in (negbig, big, two)
      n_Delete(p, zz)
   end
end

function _check_immediate_layout(qq::coeffs_ptr)
   # 0 and 1 because iszero and isone compare against their encodings, and the
   # ends of the range because rational_from_imm must not hand Singular a value
   # that overflows its own int arithmetic. Nothing checks that 2^28 is *not*
   # immediate: a Singular that widened the range would still be safe for us.
   for v in (0, 1, -1, _SR_MIN, _SR_MAX)
      p = n_Init(v, qq)
      (rational_is_imm(p) && rational_imm(p) == v) ||
         _layout_error("the immediate n_Q $v")
      rational_from_imm(v) == p ||
         _layout_error("the encoding of the immediate n_Q $v")
      n_Delete(p, qq)
   end
end

# `frac` is a normalized fraction and `bigint` an integer, so their states pin
# down where snumber keeps `s`: behind an extra `debug` field in a Singular
# built for debugging, directly after z and n otherwise.
function _find_snumber_s_offset(frac::number_ptr, bigint::number_ptr)
   base = 2 * sizeof(mpz_view)
   for offset in (base, base + sizeof(Cint))
      unsafe_load(Ptr{Cint}(raw(frac) + offset)) == RAT_NORMALIZED || continue
      unsafe_load(Ptr{Cint}(raw(bigint) + offset)) == RAT_INTEGER || continue
      return offset
   end
   _layout_error("the state of an n_Q")
end

function _check_snumber_layout(qq::coeffs_ptr)
   three = n_Init(3, qq)
   seven = n_Init(7, qq)
   two = n_Init(2, qq)
   frac = n_Div(three, seven, qq)
   # wider than any MAX_NUM_SIZE, hence never immediate
   bigint = n_Power(two, 64, qq)

   (!rational_is_imm(frac) && !rational_is_imm(bigint)) ||
      _layout_error("a large n_Q (expected a pointer)")
   (mpz_to_int(mpz_view(rational_num_mpz(frac))) == 3 &&
    mpz_to_int(mpz_view(rational_den_mpz(frac))) == 7 &&
    bigint_set!(BigInt(), rational_num_mpz(bigint)) == BigInt(2)^64) ||
      _layout_error("the numerator and denominator of an n_Q")
   global _snumber_s_offset = _find_snumber_s_offset(frac, bigint)

   for p in (frac, bigint, three, seven, two)
      n_Delete(p, qq)
   end
end
