###############################################################################
#
#   Decoding n_Z and n_Q coefficients in place
#
#   Singular hands each coefficient around as a single machine word: either a
#   pointer, or a small value stored in the word itself and marked by a tag in
#   its lowest bits. Which of the two depends on the coefficient domain.
#
#     n_Z   always a pointer to a GMP integer (mpz_t, see `mpz_view` below).
#           SI_INTEGER_VARIANT 3 would store small values in the word as n_Q
#           does, but Singular is built with variant 2, which never does.
#
#     n_Q   (SI_BIGINT_VARIANT 1, implemented in Singular's longrat.cc) either
#           an integer stored in the word itself, or a pointer to an snumber.
#           A fraction is always an snumber.
#
#   An n_Q integer stored in the word looks like this, bit 63 on the left:
#
#      63                              2   1   0
#     +---------------------------------+---+---+
#     |             value               | 0 | 1 |
#     +---------------------------------+---+---+
#
#   Bit 0 marks it; a pointer never has it set. The value is the word shifted
#   right by two, arithmetically, so a negative value comes back sign-extended.
#   longrat.cc only ever stores values in -2^28 .. 2^28-1 this way, because
#   its arithmetic adds and subtracts these words as C ints (MAX_NUM_SIZE 28;
#   its 64-bit variant is compiled out): a wider value stored here would
#   overflow and silently give wrong results. We decode any width but only
#   ever create values in that range.
#
#   longrat.cc also turns every integer result in that range into this form,
#   so 0 and 1 never occur as an snumber. iszero and isone below rely on that,
#   as Singular's own nlIsZero and nlIsOne do.
#
#   An snumber has the following layout, the labels now giving byte offsets:
#
#      0         16         32  36
#     +----------+----------+---+
#     | mpz_t z  | mpz_t n  | s |
#     +----------+----------+---+
#
#   Each mpz_t takes 16 bytes and s, a C int, 4. z is the numerator (German
#   "Zähler") and n the denominator ("Nenner"); s says how to read them:
#
#     3   an integer: the value is z, and n is undefined
#     1   the fraction z/n in lowest terms, with n > 1
#     0   the fraction z/n, possibly not in lowest terms, with n > 1
#
#   longrat.cc never sets s to 2; it tests s < 2 for "has a denominator".
#
#   The functions in this file read these words directly instead of calling
#   into Singular for each coefficient.
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
    mpz_view

Copy of the header of a GMP integer (`mpz_t`, see
https://gmplib.org/manual/Integer-Internals): `size` is the number of limbs
in use, negated for a negative value, and `d` points to the limbs, least
significant first. The limbs are not copied, so a view is only valid while
the Singular object it came from is alive: every use needs a `GC.@preserve`
of that object.
"""
struct mpz_view
   alloc::Cint
   size::Cint
   d::Ptr{Limb}
end

const mpz_view_ptr = Ptr{mpz_view}

@inline mpz_view(p::mpz_view_ptr) = unsafe_load(p)

@inline Base.iszero(m::mpz_view) = m.size == 0

@inline Base.isone(m::mpz_view) = m.size == 1 && unsafe_load(m.d) == 1

@inline AbstractAlgebra.is_unit(m::mpz_view) =
   abs(m.size) == 1 && unsafe_load(m.d) == 1

@inline Base.sign(m::mpz_view) = Int(sign(m.size))

# The magnitude of a one-limb mpz_t is representable as an Int unless it is
# larger than typemax(Int); -2^63 is the one negative value that still fits.
@inline function Nemo.fits(::Type{Int}, m::mpz_view)
   m.size == 0 && return true
   abs(m.size) == 1 || return false
   limb = unsafe_load(m.d)
   return limb <= (m.size > 0 ? Limb(typemax(Int)) : Limb(typemin(Int) % UInt))
end

# the value of `m`, which the caller has checked fits an Int
@inline function _int(m::mpz_view)
   m.size == 0 && return 0
   limb = unsafe_load(m.d) % Int
   return m.size > 0 ? limb : -limb
end

@inline (::Type{T})(m::mpz_view) where {T <: Integer} =
   Nemo.fits(Int, m) ? T(_int(m)) : T(Nemo.set!(BigInt(), m))

@inline Nemo.ZZRingElem(m::mpz_view) =
   Nemo.fits(Int, m) ? Nemo.ZZRingElem(_int(m)) : Nemo.set!(Nemo.ZZRingElem(), m)

@inline function Base.cmp(m::mpz_view, y::Int)
   # anything too large for an Int has a bigger magnitude than y
   Nemo.fits(Int, m) || return sign(m)
   return cmp(_int(m), y)
end

@inline Base.cmp(a::mpz_view, b::mpz_view) =
   Int(@ccall Base.GMP.libgmp.__gmpz_cmp(Ref(a)::Ref{mpz_view},
                                          Ref(b)::Ref{mpz_view})::Cint)

# Base hashes a one-limb BigInt as the corresponding Int, so this agrees with
# hashing the BigInt without copying it where possible
Base.hash(m::mpz_view, h::UInt) =
   Nemo.fits(Int, m) ? hash(_int(m), h) : hash(BigInt(m), h)

@inline function Nemo.set!(z::BigInt, m::mpz_view)
   @ccall Base.GMP.libgmp.__gmpz_set(z::Ref{BigInt}, Ref(m)::Ref{mpz_view})::Nothing
   return z
end

# `z` may also point into a larger FLINT object, such as half of a QQFieldElem
@inline function Nemo.set!(z::Nemo.ZZRingElem, m::mpz_view)
   @ccall libflint.fmpz_set_mpz(z::Ref{Nemo.ZZRingElem}, Ref(m)::Ref{mpz_view})::Nothing
   return z
end

@inline function Nemo.set!(z::Ptr{Nemo.ZZRingElem}, m::mpz_view)
   @ccall libflint.fmpz_set_mpz(z::Ptr{Nemo.ZZRingElem}, Ref(m)::Ref{mpz_view})::Nothing
   return z
end

# Overwrite the mpz_t at `z`, which must belong to a Singular number no one
# else holds yet.
@inline function Nemo.set!(z::mpz_view_ptr, x::Int)
   @ccall Base.GMP.libgmp.__gmpz_set_si(z::mpz_view_ptr, x::Clong)::Nothing
   return z
end

@inline function Nemo.set!(z::mpz_view_ptr, x::BigInt)
   @ccall Base.GMP.libgmp.__gmpz_set(z::mpz_view_ptr, x::Ref{BigInt})::Nothing
   return z
end

@inline function Nemo.set!(z::mpz_view_ptr, m::mpz_view)
   @ccall Base.GMP.libgmp.__gmpz_set(z::mpz_view_ptr, Ref(m)::Ref{mpz_view})::Nothing
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

# values of snumber.s, see the top of this file
const RAT_UNNORMALIZED = Cint(0)
const RAT_NORMALIZED   = Cint(1)
const RAT_INTEGER      = Cint(3)

const _SR_INT = 1
const _SNUMBER_N_OFFSET = sizeof(mpz_view)

# a Singular built without SING_NDEBUG puts an extra `debug` field before `s`
const _SNUMBER_S_OFFSET = 2 * sizeof(mpz_view)

# the range longrat.cc stores in the word, see the top of this file
const _SR_MAX = (1 << 28) - 1
const _SR_MIN = -(1 << 28)

@inline rational_is_imm(p::number_ptr) = raw(p) & _SR_INT != 0

@inline rational_imm(p::number_ptr) = raw(p) >> 2

@inline rational_imm_fits(y::Int) = _SR_MIN <= y <= _SR_MAX

@inline rational_from_imm(y::Int) = number_ptr(Ptr{number}((y << 2) + _SR_INT))

@inline rational_state(p::number_ptr) =
   unsafe_load(Ptr{Cint}(raw(p) + _SNUMBER_S_OFFSET))

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
   return BigInt(mpz_view(rational_num_mpz(p)))
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

"""
    rational_ptr(y::Integer, qq::coeffs_ptr)

An `n_Q` in `qq` equal to `y`, stored in the word itself where possible.
"""
function rational_ptr(y::Integer, qq::coeffs_ptr)
   v = machine_int(y)
   v !== nothing && rational_imm_fits(v) && return rational_from_imm(v)
   return number_ptr(y, qq)
end

"""
    rational_fraction(num, den::number_ptr, qq::coeffs_ptr)

The `n_Q` `num/den`, for an `n_Q` integer `den > 1` coprime to `num`, which
is an `Int`, a `BigInt` or an `mpz_view`. Takes ownership of `den`.

Dividing would store num/den unreduced and then spend a gcd finding nothing
to cancel. Inverting `den` gives the snumber {1, den, s = 1} without one
(nlInvers in longrat.cc), and only its numerator is left to overwrite.
"""
function rational_fraction(num, den::number_ptr, qq::coeffs_ptr)
   p = n_Invers(den, qq)
   rational_is_imm(den) || n_Delete(den, qq)
   Nemo.set!(rational_num_mpz(p), num)
   return p
end
