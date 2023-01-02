export n_Zn, N_ZnRing, characteristic

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

elem_type(::Type{N_ZnRing}) = n_Zn

parent(a::n_Zn) = a.parent

parent_type(::Type{n_Zn}) = N_ZnRing

base_ring(a::n_Zn) = ZZ

base_ring(a::N_ZnRing) = ZZ

function characteristic(R::N_ZnRing)
   return ZZ(libSingular.n_GetChar(R.ptr))
end

function deepcopy_internal(a::n_Zn, dict::IdDict)
   return parent(a)(libSingular.n_Copy(a.ptr, parent(a).ptr))
end

function hash(a::n_Zn, h::UInt)
   # get the number of limbs
   sptr = convert(Ptr{Cint}, a.ptr.cpp_object) + sizeof(Cint)
   s = unsafe_load(sptr)
   # get the pointer after the first two Cints
   d = convert(Ptr{Ptr{UInt}}, a.ptr.cpp_object) + 2*sizeof(Cint)
   p = unsafe_load(d)
   b = unsafe_load(p)
   h = xor(Base.hash_uint(xor(ifelse(s < 0, -b, b), h)), h)
   for k = 2:abs(s)
      h = xor(Base.hash_uint(xor(unsafe_load(p, k), h)), h)
   end
   return xor(h, 0xe6ebab8a56a5461b%UInt)
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

one(R::N_ZnRing) = R(1)

zero(R::N_ZnRing) = R(0)

function isone(n::n_Zn)
   c = parent(n)
   GC.@preserve n c return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_Zn)
   c = parent(n)
   GC.@preserve n c return libSingular.n_IsZero(n.ptr, c.ptr)
end

is_unit(n::n_Zn) = gcd(n, parent(n)(characteristic(parent(n)))) == 1

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(x::n_Zn) = x

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, c::N_ZnRing)
   print(io, "Residue Ring of Integer Ring modulo ", characteristic(c))
end

function BigInt(n::n_Zn)
   z = BigInt()
   GC.@preserve n begin
      ccall((:__gmpz_set, :libgmp), Cvoid,
            (Ref{BigInt}, Ptr{Any}),
            z, n.ptr.cpp_object)
   end
   return z
end

function expressify(n::n_Zn; context = nothing)
   return BigInt(n)
end

AbstractAlgebra.@enable_all_show_via_expressify n_Zn

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(x::n_Zn)
   c = parent(x)
   p = GC.@preserve x c libSingular.n_Neg(x.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(x::n_Zn, y::n_Zn)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Add(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function -(x::n_Zn, y::n_Zn)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Sub(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function *(x::n_Zn, y::n_Zn)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Mult(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(x::n_Zn, y::n_Zn)
   c = parent(x)
   GC.@preserve x y c return libSingular.n_Equal(x.ptr, y.ptr, c.ptr)
end

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

==(x::n_Zn, y::Integer) = (x ==  parent(x)(y))

==(x::Integer, y::n_Zn) = (parent(y)(x) == y)

==(x::n_Zn, y::n_Z) = (x ==  parent(x)(y))

==(x::n_Z, y::n_Zn) = (parent(y)(x) == y)

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::n_Zn, y::Int)
   y < 0 && throw(DomainError(y, "exponent must be non-negative"))
   if isone(x)
      return x
   elseif y == 0
      return one(parent(x))
   elseif y == 1
      return x
   else
      c = parent(x)
      p = GC.@preserve x y c libSingular.n_Power(x.ptr, y, c.ptr)
      return c(p)
   end
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(x::n_Zn, y::n_Zn; check::Bool=true)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_ExactDiv(x.ptr, y.ptr, c.ptr)
   z = c(p)
   libSingular.check_error()
   return z
end

###############################################################################
#
#   GCD and LCM
#
###############################################################################

function gcd(x::n_Zn, y::n_Zn)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Gcd(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Extended GCD
#
###############################################################################

function gcdx(x::n_Zn, y::n_Zn)
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
#   Unsafe functions
#
###############################################################################

function addeq!(x::n_Zn, y::n_Zn)
   c = parent(x)
   x.ptr = GC.@preserve x y c libSingular.n_InpAdd(x.ptr, y.ptr, c.ptr)
   return x
end

function mul!(x::n_Zn, y::n_Zn, z::n_Zn)
   c = parent(x)
   GC.@preserve x y z c begin
      ptr = libSingular.n_Mult(y.ptr, z.ptr, c.ptr)
      libSingular.n_Delete(x.ptr, c.ptr)
      x.ptr = ptr
      return x
   end
end

function add!(x::n_Zn, y::n_Zn, z::n_Zn)
   c = parent(x)
   GC.@preserve x y z c begin
      ptr = libSingular.n_Add(y.ptr, z.ptr, c.ptr)
      libSingular.n_Delete(x.ptr, c.ptr)
      x.ptr = ptr
      return x
   end
end

function zero!(x::n_Zn)
   c = parent(x)
   GC.@preserve x begin
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

promote_rule(C::Type{n_Zn}, ::Type{T}) where {T <: Integer} = n_Zn

promote_rule(C::Type{n_Zn}, ::Type{Nemo.fmpz}) = n_Zn

promote_rule(C::Type{n_Zn}, ::Type{n_Z}) = n_Zn

###############################################################################
#
#   Parent call functions
#
###############################################################################

(R::N_ZnRing)(n::IntegerLikeTypes = 0) = n_Zn(R, n)

(R::N_ZnRing)(n::n_Zn) = n

# take ownership of the pointer - not for general users
(R::N_ZnRing)(n::libSingular.number_ptr) = n_Zn(R, n)

(R::N_ZnRing)(n::Nemo.nmod) = n_Zn(R, lift(n))

(R::N_ZnRing)(n::Nemo.fmpz_mod) = n_Zn(R, lift(n))

(R::Nemo.NmodRing)(n::n_Zn) = R(BigInt(n))

(R::Nemo.FmpzModRing)(n::n_Zn) = R(BigInt(n))

###############################################################################
#
#   SingularResidueRing constructor
#
###############################################################################

ResidueRing(R::Integers, a::Int; cached=true) = N_ZnRing(BigInt(a), cached)

ResidueRing(R::Integers, a::BigInt; cached=true) = N_ZnRing(a, cached)
