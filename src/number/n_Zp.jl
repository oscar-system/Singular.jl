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

@doc Markdown.doc"""
    characteristic(R::N_ZpField)
> Return the characteristic of the field.
"""
function characteristic(R::N_ZpField)
   return ZZ(libSingular.n_GetChar(R.ptr))
end

function deepcopy_internal(a::n_Zp, dict::IdDict)
   return parent(a)(libSingular.n_Copy(a.ptr, parent(a).ptr))
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
   return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_Zp)
   c = parent(n)
   return libSingular.n_IsZero(n.ptr, c.ptr)
end

@doc Markdown.doc"""
    isunit(n::n_Zp)
> Return `true` if $n$ is a unit in the field, i.e. nonzero.
"""
isunit(n::n_Zp) = !iszero(n)

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

function show(io::IO, n::n_Zp)
   libSingular.StringSetS("")
   libSingular.n_Write(n.ptr, parent(n).ptr, false)
   m = libSingular.StringEndS()
   print(io, m)
end

needs_parentheses(x::n_Zp) = false

function isnegative(x::n_Zp)
   return x > parent(x)(div(characteristic(parent(x)), ZZ(2)))
end

show_minus_one(::Type{n_Zp}) = false

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(x::n_Zp) 
    C = parent(x)
    ptr = libSingular.n_Neg(x.ptr, C.ptr)
    return C(ptr) 
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(x::n_Zp, y::n_Zp)
   c = parent(x)
   p = libSingular.n_Add(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function -(x::n_Zp, y::n_Zp)
   c = parent(x)
   p = libSingular.n_Sub(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function *(x::n_Zp, y::n_Zp)
   c = parent(x)
   p = libSingular.n_Mult(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

+(x::n_Zp, y::Integer) = x + parent(x)(y)

+(x::Integer, y::n_Zp) = parent(y)(x) + y

-(x::n_Zp, y::Integer) = x - parent(x)(y)

-(x::Integer, y::n_Zp) = parent(y)(x) - y

*(x::n_Zp, y::Integer) = x*parent(x)(y)

*(x::Integer, y::n_Zp) = parent(y)(x)*y

+(x::n_Zp, y::n_Z) = x + parent(x)(y)

+(x::n_Z, y::n_Zp) = parent(y)(x) + y

-(x::n_Zp, y::n_Z) = x - parent(x)(y)

-(x::n_Z, y::n_Zp) = parent(y)(x) - y

*(x::n_Zp, y::n_Z) = x*parent(x)(y)

*(x::n_Z, y::n_Zp) = parent(y)(x)*y

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(x::n_Zp, y::n_Zp)
    return libSingular.n_Equal(x.ptr, y.ptr, parent(x).ptr)
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

function inv(x::n_Zp)
   c = parent(x)
   p = libSingular.n_Invers(x.ptr, c.ptr)
   return c(p)
end

function divexact(x::n_Zp, y::n_Zp)
   c = parent(x)
   p = libSingular.n_Div(x.ptr, y.ptr, c.ptr)
   return c(p)
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
   par = parent(x)
   p = libSingular.n_Gcd(x.ptr, y.ptr, par.ptr)
   return par(p)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addeq!(x::n_Zp, y::n_Zp)
   x.ptr = libSingular.n_InpAdd(x.ptr, y.ptr, parent(x).ptr)
   return x
end

function mul!(x::n_Zp, y::n_Zp, z::n_Zp)
   ptr = libSingular.n_Mult(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   return x
end

function add!(x::n_Zp, y::n_Zp, z::n_Zp)
   ptr = libSingular.n_Add(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   return x
end

function zero!(x::n_Zp)
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

promote_rule(C::Type{n_Zp}, ::Type{T}) where {T <: Integer} = n_Zp

promote_rule(C::Type{n_Zp}, ::Type{n_Z}) = n_Zp

###############################################################################
#
#   Parent call functions
#
###############################################################################

function (R::N_ZpField)()
   z = n_Zp(R)
   z.parent = R
   return z
end

function (R::N_ZpField)(x::Integer)
   z = R(libSingular.n_InitMPZ(BigInt(x), R.ptr)) 
   z.parent = R
   return z
end

function (R::N_ZpField)(n::Int)
   z = n_Zp(R, n)
   z.parent = R
   return z
end

function (R::N_ZpField)(n::n_Z)
   m = libSingular.nApplyMapFunc(R.from_n_Z, n.ptr, parent(n).ptr, R.ptr)
   z = n_Zp(R, m)
   z.parent = R
   return z
end

(R::N_ZpField)(n::n_Zp) = n

function (R::N_ZpField)(n::libSingular.number_ptr)
   z = n_Zp(R, n) 
   z.parent = R
   return z
end

function (R::N_ZpField)(x::Nemo.fmpz)
   a = BigInt()
   ccall((:flint_mpz_init_set_readonly, libflint), Nothing,
         (Ptr{BigInt}, Ptr{fmpz}), Ref(a), Ref(x))
   z = R(libSingular.n_InitMPZ(a, R.ptr))
      z.parent = R
   return z
end

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

   return N_ZpField(a)
end

function Base.Int(a::n_Zp)
  return reinterpret(Int, a.ptr.cpp_object)
end
