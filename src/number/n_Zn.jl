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

@doc Markdown.doc"""
    characteristic(R::N_ZnRing)

Return the characteristic $n$ of the ring.
"""
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
   return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_Zn)
   c = parent(n)
   return libSingular.n_IsZero(n.ptr, c.ptr)
end

@doc Markdown.doc"""
    isunit(n::n_Zn)

Return `true` if the given value is a unit in the integers modulo $n$.
"""
isunit(n::n_Zn) = gcd(n, parent(n)(characteristic(parent(n)))) == 1

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

function show(io::IO, n::n_Zn)
   libSingular.StringSetS("")
   libSingular.n_Write(n.ptr, parent(n).ptr, false)
   m = libSingular.StringEndS()
   print(io, m)
end

needs_parentheses(x::n_Zn) = false

isnegative(x::n_Zn) = false

show_minus_one(::Type{n_Zn}) = true

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(x::n_Zn) 
    C = parent(x)
    ptr = libSingular.n_Neg(x.ptr, C.ptr)
    return C(ptr) 
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(x::n_Zn, y::n_Zn)
   c = parent(x)
   p = libSingular.n_Add(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function -(x::n_Zn, y::n_Zn)
   c = parent(x)
   p = libSingular.n_Sub(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function *(x::n_Zn, y::n_Zn)
   c = parent(x)
   p = libSingular.n_Mult(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

+(x::n_Zn, y::Integer) = x + parent(x)(y)

+(x::Integer, y::n_Zn) = parent(y)(x) + y

-(x::n_Zn, y::Integer) = x - parent(x)(y)

-(x::Integer, y::n_Zn) = parent(y)(x) - y

*(x::n_Zn, y::Integer) = x*parent(x)(y)

*(x::Integer, y::n_Zn) = parent(y)(x)*y

+(x::n_Zn, y::n_Z) = x + parent(x)(y)

+(x::n_Z, y::n_Zn) = parent(y)(x) + y

-(x::n_Zn, y::n_Z) = x - parent(x)(y)

-(x::n_Z, y::n_Zn) = parent(y)(x) - y

*(x::n_Zn, y::n_Z) = x*parent(x)(y)

*(x::n_Z, y::n_Zn) = parent(y)(x)*y

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(x::n_Zn, y::n_Zn)
    return libSingular.n_Equal(x.ptr, y.ptr, parent(x).ptr)
end

isequal(x::n_Zn, y::n_Zn) = (x == y)

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

function divexact(x::n_Zn, y::n_Zn)
   c = parent(x)
   p = libSingular.n_ExactDiv(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   GCD and LCM
#
###############################################################################

function gcd(x::n_Zn, y::n_Zn)
   par = parent(x)
   p = libSingular.n_Gcd(x.ptr, y.ptr, par.ptr)
   return par(p)
end

###############################################################################
#
#   Extended GCD
#
###############################################################################

function gcdx(x::n_Zn, y::n_Zn)
   par = parent(x)
   s = [libSingular.n_Init(0, par.ptr)]
   t = [libSingular.n_Init(0, par.ptr)]
   g = libSingular.n_ExtGcd(x.ptr, y.ptr, pointer(s), pointer(t), par.ptr)
   return par(g), par(s[]), par(t[])
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addeq!(x::n_Zn, y::n_Zn)
   x.ptr = libSingular.n_InpAdd(x.ptr, y.ptr, parent(x).ptr)
   return x
end

function mul!(x::n_Zn, y::n_Zn, z::n_Zn)
   ptr = libSingular.n_Mult(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   return x
end

function add!(x::n_Zn, y::n_Zn, z::n_Zn)
   ptr = libSingular.n_Add(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   return x
end

function zero!(x::n_Zn)
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

promote_rule(C::Type{n_Zn}, ::Type{T}) where {T <: Integer} = n_Zn

promote_rule(C::Type{n_Zn}, ::Type{n_Z}) = n_Zn

###############################################################################
#
#   Parent call functions
#
###############################################################################

function (R::N_ZnRing)()
   z = n_Zn(R)
   z.parent = R
   return z
end

function (R::N_ZnRing)(x::Integer)
   z = R(libSingular.n_InitMPZ(BigInt(x), R.ptr)) 
   z.parent = R
   return z
end

function (R::N_ZnRing)(n::Int)
   z = n_Zn(R, n)
   z.parent = R
   return z
end

function (R::N_ZnRing)(n::n_Z)
   m = libSingular.nApplyMapFunc(R.from_n_Z, n.ptr, parent(n).ptr, R.ptr)
   z = n_Zn(R, m)
   z.parent = R
   return z
end

(R::N_ZnRing)(n::n_Zn) = n

function (R::N_ZnRing)(n::libSingular.number_ptr)
   z = n_Zn(R, n)
   z.parent = R
   return z
end

function (R::N_ZnRing)(x::Nemo.fmpz)
   z = convert_from_fmpz(R, x)
   z.parent = R
   return z
end

###############################################################################
#
#   SingularResidueRing constructor
#
###############################################################################

function ResidueRing(R::Integers, a::Int; cached=true)
   a == 0 && throw(DivideError())
   a < 0 && throw(DomainError())

   return N_ZnRing(a)
end
