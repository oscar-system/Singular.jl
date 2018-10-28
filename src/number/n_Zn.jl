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

doc"""
    characteristic(R::N_ZnRing)
> Return the characteristic $n$ of the ring.
"""
function characteristic(R::N_ZnRing)
   return ZZ(libSingular.n_GetChar(R.ptr))
end

function deepcopy_internal(a::n_Zn, dict::ObjectIdDict)
   return parent(a)(libSingular.n_Copy(a.ptr, parent(a).ptr))
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

doc"""
    isunit(n::n_Zn)
> Return `true` if the given value is a unit in the integers modulo $n$.
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
   nn = libSingular.number_ref(n.ptr)	
   libSingular.n_Write(nn, parent(n).ptr, false)
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

promote_rule{T <: Integer}(C::Type{n_Zn}, ::Type{T}) = n_Zn

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

function (R::N_ZnRing)(n::libSingular.number)
   z = n_Zn(R, n) 
   z.parent = R
   return z
end

function (R::N_ZnRing)(x::Nemo.fmpz)
   a = BigInt()
   ccall((:flint_mpz_init_set_readonly, :libflint), Void,
         (Ptr{BigInt}, Ptr{fmpz}), &a, &x)
   z = R(libSingular.n_InitMPZ(a, R.ptr))
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
