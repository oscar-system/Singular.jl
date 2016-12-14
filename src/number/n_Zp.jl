export SingularFp

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

elem_type(::SingularN_ZpField) = n_Zp

parent(a::n_Zp) = a.parent

parent_type(::Type{n_Zp}) = SingularN_ZpField

base_ring(a::n_Zp) = Union{}

base_ring(a::SingularN_ZpField) = Union{}

function characteristic(R::SingularN_ZpField)
   return SingularZZ(libSingular.n_GetChar(R.ptr))
end

function deepcopy(a::n_Zp)
   return parent(a)(libSingular.n_Copy(a.ptr, parent(a).ptr))
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

one(R::SingularN_ZpField) = R(1)

zero(R::SingularN_ZpField) = R(0)

function isone(n::n_Zp)
   c = parent(n)
   return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_Zp)
   c = parent(n)
   return libSingular.n_IsZero(n.ptr, c.ptr)
end

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

function show(io::IO, c::SingularN_ZpField)
   print(io, "Finite Field of Characteristic ", characteristic(c))
end

function show(io::IO, n::n_Zp)
   libSingular.StringSetS("")

   nn = libSingular.number_ref(n.ptr)	
   libSingular.n_Write(nn, parent(n).ptr, false)
   n.ptr = nn[]

   m = libSingular.StringEndS()
   s = unsafe_string(m) 
   libSingular.omFree(Ptr{Void}(m))

   print(io, s)
end

needs_parentheses(x::n_Zp) = false

function is_negative(x::n_Zp)
   return x > parent(x)(div(characteristic(parent(x)), SingularZZ(2)))
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

function isless(x::n_Zp, y::n_Zp)
    libSingular.n_Greater(y.ptr, x.ptr, parent(x).ptr)
end

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

function div(x::n_Zp, y::n_Zp)
   c = parent(x)
   p = libSingular.n_Div(x.ptr, y.ptr, c.ptr)
   return c(p)
end

divexact(x::n_Zp, y::n_Zp) = div(x, y)

###############################################################################
#
#   Euclidean division
#
###############################################################################

function divrem(x::n_Zp, y::n_Zp)
   par = parent(x)
   r = [libSingular.n_Init(0, par.ptr)]
   q = libSingular.n_QuotRem(x.ptr, y.ptr, pointer(r), par.ptr)
   return par(q), par(r[])
end

function rem(x::n_Zp, y::n_Zp)
   return parent(x)()
end

mod(x::n_Zp, y::n_Zp) = rem(x, y)

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

function lcm(x::n_Zp, y::n_Zp)
   if x == 0 && y == 0
      return zero(parent(x))
   end
   return div(x*y, gcd(x, y))
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addeq!(x::n_Zp, y::n_Zp)
    xx = libSingular.number_ref(x.ptr)
    libSingular.n_InpAdd(xx, y.ptr, parent(x).ptr)
    x.ptr = xx[]
    nothing
end

function mul!(x::n_Zp, y::n_Zp, z::n_Zp)
   ptr = libSingular.n_Mult(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   nothing
end

function add!(x::n_Zp, y::n_Zp, z::n_Zp)
   ptr = libSingular.n_Add(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   nothing
end

function zero!(x::n_Zp)
   ptr = libSingular.n_Init(0, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   nothing
end

###############################################################################
#
#   Conversions and promotions
#
###############################################################################

Base.promote_rule{T <: Integer}(C::Type{n_Zp}, ::Type{T}) = n_Zp

Base.promote_rule(C::Type{n_Zp}, ::Type{n_Z}) = n_Zp

###############################################################################
#
#   Parent call functions
#
###############################################################################

function (R::SingularN_ZpField)()
   z = n_Zp(R)
   z.parent = R
   return z
end

function (R::SingularN_ZpField)(x::Integer)
   z = R(libSingular.n_InitMPZ(BigInt(x), R.ptr)) 
   z.parent = R
   return z
end

function (R::SingularN_ZpField)(n::Int)
   z = n_Zp(R, n)
   z.parent = R
   return z
end

function (R::SingularN_ZpField)(n::n_Z)
   m = libSingular.nApplyMapFunc(R.from_n_Z, n.ptr, parent(n).ptr, R.ptr)
   z = n_Zp(R, m)
   z.parent = R
   return z
end

(R::SingularN_ZpField)(n::n_Zp) = n

function (R::SingularN_ZpField)(n::libSingular.number)
   z = n_Zp(R, n) 
   z.parent = R
   return z
end

function (R::SingularN_ZpField)(x::Nemo.fmpz)
   a = BigInt()
   ccall((:flint_mpz_init_set_readonly, :libflint), Void,
         (Ptr{BigInt}, Ptr{fmpz}), &a, &x)
   z = R(libSingular.n_InitMPZ(a, R.ptr))
      z.parent = R
   return z
end

###############################################################################
#
#   ResidueRing constructor
#
###############################################################################

function SingularFp(a::Int; cached=true)
   a == 0 && throw(DivideError())
   a < 0 && throw(DomainError())
   !Nemo.is_prime(UInt(a)) && throw(DomainError())

   return SingularN_ZpField(a)
end
