export SingularFiniteField

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

elem_type(::SingularN_GFField) = n_GF

parent(a::n_GF) = a.parent

parent_type(::Type{n_GF}) = SingularN_GFField

base_ring(a::n_GF) = Union{}

base_ring(a::SingularN_GFField) = Union{}

function characteristic(R::SingularN_GFField)
   return SingularZZ(libSingular.n_GetChar(R.ptr))
end

function degree(R::SingularN_GFField)
   return R.deg
end

function deepcopy(a::n_GF)
   return parent(a)(libSingular.n_Copy(a.ptr, parent(a).ptr))
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

one(R::SingularN_GFField) = R(1)

zero(R::SingularN_GFField) = R(0)

function isone(n::n_GF)
   c = parent(n)
   return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_GF)
   c = parent(n)
   return libSingular.n_IsZero(n.ptr, c.ptr)
end

isunit(n::n_GF) = !iszero(n)

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(x::n_GF) = x

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, c::SingularN_GFField)
   print(io, "Finite Field of Characteristic ", characteristic(c), " and degree ", degree(c))
end

function show(io::IO, n::n_GF)
   libSingular.StringSetS("")

   nn = libSingular.number_ref(n.ptr)	
   libSingular.n_Write(nn, parent(n).ptr, false)
   n.ptr = nn[]

   m = libSingular.StringEndS()
   s = unsafe_string(m) 
   libSingular.omFree(Ptr{Void}(m))

   print(io, s)
end

needs_parentheses(x::n_GF) = false

isnegative(x::n_GF) = x == -1

show_minus_one(::Type{n_GF}) = false

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(x::n_GF) 
    C = parent(x)
    ptr = libSingular.n_Neg(x.ptr, C.ptr)
    return C(ptr) 
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(x::n_GF, y::n_GF)
   c = parent(x)
   p = libSingular.n_Add(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function -(x::n_GF, y::n_GF)
   c = parent(x)
   p = libSingular.n_Sub(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function *(x::n_GF, y::n_GF)
   c = parent(x)
   p = libSingular.n_Mult(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

+(x::n_GF, y::Integer) = x + parent(x)(y)

+(x::Integer, y::n_GF) = parent(y)(x) + y

-(x::n_GF, y::Integer) = x - parent(x)(y)

-(x::Integer, y::n_GF) = parent(y)(x) - y

*(x::n_GF, y::Integer) = x*parent(x)(y)

*(x::Integer, y::n_GF) = parent(y)(x)*y

+(x::n_GF, y::n_Z) = x + parent(x)(y)

+(x::n_Z, y::n_GF) = parent(y)(x) + y

-(x::n_GF, y::n_Z) = x - parent(x)(y)

-(x::n_Z, y::n_GF) = parent(y)(x) - y

*(x::n_GF, y::n_Z) = x*parent(x)(y)

*(x::n_Z, y::n_GF) = parent(y)(x)*y

###############################################################################
#
#   Comparison
#
###############################################################################

function isless(x::n_GF, y::n_GF)
    libSingular.n_Greater(y.ptr, x.ptr, parent(x).ptr)
end

function ==(x::n_GF, y::n_GF)
    return libSingular.n_Equal(x.ptr, y.ptr, parent(x).ptr)
end


isequal(x::n_GF, y::n_GF) = (x == y)

###############################################################################
#
#   Ad hoc comparison
#
###############################################################################

==(x::n_GF, y::Integer) = (x ==  parent(x)(y))

==(x::Integer, y::n_GF) = (parent(y)(x) == y)

==(x::n_GF, y::n_Z) = (x ==  parent(x)(y))

==(x::n_Z, y::n_GF) = (parent(y)(x) == y)

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::n_GF, y::Int)
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

function inv(x::n_GF)
   c = parent(x)
   p = libSingular.n_Invers(x.ptr, c.ptr)
   return c(p)
end

function div(x::n_GF, y::n_GF)
   c = parent(x)
   p = libSingular.n_Div(x.ptr, y.ptr, c.ptr)
   return c(p)
end

divexact(x::n_GF, y::n_GF) = div(x, y)

###############################################################################
#
#   Euclidean division
#
###############################################################################

function divrem(x::n_GF, y::n_GF)
   par = parent(x)
   r = [libSingular.n_Init(0, par.ptr)]
   q = libSingular.n_QuotRem(x.ptr, y.ptr, pointer(r), par.ptr)
   return par(q), par(r[])
end

function rem(x::n_GF, y::n_GF)
   return parent(x)()
end

mod(x::n_GF, y::n_GF) = rem(x, y)

###############################################################################
#
#   GCD and LCM
#
###############################################################################

function gcd(x::n_GF, y::n_GF)
   if x == 0 && y == 0
      return zero(parent(x))
   end
   par = parent(x)
   p = libSingular.n_Gcd(x.ptr, y.ptr, par.ptr)
   return par(p)
end

function lcm(x::n_GF, y::n_GF)
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

function addeq!(x::n_GF, y::n_GF)
    xx = libSingular.number_ref(x.ptr)
    libSingular.n_InpAdd(xx, y.ptr, parent(x).ptr)
    x.ptr = xx[]
    nothing
end

function mul!(x::n_GF, y::n_GF, z::n_GF)
   ptr = libSingular.n_Mult(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   nothing
end

function add!(x::n_GF, y::n_GF, z::n_GF)
   ptr = libSingular.n_Add(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   nothing
end

function zero!(x::n_GF)
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

promote_rule{T <: Integer}(C::Type{n_GF}, ::Type{T}) = n_GF

promote_rule(C::Type{n_GF}, ::Type{n_Z}) = n_GF

###############################################################################
#
#   Parent call functions
#
###############################################################################

function (R::SingularN_GFField)()
   z = n_GF(R)
   z.parent = R
   return z
end

function (R::SingularN_GFField)(x::Integer)
   z = R(libSingular.n_InitMPZ(BigInt(x), R.ptr)) 
   z.parent = R
   return z
end

function (R::SingularN_GFField)(n::Int)
   z = n_GF(R, n)
   z.parent = R
   return z
end

function (R::SingularN_GFField)(n::n_Z)
   m = libSingular.nApplyMapFunc(R.from_n_Z, n.ptr, parent(n).ptr, R.ptr)
   z = n_GF(R, m)
   z.parent = R
   return z
end

(R::SingularN_GFField)(n::n_GF) = n

function (R::SingularN_GFField)(n::libSingular.number)
   z = n_GF(R, n) 
   z.parent = R
   return z
end

function (R::SingularN_GFField)(x::Nemo.fmpz)
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

function SingularFiniteField(p::Int, n::Int, S::String; cached=true)
   n >= 16 || p >= 2^16 && throw(DomainError())
   !Nemo.isprime(UInt(p)) && throw(DomainError())
   n*log(p) >= 20*log(2) && throw(DomainError())
   p^n >= 2^16 && throw(DomainError())
   n == 1 && error("Degree one finite fields not supported. Please use SingularFp")
   par = SingularN_GFField(p, n, symbol(S))
   return par, par(libSingular.n_Param(Cint(1), par.ptr))
end
