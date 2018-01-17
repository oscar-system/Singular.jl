export n_GF, N_GField

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

elem_type(::Type{N_GField}) = n_GF

parent(a::n_GF) = a.parent

parent_type(::Type{n_GF}) = N_GField

base_ring(a::n_GF) = Union{}

base_ring(a::N_GField) = Union{}

function characteristic(R::N_GField)
   return ZZ(libSingular.n_GetChar(R.ptr))
end

function degree(R::N_GField)
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

one(R::N_GField) = R(1)

zero(R::N_GField) = R(0)

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

function show(io::IO, c::N_GField)
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

function divexact(x::n_GF, y::n_GF)
   c = parent(x)
   p = libSingular.n_Div(x.ptr, y.ptr, c.ptr)
   return c(p)
end

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

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addeq!(x::n_GF, y::n_GF)
    xx = libSingular.number_ref(x.ptr)
    libSingular.n_InpAdd(xx, y.ptr, parent(x).ptr)
    x.ptr = xx[]
    return x
end

function mul!(x::n_GF, y::n_GF, z::n_GF)
   ptr = libSingular.n_Mult(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   return x
end

function add!(x::n_GF, y::n_GF, z::n_GF)
   ptr = libSingular.n_Add(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   return x
end

function zero!(x::n_GF)
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

promote_rule{T <: Integer}(C::Type{n_GF}, ::Type{T}) = n_GF

promote_rule(C::Type{n_GF}, ::Type{n_Z}) = n_GF

###############################################################################
#
#   Parent call functions
#
###############################################################################

function (R::N_GField)()
   z = n_GF(R)
   z.parent = R
   return z
end

function (R::N_GField)(x::Integer)
   z = R(libSingular.n_InitMPZ(BigInt(x), R.ptr)) 
   z.parent = R
   return z
end

function (R::N_GField)(n::Int)
   z = n_GF(R, n)
   z.parent = R
   return z
end

function (R::N_GField)(n::n_Z)
   m = libSingular.nApplyMapFunc(R.from_n_Z, n.ptr, parent(n).ptr, R.ptr)
   z = n_GF(R, m)
   z.parent = R
   return z
end

(R::N_GField)(n::n_GF) = n

function (R::N_GField)(n::libSingular.number)
   z = n_GF(R, n) 
   z.parent = R
   return z
end

function (R::N_GField)(x::Nemo.fmpz)
   a = BigInt()
   ccall((:flint_mpz_init_set_readonly, :libflint), Void,
         (Ptr{BigInt}, Ptr{fmpz}), &a, &x)
   z = R(libSingular.n_InitMPZ(a, R.ptr))
      z.parent = R
   return z
end

###############################################################################
#
#   FiniteField constructor
#
###############################################################################

function FiniteField(p::Int, n::Int, S::String; cached=true)
   n >= 16 || p >= 2^16 && throw(DomainError())
   !Nemo.isprime(Nemo.fmpz(p)) && throw(DomainError())
   n*log(p) >= 20*log(2) && throw(DomainError())
   p^n >= 2^16 && throw(DomainError())
   n == 1 && error("Degree one finite fields not supported. Please use SingularFp")
   par = N_GField(p, n, Symbol(S))
   return par, par(libSingular.n_Param(Cint(1), par.ptr))
end
