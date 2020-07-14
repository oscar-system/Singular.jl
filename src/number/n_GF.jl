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

@doc Markdown.doc"""
   characteristic(R::N_GField)
> Return the characteristic of the field.
"""
function characteristic(R::N_GField)
   return ZZ(libSingular.n_GetChar(R.ptr))
end

@doc Markdown.doc"""
    degree(R::N_GField)
> Return the degree of the field as an extension of $\mathbb{F}_p$.
"""
function degree(R::N_GField)
   return R.deg
end

function deepcopy_internal(a::n_GF, dict::IdDict)
   return parent(a)(libSingular.n_Copy(a.ptr, parent(a).ptr))
end

function hash(a::n_GF, h::UInt)
   i = reinterpret(Int, a.ptr.cpp_object)
   chash = hash(characteristic(parent(a)), h)
   ihash = hash(i, h)
   return xor(xor(chash, ihash), 0x2c42e12d0c837511%UInt)
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

@doc Markdown.doc"""
   isunit(n::n_GF)
> Return `true` if $n$ is a unit in the field, i.e. nonzero.
"""
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

   libSingular.n_Write(n.ptr, parent(n).ptr, false)

   m = libSingular.StringEndS()
   
   print(io, m)
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
   x.ptr = libSingular.n_InpAdd(x.ptr, y.ptr, parent(x).ptr)
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

promote_rule(C::Type{n_GF}, ::Type{T}) where T <: Integer = n_GF

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

# not currently implemented
# function (R::N_GField)(n::n_Z)
#    m = libSingular.nApplyMapFunc(R.from_n_Z, n.ptr, parent(n).ptr, R.ptr)
#    z = n_GF(R, m)
#    z.parent = R
#    return z
# end

(R::N_GField)(n::n_GF) = n

function (R::N_GField)(n::libSingular.number_ptr)
   z = n_GF(R, n) 
   z.parent = R
   return z
end

function (R::N_GField)(x::Nemo.fmpz)
   a = BigInt()
   ccall((:flint_mpz_init_set_readonly, libflint), Nothing,
         (Ptr{BigInt}, Ptr{fmpz}), Ref(a), Ref(x))
   z = R(libSingular.n_InitMPZ(a, R.ptr))
      z.parent = R
   return z
end

###############################################################################
#
#   FiniteField constructor
#
###############################################################################

@doc Markdown.doc"""
    FiniteField(p::Int, n::Int, S::String; cached=true)
> Returns a tuple `K, a` consisting of a finite field `K` of characteristic $p$
> and degree $n$, and its generator `a`. The string used to print the
> generator is given by `S`. If the finite field is not listed in the Conway
> tables included in Singular, an error will be raised. By default, finite
> fields are cached globally, so that there is only one finite field in the
> system with given characteristic, degree and string. If this is not the
> desired behaviour, one can pass `false` for the optional `cached` parameter.
"""
function FiniteField(p::Int, n::Int, S::String; cached=true)
   n >= 16 || p >= 2^8 && throw(DomainError())
   !Nemo.isprime(Nemo.fmpz(p)) && throw(DomainError())
   n*log(p) >= 20*log(2) && throw(DomainError())
   p^n >= 2^16 && throw(DomainError())
   par = N_GField(p, n, Symbol(S))
   return par, par(libSingular.n_Param(Cint(1), par.ptr))
end

