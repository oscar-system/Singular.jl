export n_FF, N_FField, transcendence_degree, transcendence_basis

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

elem_type(::Type{N_FField}) = n_FF

parent(a::n_FF) = a.parent

parent_type(::Type{n_FF}) = N_FField

base_ring(a::n_FF) = Union{}

base_ring(a::N_FField) = a.base_ring

@doc Markdown.doc"""
    transcendence_degree(F::N_FField)
> Return the trancendence degree of the given function field.
"""
transcendence_degree(F::N_FField) = Int(libSingular.rPar(F.ptr))

@doc Markdown.doc"""
    basis(F::N_FField)
> Return the trancendence basis of the given function field.
"""
function transcendence_basis(F::N_FField)
   n = transcendence_degree(F)
   return [F(libSingular.n_Param(Cint(i), F.ptr)) for i = 1:n]
end

@doc Markdown.doc"""
   characteristic(R::N_FField)
> Return the characteristic of the field.
"""
function characteristic(R::N_FField)
   return ZZ(libSingular.n_GetChar(R.ptr))
end

function deepcopy_internal(a::n_FF, dict::IdDict)
   return parent(a)(libSingular.n_Copy(a.ptr, parent(a).ptr))
end

#TODO: rewrite hash function
function hash(a::n_FF, h::UInt)
   return 0x2c42e12d0c837511%UInt
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

one(R::N_FField) = R(1)

zero(R::N_FField) = R(0)

function isone(n::n_FF)
   c = parent(n)
   return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_FF)
   c = parent(n)
   return libSingular.n_IsZero(n.ptr, c.ptr)
end

@doc Markdown.doc"""
   isunit(n::n_FF)
> Return `true` if $n$ is a unit in the field, i.e. nonzero.
"""
isunit(n::n_FF) = !iszero(n)

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(x::n_FF) = x

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, F::N_FField)
   print(io, "Function Field over ", base_ring(F), " with transcendence basis ", transcendence_basis(F))
end

function show(io::IO, n::n_FF)
   libSingular.StringSetS("")

   libSingular.n_Write(n.ptr, parent(n).ptr, false)

   m = libSingular.StringEndS()
   
   print(io, m)
end

needs_parentheses(x::n_FF) = false

isnegative(x::n_FF) = x == -1

show_minus_one(::Type{n_FF}) = false

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(x::n_FF) 
    C = parent(x)
    ptr = libSingular.n_Neg(x.ptr, C.ptr)
    return C(ptr) 
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(x::n_FF, y::n_FF)
   c = parent(x)
   p = libSingular.n_Add(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function -(x::n_FF, y::n_FF)
   c = parent(x)
   p = libSingular.n_Sub(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function *(x::n_FF, y::n_FF)
   c = parent(x)
   p = libSingular.n_Mult(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

+(x::n_FF, y::Integer) = x + parent(x)(y)

+(x::Integer, y::n_FF) = parent(y)(x) + y

-(x::n_FF, y::Integer) = x - parent(x)(y)

-(x::Integer, y::n_FF) = parent(y)(x) - y

*(x::n_FF, y::Integer) = x*parent(x)(y)

*(x::Integer, y::n_FF) = parent(y)(x)*y

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(x::n_FF, y::n_FF)
    return libSingular.n_Equal(x.ptr, y.ptr, parent(x).ptr)
end


isequal(x::n_FF, y::n_FF) = (x == y)

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::n_FF, y::Int)
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

function inv(x::n_FF)
   c = parent(x)
   p = libSingular.n_Invers(x.ptr, c.ptr)
   return c(p)
end

function divexact(x::n_FF, y::n_FF)
   c = parent(x)
   p = libSingular.n_Div(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   GCD and LCM
#
###############################################################################

function gcd(x::n_FF, y::n_FF)
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

function addeq!(x::n_FF, y::n_FF)
   x.ptr = libSingular.n_InpAdd(x.ptr, y.ptr, parent(x).ptr)
   return x
end

function mul!(x::n_FF, y::n_FF, z::n_FF)
   ptr = libSingular.n_Mult(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   return x
end

function add!(x::n_FF, y::n_FF, z::n_FF)
   ptr = libSingular.n_Add(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   return x
end

function zero!(x::n_FF)
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

promote_rule(C::Type{n_FF}, ::Type{T}) where T <: Integer = n_FF

promote_rule(C::Type{n_FF}, ::Type{n_Z}) = n_FF

###############################################################################
#
#   Parent call functions
#
###############################################################################

function (R::N_FField)()
   z = n_FF(R)
   z.parent = R
   return z
end

function (R::N_FField)(x::Integer)
   z = R(libSingular.n_InitMPZ(BigInt(x), R.ptr)) 
   z.parent = R
   return z
end

function (R::N_FField)(n::Int)
   z = n_FF(R, n)
   z.parent = R
   return z
end

(R::N_FField)(n::n_FF) = n

function (R::N_FField)(n::libSingular.number)
   z = n_FF(R, n) 
   z.parent = R
   return z
end

function (R::N_FField)(x::Nemo.fmpz)
   a = BigInt()
   ccall((:flint_mpz_init_set_readonly, :libflint), Nothing,
         (Ptr{BigInt}, Ptr{fmpz}), Ref(a), Ref(x))
   z = R(libSingular.n_InitMPZ(a, R.ptr))
      z.parent = R
   return z
end

###############################################################################
#
#   FunctionField constructor
#
###############################################################################

@doc Markdown.doc"""
    FunctionField(F::Field, S::Array{String, 1})
> Returns a tuple $K, a$ consisting of a function field $K$ over the field $F$
> with transcendence basis stored in the array $S$.
"""
function FunctionField(F::Field, S::Array{String, 1}; cached::Bool=true)
   typeof(F) == N_GField && error("Finite field extensions of Fp not supported.")
   l = length(S)
   isempty(l) && throw(DomainError())
   R = N_FField(F, Symbol.(S))
   return tuple(R, transcendence_basis(R))
end

@doc Markdown.doc"""
    FunctionField(F::Field, n::Int)
> Returns a tuple $K, a$ consisting of a function field $K$ over the field $F$
> with transcendence degree $n$ and transcendence basis $a1, ..., an$.
"""
function FunctionField(F::Field, n::Int; cached::Bool=true)
   S = ["a$i" for i in 1:n]
   return FunctionField(F, S)
end

