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

characteristic(R::N_GField) = ZZ(_characteristic(R))

_characteristic(R::N_GField) = Int(libSingular.n_GetChar(R.ptr))

@doc Markdown.doc"""
    degree(R::N_GField)

Return the degree of the field as an extension of $\mathbb{F}_p$.
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

gen(R::N_GField) = R(libSingular.n_Param(Cint(1), R.ptr))

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
Return `true` if $n$ is a unit in the field, i.e. nonzero.
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

function Base.show(io::IO, ::MIME"text/plain", a::n_GF)
  print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

function AbstractAlgebra.expressify(a::n_GF; context = nothing)::Any
  F = parent(a)
  i = reinterpret(Int, a.ptr.cpp_object)
  if 1 < i < characteristic(F)^degree(F)
    return Expr(:call, :^, F.S, i)
  elseif i == 1
    return F.S
  elseif i == 0
    return 1
  else
    return 0
  end
end

function show(io::IO, n::n_GF)
   libSingular.StringSetS("")
   libSingular.n_Write(n.ptr, parent(n).ptr, false)
   m = libSingular.StringEndS()
   print(io, m)
end

isnegative(x::n_GF) = false

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
    y < 0 && throw(DomainError(y, "exponent must be non-negative"))
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
#   Random functions
#
###############################################################################

Random.Sampler(::Type{RNG}, K::N_GField, n::Random.Repetition) where {RNG<:AbstractRNG} =
   SamplerSimple(K, Random.Sampler(RNG, 0:_characteristic(K) - 1, Val(Inf)))

function rand(rng::AbstractRNG, Ksp::SamplerSimple{N_GField})
   K = Ksp[]
   r = degree(K)
   alpha = gen(K)
   res = zero(K)
   for i = 0 : (r-1)
      c = rand(rng, Ksp.data)
      res += c * alpha^i
   end
   return res
end

Random.gentype(::Type{N_GField}) = elem_type(N_GField)

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

(R::N_GField)(n::IntegerLikeTypes = 0) = n_GF(R, n)

(R::N_GField)(n::n_GF) = n

(R::N_GField)(n::libSingular.number_ptr) = n_GF(R, n)

###############################################################################
#
#   FiniteField constructor
#
###############################################################################

@doc Markdown.doc"""
    FiniteField(p::Int, n::Int, S::String; cached=true)

Returns a tuple `K, a` consisting of a finite field `K` of characteristic $p$
and degree $n$, and its generator `a`. The string used to print the
generator is given by `S`. If the finite field is not listed in the Conway
tables included in Singular, an error will be raised. By default, finite
fields are cached globally, so that there is only one finite field in the
system with given characteristic, degree and string. If this is not the
desired behaviour, one can pass `false` for the optional `cached` parameter.
"""
function FiniteField(p::Int, n::Int, S::String; cached=true)
   p >= 2^8 && throw(DomainError(p, "p must be < 256"))
   !Nemo.isprime(Nemo.fmpz(p)) && throw(DomainError(p, "p must be prime"))
   (n*log(p) >= 20*log(2) || p^n >= 2^16) &&
      throw(DomainError("p^n must be < 2^16"))
   par = N_GField(p, n, Symbol(S))
   return par, gen(par)
end
