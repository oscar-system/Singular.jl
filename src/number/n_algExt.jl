export n_algExt, N_AlgExtField, AlgebraicExtensionField, modulus

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

elem_type(::Type{N_AlgExtField}) = n_algExt

parent(a::n_algExt) = a.parent

parent_type(::Type{n_algExt}) = N_AlgExtField

base_ring(a::n_algExt) = base_ring(parent(a))

base_ring(a::N_AlgExtField) = base_ring(a.minpoly)

@doc Markdown.doc"""
    gen(F::N_AlgExtField)

Return the generator the given algebraic extension.
"""
function gen(F::N_AlgExtField)
   GC.@preserve F return F(libSingular.n_Param(Cint(1), F.ptr))
end

function characteristic(R::N_AlgExtField)
   GC.@preserve R return ZZ(libSingular.n_GetChar(R.ptr))
end

function deepcopy_internal(a::n_algExt, dict::IdDict)
   R = parent(a)
   GC.@preserve a R return R(libSingular.n_Copy(a.ptr, R.ptr))
end

function hash(a::n_algExt, h::UInt)
   K = parent(a)
   F = parent(modulus(K))
   phash = hash(F(a))
   chash = hash(characteristic(K), h)
   return xor(xor(chash, phash), 0x2c42e12d0c837511%UInt)
end

function check_parent(x::n_algExt, y::n_algExt)
   parent(x) != parent(y) && error("Parents must coincide")
   nothing
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

one(R::N_AlgExtField) = R(1)

zero(R::N_AlgExtField) = R(0)

function isone(n::n_algExt)
   c = parent(n)
   GC.@preserve n c return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_algExt)
   c = parent(n)
   GC.@preserve n c return libSingular.n_IsZero(n.ptr, c.ptr)
end

isunit(n::n_algExt) = !iszero(n)

function modulus(a::N_AlgExtField)
   return a.minpoly
end

###############################################################################
#
#   Canonicalisation
#
###############################################################################

canonical_unit(x::n_algExt) = x

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, F::N_AlgExtField)
   print(IOContext(io, :compact => true), "Algebraic Extension of ",
         base_ring(F), " with defining equation ", modulus(F))
end

function expressify(a::n_algExt; context = nothing)::Any
   F = parent(modulus(parent(a)))
   return expressify(F(a), context = context)
end

AbstractAlgebra.@enable_all_show_via_expressify n_algExt

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(x::n_algExt)
    c = parent(x)
    p = GC.@preserve x c libSingular.n_Neg(x.ptr, c.ptr)
    return c(p)
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(x::n_algExt, y::n_algExt)
   check_parent(x, y)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Add(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function -(x::n_algExt, y::n_algExt)
   check_parent(x, y)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Sub(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function *(x::n_algExt, y::n_algExt)
   check_parent(x, y)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Mult(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Ad hoc arithmetic functions
#
###############################################################################

+(x::n_algExt, y::Integer) = x + parent(x)(y)

+(x::Integer, y::n_algExt) = parent(y)(x) + y

-(x::n_algExt, y::Integer) = x - parent(x)(y)

-(x::Integer, y::n_algExt) = parent(y)(x) - y

*(x::n_algExt, y::Integer) = x*parent(x)(y)

*(x::Integer, y::n_algExt) = parent(y)(x)*y

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(x::n_algExt, y::n_algExt)
   check_parent(x, y)
   c = parent(x)
   GC.@preserve x y c return libSingular.n_Equal(x.ptr, y.ptr, c.ptr)
end

isequal(x::n_algExt, y::n_algExt) = (x == y)

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::n_algExt, y::Int)
   if isone(x)
      return x
   elseif y == 0
      return one(parent(x))
   elseif y == 1
      return x
   else
      # low priority FIXME Singular still crashes upon division by zero here
      y < 0 && iszero(x) && error("division by zero")
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

function inv(x::n_algExt)
   c = parent(x)
   p = GC.@preserve x c libSingular.n_Invers(x.ptr, c.ptr)
   z = c(p)
   libSingular.check_error()
   return z
end

function divexact(x::n_algExt, y::n_algExt; check::Bool=true)
   check_parent(x, y)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Div(x.ptr, y.ptr, c.ptr)
   z = c(p)
   libSingular.check_error()
   return z
end

###############################################################################
#
#   GCD and LCM
#
###############################################################################

function gcd(x::n_algExt, y::n_algExt)
   check_parent(x, y)
   if iszero(x) && iszero(y)
      return zero(parent(x))
   end
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Gcd(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   Unsafe functions
#
###############################################################################

function addeq!(x::n_algExt, y::n_algExt)
   R = parent(x)
   x.ptr = GC.@preserve x y R libSingular.n_InpAdd(x.ptr, y.ptr, R.ptr)
   return x
end

function mul!(x::n_algExt, y::n_algExt, z::n_algExt)
   R = parent(x)
   GC.@preserve x y z R begin
      ptr = libSingular.n_Mult(y.ptr, z.ptr, R.ptr)
      libSingular.n_Delete(x.ptr, R.ptr)
      x.ptr = ptr
      return x
   end
end

function add!(x::n_algExt, y::n_algExt, z::n_algExt)
   R = parent(x)
   GC.@preserve x y z R begin
      ptr = libSingular.n_Add(y.ptr, z.ptr, R.ptr)
      libSingular.n_Delete(x.ptr, R.ptr)
      x.ptr = ptr
      return x
   end
end

function zero!(x::n_algExt)
   R = parent(x)
   GC.@preserve x R begin
      ptr = libSingular.n_Init(0, R.ptr)
      libSingular.n_Delete(x.ptr, R.ptr)
      x.ptr = ptr
      return x
   end
end

###############################################################################
#
#   Conversions and promotions
#
###############################################################################

promote_rule(C::Type{n_algExt}, ::Type{T}) where T <: Integer = n_algExt

promote_rule(C::Type{n_algExt}, ::Type{T}) where T <: Rational = n_algExt

promote_rule(C::Type{n_algExt}, ::Type{Nemo.fmpz}) = n_algExt

promote_rule(C::Type{n_algExt}, ::Type{Nemo.fmpq}) = n_algExt

promote_rule(C::Type{n_algExt}, ::Type{n_Z}) = n_algExt

promote_rule(C::Type{n_algExt}, ::Type{n_Q}) = n_algExt

# TODO really need a hand-crafted nf_elem <-> n_algExt
function (SK::Singular.N_AlgExtField)(a::Singular.Nemo.nf_elem)
  K = parent(a)
  SKa = gen(SK)
  res = SK(coeff(a, 0))
  for i in 1:degree(K)-1
    res += SK(coeff(a, i))*SKa^i
  end
  return res
end

# this is going to be dreadfully slow
function (K::Singular.Nemo.AnticNumberField)(a::Singular.n_algExt)
  SK = parent(a)
  SF = parent(modulus(SK))
  # it gets even worse: n_transExt_to_spoly only converts the "numerator"
  #                     which doesn't include the "denominator" we need
  SFa = SF(a)
  numSa = n_transExt_to_spoly(numerator(SFa))
  denSa = first(coefficients(n_transExt_to_spoly(denominator(SFa))))
  res = zero(K)
  Ka = gen(K)
  for (c, e) in zip(coefficients(numSa), exponent_vectors(numSa))
    res += Singular.Nemo.fmpq(c//denSa)*Ka^e[1]
  end
  return res
end

function (F::N_FField)(a::n_algExt)
   K = parent(a)
   F == parent(modulus(K)) || error("Parents must coincide")
   ptr = GC.@preserve a K F libSingular.algExt_to_transExt(a.ptr, K.ptr, F.ptr)
   return F(ptr)
end

###############################################################################
#
#   Parent call functions
#
###############################################################################

function (K::N_AlgExtField)(a::n_transExt)
   F = parent(a)
   F == parent(modulus(K)) || error("Parents must coincide")
   ptr = GC.@preserve a K F libSingular.transExt_to_algExt(a.ptr, K.ptr, F.ptr)
   return K(ptr)
end

function (K::N_AlgExtField)(a::n_algExt)
   K == parent(a) || error("Parents must coincide")
   return a
end

function (K::N_AlgExtField)(a::IntegerLikeTypes = 0)
  return n_algExt(K, a)
end

function (K::Singular.N_AlgExtField)(a::Union{n_Q, Nemo.fmpq, Rational})
  return K(numerator(a))//K(denominator(a))
end

# take ownership of the pointer - not for general users
function (K::N_AlgExtField)(n::libSingular.number_ptr)
  return n_algExt(K, n)
end

###############################################################################
#
#   AlgebraicExtensionField constructor
#
###############################################################################

@doc Markdown.doc"""
    AlgebraicExtensionField(F::Singular.N_FField, a::n_transExt)

Given a function field F = R(x) in one variable and a defining polynomial a
in R[x], return a tuple consisting of the new field R[x]/(a) and the class of x.
"""
function AlgebraicExtensionField(F::Singular.N_FField, a::n_transExt;
                                                                 cached = true)
   parent(a) == F || error("Parents must coincide")
   transcendence_degree(F) == 1 || throw(ArgumentError("Only algebraic extensions in one variable are supported."))
   r = N_AlgExtField(F, a, cached)
   return (r, gen(r))
end
