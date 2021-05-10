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

@doc Markdown.doc"""
    characteristic(R::N_AlgExtField)

Return the characteristic of the field.
"""
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

@doc Markdown.doc"""
   isunit(n::n_algExt)
Return `true` if $n$ is a unit in the field, i.e. nonzero.
"""
isunit(n::n_algExt) = !iszero(n)

function modulus(a::N_AlgExtField)
   return a.minpoly
end

function (F::N_FField)(a::n_algExt)
   K = parent(a)
   F == parent(modulus(K)) || error("Parents must coincide")
   ptr = GC.@preserve a K F libSingular.algExt_to_transExt(a.ptr, K.ptr, F.ptr)
   return F(ptr)
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

function AbstractAlgebra.expressify(a::n_algExt; context = nothing)::Any
   F = parent(modulus(parent(a)))
   return AbstractAlgebra.expressify(F(a), context = context)
end

function Base.show(io::IO, a::n_algExt)
   print(io, AbstractAlgebra.obj_to_string(a, context = io))
end

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
      y < 0 && iszero(x) && throw(DivideError())
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
   return c(p)
end

function divexact(x::n_algExt, y::n_algExt)
   check_parent(x, y)
   c = parent(x)
   p = GC.@preserve x y c libSingular.n_Div(x.ptr, y.ptr, c.ptr)
   return c(p)
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

promote_rule(C::Type{n_algExt}, ::Type{Nemo.fmpz}) = n_algExt

promote_rule(C::Type{n_algExt}, ::Type{n_Z}) = n_algExt

###############################################################################
#
#   Parent call functions
#
###############################################################################

(R::N_AlgExtField)(n::IntegerLikeTypes = 0) = n_algExt(R, n)

function (R::N_AlgExtField)(n::n_algExt)
   R != parent(n) && error("Parent ring does not match.")
   return n
end

(R::N_AlgExtField)(n::libSingular.number_ptr) = n_algExt(R, n)

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
