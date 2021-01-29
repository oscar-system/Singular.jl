export n_algExt, N_AlgExtField, ResidueField

###############################################################################
#
#   Data type and parent methods
#
###############################################################################

elem_type(::Type{N_AlgExtField}) = n_algExt

parent(a::n_algExt) = a.parent

parent_type(::Type{n_algExt}) = N_AlgExtField

base_ring(a::n_algExt) = base_ring(parent(a))

base_ring(a::N_AlgExtField) = a.base_ring

@doc Markdown.doc"""
    gen(F::N_AlgExtField)

Return the generator the given algebraic extension.
"""
function gen(F::N_AlgExtField)
   return F(libSingular.n_Param(Cint(1), F.ptr))
end

@doc Markdown.doc"""
    characteristic(R::N_AlgExtField)

Return the characteristic of the field.
"""
function characteristic(R::N_AlgExtField)
   return ZZ(libSingular.n_GetChar(R.ptr))
end

function deepcopy_internal(a::n_algExt, dict::IdDict)
   return parent(a)(libSingular.n_Copy(a.ptr, parent(a).ptr))
end

function hash(a::n_algExt, h::UInt)
   error("hash not implemented")
   return 0x24f9836cd67edfd9%UInt
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
   return libSingular.n_IsOne(n.ptr, c.ptr)
end

function iszero(n::n_algExt)
   c = parent(n)
   return libSingular.n_IsZero(n.ptr, c.ptr)
end

@doc Markdown.doc"""
   isunit(n::n_algExt)
Return `true` if $n$ is a unit in the field, i.e. nonzero.
"""
isunit(n::n_algExt) = !iszero(n)

function modulus(a::N_AlgExtField)
   return a(libSingular.transExt_GetMinpoly(a.ptr))
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
   print(IOContext(io, :compact => true), "Residue field of ", base_ring(F), " modulo ", modulus(F))
end

function AbstractAlgebra.expressify(n::n_algExt; context = nothing)::Any
   # TODO this easy method might not be the best
   libSingular.StringSetS("")
   libSingular.n_Write(n.ptr, parent(n).ptr, false)
   s = libSingular.StringEndS()
   e = Meta.parse(s)
   if !isa(e, Expr)
      return e
   elseif e.head == :incomplete
      return s
   elseif e.head == :call && length(e.args) == 3 && e.args[1] == :/
      e.args[1] = ://
      return e
   else
      return e
   end
end

function show(io::IO, n::n_algExt)
   libSingular.StringSetS("")
   libSingular.n_Write(n.ptr, parent(n).ptr, false)
   m = libSingular.StringEndS()
   print(io, m)
end

###############################################################################
#
#   Unary functions
#
###############################################################################

function -(x::n_algExt)
    C = parent(x)
    ptr = libSingular.n_Neg(x.ptr, C.ptr)
    return C(ptr)
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function +(x::n_algExt, y::n_algExt)
   c = parent(x)
   p = libSingular.n_Add(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function -(x::n_algExt, y::n_algExt)
   c = parent(x)
   p = libSingular.n_Sub(x.ptr, y.ptr, c.ptr)
   return c(p)
end

function *(x::n_algExt, y::n_algExt)
   c = parent(x)
   p = libSingular.n_Mult(x.ptr, y.ptr, c.ptr)
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
    return libSingular.n_Equal(x.ptr, y.ptr, parent(x).ptr)
end


isequal(x::n_algExt, y::n_algExt) = (x == y)

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::n_algExt, y::Int)
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

function inv(x::n_algExt)
   c = parent(x)
   p = libSingular.n_Invers(x.ptr, c.ptr)
   return c(p)
end

function divexact(x::n_algExt, y::n_algExt)
   c = parent(x)
   p = libSingular.n_Div(x.ptr, y.ptr, c.ptr)
   return c(p)
end

###############################################################################
#
#   GCD and LCM
#
###############################################################################

function gcd(x::n_algExt, y::n_algExt)
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

function addeq!(x::n_algExt, y::n_algExt)
   x.ptr = libSingular.n_InpAdd(x.ptr, y.ptr, parent(x).ptr)
   return x
end

function mul!(x::n_algExt, y::n_algExt, z::n_algExt)
   ptr = libSingular.n_Mult(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   return x
end

function add!(x::n_algExt, y::n_algExt, z::n_algExt)
   ptr = libSingular.n_Add(y.ptr, z.ptr, parent(x).ptr)
   libSingular.n_Delete(x.ptr, parent(x).ptr)
   x.ptr = ptr
   return x
end

function zero!(x::n_algExt)
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
    ResidueField(F::Singular.N_FField, a::n_transExt)

Return the residue field F/(a).
"""
function ResidueField(F::Singular.N_FField, a::n_transExt)
   return N_AlgExtField(F, a)
end

