###############################################################################
#
#   Data type and parent methods
#
###############################################################################

elem_type(::Type{N_UnknownSingularCoefficientRing}) = n_unknownsingularcoefficient

parent(a::n_unknownsingularcoefficient) = a.parent

parent_type(::Type{n_unknownsingularcoefficient}) = N_UnknownSingularCoefficientRing

base_ring(a::n_unknownsingularcoefficient) = Union{}

base_ring(R::N_UnknownSingularCoefficientRing) = Union{}

function characteristic(R::N_UnknownSingularCoefficientRing)
   return ZZ(_characteristic(R))
end

function _characteristic(R::N_UnknownSingularCoefficientRing)
   GC.@preserve R Int(libSingular.n_GetChar(R.ptr))
end

function check_parent(a::n_unknownsingularcoefficient, b::n_unknownsingularcoefficient)
   parent(a) != parent(b) && error("Incompatible parents")
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

function canonical_unit(a::n_unknownsingularcoefficient)
   error("canonical_unit not implemented by singular")
end

function deepcopy_internal(a::n_unknownsingularcoefficient, dict::IdDict)
   R = parent(a)
   GC.@preserve a R return R(libSingular.n_Copy(a.ptr, R.ptr))
end

function one(R::N_UnknownSingularCoefficientRing)
   GC.@preserve R return R(libSingular.n_Init(1, R.ptr))
end

function zero(R::N_UnknownSingularCoefficientRing)
   GC.@preserve R return R(libSingular.n_Init(0, R.ptr))
end

one(a::n_unknownsingularcoefficient) = one(parent(a))

zero(a::n_unknownsingularcoefficient) = zero(parent(a))

function isone(a::n_unknownsingularcoefficient)
   R = parent(a)
   GC.@preserve a R return Bool(libSingular.n_IsOne(a.ptr, R.ptr))
end

function iszero(a::n_unknownsingularcoefficient)
   R = parent(a)
   GC.@preserve a R return Bool(libSingular.n_IsZero(a.ptr, R.ptr))
end

function hash(a::n_unknownsingularcoefficient, h::UInt)
   error("hash not implemented by singular")
end

###############################################################################
#
#   String I/O
#
###############################################################################

function AbstractAlgebra.expressify(a::n_unknownsingularcoefficient; context = nothing)::Any
   # TODO this easy method might not be the best
   libSingular.StringSetS("")
   R = parent(a)
   GC.@preserve a R libSingular.n_Write(a.ptr, R.ptr, false)
   s = libSingular.StringEndS()
   e = Meta.parse(s)
   if isa(e, Expr) && e.head == :incomplete
      return s
   else
      return e
   end
end

function Base.show(io::IO, ::MIME"text/plain", a::n_unknownsingularcoefficient)
   libSingular.StringSetS("")
   R = parent(a)
   GC.@preserve a R libSingular.n_Write(a.ptr, R.ptr, false)
   print(io, libSingular.StringEndS())
end

function show(io::IO, a::n_unknownsingularcoefficient)
   libSingular.StringSetS("")
   R = parent(a)
   GC.@preserve a R libSingular.n_Write(a.ptr, R.ptr, false)
   print(io, libSingular.StringEndS())
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::n_unknownsingularcoefficient)
   R = parent(a)
   n = GC.@preserve a b R libSingular.n_Neg(a.ptr, R.ptr)
   return R(n)
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::n_unknownsingularcoefficient, b::n_unknownsingularcoefficient)
   check_parent(a, b)
   R = parent(a)
   n = GC.@preserve a b R libSingular.n_Add(a.ptr, b.ptr, R.ptr)
   return R(n)
end

function -(a::n_unknownsingularcoefficient, b::n_unknownsingularcoefficient)
   check_parent(a, b)
   R = parent(a)
   n = GC.@preserve a b R libSingular.n_Sub(a.ptr, b.ptr, R.ptr)
   return R(n)
end

function *(a::n_unknownsingularcoefficient, b::n_unknownsingularcoefficient)
   check_parent(a, b)
   R = parent(a)
   n = GC.@preserve a b R libSingular.n_Mult(a.ptr, b.ptr, R.ptr)
   return R(n)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(a::n_unknownsingularcoefficient, b::n_unknownsingularcoefficient)
   check_parent(a, b)
   R = parent(a)
   GC.@preserve a b R return libSingular.n_Equal(a.ptr, b.ptr, R.ptr)
end

###############################################################################
#
#   Powering
#
###############################################################################

function ^(x::n_unknownsingularcoefficient, y::Int)
   y < 0 && throw(DomainError(y, "exponent must be non-negative"))
   if isone(x)
      return x
   elseif y == 0
      return one(parent(x))
   elseif y == 1
      return x
   else
      R = parent(x)
      p = GC.@preserve x R libSingular.n_Power(x.ptr, y, R.ptr)
      return R(p)
   end
end

###############################################################################
#
#   Parent call functions
#
###############################################################################

function (R::N_UnknownSingularCoefficientRing)(x::Integer)
   GC.@preserve R return n_unknownsingularcoefficient(R, libSingular.n_InitMPZ(BigInt(x), R.ptr))
end

function (R::N_UnknownSingularCoefficientRing)(n::libSingular.number_ptr)
   return n_unknownsingularcoefficient(R, n)
end

function (R::N_UnknownSingularCoefficientRing)(a::n_unknownsingularcoefficient)
   parent(a) === R || error("Unable to coerce element")
   return a
end

