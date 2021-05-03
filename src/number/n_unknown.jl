###############################################################################
#
#   Basic manipulation
#
###############################################################################

elem_type(::Type{CoefficientRing{T}}) where T <: Nemo.RingElem = n_unknown{T}

parent_type(::Type{n_unknown{T}}) where T <: Nemo.RingElem = CoefficientRing{T}

parent(a::n_unknown) = a.parent

base_ring(R::CoefficientRing) = R.base_ring

base_ring(a::n_unknown) = base_ring(parent(a))

function check_parent(a::n_unknown, b::n_unknown)
   parent(a) != parent(b) && error("Incompatible parents")
end

function canonical_unit(a::n_unknown)
   R = parent(a)
   n = GC.@preserve a libSingular.julia(libSingular.cast_number_to_void(a.ptr))
   return R(canonical_unit(n))
end

function deepcopy_internal(a::n_unknown, dict::IdDict)
   R = parent(a)
   n = GC.@preserve a libSingular.julia(libSingular.cast_number_to_void(a.ptr))
   return R(deepcopy(n))
end

function one(R::CoefficientRing)
   return R(libSingular.n_Init(1, R.ptr))
end

function zero(R::CoefficientRing)
   return R(libSingular.n_Init(0, R.ptr))
end

function hash(a::n_unknown, h::UInt)
   n = GC.@preserve a libSingular.julia(libSingular.cast_number_to_void(a.ptr))
   return xor(hash(n, h), 0x664e59de562461fe%UInt)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, R::CoefficientRing)
   print(io, "Singular coefficient ring over ")
   show(io, R.base_ring)
end

function show(io::IO, a::n_unknown)
   libSingular.StringSetS("")
   R = parent(a)
   GC.@preserve a R libSingular.n_Write(a.ptr, R.ptr, false)
   print(io, libSingular.StringEndS())
end

function AbstractAlgebra.expressify(a::n_unknown; context = nothing)
   n = GC.@preserve a libSingular.julia(libSingular.cast_number_to_void(a.ptr))
   return AbstractAlgebra.expressify(n; context = context)
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::n_unknown)
   R = parent(a)
   n = GC.@preserve a R libSingular.n_Neg(a.ptr, R.ptr)
   return R(n)
end

###############################################################################
#
#   Binary operators
#
###############################################################################

function +(a::n_unknown{T}, b::n_unknown{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   n = GC.@preserve a b R libSingular.n_Add(a.ptr, b.ptr, R.ptr)
   return R(n)
end

function -(a::n_unknown{T}, b::n_unknown{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   n = GC.@preserve a b R libSingular.n_Sub(a.ptr, b.ptr, R.ptr)
   return R(n)
end

function *(a::n_unknown{T}, b::n_unknown{T}) where T <: Nemo.RingElem
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

function ==(a::n_unknown{T}, b::n_unknown{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   GC.@preserve a b R return libSingular.n_Equal(a.ptr, b.ptr, R.ptr)
end

###############################################################################
#
#   Inversion
#
###############################################################################

function inv(a::n_unknown{T}) where T <: Nemo.FieldElem
   R = parent(a)
   n = GC.@preserve a R libSingular.n_Invers(a.ptr, R.ptr)
   return R(n)
end

###############################################################################
#
#   Exact division
#
###############################################################################

function divexact(a::n_unknown{T}, b::n_unknown{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   n = GC.@preserve a b R libSingular.n_ExactDiv(a.ptr, b.ptr, R.ptr)
   return R(n)
end

###############################################################################
#
#   GCD
#
###############################################################################

function gcd(a::n_unknown{T}, b::n_unknown{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   n = GC.@preserve a b R libSingular.n_Gcd(a.ptr, b.ptr, R.ptr)
   return R(n)
end

###############################################################################
#
#   Extended GCD
#
###############################################################################

function gcdx(a::n_unknown{T}, b::n_unknown{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   s = [libSingular.n_Init(0, R.ptr)]
   t = [libSingular.n_Init(0, R.ptr)]
   n = GC.@preserve a b s t libSingular.n_ExtGcd(a.ptr, b.ptr, pointer(s), pointer(t), R.ptr)
   return R(n), R(s[]), R(t[])
end

###############################################################################
#
#   Conversions and promotions
#
###############################################################################

promote_rule(C::Type{n_unknown{S}}, ::Type{T}) where {S <: Nemo.RingElem, T <: Integer} = n_unknown{S}

promote_rule(::Type{n_unknown{T}}, ::Type{T}) where {T <: Nemo.RingElem} = n_unknown{T}

promote_rule1(::Type{U}, ::Type{n_unknown{T}}) where {T <: Nemo.RingElem, U <: Nemo.RingElem} = promote_rule(U, n_unknown{T})

function promote_rule1(::Type{n_unknown{T}}, ::Type{n_unknown{U}}) where {T <: Nemo.RingElem, U <: Nemo.RingElem}
   promote_rule(T, n_unknown{U}) == T ? n_unknown{T} : Union{}
end

function promote_rule(::Type{n_unknown{T}}, ::Type{U}) where {T <: Nemo.RingElem, U <: Nemo.RingElem}
   promote_rule(T, U) == T ? n_unknown{T} : promote_rule1(U, n_unknown{T})
end

promote_rule(::Type{n_unknown{S}}, ::Type{Nemo.fmpz}) where {S <: Nemo.RingElem} = n_unknown{S}

promote_rule(::Type{n_unknown{S}}, ::Type{Nemo.fmpq}) where {S <: Nemo.FieldElem} = n_unknown{S}

###############################################################################
#
#   Parent call overloads
#
###############################################################################

function (R::CoefficientRing{T})(a::T) where T <: Nemo.RingElem
   return n_unknown(libSingular.cast_void_to_number(libSingular.number(a)), R)
end

function (R::CoefficientRing{T})(a::Integer) where T <: Nemo.RingElem
   return R(base_ring(R)(a))
end

function (R::CoefficientRing{T})(a::Rational) where T <: Nemo.RingElem
   return R(base_ring(R)(a))
end

function (R::CoefficientRing{T})(a::n_Z) where T <: Nemo.RingElem
   return R(base_ring(R)(convert(BigInt, a)))
end

function (R::CoefficientRing{T})(a::n_Q) where T <: Nemo.RingElem
   return R(base_ring(R)(Rational{BigInt}(a)))
end

function (R::CoefficientRing{T})(a::n_unknown{T}) where T <: Nemo.RingElem
   base_ring(a) != base_ring(R) && error("Unable to coerce element")
   return a
end

function (R::CoefficientRing{T})(ptr::libSingular.number_ptr) where {T <: Nemo.RingElem}
   return n_unknown{T}(ptr, R)
end

function (R::CoefficientRing{T})(a::S) where {T <: Nemo.RingElem, S <: Nemo.RingElem}
   U = base_ring(R)
   return R(U(a))
end

###############################################################################
#
#   CoefficientRing constructor
#
###############################################################################

function CoefficientRing(R::Nemo.Ring)
   T = elem_type(R)
   return CoefficientRing{T}(R)
end

