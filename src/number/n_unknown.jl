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
   n = libSingular.julia(a.ptr)
   return R(libSingular.number(canonical_unit(n)))
end

function deepcopy_internal(a::n_unknown, dict::IdDict)
   R = parent(a)
   n = libSingular.julia(a.ptr)
   return R(libSingular.number(deepcopy(n)))
end

function one(R::CoefficientRing)
   return R(libSingular.n_Init(1, R.ptr))
end

function zero(R::CoefficientRing)
   return R(libSingular.n_Init(0, R.ptr))
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

   libSingular.n_Write(libSingular.cast_void_to_number(a.ptr), parent(a).ptr, false)

   m = libSingular.StringEndS()

   print(io, m)
end

###############################################################################
#
#   Unary operators
#
###############################################################################

function -(a::n_unknown)
   R = parent(a)
   n = libSingular.n_Neg(libSingular.cast_void_to_number(a.ptr), R.ptr)
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
   n = libSingular.n_Add(libSingular.cast_void_to_number(a.ptr), libSingular.cast_void_to_number(b.ptr), R.ptr)
   return R(n)
end

function -(a::n_unknown{T}, b::n_unknown{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   n = libSingular.n_Sub(libSingular.cast_void_to_number(a.ptr), libSingular.cast_void_to_number(b.ptr), R.ptr)
   return R(n)
end

function *(a::n_unknown{T}, b::n_unknown{T}) where T <: Nemo.RingElem
   check_parent(a, b)
   R = parent(a)
   n = libSingular.n_Mult(libSingular.cast_void_to_number(a.ptr), libSingular.cast_void_to_number(b.ptr), R.ptr)
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
   return libSingular.n_Equal(libSingular.cast_void_to_number(a.ptr), libSingular.cast_void_to_number(b.ptr), R.ptr)
end

###############################################################################
#
#   Inversion
#
###############################################################################

function inv(a::n_unknown{T}) where T <: Nemo.FieldElem
   R = parent(a)
   n = libSingular.n_Invers(libSingular.cast_void_to_number(a.ptr), R.ptr)
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
   n = libSingular.n_ExactDiv(libSingular.cast_void_to_number(a.ptr), libSingular.cast_void_to_number(b.ptr), R.ptr)
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
   n = libSingular.n_Gcd(libSingular.cast_void_to_number(a.ptr), libSingular.cast_void_to_number(b.ptr), R.ptr)
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
   n = libSingular.n_ExtGcd(libSingular.cast_void_to_number(a.ptr), libSingular.cast_void_to_number(b.ptr), pointer(s), pointer(t), R.ptr)
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
   return n_unknown(libSingular.number(a), R)
end

function (R::CoefficientRing{T})(a::Integer) where T <: Nemo.RingElem
   return R(base_ring(R)(a))
end

function (R::CoefficientRing{T})(a::n_unknown{T}) where T <: Nemo.RingElem
   base_ring(a) != base_ring(R) && error("Unable to coerce element")
   return a
end

function (R::CoefficientRing{T})(ptr::libSingular.number) where {T <: Nemo.RingElem}
   return n_unknown{T}(libSingular.cast_number_to_void(ptr), R)
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

